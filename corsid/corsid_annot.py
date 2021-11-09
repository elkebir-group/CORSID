import sys
import pysam
import argparse
import numpy as np
from .util import (
    get_name,
    make_score_func,
    get_description,
    timing,
)
from .solution import Solution
from .heuristic import (
    guess_orf1ab,
    predict_ORFs,
)
from .MWIS import (
    Interval,
)
from .annotation import get_annotation_region
from tqdm import tqdm

def get_candidate_region(fasta_file: str, annotation_file: str):
    regions = get_annotation_region(annotation_file)
    fasta = pysam.FastaFile(fasta_file)
    contig = fasta.references[0]
    ref = fasta.fetch(contig)
    _, possible_trs = predict_ORFs(ref)
    candidates = {}
    old_start = 0
    for name, (start, end) in regions.items():
        previous_ATG = ref[:start].rfind("ATG")
        previous_2ATG = ref[:previous_ATG].rfind("ATG")
        if possible_trs[previous_ATG] - previous_ATG < 100:
            candidates[name] = (max(old_start, previous_2ATG), start)
        else:
            candidates[name] = (max(old_start, previous_ATG), start)
        old_start = start
    return candidates

def get_orf_length(seq, starts):
    stop_codons = ["TAA", "TAG", "TGA"]
    stop_pos = []
    for i in range(len(seq)-3):
        if seq[i:i+3] in stop_codons:
            stop_pos.append(i)
    possible_orf = {}
    for start in starts:
        for stop in stop_pos:
            if stop > start and (stop - start) % 3 == 0:
                possible_orf[start] = stop
                break
        else:
            possible_orf[start] = len(seq)
    return possible_orf

def semi_smith_waterman(s1: str, s2: str, genes, window, match=1, mismatch=-1, indel=-1):
    s1 = s1.upper()
    s2 = s2.upper()
    len1 = len(s1)
    len2 = len(s2)
    score = make_score_func(match, mismatch, indel)
    possible_orf = get_orf_length(s2, [start for _, start in genes])

    # Prepare diagonal scores for i >= window start
    diag = np.zeros((len1+1, len2+1), dtype=int)
    for start, end in genes:
        for j in range(start, end):
            for i in range(0, len1):
                diag[i, j] = diag[i-1, j-1] + score(s1[i], s2[j])

    # Prepare local alignment scores for i < window start
    local = np.zeros((len1+1, len2+1), dtype=int)
    origin = np.zeros((len1+1, len2+1, 2), dtype=int)
    origin[0, :, 0] = 0
    origin[0, :, 1] = np.arange(0, len2 + 1)
    origin[-1, :, 0] = 0
    origin[-1, :, 1] = np.arange(1, len2 + 2)
    origin[:, 0, 0] = np.arange(0, len1 + 1)
    origin[:, 0, 1] = 0
    origin[:, -1, 0] = np.arange(1, len1 + 2)
    origin[:, -1, 1] = 0
    for start, end in genes:
        origin[:, start, 0] = np.arange(0, len1 + 1)
        origin[:, start, 1] = start
        origin[:, start-1, 0] = np.arange(0, len1 + 1)
        origin[:, start-1, 1] = start-1
        for j in range(start, end):
            for i in range(0, len1):
                new_score = local[i-1, j-1] + score(s1[i], s2[j])
                if new_score > 0:
                    local[i, j] = new_score
                    origin[i, j, :] = origin[i-1, j-1, :]
                else:
                    local[i, j] = 0
                    origin[i, j, :] = [i, j]

    # dynamic programming
    score_sweep = []
    intervals_sweep = []
    for window_start in tqdm(range(0, len1-window)):
        max_scores = []
        max_intervals = []
        for idx, (start, end) in enumerate(genes):
            max_score = float("-inf")
            max_intv = None
            for i in range(window_start + window - 1, len1):
                # Skip lower left triangular region,
                # since it won't contain full window
                for j in range(start+window+i-window_start-window+1, end):
                    idx_diag = j - i
                    prev = window_start - 1
                    # delta = the shift of the diagonal line caused by stitching with
                    # the local alignment table
                    if window_start == 0:
                        delta = 0
                    else:
                        delta = local[prev, idx_diag + prev] - diag[prev, idx_diag + prev]
                    cur_score = diag[i, j] + delta
                    if cur_score >= 2 and cur_score >= max_score:
                        max_score = cur_score
                        len_orf = possible_orf[end] - end
                        x, y = origin[prev, idx_diag + prev, :]
                        max_intv = Interval(0, y, end, cur_score,
                                            x+1, i, y+1, j, cur_score, len_orf, end, possible_orf[end])
            if max_intv is not None:
                max_scores.append(max_score)
                max_intervals.append(max_intv)
        
        score_sweep.append(sum(max_scores))
        intervals_sweep.append(max_intervals[::-1])
    return score_sweep, intervals_sweep


WINDOW = 7
SCORE = 2

@timing
def corid_annot(ref, regions, annotation, name, description, window, mismatch):
    leader_end, _, orf1ab_end = guess_orf1ab(ref)
    leader_end = min(leader_end, 500)
    genes = list((s-leader_end, e-leader_end) for s, e in regions.values())
    genes = [gene for gene in genes if gene[0] > 0 and gene[1] > 0]
    leader = ref[:leader_end].replace('N', '-')

    score_sweep, intervals_sweep = semi_smith_waterman(leader,
                                                       ref[leader_end:],
                                                       genes,
                                                       window=window,
                                                       mismatch=mismatch)

    solution = Solution(score_sweep, intervals_sweep, leader_end)
    compact_score = solution.get_compact_score(ref, orf1ab_end, offset=leader_end-orf1ab_end)
    result = solution.serialize_results(ref, window, name, description, annotation, compact=compact_score, is_corsid_a=True)

    return result


def main():
    parser = argparse.ArgumentParser()

    required = parser.add_argument_group('required arguments')
    required.add_argument("-f", "--fasta", required=True, type=str, help="FASTA genome file")
    required.add_argument("-g", "--gff", required=True, type=str, help="GFF annotation file")
    
    parser.add_argument("-n", "--name", type=str,
                        help="sample name",
                        default=None)
    parser.add_argument("-o", "--output", type=str, help="output file name")
    parser.add_argument("-w", "--window", type=int,
                        help=f"length of sliding window [{WINDOW}]",
                        default=WINDOW)
    parser.add_argument("-x", "--mismatch", type=int,
                        help=f"mismatch score [-2]",
                        default=-2)
    parser.add_argument("-s", "--score", type=int,
                        help=f"minimum alignment score threshold [{SCORE}]",
                        default=SCORE)
    args = parser.parse_args(None if sys.argv[1:] else ['-h'])

    print(' '.join(sys.argv))

    fasta = pysam.Fastafile(args.fasta)
    ref = fasta.fetch(fasta.references[0])
    regions = get_candidate_region(args.fasta, args.gff)
    annotation = get_annotation_region(args.gff)
    description = get_description(args.fasta)

    if args.name is None:
        name = get_name(args.fasta)
    else:
        name = args.name

    results = corid_annot(ref,
                            regions,
                            annotation,
                            name,
                            description,
                            args.window,
                            args.mismatch)
    results.write_result()

    if args.output:
        with open(args.output, "w") as ofile:
            ofile.write(results.to_json())

if __name__ == "__main__":
    main()
