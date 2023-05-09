import sys
import pysam
import argparse
import numpy as np
from tqdm import tqdm

from solution import Solution
from heuristic import guess_orf1ab, predict_ORFs
from MWIS import Interval, Dual_Interval
from annotation import get_annotation_region
from util import (
    get_name,
    make_score_func,
    get_description,
    timing
)

# Returns candidate dict w/ region names as keys for (start, stop) tuples
def get_candidate_region(fasta_file: str, annotation_file: str, window_size: int):
    regions = get_annotation_region(annotation_file)
    fasta = pysam.FastaFile(fasta_file)
    contig = fasta.references[0]
    ref = fasta.fetch(contig)
    _, possible_trs = predict_ORFs(ref)
    candidates, old_start = {}, 0
    for name, (start, end) in regions.items():
        previous_ATG = ref[:start].rfind("ATG")
        previous_2ATG = ref[:previous_ATG].rfind("ATG")
        if possible_trs[previous_ATG] - previous_ATG < 100:
            candidates[name] = (max(old_start, previous_2ATG), start + window_size)
        else:
            candidates[name] = (max(old_start, previous_ATG), start + window_size)
        old_start = start
    return candidates

# Returns possible_orf dict with start_pos keys and end_pos values
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

def semi_smith_waterman(s1: str, s2: str, candidates, window, tau, is_dual, match=1, mismatch=-2, indel=-1):
    s1 = s1.upper()
    s2 = s2.upper()
    len1 = len(s1)
    len2 = len(s2)
    score_f = make_score_func(match, mismatch, indel)
    possible_orf = get_orf_length(s2, [start-window for _, start in candidates])

    # Prepare diagonal scores for i >= window start
    diag = np.zeros((len1+1, len2+1), dtype=int)
    for start, end in candidates:
        for j in range(start, end):
            for i in range(0, len1):
                diag[i, j] = diag[i-1, j-1] + score_f(s1[i], s2[j])

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

    for start, end in candidates:
        origin[:, start, 0] = np.arange(0, len1 + 1)
        origin[:, start, 1] = start
        origin[:, start-1, 0] = np.arange(0, len1 + 1)
        origin[:, start-1, 1] = start-1
        for j in range(start, end):
            for i in range(0, len1):
                new_score = local[i-1, j-1] + score_f(s1[i], s2[j])

                # Produce longest possible alignment of equal score by allowing new_score = 0
                if new_score >= 0:
                    local[i, j] = new_score
                    origin[i, j, :] = origin[i-1, j-1, :]
                else:
                    local[i, j] = 0
                    origin[i, j, :] = [i, j]

    # Dynamic Programming
    score_sweep = []
    intervals_sweep = []
    for window_start in tqdm(range(0, len1-window)):
        max_scores = []
        max_intervals = []
        for start, end in candidates:
            max_score = float("-inf")
            max_intv = None

            # Skip lower left triangular region since it won't contain full window
            for i in range(window_start + window - 1, len1):
                for j in range(start + i + 1 - window_start, end):
                    idx_diag = j - i
                    prev = window_start - 1

                    # delta = the shift of the diagonal line caused by stitching with the local alignment table
                    if window_start == 0:
                        delta = 0
                    else:
                        delta = local[prev, idx_diag + prev] - diag[prev, idx_diag + prev]

                    cur_score = diag[i, j] + delta
                    if cur_score >= tau and cur_score >= max_score:
                        max_score = cur_score
                        x, y = origin[prev, idx_diag + prev, :]
                        orf_start = end - window
                        len_orf = possible_orf[orf_start] - orf_start + 2
                        
                        # idx set to window start
                        max_intv = Interval(window_start, start, end, cur_score, 
                                            x+1, i, y+1, j, cur_score,
                                            len_orf, orf_start, possible_orf[orf_start])
            
            if max_intv is not None:
                max_scores.append(max_score)
                max_intervals.append(max_intv)
            elif is_dual:
                max_scores.append(0)
                max_intervals.append(None)
        
        score_sweep.append(sum(max_scores))
        intervals_sweep.append(max_intervals)

    return score_sweep, intervals_sweep

def dualize(scores, intervals):
    """Find all admissible pairs and convert to Dual_Intervals"""
    dual_intvs, dual_weights = [], []
    for i, intv1 in enumerate(intervals):
        for j, intv2 in enumerate(intervals):
            if i < j and any(intv1) and any(intv2):
                max_core1_e = max(intv.yb for intv in intv1 if intv)
                min_core2_s = min(intv.ya for intv in intv2 if intv)

                # Admissibility Constraint
                if all(a or b for a, b in zip(intv1, intv2)) and max_core1_e < min_core2_s:
                    # The TRS-L further downstream is placed first, per Zirkel et al
                    dual_intvs.append([Dual_Interval(b, a) for a, b in zip(intv1, intv2)])
                    dual_weights.append(scores[i] + scores[j])

    return np.array(dual_weights), dual_intvs,


WINDOW = 7
TAU = 2

@timing
def corid_annot(ref, regions, annotation, name, description, window, tau, mismatch, is_dual):
    leader_end, _, orf1ab_end = guess_orf1ab(ref)
    leader_end = min(leader_end, 500)
    candidates = list((s-leader_end, e-leader_end) for s, e in regions.values())
    candidates = [gene for gene in candidates if gene[0] > 0 and gene[1] > 0]
    leader = ref[:leader_end].replace('N', '-')

    score_sweep, intervals_sweep = semi_smith_waterman(leader,
                                                       ref[leader_end:],
                                                       candidates,
                                                       window=window,
                                                       tau=tau,
                                                       mismatch=mismatch,
                                                       is_dual=is_dual)
    
    if is_dual:
        score_sweep, intervals_sweep = dualize(score_sweep, intervals_sweep)

    solution = Solution(score_sweep, intervals_sweep, leader_end, window, is_corsid_a=True, is_dual=is_dual)
    compact_score = solution.get_compact_score(ref, orf1ab_end, offset=leader_end-orf1ab_end)
    result = solution.serialize_results(ref, window, name, description, annotation, compact=compact_score)

    return result


def main():
    # Sample Commands
    # corsid_a -f test/NC_045512.fasta -g test/NC_045512.gff -o test/NC_045512.corsid_a.json > test/NC_045512.corsid_a.txt
    # corsid -f test/NC_045512.fasta -o test/NC_045512.corsid.json > test/NC_045512.corsid.txt

    parser = argparse.ArgumentParser()

    required = parser.add_argument_group('required arguments')
    required.add_argument("-f", "--fasta", required=True, type=str, help="FASTA genome file")
    required.add_argument("-g", "--gff", required=True, type=str, help="GFF annotation file")
    
    parser.add_argument("-n", "--name", type=str,
                        help="sample name",
                        default=None)
    parser.add_argument("-o", "--output", type=str,
                        help="output file name")
    parser.add_argument("-w", "--window", type=int,
                        help=f"length of sliding window [{WINDOW}]",
                        default=WINDOW)
    parser.add_argument("-x", "--mismatch", type=int,
                        help=f"mismatch score [-2]",
                        default=-2)
    parser.add_argument("-t", "--tau", type=int,
                        help=f"minimum alignment score threshold [{TAU}]",
                        default=TAU)
    parser.add_argument("-d", "--is_dual", type=bool,
                        help="enable dual alignment",
                        default=False)
    args = parser.parse_args(None if sys.argv[1:] else ['-h'])

    print(' '.join(sys.argv))

    fasta = pysam.Fastafile(args.fasta)
    ref = fasta.fetch(fasta.references[0])
    regions = get_candidate_region(args.fasta, args.gff, args.window)
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
                          args.tau,
                          args.mismatch,
                          args.is_dual)
    results.write_result()

    if args.output:
        with open(args.output, "w") as ofile:
            ofile.write(results.to_json())

if __name__ == "__main__":
    main()
