import sys
import pysam
import argparse
import numpy as np
from typing import List, Dict, Tuple
from tqdm import tqdm
import multiprocessing
import gzip
import xgboost as xgb
import pathlib
import warnings

from solution import Solution
from MWIS import Interval, Dual_Interval, MWIS
from heuristic import guess_orf1ab, predict_ORFs
from annotation import get_annotation_region
from util import make_score_func, get_description, get_name

# Default Parameters
WINDOW = 7
TAU_MAX = 7
TAU_MIN = 2
SHRINK = 0.05


def warning_on_one_line(message, category, filename, lineno, file=None, line=None):
    return f'{filename}:{lineno}: {category.__name__}:{message}\n'

warnings.formatwarning = warning_on_one_line

def information_content(seqs):
    n = len(seqs)
    seqs = [seq.upper() for seq in seqs]
    cols = [''.join(col) for col in zip(*seqs)]
    counts = [[col.count('A'), col.count('T'), col.count('C'), col.count('G'), col.count('-')] for col in cols]
    return [2+sum([c/n * np.log2(c/n) for c in count if c > 0]) for count in counts]


def remove_partial(lefts: List[List[Interval]],
                   rights: List[List[Interval]],
                   partial: List[Interval]):
    """Remove intervals overlap with partial solution from lefts and rights

    Args:
        lefts (List[List[Interval]]): intervals indexed by left endpoints
        rights (List[List[Interval]]): intervals indexed by right endpoints
        partial (List[Interval]): partial solution
    """
    for i in range(len(lefts)):
        for intv in partial:
            if intv.xa <= i <= intv.xb:
                lefts[i] = []
            else:
                lefts[i] = [l for l in lefts[i] if (l.len_x_overlap(intv.xa, intv.xb) == 0 and
                                                    l.len_x_overlap_orf(0, 200) < l.orf_length and
                                                    l.len_x_overlap_orf(intv.orf_start, intv.orf_end) < l.orf_length)]
    for i in range(len(rights)):
        for intv in partial:
            if intv.xa <= i <= intv.xb:
                rights[i] = []
            else:
                rights[i] = [r for r in rights[i] if (r.len_x_overlap(intv.xa, intv.xb) == 0 and
                                                      r.len_x_overlap_orf(0, 200) < r.orf_length and
                                                      r.len_x_overlap_orf(intv.orf_start, intv.orf_end) < r.orf_length)]
                

def merge_solutions(sol1: List[Interval], sol2: List[Interval]):
    """Merge 2 solutions and sort them

    Args:
        sol1 (List[Interval]): Solution 1
        sol2 (List[Interval]): Solution 2

    Returns:
        List[Interval]: Merged sorted solution
    """
    assert len(sol1) > 0
    assert len(sol2) > 0
    i, j = 0, 0
    sol = []
    while True:
        if sol2[j].xa <= sol1[i].xa: # Should this be flipped?
            sol.append(sol1[i])
            i += 1
        else:
            sol.append(sol2[j])
            j += 1
        if i >= len(sol1):
            sol.extend(sol2[j:])
            break
        if j >= len(sol2):
            sol.extend(sol1[i:])
            break
    return sol


def semi_smith_waterman(s1: str,
                        s2: str,
                        window: int,
                        match=1,
                        mismatch=-2,
                        indel=-1,
                        orf_thr=100,
                        shrink=0.05,
                        is_dual=False):
    """Modified Smith-Waterman

    The DP table for each window start is stitched together by the local alignment table and the 
    diagonal table at the row `window_start`.

    Args:
        s1 (str): leader sequence
        s2 (str): body sequence
        window (int): window length
        match (int, optional): matching score. Defaults to 1.
        mismatch (int, optional): mismatch score. Defaults to -2.
        indel (int, optional): indel score. Defaults to -1.
        orf_thr (int, optional): ORF length threshold. Defaults to 100.
        shrink (float, optional): ORF shrinking percentage. Defaults to 0.05.
        is_dual (bool, optional): solve for optimal pair of TRS alignments. Defaults to False.

    Returns:
        Tuple: maximum scores, and intervals
    """
    s1 = s1.upper()
    s2 = s2.upper()
    len1 = len(s1)
    len2 = len(s2)
    score_f = make_score_func(match, mismatch, indel)
    next_start, possible_trs = predict_ORFs(s2)

    # Prepare diagonal scores for i >= window start
    diag = np.zeros((len1 + 1, len2 + 1), dtype=int)
    for j in range(len2):
        for i in range(len1):
            diag[i, j] = diag[i - 1, j - 1] + score_f(s1[i], s2[j])

    # Prepare local alignment scores for i < window start
    local = np.zeros((len1 + 1, len2 + 1), dtype=int)
    origin = np.zeros((len1 + 1, len2 + 1, 2), dtype=int)
    origin[0, :, 0] = 0
    origin[0, :, 1] = np.arange(0, len2 + 1)
    origin[-1, :, 0] = 0
    origin[-1, :, 1] = np.arange(1, len2 + 2)
    origin[:, 0, 0] = np.arange(0, len1 + 1)
    origin[:, 0, 1] = 0
    origin[:, -1, 0] = np.arange(1, len1 + 2)
    origin[:, -1, 1] = 0

    for j in range(len2):
        for i in range(len1):
            new_score = local[i - 1, j - 1] + score_f(s1[i], s2[j])

            # Produce longest possible alignment of equal score by allowing new_score = 0
            if new_score >= 0:
                local[i, j] = new_score
                origin[i, j, :] = origin[i - 1, j - 1, :]
            else:
                local[i, j] = 0
                origin[i, j, :] = [i, j]

    # Dynamic Programming
    d_intervals = [{} for _ in range(len1 - window)]
    for window_start in tqdm(range(0, len1 - window)):

        # Iterate along diagonal lines
        for idx_diag in range(1 - window_start, len2 - window_start - window + 1):
            prev = window_start - 1

            # delta = the shift of the diagonal line caused by stitching with the local alignment table
            if window_start == 0:
                delta = 0
            else:
                delta = local[prev, idx_diag + prev] - diag[prev, idx_diag + prev]

            # Get closest downstream ORF
            next_orf_start = next_start[idx_diag + window_start + 1]
            if next_orf_start is None:
                continue
            next_orf_end = possible_trs[next_orf_start]
            len_orf = next_orf_end - next_orf_start
            
            # Try once more if length is not enough
            if len_orf < orf_thr:
                try:
                    next_orf_start = next_start[next_orf_start+3]
                except IndexError:
                    continue
                if next_orf_start is None:
                    continue
                next_orf_end = possible_trs[next_orf_start]
                len_orf = next_orf_end - next_orf_start
                if len_orf < orf_thr:
                    continue

            for i in range(window_start + window - 1, len1):
                # Skip lower left triangular region since it won't contain full window
                j = i + idx_diag
                if j >= len2:
                    break
                
                cur_score = diag[i, j] + delta
                if cur_score >= TAU_MIN and next_start[j] is not None:
                    shrink_start = int(next_orf_start + shrink * len_orf)
                    shrink_end = int(next_orf_start + (1-shrink) * len_orf)
                    weight = float(len_orf) + cur_score / 1000
                    x, y = origin[prev, idx_diag + prev, :]

                    if (shrink_start, shrink_end) not in d_intervals[window_start]:
                        d_intervals[window_start][shrink_start, shrink_end] = \
                                Interval(window_start, shrink_start, shrink_end, weight,
                                         x + 1, i, y + 1, j, cur_score,
                                         len_orf, next_orf_start, next_orf_end)
                    else:
                        old_intv = d_intervals[window_start][shrink_start, shrink_end]
                        if cur_score > old_intv.score_total:
                            d_intervals[window_start][shrink_start, shrink_end] = \
                                Interval(window_start, shrink_start, shrink_end, weight,
                                         x + 1, i, y + 1, j, cur_score,
                                         len_orf, next_orf_start, next_orf_end)

    # Solve MWIS
    score_sweep = {i: [] for i in range(TAU_MIN, TAU_MAX+1)}
    intervals_sweep = {i: [] for i in range(TAU_MIN, TAU_MAX+1)}

    if is_dual:
        items = [((d_intervals[i], d_intervals[j]), len2)
                for i in range(0, len1 - window) for j in range(0, len1 - window)
                if max(intv.yb for intv in filter(None, d_intervals[i].values())) < min(intv.ya for intv in filter(None, d_intervals[j].values()))]
    else:
        items = [(d_intv, len2) for d_intv in d_intervals]

    with multiprocessing.Pool() as pool:
        results = pool.starmap(solve_MWIS, tqdm(items))

    for scores, intvs in filter(None, results): 
        for k in scores.keys():
            score_sweep[k].append(scores[k])
            intervals_sweep[k].append(intvs[k])

    return score_sweep, intervals_sweep


def solve_MWIS(d_intv, len2):
    # Dualize if more than one intv provided
    if isinstance(d_intv, tuple):
        (d_intv1, d_intv2), d_intv = d_intv, {}
        for (l, r), intv in d_intv1.items():
            if (l, r) in d_intv2:
                d_intv[l, r] = Dual_Interval(d_intv2[l, r], intv)
            else:
                d_intv[l, r] = Dual_Interval(None, intv)
        for (l, r), intv in d_intv2.items():
            if (l, r) not in d_intv1:
                d_intv[l, r] = Dual_Interval(intv, None)

    lefts = [[] for _ in range(len2)]
    rights = [[] for _ in range(len2)]
    for (l, r), intv in d_intv.items():
        lefts[l].append(intv)
        rights[r].append(intv)
    idx_intv = 0
    for right in rights:
        for r in right:
            r.idx = idx_intv
            idx_intv += 1

    score_sweep = {i: 0 for i in range(TAU_MIN, TAU_MAX+1)}
    intervals_sweep = {i: None for i in range(TAU_MIN, TAU_MAX+1)}
    score, intervals = MWIS(lefts, rights, thr=TAU_MAX)
    score_sweep[TAU_MAX] = score
    intervals_sweep[TAU_MAX] = intervals[::-1]

    if len(intervals) == 0:
        return None
    remove_partial(lefts, rights, intervals)

    for thr in range(TAU_MAX-1, TAU_MIN-1, -1):
        score, intervals = MWIS(lefts, rights, thr=thr)
        if len(intervals) > 0:
            intervals_sweep[thr] = merge_solutions(intervals_sweep[thr + 1], intervals)
            score_sweep[thr] = sum(intv.w for intv in intervals_sweep[thr])
            remove_partial(lefts, rights, intervals_sweep[thr])
        else:
            intervals_sweep[thr] = intervals_sweep[thr + 1]
            score_sweep[thr] = score_sweep[thr + 1]

    return score_sweep, intervals_sweep


def corsid(ref: str,
           annotation: Dict[str, Tuple],
           description: str,
           name: str,
           window: int,
           mismatch: int=-1,
           shrink: float=0.05,
           is_dual: bool=False):
    """Entrypoint of CORSID

    Args:
        ref (str): reference genome
        annotation (Dict[str, Tuple]): genes in the annotation file
        description ([type]): description string in the annotation file
        name ([type]): name of the genome
        window ([type]): sweeping window length
        mismatch (int, optional): mismatch score. Defaults to -1.
        shrink (float, optional): fraction of positions that may overlap between consecutive genes. Defaults to 0.05.
        is_dual (bool, optional): solve for optimal pair of TRS alignments. Defaults to False.

    Returns:
        Tuple[Results, List[float]]: Results and compact scores for each result.
    """
    leader_end, _, orf1ab_end = guess_orf1ab(ref)
    leader_end = min(leader_end, 500)
    offset = orf1ab_end - 200
    leader = ref[:leader_end].replace('N', '-')
    max_weight_sweep, intervals_sweep = semi_smith_waterman(leader,
                                                            ref[offset:],
                                                            window,
                                                            mismatch=mismatch, 
                                                            orf_thr=100,
                                                            shrink=shrink,
                                                            is_dual=is_dual)
    
    solution = {s: Solution(max_weight_sweep[s], intervals_sweep[s], offset, window, is_corsid_a=False, is_dual=is_dual) for s in max_weight_sweep}
    compact_score = solution[TAU_MIN].get_compact_score(ref, solution[TAU_MIN].offset)
    results = solution[TAU_MIN].serialize_results(ref,
                                                  window,
                                                  name,
                                                  description,
                                                  annotation,
                                                  is_lex=True,
                                                  compact=compact_score)

    return results, compact_score


def main():
    parser = argparse.ArgumentParser()

    required = parser.add_argument_group('required arguments')
    required.add_argument("-f", "--fasta", required=True, type=str,
                          help="FASTA genome file")
    
    parser.add_argument("-g", "--gff", type=str,
                        help="GFF annotation file")
    parser.add_argument("-n", "--name", type=str,
                        help="sample name",
                        default=None)
    parser.add_argument("-o", "--output", type=str, help="output json file name")
    parser.add_argument("-r", "--output-orf", type=str,
                        help="output identified ORFs (FASTA), only contains the first solution")
    parser.add_argument("-3", "--output-gff3", type=str,
                        help="output identified ORFs (FASTA), only contains the first solution")
    parser.add_argument("-w", "--window", type=int,
                        help=f"length of sliding window [{WINDOW}]",
                        default=WINDOW)
    parser.add_argument("-x", "--mismatch", type=int,
                        help=f"mismatch score [-2]",
                        default=-2)
    parser.add_argument("-t", "--tau_min", type=int,
                        help=f"minimum matching score threshold [{TAU_MIN}]",
                        default=TAU_MIN)
    parser.add_argument("-T", "--tau_max", type=int,
                        help=f"maximum matching score threshold [{TAU_MAX}]",
                        default=TAU_MAX)
    parser.add_argument("--shrink", type=float,
                        help=f"fraction of positions that may overlap between consecutive genes [{SHRINK}]",
                        default=SHRINK)
    parser.add_argument("--no-missing-classifier", action='store_true',
                        help="set flag to disable missing TRS-L classifier")
    parser.add_argument("-d", "--is_dual", type=bool,
                        help="enable dual alignment",
                        default=False)
    
    args = parser.parse_args(None if sys.argv[1:] else ['-h'])
    print(' '.join(sys.argv))

    fasta = pysam.Fastafile(args.fasta)
    ref = fasta.fetch(fasta.references[0])
    description = get_description(args.fasta)

    if args.name is None:
        name = get_name(args.fasta)
    else:
        name = args.name

    if args.gff:
        annotation = get_annotation_region(args.gff)
    else:
        annotation = None

    result, compact_score = corsid(ref,
                                   annotation,
                                   description,
                                   name,
                                   args.window,
                                   args.mismatch,
                                   args.shrink,
                                   args.is_dual)

    if len(result.results) == 0:
        warnings.warn("CORSID finds no solution. The input genome is possibly incomplete and missing a TRS-L.")
        return -1

    # # Fall back to smaller window if not compact enough
    # opt_pos = result.results[0].leader_core_start
    # if compact_score[opt_pos] >= 0.2:
    #     result, _ = corsid(ref,
    #                        annotation,
    #                        description,
    #                        name,
    #                        args.window - 1,
    #                        args.tau_min,
    #                        args.tau_max,
    #                        args.mismatch,
    #                        args.shrink)

    # if not args.no_missing_classifier:
    #     clf = xgb.XGBClassifier(n_jobs=1)
    #     model = pathlib.Path(__file__).parent.resolve() / "xgboost_model.json"
    #     clf.load_model(model)
    #     # leader, score, dist, information, mean score, compact, ORF1ab start
    #     scores = [int(x.score) for x in result.results[0].bodys]
    #     start = result.results[0].leader_core_start - result.results[0].TRS_L_start
    #     end = start + result.results[0].leader_core_len
    #     seqs = [x.align[start:end] for x in result.results[0].bodys]
    #     features = [[
    #         result.results[0].leader_core_start,
    #         sum(scores),
    #         result.ORF1ab[0] - result.results[0].leader_core_start,
    #         np.mean(information_content(seqs)),
    #         np.mean(scores),
    #         1 - result.results[0].compact,
    #         result.ORF1ab[0]
    #     ]]
    #     decision = clf.predict(features)
    #     if decision[0] == 1:
    #         warnings.warn("The input genome is possibly incomplete and missing a TRS-L.")

    result.write_result()

    if args.output:
        if args.output.endswith(".gz") or args.output.endswith(".gzip"):
            with gzip.open(args.output, "wt") as ofile:
                ofile.write(result.to_json())
        else:
            with open(args.output, "w") as ofile:
                ofile.write(result.to_json())
        # write GFF3 file
        gff_name = '.'.join(args.output.split('.')[:-1] + ["gff"])
        with open(gff_name, "w") as ofile:
            ofile.write(result.to_gff())

    if args.output_orf:
        with open(args.output_orf, "w") as ofile:
            content = []
            ORFs = result.results[0].bodys
            for i, body in enumerate(ORFs):
                if body.ORF:
                    content.append(f">{body.ORF}")
                else:
                    content.append(f">putative_ORF_{i+1}")
                content.append(ref[body.ORF_start:body.ORF_start + body.ORF_len + 3])
            ofile.write('\n'.join(content))

    return 0


if __name__ == "__main__":
    sys.exit(main())
