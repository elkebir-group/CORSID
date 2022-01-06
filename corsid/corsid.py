import sys
import pysam
import argparse
import numpy as np
from .util import (
    make_score_func,
    get_description,
    get_name,
)
from .solution import Solution
from .heuristic import (
    guess_orf1ab,
    predict_ORFs,
)
from .MWIS import (
    Interval,
    MWIS,
)
from .annotation import get_annotation_region
from typing import List, Dict, Tuple
from tqdm import tqdm

# Default values
WINDOW = 7
TAU_MAX = 7
TAU_MIN = 2
SHRINK = 0.05


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
                                                    l.len_x_overlap_orf(0, 200) < l.orf_len and
                                                    l.len_x_overlap_orf(intv.orf_start, intv.orf_end) < l.orf_len)]
    for i in range(len(rights)):
        for intv in partial:
            if intv.xa <= i <= intv.xb:
                rights[i] = []
            else:
                rights[i] = [r for r in rights[i] if (r.len_x_overlap(intv.xa, intv.xb) == 0 and
                                                      r.len_x_overlap_orf(0, 200) < r.orf_len and
                                                      r.len_x_overlap_orf(intv.orf_start, intv.orf_end) < r.orf_len)]


def semi_smith_waterman(s1: str,
                        s2: str,
                        window: int,
                        match=1,
                        mismatch=-1,
                        tau_min=TAU_MIN,
                        tau_max=TAU_MAX,
                        orf_thr=100,
                        shrink=0.05):
    """Modified Smith-Waterman

    The DP table for each window start is stitched together by the local alignment table and the 
    diagonal table at the row `window_start`.

    Args:
        s1 (str): leader sequence
        s2 (str): body sequence
        window (int): window length
        match (int, optional): matching score. Defaults to 1.
        mismatch (int, optional): mismatch score. Defaults to -1.
        tau_min (int, optional): minimum tau (min score). Defaults to TAU_MIN.
        tau_max (int, optional): maximum tau (min score). Defaults to TAU_MAX.
        orf_thr (int, optional): ORF length threshold. Defaults to 100.
        shrink (float, optional): ORF shrinking percentage. Defaults to 0.05.

    Returns:
        Tuple: maximum scores, and intervals
    """
    s1 = s1.upper()
    s2 = s2.upper()
    len1 = len(s1)
    len2 = len(s2)
    score_f = make_score_func(match, mismatch)
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
            if new_score > 0:
                local[i, j] = new_score
                origin[i, j, :] = origin[i - 1, j - 1, :]
            else:
                local[i, j] = 0
                origin[i, j, :] = [i, j]

    # MWIS
    score_sweep = {i: [0] * (len1 - window) for i in range(tau_min, tau_max+1)}
    intervals_sweep = {i: [None] * (len1 - window) for i in range(tau_min, tau_max+1)}
    for window_start in tqdm(range(0, len1 - window)):
        d_intervals = {}
        # Iterate along diagonal lines
        for idx_diag in range(1 - window_start, len2 - window_start - window + 1):
            prev = window_start - 1
            # delta = the shift of the diagonal line caused by stitching with
            # the local alignment table
            if window_start == 0:
                delta = 0
            else:
                delta = local[prev, idx_diag + prev] - diag[prev, idx_diag + prev]
            next_orf_start = next_start[idx_diag + window_start + window - 3]
            if next_orf_start is None:
                continue
            next_orf_end = possible_trs[next_orf_start]
            len_orf = next_orf_end - next_orf_start
            
            # try once more if length is not enough
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
                j = i + idx_diag
                # Skip lower left triangular region,
                # since it won't contain full window
                if j >= len2:
                    break
                cur_score = diag[i, j] + delta
                if (
                    (cur_score >= 2) and
                    next_start[j] is not None
                ):
                    weight = float(len_orf) + cur_score / 1000
                    if len_orf > orf_thr:
                        y = origin[prev, idx_diag + prev, 0]
                        x = origin[prev, idx_diag + prev, 1]
                        shrink_start = int(next_orf_start + shrink * len_orf)
                        shrink_end = int(next_orf_start + (1-shrink) * len_orf)
                        if (shrink_start, shrink_end) not in d_intervals:
                            d_intervals[shrink_start, shrink_end] = Interval(0, shrink_start, shrink_end, weight,
                                        y + 1, i, x + 1, j, cur_score, len_orf, next_orf_start, next_orf_end)
                        else:
                            old_intv = d_intervals[shrink_start, shrink_end]
                            if cur_score > old_intv.score:
                                d_intervals[shrink_start, shrink_end] = Interval(0, shrink_start, shrink_end, weight,
                                        y + 1, i, x + 1, j, cur_score, len_orf, next_orf_start, next_orf_end)

        lefts = [[] for _ in range(len2)]
        rights = [[] for _ in range(len2)]
        # print(f"m={len(d_intervals)}", end=' ')
        for (l, r), intv in d_intervals.items():
            lefts[l].append(intv)
            rights[r].append(intv)
        idx_intv = 0
        for right in rights:
            for r in right:
                r.idx = idx_intv
                idx_intv += 1
        score, intervals = MWIS(lefts, rights, thr=tau_max)
        score_sweep[tau_max][window_start] = score
        intervals_sweep[tau_max][window_start] = intervals

        if len(intervals) == 0:
            continue
        remove_partial(lefts, rights, intervals)
        for thr in range(tau_max - 1, tau_min-1, -1):
            score, intervals = MWIS(lefts, rights, thr=thr)
            if len(intervals) > 0:
                intervals_sweep[thr][window_start] = merge_solutions(intervals_sweep[thr+1][window_start], intervals)
                score_sweep[thr][window_start] = sum(i.w for i in intervals_sweep[thr][window_start])
                remove_partial(lefts, rights, intervals_sweep[thr][window_start])
            else:
                intervals_sweep[thr][window_start] = intervals_sweep[thr+1][window_start]
                score_sweep[thr][window_start] = score_sweep[thr+1][window_start]
    return score_sweep, intervals_sweep


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
    i = 0
    j = 0
    sol = []
    while True:
        if sol2[j].xa <= sol1[i].xa:
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


def corsid(ref: str,
           annotation: Dict[str, Tuple],
           description: str,
           name: str,
           window: int,
           tau_min: int=TAU_MIN,
           tau_max: int=TAU_MAX,
           mismatch: int=-1,
           shrink: float=0.05):
    """Entrypoint of CORSID

    Args:
        ref (str): reference genome
        annotation (Dict[str, Tuple]): genes in the annotation file
        description ([type]): description string in the annotation file
        name ([type]): name of the genome
        window ([type]): sweeping window length
        tau_min (int, optional): minimum tau (min score). Defaults to TAU_MIN.
        tau_max (int, optional): maximum tau (min score). Defaults to TAU_MAX.
        mismatch (int, optional): mismatch score. Defaults to -1.
        shrink (float, optional): fraction of positions that may overlap between consecutive genes. Defaults to 0.05.

    Returns:
        Tuple[Results, List[float]]: Results and compact scores for each result.
    """
    leader_end, _, orf1ab_end = guess_orf1ab(ref)
    leader_end = min(leader_end, 500)
    offset = orf1ab_end - 200
    leader = ref[:leader_end].replace('N', '-')
    max_weight_sweep, intervals_sweep = semi_smith_waterman(leader, ref[offset:],
                                                            window, mismatch=mismatch, 
                                                            tau_min=tau_min, tau_max=tau_max,
                                                            orf_thr=100, shrink=shrink)
    solution1 = {s: Solution(max_weight_sweep[s], intervals_sweep[s], offset) for s in max_weight_sweep}
    compact_score = solution1[tau_min].get_compact_score(ref, solution1[tau_min].offset)
    result1 = solution1[tau_min].serialize_results(ref,
                                          window,
                                          name,
                                          description,
                                          annotation,
                                          is_lex=True,
                                          compact=compact_score)

    return result1, compact_score


def main():
    parser = argparse.ArgumentParser()

    required = parser.add_argument_group('required arguments')
    required.add_argument("-f", "--fasta", required=True,
                          type=str, help="FASTA genome file")

    parser.add_argument("-g", "--gff", type=str, help="GFF annotation file")
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
    args = parser.parse_args(None if sys.argv[1:] else ['-h'])

    print(' '.join(sys.argv))

    fasta = pysam.Fastafile(args.fasta)
    ref = fasta.fetch(fasta.references[0])
    if args.gff:
        annotation = get_annotation_region(args.gff)
    else:
        annotation = None
    description = get_description(args.fasta)

    if args.name is None:
        name = get_name(args.fasta)
    else:
        name = args.name

    result, compact_score = corsid(ref,
                                       annotation,
                                       description,
                                       name,
                                       args.window,
                                       args.tau_min,
                                       args.tau_max,
                                       args.mismatch,
                                       args.shrink)

    # Fall back to smaller window if not compact enough
    opt_pos = result.results[0].leader_core_start
    if compact_score[opt_pos] >= 0.2:
        result, _ = corsid(ref,
                               annotation,
                               description,
                               name,
                               args.window - 1,
                               args.tau_min,
                               args.tau_max,
                               args.mismatch,
                               args.shrink)

    result.write_result()

    if args.output:
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


if __name__ == "__main__":
    main()
