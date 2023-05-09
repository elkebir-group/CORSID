from typing import List, Tuple
from dataclasses import dataclass


@dataclass
class Interval:
    """Closed interval [a, b] on both x and y."""
    idx: int
    xa: int
    xb: int
    w: float
    ya: int
    yb: int
    start: int
    end: int
    score_total: int
    orf_length: int
    orf_start: int
    orf_end: int

    def __post_init__(self):
        self.orf_len = self.orf_end - self.orf_start + 1

    def __str__(self):
        return f"{self.idx}\t{self.xa}-{self.xb}\t{self.w}\t{self.ya}-{self.yb}\t{self.start}-{self.end}\t{self.orf_start}-{self.orf_end}\t{self.score_total}\t{self.orf_length}"

    def __repr__(self):
        return f"<{self.idx}|xa={self.xa}, xb={self.xb}, w={self.w}, ya={self.ya}, yb={self.yb}, start={self.start}, end={self.end}, start={self.orf_start}, end={self.orf_end}, score={self.score_total}, len={self.orf_length}>"

    def __contains__(self, item):
        return self.orf_start <= item.orf_start and item.orf_end <= self.orf_end

    def len_x_overlap(self, start, end):
        """[x_a, x_b] overlap with [start, end]"""
        if start >= end:
            return 0
        elif start > self.xb or self.xa > end:
            return 0
        else:
            return min(self.xb, end) - max(self.xa, start) + 1

    def len_x_overlap_orf(self, start, end):
        """[orf_start, orf_end] overlap with [start, end]"""
        if start >= end:
            return 0
        elif start > self.orf_end or self.orf_start > end:
            return 0
        else:
            return min(self.orf_end, end) - max(self.orf_start, start) + 1


class Dual_Interval(Interval):
    """Combined closed intervals [a, b] on both x and y."""
    idx1: int
    idx2: int
    xa: int
    xb: int
    w: float
    ya1: int
    yb1: int
    ya2: int
    yb2: int
    start1: int
    end1: int
    start2: int
    end2: int
    score1: int
    score2: int
    score_total: int
    orf_length: int
    orf_start: int
    orf_end: int

    def __init__(self, intv1, intv2):
        # Shared Attributes
        if intv1:
            self.xa = intv1.xa
            self.xb = intv1.xb
            self.orf_length = intv1.orf_length
            self.orf_start = intv1.orf_start
            self.orf_end = intv1.orf_end
        else:
            self.xa = intv2.xa
            self.xb = intv2.xb
            self.orf_length = intv2.orf_length
            self.orf_start = intv2.orf_start
            self.orf_end = intv2.orf_end

        # Individual Attributes
        if intv1:
            self.idx1 = intv1.idx
            self.ya1 = intv1.ya
            self.yb1 = intv1.yb
            self.start1 = intv1.start
            self.end1 = intv1.end
            self.score1 = intv1.score_total
        else:
            self.idx1 = None
            self.ya1 = None
            self.yb1 = None
            self.start1 = None
            self.end1 = None
            self.score1 = None

        if intv2:
            self.idx2 = intv2.idx
            self.ya2 = intv2.ya
            self.yb2 = intv2.yb
            self.start2 = intv2.start
            self.end2 = intv2.end
            self.score2 = intv2.score_total
        else:
            self.idx2 = None
            self.ya2 = None
            self.yb2 = None
            self.start2 = None
            self.end2 = None
            self.score2 = None

        # Combined Score
        if intv1 and intv2:
            self.w = intv1.w + intv2.w
            self.score_total = intv1.score_total + intv2.score_total
        elif intv1:
            self.w = intv1.w
            self.score_total = intv1.score_total
        else:
            self.w = intv2.w
            self.score_total = intv2.score_total

    def __str__(self):
        return f"{self.idx1}\t{self.idx2}\t{self.xa}-{self.xb}\t{self.w}\t{self.ya1}-{self.yb1}\t{self.ya2}-{self.yb2}\t{self.start1}-{self.end1}\t{self.start2}-{self.end2}\t{self.orf_start}-{self.orf_end}\t{self.score1}\t{self.score1}\t{self.orf_length}"

    def __repr__(self):
        return f"<{self.idx1}, {self.idx2}|xa={self.xa}, xb={self.xb}, w={self.w}, ya1={self.ya1}, yb1={self.yb1}, ya2={self.ya2}, yb2={self.yb2}, start1={self.start1}, end1={self.end1}, start2={self.start2}, end2={self.end2}, start={self.orf_start}, end={self.orf_end}, score={self.score1}, score={self.score2}, len={self.orf_length}>"


def MWIS(lefts: List[List[Interval]], rights: List[List[Interval]], thr: float=0) -> Tuple[float, List[Interval]]:
    """Maximum weight independent set.

    Intervals should be sorted by their right endpoints.

    Args:
        lefts (List[List[Interval]]): intervals indexed by left endpoints
        rights (List[List[Interval]]): intervals indexed by right endpoints
        thr (int, optional): minimum score threshold. Defaults to 0.

    Returns:
        Tuple[float, List[Interval]]: max weight and the optimal set.
    """
    n_left = sum(len(l) for l in lefts)
    n_right = sum(len(r) for r in rights)
    assert n_left == n_right
    if n_left == 0:
        return 0, []
    
    n_intv = max([max(l.idx for l in left) for left in lefts if len(left) > 0]) + 1
    sofar = [-1] * n_intv
    map_intervals = {}
    max_sofar = 0
    max_interval = None
    for j in range(len(lefts)):
        left = lefts[j]
        right = rights[j]

        for l in left:
            if l.score_total >= thr:
                sofar[l.idx] = max_sofar + l.w

        for r in right:
            if r.score_total >= thr:
                map_intervals[r.idx] = r
                if sofar[r.idx] > max_sofar:
                    max_sofar = sofar[r.idx]
                    max_interval = r

    if max_interval is None:
        return 0, []

    intervals = [max_interval]
    max_w = max_sofar - max_interval.w
    while True:
        for j in range(intervals[-1].idx-1, -1, -1):
            if abs(sofar[j] - max_w) <= 1e-9 and j in map_intervals and map_intervals[j].xb < intervals[-1].xa:
                intervals.append(map_intervals[j])
                max_w -= map_intervals[j].w
                break
        else:
            break
        if j == 0:
            break

    return max_sofar, intervals
