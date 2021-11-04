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
    score: int
    length: int
    orf_start: int
    orf_end: int

    def __post_init__(self):
        self.orf_len = self.orf_end - self.orf_start + 1

    def __str__(self):
        return f"{self.idx}\t{self.xa}-{self.xb}\t{self.w}\t{self.ya}-{self.yb}\t{self.start}-{self.end}\t{self.orf_start}-{self.orf_end}\t{self.score}\t{self.length}"

    def __repr__(self):
        return f"<{self.idx}|xa={self.xa}, xb={self.xb}, w={self.w}, ya={self.ya}, yb={self.yb}, start={self.start}, end={self.end}, start={self.orf_start}, end={self.orf_end}, score={self.score}, len={self.length}>"

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


def MWIS(lefts: List[List[Interval]], rights: List[List[Interval]], thr: float=0) -> Tuple[float, List[Interval]]:
    """Maximum weight independent set.

    Intervals should be sorted by their right endpoints.

    Args:
        lefts (List[List[Interval]]): intervals indexed by left endpoints
        rights (List[List[Interval]]): intervals indexed by left endpoints
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
            if l.score >= thr:
                sofar[l.idx] = max_sofar + l.w
        for r in right:
            if r.score >= thr:
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
