from dataclasses import dataclass
from typing import List, Tuple
import json
from util import render_color_seq
from heuristic import guess_orf1ab
from MWIS import Interval
import numpy as np
from pytablewriter import MarkdownTableWriter
from pytablewriter.style import Style
import sys


@dataclass
class BodyCore:
    ORF: str
    interval_start: int
    interval_len: int
    weight: float
    core_start: int
    core_len:int
    leader_start: int
    leader_end: int
    ORF_start: int
    ORF_len: int
    score: float
    index: int
    align: str


class Result:
    bodys: List[BodyCore]

    def __init__(self, leader_core_start, window, intervals: List[Interval], offset, max_weight, compact, ref, regions):
        self.leader_core_start = int(leader_core_start)
        self.leader_core_len = int(window)
        self.leader_core_seq = ref[leader_core_start:leader_core_start+window]
        min_core_s = min(i.ya for i in intervals)
        max_core_e = max(i.yb for i in intervals)
        self.TRS_L_start = min_core_s
        self.TRS_L_len = max_core_e - min_core_s + 1
        self.body_range_start = int(offset)
        self.body_range_len = int(len(ref) - offset) # FIXME: is ref necessary?
        self.n_intervals = len(intervals)
        self.weight = float(max_weight)
        self.compact = float(compact)
        if regions is None:
            regions = dict()

        self.bodys = []
        canon = set(regions)
        recall = set()
        for i in intervals:
            for orf, (start, end) in regions.items():
                # if not (i.end + offset) < start and not end < (i.start + offset):
                if int(i.orf_start + offset) == end:
                    assign = orf
                    recall.add(orf)
                    break
            else:
                assign = ""
            align = '-'*(i.ya - min_core_s) + ref[i.start+offset:i.end+offset+1] + '-'*(max_core_e - i.yb)
            self.bodys.append(BodyCore(
                assign,
                int(i.xa + offset),
                int(i.xb - i.xa + 1),
                i.w,
                int(i.start+offset),
                int(i.end - i.start + 1),
                int(i.ya),
                int(i.yb),
                int(i.orf_start + offset),
                int(i.length),
                i.score,
                i.idx,
                align
            ))
        self.total_ORF = list(canon)
        self.recall_ORF = list(recall)
        self.missing_ORF = list(canon - recall)
        self.TRS_L_seq = ref[min_core_s:max_core_e+1]

    def write_result(self, file=sys.stdout):
        print(f"Leader core seq: {self.leader_core_seq}, "
              f"{self.leader_core_start}-{self.leader_core_start + self.leader_core_len}, "
              f"TRS-L: {self.TRS_L_start}-{self.TRS_L_start + self.TRS_L_len}\n"
              f"max weight: {self.weight}, compact: {self.compact:.5f}, Missing: {self.missing_ORF}",
              file=file)
        header = ["ORF", "interval start-end", "interval len", "weight",
                "core start-end", "leader start-end", "ORF start-end", "score", "ORF len", "index", render_color_seq(self.TRS_L_seq)]
        table = []
        for body in self.bodys:
            # TODO: fix hard-coded index
            table.append([
                body.ORF,
                f"{body.interval_start + 1}-{body.interval_start + body.interval_len + 1}",
                body.interval_len,
                body.weight,
                f"{body.core_start + 1}-{body.core_start + body.core_len + 1}",
                f"{body.leader_start + 1}-{body.leader_end + 1}",
                f"{body.ORF_start + 1}-{body.ORF_start + body.ORF_len + 1}",
                body.score, body.ORF_len, body.index,
                render_color_seq(body.align)
            ])
        writer = MarkdownTableWriter(
            headers=header,
            value_matrix=table,
            margin=1,
        )
        writer.column_styles = [Style(align="center"), Style(align="center"), Style(),  Style(),  Style(align="center"),
                                Style(align="center"),  Style(),  Style(),  Style(),  Style()]
        writer.stream = file
        writer.write_table()

    def to_json(self):
        return json.dumps(self, default=lambda o: o.__dict__, sort_keys=True, indent=4)


class ResultEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, np.integer):
            return int(obj)
        elif isinstance(obj, np.floating):
            return float(obj)
        elif isinstance(obj, np.ndarray):
            return obj.tolist()
        elif isinstance(obj, Result) or isinstance(obj, Results) or isinstance(obj, BodyCore):
            return obj.__dict__
        else:
            return super().default(obj)


class Results:
    results: List[Result]

    def __init__(self, results: List[Result], name: str, description: str, annotation, ORF1ab: Tuple, reference: str):
        self.results = results
        self.name = name
        self.description = description
        self.genbank = annotation
        self.ORF1ab = ORF1ab
        self.sequence = reference
    
    def write_result(self, file=sys.stdout):
        for res in self.results:
            res.write_result(file=file)

    def to_json(self):
        return json.dumps(self, cls=ResultEncoder)

@dataclass
class Solution:
    """
    A solution for the MWIS problem, including a TRS-L, and several TRS-Bs.
    """
    weights: List[int]
    intervals: List[list]
    offset: int

    def __post_init__(self):
        if len(self.weights) != len(self.intervals):
            raise RuntimeError(f"Weights and intervals are not of the same length: {len(self.weights)} != {len(self.intervals)}") # TODO: fix the wording
        self.weights = np.array(self.weights)

    def argsort(self, is_lex=False, compact=None) -> np.ndarray:
        if is_lex and compact is not None:
            min_score = np.array([min(i.score for i in sol) if sol is not None and len(sol) > 0 else 0 for sol in self.intervals])
            return np.lexsort((-min_score, -self.weights, compact))
        else:
            return np.argsort(self.weights)[::-1]

    def serialize_results(self, ref, window, name, description, regions=None, annotation=None, is_lex=False, compact=None) -> Results:
        results = []
        idx = self.argsort(is_lex, compact)

        for pos in idx:
            if self.intervals[pos] is None or len(self.intervals[pos]) == 0:
                continue
            results.append(Result(pos, window, self.intervals[pos], self.offset, self.weights[pos], compact[pos], ref, regions))
        start_1a, end_1a, end_1b = guess_orf1ab(ref)
        return Results(results, name, description, annotation, (start_1a, end_1a, end_1a, end_1b), ref)

    def get_compact_score(self, ref: str, orf1_end: int):
        compact_score = []
        total_length = len(ref) - orf1_end
        for i, sol in enumerate(self.intervals):
            if sol is None:
                compact_score.append(1)
            else:
                cov = np.zeros(total_length)
                for intv in sol:
                    cov[intv.orf_start: intv.orf_end + 1] += 1
                compact_score.append((cov==0).mean())
        return np.array(compact_score)
