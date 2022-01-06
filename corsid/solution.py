from dataclasses import dataclass
from typing import List, Sequence, Tuple
import json
from .util import render_color_seq
from .heuristic import guess_orf1ab
from .MWIS import Interval
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
    """Single result, contain a TRS alignment and induced genes"""
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
        self.body_range_len = int(len(ref) - offset)
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
                if (abs(i.orf_start + offset - start) <= 3 and
                    abs(i.orf_end + offset - end) <= 3
                ):
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
        """Write result to stdout (default) or a file"""
        print(f"Leader core seq: {self.leader_core_seq}, "
              f"{self.leader_core_start}-{self.leader_core_start + self.leader_core_len}, "
              f"TRS-L: {self.TRS_L_start}-{self.TRS_L_start + self.TRS_L_len}\n"
              f"weight: {self.weight:.3f}, coverage: {1-self.compact:.5f}, Missing: {self.missing_ORF}, "
              f"min score: {min(body.score for body in self.bodys)}",
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

    def to_json(self) -> str:
        """Dump result to a json string"""
        return json.dumps(self, default=lambda o: o.__dict__, sort_keys=True, indent=4)


class ResultEncoder(json.JSONEncoder):
    """Encoder for class `Results`"""
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


@dataclass
class Results:
    """
    Collection of results from a genome, with additional meta info
    """
    results: List[Result]
    name: str
    description: str
    genbank: dict
    ORF1ab: Tuple[int]
    sequence: str
    is_corsid_a: bool
    
    def write_result(self, file=sys.stdout):
        for res in self.results:
            res.write_result(file=file)

    def to_json(self):
        return json.dumps(self, cls=ResultEncoder)

    def to_gff(self) -> str:
        """Convert result to GFF3 format"""
        gff = ["##gff-version 3.1.26"]
        gff.append(f"##sequence-region {self.name} {1} {len(self.sequence)}")
        gff.append(f"{self.name}\t.\tgene\t{self.ORF1ab[0] + 1}\t{self.ORF1ab[3] + 3}\t.\t"
                   f"+\t.\tID=gene-1;Name=ORF1;gene_biotype=protein_coding")
        gff.append(f"{self.name}\t.\tCDS\t{self.ORF1ab[0] + 1}\t{self.ORF1ab[1] + 1}\t.\t"
                   f"+\t.\tID=cds-1a;Parent=gene-1;Name=ORF1a")
        gff.append(f"{self.name}\t.\tCDS\t{self.ORF1ab[2] + 1}\t{self.ORF1ab[3] + 3}\t.\t"
                   f"+\t.\tID=cds-1b;Parent=gene-1;Name=ORF1b")
        for i, body in enumerate(self.results[0].bodys):
            if body.ORF:
                gene_id = f"gene-{body.ORF}"
            else:
                gene_id = f"putative-gene-{i+2}"
            attributes = ["ID="+gene_id, "gene_biotype=protein_coding"]
            gff.append(f"{self.name}\t.\tgene\t{body.ORF_start + 1}\t{body.ORF_start + body.ORF_len + 1}\t.\t"
                       f"+\t.\t{';'.join(attributes)}")

            attributes = [f"ID=cds-{i+2}", f"Parent={gene_id}"]
            if body.ORF:
                attributes.append(f"Name={body.ORF}")
                attributes.append(f"gene={body.ORF}")
            else:
                attributes.append(f"Name=protein-{i+2}")
                attributes.append(f"gene=ORF{i+2}")
            gff.append(f"{self.name}\t.\tCDS\t{body.ORF_start + 1}\t{body.ORF_start + body.ORF_len + 1}\t.\t"
                       f"+\t0\t{';'.join(attributes)}")
        return '\n'.join(gff)

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
        """Sort results

        Args:
            is_lex (bool, optional): set to True to sort in lexicographical order of (compact, weights, min alignment score). Defaults to False.
            compact (list, optional): list of compact score. Defaults to None.

        Returns:
            np.ndarray: index in sorted order
        """
        if is_lex and compact is not None:
            min_score = np.array([min(i.score for i in sol) if sol is not None and len(sol) > 0 else 0 for sol in self.intervals])
            return np.lexsort((-min_score, -self.weights, compact))
        else:
            return np.argsort(self.weights)[::-1]

    def serialize_results(self, ref: str, window: int, name: str, description: str,
                          annotation=None, is_lex=False, compact=None, is_corsid_a=False) -> Results:
        """Convert solution to `Results` in order to save

        Args:
            ref (str): reference genome
            window (int): window length
            name (str): name of the genome
            description (str): description of the virus
            annotation (dict, optional): genes in the annotation file. Defaults to None.
            is_lex (bool, optional): set to True to enable lexicographical sorting, need compact as well. Defaults to False.
            compact (ndarray, optional): list of genome coverage. Defaults to None.
            is_corsid_a (bool, optional): is the solution created using CORSID-A. Defaults to False.

        Returns:
            Results: collection of TRS alignments
        """
        results = []
        idx = self.argsort(is_lex, compact)

        for pos in idx:
            if self.intervals[pos] is None or len(self.intervals[pos]) == 0:
                continue
            results.append(Result(pos, window, self.intervals[pos], self.offset, self.weights[pos], compact[pos], ref, annotation))
        start_1a, end_1a, end_1b = guess_orf1ab(ref)
        return Results(results, name, description, annotation, (start_1a, end_1a, end_1a, end_1b), ref, is_corsid_a)

    def get_compact_score(self, ref: str, orf1_end: int, offset: int=0) -> List[float]:
        """Return genome coverage of solutions

        Args:
            ref (str): reference genome
            orf1_end (int): the position of ORF1ab end

        Returns:
            List[float]: array of genome coverage
        """
        compact_score = []
        total_length = len(ref) - orf1_end
        for i, sol in enumerate(self.intervals):
            if sol is None:
                compact_score.append(1)
            else:
                cov = np.zeros(total_length)
                for intv in sol:
                    cov[max(0, intv.orf_start + offset): intv.orf_end + offset + 1] += 1
                compact_score.append((cov==0).mean())
        return np.array(compact_score)
