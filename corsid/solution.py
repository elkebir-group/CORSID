from dataclasses import dataclass
from typing import List, Tuple
import json
from util import render_color_seq
from heuristic import guess_orf1ab
from MWIS import Interval, Dual_Interval
import numpy as np
from pytablewriter import MarkdownTableWriter
from pytablewriter.style import Style
import sys


@dataclass
class BodyCore:
    ORF_name: str
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
                int(i.orf_length),
                i.score_total,
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
                body.ORF_name,
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
        writer.column_styles = [Style(align="center")] * len(header)
        writer.stream = file
        writer.write_table()

    def to_json(self) -> str:
        """Dump result to a json string"""
        return json.dumps(self, default=lambda o: o.__dict__, sort_keys=True, indent=4)


@dataclass
class Dual_BodyCore:
    ORF_name: str
    ORF_start: int
    ORF_len: int
    candidate_start: int
    candidate_len: int
    core_start1: int
    core_len1: int
    core_start2: int
    core_len2: int
    leader_start1: int
    leader_end1: int
    leader_start2: int
    leader_end2: int
    score1: float
    score2: float
    weight: float
    align1: str
    align2: str


class Dual_Result(Result):
    """Single result, contain a TRS alignment and induced genes"""
    bodys: List[Dual_BodyCore]

    def __init__(self, window, intervals: List[Dual_Interval], offset, max_weight, compact, ref, regions):
        self.leader_core_start1 = int([i.idx1 for i in intervals if i.idx1][0])
        self.leader_core_seq1 = ref[self.leader_core_start1:self.leader_core_start1+window]
        self.leader_core_start2 = int([i.idx2 for i in intervals if i.idx2][0])
        self.leader_core_seq2 = ref[self.leader_core_start2:self.leader_core_start2+window]
        self.leader_core_len = int(window)

        min_core1_s = min(i.ya1 for i in intervals if i.ya1)
        max_core1_e = max(i.yb1 for i in intervals if i.yb1)
        self.TRS1_L_start = min_core1_s
        self.TRS1_L_len = max_core1_e - min_core1_s + 1
        self.TRS1_L_seq = ref[min_core1_s:max_core1_e+1]

        min_core2_s = min(i.ya2 for i in intervals if i.ya2)
        max_core2_e = max(i.yb2 for i in intervals if i.yb2)
        self.TRS2_L_start = min_core2_s
        self.TRS2_L_len = max_core2_e - min_core2_s + 1
        self.TRS2_L_seq = ref[min_core2_s:max_core2_e+1]
        
        self.body_range_start = int(offset)
        self.body_range_len = int(len(ref) - offset)

        intv_scores = [i.score1 for i in intervals if i.score1] + [i.score2 for i in intervals if i.score2]
        self.min_score = float(min(intv_scores))
        self.weight = float(max_weight)
        self.compact = float(compact)
        if regions is None:
            regions = dict()
        self.bodys = []
        canon = set(regions)
        recall = set()

        for i in intervals:
            for orf, (start, end) in regions.items():
                if (abs(i.orf_start + offset - start) <= 3 and abs(i.orf_end + offset - end) <= 3):
                    assign = orf
                    recall.add(orf)
                    break
            else:
                assign = ''

            if i.start1:
                TRS1_B_start = int(i.start1 + offset)
                TRS1_B_len = int(i.end1 - i.start1)
                align1 = '-'*(i.ya1 - min_core1_s) + ref[i.start1+offset:i.end1+offset+1] + '-'*(max_core1_e - i.yb1)
                # align1 = ref[i.start1+offset - (i.ya1 - min_core1_s):i.end1+offset + (max_core1_e - i.yb1) + 1]
            else:
                TRS1_B_start = None
                TRS1_B_len = None
                align1 = ''

            if i.start2:
                TRS2_B_start = int(i.start2 + offset)
                TRS2_B_len = int(i.end2 - i.start2)
                align2 = '-'*(i.ya2 - min_core2_s) + ref[i.start2+offset:i.end2+offset+1] + '-'*(max_core2_e - i.yb2)
                # align2 = ref[i.start2+offset - (i.ya2 - min_core2_s):i.end2+offset + (max_core2_e - i.yb2) + 1]
            else:
                TRS2_B_start = None
                TRS2_B_len = None
                align2 = ''

            self.bodys.append(Dual_BodyCore(
                assign,
                int(i.orf_start + offset),
                int(i.orf_length),
                int(i.xa + offset),
                int(i.xb - i.xa),
                TRS1_B_start,
                TRS1_B_len,
                TRS2_B_start,
                TRS2_B_len,
                i.ya1,
                i.yb1,
                i.ya2,
                i.yb2,
                i.score1,
                i.score2,
                i.w,
                align1,
                align2
            ))

        self.total_ORF = list(canon)
        self.recall_ORF = list(recall)
        self.missing_ORF = list(canon - recall)

    def write_result(self, file=sys.stdout):
        """Write result to stdout (default) or a file"""
        print(f"\nLeader core seq: {self.leader_core_seq1}, "
              f"Core1: {self.leader_core_start1}-{self.leader_core_start1 + self.leader_core_len - 1}, "
              f"TRS1-L: {self.TRS1_L_start}-{self.TRS1_L_start + self.TRS1_L_len - 1}\n"
              f"Leader core seq: {self.leader_core_seq2}, "
              f"Core2: {self.leader_core_start2}-{self.leader_core_start2 + self.leader_core_len - 1}, "
              f"TRS2-L: {self.TRS2_L_start}-{self.TRS2_L_start + self.TRS2_L_len - 1}\n"
              f"Weight: {self.weight:.3f}, Coverage: {1-self.compact:.5f},  "
              f"Min Score: {self.min_score}, Missing: {self.missing_ORF}",
              file=file)
        
        header = ["Gene", "ORF", "Candidate Region",
                  "TRS1-B", "Leader Match", "Score",
                  "TRS2-B", "Leader Match", "Score", 
                  "Total", render_color_seq(self.TRS1_L_seq), render_color_seq(self.TRS2_L_seq)]
        table = []

        for body in self.bodys:

            if body.core_start1:
                core_range1 = f"{body.core_start1}-{body.core_start1 + body.core_len1}"
                leader_range1 = f"{body.leader_start1}-{body.leader_end1}"
            else:
                core_range1 = ""
                leader_range1 = ""

            if body.core_start2:
                core_range2 = f"{body.core_start2}-{body.core_start2 + body.core_len2}"
                leader_range2 = f"{body.leader_start2}-{body.leader_end2}"
            else:
                core_range2 = ""
                leader_range2 = ""

            table.append([
                body.ORF_name,
                f"{body.ORF_start + 1}-{body.ORF_start + body.ORF_len + 1}",
                # body.ORF_len,
                f"{body.candidate_start}-{body.candidate_start + body.candidate_len}",
                # body.candidate_len,
                core_range1,
                leader_range1,
                body.score1,
                core_range2,
                leader_range2,
                body.score2,
                f"{body.weight:.5f}", 
                render_color_seq(body.align1),
                render_color_seq(body.align2)
            ])

        writer = MarkdownTableWriter(headers=header, value_matrix=table, margin=1)
        writer.column_styles = [Style(align="center")] * len(header)
        writer.stream = file
        writer.write_table()


class ResultEncoder(json.JSONEncoder):
    """Encoder for class `Results`"""
    def default(self, obj):
        if isinstance(obj, np.integer):
            return int(obj)
        elif isinstance(obj, np.floating):
            return float(obj)
        elif isinstance(obj, np.ndarray):
            return obj.tolist()
        elif (isinstance(obj, Result) or isinstance(obj, Results) or isinstance(obj, BodyCore)
              or isinstance(obj, Dual_Result) or isinstance(obj, Dual_BodyCore)):
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
        for i, body in enumerate(self.results[0].bodys[::-1]):
            if body.ORF_name:
                gene_id = f"gene-{body.ORF_name}"
            else:
                gene_id = f"putative-gene-{i+2}"
            attributes = ["ID="+gene_id, "gene_biotype=protein_coding"]
            gff.append(f"{self.name}\t.\tgene\t{body.ORF_start + 1}\t{body.ORF_start + body.ORF_len + 3}\t.\t"
                       f"+\t.\t{';'.join(attributes)}")

            attributes = [f"ID=cds-{i+2}", f"Parent={gene_id}"]
            if body.ORF_name:
                attributes.append(f"Name={body.ORF_name}")
                attributes.append(f"gene={body.ORF_name}")
            else:
                attributes.append(f"Name=protein-{i+2}")
                attributes.append(f"gene=ORF{i+2}")
            gff.append(f"{self.name}\t.\tCDS\t{body.ORF_start + 1}\t{body.ORF_start + body.ORF_len + 3}\t.\t"
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
    window: int
    is_corsid_a: bool
    is_dual: bool

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
            min_score = np.array([min(i.score_total for i in sol) if sol is not None and len(sol) > 0 else 0 for sol in self.intervals])
            align_score = (self.weights - self.weights.astype(int)) * 1000
            return np.lexsort((-min_score, -align_score, compact))
        else:
            return np.argsort(self.weights)[::-1]

    def serialize_results(self, ref: str, window: int, name: str, description: str,
                          annotation=None, is_lex=False, compact=None) -> Results:
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

        for pos in idx[:100]:
            if self.intervals[pos] is None or len(self.intervals[pos]) == 0:
                continue
            if self.is_dual:
                results.append(Dual_Result(window, self.intervals[pos], self.offset, self.weights[pos], compact[pos], ref, annotation))
            else:
                results.append(Result(pos, window, self.intervals[pos], self.offset, self.weights[pos], compact[pos], ref, annotation))

        start_1a, end_1a, end_1b = guess_orf1ab(ref)
        return Results(results, name, description, annotation, (start_1a, end_1a, end_1a, end_1b), ref, self.is_corsid_a)

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
        for sol in self.intervals:
            if sol is None:
                compact_score.append(1)
            else:
                cov = np.zeros(total_length)
                for intv in sol:
                    cov[max(0, intv.orf_start + offset): intv.orf_end + offset + 1] += 1
                compact_score.append((cov==0).mean())
        return np.array(compact_score)
