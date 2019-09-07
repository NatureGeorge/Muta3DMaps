# @Date:   2019-09-03T16:37:05+08:00
# @Email:  1730416009@stu.suda.edu.cn
# @Filename: PdSeq_unit.py
# @Last modified time: 2019-09-07T23:27:29+08:00
from Unit import Unit
import sys
from Bio import pairwise2
from Bio.SubsMat import MatrixInfo as matlist
from Bio import Align
import json
sys.path.append('./')


class RowChecker:
    def __init__(self, colName):
        self.colName = colName
        self.content = ''

    def check(self, row):
        temp = row[self.colName]
        if temp != self.content:
            self.content = temp
            return 1
        else:
            return 0


class IndexRecorder:
    def __init__(self):
        self.count = 0

    def check(self, seq):
        if seq != '-':
            self.count += 1
            return self.count
        else:
            return -1

    def reset(self):
        self.count = 0


class PdSeqAlign:
    def getAlignmentSegment(alignment):
        segments1 = []
        segments2 = []
        i1, i2 = alignment.path[0]
        for node in alignment.path[1:]:
            j1, j2 = node
            if j1 > i1 and j2 > i2:
                segment1 = [i1+1, j1]
                segment2 = [i2+1, j2]
                segments1.append(segment1)
                segments2.append(segment2)
            i1, i2 = j1, j2
        return segments1, segments2

    def __init__(self):
        self.seqa = ''
        self.seqb = ''
        self.aligner = Align.PairwiseAligner()
        self.aligner.mode = 'global'
        self.aligner.open_gap_score = -10
        self.aligner.extend_gap_score = -0.5
        self.aligner.substitution_matrix = matlist.blosum62

    def makeAlignment_pairwise2(self, seqa, seqb):
        if seqa == self.seqa and seqb == self.seqb:
            return self.range_a, self.range_b
        else:
            self.seqa = seqa
            self.seqb = seqb
            self.alns = pairwise2.align.globalds(seqa, seqb, matlist.blosum62, -10, -0.5)
            # aln_seqa, aln_seqb, score, begin, end = alns[0]
            indexRecord = IndexRecorder()
            a = list(map(indexRecord.check, self.alns[0][0]))
            indexRecord.reset()
            b = list(map(indexRecord.check, self.alns[0][1]))
            indexRecord.reset()
            mapped = list(filter(lambda x: x[0] != -1 and x[1] != -1, zip(a, b)))
            self.range_a = Unit.getInterval([i[0] for i in mapped])
            self.range_b = Unit.getInterval([i[1] for i in mapped])
            return self.range_a, self.range_b

    def makeAlignment_align(self, seqa, seqb):
        if not (seqa == self.seqa and seqb == self.seqb):
            self.seqa = seqa
            self.seqb = seqb
            alignments = self.aligner.align(seqa, seqb)
            result = PdSeqAlign.getAlignmentSegment(alignments[0])
            self.range_a, self.range_b = json.dumps(result[0]), json.dumps(result[1])

        return self.range_a, self.range_b
