# @Date:   2019-09-03T16:37:05+08:00
# @Email:  1730416009@stu.suda.edu.cn
# @Filename: PdSeq_unit.py
# @Last modified time: 2019-09-03T20:41:50+08:00
import json
import pandas as pd
import numpy as np
from Bio import pairwise2
from Bio.SubsMat import MatrixInfo as matlist
sys.path.append('./')
from Unit import Unit


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


class PdSeq:
    def __init__(self, seqCols):
        self.seqCola, self.seqColb = seqCols
        self.seqa = ''
        self.seqb = ''

    def makeAlignment(self, row):
        seqa = row[self.seqCola]
        seqb = row[self.seqColb]
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
