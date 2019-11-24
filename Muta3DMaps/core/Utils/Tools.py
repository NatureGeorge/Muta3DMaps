# @Date:   2019-08-16T23:19:34+08:00
# @Email:  1730416009@stu.suda.edu.cn
# @Filename: Tools.py
# @Last modified time: 2019-11-20T19:31:42+08:00
import pandas as pd
import numpy as np
import json


class Gadget:

    SEQ_DICT = {
        "GLY": "G", "ALA": "A", "SER": "S", "THR": "T", "CYS": "C", "VAL": "V", "LEU": "L",
        "ILE": "I", "MET": "M", "PRO": "P", "PHE": "F", "TYR": "Y", "TRP": "W", "ASP": "D",
        "GLU": "E", "ASN": "N", "GLN": "Q", "HIS": "H", "LYS": "K", "ARG": "R", "HSD": "H",
        "A": "A", "G": "G", "C": "C", "T": "T", "U": "U", "DA": "A", "DT": "T",
        "DU": "U", "DC": "C", "DG": "G", "DI": "I", "?": "?", "UNK": "!"}

    class GroupER:
        def __init__(self, name_tp, filter_ele, dfrm, new_col):
            self.name = name_tp
            self.filter_ele = filter_ele
            self.content = []
            self.index = []
            self.dfrm = dfrm
            self.new_col = new_col

        def output(self):
            if self.index:
                self.dfrm.loc[self.index, self.new_col] = str(
                    self.content).replace('\'', '"')

        def check(self, tp, data):  # (index, type, len, chain_id)
            if self.name == tp:
                self.index.append(data[0])
                if data[1] in self.filter_ele:
                    self.content.append(data[2:])
            else:
                self.output()
                self.name = tp

                if data[1] in self.filter_ele:
                    self.content = [data[2:]]
                else:
                    self.content = []

                self.index = [data[0]]

    class RangeSetER:
        def __init__(self, name_tp):
            self.name = name_tp  # ('pdb_id', 'chain_id', 'UniProt')
            self.pdb_range = []
            self.unp_range = []

        def output(self):
            if self.pdb_range:
                pdbRange = json.dumps(self.pdb_range)
                unpRange = json.dumps(self.unp_range)
                return '%s|%s' % (pdbRange, unpRange)
            else:
                return '%s|%s' % (self.temp1, self.temp2)

        def check(self, tp_1, tp_2):
            self.temp1 = '[[%s, %s]]' % tp_2[:2]
            self.temp2 = '[[%s, %s]]' % tp_2[2:4]

            if self.name == tp_1:
                self.pdb_range.append([int(tp_2[0]), int(tp_2[1])])
                self.unp_range.append([int(tp_2[2]), int(tp_2[3])])
                out = self.output()
            else:
                self.name = tp_1
                self.pdb_range = [[int(tp_2[0]), int(tp_2[1])]]
                self.unp_range = [[int(tp_2[2]), int(tp_2[3])]]
                out = self.output()

            return out

    class ConstraintDict:
        """
        This Class is a new kind of dictionary.
        Initial Dict: Stable
        Input Dict: Mutal
        """

        def __init__(self, iniDict):
            self.iniDict = iniDict
            self.xkeys = self.iniDict.keys()
            self.xitems = self.iniDict.items()
            self.mutalDict = {}

        def setMutalDict(self, inputDict):
            self.mutalDict = inputDict
            self.xkeys = list(self.iniDict.keys())
            self.xkeys.extend(list(self.mutalDict.keys()))
            self.xitems = list(self.iniDict.items())
            self.xitems.extend(list(self.mutalDict.items()))

        def __getitem__(self, key):
            self.result = self.iniDict.get(key, False)
            if self.result:
                return self.result
            else:
                return self.mutalDict[key]

        def keys(self):
            return self.xkeys

        def items(self):
            return self.xitems

        def addConstraintToDf(df, constraint_dict):
            for i, j in constraint_dict.items():
                j1, j2 = j
                if j2 == 'eq':
                    df = df[df[i] == j1]
                elif j2 == 'ne':
                    df = df[df[i] != j1]
                elif j2 == 'gt':
                    df = df[df[i] > j1]
                elif j2 == 'lt':
                    df = df[df[i] < j1]
                elif j2 == 'ge':
                    df = df[df[i] >= j1]
                elif j2 == 'le':
                    df = df[df[i] <= j1]
            return df

    class MultiToOne:
        def __init__(self):
            self.aa_map = Gadget.SEQ_DICT

        def multi_letter_convert_to_one_letter(self, a):
            return self.aa_map.get(a, 'X')

    def handleResolution(resolution):
        def float_fun(x): return float(x) if x not in '?.' else 1000
        if pd.isna(resolution):
            return 1000
        elif isinstance(resolution, str):
            if ',' in resolution:
                reso_li = map(float_fun, resolution.split(','))
                return min(reso_li)
            elif resolution in '?.':
                return 1000
            else:
                return float(resolution)
        else:
            return resolution

    def selectChain(grouped_df, df, rank_list, rankName, rankFormat, rangeName, selectName, r1_cutoff, r2_cutoff):
        def getRange(li):
            rangeSet = set()
            for ran in li:
                rangeSet = rangeSet | set(range(ran[0], ran[1] + 1))
            return rangeSet
        '''
        rankName: 'rank'
        rankFormat: '%d-%d-%d-%d'
        rangeName: seg_unp_range
        selectName: pdb_chain_select
        r1_cutoff: 0.3
        r2_cutoff: 0.2
        # overlap_type: ['shorter']
        '''
        # overlap_list = []
        if len(grouped_df) > 1:
            rank_df = pd.DataFrame([grouped_df[ele].rank(ascending=0, method='dense')
                                    for ele in rank_list]).T
            if rankName:
                df.loc[rank_df.index, rankName] = rank_df.apply(
                    lambda x: rankFormat % tuple(x[ele] for ele in rank_list), axis=1)
            index_list = grouped_df.sort_values(
                by=rank_list, ascending=False).index
            # FIRST
            repreSet = getRange(json.loads(
                grouped_df.loc[index_list[0], rangeName]))
            df.loc[index_list[0], selectName] = True
            # LATTER
            for tr in index_list[1:]:
                temp_range = getRange(json.loads(
                    grouped_df.loc[tr, rangeName]))
                if temp_range <= repreSet:
                    continue
                else:
                    # overlap = getOverlap(temp_range, repreSet, overlap_type)
                    # overlap_list.append(overlap)
                    overlap = len(temp_range & repreSet)
                    temp_range_len = len(temp_range)
                    # ---------------------------------------------------------------
                    if (overlap / temp_range_len <= r1_cutoff) and ((temp_range_len - overlap) / len(repreSet) >= r2_cutoff):
                        # ---------------------------------------------------------------
                        repreSet = repreSet | temp_range
                        df.loc[tr, selectName] = True
        else:
            df.loc[grouped_df.index, selectName] = True

    def getInterval(rangeSet):
        if rangeSet == '' or rangeSet == set() or rangeSet == [] or isinstance(rangeSet, float):
            return np.nan
        else:
            start = []
            interval_list = []
            true_interval_list = []
            maxRange = max(rangeSet)
            minRange = min(rangeSet)
            if len(rangeSet) == (maxRange + 1 - minRange):
                true_interval_list.append([minRange, maxRange])
            else:
                rangeSet_list = list(rangeSet)
                rangeSet_list.sort()
                for j in rangeSet_list:
                    if not start:
                        i = j
                        start.append(j)
                        i += 1
                    else:
                        if (i != j) or (j == max(rangeSet_list)):
                            if j == max(rangeSet_list):
                                if (i != j):
                                    interval_list.append(start)
                                    interval_list.append([j])
                                    break
                                else:
                                    start.append(j)
                            interval_list.append(start)
                            start = [j]
                            i = j + 1
                        else:
                            start.append(j)
                            i += 1
                for li in interval_list:
                    maxRange = max(li)
                    minRange = min(li)
                    true_interval_list.append([minRange, maxRange])
            return true_interval_list

    def set_lists(self, pdb_list, unp_list):
        self.unp_list = unp_list
        self.pdb_list = pdb_list

    def __init__(self):
        self.file_list = []
