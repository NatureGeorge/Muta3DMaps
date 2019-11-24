# @Date:   2019-08-13T16:19:30+08:00
# @Email:  1730416009@stu.suda.edu.cn
# @Filename: SIFTS_Utils.py
# @Last modified time: 2019-08-16T23:18:29+08:00
import pandas as pd
import numpy as np
import json, os, re, wget, gzip, time, requests
from urllib import request, error
from retrying import retry
from multiprocessing.dummy import Pool
from bs4 import BeautifulSoup
from Bio.PDB import PDBList
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from itertools import combinations


class Unit:

    class RangeSetER:
        def __init__(self, name_tp):
            self.name = name_tp# ('pdb_id', 'chain_id', 'UniProt')
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
                out = self.output()
                self.name = tp_1
                self.pdb_range = [[int(tp_2[0]), int(tp_2[1])]]
                self.unp_range = [[int(tp_2[2]), int(tp_2[3])]]

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
            for i,j in constraint_dict.items():
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
            self.aa_map = {
                "GLY": "G", "ALA": "A", "SER": "S", "THR": "T", "CYS": "C", "VAL": "V", "LEU": "L",
                "ILE": "I", "MET": "M", "PRO": "P","PHE": "F", "TYR": "Y", "TRP": "W", "ASP": "D",
                "GLU": "E", "ASN": "N", "GLN": "Q", "HIS": "H", "LYS": "K", "ARG": "R", "HSD": "H",
                "A": "A", "G": "G", "C": "C", "T": "T", "U": "U", "DA": "DA", "DT": "DT",
                "DU": "DU", "DC": "DC", "DG": "DG","DI":"DI","?":"?","UNK":"!"}

        def multi_letter_convert_to_one_letter(self, a):
            return self.aa_map.get(a, 'X')

    def selectChain(grouped_df, df, rank_list, rankName, rankFormat, rangeName, selectName, r1_cutoff, r2_cutoff):
        def getRange(li):
            rangeSet = set()
            for ran in li:
                rangeSet = rangeSet | set(range(ran[0], ran[1]+1))
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
        #overlap_list = []
        if len(grouped_df) > 1:
            rank_df = pd.DataFrame([grouped_df[ele].rank(ascending=0, method='dense')
                                    for ele in rank_list]).T
            if rankName:
                df.loc[rank_df.index, rankName] = rank_df.apply(
                    lambda x: rankFormat % tuple(x[ele] for ele in rank_list), axis=1)
            index_list = grouped_df.sort_values(by=rank_list, ascending=False).index
            # FIRST
            repreSet = getRange(json.loads(grouped_df.loc[index_list[0], rangeName]))
            df.loc[index_list[0], selectName] = True
            # LATTER
            for tr in index_list[1:]:
                temp_range = getRange(json.loads(grouped_df.loc[tr, rangeName]))
                if temp_range <= repreSet:
                    continue
                else:
                    #overlap = getOverlap(temp_range, repreSet, overlap_type)
                    #overlap_list.append(overlap)
                    overlap = len(temp_range & repreSet)
                    temp_range_len = len(temp_range)
                    #---------------------------------------------------------------
                    if (overlap/temp_range_len <= r1_cutoff) and ((temp_range_len - overlap)/len(repreSet) >= r2_cutoff):
                    #---------------------------------------------------------------
                        repreSet = repreSet | temp_range
                        df.loc[tr, selectName] = True
        else:
            df.loc[grouped_df.index, selectName] = True

    def getInterval(rangeSet):
        if rangeSet == '' or rangeSet == set():
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

    def file_i(self, path, df, va_tp, sep='\t'):
        try:
            dfrm = pd.read_csv(path, sep=sep)
        except Exception:
            dfrm = df
            if not isinstance(dfrm, pd.DataFrame):
                raise Exception('Input: %s=pd.DataFrame() or %s=str' % va_tp)
        return dfrm

    def file_o(self, path, df):
        if path:
            df.to_csv(path, sep='\t', index=False)
            self.file_list.append(path)

    def __init__(self):
        self.file_list = []


class SIFTS_unit(Unit):
    CONFIG = {
        'PDB_ID': 'pdb_id',
        'RAW_SIFTS_COLUMNS':[
            'pdb_id', 'chain_id', 'UniProt', 'identity', 'identifier',
            'pdb_start', 'pdb_end', 'unp_start', 'unp_end',
            'is_canonical', 'start', 'end', 'entity_id', 'struct_asym_id'],
        'DOWNLOAD_FOLDER': '/home/zzf/Work/StructDDG_0427/data/Mapping_Pipeline/sifts_files/',
        'UNP_LIST_PATH': '/home/zzf/Work/StructDDG_0427/data/Mapping_Pipeline/sifts_uniprot_list.tsv',
        'ELE_LIST': ["coordinates_len", "mappedOut", "metal_count", "delHT_MissingNum", "if_2", "if_1"],
        'INI_LIST': [1, 1, 2, 3, 4, 1, 2, 3, 4, 2, 3, 4, 2, 3, 2],
        'SEG_SIFTS_MAP': {
              "SP_PRIMARY":"UniProt",# (Canonical)
              "RES_BEG":"pdb_start",
              "RES_END":"pdb_end",
              "PDB_BEG":"residue index(start) in pdb",
              "PDB_END":"residue index(end) in pdb",
              "SP_BEG":"unp_start",
              "SP_END":"unp_end",
              "PDB": 'pdb_id',
              "CHAIN": 'chain_id',
              },
        'SCORE_COL': 'BS',
        'PDB_SELECT_COL': 'pdb_chain_select',
        'PDB_RANK_COL': 'pdb_chain_rank',
        'PDB_RANK_LIST':['BS', 'ne_resolution_score', 'initial_version_time'],
        'PDB_RANK_FORMAT': '%d-%d-%d',
        'PDB_INIT_GROUPBY_LIST': ['UniProt'],
        'PDB_SELECT_RANGE_NAME': 'seg_unp_range',

    }

    def get_raw_SIFTS(self, outputPath):
        # Have to use set_lists() before using this fun.
        def order_SIFTS_info(pdbId, info):
            first = True
            order = SIFTS_unit.CONFIG['RAW_SIFTS_COLUMNS']
            for uniprot in info.keys():
                gene = info[uniprot]['identifier']
                pdbChain = info[uniprot]['mappings']
                for chain in pdbChain:
                    chain['start'] = json.dumps(chain['start'])
                    chain['end'] = json.dumps(chain['end'])
                    chain[SIFTS_unit.CONFIG['PDB_ID']] = pdbId
                    chain['UniProt'] = uniprot
                    chain['identifier'] = gene
                    if first:
                        df = pd.DataFrame(chain, index=[0])
                        df = df[order]
                        first = False
                    else:
                        temp = pd.DataFrame(chain, index=[0])
                        df = df.append(temp[order])
            if 'df' in locals().keys():
                return df.reset_index().drop(['index'], axis=1)
            else:
                return False

        fail_list = []
        self.raw_SIFTS_filePath = outputPath
        allPDB, current = len(self.pdb_list), 0
        for pdbId in self.pdb_list:
            time.sleep(1.2)
            print('getSiftsInfo(): Start to get the pdb info from SIFTS.',pdbId)
            url = ('http://www.ebi.ac.uk/pdbe/api/mappings/all_isoforms/%s') % (pdbId)
            try:
                req = request.Request(url)
                page = request.urlopen(req).read()
                page = page.decode('utf-8')
                info = json.loads(page)[pdbId]['UniProt']
                if info:
                    df = order_SIFTS_info(pdbId, info)
                    if not isinstance(df, bool):
                        df.to_csv(outputPath, index=False, sep='\t', mode='a+')
            except Exception as e:
                print('PDB: [', pdbId, ']\n', e)
                fail_list.append(pdbId)
            current += 1
            print('getSiftsInfo(): End a circle.[',pdbId,'] current:',current,'ALL:',allPDB)
        return fail_list

    def get_info_from_uniprot_pdb_file(self, filePath=False, related_unp=False, related_pdb=False):
        '''
        Reference: (http://www.ebi.ac.uk/pdbe/docs/sifts/quick.html)
        "A summary of the UniProt to PDB mappings showing the UniProt accession
        followed by a semicolon-separated list of PDB four letter codes."
        '''
        if not filePath:
            url = "ftp://ftp.ebi.ac.uk/pub/databases/msd/sifts/flatfiles/csv/uniprot_pdb.csv.gz"
            filePath = SIFTS_unit.CONFIG['DOWNLOAD_FOLDER'] + 'uniprot_pdb_%s.csv.gz' % (time.strftime("%Y_%m_%d", time.localtime()))
            wget.download(url, out=filePath)
            g = gzip.GzipFile(mode="rb", fileobj=open(filePath, 'rb'))
            filePath = filePath[:-3]
            open(filePath, "wb").write(g.read())
            self.file_list.append(filePath)

        dfrm = pd.read_csv(filePath, sep=',', header=1)
        dfrm['PDB'] = dfrm.apply(lambda x: x['PDB'].upper(), axis=1)
        pdb_list = []
        if related_unp:
            dfrm = dfrm[dfrm['SP_PRIMARY'].isin(related_unp)]
        for i in dfrm.index:
            pdb_list.extend(dfrm.loc[i, 'PDB'].split(';'))
        if related_pdb:
            return {'pdb_set':set(pdb_list) & related_pdb, 'unp_set':set(dfrm['SP_PRIMARY'])}
        else:
            return {'pdb_set':set(pdb_list), 'unp_set':set(dfrm['SP_PRIMARY'])}

    def handle_SIFTS(self, sifts_filePath=False, outputPath=False):
        def addSiftsRange(in_df):
            group_info_col = SIFTS_unit.CONFIG['RAW_SIFTS_COLUMNS'][:3]
            range_info_col = SIFTS_unit.CONFIG['RAW_SIFTS_COLUMNS'][5:9]
            rangeSetER = SIFTS_unit.RangeSetER(group_info_col)
            in_df['rangeInfo'] = in_df.apply(lambda x: rangeSetER.check(
                tuple(x[i] for i in group_info_col),
                tuple(x[i] for i in range_info_col)),
                axis=1)
            return in_df.drop(columns=range_info_col).drop_duplicates(subset=group_info_col, keep='last')

        filePath = sifts_filePath or self.raw_SIFTS_filePath
        sifts_df = pd.read_csv(filePath, sep='\t', na_values=[SIFTS_unit.CONFIG['PDB_ID']])
        sifts_df.dropna(subset=[SIFTS_unit.CONFIG['PDB_ID']], inplace=True)
        sifts_df[SIFTS_unit.CONFIG['PDB_ID']] = sifts_df.apply(lambda x: x[SIFTS_unit.CONFIG['PDB_ID']].upper(), axis=1)

        new_sifts_df = addSiftsRange(sifts_df)
        new_sifts_df['sifts_pdb_range'] = new_sifts_df.apply(lambda x: x['rangeInfo'].split('|')[0], axis=1)
        new_sifts_df['sifts_unp_range'] = new_sifts_df.apply(lambda x: x['rangeInfo'].split('|')[1], axis=1)
        new_sifts_df.drop(columns=['rangeInfo'], inplace=True)
        new_sifts_df[new_sifts_df['is_canonical'] == True][['UniProt']].drop_duplicates().to_csv(SIFTS_unit.CONFIG['UNP_LIST_PATH'], sep='\t', index=False)

        self.file_o(outputPath, new_sifts_df)
        return new_sifts_df

    def deal_with_insertionDeletion_SIFTS(self, sifts_df=False, sifts_filePath=False, outputPath=False):
        def get_gap_list(li):
            gap_num = len(li) - 1
            gap_list = []
            for i in range(gap_num):
                gap_list.append(li[i+1][0] - li[i][1] - 1)
            return gap_list

        def get_ran_var(li_a, li_b):
            array_a = np.array([ran[1] - ran[0] + 1 for ran in li_a])
            array_b = np.array([ran[1] - ran[0] + 1 for ran in li_b])
            return (array_a - array_b).tolist()

        def add_tage_to_range(df, tage_name):
            # ADD TAGE FOR SIFTS
            df[tage_name] = 'Safe'
            # No Insertion But Deletion[Pure Deletion]
            df.loc[(df[(df['group_info'] == 1) & (df['sifts_unp_pdb_var'] > 0)].index), tage_name] = 'Deletion'
            # Insertion & No Deletion
            df.loc[df[
                (df['group_info'] != 1) &
                (df['var_0_count'] == df['group_info']) &
                (df['unp_GAP_0_count'] == (df['group_info'] -1))].index, tage_name] = 'Insertion'
            # Insertion & Deletion
            df.loc[df[
                (df['group_info'] != 1) &
                ((df['var_0_count'] != df['group_info']) |
                (df['unp_GAP_0_count'] != (df['group_info'] -1)))].index, tage_name] = 'Insertion & Deletion'

        dfrm = self.file_i(sifts_filePath, sifts_df, ('sifts_filePath', 'sifts_df'))
        dfrm['pdb_GAP_list'] = dfrm.apply(lambda x: json.dumps(get_gap_list(json.loads(x['sifts_pdb_range']))), axis=1)
        dfrm['unp_GAP_list'] = dfrm.apply(lambda x: json.dumps(get_gap_list(json.loads(x['sifts_unp_range']))), axis=1)
        dfrm['var_list'] = dfrm.apply(lambda x: json.dumps(get_ran_var(json.loads(x['sifts_unp_range']), json.loads(x['sifts_pdb_range']))), axis=1)
        dfrm['delete'] = dfrm.apply(lambda x: x['var_list'].find('-') != -1, axis=1)
        dfrm['delete'] = dfrm.apply(lambda x: True if x['unp_GAP_list'].find('-') != -1 else x['delete'], axis=1)
        dfrm['var_0_count'] = dfrm.apply(lambda x: json.loads(x['var_list']).count(0), axis=1)
        dfrm['unp_GAP_0_count'] = dfrm.apply(lambda x: json.loads(x['unp_GAP_list']).count(0), axis=1)
        dfrm['group_info'] = dfrm.apply(lambda x: len(json.loads(x['sifts_pdb_range'])), axis=1)
        dfrm['sifts_unp_pdb_var'] = dfrm.apply(lambda x: json.loads(x['var_list'])[0], axis=1)

        add_tage_to_range(dfrm, tage_name='sifts_range_tage')
        self.file_o(outputPath, dfrm)
        return dfrm

    def add_mmcif_info_SIFTS(self, sifts_df=False, sifts_filePath=False, mmcif_df=False, mmcif_filePath=False, outputPath=False):
        sifts_dfrm = self.file_i(sifts_filePath, sifts_df, ('sifts_filePath', 'sifts_df'))
        mmcif_dfrm = self.file_i(mmcif_filePath, mmcif_df, ('mmcif_filePath', 'mmcif_df'))
        dfrm = pd.merge(sifts_dfrm, mmcif_dfrm, on=['pdb_id', 'chain_id'], how='left')
        self.file_o(outputPath, dfrm)
        return dfrm

    def add_unp_len_SIFTS(self, sifts_df=False, sifts_filePath=False, unpLen_df=False, unpLen_filePath=False, outputPath=False, sep=','):
        # unpLen_filePath = '/home/zzf/Work/StructDDG_0427/data/Mapping_Pipeline/sifts_unp_len_list.csv'
        sifts_dfrm = self.file_i(sifts_filePath, sifts_df, ('sifts_filePath', 'sifts_df'))
        sifts_dfrm['Entry'] = sifts_dfrm.apply(lambda x: x['UniProt'].split('-')[0], axis=1)
        unpLen_dfrm = self.file_i(unpLen_filePath, unpLen_df, ('unpLen_filePath', 'unpLen_df'), sep=sep)
        unpLen_dfrm.drop(columns=['yourlist'], inplace=True)
        unpLen_dfrm.rename(columns={'Length':'UNP_len'}, inplace=True)
        dfrm = pd.merge(sifts_dfrm, unpLen_dfrm, on=['Entry'], how='left')
        self.file_o(outputPath, dfrm)
        return dfrm

    def score_SIFTS(self, sifts_df=False, sifts_filePath=False, outputPath=False):
        def getMappedRange(pdb_range, mis_range):
            try:
                pdb_range = json.loads(pdb_range)
            except Exception:
                return np.nan
            try:
                mis_range = json.loads(mis_range)
            except Exception:
                mis_range = False
            pdb_range_set = set()
            mis_range_set = set()
            for ran in pdb_range:
                pdb_range_set = pdb_range_set | set(range(ran[0], ran[1]+1))
            if mis_range:
                for ran in mis_range:
                    mis_range_set = mis_range_set | set(range(ran[0], ran[1]+1))
            return MMCIF_unit.getInterval(pdb_range_set - mis_range_set)

        def getMappedOut(pdb_range, seqres_len):
            try:
                li = json.loads(pdb_range)
            except Exception:
                return np.nan
            out_head = li[0][0]-1
            out_tail = seqres_len - li[-1][-1]
            if out_head <= 5:
                out_head = 0
            else:
                out_head -= 5
            if out_tail <= 5:
                out_tail = 0
            else:
                out_tail -= 5
            return out_head + out_tail

        def getHeadTailMisNum(head, tail, li):
            try:
                mis_li = json.loads(li)
            except Exception:
                return 0
            if not isinstance(head, float) and not isinstance(tail, float): # WARNING
                return len(list(filter(lambda x: (x >= head+5) & (x <= tail-5), mis_li)))
            else:
                return 0

        def getRangeLen(li):
            length = 0
            for ran in li:
                length += len(range(ran[0], ran[1]+1))
            return length

        def get_weight():
            ele_list = SIFTS_unit.CONFIG['ELE_LIST']
            ini_list = SIFTS_unit.CONFIG['INI_LIST']
            com_list = list(combinations(ele_list, 2))
            com_dict = {}
            count = 0
            for i in com_list:
                com_dict[i] = ini_list[count]
                count += 1
            mat = pd.DataFrame([[1]*len(ele_list)]*len(ele_list))
            mat.columns = ele_list
            mat.index = ele_list
            for i in ele_list:
                for j in ele_list:
                    mat[i][j] = com_dict.get((i,j), 1)
            x = mat.values
            y = 1/x
            value, vector = np.linalg.eig(np.triu(x.T).T + np.triu(y.T) - np.diag(x.diagonal()))
            max_val = np.max(value)
            index = list(value).index(max_val)
            max_vector = vector[:,index]
            select_vector = np.real(max_vector/np.linalg.norm(max_vector,ord=1))
            return -select_vector*50

        def calScore(df, colName, ele_list, weight):
            def cal(data, ele_list, weight):
                score = data[ele_list[0]] * weight[0]
                count = 1
                for ele in ele_list[1:]:
                    score -= data[ele] * weight[count]
                    count += 1
                return score

            df[colName] = df.apply(lambda x: cal(x, ele_list, weight) / x['UNP_len'] if not isinstance(x['sifts_pdb_range'], float) else np.nan, axis=1)

        def getUsefulChainNum(s, cutoff):
            # useful = []
            count = 0
            li = s.split(',')
            for i in li:
                temp = i.split(':')
                if int(temp[1]) > cutoff:
                    # useful.append(temp[0])
                    count += 1
            # return useful
            return count

        sifts_dfrm = self.file_i(sifts_filePath, sifts_df, ('sifts_filePath', 'sifts_df'))
        # Metal ligand Info
        sifts_dfrm['metal_list'] = sifts_dfrm.apply(lambda x: list(map(int, x['ligand_position_in_seqres'].split(
            ';'))) if not isinstance(x['ligand_position_in_seqres'], float) else np.nan, axis=1)
        # MappedRange Info
        sifts_dfrm['pdb_mapped_range'] = sifts_dfrm.apply(
            lambda x: getMappedRange(x['sifts_pdb_range'], x['mis_range']), axis=1)
        sifts_dfrm['pdb_mappedRange_head'] = sifts_dfrm.apply(
            lambda x: x['pdb_mapped_range'][0][0] if not isinstance(x['pdb_mapped_range'], float) else np.nan, axis=1)
        sifts_dfrm['pdb_mappedRange_tail'] = sifts_dfrm.apply(
            lambda x: x['pdb_mapped_range'][-1][-1] if not isinstance(x['pdb_mapped_range'], float) else np.nan, axis=1)
        # MappedOut Info
        sifts_dfrm['mappedOut'] = sifts_dfrm.apply(
            lambda x: getMappedOut(x['sifts_pdb_range'], x['seqres_len']), axis=1)
        # Del HeadTail Missing
        sifts_dfrm['delHT_MissingNum'] = sifts_dfrm.apply(lambda x: getHeadTailMisNum(
            x['pdb_mappedRange_head'], x['pdb_mappedRange_tail'], x['mis_index']), axis=1)
        # Get MappedRange
        sifts_dfrm['pdb_mapped_range_len'] = sifts_dfrm.apply(lambda x: getRangeLen(
            x['pdb_mapped_range']) if not isinstance(x['pdb_mapped_range'], float) else 0, axis=1)
        # Get BasicScore
        sifts_dfrm['if_1'] = sifts_dfrm.apply(
            lambda x: sum(json.loads(x['pdb_GAP_list'])) + sum(
            json.loads(x['unp_GAP_list'])) + sum(
            json.loads(x['var_list'])), axis=1)
        sifts_dfrm['if_2'] = sifts_dfrm.apply(
            lambda x: len(x['Modification_num']) + x['mutation_num'] if not isinstance(x['Modification_num'], float) else x['mutation_num'], axis=1)
        sifts_dfrm['metal_count'] = sifts_dfrm.apply(lambda x: len(x['metal_list']) if not isinstance(x['metal_list'], float) else 0, axis=1)

        select_vector = get_weight()
        calScore(sifts_dfrm, SIFTS_unit.CONFIG['SCORE_COL'], SIFTS_unit.CONFIG['ELE_LIST'], select_vector)

        # Add resolution-related socre
        sifts_dfrm['ne_resolution_score'] = sifts_dfrm.apply(
            lambda x: -x['resolution_score'], axis=1)

        # Find Useful chain
        sifts_dfrm['pdb_SIFTS_useful_chain_num'] = np.nan
        for i,j in sifts_dfrm.groupby([SIFTS_unit.CONFIG['PDB_ID']]):
            sifts_dfrm.loc[j.index, 'pdb_SIFTS_useful_chain_num'] = len(set(j[j['coordinates_len'] > 20]['chain_id']))
        # ---
        sifts_dfrm['pdb_useful_chain_num'] = sifts_dfrm.apply(
            lambda x: getUsefulChainNum(x['protein_chain_and_length'], 20), axis=1)


        self.file_o(outputPath, sifts_dfrm)
        return sifts_dfrm

    def get_seg_info_from_uniprot_segments_file(self, filePath=False, outputPath=False):
        '''
        Reference: (http://www.ebi.ac.uk/pdbe/docs/sifts/quick.html)
        "A summary of the UniProt to PDBe residue level mapping
        (observed residues only), showing the start and end residues of the
        mapping using SEQRES, PDB sequence and UniProt numbering."
        '''
        def add_seg_range(in_df):
            group_info_col = ['PDB', 'CHAIN', 'SP_PRIMARY']
            range_info_col = ['RES_BEG', 'RES_END', 'SP_BEG', 'SP_END']
            rangeSetER = SIFTS_unit.RangeSetER(group_info_col)
            in_df['rangeInfo'] = in_df.apply(lambda x: rangeSetER.check(
                tuple(x[i] for i in group_info_col),
                tuple(x[i] for i in range_info_col)),
                axis=1)
            return in_df.drop(columns=range_info_col+['PDB_BEG','PDB_END']).drop_duplicates(subset=group_info_col, keep='last')

        if not filePath:
            url = "ftp://ftp.ebi.ac.uk/pub/databases/msd/sifts/flatfiles/csv/uniprot_segments_observed.csv.gz"
            filePath = SIFTS_unit.CONFIG['DOWNLOAD_FOLDER'] + 'uniprot_segments_observed_%s.csv.gz' % (time.strftime("%Y_%m_%d", time.localtime()))
            wget.download(url, out=filePath)
            g = gzip.GzipFile(mode="rb", fileobj=open(filePath, 'rb'))
            filePath = filePath[:-3]
            open(filePath, "wb").write(g.read())
            self.file_list.append(filePath)

        dfrm = pd.read_csv(filePath, sep=',', header=1)
        new_dfrm = add_seg_range(dfrm)
        new_dfrm['seg_pdb_range'] = new_dfrm.apply(lambda x: x['rangeInfo'].split('|')[0], axis=1)
        new_dfrm['seg_unp_range'] = new_dfrm.apply(lambda x: x['rangeInfo'].split('|')[1], axis=1)
        new_dfrm.drop(columns=['rangeInfo'], inplace=True)
        new_dfrm.columns = [SIFTS_unit.CONFIG['SEG_SIFTS_MAP'].get(i, i) for i in new_dfrm.columns]
        new_dfrm[SIFTS_unit.CONFIG['PDB_ID']] = new_dfrm.apply(lambda x: x[SIFTS_unit.CONFIG['PDB_ID']].upper(), axis=1)
        self.file_o(outputPath, new_dfrm)
        return new_dfrm

    def add_seg_info_to_SIFTS(self, sifts_df=False, sifts_filePath=False, seg_df=False, seg_filePath=False, outputPath=False):
        sifts_dfrm = self.file_i(sifts_filePath, sifts_df, ('sifts_filePath', 'sifts_df'))
        seg_dfrm = self.file_i(seg_filePath, seg_df, ('seg_filePath', 'seg_df'))
        dfrm = pd.merge(sifts_dfrm, seg_dfrm, on=['pdb_id', 'chain_id', 'UniProt'], how='left')
        self.file_o(outputPath, dfrm)
        return dfrm

    def select_PDB_SIFTS(self, groupby_list, sifts_df=False, sifts_filePath=False, outputPath=False):
        sifts_dfrm = self.file_i(sifts_filePath, sifts_df, ('sifts_filePath', 'sifts_df'))
        select_col = SIFTS_unit.CONFIG['PDB_SELECT_COL']
        rank_col = SIFTS_unit.CONFIG['PDB_RANK_COL']
        rank_list = SIFTS_unit.CONFIG['PDB_RANK_LIST']
        rank_format = SIFTS_unit.CONFIG['PDB_RANK_FORMAT']
        # groupby_list = SIFTS_unit.CONFIG['PDB_INIT_GROUPBY_LIST']
        range_name = SIFTS_unit.CONFIG['PDB_SELECT_RANGE_NAME']

        sifts_dfrm[select_col] = False
        sifts_dfrm[rank_col] = np.nan

        for i, j in sifts_dfrm.groupby(groupby_list):
            SIFTS_unit.selectChain(
                j, sifts_dfrm, rank_list, rank_col, rank_format, range_name, select_col, 0.3, 0.2)

        self.file_o(outputPath, sifts_dfrm)
        return sifts_dfrm

    def find_mo_SIFTS(self, groupby_list, sifts_df=False, sifts_filePath=False, outputPath=False):
        sifts_dfrm = self.file_i(sifts_filePath, sifts_df, ('sifts_filePath', 'sifts_df'))
        constraint_dict = ConstraintDict(
            {
                'only_contains_unk_in_chain_pdb': ('no', 'eq'),
                'UNK_ALL_IN_CHAIN': ('no', 'eq'),
                'pdb_contain_chain_type': ('protein', 'eq'),
                'pdb_SIFTS_useful_chain_num': (0, 'gt'),
                'coordinates_len': (20, 'gt'),
                'pdb_mapped_range_len': (20, 'gt'),
                'delete': (False, 'eq'),
                'identity': (0.9, 'ge'),
                'is_canonical': (True, 'eq'),

            }
        )
        constraint_dict.setMutalDict(
            {
                'pdb_useful_chain_num': (1, 'eq'),
                'pdb_SIFTS_useful_chain_num': (1, 'eq'),
                'resolution_score': (3, 'le')
            }
        )

        mo_fakeHoHe_df = ConstraintDict.addConstraintToDf(sifts_dfrm, constraint_dict)
        file_o(outputPath, mo_fakeHoHe_df)
        return mo_fakeHoHe_df

    def map_muta_to_PDB_SIFTS(self, sifts_df=False, sifts_filePath=False, outputPath=False):
        def getUniprotMutaSite(pdb_position, coor_list, sifts_unp_range, sifts_pdb_range):
            position = pdb_position[1:-1]
            seqresMutaSite = coor_list.index(position)+1
            try:
                pdb_range = json.loads(sifts_pdb_range)
                unp_range = json.loads(sifts_unp_range)
            except Exception:
                print('getUniprotMutaSite(): Fail to load range')
                return 0
            pdb_li = []
            unp_li = []
            for ran in pdb_range:
                pdb_li.extend(list(range(ran[0], ran[1]+1)))
            for ran in unp_range:
                unp_li.extend(list(range(ran[0], ran[1]+1)))
            try:
                r = unp_li[pdb_li.index(int(seqresMutaSite))]
                return r
            except Exception:
                print('getUniprotMutaSite(): Fail to get Site')
                return 0

        def mapMutaFromPDBToUniprot(dfrm, groupby_list, pdb_muta_col, unp_muta_col):
            dfrm[unp_muta_col] = np.nan # 'Mutation_Uniprot'
            for _, groupData in dfrm.groupby(groupby_list): # ['pdb_id', 'chain_id', 'iso_id']
                coor_list = groupData.loc[groupData.index[0], 'pdb_ins_position'].split(';') # pdb_ins_seqres_position
                sifts_unp_range = groupData.loc[groupData.index[0], 'sifts_unp_range']
                sifts_pdb_range = groupData.loc[groupData.index[0], 'sifts_pdb_range']
                dfrm.loc[groupData.index, unp_muta_col] = groupData.apply(
                    lambda x: getUniprotMutaSite(
                        x[pdb_muta_col], coor_list, sifts_unp_range, sifts_pdb_range
                        ), axis=1)

        sifts_dfrm = self.file_i(sifts_filePath, sifts_df, ('sifts_filePath', 'sifts_df'))
        mapMutaFromPDBToUniprot(sifts_dfrm)
        file_o(outputPath, sifts_dfrm)
        return sifts_dfrm


class MMCIF_unit(Unit):
    CONFIG = {
        'PDB_ID': 'pdb_id',
        'MMCIF_OLD_FOLDER': ['/data1/suntt/process0606/cgc_mmcif_file/', '/data1/suntt/CanDriver/Data/PDB_cgc/cgc_mmcif_file/', '/data1/suntt/CanDriver/Data/PDB_NEW/mmcif_file/'],
        'MMCIF_FOLDER': '../../data/Mapping_Pipeline/mmcif_file/',
        'OUTPUT_FOLDER': '../../data/Mapping_Pipeline/output_files/',
        'HEADERS': {'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/64.0.3282.140 Safari/537.36 Edge/17.17134'},
        'CHAIN_TYPE_FILE_LIST': ('cgc_pdb_chain_type.txt', 'cgc_pdb_chain_type_extra', 'cgc_chain_type_error_pdb.txt', 'output/cgc_pdb_chain_type_all_new.txt'),
        'ATOM_SEQRES_FILE_LIST': ('cgc_pdb_atom_seqres.txt', 'cgc_seqres_atom_error_pdb.txt', 'output/cgc_pdb_atom_seqres_add_chain_type_new.txt', 'output/cgc_pdb_coordinates_site_range.txt'),
        'LIGAND_FILE_LIST': ('ligand_info0605.txt', 'ligand_info_extra0605.txt', 'ligand_info_error.txt', 'output/ligand_info_final1.txt', 'output/ligand_info_final2.txt'),
        'CHAIN_AND_SEQRES_FILE_LIST': ('output/cgc_pdb_atom_seqres_protein_chain_info.txt', 'cgc_protein_chain_id_in_pdb.txt', 'output/cgc_pdb_atom_seqres_info_integration.txt'),
        'ADD_MODIFICATION_FILE': 'output/cgc_pdb_atom_seqres_info_integration_new.txt',
        'ADD_MISSING_FILE': ('output/cgc_information_statistics1.txt', 'output/cgc_pdb_atom_seqres_info_integration_new_add_coor_site.txt'),
        'PDB_MUTATION_FILE': 'output/pdb_mutation_info.txt',
        'RESOLUTION_FILE': ('output/resoluton_error.txt', 'output/pdb_resoluton_info.txt'),
        'YEAR_INFO_1_FILE': ('../../data/Mapping_Pipeline/output_files/pdb_date_info_newest.txt', '../../data/Mapping_Pipeline/output_files/pdb_date_error_newest.txt'),
        'YEAR_INFO_2_FILE': ('../../data/Mapping_Pipeline/output_files/pdb_date_supp_newest.txt', '../../data/Mapping_Pipeline/output_files/pdb_date_error_supp_newest.txt'),
        'YEAR_INFO_ALL': '../../data/Mapping_Pipeline/output_files/pdb_date_info_newest_all.txt',
        'FINAL_FILE': ('output/cgc_pdb_atom_seqres_info_integration_final.txt', 'PDB_cgc/output/cgc_pdb_atom_seqres_info_integration_final0614.txt'),
        'MMICF_USECOLS': ['pdb_id', 'chain_id', 'seqres_len', 'coordinates_len', 'Modification_position', 'ligand_position_in_seqres', 'mis_range', 'mis_index', 'Modification_num', 'mutation_num'],

        }

    def set_output_folder(self, path):
        self.CONFIG['OUTPUT_FOLDER'] = path

    def download_cif_file(pdbId, path):
        url = 'https://files.rcsb.org/view/%s.cif' % pdbId
        html = request.urlopen(url).read()
        html = html.decode('utf-8')
        with open(path, 'w') as fw:
            fw.write(html)
            time.sleep(2)

    def get_mmcif_file_path(self, pdbId, download=False):
        print('get_mmcif_file_path(): Working on [%s]' % pdbId)
        new_path  = '%s%s.cif' % (self.CONFIG['MMCIF_FOLDER'], pdbId)

        for path in self.CONFIG['MMCIF_OLD_FOLDER']:
            old_path = '%s%s.cif' % (path, pdbId)
            if os.path.exists(old_path):
                return old_path

        if os.path.exists(new_path):
            return new_path
        else:
            if download:
                MMCIF_unit.download_cif_file(pdbId, new_path)
                return False
            else:
                return new_path


    '''
    def download_mmcif_file(self):
        # 暂时不用(zzf)
        @retry(stop_max_attempt_number=3, wait_fixed=1000)
        pool = Pool(processes=20)
        pool.map(download_cif_file, self.pdb_list)
    '''

    def extract_chain_type_info(self):
        chain_type_file1, chain_type_file2, chain_type_error, chain_type_file_all = MMCIF_unit.CONFIG['CHAIN_TYPE_FILE_LIST']
        outpath = self.CONFIG['OUTPUT_FOLDER']

        demo_dict_df_list = []
        fw = open(outpath + chain_type_file2, 'w')
        error_pdb_file = open(outpath + chain_type_error, 'w')
        for pdbId in self.pdb_list:
            # pdbFileSavePath = '%s%s.cif' % (MMCIF_unit.CONFIG['MMCIF_FOLDER'], pdbId)
            pdbFileSavePath = self.get_mmcif_file_path(pdbId, True)
            if not pdbFileSavePath:
                continue

            mmcif_dict = MMCIF2Dict(pdbFileSavePath)
            demo_dict = {}
            index = ['_entity_poly.type', '_entity_poly.pdbx_strand_id']
            try:
                for i in index:
                    demo_dict[i] = mmcif_dict[i]
                df = pd.DataFrame(demo_dict)
                df['pdb_id'] = pdbId
                demo_dict_df_list.append(df)
            except:
                try:
                    fw.write('%s\t%s\t%s\n' % (
                    pdbId, mmcif_dict['_entity_poly.type'], mmcif_dict['_entity_poly.pdbx_strand_id']))
                except:
                    error_pdb_file.write(pdbId + '\n')

        demo_df = pd.concat(demo_dict_df_list)
        demo_df.to_csv(outpath + chain_type_file1, sep='\t', index=False)
        fw.close()
        error_pdb_file.close()
        # 将chain_type的信息合并到一起
        info = pd.read_csv(outpath + chain_type_file1, sep='\t', dtype=str)
        info1 = pd.read_csv(outpath + chain_type_file2, sep='\t', dtype=str,
                            names=['pdb_id', '_entity_poly.type', '_entity_poly.pdbx_strand_id'])
        info2 = pd.concat([info, info1], axis=0)
        info2.rename(columns={'_entity_poly.pdbx_strand_id': 'chain_id', '_entity_poly.type': 'chain_type_details'},
                     inplace=True)
        info2['chain_type'] = info2['chain_type_details'].replace('polypeptide(L)', 'protein').replace('polypeptide(D)',
                                                                                                       'protein').replace(
            'polydeoxyribonucleotide', 'DNA').replace('polyribonucleotide', 'RNA').replace(
            'polydeoxyribonucleotide/polyribonucleotide hybrid', 'RNA+DNA')
        # info2.to_csv(outpath+'PDB_cgc/cgc_pdb_chain_type_all.txt',sep='\t',index=False)

        #重新设置索引号，避免同一索引对应不同行，因为两个数据concat时各自文件的索引仍是之前的
        info2.index = range(len(info2))

        #由于重新设置了索引不会造成混淆，所以可以使用以下方法，比较快
        result = info2.drop('chain_id', axis=1).join(
            info2['chain_id'].str.split(',', expand=True).stack().reset_index(level=1, drop=True).rename('chain_id_new'))
        info3 = result[['pdb_id', 'chain_type']].drop_duplicates()
        info4 = info3.sort_values(by=['chain_type']).groupby(['pdb_id'], as_index=False).agg(lambda x: ','.join(x))
        info4.rename(columns={'chain_type': 'pdb_contain_chain_type'}, inplace=True)
        info5 = pd.merge(result, info4, on=['pdb_id'], how='left')
        info5.to_csv(outpath + chain_type_file_all, sep='\t', index=False)

    def extract_seqres_and_atom_info(self):
        atom_seqres_file, atom_seqres_error, atom_seqres_chain_type_oringnal, coordinates_file = MMCIF_unit.CONFIG['ATOM_SEQRES_FILE_LIST']
        outpath = self.CONFIG['OUTPUT_FOLDER']
        chain_type_file_all = MMCIF_unit.CONFIG['CHAIN_TYPE_FILE_LIST'][3]

        demo_dict_df_list = []
        error_pdb_file = open(outpath + atom_seqres_error, 'w')
        for pdbId in self.pdb_list:
            # pdbFileSavePath = '%s%s.cif' % (MMCIF_unit.CONFIG['MMCIF_FOLDER'], pdbId)
            pdbFileSavePath = self.get_mmcif_file_path(pdbId)
            if not pdbFileSavePath:
                continue

            mmcif_dict = MMCIF2Dict(pdbFileSavePath)
            demo_dict = {}
            index = ['_pdbx_poly_seq_scheme.mon_id', '_pdbx_poly_seq_scheme.ndb_seq_num',
                     '_pdbx_poly_seq_scheme.pdb_seq_num', '_pdbx_poly_seq_scheme.auth_seq_num',
                     '_pdbx_poly_seq_scheme.pdb_mon_id', '_pdbx_poly_seq_scheme.auth_mon_id',
                     '_pdbx_poly_seq_scheme.pdb_strand_id', '_pdbx_poly_seq_scheme.pdb_ins_code']
            try:
                for i in index:
                    demo_dict[i] = mmcif_dict[i]
                df = pd.DataFrame(demo_dict)
                df['pdb_id'] = pdbId
                demo_dict_df_list.append(df)
            except:
                error_pdb_file.write(pdbId + '\n')

        demo_df1 = pd.concat(demo_dict_df_list)
        demo_df1.to_csv(outpath + atom_seqres_file, sep='\t', index=False)
        error_pdb_file.close()


        # 将chain_type信息加入seqres和atom部分
        file1 = pd.read_csv(outpath + atom_seqres_file, sep='\t', dtype=str)
        file2 = pd.read_csv(outpath + chain_type_file_all, sep='\t', dtype=str)
        file2.rename(columns={'chain_id_new': 'chain_id'}, inplace=True) # ?
        file3 = pd.merge(file1, file2, left_on=['pdb_id', '_pdbx_poly_seq_scheme.pdb_strand_id'],
                         right_on=['pdb_id', 'chain_id'], how='left')
        # file3.to_csv(outpath+'PDB_cgc/cgc_pdb_atom_seqres_add_chain_type.txt',sep='\t',index=False)
        file3.rename(columns={'_pdbx_poly_seq_scheme.mon_id': 'SEQRES', '_pdbx_poly_seq_scheme.pdb_mon_id': 'Coordinates',
                              '_pdbx_poly_seq_scheme.ndb_seq_num': 'pdb_index',
                              '_pdbx_poly_seq_scheme.pdb_seq_num': 'position_in_seqres',
                              '_pdbx_poly_seq_scheme.auth_seq_num': 'position_in_coordinates',
                              '_pdbx_poly_seq_scheme.pdb_ins_code': 'inside_code'}, inplace=True)
        file4 = file3.drop(['_pdbx_poly_seq_scheme.auth_mon_id', '_pdbx_poly_seq_scheme.pdb_strand_id'], axis=1)
        file4.to_csv(outpath + atom_seqres_chain_type_oringnal, sep='\t', index=False)
        # 加入coordinates_start和coordinates_end信息
        coordinates_range = file4[file4['pdb_contain_chain_type'].notna() & file4['pdb_contain_chain_type'].str.contains('protein')]
        coordinates_range['pdb_ins_position'] = coordinates_range['position_in_seqres'] + coordinates_range['inside_code']
        coordinates_range['pdb_ins_position'] = coordinates_range['pdb_ins_position'].str.replace('.', '')
        coordinates_range1 = coordinates_range.groupby(['pdb_id', 'chain_id'], as_index=False)['pdb_ins_position'].agg(
            lambda x: ';'.join(x))
        coordinates_range1.to_csv(outpath + coordinates_file, sep='\t', index=False)

    def extract_pdb_ligand_info(self):
        outpath = self.CONFIG['OUTPUT_FOLDER']
        ligand_file1, ligand_file2, ligand_file_error, ligand_file_final1, ligand_file_final2 = MMCIF_unit.CONFIG['LIGAND_FILE_LIST']
        atom_seqres_chain_type_oringnal = MMCIF_unit.CONFIG['ATOM_SEQRES_FILE_LIST'][2]

        demo_dict_df_list = []
        fw = open(outpath + ligand_file2,'w')
        fp = open(outpath + ligand_file_error,'w')
        for pdbId in self.pdb_list:
            # pdbFileSavePath = '%s%s.cif' % (MMCIF_unit.CONFIG['MMCIF_FOLDER'], pdbId)
            pdbFileSavePath = self.get_mmcif_file_path(pdbId)
            if not pdbFileSavePath:
                continue

            mmcif_dict = MMCIF2Dict(pdbFileSavePath)
            demo_dict = {}
            index = ['_struct_conn.conn_type_id','_struct_conn.ptnr1_auth_asym_id','_struct_conn.ptnr1_auth_comp_id','_struct_conn.ptnr1_auth_seq_id',
                     '_struct_conn.ptnr2_auth_asym_id','_struct_conn.ptnr2_auth_comp_id','_struct_conn.ptnr2_auth_seq_id']
            try:
                for i in index:
                    demo_dict[i] = mmcif_dict[i]
                df = pd.DataFrame(demo_dict)
                df['pdb_id'] = pdbId
                demo_dict_df_list.append(df)
            except:
                try:
                    fw.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n'%(pdbId,mmcif_dict['_struct_conn.conn_type_id'],mmcif_dict['_struct_conn.ptnr1_auth_asym_id'],mmcif_dict['_struct_conn.ptnr1_auth_comp_id'],
                                             mmcif_dict['_struct_conn.ptnr1_auth_seq_id'],mmcif_dict['_struct_conn.ptnr2_auth_asym_id'],mmcif_dict['_struct_conn.ptnr2_auth_comp_id'],
                                             mmcif_dict['_struct_conn.ptnr2_auth_seq_id']))
                except:
                    fp.write(pdbId+'\n')
        demo_df = pd.concat(demo_dict_df_list)
        demo_df.to_csv(outpath + ligand_file1,sep='\t',index=False)
        fw.close()
        fp.close()

        def ligand_count(ligand_ptnr_seq_id):
            a = len(ligand_ptnr_seq_id.split(';'))
            return a

        def metal_check(connection_type):
            if pd.isnull(connection_type):
                return '0'
            else:
                if connection_type == 'metalc':
                    return '1'
                else:
                    return '0'

        # atom_seqres_chain_type_oringnal = 'output/cgc_pdb_atom_seqres_add_chain_type_new.txt'
        ligand_info1 = pd.read_csv(outpath + ligand_file1, sep='\t', dtype=str, keep_default_na=False)
        ligand_info2 = pd.read_csv(outpath + ligand_file2, sep='\t', dtype=str, keep_default_na=False,
                                   names=['pdb_id', '_struct_conn.conn_type_id', '_struct_conn.ptnr1_auth_asym_id',
                                          '_struct_conn.ptnr1_auth_comp_id', '_struct_conn.ptnr1_auth_seq_id',
                                          '_struct_conn.ptnr2_auth_asym_id', '_struct_conn.ptnr2_auth_comp_id',
                                          '_struct_conn.ptnr2_auth_seq_id'])
        ligand_info_all = pd.concat([ligand_info1, ligand_info2], axis=0)
        metal_ligand = ['ZN', 'MG', 'CA', 'FE', 'NA', 'MN', 'K', 'NI', 'CU', 'CO', 'CD', 'HG', 'PT', 'MO', 'BE', 'AL', 'BA',
                        'RU', 'SR', 'V', 'CS', 'W', 'AU', 'YB', 'LI', 'GD', 'PB', 'Y', 'TL', 'IR', 'RB', 'SM', 'AG',
                        'OS', 'PR', 'PD', 'EU', 'RH', 'RE', 'TB', 'TA', 'LU', 'HO', 'CR', 'GA', 'LA', 'SN', 'SB', 'CE',
                        'ZR',
                        'ER', 'TH', 'TI', 'IN', 'HF', 'SC', 'DY', 'BI', 'PA', 'PU', 'AM', 'CM', 'CF', 'GE', 'NB', 'TC',
                        'ND',
                        'PM', 'TM', 'PO', 'FR', 'RA', 'AC', 'NP', 'BK', 'ES', 'FM', 'MD', 'NO', 'LR', 'RF', 'DB', 'SG']

        ligand_info_all1 = ligand_info_all[(ligand_info_all['_struct_conn.conn_type_id'] == 'metalc') & (
            ligand_info_all['_struct_conn.ptnr1_auth_comp_id'].isin(metal_ligand))]
        ligand_info_all1.rename(
            columns={'_struct_conn.conn_type_id': 'connection_type', '_struct_conn.ptnr1_auth_asym_id': 'ligand_chain',
                     '_struct_conn.ptnr1_auth_comp_id': 'ligand_comp',
                     '_struct_conn.ptnr1_auth_seq_id': 'ligand_seq_id', '_struct_conn.ptnr2_auth_asym_id': 'chain_id',
                     '_struct_conn.ptnr2_auth_comp_id': 'ligand_ptnr_comp',
                     '_struct_conn.ptnr2_auth_seq_id': 'position_in_seqres'}, inplace=True)

        ligand_info_all2 = ligand_info_all[(ligand_info_all['_struct_conn.conn_type_id'] == 'metalc') & (
            ligand_info_all['_struct_conn.ptnr2_auth_comp_id'].isin(metal_ligand))]

        ligand_info_all2.rename(
            columns={'_struct_conn.conn_type_id': 'connection_type', '_struct_conn.ptnr2_auth_asym_id': 'ligand_chain',
                     '_struct_conn.ptnr2_auth_comp_id': 'ligand_comp',
                     '_struct_conn.ptnr2_auth_seq_id': 'ligand_seq_id', '_struct_conn.ptnr1_auth_asym_id': 'chain_id',
                     '_struct_conn.ptnr1_auth_comp_id': 'ligand_ptnr_comp',
                     '_struct_conn.ptnr1_auth_seq_id': 'position_in_seqres'}, inplace=True)
        ligand_info_all3 = pd.concat([ligand_info_all1, ligand_info_all2], axis=0)
        ligand_info_all3.reset_index(drop=True)
        ligand_info_all3['ismetal'] = ligand_info_all3.apply(lambda x: metal_check(x.connection_type), axis=1)
        ligand_info_all4 = ligand_info_all3.drop(['ligand_chain', 'ligand_seq_id'], axis=1)
        ligand_info_all4 = ligand_info_all4.drop_duplicates()

        data_index_position = pd.read_csv(outpath + atom_seqres_chain_type_oringnal, sep='\t', dtype=str)
        data_index_position2 = data_index_position[
            ['pdb_id', 'chain_id', 'pdb_index', 'position_in_seqres']].drop_duplicates()
        infomerge = pd.merge(ligand_info_all4, data_index_position2, how='left',
                             on=['pdb_id', 'chain_id', 'position_in_seqres'])
        infomerge1 = infomerge[infomerge['pdb_index'].notna()]
        infomerge1.to_csv(outpath + ligand_file_final1, sep='\t', index=False)
        infomerge2 = infomerge1.groupby(['pdb_id', 'chain_id'], as_index=False)[
            'ligand_comp', 'ligand_ptnr_comp', 'position_in_seqres', 'pdb_index', 'ismetal'].agg(lambda x: ';'.join(x))
        infomerge2['ligand_count'] = infomerge2.apply(lambda x: ligand_count(x.position_in_seqres), axis=1)
        infomerge2.to_csv(outpath + ligand_file_final2, sep='\t', index=False)

    def deal_with_chain_and_seqres_atom(self):
        def get_modification(m, k):
            m1 = str(m).replace('?', '')
            yes = [i + 1 for i, v in enumerate(m1) if v == 'X']
            if yes != []:
                length = int(len(yes))
                if k == '1':
                    yes1 = str(yes).replace('[', '').replace(']', '').replace(' ','')
                    return yes1
                elif k == '2':
                    return length

        def get_modification_seqres_index(m):
            #     print m
            yes = [i + 1 for i, v in enumerate(m) if v == 'X']
            if yes != []:
                yes1 = str(yes).replace('[', '').replace(']', '').replace(' ','')
                return yes1

        atom_seqres_protein_chain, cgc_protein_chain_id, integration_file = MMCIF_unit.CONFIG['CHAIN_AND_SEQRES_FILE_LIST']
        outpath = self.CONFIG['OUTPUT_FOLDER']
        atom_seqres_chain_type_oringnal = MMCIF_unit.CONFIG['ATOM_SEQRES_FILE_LIST'][2]
        multiToOne = MMCIF_unit.MultiToOne()

        f = pd.read_csv(outpath + atom_seqres_chain_type_oringnal, sep='\t', dtype=str)
        f['SEQRES'] = f.apply(lambda x: multiToOne.multi_letter_convert_to_one_letter(x.SEQRES), axis=1)
        f['Coordinates'] = f.apply(lambda x: multiToOne.multi_letter_convert_to_one_letter(x.Coordinates), axis=1)
        ##以下仅针对蛋白链
        f1 = f[f['chain_type'] == 'protein']# .reset_index(drop=True)
        f1.to_csv(outpath + atom_seqres_protein_chain, sep='\t', index=False)
        # 将seqres信息放入一行
        f2 = f1[['pdb_id', 'chain_id', 'SEQRES', 'inside_code']]  # 此处不能去重

        f3 = f2.groupby(['pdb_id', 'chain_id'], as_index=False)['SEQRES'].agg(lambda x: ''.join(x))
        f3['seqres_len'] = f3['SEQRES'].str.len()
        # 将coordinates信息放入一行
        f4 = f1[['pdb_id', 'chain_id', 'Coordinates', 'inside_code']]  # 此处不能去重
        f5 = f4.groupby(['pdb_id', 'chain_id'], as_index=False)['Coordinates'].agg(lambda x: ''.join(x))
        f5['coordinates_len'] = f5['Coordinates'].str.replace('?', '').str.len()
        # 合并两部分信息
        f6 = pd.merge(f3, f5, on=['pdb_id', 'chain_id'], how='left')

        # 提取所有蛋白链的chain_id
        allchain = f6[['pdb_id', 'chain_id']].drop_duplicates()
        allchain1 = allchain.sort_values(by=['chain_id']).groupby(['pdb_id'], as_index=False).agg(lambda x: ','.join(x))
        allchain1.rename(columns={'chain_id': 'pdb_protein_chain_id'}, inplace=True)
        allchain1.to_csv(outpath + cgc_protein_chain_id, sep='\t', index=False)
        ##以下针对非蛋白链
        ff1 = f[f['chain_type'] != 'protein']# .reset_index(drop=True)
        if len(ff1) != 0:
            ff2 = ff1[['pdb_id', 'chain_id', 'SEQRES', 'inside_code']]  # 此处不能去重
            ff3 = ff2.groupby(['pdb_id', 'chain_id'], as_index=False)['SEQRES'].agg(lambda x: ''.join(x))
            ff3['seqres_len'] = ff3['SEQRES'].str.replace('D', '').str.len()
            # 将coordinates信息放入一行
            ff4 = ff1[['pdb_id', 'chain_id', 'Coordinates', 'inside_code']]  # 此处不能去重
            ff5 = ff4.groupby(['pdb_id', 'chain_id'], as_index=False)['Coordinates'].agg(lambda x: ''.join(x))
            ff5['coordinates_len'] = ff5['Coordinates'].str.replace('D', '').str.replace('?', '').str.len()
            # 合并两部分信息
            ff6 = pd.merge(ff3, ff5, on=['pdb_id', 'chain_id'], how='left')
            ff6['Coordinates'] = ff6['Coordinates'].str.replace('D', '')
            full = pd.concat([f6, ff6], axis=0)
        else:
            full = f6
        #将修饰信息加入文件中
        full['Modification_position'] = full.apply(lambda x: get_modification(x.Coordinates, '1'), axis=1)
        full['Modification_num'] = full.apply(lambda x: get_modification(x.Coordinates, '2'), axis=1)
        full['Modification_position_seqres_index'] = full.apply(lambda x: get_modification_seqres_index(x.Coordinates),axis=1)
        full.to_csv(outpath + integration_file, sep='\t', index=False)

    def add_modification_pdb_type_to_integraton_file(self):
        def modification_type(Modification_position,coordinates_len):
            if pd.isnull(Modification_position):
                modification_site='no_modification'
            else:
                modification_site=[]
                modification_list = str(Modification_position).split(',')
                for i in modification_list:
                    #判断modify的类型，i的类型定义为整型
                    if int(i) <= 5:
                        modify='start'
                    elif int(i) >= int(coordinates_len)-5:
                        modify='end'
                    else:
                        modify='middle'
                    #如果出现过就不放入list中
                    if modify not in modification_site:
                        modification_site.append(modify)
                modification_site = ','.join(modification_site)
            return modification_site

        outpath = self.CONFIG['OUTPUT_FOLDER']
        integration_file = MMCIF_unit.CONFIG['CHAIN_AND_SEQRES_FILE_LIST'][-1]
        chain_type_file_all = MMCIF_unit.CONFIG['CHAIN_TYPE_FILE_LIST'][-1]
        integration_file_new = MMCIF_unit.CONFIG['ADD_MODIFICATION_FILE']

        # 将modification_type、chain_type、pdb_type信息加入到integration文件中
        ff = pd.read_csv(outpath + integration_file, sep='\t')
        ff['Modification_position'] = ff.apply(lambda x: x['Modification_position'].replace('[', '').replace(']', '').replace(' ', '') if isinstance(x['Modification_position'], str) else np.nan, axis=1)
        ff['modification_site'] = ff.apply(lambda x: modification_type(x.Modification_position, x.coordinates_len), axis=1)
        '''
        pdb_type = pd.read_csv(outpath + pdb_and_sifts_protein_chain, sep='\t')
        pdb_type1 = pdb_type[['pdb_id', 'pdb_type', 'pdb_protein_chain_id']].drop_duplicates()
        ff1 = pd.merge(ff, pdb_type1, on=['pdb_id'], how='left')
        '''
        chain_type = pd.read_csv(outpath + chain_type_file_all, sep='\t', dtype=str)
        chain_type.rename(columns={'chain_id_new': 'chain_id'}, inplace=True)
        ff2 = pd.merge(ff, chain_type, on=['pdb_id', 'chain_id'], how='left') # ff1(before)
        ff2.to_csv(outpath + integration_file_new, sep='\t', index=False)

    def add_missing_coordinates_start_end(self):
        def getmisindex(a):
            str = a
            word = '\\?'
            b = [m.start() + 1 for m in re.finditer(word, str)]
            if b != []:
                return b
            else:
                return ''

        def select_UNK(m):# m is the Seqres'content for one line
            if len(set(m))==1 and '!' in list(set(m)):
                return 'yes'
            else:
                return 'no'

        def mis_or_not(a):
            if '?' in a:
                return 'yes'
            else:
                return 'no'

        outpath = self.CONFIG['OUTPUT_FOLDER']
        integration_file_new = MMCIF_unit.CONFIG['ADD_MODIFICATION_FILE']
        all_chain_and_length, integration_new_missing_range = MMCIF_unit.CONFIG['ADD_MISSING_FILE']
        coordinates_file = MMCIF_unit.CONFIG['ATOM_SEQRES_FILE_LIST'][-1]

        ff = pd.read_csv(outpath + integration_file_new, sep='\t', dtype=str)
        ff['chain_and_length'] = ff['chain_id'] + ':' + ff['coordinates_len']
        # protein链和长度合并
        ff1 = ff[ff['chain_type'] == 'protein'][['pdb_id', 'chain_and_length']]
        ff1.rename(columns={'chain_and_length': 'protein_chain_and_length'}, inplace=True)
        protein = ff1.sort_values(by=['protein_chain_and_length']).groupby(['pdb_id'], as_index=False).agg(
            lambda x: ','.join(x))
        # DNA链和长度合并
        ff2 = ff[ff['chain_type'] == 'DNA'][['pdb_id', 'chain_and_length']]
        ff2.rename(columns={'chain_and_length': 'DNA_chain_and_length'}, inplace=True)
        DNA = ff2.sort_values(by=['DNA_chain_and_length']).groupby(['pdb_id'], as_index=False).agg(lambda x: ','.join(x))
        # RNA链和长度合并
        ff3 = ff[ff['chain_type'] == 'RNA'][['pdb_id', 'chain_and_length']]
        ff3.rename(columns={'chain_and_length': 'RNA_chain_and_length'}, inplace=True)
        RNA = ff3.sort_values(by=['RNA_chain_and_length']).groupby(['pdb_id'], as_index=False).agg(lambda x: ','.join(x))
        # mix链和长度合并
        ff4 = ff[ff['chain_type'] == 'RNA+DNA'][['pdb_id', 'chain_and_length']]
        ff4.rename(columns={'chain_and_length': 'mix_chain_and_length'}, inplace=True)
        mix = ff4.sort_values(by=['mix_chain_and_length']).groupby(['pdb_id'], as_index=False).agg(lambda x: ','.join(x))
        # 合并上面四类
        all1 = pd.merge(protein, DNA, on=['pdb_id'], how='left')
        all2 = pd.merge(all1, RNA, on=['pdb_id'], how='left')
        all3 = pd.merge(all2, mix, on=['pdb_id'], how='left')
        '''
        pdb_type = pd.read_csv(outpath + pdb_and_sifts_protein_chain, sep='\t', dtype=str)
        pdb_type1 = pdb_type[['pdb_id', 'pdb_type']].drop_duplicates()
        all4 = pd.merge(all3, pdb_type1, on=['pdb_id'], how='left')
        '''
        all3.to_csv(outpath + all_chain_and_length, sep='\t', index=False) # all4(before)

        # 添加整条链均为空的注释
        ff['UNK_ALL_IN_CHAIN'] = ff.apply(lambda x: select_UNK(x.SEQRES), axis=1)
        unk_file = ff[ff['UNK_ALL_IN_CHAIN'] == 'yes']
        unk_pdb = set(unk_file['pdb_id'])
        integration_file1 = ff[ff['pdb_id'].isin(unk_pdb)]
        integration_file1['only_contains_unk_in_chain_pdb'] = 'yes'
        integration_file2 = ff[~ff['pdb_id'].isin(unk_pdb)]
        integration_file2['only_contains_unk_in_chain_pdb'] = 'no'
        integration_file3 = pd.concat([integration_file1, integration_file2], axis=0)
        # 将统计的信息与整合后的信息合并
        statistics = pd.merge(integration_file3, all3, on=['pdb_id'], how='left') # , 'pdb_type'], how='left') # all4(before)
        # 增加丢失信息
        statistics['mis_index'] = statistics.apply(lambda x: getmisindex(x.Coordinates), axis=1)
        statistics['mis_range'] = statistics.apply(lambda x: MMCIF_unit.getInterval(x.mis_index), axis=1)
        '''
        # statistics['mis_index1'] = statistics['mis_index'].astype(str).str.replace('[','').str.replace('','')
        statistics['ismis'] = statistics.apply(lambda x: mis_or_not(x.Coordinates), axis=1)
        statistics['mis_each_len'] = statistics.apply(lambda x: geteachmis_len(x.mis_range), axis=1)
        statistics['mis_distance_each'] = statistics.apply(lambda x: getdistance_eachmis(x.mis_range), axis=1)
        statistics['mis_distance_judgeup5'] = statistics.apply(lambda x: judge_distance5(x.mis_distance_each), axis=1)
        statistics['mis_count'] = statistics.apply(lambda x: miscount(x.Coordinates), axis=1)
        '''
        coordinates_range = pd.read_csv(outpath + coordinates_file, sep='\t', dtype=str)
        statistics1 = pd.merge(statistics, coordinates_range, on=['pdb_id', 'chain_id'], how='left')
        statistics1.to_csv(outpath + integration_new_missing_range, sep='\t', index=False)

    def get_pdb_muta_info(self):
        @retry(stop_max_attempt_number=3, wait_fixed=1000)
        def getmutation_poly(pdbid):
            fw = open(self.CONFIG['OUTPUT_FOLDER'] + MMCIF_unit.CONFIG['PDB_MUTATION_FILE'], 'a')
            url = 'https://www.rcsb.org/structure/{}'.format(pdbid)
            r = requests.get(url, headers=MMCIF_unit.CONFIG['HEADERS'])
            soup = BeautifulSoup(r.text, 'html.parser')
            s1 = soup.find(id='MacromoleculeTable')
            table = s1.find_all(class_='table-responsive')
            table_all = []
            for i in table:
                table2 = i.find(class_='table table-bordered table-condensed')
                table3 = table2.find('tbody')
                table4 = table3.find(id=re.compile(r'macromolecule-entityId-'))
                table5_chain = table4.find(class_='ellipsisToolTip').text
                table5_mutation_0 = table4.find_all('td')[4]
                table5_mutation_1 = re.split(r'\xa0', table5_mutation_0.text)[0].split(': ')[1]
                table_all.append([pdbid, str(table5_chain), str(table5_mutation_1)])

            fw.write('\n'.join([str('\t'.join(x)) for x in table_all]))
            fw.write('\n')
            fw.close()
            time.sleep(2)

        # div id="MacromoleculeTable"
        # div class="table-responsive"
        # table class="table table-bordered table-condensed"
        # tbody
        # tr id="macromolecule-entityId-3-rowDescription"
        # 第五个td的第一个strong后面

        fw = open(self.CONFIG['OUTPUT_FOLDER'] + MMCIF_unit.CONFIG['PDB_MUTATION_FILE'], 'w')
        fw.write('pdb_id\tchain_id\tmutation_num\n')
        fw.close()
        # 测试
        # pdbid = '1a02'
        # getmutation_poly(pdbid)
        # 读取我们的突变数据集，获取当前的pdb信息
        pool = Pool(processes=10)
        pool.map(getmutation_poly, self.pdb_list)

    def get_resolution(self):
        @retry(stop_max_attempt_number=3, wait_fixed=1000)
        def getResolution(pdbid):
        #     output infomation of error
            fw=open(self.CONFIG['OUTPUT_FOLDER'] + MMCIF_unit.CONFIG['RESOLUTION_FILE'][0],'a')
        #     output resolution of pdb
            fw2=open(self.CONFIG['OUTPUT_FOLDER'] + MMCIF_unit.CONFIG['RESOLUTION_FILE'][1],'a')
            url = 'https://www.rcsb.org/structure/{}'.format(pdbid)
            r = requests.get(url, headers=MMCIF_unit.CONFIG['HEADERS'])
            soup = BeautifulSoup(r.text, 'html.parser')
            try:
                s1 = soup.find(id='exp_header_0_diffraction_resolution')  # x-ray
                s2 = soup.find(id='exp_header_0_em_resolution')   # ELECTRON MICROSCOPY
                s3 = soup.find(id='exp_header_0_method')          # NMR
                if s1:
                    resolution = re.split(r'\xa0',s1.text)[1]
                    method = 'x-ray'
                elif s2:
                    resolution = re.split(r'\xa0',s2.text)[1]
                    method = 'electron'
                elif s3:
                    resolution = 'none'
                    method = 'nmr'
                fw2.write(pdbid+'\t'+str(resolution)+'_'+method+'\n')
                return [pdbid,str(resolution)+'_'+method]#②
        #         return pdbid+'\t'+str(resolution)+'_'+method+'\n'①
            except Exception as e:  # 不属于以上三种情况
                fw.write(e)
            fw.close()
            fw2.close()
            time.sleep(2)
        # use my pdblist to running function
        fw=open(self.CONFIG['OUTPUT_FOLDER'] + MMCIF_unit.CONFIG['RESOLUTION_FILE'][1],'w')
        fw.write('pdb_id\tresolution\n')
        fw.close()
        # give 10 processes
        pool=Pool(processes=10)
        pool.map(getResolution, self.pdb_list)

    def get_year_info1(self):
        @retry(stop_max_attempt_number=3, wait_fixed=1000)
        def get_date1(pdbid):
            print(pdbid)
            fw = open(MMCIF_unit.CONFIG['YEAR_INFO_1_FILE'][0], 'a')
            fp = open(MMCIF_unit.CONFIG['YEAR_INFO_1_FILE'][1], 'a')
            url = 'https://www.rcsb.org/structure/' + pdbid
            r = requests.get(url, headers=MMCIF_unit.CONFIG['HEADERS'])
            soup = BeautifulSoup(r.text, 'html.parser')
            s1 = soup.findAll("div", {"class": re.compile("col-md-6 col-sm-6 col-xs-12 col-xs-12")})
            try:
                a = str(s1[-1].text.split(':')[1].split('Type')[0].strip(' '))
                b = str(s1[-1].text.split(':')[-2].split('Type')[0].strip(' '))
                fw.write('%s\t%s\t%s\n' % (pdbid, a, b))
            except:
                fp.write(pdbid + '\n')
            fw.close()
            fp.close()
            time.sleep(2)

        fw = open(MMCIF_unit.CONFIG['YEAR_INFO_1_FILE'][0], 'w')
        fw.write('pdb_id\tinitial_version_time\tnewest_version_time\n')
        fw.close()
        pool = Pool(processes=20)
        pool.map(get_date1, self.pdb_list)

    def get_year_info2(self):
        @retry(stop_max_attempt_number=3, wait_fixed=1000)
        def get_date2(pdbid):
            print(pdbid)
            fw = open(MMCIF_unit.CONFIG['YEAR_INFO_2_FILE'][0], 'a')
            fp = open(MMCIF_unit.CONFIG['YEAR_INFO_2_FILE'][1], 'a')
            url = 'https://www.rcsb.org/structure/' + pdbid
            r = requests.get(url, headers=MMCIF_unit.CONFIG['HEADERS'])
            soup = BeautifulSoup(r.text, 'html.parser')
            s1 = soup.findAll("div", {"class": re.compile("col-md-6 col-sm-6 col-xs-12 col-xs-12")})
            try:
                a = str(s1[-1].text.split(':')[1].split('Type')[0].strip(' '))
                b = str(s1[-1].text.split(':')[-2].split('Type')[0].strip(' '))
                fw.write('%s\t%s\t%s\n' % (pdbid, a, b))
            except:
                fp.write(pdbid + '\n')
            fw.close()
            fp.close()
            time.sleep(2)

        fw = open(MMCIF_unit.CONFIG['YEAR_INFO_2_FILE'][0], 'w')
        fw.write('pdb_id\tinitial_version_time\tnewest_version_time\n')
        fw.close()
        pool = Pool(processes=5)
        data_temp = pd.read_csv(MMCIF_unit.CONFIG['YEAR_INFO_1_FILE'][0], sep='\t', names=['pdb_id'])
        pdbid = list(set(data_temp['pdb_id']))
        pool.map(get_date2, pdbid)

    def get_year_info_all(self):
        year1 = pd.read_csv(MMCIF_unit.CONFIG['YEAR_INFO_1_FILE'][0], sep='\t', dtype=str)
        year2 = pd.read_csv(MMCIF_unit.CONFIG['YEAR_INFO_2_FILE'][0], sep='\t', dtype=str)
        yearall = pd.concat([year1, year2], axis=0)
        yearall.to_csv(MMCIF_unit.CONFIG['YEAR_INFO_ALL'], sep='\t', index=False)

    def get_final_file_before_mapmuta(self):
        integration_new_missing_range = MMCIF_unit.CONFIG['ADD_MISSING_FILE'][-1]
        outpath = self.CONFIG['OUTPUT_FOLDER']
        resolution_file = MMCIF_unit.CONFIG['RESOLUTION_FILE'][1]
        pdb_mutation_file = MMCIF_unit.CONFIG['PDB_MUTATION_FILE']
        year_all_info = MMCIF_unit.CONFIG['YEAR_INFO_ALL']
        integration_final, final_file_before_mapping = MMCIF_unit.CONFIG['FINAL_FILE']
        ligand_file_final2 = MMCIF_unit.CONFIG['LIGAND_FILE_LIST'][-1]

        statistics = pd.read_csv(outpath + integration_new_missing_range,sep='\t',dtype=str)
        '''
        statistics['mistype'] = statistics.apply(lambda x:missing_type(x.mis_range,x.seqres_len),axis=1)
        statistics.to_csv(outpath + integration_new_missing_range1, sep='\t',index=False)
        '''
        #文件中加入resolution和pdb_mutation信息
        resolution = pd.read_csv(outpath + resolution_file,sep='\t',dtype=str)
        resolution['resolution_score'] = resolution['resolution'].str.split('_').str[0].str.replace('none','1000')
        resolution['resolution_method'] = resolution['resolution'].str.split('_').str[1]
        resolution1 = resolution.drop(['resolution'],axis=1)
        statistics1 = pd.merge(statistics,resolution1,on=['pdb_id'],how='left')
        pdb_muta = pd.read_csv(outpath + pdb_mutation_file,sep='\t',dtype=str)
        #文件内有空白行
        pdb_muta = pdb_muta[pdb_muta['pdb_id'].notna()]

        # 重新设置索引号，避免同一索引对应不同行
        pdb_muta.index = range(len(pdb_muta))

        # 由于重新设置了索引不会造成混淆，所以可以使用以下方法，比较快
        pdb_muta1 = pdb_muta.drop('chain_id', axis=1).join(
            pdb_muta['chain_id'].str.split(', ', expand=True).stack().reset_index(level=1, drop=True).rename('chain_id'))

        pdb_muta1['WT'] = ['1'if int(x)==0 else '0' for x in pdb_muta1['mutation_num']]
        statistics2 = pd.merge(statistics1,pdb_muta1,on=['pdb_id','chain_id'],how='left')
        # 增加金属配体的信息
        ligand_file = pd.read_csv(outpath + ligand_file_final2, sep='\t', dtype=str)
        statistics3 = pd.merge(statistics2, ligand_file, on=['pdb_id', 'chain_id'], how='left')
        statistics3.rename(columns={'pdb_index': 'ligand_pdb_index', 'position_in_seqres': 'ligand_position_in_seqres'},
                           inplace=True)

        # 增加年份信息
        yearall = pd.read_csv(year_all_info, sep='\t', dtype=str)
        mergeall = pd.merge(statistics3, yearall, on=['pdb_id'], how='left')

        '''
        statistics3.to_csv(outpath + integration_final, sep='\t', index=False)
        sifts_file_zzf = pd.read_csv(outpath + sifts_cgc, sep='\t', dtype=str)
        mergeall1 = pd.merge(sifts_file_zzf, mergeall, on=['pdb_id', 'chain_id'], how='left')
        mergeall2 = mergeall1[mergeall1['SEQRES'].notna()]
        mergeall2.to_csv(outpath + final_file_before_mapping, sep='\t',
                         index=False)
        '''
        mergeall['pdb_id'] = mergeall.apply(lambda x: x['pdb_id'].upper(), axis=1)
        mergeall.to_csv(outpath + final_file_before_mapping, sep='\t',
                         index=False)


class Interactome3D(Unit):
    CONFIG = {
        'DOWNLOAD_FOLDER': '/home/zzf/Work/StructDDG_0427/data/Mapping_Pipeline/Interactome3D_files/',
        'INTERACTION_SIFTS_COL': ['interac_TYPE', 'pdb_id', 'interac_BIO_UNIT',
                                  'interac_FILENAME', 'interac_group_compo',
                                  'UniProt', 'chain_id', 'interac_MODEL',
                                  'interac_SEQ_IDENT', 'interac_COVERAGE',
                                  'interac_DOMAIN', 'interac_model_len',
                                  'interac_model_range'],

    }

    def get_interactions_meta(self, species='human', struct_type=False, filePath=False, related_unp=False, related_pdb=False, outputPath=False):
        if not filePath:
            url = "https://interactome3d.irbbarcelona.org/user_data/%s/download/complete/interactions.dat" % species
            filePath = self.CONFIG['DOWNLOAD_FOLDER'] + 'interactions_%s_%s.dat' % (species, time.strftime("%Y_%m_%d", time.localtime()))
            wget.download(url, out=filePath)
            self.file_list.append(filePath)

        dfrm = pd.read_csv(filePath, sep='\t')
        if struct_type:
            dfrm = dfrm[dfrm['TYPE'] == struct_type].reset_index(drop=True)
        dfrm['PDB_ID'] = dfrm.apply(lambda x: x['PDB_ID'].upper(), axis=1)
        dfrm['group_compo'] = dfrm.apply(lambda x: '%s_%s' % tuple(sorted([x['PROT1'], x['PROT2']])), axis=1)
        common_cols = ['TYPE','PDB_ID','BIO_UNIT','FILENAME', 'group_compo']
        s_cols = ['PROT', 'CHAIN', 'MODEL', 'SEQ_IDENT', 'COVERAGE', 'SEQ_BEGIN', 'SEQ_END', 'DOMAIN']
        get_s_cols = lambda num: ['%s%s' % (i, num) for i in s_cols]

        df1, df2 = dfrm[common_cols+get_s_cols(1)].copy(), dfrm[common_cols+get_s_cols(2)].copy()
        df1.columns, df2.columns = common_cols+s_cols, common_cols+s_cols
        df12 = pd.concat([df1,df2]).reset_index(drop=True)
        df12['model_len'] = df12.apply(lambda x: x['SEQ_END'] - x['SEQ_BEGIN'] + 1, axis=1)
        df12['model_range'] = df12.apply(lambda x: '[[%d, %d]]'%(x['SEQ_BEGIN'],x['SEQ_END']), axis=1)

        if related_unp:
            df12 = df12[df12['PROT'].isin(related_unp)]

        if related_pdb:
            df12 = df12[df12['PDB_ID'].isin(related_pdb)]

        self.file_o(outputPath, df12)
        return df12

    def add_SIFTS_to_interactions_meta(self, sifts_df=False, sifts_filePath=False, interac_df=False, interac_filePath=False, outputPath=False):
        sifts_dfrm = self.file_i(sifts_filePath, sifts_df, ('sifts_filePath', 'sifts_df'))
        interac_dfrm = self.file_i(interac_filePath, interac_df, ('interac_filePath', 'interac_df'))
        interac_dfrm = interac_dfrm[interac_dfrm['TYPE'] == 'Structure'].reset_index(drop=True)
        interac_dfrm.drop(columns=['SEQ_BEGIN', 'SEQ_END'], inplace=True)
        interac_dfrm.columns = self.CONFIG['INTERACTION_SIFTS_COL']
        dfrm = pd.merge(interac_dfrm, sifts_dfrm, on=['pdb_id', 'chain_id', 'UniProt'], how='left')
        self.file_o(outputPath, dfrm)
        return dfrm



if __name__ == '__main__':
    # interact_df = pd.read_csv('../../data/Mapping_Pipeline/Interactome3D_files/interactions_human_0814.tsv', sep='\t')
    # pdb_set = set(interact_df['PDB_ID'])
    # pdb_list = ['2MSE', '4OV6', '5HHM', '5HHO', '5NQK', '1AO7', '4ZDH', '3PL6']
    demo_sifts = SIFTS_unit()
    demo_sifts.deal_with_insertionDeletion_SIFTS(sifts_filePath='../../data/pdb_uniprot_SIFTS_NEW0813.tsv', outputPath='../../data/Mapping_Pipeline/sifts_files/pdb_uniprot_SIFTS_NEW0813.tsv')
