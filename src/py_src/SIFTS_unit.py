# @Date:   2019-08-16T23:24:17+08:00
# @Email:  1730416009@stu.suda.edu.cn
# @Filename: SIFTS_unit.py
# @Last modified time: 2019-08-28T22:13:19+08:00
import pandas as pd
import numpy as np
import json, wget, gzip, time, sys
from urllib import request
from itertools import combinations
sys.path.append('./')
from Unit import Unit


class SIFTS_unit(Unit):
    CONFIG = {
        'PDB_ID': 'pdb_id',
        'RAW_SIFTS_COLUMNS':[
            'pdb_id', 'chain_id', 'UniProt', 'identity', 'identifier',
            'pdb_start', 'pdb_end', 'unp_start', 'unp_end',
            'is_canonical', 'start', 'end', 'entity_id', 'struct_asym_id'],
        'DOWNLOAD_FOLDER': '../../data/sifts_files/',
        'UNP_LIST_PATH': '../../data/sifts_files/sifts_uniprot_list.tsv',
        'ELE_LIST': ["coordinates_len", "mappedOut", "metal_ligand_num", "delHT_MissingNum", "if_2", "if_1"],
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

        try:
            rows = len(pd.read_csv(outputPath, sep='\t', usecols=['UniProt']))
        except FileNotFoundError:
            rows = 0
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
                        df.to_csv(outputPath, index=False, sep='\t', mode='a+', header=False)# Attention: 初始Build dataset时需指定header
            except Exception as e:
                print('PDB: [', pdbId, ']\n', e)
                fail_list.append(pdbId)
            current += 1
            print('getSiftsInfo(): End a circle.[',pdbId,'] current:',current,'ALL:',allPDB)
        return rows, fail_list

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

    def handle_SIFTS(self, sifts_filePath=False, skiprows=0, outputPath=False):
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
        sifts_df = pd.read_csv(filePath, sep='\t', na_values=[SIFTS_unit.CONFIG['PDB_ID'], '', None], keep_default_na=False, skiprows=skiprows, names=SIFTS_unit.CONFIG['RAW_SIFTS_COLUMNS'])
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
        dfrm = pd.merge(sifts_dfrm, mmcif_dfrm, on=['pdb_id', 'chain_id', 'entity_id'], how='left')
        self.file_o(outputPath, dfrm)
        return dfrm

    def add_unp_len_SIFTS(self, sifts_df=False, sifts_filePath=False, unpLen_df=False, unpLen_filePath=False, outputPath=False, sep='\t'):
        # unpLen_filePath = '/home/zzf/Work/StructDDG_0427/data/Mapping_Pipeline/sifts_unp_len_list.csv'
        sifts_dfrm = self.file_i(sifts_filePath, sifts_df, ('sifts_filePath', 'sifts_df'))
        sifts_dfrm['Entry'] = sifts_dfrm.apply(lambda x: x['UniProt'].split('-')[0], axis=1)
        unpLen_dfrm = self.file_i(unpLen_filePath, unpLen_df, ('unpLen_filePath', 'unpLen_df'), sep=sep)
        # unpLen_dfrm.drop(columns=['yourlist'], inplace=True)
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
            return SIFTS_unit.getInterval(pdb_range_set - mis_range_set)

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
                mis_li = json.loads(li.replace('\'', '"'))
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

        def getUsefulChainNum(sli, cutoff):
            li = json.loads(sli.replace('(','[').replace(')',']'))
            return len(list(filter(lambda x: int(x[0]) > cutoff, li)))
            '''
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
            '''

        sifts_dfrm = self.file_i(sifts_filePath, sifts_df, ('sifts_filePath', 'sifts_df'))
        '''
        # Metal ligand Info
        sifts_dfrm['metal_list'] = sifts_dfrm.apply(lambda x: list(map(int, x['ligand_position_in_seqres'].split(
            ';'))) if not isinstance(x['ligand_position_in_seqres'], float) else np.nan, axis=1)
        '''
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

        '''sifts_dfrm['metal_count'] = sifts_dfrm.apply(lambda x: len(x['metal_list']) if not isinstance(x['metal_list'], float) else 0, axis=1)'''

        select_vector = get_weight()
        calScore(sifts_dfrm, SIFTS_unit.CONFIG['SCORE_COL'], SIFTS_unit.CONFIG['ELE_LIST'], select_vector)

        # Add resolution-related socre
        deal_rs = lambda x: -float(x.split(',')[0]) if ',' in x else -1200
        def deal_reso_score(score):
            if isinstance(score, str):
                if len(set(score) & set('?,')) == 0:
                    return -float(score)
                else:
                    return deal_rs(score)
            else:
                return -score
        sifts_dfrm['ne_resolution_score'] = sifts_dfrm.apply(lambda x: deal_reso_score(x['resolution_score']), axis=1)

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
                'contains_unk_in_chain_pdb': (False, 'eq'),
                'UNK_ALL_IN_CHAIN': (False, 'eq'),
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

    def map_muta_from_PDB_to_UNP(self, sifts_df=False, sifts_filePath=False, outputPath=False):
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

    def map_muta_from_UNP_to_PDB(self, mutaSiteCol, sifts_df=False, sifts_filePath=False, outputPath=False):
        sifts_dfrm = self.file_i(sifts_filePath, sifts_df, ('sifts_filePath', 'sifts_df'))



        file_o(outputPath, sifts_dfrm)
        return sifts_dfrm



if __name__ == '__main__':
    demo = SIFTS_unit()
    demo.set_lists(['Test is success.'], [])
    print(demo.pdb_list)
