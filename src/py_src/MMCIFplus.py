# @Date:   2019-09-09T16:32:43+08:00
# @Email:  1730416009@stu.suda.edu.cn
# @Filename: MMCIFplus.py
# @Last modified time: 2019-09-10T19:45:56+08:00
import os
import time
import requests
import sys
import json
import pandas as pd
import numpy as np
from urllib import request
from retrying import retry
from multiprocessing.dummy import Pool
from collections import defaultdict, Iterable, Iterator
from Bio.File import as_handle
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from Unit import Unit
sys.path.append('./')

MMCIF_FILE_FOLDER = {
    'MMCIF_OLD_FOLDER': [
        '/data1/suntt/process0606/cgc_mmcif_file/',
        '/data1/suntt/CanDriver/Data/PDB_cgc/cgc_mmcif_file/',
        '/data1/suntt/CanDriver/Data/PDB_NEW/mmcif_file/'
        ],
    'MMCIF_NEW_FOLDER': '/home/zzf/Work/SIFTS_Plus_Muta_Maps/data/mmcif_file/',
}

CONFIG = {
    'LIGAND_LIST': [
        'ZN', 'MG', 'CA', 'FE', 'NA', 'MN', 'K', 'NI', 'CU', 'CO', 'CD', 'HG', 'PT', 'MO', 'BE', 'AL', 'BA',
        'RU', 'SR', 'V', 'CS', 'W', 'AU', 'YB', 'LI', 'GD', 'PB', 'Y', 'TL', 'IR', 'RB', 'SM', 'AG',
        'OS', 'PR', 'PD', 'EU', 'RH', 'RE', 'TB', 'TA', 'LU', 'HO', 'CR', 'GA', 'LA', 'SN', 'SB', 'CE',
        'ZR', 'ER', 'TH', 'TI', 'IN', 'HF', 'SC', 'DY', 'BI', 'PA', 'PU', 'AM', 'CM', 'CF', 'GE', 'NB', 'TC',
        'ND', 'PM', 'TM', 'PO', 'FR', 'RA', 'AC', 'NP', 'BK', 'ES', 'FM', 'MD', 'NO', 'LR', 'RF', 'DB', 'SG'],
    'HEADERS': {'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/64.0.3282.140 Safari/537.36 Edge/17.17134'},
    'CHAIN_TYPE_DICT': {
        'polypeptide(L)': 'protein', 'polypeptide(D)': 'protein', 'polydeoxyribonucleotide': 'DNA',
        'polyribonucleotide': 'RNA', 'polydeoxyribonucleotide/polyribonucleotide hybrid': 'DNA+RNA'},
}

DEFAULT_COLS = {
    'COMMON_COL': ['data_', '_pdbx_audit_revision_history.revision_date', '_exptl.method', '_em_3d_reconstruction.resolution', '_refine.ls_d_res_high'],
    'ENTITY_COL': ['_entity.pdbx_mutation', '_entity.id'],
    'TYPE_COL': ['_entity_poly.entity_id', '_entity_poly.pdbx_strand_id', '_entity_poly.type'],
    'SEQRES_COL': ['_pdbx_poly_seq_scheme.pdb_strand_id',
                   '_pdbx_poly_seq_scheme.mon_id', '_pdbx_poly_seq_scheme.pdb_mon_id', '_pdbx_poly_seq_scheme.auth_mon_id',
                   '_pdbx_poly_seq_scheme.ndb_seq_num', '_pdbx_poly_seq_scheme.pdb_seq_num',
                   '_pdbx_poly_seq_scheme.auth_seq_num', '_pdbx_poly_seq_scheme.pdb_ins_code'],
    'LIGAND_COL': [
        '_struct_conn.conn_type_id', '_struct_conn.ptnr1_auth_comp_id', '_struct_conn.ptnr2_auth_comp_id',
        '_struct_conn.ptnr1_auth_asym_id', '_struct_conn.ptnr2_auth_asym_id',
        '_struct_conn.ptnr1_auth_seq_id', '_struct_conn.ptnr2_auth_seq_id'],
    'BIOASS_COL': ['_pdbx_struct_assembly_gen.assembly_id',
                   '_pdbx_struct_assembly_gen.oper_expression',
                   '_pdbx_struct_assembly_gen.asym_id_list',
                   '_pdbx_struct_assembly.oligomeric_count'],
    'METAL_LIGAND_COL': ['metal_ligand_chain_id', 'metal_ligand_content'],
}

class MMCIF2DictPlus(MMCIF2Dict):
    """Parse a MMCIF file and return a dictionary"""

    def __init__(self, filename, useKeyList):
        """
        <h1>Parse a mmCIF file and return a dictionary</h1>
        <h3>Override ```__init__``` method of ```Bio.PDB.MMCIF2Dict.MMCIF2Dict```</h3>
        <ul>
            <li>```filename```: "Name of the PDB file OR an open filehandle."</li>
            <li>```useKeyList```: "Those MMCIF-keys that user want to collect data from."</li>
        </ul>
        """
        self.quote_chars = ["'", '"']
        self.whitespace_chars = [" ", "\t"]
        with as_handle(filename) as handle:
            loop_flag = False
            key = None
            tokens = self._tokenize(handle)
            try:
                token = next(tokens)
            except StopIteration:
                return  # for Python 3.7 and PEP 479
            self[token[0:5]] = token[5:]
            i = 0
            n = 0
            use = []
            for token in tokens:
                if token.lower() == "loop_":
                    loop_flag = True
                    keys = []
                    i = 0
                    n = 0
                    use = []
                    continue
                elif loop_flag:
                    # The second condition checks we are in the first column
                    # Some mmCIF files (e.g. 4q9r) have values in later columns
                    # starting with an underscore and we don't want to read
                    # these as keys
                    if token.startswith("_") and (n == 0 or i % n == 0):
                        if i > 0:
                            loop_flag = False
                        else:
                            # Additional statement:focus on target data
                            if token in useKeyList:
                                use.append(n)
                                self[token] = []
                            keys.append(token)
                            n += 1
                            continue
                    else:
                        key_index = i % n
                        try:
                            if key_index in use:
                                self[keys[key_index]].append(token)
                        except Exception as e:
                            print(keys, key_index, use)
                            raise Exception(e)
                        i += 1
                        continue
                if key is None:
                    # Additional statement:focus on target data
                    if token in useKeyList:
                        key = token
                else:
                    # Always returns a list
                    self[key] = [token]
                    key = None


class MMCIF2Dfrm(Unit):
    """Convert MMCIF data into a DataFrame"""

    FUNC_LI_DI = []
    FUNC_LI_DF = []

    @property
    def default_use_keys(self):
        use_keys = []
        use_col_li = ['COMMON_COL', 'ENTITY_COL', 'TYPE_COL', 'SEQRES_COL', 'LIGAND_COL', 'BIOASS_COL']
        for col in use_col_li:
            use_keys.extend(DEFAULT_COLS[col])
        return use_keys

    def checkEntityType(sli, cli):
        if isinstance(sli, str):
            sli = json.loads(sli.replace('\'', '"'))
        if isinstance(cli, str):
            cli = json.loads(cli.replace('\'', '"'))
        li = list(zip(sli, cli))
        pli = list(filter(lambda x: 'polypeptide' in x[0], li))
        pli_len = len(pli)
        if pli_len == 1:
            count = pli[0][1].count(',')
            if count == 0:
                return 'mo'
            elif count > 0:
                return 'ho:%d' % (count + 1)
        elif pli_len > 1:
            return 'he:' + ';'.join(i[1] for i in pli)

    def get_mmcif_dict(info_key, info_dict, path):
        '''
        <h1>Get MMCIF info in dict-format with the help of ```MMCIF2DictPlus```</h1>
        <b>Creating lists of lists stored in a ```defaultdict``` and prepare for converting it into a ```DataFrame```.</b>
        <h3>Param</h3>
        <ul>
            <li>```info_key```: "Those MMCIF-keys with multiple data which are converted into a list by ```MMCIF2DictPlus```."</li>
            <li>```info_dict```: "The defaultdict that stores data."</li>
            <li>```path```: "The file path of MMCIF-file."</li>
        </ul>
        '''
        assert isinstance(info_key, (Iterable, Iterator)), "Invalid Input, info_key should be iterable"
        assert isinstance(info_dict, defaultdict), "Invalid Input, info_dict should be a defaultdict"

        # info_dict = defaultdict(list)
        mmcif_dict = MMCIF2DictPlus(path, info_key)
        for key in info_key:
            data = mmcif_dict.get(key, np.nan)
            info_dict[key].append(data)

    def dispatch_on_set(keys, func_li):
        """Decorator to add new dispatch functions."""
        def register(func):
            func_li.append((func, set(keys)))
            return func
        return register

    def handle_mmcif_data(query, data, fun_li):
        use = False
        for func, keySet in fun_li:
            if set(query) >= keySet:
                func(data)
                use = True
        return use

    def get_index(x, y, z): return y[x[z]:x[z + 1]] if len(x) != 1 and z + 1 < len(x) else y[x[z]:]

    @dispatch_on_set(DEFAULT_COLS['SEQRES_COL'], func_li=FUNC_LI_DI)
    def handle_seqres_di(info_dict):
        # Deal with SEQRES_COL
        resides_col_li = DEFAULT_COLS['SEQRES_COL'][1:4]
        mtoTool = Unit.MultiToOne()
        for i in range(len(info_dict[resides_col_li[0]])):
            for resides_col in resides_col_li:
                info_dict[resides_col][i] = ''.join(
                    mtoTool.multi_letter_convert_to_one_letter(j) for j in info_dict[resides_col][i])

        pdbx_poly_key = DEFAULT_COLS['SEQRES_COL'][0]
        for i in range(len(info_dict[pdbx_poly_key])):
            strand_id_index = [0]
            li = info_dict[pdbx_poly_key][i]
            save_id = li[0]
            strand_id_li = [save_id]
            for j in range(len(li)):
                if li[j] != save_id:
                    save_id = li[j]
                    strand_id_index.append(j)
                    strand_id_li.append(save_id)
            info_dict[pdbx_poly_key][i] = strand_id_li

            for col in DEFAULT_COLS['SEQRES_COL'][1:4]:
                info_dict[col][i] = [
                    MMCIF2Dfrm.get_index(strand_id_index, info_dict[col][i], j)
                    for j in range(len(strand_id_index))]

            for col in DEFAULT_COLS['SEQRES_COL'][4:]:
                info_dict[col][i] = [';'.join(
                    MMCIF2Dfrm.get_index(strand_id_index, info_dict[col][i], j))
                    for j in range(len(strand_id_index))]

    @dispatch_on_set(DEFAULT_COLS['LIGAND_COL'], func_li=FUNC_LI_DI)
    def handle_ligand_di(info_dict):
        ligand_col_list = DEFAULT_COLS['LIGAND_COL']
        metal_li = CONFIG['LIGAND_LIST']

        for i in range(len(info_dict[ligand_col_list[0]])):
            temp = pd.isna(info_dict[ligand_col_list[0]][i])
            if temp is True:
                info_dict[DEFAULT_COLS['METAL_LIGAND_COL']
                          [0]].append(np.nan)
                info_dict[DEFAULT_COLS['METAL_LIGAND_COL']
                          [1]].append(np.nan)
                continue
            ligand_col_tp = tuple(info_dict[col][i] for col in ligand_col_list)
            ligand_col_zip_li = list(zip(*ligand_col_tp))

            aa_li = list(MMCIF2Dfrm.SEQ_DICT.keys())[:21]
            metal_ligand_info = list(
                filter(lambda x: x[0] == 'metalc', ligand_col_zip_li))
            # chain_id: _struct_conn.ptnr2_auth_asym_id [4]
            sub_metal_ligand_info_1 = filter(
                lambda x: x[1] in metal_li and x[2] in aa_li, metal_ligand_info)
            # chain_id: _struct_conn.ptnr1_auth_asym_id [3]
            sub_metal_ligand_info_2 = filter(
                lambda x: x[2] in metal_li and x[1] in aa_li, metal_ligand_info)

            new_metal_ligand_info = []
            for tp in sub_metal_ligand_info_1:
                new_metal_ligand_info.append(
                    (tp[4], tp[1], tp[5], tp[2], tp[6]))
            for tp in sub_metal_ligand_info_2:
                new_metal_ligand_info.append(
                    (tp[3], tp[2], tp[6], tp[1], tp[5]))

            new_metal_ligand_info.sort(key=lambda x: x[0])
            try:
                save_id = new_metal_ligand_info[0][0]
                # print(new_metal_ligand_info)
            except IndexError:
                info_dict[DEFAULT_COLS['METAL_LIGAND_COL']
                          [0]].append(np.nan)
                info_dict[DEFAULT_COLS['METAL_LIGAND_COL']
                          [1]].append(np.nan)
                continue

            strand_id_li = [save_id]
            strand_id_index = [0]
            for j in range(len(new_metal_ligand_info)):
                if new_metal_ligand_info[j][0] != save_id:
                    save_id = new_metal_ligand_info[j][0]
                    strand_id_index.append(j)
                    strand_id_li.append(save_id)

            info_dict[DEFAULT_COLS['METAL_LIGAND_COL'][0]].append(strand_id_li)
            info_dict[DEFAULT_COLS['METAL_LIGAND_COL'][1]].append(
                [MMCIF2Dfrm.get_index(strand_id_index, [ele[1:] for ele in new_metal_ligand_info], j) for j in range(len(strand_id_index))]
            )

    @dispatch_on_set(DEFAULT_COLS['COMMON_COL'], func_li=FUNC_LI_DF)
    def handle_common_df(df):
        # Deal with the date of structure
        df['initial_version_time'] = df.apply(
            lambda x: x[DEFAULT_COLS['COMMON_COL'][1]][0], axis=1)
        df['newest_version_time'] = df.apply(
            lambda x: x[DEFAULT_COLS['COMMON_COL'][1]][-1], axis=1)
        # Deal with the resolution
        df['resolution'] = df.apply(lambda x: x[DEFAULT_COLS['COMMON_COL'][4]], axis=1)
        df['resolution'] = df.apply(lambda x: x[DEFAULT_COLS['COMMON_COL'][3]] if isinstance(x['resolution'], float) else x['resolution'], axis=1)

    @dispatch_on_set(DEFAULT_COLS['ENTITY_COL'], func_li=FUNC_LI_DF)
    def handle_entity_df(df):
        # Deal with the mutations
        def muta_count(x): return x.count(',') + 1 if x != '?' else 0
        df['mutation_num'] = df.apply(lambda x: [muta_count(i) for i in x['_entity.pdbx_mutation']], axis=1)

    @dispatch_on_set(DEFAULT_COLS['TYPE_COL'], func_li=FUNC_LI_DF)
    def handle_type_df(df):
        # Deal with chain type
        def get_chainType_fun(ele): return CONFIG['CHAIN_TYPE_DICT'].get(ele, 'other')
        df['pdb_contain_chain_type'] = df.apply(lambda x: ','.join(sorted(set(map(get_chainType_fun, json.loads(x['_entity_poly.type'].replace('\'', '"'))))))
                                                if isinstance(x['_entity_poly.type'], str)
                                                else ','.join(sorted(set(map(get_chainType_fun, x['_entity_poly.type'])))), axis=1)
        # Add Info about pdb_type
        df['pdb_type_MMCIF'] = df.apply(lambda x: MMCIF2Dfrm.checkEntityType(
            x['_entity_poly.type'], x['_entity_poly.pdbx_strand_id']), axis=1)

    @dispatch_on_set(['_pdbx_poly_seq_scheme.mon_id', '_entity_poly.type'], func_li=FUNC_LI_DF)
    def handle_unk_df(df):
        # Deal with UNK_ALL in chain
        def get_unk_fun(ele): return len(ele) == ele.count('!')
        df['UNK_ALL_IN_CHAIN'] = df.apply(lambda x: list(map(get_unk_fun, json.loads(x['_pdbx_poly_seq_scheme.mon_id'].replace('\'', '"'))))
                                          if isinstance(x['_entity_poly.type'], str)
                                          else list(map(get_unk_fun, x['_pdbx_poly_seq_scheme.mon_id'])), axis=1)
        # Deal with UNK_ALL in chains of a pdb
        df['contains_unk_in_chain_pdb'] = df.apply(
            lambda x: len(set(x['UNK_ALL_IN_CHAIN'])) == 2, axis=1)

    def mmcif_dict2dfrm(self, path_list, useKeyList=False, outputPath=False):
        '''
        <h1>Convert MMCIF data into a DataFrame</h1>
        <b>Creating a DataFrame stores all the target data of related MMCIF file in ```path_list```.</b>
        <h3>Param</h3>
        <ul>
            <li>```path_list```: "The list of MMCIF file path OR open filehandle"</li>
            <li>```useKeyList```: "The MMCIF-keys that you want to collect data from."</li>
            <li>```outputPath```: "The savepath of the final DataFrame. Default value:```False```"</li>
        </ul>
        <h3>Return</h3>
        <ul><li>```pandas.DataFrame```</li></ul>
        '''
        info_dict = defaultdict(list)
        for path in path_list:
            if path[-3:] != 'cif':
                print('Not a valid MMCIF file path: %s' % path)
            else:
                if not useKeyList:
                    useKeyList = self.default_use_keys
                MMCIF2Dfrm.get_mmcif_dict(useKeyList, info_dict, path)
        # Modify the dict
        modified_di = MMCIF2Dfrm.handle_mmcif_data(useKeyList, info_dict, MMCIF2Dfrm.FUNC_LI_DI)
        if modified_di:
            print('handle_mmcif_data(): Modified Dict')
        # Transform dict into dfrm
        df = pd.DataFrame(info_dict)
        # Modify the dfrm
        modified_df = MMCIF2Dfrm.handle_mmcif_data(useKeyList, df, MMCIF2Dfrm.FUNC_LI_DF)
        if modified_df:
            print('handle_mmcif_data(): Modified Dfrm')
        # Change the columns
        df.rename(
            columns={DEFAULT_COLS['COMMON_COL'][0]: 'pdb_id', DEFAULT_COLS['COMMON_COL'][2]: 'method'}, inplace=True)

        if os.path.exists(outputPath):
            self.file_o(outputPath, df, mode='a+', header=False)
        else:
            self.file_o(outputPath, df)
        return df


if __name__ == '__main__':
    route = 'C:\\Users\\Nature\\Desktop\\LiGroup\\Filter_new_20190123\\doc_in\\spe\\'
    file_list = os.listdir(route)
    file_p_list = [route + i for i in file_list]
    mmcif_demo = MMCIF2Dfrm()
    df = mmcif_demo.mmcif_dict2dfrm(file_p_list, useKeyList=['data_', '_entity_poly.type'])

    for i in df.index:
        print(df.loc[i, ])
