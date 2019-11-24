---
title: "CODE"
author: ZeFeng Zhu
date: Oct 30, 2019
output:
  word_document:
    path: C:/Users/Nature/Desktop/MySoftWare/code2.docx
export_on_save:
pandoc: true
---

```py
# @Date:   2019-09-09T16:32:43+08:00
# @Email:  1730416009@stu.suda.edu.cn
# @Filename: MMCIFplus.py
# @Last modified time: 2019-10-29T17:37:00+08:00
import os
import sys
import json
import pandas as pd
import numpy as np
from collections import defaultdict, Iterable, Iterator
from Bio.File import as_handle
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from Unit import Unit
from RetrievePDB import MPWrapper
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
                   '_pdbx_poly_seq_scheme.auth_seq_num', '_pdbx_poly_seq_scheme.pdb_ins_code',
                   '_pdbx_poly_seq_scheme.asym_id'],
    'LIGAND_COL': [
        '_struct_conn.conn_type_id', '_struct_conn.ptnr1_auth_comp_id', '_struct_conn.ptnr2_auth_comp_id',
        '_struct_conn.ptnr1_auth_asym_id', '_struct_conn.ptnr2_auth_asym_id',
        '_struct_conn.ptnr1_auth_seq_id', '_struct_conn.ptnr2_auth_seq_id'],
    'BIOASS_COL': ['_pdbx_struct_assembly_gen.assembly_id',
                   '_pdbx_struct_assembly_gen.oper_expression',
                   '_pdbx_struct_assembly_gen.asym_id_list',
                   '_pdbx_struct_assembly.oligomeric_count'],
    'NON_POLY_COL': ['_pdbx_entity_nonpoly.entity_id', '_pdbx_entity_nonpoly.name', '_pdbx_entity_nonpoly.comp_id'],
    'COORDINATE_MODEL_COL': ['_pdbx_coordinate_model.asym_id', '_pdbx_coordinate_model.type'],
    'METAL_LIGAND_COL': ['metal_ligand_chain_id', 'metal_ligand_content'],
}


class MMCIF2DictPlus(MMCIF2Dict):
    """# Parse a MMCIF file and return a dictionary"""

    def __init__(self, filename, useKeyList):
        """
        <h1>Parse a mmCIF file and return a dictionary</h1>
        <h3>Override methods of ```Bio.PDB.MMCIF2Dict.MMCIF2Dict```</h3>
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
                print(token[1])
            except StopIteration:
                return  # for Python 3.7 and PEP 479
            self[token[1][0:5]] = token[1][5:]
            i = 0
            n = 0
            use = []
            for token in tokens:
                if token[1].lower() == "loop_":
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
                    if token[1].startswith("_") and (n == 0 or i % n == 0):
                        if i > 0:
                            loop_flag = False
                        else:
                            # Additional statement:focus on target data
                            if token[1] in useKeyList and token[0] == 0:
                                use.append(n)
                                self[token[1]] = []
                            keys.append(token[1])
                            n += 1
                            continue
                    else:
                        key_index = i % n
                        try:
                            if key_index in use:
                                self[keys[key_index]].append(token[1])
                        except Exception as e:
                            print(keys, key_index, use)
                            raise Exception(e)
                        i += 1
                        continue
                if key is None:
                    # Additional statement:focus on target data
                    if token[1] in useKeyList and token[0] == 0:
                        key = token[1]
                else:
                    # Always returns a list
                    self[key] = [token[1]]
                    key = None

    def _splitline(self, line):
        # See https://www.iucr.org/resources/cif/spec/version1.1/cifsyntax for the syntax
        in_token = False
        # quote character of the currently open quote, or None if no quote open
        quote_open_char = None
        start_i = 0
        for (i, c) in enumerate(line):
            if c in self.whitespace_chars:
                if in_token and not quote_open_char:
                    in_token = False
                    # yield line[start_i:i]
                    yield start_i, line[start_i:i]
            elif c in self.quote_chars:
                if not quote_open_char:
                    if in_token:
                        raise ValueError("Opening quote in middle of word: " + line)
                    quote_open_char = c
                    in_token = True
                    start_i = i + 1
                elif c == quote_open_char and (i + 1 == len(line) or line[i + 1] in self.whitespace_chars):
                    quote_open_char = None
                    in_token = False
                    # yield line[start_i:i]
                    yield start_i, line[start_i:i]
            elif c == "#" and not in_token:
                # Skip comments. "#" is a valid non-comment char inside of a
                # quote and inside of an unquoted token (!?!?), so we need to
                # check that the current char is not in a token.
                return
            elif not in_token:
                in_token = True
                start_i = i
        if in_token:
            # yield line[start_i:]
            yield start_i, line[start_i:]
        if quote_open_char:
            raise ValueError("Line ended with quote open: " + line)

    def _tokenize(self, handle):
        for line in handle:
            if line.startswith("#"):
                continue
            elif line.startswith(";"):
                # The spec says that leading whitespace on each line must be
                # preserved while trailing whitespace may be stripped.  The
                # trailing newline must be stripped.
                token_buffer = [line[1:].rstrip()]
                for line in handle:
                    line = line.rstrip()
                    if line == ";":
                        break
                    token_buffer.append(line)
                # yield "\n".join(token_buffer)
                yield 1, "\n".join(token_buffer)
            else:
                for token in self._splitline(line.strip()):
                    yield token


class MMCIF2Dfrm(Unit):
    """
    # Convert MMCIF data into a DataFrame

    ## Key Feactures of PDB

    * Which method & resolution (Filter)
    * Date (Filter)
    * Whether ATOM ONLY (Filter)
    * Whether contain UNK residue (Filter)
    * Whether contain Nucleotide (Filter)
    * Mo/Ho/He (Filter)
    * Protein Chains that have length more than 20 aa (Filter)
    * Nucleotide Chains that have length more than 5 bases (Filter)
    * BioUnit Component
    * Mutations
    * Ligands
    * Missing Region of each chain
    * ...

    ## Ability of Building/Updating DataSet
    """

    FUNC_LI_DI = []
    FUNC_LI_DF = []
    pdb_path_li = []
    protein_filter = {
        "coordinates_len": (20, "gt"),
        "_pdbx_coordinate_model.type": (None, "notNull"),
        "UNK_ALL_IN_CHAIN": (False, "eq"),
        "contains_unk_in_chain_pdb": (False, "eq"),
        "proteinChainCount": (26, "le"),
        "method": (['SOLUTION NMR', 'X-RAY DIFFRACTION'], "isIn")
    }

    @property
    def default_use_keys(self):
        use_keys = []
        use_col_li = ['COMMON_COL', 'ENTITY_COL', 'TYPE_COL', 'SEQRES_COL', 'LIGAND_COL', 'BIOASS_COL', 'NON_POLY_COL', 'COORDINATE_MODEL_COL']
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

    def check_mmcif_file(self, pdb_list, processes=4, maxSleep=3):
        def find_unDownloaded_file(pdbId):
            for path in MMCIF_FILE_FOLDER['MMCIF_OLD_FOLDER'] + [MMCIF_FILE_FOLDER['MMCIF_NEW_FOLDER']]:
                old_path = '%s%s.cif' % (path, pdbId)
                if os.path.exists(old_path):
                    MMCIF2Dfrm.pdb_path_li.append(old_path)
                    return False
            MMCIF2Dfrm.pdb_path_li.append(old_path)
            return True

        unDownload = list(filter(find_unDownloaded_file, pdb_list))
        mpw = MPWrapper(MMCIF_FILE_FOLDER['MMCIF_NEW_FOLDER'], processes=processes, maxSleep=maxSleep)
        # @retry(stop_max_attempt_number=3, wait_fixed=1000)
        mpw.http_retrive(unDownload)

    def update_mmcif_result(self, rawOutputPath, handledOutputPath, chunksize=100, finished=[]):
        mmcif_file_li = []
        for path in self.pdb_path_li:
            if path[-8:-4] not in finished:
                mmcif_file_li.append(path)
        for i in range(0, len(mmcif_file_li), chunksize):
            chunk_li = mmcif_file_li[i:i+chunksize]
            chunk_df = self.mmcif_dict2dfrm(chunk_li, outputPath=rawOutputPath)
            self.handle_mmcif_dfrm(chunk_df, outputPath=handledOutputPath)

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
        """# Decorator to add new dispatch functions."""
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
        coordinates_model_key = DEFAULT_COLS['SEQRES_COL'][8]
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

            for col in DEFAULT_COLS['SEQRES_COL'][4:8]:
                info_dict[col][i] = [';'.join(
                    MMCIF2Dfrm.get_index(strand_id_index, info_dict[col][i], j))
                    for j in range(len(strand_id_index))]

            new_comodel_li = []
            for ele in info_dict[coordinates_model_key][i]:
                if ele not in new_comodel_li:
                    new_comodel_li.append(ele)
            info_dict[coordinates_model_key][i] = new_comodel_li

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

    @dispatch_on_set(DEFAULT_COLS['BIOASS_COL']+['_pdbx_poly_seq_scheme.asym_id'], func_li=FUNC_LI_DF)
    def handle_bioas_df(df):
        '''
        ### Deal with Biological Assemblies in PDB
        #### Example
        Key|Value
        -|-
        pdb_id|3A6P
        _entity.id|[1, 2, 3, 4, 5, 6, 7]
        _entity_poly.entity_id|[1, 2, 3, 4, 5]
        _entity_poly.pdbx_strand_id|[A,F, B,G, C,H, D,I, E,J]
        _entity_poly.type|[polypeptide(L), polypeptide(L), polypeptide(L)...
        _pdbx_struct_assembly_gen.assembly_id|[1, 2]
        _pdbx_struct_assembly_gen.oper_expression|[1, 1]
        _pdbx_struct_assembly_gen.asym_id_list|[A,B,C,D,E,K,L, F,G,H,I,J,M,N]
        _pdbx_struct_assembly.oligomeric_count|[5, 5]
        ---
        ```_pdbx_struct_assembly_gen.asym_id_list``` -> [A,B,C,D,E,~~K,L,~~ F,G,H,I,J,~~M,N~~]
        ```~~_pdbx_poly_seq_scheme.pdb_strand_id~~ [wrong], _pdbx_poly_seq_scheme.asym_id [correct]```
        > For more examples, please goto [here](https://naturegeorge.github.io/BioinforResearch/md_for_MMCIF2Dfrm.html "Link")
        '''

        def getBioUnitInfo(au_id_li, au_asym_id_li, oli, pl_asym_id_li):
            au_dict = defaultdict(list)
            for index in range(len(au_id_li)):
                au_dict[au_id_li[index]].extend([chain for chain in au_asym_id_li[index].split(',') if chain in pl_asym_id_li])
            for index, key in enumerate(sorted(au_dict.keys())):
                au_dict[key] = (int(oli[index]), au_dict[key])
            return json.dumps(au_dict)

        df['bioUnit'] = df.apply(lambda x: getBioUnitInfo(
                x['_pdbx_struct_assembly_gen.assembly_id'],
                x['_pdbx_struct_assembly_gen.asym_id_list'],
                x['_pdbx_struct_assembly.oligomeric_count'],
                x['_pdbx_poly_seq_scheme.asym_id']
                ) if not isinstance(x['_pdbx_struct_assembly_gen.assembly_id'], float) else np.nan, axis=1)

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

    def handle_mmcif_dfrm(self, dfrm, outputPath=False):

        def get_sub_df(df, i, spe_col_li, common_col_li):
            try:
                a = pd.DataFrame({key: df.loc[i, key] for key in spe_col_li})
            except Exception:
                try:
                    a = pd.DataFrame({key: json.loads(df.loc[i, key].replace('\'', '"').replace('False', 'false').replace('True', 'true')) for key in spe_col_li})
                except Exception:
                    a = pd.DataFrame({key: [df.loc[i, key]] for key in spe_col_li})

            for common_col in common_col_li:
                da = df.loc[i, common_col]
                if isinstance(da, list):
                    a[common_col] = ','.join(df.loc[i, common_col])
                elif isinstance(da, (str, bool, np.bool_)) or pd.isna(da):
                    a[common_col] = da
                else:
                    print('get_sub_df(): WARNING: %s -> %s(%s)' % (common_col, da, type(da)))
                    a[common_col] = da
            return a

        def sub_handle_df(df, spe_col_li, common_col_li):
            df_li = []
            for i in df.index:
                df_li.append(get_sub_df(df, i, spe_col_li, common_col_li))
            return pd.concat(df_li, ignore_index=True)

        entity_poly_df = sub_handle_df(
            dfrm, DEFAULT_COLS['ENTITY_COL'] + ['mutation_num'], ['pdb_id'])
        type_poly_df = sub_handle_df(
            dfrm, DEFAULT_COLS['TYPE_COL'], ['pdb_id'])
        basic_df = sub_handle_df(dfrm, DEFAULT_COLS['SEQRES_COL'] + ['UNK_ALL_IN_CHAIN'],
            [
                'pdb_id', 'method', 'initial_version_time',
                'newest_version_time', 'resolution', 'pdb_contain_chain_type',
                'contains_unk_in_chain_pdb', 'pdb_type_MMCIF', 'bioUnit'
            ]
            + DEFAULT_COLS['COMMON_COL'][3:5]
            + DEFAULT_COLS['BIOASS_COL']
            + DEFAULT_COLS['NON_POLY_COL'])
        ligand_df = sub_handle_df(
            dfrm, DEFAULT_COLS['METAL_LIGAND_COL'], ['pdb_id'])

        new_type_poly_df = type_poly_df.drop(DEFAULT_COLS['TYPE_COL'][1], axis=1).join(
            type_poly_df[DEFAULT_COLS['TYPE_COL'][1]].str.split(',', expand=True).stack().reset_index(level=1, drop=True).rename('chain_id'))

        coordinates_model_df = sub_handle_df(dfrm, DEFAULT_COLS['COORDINATE_MODEL_COL'], ['pdb_id'])
        coordinates_model_df.rename(columns={'_pdbx_coordinate_model.asym_id': 'asym_id'}, inplace=True)

        entity_poly_df.rename(columns={
                              '_entity.pdbx_mutation': 'mutation_content', '_entity.id': 'entity_id'}, inplace=True)
        new_type_poly_df.rename(columns={
                                '_entity_poly.entity_id': 'entity_id', '_entity_poly.type': 'protein_type'}, inplace=True)
        basic_df.rename(
            columns={'_pdbx_poly_seq_scheme.pdb_strand_id': 'chain_id', '_pdbx_poly_seq_scheme.asym_id': 'asym_id'}, inplace=True)

        ligand_df.rename(
            columns={DEFAULT_COLS['METAL_LIGAND_COL'][0]: 'asym_id'}, inplace=True)

        if ligand_df['asym_id'].isnull().sum() == len(ligand_df):
            basic_df[DEFAULT_COLS['METAL_LIGAND_COL'][1]] = np.nan
            df_1 = basic_df
        else:
            df_1 = pd.merge(basic_df, ligand_df, how='left')

        if coordinates_model_df['asym_id'].isnull().sum() == len(coordinates_model_df):
            df_1[DEFAULT_COLS['COORDINATE_MODEL_COL'][1]] = np.nan
            df_1_1 = df_1
        else:
            df_1_1 = pd.merge(df_1, coordinates_model_df, how='left')

        df_2 = pd.merge(new_type_poly_df, df_1_1, how='left')
        df_3 = pd.merge(df_2, entity_poly_df, how='left')

        df_3['metal_ligand_num'] = df_3.apply(lambda x: str(x['metal_ligand_content']).count(
            '),') + 1 if not isinstance(x['metal_ligand_content'], float) else 0, axis=1)
        df_3['Modification_num'] = df_3.apply(lambda x: x['_pdbx_poly_seq_scheme.mon_id'].count(
            'X') if not isinstance(x['_pdbx_poly_seq_scheme.mon_id'], float) else np.nan, axis=1)
        df_3['seqres_len'] = df_3.apply(lambda x: len(x['_pdbx_poly_seq_scheme.mon_id']) if not isinstance(
            x['_pdbx_poly_seq_scheme.mon_id'], float) else np.nan, axis=1)
        df_3['coordinates_len'] = df_3.apply(lambda x: len(x['_pdbx_poly_seq_scheme.pdb_mon_id'].replace(
            '?', '')) if not isinstance(x['_pdbx_poly_seq_scheme.pdb_mon_id'], float) else np.nan, axis=1)

        def find_charIndex_fun(s, char): return [x for x in range(
            s.find(char), len(s)) if s[x] == char]
        df_3['Modification_index'] = df_3.apply(lambda x: find_charIndex_fun(
            x['_pdbx_poly_seq_scheme.mon_id'], 'X') if not isinstance(x['_pdbx_poly_seq_scheme.mon_id'], float) else np.nan, axis=1)
        df_3['mis_index'] = df_3.apply(lambda x: find_charIndex_fun(x['_pdbx_poly_seq_scheme.pdb_mon_id'], '?') if not isinstance(
            x['_pdbx_poly_seq_scheme.pdb_mon_id'], float) else np.nan, axis=1)

        df_3['mis_range'] = df_3.apply(lambda x: MMCIF2Dfrm.getInterval(
            x['mis_index']) if not isinstance(x['mis_index'], float) else np.nan, axis=1)
        df_3['resolution_score'] = df_3.apply(lambda x: MMCIF2Dfrm.handleResolution(x['resolution']), axis=1)

        col_list = ['pdb_id', 'chain_id', 'protein_type', 'coordinates_len']
        pro_chain_grouper = MMCIF2Dfrm.GroupER(
            col_list[0], ['polypeptide(L)', 'polypeptide(D)'], df_3, 'protein_chain_and_length')
        df_3['protein_chain_and_length'] = np.nan
        for i in df_3.index:
            pro_chain_grouper.check(df_3.loc[i, col_list[0]], (
                i, df_3.loc[i, col_list[-2]], df_3.loc[i, col_list[-1]], df_3.loc[i, col_list[1]]))
        pro_chain_grouper.output()

        if os.path.exists(outputPath):
            self.file_o(outputPath, df_3, mode='a+', header=False)
        else:
            self.file_o(outputPath, df_3)
        return df_3
```

```py
# @Date:   2019-08-16T23:24:17+08:00
# @Email:  1730416009@stu.suda.edu.cn
# @Filename: SIFTS_unit.py
# @Last modified time: 2019-10-19T19:28:02+08:00
import pandas as pd
import numpy as np
import json, wget, gzip, time, sys
from urllib import request
from itertools import combinations
from Bio import SeqIO
sys.path.append('./')
from Unit import Unit
from PdSeq_unit import PdSeqAlign


class SIFTS_unit(Unit):
    CONFIG = {
        'PDB_ID': 'pdb_id',
        'RAW_SIFTS_COLUMNS': [
            'pdb_id', 'chain_id', 'UniProt', 'identity', 'identifier',
            'pdb_start', 'pdb_end', 'unp_start', 'unp_end',
            'is_canonical', 'start', 'end', 'entity_id', 'struct_asym_id'],
        'DOWNLOAD_FOLDER': '../../data/sifts_files/',
        'UNP_LIST_PATH': '../../data/sifts_files/sifts_uniprot_list.tsv',
        'ELE_LIST': ["coordinates_len", "mappedOut", "metal_ligand_num", "delHT_MissingNum", "if_2", "if_1"],
        'INI_LIST': [1, 1, 2, 3, 4, 1, 2, 3, 4, 2, 3, 4, 2, 3, 2],
        'SEG_SIFTS_MAP': {
              "SP_PRIMARY": "UniProt",  # (Canonical)
              "RES_BEG": "pdb_start",
              "RES_END": "pdb_end",
              "PDB_BEG": "residue index(start) in pdb",
              "PDB_END": "residue index(end) in pdb",
              "SP_BEG": "unp_start",
              "SP_END": "unp_end",
              "PDB": 'pdb_id',
              "CHAIN": 'chain_id',
              },
        'SCORE_COL': 'BS',
        'PDB_SELECT_COL': 'pdb_chain_select',
        'PDB_RANK_COL': 'pdb_chain_rank',
        'PDB_RANK_LIST': ['BS', 'ne_resolution_score', 'initial_version_time'],
        'PDB_RANK_FORMAT': '%d-%d-%d',
        'PDB_INIT_GROUPBY_LIST': ['UniProt'],
        'PDB_SELECT_RANGE_NAME': 'seg_unp_range',
        'PDB_SEQRES_COL': '_pdbx_poly_seq_scheme.mon_id',

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
            print('getSiftsInfo(): Start to get the pdb info from SIFTS.', pdbId)
            url = ('http://www.ebi.ac.uk/pdbe/api/mappings/all_isoforms/%s') % (pdbId)
            try:
                req = request.Request(url)
                page = request.urlopen(req).read()
                page = page.decode('utf-8')
                info = json.loads(page)[pdbId]['UniProt']
                if info:
                    df = order_SIFTS_info(pdbId, info)
                    if not isinstance(df, bool):
                        df.to_csv(outputPath, index=False, sep='\t', mode='a+', header=False)  # Attention: 初始Build dataset时需指定header
            except Exception as e:
                print('PDB: [', pdbId, ']\n', e)
                fail_list.append(pdbId)
            current += 1
            print('getSiftsInfo(): End a circle.[', pdbId, '] current:', current, 'ALL:', allPDB)
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
            return {'pdb_set': set(pdb_list) & related_pdb, 'unp_set': set(dfrm['SP_PRIMARY'])}
        else:
            return {'pdb_set': set(pdb_list), 'unp_set': set(dfrm['SP_PRIMARY'])}

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
        unpLen_dfrm.rename(columns={'Length': 'UNP_len'}, inplace=True)
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
            if not isinstance(head, float) and not isinstance(tail, float):  # WARNING
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
                    mat[i][j] = com_dict.get((i, j), 1)
            x = mat.values
            y = 1/x
            value, vector = np.linalg.eig(np.triu(x.T).T + np.triu(y.T) - np.diag(x.diagonal()))
            max_val = np.max(value)
            index = list(value).index(max_val)
            max_vector = vector[:, index]
            select_vector = np.real(max_vector/np.linalg.norm(max_vector, ord=1))
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
            li = json.loads(sli.replace('(', '[').replace(')', ']'))
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
            lambda x: x['Modification_num'] + x['mutation_num'] if not isinstance(x['Modification_num'], float) else x['mutation_num'], axis=1)

        '''sifts_dfrm['metal_count'] = sifts_dfrm.apply(lambda x: len(x['metal_list']) if not isinstance(x['metal_list'], float) else 0, axis=1)'''

        select_vector = get_weight()
        calScore(sifts_dfrm, SIFTS_unit.CONFIG['SCORE_COL'], SIFTS_unit.CONFIG['ELE_LIST'], select_vector)

        # Add resolution-related socre
        deal_rs = lambda x: -float(x.split(',')[0]) if ',' in x else -1200

        def deal_reso_score(score):
            try:
                if isinstance(score, str):
                    if len(set(score) & set('?,')) == 0:
                        return -float(score)
                    else:
                        return deal_rs(score)
                else:
                    return -score
            except Exception:
                return -1000

        sifts_dfrm['ne_resolution_score'] = sifts_dfrm.apply(lambda x: deal_reso_score(x['resolution_score']), axis=1)

        # Find Useful chain
        sifts_dfrm['pdb_SIFTS_useful_chain_num'] = np.nan
        for i, j in sifts_dfrm.groupby([SIFTS_unit.CONFIG['PDB_ID']]):
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
            return in_df.drop(columns=range_info_col+['PDB_BEG', 'PDB_END']).drop_duplicates(subset=group_info_col, keep='last')

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

    def update_range_info_SIFTS(self, unp_fasta_files_path, new_range_cols=('new_sifts_unp_range', 'new_sifts_pdb_range'), sifts_df=False, sifts_filePath=False, outputPath=False):
        sifts_dfrm = self.file_i(sifts_filePath, sifts_df, ('sifts_filePath', 'sifts_df'))
        focus = ['Deletion', 'Insertion & Deletion']
        focus_index = sifts_dfrm[sifts_dfrm['sifts_range_tage'].isin(focus)].index

        da1 = []
        da2 = []
        pdSeqAligner = PdSeqAlign()
        for index in focus_index:
            pdbSeq = sifts_dfrm.loc[index, self.CONFIG['PDB_SEQRES_COL']]
            unpSeqOb = SeqIO.read(unp_fasta_files_path % sifts_dfrm.loc[index, 'UniProt'], "fasta")
            da = pdSeqAligner.makeAlignment_align(unpSeqOb.seq, pdbSeq)
            da1.append(da[0])
            da2.append(da[1])

        df = pd.DataFrame({new_range_cols[0]: da1, new_range_cols[1]: da2}, index=focus_index)
        new_sifts_df = pd.merge(sifts_dfrm, df, left_index=True, right_index=True, how='left')
        new_sifts_df[new_range_cols[0]] = new_sifts_df.apply(lambda x: x['sifts_unp_range'] if isinstance(x[new_range_cols[0]], float) else x[new_range_cols[0]], axis=1)
        new_sifts_df[new_range_cols[1]] = new_sifts_df.apply(lambda x: x['sifts_pdb_range'] if isinstance(x[new_range_cols[1]], float) else x[new_range_cols[1]], axis=1)

        self.file_o(outputPath, new_sifts_df)
        return new_sifts_df

    def select_PDB_SIFTS(self, groupby_list,
                               select_col=CONFIG['PDB_SELECT_COL'],
                               rank_col=CONFIG['PDB_RANK_COL'],
                               rank_list=CONFIG['PDB_RANK_LIST'],
                               rank_format=CONFIG['PDB_RANK_FORMAT'],
                               range_name=CONFIG['PDB_SELECT_RANGE_NAME'],
                               sifts_df=False, sifts_filePath=False, outputPath=False):
        sifts_dfrm = self.file_i(sifts_filePath, sifts_df, ('sifts_filePath', 'sifts_df'))
        sifts_dfrm[select_col] = False
        sifts_dfrm[rank_col] = np.nan

        for i, j in sifts_dfrm.groupby(groupby_list):
            SIFTS_unit.selectChain(
                j, sifts_dfrm, rank_list, rank_col, rank_format, range_name, select_col, 0.3, 0.2)

        self.file_o(outputPath, sifts_dfrm)
        return sifts_dfrm

    def find_mo_SIFTS(self, groupby_list, sifts_df=False, sifts_filePath=False, outputPath=False):
        sifts_dfrm = self.file_i(sifts_filePath, sifts_df, ('sifts_filePath', 'sifts_df'))
        constraint_dict = self.ConstraintDict(
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

        mo_fakeHoHe_df = self.ConstraintDict.addConstraintToDf(sifts_dfrm, constraint_dict)
        self.file_o(outputPath, mo_fakeHoHe_df)
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
            dfrm[unp_muta_col] = np.nan  # 'Mutation_Uniprot'
            for _, groupData in dfrm.groupby(groupby_list):  # ['pdb_id', 'chain_id', 'iso_id']
                coor_list = groupData.loc[groupData.index[0], 'pdb_ins_position'].split(';')  # pdb_ins_seqres_position
                sifts_unp_range = groupData.loc[groupData.index[0], 'sifts_unp_range']
                sifts_pdb_range = groupData.loc[groupData.index[0], 'sifts_pdb_range']
                dfrm.loc[groupData.index, unp_muta_col] = groupData.apply(
                    lambda x: getUniprotMutaSite(
                        x[pdb_muta_col], coor_list, sifts_unp_range, sifts_pdb_range
                        ), axis=1)

        sifts_dfrm = self.file_i(sifts_filePath, sifts_df, ('sifts_filePath', 'sifts_df'))
        mapMutaFromPDBToUniprot(sifts_dfrm)
        self.file_o(outputPath, sifts_dfrm)
        return sifts_dfrm

    def map_muta_from_unp_to_pdb(x, muta_col, unp_range_col, pdb_range_col, error_li, addInscode=True):
        sub_error_li = []
        muta_li = x[muta_col]
        if isinstance(muta_li, str):
            muta_li = json.loads(muta_li.replace('\'', '"'))
        unp_range = json.loads(x[unp_range_col])
        pdb_range = json.loads(x[pdb_range_col])

        auth_seq_li = x['_pdbx_poly_seq_scheme.auth_seq_num'].split(';')
        com_seq_li = x['_pdbx_poly_seq_scheme.pdb_seq_num'].split(';')
        inscode_seq_li = x['_pdbx_poly_seq_scheme.pdb_ins_code'].split(';')
        found_muta_1 = x['mutation_content'].split(',')
        found_muta_2 = [i[1:-1] for i in found_muta_1]

        pdb_li = []
        unp_li = []
        new_muta_site = []
        for ran in pdb_range:
            pdb_li.extend(list(range(ran[0], ran[1]+1)))
        for ran in unp_range:
            unp_li.extend(list(range(ran[0], ran[1]+1)))
        if addInscode:
            for muta in muta_li:
                try:
                    seqresSite = pdb_li[unp_li.index(int(muta[1:-1]))]
                except ValueError:
                    new_muta_site.append('#')
                    sub_error_li.append('Unmapped #: %s' % muta)
                    continue
                except IndexError:
                    new_muta_site.append('$')
                    sub_error_li.append('Unmapped $: %s' % muta)
                    continue
                try:
                    seq_aa = x['_pdbx_poly_seq_scheme.mon_id'][seqresSite-1]
                except IndexError:
                    print(x['pdb_id'], x['UniProt'], x['new_sifts_pdb_range'], x['new_sifts_unp_range'])
                    print(x['_pdbx_poly_seq_scheme.mon_id'], seqresSite)
                    raise IndexError
                ref_aa = muta[0]
                inscode = inscode_seq_li[seqresSite-1]

                # Analyse the Situation of unsuccessful mapping
                if seq_aa != ref_aa:
                    try:
                        found_muta = found_muta_1[found_muta_2.index(com_seq_li[seqresSite-1])]
                    except ValueError:
                        found_muta = False

                    error_info = []
                    if found_muta:
                        # May because the mutation that already found, check the muta aa
                        if found_muta[0] == ref_aa:
                            error_info.append('EntityMutation: %s' % found_muta)
                    if seq_aa == 'X':
                        error_info.append('ModifiedResidue: %s%s%s' % (seq_aa, com_seq_li[seqresSite-1], inscode))
                    else:
                        error_info.append('PossibleMutation: %s%s%s' % (seq_aa, com_seq_li[seqresSite-1], inscode))
                    sub_error_li.append('' + ','.join(error_info))
                else:
                    sub_error_li.append('Safe')

                # Check inscode
                if inscode != '.':
                    new_muta_site.append('%s%s' % (auth_seq_li[seqresSite-1], inscode))
                else:
                    new_muta_site.append(auth_seq_li[seqresSite-1])
        else:
            for muta in muta_li:
                try:
                    seqresSite = pdb_li[unp_li.index(int(muta[1:-1]))]
                except ValueError:
                    new_muta_site.append('#')
                    continue
                except IndexError:
                    new_muta_site.append('$')
                    continue

                seq_aa = x['_pdbx_poly_seq_scheme.mon_id'][seqresSite-1]
                ref_aa = muta[0]

                # Analyse the Situation of unsuccessful mapping
                if seq_aa != ref_aa:
                    try:
                        found_muta = found_muta_1[found_muta_2.index(com_seq_li[seqresSite-1])]
                    except ValueError:
                        found_muta = False

                    error_info = []
                    if found_muta:
                        # May because the mutation that already found, check the muta aa
                        if found_muta[0] == ref_aa:
                            error_info.append('EntityMutation: %s' % found_muta)
                    if seq_aa == 'X':
                        error_info.append('ModifiedResidue: %s%s%s' % (seq_aa, com_seq_li[seqresSite-1], inscode))
                    else:
                        error_info.append('PossibleMutation: %s%s%s' % (seq_aa, com_seq_li[seqresSite-1], inscode))
                    error_li.append('' + ','.join(error_info))
                else:
                    error_li.append('Safe')

                new_muta_site.append(auth_seq_li[seqresSite-1])

        error_li.append(sub_error_li)
        return new_muta_site

    def map_muta_from_pdb_to_unp(x, muta_col, unp_range_col, pdb_range_col, error_li, unp_fasta_files_path):
        sub_error_li = []
        muta_li = x[muta_col]
        if isinstance(muta_li, str):
            muta_li = json.loads(muta_li.replace('\'', '"'))
        unp_range = json.loads(x[unp_range_col])
        pdb_range = json.loads(x[pdb_range_col])

        unpSeqOb = SeqIO.read(unp_fasta_files_path % x['UniProt'], "fasta")
        unpSeq = unpSeqOb.seq

        # auth_seq_li = x['_pdbx_poly_seq_scheme.auth_seq_num'].split(';')
        com_seq_li = x['_pdbx_poly_seq_scheme.pdb_seq_num'].split(';')
        inscode_seq_li = x['_pdbx_poly_seq_scheme.pdb_ins_code'].replace('.', '').split(';')
        res_li = ['%s%s' % x for x in zip(com_seq_li, inscode_seq_li)]
        found_muta_1 = x['mutation_content'].split(',')
        found_muta_2 = [i[1:-1] for i in found_muta_1]

        pdb_li = []
        unp_li = []
        new_muta_site = []
        for ran in pdb_range:
            pdb_li.extend(list(range(ran[0], ran[1]+1)))
        for ran in unp_range:
            unp_li.extend(list(range(ran[0], ran[1]+1)))
        for muta in muta_li:
            # muta = auth_seq_li[muta-1]
            muta_seqres_index = res_li.index(muta[1:-1])
            try:
                # seqresSite = pdb_li[unp_li.index(int(muta[1:-1]))]
                uniprotSite = unp_li[pdb_li.index(muta_seqres_index+1)]
            except ValueError:
                new_muta_site.append('#')
                sub_error_li.append('Unmapped #: %s' % muta)
                continue
            except IndexError:
                new_muta_site.append('$')
                sub_error_li.append('Unmapped $: %s' % muta)
                continue

            seq_aa = unpSeq[uniprotSite-1]
            ref_aa = muta[0]
            # inscode = inscode_seq_li[seqresSite-1]

            # Analyse the Situation of unsuccessful mapping
            if seq_aa != ref_aa:
                try:
                    found_muta = found_muta_1[found_muta_2.index(muta[1:-1])]
                except ValueError:
                    found_muta = False

                error_info = []
                if found_muta:
                    # May because the mutation that already found, check the muta aa
                    if found_muta[0] == ref_aa:
                        error_info.append('EntityMutation: %s' % found_muta)
                if seq_aa == 'X':
                    error_info.append('ModifiedResidue: %s' % (seq_aa))
                else:
                    error_info.append('PossibleMutation: %s' % (seq_aa))
                sub_error_li.append('' + ','.join(error_info))
            else:
                sub_error_li.append('Safe')

            new_muta_site.append(uniprotSite)

        error_li.append(sub_error_li)
        return new_muta_site

    def get_muta_interactions_from_sifts(self, muta_li, score_df, handle_mmcif_df, target_col='Target_Mutation_unp'):
        def getCompoPDBIndex(dfrm, target, interactor):
            target_df = dfrm[dfrm['UniProt'] == target]
            target_set = target_df['pdb_id'].drop_duplicates()
            if target != interactor:
                interactor_df = dfrm[dfrm['UniProt'] == interactor]
                interactor_set = interactor_df['pdb_id'].drop_duplicates()
                result_set = pd.merge(target_set, interactor_set)
                return result_set['pdb_id'], target_df[target_df['pdb_id'].isin(result_set['pdb_id'])].index | interactor_df[interactor_df['pdb_id'].isin(result_set['pdb_id'])].index
            else:
                return target_set, target_df[target_df['pdb_id'].isin(target_set)].index

        def is_PDB_pure(sifts_df, mmcif_df, pdb_li):
            # nucle_key=['chain_id', 'protein_type', 'coordinates_len', '_pdbx_poly_seq_scheme.auth_mon_id']
            if isinstance(pdb_li, pd.DataFrame):
                pdb_li = pdb_li['pdb_id']
            info_di = {}
            nucle_di = {}
            for pdb in pdb_li:
                sifts_entityIdSet = set(sifts_df[sifts_df['pdb_id'] == pdb]['entity_id'])
                temp = mmcif_df[mmcif_df['pdb_id'] == pdb]
                mmcif_entityIdSet = set(temp['entity_id'])
                if mmcif_entityIdSet > sifts_entityIdSet:
                    info_di[pdb] = False
                    # nucle_df = temp[~temp['protein_type'].isin(['polypeptide(L)', 'polypeptide(D)'])][nucle_key]
                    nucle_se = temp[~temp['protein_type'].isin(['polypeptide(L)', 'polypeptide(D)'])]['chain_id']
                    # nucle_di[pdb] = nucle_df.to_dict('records')
                    nucle_di[pdb] = nucle_se.to_dict()
                else:
                    info_di[pdb] = True
                    nucle_di[pdb] = {}
            return info_di, nucle_di

        usecols = ['pdb_id', 'chain_id', 'UniProt', 'identity', 'identifier',
                    'BS', 'ne_resolution_score', 'initial_version_time',
                    'new_sifts_unp_range','new_sifts_pdb_range',
                    'mutation_content', 'method', 'pdb_type_MMCIF',
                    'pdb_contain_chain_type','bioUnit',
                    '_pdbx_poly_seq_scheme.auth_seq_num',
                    '_pdbx_poly_seq_scheme.pdb_ins_code',
                    '_pdbx_poly_seq_scheme.mon_id',
                    '_pdbx_poly_seq_scheme.pdb_seq_num']
        copo_df_li = []
        warn_li = []

        for target, interactor in muta_li.index:
            pdb_se, index = getCompoPDBIndex(score_df, target, interactor)
            if len(index) == 0:
                warn_li.append((target, interactor))
                continue

            info_di, nucle_di = is_PDB_pure(score_df.loc[index], handle_mmcif_df, pdb_se)
            copo_df = score_df.loc[index, usecols].copy()
            copo_df['pure'] = copo_df.apply(lambda x: info_di[x['pdb_id']], axis=1)
            copo_df['nonProtein'] = copo_df.apply(lambda x: nucle_di[x['pdb_id']], axis=1)

            muta_content = muta_li[target][interactor]
            copo_df[target_col] = copo_df.apply(lambda x: muta_content if x['UniProt'] == target else np.nan, axis=1)
            muta_info_li = []
            copo_df[target_col] = copo_df.apply(lambda x: SIFTS_unit.map_muta_from_unp_to_pdb(x, target_col, 'new_sifts_unp_range', 'new_sifts_pdb_range', muta_info_li) if not isinstance(x['new_sifts_pdb_range'], float) and not isinstance(x[target_col], float) else np.nan, axis=1)
            copo_df['muta_map_info'] = pd.Series(muta_info_li, index=copo_df.dropna(subset=[target_col]).index)
            copo_df['compo'] = "%s, %s" % (target, interactor)
            copo_df_li.append(copo_df)

        return pd.concat(copo_df_li), pd.DataFrame(warn_li, columns=['Target_UniprotID', 'Interactor_UniprotID'])

    def get_unp_len_from_fasta(self, unp, unp_fasta_files_path):
        unpSeqOb = SeqIO.read(unp_fasta_files_path % unp, "fasta")
        return len(unpSeqOb.seq)
```

```py
# @Date:   2019-08-16T23:19:34+08:00
# @Email:  1730416009@stu.suda.edu.cn
# @Filename: Unit.py
# @Last modified time: 2019-09-12T10:37:28+08:00
import pandas as pd
import numpy as np
import json



class Unit:

    SEQ_DICT = {
        "GLY": "G", "ALA": "A", "SER": "S", "THR": "T", "CYS": "C", "VAL": "V", "LEU": "L",
        "ILE": "I", "MET": "M", "PRO": "P","PHE": "F", "TYR": "Y", "TRP": "W", "ASP": "D",
        "GLU": "E", "ASN": "N", "GLN": "Q", "HIS": "H", "LYS": "K", "ARG": "R", "HSD": "H",
        "A": "A", "G": "G", "C": "C", "T": "T", "U": "U", "DA": "A", "DT": "T",
        "DU": "U", "DC": "C", "DG": "G","DI":"I","?":"?","UNK":"!"}

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
                self.dfrm.loc[self.index, self.new_col] = str(self.content).replace('\'', '"')

        def check(self, tp, data): # (index, type, len, chain_id)
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
            self.aa_map = Unit.SEQ_DICT

        def multi_letter_convert_to_one_letter(self, a):
            return self.aa_map.get(a, 'X')

    def handleResolution(resolution):
        float_fun = lambda x: float(x) if x not in '?.' else 1000
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
        # overlap_list = []
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
                    # overlap = getOverlap(temp_range, repreSet, overlap_type)
                    # overlap_list.append(overlap)
                    overlap = len(temp_range & repreSet)
                    temp_range_len = len(temp_range)
                    # ---------------------------------------------------------------
                    if (overlap/temp_range_len <= r1_cutoff) and ((temp_range_len - overlap)/len(repreSet) >= r2_cutoff):
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

    def file_i(self, path, df, va_tp, sep='\t'):
        try:
            dfrm = pd.read_csv(path, sep=sep, na_values=['', None], keep_default_na=False)
        except Exception:
            dfrm = df
            if not isinstance(dfrm, pd.DataFrame):
                raise Exception('Input: %s=pd.DataFrame() or %s=str' % va_tp)
        return dfrm

    def file_o(self, path, df, mode='w+', header=True):
        if path:
            df.to_csv(path, sep='\t', index=False, mode=mode, header=header)
            self.file_list.append(path)

    def __init__(self):
        self.file_list = []
```

```py
# @Date:   2019-10-24T23:35:42+08:00
# @Email:  1730416009@stu.suda.edu.cn
# @Filename: RetrievePDB.py
# @Last modified time: 2019-10-28T11:39:11+08:00
import wget
import gzip
import urllib
import ftplib
import shutil
import os
from collections import Iterable, Iterator
from multiprocessing.dummy import Pool
from time import sleep
from random import uniform
"""
Requirement:
1. Fast and Legal
2. Keep Update

Reference:
http://www.wwpdb.org/ftp/pdb-ftp-sites
http://zetcode.com/python/ftp/
"""
_FTP_SITE = ["RCSB", "PDBE", "PDBJ"]
_FTP_HEADER = "ftp://"
_RCSB_FTP = "ftp.rcsb.org"
_PDBE_FTP = "ftp.ebi.ac.uk"
_PDBJ_FTP = "ftp.pdbj.org"
# _FTP_HOST = {"RCSB": _RCSB_FTP, "PDBE": _PDBE_FTP, "PDBJ": _PDBJ_FTP}
_FTP_HOST = dict(zip(_FTP_SITE, [_RCSB_FTP, _PDBE_FTP, _PDBJ_FTP]))
_RCSB_DIVIDED = "pub/pdb/data/structures/divided"
_PDBE_DIVIDED = "pub/databases/pdb/data/structures/divided"
_PDBJ_DIVIDED = "pub/pdb/data/structures/divided"
# _DIVIDED_PATH = {"RCSB": _RCSB_DIVIDED, "PDBE": _PDBE_DIVIDED, "PDBJ": _PDBJ_DIVIDED}
_DIVIDED_PATH = dict(
    zip(_FTP_SITE, [_RCSB_DIVIDED, _PDBE_DIVIDED, _PDBJ_DIVIDED]))
_COMPLETE_TAGE = "226"  # "226-File successfully transferred", "226 Transfer complete"
_FORMAT_DICT = {
    "XML": ".xml",
    "mmCIF": ".cif",
    "pdb": ".pdb",
}
_RAW_FORMAT = {".pdb": ".ent"}
_RCSB_HTTP_VIEW = "https://files.rcsb.org/view/"
_RCSB_HTTP = "https://files.rcsb.org/download/"


def printList(list):
    string = ""
    for i in list:
        string += "{0}\n".format(i)
    print(string)


class RetrievePDB:
    """
    Retrieve PDB File

    PDB Archive Reference: wwPDB [1]_ [2]_ [3]_ [4]_, RCSB [5]_ [6]_, PDBe [7]_, PDBj [8]_
    FTP Site Reference: [9]_
    Code Reference: [10]_

    .. [1] http://www.wwpdb.org/
    .. [2] H.M. Berman, K. Henrick, H. Nakamura (2003) Announcing the worldwide Protein Data Bank Nature Structural Biology 10 (12): 980.
    .. [3] H.M. Berman, K. Henrick, H.Nakamura, J.L. Markley (2007) The Worldwide Protein Data Bank (wwPDB): Ensuring a single, uniform archive of PDB data Nucleic Acids Res. 35 (Database issue): D301-3.
    .. [4] wwPDB consortium. (2019) Protein Data Bank: the single global archive for 3D macromolecular structure data. Nucleic Acids Res 47: D520-D528 doi: 10.1093/nar/gky949.
    .. [5] http://www.rcsb.org/
    .. [6] H.M. Berman, J. Westbrook, Z. Feng, G. Gilliland, T.N. Bhat, H. Weissig, I.N. Shindyalov, P.E. Bourne. (2000) The Protein Data Bank Nucleic Acids Research, 28: 235-242.
    .. [7] https://www.ebi.ac.uk/pdbe/
    .. [8] https://pdbj.org/
    .. [9] http://www.wwpdb.org/ftp/pdb-ftp-sites
    .. [10] http://zetcode.com/python/ftp/

    Script:

    .. code-block:: python
        :linenos:

        ftpPDB = RetrievePDB("C:/Users/Nature/Downloads/PDBE/", ftpSite="PDBE", format="pdb")
        ftpPDB.ftp_retrieve(pdbs=['10z1', '10z2', '2xyn', '10z3'], remove=True)
        print(ftpPDB)
        printList(ftpPDB.getFail())
        ftpPDB.quick_ftp_retrieve('5js8', remove=False)
        ftpPDB.quick_http_retrieve('5js1', module="urllib", view=False, bioAssembly=1, remove=True)

    """

    class HandleIO:
        def __init__(self, handle):
            self.handle = handle

        def append(self, block):
            self.handle.write(block)

        def close(self):
            self.handle.close()

    def __init__(self, downloadPath, pdbs=None, ftpSite="PDBE", format="mmCIF"):
        self.setDownloadPath(downloadPath)
        self.setFormat(format)
        self.setFTPSite(ftpSite)
        self.pdbs = pdbs
        self.fail = []
        print(self)

    def __len__(self):
        if self.pdbs is None:
            return 0
        else:
            return len(self.pdbs)

    def __getitem__(self, x):
        path = os.path.join(self.downloadPath, "%s%s" % (x, self.tail))
        if os.path.exists(path):
            return path
        else:
            return None

    def __repr__(self):
        format = "%s: %s, "
        string = "RetrievePDB: {%s}"
        content = ""
        for key, value in self.__dict__.items():
            if key not in ['pdbs', 'fail']:
                content += format % (key, value)
        content += format % ("len(pdbs)", len(self))
        content += format % ("len(fail)", len(self.fail))
        return string % content

    def getFail(self):
        """
        Get the list of PDB ids that failed to fetch files
        """
        return self.fail

    def setFTPSite(self, ftpSite):
        """
        Set the value of PDB FTP site

        :param str ftpSite: The FTP site of retriving PDB files, default value: `RCSB`, {RCSB, PDBE, PDBJ}

        """
        if ftpSite not in _FTP_SITE:
            raise ValueError(
                "Illegal site name. Please select from %s" % _FTP_SITE)
        else:
            self.ftpSite = ftpSite
            self.host = _FTP_HOST[ftpSite]
            self.dividedPath = _DIVIDED_PATH[ftpSite]

    def setFormat(self, format):
        """
        Set the format of PDB file

        :param str format: The file format of PDB file, default value: `mmCIF`, {mmCIF, pdb, XML}

        """
        if format not in _FORMAT_DICT.keys():
            raise ValueError(
                "Illegal format name. Please select from mmCIF, pdb, XML")
        else:
            self.format = format
            self.tail = _FORMAT_DICT[format]
            self.raw_tail = _RAW_FORMAT.get(self.tail, self.tail)
            if format == "pdb":
                self.prefix = format
            else:
                self.prefix = ""

    def setDownloadPath(self, path):
        """Set the value of PDB file download path"""
        if os.path.isdir(path):
            self.downloadPath = path
        else:
            raise ValueError("Not a valid download path")

    def setPDBs(self, pdbs):
        """
        Set the value of PDB id(s)

        :param pdbs: single PDB id or PDB ids
        :type pdbs: Iterable or Iterator or str

        * ``str``: single PDB id
        * ``Iterable, Iterator``: PDB ids

        """
        if isinstance(pdbs, str):
            self.pdbs = [pdbs]
        elif isinstance(pdbs, Iterable):
            self.pdbs = sorted(pdbs, key=lambda x: x[1:3] + x[0] + x[3])
        elif isinstance(pdbs, Iterator):
            self.pdbs = pdbs
        elif pdbs is None:
            raise ValueError("pdbs should not be None. Please Specify pdbs!")
        else:
            raise ValueError("Invalid Input")

    def ftp_retrieve(self, **kwargs):
        """
        Retrieve PDB files via FTP Connection

        :param pdbs: single PDB id or PDB ids
        :param bool remove: whether remove the compressed file, default value: `True`
        :type pdbs: Iterable or Iterator or str

        """
        # Check PDB List
        pdbs = kwargs.get('pdbs', False)
        remove = kwargs.get('remove', True)
        if pdbs is not False:
            self.setPDBs(pdbs)
        elif self.pdbs is None:
            raise ValueError("Please Specify pdbs!")
        # Connect
        with ftplib.FTP(self.host) as ftp:
            print(ftp.getwelcome())
            ftp.login()  # anonymous account
            ftp.cwd("%s/%s" % (self.dividedPath, self.format))
            # Start to retrieve
            cur = ""
            for pdb in self.pdbs:
                pdb = pdb.lower()
                subPath = pdb[1:3]
                try:
                    if cur != subPath:
                        if cur == "":
                            ftp.cwd(subPath)
                        else:
                            ftp.cwd("../%s" % subPath)
                        cur = subPath
                    file_orig = "%s%s%s.gz" % (self.prefix, pdb, self.raw_tail)
                    # Start to download
                    filename = os.path.join(
                        self.downloadPath, "%s%s.gz" % (pdb, self.tail))
                    print("Downloading File: %s" % filename)
                    data = self.HandleIO(open(filename, 'w+b'))
                    res = ftp.retrbinary('RETR ' + file_orig, data.append)
                    data.close()
                    # Check Data Completeness
                    if not res.startswith(_COMPLETE_TAGE):
                        print('Download failed', res)
                        if os.path.isfile(filename):
                            os.remove(filename)
                        self.fail.append(pdb)
                        continue
                    self.decompression(filename, remove=remove)
                except ftplib.error_perm as e:
                    print('FTP error:', e)
                    if 'filename' in locals().keys():
                        data.close()
                        os.remove(filename)
                    self.fail.append(pdb)
                    continue
        print("FTP Closed")

    def decompression(self, path, extension=".gz", remove=True, outputPath=None):
        """
        Decompress gz file

        :param str path: The file path of the compressed file
        :param str extension: Compressed file tail, default value: `.gz`
        :param bool remove: Whether remove the compressed file, default value: `True`
        :param str outputPath: File safe path, default value: `None`
        """

        """
        with gzip.GzipFile(mode="rb", fileobj=open(path, 'rb')) as raw:
            with open(path[:-len(extension)], "wb") as file:
                file.write(raw.read())
        """
        if outputPath is None:
            outputPath = path[:-len(extension)]

        with gzip.open(path, 'rb') as raw:
            with open(outputPath, 'wb') as file:
                shutil.copyfileobj(raw, file)
        try:
            if remove:
                os.remove(path)
        except Exception as e:
            print(e)

    def quick_ftp_retrieve(self, pdb, remove=True):
        """
        Download PDB file via FTP with ``wget``

        *Download only one file at a time*

        :param str pdb: single PDB id
        :param bool remove: whether remove the compressed file, default value: `True`
        """
        pdb = pdb.lower()
        file_orig = "%s%s%s.gz" % (self.prefix, pdb, self.raw_tail)
        site = "%s%s/%s/%s/%s/%s" % (_FTP_HEADER, self.host,
                                     self.dividedPath, self.format, pdb[1:3], file_orig)
        path = os.path.join(self.downloadPath, "%s%s.gz" % (pdb, self.tail))
        print("Downloading File: %s" % path)
        try:
            wget.download(site, out=path)
        except urllib.error.URLError:
            print("Download failed")
            self.fail.append(pdb)
        self.decompression(path, remove=remove)

    def quick_http_retrieve(self, pdb, module="wget", view=False, bioAssembly="", extension=".gz", remove=True):
        """
        Download PDB file via HTTP with ``wget.download`` or ``urllib.request.urlopen``

        *RCSB Only*

        :param str pdb: single PDB id
        :param str module: For user to select the module, default value: `wget`, {"wget", "urllib"}
        :param bool view: Whether get the PDB data from the web pages
        :param str bioAssembly: The bioAssembly id, default value: `''`
        :param str extension: Compressed file tail, default value: `.gz`
        :param bool remove: whether remove the compressed file, default value: `True`

        """
        # Whether download from web page
        if view:
            fileName = "%s%s" % (pdb, self.tail)
            url = "%s%s" % (_RCSB_HTTP_VIEW, fileName)
        else:
            fileName = "{pdb}{tail}{bioAssembly}{extension}".format(
                pdb=pdb, tail=self.tail, bioAssembly=bioAssembly, extension=extension)
            url = "{site}{fileName}".format(site=_RCSB_HTTP, fileName=fileName)
        # r = requests.get(url)
        path = os.path.join(self.downloadPath, fileName)
        print("Downloading File: %s" % path)
        try:
            # Select module method: urllib
            if module != "wget":
                with open(path, 'w+b') as fw:
                    # fw.write(r.content)
                    text = urllib.request.urlopen(url).read()
                    fw.write(text)  # .decode('utf-8')
            # Select module method: wget
            else:
                wget.download(url, out=path)
        except urllib.error.URLError:
            print("Download failed")
            self.fail.append(pdb)
            return
        # Whether to decompress
        if not view and extension == '.gz':
            self.decompression(path, remove=remove)


class MPWrapper:
    """
    Multiprocessing wrapper for ``RetrievePDB``

    When there is a large number of PDB files to download, this class is helpful.
    But Need to be careful with the numbers of processes and the time of sleep.

    Script:

    .. code-block:: python
        :linenos:

        pdbs = ['1A02', '3KBZ', '3KC0', '3KC1', '3KMU', '3KMW', '3KYC', '3KYD', ...]
        mpw = MPWrapper("C:/Users/Nature/Downloads/")
        # fail = mpw.http_retrieve(pdbs)
        # fail = mpw.http_retrieve(pdbs, module="urllib")
        # fail = mpw.ftp_retrieve_wget(pdbs)
        fail = mpw.ftp_retrieve_batch(pdbs)
        printList(fail)

    :param str downloadPath: File folder of Downloaded PDB files
    :param int processes: Number of processes, default value: `3`
    :param int maxSleep: Max sleep time, default value: `3`
    :param str ftpSite: The FTP site of retriving PDB files, default value: `RCSB`, {RCSB, PDBE, PDBJ}
    :param str format: The file format of PDB file, default value: `mmCIF`, {mmCIF, pdb, XML}

    """

    def __init__(self, downloadPath, processes=3, maxSleep=3, ftpSite="RCSB", format="mmCIF"):
        self.setProcesses(processes, maxSleep)
        self.retrievePDB = RetrievePDB(
            downloadPath, ftpSite=ftpSite, format=format)

    def setProcesses(self, processes, maxSleep):
        """
        Set the value of ``processes`` and ``maxSleep``

        :param int processes: Number of processes
        :param int maxSleep: Max sleep time
        """
        if processes > 20:
            print("MPWrapper: Too many processes. Be careful !")
        self.processes = processes
        self.maxSleep = maxSleep

    def http_retrieve(self, pdbs, module="wget", view=False, bioAssembly="", extension=".gz", remove=True):
        """
        Retrieve PDB file via http with ``wget.download`` or ``urllib.request.urlopen``

        **HTTP SITE: RCSB ONLY**

        :param pdbs: An object containing the PDB ids that need to be download
        :param str module: For user to select the module, default value: `wget`, {"wget", "urllib"}
        :param bool view: Whether get the PDB data from the web pages, default value: `False`
        :param str bioAssembly: The bioAssembly id, default value: `''`
        :param str extension: Compressed file tail, default value: `.gz`
        :param bool remove: whether remove the compressed file, default value: `True`
        :type pdbs: Iterable or Iterator
        :return: A fail list that contains the PDB ids falied to download
        :rtype: list(str)
        """
        def register(pdb):
            stop = uniform(0, self.maxSleep)
            sleep(stop)
            self.retrievePDB.quick_http_retrieve(pdb, module=module, view=view, bioAssembly=bioAssembly, extension=extension, remove=remove)
            # print(pdb, stop)

        pool = Pool(processes=self.processes)
        pool.map(register, pdbs)
        return self.retrievePDB.getFail()

    def ftp_retrieve_wget(self, pdbs, remove=True):
        """
        Download PDB file via FTP with ``wget``

        :param pdbs: An object containing the PDB ids that need to be download
        :param bool remove: whether remove the compressed file, default value: `True`
        :type pdbs: Iterable or Iterator
        :return: A fail list that contains the PDB ids falied to download
        :rtype: list(str)
        """
        def register(pdb):
            stop = uniform(0, self.maxSleep)
            sleep(stop)
            self.retrievePDB.quick_ftp_retrieve(pdb, remove=remove)
            # print(pdb, stop)

        pool = Pool(processes=self.processes)
        pool.map(register, pdbs)
        return self.retrievePDB.getFail()

    def ftp_retrieve_batch(self, pdbs, remove=True, chunksize=100):
        """
        Retrieve PDB files via FTP Connection

        :param pdbs: An object containing the PDB ids that need to be download
        :param bool remove: whether remove the compressed file, default value: `True`
        :param int chunksize: the size of PDBs that query during a single FTP connection, default value: `100`
        :type pdbs: Iterable or Iterator
        :return: A fail list that contains the PDB ids falied to download
        :rtype: list(str)

        """
        assert isinstance(pdbs, Iterable), "pdbs should be an Iterable object in this function!"

        def register(chunk):
            sleep(uniform(0, self.maxSleep))
            self.retrievePDB.ftp_retrieve(pdbs=chunk, remove=remove)
            # print(chunk)

        chunks = [pdbs[i:i + chunksize]
                  for i in range(0, len(pdbs), chunksize)]
        pool = Pool(processes=self.processes)
        pool.map(register, chunks)
        return self.retrievePDB.getFail()

```

```py
# @Date:   2019-09-03T16:37:05+08:00
# @Email:  1730416009@stu.suda.edu.cn
# @Filename: PdSeq_unit.py
# @Last modified time: 2019-09-10T15:39:57+08:00
from Unit import Unit
import sys
from Bio import pairwise2
from Bio.SubsMat import MatrixInfo as matlist
from Bio import Align
import json
import functools
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
                segment1 = [i1 + 1, j1]
                segment2 = [i2 + 1, j2]
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
        self.alignment_count = 0

    def makeAlignment_pairwise2(self, seqa, seqb):
        if seqa == self.seqa and seqb == self.seqb:
            return self.range_a, self.range_b
        else:
            self.seqa = seqa
            self.seqb = seqb
            self.alns = pairwise2.align.globalds(
                seqa, seqb, matlist.blosum62, -10, -0.5)
            # aln_seqa, aln_seqb, score, begin, end = alns[0]
            indexRecord = IndexRecorder()
            a = list(map(indexRecord.check, self.alns[0][0]))
            indexRecord.reset()
            b = list(map(indexRecord.check, self.alns[0][1]))
            indexRecord.reset()
            mapped = list(
                filter(lambda x: x[0] != -1 and x[1] != -1, zip(a, b)))
            self.range_a = Unit.getInterval([i[0] for i in mapped])
            self.range_b = Unit.getInterval([i[1] for i in mapped])
            return self.range_a, self.range_b

    def makeAlignment_align(self, seqa, seqb):
        if not (seqa == self.seqa and seqb == self.seqb):
            self.seqa = seqa
            self.seqb = seqb
            alignments = self.aligner.align(seqa, seqb)
            result = PdSeqAlign.getAlignmentSegment(alignments[0])
            self.range_a, self.range_b = json.dumps(
                result[0]), json.dumps(result[1])

        return self.range_a, self.range_b

    @functools.lru_cache()
    def new_makeAlignment_align(self, seqa, seqb):
        alignments = self.aligner.align(seqa, seqb)
        result = PdSeqAlign.getAlignmentSegment(alignments[0])
        self.alignment_count += 1
        print("MAKE ALIGNMENT: %s" % self.alignment_count)
        return json.dumps(result[0]), json.dumps(result[1])

```

```py
# @Date:   2019-08-16T20:26:58+08:00
# @Email:  1730416009@stu.suda.edu.cn
# @Filename: UniProt_unit.py
# @Last modified time: 2019-10-27T01:38:26+08:00
import urllib.parse
import urllib.request
from retrying import retry
import pandas as pd
import numpy as np
from random import uniform
from time import sleep
import os, sys, re, json
from collections import Counter
sys.path.append('./')
from Unit import Unit


class UniProt_unit(Unit):
    '''
    Use UniProt ID Mapping API

    Key parameters in using UniProt ID Mapping API:

    .. code-block:: python
        :linenos:

        params = {
            'from': 'ACC+ID',
            'to': 'ACC',
            'format': 'tab',
            'columns': 'id,length,reviewed,comment(ALTERNATIVE%20PRODUCTS),feature(ALTERNATIVE%20SEQUENCE)'...
            'query': list_str,
        }

    .. csv-table:: Details of ID Abbreviation
        :header: "Abbreviation", "Name", "Direction"

        "ACC+ID", "UniProtKB AC/ID", "from"
        "ACC", "UniProtKB AC", "both"
        "ID", "UniProtKB ID", "both"
        "EMBL_ID", "EMBL/GenBank/DDBJ", "both"
        "REFSEQ_NT_ID", "RefSeq Nucleotide", "both"
        "P_REFSEQ_AC", "RefSeq Protein", "both"
        "PDB_ID", "PDB", "both"
        "ENSEMBL_TRS_ID", "Ensembl Transcript", "both"
        "ENSEMBL_ID", "Ensembl", "both"
        "...", "...", "..."

    Here is the `Reference`_.

    .. _Reference: https://www.uniprot.org/help/api_idmapping

    Following table show the details of rest params:

    .. csv-table:: Details of other parameters
        :header: "param", "Explanation", "Reference"

        "columns", "comma-separated list of column names", "https://www.uniprot.org/help/api_queries"
        "columns", "Lists the column names for programmatic (RESTful) access to tab-separated or Excel downloads of UniProtKB search results", "https://www.uniprot.org/help/uniprotkb_column_names"
        "format", "html | tab | xls | fasta | gff | txt | xml | rdf | list | rss", "https://www.uniprot.org/help/api_queries"
        "query", "query text/id(s)", "https://www.uniprot.org/help/text-search"

    '''

    URL = 'https://www.uniprot.org/uploadlists/'
    COLUMNS = ['id', 'length', 'reviewed', 'comment(ALTERNATIVE%20PRODUCTS)', 'feature(ALTERNATIVE%20SEQUENCE)', 'genes', 'organism', 'sequence', 'protein%20names']
    COLUMN_DICT = {
                    'id': 'Entry', 'length': 'Length', 'reviewed': 'Status',
                    'comment(ALTERNATIVE%20PRODUCTS)': 'Alternative products (isoforms)',
                    'feature(ALTERNATIVE%20SEQUENCE)': 'Alternative sequence (isoforms)',
                    'genes': 'Gene names', 'organism': 'Organism', 'sequence': 'Sequence',
                    'protein%20names': 'Protein names'}
    params = {
        'from': 'ACC+ID',
        'to': 'ACC',
        'format': 'tab',
    }

    @retry(stop_max_attempt_number=3, wait_fixed=1000)
    def go_to_uniprot(url, params, code='utf-8'):
        sleep(uniform(0.99, 5))
        data = urllib.parse.urlencode(params)
        data = data.encode('utf-8')
        req = urllib.request.Request(url, data)
        with urllib.request.urlopen(req) as f:
            response = f.read()
        return response.decode(code)

    def get_info_from_uniprot(self, usecols, outputPath, from_list=None, from_list_file_path=None, sep='\t', chunksize=100, header=None):

        def iter_io(iter_object, params, url, outputPath):
            params['query'] = ','.join(iter_object)  # list_str
            result = UniProt_unit.go_to_uniprot(url, params)
            with open(outputPath, 'a+') as outputFile:
                outputFile.write(result)

        def tidy_result(path, colName='Entry', sep='\t'):
            df = pd.read_csv(path, sep=sep, na_values=[colName])
            df.dropna(subset=[colName], inplace=True)
            df.to_csv(path, sep=sep, index=False)

        if usecols != 'all':
            if not set(self.COLUMNS) >= set(usecols):
                print('get_info_from_uniprot(): please specified usecols with elements in %s' % self.COLUMNS)
                return False
            else:
                self.params['columns'] = ','.join(usecols)
        else:
            usecols = self.COLUMNS
            self.params['columns'] = ','.join(self.COLUMNS)

        if from_list_file_path:
            df = pd.read_csv(from_list_file_path, header=header, chunksize=chunksize, sep=sep)
            for chunk in df:
                iter_io(chunk[0], self.params, self.URL, outputPath)
        else:
            if os.path.exists(outputPath):
                new_colNames = [self.CONFIG['COLUMN_DICT'].get(i, i) for i in usecols] + ['yourlist', 'isomap']
                finish = pd.read_csv(outputPath, sep='\t', usecols=['yourlist'], names=new_colNames, skiprows=1, header=None)['yourlist']
                finish_li = []
                for i in finish:
                    if i.count(',') > 0:
                        finish_li.extend(i.split(','))
                    else:
                        finish_li.append(i)
            else:
                finish_li = []

            if finish_li:
                new_li = list(set(from_list) - set(finish_li))
            else:
                new_li = from_list

            for i in range(0, len(new_li), chunksize):
                iter_io(new_li[i:i+chunksize], self.params, self.URL, outputPath)

        tidy_result(outputPath)
        return True

    def __init__(self, dfrm, id_col, id_type, usecols, reportPath, muta_col=None, muta_type=None, gene_col=None):
        """
        Prerequisite:
        * Assume that the parameters are all legal
        * Never Change the input dataFrame
        * Data has been filtered
        """
        self.dfrm = dfrm
        self.index = dfrm.index
        self.id_col = id_col
        self.id_type = id_type
        self.muta_col = muta_col
        self.muta_type = muta_type
        self.usecols = usecols
        self.gene_col = gene_col
        if muta_col is not None:
            self.muta_li = dfrm.groupby(by=[id_col]).apply(lambda x: [i for i in x[muta_col]])
        self.report = open(reportPath, 'w+')

    def split_df(dfrm, colName, sep):
        """Split DataFrame"""
        df = dfrm.copy()
        return df.drop([colName], axis=1).join(df[colName].str.split(sep, expand=True).stack().reset_index(level=1, drop=True).rename(colName))

    def getCanonicalInfo(self, dfrm):
        """
        Will Change the dfrm

        * Add new column (canonical_isoform)
        """
        canonical_pattern = re.compile(r'IsoId=([0-9A-Z-]+); Sequence=Displayed')
        dfrm['canonical_isoform'] = dfrm.apply(lambda x: ','.join(canonical_pattern.findall(x['Alternative products (isoforms)'])) if not isinstance(x['Alternative products (isoforms)'], float) else np.nan, axis=1)
        special_case = dfrm[dfrm['canonical_isoform'] == ''].index
        if len(special_case) > 0:
            canonical_pattern = re.compile(r'IsoId=([0-9A-Z-,\s]+); Sequence=Displayed')
            special_se = dfrm.loc[special_case].apply(lambda x: ','.join(canonical_pattern.findall(x['Alternative products (isoforms)'])), axis=1)
            dfrm.loc[special_case, 'canonical_isoform'] = special_se
            self.report.write("# Special Cases of Canonical Info\n")
            self.report.write(str(special_se))
        return special_se

    def get_raw_ID_Mapping(self, outputPath):
        """
        Get Raw ID MApping Result
        """
        self.params['from'] = self.id_type
        status = self.get_info_from_uniprot(
            self.usecols,
            outputPath,
            self.dfrm[self.id_col].drop_duplicates())

        if not status:
            return False

        new_colNames = [self.CONFIG['COLUMN_DICT'].get(i, i) for i in self.usecols] + ['yourlist', 'isomap']
        dfrm = pd.read_csv(outputPath, sep='\t', names=new_colNames, skiprows=1, header=None)
        # WRITE REPORT
        self.report.write("# RAW ID MAPPING FILE RESULT\n# %s\n" % (outputPath))
        self.report.write(str(dfrm.isnull().sum()))
        self.raw_id_mapping = dfrm
        return True

    def set_raw_id_mapping(self, raw_id_mapping):
        """
        Set the ``pandas.DataFrame`` object that will be deal with in ``handle_ID_Mapping()``
        """
        # assert isinstance(raw_id_mapping, pd.DataFrame), "WARNING, raw_id_mapping should a pandas.DataFrame"
        self.raw_id_mapping = raw_id_mapping

    def handle_ID_Mapping(self):
        """
        Deal with different situations

        1. Add New Columns: ``unp_map_tage`` -> Description for the reliability of mapping, ``canonical_isoform``
        2. Split DataFrame
        3. Classification

        **About unp_map_tage**

        * ``Untrusted & No Isoform``
            * 是指UniProt存在Isoform但是Mapping结果没有明确给出是Map上哪条Isoform,转录本序列与蛋白序列不一致
            * It means that there are isoforms in UniProt, but the mapping result does not clearly indicate which isoform is correspond with the transcript(e.g), and the transcript sequence is inconsistent with the protein sequence.
        * ``Trusted & No Isoform``
            * 是指UniProt不存在Isoform,Mapping结果没问题
        * ``Trusted & Isoform``
            * 是指UniProt存在Isoform,Mapping结果没问题


        .. csv-table:: Split DataFrame Example 1
            :header: "Entry", "Gene names", "Status", "Alternative products (isoforms)", "Organism", "Protein names", "yourlist", "isomap"

            "O94827", "PLEKHG5 KIAA0720", "reviewed", "ALTERNATIVE PRODUCTS:  Event=Alternative splic...", "Homo sapiens (Human)", "Pleckstrin homology domain-containing family G...", "NP_065682,NP_001252521", "NP_001252521 -> O94827-6,NP_065682 -> O94827-5"

        .. csv-table:: Splited DataFrame Example 1
            :header: "Entry", "Gene names", "Status", "Alternative products (isoforms)", "Organism", "Protein names", "yourlist", "UniProt"

            "O94827", "PLEKHG5 KIAA0720", "reviewed", "ALTERNATIVE PRODUCTS:  Event=Alternative splic...", "Homo sapiens (Human)", "Pleckstrin homology domain-containing family G...", "NP_065682", "O94827-5"
            "O94827", "PLEKHG5 KIAA0720", "reviewed", "ALTERNATIVE PRODUCTS:  Event=Alternative splic...", "Homo sapiens (Human)", "Pleckstrin homology domain-containing family G...", "NP_001252521", "O94827-6"

        ...

        For more details, please go to https://github.com/NatureGeorge/SIFTS_Plus_Muta_Maps/blob/master/md/UniProt_ID_Mapping.md
        """

        df = self.raw_id_mapping
        # Add New Column: canonical_isoform
        canonicalInfo_special_se = self.getCanonicalInfo(df)
        # Add New Column: unp_map_tage
        df['unp_map_tage'] = np.nan
        # Classification
        df_with_isomap = df[~df['isomap'].isnull()]  # Class B
        df_with_no_isomap = df[df['isomap'].isnull()]  # Class A
        # ----------------------------------------------------------------------
        # In Class A
        # ----------------------------------------------------------------------
        df_wni_split = UniProt_unit.split_df(df_with_no_isomap, 'yourlist', ',')
        df_wni_split.drop(columns=['isomap'], inplace=True)
        df_wni_split['UniProt'] = df_wni_split['Entry']  # [yourlist <-> UniProt]
        df_wni_split['unp_map_tage'] = 'Trusted & No Isoform'
        # Find out special cases 1
        df_wni_split_warn = df_wni_split[~df_wni_split['Alternative products (isoforms)'].isnull()].index
        df_wni_split.loc[df_wni_split_warn, 'unp_map_tage'] = 'Untrusted & No Isoform'
        # 'Entry', 'Gene names', 'Status', 'Alternative products (isoforms)', 'Organism', 'yourlist', 'UniProt'
        # ----------------------------------------------------------------------
        # In Class B
        # ----------------------------------------------------------------------
        wi_yourlist_count = df_with_isomap.apply(lambda x: x['yourlist'].count(','), axis=1)
        wi_isomap_count = df_with_isomap.apply(lambda x: x['isomap'].count(','), axis=1)
        # In subClass 1
        df_wi_eq = df_with_isomap.loc[wi_yourlist_count[wi_yourlist_count == wi_isomap_count].index]
        df_wi_eq_split = UniProt_unit.split_df(df_wi_eq.drop(columns=['yourlist']), 'isomap', ',')
        df_wi_eq_split['yourlist'], df_wi_eq_split['UniProt'] = df_wi_eq_split['isomap'].str.split(' -> ', 1).str
        df_wi_eq_split.drop(columns=['isomap'], inplace=True)  # [yourlist <-> UniProt]
        df_wi_eq_split['unp_map_tage'] = 'Trusted & Isoform'
        # # 'Entry', 'Gene names', 'Status', 'Alternative products (isoforms)', 'Organism', 'yourlist', 'UniProt'

        # In subClass 2
        df_wi_ne = df_with_isomap.loc[wi_yourlist_count[wi_yourlist_count != wi_isomap_count].index]
        df_wi_ne_split = UniProt_unit.split_df(df_wi_ne, 'isomap', ',')
        df_wi_ne_split.rename(columns={'yourlist': 'checkinglist'}, inplace=True)
        df_wi_ne_split['yourlist'], df_wi_ne_split['UniProt'] = df_wi_ne_split['isomap'].str.split(' -> ', 1).str
        df_wi_ne_split.drop(columns=['isomap'], inplace=True)
        df_wi_ne_split['unp_map_tage'] = 'Trusted & Isoform & Contain Warnings'
        # 'Entry', 'Gene names', 'Status', 'Alternative products (isoforms)', 'Organism', 'yourlist', 'UniProt', 'checkinglist'
        # Find out special cases 2
        usecols = pd.Index(set(df_wi_ne_split.columns) - {'yourlist', 'UniProt'})
        df_wi_ne_warn = UniProt_unit.split_df(df_wi_ne_split[usecols].drop_duplicates(), 'checkinglist', ',')
        df_wi_ne_warn = df_wi_ne_warn[~df_wi_ne_warn['checkinglist'].isin(df_wi_ne_split['yourlist'])].rename(columns={'checkinglist': 'yourlist'})
        df_wi_ne_warn['UniProt'] = df_wi_ne_warn['Entry']
        df_wi_ne_warn['unp_map_tage'] = 'Untrusted & No Isoform'

        # Update UniProt
        final_df = pd.concat((df_wni_split, df_wi_eq_split, df_wi_ne_split.drop(columns=['checkinglist']), df_wi_ne_warn), sort=False).reset_index(drop=True)
        final_df['UniProt'] = final_df.apply(lambda x: x['Entry'] if x['UniProt'] == x['canonical_isoform'] else x['UniProt'], axis=1)
        if len(canonicalInfo_special_se) > 0:
            canonicalInfo_special_case = final_df[final_df['canonical_isoform'].isin(canonicalInfo_special_se)].index
            final_df.loc[canonicalInfo_special_case, 'UniProt'] = final_df.loc[canonicalInfo_special_case].apply(lambda x: x['Entry'] if x['UniProt'] in x['canonical_isoform'] else x['UniProt'], axis=1)
        return final_df

    def getGeneStatus(self, handled_df):
        """
        Will Change the dfrm, add Gene Status

        * Add new column (GENE)
        * Add new column (GENE_status)

        **About GENE_status**

        * ``False`` : First element of Gene names is not correspond with refSeq's GENE (e.g)
        * others(corresponding GENE)

        """
        if not self.gene_col:
            return None
        gene_map = self.dfrm[[self.id_col, self.gene_col]].drop_duplicates()
        gene_map.index = gene_map[self.id_col]
        gene_map.drop(columns=[self.id_col], inplace=True)
        handled_df['GENE'] = handled_df.apply(lambda z: gene_map['GENE'][z['yourlist']], axis=1)
        handled_df['GENE_status'] = handled_df.apply(lambda x: x['GENE'] == x['Gene names'].split(' ')[0] if not isinstance(x['Gene names'], float) else False, axis=1)

    def label_mapping_status(self, dfrm, constraint_dict):
        """
        Will Change the dfrm, label Mapping Status

        * Add new column (Mapping_status)

        **About Mapping_status**

        * ``Yes``: 可信的结果，进行后续的PDB Mapping; 通过 ``constraint_dict`` 的限制
        * ``Error``: 一个id对应多个UniProt; 通过 ``constraint_dict`` 的限制
        * ``No``: 不可信的结果; 未通过 ``constraint_dict`` 的限制

        Example of ``constraint_dict``:

        .. code-block:: python
            :linenos:

            constraint_dict = {
                "GENE_status": (False, "ne"),  # 'ne' for !=
                "Status": ("reviewed", "eq"),  # 'eq' for ==
                "unp_map_tage": ("Untrusted & No Isoform", "ne")
            }

        """
        dfrm['Mapping_status'] = 'No'
        # pass_index = dfrm[(dfrm['GENE_status'] != False) & (dfrm['Status'] == 'reviewed') & (dfrm['unp_map_tage'] != 'Untrusted & No Isoform')].index
        pass_df = Unit.ConstraintDict.addConstraintToDf(dfrm, constraint_dict)
        pass_index = pass_df.index
        dfrm.loc[pass_index, 'Mapping_status'] = 'Yes'

        multipleCounter, err_li = Counter(dfrm.loc[pass_index, 'yourlist']), []
        for i, j in multipleCounter.items():
            if j > 1:
                err_li.append(i)

        err_index = pass_df[pass_df['yourlist'].isin(err_li)].index
        dfrm.loc[err_index, 'Mapping_status'] = 'Error'

        # Write Report
        all_id = set(self.dfrm[self.id_col])
        self.report.write("\n# All id: %s\n" % len(all_id))
        unmapped_id = all_id - set(dfrm['yourlist'])
        self.report.write("# Unmapped id: %s\n" % len(unmapped_id))
        for i in sorted(unmapped_id):
            self.report.write("%s\n" % i)
        untrusted_id = set(dfrm[dfrm['Mapping_status'] != 'Yes']['yourlist']) - set(dfrm[dfrm['Mapping_status'] == 'Yes']['yourlist'])
        self.report.write("# Untrusted id: %s\n" % len(untrusted_id))
        for i in sorted(untrusted_id):
            self.report.write("%s\n" % i)
        self.report.write("# Error id: %s\n" % len(dfrm.loc[err_index, 'yourlist'].drop_duplicates()))
        self.report.write(str(dfrm.loc[err_index]))

    def script(self, rawOutputPath, handledOutputPath):
        """
        .. code-block:: python
            :linenos:

            from UniProt_unit import UniProt_unit
            id_col = 'RefSeq_protein'
            id_type = 'P_REFSEQ_AC'
            muta_col = 'mutation_unp'
            gene_col = 'GENE'
            usecols = ['id', 'genes', 'reviewed', 'comment(ALTERNATIVE%20PRODUCTS)', 'organism', 'protein%20names']  # Necessary Columns
            reportPath = '/data/zzf/UniProt_files/id_mapping_files/HGMD_RefSeq_protein_mapping_Report_1021.txt'  # Report of Data Processing
            rawOutputPath = '/data/zzf/UniProt_files/id_mapping_files/HGMD_RefSeq_protein_mapping_1021.tsv' # OutPut File of ID Mapping (RAW)
            handledOutputPath = '/data/zzf/groupWorks/HGMD_RefSeq_protein_mapping_modified_1024.tsv'  # OutPut File of ID Mapping (Final Result)
            # Add constraint/filter to data
            constraint_dict = {
                "GENE_status": (False, "ne"),  # 'ne' for !=
                "Status": ("reviewed", "eq"),  # 'eq' for ==
                "unp_map_tage": ("Untrusted & No Isoform", "ne")
            }

            # Initial
            unp_demo = UniProt_unit(group_df, id_col, id_type, usecols, reportPath, muta_col=muta_col, gene_col=gene_col)
            # Return True if get RAW Result Successfully
            unp_demo.get_raw_ID_Mapping(rawOutputPath)
            # Deal with different situations
            handled_df = unp_demo.handle_ID_Mapping()
            # Add Gene Status
            unp_demo.getGeneStatus(handled_df)
            # Label Mapping Status
            unp_demo.label_mapping_status(handled_df, constraint_dict)
            # close the file-handle of report
            unp_demo.report.close()
            # Output the final result
            handled_df.to_csv(handledOutputPath, sep='\t', index=False)

        """
        step1 = self.get_raw_ID_Mapping(rawOutputPath)
        if not step1:
            return None
        handled_df = self.handle_ID_Mapping()
        self.getGeneStatus(handled_df)
        constraint_dict = {
            "GENE_status": (False, "ne"),
            "Status": ("reviewed", "eq"),
            "unp_map_tage": ("Untrusted & No Isoform", "ne")
        }
        self.report.write("Constraint Dict\n"+json.dumps(constraint_dict)+"\n")
        self.label_mapping_status(handled_df, constraint_dict)
        handled_df.to_csv(handledOutputPath, sep='\t', index=False)
        self.report.close()

```
