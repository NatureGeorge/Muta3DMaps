# @Date:   2019-08-19T19:29:29+08:00
# @Email:  1730416009@stu.suda.edu.cn
# @Filename: MMCIFplus_unit.py
# @Last modified time: 2019-08-21T00:16:09+08:00
from collections import defaultdict
import pandas as pd
import numpy as np
import os, re, time, requests, sys
from urllib import request, error
from retrying import retry
from multiprocessing.dummy import Pool
from bs4 import BeautifulSoup
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
sys.path.append('./')
from Unit import Unit

class MMCIF_unit(Unit):
    CONFIG = {
        'MMCIF_OLD_FOLDER': [
            '/data1/suntt/process0606/cgc_mmcif_file/',
            '/data1/suntt/CanDriver/Data/PDB_cgc/cgc_mmcif_file/',
            '/data1/suntt/CanDriver/Data/PDB_NEW/mmcif_file/'
            ],
        'MMCIF_FOLDER': '../../data/Mapping_Pipeline/mmcif_file/',
        'COMMON_COL': ['_pdbx_audit_revision_history.revision_date', '_exptl.method', '_em_3d_reconstruction.resolution', '_refine.ls_d_res_high'],
        'ENTITY_COL': ['_entity.pdbx_mutation', '_entity.id'],
        'TYPE_COL':['_entity_poly.entity_id', '_entity_poly.pdbx_strand_id', '_entity_poly.type'],
        'SEQRES_COL':['_pdbx_poly_seq_scheme.pdb_strand_id',
                 '_pdbx_poly_seq_scheme.mon_id','_pdbx_poly_seq_scheme.pdb_mon_id', '_pdbx_poly_seq_scheme.auth_mon_id',
                 '_pdbx_poly_seq_scheme.ndb_seq_num', '_pdbx_poly_seq_scheme.pdb_seq_num',
                 '_pdbx_poly_seq_scheme.auth_seq_num', '_pdbx_poly_seq_scheme.pdb_ins_code'],
        'LIGAND_COL': [
                 '_struct_conn.ptnr2_auth_asym_id','_struct_conn.ptnr2_auth_comp_id',
                 '_struct_conn.ptnr2_auth_seq_id',
                 '_struct_conn.conn_type_id',
                 '_struct_conn.ptnr1_auth_asym_id', '_struct_conn.ptnr1_auth_comp_id',
                 '_struct_conn.ptnr1_auth_seq_id'],
        'LIGAND_LIST': [
                        'ZN', 'MG', 'CA', 'FE', 'NA', 'MN', 'K', 'NI', 'CU', 'CO', 'CD', 'HG', 'PT', 'MO', 'BE', 'AL', 'BA',
                        'RU', 'SR', 'V', 'CS', 'W', 'AU', 'YB', 'LI', 'GD', 'PB', 'Y', 'TL', 'IR', 'RB', 'SM', 'AG',
                        'OS', 'PR', 'PD', 'EU', 'RH', 'RE', 'TB', 'TA', 'LU', 'HO', 'CR', 'GA', 'LA', 'SN', 'SB', 'CE',
                        'ZR', 'ER', 'TH', 'TI', 'IN', 'HF', 'SC', 'DY', 'BI', 'PA', 'PU', 'AM', 'CM', 'CF', 'GE', 'NB', 'TC',
                        'ND', 'PM', 'TM', 'PO', 'FR', 'RA', 'AC', 'NP', 'BK', 'ES', 'FM', 'MD', 'NO', 'LR', 'RF', 'DB', 'SG']
    }

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

    def get_mmcif_info_old(info_key, info_dict, path):
        mmcif_dict = MMCIF2Dict(path)
        for key in info_key:
            data = mmcif_dict.get(key,[])
            if isinstance(data, str):
                data = data.split(',')
            else:
                print(data)
            info_dict[key].append(list(filter(lambda x :x != '?' and x != '.', data)))

    def get_mmcif_info(info_key, info_key_nli, info_dict, path):
        mmcif_dict = MMCIF2Dict(path)
        for key in info_key:

            if key in info_key_nli:
                data = mmcif_dict.get(key,np.nan)
                info_dict[key].append(data)
            else:
                data = mmcif_dict.get(key,[])
                if isinstance(data, str):
                    info_dict[key].append([data])
                    ## data = data.split(',')
                    # print(key)
                # info_dict[key].append(list(filter(lambda x :x not in '?.', data)))
                else:
                    info_dict[key].append(data)

    def get_data_from_mmcif(self, path_list):
        '''
        {
            '_pdbx_audit_revision_history.revision_date': ['initial_version_time', 'newest_version_time'], # sometimes not a list
            '_entity.pdbx_mutation': ['mutation_num', 'mutation_content'], # sometimes not a list
            '_entity.id': ['entity_id_aidMuta'], # sometimes not a list
            ['_em_3d_reconstruction.resolution','_refine.ls_d_res_high']: ['resolution'], # not a lists
            '_exptl.method': ['method'], # not a list
        }

        '''
        info_dict = defaultdict(list)
        for path in path_list:
            if path[-3:] == 'cif':
                print(path)
                info_dict['pdb_id'].append(path[-8:-4])
                MMCIF_unit.get_mmcif_info(
                    MMCIF_unit.CONFIG['COMMON_COL'] + \
                    MMCIF_unit.CONFIG['ENTITY_COL'] + \
                    MMCIF_unit.CONFIG['TYPE_COL'] + \
                    MMCIF_unit.CONFIG['SEQRES_COL'] + \
                    MMCIF_unit.CONFIG['LIGAND_COL'],
                    MMCIF_unit.CONFIG['COMMON_COL'][1:],
                    info_dict,
                    path)

        # Deal with Residues in SEQRES_COL
        resides_col_li = MMCIF_unit.CONFIG['SEQRES_COL'][1:4]
        mtoTool = Unit.MultiToOne()
        for i in range(len(resides_col_li[0])):
            for resides_col in resides_col_li:
                info_dict[resides_col][i] = ''.join([mtoTool.multi_letter_convert_to_one_letter(j) for j in info_dict[resides_col][i]])

        get_index = lambda x, y, z: y[x[z]:x[z+1]] if len(x) != 1 and z+1 < len(x) else y[x[z]:]
        # Deal with SEQRES_COL
        pdbx_poly_key = MMCIF_unit.CONFIG['SEQRES_COL'][0]
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

            info_dict['_pdbx_poly_seq_scheme.mon_id'][i] = [
                get_index(strand_id_index, info_dict['_pdbx_poly_seq_scheme.mon_id'][i], j)
                for j in range(len(strand_id_index))]

            info_dict['_pdbx_poly_seq_scheme.auth_mon_id'][i] = [
                get_index(strand_id_index, info_dict['_pdbx_poly_seq_scheme.auth_mon_id'][i], j)
                for j in range(len(strand_id_index))]

            info_dict['_pdbx_poly_seq_scheme.pdb_mon_id'][i] = [
                get_index(strand_id_index, info_dict['_pdbx_poly_seq_scheme.pdb_mon_id'][i], j)
                for j in range(len(strand_id_index))]

            info_dict['_pdbx_poly_seq_scheme.ndb_seq_num'][i] = [';'.join(
                get_index(strand_id_index, info_dict['_pdbx_poly_seq_scheme.ndb_seq_num'][i], j))
                for j in range(len(strand_id_index))]

            info_dict['_pdbx_poly_seq_scheme.pdb_seq_num'][i] = [';'.join(
                get_index(strand_id_index, info_dict['_pdbx_poly_seq_scheme.pdb_seq_num'][i], j))
                for j in range(len(strand_id_index))]

            info_dict['_pdbx_poly_seq_scheme.auth_seq_num'][i] = [';'.join(
                get_index(strand_id_index, info_dict['_pdbx_poly_seq_scheme.auth_seq_num'][i], j))
                for j in range(len(strand_id_index))]

            info_dict['_pdbx_poly_seq_scheme.pdb_ins_code'][i] = [';'.join(
                get_index(strand_id_index, info_dict['_pdbx_poly_seq_scheme.pdb_ins_code'][i], j))
                for j in range(len(strand_id_index))]
        # Deal with LIGAND_COL: Sort the data
        ligand_col_list = MMCIF_unit.CONFIG['LIGAND_COL']
        for i in range(len(info_dict[ligand_col_list[0]])):
            ligand_col_tp = tuple(info_dict[col][i] for col in ligand_col_list)
            ligand_col_zip_li = list(zip(*ligand_col_tp))
            ligand_col_zip_li.sort()
            for col_index in range(len(ligand_col_list)):
                info_dict[ligand_col_list[col_index]][i] = [tp[col_index] for tp in ligand_col_zip_li]
        # Deal with LIGAND_COL: Group the data
        ligand_group_col = ligand_col_list[0]
        for i in range(len(info_dict[ligand_group_col])):
            strand_id_index = [0]
            li = info_dict[ligand_group_col][i]
            if not li:
                info_dict['%s_index'%ligand_group_col].append([])
                info_dict['%s_li'%ligand_group_col].append([])
                continue
            save_id = li[0]
            strand_id_li = [save_id]
            for j in range(len(li)):
                if li[j] != save_id:
                    save_id = li[j]
                    strand_id_index.append(j)
                    strand_id_li.append(save_id)
            info_dict['%s_li'%ligand_group_col].append(strand_id_li)
            info_dict['%s_index'%ligand_group_col].append(strand_id_index)

            for col in ligand_col_list:
                info_dict[col][i] = [
                    get_index(strand_id_index, info_dict[col][i], j)
                    for j in range(len(strand_id_index))]

        df = pd.DataFrame(info_dict)
        # Deal with the date of structure
        df['initial_version_time'] = df.apply(lambda x: x['_pdbx_audit_revision_history.revision_date'][0], axis=1)
        df['newest_version_time'] = df.apply(lambda x: x['_pdbx_audit_revision_history.revision_date'][-1], axis=1)
        # Deal with the mutations
        muta_count = lambda x: x.count(',')+1 if x!= '?' else 0
        df['mutation_num'] = df.apply(lambda x: [muta_count(i) for i in x['_entity.pdbx_mutation']], axis=1)
        # Deal with the resolution
        df['resolution'] = df.apply(lambda x: x['_refine.ls_d_res_high'] if x['_exptl.method']=='X-RAY DIFFRACTION' else x['_em_3d_reconstruction.resolution'], axis=1)
        # Change the columns
        df.rename(columns={'_exptl.method':'method'},inplace=True)
        df.drop(columns=['_pdbx_audit_revision_history.revision_date','_em_3d_reconstruction.resolution','_refine.ls_d_res_high'],inplace=True)

        return df



if __name__ == '__main__':
    route = 'C:\\Users\\Nature\\Desktop\\LiGroup\\Filter_new_20190123\\doc_in\\'
    file_list = os.listdir(route)
    file_p_list = [route+i for i in file_list]
    mmcif_demo = MMCIF_unit()
    df = mmcif_demo.get_data_from_mmcif(file_p_list)
    for i in df.index:
        print(df.loc[i,])
