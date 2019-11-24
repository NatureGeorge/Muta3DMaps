# @Date:   2019-08-16T23:31:03+08:00
# @Email:  1730416009@stu.suda.edu.cn
# @Filename: MMCIF_unit.py
# @Last modified time: 2019-08-21T16:02:36+08:00
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
        'PDB_ID': 'pdb_id',
        'MMCIF_OLD_FOLDER': ['/data1/suntt/process0606/cgc_mmcif_file/', '/data1/suntt/CanDriver/Data/PDB_cgc/cgc_mmcif_file/', '/data1/suntt/CanDriver/Data/PDB_NEW/mmcif_file/'],
        'MMCIF_FOLDER': '/home/zzf/Work/SIFTS_Plus_Muta_Maps/data/mmcif_file/',
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


if __name__ == '__main__':
    demo = MMCIF_unit()
    demo.set_lists(['Test is success.'], [])
    print(demo.pdb_list)
