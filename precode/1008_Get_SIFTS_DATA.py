# @Date:   2019-10-08T15:14:21+08:00
# @Email:  1730416009@stu.suda.edu.cn
# @Filename: 1008_Get_SIFTS_DATA.py
# @Last modified time: 2019-10-08T15:34:48+08:00
import sys
sys.path.append('/home/zzf/Work/SIFTS_Plus_Muta_Maps/src/py_src')
from SIFTS_unit import SIFTS_unit

raw_sifts_file_path = 'pdb_uniprot_SIFTS_raw_demo.csv'
add_rangeInfo_sifts_file_path = 'pdb_uniprot_SIFTS_addRangeInfo_demo.tsv'
add_InDe_sifts_file_path = 'pdb_uniprot_SIFTS_delwithInDe_demo.tsv'
SIFTS_unit.CONFIG['DOWNLOAD_FOLDER'] = '/data/zzf/SIFTS_files/'
sifts_demo = SIFTS_unit()
info_dict = sifts_demo.get_info_from_uniprot_pdb_file()
sifts_demo.pdb_list = sorted(info_dict['pdb_set'])
sifts_demo.get_raw_SIFTS(outputPath='/data/zzf/SIFTS_files/%s' % raw_sifts_file_path)
handle_df = sifts_demo.handle_SIFTS(outputPath='/data/zzf/SIFTS_files/%s' % add_rangeInfo_sifts_file_path)
sifts_demo.deal_with_insertionDeletion_SIFTS(sifts_df=handle_df, outputPath='/data/zzf/SIFTS_files/%s' % add_InDe_sifts_file_path)
