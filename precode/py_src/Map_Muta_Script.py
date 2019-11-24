# @Date:   2019-10-04T19:20:23+08:00
# @Email:  1730416009@stu.suda.edu.cn
# @Filename: Map_Muta_Script.py
# @Last modified time: 2019-10-04T19:29:44+08:00
import pandas as pd
import numpy as np
import sys, os, re
sys.path.append('/home/zzf/Work/SIFTS_Plus_Muta_Maps/src/py_src')
from MMCIFplus import *
from SIFTS_unit import SIFTS_unit

mmcif_demo = MMCIF2Dfrm()
sifts_demo = SIFTS_unit()

group_df = pd.read_csv('/data/cheny/MapTool_zzf/SNV_2009_uniprot_unique.txt', sep='\t')
unp_li = group_df['UniprotID'].drop_duplicates()
muta_li = group_df[['UniprotID','Mutation_unp']].groupby(by='UniprotID').apply(lambda x:[i for i in x['Mutation_unp']])

demo_file_route = '/home/zzf/Work/SIFTS_Plus_Muta_Maps/data/demo_files/'
add_unpLen_segInfo_file_path = demo_file_route + 'pdb_uniprot_SIFTS_addunpLenSegInfo_demo.tsv'
sifts_df = pd.read_csv(add_unpLen_segInfo_file_path, sep='\t', converters={'chain_id':str, 'struct_asym_id':str, 'entity_id':str})
use_sifts_df = sifts_df[sifts_df['UniProt'].isin(unp_li)].drop_duplicates().copy()
pdb_li = use_sifts_df['pdb_id'].drop_duplicates()

raw_mmcif_path = '/data/zzf/MMCIF_files/new_rawMMCIF2Dfrm.tsv'
handled_mmcif_path = '/data/zzf/MMCIF_files/new_handledMMCIF2Dfrm.tsv'
mmcif_demo.check_mmcif_file(pdb_li)

if os.path.exists(handled_mmcif_path):
    finished_li = set(pd.read_csv(handled_mmcif_path, sep='\t', usecols=['pdb_id'])['pdb_id'])
    mmcif_file_li = []
    chunksize = 100
    for path in mmcif_demo.pdb_path_li:
        if path[-8:-4] not in finished_li:
            mmcif_file_li.append(path)
        else:
            print(path[-8:-4],'FINISHED')
else:
    mmcif_file_li = mmcif_demo.pdb_path_li

chunksize = 100
for i in range(0, len(mmcif_file_li), chunksize):
    chunk_li = mmcif_file_li[i:i+chunksize]
    chunk_df = mmcif_demo.mmcif_dict2dfrm(chunk_li, outputPath=raw_mmcif_path)
    mmcif_demo.handle_mmcif_dfrm(chunk_df, outputPath=handled_mmcif_path)
