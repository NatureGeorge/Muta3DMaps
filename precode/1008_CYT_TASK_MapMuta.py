# @Date:   2019-10-08T11:04:52+08:00
# @Email:  1730416009@stu.suda.edu.cn
# @Filename: 1008_CYT_TASK_MapMuta.py
# @Last modified time: 2019-10-08T11:13:54+08:00
import pandas as pd
import numpy as np
import sys, os, re
sys.path.append('/home/zzf/Work/SIFTS_Plus_Muta_Maps/src/py_src')
from MMCIFplus import *
from SIFTS_unit import SIFTS_unit


def setCompo(a, b):
    li = sorted([a, b])
    return '%s_%s' % (li[0], li[1])

def manual(groupFilePath, usecols=None):
    mmcif_demo = MMCIF2Dfrm()
    sifts_demo = SIFTS_unit()
    group_df = pd.read_csv(groupFilePath, sep='\t', usecols=usecols)
    group_df['compo'] = group_df.apply(lambda x: setCompo(x['Target_UniprotID'],x['Interactor_UniprotID']), axis=1)
    unp_li = group_df['Target_UniprotID'].drop_duplicates()
    unp_set = set(group_df['Target_UniprotID']) | set(group_df['Interactor_UniprotID'])

def manual_gotoUniProt(mapping_file_path):
    from UniProt_unit import UniProt_unit
    unp_demo = UniProt_unit()
    usecols = ['id', 'comment(ALTERNATIVE%20PRODUCTS)']
    new_colNames = ['Entry', 'Alternative products (isoforms)', 'yourlist','isomap']

    if os.path.exists(mapping_file_path):
        finish = pd.read_csv(mapping_file_path,
                             sep='\t', usecols=['yourlist'],
                             names=new_colNames,
                             skiprows=1, header=None)['yourlist']
        finish_li = []
        for i in finish:
            if i.count(',') > 0:
                finish_li.extend(i.split(','))
            else:
                finish_li.append(i)
    else:
        finish_li = []
    unp_demo.get_info_from_uniprot(usecols, mapping_file_path, unp_li)

def updateUniProtID(mapping_file_path, new_colNames):
    unp_df = pd.read_csv(mapping_file_path, sep='\t', names=new_colNames, skiprows=1, header=None)
    find_canonical_pattern = re.compile(r'IsoId=([0-9A-Z-]+); Sequence=Displayed')
    unp_df['canonical_isoform'] = unp_df.apply(lambda x: ','.join(find_canonical_pattern.findall(x['Alternative products (isoforms)'])) if not isinstance(x['Alternative products (isoforms)'], float) else np.nan, axis=1)
    canonical_se = set(unp_df['canonical_isoform'].dropna())
    group_df['Target_UniprotID'] = group_df.apply(lambda x: x['Target_UniprotID'].split('-')[0] if x['Target_UniprotID'] in canonical_se else x['Target_UniprotID'], axis=1)
    unp_li = group_df['Target_UniprotID'].drop_duplicates()



if __name__ == "__main__":
    mapping_file_path = '/home/zzf/Work/SIFTS_Plus_Muta_Maps/data/PremPS/data/findUniProtCanonical_1004.tsv'
    manual('/data/cheny/yuhaiyuan/SNV_2009_addUniprot_del.txt', ['Target_UniprotID','Interactor_UniprotID','Target_Mutation_unp'])
    # manual_gotoUniProt(mapping_file_path)
