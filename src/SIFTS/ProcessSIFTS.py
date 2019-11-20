# @Date:   2019-11-20T16:59:18+08:00
# @Email:  1730416009@stu.suda.edu.cn
# @Filename: ProcessSIFTS.py
# @Last modified time: 2019-11-20T19:43:45+08:00
import pandas as pd
import numpy as np
import json
import time
import wget
import os
from urllib import request
from Utils.Logger import RunningLogger
from Utils.FileIO import decompression, file_i, file_o
from Utils.Tools import Gadget


SIFTS_URL = 'http://www.ebi.ac.uk/pdbe/api/mappings/all_isoforms/%s'
SIFTS_FILE_URL = {
    "uniprot_pdb": "ftp://ftp.ebi.ac.uk/pub/databases/msd/sifts/flatfiles/csv/uniprot_pdb.csv.gz"
}
PDB_ID = "pdb_id"
RAW_SIFTS_COLUMNS = [
    'pdb_id', 'chain_id', 'UniProt', 'identity', 'identifier', 'pdb_start',
    'pdb_end', 'unp_start', 'unp_end', 'is_canonical', 'start', 'end',
    'entity_id', 'struct_asym_id'
]
SEG_SIFTS_MAP = {
      "SP_PRIMARY": "UniProt",  # (Canonical)
      "RES_BEG": "pdb_start",
      "RES_END": "pdb_end",
      "PDB_BEG": "residue index(start) in pdb",
      "PDB_END": "residue index(end) in pdb",
      "SP_BEG": "unp_start",
      "SP_END": "unp_end",
      "PDB": 'pdb_id',
      "CHAIN": 'chain_id',
}


class RetrieveSIFTS:
    def __init__(self, **kwargs):
        # Initialize the logger
        loggingPath = kwargs.get("loggingPath", None)
        self.Logger = RunningLogger("RetrieveSIFTS", loggingPath)
        # Need to be checked
        self.rawSIFTSpath = kwargs.get("rawSIFTSpath", None)
        self.downloadFolder = kwargs.get("downloadFolder", None)

    def sort_SIFTS_info(self, pdbId, info):
        first = True
        order = RAW_SIFTS_COLUMNS
        for uniprot in info.keys():
            gene = info[uniprot]['identifier']
            pdbChain = info[uniprot]['mappings']
            for chain in pdbChain:
                chain['start'] = json.dumps(chain['start'])
                chain['end'] = json.dumps(chain['end'])
                chain[PDB_ID] = pdbId
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
            return None

    def retrieve_raw_SIFTS(self, pdbs, outputPath=None, sleepTime=1.5):
        def addHeader(rows):
            if rows == 0:
                return RAW_SIFTS_COLUMNS
            else:
                return False

        if outputPath is None:
            outputPath = self.rawSIFTSpath

        try:
            rows = len(pd.read_csv(outputPath, sep='\t', usecols=['UniProt']))
        except FileNotFoundError:
            rows = 0

        fail_list = []
        allPDB, current = len(pdbs), 0
        for pdbId in pdbs:
            time.sleep(sleepTime)
            url = SIFTS_URL % (pdbId)
            try:
                req = request.Request(url)
                page = request.urlopen(req).read()
                page = page.decode('utf-8')
                info = json.loads(page)[pdbId]['UniProt']
                if info:
                    df = self.sort_SIFTS_info(pdbId, info)
                    if df is not None:
                        df.to_csv(outputPath, index=False, sep='\t', mode='a+', header=addHeader(rows))
            except Exception as e:
                self.Logger.logger.error('%s, %s' % (pdbId, e))
                fail_list.append(pdbId)
            current += 1
            self.Logger.logger.info('retrieve_raw_SIFTS: End a circle. %s, cur:%s, sum:%s' % (pdbId, current, allPDB))
        return rows, fail_list

    def get_info_from_uniprot_pdb_file(self, filePath=None, related_unp=None, related_pdb=None):
        '''
            Reference: (http://www.ebi.ac.uk/pdbe/docs/sifts/quick.html)

                A summary of the UniProt to PDB mappings showing the UniProt accession
                followed by a semicolon-separated list of PDB four letter codes.
        '''
        if filePath is None:
            url = SIFTS_FILE_URL["uniprot_pdb"]
            filePath = os.path.join(self.downloadFolder, 'uniprot_pdb_%s.csv.gz' % time.strftime("%Y_%m_%d", time.localtime()))

            try:
                wget.download(url, out=filePath)
                decompression(filePath, remove=True, logger=self.Logger.logger)
            except Exception:
                self.Logger.logger.error("Download failed: %s" % filePath)

            self.Logger.logger.info("Download File: %s" % filePath[:-3])

        dfrm = pd.read_csv(filePath, sep=',', header=1)
        pdb_list = []
        if related_unp is not None:
            dfrm = dfrm[dfrm['SP_PRIMARY'].isin(related_unp)]
        for i in dfrm.index:
            pdb_list.extend(dfrm.loc[i, 'PDB'].split(';'))
        if related_pdb is not None:
            return {'pdb_set': set(pdb_list) & set(related_pdb), 'unp_set': set(dfrm['SP_PRIMARY'])}
        else:
            return {'pdb_set': set(pdb_list), 'unp_set': set(dfrm['SP_PRIMARY'])}


def handle_SIFTS(filePath, skiprows=0, outputPath=None):
    def addSiftsRange(in_df):
        group_info_col = RAW_SIFTS_COLUMNS[:3]
        range_info_col = RAW_SIFTS_COLUMNS[5:9]
        rangeSetER = Gadget.RangeSetER(group_info_col)
        in_df['rangeInfo'] = in_df.apply(lambda x: rangeSetER.check(
            tuple(x[i] for i in group_info_col),
            tuple(x[i] for i in range_info_col)),
            axis=1)
        return in_df.drop(columns=range_info_col).drop_duplicates(subset=group_info_col, keep='last')

    sifts_df = pd.read_csv(filePath, sep='\t', na_values=[PDB_ID, '', None], keep_default_na=False, skiprows=skiprows, names=RAW_SIFTS_COLUMNS)
    sifts_df.dropna(subset=[PDB_ID], inplace=True)
    sifts_df[PDB_ID] = sifts_df[PDB_ID].str.upper()

    new_sifts_df = addSiftsRange(sifts_df)
    new_sifts_df['sifts_pdb_range'] = new_sifts_df.apply(lambda x: x['rangeInfo'].split('|')[0], axis=1)
    new_sifts_df['sifts_unp_range'] = new_sifts_df.apply(lambda x: x['rangeInfo'].split('|')[1], axis=1)
    new_sifts_df.drop(columns=['rangeInfo'], inplace=True)

    file_o(outputPath, new_sifts_df)
    return new_sifts_df


def deal_with_insertionDeletion_SIFTS(self, sifts_df=None, sifts_filePath=None, outputPath=None):
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

    dfrm = file_i(sifts_filePath, sifts_df, ('sifts_filePath', 'sifts_df'))
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
    file_o(outputPath, dfrm)
    return dfrm
