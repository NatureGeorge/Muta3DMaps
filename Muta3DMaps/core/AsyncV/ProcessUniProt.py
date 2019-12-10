# @Date:   2019-11-22T15:19:51+08:00
# @Email:  1730416009@stu.suda.edu.cn
# @Filename: ProcessUniProt.py
# @Last modified time: 2019-12-08T18:47:01+08:00
import os
import time
import sys
import asyncio
import aiohttp
from aiohttp import web
import pandas as pd
import numpy as np
import argparse
from random import uniform
from collections import Counter
from collections.abc import Iterable
import re
from Logger import RunningLogger


DEMO_ID_LYST = (
    'NP_003318',
    'NP_542172',
    'NP_689414',
    'NP_001034300',
    'NP_001108220',
    'NP_000806',
    'NP_003027',
    'NP_689705',
    'NP_001157196',
    'NP_056372',
    'NP_997253',
    'NP_009193',
    'NP_001073866',
    'NP_004276',
    'NP_079382',
    'NP_073624',
    'NP_004949',
    'NP_036300',
    'NP_001120797',
    'NP_006163')

COLUMNS = [
    'id', 'length', 'reviewed', 'comment(ALTERNATIVE%20PRODUCTS)',
    'feature(ALTERNATIVE%20SEQUENCE)', 'genes', 'organism', 'protein%20names'
    ]

COLUMN_DICT = {
    'id': 'Entry', 'length': 'Length', 'reviewed': 'Status',
    'comment(ALTERNATIVE%20PRODUCTS)': 'Alternative products (isoforms)',
    'feature(ALTERNATIVE%20SEQUENCE)': 'Alternative sequence (isoforms)',
    'genes': 'Gene names', 'organism': 'Organism', 'sequence': 'Sequence',
    'protein%20names': 'Protein names'}

BASE_URL = 'https://www.uniprot.org/uploadlists/'

PARAMS = {
    'from': 'P_REFSEQ_AC',
    'to': 'ACC',
    'format': 'tab',
    }

LOGGER_PATH = 'AsyncProcessUnp.log'
UniProt_ID_Mapping_RAW_PATH = 'UniProt_ID_Mapping_raw.tsv'
UniProt_ID_Mapping_MODIFIED_PATH = 'UniProt_ID_Mapping_modified.tsv'


class MapUniProtID:

    def __init__(self, dfrm, id_col, id_type, usecols, loggingPath, site_col=None, gene_col=None):
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
        self.site_col = site_col
        self.usecols = usecols
        self.gene_col = gene_col
        if site_col is not None:
            self.site_li = dfrm.groupby(by=[id_col]).apply(
                lambda x: [i for i in x[site_col]])
        self.Logger = RunningLogger("MapUniProtID", loggingPath)

    @staticmethod
    def createLystStr(lyst, chunksize=100):
        strLyst = [','.join(lyst[i:i+chunksize])
                   for i in range(0, len(lyst), chunksize)]
        return strLyst

    @staticmethod
    def save_data(path, data):
        with open(path, 'ab') as fp:
            fp.write(data)

    @staticmethod
    async def get_data(session, chunkLystStr):
        PARAMS['query'] = chunkLystStr
        async with session.get(BASE_URL, params=PARAMS) as resp:
            if resp.status == 200:
                return await resp.read()
            elif resp.status == 404:
                raise web.HTTPNotFound()
            else:
                mes = "code={resp.status}, message={resp.reason}, headers={resp.headers}"
                raise Exception(mes.format(resp=resp))

    async def download_one_chunk(self, session, chunkLystStr, path, semaphore, index):
        async with semaphore:
            status = index % 2
            if status:
                # sec = uniform(5, 10)
                sec = (index/2)**1.2
                self.Logger.logger.info("Start to sleep {:.2f}s".format(sec))
                await asyncio.sleep(sec)
            
            
            self.Logger.logger.info("Start to get data")
            rawData = await self.get_data(session, chunkLystStr)
            
            self.Logger.logger.info(chunkLystStr[:50] + '...')
            self.save_data(path, rawData)
        return chunkLystStr

    async def download_many(self, chunkLystStr_list, path, concur_req):
        semaphore = asyncio.Semaphore(concur_req)
        async with aiohttp.ClientSession() as session:
            res = await asyncio.gather(
                *[asyncio.create_task(self.download_one_chunk(session, chunkLystStr, path, semaphore, index))
                    for index, chunkLystStr in enumerate(chunkLystStr_list)])

        return len(res)

    def get_info_from_uniprot(self, outputPath, from_list=None, sep='\t', chunksize=100, concur_req=10, header=None):

        def tidy_result(path, colName='Entry', sep='\t'):
            df = pd.read_csv(path, sep=sep, na_values=[colName])
            df.dropna(subset=[colName], inplace=True)
            df.to_csv(path, sep=sep, index=False)

        if self.usecols != 'default':
            if not set(COLUMNS) >= set(self.usecols):
                self.Logger.logger.error(
                    'get_info_from_uniprot(): please specified usecols with elements in %s' % COLUMNS)
                return False
        else:
            self.usecols = COLUMNS
        
        PARAMS['columns'] = ','.join(self.usecols)

        if os.path.exists(outputPath):
            new_colNames = [COLUMN_DICT.get(i, i)
                            for i in self.usecols] + ['yourlist', 'isomap']
            try:
                finish = pd.read_csv(outputPath, sep='\t', usecols=[
                                     'yourlist'], names=new_colNames, skiprows=1, header=None)['yourlist']
            except Exception:
                self.Logger.logger.warning(
                    "Something wrong with finished raw file, probably without 'isomap' column.")
                finish = pd.read_csv(
                    outputPath, sep='\t', names=new_colNames[:-1], skiprows=1, header=None)
                finish['isomap'] = np.nan
                finish.to_csv(outputPath, sep="\t", index=False)
                finish = finish['yourlist']

            finish_li = []
            for i in finish:
                if i.count(',') > 0:
                    finish_li.extend(i.split(','))
                else:
                    finish_li.append(i)

        else:
            finish_li = []

        self.Logger.logger.info("Have finished {} ids".format(len(finish_li)))

        if finish_li:
            new_li = list(set(from_list) - set(finish_li))
        else:
            new_li = from_list

        t0 = time.perf_counter()
        count = asyncio.run(
            self.download_many(
                self.createLystStr(new_li, chunksize), outputPath, concur_req))
        elapsed = time.perf_counter() - t0
        self.Logger.logger.info(
            '\n{} chunks downloaded in {:.2f}s'.format(count, elapsed))

        tidy_result(outputPath)
        return True

    def get_raw_ID_Mapping(self, outputPath, chunksize=100, concur_req=10):
        """
        Get Raw ID MApping Result
        """
        PARAMS['from'] = self.id_type
        status = self.get_info_from_uniprot(
            outputPath,
            self.dfrm[self.id_col].drop_duplicates(),
            chunksize=chunksize,
            concur_req=concur_req
        )

        if not status:
            return False

        new_colNames = [COLUMN_DICT.get(i, i) for i in self.usecols] + ['yourlist', 'isomap']
        dfrm = pd.read_csv(outputPath, sep='\t', names=new_colNames, skiprows=1, header=None)
        # WRITE REPORT
        self.Logger.logger.warning("\n%s\n" % str(dfrm.isnull().sum()))
        self.raw_id_mapping = dfrm
        return True

    @staticmethod
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
        dfrm['canonical_isoform'] = dfrm.apply(lambda x: ','.join(canonical_pattern.findall(
            x['Alternative products (isoforms)'])) if not isinstance(x['Alternative products (isoforms)'], float) else np.nan, axis=1)
        special_case = dfrm[dfrm['canonical_isoform'] == ''].index
        if len(special_case) > 0:
            canonical_pattern = re.compile(r'IsoId=([0-9A-Z-,\s]+); Sequence=Displayed')
            special_se = dfrm.loc[special_case].apply(lambda x: ','.join(
                canonical_pattern.findall(x['Alternative products (isoforms)'])), axis=1)
            dfrm.loc[special_case, 'canonical_isoform'] = special_se
            self.Logger.logger.warning("Special Cases of Canonical Info:")
            self.Logger.logger.warning("\n%s\n" % str(special_se))
        else:
            special_se = pd.Series([])

        return special_se

    def handle_ID_Mapping(self):
        self.Logger.logger.info("Start to handle id mapping result")
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
        df_wni_split = MapUniProtID.split_df(df_with_no_isomap, 'yourlist', ',')
        df_wni_split.drop(columns=['isomap'], inplace=True)
        # [yourlist <-> UniProt]
        df_wni_split['UniProt'] = df_wni_split['Entry']
        df_wni_split['unp_map_tage'] = 'Trusted & No Isoform'
        # Find out special cases 1
        df_wni_split_warn = df_wni_split[~df_wni_split['Alternative products (isoforms)'].isnull()].index
        df_wni_split.loc[df_wni_split_warn, 'unp_map_tage'] = 'Untrusted & No Isoform'
        # 'Entry', 'Gene names', 'Status', 'Alternative products (isoforms)', 'Organism', 'yourlist', 'UniProt'
        # ----------------------------------------------------------------------
        # In Class B
        # ----------------------------------------------------------------------
        wi_yourlist_count = df_with_isomap.apply(
            lambda x: x['yourlist'].count(','), axis=1)
        wi_isomap_count = df_with_isomap.apply(
            lambda x: x['isomap'].count(','), axis=1)
        # In subClass 1
        df_wi_eq = df_with_isomap.loc[wi_yourlist_count[wi_yourlist_count ==
                                                        wi_isomap_count].index]
        df_wi_eq_split = MapUniProtID.split_df(
            df_wi_eq.drop(columns=['yourlist']), 'isomap', ',')
        df_wi_eq_split['yourlist'], df_wi_eq_split['UniProt'] = df_wi_eq_split['isomap'].str.split(
            ' -> ', 1).str
        # [yourlist <-> UniProt]
        df_wi_eq_split.drop(columns=['isomap'], inplace=True)
        df_wi_eq_split['unp_map_tage'] = 'Trusted & Isoform'
        # # 'Entry', 'Gene names', 'Status', 'Alternative products (isoforms)', 'Organism', 'yourlist', 'UniProt'

        # In subClass 2
        df_wi_ne = df_with_isomap.loc[wi_yourlist_count[wi_yourlist_count !=
                                                        wi_isomap_count].index]
        if len(df_wi_ne) > 0:
            df_wi_ne_split = MapUniProtID.split_df(df_wi_ne, 'isomap', ',')
            df_wi_ne_split.rename(
                columns={'yourlist': 'checkinglist'}, inplace=True)
            df_wi_ne_split['yourlist'], df_wi_ne_split['UniProt'] = df_wi_ne_split['isomap'].str.split(
                ' -> ', 1).str
            df_wi_ne_split.drop(columns=['isomap'], inplace=True)
            df_wi_ne_split['unp_map_tage'] = 'Trusted & Isoform & Contain Warnings'
            # 'Entry', 'Gene names', 'Status', 'Alternative products (isoforms)', 'Organism', 'yourlist', 'UniProt', 'checkinglist'
            # Find out special cases 2
            usecols = pd.Index(set(df_wi_ne_split.columns) -
                               {'yourlist', 'UniProt'})
            df_wi_ne_warn = MapUniProtID.split_df(
                df_wi_ne_split[usecols].drop_duplicates(), 'checkinglist', ',')
            df_wi_ne_warn = df_wi_ne_warn[~df_wi_ne_warn['checkinglist'].isin(
                df_wi_ne_split['yourlist'])].rename(columns={'checkinglist': 'yourlist'})
            df_wi_ne_warn['UniProt'] = df_wi_ne_warn['Entry']
            df_wi_ne_warn['unp_map_tage'] = 'Untrusted & No Isoform'

            # Update UniProt
            final_df = pd.concat((df_wni_split, df_wi_eq_split, df_wi_ne_split.drop(
                columns=['checkinglist']), df_wi_ne_warn), sort=False).reset_index(drop=True)
        else:
            final_df = pd.concat((df_wni_split, df_wi_eq_split),
                                 sort=False).reset_index(drop=True)

        final_df['UniProt'] = final_df.apply(
            lambda x: x['Entry'] if x['UniProt'] == x['canonical_isoform'] else x['UniProt'], axis=1)
        if len(canonicalInfo_special_se) > 0:
            canonicalInfo_special_case = final_df[final_df['canonical_isoform'].isin(
                canonicalInfo_special_se)].index
            final_df.loc[canonicalInfo_special_case, 'UniProt'] = final_df.loc[canonicalInfo_special_case].apply(
                lambda x: x['Entry'] if x['UniProt'] in x['canonical_isoform'] else x['UniProt'], axis=1)
        return final_df

    def getGeneStatus(self, handled_df):
        """
        Will Change the dfrm, add Gene Status

        * Add new column (GENE) # if id_col != gene_col
        * Add new column (GENE_status)

        **About GENE_status**

        * ``False`` : First element of Gene names is not correspond with refSeq's GENE (e.g)
        * others(corresponding GENE)

        """
        if self.gene_col is None:
            handled_df['GENE_status'] = np.nan
            return None
        if self.id_col != self.gene_col:
            gene_map = self.dfrm[[self.id_col,
                                  self.gene_col]].drop_duplicates()
            gene_map.index = gene_map[self.id_col]
            gene_map.drop(columns=[self.id_col], inplace=True)
            handled_df['GENE'] = handled_df.apply(
                lambda z: gene_map[self.gene_col][z['yourlist']], axis=1)
            handled_df['GENE_status'] = handled_df.apply(lambda x: x['GENE'] == x['Gene names'].split(
                ' ')[0] if not isinstance(x['Gene names'], float) else False, axis=1)
        else:
            handled_df['GENE_status'] = handled_df.apply(lambda x: x['yourlist'] == x['Gene names'].split(
                ' ')[0] if not isinstance(x['Gene names'], float) else False, axis=1)

    def label_mapping_status(self, dfrm):
        dfrm['Mapping_status'] = 'No'        
        dfrm['GENE_status'] = dfrm['GENE_status'].apply(
            lambda x: x.any() if isinstance(x, Iterable) else x)
        pass_df = dfrm[
            (dfrm['GENE_status'] == True) &
            (dfrm['Status'] == 'reviewed') &
            (dfrm['unp_map_tage'] != 'Untrusted & No Isoform')]
        pass_index = pass_df.index
        dfrm.loc[pass_index, 'Mapping_status'] = 'Yes'

        # Deal with 'one to many' situation
        multipleCounter = Counter(dfrm.loc[pass_index, 'yourlist'])
        err_li = [i for i, j in multipleCounter.items() if j > 1]
        err_index = pass_df[pass_df['yourlist'].isin(err_li)].index
        dfrm.loc[err_index, 'Mapping_status'] = 'Error'

        # Write Report
        all_id = set(self.dfrm[self.id_col])
        self.Logger.logger.warning("All id: %s" % len(all_id))
        unmapped_id = all_id - set(dfrm['yourlist'])
        self.Logger.logger.warning("Unmapped id: %s" % len(unmapped_id))
        for i in sorted(unmapped_id):
            self.Logger.logger.warning("Unmapped: %s" % i)
        untrusted_id = set(dfrm[dfrm['Mapping_status'] != 'Yes']['yourlist']) - \
            set(dfrm[dfrm['Mapping_status'] == 'Yes']['yourlist'])
        self.Logger.logger.warning("Untrusted id: %s" % len(untrusted_id))
        for i in sorted(untrusted_id):
            self.Logger.logger.warning("Untrusted: %s" % i)
        self.Logger.logger.warning("Error id: %s" % len(
            dfrm.loc[err_index, 'yourlist'].drop_duplicates()))
        self.Logger.logger.warning("\n%s\n" % str(dfrm.loc[err_index]))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Commandline mode for ProcessUniProt')
    parser.add_argument('-r', '--referenceFile', type=str,
                        help='The reference file of IDs(with mutation Site) that need to map via UniProt RESTful API')
    parser.add_argument('-s', '--sep', type=str,
                        default='\t',
                        help='The seperator of referenceFile')
    parser.add_argument('-i', '--idCol', type=str, default='RefSeq_protein',
                        help='The column name of IDs in referenceFile.')
    parser.add_argument('-t', '--idType', type=str, default='P_REFSEQ_AC',
                        help='ID Abbreviation that stands for the type of ID.')
    parser.add_argument('-o', '--outputFolder', type=str,
                        help='Output file path.')
    parser.add_argument('-c', '--chunksize', type=int, default=100)
    parser.add_argument('-u', '--concurReq', type=int, default=10)
    parser.add_argument('-p', '--procceed', type=bool, default=True)
    parser.add_argument('-n', '--nrows', type=int, default=None)
    parser.add_argument('-f', '--finishedRaw', type=str, default='')
    args = parser.parse_args()

    rawPath = os.path.join(args.outputFolder, UniProt_ID_Mapping_RAW_PATH)
    modPath = os.path.join(args.outputFolder, UniProt_ID_Mapping_MODIFIED_PATH)
    logPath = os.path.join(args.outputFolder, LOGGER_PATH)

    demo = MapUniProtID(
        dfrm=pd.read_csv(args.referenceFile, sep=args.sep, nrows=args.nrows),
        id_col=args.idCol,
        id_type=args.idType,
        usecols='default',
        loggingPath=logPath,
        site_col='mutation_unp',
        gene_col='GENE'
    )

    if args.finishedRaw:
        new_colNames = [COLUMN_DICT.get(i, i) for i in COLUMNS] + ['yourlist', 'isomap']
        dfrm = pd.read_csv(args.finishedRaw, sep='\t', names=new_colNames, skiprows=1, header=None)
        demo.raw_id_mapping = dfrm
        status = True
    else:
        for index in range(2):
            status = demo.get_raw_ID_Mapping(rawPath, args.chunksize, args.concurReq)
            if not index:
                demo.Logger.logger.info(
                    "Will restart to fix unmapped id mapping result after 8 seconds...")
                time.sleep(8)
            else:
                demo.Logger.logger.info("Finish raw id mapping")
    
    if status and args.procceed:
        # Deal with different situations
        handled_ID_Mapping = demo.handle_ID_Mapping()
        # Add Gene Status
        demo.getGeneStatus(handled_ID_Mapping)
        # Label Mapping Status
        demo.label_mapping_status(handled_ID_Mapping)
        # Output the final result
        handled_ID_Mapping.to_csv(modPath, sep='\t', index=False)
    else:
        demo.Logger.logger.error("There is something wrong.")
