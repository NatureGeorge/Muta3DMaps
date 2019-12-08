# @Date:   2019-11-22T15:19:51+08:00
# @Email:  1730416009@stu.suda.edu.cn
# @Filename: ProcessUniProt.py
# @Last modified time: 2019-12-08T18:47:01+08:00
import os
import time
import sys
import asyncio
import aiohttp
import pandas as pd
from numpy import nan
from Logger import RunningLogger
import argparse
from random import uniform


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

LOGGER_PATH = 'downloads/async_processUnp.log'


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
            return await resp.read()

    async def download_one_chunk(self, session, chunkLystStr, path):
        rawData = await self.get_data(session, chunkLystStr)
        await asyncio.sleep(uniform(1, 5))
        self.Logger.logger.info(chunkLystStr[:50] + '...')
        self.save_data(path, rawData)
        return chunkLystStr

    async def download_many(self, chunkLystStr_list, path):
        async with aiohttp.ClientSession() as session:
            res = await asyncio.gather(
                *[asyncio.create_task(self.download_one_chunk(session, chunkLystStr, path))
                    for chunkLystStr in chunkLystStr_list])

        return len(res)

    def get_info_from_uniprot(self, outputPath, from_list=None, sep='\t', chunksize=100, header=None):

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
                finish['isomap'] = nan
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
                self.createLystStr(new_li, chunksize), outputPath))
        elapsed = time.perf_counter() - t0
        self.Logger.logger.info(
            '\n{} chunks downloaded in {:.2f}s'.format(count, elapsed))

        tidy_result(outputPath)
        return True

    def get_raw_ID_Mapping(self, outputPath, chunksize=100):
        """
        Get Raw ID MApping Result
        """
        PARAMS['from'] = self.id_type
        status = self.get_info_from_uniprot(
            outputPath,
            self.dfrm[self.id_col].drop_duplicates(),
            chunksize=chunksize
        )

        if not status:
            return False

        new_colNames = [COLUMN_DICT.get(i, i) for i in self.usecols] + ['yourlist', 'isomap']
        dfrm = pd.read_csv(outputPath, sep='\t', names=new_colNames, skiprows=1, header=None)
        # WRITE REPORT
        self.Logger.logger.warning("\n%s\n" % str(dfrm.isnull().sum()))
        self.raw_id_mapping = dfrm
        return True


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
    parser.add_argument('-o', '--output', type=str,
                        help='Output file path.')
    parser.add_argument('-c', '--chunksize', type=int, default=100)
    args = parser.parse_args()
    demo = MapUniProtID(
        dfrm=pd.read_csv(args.referenceFile, sep=args.sep),
        id_col=args.idCol,
        id_type=args.idType,
        usecols='default',
        loggingPath=LOGGER_PATH,
        site_col='mutation_unp',
        gene_col='GENE'
    )
    demo.get_raw_ID_Mapping(args.output, args.chunksize)
