# @Created Date: 2019-12-08 06:46:49 pm
# @Filename: ProcessUniProt.py
# @Email:  1730416009@stu.suda.edu.cn
# @Author: ZeFeng Zhu
# @Last Modified: 2019-12-22 04:37:24 pm
# @Copyright (c) 2019 MinghuiGroup, Soochow University
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
from collections.abc import Iterable, Iterator
import re
from Logger import RunningLogger
import json


COLUMNS = [
    'id', 'length', 'reviewed', 'comment(ALTERNATIVE%20PRODUCTS)',
    'feature(ALTERNATIVE%20SEQUENCE)', 'genes', 'organism', 'protein%20names'
    ]

COLUMN_DICT = {
    'id': 'Entry', 'length': 'Length', 'reviewed': 'Status',
    'comment(ALTERNATIVE%20PRODUCTS)': 'Alternative products (isoforms)',
    'feature(ALTERNATIVE%20SEQUENCE)': 'Alternative sequence',
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
SITE_INFO_PATH = "Site_info.tsv"

class ExtractIsoAlt:
    # pattern_all = re.compile(r"([A-z0-9]+):VAR_SEQ\s([0-9\s]+)\s([A-z\s\->]+)\s\(in\s([^\)]+)\)[^=]+FTId=([A-z0-9_,\s]+)")
    pattern_sep = re.compile(r"([A-z0-9]+):VAR_SEQ\s([0-9\.]+);\s+")
    pattern_iso_Value = re.compile(r'/([^=]+)="([^;]+)"')
    pattern_inIso = re.compile(r"([A-z\s\->]+)\s\(in\s([^\)]+)\)")
    pattern_iso_keyValue = re.compile(r"([^=]+)=([^;]+);\s+")
    pattern_inDe = re.compile(r"([A-z]+)\s->\s([A-z]+)")
    usecols = ["Entry", "Alternative sequence",
               "Alternative products (isoforms)"]

    def __init__(self, path, outFolder, usecols=usecols):
        self.dfrm = pd.read_csv(path, sep="\t", usecols=usecols)
        self.dfrm.dropna(inplace=True)
        self.dfrm.drop_duplicates(inplace=True)
        self.dfrm.reset_index(drop=True, inplace=True)
        self.altSeqPath = os.path.join(outFolder, "UniProt_AltSeq_info.tsv")
        self.altProPath = os.path.join(outFolder, "UniProt_AltPro_info.tsv")

    def __len__(self):
        return len(self.dfrm)

    def extractAltSeq(self):
        target_col = "Alternative sequence"
        self.dfrm[target_col] = self.dfrm.apply(lambda x: x[target_col].replace(
            "VAR_SEQ", "{}:VAR_SEQ".format(x["Entry"])), axis=1)
        target_col = "Alternative products (isoforms)"
        self.dfrm[target_col] = self.dfrm.apply(
            lambda x: "Entry={}; {}".format(x["Entry"], x[target_col]), axis=1)

        '''
        altSeq_li = []
        for i in self.dfrm["Alternative sequence"].dropna():
            find = self.pattern_all.findall(i)
            altSeq_li.extend(find)
        
        altSeq_df = pd.DataFrame(
            altSeq_li, columns=["Entry", "AltRange", "AltInfo", "AltIso", "AltID"])
        '''
        altSeq_li = []
        for content in self.dfrm["Alternative sequence"].dropna():
            result = self.pattern_sep.split(content)
            for i in range(1, len(result)-1, 3):
                altSeq_li.append(
                    result[i+2] + '; /Entry="%s"; /AltRange="%s"; ' % tuple(result[i:i+2]))
        
        altSeq_df = pd.DataFrame(dict(i.groups() for i in self.pattern_iso_Value.finditer(content)) for content in altSeq_li)
        altSeq_df.rename(columns={"id": "AltID"}, inplace=True)
        altSeq_df["AltInfo"], altSeq_df["AltIso"] = altSeq_df["note"].apply(lambda x: self.pattern_inIso.search(x).groups()).str
        '''
            for i in range(0, len(result)-1, 3):
                temp = result[1:][i:i+3]
                cur = temp.pop(2)
                temp.extend(self.pattern_iso_Value.findall(cur))
                cur = temp.pop(2)
                temp.extend(list(self.pattern_inIso.search(cur).groups()))
                if len(temp) == 6:
                    altSeq_li.append(temp)
                elif len(temp) == 5:
                    temp.insert(2, '')
                    altSeq_li.append(temp)
                else:
                    raise Exception(content)
     
        altSeq_df = pd.DataFrame(
            altSeq_li, columns=["Entry", "AltRange", "AltNote", "AltID", "AltInfo", "AltIso"])
        '''
        altSeq_df["AltRange"] = altSeq_df["AltRange"].apply(
            lambda x: [int(i) for i in x.split("..")])

        altSeq_df["AltLen"] = altSeq_df.apply(lambda x: [len(i) for i in self.pattern_inDe.search(
            x["AltInfo"]).groups()] if x["AltInfo"] != "Missing" else len(range(*x["AltRange"]))+1, axis=1)
        altSeq_df.to_csv(self.altSeqPath, sep="\t", index=False)
        self.altSeq_dict = altSeq_df[[
            "AltID", "AltLen", "AltRange"]].to_dict("list")

    @staticmethod
    def getAltProInfo(inputLyst, groupCol="isoforms", flagCol="Name"):
        def toFrame(dict):
            df = pd.DataFrame(dict[groupCol])
            for key, value in dict.items():
                if key == groupCol:
                    continue
                else:
                    df[key] = value
            return df

        outputDict = {groupCol: []}
        flag = 0
        for key, value in inputLyst:
            if not flag and key != flagCol:
                outputDict[key] = value
            if key == flagCol:
                outputDict[groupCol].append({key: value})
                flag += 1
            else:
                if flag:
                    outputDict[groupCol][flag-1][key] = value

        return toFrame(outputDict)

    @staticmethod
    def getAltInterval(alt_id, altSeq_dict):
        if alt_id in ["Displayed", "External", "Not described"]:
            return np.nan, np.nan
        elif not alt_id.startswith("VSP"):
            raise ValueError("Unexcepted alt_id: %s" % alt_id)
        else:
            mis_info = []
            inDe_info = []
            alt_id_lyst = alt_id.split(", ")
            for altID in alt_id_lyst:
                index = altSeq_dict["AltID"].index(altID)
                len_info = altSeq_dict["AltLen"][index]
                range_info = altSeq_dict["AltRange"][index]
                if isinstance(len_info, int):
                    mis_info.append((range_info, len_info))
                else:
                    inDe_info.append((range_info, len_info))
                    if len_info[0] > len_info[1]:
                        mis_info.append(
                            ([range_info[0], range_info[1]-1], len_info[0]-len_info[1]))

        return mis_info, inDe_info

    @staticmethod
    def getAffectedInterval(mis_info, inDe_info):
        if isinstance(mis_info, float):
            return np.nan
        affect_Interval = []
        for inDeRange, inDeLen in inDe_info:
            start = inDeRange[0]
            affect_Interval.append((start, start+inDeLen[1]-1))

        for index in range(len(affect_Interval)):
            raw = affect_Interval[index]
            interval = list(raw)
            for misRange, misLen in mis_info:
                if misRange[-1] < raw[0]:
                    interval = [x-misLen for x in interval]
            affect_Interval[index] = interval

        if affect_Interval:
            return affect_Interval
        else:
            return np.nan

    def extractAltPro(self):
        altPro_df = pd.concat((self.getAltProInfo(self.pattern_iso_keyValue.findall(
            i)) for i in self.dfrm["Alternative products (isoforms)"].dropna()), sort=False)
        altPro_df["AltInterval"] = altPro_df["Sequence"].apply(
            lambda x: self.getAffectedInterval(*self.getAltInterval(x, self.altSeq_dict)))
        altPro_df["UniProt"] = altPro_df.apply(lambda x: x["IsoId"].split(
            ",")[0] if x["Sequence"] != "Displayed" else x["Entry"], axis=1)
        altPro_df.to_csv(self.altProPath, sep="\t", index=False)

    @classmethod
    def main(cls, **kwargs):
        demo = cls(**kwargs)
        demo.extractAltSeq()
        demo.extractAltPro()


def asIter(func, typeTp=[str]):
    """
    Always let the args are all Iterable/Iterator object
    """
    typeTp = tuple(typeTp)

    def register(*args):
        new_args = []
        for arg in args:
            if isinstance(arg, (Iterable, Iterator)) or arg is None:
                if isinstance(arg, typeTp):
                    new_args.append([arg])
                else:
                    new_args.append(arg)
            else:
                new_args.append([arg])
        new_args = tuple(new_args)
        return func(*new_args)

    return register


class MapUniProtID:

    def __init__(self, dfrm, id_col, id_type, loggingPath, site_col=None, site_type=None, gene_col=None, usecols=None):
        self.dfrm = dfrm
        self.index = dfrm.index
        self.id_col = id_col
        self.id_type = id_type
        self.site_col = site_col
        self.site_type = site_type
        self.gene_col = gene_col
        self.usecols = usecols
        self.Logger = RunningLogger("MapUniProtID", loggingPath)

    @classmethod
    def fromDfrm(cls, **kwargs):
        return cls(**kwargs)

    @classmethod
    def fromIter(cls, ids, id_type, loggingPath, sites=None, site_type=None, genes=None, usecols=None, **kwargs):
        @asIter
        def forIter(*args):
            return args

        ids, sites, genes = forIter(ids, sites, genes)
        id_col = kwargs.get("id_col", "id")
        site_col = kwargs.get("site_col", "site")
        gene_col = kwargs.get("gene_col", "gene")
        dfrm = pd.DataFrame({
            id_col: ids,
            site_col: sites,
            gene_col: genes
        })
        return cls(dfrm, id_col, id_type, loggingPath, site_col, gene_col, usecols)

    @property
    def site_li(self):
        if self.site_col is not None:
            # return self.dfrm.groupby(by=[self.id_col]).apply(lambda x: [i for i in x[self.site_col]])
            for name, group in self.dfrm.groupby(by=[self.id_col]):
                yield name, [i for i in group[self.site_col]]
        else:
            return None

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
            status = index % 5
            if status:
                # sec = uniform(5, 10)
                sec = (index/5)**1.2
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

        if self.usecols is not None:
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

    def getCanonicalInfo(self, dfrm, folder):
        """
        Will Change the dfrm

        * Add new column (canonical_isoform)
        * Change the content of column (UniProt)
        """
        # Get info from Alt Product file
        usecols = ["IsoId", "Sequence", "Entry", "UniProt"]
        altPro_path = os.path.join(folder, "UniProt_AltPro_info.tsv")
        altPro_df = pd.read_csv(altPro_path, sep="\t", usecols=usecols)
        altPro_df = altPro_df[altPro_df["Sequence"] == "Displayed"].reset_index(drop=True)
        altPro_df.rename(columns={"IsoId": "canonical_isoform"}, inplace=True)
        # altPro_df["canonical_isoform"] = altPro_df.apply(lambda x: x.split(", ")[0] if "," in x else x)
        # Modify dfrm
        dfrm = pd.merge(dfrm, altPro_df[["canonical_isoform", "Entry"]], how="left")
        '''
        cano_index = dfrm[dfrm["canonical_isoform"].notnull()].index
        dfrm.loc[cano_index, "UniProt"] = dfrm.loc[cano_index, ].apply(
            lambda x: x["Entry"] if x["UniProt"] in x["canonical_isoform"] else x["UniProt"], axis=1)
        '''
        return dfrm

    def handle_ID_Mapping(self, folder):
        self.Logger.logger.info("Start to handle id mapping result")
        df = self.raw_id_mapping
        # Add New Column: canonical_isoform
        df = self.getCanonicalInfo(df, folder)
        # Add New Column: unp_map_tage
        df['unp_map_tage'] = np.nan
        # Classification
        df_with_no_isomap = df[df['isomap'].isnull()]  # Class A
        df_with_isomap = df[df['isomap'].notnull()]  # Class B
        # ----------------------------------------------------------------------
        # In Class A
        # ----------------------------------------------------------------------
        if len(df_with_no_isomap) > 0:
            df_wni_split = self.split_df(df_with_no_isomap, 'yourlist', ',')
            df_wni_split.drop(columns=['isomap'], inplace=True)
            # [yourlist <-> UniProt]
            df_wni_split['UniProt'] = df_wni_split['Entry']
            df_wni_split['unp_map_tage'] = 'Trusted & No Isoform'
            # Find out special cases 1
            df_wni_split_warn = df_wni_split[df_wni_split['Alternative products (isoforms)'].notnull()].index
            df_wni_split.loc[df_wni_split_warn, 'unp_map_tage'] = 'Untrusted & No Isoform'
            # 'Entry', 'Gene names', 'Status', 'Alternative products (isoforms)', 'Organism', 'yourlist', 'UniProt'
        # ----------------------------------------------------------------------
        # In Class B
        # ----------------------------------------------------------------------
        if len(df_with_isomap) > 0:
            wi_yourlist_count = df_with_isomap.apply(
                lambda x: x['yourlist'].count(','), axis=1)
            wi_isomap_count = df_with_isomap.apply(
                lambda x: x['isomap'].count(','), axis=1)
            # In subClass 1
            df_wi_eq = df_with_isomap.loc[wi_yourlist_count[wi_yourlist_count == wi_isomap_count].index]
            if len(df_wi_eq) > 0:
                df_wi_eq_split = self.split_df(
                    df_wi_eq.drop(columns=['yourlist']), 'isomap', ',')
                df_wi_eq_split['yourlist'], df_wi_eq_split['UniProt'] = df_wi_eq_split['isomap'].str.split(
                    ' -> ', 1).str
                # [yourlist <-> UniProt]
                df_wi_eq_split.drop(columns=['isomap'], inplace=True)
                df_wi_eq_split['unp_map_tage'] = 'Trusted & Isoform'
                # # 'Entry', 'Gene names', 'Status', 'Alternative products (isoforms)', 'Organism', 'yourlist', 'UniProt'

            # In subClass 2
            df_wi_ne = df_with_isomap.loc[wi_yourlist_count[wi_yourlist_count != wi_isomap_count].index]
            if len(df_wi_ne) > 0:
                df_wi_ne_split = self.split_df(df_wi_ne, 'isomap', ',')
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
                df_wi_ne_warn = self.split_df(
                    df_wi_ne_split[usecols].drop_duplicates(), 'checkinglist', ',')
                df_wi_ne_warn = df_wi_ne_warn[~df_wi_ne_warn['checkinglist'].isin(
                    df_wi_ne_split['yourlist'])].rename(columns={'checkinglist': 'yourlist'})
                df_wi_ne_warn['UniProt'] = df_wi_ne_warn['Entry']
                # sequence conflict
                df_wi_ne_warn['unp_map_tage'] = 'Untrusted & No Isoform'
                df_wi_ne_split.drop(columns=['checkinglist'], inplace=True)

        # Concat Dfrm
        variables = ["df_wni_split", "df_wi_eq_split", "df_wi_ne_split", "df_wi_ne_warn"]
        lvs = locals()
        varLyst = [lvs[variable] for variable in variables if variable in lvs]
        final_df = pd.concat(varLyst, sort=False).reset_index(drop=True)
        cano_index = final_df[final_df["canonical_isoform"].notnull()].index
        final_df.loc[cano_index, "UniProt"] = final_df.loc[cano_index, ].apply(
            lambda x: x["Entry"] if x["UniProt"] in x["canonical_isoform"] else x["UniProt"], axis=1)
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
        if self.id_col != 'GENENAME':
            
            if self.gene_col is None:
                handled_df['GENE_status'] = True
                return None

            gene_map = self.dfrm[[self.id_col,
                                  self.gene_col]].drop_duplicates()
            gene_map = gene_map.groupby(self.id_col)[self.gene_col].apply(
                lambda x: np.array(x) if len(x) > 1 else list(x)[0])
            handled_df['GENE'] = handled_df.apply(lambda z: gene_map[z['yourlist']], axis=1)
            handled_df['GENE_status'] = handled_df.apply(lambda x: x['GENE'] == x['Gene names'].split(
                ' ')[0] if not isinstance(x['Gene names'], float) else False, axis=1)
            handled_df['GENE'] = handled_df['GENE'].apply(lambda x: ','.join(x) if not isinstance(x, str) else x)
        else:
            handled_df['GENE_status'] = handled_df.apply(lambda x: x['yourlist'] == x['Gene names'].split(
                ' ')[0] if not isinstance(x['Gene names'], float) else False, axis=1)

    def label_mapping_status(self, dfrm):
        dfrm['Mapping_status'] = 'No'        
        dfrm['GENE_status'] = dfrm['GENE_status'].apply(
            lambda x: x.any() if isinstance(x, Iterable) else x)
        
        if self.id_col == 'GENENAME':
            pass_df = dfrm[
                (dfrm['GENE_status'] == True) &
                (dfrm['Status'] == 'reviewed') &
                (dfrm['unp_map_tage'] != 'Untrusted & No Isoform')]
        else:
            pass_df = dfrm[
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

    @classmethod
    def main(cls, **kwargs):
        defaultMode = "fromDfrm"
        mode = kwargs.get("mode", defaultMode)
        finishedRaw = kwargs.get("finishedRaw", "")
        outFolder = kwargs["outputFolder"]
        rawPath = os.path.join(outFolder, UniProt_ID_Mapping_RAW_PATH)
        modPath = os.path.join(outFolder, UniProt_ID_Mapping_MODIFIED_PATH)
        logPath = os.path.join(outFolder, LOGGER_PATH)
        id_col = kwargs.get("id_col", "id")
        site_col = kwargs.get("site_col", "site")
        gene_col = kwargs.get("gene_col", "gene")

        usecols = [col for col in [id_col, site_col, gene_col] if col is not None]
        
        if mode == defaultMode:
            demo = cls(
                dfrm=pd.read_csv(kwargs["referenceFile"], sep=kwargs["sep"], nrows=kwargs["nrows"], usecols=usecols),
                id_col=id_col,
                id_type=kwargs["id_type"],
                loggingPath=logPath,
                site_col=site_col,
                gene_col=gene_col,
                usecols=kwargs.get("usecols", None),
            )
        elif mode == "fromIter":
            kwargs["loggingPath"] = logPath
            demo = cls.fromIter(**kwargs)
        
        if finishedRaw:
            new_colNames = [COLUMN_DICT.get(i, i)
                            for i in COLUMNS] + ['yourlist', 'isomap']
            dfrm = pd.read_csv(finishedRaw, sep='\t',
                            names=new_colNames, skiprows=1, header=None)
            demo.raw_id_mapping = dfrm
            status = True
        else:
            chunksize = kwargs.get("chunksize", 100)
            concurReq = kwargs.get("concurReq", 20)
            # I am not sure whether this step should exist?
            for index in range(2):
                status = demo.get_raw_ID_Mapping(
                    rawPath, chunksize, concurReq)
                if not index:
                    demo.Logger.logger.info(
                        "Will restart to fix unmapped id mapping result after 8 seconds...")
                    time.sleep(8)
                else:
                    demo.Logger.logger.info("Finish raw id mapping")
            # [Can be done parallelly]
            if site_col is not None:
                siteInfoPath = os.path.join(outFolder, SITE_INFO_PATH)
                with open(siteInfoPath, "w+") as siteOutput:
                    siteOutput.write("id\tsite\n")
                    formatStr = "%s\t%s\n"
                    for name, lyst in demo.site_li:
                        siteOutput.write(formatStr % (name, json.dumps(lyst)))
                demo.Logger.logger.info(
                    "Site Info has been safed in %s" % siteInfoPath)

        
        # Extract Alt Seq/Pro Info [Can be done parallelly]
        ExtractIsoAlt.main(path=rawPath, outFolder=outFolder)

        if status and kwargs.get("procceed", True):
            # Deal with different situations
            handled_ID_Mapping = demo.handle_ID_Mapping(outFolder)
            # Add Gene Status
            demo.getGeneStatus(handled_ID_Mapping)
            # Label Mapping Status
            demo.label_mapping_status(handled_ID_Mapping)
            # Output the final result
            handled_ID_Mapping.to_csv(modPath, sep='\t', index=False)
        else:
            demo.Logger.logger.error("There is something wrong.")
        


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Commandline mode for ProcessUniProt')
    parser.add_argument('-r', '--referenceFile', type=str, help='The reference file of IDs(with mutation Site) that need to map via UniProt RESTful API')
    parser.add_argument('-s', '--sep', type=str, default='\t', help='The seperator of referenceFile')
    parser.add_argument('-i', '--idCol', type=str, default='RefSeq_protein', help='The column name of IDs in referenceFile.')
    parser.add_argument('-t', '--idType', type=str, default='P_REFSEQ_AC', help='ID Abbreviation that stands for the type of ID.')
    parser.add_argument('-o', '--outputFolder', type=str, help='Output file path.')
    parser.add_argument('-c', '--chunksize', type=int, default=100)
    parser.add_argument('-u', '--concurReq', type=int, default=20)
    parser.add_argument('-p', '--procceed', type=bool, default=True)
    parser.add_argument('-n', '--nrows', type=int, default=None)
    parser.add_argument('-f', '--finishedRaw', type=str, default='')
    parser.add_argument('-l', '--siteCol', type=str, default=None)
    parser.add_argument('-g', '--geneCol', type=str, default=None)
    # parser.add_argument('-m', '--mode', type=str, default='fromDfrm')
    args = parser.parse_args()
    
    MapUniProtID.main(finishedRaw=args.finishedRaw, outputFolder=args.outputFolder,
                      referenceFile=args.referenceFile, sep=args.sep,
                      nrows=args.nrows, id_col=args.idCol,
                      id_type=args.idType, site_col=args.siteCol,
                      gene_col=args.geneCol, procceed=args.procceed,
                      chunksize=args.chunksize, concurReq=args.concurReq)
    '''
    MapUniProtID.main(
        mode="fromIter", 
        outputFolder="C:\\OmicData\\LiGroupWork\\compareTool\\1222\\fromIter",
        ids=["NP_689699", "NP_689699", "NP_940978"],
        genes=["SAMD11", "SAMD11", "AGRN"], 
        sites=["K45E", "P293A", "G76S"], 
        id_type="P_REFSEQ_AC")
    '''
