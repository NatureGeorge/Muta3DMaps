# @Created Date: 2019-12-08 06:46:49 pm
# @Filename: decode.py
# @Email:  1730416009@stu.suda.edu.cn
# @Author: ZeFeng Zhu
# @Last Modified: 2020-02-16 10:54:32 am
# @Copyright (c) 2020 MinghuiGroup, Soochow University
from typing import Iterable, Iterator, Optional, Union, Generator, Dict, List
import re, time, sys
import ujson as json
import pandas as pd
import numpy as np
import logging
from pathlib import Path
from unsync import unsync, Unfuture
from collections import Counter
sys.path.append(r'C:\GitWorks\Muta3DMaps')
from Muta3DMaps.core.retrieve.fetchFiles import UnsyncFetch

QUERY_COLUMNS: List[str] = [
    'id', 'length', 'reviewed', 
    'comment(ALTERNATIVE%20PRODUCTS)',
    'feature(ALTERNATIVE%20SEQUENCE)', 
    'genes', 'organism', 'protein%20names']

RESULT_COLUMNS: List[str] = [
    'Entry', 'Length', 'Status',
    'Alternative products (isoforms)',
    'Alternative sequence',
    'Gene names', 'Organism', 'Protein names']

COLUMNS_DICT: Dict = dict(zip(QUERY_COLUMNS, RESULT_COLUMNS))

RESULT_NEW_COLUMN: List[str] = ['yourlist', 'isomap']

BASE_URL: str = 'https://www.uniprot.org/uploadlists/'

PARAMS: Dict = {
    'columns': None,
    'query': None,
    'from': None,
    'to': 'ACC',
    'format': 'tab'}


class ExtractIsoAlt:
    # pattern_all = re.compile(r"([A-z0-9]+):VAR_SEQ\s([0-9\s]+)\s([A-z\s\->]+)\s\(in\s([^\)]+)\)[^=]+FTId=([A-z0-9_,\s]+)")
    pattern_sep = re.compile(r"([A-z0-9]+):VAR_SEQ\s([0-9\.]+);\s+")
    pattern_iso_Value = re.compile(r'/([^=]+)="([^;]+)"')
    pattern_inIso = re.compile(r"([A-z\s\->]+)\s\(in\s([^\)]+)\)")
    pattern_iso_keyValue = re.compile(r"([^=]+)=([^;]+);\s+")
    pattern_inDe = re.compile(r"([A-z]+)\s->\s([A-z]+)")
    usecols = ["Entry", "Alternative sequence", "Alternative products (isoforms)"]

    def __init__(self, path: str, usecols=usecols, sep: str = '\t'):
        pathOb = Path(path)
        prefix = pathOb.stem
        self.dfrm = pd.read_csv(path, sep=sep, usecols=usecols)
        self.dfrm.dropna(inplace=True)
        self.dfrm.drop_duplicates(inplace=True)
        self.dfrm.reset_index(drop=True, inplace=True)
        self.altSeqPath = str(Path(pathOb.parent, f"{prefix}_UNP_AltSeq.tsv"))
        self.altProPath = str(Path(pathOb.parent, f"{prefix}_UNP_AltPro.tsv"))

    def __len__(self):
        return len(self.dfrm)

    def extractAltSeq(self):
        target_col_1 = "Alternative sequence"
        self.dfrm[target_col_1] = self.dfrm.apply(lambda x: x[target_col_1].replace(
            "VAR_SEQ", "{}:VAR_SEQ".format(x["Entry"])), axis=1)
        target_col_2 = "Alternative products (isoforms)"
        self.dfrm[target_col_2] = self.dfrm.apply(
            lambda x: "Entry={}; {}".format(x["Entry"], x[target_col_2]), axis=1)

        altSeq_li = []
        for content in self.dfrm[target_col_1].dropna():
            result = self.pattern_sep.split(content)
            for i in range(1, len(result)-1, 3):
                altSeq_li.append(
                    result[i+2] + '; /Entry="%s"; /AltRange="%s"; ' % tuple(result[i:i+2]))

        altSeq_df = pd.DataFrame(dict(i.groups() for i in self.pattern_iso_Value.finditer(
            content)) for content in altSeq_li)
        altSeq_df.rename(columns={"id": "AltID"}, inplace=True)
        altSeq_df["AltInfo"], altSeq_df["AltIso"] = altSeq_df["note"].apply(
            lambda x: self.pattern_inIso.search(x).groups()).str
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
        if len(demo.dfrm) > 0:
            demo.extractAltSeq()
            demo.extractAltPro()
            return demo.altSeqPath, demo.altProPath
        else:
            return None, None

class MapUniProtID(object):
    '''
    Implement UniProt Retrieve/ID Mapping API
    '''
    logger = logging.getLogger("MapUniProtID")
    logger.setLevel(logging.DEBUG)
    streamHandler = logging.StreamHandler(); streamHandler.setLevel(logging.INFO)
    formatter = logging.Formatter("%(asctime)s %(name)s %(levelname)s %(message)s"); streamHandler.setFormatter(formatter)
    logger.addHandler(streamHandler)

    def __init__(self, id_col: str, id_type: str, 
        dfrm: Optional[pd.DataFrame], 
        ids: Optional[Iterable] = None,
        sites: Optional[Iterable] = None,
        genes: Optional[Iterable] = None,
        usecols: Optional[Iterable] = QUERY_COLUMNS,
        site_col: Optional[str] = None, 
        site_type: Optional[str] = None, 
        gene_col: Optional[str] = None,  
        loggingPath: Optional[str] = None):
        
        if dfrm is not None:
            self.dfrm = dfrm
        else:
            '''
            the length of dataframe is based on:
            
                * the num of `ids` if there is more than one id
                * the num of `sites` if there is just one id with specified `sites`
            '''
            if isinstance(ids, str):
                if sites is not None and not isinstance(sites, str):
                    index_len = len(sites)
                else:
                    index_len = 1
            else:
                index_len = len(ids)
            
            self.dfrm = pd.DataFrame(dict(zip(
                (col for col in (id_col, site_col, gene_col) if col is not None),
                (value for value in (ids, sites, genes) if value is not None))), 
                index=list(range(index_len)))
        
        self.index = dfrm.index
        self.id_col = id_col
        self.id_type = id_type
        self.site_col = site_col
        self.site_type = site_type
        self.gene_col = gene_col
        self.usecols = usecols
        self.loggingPath = loggingPath
        PARAMS['columns'] = ','.join(self.usecols)
        PARAMS['from'] = id_type
        if isinstance(loggingPath, str):
            self.set_logging_fileHandler(loggingPath)

    @classmethod
    def set_logging_fileHandler(cls, path: str, level: int = logging.DEBUG, formatter=formatter):
        try:
            fileHandler = logging.FileHandler(filename=path)
            fileHandler.setLevel(level)
            fileHandler.setFormatter(formatter)
            cls.logger.addHandler(fileHandler)
            cls.logger.info(f"Logging file in {path}")
        except Exception:
            cls.logger.warning(
                "Invalid file path for logging file ! Please specifiy path=...")

    @property
    def sites(self) -> Generator:
        if self.site_col is not None:
            for name, group in self.dfrm.groupby(by=self.id_col, sort=False):
                yield name, group[self.site_col]
        else:
            yield None

    @staticmethod
    def split_df(dfrm, colName, sep):
        """Split DataFrame"""
        df = dfrm.copy()
        return df.drop([colName], axis=1).join(df[colName].str.split(sep, expand=True).stack().reset_index(level=1, drop=True).rename(colName))

    def yieldTasks(self, lyst: Iterable, chunksize: int = 100, sep: str = ',') -> Generator:
        fileName = self.outputPath.stem
        for i in range(0, len(lyst), chunksize):
            cur_fileName = f'{fileName}_{i}'
            cur_params = PARAMS.copy()
            cur_params['query'] = sep.join(lyst[i:i+chunksize])
            yield ('get', {'url': BASE_URL, 'params': cur_params}, str(Path(self.outputPath.parent, cur_fileName+self.outputPath.suffix)))

    def retrieve(self, outputPath: str, finishedPath: Optional[str] = None, sep: str = '\t', chunksize: int = 100, concur_req: int = 20, rate: float = 1.5):
        finish_id = list()
        self.outputPath = Path(outputPath)
        self.result_cols = [COLUMNS_DICT.get(i, i) for i in self.usecols] + RESULT_NEW_COLUMN
        if finishedPath is not None:
            try:
                target_col = RESULT_NEW_COLUMN[0]
                finish: pd.Series = pd.read_csv(
                    finishedPath,
                    sep=sep, 
                    usecols=[target_col],
                    names=self.result_cols, 
                    skiprows=1, 
                    header=None)[target_col]
            except Exception as e:
                col_to_add = RESULT_NEW_COLUMN[1]
                self.logger.warning(f"{e}\nSomething wrong with finished raw file, probably without '{col_to_add}' column.")
                finish_df = pd.read_csv(
                    finishedPath, sep=sep, names=self.result_cols[:-1], skiprows=1, header=None)
                finish_df[col_to_add] = np.nan
                finish_df.to_csv(finishedPath, sep=sep, index=False)
                finish: pd.Series = finish_df[target_col]

            for query_id in finish:
                if ',' in query_id:
                    finish_id.extend(query_id.split(','))
                else:
                    finish_id.append(query_id)

        query_id: pd.Series = self.dfrm[self.id_col]
        if finish_id:
            rest_id = list(set(query_id) - set(finish_id))
        else:
            rest_id = query_id.unique()
        
        self.logger.info(f"Have finished {len(finish_id)} ids, {len(rest_id)} ids left.")
        t0 = time.perf_counter()
        res = UnsyncFetch.multi_tasks(self.yieldTasks(rest_id, chunksize), self.process, concur_req, rate).result()
        elapsed = time.perf_counter() - t0
        self.logger.info('{} chunks downloaded in {:.2f}s'.format(len(res), elapsed))

    def getCanonicalInfo(self, dfrm: pd.DataFrame):
        """
        Will Change the dfrm

        * Add new column (canonical_isoform)
        * Change the content of column (UniProt)
        """
        # Get info from Alt Product file
        if self.altProPath is None:
            dfrm['canonical_isoform'] = np.nan
            return dfrm
        else:
            usecols = ["IsoId", "Sequence", "Entry", "UniProt"]
            altPro_df = pd.read_csv(self.altProPath, sep="\t", usecols=usecols)
            altPro_df = altPro_df[altPro_df["Sequence"] == "Displayed"].reset_index(drop=True)
            altPro_df.rename(columns={"IsoId": "canonical_isoform"}, inplace=True)
            # Modify dfrm
            dfrm = pd.merge(dfrm, altPro_df[["canonical_isoform", "Entry"]], how="left")
            return dfrm

    def getGeneStatus(self, handled_df: pd.DataFrame, colName: str = 'GENE_status'):
        """
        Will Change the dfrm, add Gene Status

        * Add new column (GENE) # if id_col != gene_col
        * Add new column (GENE_status)

        **About GENE_status**

        * ``False`` : First element of Gene names is not correspond with refSeq's GENE (e.g)
        * others(corresponding GENE)

        """
        self.gene_status_col = colName
        if self.id_col != 'GENENAME':

            if self.gene_col is None:
                handled_df[colName] = True
                return None

            gene_map = self.dfrm[[self.id_col,
                                  self.gene_col]].drop_duplicates()
            gene_map = gene_map.groupby(self.id_col)[self.gene_col].apply(
                lambda x: np.array(x) if len(x) > 1 else list(x)[0])
            handled_df['GENE'] = handled_df.apply(
                lambda z: gene_map[z['yourlist']], axis=1)
            handled_df[colName] = handled_df.apply(lambda x: x['GENE'] == x['Gene names'].split(
                ' ')[0] if not isinstance(x['Gene names'], float) else False, axis=1)
            handled_df['GENE'] = handled_df['GENE'].apply(
                lambda x: ','.join(x) if not isinstance(x, str) else x)
        else:
            handled_df[colName] = handled_df.apply(lambda x: x['yourlist'] == x['Gene names'].split(
                ' ')[0] if not isinstance(x['Gene names'], float) else False, axis=1)

    def label_mapping_status(self, dfrm: pd.DataFrame, colName: str = 'Mapping_status'):
        self.mapping_status_col = colName
        gene_status_col = self.gene_status_col
        dfrm[colName] = 'No'
        dfrm[gene_status_col] = dfrm[gene_status_col].apply(
            lambda x: x.any() if isinstance(x, Iterable) else x)

        if self.id_col == 'GENENAME':
            pass_df = dfrm[
                (dfrm[gene_status_col] == True) &
                (dfrm['Status'] == 'reviewed') &
                (dfrm['unp_map_tage'] != 'Untrusted & No Isoform')]
        else:
            pass_df = dfrm[
                (dfrm['Status'] == 'reviewed') &
                (dfrm['unp_map_tage'] != 'Untrusted & No Isoform')]
        pass_index = pass_df.index
        dfrm.loc[pass_index, colName] = 'Yes'

        # Deal with 'one to many' situation
        multipleCounter = Counter(dfrm.loc[pass_index, 'yourlist'])
        err_li = [i for i, j in multipleCounter.items() if j > 1]
        err_index = pass_df[pass_df['yourlist'].isin(err_li)].index
        dfrm.loc[err_index, colName] = 'Error'

    @unsync
    def process(self, path: Union[str, Unfuture], sep: str = '\t'):
        self.logger.debug("Start to handle id mapping result")
        if not isinstance(path, str):
            path = path.result()
        self.altSeqPath, self.altProPath = ExtractIsoAlt.main(path=path)
        try:
            df = pd.read_csv(path, sep='\t', names=self.result_cols, skiprows=1, header=None)
        except ValueError:
            df = pd.read_csv(path, sep='\t', names=self.result_cols[:-1], skiprows=1, header=None)

        # Add New Column: canonical_isoform
        df = self.getCanonicalInfo(df)
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
            df_wni_split_warn = df_wni_split[df_wni_split['Alternative products (isoforms)'].notnull(
            )].index
            df_wni_split.loc[df_wni_split_warn,
                             'unp_map_tage'] = 'Untrusted & No Isoform'
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
            df_wi_eq = df_with_isomap.loc[wi_yourlist_count[wi_yourlist_count ==
                                                            wi_isomap_count].index]
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
            df_wi_ne = df_with_isomap.loc[wi_yourlist_count[wi_yourlist_count !=
                                                            wi_isomap_count].index]
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
        variables = ["df_wni_split", "df_wi_eq_split",
                     "df_wi_ne_split", "df_wi_ne_warn"]
        lvs = locals()
        varLyst = [lvs[variable] for variable in variables if variable in lvs]
        final_df = pd.concat(varLyst, sort=False).reset_index(drop=True)
        cano_index = final_df[final_df["canonical_isoform"].notnull()].index
        final_df.loc[cano_index, "UniProt"] = final_df.loc[cano_index, ].apply(
            lambda x: x["Entry"] if x["UniProt"] in x["canonical_isoform"] else x["UniProt"], axis=1)
        
        # Add Gene Status
        self.getGeneStatus(final_df)
        # Label Mapping Status
        self.label_mapping_status(final_df)

        pathOb = Path(path)
        edPath = str(Path(pathOb.parent, f'{pathOb.stem}_ed{pathOb.suffix}'))
        final_df.to_csv(edPath, sep=sep)
        return edPath


if __name__ == "__main__":
    dfrm = pd.read_csv(r'C:\OmicData\LiGroupWork\compareTool\1211\exac_missense_mo.txt', sep='\t', nrows=25000)
    demo = MapUniProtID(id_col='ENST', id_type='ENSEMBL_TRS_ID', dfrm=dfrm, site_col='aa_pos', gene_col='GENE')
    demo.set_logging_fileHandler(r'C:\GitWorks\Muta3DMaps\Muta3DMaps\test\data\uniprot\id_map.log')
    demo.retrieve(r'C:\GitWorks\Muta3DMaps\Muta3DMaps\test\data\uniprot\uniprot_id_mapping_test.tsv')
