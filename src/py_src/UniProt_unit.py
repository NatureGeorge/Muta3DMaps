# @Date:   2019-08-16T20:26:58+08:00
# @Email:  1730416009@stu.suda.edu.cn
# @Filename: UniProt_unit.py
# @Last modified time: 2019-11-14T20:30:58+08:00
# import urllib.parse
# import urllib.request
import requests
import pandas as pd
import numpy as np
from random import uniform
from time import sleep
import os, sys, re, json
from collections import Counter
from multiprocessing import Pool
sys.path.append('./')
from Unit import Unit


class UniProt_unit(Unit):
    '''
    Use UniProt ID Mapping API

    Key parameters in using UniProt ID Mapping API:

    .. code-block:: python
        :linenos:

        params = {
            'from': 'ACC+ID',
            'to': 'ACC',
            'format': 'tab',
            'columns': 'id,length,reviewed,comment(ALTERNATIVE%20PRODUCTS),feature(ALTERNATIVE%20SEQUENCE)'...
            'query': list_str,
        }

    .. csv-table:: Details of ID Abbreviation
        :header: "Abbreviation", "Name", "Direction"

        "ACC+ID", "UniProtKB AC/ID", "from"
        "ACC", "UniProtKB AC", "both"
        "ID", "UniProtKB ID", "both"
        "EMBL_ID", "EMBL/GenBank/DDBJ", "both"
        "REFSEQ_NT_ID", "RefSeq Nucleotide", "both"
        "P_REFSEQ_AC", "RefSeq Protein", "both"
        "PDB_ID", "PDB", "both"
        "ENSEMBL_TRS_ID", "Ensembl Transcript", "both"
        "ENSEMBL_ID", "Ensembl", "both"
        "...", "...", "..."

    Here is the `Reference`_.

    .. _Reference: https://www.uniprot.org/help/api_idmapping

    Following table show the details of rest params:

    .. csv-table:: Details of other parameters
        :header: "param", "Explanation", "Reference"

        "columns", "comma-separated list of column names", "https://www.uniprot.org/help/api_queries"
        "columns", "Lists the column names for programmatic (RESTful) access to tab-separated or Excel downloads of UniProtKB search results", "https://www.uniprot.org/help/uniprotkb_column_names"
        "format", "html | tab | xls | fasta | gff | txt | xml | rdf | list | rss", "https://www.uniprot.org/help/api_queries"
        "query", "query text/id(s)", "https://www.uniprot.org/help/text-search"

    '''

    URL = 'https://www.uniprot.org/uploadlists/'
    USER_AGENT = 'Mozilla/5.0 (Windows NT 6.1; WOW64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/63.0.3239.132 Safari/537.36 QIHU 360SE'
    COLUMNS = ['id', 'length', 'reviewed', 'comment(ALTERNATIVE%20PRODUCTS)', 'feature(ALTERNATIVE%20SEQUENCE)', 'genes', 'organism', 'sequence', 'protein%20names']
    COLUMN_DICT = {
                    'id': 'Entry', 'length': 'Length', 'reviewed': 'Status',
                    'comment(ALTERNATIVE%20PRODUCTS)': 'Alternative products (isoforms)',
                    'feature(ALTERNATIVE%20SEQUENCE)': 'Alternative sequence (isoforms)',
                    'genes': 'Gene names', 'organism': 'Organism', 'sequence': 'Sequence',
                    'protein%20names': 'Protein names'}
    params = {
        'from': 'ACC+ID',
        'to': 'ACC',
        'format': 'tab',
    }

    # @retry(stop_max_attempt_number=3, wait_fixed=1000)
    def go_to_uniprot(url, params, code='utf-8'):
        sleep(uniform(0.99, 5))
        data = urllib.parse.urlencode(params)
        data = data.encode('utf-8')
        req = urllib.request.Request(url, data, headers={'User-Agent': 'User-Agent:Mozilla/5.0'})
        with urllib.request.urlopen(req) as f:
            response = f.read()
        return response.decode(code)
        '''
        session = requests.Session()
        adapter = requests.adapters.HTTPAdapter(max_retries=2)
        session.mount('https://', adapter)
        session.mount('http://', adapter)
        session.headers.update({'User-Agent': UniProt_unit.USER_AGENT})

        with session.get(url, params=params) as r:
            result = r.text
        return result
        # with requests.get(url, params=params, headers={'User-Agent': UniProt_unit.USER_AGENT}, stream=True) as r:
        '''

    def get_info_from_uniprot(self, usecols, outputPath, from_list=None, from_list_file_path=None, sep='\t', chunksize=100, header=None):

        def iter_io(iter_object, params, url, outputPath):
            params['query'] = ','.join(iter_object)  # list_str
            result = UniProt_unit.go_to_uniprot(url, params)
            with open(outputPath, 'a+') as outputFile:
                outputFile.write(result)

        def tidy_result(path, colName='Entry', sep='\t'):
            df = pd.read_csv(path, sep=sep, na_values=[colName])
            df.dropna(subset=[colName], inplace=True)
            df.to_csv(path, sep=sep, index=False)

        if usecols != 'all':
            if not set(self.COLUMNS) >= set(usecols):
                print('get_info_from_uniprot(): please specified usecols with elements in %s' % self.COLUMNS)
                return False
            else:
                self.params['columns'] = ','.join(usecols)
        else:
            usecols = self.COLUMNS
            self.params['columns'] = ','.join(self.COLUMNS)

        if from_list_file_path:
            df = pd.read_csv(from_list_file_path, header=header, chunksize=chunksize, sep=sep)
            for chunk in df:
                iter_io(chunk[0], self.params, self.URL, outputPath)
        else:
            if os.path.exists(outputPath):
                new_colNames = [self.COLUMN_DICT.get(i, i) for i in usecols] + ['yourlist', 'isomap']
                try:
                    finish = pd.read_csv(outputPath, sep='\t', usecols=['yourlist'], names=new_colNames, skiprows=1, header=None)['yourlist']
                except Exception:
                    finish = pd.read_csv(outputPath, sep='\t', names=new_colNames[:-1], skiprows=1, header=None)
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

            if finish_li:
                new_li = list(set(from_list) - set(finish_li))
            else:
                new_li = from_list

            for i in range(0, len(new_li), chunksize):
                iter_io(new_li[i:i+chunksize], self.params, self.URL, outputPath)

        tidy_result(outputPath)
        return True

    def __init__(self, dfrm, id_col, id_type, usecols, reportPath, muta_col=None, muta_type=None, gene_col=None):
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
        self.muta_col = muta_col
        self.muta_type = muta_type
        self.usecols = usecols
        self.gene_col = gene_col
        if muta_col is not None:
            self.muta_li = dfrm.groupby(by=[id_col]).apply(lambda x: [i for i in x[muta_col]])
        self.report = open(reportPath, 'a+')

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
        dfrm['canonical_isoform'] = dfrm.apply(lambda x: ','.join(canonical_pattern.findall(x['Alternative products (isoforms)'])) if not isinstance(x['Alternative products (isoforms)'], float) else np.nan, axis=1)
        special_case = dfrm[dfrm['canonical_isoform'] == ''].index
        if len(special_case) > 0:
            canonical_pattern = re.compile(r'IsoId=([0-9A-Z-,\s]+); Sequence=Displayed')
            special_se = dfrm.loc[special_case].apply(lambda x: ','.join(canonical_pattern.findall(x['Alternative products (isoforms)'])), axis=1)
            dfrm.loc[special_case, 'canonical_isoform'] = special_se
            self.report.write("# Special Cases of Canonical Info\n")
            self.report.write(str(special_se))
        else:
            special_se = pd.Series([])

        return special_se

    def get_raw_ID_Mapping(self, outputPath, chunksize=100):
        """
        Get Raw ID MApping Result
        """
        self.params['from'] = self.id_type
        status = self.get_info_from_uniprot(
            self.usecols,
            outputPath,
            self.dfrm[self.id_col].drop_duplicates(),
            chunksize=chunksize
        )

        if not status:
            return False

        new_colNames = [self.COLUMN_DICT.get(i, i) for i in self.usecols] + ['yourlist', 'isomap']
        dfrm = pd.read_csv(outputPath, sep='\t', names=new_colNames, skiprows=1, header=None)
        # WRITE REPORT
        self.report.write("# RAW ID MAPPING FILE RESULT\n# %s\n" % (outputPath))
        self.report.write(str(dfrm.isnull().sum()))
        self.raw_id_mapping = dfrm
        return True

    def set_raw_id_mapping(self, raw_id_mapping):
        """
        Set the ``pandas.DataFrame`` object that will be deal with in ``handle_ID_Mapping()``
        """
        # assert isinstance(raw_id_mapping, pd.DataFrame), "WARNING, raw_id_mapping should a pandas.DataFrame"
        self.raw_id_mapping = raw_id_mapping

    def handle_ID_Mapping(self):
        """
        Deal with different situations

        1. Add New Columns: ``unp_map_tage`` -> Description for the reliability of mapping, ``canonical_isoform``
        2. Split DataFrame
        3. Classification

        **About unp_map_tage**

        * ``Untrusted & No Isoform``
            * 是指UniProt存在Isoform但是Mapping结果没有明确给出是Map上哪条Isoform,转录本序列与蛋白序列不一致
            * It means that there are isoforms in UniProt, but the mapping result does not clearly indicate which isoform is correspond with the transcript(e.g), and the transcript sequence is inconsistent with the protein sequence.
        * ``Trusted & No Isoform``
            * 是指UniProt不存在Isoform,Mapping结果没问题
        * ``Trusted & Isoform``
            * 是指UniProt存在Isoform,Mapping结果没问题


        .. csv-table:: Split DataFrame Example 1
            :header: "Entry", "Gene names", "Status", "Alternative products (isoforms)", "Organism", "Protein names", "yourlist", "isomap"

            "O94827", "PLEKHG5 KIAA0720", "reviewed", "ALTERNATIVE PRODUCTS:  Event=Alternative splic...", "Homo sapiens (Human)", "Pleckstrin homology domain-containing family G...", "NP_065682,NP_001252521", "NP_001252521 -> O94827-6,NP_065682 -> O94827-5"

        .. csv-table:: Splited DataFrame Example 1
            :header: "Entry", "Gene names", "Status", "Alternative products (isoforms)", "Organism", "Protein names", "yourlist", "UniProt"

            "O94827", "PLEKHG5 KIAA0720", "reviewed", "ALTERNATIVE PRODUCTS:  Event=Alternative splic...", "Homo sapiens (Human)", "Pleckstrin homology domain-containing family G...", "NP_065682", "O94827-5"
            "O94827", "PLEKHG5 KIAA0720", "reviewed", "ALTERNATIVE PRODUCTS:  Event=Alternative splic...", "Homo sapiens (Human)", "Pleckstrin homology domain-containing family G...", "NP_001252521", "O94827-6"

        ...

        For more details, please go to https://github.com/NatureGeorge/SIFTS_Plus_Muta_Maps/blob/master/md/UniProt_ID_Mapping.md
        """

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
        df_wni_split = UniProt_unit.split_df(df_with_no_isomap, 'yourlist', ',')
        df_wni_split.drop(columns=['isomap'], inplace=True)
        df_wni_split['UniProt'] = df_wni_split['Entry']  # [yourlist <-> UniProt]
        df_wni_split['unp_map_tage'] = 'Trusted & No Isoform'
        # Find out special cases 1
        df_wni_split_warn = df_wni_split[~df_wni_split['Alternative products (isoforms)'].isnull()].index
        df_wni_split.loc[df_wni_split_warn, 'unp_map_tage'] = 'Untrusted & No Isoform'
        # 'Entry', 'Gene names', 'Status', 'Alternative products (isoforms)', 'Organism', 'yourlist', 'UniProt'
        # ----------------------------------------------------------------------
        # In Class B
        # ----------------------------------------------------------------------
        wi_yourlist_count = df_with_isomap.apply(lambda x: x['yourlist'].count(','), axis=1)
        wi_isomap_count = df_with_isomap.apply(lambda x: x['isomap'].count(','), axis=1)
        # In subClass 1
        df_wi_eq = df_with_isomap.loc[wi_yourlist_count[wi_yourlist_count == wi_isomap_count].index]
        df_wi_eq_split = UniProt_unit.split_df(df_wi_eq.drop(columns=['yourlist']), 'isomap', ',')
        df_wi_eq_split['yourlist'], df_wi_eq_split['UniProt'] = df_wi_eq_split['isomap'].str.split(' -> ', 1).str
        df_wi_eq_split.drop(columns=['isomap'], inplace=True)  # [yourlist <-> UniProt]
        df_wi_eq_split['unp_map_tage'] = 'Trusted & Isoform'
        # # 'Entry', 'Gene names', 'Status', 'Alternative products (isoforms)', 'Organism', 'yourlist', 'UniProt'

        # In subClass 2
        df_wi_ne = df_with_isomap.loc[wi_yourlist_count[wi_yourlist_count != wi_isomap_count].index]
        if len(df_wi_ne) > 0:
            df_wi_ne_split = UniProt_unit.split_df(df_wi_ne, 'isomap', ',')
            df_wi_ne_split.rename(columns={'yourlist': 'checkinglist'}, inplace=True)
            df_wi_ne_split['yourlist'], df_wi_ne_split['UniProt'] = df_wi_ne_split['isomap'].str.split(' -> ', 1).str
            df_wi_ne_split.drop(columns=['isomap'], inplace=True)
            df_wi_ne_split['unp_map_tage'] = 'Trusted & Isoform & Contain Warnings'
            # 'Entry', 'Gene names', 'Status', 'Alternative products (isoforms)', 'Organism', 'yourlist', 'UniProt', 'checkinglist'
            # Find out special cases 2
            usecols = pd.Index(set(df_wi_ne_split.columns) - {'yourlist', 'UniProt'})
            df_wi_ne_warn = UniProt_unit.split_df(df_wi_ne_split[usecols].drop_duplicates(), 'checkinglist', ',')
            df_wi_ne_warn = df_wi_ne_warn[~df_wi_ne_warn['checkinglist'].isin(df_wi_ne_split['yourlist'])].rename(columns={'checkinglist': 'yourlist'})
            df_wi_ne_warn['UniProt'] = df_wi_ne_warn['Entry']
            df_wi_ne_warn['unp_map_tage'] = 'Untrusted & No Isoform'

            # Update UniProt
            final_df = pd.concat((df_wni_split, df_wi_eq_split, df_wi_ne_split.drop(columns=['checkinglist']), df_wi_ne_warn), sort=False).reset_index(drop=True)
        else:
            final_df = pd.concat((df_wni_split, df_wi_eq_split), sort=False).reset_index(drop=True)

        final_df['UniProt'] = final_df.apply(lambda x: x['Entry'] if x['UniProt'] == x['canonical_isoform'] else x['UniProt'], axis=1)
        if len(canonicalInfo_special_se) > 0:
            canonicalInfo_special_case = final_df[final_df['canonical_isoform'].isin(canonicalInfo_special_se)].index
            final_df.loc[canonicalInfo_special_case, 'UniProt'] = final_df.loc[canonicalInfo_special_case].apply(lambda x: x['Entry'] if x['UniProt'] in x['canonical_isoform'] else x['UniProt'], axis=1)
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
            return None
        if self.id_col != self.gene_col:
            gene_map = self.dfrm[[self.id_col, self.gene_col]].drop_duplicates()
            gene_map.index = gene_map[self.id_col]
            gene_map.drop(columns=[self.id_col], inplace=True)
            handled_df['GENE'] = handled_df.apply(lambda z: gene_map[self.gene_col][z['yourlist']], axis=1)
            handled_df['GENE_status'] = handled_df.apply(lambda x: x['GENE'] == x['Gene names'].split(' ')[0] if not isinstance(x['Gene names'], float) else False, axis=1)
        else:
            handled_df['GENE_status'] = handled_df.apply(lambda x: x['yourlist'] == x['Gene names'].split(' ')[0] if not isinstance(x['Gene names'], float) else False, axis=1)

    def label_mapping_status(self, dfrm, constraint_dict):
        """
        Will Change the dfrm, label Mapping Status

        * Add new column (Mapping_status)

        **About Mapping_status**

        * ``Yes``: 可信的结果，进行后续的PDB Mapping; 通过 ``constraint_dict`` 的限制
        * ``Error``: 一个id对应多个UniProt; 通过 ``constraint_dict`` 的限制
        * ``No``: 不可信的结果; 未通过 ``constraint_dict`` 的限制

        Example of ``constraint_dict``:

        .. code-block:: python
            :linenos:

            constraint_dict = {
                "GENE_status": (False, "ne"),  # 'ne' for !=
                "Status": ("reviewed", "eq"),  # 'eq' for ==
                "unp_map_tage": ("Untrusted & No Isoform", "ne")
            }

        """
        dfrm['Mapping_status'] = 'No'
        # pass_index = dfrm[(dfrm['GENE_status'] != False) & (dfrm['Status'] == 'reviewed') & (dfrm['unp_map_tage'] != 'Untrusted & No Isoform')].index
        pass_df = Unit.ConstraintDict.addConstraintToDf(dfrm, constraint_dict)
        pass_index = pass_df.index
        dfrm.loc[pass_index, 'Mapping_status'] = 'Yes'

        multipleCounter, err_li = Counter(dfrm.loc[pass_index, 'yourlist']), []
        for i, j in multipleCounter.items():
            if j > 1:
                err_li.append(i)

        err_index = pass_df[pass_df['yourlist'].isin(err_li)].index
        dfrm.loc[err_index, 'Mapping_status'] = 'Error'

        # Write Report
        all_id = set(self.dfrm[self.id_col])
        self.report.write("\n# All id: %s\n" % len(all_id))
        unmapped_id = all_id - set(dfrm['yourlist'])
        self.report.write("# Unmapped id: %s\n" % len(unmapped_id))
        for i in sorted(unmapped_id):
            self.report.write("%s\n" % i)
        untrusted_id = set(dfrm[dfrm['Mapping_status'] != 'Yes']['yourlist']) - set(dfrm[dfrm['Mapping_status'] == 'Yes']['yourlist'])
        self.report.write("# Untrusted id: %s\n" % len(untrusted_id))
        for i in sorted(untrusted_id):
            self.report.write("%s\n" % i)
        self.report.write("# Error id: %s\n" % len(dfrm.loc[err_index, 'yourlist'].drop_duplicates()))
        self.report.write(str(dfrm.loc[err_index]))

    def script(self, rawOutputPath, handledOutputPath):
        """
        .. code-block:: python
            :linenos:

            from UniProt_unit import UniProt_unit
            id_col = 'RefSeq_protein'
            id_type = 'P_REFSEQ_AC'
            muta_col = 'mutation_unp'
            gene_col = 'GENE'
            usecols = ['id', 'genes', 'reviewed', 'comment(ALTERNATIVE%20PRODUCTS)', 'organism', 'protein%20names']  # Necessary Columns
            reportPath = '/data/zzf/UniProt_files/id_mapping_files/HGMD_RefSeq_protein_mapping_Report_1021.txt'  # Report of Data Processing
            rawOutputPath = '/data/zzf/UniProt_files/id_mapping_files/HGMD_RefSeq_protein_mapping_1021.tsv' # OutPut File of ID Mapping (RAW)
            handledOutputPath = '/data/zzf/groupWorks/HGMD_RefSeq_protein_mapping_modified_1024.tsv'  # OutPut File of ID Mapping (Final Result)
            # Add constraint/filter to data
            constraint_dict = {
                "GENE_status": (False, "ne"),  # 'ne' for !=
                "Status": ("reviewed", "eq"),  # 'eq' for ==
                "unp_map_tage": ("Untrusted & No Isoform", "ne")
            }

            # Initial
            unp_demo = UniProt_unit(group_df, id_col, id_type, usecols, reportPath, muta_col=muta_col, gene_col=gene_col)
            # Return True if get RAW Result Successfully
            unp_demo.get_raw_ID_Mapping(rawOutputPath)
            # Deal with different situations
            handled_df = unp_demo.handle_ID_Mapping()
            # Add Gene Status
            unp_demo.getGeneStatus(handled_df)
            # Label Mapping Status
            unp_demo.label_mapping_status(handled_df, constraint_dict)
            # close the file-handle of report
            unp_demo.report.close()
            # Output the final result
            handled_df.to_csv(handledOutputPath, sep='\t', index=False)

        """
        step1 = self.get_raw_ID_Mapping(rawOutputPath)
        if not step1:
            return None
        handled_df = self.handle_ID_Mapping()
        self.getGeneStatus(handled_df)
        constraint_dict = {
            "GENE_status": (False, "ne"),
            "Status": ("reviewed", "eq"),
            "unp_map_tage": ("Untrusted & No Isoform", "ne")
        }
        self.report.write("Constraint Dict\n"+json.dumps(constraint_dict)+"\n")
        self.label_mapping_status(handled_df, constraint_dict)
        handled_df.to_csv(handledOutputPath, sep='\t', index=False)
        self.report.close()
