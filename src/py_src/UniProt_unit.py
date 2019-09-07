# @Date:   2019-08-16T20:26:58+08:00
# @Email:  1730416009@stu.suda.edu.cn
# @Filename: UniProt_unit.py
# @Last modified time: 2019-09-06T23:01:26+08:00
import urllib.parse
import urllib.request
from retrying import retry
import pandas as pd
from random import uniform
from time import sleep
import os, sys
sys.path.append('./')
from Unit import Unit


class UniProt_unit(Unit):
    '''
        params = {
            'from': 'ACC+ID',
            'to': 'ACC',
            'format': 'tab',
            'columns': 'id,length,reviewed,comment(ALTERNATIVE%20PRODUCTS),feature(ALTERNATIVE%20SEQUENCE)'...
            'query': list_str,
        }

       self.params['from'], params['to'] (Particial)
       ACC+ID, ACC, ID, EMBL_ID, P_REFSEQ_AC, PDB_ID, ENSEMBL_TRS_ID, ENSEMBL_ID, ...
       (Reference: https://www.uniprot.org/help/api_idmapping)

       ---

       self.params['columns'] (Particial)
       'comma-separated list of column names'
       (Reference: https://www.uniprot.org/help/api_queries)

       Lists the column names for programmatic (RESTful) access to tab-separated
       or Excel downloads of UniProtKB search results.
       (Reference: https://www.uniprot.org/help/uniprotkb_column_names)

       ---

       self.params['format'] (All)
       html | tab | xls | fasta | gff | txt | xml | rdf | list | rss
       (Reference: https://www.uniprot.org/help/api_queries)

       ---

       self.params['query']
       (Reference: https://www.uniprot.org/help/text-search)
    '''

    URL = 'https://www.uniprot.org/uploadlists/'
    COLUMNS = ['id', 'length', 'reviewed', 'comment(ALTERNATIVE%20PRODUCTS)', 'feature(ALTERNATIVE%20SEQUENCE)', 'genes', 'organism', 'sequence']
    COLUMN_DICT = {'id': 'Entry', 'length': 'Length', 'reviewed': 'Status',
                    'comment(ALTERNATIVE%20PRODUCTS)': 'Alternative products (isoforms)',
                    'feature(ALTERNATIVE%20SEQUENCE)': 'Alternative sequence (isoforms)',
                    'genes': 'Gene names', 'organism': 'Organism', 'sequence': 'Sequence'}
    params = {
        'from': 'ACC+ID',
        'to': 'ACC',
        'format': 'tab',
    }

    @retry(stop_max_attempt_number=3, wait_fixed=1000)
    def go_to_uniprot(url, params, code='utf-8'):
        sleep(uniform(0.5, 5))
        data = urllib.parse.urlencode(params)
        data = data.encode('utf-8')
        req = urllib.request.Request(url, data)
        with urllib.request.urlopen(req) as f:
            response = f.read()
        return response.decode(code)

    def get_info_from_uniprot(self, usecols, outputPath, from_list=False, from_list_file_path=False, sep='\t', chunksize=100, header=None):

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
                new_colNames = [self.CONFIG['COLUMN_DICT'].get(i, i) for i in usecols] + ['yourlist', 'isomap']
                finish = pd.read_csv(outputPath, sep='\t', usecols=['yourlist'], names=new_colNames, skiprows=1, header=None)['yourlist']
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

    def split_df(dfrm, colName, sep):
        df = dfrm.copy()
        return df.drop([colName], axis=1).join(df[colName].str.split(sep, expand=True).stack().reset_index(level=1, drop=True).rename(colName))

    def handle_unp_info(self, unp_df=False, unp_filPath=False, usecols=COLUMNS, outputPath=False):
        unp_df = self.file_i(unp_filPath, unp_df, ('unp_filPath', 'unp_df'))
        # Set the columns
        col_1 = usecols.copy()
        col_1.remove('isomap')
        col_2 = usecols.copy()
        col_2.remove('yourlist')
        # Split the DataFrame
        sub_df_1 = UniProt_unit.split_df(unp_df[col_1], 'yourlist', ',')
        sub_df_2 = UniProt_unit.split_df(unp_df[col_2], 'isomap', ',')
        # Merge Back
        sub_df = pd.merge(sub_df_1, sub_df_2, how='left')
        sub_df['drop'] = sub_df.apply(lambda x: x['yourlist'] != x['isomap'].split(' -> ')[0] if not isinstance(x['isomap'], float) else False, axis=1)
        new_unp_df = sub_df[sub_df['drop'] != True].drop(columns=['drop']).reset_index(drop=True)
        # Output
        self.file_o(outputPath, new_unp_df)
        return new_unp_df


if __name__ == '__main__':
    from_list_file_path = '../../data/demo_files/unp_list.txt'
    outputPath = '../../data/demo_files/info_of_unp.tsv'
    uniprot_demo = UniProt_unit()
    uniprot_demo.get_info_from_uniprot(['id', 'length'], outputPath, from_list_file_path=from_list_file_path)

    df = pd.read_csv(outputPath, sep='\t')
    unp_len_df = df[['Entry', 'Length']]
