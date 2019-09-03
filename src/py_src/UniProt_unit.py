# @Date:   2019-08-16T20:26:58+08:00
# @Email:  1730416009@stu.suda.edu.cn
# @Filename: UniProt_unit.py
# @Last modified time: 2019-09-03T21:36:02+08:00
import urllib.parse
import urllib.request
import pandas as pd
from random import uniform
from time import sleep

class UniProt_unit:
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
    COLUMNS = ['id','length','reviewed','comment(ALTERNATIVE%20PRODUCTS)','feature(ALTERNATIVE%20SEQUENCE)', 'genes', 'organism', 'sequence']
    params = {
        'from': 'ACC+ID',
        'to': 'ACC',
        'format': 'tab',
    }

    def go_to_uniprot(url, params, code='utf-8'):
        sleep(uniform(0.5,5))
        data = urllib.parse.urlencode(params)
        data = data.encode('utf-8')
        req = urllib.request.Request(url, data)
        with urllib.request.urlopen(req) as f:
           response = f.read()
        return response.decode(code)

    def get_info_from_uniprot(self, usecols, outputPath, from_list=False, from_list_file_path=False, sep='\t', chunksize=100, header=None):

        def iter_io(iter_object, params, url, outputPath):
            params['query'] = ','.join(iter_object) # list_str
            result = UniProt_unit.go_to_uniprot(url, params)
            with open(outputPath, 'a+') as outputFile:
                outputFile.write(result)

        def tidy_result(path, colName='Entry', sep='\t'):
            df = pd.read_csv(path, sep=sep, na_values=[colName])
            df.dropna(subset=[colName], inplace=True)
            df.to_csv(path, sep=sep, index=False)

        if usecols != 'all':
            if not set(self.COLUMNS) >= set(usecols):
                print('get_info_from_uniprot(): please specified usecols with elements in %s'%self.COLUMNS)
                return False
            else:
                self.params['columns'] = ','.join(usecols)
        else:
            self.params['columns'] = ','.join(self.COLUMNS)

        if from_list_file_path:
            df = pd.read_csv(from_list_file_path, header=header, chunksize=chunksize, sep=sep)
            for chunk in df:
                iter_io(chunk[0], self.params, self.URL, outputPath)
        else:
            for i in range(0, len(from_list), chunksize):
                iter_io(from_list[i:i+chunksize], self.params, self.URL, outputPath)

        tidy_result(outputPath)
        return True


if __name__ == '__main__':
    from_list_file_path = '../../data/demo_files/unp_list.txt'
    outputPath = '../../data/demo_files/info_of_unp.tsv'
    uniprot_demo = UniProt_unit()
    uniprot_demo.get_info_from_uniprot(['id','length'], outputPath, from_list_file_path=unp_list_file_path)

    df = pd.read_csv(outputPath, sep='\t')
    unp_len_df = df[['Entry', 'Length']]
