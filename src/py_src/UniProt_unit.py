# @Date:   2019-08-16T20:26:58+08:00
# @Email:  1730416009@stu.suda.edu.cn
# @Filename: get_unp_len_from_uniprot.py
# @Last modified time: 2019-08-16T23:07:24+08:00
import urllib.parse
import urllib.request
import pandas as pd

class UniProt_unit:
    '''
        params = {
            'from': 'ACC+ID',
            'to': 'ACC',
            'format': 'tab',
            'columns': 'id,length,reviewed,comment(ALTERNATIVE%20PRODUCTS),feature(ALTERNATIVE%20SEQUENCE)', # (Particial)
            'query': pdb_list_str,
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
    params = {
        'from': 'ACC+ID',
        'to': 'ACC',
        'format': 'tab',
        'columns': 'id,length,reviewed,comment(ALTERNATIVE%20PRODUCTS),feature(ALTERNATIVE%20SEQUENCE)',
    }

    def go_to_uniprot(url, params, code='utf-8'):
        data = urllib.parse.urlencode(params)
        data = data.encode('utf-8')
        req = urllib.request.Request(url, data)
        with urllib.request.urlopen(req) as f:
           response = f.read()
        return response.decode(code)

    def get_info_from_uniprot(self, outputPath, unp_list=False, unp_list_file_path=False, sep='\t', chunksize=100, header=None, unp_col=0):

        def iter_io(iter_object, params, url, outputPath):
            params['query'] = ','.join(iter_object) # pdb_list_str
            result = UniProt_unit.go_to_uniprot(url, params)
            with open(outputPath, 'a+') as outputFile:
                outputFile.write(result)

        if unp_list_file_path:
            df = pd.read_csv(unp_list_file_path, header=header, chunksize=chunksize, sep=sep)
            for chunk in df:
                iter_io(chunk[0], self.params, self.URL, outputPath)
        else:
            for i in range(0, len(unp_list), chunksize):
                iter_io(unp_list[i:i+chunksize], self.params, self.URL, outputPath)


if __name__ == '__main__':
    # print(help(UniProt_unit))
    # '''
    unp_list_file_path = 'C:/Users/Nature/Desktop/LiGroup/Work/StructDDG_0427/data/0816_unp_list.txt'
    outputPath = '../../data/info_of_UniprotID_list.tab'
    uniprot_demo = UniProt_unit()
    uniprot_demo.get_info_from_uniprot(outputPath, unp_list_file_path=unp_list_file_path)
    # '''
