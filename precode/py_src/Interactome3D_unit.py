# @Date:   2019-08-16T23:34:20+08:00
# @Email:  1730416009@stu.suda.edu.cn
# @Filename: Interactome3D_unit.py
# @Last modified time: 2019-11-07T21:56:39+08:00
import pandas as pd
import wget, time, sys
from urllib import request
from retrying import retry
from multiprocessing.dummy import Pool
from xml.etree import ElementTree
sys.path.append('./')
from Unit import Unit


class Interactome3D_unit(Unit):
    # Reference: https://interactome3d.irbbarcelona.org/help.php#restful
    CONFIG = {
        'DOWNLOAD_FOLDER': '/data/zzf/Interactome3D_files/',
        'INTERACTION_SIFTS_COL': ['interac_TYPE', 'pdb_id', 'interac_BIO_UNIT',
                                  'interac_FILENAME', 'interac_group_compo',
                                  'UniProt', 'chain_id', 'interac_MODEL',
                                  'interac_SEQ_IDENT', 'interac_COVERAGE',
                                  'interac_DOMAIN', 'interac_model_len',
                                  'interac_model_range'],
        'DOWNLOAD_URL': 'https://interactome3d.irbbarcelona.org/api/%s?',  # %s=%s&%s=%s

    }

    def get_interactions_meta(self, species='human', struct_type=False, filePath=False, related_unp=False, related_pdb=False, outputPath=False):
        if not filePath:
            url = "https://interactome3d.irbbarcelona.org/user_data/%s/download/complete/interactions.dat" % species
            filePath = self.CONFIG['DOWNLOAD_FOLDER'] + 'interactions_%s_%s.dat' % (species, time.strftime("%Y_%m_%d", time.localtime()))
            wget.download(url, out=filePath)
            self.file_list.append(filePath)

        dfrm = pd.read_csv(filePath, sep='\t')
        if struct_type:
            dfrm = dfrm[dfrm['TYPE'] == struct_type].reset_index(drop=True)
        dfrm['PDB_ID'] = dfrm.apply(lambda x: x['PDB_ID'].upper(), axis=1)
        dfrm['group_compo'] = dfrm.apply(lambda x: '%s_%s' % tuple(sorted([x['PROT1'], x['PROT2']])), axis=1)
        common_cols = ['TYPE', 'PDB_ID', 'BIO_UNIT', 'FILENAME', 'group_compo']
        s_cols = ['PROT', 'CHAIN', 'MODEL', 'SEQ_IDENT', 'COVERAGE', 'SEQ_BEGIN', 'SEQ_END', 'DOMAIN']
        get_s_cols = lambda num: ['%s%s' % (i, num) for i in s_cols]

        df1, df2 = dfrm[common_cols+get_s_cols(1)].copy(), dfrm[common_cols+get_s_cols(2)].copy()
        df1.columns, df2.columns = common_cols+s_cols, common_cols+s_cols
        df12 = pd.concat([df1, df2]).reset_index(drop=True)
        df12['model_len'] = df12.apply(lambda x: x['SEQ_END'] - x['SEQ_BEGIN'] + 1, axis=1)
        df12['model_range'] = df12.apply(lambda x: '[[%d, %d]]' % (x['SEQ_BEGIN'], x['SEQ_END']), axis=1)

        if related_unp:
            df12 = df12[df12['PROT'].isin(related_unp)]

        if not isinstance(related_pdb, bool):
            df12 = df12[df12['PDB_ID'].isin(related_pdb)]

        self.file_o(outputPath, df12)
        return df12

    def add_SIFTS_to_interactions_meta(self, sifts_df=False, sifts_filePath=False, interac_df=False, interac_filePath=False, outputPath=False):
        sifts_dfrm = self.file_i(sifts_filePath, sifts_df, ('sifts_filePath', 'sifts_df'))
        interac_dfrm = self.file_i(interac_filePath, interac_df, ('interac_filePath', 'interac_df'))
        interac_dfrm = interac_dfrm[interac_dfrm['TYPE'] == 'Structure'].reset_index(drop=True)
        interac_dfrm.drop(columns=['SEQ_BEGIN', 'SEQ_END'], inplace=True)
        interac_dfrm.columns = self.CONFIG['INTERACTION_SIFTS_COL']
        dfrm = pd.merge(interac_dfrm, sifts_dfrm, on=['pdb_id', 'chain_id', 'UniProt'], how='left')
        self.file_o(outputPath, dfrm)
        return dfrm

    @retry(stop_max_attempt_number=3, wait_fixed=1000)
    def download_pdb_from_Interactome3D(filename, type='interaction'):
        url = Interactome3D_unit.CONFIG['DOWNLOAD_URL'] % 'getPdbFile' + 'filename=%s&type=%s' % (filename, type)
        xmlPage = request.urlopen(url).read()
        xmlPage = xmlPage.decode('utf-8')
        node = ElementTree.XML(xmlPage)
        with open(Interactome3D_unit.CONFIG['DOWNLOAD_FOLDER']+filename, 'w') as fw:
            fw.write(node[0][0][1].text)
            time.sleep(2)

    def download_model_script(self, fileName_list, chunksize=100):
        for i in range(0, len(fileName_list), chunksize):
            chunk_li = fileName_list[i:i+chunksize]
            pool = Pool(processes=20)
            pool.map(Interactome3D_unit.download_pdb_from_Interactome3D, chunk_li)


if __name__ == '__main__':
    demo = Interactome3D_unit()
    demo.set_lists(['Test is success.'], [])
    print(demo.pdb_list)
