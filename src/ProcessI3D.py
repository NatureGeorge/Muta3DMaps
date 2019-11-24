# @Date:   2019-08-16T23:34:20+08:00
# @Email:  1730416009@stu.suda.edu.cn
# @Filename: ProcessI3D.py
# @Last modified time: 2019-11-24T15:06:15+08:00
import pandas as pd
import wget, time, os
from urllib import request
from retrying import retry
from multiprocessing.dummy import Pool
from xml.etree import ElementTree
from Utils.Logger import RunningLogger
from Utils.FileIO import file_o


class RetrieveI3D:
    # Reference: https://interactome3d.irbbarcelona.org/help.php#restful
    CONFIG = {
        'INTERACTION_SIFTS_COL': ['i3d_TYPE', 'pdb_id', 'i3d_BIO_UNIT',
                                  'i3d_FILENAME', 'i3d_INTERACT_COMPO',
                                  'Entry', 'chain_id', 'i3d_MODEL',
                                  'i3d_SEQ_IDENT', 'i3d_COVERAGE',
                                  'i3d_DOMAIN',
                                  'i3d_SAME_MODEL', 'i3d_CHAIN_COMPO',
                                  'i3d_pdb_type', 'i3d_INTERACT_COMPO'
                                  'i3d_model_len', 'i3d_model_range'],
        'DOWNLOAD_URL': 'https://interactome3d.irbbarcelona.org/api/%s?',  # %s=%s&%s=%s
    }

    def __init__(self, **kwargs):
        self.downloadFolder = kwargs["downloadFolder"]
        self.Logger = RunningLogger("RetrieveI3D", kwargs["loggingPath"])

    def get_interactions_meta(self, species='human', struct_type=False, filePath=False, related_unp=False, related_pdb=False, outputPath=False):
        if not filePath:
            url = "https://interactome3d.irbbarcelona.org/user_data/%s/download/complete/interactions.dat" % species
            filePath = os.path.join(self.downloadFolder, 'I3D_META_interactions_%s_%s.dat' % (species, time.strftime("%Y_%m_%d", time.localtime())))
            self.Logger.logger.info("Downloading File: %s" % filePath)
            wget.download(url, out=filePath)

        dfrm = pd.read_csv(filePath, sep='\t', na_values=["-"])
        if struct_type:
            dfrm = dfrm[dfrm['TYPE'] == struct_type].reset_index(drop=True)

        dfrm["i3d_SAME_MODEL"] = dfrm.apply(lambda x: x["MODEL1"] == x["MODEL2"], axis=1)
        dfrm['i3d_CHAIN_COMPO'] = dfrm.apply(lambda x: '%s_%s' % (x['CHAIN1'], x['CHAIN2']), axis=1)
        dfrm["i3d_pdb_type"] = dfrm.apply(lambda x: "ho" if x['PROT1'] == x['PROT2'] else "he", axis=1)
        dfrm['i3d_INTERACT_COMPO'] = dfrm.apply(lambda x: '%s_%s' % tuple(sorted([x['PROT1'], x['PROT2']])), axis=1)  # Need to check

        dfrm['PDB_ID'] = dfrm.apply(lambda x: x['PDB_ID'].upper(), axis=1)

        common_cols = ['TYPE', 'PDB_ID', 'BIO_UNIT', 'FILENAME', 'i3d_SAME_MODEL', 'i3d_CHAIN_COMPO', 'i3d_pdb_type', 'i3d_INTERACT_COMPO']
        s_cols = ['PROT', 'CHAIN', 'MODEL', 'SEQ_IDENT', 'COVERAGE', 'SEQ_BEGIN', 'SEQ_END', 'DOMAIN']

        get_s_cols = lambda num: ['%s%s' % (i, num) for i in s_cols]

        df1, df2 = dfrm[common_cols+get_s_cols(1)].copy(), dfrm[common_cols+get_s_cols(2)].copy()
        df1.columns, df2.columns = common_cols+s_cols, common_cols+s_cols
        df12 = pd.concat([df1, df2]).reset_index(drop=True)

        df12['i3d_model_len'] = df12.apply(lambda x: x['SEQ_END'] - x['SEQ_BEGIN'] + 1, axis=1)
        df12['i3d_model_range'] = df12.apply(lambda x: '[[%d, %d]]' % (x['SEQ_BEGIN'], x['SEQ_END']), axis=1)

        if related_unp:
            df12 = df12[df12['PROT'].isin(related_unp)]

        if not isinstance(related_pdb, bool):
            df12 = df12[df12['PDB_ID'].isin(related_pdb)]

        df12.drop(columns=['SEQ_BEGIN', 'SEQ_END'], inplace=True)
        df12.columns = self.CONFIG['INTERACTION_SIFTS_COL']

        file_o(outputPath, df12)
        return df12

    @retry(stop_max_attempt_number=3, wait_fixed=1000)
    def download_pdb_from_Interactome3D(self, filename, type='interaction'):
        url = RetrieveI3D.CONFIG['DOWNLOAD_URL'] % 'getPdbFile' + 'filename=%s&type=%s' % (filename, type)
        xmlPage = request.urlopen(url).read()
        xmlPage = xmlPage.decode('utf-8')
        node = ElementTree.XML(xmlPage)
        with open(os.path.join(self.downloadFolder, filename), 'w') as fw:
            fw.write(node[0][0][1].text)
            time.sleep(2)

    def download_model_script(self, fileName_list, chunksize=100):
        for i in range(0, len(fileName_list), chunksize):
            chunk_li = fileName_list[i:i+chunksize]
            pool = Pool(processes=20)
            pool.map(self.download_pdb_from_Interactome3D, chunk_li)
