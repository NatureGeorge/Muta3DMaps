# @Date:   2019-08-16T23:34:20+08:00
# @Email:  1730416009@stu.suda.edu.cn
# @Filename: Interactome3D_unit.py
# @Last modified time: 2019-08-16T23:35:51+08:00
import pandas as pd
import numpy as np
import wget, time, sys
sys.path.append('./')
from Unit import Unit


class Interactome3D(Unit):
    CONFIG = {
        'DOWNLOAD_FOLDER': '/home/zzf/Work/StructDDG_0427/data/Mapping_Pipeline/Interactome3D_files/',
        'INTERACTION_SIFTS_COL': ['interac_TYPE', 'pdb_id', 'interac_BIO_UNIT',
                                  'interac_FILENAME', 'interac_group_compo',
                                  'UniProt', 'chain_id', 'interac_MODEL',
                                  'interac_SEQ_IDENT', 'interac_COVERAGE',
                                  'interac_DOMAIN', 'interac_model_len',
                                  'interac_model_range'],

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
        common_cols = ['TYPE','PDB_ID','BIO_UNIT','FILENAME', 'group_compo']
        s_cols = ['PROT', 'CHAIN', 'MODEL', 'SEQ_IDENT', 'COVERAGE', 'SEQ_BEGIN', 'SEQ_END', 'DOMAIN']
        get_s_cols = lambda num: ['%s%s' % (i, num) for i in s_cols]

        df1, df2 = dfrm[common_cols+get_s_cols(1)].copy(), dfrm[common_cols+get_s_cols(2)].copy()
        df1.columns, df2.columns = common_cols+s_cols, common_cols+s_cols
        df12 = pd.concat([df1,df2]).reset_index(drop=True)
        df12['model_len'] = df12.apply(lambda x: x['SEQ_END'] - x['SEQ_BEGIN'] + 1, axis=1)
        df12['model_range'] = df12.apply(lambda x: '[[%d, %d]]'%(x['SEQ_BEGIN'],x['SEQ_END']), axis=1)

        if related_unp:
            df12 = df12[df12['PROT'].isin(related_unp)]

        if related_pdb:
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


if __name__ == '__main__':
    demo = Interactome3D()
    demo.set_lists(['Test is success.'], [])
    print(demo.pdb_list)
