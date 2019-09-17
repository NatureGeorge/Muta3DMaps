# @Date:   2019-09-17T18:14:49+08:00
# @Email:  1730416009@stu.suda.edu.cn
# @Filename: SMR_unit.py
# @Last modified time: 2019-09-17T20:45:17+08:00
import sys
import pandas as pd
import numpy as np
import json
import os
import wget
import tarfile
import time
from collections import Iterable, Iterator
sys.path.append('./')
from Unit import Unit


class SMR_unit(Unit):

    CONFIG = {
            'DOWNLOAD_FOLDER': '/data/zzf/SMR_files/',
            'SPECIES_DICT': {'Human': '9606', 'Mouse': '10090'},
            'SMR_INDEX_FILE_PATH': '/data/zzf/SMR_files/SWISS-MODEL_Repository/INDEX',
        }

    def get_SMR_meta(self, species_code='9606', provider='SWISSMODEL', filePath=False, related_unp=False, outputPath=False):
        if not filePath:
            url = "https://swissmodel.expasy.org/repository/download/core_species/%s_meta.tar.gz" % species_code
            filePath = self.CONFIG['DOWNLOAD_FOLDER'] + 'SMR_Meta_%s_%s.tar.gz' % (species_code, time.strftime("%Y_%m_%d", time.localtime()))
            wget.download(url, out=filePath)
            self.file_list.append(filePath)
            tar = tarfile.open(filePath)
            for name in tar:
                tar.extract(name, self.CONFIG['DOWNLOAD_FOLDER'])
                print(name)
                self.file_list.append(name)
            dfrm = pd.read_csv(self.CONFIG['DOWNLOAD_FOLDER']+'SWISS-MODEL_Repository/INDEX', sep='\t', skiprows=6)
        else:
            dfrm = pd.read_csv(filePath, sep='\t', skiprows=6)

        if provider:
            dfrm = dfrm[dfrm['provider'] == provider].reset_index(drop=True)

        dfrm.rename(columns={'UniProtKB_ac': 'UniProt', 'seqid':'smr_meta_seqid','url':'smr_url', 'template':'smr_template'}, inplace=True)

        if isinstance(related_unp, (Iterable, Iterator)):
            dfrm = dfrm[dfrm['UniProt'].isin(related_unp)]

        dfrm['smr_meta_range'] = dfrm.apply(lambda x: '[[%d, %d]]' % (x['from'], x['to']), axis=1)

        self.file_o(outputPath, dfrm)
        return dfrm

    def webDownload(out_fname, url):
        if os.path.exists(out_fname):
            return
        try:
            wget.download(url, out=out_fname)
        except Exception as e:
            if hasattr(e, 'code'):
                print(e.code, 1)
            if hasattr(e, 'reason'):
                print(e.reason, 2)
            out_fname = 'null'
        print('webDownload():', out_fname)
