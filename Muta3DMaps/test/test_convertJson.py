import pytest
from Muta3DMaps.core.pdbe.decode import convertJson2other
import ujson
import tablib
from typing import Union, Optional, Iterator, Iterable, Set, Dict, List, Any
from time import perf_counter
import pandas as pd
import matplotlib.pyplot as plt

infile = r'./data/1a01,2xyn,1miu_summary.json'
outfile = r'./data/1a01,2xyn,1miu_summary.tsv'

def task_convertJson_ujson(infile, outfile):
    t0 = perf_counter()
    ob = tablib.Dataset()
    with open(outfile, 'w') as outOb:
        with open(infile, 'r') as inOb:
            data = ujson.load(inOb)
        for index, pdb in enumerate(data):
            sub = ujson.dumps(data[pdb])
            outOb.write(convertJson2other(sub, ob, 'tsv', index, append_data=[pdb], append_header='pdb_id'))
    return perf_counter() - t0


def task_convertJson_pandas(infile, outfile):
    def func(data):
        for pdb in data:
            dfrm = pd.DataFrame(data[pdb])
            dfrm['pdb_id'] = pdb
            yield dfrm
    
    t0 = perf_counter()
    with open(infile, 'r') as inOb:
        data = ujson.load(inOb)
    dfrm = pd.concat((df for df in func(data)),
              ignore_index=True, sort=False)
    dfrm.to_csv(outfile, sep='\t')
    return perf_counter() - t0


def test_convertJson():
    logPath = r'./data/test_convertJson.log'
    with open(logPath, 'a+') as logFile:
        for _ in range(1):
            ujson_res = task_convertJson_ujson(infile, outfile)
            pandas_res = task_convertJson_pandas(infile, outfile)
            logFile.write(f'ujson, {ujson_res}, pandas, {pandas_res}\n')
            assert ujson_res < pandas_res
    dfrm = pd.read_csv(logPath, sep=',', names=['ujson_tag', 'ujson', 'pandas_tag', 'pandas'])
    plt.plot(list(range(len(dfrm))), dfrm['ujson'], label='ujson')
    plt.plot(list(range(len(dfrm))), dfrm['pandas'], label='pandas')
    plt.legend()
    plt.savefig('./data/test_convertJson.png')

