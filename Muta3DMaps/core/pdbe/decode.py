# @Created Date: 2020-01-12 01:27:18 pm
# @Filename: decode.py
# @Email:  1730416009@stu.suda.edu.cn
# @Author: ZeFeng Zhu
# @Last Modified: 2020-02-11 04:22:22 pm
# @Copyright (c) 2020 MinghuiGroup, Soochow University
import os
import pandas as pd
import pyexcel as pe
import tablib
from tablib import InvalidDimensions, UnsupportedFormat
from typing import Union, Optional, Iterator, Iterable, Set, Dict, List, Any, Generator, Callable
from json import JSONDecodeError
import ujson as json
import time
from pathlib import Path
from logging import Logger
from collections import OrderedDict
from Muta3DMaps.core.utils import decompression
from Muta3DMaps.core.log import Abclog
from Muta3DMaps.core.retrieve.fetchFiles import UnsyncFetch

API_LYST: List = sorted(['summary', 'molecules', 'experiment', 'ligand_monomers',
                   'modified_AA_or_NA', 'mutated_AA_or_NA', 'status',
                   'polymer_coverage', 'secondary_structure',
                   'residue_listing', 'binding_sites', 'files', 'observed_residues_ratio',
                   'assembly', 'electron_density_statistics',
                   'cofactor', 'drugbank', 'related_experiment_data'])

BASE_URL: str = 'http://www.ebi.ac.uk/pdbe/api'

FTP_URL: str = 'ftp://ftp.ebi.ac.uk'

FTP_DEFAULT_PATH: str = '/pub/databases/msd/sifts/flatfiles/tsv/uniprot_pdb.tsv.gz'

FUNCS = list()

def dispatch_on_set(keys: Set):
    '''
    Decorator to add new dispatch functions
    '''
    def register(func):
        FUNCS.append((func, set(keys)))
        return func
    return register


def traversePDBeData(query: Any, *args):
    for func, keySet in FUNCS:
        if query in keySet:
            return func(*args)
    else:
        raise ValueError('Invalid query')


def convertJson2other(
        data: Union[List, str, None], 
        append_data: Union[Iterable, Iterator],
        converter: Optional[tablib.Dataset] = None, 
        export_format: str = 'tsv', 
        ignore_headers: Union[bool, int] = False, 
        log_func=print) -> Any:
    '''
    Convert valid json-string/dict into specified format via `tablib.Dataset.export`
    '''
    if converter is None:
        converter = tablib.Dataset()
    try:
        if isinstance(data, str):
            converter.json = data
        elif isinstance(data, List):
            converter.dict = data
        elif data is None:
            pass
        else:
            log_func(f'Invalid type for data`: {type(data)}')
            return None
        for data_to_append in append_data:
            converter.append_col(*data_to_append)
        if ignore_headers:
            converter.headers = None
        return converter.export(export_format)
    except KeyError:
        log_func('Not a valid json-string/dict to convert format via `tablib.Dataset.export`')
    except JSONDecodeError:
        log_func(f'Invalid json string')
    except InvalidDimensions:
        log_func('Invalid data or append_data')
    except UnsupportedFormat:
        log_func(f'Invalid export_format: {export_format}')

class ProcessSIFTS(Abclog):
    @classmethod
    def related_UNP_PDB(cls, filePath: Union[str, Path], related_unp: Optional[Iterable] = None, related_pdb: Optional[Iterable] = None):
        '''
        Reference
        
            * http://www.ebi.ac.uk/pdbe/docs/sifts/quick.html
            * A summary of the UniProt to PDB mappings showing the UniProt accession
            followed by a semicolon-separated list of PDB four letter codes.
        '''
        filePath = Path(filePath)
        if filePath.is_dir():
            url = FTP_URL+FTP_DEFAULT_PATH
            task = ('ftp', {'url': url}, str(filePath))
            res = UnsyncFetch.multi_tasks([task]).result()
            filePath = decompression(res[0], remove=False, logger=cls.logger)
        elif filePath.is_file():
            filePath = str(filePath)
        else:
            raise ValueError('Valid value for filePath')

        dfrm = pd.read_csv(filePath, sep='\t', header=1)
        pdb_list = list()
        if related_unp is not None:
            dfrm = dfrm[dfrm['SP_PRIMARY'].isin(related_unp)]
        for i in dfrm.index:
            pdb_list.extend(dfrm.loc[i, 'PDB'].split(';'))
        if related_pdb is not None:
            return set(pdb_list) & set(related_pdb), set(dfrm['SP_PRIMARY'])
        else:
            return set(pdb_list), set(dfrm['SP_PRIMARY'])

    @staticmethod
    def yieldTasks(pdbs: Union[Iterable, Iterator], folder: Union[str, Path]) -> Generator:
        for pdb in pdbs:
            pdb = pdb.lower()
            yield 'get', {'url': f'{BASE_URL}/mappings/all_isoforms/{pdb}'}, str(Path(folder, f'{pdb}_sifts.json'))

    @classmethod
    def process(cls):
        pass

    @classmethod
    def retrieve(cls, pdbs: Union[Iterable, Iterator], outputFolder: Union[str, Path], concur_req: int = 20, rate: float = 1.5):
        outputFolder = Path(outputFolder)
        if outputFolder.is_file():
            folder = outputFolder.stem
        elif outputFolder.is_dir():
            folder = outputFolder
        else:
            raise ValueError('Invalid value for outputFolder')

        t0 = time.perf_counter()
        res = UnsyncFetch.multi_tasks(cls.yieldTasks(
            pdbs, folder), cls.process, concur_req, rate, cls.logger).result()
        elapsed = time.perf_counter() - t0
        cls.logger.info(f'{len(pdbs)} pdbs downloaded in {elapsed}s')
        return res


class ProcessPDBe(Abclog):
    pass


class PDBeDecoder(object):
    @staticmethod
    def sync_with_pyexcel(*args) -> pe.Sheet:
        records, *remain = args
        sheet = pe.get_sheet(records=records, name_columns_by_row=0)
        if len(remain) > 1:
            append_header, append_value = remain
            append_data = [append_header] + [append_value]*len(records)
            sheet.column += pe.Sheet(append_data)
        return sheet

    @classmethod
    def pyexcel_io(cls, suffix: str, data: Dict, **kwargs) -> pe.Sheet:
        cur_sheet = None
        for res in traversePDBeData(suffix, data):
            try:
                cur_sheet.row += cls.sync_with_pyexcel(*res)
            except AttributeError:
                cur_sheet = cls.sync_with_pyexcel(*res)
        if kwargs:
            cur_sheet.save_as(**kwargs)
        return cur_sheet

    @staticmethod
    def sync_with_tablib(*args) -> tablib.Dataset:
        records, *remain = args
        ob = tablib.Dataset()
        for index in range(len(records)):
            records[index] = OrderedDict(sorted(records[index].items()))
        ob.dict = records
        if len(remain) > 1:
            append_header, append_value = remain
            for i in range(len(append_value)):
                ob.append_col([append_value[i]]*len(records), append_header[i])
        return ob

    @classmethod
    def tablib_io(cls, suffix: str, data: Dict, **kwargs) -> tablib.Dataset:
        cur_ob = None
        for res in traversePDBeData(suffix, data):
            try:
                cur_ob.dict += cls.sync_with_tablib(*res).dict
            except AttributeError:
                cur_ob = cls.sync_with_tablib(*res)
        if kwargs:
            with open(file=kwargs['file'], mode=kwargs.get('mode', 'w+')) as outputFile:
                outputFile.write(cur_ob.export(kwargs['format']))
        return cur_ob

    @staticmethod
    @dispatch_on_set({'/pdb/entry/status/', '/pdb/entry/summary/', '/pdb/entry/modified_AA_or_NA/',
                      '/pdb/entry/mutated_AA_or_NA/', '/pdb/entry/cofactor/', '/pdb/entry/molecules/',
                      '/pdb/entry/ligand_monomers/', '/pdb/entry/experiment/',
                      '/pdb/entry/electron_density_statistics/',
                      '/pdb/entry/related_experiment_data/', '/pdb/entry/drugbank/'})
    def yieldCommon(data: Dict) -> Generator:
        for pdb in data:
            values = data[pdb]
            for value in values:
                for key in value:
                    if isinstance(value[key], (Dict, List)):
                        value[key] = json.dumps(value[key])
            yield values, ('pdb_id',), (pdb,)

    @staticmethod
    @dispatch_on_set({'/pdb/entry/polymer_coverage/'})
    def yieldPolymerCoverage(data: Dict) -> Generator:
        for pdb in data:
            molecules = data[pdb]['molecules']
            for entity in molecules:
                chains = entity['chains']
                for chain in chains:
                    observed = chain['observed']
                    for fragement in observed:
                        for key in ('start', 'end'):
                            fragement[key] = json.dumps(fragement[key])
                    yield observed, ('chain_id', 'struct_asym_id', 'entity_id', 'pdb_id'), (chain['chain_id'], chain['struct_asym_id'], entity['entity_id'], pdb)

    @staticmethod
    @dispatch_on_set({'/pdb/entry/observed_residues_ratio/'})
    def yieldObservedResiduesRatio(data: Dict) -> Generator:
        for pdb in data:
            for entity_id, entity in data[pdb].items():
                yield entity, ('entity_id', 'pdb_id'), (entity_id, pdb)

    @staticmethod
    @dispatch_on_set({'/pdb/entry/residue_listing/'})
    def yieldResidues(data: Dict) -> Generator:
        for pdb in data:
            molecules = data[pdb]['molecules']
            for entity in molecules:
                chains = entity['chains']
                for chain in chains:
                    residues = chain['residues']
                    for res in residues:
                        if 'multiple_conformers' not in res:
                            res['multiple_conformers'] = None
                        else:
                            res['multiple_conformers'] = json.dumps(res['multiple_conformers'])
                    yield residues, ('chain_id', 'struct_asym_id', 'entity_id', 'pdb_id'), (chain['chain_id'], chain['struct_asym_id'], entity['entity_id'], pdb)

    @staticmethod
    @dispatch_on_set({'/pdb/entry/secondary_structure/'})
    def yieldSecondaryStructure(data: Dict):
        for pdb in data:
            molecules = data[pdb]['molecules']
            for entity in molecules:
                chains = entity['chains']
                for chain in chains:
                    secondary_structure = chain['secondary_structure']
                    for name in secondary_structure:
                        fragment = secondary_structure[name]
                        for record in fragment:
                            for key in record:
                                if isinstance(record[key], (Dict, List)):
                                    record[key] = json.dumps(record[key])
                            if 'sheet_id' not in record:
                                record['sheet_id'] = None
                        yield fragment, ('secondary_structure', 'chain_id', 'struct_asym_id', 'entity_id', 'pdb_id'), (name, chain['chain_id'], chain['struct_asym_id'], entity['entity_id'], pdb)

    @staticmethod
    @dispatch_on_set({'/pdb/entry/binding_sites/'})
    def yieldBindingSites(data: Dict) -> Generator:
        for pdb in data:
            for site in data[pdb]:
                for tage in ('site_residues', 'ligand_residues'):
                    residues = site[tage]
                    for res in residues:
                        if 'symmetry_symbol' not in res:
                            res['symmetry_symbol'] = None
                    yield residues, ('residues_type', 'details', 'evidence_code', 'site_id', 'pdb_id'), (tage, site['details'], site['evidence_code'], site['site_id'], pdb)

    @staticmethod
    @dispatch_on_set({'/pdb/entry/assembly/'})
    def yieldAssembly(data: Dict):
        for pdb in data:
            for biounit in data[pdb]:
                entities = biounit['entities']
                for entity in entities:
                    for key in entity:
                        if isinstance(entity[key], (Dict, List)):
                            entity[key] = json.dumps(entity[key])
                keys = list(biounit)
                keys.remove('entities')
                yield entities, tuple(keys)+('pdb_id',), tuple(biounit[key] for key in keys)+(pdb, )

    @staticmethod
    @dispatch_on_set({'/pdb/entry/files/'})
    def yieldAssociatedFiles(data: Dict):
        for pdb in data:
            for key in data[pdb]:
                for innerKey in data[pdb][key]:
                    record = data[pdb][key][innerKey]
                    if record:
                        yield record, ('innerKey', 'key', 'pdb_id'), (innerKey, key, pdb)
                    else:
                        continue

    @staticmethod
    @dispatch_on_set({'/mappings/all_isoforms/'})
    def yieldSIFTSRange(data: Dict) -> Generator:
        top_root = next(iter(data))  # PDB_ID or UniProt Isoform ID
        sec_root = next(iter(data[top_root]))  # 'UniProt' or 'PDB'
        child = data[top_root][sec_root]
        thi_root = next(iter(child))
        test_value = child[thi_root]
        # from PDB to UniProt
        if isinstance(test_value, Dict) and sec_root == 'UniProt':
            for uniprot in child:
                name = child[uniprot]['name']
                identifier = child[uniprot]['identifier']
                chains = child[uniprot]['mappings']
                for chain in chains:
                    chain['start'] = json.dumps(chain['start'])
                    chain['end'] = json.dumps(chain['end'])
                    chain['pdb_id'] = top_root
                    chain[sec_root] = uniprot
                    chain['identifier'] = identifier
                    chain['name'] = name
                yield chains, None
        # from UniProt to PDB
        elif isinstance(test_value, List) and sec_root == 'PDB':
            for pdb in child:
                chains = child[pdb]
                for chain in chains:
                    chain['start'] = json.dumps(chain['start'])
                    chain['end'] = json.dumps(chain['end'])
                yield chains, ('pdb_id', 'UniProt'), (pdb, top_root)
        else:
            raise ValueError(f'Unexpected data structure for inputted data: {data}')


# TODO: Chain UniProt ID Mapping -> ProcessSIFTS -> ProcessPDBe
# TODO: Deal with oligomeric PDB
