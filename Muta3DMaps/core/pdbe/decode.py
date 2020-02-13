# @Created Date: 2020-01-12 01:27:18 pm
# @Filename: decode.py
# @Email:  1730416009@stu.suda.edu.cn
# @Author: ZeFeng Zhu
# @Last Modified: 2020-02-11 04:22:22 pm
# @Copyright (c) 2020 MinghuiGroup, Soochow University
import os
import pandas as pd
import tablib
from tablib import InvalidDimensions, UnsupportedFormat
from typing import Union, Optional, Iterator, Iterable, Set, Dict, List, Any, Generator
from json import JSONDecodeError


try:
    import ujson as json
except ImportError:
    import json

API_LYST = sorted(["summary", "molecules", "experiment", "ligand_monomers",
                   "modified_AA_or_NA", "mutated_AA_or_NA", "status",
                   "polymer_coverage", "secondary_structure",
                   "residue_listing", "binding_sites", "files", "observed_residues_ratio",
                   "assembly", "electron_density_statistics",
                   "cofactor", "drugbank", "related_experiment_data"])


FUNCS = list()


def getFiles(workdir: str, suffix: str):
    for file in os.listdir(workdir):
        if suffix in file:
            yield os.path.join(workdir, file)


def dispatch_on_set(keys: Set):
    """
    Decorator to add new dispatch functions
    """
    def register(func):
        FUNCS.append((func, set(keys)))
        return func
    return register


def traversePDBeData(query: Any, *args):
    for func, keySet in FUNCS:
        if query in keySet:
            return func(*args)


class PDBeJsonDecoder(object):
    @staticmethod
    def yieldJsonDataFromFiles(files: Union[Iterable[str], Iterator[str]]):
        for file in files:
            with open(file, 'rt') as handle:
                data = json.load(handle)
            yield data

    @classmethod
    @dispatch_on_set({"status", "summary", "modified_AA_or_NA",
                      "mutated_AA_or_NA", "cofactor", "molecules",
                      "ligand_monomers", "experiment",
                      "electron_density_statistics",
                      "related_experiment_data", "drugbank"})
    def pdb_common(cls, files: Union[Iterable[str], Iterator[str]], export_format: str = 'tsv') -> Generator:
        jsonDataGenerator = cls.yieldJsonDataFromFiles(files)
        ob = tablib.Dataset()
        for index, data in enumerate(jsonDataGenerator):
            for pdb in data:
                yield convertJson2other(data[pdb], ([pdb], 'pdb_id'), ob, export_format, index)

    # TODO: integrate 'residue_listing' and 'secondary_structure'
    @classmethod
    @dispatch_on_set({'polymer_coverage'})
    def pdb_polymerCoverage(cls, files: Union[Iterable[str], Iterator[str]], export_format: str = 'tsv') -> Generator:
        jsonDataGenerator = cls.yieldJsonDataFromFiles(files)
        ob = tablib.Dataset()
        count = -1
        for data in jsonDataGenerator:
            for pdb in data:
                molecules = data[pdb]['molecules']
                for entity in molecules:
                    chains = entity['chains']
                    for chain in chains:
                        observed = chain['observed']
                        for fragement in observed:
                            for key in ('start', 'end'):
                                fragement[key] = json.dumps(fragement[key])
                        observed_len = len(observed)
                        append_data = zip(
                            (
                                [chain['chain_id']] * observed_len,
                                [chain['struct_asym_id']] * observed_len, 
                                [entity['entity_id']] * observed_len, 
                                [pdb] * observed_len),
                            (
                                'chain_id', 
                                'struct_asym_id', 
                                'entity_id', 
                                'pdb'))
                        count += 1
                        yield convertJson2other(observed, append_data, ob, export_format, count)

    @classmethod
    @dispatch_on_set({'observed_residues_ratio'})
    def pdb_observedResiduesRatio(cls, files: Union[Iterable[str], Iterator[str]], export_format: str = 'tsv') -> Generator:
        jsonDataGenerator = cls.yieldJsonDataFromFiles(files)
        ob = tablib.Dataset()
        count = -1
        for data in jsonDataGenerator:
            for pdb in data:
                for entity_id, entity in data[pdb].items():
                    entity_len = len(entity)
                    append_data = zip(
                        ([entity_id]*entity_len, [pdb]*entity_len),
                        ('entity_id', 'pdb_id'))
                    count += 1
                    yield convertJson2other(entity, append_data, ob, export_format, count)

    @classmethod
    @dispatch_on_set({'residue_listing'})
    def pdb_residueListing(cls, files: Union[Iterable[str], Iterator[str]], export_format: str = 'tsv') -> Generator:
        jsonDataGenerator = cls.yieldJsonDataFromFiles(files)
        ob = tablib.Dataset()
        count = -1
        for data in jsonDataGenerator:
            for pdb in data:
                molecules = data[pdb]['molecules']
                for entity in molecules:
                    chains = entity['chains']
                    for chain in chains:
                        residues = chain['residues']
                        residues_len = len(residues)
                        append_data = zip(
                                (
                                    [chain['chain_id']]*residues_len, 
                                    [chain['struct_asym_id']]*residues_len, 
                                    [entity['entity_id']]*residues_len, 
                                    [pdb]*residues_len),
                                (
                                    'chain_id', 
                                    'struct_asym_id', 
                                    'entity_id', 
                                    'pdb_id'))
                        count += 1
                        yield convertJson2other(residues, append_data, ob, export_format, count)

    @classmethod
    @dispatch_on_set({'secondary_structure'})
    def pdb_secondaryStructure(cls, files: Union[Iterable[str], Iterator[str]]):
        def yieldDfrm(jsonDataGenerator):
            for data in jsonDataGenerator:
                for pdb in data:
                    molecules = data[pdb]['molecules']
                    for entity in molecules:
                        chains = entity['chains']
                        for chain in chains:
                            secondary_structure = chain['secondary_structure']
                            for name, fragment in secondary_structure.items():
                                dfrm = pd.DataFrame(fragment)
                                dfrm['pdb_id'] = pdb
                                dfrm['entity_id'] = entity['entity_id']
                                dfrm['struct_asym_id'] = chain['struct_asym_id']
                                dfrm['chain_id'] = chain['chain_id']
                                dfrm['secondary_structure'] = name
                                yield dfrm

        jsonDataGenerator = cls.yieldJsonDataFromFiles(files)
        return pd.concat((df for df in yieldDfrm(jsonDataGenerator)), ignore_index=True, sort=False)

    @classmethod
    @dispatch_on_set({'binding_sites'})
    def pdb_bindingSites(cls, files: Union[Iterable[str], Iterator[str]]):
        def yieldDfrm(jsonDataGenerator):
            for data in jsonDataGenerator:
                for pdb in data:
                    for site in data[pdb]:
                        site_dfrm = pd.DataFrame(site['site_residues'])
                        site_dfrm['residues_type'] = 'site_residues'
                        ligand_dfrm = pd.DataFrame(site['ligand_residues'])
                        ligand_dfrm['residues_type'] = 'ligand_residues'
                        dfrm = pd.concat([site_dfrm, ligand_dfrm],
                                         ignore_index=True, sort=False)
                        dfrm['pdb_id'] = pdb
                        dfrm['site_id'] = site['site_id']
                        dfrm['evidence_code'] = site['evidence_code']
                        dfrm['details'] = site['details']
                        yield dfrm

        jsonDataGenerator = cls.yieldJsonDataFromFiles(files)
        return pd.concat((df for df in yieldDfrm(jsonDataGenerator)), ignore_index=True, sort=False)

    @classmethod
    @dispatch_on_set({'assembly'})
    def pdb_assembly(cls, files: Union[Iterable[str], Iterator[str]]):
        def yieldDfrm(jsonDataGenerator):
            for data in jsonDataGenerator:
                for pdb in data:
                    for biounit in data[pdb]:
                        dfrm = pd.DataFrame(biounit['entities'])
                        dfrm['pdb_id'] = pdb
                        for key in biounit:
                            cur = biounit[key]
                            if not isinstance(cur, list):
                                dfrm[key] = biounit[key]
                        yield dfrm

        jsonDataGenerator = cls.yieldJsonDataFromFiles(files)
        return pd.concat((df for df in yieldDfrm(jsonDataGenerator)), ignore_index=True, sort=False)

    @classmethod
    @dispatch_on_set({'files'})
    def pdb_files(cls, files: Union[Iterable[str], Iterator[str]]):
        def yieldDfrm(jsonDataGenerator):
            for data in jsonDataGenerator:
                for pdb in data:
                    for key in data[pdb]:
                        for innerKey in data[pdb][key]:
                            record = data[pdb][key][innerKey]
                            if record:
                                dfrm = pd.DataFrame(record)
                                dfrm['pdb_id'] = pdb
                                dfrm['innerKey'] = innerKey
                                dfrm['key'] = key
                                yield dfrm
                            else:
                                continue

        jsonDataGenerator = cls.yieldJsonDataFromFiles(files)
        return pd.concat((df for df in yieldDfrm(jsonDataGenerator)), ignore_index=True, sort=False)


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
        
        '''
        TODO: something need to fix
            1. a more pretty way to append data
            2. export correct json-format dict in tsv/csv format
            3. current not a universal function
        '''
    except KeyError:
        log_func(f'Not a valid json-string/dict to convert format via `tablib.Dataset.export`: {data}')
    except JSONDecodeError:
        log_func(f'Invalid json string: {data}')
    except InvalidDimensions:
        log_func(f'Invalid append_data: {append_data}')
    except UnsupportedFormat:
        log_func(f'Invalid export_format: {export_format}')


if __name__ == "__main__":
    files = [
        r'C:\OmicData\LiGroupWork\PDBeAPI\0117\4w9p_residue_listing.json',
        r'C:\OmicData\LiGroupWork\PDBeAPI\0117\1fm9_residue_listing.json']
    demo = PDBeJsonDecoder.pdb_residueListing(files, 'csv')
    for i in demo:
        print(i)
