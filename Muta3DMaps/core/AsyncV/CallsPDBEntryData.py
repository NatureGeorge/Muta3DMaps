# @Created Date: 2020-01-12 01:27:18 pm
# @Filename: CallsPDBEntryData.py
# @Email:  1730416009@stu.suda.edu.cn
# @Author: ZeFeng Zhu
# @Last Modified: 2020-01-12 07:29:56 pm
# @Copyright (c) 2020 MinghuiGroup, Soochow University

'''
[REST calls based on PDB entry data](https://www.ebi.ac.uk/pdbe/api/doc/pdb.html)

Note that at present `mutated_AA_or_NA` does not provide information
about mutated nucleotides in RNA or DNA chains, but it would do so
in near future.(2020-01-12)
'''

import asyncio
import aiohttp
from aiohttp import web
import pandas as pd
import numpy as np
from collections import Counter
from collections.abc import Iterable, Iterator
import re
import json
import time
import os

try:
    from .Logger import RunningLogger
except Exception:
    from Logger import RunningLogger

BASE_URL = "https://www.ebi.ac.uk/pdbe/api/pdb/entry/"
CHAIN_suffix = "/chain"
'''
SUMMARY = "summary/{}"
MOLECULES = "molecules/{}"
EXPERIMENT = "experiment/{}"
# PUBLICATIONS = "publications/{}"
# RELATED_PUBLICATIONS"related_publications/{}"
# NMR_RESOURCES = "nmr_resources/{}"
LIGAND_MONOMERS = "ligand_monomers/{}"
MODIFIED_AA_or_NA = "modified_AA_or_NA/{}"
MUTATED_AA_or_NA = "mutated_AA_or_NA/{}"
RELEASE_STATUS = "status/{}"
OBSERVED_RANGES = "polymer_coverage/{}"
OBSERVED_RANGES_CHAIN = OBSERVED_RANGES + CHAIN_suffix
SECONDARY_STRUCTURE = "secondary_structure/{}"
RESIDUE_LISTING = "residue_listing/{}"
RESIDUE_LISTING_CHAIN = RESIDUE_LISTING + CHAIN_suffix
BINDING_SITES = "binding_sites/{}"
FILES = "files/{}"
OBESRVED_RESIDUES_RATIO = "observed_residues_ratio/{}"
ASSEMBLY = "assembly{}"
ELECTRON_DENSITY_STATISTICS = "electron_density_statistics/{}"
FUNCTIONAL_ANNOT_COFACTOR = "cofactor/{}"
DRUGBANK_ANNOT = "drugbank/{}"
RELATED_EXPERIMENT_DATA = "related_experiment_data/{}"
'''
API_LYST = sorted(["summary", "molecules", "experiment", "ligand_monomers",
            "modified_AA_or_NA", "mutated_AA_or_NA", "status",
            "polymer_coverage", "secondary_structure",
            "residue_listing", "binding_sites", "files", "observed_residues_ratio",
            "assembly", "electron_density_statistics",
            "cofactor", "drugbank", "related_experiment_data"])


class PathNode(object):
    '''
    features:
        * Provide a geneator that yield particular (url, filePath) for a pdb with specified api
        * Record files' status of the pdb

    :param pdb: pdb id
    :type pdb: str
    '''
    logpath = None

    def __init__(self, pdb):
        self.id = pdb
        self.api_status = {}

    def __repr__(self):
        return "PathNode<id:{}, api_status:{}>".format(self.id, self.api_status)
    
    @property
    def paths(self):
        for api in API_LYST:
            yield "{}{}/{}".format(BASE_URL, api, self.id), "{}_{}.json".format(self.id, api)
    
    def log(self):
        with open(self.logpath, 'a+') as outFile:
            outFile.write("{}\t{}\n".format(
                self.id,
                "\t".join(str(self.api_status[key]) for key in API_LYST)
            ))

    @classmethod
    def main(cls, pdbs, logpath):
        cls.logpath = logpath
        for pdb in pdbs:
            cur = cls(pdb)
            yield cur
        

class RetrieveEntryData(object):
    '''
    Retrieve PDB entry data
    -----------------------

    reference: https://www.ebi.ac.uk/pdbe/api/doc/pdb.html

    feature:
        * async (base on event loop)
        * organize data within one entry into a tree structure
    '''
    
    def __init__(self, pdbs=None, loggingPath=None, workdir=None, overviewFileName=None):
        if isinstance(pdbs, str):
            pdbs = [pdbs]
        self.pdbs = pdbs
        self.Logger = RunningLogger("RetrieveEntryData", loggingPath)
        self.workdir = workdir
        self.overview = os.path.join(workdir, overviewFileName)
    
    @staticmethod
    async def retrieve(session, url, pdbNode):
        # url = BASE_URL + sub_url.format(pdb)
        async with session.get(url) as resp:
            token = url[len(BASE_URL):-5]
            if resp.status == 200:
                pdbNode.api_status[token] = 1
                return await resp.read()
            elif resp.status == 404:
                # raise web.HTTPNotFound()
                pdbNode.api_status[token] = 0
                return None
            else:
                mes = "code={resp.status}, message={resp.reason}, headers={resp.headers}"
                raise Exception(mes.format(resp=resp))
    
    @staticmethod
    def save(path, data):
        with open(path, 'wb') as fp:
            fp.write(data)

    async def unit(self, pdbNode, url, path, session, semaphore):
        async with semaphore:
            self.Logger.logger.info("Start to get data")
            rawData = await self.retrieve(session, url, pdbNode)
            if rawData is not None:
                self.save(os.path.join(self.workdir, path), rawData)
            else:
                self.Logger.logger.info("{} Not Found {}".format(pdbNode.id, path[4:]))
        if len(pdbNode.api_status) == len(API_LYST):
            # self.Logger.logger.info(str(pdbNode))
            pdbNode.log()
        return url

    async def multi(self, nodes, concur_req):
        semaphore = asyncio.Semaphore(concur_req)
        async with aiohttp.ClientSession() as session:
            res = await asyncio.gather(
                *[asyncio.create_task(self.unit(pdbNode, url, path, session, semaphore))
                    for pdbNode in nodes for url, path in pdbNode.paths])

        return len(res)

    def main(self, concur_req=100):
        nodes = PathNode.main(self.pdbs, self.overview)
        t0 = time.perf_counter()
        count = asyncio.run(self.multi(nodes, concur_req))
        elapsed = time.perf_counter() - t0
        self.Logger.logger.info(
            '\n{} entries downloaded in {:.2f}s'.format(count, elapsed))


class PDBeJsonDecoder(object):
    @staticmethod
    def yieldJsonDataFromFiles(files):
        for file in files:
            yield json.load(open(file, 'rt'))

    @classmethod
    def pdb_Root(cls, files, fileType):
        '''
        def wrapper(jsonDataGenerator):
            for data in jsonDataGenerator:
                pdb = next(iter(data))
                data[pdb][0]['pdb_id'] = pdb
                yield pdb, data
        '''

        def yieldDfrm(jsonDataGenerator):
            for data in jsonDataGenerator:
                pdb = next(iter(data))
                dfrm = pd.DataFrame(data[pdb])
                dfrm['pdb_id'] = pdb
                yield dfrm
        
        if not set(fileType) <= {"status", "summary", "modified_AA_or_NA",
                                 "cofactor", "molecules", "ligand_monomers", 
                                 "experiment", "electron_density_statistics", 
                                 "related_experiment_data", "drugbank"}:
            return None
        jsonDataGenerator = cls.yieldJsonDataFromFiles(files)
        # dfrm = pd.DataFrame(data[key][0] for key, data in wrapper(jsonDataGenerator))
        return pd.concat((df for df in yieldDfrm(jsonDataGenerator)), ignore_index=True, sort=False)

    @classmethod
    def pdb_polymerCoverage(cls, files):
        # if not set(fileType) <= {"polymer_coverage", "observed_residues_ratio"}:
        # return None
        def yieldDfrm(jsonDataGenerator):
            for data in jsonDataGenerator:
                pdb = next(iter(data))
                molecules = data[pdb]['molecules']
                for entity in molecules:
                    dfrm = pd.DataFrame(entity['chains'])
                    dfrm['entity_id'] = entity['entity_id']
                    dfrm['pdb_id'] = pdb
                    yield dfrm
        
        jsonDataGenerator = cls.yieldJsonDataFromFiles(files)
        return pd.concat((df for df in yieldDfrm(jsonDataGenerator)), ignore_index=True, sort=False)
    
    @classmethod
    def pdb_observedResiduesRatio(cls, files):
        def yieldDfrm(jsonDataGenerator):
            for data in jsonDataGenerator:
                pdb = next(iter(data))
                for entity_id, entity in data[pdb].items():
                    dfrm = pd.DataFrame(entity)
                    dfrm['entity_id'] = entity_id
                    dfrm['pdb_id'] = pdb
                    yield dfrm
        
        jsonDataGenerator = cls.yieldJsonDataFromFiles(files)
        return pd.concat((df for df in yieldDfrm(jsonDataGenerator)), ignore_index=True, sort=False)

    @classmethod
    def pdb_mutatedAAorNA(cls, files):
        def yieldDfrm(jsonDataGenerator):
            for data in jsonDataGenerator:
                pdb = next(iter(data))
                for entity in data[pdb]:
                    dfrm = pd.DataFrame(entity)
                    dfrm['pdb_id'] = pdb
                    dfrm['mutation_from'], dfrm['mutation_to'], dfrm['mutation_type'] = dfrm['mutation_details'].apply(lambda info: (info['from'], info['to'], info['type'])).str
                    yield dfrm
        
        jsonDataGenerator = cls.yieldJsonDataFromFiles(files)
        return pd.concat((df for df in yieldDfrm(jsonDataGenerator)), ignore_index=True, sort=False)

    @classmethod
    def pdb_residueListing(cls, files):
        def yieldDfrm(jsonDataGenerator):
            for data in jsonDataGenerator:
                pdb = next(iter(data))
                molecules = data[pdb]['molecules']
                for entity in molecules:
                    chains = entity['chains']
                    for chain in chains:
                        dfrm = pd.DataFrame(chain['residues'])
                        dfrm['pdb_id'] = pdb
                        dfrm['entity_id'] = entity['entity_id']
                        dfrm['struct_asym_id'] = chain['struct_asym_id']
                        dfrm['chain_id'] = chain['chain_id']
                        yield dfrm

        jsonDataGenerator = cls.yieldJsonDataFromFiles(files)
        return pd.concat((df for df in yieldDfrm(jsonDataGenerator)), ignore_index=True, sort=False)

    '''
    @staticmethod
    def splitDictInDfrm(dfrm, dictColName):
        newColumns = dfr
        dfrm[] = dfrm[dictColName].apply(lambda x: tuple(x.values())
    '''

    @classmethod
    def pdb_secondaryStructure(cls, files):
        def yieldDfrm(jsonDataGenerator):
            for data in jsonDataGenerator:
                pdb = next(iter(data))
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
    def pdb_bindingSites(cls, files):
        def yieldDfrm(jsonDataGenerator):
            for data in jsonDataGenerator:
                pdb = next(iter(data))
                sites = data[pdb]
                for site in sites:
                    site_dfrm = pd.DataFrame(site['site_residues'])
                    site_dfrm['residues_type'] = 'site_residues'
                    ligand_dfrm = pd.DataFrame(site['ligand_residues'])
                    ligand_dfrm['residues_type'] = 'ligand_residues'
                    dfrm = pd.concat([site_dfrm, ligand_dfrm], ignore_index=True, sort=False)
                    dfrm['pdb_id'] = pdb
                    dfrm['site_id'] = site['site_id']
                    dfrm['evidence_code'] = site['evidence_code']
                    dfrm['details'] = site['details']
                    yield dfrm

        jsonDataGenerator = cls.yieldJsonDataFromFiles(files)
        return pd.concat((df for df in yieldDfrm(jsonDataGenerator)), ignore_index=True, sort=False)

    @classmethod
    def pdb_assembly(cls, files):
        def yieldDfrm(jsonDataGenerator):
            for data in jsonDataGenerator:
                pdb = next(iter(data))
                for biounit in data[pdb]:
                    dfrm = pd.DataFrame(biounit['entities'])
                    dfrm['pdb_id'] = pdb
                    dfrm['polymeric_count'] = biounit['polymeric_count']
                    dfrm['assembly_composition'] = biounit['assembly_composition']
                    dfrm['molecular_weight'] = biounit['molecular_weight']
                    dfrm['details'] = biounit['details']
                    dfrm['assembly_id'] = biounit['assembly_id']
                    yield dfrm

        jsonDataGenerator = cls.yieldJsonDataFromFiles(files)
        return pd.concat((df for df in yieldDfrm(jsonDataGenerator)), ignore_index=True, sort=False)

    @classmethod
    def pdb_files(cls, files):
        jsonDataGenerator = cls.yieldJsonDataFromFiles(files)
        def yieldDfrm(jsonDataGenerator):
            for data in jsonDataGenerator:
                pdb = next(iter(data))
                info = data[pdb]
                for key in info:
                    for innerKey in info[key]:
                        record = info[key][innerKey]
                        if record:
                            dfrm = pd.DataFrame(record)
                            dfrm['innerKey'] = innerKey
                            dfrm['key'] = key
                            yield dfrm
                        else:
                            continue
        
        jsonDataGenerator = cls.yieldJsonDataFromFiles(files)
        return pd.concat((df for df in yieldDfrm(jsonDataGenerator)), ignore_index=True, sort=False)


if __name__ == "__main__":
    demo = RetrieveEntryData(pdbs=["5o8b"],  # ["1a01", "5dfz", "2hev", "2hey", "5hkr", "4ddg", "4v5j", "3g96", "5hht"]
                             loggingPath="C:/OmicData/LiGroupWork/PDBeAPI/log.log",
                             workdir="C:/OmicData/LiGroupWork/PDBeAPI/test/",
                             overviewFileName="PDB_Entry_File_overview.tsv")
    demo.main()
    
