# @Created Date: 2020-02-29 04:17:36 pm
# @Filename: oligomer.py
# @Email:  1730416009@stu.suda.edu.cn
# @Author: ZeFeng Zhu
# @Last Modified: 2020-02-29 10:37:06 pm
# @Copyright (c) 2020 MinghuiGroup, Soochow University
import ujson as json
from collections import defaultdict
from itertools import product, combinations, combinations_with_replacement
from typing import Dict, Iterable, Union, Set, List


def validateSIFTS(data: Dict):
    '''
    验证同一Entry下, 每个Entity每条链下的Isoform情况相同
    '''
    for pdb_id, entries in data.items():
        for entry_id, entities in entries.items():
            for entity_id, chains in entities.items():
                isoformSet = set(tuple(isoforms.keys()) for chain_id, isoforms in chains.items())
                if len(isoformSet) > 1:
                    raise ValueError(f'{pdb_id}_{entry_id}_{entity_id}: {isoformSet}')


def range2Set(rangeLyst: Iterable):
    rangeSet = set()
    for res in rangeLyst:
        rangeSet = rangeSet | set(range(res[0], res[1]+1))
    return rangeSet


def jaccardIndex(range_a: Union[str, Iterable], range_b: [str, Iterable]):
    if isinstance(range_a, str):
        range_a = range2Set(json.loads(range_a))
    if isinstance(range_b, str):
        range_b = range2Set(json.loads(range_b))
    elif isinstance(range_a, Set):
        pass
    elif isinstance(range_b, Set):
        pass
    elif isinstance(range_a, Iterable):
        range_a = range2Set(range_a)
    elif isinstance(range_b, Iterable):
        range_b = range2Set(range_b)
    
    return len(range_a & range_b)/len(range_a | range_b)


def yieldHe(data: Dict):
    for pdb_id, entries in data.items():
        if len(entries) == 1:
            continue
        else:
            for entry_l, entry_r in combinations(entries, 2):
                eitities_l, eitities_r = entries[entry_l], entries[entry_r]
                for entity_id_l, entity_id_r in product(eitities_l, eitities_r):
                    chains_l, chains_r = eitities_l[entity_id_l], eitities_r[entity_id_r]
                    for chain_id_l, chain_id_r in product(chains_l, chains_r):
                        isoforms_l, isoforms_r = chains_l[chain_id_l], chains_r[chain_id_r]
                        for isoform_id_l, isoform_id_r in product(isoforms_l, isoforms_r):
                            unit_l, unit_r = isoforms_l[isoform_id_l], isoforms_r[isoform_id_r]
                            if entity_id_l == entity_id_r:
                                jaccard = jaccardIndex(
                                    unit_l['sifts_unp_range'], unit_r['sifts_unp_range'])
                            else:
                                jaccard = 0
                            yield pdb_id, (entry_l, entry_r), (entity_id_l, entity_id_r), (chain_id_l, chain_id_r), (isoform_id_l, isoform_id_r), jaccard


def yieldHo(data):
    '''
    同一PDB下, 内任意两条链是否属于同一蛋白, 且覆盖范围是否相似
    还需验证，不同Entry下对应的entity不同
    '''
    for pdb_id, entries in data.items():
        for entry_id, entities in entries.items():
            for entity_id_l, entity_id_r in combinations_with_replacement(entities.keys(), 2):
                if entity_id_l != entity_id_r:
                    chains_l, chains_r = entities[entity_id_l], entities[entity_id_r]
                    chain_ids_l, chain_ids_r = chains_l.keys(), chains_r.keys()
                    for chain_id_l, chain_id_r in product(chain_ids_l, chain_ids_r):
                        chain_l, chain_r = chains_l[chain_id_l], chains_r[chain_id_r]
                        isoform_ids_l, isoform_ids_r = chain_l.keys(), chain_r.keys()
                        assert isoform_ids_l == isoform_ids_r
                        for isoform_id in isoform_ids_l:
                            yield pdb_id, entry_id, 'ho', False, (chain_id_l, chain_id_r), jaccardIndex(chain_l[isoform_id]['sifts_unp_range'], chain_r[isoform_id]['sifts_unp_range'])
                else:
                    chain_ids = entities[entity_id_l].keys()
                    for res in combinations(chain_ids, 2):
                        yield pdb_id, entry_id, 'ho', True, res, 1
