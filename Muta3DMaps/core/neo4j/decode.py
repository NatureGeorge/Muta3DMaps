# @Created Date: 2020-04-08 09:44:01 am
# @Filename: decode.py
# @Email:  1730416009@stu.suda.edu.cn
# @Author: ZeFeng Zhu
# @Last Modified: 2020-04-08 09:44:05 am
# @Copyright (c) 2020 MinghuiGroup, Soochow University
from typing import List, Iterable, Iterator, Union, Dict
from neo4j import GraphDatabase
import pandas as pd
from pathlib import Path
from collections import defaultdict
from itertools import combinations, product
from numpy import nan
import numpy as np
import ujson as json
from functools import lru_cache
from textdistance import jaccard, overlap
from Bio import Align
from Bio.SubsMat import MatrixInfo as matlist

SEQ_DICT = {
    "GLY": "G", "ALA": "A", "SER": "S", "THR": "T", "CYS": "C", "VAL": "V", "LEU": "L",
    "ILE": "I", "MET": "M", "PRO": "P", "PHE": "F", "TYR": "Y", "TRP": "W", "ASP": "D",
    "GLU": "E", "ASN": "N", "GLN": "Q", "HIS": "H", "LYS": "K", "ARG": "R"}

standardAA = list(SEQ_DICT.keys())

standardNu = ['DA', 'DT', 'DC', 'DG', 'A', 'U', 'C', 'G']


def to_interval(lyst: Union[Iterable, Iterator]) -> Union[float, List]:
    def pass_check(lyst):
        if not lyst or pd.isna(lyst):
            return False
        else:
            return True
    if not pass_check(lyst): return nan
    else:
        lyst = set(lyst)
        if not pass_check(lyst): return nan
        start = []
        interval_lyst = []
        true_interval_lyst = []
        max_edge = max(lyst)
        min_edge = min(lyst)
        if len(lyst) == (max_edge + 1 - min_edge):
            true_interval_lyst.append((min_edge, max_edge))
        else:
            lyst_list = sorted(lyst)
            for j in lyst_list:
                if not start:
                    i = j
                    start.append(j)
                    i += 1
                else:
                    if (i != j) or (j == max(lyst_list)):
                        if j == max(lyst_list):
                            if (i != j):
                                interval_lyst.append(start)
                                interval_lyst.append([j])
                                break
                            else:
                                start.append(j)
                        interval_lyst.append(start)
                        start = [j]
                        i = j + 1
                    else:
                        start.append(j)
                        i += 1
            for li in interval_lyst:
                max_edge = max(li)
                min_edge = min(li)
                true_interval_lyst.append((min_edge, max_edge))
        return true_interval_lyst


def interval2set(lyst: Union[Iterable, Iterator, str]):
    if isinstance(lyst, str):
        lyst = json.loads(lyst)
    range_set = set()
    for left, right in lyst:
        range_set = range_set | set(range(left, right+1))
    return range_set


def subtract_range(pdb_range: Union[str, Iterable], mis_range: Union[str, Iterable]) -> List:
    if isinstance(mis_range, float):
        return pdb_range
    pdb_range_set = interval2set(pdb_range)
    mis_range_set = interval2set(mis_range)
    return to_interval(pdb_range_set - mis_range_set)


def overlap_range(obs_range:Union[str, Iterable], unk_range: Union[str, Iterable]) -> List:
    if isinstance(unk_range, float):
        return nan
    obs_range_set = interval2set(obs_range)
    unk_range_set = interval2set(unk_range)
    return to_interval(obs_range_set & unk_range_set)


def outside_range_len(pdb_range: str, seqres_len: int, omit: int = 5) -> int:
    lyst = json.loads(pdb_range)
    out_head = lyst[0][0]-1
    out_tail = seqres_len - lyst[-1][-1]
    if out_head <= omit:
        out_head = 0
    else:
        out_head -= omit
    if out_tail <= omit:
        out_tail = 0
    else:
        out_tail -= omit
    return out_head + out_tail


def range_len(lyst: Union[List, str, float]) -> int:
    if isinstance(lyst, float):
        return 0
    elif isinstance(lyst, str):
        lyst = json.loads(lyst)
    length = 0
    for left, right in lyst:
        length += right - left + 1
    return length


def lyst2dict(lyst: List) -> Dict:
    res = dict(lyst)
    for key in res.keys():
        res[key] = dict(res[key])
    return res


def sub_index(init_index, subtract_index) -> pd.Index:
    if len(subtract_index) == 0:
        return init_index
    else:
        return pd.Index(set(init_index)-set(subtract_index))


class Entry(object):
    @classmethod
    def set_session(cls, session):
        cls.session = session
    
    @staticmethod
    def to_data_frame(res):
        return pd.DataFrame(dict(i) for i in res)
    
    @classmethod
    def summary_method(cls, pdbs, session=None):
        if session is None:
            session = cls.session
        query = '''
            MATCH (entry:Entry)-[:EXPERIMENT]->(method:Method)
            WHERE entry.ID in $lyst
            RETURN entry.ID as pdb_id, entry.PDB_REV_DATE_ORIGINAL, entry.FIRST_REV_DATE, entry.PDB_REV_DATE, entry.REVISION_DATE, entry.RESOLUTION, method.METHOD_CLASS
        '''
        return session.run(query, lyst=list(pdbs))
    
    @classmethod
    def summary_ligand(cls, pdbs, session=None):
        if session is None:
            session = cls.session
        query = '''
            MATCH (entry:Entry)-[:HAS_BOUND_MOLECULE]->(bmol:BoundMolecule)<-[:IS_PART_OF]-(bli:BoundLigand) 
            WHERE entry.ID in $lyst
            RETURN entry.ID as pdb_id, COUNT(bli) as BOUND_LIGAND_COUNT, COUNT(distinct bmol) as BOUND_MOL_COUNT
        '''
        return session.run(query, lyst=list(pdbs))
    
    @classmethod
    def summary_nucleotides(cls, pdbs, session=None):
        if session is None:
            session = cls.session
        query = '''
            MATCH (entry:Entry)-[:HAS_ENTITY]->(entity:Entity)-[:HAS_PDB_RESIDUE]->(pdbRes:PDBResidue)
            WHERE entity.POLYMER_TYPE IN ['D','R','D/R'] AND entry.ID in $lyst
            RETURN 
                entry.ID as pdb_id, 
                entity.ID as entity_id, 
                entity.POLYMER_TYPE as polymer_type,
                size([x in collect(pdbRes.CHEM_COMP_ID) where x in ['DA','DT','DC','DG']]) as DNA_COUNT,
                size([x in collect(pdbRes.CHEM_COMP_ID) where x in ['A','U','C','G']]) as RNA_COUNT,
                size([x in collect(pdbRes.CHEM_COMP_ID) where not x in $standardNu]) as OTHER_COUNT
        '''
        return session.run(query, lyst=list(pdbs), standardNu=standardNu)
    
    @staticmethod
    def get_polymer_type(dna, rna):
        if dna and rna:
            return 'D/R'
        elif dna and not rna:
            return 'D'
        elif not dna and rna:
            return 'R'
        elif not dna and not rna:
            return nan

    @classmethod
    def deal_nucleotides(cls, res):
        dfrm = cls.to_data_frame(res)
        dfrm.polymer_type = dfrm.apply(lambda x: cls.get_polymer_type(x['DNA_COUNT'], x['RNA_COUNT']), axis=1)
        return pd.DataFrame(
            ((pdb_id, json.dumps(dict(zip(data.entity_id, data.polymer_type)))) for pdb_id, data in res.groupby('pdb_id')),
            columns=('pdb_id', 'nucleotides_entity_type'))

    @classmethod
    def summary_seq(cls, pdbs, session=None):
        if session is None:
            session = cls.session
        query = '''
            MATCH (entry:Entry)-[:HAS_ENTITY]->(entity:Entity{POLYMER_TYPE:'P'})-[:CONTAINS_CHAIN]->(chain:Chain)-[inChain:IS_IN_CHAIN]-(res:PDBResidue)
            WHERE entry.ID in $lyst
            RETURN 
                entry.ID as pdb_id,
                entity.ID as entity_id,
                chain.AUTH_ASYM_ID as chain_id,
                COUNT(res) as SEQRES_COUNT,
                avg(tofloat(inChain.OBSERVED_RATIO)) as AVG_OBS_RATIO,
                [x in COLLECT([res.ID, inChain.AUTH_COMP_ID]) WHERE NOT x[1] IN $standardAA] as NON_INDEX,
                [x in COLLECT([res.ID, inChain.AUTH_COMP_ID]) WHERE x[1] = 'UNK'] as UNK_INDEX,
                [x in COLLECT([res.ID, inChain.OBSERVED]) WHERE x[1]='N'] as MIS_INDEX
        '''
        return session.run(query, lyst=list(pdbs), standardAA=standardAA)

    @staticmethod
    def deal_seq_index(dfrm: pd.DataFrame) -> pd.DataFrame:
        def index2range(lyst):
            res = to_interval(int(i[0]) for i in lyst)
            if isinstance(res, List):
                return json.dumps(res)
            else:
                return res
        
        def get_obs_rel_count(seqres_count, mis_range, unk_range, non_range):
            seq_range = [(1, seqres_count)]
            obs_range = subtract_range(seq_range, mis_range)
            obs_unk_range = overlap_range(obs_range, unk_range)
            obs_pure_range = subtract_range(obs_range, non_range)
            return range_len(obs_unk_range), range_len(obs_pure_range)

        dfrm['UNK_COUNT'] = dfrm.UNK_INDEX.apply(len)
        dfrm['PURE_SEQRES_COUNT'] = dfrm.SEQRES_COUNT - dfrm.NON_INDEX.apply(len)
        dfrm['OBS_RECORD_COUNT'] = dfrm.SEQRES_COUNT - dfrm.MIS_INDEX.apply(len)
        dfrm.NON_INDEX = dfrm.NON_INDEX.apply(index2range)
        dfrm.UNK_INDEX = dfrm.UNK_INDEX.apply(index2range)
        dfrm.MIS_INDEX = dfrm.MIS_INDEX.apply(index2range)
        dfrm[['OBS_UNK_COUNT', 'ATOM_RECORD_COUNT']] = dfrm.apply(lambda x: get_obs_rel_count(
            x['SEQRES_COUNT'], x['MIS_INDEX'], x['UNK_INDEX'], x['NON_INDEX']), axis=1, result_type='expand')
        return dfrm

    @classmethod
    @lru_cache()
    def get_seqres(cls, pdb_id: str, entity_id: str = None, chain_id: str = None, three2one: bool = True, session=None):
        if session is None:
            session = cls.session
        if entity_id is not None:
            query = '''
                MATCH (entity:Entity)-[:HAS_PDB_RESIDUE]->(pdbRes:PDBResidue)
                WHERE entity.UNIQID = $uniqid
                RETURN pdbRes.CHEM_COMP_ID AS residue_name ORDER BY toInteger(pdbRes.ID)
            '''
            res = session.run(query, uniqid=f'{pdb_id}_{entity_id}')
        elif chain_id is not None:
            query = '''
                MATCH (entry:Entry)-[:HAS_ENTITY]-(:Entity)-[:CONTAINS_CHAIN]-(chain:Chain)-[:IS_IN_CHAIN]-(pdbRes:PDBResidue)
                WHERE entry.ID = $pdb_id AND chain.AUTH_ASYM_ID = $chain_id
                RETURN pdbRes.CHEM_COMP_ID AS residue_name ORDER BY toInteger(pdbRes.ID)
            '''
            res = session.run(query, pdb_id=pdb_id, chain_id=chain_id)
        else:
            raise ValueError('please specify entity_id or chain_id')
        
        if three2one:
            return ''.join(SEQ_DICT.get(r['residue_name'], 'X') for r in res)
        else:
            return res

    @classmethod
    def get_residues(cls, pdb_id: str, entity_id: str, chain_id: str, res_ids=None, observed_only:bool=False, session=None):
        if session is None:
            session = cls.session
        query = '''
            MATCH (entry:Entry)-[:HAS_ENTITY]->(entity:Entity)-[:CONTAINS_CHAIN]->(chain:Chain)-[inChain:IS_IN_CHAIN]-(res:PDBResidue)
            WHERE entry.ID = $pdb_id AND entity.ID = $entity_id AND chain.AUTH_ASYM_ID = $chain_id {}
            RETURN entry.ID as pdb_id, entity.ID as entity_id, chain.AUTH_ASYM_ID as chain_id, res.CHEM_COMP_ID as residue_name, toInteger(res.ID) as residue_number, tofloat(inChain.OBSERVED_RATIO) as obs_ratio, inChain.AUTH_SEQ_ID as author_residue_number, inChain.PDB_INS_CODE as author_insertion_code ORDER BY residue_number
        '''
        if res_ids is None:
            if observed_only:
                query = query.format("AND inChain.OBSERVED = 'Y'")
            else:
                query = query.format('')
            return session.run(query, pdb_id=pdb_id, entity_id=str(entity_id), chain_id=chain_id)
        else:
            if observed_only:
                query = query.format("AND res.ID IN $res_ids AND inChain.OBSERVED = 'Y'")
            else:
                query = query.format('AND res.ID IN $res_ids')
            res_ids = [str(i) for i in res_ids]
            return session.run(query, pdb_id=pdb_id, entity_id=str(entity_id), chain_id=chain_id, res_ids=res_ids)

    @classmethod
    def summary_entity_chain(cls, pdbs, session=None):
        if session is None:
            session = cls.session
        query = """
                MATCH (entry:Entry)-[:HAS_ENTITY]-(entity:Entity{POLYMER_TYPE:'P'})-[:CONTAINS_CHAIN]-(chain:Chain)
                WHERE entry.ID IN $pdbs
                WITH entry, entity, COLLECT(chain.AUTH_ASYM_ID) as chain_ids
                RETURN entry.ID as pdb_id, COLLECT([entity.ID, chain_ids]) as entity_chain_map
            """
        return session.run(query, pdbs=list(pdbs))
    
    @classmethod
    def deal_entity_chain(cls, res):
        dfrm = cls.to_data_frame(res)
        dfrm.entity_chain_map = dfrm.entity_chain_map.apply(dict)
        dfrm['entity_count'] = dfrm.entity_chain_map.apply(len)
        dfrm['chain_count'] = dfrm.entity_chain_map.apply(lambda x: sum(len(i) for i in x.values()))
        return dfrm
 

class SIFTS(Entry):
    @classmethod
    def summary_entity_unp(cls, pdbs, session=None):
        if session is None:
            session = cls.session
        query = """
                MATCH (entry:Entry)-[:HAS_ENTITY]-(entity:Entity{POLYMER_TYPE:'P'})-[unp_rel:HAS_UNIPROT{BEST_MAPPING:"1"}]-(unp:UniProt)
                WHERE entry.ID in $pdbs
                WITH entry, entity, unp
                MATCH (unp)-[seg:HAS_UNIPROT_SEGMENT]-(entity)
                WITH COLLECT(DISTINCT [seg.UNP_START, seg.UNP_END]) as range_info, unp, entry, entity
                WITH COLLECT([entity.ID, range_info]) as eneitiy_unp_info, entry, unp
                RETURN entry.ID as pdb_id, COLLECT([unp.ACCESSION, eneitiy_unp_info]) as unp_entity_info
            """
        return session.run(query, pdbs=list(pdbs))
    
    @classmethod
    def deal_entity_unp(cls, res):
        dfrm = cls.to_data_frame(res)
        dfrm.unp_entity_info = dfrm.unp_entity_info.apply(lyst2dict)
        dfrm['unp_count'] = dfrm.unp_entity_info.apply(len)
        dfrm['unp_rel_count'] = dfrm.unp_entity_info.apply(lambda x: sum(len(i) for i in x.values()))
        return dfrm

    @staticmethod
    def overall_check(df: pd.DataFrame):
        error_index = df[df.entity_count < df.unp_rel_count].index
        pass_index = sub_index(df.index, error_index)
        return df.loc[pass_index], df.loc[error_index]
    
    @staticmethod
    def get_raw_oligo_state(df: pd.DataFrame):
        state_name, unmapped_name = 'raw_oligo_state', 'has_unmapped_protein'
        cur_df = df
        # He
        he = cur_df[
            ((cur_df.entity_count > cur_df.unp_rel_count) &
            (cur_df.unp_count.eq(1))) | (cur_df.unp_count.gt(1))
        ].index
        cur_df = cur_df.loc[sub_index(cur_df.index, he)]
        # Ho
        ho = cur_df[
            (cur_df.entity_count == cur_df.unp_rel_count)  # already satisfied
            & (cur_df.unp_count.eq(1))  # already satisfied
            & (cur_df.chain_count.gt(1))
        ].index
        cur_df = cur_df.loc[sub_index(cur_df.index, ho)]
        # Mo
        mo = df[df.chain_count.eq(1)].index
        cur_df_index = sub_index(cur_df.index, mo)
        # Tag
        df[state_name] = nan
        df.loc[mo, state_name] = 'mo'
        df.loc[ho, state_name] = 'ho'
        df.loc[he, state_name] = 'he'
        # Tag
        non_ho = df.loc[sub_index(df.index, ho)]
        unmapped = non_ho[non_ho.entity_count > non_ho.unp_rel_count].index
        df[unmapped_name] = False
        df.loc[unmapped, unmapped_name] = True
        # return
        return df, (mo, ho, he, cur_df_index)

    @classmethod
    def summary_oligo_state(cls, pdbs):
        entity_chain_info = cls.deal_entity_chain(cls.summary_entity_chain(pdbs))
        unp_entity_info = cls.deal_entity_unp(cls.summary_entity_unp(pdbs))
        unp_entity_chain_info = pd.merge(entity_chain_info, unp_entity_info, how='left')
        unp_entity_chain_info.unp_count.fillna(0, inplace=True)
        unp_entity_chain_info.unp_rel_count.fillna(0, inplace=True)
        unp_entity_chain_info.unp_count = unp_entity_chain_info.unp_count.apply(int)
        unp_entity_chain_info.unp_rel_count = unp_entity_chain_info.unp_rel_count.apply(int)
        pass_result, error_result = cls.overall_check(unp_entity_chain_info)
        pass_result_oli, indexes_oli = cls.get_raw_oligo_state(pass_result)

    @classmethod
    def set_from(cls, id_type: str):
        id_type = id_type.lower()
        if id_type in ('unp', 'uniprot', 'unp.ACCESSION'):
            cls.from_str = 'unp.ACCESSION'
        elif id_type in ('pdb', 'pdb_id', 'entry.ID'):
            cls.from_str = 'entry.ID'
        else:
            raise ValueError('unknown id type')

    @classmethod
    def summary_mapping(cls, lyst, id_type: str, session=None):
        if session is None:
            session = cls.session
        cls.set_from(id_type)
        query = '''
            MATCH (entry:Entry)-[:HAS_ENTITY]-(entity:Entity)-[seg:HAS_UNIPROT_SEGMENT]-(unp:UniProt)
            WHERE {} in $lyst
            RETURN unp.ACCESSION as UniProt, entry.ID as pdb_id, entity.ID as entity_id, seg.AUTH_ASYM_ID as chain_id, tofloat(seg.IDENTITY) as identity, COLLECT([toInteger(seg.PDB_START), toInteger(seg.PDB_END)]) as pdb_range, COLLECT([toInteger(seg.UNP_START), toInteger(seg.UNP_END)]) as unp_range
        '''.format(cls.from_str)
        return session.run(query, lyst=list(lyst))
    
    @classmethod
    def deal_mapping(cls, res):
        dfrm = cls.to_data_frame(res)
        dfrm = cls.deal_InDe(dfrm)
        dfrm.pdb_range = dfrm.pdb_range.apply(json.dumps)
        dfrm.unp_range = dfrm.unp_range.apply(json.dumps)
        dfrm = cls.update_range(dfrm)
        return dfrm
    
    @staticmethod
    def sort_2_range(unp_range: List, pdb_range: List):
        unp_range, pdb_range = zip(*sorted(zip(unp_range, pdb_range), key=lambda x: x[0][0]))
        return unp_range, pdb_range

    @classmethod
    def deal_InDe(cls, dfrm: pd.DataFrame) -> pd.DataFrame:
        def get_gap_list(li: List):
            return [li[i+1][0] - li[i][1] - 1 for i in range(len(li)-1)]

        def get_range_diff(lyst_a: List, lyst_b: List):
            array_a = np.array([right - left + 1 for left, right in lyst_a])
            array_b = np.array([right - left + 1 for left, right in lyst_b])
            return (array_a - array_b).tolist()

        def add_tage_to_range(df: pd.DataFrame, tage_name: str):
            # ADD TAGE FOR SIFTS
            df[tage_name] = 'Safe'
            # No Insertion But Deletion[Pure Deletion]
            df.loc[df[(df['group_info'] == 1) & (
            df['unp_pdb_var'] > 0)].index, tage_name] = 'Deletion'
            # Insertion & No Deletion
            df.loc[df[
                (df['group_info'] != 1) &
                (df['var_0_count'] == df['group_info']) &
                (df['unp_gap_0_count'] == (df['group_info'] - 1))].index, tage_name] = 'Insertion'
            # Insertion & Deletion
            df.loc[df[
                (df['group_info'] != 1) &
                ((df['var_0_count'] != df['group_info']) |
                (df['unp_gap_0_count'] != (df['group_info'] - 1)))].index, tage_name] = 'Insertion & Deletion'

        dfrm['group_info'] = dfrm.apply(lambda x: len(
            x['pdb_range']), axis=1)
        
        focus_index = dfrm[dfrm.group_info.gt(1)].index
        focus_df = dfrm.loc[focus_index].apply(lambda x: cls.sort_2_range(x['unp_range'], x['pdb_range']), axis=1, result_type='expand')
        focus_df.index = focus_index
        focus_df.columns = ['unp_range', 'pdb_range']
        dfrm.loc[focus_index, ['unp_range', 'pdb_range']] = focus_df
        
        dfrm['pdb_gap_list'] = dfrm.apply(lambda x: json.dumps(
            get_gap_list(x['pdb_range'])), axis=1)
        dfrm['unp_gap_list'] = dfrm.apply(lambda x: json.dumps(
            get_gap_list(x['unp_range'])), axis=1)
        dfrm['var_list'] = dfrm.apply(lambda x: json.dumps(get_range_diff(
            x['unp_range'], x['pdb_range'])), axis=1)
        dfrm['repeated'] = dfrm.apply(
            lambda x: '-' in x['var_list'], axis=1)
        dfrm['repeated'] = dfrm.apply(
            lambda x: True if '-' in x['unp_gap_list'] else x['repeated'], axis=1)
        dfrm['var_0_count'] = dfrm.apply(
            lambda x: json.loads(x['var_list']).count(0), axis=1)
        dfrm['unp_gap_0_count'] = dfrm.apply(
            lambda x: json.loads(x['unp_gap_list']).count(0), axis=1)
        dfrm['unp_pdb_var'] = dfrm.apply(
            lambda x: json.loads(x['var_list'])[0], axis=1)
        add_tage_to_range(dfrm, tage_name='sifts_range_tag')
        return dfrm

    @classmethod
    def update_range(cls, dfrm: pd.DataFrame, new_range_cols=('new_unp_range', 'new_pdb_range')) -> pd.DataFrame:
        focus_index = dfrm[
            (dfrm.sifts_range_tag.isin(('Deletion', 'Insertion & Deletion')))
            & (dfrm.repeated.eq(False))].index
        updated_pdb_range, updated_unp_range = list(), list()
        seqAligner = SeqPairwiseAlign()
        for index in focus_index:
            record = dfrm.loc[index]
            pdbSeq = cls.get_seqres(record['pdb_id'], record['entity_id'])
            unpSeq = cls.get_unp_seq(record["UniProt"])
            res = seqAligner.makeAlignment(unpSeq, pdbSeq)
            updated_unp_range.append(res[0])
            updated_pdb_range.append(res[1])

        updated_range_df = pd.DataFrame(
            {new_range_cols[0]: updated_unp_range, new_range_cols[1]: updated_pdb_range}, index=focus_index)
        dfrm = pd.merge(dfrm, updated_range_df, left_index=True,
                        right_index=True, how='left')
        dfrm[new_range_cols[0]] = dfrm.apply(lambda x: x['unp_range'] if pd.isna(
            x[new_range_cols[0]]) else x[new_range_cols[0]], axis=1)
        dfrm[new_range_cols[1]] = dfrm.apply(lambda x: x['pdb_range'] if pd.isna(
            x[new_range_cols[1]]) else x[new_range_cols[1]], axis=1)
        return dfrm

    @classmethod
    @lru_cache()
    def get_unp_seq(cls, unp: str, session=None):
        if session is None:
            session = cls.session
        query = '''
            MATCH (unp:UniProt{ACCESSION: "%s"})-[:HAS_UNP_RESIDUE]-(unpRes:UNPResidue)
            RETURN unpRes.ONE_LETTER_CODE AS aa ORDER BY toInteger(unpRes.ID)
        ''' % unp
        return ''.join(r['aa'] for r in session.run(query))
        
    @classmethod
    def summary_muta(cls, lyst, session=None):
        def get_mutaRes(dfrm: pd.DataFrame):
            sites = defaultdict(list)
            def storeSites(lyst, sitelyst):
                sitelyst.extend([i.split('|')[1] for i in lyst])
            dfrm.apply(lambda x: storeSites(x[2], sites[(x[0], x[1], x[3])]), axis=1)
            return sites

        def yield_dfrm_of_mutaRes(sites: Dict):
            for (pdb, entity, tag), value in sites.items():
                cur_df = pd.DataFrame(set(value), columns=['residue_number'])
                cur_df['pdb_id'] = pdb
                cur_df['entity_id'] = entity
                cur_df['tag'] = tag
                yield cur_df
        
        if session is None:
            session = cls.session
        muta_query = '''
            MATCH (entry:Entry)-[:HAS_ENTITY]->(entity:Entity{POLYMER_TYPE:'P'})-[:HAS_RESIDUE_CONFLICT]->(resCon:ResidueConflict)
            WHERE entry.ID in $lyst
            RETURN entry.ID as pdb_id, entity.ID as entity_id, resCon.DETAILS, resCon.ID
        '''
        result = cls.to_data_frame(session.run(muta_query, lyst=lyst))
        sites = get_mutaRes(result)
        mutaRes_df = pd.concat(yield_dfrm_of_mutaRes(sites), sort=False, ignore_index=True)
        focusRes = (mutaRes_df.pdb_id+'_'+mutaRes_df.entity_id+'_' + mutaRes_df.residue_number).drop_duplicates().to_list()
        query = """
                MATCH (entry:Entry)-[:HAS_ENTITY]->(entity:Entity{POLYMER_TYPE:'P'})-[:CONTAINS_CHAIN]->(chain:Chain)-[inChain:IS_IN_CHAIN]-(res:PDBResidue)
                WHERE res.UNIQID IN $focusRes AND entry.ID in $lyst
                RETURN entry.ID as pdb_id, entity.ID as entity_id, chain.AUTH_ASYM_ID as chain_id, res.ID as residue_number, inChain.OBSERVED as OBSERVED
            """
        result = cls.to_data_frame(session.run(query, focusRes=focusRes, lyst=lyst))
        mutaRes_df = pd.merge(mutaRes_df, result)
        return mutaRes_df


class SeqPairwiseAlign(object):
    def __init__(self):
        self.seqa = None
        self.seqb = None
        self.aligner = Align.PairwiseAligner()
        self.aligner.mode = 'global'
        self.aligner.open_gap_score = -10
        self.aligner.extend_gap_score = -0.5
        self.aligner.substitution_matrix = matlist.blosum62
        self.alignment_count = 0

    @lru_cache()
    def makeAlignment(self, seqa, seqb):
        if seqa is None or seqb is None:
            return np.nan, np.nan
        self.seqa = seqa
        self.seqb = seqb
        alignments = self.aligner.align(seqa, seqb)
        for alignment in alignments:
            result = self.getAlignmentSegment(alignment)
            self.alignment_count += 1
            return json.dumps(result[0]), json.dumps(result[1])

    @staticmethod
    def getAlignmentSegment(alignment):
        segments1 = []
        segments2 = []
        i1, i2 = alignment.path[0]
        for node in alignment.path[1:]:
            j1, j2 = node
            if j1 > i1 and j2 > i2:
                segment1 = (i1 + 1, j1)
                segment2 = (i2 + 1, j2)
                segments1.append(segment1)
                segments2.append(segment2)
            i1, i2 = j1, j2
        return segments1, segments2
