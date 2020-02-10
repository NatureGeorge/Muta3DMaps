import pytest
from Muta3DMaps.core.retrieve.fetchFiles import UnsyncFetch


UNP_ID_MAP = 'https://www.uniprot.org/uploadlists/'
PDBE_URL = 'https://www.ebi.ac.uk/pdbe/api'
SMR_URL = 'https://swissmodel.expasy.org/repository/uniprot'
MODB_URL = 'http://salilab.org/modbase-cgi/model_search.cgi'
I3D_URL = 'https://interactome3d.irbbarcelona.org/api'

UNP_ID_MAP_COLUMNS = [
    'id', 'length', 'reviewed', 'comment(ALTERNATIVE%20PRODUCTS)',
    'feature(ALTERNATIVE%20SEQUENCE)', 'genes','organism','protein%20names']
ID_DEMO = [
    'ENST00000335137', 'ENST00000420190', 'ENST00000434641',
    'ENST00000379407', 'ENST00000379409']

TASKS = [
    ('get', {'url': UNP_ID_MAP, 'params': {'from': 'ENSEMBL_TRS_ID', 'to': 'ACC', 'format': 'tab', 'query': ','.join(ID_DEMO), 'columns': ','.join(UNP_ID_MAP_COLUMNS)}}, 'uniprot_id_mapping_test.tsv'),
    ('get', {'url': f'{PDBE_URL}/pdb/entry/molecules/', 'params': '1a01'}, '1a01_molecules.json'),  # 405
    ('get', {'url': f'{PDBE_URL}/pdb/entry/molecules/1a01'}, '1a01_molecules.json'),
    ('get', {'url': f'{PDBE_URL}/pdb/entry/molecules/2xyn'}, '2xyn_molecules.json'),
    ('get', {'url': f'{PDBE_URL}/pdb/entry/molecules/1miu'}, '1miu_molecules.json'),
    ('get', {'url': f'{PDBE_URL}/pdb/entry/molecules/2hev'}, '2hev_molecules.json'),
    ('get', {'url': f'{PDBE_URL}/pdb/entry/molecules/3g96'}, '3g96_molecules.json'),
    ('get', {'url': f'{PDBE_URL}/pdb/entry/molecules/6lu7'}, '6lu7_molecules.json'),
    ('post', {'url': f'{PDBE_URL}/pdb/entry/summary/', 'data': '6lu7,3g96,2hev'}, '6lu7,3g96,2hev_summary.json'),
    ('post', {'url': f'{PDBE_URL}/pdb/entry/summary/', 'data': '1a01,2xyn,1miu'}, '1a01,2xyn,1miu_summary.json'),
    ('get', {'url': f'{PDBE_URL}/pdb/entry/residue_listing/1a01'}, '1a01_residue_listing.json'),
    ('post', {'url': f'{PDBE_URL}/pdb/entry/status/', 'data': '1a01,pppp,2xyn'}, '1a01,pppp,2xyn_status.json'),  # contail 404, the server would ignore pppp
    ('get', {'url': f'{PDBE_URL}/mappings/all_isoforms/1a01'}, '1a01_sifts.json'),
    ('get', {'url': f'{SMR_URL}/P07900.json', 'params': {'provider': 'swissmodel'}}, 'P07900_SMR.json'),
    ('get', {'url': f'{SMR_URL}/P0DP24.json', 'params': {'provider': 'swissmodel'}}, 'P0DP24_SMR.json'),
    ('get', {'url': f'{SMR_URL}/P07900-2.json', 'params': {'provider': 'swissmodel'}}, 'P07900-2_SMR.json'),
    ('get', {'url': f'{SMR_URL}/P07900-2.pdb', 'params': {'provider': 'swissmodel', 'sort': 'seqid'}}, 'P07900-2_SMR.pdb'),
    ('get', {'url': f'{I3D_URL}/getInteractionStructures', 'params': {'queryProt1': 'A0A5B9', 'queryProt2': 'P01848'}}, 'A0A5B9-P01848_I3D.xml'),
    ('get', {'url': f'{I3D_URL}/getInteractionStructures', 'params': {'queryProt1': 'Q16543', 'queryProt2': 'P07900'}}, 'Q16543-P07900_I3D.xml'),
    ('get', {'url': f'{I3D_URL}/getPdbFile', 'params': {'filename': 'P07900-Q16543-EXP-pdb2k5b.ent-A-0-B-0.pdb', 'type': 'interaction'}}, 'P07900-Q16543_I3D.pdb'),
    ('get', {'url': f'{MODB_URL}', 'params': {'searchkw': 'name', 'kword': 'Q99777'}}, 'ModB_test.html'),
]
def task():
    try:
        return len(UnsyncFetch.main(r'./data/', TASKS, 8))
    except Exception:
        return -1

def test_task():
    res = task()
    assert res == len(TASKS)  # and len([i for i in res if i is None]) == 2
