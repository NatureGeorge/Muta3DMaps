import pytest
from Muta3DMaps.core.retrieve.fetchFiles import UnsyncFetch


BASE_URL = 'https://www.ebi.ac.uk/pdbe/api/pdb/entry'

TASKS = [
    ('get', {'url': f'{BASE_URL}/molecules/', 'params': '1a01'}, '1a01_molecules.json'), # 405
    ('get', {'url': f'{BASE_URL}/molecules/1a01'}, '1a01_molecules.json'),
    ('get', {'url': f'{BASE_URL}/molecules/2xyn'}, '2xyn_molecules.json'),
    ('get', {'url': f'{BASE_URL}/molecules/1miu'}, '1miu_molecules.json'),
    ('get', {'url': f'{BASE_URL}/molecules/2hev'}, '2hev_molecules.json'),
    ('get', {'url': f'{BASE_URL}/molecules/3g96'}, '3g96_molecules.json'),
    ('get', {'url': f'{BASE_URL}/molecules/6lu7'}, '6lu7_molecules.json'),
    ('post', {'url': f'{BASE_URL}/summary/', 'data': '6lu7,3g96,2hev'}, '6lu7,3g96,2hev_summary.json'),
    ('post', {'url': f'{BASE_URL}/summary/', 'data': '1a01,2xyn,1miu'}, '1a01,2xyn,1miu_summary.json'),
    ('get', {'url': f'{BASE_URL}/residue_listing/1a01'}, '1a01_residue_listing.json'),
    ('post', {'url': f'{BASE_URL}/status/', 'data': '1a01,pppp,2xyn'}, '1a01,pppp,2xyn_status.json') # contail 404, the server would ignore pppp
]

def test_main():
    try:
        return UnsyncFetch.main(r'./data/', TASKS, 8)
    except Exception:
        return -1

def answer_main():
    res = test_main()
    assert res == len(TASKS)  # and len([i for i in res if i is None]) == 2