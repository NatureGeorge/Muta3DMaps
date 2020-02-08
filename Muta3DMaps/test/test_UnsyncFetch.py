import pytest
from Muta3DMaps.core.retrieve.fetchFiles import UnsyncFetch

DEMO = [
    ('https://www.ebi.ac.uk/pdbe/api/pdb/entry/molecules/1a01',
     '1a01_molecules.json'),
    ('https://www.ebi.ac.uk/pdbe/api/pdb/entry/molecules/2xyn',
     '2xyn_molecules.json'),
    ('https://www.ebi.ac.uk/pdbe/api/pdb/entry/molecules/1miu',
     '1miu_molecules.json'),
    ('https://www.ebi.ac.uk/pdbe/api/pdb/entry/molecules/2hev',
     '2hev_molecules.json'),
    ('https://www.ebi.ac.uk/pdbe/api/pdb/entry/molecules/3g96',
     '3g96_molecules.json'),
    ('https://www.ebi.ac.uk/pdbe/api/pdb/entry/residue_listing/1a01',
     '1a01_residue_listing.json'),
    ('https://www.ebi.ac.uk/pdbe/api/pdb/entry/residue_listing/2xyn',
     '2xyn_residue_listing.json'),
    ('https://www.ebi.ac.uk/pdbe/api/pdb/entry/residue_listing/1miu',
     '1miu_residue_listing.json'),
    ('https://www.ebi.ac.uk/pdbe/api/pdb/entry/residue_listing/2hev',
     '2hev_residue_listing.json'),
    ('https://www.ebi.ac.uk/pdbe/api/pdb/entry/residue_listing/3g96',
     '3g96_residue_listing.json')]

def test_main():
    try:
        return UnsyncFetch.main(r'./data/', DEMO)
    except Exception:
        return -1

def answer_main():
    assert test_main() == len(DEMO)