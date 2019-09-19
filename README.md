# SIFTS_Plus_Muta_Maps
> Last modified time: 2019-08-16T23:07:24+08:00
A tool developed by Minghui Li Group.

## Dependent Tool & DB
* SIFTS
* UniProt
* PDB(MMCIF)
* SWISS-MODEL Repository
* ModBase
* Interactome3D

## Function
* Collect data from SIFTS, UniProt, PDB(MMCIF), SMR, ModBase, Interactome3D
* Define representative structures(PDB or Model) of a uniprot(Canonical) with a score-based approach
* Map mutation from (transcript/UniProt) isoform to (Canonical UniProt/PDB/Model) or from PDB  to UniProt Isoform.

## Correspondence Between Function and ```PYTHON``` Files
* ```MMCIFplus.py```__\[CURRENTLY WORKING ON]__
  * Extract Info From mmCIF-format PDB Files
* ```SIFTS_units.py```__\[CURRENTLY WORKING ON]__
  * Map PDB to UniProt and vice versa
  * Map PDB Mutation Site to UniProt and vice versa
* ```Interactome3D_unit.py``` __\[NEED TO BE FIXED]__
  * Get Interactome3D Info
* ```SMR_unit.py``` __\[NEED TO BE FIXED]__
  * Get SWISS-MODEL Repository Info
* ```UniProt_unit.py``` __\[FINISHED]__
  * UniProt ID Mapping API

<img src="./docs/figs/code_fig_1.bmp"></img>
</br>
<img src="./docs/figs/code_fig_2.bmp"></img>
