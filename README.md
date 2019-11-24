# Muta3DMaps
> Date: 2019-08-16T23:07:24+08:00
A tool developed by Minghui Li Group and maintained by ZeFeng Zhu.

* 🔨 stands for __\[under construction]__

`Muta3DMaps` is Python Package that retrieves, filtering and organizes data from various databases or tools via API to process the residue-level mapping between protein sequences and protein 3D structures.

## Dependent Tool & DB
* SIFTS
* UniProt
* wwPDB
* Interactome3D
* SWISS-MODEL Repository 🔨
* ModBase 🔨




## Function
* Collect data from SIFTS, UniProt, wwPDB(MMCIF), Interactome3D, SMR, ModBase
* Define representative structures(PDB or Model) of a uniprot(Canonical) with a score-based approach 🔨
* Map mutation from (transcript/UniProt) isoform to (Canonical UniProt/PDB/Model) or from PDB  to UniProt Isoform.

## Regarded as Modules
> Correspondence Between Core Functions and ```PYTHON``` Files

* ```ProcessUniProt.py```
  * UniProt ID Mapping API
* ```ProcessSIFTS.py```
  * Map PDB to UniProt and vice versa
  * Map PDB Mutation Site to UniProt and vice versa
* ```RetrievePDB.py```
  * Download PDB Files
* ```ProcessMMCIF.py```
  * Extract Info From mmCIF-format PDB Files
* ```ProcessI3D.py```
  * Retrieve and modify interactions info of Interactome3D
* ```Representative.py``` 🔨
  * Generate Representative Structure Dataset
* ```ProcessSMR.py``` 🔨
  * Retrieve and modify SWISS-MODEL Repository Info
* ```ProcessModB.py``` 🔨
  * Retrieve and modify ModBase Info

## Regarded as a Command Line Tool

```bash
Usage: Run.py [OPTIONS] COMMAND [ARGS]...

Options:
  --folder PATH  The file folder of new files.
  --help         Show this message and exit.

Commands:
  i3dmap
  initmmcif
  initsifts
  inituniprot
  initunpfasta
  unp2pdb
```

## Install with `setuptools`

```bash
python setup.py install --record record.txt
```

### Uninstall

```bash
FOR /F "delims=" %f in (record.txt) DO del "%f"
```
