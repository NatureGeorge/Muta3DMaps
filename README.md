# Muta3DMaps

> Date: 2019-08-16T23:07:24+08:00
A tool developed by Minghui Li Group and maintained by ZeFeng Zhu.

* 🔨 stands for __\[under construction]__

`Muta3DMaps` is Python Package that retrieves, filtering and organizes data from various databases or tools via API to process the residue-level mapping between protein sequences and protein 3D structures.

[toc]

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

```py
# Example
>> from Muta3DMaps.core.Mods import ProcessUniProt, RetrievePDB
>> help(RetrievePDB)
```

```txt
Help on module Muta3DMaps.core.Mods.RetrievePDB in Muta3DMaps.core.Mods:

NAME
    Muta3DMaps.core.Mods.RetrievePDB

DESCRIPTION
    …

CLASSES
    builtins.object
        MPWrapper
        RetrievePDB

    class MPWrapper(builtins.object)
     |  MPWrapper(downloadPath, loggingPath, processes=3, maxSleep=3, ftpSite='RCSB', format='mmCIF')
     |
     |  Multiprocessing wrapper for ``RetrievePDB``
     |
     |  When there is a large number of PDB files to download, this class is helpful.
     |  But Need to be careful with the numbers of processes and the time of sleep.
…
```

## Regarded as a Command Line Tool

```bash
Usage: Muta3DMaps [OPTIONS] COMMAND [ARGS]...

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

### `inituniprot`

```bash
>Muta3DMaps --folder ./ inituniprot --help
Usage: Muta3DMaps inituniprot [OPTIONS]

Options:
  --referenceFile PATH      The reference file of IDs(with mutation Site) that
                            need to map via UniProt RESTful API.
  --sep TEXT                The seperator of referenceFile.
  --idCol TEXT              The column name of IDs in referenceFile.
  --idType TEXT             ID Abbreviation that stands for the type of ID.
  --addUseCols TEXT         Comma-separated list of the column names for
                            programmatic access to the UniProtKB search
                            results.
  --siteCol TEXT            The column name of aa site in referenceFile.
  --geneCol TEXT            The column name of gene info in referenceFile.
  --procced / --no-procced  Whether to procced after saving the site info.
  --help                    Show this message and exit.
```

### `initunpfasta`

```bash
>Muta3DMaps --folder ./ initunpfasta --help
Usage: Muta3DMaps initunpfasta [OPTIONS]

Options:
  --fastaFolder PATH    The file folder of UniProt FASTA Seq repository.
  --unreviewed BOOLEAN  Whethter to include FASTA Seq of unreviewed UniProt
                        Entry.
  --isoform BOOLEAN     Whethter to include isoform Seq.
  --split BOOLEAN       Whethter to split FASTA files.
  --mode [wget|ftplib]  Retrieve mode.
  --fastaPath PATH      The file path of downloaded fasta file.
  --referenceFile PATH  The file path of reference file that contains target
                        UniProt ID.
  --sep TEXT            The seperator of referenceFile.
  --colName TEXT        The column name of UniProt IDs in referenceFile.
  --help                Show this message and exit.
```

### `initsifts`

```bash
>Muta3DMaps --folder ./ initsifts --help
Usage: Muta3DMaps initsifts [OPTIONS]

Options:
  --test INTEGER              Num of PDB IDs to test the program. Only for
                              test.
  --unpFile PATH              The file that comtains Target UniProt IDs.
  --unpCol TEXT               The column of UniProt IDs in unpFile.
  --sep TEXT                  The seperator of unpFile.
  --filtering <TEXT TEXT>...  [filterColumn filterValue]: The filter of
                              unpFile. Keep the rows that have equal value in
                              filter column.
  --useInitizedUnp BOOLEAN    Whether to set the initialized result as the
                              unpFile.
  --help                      Show this message and exit.
```

### `initmmcif`

```bash
>Muta3DMaps --folder ./ initmmcif --help
Usage: Muta3DMaps initmmcif [OPTIONS]

Options:
  --pdbFolder PATH  The file folder of PDB repository.
  --pdbsFile PATH   The file that comtains PDB IDs.
  --pdbCol TEXT     The column of PDB IDs in pdbsFile.
  --sep TEXT        The seperator of pdbsFile.
  --help            Show this message and exit.
```

### `unp2pdb`

```bash
>Muta3DMaps --folder ./ unp2pdb --help
Usage: Muta3DMaps unp2pdb [OPTIONS]

Options:
  --fastaFolder PATH   The file folder of UniProt FASTA Seq repository.
  --siteInfoFile PATH  The file that comtains site info.
  --help               Show this message and exit.
```

### `i3dmap`

```bash
>Muta3DMaps --folder ./ i3dmap --help
Usage: Muta3DMaps i3dmap [OPTIONS]

Options:
  --i3dPath PATH  The downloaded file path of Interactome3D Meta file.
  --help          Show this message and exit.
```

## Install with `setuptools`

```bash
git clone https://github.com/NatureGeorge/Muta3DMaps
python setup.py install --record record.txt
```

### Some bugs about Python 3.7

if you encounter this:

```bash
AttributeError: type object 'Callable' has no attribute '_abc_registry'
```

when you use `Muta3DMaps`, you can either downgrad your Python to 3.6 or `pip uninstall typing`

> See [here](https://stackoverflow.com/questions/55833509/attributeerror-type-object-callable-has-no-attribute-abc-registry "link") for more information.

### Something wrong about `biopython 1.75`

For the reason that `biopython 1.75` has changed the way to set sub-matrix in `Bio.Align` and makes it becomes an unstable module to align sequences, `1.73, 1.74` versions are recommended.

### Uninstall

```bash
FOR /F "delims=" %f in (record.txt) DO del "%f"
```

## Install with `pip` (Need to be fixed)

```bash
pip install -i https://test.pypi.org/simple/ Muta3DMaps
```
