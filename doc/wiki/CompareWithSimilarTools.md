# Compare with some similar tools

## RCSB: Map Genomic Position to Protein Sequence and 3D Structure (2016)

### Reference

> __Integrating genomic information with protein sequence and 3D atomic level structure at the RCSB protein data bank.__
Prlic A1, Kalro T1, Bhattacharya R2, Christie C1, Burley SK1,3, Rose PW1.
> _Bioinformatics. 2016 Dec 15;32(24):3833-3835. doi: 10.1093/bioinformatics/btw547. Epub 2016 Aug 22._

#### Link

* <http://www.rcsb.org/pdb/chromosome.do>
* <http://www.rcsb.org/pdb/browse/homo_sapiens.do>

### Limitation

1. Can only input chromosome position
2. Only RefSeq, Ensembl Gene with Canonical UniProt
3. Isoform mutations are mapped to Canoncial mutations
4. Not a bidirectional tool?
5. Can not batch retrieve
6. No sequence identity info

### Question

1. whether mapping mutations of isoform `UniProt` to canonical sequence is trusted?
 * detect those mutations located in the missing/intertion/deletion/varience sites and do some statistic

2. how to get chromosome position so that I can make comparison between it and my tool

3. ...



## GS2 (2018 June)

![fig](https://g2s.genomenexus.org/images/workflow.png)

### Reference

> __G2S: a web-service for annotating genomic variants on 3D protein structures.__
Wang J1, Sheridan R2, Sumer SO2, Schultz N2,3, Xu D1, Gao J2.
> _Bioinformatics. 2018 Jun 1;34(11):1949-1950. doi: 10.1093/bioinformatics/bty047._

#### Link

* <https://g2s.genomenexus.org/>

### Limitation

1. Only RefSeq, Ensembl Gene
2. __Surpassed__ by `SIFTS(2019)` + `UniProt ID Mapping`
3. No model source

## PhyreRisk (2019 June)

![fig](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6597944/bin/gr2.jpg)

### Reference

> __PhyreRisk: A Dynamic Web Application to Bridge Genomics, Proteomics and 3D Structural Data to Guide Interpretation of Human Genetic Variants.__
Ofoegbu TC, David A, Kelley LA, Mezulis S, Islam SA, Mersmann SF, Strömich L, Vakser IA, Houlston RS, Sternberg MJE.
> _J Mol Biol, 2019 Jun 14;431(13):2460-2466. doi: 10.1016/j.jmb.2019.04.043. Epub 2019 May 7. PMID: 31075275_

#### Link

* <http://phyrerisk.bc.ic.ac.uk/>

### Limitation

1. Only Canonical UniProt are attached with PDB structures
    * For some transcripts that mapped to UniProt Isoforms, related Crystal Structures will be lost 
    * Furthermore, without sequence identity info (Between PDB Chain and UniProt Canonical/Isoform Sequence)

2. Not a bidirectional mapping between PDB residues and UniProt seq-site

3. Some parts of it are based on some databases, need to keep updating (`UniProt` e.g)

4. Limited Protein Structure Model resources
   * Only `Phyre2` and has no interface to `Interactome3D, ModBase, SWISS-MODEL Repository`
  
5. Without some info contained in MMCIF
   * Coordinates Length
   * Ligand Info
   * Without info of EntityMutation in PDB Chain
   * UNK residues in PDB Chain
   * CA-ONLY Chain in PDB
   * There is no judgment on the polymer form and whether it contains nucleic acid chain
   * ...

6. Only wild type structure

7. Without info of Insertion or Deletion between PDB Chain and UniProt

8. No interface/API to retrieve queried results (except for `.pdb` files)
   * Furthermore, results are split into pages for each query mutation
   * In concluding, low data accessibility

9. Without feedback on unmapped query
