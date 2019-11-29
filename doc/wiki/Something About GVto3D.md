# Some thing About GVto3D

This is a personal summary of the following article.

> __Mapping genetic variations to three-dimensional protein structures to enhance variant interpretation: a proposed framework.__
Glusman G, Rose PW, ... [Link](https://genomemedicine.biomedcentral.com/articles/10.1186/s13073-017-0509-y "Link")

* `GVto3D`: Gene Variation to 3D (workshop)

## Written in front (Overall Summmary)

Introduce a workshop that proposed a framework/portal (GVto3D) aimed to promote the integration of some data resources. It is just a proposed portal of the workshop, not a tool or database.

## Goal

 > The overarching goal of the workshop was to address the question: what can be done together as a community to __advance the integration of genetic variants and 3D protein structures__ that could not be done by a single investigator or laboratory?

## Outcomes

 > ..propose the development of a __framework__ with which to promote progress in this arena. The framework will include a set of standard formats, common ontologies, a common application programming interface to enable interoperation of the resources, and a Tool Registry to make it easy to find and apply the tools to specific analysis problems. Interoperability will enable integration of diverse data sources and tools and collaborative development of variant effect prediction methods.

## Background: What has been done

Note: *Only A Brief Summary*

* Various Mutation Resources (SNVs e.g)
  * dbSNP, dbVar, COSMIC, cBioPortal, UniProt, Kaviar, Clinvar, HGMD, ExAC, and gnomAD
  * ...
* Various Protein Structure Resources
  * PDB, ~~Protein Model Portal~~__(Deprecated)__, I-TASSER，ModWeb, Phyre2，HHpred and SWISS-MODEL
  * ...
* Some Attempts To Find Connections Between Genetic Variations and Protein Structure
  * cBioPortal, COSMIC-3D, CRAVAT, Jalview, MuPIT, MutDB, STRUM, Cancer3D
  * need to do sequence alignments _(Note: some kind of extra work)_
* ...

## What will be done (from 2017)

In the context of what I concentrate on:

> Mapping between genomic and protein coordinate systems is error prone due to, for example, different genome assembly versions and alternative splicing. Where a mapping from genome to UniProt is possible, `SIFTS` and CRAVAT provide consistent residue-level mapping to and from PDB structures and other resources.

## Details of GVto3D

### Components of the GVto3D portal

> The Tools Registry contains a searchable description and metadata for tools, resources, and reference data sets for __third-party variant effect prediction and annotation services__. Standardized application programming interfaces (APIs) provide interoperability for data input and output of these third-party tools. Custom adapters can provide limited interoperability for tools that cannot adopt the API. __A mapping service provides bidirectional mappings from reference genome coordinates to UniProt protein positions and to Protein Data Bank (PDB) residue positions.__ The tools can use the mapping service to accept variant positions in any of the three coordinate systems. A beacon system enables queries about variant positions where three-dimensional (3D) structural information and annotation are available

![fig](https://media.springernature.com/full/springer-static/image/art%3A10.1186%2Fs13073-017-0509-y/MediaObjects/13073_2017_509_Fig1_HTML.gif?as=webp "Link")

#### What I Focus on

> __A mapping service provides bidirectional mappings from reference genome coordinates to UniProt protein positions and to Protein Data Bank (PDB) residue positions.__

1. bidirectional mappings from reference genome coordinates to UniProt protein positions __(Via UniProt ID Mapping ?)__
2. (bidirectional mappings) from UniProt protein positions to to Protein Data Bank (PDB) residue positions. __(Via SIFTS ?)__

## Any works have been done since this workshop(?)

> <https://www.ncbi.nlm.nih.gov/pubmed?linkname=pubmed_pubmed_citedin&from_uid=29254494>

### PhyreRisk (2019 June)

> __PhyreRisk: A Dynamic Web Application to Bridge Genomics, Proteomics and 3D Structural Data to Guide Interpretation of Human Genetic Variants.__
Ofoegbu TC, David A, Kelley LA, Mezulis S, Islam SA, Mersmann SF, Strömich L, Vakser IA, Houlston RS, Sternberg MJE.
> _J Mol Biol, 2019 Jun 14;431(13):2460-2466. doi: 10.1016/j.jmb.2019.04.043. Epub 2019 May 7. PMID: 31075275_

PhyreRisk is an open-access, publicly accessible web application for interactively bridging genomic, proteomic and structural data facilitating the mapping of human variants onto protein structures. A major advance over other tools for sequence-structure variant mapping is that PhyreRisk provides information on 20,214 human canonical proteins and an additional 22,271 alternative protein sequences (isoforms)...

Note: ***Very similar to my work.***

* [Link of `phyrerisk`](<http://phyrerisk.bc.ic.ac.uk/home> "Link")

#### Features

* From Ensembl, VCF, reference SNP ID and HGVS notations (Human Build GRCh37 and GRCh38) to PDB residue, Via [`VEP`](https://www.ensembl.org/info/docs/tools/vep/index.html "Link")

* From UniProt (Canonical & Isoform) to PDB, Via `?`
* Provides structural coverage (partial or complete) for 70% (14,035 of 20,214 __canonical__ proteins) of the human proteome, by storing 18,874 experimental structures and 84,818 pre-built models of __canonical__ proteins and their __isoforms__ generated using our in house `Phyre2`.

* PhyreRisk reports 55,732 experimentally, multi-validated protein __interactions__ from `IntAct` (as per UniProt filtering) and 24,260 experimental structures of protein complexes.

* Currently, 163,286 variants from UniProt Humsavar database have been curated. In the current version of PhyreRisk, we have chosen not to store all known human variants catalogued by other databases, such as ExAC and dbSNP.

![fig](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6597944/bin/gr2.jpg)

#### Limitation (Personal Summary)

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

### GS2 (2018 June)

> __G2S: a web-service for annotating genomic variants on 3D protein structures.__
Wang J1, Sheridan R2, Sumer SO2, Schultz N2,3, Xu D1, Gao J2.
> _Bioinformatics. 2018 Jun 1;34(11):1949-1950. doi: 10.1093/bioinformatics/bty047._

Note/Personal Conclusion: __Surpassed__ by `SIFTS(2019)` + `UniProt ID Mapping`

* [Link of `G2S`](<https://g2s.genomenexus.org/> "Link")

![fig](https://g2s.genomenexus.org/images/workflow.png)

## Some noteworthy things

> <http://www.rcsb.org/pdb/chromosome.do> (Since 2016)

> __Integrating genomic information with protein sequence and 3D atomic level structure at the RCSB protein data bank.__
Prlic A1, Kalro T1, Bhattacharya R2, Christie C1, Burley SK1,3, Rose PW1.
> _Bioinformatics. 2016 Dec 15;32(24):3833-3835. doi: 10.1093/bioinformatics/btw547. Epub 2016 Aug 22._

Note: My work has a strong similarity with the job that RCSB has done, however, their approach doesn't provide the users with enough accessibilities to raw data. For instance, the users can not upload a muta-site of UniProt Sequence and then retrieve all possible results via API

### Data Provided By RCSB

> <http://www.rcsb.org/pdb/browse/homo_sapiens.do>

#### Limitation

1. Only RefSeq, Ensembl Gene with Canonical UniProt
2. Not a bidirectional tool
3. Can not batch retrieve?
4. Isoform mutations are mapped to Canoncial mutations
5. No sequence identity info