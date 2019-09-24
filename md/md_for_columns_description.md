# Description For Columns
> Last Modified Time: 20190924

## UniProt Related Columns

```json
{
  "Entry": "UniProt Entry",
  "Gene names": "A list of this Entry's gene names",
  "Status": "Reviewed or Unreciewed",
  "Alternative products (isoforms)": "Details of the difference between isoforms",
  "Organism": "Species",
  "GENE": "Query Gene Name",
  "refSeq_protein": "Quert RefSeq Protein Accession Number",
  "UniProt": "Isoform ID",
  "canonical_isoform": "Isoform of Canonical Sequence",
}
```

## SIFTS Related Columns

```json
{
  "pdb_id": "Upper Form of PDB ID",
  "chain_id": "pdb_strand_id",
  "identity": "Sequence Identity between PDB_SEQRES and UniProt Isoform Sequence",
  "identifier": "GeneName_Species",
  "is_canonical": "boolean",
  "start": "Details of mapped start site",
  "end": "Details of mapped end site",
  "entity_id": "same as entity_id in mmCIF file",
  "struct_asym_id": "same as asym_id in mmCIF file",
  "sifts_pdb_range": "mapped range of pdb_chain",
  "sifts_unp_range": "mapped range of uniprot_isoform",
  "pdb_GAP_list": "Gap between ranges of sifts_pdb_range",
  "unp_GAP_list": "Gap between ranges of sifts_unp_range",
  "var_list": "Difference between ranges of sifts_pdb_range and that of sifts_unp_range",
  "delete": "whether delete or not",
  "var_0_count": "count 0 in var_list",
  "unp_GAP_0_count": "count 0 in unp_GAP_list",
  "group_info": "",
  "sifts_unp_pdb_var": "",
  "sifts_range_tage": "Safe or Deletion or Insertion or Insertion&Deletion",
  "UNP_len": "the length of UniProt Isoform Sequence",
  "seg_pdb_range": "",
  "seg_unp_range": "",
}

```

## MMCIF Related Columns

```json
{
  "protein_type": "polypeptide(L) or polypeptide(D) or ...",
  "_pdbx_poly_seq_scheme.mon_id": "SEQRES",
  "_pdbx_poly_seq_scheme.pdb_mon_id": "SEQRES: MISSING->?",
  "_pdbx_poly_seq_scheme.auth_mon_id": "SEQRES: MISSING->?",
  "_pdbx_poly_seq_scheme.ndb_seq_num": "SEQRES INDEX: Start at 1",
  "_pdbx_poly_seq_scheme.pdb_seq_num": "SEQRES INDEX",
  "_pdbx_poly_seq_scheme.auth_seq_num": "SEQRES INDEX: MISSING->?",
  "_pdbx_poly_seq_scheme.pdb_ins_code": "SEQRES INSIDE CODE",
  "asym_id": "",
  "UNK_ALL_IN_CHAIN": "Whether all the chains contains UNK",
  "method": "Resolution Method",
  "initial_version_time": "",
  "newest_version_time": "",
  "resolution": "",
  "pdb_contain_chain_type": "",
  "contains_unk_in_chain_pdb": "Whether the PDB contains chains that contain UNK",
  "pdb_type_MMCIF": "mo or ho:num or he:chain_group1;chain_group2'...",
  "bioUnit": "{'bioUnit_1':[num, chain_list], 'bioUnit_2':...}",
  "_em_3d_reconstruction.resolution": "ELECTRON MICROSCOPY Resolution ...",
  "_refine.ls_d_res_high": "X-RAY Resolution ...",
  "_pdbx_struct_assembly_gen.assembly_id": "bioUnit ID",
  "_pdbx_struct_assembly_gen.oper_expression": "",
  "_pdbx_struct_assembly_gen.asym_id_list": "",
  "_pdbx_struct_assembly.oligomeric_count": "",
  "_pdbx_entity_nonpoly.entity_id": "",
  "_pdbx_entity_nonpoly.name": "",
  "_pdbx_entity_nonpoly.comp_id": "",
  "metal_ligand_content": "",
  "_pdbx_coordinate_model.type": "The Info of ATOM-ONLY",
  "mutation_content": "",
  "mutation_num": "",
  "metal_ligand_num": "",
  "Modification_num": "",
  "seqres_len": "",
  "coordinates_len": "Length of ATOM Residues",
  "Modification_index": "Index of Modified Residues(Start at 1)",
  "mis_index": "Index of Missing Residues(Start at 1)",
  "mis_range": "",
  "resolution_score": "resolution, if nan, 1000",
  "protein_chain_and_length": "A list of Length of all Protein chains",
}

```

## Other New Added Columns

```json
{
  "new_sifts_unp_range": "New Range That Should use",
  "new_sifts_pdb_range": "New Range That Should use",
  "mutation_unp": "Mutation list of uniprot",
  "mutation_pdb": "MutaSite list of PDB_Chain",
  "muta_map_info": "Mutation Mapping Info List",
  "pdb_mapped_range",
  "pdb_mappedRange_head",
  "pdb_mappedRange_tail",
  "mappedOut",
  "delHT_MissingNum",
  "pdb_mapped_range_len",
  "if_1",
  "if_2",
  "BS": "Score",
  "ne_resolution_score": "Negative Resolution Score",
  "pdb_SIFTS_useful_chain_num",
  "pdb_useful_chain_num",
  "ab",
  "pdb_chain_select": "Whether select or not",
  "pdb_chain_rank": "Rank of BS - ne_resolution_score - initial_version_time"
}
```

## Columns That You Should Focus On

|col|description|
-|-|-
|GENE|Query Gene Name|
|refSeq_protein|Quert RefSeq Protein Accession Number|
|UniProt|Isoform ID|
|pdb_id|Upper Form of PDB ID|
|chain_id|pdb_strand_id|
|identity|Sequence Identity between PDB_SEQRES and UniProt Isoform Sequence|
|sifts_range_tage|```Safe``` or ```Deletion``` or ```Insertion``` or ```Insertion & Deletion```|
|UNK_ALL_IN_CHAIN|Whether all the chains contains UNK|
|contains_unk_in_chain_pdb|Whether the PDB contains chains that contain UNK|
|method|Resolution Method|
|ne_resolution_score|Negative Resolution Score|
|pdb_contain_chain_type|Protein, DNA, RNA|
|pdb_type_MMCIF|mo or ho:num or he:chain_group1;chain_group2'...|
|bioUnit|{'bioUnit_1':[num, chain_list], 'bioUnit_2':...}|
|new_sifts_unp_range|New Range That Should use|
|new_sifts_pdb_range|New Range That Should use|
|mutation_unp|Mutation list of uniprot|
|mutation_pdb|MutaSite list of PDB_Chain|
|muta_map_info|Mutation Mapping Info List|
|pdb_chain_select|Whether select or not"|
