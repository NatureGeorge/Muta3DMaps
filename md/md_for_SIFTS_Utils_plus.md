# FILENAME: SIFTS_Unit.py
> Code By ZZF, 2019-08-19

## Related Files
* ```raw_sifts_file_path = 'pdb_uniprot_SIFTS_raw_demo.tsv'```
* ```add_rangeInfo_sifts_file_path = 'pdb_uniprot_SIFTS_addRangeInfo_demo.tsv'```
* ```add_InDe_sifts_file_path = 'pdb_uniprot_SIFTS_delwithInDe_demo.tsv'```
* ```unp_len_file_path = 'unp_len_list_demo.tsv'```
* ```new_seg_edited_file_path = 'uniprot_segments_observed_edited_demo.tsv'```
* ```add_unpLen_segInfo_file_path = 'pdb_uniprot_SIFTS_addunpLenSegInfo_demo.tsv'```

### Raw SIFTS File: ```raw_sifts_file```
* File Format: ```tsv```
* sep: ```\t```

#### Columns
pdb_id	|	chain_id	|	UniProt	|	identity	|	identifier | pdb_start	|	pdb_end	|	unp_start	|	unp_end | is_canonical	|	start	|	end	|	entity_id	|	struct_asym_id
-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-
str | str | str | float | str | int| int | int | int | bool | str | str | str | str

##### Note
* No ```NAN``` value is permitted
* Some values of ```chain_id```or```struct_asym_id``` may be ```'NA'``` or something else which will be recognized as ```NaN``` by ```pandas.read_csv()```. To solve this problem, there are 2 ways when parsing the file with ```pandas.read_csv()```:
  * Let ```keep_default_na=False``` and specified ```na_values=['', None]```
  * specified ```converters={'chain_id':str, 'struct_asym_id':str}```
> [Reference 1](https://stackoverflow.com/questions/16596188/pandas-convert-na-to-nan "stackoverflow.com"), [Reference 2](https://pandas.pydata.org/pandas-docs/stable/reference/api/pandas.read_csv.html "pandas.pydata.org")

---

#### Example
#### JSON-format results from SIFTS API
```json
{
  "2znj": {
    "UniProt": {
      "B0S4P3": {
        "identifier": "B0S4P3_DESHA",
        "name": "B0S4P3_DESHA",
        "mappings": [
          {
            "entity_id": 1,
            "chain_id": "A",
            "pdb_start": 21,
            "start": {...},
            "unp_end": 288,
            "pdb_end": 308,
            "identity": 1,
            "unp_start": 1,
            "end": {...},
            "is_canonical": true,
            "struct_asym_id": "A"
          },
          {...},
          {...}
        ]
      }
    }
  }
}
```


#### Results of ```SIFTS_Unit.get_raw_SIFTS()```
pdb_id	|	chain_id	|	UniProt	|	identity	|	identifier | pdb_start	|	pdb_end	|	unp_start	|	unp_end | is_canonical	|	start	|	end	|	entity_id	|	struct_asym_id
-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-
2znj	|	A	|	B0S4P3	|	1	|	B0S4P3_DESHA	|	21	|	308	|	1	|	288	|	True	|	{author_residue_number: null, author_insertion_code: "", residue_number: 21}	|	{author_residue_number: 288, author_insertion_code: "", residue_number: 308}	|	1	|	A
2znj	|	B	|	B0S4P3	|	1	|	B0S4P3_DESHA	|	21	|	308	|	1	|	288	|	True	|	{author_residue_number: null, author_insertion_code: "", residue_number: 21}	|	{author_residue_number: 288, author_insertion_code: "", residue_number: 308}	|	1	|	B

### SIFTS File With Range Info: ```add_rangeInfo_sifts_file```
* File Format: ```tsv```
* sep: ```\t```

#### Columns
> Difference with ```raw_sifts_file```
##### Add
sifts_pdb_range | sifts_unp_range
-|-
str (JSON-format) | str (JSON-format)
##### Del
pdb_start | pdb_end | unp_start | unp_end
-|-|-|-
str | str | str | str

##### Note
* No ```NAN``` value is permitted in these new columns
* ```pdb_id``` will be transform into capital letters

#### Example
pdb_id	|	chain_id	|	UniProt	| ... |is_canonical | ... | sifts_pdb_range | sifts_unp_range
-|-|-|-|-|-|-|-
2ZNJ	|A	|B0S4P3	| ... | True | ... |\[\[21, 308]]	|\[\[1, 288]]
2ZNJ	|B	|B0S4P3	| ... | True | ... |\[\[21, 308]]	|\[\[1, 288]]
1XZ9	|A	|P55196-1	| ... | False | ... |\[\[6, 52], \[54, 101]]	|\[\[985, 1031], \[1032, 1079]]
1XZ9	|A	|P55196-3	| ... | False | ... |\[\[6, 52], \[54, 101]]	|\[\[985, 1031], \[1032, 1079]]

### SIFTS File With Insertion & Deletion Info: ```add_InDe_sifts_file```
* File Format: ```tsv```
* sep: ```\t```

#### Columns
> Difference with ```add_rangeInfo_sifts_file```

##### Add
pdb_GAP_list | unp_GAP_list | var_list | delete | var_0_count | unp_GAP_0_count | group_info | sifts_unp_pdb_var | sifts_range_tage
-|-|-|-|-|-|-|-|-
str (JSON-format) | str (JSON-format) | str (JSON-format) | bool | int | int | int | int | str

##### Note
* No ```NAN``` value is permitted in these new columns

#### Example
pdb_id	|	chain_id	|	UniProt	| identity |... |is_canonical | ... | sifts_pdb_range | sifts_unp_range | pdb_GAP_list | unp_GAP_list | var_list | delete | var_0_count | unp_GAP_0_count | group_info | sifts_unp_pdb_var | sifts_range_tage
-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-
3WZU	|A	|O14733-4	|0.978	|...|False|...	|\[\[2, 318]]	|\[\[103, 426]]	|\[]	|\[]	|\[7]	|False	|0	|0	|1	|7	|Deletion
3WZU	|A	|O14733	|1.000	|...|True|...	|\[\[2, 318]]	|\[\[103, 419]]	|\[]	|\[]	|\[0]	|False	|1	|0	|1	|0	|Safe
1XZ9	|A	|P55196-1	|0.990 |	...|False|...	|\[\[6, 52], \[54, 101]]	|\[\[985, 1031], \[1032, 1079]]	|\[1]	|\[0]	|\[0, 0]	|False	|2	|1	|2	|0	|Insertion
1XZ9	|A	|P55196	|1.000	|...|True|...	|\[\[6, 101]]	|\[\[1001, 1096]]	|\[]	|\[]	|\[0]	|False	|1	|0	|1	|0	|Safe
2WTD	|A	|O96017	|0.997	|...|True|...	|\[\[7, 329]]	|\[\[209, 531]]	|\[]	|\[]	|\[0]	|False	|1	|0	|1	|0	|Safe
2WTD	|A	|O96017-2	|0.210 |	...|False|...	|\[\[1, 39], \[42, 55]]	|\[\[151, 196], \[197, 210]]	|\[2]	|\[0]	|\[7, 0]	|False	|1	|1	|2	|7	|Insertion & Deletion
3CMU	|A	|P0A7G6	|0.987	|...|True|...	|\[\[5, 309], \[325, 658], \[674, 1007], \[1023, 135...	|\[\[31, 335], \[2, 335], \[2, 335], \[2, 335], \[2, ...	|\[15, 15, 15, 12, 14]	|\[-334, -334, -334, -334, -334]	|\[0, 0, 0, 0, 0, 0]	|True	|6	|0	|6	|0	|Insertion & Deletion
5LOH	|A	|Q96GX5-2	|0.718 |...|False|...	|\[\[4, 197], \[202, 218], \[256, 341]]	|\[\[1, 194], \[739, 209], \[755, 840]]	|\[4, 37]	|\[544, 545]	|\[0, -546, 0]	|True	|2	|0	|3	|0	|Insertion & Deletion
5LOH	|A	|Q96GX5	|0.986	|...|True|...	|\[\[4, 197], \[202, 341]]	|\[\[1, 194], \[740, 879]]	|\[4]	|\[545]	|\[0, 0]	|False	|2	|0	|2	|0	|Insertion & Deletion

### SIFTS File With uniprot_length and SIFTS_segment Info: ```add_unpLen_segInfo_file```
* File Format: ```tsv```
* sep: ```\t```

#### Columns
> Difference with ```add_InDe_sifts_file```
##### Add
UNP_len | seg_pdb_range | seg_unp_range
-|-|-
int | str (JSON-format) | str (JSON-format)

##### Note
* No ```NAN``` value is permitted in these new columns

#### Example
pdb_id | chain_id | UniProt | identity |...| is_canonical |...| sifts_pdb_range | sifts_unp_range | ...|sifts_range_tage | UNP_len | seg_pdb_range| seg_unp_range
-|-|-|-|-|-|-|-|-|-|-|-|-|-
2ZNJ	|A	|B0S4P3	|1.000	|...|True|...	|\[\[21, 308]]	|[\[1, 288]]	|...|Safe	|288.0	|[\[31, 308]]	|[\[11, 288]]
