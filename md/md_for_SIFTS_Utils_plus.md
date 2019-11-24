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
-|-|-|-|-|-|-|-|-|-|-|-|-|-
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
-|-|-|-|-|-|-|-|-|-|-|-|-|-
2znj	|	A	|	B0S4P3	|	1	|	B0S4P3_DESHA	|	21	|	308	|	1	|	288	|	True	|	{author_residue_number: null, author_insertion_code: "", residue_number: 21}	|	{author_residue_number: 288, author_insertion_code: "", residue_number: 308}	|	1	|	A
2znj	|	B	|	B0S4P3	|	1	|	B0S4P3_DESHA	|	21	|	308	|	1	|	288	|	True	|	{author_residue_number: null, author_insertion_code: "", residue_number: 21}	|	{author_residue_number: 288, author_insertion_code: "", residue_number: 308}	|	1	|	B

### self.get_raw_SIFTS()

#### Parameters

```py
def get_raw_SIFTS(self, outputPath):[...]
```

#### Usage
Get raw SIFTS data via SIFTS API: [SIFTS Mappings (PDB <-> UniProt all isoforms)](http://www.ebi.ac.uk/pdbe/api/doc/sifts.html "Web Link").
* Get data of ```pdb_id``` in ```self.pdb_list``` which is set by ```self.set_lists(pdb_list, [])```.
* The JSON-format data will be transformed into a dataframe.
* The dataframe will be saved in ```outputPath```.
* The ```outputPath``` is save by ```self.raw_SIFTS_filePath = outputPath```, and can be use in ```self.handle_SIFTS()```

#### Columns in DataFrame
Set by ```SIFTS_unit.CONFIG['RAW_SIFTS_COLUMNS']```

```py
[
        'pdb_id', 'chain_id', 'UniProt', 'identity', 'identifier',
        'pdb_start', 'pdb_end', 'unp_start', 'unp_end',
        'is_canonical', 'start', 'end', 'entity_id', 'struct_asym_id'
]
```

#### Return
This function returns a ```pdb_list``` that did not get data successfully.

```py
return fail_list
```

---


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

### self.handle_SIFTS()

#### Parameters

```py
def handle_SIFTS(self, sifts_filePath=False, outputPath=False):[...]
```

#### Usage
Handle the range of pdb and unp (_the info convey by ```['pdb_start', 'pdb_end', 'unp_start', 'unp_end']```, belongs to a particular PDB chain with corresponding UniProt isoform_) in the ___raw SIFTS file___ and then transform the range info into a interval-format.
* Create two new columns: ```sifts_pdb_range```  and ```sifts_unp_range```.
* With no specified ```sifts_filePath```, the function will take ```self.raw_SIFTS_filePath``` as the path of ___raw SIFTS file___.
* With specified ```outputPath```, the function will save the new DataFrame in it.
>  ___raw SIFTS file___ : File created by ```self.get_raw_SIFTS()```

#### Return
This function returns the new DataFrame.

```py
return new_sifts_df
```

---


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

### self.deal_with_insertionDeletion_SIFTS()

#### Parameters

```py
def deal_with_insertionDeletion_SIFTS(self, sifts_df=False, sifts_filePath=False, outputPath=False):[...]
```

#### Usage
Find possible Insertion and Deletion of a particular PDB chain, and the chain belongs to a corresponding UniProt isoform.
* Create New columns: ```['pdb_GAP_list', 'unp_GAP_list', 'var_list', 'delete', 'var_0_count', 'unp_GAP_0_count', 'group_info', 'sifts_unp_pdb_var', 'sifts_range_tage']```
* Value of column ```sifts_range_tage``` can be:
    * 'Safe': 'No Insertion & No Deletion'
    * 'Deletion': 'No Insertion But Deletion'
    * 'Insertion': 'Insertion, No Deletion'
    * 'Insertion & Deletion': 'Both Insertion and Deletion have been found'
* Value of column ```delete``` can be:
    * False: 'This correspondence of PDB chain and UniProt isoform is correct.'
    * True: 'This correspondence of PDB chain and UniProt isoform contains error.'
* With specified ```sifts_df```, the function will add columns in this DataFrame. But when specified ```sifts_filePath```, the function will open the file then add columns in it and ignore ```sifts_df```.


#### Return
This function returns the new DataFrame.

```py
return dfrm
```

---


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

---

### self.score_SIFTS()

#### Parameters

```py
def score_SIFTS(self, sifts_df=False, sifts_filePath=False, outputPath=False):[...]
```

#### Usage
Create columns that help examine the quality of pdb chain with corresponding uniprot and calculate the score.
* Create New columns: ```['pdb_mapped_range', 'pdb_mappedRange_head', 'pdb_mappedRange_tail', 'mappedOut', 'delHT_MissingNum', 'pdb_mapped_range_len', 'if_1', 'if_2', 'BS','ne_resolution_score', 'pdb_useful_chain_num']```
  * ```pdb_mapped_range```: interval-format data that represent the mapped range between a pdb chain sequence that have coordinates with its corresponding uniprot sequence
  * ```pdb_mappedRange_head```: the index of the first mapped residue in ```pdb_mapped_range```
  * ```pdb_mappedRange_tail```: the index of the last mapped residue in ```pdb_mapped_range```
  * ```mappedOut```: the length of those unmapped pdb sequences in the head and tail of a pdb chain sequence. Counts the length when the length is > 5 either in head or tail.
  * ```delHT_MissingNum```: the counts of missing residues in a pdb chain. Ignore those missing which locate in the head(<=5) and tail(<=5)
  * ```pdb_mapped_range_len```: the length of ```pdb_mapped_range```
  * ```if_1```: the sum of ```['pdb_GAP_list','unp_GAP_list','var_list']```, which represents the num of Insertion and Deletion related residues.
  * ```if_2```: the sum of ```['Modification_num', 'mutation_num']```, which represents the num of Modification and Mutatuin related residues
  * ```BS```: the basic score of a pdb chain, calculate by ```["coordinates_len", "mappedOut", "metal_ligand_num", "delHT_MissingNum", "if_2", "if_1", "UNP_len"]```
  * ```ne_resolution_score```: negative resolution_score, for those pdb without resolution data, we set it with -1000
  * ```pdb_useful_chain_num```: the num of protein chains that have coordinates_len more than 20.


  #### Return
  This function returns the new DataFrame.

  ```py
  return sifts_dfrm
  ```

---

### self.select_PDB_SIFTS()

#### Parameters

```py
def select_PDB_SIFTS(self, groupby_list, sifts_df=False, sifts_filePath=False, outputPath=False):[...]
```

#### Usage
Select the represent set of the groupby_list, composed of pdb chains
* ```groupby_list``` can be ```['UniProt']``` or ```['UniProt', 'pdb_type_MMCIF']``` or ```['UniProt', 'interaction_compo']``` .etc
