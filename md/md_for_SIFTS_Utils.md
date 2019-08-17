# FILENAME: SIFTS_Utils.py
> Code By ZZF, 2019-08-13

## Requirement
1. SIFTS_unit 能够补充新更新的PDB
2. SIFTS_unit 各列缺失值的定义与优化处理
3. SIFTS_unit 清晰的文件路径配置
4. MMCIF_unit 功能函数优化
5. MMCIF_unit 各列缺失值的定义与优化处理
5. MMCIF_unit 清晰的文件路径配置

## SIFTS_unit 流程
### Build Complete Basic DataSet
```py
import sys
import pandas as pd
sys.path.append('./')
from SIFTS_unit import SIFTS_unit
from UniProt_unit import UniProt_unit

raw_sifts_file_path = '../../data/pdb_uniprot_SIFTS_0714.csv'
add_rangeInfo_sifts_file_path = '../../data/pdb_uniprot_SIFTS_NEW0815.tsv'
add_InDe_sifts_file_path = '../../data/Mapping_Pipeline/sifts_files/pdb_uniprot_SIFTS_delwithInDe_0815.tsv'
unp_len_file_path = '../../data/Mapping_Pipeline/unp_len_list_0817.tsv'
seg_edited_file_path = '../../data/Mapping_Pipeline/sifts_files/uniprot_segments_observed_edited_0815.tsv'
add_unpLen_segInfo_file_path = '../../data/Mapping_Pipeline/sifts_files/???'

sifts_demo = SIFTS_unit()
info_dict = sifts_demo.get_info_from_uniprot_pdb_file() # Download a new file
sifts_demo.set_lists(info_dict['pdb_set'], [])

fail_list = sifts_demo.get_raw_SIFTS(raw_sifts_file_path)
sifts_df_1 = sifts_demo.handle_SIFTS(outputPath=add_rangeInfo_sifts_file_path)
sifts_df_2 = sifts_demo.deal_with_insertionDeletion_SIFTS(sifts_df=sifts_df_1, outputPath=add_InDe_sifts_file_path)

# Add UniProt Length Info
uniprot_demo = UniProt_unit()
uniprot_demo.get_info_from_uniprot(['id', 'length'], unp_len_file_path, unp_list=list(info_dict['unp_set']))
unp_len_df = pd.read_csv(unp_len_file_path, sep='\t', usecols=['Entry', 'Length'])
sifts_df_3 = sifts_demo.add_unp_len_SIFTS(sifts_df=sifts_df_2, unpLen_df=unp_len_df)

# Add Segment Info
seg_df = sifts_demo.get_seg_info_from_uniprot_segments_file(outputPath=seg_edited_file_path) # Download a new file
sifts_demo.add_seg_info_to_SIFTS(sifts_df=sifts_df_3, seg_df=seg_df, outputPath=add_unpLen_segInfo_file_path)
```

### Update Basic DataSet
```py
import sys
import pandas as pd
sys.path.append('./')
from SIFTS_unit import SIFTS_unit
from UniProt_unit import UniProt_unit

demo_pdb_list_path = '../../data/demo_pdb_list.tsv'
raw_sifts_file_path = '../../data/pdb_uniprot_SIFTS_demo.csv'
add_rangeInfo_sifts_file_path = '../../data/pdb_uniprot_SIFTS_NEW0815.tsv'
add_InDe_sifts_file_path = '../../data/Mapping_Pipeline/sifts_files/pdb_uniprot_SIFTS_delwithInDe_0815.tsv'
unp_len_file_path = '../../data/Mapping_Pipeline/unp_len_list_0817.tsv'
new_seg_edited_file_path = '../../data/Mapping_Pipeline/sifts_files/uniprot_segments_observed_edited_new.tsv' # Has Not exit
add_unpLen_segInfo_file_path = '../../data/Mapping_Pipeline/sifts_files/???'

demo_pdb_list = pd.read_csv(demo_pdb_list_path, sep='\t', usecols=['pdb_id'])['pdb_id']
old_pdb_list = pd.read_csv(add_InDe_sifts_file_path, sep='\t', usecols=['pdb_id'])['pdb_id']

new_pdb_list = list(set(demo_pdb_list) - set(old_pdb_list))

sifts_demo = SIFTS_unit()
sifts_demo.set_lists(new_pdb_list, [])

fail_list = sifts_demo.get_raw_SIFTS(raw_sifts_file_path)
sifts_df_1 = sifts_demo.handle_SIFTS()
sifts_df_2 = sifts_demo.deal_with_insertionDeletion_SIFTS(sifts_df=sifts_df_1)

sifts_df_1.to_csv(add_rangeInfo_sifts_file_path, sep='\t', header=False, index=False, mode='a+')
sifts_df_2.to_csv(add_InDe_sifts_file_path, sep='\t', header=False, index=False, mode='a+')

# Update(add) UniProt Length Info
unp_list = sifts_df_2.apply(lambda x: x['UniProt'].split('-')[0], axis=1).drop_duplicates().tolist()
uniprot_demo = UniProt_unit()
uniprot_demo.get_info_from_uniprot(['id', 'length'], unp_len_file_path, unp_list=unp_list)
unp_len_df = pd.read_csv(unp_len_file_path, sep='\t', usecols=['Entry', 'Length'])
sifts_df_3 = sifts_demo.add_unp_len_SIFTS(sifts_df=sifts_df_2, unpLen_df=unp_len_df)

# Update(add) Segment info
seg_df = sifts_demo.get_seg_info_from_uniprot_segments_file(outputPath=new_seg_edited_file_path) # Download a new file
sifts_df_4 = sifts_demo.add_seg_info_to_SIFTS(sifts_df=sifts_df_3, seg_df=seg_df)

sifts_df_4.to_csv(add_unpLen_segInfo_file_path, sep='\t', mode='a+', index=False, header=False)
```

### self.get_info_from_uniprot_pdb_file()

#### Parameters

```py
def get_info_from_uniprot_pdb_file(self, filePath=False, related_unp=False, related_pdb=False):[...]
```

#### Usage
Get pdb_list and unp_list from uniprot_pdb.csv.gz from [SIFTS](http://www.ebi.ac.uk/pdbe/docs/sifts/quick.html "Reference"). ([Download Link](ftp://ftp.ebi.ac.uk/pub/databases/msd/sifts/flatfiles/csv/uniprot_pdb.csv.gz "Click to download the file"))
* To get pdb_list and unp_list via downloaded file(uniprot_pdb.csv), specify a path: ```filePath = '../../data/uniprot_pdb.csv'``` e.g
* With no specified path, the function will download the file automatically and save the file in ```filePath = SIFTS_unit.CONFIG['DOWNLOAD_FOLDER'] + 'uniprot_pdb_%s.csv.gz' % (time.strftime("%Y_%m_%d", time.localtime()))```.
* With specified ```related_unp=cgc_unp_list```, the function will return a dict with related pdb_set and unp_set.
* With specified ```related_pdb=cgc_pdb_list```, the function will return a dict with related pdb_set and unp_set.

#### Return

```py
return {'pdb_set':set(...), 'unp_set':set(...)}
```

---

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

### self.add_mmcif_info_SIFTS()

#### Parameters

```py
def :[...]
```

#### Usage
Get


#### Return

```py
return
```














```mermaid
graph TB
  Z_0("self.set_lists(pdb_list, [])")
  A_0("get_info_from_uniprot_pdb_file()")
  A_1("get_raw_SIFTS()")
  C("")
```
