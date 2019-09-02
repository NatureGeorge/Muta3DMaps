# FILENAME: UniProt_Unit.py
> Last modified time: 2019-09-02 By ZeFeng Zhu

### self.get_info_from_uniprot()

#### Parameters

```py
def get_info_from_uniprot(self, usecols, outputPath, from_list=False, from_list_file_path=False, sep='\t', chunksize=100, header=None):[...]
```

#### Usage
Get Data Via UniProt RESTful API:
* ID Mapping $\rightarrow$ Set 'from' variable in ```params```
* Collect Info $\rightarrow$ Set 'columns' variable in ```params```
* ...

```py
unp_list_file_path = '../demo_files/unp_list.txt'
outputPath = '../demo_files/info_of_unp.tsv'
uniprot_demo = UniProt_unit()
uniprot_demo.get_info_from_uniprot(['id','length'], outputPath, unp_list_file_path=unp_list_file_path)
```

### Details about UniProt RESTful API

```URL = 'https://www.uniprot.org/uploadlists/'```

#### Reference
##### Values and Description of Available params
> https://www.uniprot.org/help/api_queries

```json
params = {
        "from": "Abbreviations of database name",
        "to": "Abbreviations of database name",
        "format": "Format in which to return results",
        "columns": "Comma-separated list of column names",
        "query": "Query String",
    }

```

##### Abbreviations of database name
> https://www.uniprot.org/help/api_idmapping

```ACC+ID, ACC, ID, EMBL_ID, P_REFSEQ_AC, PDB_ID, ENSEMBL_TRS_ID, ENSEMBL_ID, ...```

##### Columns Names
> https://www.uniprot.org/help/uniprotkb_column_names

```id,length,reviewed,comment(ALTERNATIVE%20PRODUCTS),feature(ALTERNATIVE%20SEQUENCE)'...```

### Return
This function returns a ```bool``` variable that represent whether the process is function well or not

---
