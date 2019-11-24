---
title: "Guidence"
author: ZeFeng Zhu
date: Oct 30, 2019
output:
  word_document:
    path: C:/Users/Nature/Desktop/MySoftWare/Guidence.docx
export_on_save:
pandoc: true
---

# Guidence

[toc]

## 1. 建立数据集

### 1.1 UniProt ID Mapping数据集的创建

可以通过 UniProt RESTful API(https://www.uniprot.org/uploadlists/)来获取ID映射数据

#### 1.1.1 在使用 UniProt ID Mapping API过程中的关键参数

```py
params = {
    'from': 'ACC+ID',
    'to': 'ACC',
    'format': 'tab',
    'columns': 'id,length,reviewed,comment(ALTERNATIVE%20PRODUCTS),feature(ALTERNATIVE%20SEQUENCE)'...
    'query': list_str, # "comma-separated id list"
      }
```

```rst
  .. csv-table:: Details of ID Abbreviation
      :header: "Abbreviation", "Name", "Direction"

      "ACC+ID", "UniProtKB AC/ID", "from"
      "ACC", "UniProtKB AC", "both"
      "ID", "UniProtKB ID", "both"
      "EMBL_ID", "EMBL/GenBank/DDBJ", "both"
      "REFSEQ_NT_ID", "RefSeq Nucleotide", "both"
      "P_REFSEQ_AC", "RefSeq Protein", "both"
      "PDB_ID", "PDB", "both"
      "ENSEMBL_TRS_ID", "Ensembl Transcript", "both"
      "ENSEMBL_ID", "Ensembl", "both"
      "...", "...", "..."
```

详情可参考 https://www.uniprot.org/help/api_idmapping


#### 1.1.2 输入数据细节

```rst
  .. csv-table:: Details of other parameters
      :header: "param", "Explanation", "Reference"

      "columns", "comma-separated list of column names", "https://www.uniprot.org/help/api_queries"
      "columns", "Lists the column names for programmatic (RESTful) access to tab-separated or Excel downloads of UniProtKB search results", "https://www.uniprot.org/help/uniprotkb_column_names"
      "format", "html | tab | xls | fasta | gff | txt | xml | rdf | list | rss", "https://www.uniprot.org/help/api_queries"
      "query", "query text/id(s)", "https://www.uniprot.org/help/text-search"
```

#### 1.1.3 输入数据样例

* ```group_df```

| mutation_unp  |   GENE    |RefSeq_protein |RefSeq_nucleotide|
|:-------------:|:-------------:|:-------------:|:-------------:|
|     K45E      |      SAMD11      |   NP_689699   |   NM_152486   |
|     P293A     |      SAMD11      |   NP_689699   |   NM_152486   |
|     G76S      |      AGRN      |   NP_940978   |   NM_198576   |
|     N105I     |      AGRN      |   NP_940978   |   NM_198576   |
|     A375S     |      AGRN      |   NP_940978   |   NM_198576   |
|     ...     |      ...      |   ...   |   ...   |

#### 1.1.4 目标列及对应信息

* id_col: ```RefSeq_protein```
* muta_col: ```mutation_unp```
* gene_col: ```GENE```
* id_type: Refseq Protein -> P_REFSEQ_AC
* muta_type: mutation in UniProt Site

#### 1.1.5 样例脚本

```py
from UniProt_unit import UniProt_unit
# 设置初始参数值
id_col = 'RefSeq_protein'
id_type = 'P_REFSEQ_AC'
muta_col = 'mutation_unp'
gene_col = 'GENE'
usecols = ['id', 'genes', 'reviewed', 'comment(ALTERNATIVE%20PRODUCTS)', 'organism', 'protein%20names']  # Necessary Columns
reportPath = 'RefSeq_protein_mapping_Report.txt'  # Report of Data Processing
rawOutputPath = 'RefSeq_protein_mapping.tsv' # OutPut File of ID Mapping (RAW)
handledOutputPath = 'RefSeq_protein_mapping_modified.tsv'  # OutPut File of ID Mapping (Final Result)
# Add constraint/filter to data
constraint_dict = {
    "GENE_status": (False, "ne"),  # 'ne' for !=
    "Status": ("reviewed", "eq"),  # 'eq' for ==
    "unp_map_tage": ("Untrusted & No Isoform", "ne")
}

# 初始化
unp_demo = UniProt_unit(group_df, id_col, id_type, usecols, reportPath, muta_col=muta_col, gene_col=gene_col)
# Return True if get RAW Result Successfully
unp_demo.get_raw_ID_Mapping(rawOutputPath)
# Deal with different situations
handled_ID_Mapping = unp_demo.handle_ID_Mapping()
# Add Gene Status
unp_demo.getGeneStatus(handled_ID_Mapping)
# Label Mapping Status
unp_demo.label_mapping_status(handled_ID_Mapping, constraint_dict)
# close the file-handle of report
unp_demo.report.close()
# Output the final result
handled_ID_Mapping.to_csv(handledOutputPath, sep='\t', index=False)
```

结果存储于 ```pandas.DataFrame```类型中。

#### 1.1.6 标签以及报告的解释
##### 1.1.6.1 关于 ```unp_map_tage```

* ```Untrusted & No Isoform```

是指UniProt存在Isoform但是Mapping结果没有明确给出是Map上哪条Isoform,转录本序列与蛋白序列不一致

* ```Trusted & No Isoform```

是指UniProt不存在Isoform,Mapping结果没问题

* ```Trusted & Isoform```

是指UniProt存在Isoform,Mapping结果没问题

##### 1.1.6.2 关于 ```Mapping_status```

* ```Yes```: 可信的结果，进行后续的PDB Mapping; 通过```constraint_dict```的限制
* ```Error```: 一个id对应多个UniProt; 通过```constraint_dict```的限制
* ```No```: 不可信的结果; 未通过```constraint_dict```的限制

##### 1.1.6.3 关于 ```GENE_status```

* ```False```: ```Gene names```的第一个元素 与 refSeq的 ```GENE```不匹配
* 其他值即为对应 ```GENE```

##### 1.1.6.4 输出范例

|     Entry     |  Gene names   |    Status     |Alternative products (isoforms)|   Organism    | Protein names |canonical_isoform| unp_map_tage  |   yourlist    |    UniProt    |     GENE      |  GENE_status  |Mapping_status |
|:-------------:|:-------------:|:-------------:|:-------------:|:-------------:|:-------------:|:-------------:|:-------------:|:-------------:|:-------------:|:-------------:|:-------------:|:-------------:|
|    Q4VBY6     |    CDC2L2     |  unreviewed   |      NaN      |Homo sapiens (Human)|CDC2L2 protein (Fragment)|      NaN      |Trusted & No Isoform|   NP_076916   |    Q4VBY6     |    CDK11A     |     False     |      No       |
|    Q9C0B2     |CFAP74 C1orf222 KIAA1751|   reviewed    |ALTERNATIVE PRODUCTS:  Event=Alternative splic...|Homo sapiens (Human)|Cilia- and flagella-associated protein 74|   Q9C0B2-1    |Untrusted & No Isoform| NP_001291289  |    Q9C0B2     |   KIAA1751    |     False     |      No       |
|    Q96NU1     |    SAMD11     |   reviewed    |ALTERNATIVE PRODUCTS:  Event=Alternative promo...|Homo sapiens (Human)|Sterile alpha motif domain-containing protein ...|   Q96NU1-3    |Untrusted & No Isoform|   NP_689699   |    Q96NU1     |    SAMD11     |     True      |      No       |
|    P43489     |TNFRSF4 TXGP1L |   reviewed    |      NaN      |Homo sapiens (Human)|Tumor necrosis factor receptor superfamily mem...|      NaN      |Trusted & No Isoform|   NP_003318   |    P43489     |    TNFRSF4    |     True      |      Yes      |


##### 1.1.6.5 Report Statistis

下列统计信息将输出至样例脚本中定义的 ```reportPath``` 对应文件中

* all RefSeq Count: 8450
* Unmapped Results: 34
* Trusted Results: 8066
  * ```handled_df[handled_df['Mapping_status'] == 'Yes']```
* Untrusted Resultes: 350
  * Error Results: 1 ```{NP_002115: [P01911, Q29974]}, handled_df[handled_df['Mapping_status'] == 'Error']```
  * Others: ```handled_df[handled_df['Mapping_status'] == 'No']```

### 1.2 获取整合好的突变信息
```py
print(unp_demo.muta_li)
```

```txt
RefSeq_protein
NP_000005                              [R1297C, V1000I, C972Y]
NP_000006    [L24I, R64W, R64Q, I114T, D122N, L137F, Q145P,...
NP_000007    [R31H, E43K, Q45R, Q49E, A52V, R53C, R53H, A56...
NP_000008    [R46W, P55L, G90S, G92C, L93I, D94H, I105N, R1...
NP_000009    [S72F, G76E, P89S, P91Q, P91L, S110Y, F113L, N...
NP_000010    [I52T, Q73P, N93S, G100E, Q101K, N123K, K124E,...
NP_000011    [K8T, K8N, V18M, P30S, V32G, C34Y, C36Y, S38I,...
NP_000012    [Q15H, R35Q, N39Y, D40N, R42L, V63G, A79T, A79...
NP_000013    [A329V, M310T, L304R, P297Q, S291L, R282Q, P27...
NP_000014    [P30L, L31P, V32A, R34C, R34H, R34L, V37L, A53...
NP_000015                           [R16G, Q27E, A119D, T164I]
...
```

### 1.3 SIFTS数据集的创建
将Uniprot Sequence(包括isoform)上的突变位点信息映射到PDB结构上，就必须进行序列比对以获得残基水平的一一对应关系。而这个工作已由[SIFTS](https://www.ebi.ac.uk/pdbe/docs/sifts/index.html "Link")完成。

SIFTS已将PDB的SEQRES序列与各UniProt Isoform进行序列比对，并且经常更新。通过调用SIFTS的API接口(http://www.ebi.ac.uk/pdbe/api/mappings/all_isoforms/:accession)，解析JSON格式的返回数据，即可获取到PDB相关链与Uniprot Isoform的残基位置对应关系与范围以及序列相似性等信息。残基对应关系中的PDB残基号码是从1开始计数的SEQRES残基索引。

pdb_id	|	chain_id	|	UniProt	|	identity	|	identifier | pdb_start	|	pdb_end	|	unp_start	|	unp_end | is_canonical	|	start	|	end	|	entity_id	|	struct_asym_id
-|-|-|-|-|-|-|-|-|-|-|-|-|-
2znj	|	A	|	B0S4P3	|	1	|	B0S4P3_DESHA	|	21	|	308	|	1	|	288	|	True	|	{author_residue_number: null, author_insertion_code: "", residue_number: 21}	|	{author_residue_number: 288, author_insertion_code: "", residue_number: 308}	|	1	|	A

解析SIFTS提供的残基对应关系，可以间接判断出SEQRES序列相对于完整蛋白序列(Uniprot Sequence)的差异，包括Insertion，Deletion以及SEQRES序列头尾的差异部分 (SEQRES序列在头或尾具有而Uniprot序列不具有的部分)。我们将原本SIFTS提供的```pdb_start, pdb_end, unp_start, unp_end``` 残基水平对应关系转换为区间格式，如下表```sifts_pdb_range, sifts_unp_range```所示。

pdb_id	|	chain_id	|	UniProt	| identity |... |is_canonical | ... | sifts_pdb_range | sifts_unp_range | delete | sifts_range_tage
-|-|-|-|-|-|-|-|-|-|-
3WZU	|A	|O14733-4	|0.978	|...|False|...	|\[\[2, 318]]	|\[\[103, 426]]	|False |Deletion
3WZU	|A	|O14733	|1.000	|...|True|...	|\[\[2, 318]]	|\[\[103, 419]]	| False |Safe
1XZ9	|A	|P55196-1	|0.990 |	...|False|...	|\[\[6, 52], \[54, 101]]	|\[\[985, 1031], \[1032, 1079]]	| False |Insertion

大多数情况下，SIFTS提供的残基水平对应关系转换为区间后只会有一个子区间，pdb与unp子区间对应的残基序列长度是相等的，我们记为```Safe```。但是，存在特殊情况：
* 对应子区间的长度并不总是相等 $\rightarrow$ pdb链相对于uniprot Isoform Sequence发生了Deletion
* 有时pdb与unp的对应关系会有多个子区间 $\rightarrow$ pdb链相对于uniprot Isoform Sequence发生了Insertion
* unp子区间之间并不连续 $\rightarrow$ pdb链相对于uniprot Isoform Sequence发生了Insertion 以及 Deletion，存在较长差异序列

发现上述特殊情况中的```Insertion```以及```Insertion & Deletion```时后续可以调用以及封装好的函数进行pair-wise alignment(采用Biopyton模块的pairwise2)以确定具体的差异内容, 更正```sifts_pdb_range, sifts_unp_range```使得残基对应关系更为精确。

下面通过SIFTS API 建立UniProt与PDB残基水平对应关系的原始数据集。

```py
from SIFTS_unit import SIFTS_unit

# SIFTS 原始数据文件路径
raw_sifts_file_path = 'pdb_uniprot_SIFTS_raw_demo.csv'
# SIFTS 加上经过整合的覆盖范围区间信息后的文件路径
add_rangeInfo_sifts_file_path = 'pdb_uniprot_SIFTS_addRangeInfo_demo.tsv'
# SIFTS 判断PDB链相对UniProt Isoform的序列是否有Intertion,Deletion,并加上标签
add_InDe_sifts_file_path = 'pdb_uniprot_SIFTS_delwithInDe_demo.tsv'
# 文件保存文件夹
SIFTS_unit.CONFIG['DOWNLOAD_FOLDER'] = '/data/zzf/SIFTS_files/'
# 初始化
sifts_demo = SIFTS_unit()
# 获取SIFTS当前可提供UniProt-PDB映射信息的UniProt与PDB ID
info_dict = sifts_demo.get_info_from_uniprot_pdb_file()
# 设置预备获取信息的PDB ID集
sifts_demo.pdb_list = sorted(info_dict['pdb_set'])
# 获取SIFTS 原始数据
sifts_demo.get_raw_SIFTS(outputPath='%s%s' % (SIFTS_unit.CONFIG['DOWNLOAD_FOLDER'], raw_sifts_file_path))
# 处理SIFTS原始数据: 加上经过整合的覆盖范围区间等信息
handle_sifts_df = sifts_demo.handle_SIFTS(outputPath='%s%s' % (SIFTS_unit.CONFIG['DOWNLOAD_FOLDER'], add_rangeInfo_sifts_file_path))
# 处理SIFTS原始数据: 判断PDB链相对UniProt Isoform的序列是否有Intertion,Deletion,并加上标签
sifts_demo.deal_with_insertionDeletion_SIFTS(sifts_df=handle_sifts_df, outputPath='%s%s' % (SIFTS_unit.CONFIG['DOWNLOAD_FOLDER'], add_InDe_sifts_file_path))
```

结果存储于 ```pandas.DataFrame```类型中。

### 1.4 PDB文件的下载

#### 1.4.1 样例脚本
```py
# 预下载的PDB
pdbs = ['1A02', '3KBZ', '3KC0', '3KC1', '3KMU', '3KMW', '3KYC', '3KYD', ...]
# 下载路径
path = "./"
# 默认下载MMCIF格式
mpw = MPWrapper(path)
fail = mpw.ftp_retrieve_batch(pdbs)
# 查看下载失败的PDB ID集
print(fail)
```

### 1.5 MMCIF文件信息的获取与整合

#### 1.5.1 获取相关的PDB ID集

获取当前研究涉及到的PDB ID:

```py
dfA = handle_sifts_df
dfB = handled_ID_Mapping
pdbs = dfA[dfA['UniProt'].isin(dfB['UniProt'])]['pdb_id'].drop_duplicates()
```

#### 1.5.2 建立MMCIF信息数据集

```py
import pandas as pd
import os
from MMCIFplus import *

# 设置文件路径
raw_mmcif_path = 'rawMMCIF2Dfrm.tsv'
handled_mmcif_path = 'handledMMCIF2Dfrm.tsv'
# 初始化
mmcif_demo = MMCIF2Dfrm()
# 下载未下载过的MMCIF文件
mmcif_demo.check_mmcif_file(pdbs)
# 接续
if os.path.exists(handled_mmcif_path):
    finished_li = set(pd.read_csv(handled_mmcif_path, sep='\t', usecols=['pdb_id'])['pdb_id'])
else:
    finished_li = []

mmcif_demo.update_mmcif_result(raw_mmcif_path, handled_mmcif_path, finished=finished_li)
```


### 1.6 Interactome3D meta-data的获取

```py
from Interactome3D_unit import Interactome3D

interactDemo = Interactome3D()
# 原文件下载路径
interactDemo.CONFIG['DOWNLOAD_FOLDER'] = './'
# 获取经过处理的整合文件
interact_df = interactDemo.get_interactions_meta(outputPath=interactDemo.CONFIG['DOWNLOAD_FOLDER']+'interactions_modified.tsv')
```

## 2 整合SIFTS数据与PDB(MMCIF格式)数据

```py
mmcif_df = pd.read_csv(handled_mmcif_path, sep='\t', converters={'chain_id':str, 'asym_id':str, 'entity_id':str})
sifts_mmcif_df = sifts_demo.add_mmcif_info_SIFTS(sifts_df=handle_sifts_df, mmcif_df=mmcif_df)
```

## 3 将突变从UniProt映射至PDB链

### 3.1 对Deletion相关数据进行修正

```py
# UniProt Fasta序列文件夹，可至UniProt官网下载
unp_fasta_files_path = './fasta_files/%s.fasta'
update_sifts_mmcif_df = sifts_demo.update_range_info_SIFTS(unp_fasta_files_path, sifts_df=sifts_mmcif_df)
```

### 3.2 映射

```py
# 将突变数据整合
update_sifts_mmcif_df['mutation_unp'] = update_sifts_mmcif_df.apply(lambda x: unp_demo.muta_li[x['UniProt']], axis=1)
# PDB突变位点信息
muta_info_li = []
update_sifts_mmcif_df['mutation_pdb'] = update_sifts_mmcif_df.apply(lambda x: SIFTS_unit.map_muta_from_unp_to_pdb(x, 'mutation_unp', 'new_sifts_unp_range', 'new_sifts_pdb_range', muta_info_li, unp_fasta_files_path) if not isinstance(x['new_sifts_pdb_range'], float) and not isinstance(x['mutation_unp'], float) else np.nan, axis=1)
# 映射情况信息
update_sifts_mmcif_df['muta_map_info'] = pd.Series(muta_info_li, index=update_sifts_mmcif_df.dropna(subset=['mutation_unp']).index)
```

## 4 参考数据库与工具

[1] H.M. Berman, K. Henrick, H. Nakamura (2003) Announcing the worldwide Protein Data Bank Nature Structural Biology 10 (12): 980.

[2] H.M. Berman, K. Henrick, H.Nakamura, J.L. Markley (2007) The Worldwide Protein Data Bank (wwPDB): Ensuring a single, uniform archive of PDB data Nucleic Acids Res. 35 (Database issue): D301-3.

[3] wwPDB consortium. (2019) Protein Data Bank: the single global archive for 3D macromolecular structure data. Nucleic Acids Res 47: D520-D528 doi: 10.1093/nar/gky949.

[4] H.M. Berman, J. Westbrook, Z. Feng, G. Gilliland, T.N. Bhat, H. Weissig, I.N. Shindyalov, P.E. Bourne. (2000) The Protein Data Bank Nucleic Acids Research, 28: 235-242.

[5] http://www.wwpdb.org/

[6] http://www.rcsb.org/

[7] https://www.ebi.ac.uk/pdbe/

[8] https://pdbj.org/

[9] UniProt: a worldwide hub of protein knowledge Nucleic Acids Res. 47: D506-515 (2019)

[10] Dana, Gutmanas et al., Nucleic Acids Research 47, D482 (2019), https://www.ebi.ac.uk/pdbe/docs/sifts/index.html

[11] Mosca R, Céol A, Aloy P, Interactome3D: adding structural details to protein networks, Nature Methods (2013) 10(1):47-53, doi:10.1038/nmeth.2289, https://interactome3d.irbbarcelona.org
