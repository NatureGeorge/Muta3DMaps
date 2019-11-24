# Mapping Pipeline
> 20190830

### Qusetion
* 解析SIFTS提供的残基对应关系，可以间接判断出SEQRES序列相对于完整蛋白序列(Uniprot Sequence)的差异，包括Insertion，Deletion以及SEQRES序列头尾的差异部分 (SEQRES序列在头或尾具有而Uniprot序列不具有的部分)。发现上述差异时我们会进行pair-wise alignment以确定具体的差异内容。```这块就开始alignment了?```
* _pdbx_poly_seq_scheme下的ndb_seq_num、pdb_seq_num以及auth_seq_num ```意义？要写意思 而不是标识符```
* sifts利用的信息 ```有哪些```
* 所有因素都会除以Uniprot Sequence的长度。 权重确定方法采用基于经验的主观赋权——层次分析法(AHP)。相对重要性矩阵如下。一致性检验通过，最终得到权重如下。 多些细节，图要有解释说明
* 异聚体按蛋白组成来分组，同组内选择含有按上述原则选择的链的。```什么意思```
* Loop-content大于30%(?)的。```有没有疑问```
* 若overlap/该链长度 <= 0.3 且 (该链长度-overlap)/已加入代表集的链的非冗余长度 >= 0.2，则将该链加入代表集，循环此操作。```这句话表示什么？```

### Part1
蛋白质三维结构数据是对蛋白质空间构象的描述，其中包含一维序列信息。
通常情况，由于技术限制，单一三维结构数据往往只覆盖蛋白质的部分序列，不能表示完整的蛋白结构。在这种情况下，多个结构数据可以表示单一蛋白质，因此我们采用UniProtKB的标识符来辨别蛋白质。
#### Protein Structural Data Source
##### Crystal Structure
将蛋白序列层面的氨基酸突变位点映射到三维结构数据中相应的位点是一个必要的任务，反之亦然。这个任务关键在于知悉三维结构数据与蛋白完整序列的残基层面的一一对应关系统一
蛋白实验晶体结构数据来自Protein Data Bank。我们统一采用mmcif格式来解析，能够更高效准确地获取结构覆盖范围、缺失范围、配体结合位点、修饰位点、突变位点以及是否包含核酸链等信息。
##### Model Structure
实验晶体未能覆盖的部分由蛋白质结构模型来补充。蛋白质单体结构模型有两个来源：Swiss-Model Repository以及ModBase，相互作用模型来自于Interactome3D。
> SWISS-MODEL Repository（SMR）是由自动SWISS-MODEL同源建模管道生成的带注释的3D蛋白质结构模型的数据库。其定期更新确保目标覆盖范围完整，使用最新的序列和模板结构数据库构建模型，并充分利用底层建模管道中的改进。其在建模管道中即采用了QMEAN方法来对模型进行质量评估。只有高质量的模型才会被导入SMR：对于少于100个残基的模型，序列同一性必须超过30％；对于具有大于100个残基的模型，QMEAN4得分必须大于-5。它目前拥有> 40万个高质量模型，覆盖了近20％的Swiss-Prot / UniProtKB条目。对于智人（Homo sapiens），该数据库提供可供下载的模型结构数据以及坐标文件，采用的UniProtKB映射关系可随时更新。(https://www.swissmodel.expasy.org/repository)

> ModBase是由一个依赖于PSI-BLAST、HHSuite和MODELLER程序的自动建模管道（ModPipe）生成的蛋白质结构模型的可查询数据库。其会通过若干附加的质量评估标准来过滤模型，最终只有质量标准值高于指定阈值或E值<10 - 4的模型才包含在最终模型集中。目前包含近3000万个可靠的模型，每周更新模型的元数据。可经由UniprotKB、GI标识符、PDB的id等来搜索模型。对于智人同样也提供可供下载的模型结构数据与座标文件。(https://modbase.compbio.ucsf.edu/modbase-cgi/index.cgi)

> Interactome3D是一个用于蛋白质-蛋白质相互作用的结构注释和建模的资源平台。其提供了八种生物体模型中超过12,000个蛋白质相互作用的结构细节。其从蛋白质数据库（PDB）收集蛋白质和相互作用的结构。而对于模型，单个蛋白质的同源模型从Modbase收集；蛋白相互作用模型使用Modeller建模而成，模板来自于PDB的生物单元文件（biological unit files）或者3did的域-域结构模板（domain-domain structural templates）。对于从ModBase获取的同源模型，其附加若干质量评估标准来过滤模型。对于采用MODELLER建模的来的模型，其选择DOPE分数低于60、主链不出现结（knots）的模型，并且丢弃具有少于5对残基间相互作用的模型。(https://interactome3d.irbbarcelona.org)

~~我们分别从SMR以及ModBase抽取了部分的具有相同uniprot覆盖范围的模型，通过SAVES(http://servicesn.mbi.ucla.edu/SAVES/)进行模型质量评估，发现SMR的模型质量普遍优于Modbase的模型。因此，蛋白质单体结构模型优先选取来自于Swiss-Model Repository的而后选取来自ModBase的~~。必会获取的模型结构数据元信息包括序列相似性，模板pdb链以及相关模型质量评分。

#### ID Mapping
处理的突变数据集可以来自于COSMIC等，突变以转录本(ENST，Ref等)为载体。转录本或GeneName等与蛋白质的对应关系转换由UniProt的ID Mapping Tool (https://www.uniprot.org/uploadlists/)完成。

在利用UniProt的ID Mapping Tool后，转录本上的突变便成了对应isoform蛋白的突变，位点信息不变。Uniprot Isoform与转录本呈一对多关系(_部分存在多对多关系，需要仔细检查_)，且对于一个Uniprot Entry，可能存在多个isoform， 但是只有一个isoform为canonical。

~~我们将同一Uniprot Entry下的isoform上的突变位点，经由pair-wise alignment(采用Biopyton模块的pairwise2)，映射到Canonical Sequence上，进行后续操作，舍弃部分未成功映射的突变~~。

#### Residue Level Mapping
##### SIFTS
将Uniprot Sequence(包括isoform)上的突变位点信息映射到PDB结构上，就必须进行序列比对以获得残基水平的一一对应关系。而这个工作已由SIFTS完成。
> Structure Integration with Function, Taxonomy and Sequence (SIFTS) is a project in the PDBe-KB resource for residue-level mapping between UniProt and PDB entries. SIFTS also provides annotation from the IntEnz, GO, InterPro, Pfam, CATH, SCOP, PubMed, Ensembl and Homologene resources. The information is updated and released every week concurrently with the release of new PDB entries and is widely used by resources such as RCSB PDB, PDBj, PDBsum, Pfam, SCOP and InterPro. (https://www.ebi.ac.uk/pdbe/docs/sifts/index.html)

SIFTS已将PDB的SEQRES序列与各UniProt Isoform进行序列比对，并且经常更新。通过调用SIFTS的API接口(http://www.ebi.ac.uk/pdbe/api/mappings/all_isoforms/:accession)，解析JSON格式的返回数据，即可获取到PDB相关链与Uniprot Isoform的残基位置对应关系与范围以及序列相似性等信息。残基对应关系中的PDB残基号码是从1开始计数的SEQRES残基索引。
###### Raw Data That Collected From SIFTS (Example of 2znj_A_B0S4P3)
pdb_id	|	chain_id	|	UniProt	|	identity	|	identifier | pdb_start	|	pdb_end	|	unp_start	|	unp_end | is_canonical	|	start	|	end	|	entity_id	|	struct_asym_id
-|-|-|-|-|-|-|-|-|-|-|-|-|-
2znj	|	A	|	B0S4P3	|	1	|	B0S4P3_DESHA	|	21	|	308	|	1	|	288	|	True	|	{author_residue_number: null, author_insertion_code: "", residue_number: 21}	|	{author_residue_number: 288, author_insertion_code: "", residue_number: 308}	|	1	|	A

解析SIFTS提供的残基对应关系，可以间接判断出SEQRES序列相对于完整蛋白序列(Uniprot Sequence)的差异，包括Insertion，Deletion以及SEQRES序列头尾的差异部分 (SEQRES序列在头或尾具有而Uniprot序列不具有的部分)。我们将原本SIFTS提供的```pdb_start, pdb_end, unp_start, unp_end``` 残基水平对应关系转换为区间格式，如下表```sifts_pdb_range, sifts_unp_range```所示。
大多数情况下，SIFTS提供的残基水平对应关系转换为区间后只会有一个子区间，pdb与unp子区间对应的残基序列长度是相等的，我们记为```Safe```。但是，存在特殊情况：
* 对应子区间的长度并不总是相等 $\rightarrow$ pdb链相对于uniprot Isoform Sequence发生了Deletion
* 有时pdb与unp的对应关系会有多个子区间 $\rightarrow$ pdb链相对于uniprot Isoform Sequence发生了Insertion
* unp子区间之间并不连续 $\rightarrow$ pdb链相对于uniprot Isoform Sequence发生了Insertion 以及 Deletion，存在较长差异序列

###### Processed SIFTS Data (part)
pdb_id	|	chain_id	|	UniProt	| identity |... |is_canonical | ... | sifts_pdb_range | sifts_unp_range | delete | sifts_range_tage
-|-|-|-|-|-|-|-|-|-|-
3WZU	|A	|O14733-4	|0.978	|...|False|...	|\[\[2, 318]]	|\[\[103, 426]]	|False |Deletion
3WZU	|A	|O14733	|1.000	|...|True|...	|\[\[2, 318]]	|\[\[103, 419]]	| False |Safe
1XZ9	|A	|P55196-1	|0.990 |	...|False|...	|\[\[6, 52], \[54, 101]]	|\[\[985, 1031], \[1032, 1079]]	| False |Insertion
1XZ9	|A	|P55196	|1.000	|...|True|...	|\[\[6, 101]]	|\[\[1001, 1096]]	|False |Safe
2WTD	|A	|O96017	|0.997	|...|True|...	|\[\[7, 329]]	|\[\[209, 531]]	|False |Safe
2WTD	|A	|O96017-2	|0.210 |	...|False|...	|\[\[1, 39], \[42, 55]]	|\[\[151, 196], \[197, 210]] |False	|Insertion & Deletion
3CMU	|A	|P0A7G6	|0.987	|...|True|...	|\[\[5, 309], \[325, 658], \[674, 1007], \[1023, 135...	|\[\[31, 335], \[2, 335], \[2, 335], \[2, 335], \[2, ...|True	|Insertion & Deletion
5LOH	|A	|Q96GX5-2	|0.718 |...|False|...	|\[\[4, 197], \[202, 218], \[256, 341]]	|\[\[1, 194], \[739, 209], \[755, 840]]	|True |Insertion & Deletion
5LOH	|A	|Q96GX5	|0.986	|...|True|...	|\[\[4, 197], \[202, 341]]	|\[\[1, 194], \[740, 879]]	|False|Insertion & Deletion

发现上述特殊情况中的```Insertion```以及```Insertion & Deletion```时我们会进行pair-wise alignment(采用Biopyton模块的pairwise2)以确定具体的差异内容, 更正```sifts_pdb_range, sifts_unp_range```使得残基对应关系更为精确。

##### MMCIF
PDB的mmcif格式中，我们获取但不限于下列信息：
* ```_pdbx_poly_seq_scheme.ndb_seq_num``` $\rightarrow$ 从1开始计数的SEQRES残基索引
* ```_pdbx_poly_seq_scheme.pdb_seq_num``` $\rightarrow$ pdb文件中作者定义的残基号码索引
* ```_pdbx_poly_seq_scheme.auth_seq_num``` $\rightarrow$ 以及pdb文件中作者定义的残基号码索引但是用'?'符号来标识Missing Residue的残基号码索引
* ```_pdbx_poly_seq_scheme.pdb_ins_code``` $\rightarrow$ pdb文件中作者定义的残基号码索引的inside code
> Reference: http://mmcif.rcsb.org/dictionaries/mmcif_pdbx_v50.dic/Items/_pdbx_poly_seq_scheme.auth_seq_num.html
###### Data That Collected From MMCIF-Format PDB Files (Example of 5PAW_B)
> ```...``` stand for Ellipsis
```DataFrame
entity_id                                                                             2
protein_type                                                             polypeptide(L)
pdb_id                                                                             5PAW
chain_id                                                                              B
_pdbx_poly_seq_scheme.mon_id                                          IVG...SRKVGDSP...
_pdbx_poly_seq_scheme.pdb_mon_id                                      IVG...SR????SP...
_pdbx_poly_seq_scheme.auth_mon_id                                     IVG...SR????SP...
_pdbx_poly_seq_scheme.ndb_seq_num          1;2;3;...162;163;164;165;166;167;168;169;...
_pdbx_poly_seq_scheme.pdb_seq_num    213;214;215;...374;375;376;377;378;379;380;381;...
_pdbx_poly_seq_scheme.auth_seq_num   213;214;215;...374;375;  ?;  ?;  ?;  ?;380;381;...
_pdbx_poly_seq_scheme.pdb_ins_code    .;.;.;.;.;.;.;.;.;.;.;.;.;.;.;.;.;.;.;.;.;.;.;...
UNK_ALL_IN_CHAIN                                                                  False
method                                                                X-RAY DIFFRACTION
initial_version_time                                                         2017-06-21
newest_version_time                                                          2018-02-21
resolution                                                                          2.2
pdb_contain_chain_type                                                          protein
contains_unk_in_chain_pdb                                                         False
metal_ligand_content                  [('CA', '501', 'GLU', '270'), ('CA', '501', 'A...
mutation_content                                                                      ?
mutation_num                                                                          0
metal_ligand_num                                                                      4
Modification_num                                                                      0
seqres_len                                                                          254
coordinates_len                                                                     250
Modification_index                                                                   []
mis_index                                                          [163, 164, 165, 166]
mis_range                                                                  [[163, 166]]
resolution_score                                                                    2.2
protein_chain_and_length                                        [(56, "A"), (250, "B")]
```

##### Mutation Mapping
利用MMCIF提取的```_pdbx_poly_seq_scheme.auth_seq_num, _pdbx_poly_seq_scheme.pdb_ins_code```信息以及上述的SEQRES与Uniprot Sequence的对应关系(```sifts_pdb_range, sifts_unp_range```)即可得到特定突变位点在PDB的SEQRES中的相对位置进而得知其在pdb文件中的由作者定义的残基号码索引为何。反向操作亦然。

###### Example of 1XZ9_A_P55196-1
__Alignment__

```alignment
Query = sp|P55196-1|AFAD_HUMAN
Sbjct = 1XZ9:A|PDBID|CHAIN|SEQUENC

Query: 985  LRKEPEIITVTLKKQNGMGLSIVAAKGAGQDKLGIYVKSVVKGGAAD-DGRLAAGDQLLS 1043
            LRKEPEIITVTLKKQNGMGLSIVAAKGAGQDKLGIYVKSVVKGGAAD DGRLAAGDQLLS
Sbjct: 6    LRKEPEIITVTLKKQNGMGLSIVAAKGAGQDKLGIYVKSVVKGGAADVDGRLAAGDQLLS 65


Query: 1044 VDGRSLVGLSQERAAELMTRTSSVVTLEVAKQGAIY 1079
            VDGRSLVGLSQERAAELMTRTSSVVTLEVAKQGAIY
Sbjct: 66   VDGRSLVGLSQERAAELMTRTSSVVTLEVAKQGAIY 101
```

__SIFTS Residue level Mapping__
pdb_id	|	chain_id	|	UniProt	| identity |... |is_canonical | ... | sifts_pdb_range | sifts_unp_range | delete | sifts_range_tage
-|-|-|-|-|-|-|-|-|-|-
1XZ9	|A	|P55196-1	|0.990 |	...|False|...	|\[\[6, 52], \[54, 101]]	|\[\[985, 1031], \[1032, 1079]]	| False |Insertion

__MMCIF Data__
```DataFrame
entity_id                                                                             1
protein_type                                                             polypeptide(L)
pdb_id                                                                             1XZ9
chain_id                                                                              A
_pdbx_poly_seq_scheme.mon_id          GPLGSLRKEPEIITVTLKKQNGMGLSIVAAKGAGQDKLGIYVKSVV...
_pdbx_poly_seq_scheme.pdb_mon_id      GPLGSLRKEPEIITVTLKKQNGMGLSIVAAKGAGQDKLGIYVKSVV...
_pdbx_poly_seq_scheme.auth_mon_id     GPLGSLRKEPEIITVTLKKQNGMGLSIVAAKGAGQDKLGIYVKSVV...
_pdbx_poly_seq_scheme.ndb_seq_num     1;2;3;4;5;6;7;8;9;10;11;12;13;14;15;16;17;18;1...
_pdbx_poly_seq_scheme.pdb_seq_num     1;2;3;4;5;6;7;8;9;10;11;12;13;14;15;16;17;18;1...
_pdbx_poly_seq_scheme.auth_seq_num    1;2;3;4;5;6;7;8;9;10;11;12;13;14;15;16;17;18;1...
_pdbx_poly_seq_scheme.pdb_ins_code    .;.;.;.;.;.;.;.;.;.;.;.;.;.;.;.;.;.;.;.;.;.;.;...
```

将Uniprot Sequence上的突变位点信息映射到蛋白质结构模型上基本无太大阻力，SMR、ModBase以及Interactome3D都支持UniProt对蛋白的唯一标识。除极少数特例，结构模型内的残基编号都与Uniprot Sequence一致。Uniprot Sequence上的突变位点即为结构模型对应号码残基。

### Part2

对于一个蛋白(Uniprot Entry)，在有多个蛋白质结构数据的情况下，为减少数据冗余性，确定蛋白结构数据代表集是必要的。代表集分为实验晶体结构代表集与结构模型代表集，其中实验晶体结构代表集优先。

在确定实验晶体结构代表集之前，会过滤掉x-ray中高于2.5的结构。确定实验晶体结构代表集的原则是同时考虑覆盖范围与质量，选择最佳的。采用加权打分制，考虑的因素有(信息从mmcif格式文件中以及SIFTS提供的信息(PDB identity>90%)中获取)：
* 该PDB链具有三维坐标的残基的数目(coordinates_len)
* 该PDB链残基丢失数目(mis_num)
* 该PDB链头尾部分Uniprot Sequence不具有的特异序列长度，超过5开始计数(mappedOut_num)
* 该PDB链中金属配体结合位点数目(metalLigand_num)
* 该PDB链中Insertion与Deletion的长度(Insertion_length, Deletion_length)
* 该PDB链中修饰残基数目与突变残基数目(Modification_num, mutation_num)

所有因素都会除以Uniprot Sequence的长度。 权重确定方法采用基于经验的主观赋权——层次分析法(AHP)。相对重要性矩阵如下。一致性检验通过，最终得到权重如下。

在同一Uniprot Entry下，对于PDB链，首先选择拥有最高打分的，加入代表集；如若打分相同，选择具有更好分辨率的(x-Ray与NMR都有的时优先考虑x-ray)；仍无法选出最优，便选择最新的结构。接着看排名次位的链，若overlap/该链长度 <= 0.3 且 (该链长度-overlap)/已加入代表集的链的非冗余长度 >= 0.2，则将该链加入代表集，循环此操作。

单体、同聚体、异聚体的代表集确定工作是分开处理。单体指代PDB结构数据文件中仅有1条长度大于20氨基酸长度的链，不含核酸。同聚体指代PDB结构中所有大于20氨基酸长度的链都属于同一Uniprot Entry，同时覆盖范围基本一致（剔除极少数覆盖范围不一致的），不含核酸链。异聚体指代PDB结构中大于20氨基酸长度的链各有各的对应Uniprot Entry，不含核酸链。对于PDB中氨基酸长度小于20的链，我们将其从PDB文件中剔除。

单体直接按上述原则选取链，直接确定出选取的PDB id；同聚体选择含有按上述原则选择的链的；异聚体按蛋白组成来分组，同组内选择含有按上述原则选择的链的。


在确认模型结构代表集之前，会过滤掉模型质量评分不合格的(来自ModBase的模型MPQS需大于1.1)、序列相似性低于30%的、长度小于20的以及Loop-content大于30%(?)的。确定结构模型代表集的原则是优先选择覆盖率高的，覆盖率一致再考虑序列相似性，次之考虑模型质量评分(QMEAN,MPSQ,DOPE)，再次考虑Loop-Content，而后加入代表集。与确定晶体结构代表集时同理，接着看排名次位的链，若overlap/该链长度 <= 0.3 且 (该链长度-overlap)/已加入代表集的链的非冗余长度 >= 0.2，则将该链加入代表集，循环此操作。

单体结构模型与相互作用模型分别处理。单体结构模型中优先选择来自SMR的模型。相互作用模型中，确定好代表集后还需计算模型的范德华能量以判断模型的合理性，如若不合理则需要重新确定该Uniprot Entry下的相互作用结构模型代表集。
