# UniProt ID Mapping
> Last modified time 2019-10-24

## Example

### InputFile
* ```group_df```

|     CHROM     |      POS      |      ID       |      REF      |      ALT      |     QUAL      |    FILTER     |      MUT      |      DNA      |     PROT      |      DB       |     PHEN      |   RANKSCORE   |     GENE      |     CLASS     |    STRAND     | mutation_unp  |   missense    |RefSeq_protein |RefSeq_nucleotide|
|:-------------:|:-------------:|:-------------:|:-------------:|:-------------:|:-------------:|:-------------:|:-------------:|:-------------:|:-------------:|:-------------:|:-------------:|:-------------:|:-------------:|:-------------:|:-------------:|:-------------:|:-------------:|:-------------:|:-------------:|
|       1       |    930215     |   CM1613956   |       A       |      [G]      |      NaN      |      NaN      |      ALT      |NM_152486.2:c.133A>G|NP_689699.2:p.K45E|  rs903331232  |"Retinitis_pigmentosa"|     0.21      |    SAMD11     |      DM?      |       +       |     K45E      |      yes      |   NP_689699   |   NM_152486   |
|       1       |    942143     |   CM1511864   |       C       |      [G]      |      NaN      |      NaN      |      ALT      |NM_152486.2:c.877C>G|NP_689699.2:p.P293A|  rs200195897  |"Autism_spectrum_disorder"|      0.1      |    SAMD11     |      DM?      |       +       |     P293A     |      yes      |   NP_689699   |   NM_152486   |
|       1       |    1022225    |   CM148517    |       G       |      [A]      |      NaN      |      NaN      |      ALT      |NM_198576.3:c.226G>A|NP_940978.2:p.G76S|  rs756623659  |"Congenital_myasthenic_syndrome_with_distal_mu...|     0.15      |     AGRN      |      DM       |       +       |     G76S      |      yes      |   NP_940978   |   NM_198576   |
|       1       |    1022313    |   CM148518    |       A       |      [T]      |      NaN      |      NaN      |      ALT      |NM_198576.3:c.314A>T|NP_940978.2:p.N105I|  rs879253787  |"Congenital_myasthenic_syndrome_with_distal_mu...|     0.91      |     AGRN      |      DM       |       +       |     N105I     |      yes      |   NP_940978   |   NM_198576   |
|       1       |    1041648    |   CM1613410   |       G       |      [T]      |      NaN      |      NaN      |      ALT      |NM_198576.3:c.1123G>T|NP_940978.2:p.A375S|  rs138031468  |"Ovarian_cancer_epithelial_reduced_risk"|0.12000000000000001|     AGRN      |      DP       |       +       |     A375S     |      yes      |   NP_940978   |   NM_198576   |

#### Target Columns & Info
* id_col: ```RefSeq_protein```
* muta_col: ```mutation_unp```
* gene_col: ```GENE```
* id_type: Refseq Protein -> P_REFSEQ_AC
* muta_type: mutation in UniProt Site

### Get ID Mapping File

#### Script

```py
from UniProt_unit import UniProt_unit
id_col = 'RefSeq_protein'
id_type = 'P_REFSEQ_AC'
muta_col = 'mutation_unp'
gene_col = 'GENE'
usecols = ['id', 'genes', 'reviewed', 'comment(ALTERNATIVE%20PRODUCTS)', 'organism', 'protein%20names']  # Necessary Columns
reportPath = '/data/zzf/UniProt_files/id_mapping_files/HGMD_RefSeq_protein_mapping_Report_1021.txt'  # Report of Data Processing
rawOutputPath = '/data/zzf/UniProt_files/id_mapping_files/HGMD_RefSeq_protein_mapping_1021.tsv' # OutPut File of ID Mapping (RAW)
handledOutputPath = '/data/zzf/groupWorks/HGMD_RefSeq_protein_mapping_modified_1024.tsv'  # OutPut File of ID Mapping (Final Result)
# Add constraint/filter to data
constraint_dict = {
    "GENE_status": (False, "ne"),  # 'ne' for !=
    "Status": ("reviewed", "eq"),  # 'eq' for ==
    "unp_map_tage": ("Untrusted & No Isoform", "ne")
}

# Initial
unp_demo = UniProt_unit(group_df, id_col, id_type, usecols, reportPath, muta_col=muta_col, gene_col=gene_col)
# Return True if get RAW Result Successfully
unp_demo.get_raw_ID_Mapping(rawOutputPath)
# Deal with different situations
handled_df = unp_demo.handle_ID_Mapping()
# Add Gene Status
unp_demo.getGeneStatus(handled_df)
# Label Mapping Status
unp_demo.label_mapping_status(handled_df, constraint_dict)
# close the file-handle of report
unp_demo.report.close()
# Output the final result
handled_df.to_csv(handledOutputPath, sep='\t', index=False)
```

#### About Label/Report
##### About ```unp_map_tage```

* ```Untrusted & No Isoform```

是指UniProt存在Isoform但是Mapping结果没有明确给出是Map上哪条Isoform,转录本序列与蛋白序列不一致

It means that there are isoforms in UniProt, but the mapping result does not clearly indicate which isoform is correspond with the transcript(e.g), and the transcript sequence is inconsistent with the protein sequence.

* ```Trusted & No Isoform```
是指UniProt不存在Isoform,Mapping结果没问题
* ```Trusted & Isoform```
是指UniProt存在Isoform,Mapping结果没问题

##### About ```Mapping_status```
* ```Yes```: 可信的结果，进行后续的PDB Mapping; 通过```constraint_dict```的限制
* ```Error```: 一个id对应多个UniProt; 通过```constraint_dict```的限制
* ```No```: 不可信的结果; 未通过```constraint_dict```的限制

#### RAW ID Mapping File
##### with no isomap data
|     Entry     |  Gene names   |    Status     |Alternative products (isoforms)|   Organism    | Protein names |   yourlist    |    isomap     |
|:-------------:|:-------------:|:-------------:|:-------------:|:-------------:|:-------------:|:-------------:|:-------------:|
|    P43489     |TNFRSF4 TXGP1L |   reviewed    |      NaN      |Homo sapiens (Human)|Tumor necrosis factor receptor superfamily mem...|   NP_003318   |      NaN      |
|  A0A024R084   |SDF4 hCG_19193 |  unreviewed   |      NaN      |Homo sapiens (Human)|Stromal cell derived factor 4, isoform CRA_c|   NP_057260   |      NaN      |


##### with isomap data
|     Entry     |  Gene names   |    Status     |Alternative products (isoforms)|   Organism    | Protein names |   yourlist    |    isomap     |
|:-------------:|:-------------:|:-------------:|:-------------:|:-------------:|:-------------:|:-------------:|:-------------:|
|    O00468     |  AGRN AGRIN   |   reviewed    |ALTERNATIVE PRODUCTS:  Event=Alternative splic...|Homo sapiens (Human)|Agrin [Cleaved into: Agrin N-terminal 110 kDa ...|   NP_940978   |NP_940978 -> O00468-6|
|    Q9BRK5     |SDF4 CAB45 PSEC0034|   reviewed    |ALTERNATIVE PRODUCTS:  Event=Alternative splic...|Homo sapiens (Human)|45 kDa calcium-binding protein (Cab45) (Stroma...|   NP_057260   |NP_057260 -> Q9BRK5-1|

##### with no isomap data & need to be split
|     Entry     |  Gene names   |    Status     |Alternative products (isoforms)|   Organism    | Protein names |   yourlist    |    isomap     |
|:-------------:|:-------------:|:-------------:|:-------------:|:-------------:|:-------------:|:-------------:|:-------------:|
|    O75015     |FCGR3B CD16B FCG3 FCGR3 IGFR3|   reviewed    |      NaN      |Homo sapiens (Human)|Low affinity immunoglobulin gamma Fc region re...|NP_000561,NP_001231682|      NaN      |
|    Q53SH4     |CHRNA1 hCG_1811440|  unreviewed   |      NaN      |Homo sapiens (Human)|Cholinergic receptor, nicotinic, alpha 1 (Musc...|NP_000070,NP_001034612|      NaN      |

##### with isomap data & need to be split
|     Entry     |  Gene names   |    Status     |Alternative products (isoforms)|   Organism    | Protein names |   yourlist    |    isomap     |
|:-------------:|:-------------:|:-------------:|:-------------:|:-------------:|:-------------:|:-------------:|:-------------:|
|    O94827     |PLEKHG5 KIAA0720 |   reviewed    |ALTERNATIVE PRODUCTS:  Event=Alternative splic...|Homo sapiens (Human)|Pleckstrin homology domain-containing family G...|NP_065682,NP_001252521|NP_001252521 -> O94827-6,NP_065682 -> O94827-5|
|    O60333     |KIF1B KIAA0591 KIAA1448|   reviewed    |ALTERNATIVE PRODUCTS:  Event=Alternative splic...|Homo sapiens (Human)|Kinesin-like protein KIF1B (Klp)|NP_055889,NP_904325|NP_055889 -> O60333-2,NP_904325 -> O60333-3|

##### Untrusted Results
Example 1:

|     Entry     |  Gene names   |    Status     |Alternative products (isoforms)|   Organism    | Protein names |   yourlist    |    isomap     |
|:-------------:|:-------------:|:-------------:|:-------------:|:-------------:|:-------------:|:-------------:|:-------------:|
|    Q15025     |TNIP1 KIAA0113 NAF1|   reviewed    |ALTERNATIVE PRODUCTS:  Event=Alternative splic...|Homo sapiens (Human)|TNFAIP3-interacting protein 1 (A20-binding inh...|NP_006049,~~NP_001239314~~|NP_006049 -> Q15025-1|
|    P63092     |GNAS GNAS1 GSP |   reviewed    |ALTERNATIVE PRODUCTS:  Event=Alternative splic...|Homo sapiens (Human)|Guanine nucleotide-binding protein G(s) subuni...|~~NP_001070958~~,NP_000507,~~NP_001296769~~|NP_000507 -> P63092-1|

* In the example 1, ```NP_001239314, NP_001070958, NP_001296769``` have no corresponding isoform, but the corresponding UniProt Entry has.
* Here may be the reason, there is no identical sequence in isoforms:

```clustal
CLUSTAL O(1.2.4) multiple sequence alignment


NP_001239314.1               MEGRGPYRIYDPGGSVPSGEASAAFERLVKENSRLKEKMQGIKMLGELLEESQMEATRLR
NP_006049.3                  MEGRGPYRIYDPGGSVPSGEASAAFERLVKENSRLKEKMQGIKMLGELLEESQMEATRLR
sp|Q15025|TNIP1_HUMAN        MEGRGPYRIYDPGGSVPSGEASAAFERLVKENSRLKEKMQGIKMLGELLEESQMEATRLR
sp|Q15025-2|TNIP1_HUMAN      MEGRGPYRIYDPGGSVPSGEASAAFERLVKENSRLKEKMQGIKMLGELLEESQMEATRLR
sp|Q15025-3|TNIP1_HUMAN      -----------------------------------------------------MEATRLR
sp|Q15025-4|TNIP1_HUMAN      MEGRGPYRIYDPGGSVPSGEASAAFERLVKENSRLKEKMQGIKMLGELLEESQMEATRLR
sp|Q15025-5|TNIP1_HUMAN      MEGRGPYRIYDPGGSVPSGEASAAFERLVKENSRLKEKMQGIKMLGELLEESQMEATRLR
sp|Q15025-6|TNIP1_HUMAN      MEGRGPYRIYDPGGSVPSGEASAAFERLVKENSRLKEKMQGIKMLGELLEESQMEATRLR
sp|Q15025-7|TNIP1_HUMAN      MEGRGPYRIYDPGGSVPSGEASAAFERLVKENSRLKEKMQGIKMLGELLEESQMEATRLR
sp|Q15025-8|TNIP1_HUMAN      MEGRGPYRIYDPGGSVPSGEASAAFERLVKENSRLKEKMQGIKMLGELLEESQMEATRLR
                                                                                  *******

NP_001239314.1               QKAEELVKDNELLPPPSPSLGSFDPLAELTGKDSNVTASPTAPACPSDKPAPVQKPPSSG
NP_006049.3                  QKAEELVKDNELLPPPSPSLGSFDPLAELTGKDSNVTASPTAPACPSDKPAPVQKPPSSG
sp|Q15025|TNIP1_HUMAN        QKAEELVKDNELLPPPSPSLGSFDPLAELTGKDSNVTASPTAPACPSDKPAPVQKPPSSG
sp|Q15025-2|TNIP1_HUMAN      QKAEELVKDNELLPPPSPSLGSFDPLAELTGKDSNVTASPTAPACPSDKPAPVQKPPSSG
sp|Q15025-3|TNIP1_HUMAN      QKAEELVKDNELLPPPSPSLGSFDPLAELTGKDSNVTASPTAPACPSDKPAPVQKPPSSG
sp|Q15025-4|TNIP1_HUMAN      QKAEELVKDNELLPPPSPSLGSFDPLAELTGKDSNVTASPTAPACPSDKPAPVQKPPSSG
sp|Q15025-5|TNIP1_HUMAN      QKAEELVKDNELLPPPSPSLGSFDPLAELTGKDSNVTASPTAPACPSDKPAPVQKPPSSG
sp|Q15025-6|TNIP1_HUMAN      QKAEELVKDNELLPPPSPSLGSFDPLAELTGKDSNVTASPTAPACPSDKPAPVQKPPSSG
sp|Q15025-7|TNIP1_HUMAN      QKAEELVKDNELLPPPSPSLGSFDPLAELTGKDSNVTASPTAPACPSDKPAPVQKPPSSG
sp|Q15025-8|TNIP1_HUMAN      QKAEELVKDNELLPPPSPSLGSFDPLAELTGKDSNVTASPTAPACPSDKPAPVQKPPSSG
                             ************************************************************
........................

NP_001239314.1               ERMNEEKEELKKQVEKLQAQVTLSNAQLKAFKDEEKAREALRQQKRKAKASGERYHVEPH
NP_006049.3                  ERMNEEKEELKKQVEKLQAQVTLSNAQLKAFKDEEKAREALRQQKRKAKASGERYHVEPH
sp|Q15025|TNIP1_HUMAN        ERMNEEKEELKKQVEKLQAQVTLSNAQLKAFKDEEKAREALRQQKRKAKASGERYHVEPH
sp|Q15025-2|TNIP1_HUMAN      ERMNEEKEELKKQVEKLQAQVTLSNAQLKAFKDEEKAREALRQQKRKAKASGERYHVEPH
sp|Q15025-3|TNIP1_HUMAN      ERMNEEKEELKKQVEKLQAQVTLSNAQLKAFKDEEKAREALRQQKRKAKASGERYHVEPH
sp|Q15025-4|TNIP1_HUMAN      ERMNEEKEELKKQVEKLQAQVTLSNAQLKAFKDEEKAREALRQQKRKA------------
sp|Q15025-5|TNIP1_HUMAN      ERMNEEKEELKKQVEKLQAQVTLSNAQLKAFKDEEKAREALRQQKRKAKSLQKMTVRGLS
sp|Q15025-6|TNIP1_HUMAN      ERMNEEKEELKKQVEKLQAQVTLSNAQLKAFKDEEKAREALRQQKRKA------------
sp|Q15025-7|TNIP1_HUMAN      ERMNEEKEELKKQVEKLQAQVTLSNAQLKAFKDEEKAREALRQQKRKAKGTHRGCPRRLP
sp|Q15025-8|TNIP1_HUMAN      ERMNEEKEELKKQVEKLQAQVTLSNAQLKAFKDEEKAREALRQQKRKAKSQLISDCQ-ET
                             ************************************************

NP_001239314.1               PEHLCGAYPYAYPPMPAMVPHHGFEDWSQIRYPPPPMAMEHPPPLPNSRLFHLHQFCRSR
NP_006049.3                  PEHLCGAYPYAYPPMPAMVPHHGFEDWSQIRYPPPPMAMEHPPPLPNSRLFHLPEYTWRL
sp|Q15025|TNIP1_HUMAN        PEHLCGAYPYAYPPMPAMVPHHGFEDWSQIRYPPPPMAMEHPPPLPNSRLFHLPEYTWRL
sp|Q15025-2|TNIP1_HUMAN      PEHLCGAYPYAYPPMPAMVPHHGFEDWSQIRYPPPPMAMEHPPPLPNSRLFHLPEYTWRL
sp|Q15025-3|TNIP1_HUMAN      PEHLCGAYPYAYPPMPAMVPHHGFEDWSQIRYPPPPMAMEHPPPLPNSRLFHLPEYTWRL
sp|Q15025-4|TNIP1_HUMAN      ----------------------------------------------------KPEYTWRL
sp|Q15025-5|TNIP1_HUMAN      ETRLCHLAPPSSCRAS--------------------------------------------
sp|Q15025-6|TNIP1_HUMAN      ----------------------------------------------------KPEYTWRL
sp|Q15025-7|TNIP1_HUMAN      ERKVK-------------------------------------------------------
sp|Q15025-8|TNIP1_HUMAN      RSHLHGVARASAG-----------------------------------------------


NP_001239314.1               NTPGVYPVEGFEIQIRAPK-----------------
NP_006049.3                  PCGGVRNPNQSSQVMDPPTARPTEPESPKNDREGPQ
sp|Q15025|TNIP1_HUMAN        PCGGVRNPNQSSQVMDPPTARPTEPESPKNDREGPQ
sp|Q15025-2|TNIP1_HUMAN      PCGGVRNPNQSSQVMDPPTARPTEPEPADLRLPRN-
sp|Q15025-3|TNIP1_HUMAN      PCGGVRNPNQSSQVMDPPTARPTEPESPKNDREGPQ
sp|Q15025-4|TNIP1_HUMAN      PCGGVRNPNQSSQVMDPPTARPTEPESPKNDREGPQ
sp|Q15025-5|TNIP1_HUMAN      ------------------------------------
sp|Q15025-6|TNIP1_HUMAN      PCGGVRNPNQSSQVMDPPTARPTEPEPADLRLPRN-
sp|Q15025-7|TNIP1_HUMAN      ------------------------------------
sp|Q15025-8|TNIP1_HUMAN      ------------------------------------
```

Example 2:

|     Entry     |  Gene names   |    Status     |Alternative products (isoforms)|   Organism    | Protein names |   yourlist    |    isomap     |
|:-------------:|:-------------:|:-------------:|:-------------:|:-------------:|:-------------:|:-------------:|:-------------:|
|    Q96NU1     |    SAMD11     |   reviewed    |ALTERNATIVE PRODUCTS:  Event=Alternative promo...|Homo sapiens (Human)|Sterile alpha motif domain-containing protein ...|   ~~NP_689699~~   |      NaN      |
|    Q9C0B2     |CFAP74 C1orf222 KIAA1751|   reviewed    |ALTERNATIVE PRODUCTS:  Event=Alternative splic...|Homo sapiens (Human)|Cilia- and flagella-associated protein 74| ~~NP_001291289~~  |      NaN      |

```clustal
CLUSTAL O(1.2.4) multiple sequence alignment


NP_689699.2                  MSKGILQVHPPICDCPGCRISSPVNRGRLADKRTVALPAARNLKKERTPSFSASDGDSDG
sp|Q96NU1|SAM11_HUMAN        MSKGILQVHPPICDCPGCRISSPVNRGRLADKRTVALPAARNLKKERTPSFSASDGDSDG
sp|Q96NU1-1|SAM11_HUMAN      ------------------------------------------------------------
sp|Q96NU1-2|SAM11_HUMAN      ------------------------------------------------------------
sp|Q96NU1-4|SAM11_HUMAN      MSKGILQVHPPICDCPGCRISSPVNRGRLADKRTVALPAARNLKKERTPSFSASDGDSDG
sp|Q96NU1-5|SAM11_HUMAN      MSKGILQVHPPICDCPGCRISSPVNRGRLADKRTVALPAARNLKKERTPSFSASDGDSDG
sp|Q96NU1-6|SAM11_HUMAN      MSKGILQVHPPICDCPGCRISSPVNRGRLADKRTVALPAARNLKKERTPSFSASDGDSDG


NP_689699.2                  SGPTCGRRPGLKQEDGPHIRIMKRRVHTHWDVNISFREASCSQDGNLPTLISSVHRSRHL
sp|Q96NU1|SAM11_HUMAN        SGPTCGRRPGLKQEDGPHIRIMKRRVHTHWDVNISFREASCSQDGNLPTLISSVHRSRHL
sp|Q96NU1-1|SAM11_HUMAN      ------------------------------------------------------------
sp|Q96NU1-2|SAM11_HUMAN      ------------------------------------------------------------
sp|Q96NU1-4|SAM11_HUMAN      SGPTCGRRPGLKQEDGPHIRIMKRRVHTHWDVNISFREASCSQDGNLPTLISSVHRSRHL
sp|Q96NU1-5|SAM11_HUMAN      SGPTCGRRPGLKQEDGPHIRIMKRRVHTHWDVNISFREASCSQDGNLPTLISSVHRSRHL
sp|Q96NU1-6|SAM11_HUMAN      SGPTCGRRPGLKQEDGPHIRIMKRRVHTHWDVNISFREASCSQDGNLPTLISSVHRSRHL


NP_689699.2                  VMPEHQSRCEFQRGSLEIGLRPAGDLLGKRLGRSPRISSDCFSEKRARSESPQ-EALLLP
sp|Q96NU1|SAM11_HUMAN        VMPEHQSRCEFQRGSLEIGLRPAGDLLGKRLGRSPRISSDCFSEKRARSESPQ-EALLLP
sp|Q96NU1-1|SAM11_HUMAN      ------------------------------------------------------------
sp|Q96NU1-2|SAM11_HUMAN      ------------------------------------------------------------
sp|Q96NU1-4|SAM11_HUMAN      VMPEHQSRCEFQRGSLEIGLRPAGDLLGKRLGRSPRISSDCFSEKRARSESPQAEALLLP
sp|Q96NU1-5|SAM11_HUMAN      VMPEHQSRCEFQRGSLEIGLRPAGDLLGKRLGRSPRISSDCFSEKRARSESPQAEALLLP
sp|Q96NU1-6|SAM11_HUMAN      VMPEHQSRCEFQRGSLEIGLRPAGDLLGKRLGRSPRISSDCFSEKRARSESPQ-EALLLP


NP_689699.2                  RELGPSMAPEDHYRRLVSALSEASTFEDPQRLYHLGLPSHGEDPPWHDPPHHLPSHDLLR
sp|Q96NU1|SAM11_HUMAN        RELGPSMAPEDHYRRLVSALSEASTFEDPQRLYHLGLPSHGEDPPWHDPPHHLPSHDLLR
sp|Q96NU1-1|SAM11_HUMAN      ------MAPEDHYRRLVSALSEASTFEDPQRLYHLGLP----------------SHDLLR
sp|Q96NU1-2|SAM11_HUMAN      ------MAPEDHYRRLVSALSEASTFEDPQRLYHLGLPSHGEDPPWHDPPHHLPSHDLLR
sp|Q96NU1-4|SAM11_HUMAN      RELGPSMAPEDHYRRLVSALSEASTFEDPQRLYHLGLP----------------SHDLLR
sp|Q96NU1-5|SAM11_HUMAN      RELGPSMAPEDHYRRLVSALSEASTFEDPQRLYHLGLP----------------SHDLLR
sp|Q96NU1-6|SAM11_HUMAN      RELGPSMAPEDHYRRLVSALSEASTFEDPQRLYHLGLP----------------SHDLLR
                                   ********************************                ******

NP_689699.2                  VRQEVAAAALRGPSGLEAHLPSSTAGQRRKQGLAQHREGAAPAAAPSFSERELPQPPPLL
sp|Q96NU1|SAM11_HUMAN        VRQEVAAAALRGPSGLEAHLPSSTAGQRRKQGLAQHREGAAPAAAPSFSERELPQPPPLL
sp|Q96NU1-1|SAM11_HUMAN      VRQEVAAAALRGPSGLEAHLPSSTAGQRRKQGLAQHREGAAPAAAPSFSERELPQPPPLL
sp|Q96NU1-2|SAM11_HUMAN      VRQEVAAAALRGPSGLEAHLPSSTAGQRRKQGLAQHREGAAPAAAPSFSERELPQPPPLL
sp|Q96NU1-4|SAM11_HUMAN      VRQEVAAAALRGPSGLEAHLPSSTAGQRRKQGLAQHREGAAPAAAPSFSERELPQPPPLL
sp|Q96NU1-5|SAM11_HUMAN      VRQEVAAAALRGPSGLEAHLPSSTAGQRRKQGLAQHREGAAPAAAPSFSERELPQPPPLL
sp|Q96NU1-6|SAM11_HUMAN      VRQEVAAAALRGPSGLEAHLPSSTAGQRRKQGLAQHREGAAPAAAPSFSERELPQPPPLL
                             ************************************************************

NP_689699.2                  SPQNAPHVALGPHLRPPFLGVPSALCQTPGYGFLPPAQAEMFARQQELLRKQNLARLELP # Notice that there is little difference in residues
sp|Q96NU1|SAM11_HUMAN        SPQNAPHVALGPHLRPPFLGVPSALCQTPGYGFLPPAQAEMFAWQQELLRKQNLARLELP
sp|Q96NU1-1|SAM11_HUMAN      SPQNAPHVALGPHLRPPFLGVPSALCQTPGYGFLPPAQAEMFAWQQELLRKQNLARLELP
sp|Q96NU1-2|SAM11_HUMAN      SPQNAPHVALGPHLRPPFLGVPSALCQTPGYGFLPPAQAEMFAWQQELLRKQNLARLELP
sp|Q96NU1-4|SAM11_HUMAN      SPQNAPHVALGPHLRPPFLGVPSALCQTPGYGFLPPAQAEMFAWQQELLRKQNLARLELP
sp|Q96NU1-5|SAM11_HUMAN      SPQNAPHVALGPHLRPPFLGVPSALCQTPGYGFLPPAQAEMFAWQQELLRKQNLARLELP
sp|Q96NU1-6|SAM11_HUMAN      SPQNAPHVALGPHLRPPFLGVPSALCQTPGYGFLPPAQAEMFAWQQELLRKQNLARLELP
                             ******************************************* ****************

NP_689699.2                  ADLLRQKELESARPQLLAPETALRPNDGAEELQRRGALLVLNHGAAPLLALPPQGPPGSG
sp|Q96NU1|SAM11_HUMAN        ADLLRQKELESARPQLLAPETALRPNDGAEELQRRGALLVLNHGAAPLLALPPQGPPGSG
sp|Q96NU1-1|SAM11_HUMAN      ADLLRQKELESARPQLLAPETALRPNDGAEELQRRGALLVLNHGAAPLLALPPQGPPGSG
sp|Q96NU1-2|SAM11_HUMAN      ADLLRQKELESARPQLLAPETALRPNDGAEELQRRGALLVLNHGAAPLLALPPQGPPGSG
sp|Q96NU1-4|SAM11_HUMAN      ADLLRQKELESARPQLLAPETALRPNDGAEELQRRGALLVLNHGAAPLLALPPQGPPGSG
sp|Q96NU1-5|SAM11_HUMAN      ADLLRQKELESARPQLLAPETALRPNDGAEELQRRGALLVLNHGAAPLLALPPQGPPGSG
sp|Q96NU1-6|SAM11_HUMAN      ADLLRQKELESARPQLLAPETALRPNDGAEELQRRGALLVLNHGAAPLLALPPQGPPGSG
                             ************************************************************

NP_689699.2                  PPTPSRDSARRAPRKGGPGPASARPSESKEMTGARLWAQDGSEDEPPKDSDGEDPETAAV
sp|Q96NU1|SAM11_HUMAN        PPTPSRDSARRAPRKGGPGPASARPSESKEMTGARLWAQDGSEDEPPKDSDGEDPETAAV
sp|Q96NU1-1|SAM11_HUMAN      PPTPSRDSARRAPRKGGPGPASARPSESKEMTGARLWAQDGSEDEPPKDSDGEDPETAAV
sp|Q96NU1-2|SAM11_HUMAN      PPTPSRDSARRAPRKGGPGPASARPSESKEMTGARLWAQDGSEDEPPKDSDGEDPETAAV
sp|Q96NU1-4|SAM11_HUMAN      PPTPSRDSARRAPRKGGPGPASARPSESKEMTGARLWAQDGSEDEPPKDSDGEDPETAAV
sp|Q96NU1-5|SAM11_HUMAN      PPTPSRDSARRAPRKGGPGPASARPSESKEMTGARLWAQDGSEDEPPKDSDGEDPETAAV
sp|Q96NU1-6|SAM11_HUMAN      PPTPSRDSARRAPRKGGPGPASARPSESKEMTGARLWAQDGSEDEPPKDSDGEDPETAAV
                             ************************************************************

NP_689699.2                  GCRGPTPGQAPAGGAGAEGKGLFPGSTLPLGFPYAVSPYFHTGAVGGLSMDGEEAPAPED
sp|Q96NU1|SAM11_HUMAN        GCRGPTPGQAPAGGAGAEGKGLFPGSTLPLGFPYAVSPYFHTGAVGGLSMDGEEAPAPED
sp|Q96NU1-1|SAM11_HUMAN      GCRGPTPGQAPAGGAGAEGKGLFPGSTLPLGFPYAVSPYFHTGAVGGLSMDGEEAPAPED
sp|Q96NU1-2|SAM11_HUMAN      GCRGPTPGQAPAGGAGAEGKGLFPGSTLPLGFPYAVSPYFHTGAVGGLSMDGEEAPAPED
sp|Q96NU1-4|SAM11_HUMAN      GCRGPTPGQAPAGGAGAEGKGLFPGSTLPLGFPYAVSPYFHTGAVGGLSMDGEEAPAPED
sp|Q96NU1-5|SAM11_HUMAN      GCRGPTPGQAPAGGAGAEGKGLFPGSTLPLGFPYAVSPYFHTGAVGGLSMDGEEAPAPED
sp|Q96NU1-6|SAM11_HUMAN      GCRGPTPGQAPAGGAGAEGKGLFPGSTLPLGFPYAVSPYFHTGAVGGLSMDGEEAPAPED
                             ************************************************************

NP_689699.2                  VTKWTVDDVCSFVGGLSGCGEYTRVFREQGIDGETLPLLTEEHLLTNMGLKLGPALKIRA
sp|Q96NU1|SAM11_HUMAN        VTKWTVDDVCSFVGGLSGCGEYTRVFREQGIDGETLPLLTEEHLLTNMGLKLGPALKIRA
sp|Q96NU1-1|SAM11_HUMAN      VTKWTVDDVCSFVGGLSGCGEYTRVFREQGIDGETLPLLTEEHLLTNMGLKLGPALKIRA
sp|Q96NU1-2|SAM11_HUMAN      VTKWTVDDVCSFVGGLSGCGEYTRVFREQGIDGETLPLLTEEHLLTNMGLKLGPALKIRA
sp|Q96NU1-4|SAM11_HUMAN      VTKWTVDDVCSFVGGLSGCGEYTRVFREQGIDGETLPLLTEEHLLTNMGLKLGPALKIRA
sp|Q96NU1-5|SAM11_HUMAN      VTKWTVDDVCSFVGGLSGCGEYTRVFREQGIDGETLPLLTEEHLLTNMGLKLGPALKIRA
sp|Q96NU1-6|SAM11_HUMAN      VTKWTVDDVCSFVGGLSGCGEYTRVFREQGIDGETLPLLTEEHLLTNMGLKLGPALKIRA
                             ************************************************************

NP_689699.2                  Q---------------------------------VARRLGRVFYVASFPVALPLQPPTLR
sp|Q96NU1|SAM11_HUMAN        Q---------------------------------VARRLGRVFYVASFPVALPLQPPTLR
sp|Q96NU1-1|SAM11_HUMAN      Q---------------------------------VARRLGRVFYVASFPVALPLQPPTLR
sp|Q96NU1-2|SAM11_HUMAN      QVRRWGVRSGSPDHSWAESSGWVCDSPHQAISLQVARRLGRVFYVASFPVALPLQPPTLR
sp|Q96NU1-4|SAM11_HUMAN      Q---------------------------------VARRLGRVFYVASFPVALPLQPPTLR
sp|Q96NU1-5|SAM11_HUMAN      QVRRWGVRSGSPDHSWAESSGWVCDSPHQAISLQVARRLGRVFYVASFPVALPLQPPTLR
sp|Q96NU1-6|SAM11_HUMAN      Q---------------------------------VARRLGRVFYVASFPVALPLQPPTLR
                             *                                 **************************

NP_689699.2                  APERELGTGEQPLSPTTATSPYGGGHALAGQTSPKQENGTLALLPGAPDPSQPLC
sp|Q96NU1|SAM11_HUMAN        APERELGTGEQPLSPTTATSPYGGGHALAGQTSPKQENGTLALLPGAPDPSQPLC
sp|Q96NU1-1|SAM11_HUMAN      APERELGTGEQPLSPTTATSPYGGGHALAGQTSPKQENGTLALLPGAPDPSQPLC
sp|Q96NU1-2|SAM11_HUMAN      APERELGTGEQPLSPTTATSPYGGGHALAGQTSPKQENGTLALLPGAPDPSQPLC
sp|Q96NU1-4|SAM11_HUMAN      APERELGTGEQPLSPTTATSPYGGGHALAGQTSPKQENGTLALLPGAPDPSQPLC
sp|Q96NU1-5|SAM11_HUMAN      APERELGTGEQPLSPTTATSPYGGGHALAGQTSPKQENGTLALLPGAPDPSQPLC
sp|Q96NU1-6|SAM11_HUMAN      APERELGTGEQPLSPTTATSPYGGGHALAGQTSPKQENGTLALLPGAPDPSQPLC
                             *******************************************************
```

##### Error Results
|     Entry     |  Gene names   |    Status     |Alternative products (isoforms)|   Organism    | Protein names |canonical_isoform| unp_map_tage  |   yourlist    |    UniProt    |     GENE      |  GENE_status  |Mapping_status |
|:-------------:|:-------------:|:-------------:|:-------------:|:-------------:|:-------------:|:-------------:|:-------------:|:-------------:|:-------------:|:-------------:|:-------------:|:-------------:|
|    P01911     |HLA-DRB1 HLA-DRB2|   reviewed    |      NaN      |Homo sapiens (Human)|HLA class II histocompatibility antigen, DRB1-...|      NaN      |Trusted &amp; No Isoform|   NP_002115   |    P01911     |   HLA-DRB1    |     True      |     Error     |
|    Q29974     |   HLA-DRB1    |   reviewed    |      NaN      |Homo sapiens (Human)|HLA class II histocompatibility antigen, DRB1-...|      NaN      |Trusted &amp; No Isoform|   NP_002115   |    Q29974     |   HLA-DRB1    |     True      |     Error     |

#### Report Statistis
* all RefSeq Count: 8450
* Unmapped Results: 34
* Trusted Results: 8066
  * ```handled_df[handled_df['Mapping_status'] == 'Yes']```
* Untrusted Resultes: 350
  * Error Results: 1 ```{NP_002115: [P01911, Q29974]}, handled_df[handled_df['Mapping_status'] == 'Error']```
  * Others: ```handled_df[handled_df['Mapping_status'] == 'No']```

### Deal With Mutation Info
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
NP_000016                           [S257P, S165P, W64R, I62M]
NP_000017    [A2V, A3P, A3V, S23R, P24L, M26L, I72V, E80D, ...
NP_000018    [C306R, G302R, T257I, G252E, G252R, C163S, R16...
NP_000020    [R458C, V385E, Q285P, Y281C, M268T, L244R, T24...
NP_000021    [T9N, P11R, P11L, L18F, N22S, L25R, L26P, P28S...
NP_000022    [V275M, A274T, R240W, V153M, G133R, C132R, E89...
NP_000023    [Y586F, R572H, S568G, M567I, M567V, E565K, H56...
NP_000024    [T14R, A19S, M67V, R74W, P84S, P84L, C88W, E90...
NP_000025                         [D129G, E207K, A280V, C339Y]
NP_000026    [Y343H, A338V, N335K, L311P, R304Q, R304W, L28...
NP_000027    [M712T, P684S, S626C, P572S, R458H, R421W, M34...
NP_000028    [D1592N, A1462V, L1340P, T1075I, I1055T, L1046...
NP_000029    [N32S, E56Q, S92C, R99W, R106H, S127G, S130G, ...
NP_000030    [K262N, E222K, L202H, L202P, L202R, A199P, L19...
NP_000031           [A43T, D45V, V50M, Q58K, D65N, K78E, T94A]
NP_000032    [E21K, E31K, E37K, R43C, L46P, Q64H, K90E, A12...
NP_000033                          [W335S, C325G, L266V, K38E]
NP_000034    [L7P, L8P, A16T, T28A, G66D, D78E, C82R, C85R,...
NP_000035    [E2K, V30M, A45G, Q58L, Q71R, S176R, C177G, Q1...
                                   ...
NP_997253                                              [T102I]
NP_997254                                              [A867V]
NP_997260                                               [W84L]
NP_997279                                              [N573S]
NP_997281                                              [M315V]
NP_997297                                              [V897M]
NP_997304    [L37F, K145E, H211Q, P289L, R301Q, D349G, L375...
NP_997320    [H827L, M1518V, A1616T, V1896M, Y1989C, R2218H...
NP_997329                                              [D155N]
NP_997351                                               [R73Q]
NP_997400                                             [T1660I]
NP_997464    [W4R, I26M, R50C, G79R, R94C, G114A, R133H, Y1...
NP_997468                                   [S8F, F10V, I196L]
NP_997646                                        [R55L, G120V]
NP_997647    [V2046I, L1974P, L1974R, W1925C, W1925R, R1784...
NP_997657                                              [R492Q]
NP_997698                      [V1845M, P1635L, V1453F, M255K]
NP_997700                                        [G92V, A514S]
NP_997717                                   [R99W, H96R, I67N]
NP_998760                                         [G50E, G96V]
NP_998761                                              [D473N]
NP_998763                                              [P226A]
NP_998764    [N52S, F54S, R57W, R58W, D81G, V87I, D93E, L10...
NP_998771    [S453R, G437S, S406C, P395Q, N333D, G205V, S15...
NP_998772                                        [R35P, H154P]
NP_998778                   [S358L, L348P, T185M, M132T, A56T]
NP_998813                                        [P19T, R191M]
NP_998818    [C321W, G320V, R288W, I287S, I281T, V274M, G25...
NP_998820    [S363P, E401K, S661F, V690G, R718C, G975R, S13...
NP_998885    [E127K, C122R, Q108P, P96T, P80L, G66V, G66S, ...
```
