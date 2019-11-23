# @Date:   2019-11-20T23:30:02+08:00
# @Email:  1730416009@stu.suda.edu.cn
# @Filename: Run.py
# @Last modified time: 2019-11-23T21:23:01+08:00
import click
import configparser
import os
import json
from pandas import read_csv, merge, Series
from numpy import nan
from ProcessSIFTS import RetrieveSIFTS, handle_SIFTS, deal_with_insertionDeletion_SIFTS, update_range_SIFTS, map_muta_from_unp_to_pdb
from ProcessMMCIF import MMCIF2Dfrm
from ProcessUniProt import MapUniProtID, retrieveUniProtSeq
from Utils.Logger import RunningLogger


_FOLDER = ""
_config = configparser.RawConfigParser()
_config.read('settings.ini', encoding='utf8')
_LOGGER_PATH = _config.get("DEFAULT", "LOGGER_PATH")
_SIFTS_RAW_PATH = _config.get("DEFAULT", "SIFTS_RAW_PATH")
_SIFTS_MODIFIED_PATH = _config.get("DEFAULT", "SIFTS_MODIFIED_PATH")
_SIFTS_PDB = _config.get("DEFAULT", "SIFTS_PDB")
_MMCIF_RAW_PATH = _config.get("DEFAULT", "MMCIF_RAW_PATH")
_MMCIF_MODIFIED_PATH = _config.get("DEFAULT", "MMCIF_MODIFIED_PATH")
_UniProt_ID_Mapping_RAW_PATH = _config.get("DEFAULT", "UniProt_ID_Mapping_RAW_PATH")
_UniProt_ID_Mapping_MODIFIED_PATH = _config.get("DEFAULT", "UniProt_ID_Mapping_MODIFIED_PATH")
_UniProt_DEFAULT_COL = _config.get("DEFAULT", "UniProt_DEFAULT_COL").split(",")
_SITE_INFO_PATH = _config.get("DEFAULT", "SITE_INFO_PATH")
_INTERGRATE_PATH = _config.get("DEFAULT", "INTERGRATE_PATH")
_CONVERTER = {"chain_id": str, "struct_asym_id": str, "entity_id": int, "asym_id": str}


def colorClick(a, b="Initializing %s DataSet", fg="green"):
    return click.style(b % a, fg=fg)


@click.group()
@click.option("--folder", default="", help="The file folder of new files.", type=click.Path())
def interface(folder):
    global _FOLDER, _LOGGER_PATH, _SIFTS_RAW_PATH, _SIFTS_MODIFIED_PATH, _SIFTS_PDB, _MMCIF_RAW_PATH, _MMCIF_MODIFIED_PATH, _UniProt_ID_Mapping_RAW_PATH, _UniProt_ID_Mapping_MODIFIED_PATH, _SITE_INFO_PATH, _INTERGRATE_PATH
    _FOLDER = folder
    _LOGGER_PATH = os.path.join(_FOLDER, _LOGGER_PATH)
    _SIFTS_RAW_PATH = os.path.join(_FOLDER, _SIFTS_RAW_PATH)
    _SIFTS_MODIFIED_PATH = os.path.join(_FOLDER, _SIFTS_MODIFIED_PATH)
    _SIFTS_PDB = os.path.join(_FOLDER, _SIFTS_PDB)
    _MMCIF_RAW_PATH = os.path.join(_FOLDER, _MMCIF_RAW_PATH)
    _MMCIF_MODIFIED_PATH = os.path.join(_FOLDER, _MMCIF_MODIFIED_PATH)
    _UniProt_ID_Mapping_RAW_PATH = os.path.join(_FOLDER, _UniProt_ID_Mapping_RAW_PATH)
    _UniProt_ID_Mapping_MODIFIED_PATH = os.path.join(_FOLDER, _UniProt_ID_Mapping_MODIFIED_PATH)
    _SITE_INFO_PATH = os.path.join(_FOLDER, _SITE_INFO_PATH)
    _INTERGRATE_PATH = os.path.join(_FOLDER, _INTERGRATE_PATH)


@interface.command()
@click.option("--test", default=0, help="Num of PDB IDs to test the program. Only for test.", type=int)
@click.option("--unpFile", default="", help="The file that comtains Target UniProt IDs.", type=click.Path())
@click.option("--unpCol", default="UniProt", help="The column of UniProt IDs in unpFile.", type=str)
@click.option("--sep", default="\t", help="The seperator of unpFile.", type=str)
@click.option("--filtering", default=None, help="[filterColumn filterValue]: The filter of unpFile. Keep the rows that have equal value in filter column.", type=(str, str))
@click.option("--useInitizedUnp", default=False, help="Whether to set the initialized result as the unpFile.", type=bool)
def initSIFTS(test, unpfile, unpcol, sep, filtering, useInitizedUnp):
    click.echo(colorClick("SIFTS"))
    retrieveOb = RetrieveSIFTS(
        loggingPath=_LOGGER_PATH,
        rawSIFTSpath=_SIFTS_RAW_PATH,
        downloadFolder=_FOLDER)

    if useInitizedUnp:
        unpfile = _UniProt_ID_Mapping_MODIFIED_PATH
        filtering = ("Mapping_status", "Yes")

    if unpfile != "":
        if filtering is not None:
            filtercol, filtervalue = filtering
            filterDfrm = read_csv(unpfile, usecols=[unpcol, filtercol], sep=sep)
            related_unp = filterDfrm[filterDfrm[filtercol] == filtervalue][unpcol].drop_duplicates()
        else:
            related_unp = read_csv(unpfile, usecols=[unpcol], sep=sep)[unpcol].drop_duplicates()

        related_unp = related_unp.apply(lambda x: x.split("-")[0])
    else:
        related_unp = None

    pdbs = sorted(retrieveOb.get_info_from_uniprot_pdb_file(related_unp=related_unp)['pdb_set'])

    if test:
        pdbs = pdbs[:test]

    rows, fail_list = retrieveOb.retrieve_raw_SIFTS(pdbs)
    click.echo('\nFinished Rows: %s\nFail PDBs: %s' % (rows, json.dumps(fail_list)))
    sifts_df = deal_with_insertionDeletion_SIFTS(
        sifts_df=handle_SIFTS(_SIFTS_RAW_PATH, rows),
        outputPath=_SIFTS_MODIFIED_PATH
    )
    click.echo("The Raw SIFTS File Has Been Modified")
    sifts_df[['pdb_id']].drop_duplicates().to_csv(_SIFTS_PDB, index=False)


@interface.command()
@click.option("--pdbFolder", help="The file folder of PDB repository.", type=click.Path())
@click.option("--pdbsFile", default="", help="The file that comtains PDB IDs.", type=click.Path())
@click.option("--pdbCol", default="pdb_id", help="The column of PDB IDs in pdbsFile.", type=str)
@click.option("--sep", default="\t", help="The seperator of pdbsFile.", type=str)
def initMMCIF(pdbfolder, pdbsfile, pdbcol, sep):
    click.echo(colorClick("MMCIF"))
    if pdbsfile == "":
        pdbsfile = _SIFTS_PDB
    pdbs = read_csv(pdbsfile, usecols=[pdbcol], sep=sep)[pdbcol].drop_duplicates()
    mmcif_demo = MMCIF2Dfrm(loggingPath=_LOGGER_PATH, downloadFolder=_FOLDER)
    mmcif_demo.check_mmcif_file(pdbs)
    mmcif_demo.update_mmcif_result(_MMCIF_RAW_PATH, _MMCIF_MODIFIED_PATH)


@interface.command()
@click.option("--referenceFile", default="./", help="The reference file of IDs(with mutation Site) that need to map via UniProt RESTful API.", type=click.Path(exists=True))
@click.option("--sep", default="\t", help="The seperator of referenceFile.", type=str)
@click.option("--idCol", default="RefSeq_protein", help="The column name of IDs in referenceFile.", type=str)
@click.option("--idType", default="P_REFSEQ_AC", help="ID Abbreviation that stands for the type of ID.", type=str)
@click.option("--addUseCols", default="", help="Comma-separated list of the column names for programmatic access to the UniProtKB search results.", type=str)
@click.option("--siteCol", default="mutation_unp", help="The column name of aa site in referenceFile.", type=str)
@click.option("--geneCol", default="GENE", help="The column name of gene info in referenceFile.", type=str)
@click.option("--procced/--no-procced", default=True, help="Whether to procced after saving the site info.", is_flag=True)
def initUniProt(referencefile, sep, idcol, idtype, addusecols, sitecol, genecol, procced):
    click.echo(colorClick("UniProt"))
    if addusecols != "":
        usecols = _UniProt_DEFAULT_COL + addusecols
    else:
        usecols = _UniProt_DEFAULT_COL
    # Initializing
    unp_demo = MapUniProtID(
        dfrm=read_csv(referencefile, sep=sep),
        id_col=idcol,
        id_type=idtype,
        usecols=usecols,
        loggingPath=_LOGGER_PATH,
        site_col=sitecol,
        gene_col=genecol
    )

    # Get Site Info
    if sitecol is not None:
        unp_demo.site_li.apply(json.dumps).to_csv(_SITE_INFO_PATH, sep="\t", header=["site"])
        click.echo("Site Info has been safed in %s: \n%s\n..." % (_SITE_INFO_PATH, str(unp_demo.site_li[:5])))

    if not procced:
        return
    # Return True if get RAW Result Successfully
    if unp_demo.get_raw_ID_Mapping(_UniProt_ID_Mapping_RAW_PATH):
        # Deal with different situations
        handled_ID_Mapping = unp_demo.handle_ID_Mapping()
        # Add Gene Status
        unp_demo.getGeneStatus(handled_ID_Mapping)
        # Label Mapping Status
        unp_demo.label_mapping_status(handled_ID_Mapping)
        # Output the final result
        handled_ID_Mapping.to_csv(_UniProt_ID_Mapping_MODIFIED_PATH, sep='\t', index=False)
    else:
        click.echo("There is something wrong.")


@interface.command()
@click.option("--fastaFolder", default="./", help="The file folder of UniProt FASTA Seq repository.", type=click.Path())
@click.option("--unreviewed", default=True, help="Whethter to include FASTA Seq of unreviewed UniProt Entry.", type=bool)
@click.option("--isoform", default=True, help="Whethter to include isoform Seq.", type=bool)
@click.option("--split", default=True, help="Whethter to split FASTA files.", type=bool)
@click.option("--mode", default="wget", help="Retrieve mode.", type=click.Choice(['wget', 'ftplib'], case_sensitive=False))
def initUnpFASTA(fastafolder, unreviewed, isoform, split, mode):
    click.echo(colorClick("UniProt FASTA Seq"))
    retrieveUniProtSeq(
        fastafolder,
        logger=RunningLogger("retrieveUniProtSeq", _LOGGER_PATH).logger,
        unreviewed=unreviewed,
        isoform=isoform,
        split=split,
        mode=mode
    )


@interface.command()
@click.option("--fastaFolder", default="./", help="The file folder of UniProt FASTA Seq repository.", type=click.Path())
@click.option("--siteInfoFile", default="", help="The file that comtains site info.", type=click.Path())
def mappingFromUnpToPDB(fastafolder, siteInfoFile):
    click.echo(colorClick("Mapping"))
    logger = RunningLogger("mappingFromUnpToPDB", _LOGGER_PATH).logger
    sifts_df = read_csv(_SIFTS_MODIFIED_PATH, sep="\t", converters=_CONVERTER).drop_duplicates(subset=['UniProt', 'pdb_id', 'chain_id'], keep='last')
    mmcif_df = read_csv(_MMCIF_MODIFIED_PATH, sep="\t", converters=_CONVERTER).drop_duplicates(subset=['pdb_id', 'chain_id'], keep='last')
    sifts_mmcif_df = merge(sifts_df, mmcif_df)
    sifts_mmcif_df = update_range_SIFTS(fastafolder, sifts_df=sifts_mmcif_df)

    if os.path.exists(_UniProt_ID_Mapping_MODIFIED_PATH):  # Procced
        id_map_df = read_csv(_UniProt_ID_Mapping_MODIFIED_PATH, sep="\t")
    sifts_mmcif_df = merge(sifts_mmcif_df, id_map_df[['UniProt', 'yourlist']].drop_duplicates())

    logger.warning("\n%s\n" % sifts_mmcif_df.isnull().sum())

    if siteInfoFile == "":
        siteInfoFile = _SITE_INFO_PATH

    if os.path.exists(siteInfoFile):
        siteSe = read_csv(siteInfoFile, sep="\t", index_col=0)['site']
        siteSe = siteSe.apply(json.loads)
        sifts_mmcif_df['mutation_unp'] = sifts_mmcif_df.apply(lambda x: siteSe[x['yourlist']], axis=1)
        muta_info_li = []

        sifts_mmcif_df['mutation_pdb'] = sifts_mmcif_df.apply(lambda x: map_muta_from_unp_to_pdb(
            x, 'mutation_unp', 'new_sifts_unp_range', 'new_sifts_pdb_range', muta_info_li
            ) if not isinstance(x['new_sifts_pdb_range'], float) and not isinstance(x['mutation_unp'], float) else nan, axis=1)

        sifts_mmcif_df['muta_map_info'] = Series(muta_info_li, index=sifts_mmcif_df.dropna(subset=['mutation_unp']).index)

    sifts_mmcif_df.to_csv(_INTERGRATE_PATH, sep="\t", index=False)


@interface.command()
def constraintMapping():
    intergrate_df = read_csv(_INTERGRATE_PATH, sep="\t", converters=_CONVERTER)
    filter_df = intergrate_df[
        (intergrate_df['_pdbx_coordinate_model.type'].isnull())
        & (intergrate_df['contains_unk_in_chain_pdb'] == False)
        & (intergrate_df['coordinates_len'] > 20)  # For Chain
        & (intergrate_df['method'].isin(["X-RAY DIFFRACTION", "SOLUTION NMR"]))
        & (intergrate_df['delete'] == False)
        & (intergrate_df['identity'] >= 0.9)
        & (intergrate_df['pdb_contain_chain_type'].isin(["protein", "DNA,protein", "protein,DNA", "RNA,protein", "protein,RNA"]))
    ].reset_index(drop=True)



interface.add_command(initUniProt)
interface.add_command(initUnpFASTA)
interface.add_command(initSIFTS)
interface.add_command(initMMCIF)
interface.add_command(mappingFromUnpToPDB)
interface.add_command(constraintMapping)


if __name__ == '__main__':
    interface()
    # python Run.py --folder C:\OmicData\temp initunpfasta --fastaFolder C:\OmicData\temp\fasta --unreviewed False --isoform False --split False
