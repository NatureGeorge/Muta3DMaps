# @Date:   2019-10-24T23:35:42+08:00
# @Email:  1730416009@stu.suda.edu.cn
# @Filename: RetrievePDB.py
# @Last modified time: 2019-11-24T22:59:10+08:00
import wget
import urllib
import ftplib
import os
from collections.abc import Iterable, Iterator
from multiprocessing.dummy import Pool
from time import sleep
from random import uniform
from ..Utils.Logger import RunningLogger
from ..Utils.FileIO import decompression, HandleIO
"""
Requirement:
1. Fast and Legal
2. Keep Update

Reference:
http://www.wwpdb.org/ftp/pdb-ftp-sites
http://zetcode.com/python/ftp/
"""
_FTP_SITE = ["RCSB", "PDBE", "PDBJ"]
_FTP_HEADER = "ftp://"
_RCSB_FTP = "ftp.rcsb.org"
_PDBE_FTP = "ftp.ebi.ac.uk"
_PDBJ_FTP = "ftp.pdbj.org"
# _FTP_HOST = {"RCSB": _RCSB_FTP, "PDBE": _PDBE_FTP, "PDBJ": _PDBJ_FTP}
_FTP_HOST = dict(zip(_FTP_SITE, [_RCSB_FTP, _PDBE_FTP, _PDBJ_FTP]))
_RCSB_DIVIDED = "pub/pdb/data/structures/divided"
_PDBE_DIVIDED = "pub/databases/pdb/data/structures/divided"
_PDBJ_DIVIDED = "pub/pdb/data/structures/divided"
# _DIVIDED_PATH = {"RCSB": _RCSB_DIVIDED, "PDBE": _PDBE_DIVIDED, "PDBJ": _PDBJ_DIVIDED}
_DIVIDED_PATH = dict(
    zip(_FTP_SITE, [_RCSB_DIVIDED, _PDBE_DIVIDED, _PDBJ_DIVIDED]))
_COMPLETE_TAGE = "226"  # "226-File successfully transferred", "226 Transfer complete"
_FORMAT_DICT = {
    "XML": ".xml",
    "mmCIF": ".cif",
    "pdb": ".pdb",
}
_RAW_FORMAT = {".pdb": ".ent"}
_RCSB_HTTP_VIEW = "https://files.rcsb.org/view/"
_RCSB_HTTP = "https://files.rcsb.org/download/"
_RCSB_REPORT_HTTP = "https://files.rcsb.org/pub/pdb/validation_reports/"
_XML_REPORT_KEY = ["percent-rota-outliers", "percent-rama-outliers", "clashscore"]


class RetrievePDB:
    """
    Retrieve PDB File

    PDB Archive Reference: wwPDB [1]_ [2]_ [3]_ [4]_, RCSB [5]_ [6]_, PDBe [7]_, PDBj [8]_
    FTP Site Reference: [9]_
    Code Reference: [10]_

    .. [1] http://www.wwpdb.org/
    .. [2] H.M. Berman, K. Henrick, H. Nakamura (2003) Announcing the worldwide Protein Data Bank Nature Structural Biology 10 (12): 980.
    .. [3] H.M. Berman, K. Henrick, H.Nakamura, J.L. Markley (2007) The Worldwide Protein Data Bank (wwPDB): Ensuring a single, uniform archive of PDB data Nucleic Acids Res. 35 (Database issue): D301-3.
    .. [4] wwPDB consortium. (2019) Protein Data Bank: the single global archive for 3D macromolecular structure data. Nucleic Acids Res 47: D520-D528 doi: 10.1093/nar/gky949.
    .. [5] http://www.rcsb.org/
    .. [6] H.M. Berman, J. Westbrook, Z. Feng, G. Gilliland, T.N. Bhat, H. Weissig, I.N. Shindyalov, P.E. Bourne. (2000) The Protein Data Bank Nucleic Acids Research, 28: 235-242.
    .. [7] https://www.ebi.ac.uk/pdbe/
    .. [8] https://pdbj.org/
    .. [9] http://www.wwpdb.org/ftp/pdb-ftp-sites
    .. [10] http://zetcode.com/python/ftp/

    Script:

    .. code-block:: python
        :linenos:

        ftpPDB = RetrievePDB("C:/Users/Nature/Downloads/PDBE/", ftpSite="PDBE", format="pdb")
        ftpPDB.ftp_retrieve(pdbs=['10z1', '10z2', '2xyn', '10z3'], remove=True)
        print(ftpPDB)
        printList(ftpPDB.getFail())
        ftpPDB.quick_ftp_retrieve('5js8', remove=False)
        ftpPDB.quick_http_retrieve('5js1', module="urllib", view=False, bioAssembly=1, remove=True)

    """

    def __init__(self, downloadPath, loggingPath, ftpSite="PDBE", format="mmCIF"):
        self.setDownloadPath(downloadPath)
        self.setFormat(format)
        self.setFTPSite(ftpSite)
        self.pdbs = []
        self.fail = []
        self.Logger = RunningLogger("RetrievePDB", loggingPath)
        self.Logger.logger.info(repr(self))

    def __len__(self):
        if self.pdbs is None:
            return 0
        else:
            return len(self.pdbs)

    def __getitem__(self, x):
        path = os.path.join(self.downloadPath, "%s%s" % (x, self.tail))
        if os.path.exists(path):
            return path
        else:
            return None

    def __repr__(self):
        format = "%s: %s, "
        string = "RetrievePDB: {%s}"
        content = ""
        for key, value in self.__dict__.items():
            if key not in ['pdbs', 'fail']:
                content += format % (key, value)
        content += format % ("len(pdbs)", len(self))
        content += format % ("len(fail)", len(self.fail))
        return string % content

    def getFail(self):
        """
        Get the list of PDB ids that failed to fetch files
        """
        return self.fail

    def setFTPSite(self, ftpSite):
        """
        Set the value of PDB FTP site

        :param str ftpSite: The FTP site of retriving PDB files, default value: `RCSB`, {RCSB, PDBE, PDBJ}

        """
        if ftpSite not in _FTP_SITE:
            raise ValueError(
                "Illegal site name. Please select from %s" % _FTP_SITE)
        else:
            self.ftpSite = ftpSite
            self.host = _FTP_HOST[ftpSite]
            self.dividedPath = _DIVIDED_PATH[ftpSite]

    def setFormat(self, format):
        """
        Set the format of PDB file

        :param str format: The file format of PDB file, default value: `mmCIF`, {mmCIF, pdb, XML}

        """
        if format not in _FORMAT_DICT.keys():
            raise ValueError(
                "Illegal format name. Please select from mmCIF, pdb, XML")
        else:
            self.format = format
            self.tail = _FORMAT_DICT[format]
            self.raw_tail = _RAW_FORMAT.get(self.tail, self.tail)
            if format == "pdb":
                self.prefix = format
            else:
                self.prefix = ""

    def setDownloadPath(self, path):
        """Set the value of PDB file download path"""
        if os.path.isdir(path):
            self.downloadPath = path
        else:
            raise ValueError("Not a valid download path")

    def setPDBs(self, pdbs):
        """
        Set the value of PDB id(s)

        :param pdbs: single PDB id or PDB ids
        :type pdbs: Iterable or Iterator or str

        * ``str``: single PDB id
        * ``Iterable, Iterator``: PDB ids

        """
        if isinstance(pdbs, str):
            self.pdbs = [pdbs]
        elif isinstance(pdbs, Iterable):
            self.pdbs = sorted(pdbs, key=lambda x: x[1:3] + x[0] + x[3])
        elif isinstance(pdbs, Iterator):
            self.pdbs = pdbs
        elif pdbs is None:
            raise ValueError("pdbs should not be None. Please Specify pdbs!")
        else:
            raise ValueError("Invalid Input")

    def ftp_retrieve(self, **kwargs):
        """
        Retrieve PDB files via FTP Connection

        :param pdbs: single PDB id or PDB ids
        :param bool remove: whether remove the compressed file, default value: `True`
        :type pdbs: Iterable or Iterator or str

        """
        # Check PDB List
        pdbs = kwargs.get('pdbs', False)
        remove = kwargs.get('remove', True)
        if pdbs is not False:
            self.setPDBs(pdbs)
        elif self.pdbs is None:
            raise ValueError("Please Specify pdbs!")
        # Connect
        with ftplib.FTP(self.host) as ftp:
            self.Logger.logger.warning(ftp.getwelcome())
            ftp.login()  # anonymous account
            ftp.cwd("%s/%s" % (self.dividedPath, self.format))
            # Start to retrieve
            cur = ""
            for pdb in self.pdbs:
                pdb = pdb.lower()
                subPath = pdb[1:3]
                try:
                    if cur != subPath:
                        if cur == "":
                            ftp.cwd(subPath)
                        else:
                            ftp.cwd("../%s" % subPath)
                        cur = subPath
                    file_orig = "%s%s%s.gz" % (self.prefix, pdb, self.raw_tail)
                    # Start to download
                    filename = os.path.join(
                        self.downloadPath, "%s%s.gz" % (pdb.upper(), self.tail))
                    self.Logger.logger.info("Downloading File: %s" % filename)
                    data = HandleIO(open(filename, 'w+b'))
                    res = ftp.retrbinary('RETR ' + file_orig, data.append)
                    data.close()
                    # Check Data Completeness
                    if not res.startswith(_COMPLETE_TAGE):
                        self.Logger.logger.warning('Download failed: %s' % res)
                        if os.path.isfile(filename):
                            os.remove(filename)
                        self.fail.append(pdb)
                        continue
                    decompression(filename, remove=remove, logger=self.Logger.logger)
                except ftplib.error_perm as e:
                    self.Logger.logger.warning('FTP error: %s' % e)
                    if 'filename' in locals().keys():
                        data.close()
                        os.remove(filename)
                    self.fail.append(pdb)
                    continue
        self.Logger.logger.info("FTP Closed")

    def quick_ftp_retrieve(self, pdb, remove=True):
        """
        Download PDB file via FTP with ``wget``

        *Download only one file at a time*

        :param str pdb: single PDB id
        :param bool remove: whether remove the compressed file, default value: `True`
        """
        pdb = pdb.lower()
        file_orig = "%s%s%s.gz" % (self.prefix, pdb, self.raw_tail)
        site = "%s%s/%s/%s/%s/%s" % (_FTP_HEADER, self.host,
                                     self.dividedPath, self.format, pdb[1:3], file_orig)
        path = os.path.join(self.downloadPath, "%s%s.gz" % (pdb.upper(), self.tail))
        self.Logger.logger.info("Downloading File: %s" % path)
        try:
            wget.download(site, out=path)
            self.Logger.logger.info("\n")
        except urllib.error.URLError:
            self.Logger.logger.warning("Download failed")
            self.fail.append(pdb)
        decompression(path, remove=remove, logger=self.Logger.logger)

    def quick_http_retrieve(self, pdb, module="wget", view=False, bioAssembly="", extension=".gz", remove=True):
        """
        Download PDB file via HTTP with ``wget.download`` or ``urllib.request.urlopen``

        *RCSB Only*

        :param str pdb: single PDB id
        :param str module: For user to select the module, default value: `wget`, {"wget", "urllib"}
        :param bool view: Whether get the PDB data from the web pages
        :param str bioAssembly: The bioAssembly id, default value: `''`
        :param str extension: Compressed file tail, default value: `.gz`
        :param bool remove: whether remove the compressed file, default value: `True`

        """
        # Whether download from web page
        if view:
            fileName = "%s%s" % (pdb, self.tail)
            url = "%s%s" % (_RCSB_HTTP_VIEW, fileName)
        else:
            fileName = "{pdb}{tail}{bioAssembly}{extension}".format(
                pdb=pdb, tail=self.tail, bioAssembly=bioAssembly, extension=extension)
            url = "{site}{fileName}".format(site=_RCSB_HTTP, fileName=fileName)
        # r = requests.get(url)
        path = os.path.join(self.downloadPath, fileName)
        self.Logger.logger.info("Downloading File: %s" % path)
        try:
            # Select module method: urllib
            if module != "wget":
                with open(path, 'w+b') as fw:
                    # fw.write(r.content)
                    text = urllib.request.urlopen(url).read()
                    fw.write(text)  # .decode('utf-8')
            # Select module method: wget
            else:
                wget.download(url, out=path)
                self.Logger.logger.info("\n")
        except urllib.error.URLError:
            self.Logger.logger.warning("Download failed")
            self.fail.append(pdb)
            return
        # Whether to decompress
        if not view and extension == '.gz':
            decompression(path, remove=remove, logger=self.Logger.logger)

    def quick_http_retrieve_report(self, pdb, remove=True):
        # https://files.rcsb.org/pub/pdb/validation_reports/g9/3g96/3g96_validation.xml.gz
        pdb = pdb.lower()
        fileName = "{pdb}_validation.xml.gz".format(pdb=pdb)
        url = "{reportSite}{sub}/{pdb}/{fileName}".format(reportSite=_RCSB_REPORT_HTTP, sub=pdb[1:3], pdb=pdb, fileName=fileName)
        path = os.path.join(self.downloadPath, fileName)
        self.Logger.logger.info("Downloading File: %s" % path)
        try:
            wget.download(url, out=path)
            self.Logger.logger.info("\n")
        except urllib.error.URLError:
            self.Logger.logger.warning("Download failed")
            self.fail.append(pdb)
            return
        decompression(path, remove=remove, logger=self.Logger.logger)


class MPWrapper:
    """
    Multiprocessing wrapper for ``RetrievePDB``

    When there is a large number of PDB files to download, this class is helpful.
    But Need to be careful with the numbers of processes and the time of sleep.

    Script:

    .. code-block:: python
        :linenos:

        pdbs = ['1A02', '3KBZ', '3KC0', '3KC1', '3KMU', '3KMW', '3KYC', '3KYD', ...]
        mpw = MPWrapper("C:/Users/Nature/Downloads/")
        # fail = mpw.http_retrieve(pdbs)
        # fail = mpw.http_retrieve(pdbs, module="urllib")
        # fail = mpw.ftp_retrieve_wget(pdbs)
        fail = mpw.ftp_retrieve_batch(pdbs)
        printList(fail)

    :param str downloadPath: File folder of Downloaded PDB files
    :param int processes: Number of processes, default value: `3`
    :param int maxSleep: Max sleep time, default value: `3`
    :param str ftpSite: The FTP site of retriving PDB files, default value: `RCSB`, {RCSB, PDBE, PDBJ}
    :param str format: The file format of PDB file, default value: `mmCIF`, {mmCIF, pdb, XML}

    """

    def __init__(self, downloadPath, loggingPath, processes=3, maxSleep=3, ftpSite="RCSB", format="mmCIF"):
        self.setProcesses(processes, maxSleep)
        self.retrievePDB = RetrievePDB(
            downloadPath, loggingPath, ftpSite=ftpSite, format=format)

    def setProcesses(self, processes, maxSleep):
        """
        Set the value of ``processes`` and ``maxSleep``

        :param int processes: Number of processes
        :param int maxSleep: Max sleep time
        """
        if processes > 20:
            self.retrievePDB.Logger.logger.warning("MPWrapper: Too many processes. Be careful !")
        self.processes = processes
        self.maxSleep = maxSleep

    def http_retrieve(self, pdbs, module="wget", view=False, bioAssembly="", extension=".gz", remove=True):
        """
        Retrieve PDB file via http with ``wget.download`` or ``urllib.request.urlopen``

        **HTTP SITE: RCSB ONLY**

        :param pdbs: An object containing the PDB ids that need to be download
        :param str module: For user to select the module, default value: `wget`, {"wget", "urllib"}
        :param bool view: Whether get the PDB data from the web pages, default value: `False`
        :param str bioAssembly: The bioAssembly id, default value: `''`
        :param str extension: Compressed file tail, default value: `.gz`
        :param bool remove: whether remove the compressed file, default value: `True`
        :type pdbs: Iterable or Iterator
        :return: A fail list that contains the PDB ids falied to download
        :rtype: list(str)
        """
        def register(pdb):
            stop = uniform(0, self.maxSleep)
            sleep(stop)
            self.retrievePDB.quick_http_retrieve(pdb, module=module, view=view, bioAssembly=bioAssembly, extension=extension, remove=remove)
            # print(pdb, stop)

        pool = Pool(processes=self.processes)
        pool.map(register, pdbs)
        return self.retrievePDB.getFail()

    def http_retrieve_report(self, pdbs, remove=True):
        def register(pdb):
            stop = uniform(0, self.maxSleep)
            sleep(stop)
            self.retrievePDB.quick_http_retrieve_report(pdb, remove=remove)
            # print(pdb, stop)

        pool = Pool(processes=self.processes)
        pool.map(register, pdbs)
        return self.retrievePDB.getFail()

    def ftp_retrieve_wget(self, pdbs, remove=True):
        """
        Download PDB file via FTP with ``wget``

        :param pdbs: An object containing the PDB ids that need to be download
        :param bool remove: whether remove the compressed file, default value: `True`
        :type pdbs: Iterable or Iterator
        :return: A fail list that contains the PDB ids falied to download
        :rtype: list(str)
        """
        def register(pdb):
            stop = uniform(0, self.maxSleep)
            sleep(stop)
            self.retrievePDB.quick_ftp_retrieve(pdb, remove=remove)
            # print(pdb, stop)

        pool = Pool(processes=self.processes)
        pool.map(register, pdbs)
        return self.retrievePDB.getFail()

    def ftp_retrieve_batch(self, pdbs, remove=True, chunksize=100):
        """
        Retrieve PDB files via FTP Connection

        :param pdbs: An object containing the PDB ids that need to be download
        :param bool remove: whether remove the compressed file, default value: `True`
        :param int chunksize: the size of PDBs that query during a single FTP connection, default value: `100`
        :type pdbs: Iterable or Iterator
        :return: A fail list that contains the PDB ids falied to download
        :rtype: list(str)

        """
        assert isinstance(pdbs, Iterable), "pdbs should be an Iterable object in this function!"

        def register(chunk):
            sleep(uniform(0, self.maxSleep))
            self.retrievePDB.ftp_retrieve(pdbs=chunk, remove=remove)
            # print(chunk)

        chunks = [pdbs[i:i + chunksize]
                  for i in range(0, len(pdbs), chunksize)]
        pool = Pool(processes=self.processes)
        pool.map(register, chunks)
        return self.retrievePDB.getFail()


# if __name__ == "__main__":
    # ftpPDB = RetrievePDB("C:/Users/Nature/Downloads/PDBE/", ftpSite="PDBE", format="pdb")
    # ftpPDB.ftp_retrieve(pdbs=['10z1', '10z2', '10z3'], remove=True)
    # ftpPDB.quick_ftp_retrieve('5js8', remove=False)
    # ftpPDB.quick_http_retrieve('5js1', module="urllib", view=False, bioAssembly=1, remove=True)
    # print(ftpPDB)
    # printList(ftpPDB.getFail())
    # pdbs = ['1A02', ...]
    # mpw = MPWrapper("C:/Users/Nature/Downloads/")  # 1.1: 79s 1.2: 90s 2: 114s 3:154s
    # mpw.http_retrieve_report(pdbs[:10])  # 79s
    # mpw.http_retrieve(pdbs, module="urllib")  # 90s
    # mpw.ftp_retrieve_wget(pdbs)  # 114s
    # mpw.ftp_retrieve_batch(pdbs)  # 154s
