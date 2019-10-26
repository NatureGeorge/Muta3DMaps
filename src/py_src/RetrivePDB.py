# @Date:   2019-10-24T23:35:42+08:00
# @Email:  1730416009@stu.suda.edu.cn
# @Filename: RetrivePDB.py
# @Last modified time: 2019-10-27T01:20:38+08:00
import wget
import gzip
import urllib
import ftplib
import shutil
import os
from collections import Iterable, Iterator
from multiprocessing.dummy import Pool
from time import sleep
from random import uniform
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


def printList(list):
    string = ""
    for i in list:
        string += "{0}\n".format(i)
    print(string)


class RetrivePDB:
    """
    Retrive PDB File

    Script:

    .. code-block:: python
        :linenos:

        ftpPDB = RetrivePDB("C:/Users/Nature/Downloads/PDBE/", ftpSite="PDBE", format="pdb")
        ftpPDB.ftp_retrive(pdbs=['10z1', '10z2', '2xyn', '10z3'], remove=True)
        print(ftpPDB)
        printList(ftpPDB.getFail())
        ftpPDB.quick_ftp_retrive('5js8', remove=False)
        ftpPDB.quick_http_retrive('5js1', module="urllib", view=False, bioAssembly=1, remove=True)

    """
    class HandleIO:
        def __init__(self, handle):
            self.handle = handle

        def append(self, block):
            self.handle.write(block)

        def close(self):
            self.handle.close()

    def __init__(self, downloadPath, pdbs=None, ftpSite="PDBE", format="mmCIF"):
        self.setDownloadPath(downloadPath)
        self.setFormat(format)
        self.setFTPSite(ftpSite)
        self.pdbs = pdbs
        self.fail = []
        print(self)

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
        string = "RetrivePDB: {%s}"
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

    def ftp_retrive(self, **kwargs):
        """
        Retrive PDB files via FTP Connection

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
            print(ftp.getwelcome())
            ftp.login()  # anonymous account
            ftp.cwd("%s/%s" % (self.dividedPath, self.format))
            # Start to retrive
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
                        self.downloadPath, "%s%s.gz" % (pdb, self.tail))
                    print("Downloading File: %s" % filename)
                    data = self.HandleIO(open(filename, 'w+b'))
                    res = ftp.retrbinary('RETR ' + file_orig, data.append)
                    data.close()
                    # Check Data Completeness
                    if not res.startswith(_COMPLETE_TAGE):
                        print('Download failed', res)
                        if os.path.isfile(filename):
                            os.remove(filename)
                        self.fail.append(pdb)
                        continue
                    self.decompression(filename, remove=remove)
                except ftplib.error_perm as e:
                    print('FTP error:', e)
                    if 'filename' in locals().keys():
                        data.close()
                        os.remove(filename)
                    self.fail.append(pdb)
                    continue
        print("FTP Closed")

    def decompression(self, path, extension=".gz", remove=True, outputPath=None):
        """
        Decompress gz file

        :param str path: The file path of the compressed file
        :param str extension: Compressed file tail, default value: `.gz`
        :param bool remove: Whether remove the compressed file, default value: `True`
        :param str outputPath: File safe path, default value: `None`
        """

        """
        with gzip.GzipFile(mode="rb", fileobj=open(path, 'rb')) as raw:
            with open(path[:-len(extension)], "wb") as file:
                file.write(raw.read())
        """
        if outputPath is None:
            outputPath = path[:-len(extension)]

        with gzip.open(path, 'rb') as raw:
            with open(outputPath, 'wb') as file:
                shutil.copyfileobj(raw, file)
        try:
            if remove:
                os.remove(path)
        except Exception as e:
            print(e)

    def quick_ftp_retrive(self, pdb, remove=True):
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
        path = os.path.join(self.downloadPath, "%s%s.gz" % (pdb, self.tail))
        print("Downloading File: %s" % path)
        try:
            wget.download(site, out=path)
        except urllib.error.URLError:
            print("Download failed")
            self.fail.append(pdb)
        self.decompression(path, remove=remove)

    def quick_http_retrive(self, pdb, module="wget", view=False, bioAssembly="", extension=".gz", remove=True):
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
        print("Downloading File: %s" % path)
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
        except urllib.error.URLError:
            print("Download failed")
            self.fail.append(pdb)
            return
        # Whether to decompress
        if not view and extension == '.gz':
            self.decompression(path, remove=remove)


class MPWrapper:
    """
    Multiprocessing wrapper for ``RetrivePDB``

    When there is a large number of PDB files to download, this class is helpful.
    But Need to be careful with the numbers of processes and the time of sleep.

    Script:

    .. code-block:: python
        :linenos:

        pdbs = ['1A02', '3KBZ', '3KC0', '3KC1', '3KMU', '3KMW', '3KYC', '3KYD', ...]
        mpw = MPWrapper("C:/Users/Nature/Downloads/")
        # fail = mpw.http_retrive(pdbs)
        # fail = mpw.http_retrive(pdbs, module="urllib")
        # fail = mpw.ftp_retrive_wget(pdbs)
        fail = mpw.ftp_retrive_batch(pdbs)
        printList(fail)

    :param str downloadPath: File folder of Downloaded PDB files
    :param int processes: Number of processes, default value: `3`
    :param int maxSleep: Max sleep time, default value: `3`
    :param str ftpSite: The FTP site of retriving PDB files, default value: `RCSB`, {RCSB, PDBE, PDBJ}
    :param str format: The file format of PDB file, default value: `mmCIF`, {mmCIF, pdb, XML}

    """

    def __init__(self, downloadPath, processes=3, maxSleep=3, ftpSite="RCSB", format="mmCIF"):
        self.setProcesses(processes, maxSleep)
        self.retrivePDB = RetrivePDB(
            downloadPath, ftpSite=ftpSite, format=format)

    def setProcesses(self, processes, maxSleep):
        """
        Set the value of ``processes`` and ``maxSleep``

        :param int processes: Number of processes
        :param int maxSleep: Max sleep time
        """
        if processes > 20:
            print("MPWrapper: Too many processes. Be careful !")
        self.processes = processes
        self.maxSleep = maxSleep

    def http_retrive(self, pdbs, module="wget", view=False, bioAssembly="", extension=".gz", remove=True):
        """
        Retrive PDB file via http with ``wget.download`` or ``urllib.request.urlopen``

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
            self.retrivePDB.quick_http_retrive(pdb, module=module, view=view, bioAssembly=bioAssembly, extension=extension, remove=remove)
            # print(pdb, stop)

        pool = Pool(processes=self.processes)
        pool.map(register, pdbs)
        return self.retrivePDB.getFail()

    def ftp_retrive_wget(self, pdbs, remove=True):
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
            self.retrivePDB.quick_ftp_retrive(pdb, remove=remove)
            # print(pdb, stop)

        pool = Pool(processes=self.processes)
        pool.map(register, pdbs)
        return self.retrivePDB.getFail()

    def ftp_retrive_batch(self, pdbs, remove=True, chunksize=100):
        """
        Retrive PDB files via FTP Connection

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
            self.retrivePDB.ftp_retrive(pdbs=chunk, remove=remove)
            # print(chunk)

        chunks = [pdbs[i:i + chunksize]
                  for i in range(0, len(pdbs), chunksize)]
        pool = Pool(processes=self.processes)
        pool.map(register, chunks)
        return self.retrivePDB.getFail()


if __name__ == "__main__":
    # ftpPDB = RetrivePDB("C:/Users/Nature/Downloads/PDBE/", ftpSite="PDBE", format="pdb")
    # ftpPDB.ftp_retrive(pdbs=['10z1', '10z2', '10z3'], remove=True)
    # ftpPDB.quick_ftp_retrive('5js8', remove=False)
    # ftpPDB.quick_http_retrive('5js1', module="urllib", view=False, bioAssembly=1, remove=True)
    # print(ftpPDB)
    # printList(ftpPDB.getFail())
    pdbs = ['1A02',
            '3KBZ',
            '3KC0',
            '3KC1',
            '3KMU',
            '3KMW',
            '3KYC',
            '3KYD',
            '3L3C',
            '3L5P',
            '3L5R',
            '3L5S',
            '3L5T',
            '3L5U',
            '3L5V',
            '3L7U',
            '3LHR',
            '3M0D',
            '3M1D',
            '3MK4',
            '3MUD',
            '3MUP',
            '3NR2',
            '3OD5',
            '3OQ9',
            '3OZ1',
            '3P0E',
            '3P45',
            '3P4U',
            '3P87',
            '3PHB',
            '3Q0K',
            '3Q4F',
            '3Q84',
            '3Q91',
            '3QNI',
            '3R1H',
            '3R1L',
            '3R5J',
            '3R6G',
            '3R6L',
            '3R7B',
            '3R7N',
            '3R7S',
            '3REP',
            '3RJM']
    mpw = MPWrapper("C:/Users/Nature/Downloads/")  # 1.1: 79s 1.2: 90s 2: 114s 3:154s
    # mpw.http_retrive(pdbs)  # 79s
    # mpw.http_retrive(pdbs, module="urllib")  # 90s
    # mpw.ftp_retrive_wget(pdbs)  # 114s
    # mpw.ftp_retrive_batch(pdbs)  # 154s
