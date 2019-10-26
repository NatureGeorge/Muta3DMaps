# @Date:   2019-10-24T23:35:42+08:00
# @Email:  1730416009@stu.suda.edu.cn
# @Filename: RetrivePDB.py
# @Last modified time: 2019-10-26T12:22:18+08:00
import wget
import gzip
import urllib
import ftplib
import shutil
import os
from collections import Iterable, Iterator

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
_DIVIDED_PATH = dict(zip(_FTP_SITE, [_RCSB_DIVIDED, _PDBE_DIVIDED, _PDBJ_DIVIDED]))
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

    class HandleIO:
        def __init__(self, handle):
            self.handle = handle

        def append(self, block):
            self.handle.write(block)

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
        string = "FTP_PDB: {%s}"
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
        """Set the value of PDB FTP site
        * RCSB
        * PDBE
        * PDBJ
        """
        if ftpSite not in _FTP_SITE:
            raise ValueError("Illegal site name. Please select from %s" % _FTP_SITE)
        else:
            self.ftpSite = ftpSite
            self.host = _FTP_HOST[ftpSite]
            self.dividedPath = _DIVIDED_PATH[ftpSite]

    def setFormat(self, format):
        """Set the format of PDB file
        * mmCIF
        * pdb
        * XML
        """
        if format not in _FORMAT_DICT.keys():
            raise ValueError("Illegal format name. Please select from mmCIF, pdb, XML")
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
        * str: single PDB id
        * Iterable, Iterator: PDB ids
        """
        if isinstance(pdbs, str):
            self.pdbs = [pdbs]
        elif isinstance(pdbs, Iterable):
            self.pdbs = sorted(pdbs, key=lambda x: x[1:3]+x[0]+x[3])
        elif isinstance(pdbs, Iterator):
            self.pdbs = pdbs
        elif pdbs is None:
            raise ValueError("pdbs should not be None. Please Specify pdbs!")
        else:
            raise ValueError("Invalid Input")

    def ftp_retrive(self, **kwargs):
        """
        Retrive PDB files
        * pdbs: PDB ids
        * remove: whether remove the compressed file
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
            try:
                ftp.login()  # anonymous account
                ftp.cwd("%s/%s" % (self.dividedPath, self.format))
                # print(ftp.pwd())
                # Start to retrive
                cur = ""
                for pdb in self.pdbs:
                    pdb = pdb.lower()
                    subPath = pdb[1:3]
                    try:
                        """
                        # Check PDB id
                        if subPath not in ftp.nlst():
                            self.fail.append(pdb)
                            continue
                        """
                        if cur != subPath:
                            if cur == "":
                                ftp.cwd(subPath)
                                print(ftp.pwd())
                            else:
                                ftp.cwd("../%s" % subPath)
                                print(ftp.pwd())
                            cur = subPath
                        file_orig = "%s%s%s.gz" % (self.prefix, pdb, self.raw_tail)
                        """
                        # Check PDB id
                        if file_orig not in ftp.nlst():
                            self.fail.append(pdb)
                            continue
                        """
                        # Start to download
                        filename = os.path.join(self.downloadPath, "%s%s.gz" % (pdb, self.tail))
                        print("Downloading File: %s" % filename)
                        data = self.HandleIO(open(filename, 'w+b'))
                        res = ftp.retrbinary('RETR ' + file_orig, data.append)
                        data.handle.close()
                        # Check Data Completeness
                        if not res.startswith(_COMPLETE_TAGE):
                            print('Download failed', res)
                            if os.path.isfile(filename):
                                os.remove(filename)
                            # ftp.cwd("../")
                            self.fail.append(pdb)
                            continue
                        self.decompression(filename, remove=remove)
                        # ftp.cwd("../")
                    except ftplib.error_perm as e:
                        print('FTP error:', e)
                        if 'filename' in locals().keys():
                            data.handle.close()
                            os.remove(filename)
                        self.fail.append(pdb)
                        continue
            except ftplib.all_errors as e:
                print('FTP error:', e)
        print("Close FTP")

    def decompression(self, path, extension=".gz", remove=True, outputPath=None):
        """
        Decompression gz file
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
        Download PDB file via ftp with ```wget```
        """
        pdb = pdb.lower()
        file_orig = "%s%s%s.gz" % (self.prefix, pdb, self.raw_tail)
        site = "%s%s/%s/%s/%s/%s" % (_FTP_HEADER, self.host, self.dividedPath, self.format, pdb[1:3], file_orig)
        path = os.path.join(self.downloadPath, "%s%s.gz" % (pdb, self.tail))
        print("Downloading File: %s" % path)
        try:
            wget.download(site, out=path)
        except urllib.error.URLError:
            print("Download failed")
        self.decompression(path, remove=remove)

    def quick_http_retrive(self, pdb, module="wget", view=False, bioAssembly="", extension=".gz", remove=True):
        """
        Download PDB file via http with ```wget.download``` or ```urllib.request.urlopen```
        * RCSB Only
        """
        # Whether download from web page
        if view:
            fileName = "%s%s" % (pdb, self.tail)
            url = "%s%s" % (_RCSB_HTTP_VIEW, fileName)
        else:
            fileName = "{pdb}{tail}{bioAssembly}{extension}".format(pdb=pdb, tail=self.tail, bioAssembly=bioAssembly, extension=extension)
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
                return
        # Whether to decompress
        if not view and extension == '.gz':
            self.decompression(path, remove=remove)


if __name__ == "__main__":
    ftpPDB = RetrivePDB("C:/Users/Nature/Downloads/PDBJ/", ftpSite="PDBJ", format="pdb")
    # ftpPDB.ftp_retrive(pdbs=['5jp8', '5jp5'], remove=True)
    # ftpPDB.quick_ftp_retrive('5js8', remove=False)
    # ftpPDB.quick_http_retrive('5js1', module="urllib", view=False, bioAssembly=1, remove=True)
    # print(ftpPDB)
