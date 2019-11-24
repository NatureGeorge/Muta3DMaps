# @Date:   2019-10-25T13:16:43+08:00
# @Email:  1730416009@stu.suda.edu.cn
# @Filename: FileReport.py
# @Last modified time: 2019-10-25T13:31:09+08:00


class FileReport(dict):
    """
    class Info:
        def __init__(self, key):
            self.key
            self.content = []

        def append(ele):
            self.content.append(ele)
    """

    def __init__(self, path, mode='w+'):
        self.path = self.path
        self.mode = mode
        # self.info = [Info("message"), Info("warning"), Info("bugs")]
        self.mes = []
        self.warning = []
        self.bugs = []

    def dead(self):
        self.handle.close()

    def info(self, mes):
        self.mes.append(mes)

    def warn(self, warning):
        self.warning.append(warning)

    def bug(self, bug):
        self.bugs.append(bug)

    def output(self):
        # ob = [self.mes, self.warning, self.bugs]
        # for info in ob:
        return
