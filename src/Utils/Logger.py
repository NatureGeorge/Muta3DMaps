# @Date:   2019-11-20T18:01:14+08:00
# @Email:  1730416009@stu.suda.edu.cn
# @Filename: Logger.py
# @Last modified time: 2019-11-20T18:17:40+08:00
import logging


class RunningLogger:
    """
    Initialize the logger

    :param str loggerName: Name for RunningLogger Object
    :param str format: Format for loggering
    :param str loggingPath: FilePath for logger file

    """

    def __init__(self, loggerName="RunningLogger", loggingPath=None, format="%(asctime)s %(name)s %(levelname)s %(message)s"):
        formatter = logging.Formatter(format)
        self.logger = logging.getLogger(loggerName)
        self.logger.setLevel(logging.DEBUG)

        streamHandler = logging.StreamHandler()
        streamHandler.setLevel(logging.WARNING)
        streamHandler.setFormatter(formatter)
        self.logger.addHandler(streamHandler)

        try:
            self.loggingPath = loggingPath
            fileHandler = logging.FileHandler(filename=self.loggingPath)
            fileHandler.setLevel(logging.DEBUG)
            fileHandler.setFormatter(formatter)
            self.logger.addHandler(fileHandler)
        except Exception:
            self.logger.warning("Invalid file path for logging file ! Please specifiy loggingPath=...")
