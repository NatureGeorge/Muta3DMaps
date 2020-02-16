# @Created Date: 2020-02-16 09:58:31 pm
# @Filename: log.py
# @Email:  1730416009@stu.suda.edu.cn
# @Author: ZeFeng Zhu
# @Last Modified: 2020-02-16 09:58:41 pm
# @Copyright (c) 2020 MinghuiGroup, Soochow University
import logging
from typing import Optional

class Abclog(object):
    '''
    Abstract class of logging template
    '''
    @classmethod
    def init_logger(cls, logName: Optional[str] = None, logger: Optional[logging.Logger] = None, log_level: int = logging.DEBUG, stream_log_level: int = logging.INFO):
        if hasattr(cls, 'logger') and cls.logger is not None:
            pass
        elif logger is None:
            cls.logger = logging.getLogger(logName)
            cls.logger.setLevel(log_level)
            cls.streamHandler = logging.StreamHandler()
            cls.streamHandler.setLevel(stream_log_level)
            cls.formatter = logging.Formatter(
                "%(asctime)s %(name)s %(levelname)s %(message)s")
            cls.streamHandler.setFormatter(cls.formatter)
            cls.logger.addHandler(cls.streamHandler)
        else:
            cls.logger = logger

    @classmethod
    def set_logging_fileHandler(cls, path: str, level: int = logging.DEBUG, logName: Optional[str] = None):
        if logName is not None and not hasattr(cls, 'logger'):
            cls.init_logger(logName)
        try:
            fileHandler = logging.FileHandler(filename=path)
            fileHandler.setLevel(level)
            fileHandler.setFormatter(cls.formatter)
            cls.logger.addHandler(fileHandler)
            cls.logger.info(f"Logging file in {path}")
        except Exception:
            cls.logger.warning("Invalid file path for logging file ! Please specifiy path=...")
