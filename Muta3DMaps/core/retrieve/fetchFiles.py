# @Created Date: 2020-01-18 10:53:07 am
# @Filename: fetchFiles.py
# @Email:  1730416009@stu.suda.edu.cn
# @Author: ZeFeng Zhu
# @Last Modified: 2020-02-08 04:15:07 pm
# @Copyright (c) 2020 MinghuiGroup, Soochow University
import os
from time import perf_counter
import asyncio
import aiohttp
import aiofiles
from unsync import unsync
from tenacity import retry, wait_random_exponential, stop_after_attempt, after_log, RetryError
import logging
from tqdm import tqdm
from typing import Iterable, Iterator, Union


class UnsyncFetch(object):

    logger = logging.getLogger("UnsyncFetch"); logger.setLevel(logging.DEBUG)
    streamHandler = logging.StreamHandler(); streamHandler.setLevel(logging.INFO)
    formatter = logging.Formatter("%(asctime)s %(name)s %(levelname)s %(message)s"); streamHandler.setFormatter(formatter)
    logger.addHandler(streamHandler)

    @classmethod
    def set_logging_fileHandler(cls, path: str, level: int = logging.DEBUG, formatter=formatter):
        try:
            fileHandler = logging.FileHandler(filename=path)
            fileHandler.setLevel(level)
            fileHandler.setFormatter(formatter)
            cls.logger.addHandler(fileHandler)
            cls.logger.info(f"Logging file has been created in {path}")
        except Exception:
            cls.logger.warning("Invalid file path for logging file ! Please specifiy path=...")

    @classmethod
    @retry(wait=wait_random_exponential(multiplier=1, max=15), stop=stop_after_attempt(3), after=after_log(logger, logging.WARNING))
    async def download_file(cls, url: str):
        cls.logger.debug(f"Start to get file: {url}")
        async with aiohttp.ClientSession() as session:
            async with session.get(url) as resp:
                if resp.status == 200:
                    return await resp.read()
                elif resp.status == 404:
                    return None
                else:
                    mes = "code={resp.status}, message={resp.reason}, headers={resp.headers}"
                    raise Exception(mes.format(resp=resp))

    @classmethod
    @unsync
    async def save_file(cls, path: str, data: bytes):
        cls.logger.debug(f"Start to save file: {path}")
        async with aiofiles.open(path, 'wb') as jsonFile:
            await jsonFile.write(data)

    @classmethod
    @unsync
    async def fetch_file(cls, semaphore: asyncio.Semaphore, url: str, path: str, rate: float):
        async with semaphore:
            try:
                data = await cls.download_file(url)
                if data is not None:
                    await cls.save_file(path, data)
                    await asyncio.sleep(rate)
                    return path
                else:
                    cls.logger.warning(f"404 for: {url}")
            except RetryError:
                cls.logger.error(f"Retry failed for: {url}")

    @classmethod
    @unsync
    async def multi_tasks(cls, workdir: str, tasks, concur_req: int = 4, rate: float = 1.5):
        '''
        Template for multiTasking

        TODO
            1. asyncio.Semaphore
            2. unit func
        '''
        semaphore = asyncio.Semaphore(concur_req)
        # await asyncio.gather(*[cls.fetch_file(semaphore, url, os.path.join(workdir, path), rate) for url, path in tasks])
        tasks = [cls.fetch_file(semaphore, url, os.path.join(workdir, path), rate) for url, path in tasks]
        return [await fob for fob in tqdm(asyncio.as_completed(tasks), total=len(tasks))]

    @classmethod
    def main(cls, workdir: str, data: Union[Iterable, Iterator], logName='UnsyncFetch'):
        cls.set_logging_fileHandler(os.path.join(workdir, f'{logName}.log'))
        t0 = perf_counter()
        res = cls.multi_tasks(workdir, data, 5).result()
        elapsed = perf_counter() - t0
        cls.logger.info(f'downloaded in {elapsed}s')
        return res

