# @Created Date: 2020-01-18 10:53:07 am
# @Filename: fetchFiles.py
# @Email:  1730416009@stu.suda.edu.cn
# @Author: ZeFeng Zhu
# @Last Modified: 2020-02-09 08:50:23 pm
# @Copyright (c) 2020 MinghuiGroup, Soochow University
import os
from time import perf_counter
import asyncio
import aiohttp
import aiofiles
from unsync import unsync, Unfuture
from tenacity import retry, wait_random, stop_after_attempt, after_log, RetryError
import logging
from tqdm import tqdm
from typing import Iterable, Iterator, Union, Any, Optional, List, Dict, Coroutine, Callable
import ujson as json
from Muta3DMaps.core.log import Abclog
from Muta3DMaps.core.pdbe.decode import PDBeJsonDecoder, traversePDBeData, convertJson2other
import re


class UnsyncFetch(Abclog):
    '''
    Fetch files through API in unsync/async way

    Note:

    * Implement both GET and POST method
    * Load all the available data(stream) in memory as soon as it is received
    * Since the methods in this class would not load the entire file in memory,
      the response data could be a large file
    * Parameters of `tenacity.retry` is built-in
    
    Reference (following packages provide me with a lot of insights and inspiration)

        `aiohttp`

        * https://github.com/aio-libs/aiohttp
        * https://docs.aiohttp.org/en/stable/client_quickstart.html
        * https://docs.aiohttp.org/en/stable/streams.html

        `unsync`

        * https://github.com/alex-sherman/unsync

        `tenacity`

        * https://github.com/jd/tenacity
        * https://tenacity.readthedocs.io/en/latest/
    '''

    retry_kwargs = {
        'wait': wait_random(max=10),
        'stop': stop_after_attempt(3),
        # 'after': after_log(logger, logging.WARNING)
        }

    @classmethod
    async def download_file(cls, method: str, info: Dict, path: str):
        cls.logger.debug(f"Start to download file: {info}")
        async with aiohttp.ClientSession() as session:
            async_func = getattr(session, method)
            async with async_func(**info) as resp:
                if resp.status == 200:
                    async with aiofiles.open(path, 'wb') as fileOb:
                        # Asynchronous iterator implementation of readany()
                        async for chunk in resp.content.iter_any():
                            await fileOb.write(chunk)
                    cls.logger.debug(f"File has been saved in: {path}")
                    return path
                elif resp.status in (404, 405):
                    cls.logger.warning(f"404/405 for: {info}")
                    return None
                else:
                    mes = "code={resp.status}, message={resp.reason}, headers={resp.headers}".format(resp=resp)
                    cls.logger.error(mes)
                    raise Exception(mes)

    @classmethod
    @unsync
    async def save_file(cls, path: str, data: bytes):
        '''Deprecated'''
        async with aiofiles.open(path, 'wb') as fileOb:
            await fileOb.write(data)

    @classmethod
    @unsync
    async def fetch_file(cls, semaphore: asyncio.Semaphore, method: str, info: Dict, path: str, rate: float):
        try:
            async with semaphore:
                res = await cls.download_file(method, info, path)
                if res is not None:
                    await asyncio.sleep(rate)
                return res
        except RetryError:
            cls.logger.error(f"Retry failed for: {info}")

    @classmethod
    @unsync
    async def multi_tasks(cls, tasks: Union[Iterable, Iterator], to_do_func: Optional[Callable] = None, concur_req: int = 4, rate: float = 1.5, logger: Optional[logging.Logger] = None):
        '''
        Template for multiTasking

        TODO
            1. asyncio.Semaphore
            2. unit func
        '''
        cls.init_logger('UnsyncFetch', logger)
        cls.retry_kwargs['after'] = after_log(cls.logger, logging.WARNING)
        cls.download_file = retry(cls.download_file, **cls.retry_kwargs)
        semaphore = asyncio.Semaphore(concur_req)
        if to_do_func is None:
            # to_do_func = cls.default_func
            tasks = [cls.fetch_file(semaphore, method, info, path, rate) for method, info, path in tasks]
        else:
            tasks = [cls.fetch_file(semaphore, method, info, path, rate).then(to_do_func) for method, info, path in tasks]
        # return await asyncio.gather(*tasks)
        return [await fob for fob in tqdm(asyncio.as_completed(tasks), total=len(tasks))]

    @classmethod
    def main(cls, workdir: str, data: Union[Iterable, Iterator], concur_req: int = 4, rate: float = 1.5, logName: str = 'UnsyncFetch'):
        cls.set_logging_fileHandler(os.path.join(workdir, f'{logName}.log'), logName='UnsyncFetch')
        t0 = perf_counter()
        # res = asyncio.run(cls.multi_tasks(data, concur_req=concur_req, rate=rate))
        res = cls.multi_tasks(data, concur_req=concur_req, rate=rate).result()
        elapsed = perf_counter() - t0
        cls.logger.info(f'downloaded in {elapsed}s')
        return res
