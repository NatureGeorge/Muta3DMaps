import os
from time import perf_counter
import asyncio
import aiohttp
import aiofiles
from unsync import unsync
from tenacity import retry, wait_random_exponential, stop_after_attempt, after_log, RetryError, retry_if_exception_type
import logging
from tqdm import tqdm
from typing import Iterable, Iterator, Union


def MyException(info: str):
    raise ValueError(info)


class Test_UnsyncFetch(object):

    name = 'Test_tenacity_UnsyncFetch'
    logger = logging.getLogger(name)
    logger.setLevel(logging.DEBUG)
    streamHandler = logging.StreamHandler()
    streamHandler.setLevel(logging.INFO)
    formatter = logging.Formatter("%(asctime)s %(name)s %(levelname)s %(message)s")
    streamHandler.setFormatter(formatter)
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
            cls.logger.warning("Invalid file path! Please specifiy path=...")

    @classmethod
    @retry(retry=retry_if_exception_type(ValueError), wait=wait_random_exponential(multiplier=1, max=15), stop=stop_after_attempt(3), after=after_log(logger, logging.WARNING))
    async def download_file(cls, url: str):
        cls.logger.debug(f"Start to get file: {url}")
        '''
        async with aiohttp.ClientSession() as session:
            async with session.get(url) as resp:
                if resp.status == 200:
                    return await resp.read()
                elif resp.status == 404:
                    return None
                else:
                    mes = "code={resp.status}, message={resp.reason}, headers={resp.headers}"
                    raise Exception(mes.format(resp=resp))
        '''
        for ic in url:
            if ic == 'x':
                MyException("Invalid value")
                # pass
            elif ic == '6':
                await asyncio.sleep(.1)
                return None
            else:
                await asyncio.sleep(.05)
        return f'test={url}'

    @classmethod
    @unsync
    async def save_file(cls, path: str, data: bytes):
        cls.logger.debug(f"Start to save file: {path}")
        '''
        async with aiofiles.open(path, 'wb') as jsonFile:
            await jsonFile.write(data)
        '''
        await asyncio.sleep(.5)

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
        semaphore = asyncio.Semaphore(concur_req)
        # await asyncio.gather(*[cls.fetch_file(semaphore, url, os.path.join(workdir, path), rate) for url, path in tasks])
        tasks = [cls.fetch_file(semaphore, url, os.path.join(
            workdir, path), rate) for url, path in tasks]
        return [await fob for fob in tqdm(asyncio.as_completed(tasks), total=len(tasks))]

    @classmethod
    def main(cls, workdir: str, data: Union[Iterable, Iterator], logName=name):
        cls.set_logging_fileHandler(os.path.join(workdir, f'{logName}.log'))
        t0 = perf_counter()
        res = cls.multi_tasks(workdir, data, 5).result()
        elapsed = perf_counter() - t0
        cls.logger.info(f'downloaded in {elapsed}s')
        return res


DEMO = [
    ('https://www.ebi.ac.uk/pdbe/api/pdb/entry/molecules/1a01',
     '1a01_molecules.json'),
    ('https://www.ebi.ac.uk/pdbe/api/pdb/entry/molecules/2xyn',
     '2xyn_molecules.json'),
    ('https://www.ebi.ac.uk/pdbe/api/pdb/entry/molecules/1miu',
     '1miu_molecules.json'),
    ('https://www.ebi.ac.uk/pdbe/api/pdb/entry/molecules/2hev',
     '2hev_molecules.json'),
    ('https://www.ebi.ac.uk/pdbe/api/pdb/entry/molecules/3g96',
     '3g96_molecules.json'),
    ('https://www.ebi.ac.uk/pdbe/api/pdb/entry/residue_listing/1a01',
     '1a01_residue_listing.json'),
    ('https://www.ebi.ac.uk/pdbe/api/pdb/entry/residue_listing/2xyn',
     '2xyn_residue_listing.json'),
    ('https://www.ebi.ac.uk/pdbe/api/pdb/entry/residue_listing/1miu',
     '1miu_residue_listing.json'),
    ('https://www.ebi.ac.uk/pdbe/api/pdb/entry/residue_listing/2hev',
     '2hev_residue_listing.json'),
    ('https://www.ebi.ac.uk/pdbe/api/pdb/entry/residue_listing/3g96',
     '3g96_residue_listing.json')]


if __name__ == "__main__":
    Test_UnsyncFetch.main(r'./data/', DEMO)
