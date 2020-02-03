# @Created Date: 2020-01-18 10:53:07 am
# @Filename: unsyncTry.py
# @Email:  1730416009@stu.suda.edu.cn
# @Author: ZeFeng Zhu
# @Last Modified: 2020-02-03 06:09:14 pm
# @Copyright (c) 2020 MinghuiGroup, Soochow University
import os
from time import perf_counter
import datetime
import asyncio
import aiohttp
import aiofiles
from unsync import unsync

BASE_URL = "https://www.ebi.ac.uk/pdbe/api/pdb/entry/"
DEMO = [
    (BASE_URL+'molecules/1a01', '1a01_molecules.json'),
    (BASE_URL+'molecules/2xyn', '2xyn_molecules.json'),
    (BASE_URL+'molecules/1miu', '1miu_molecules.json'),
    (BASE_URL+'molecules/2hev', '2hev_molecules.json'),
    (BASE_URL+'molecules/3g96', '3g96_molecules.json'),
    (BASE_URL+'residue_listing/1a01', '1a01_residue_listing.json'),
    (BASE_URL+'residue_listing/2xyn', '2xyn_residue_listing.json'),
    (BASE_URL+'residue_listing/1miu', '1miu_residue_listing.json'),
    (BASE_URL+'residue_listing/2hev', '2hev_residue_listing.json'),
    (BASE_URL+'residue_listing/3g96', '3g96_residue_listing.json')]


@unsync
async def download_file(url: str):
    print(f"[{datetime.datetime.now()}] Start to get file: {url}")
    async with aiohttp.ClientSession() as session:
        async with session.get(url) as resp:
            if resp.status == 200:
                return await resp.read()
            elif resp.status == 404:
                return None
            else:
                mes = "code={resp.status}, message={resp.reason}, headers={resp.headers}"
                raise Exception(mes.format(resp=resp))


@unsync
async def save_file(path: str, data: bytes):
    print(f"[{datetime.datetime.now()}] Start to save file: {path}")
    async with aiofiles.open(path, 'wb') as jsonFile:
        await jsonFile.write(data)


@unsync
async def fetch_file(semaphore: asyncio.Semaphore, url: str, path: str, rate: float):
    async with semaphore:
        data = await download_file(url)
        await asyncio.sleep(rate)
        if data is not None:
            await save_file(path, data)
            return path


@unsync
async def multi_tasks(workdir: str, concur_req: int, rate=1.5):
    semaphore = asyncio.Semaphore(concur_req)
    await asyncio.gather(*[fetch_file(semaphore, url, os.path.join(workdir, path), rate) for url, path in DEMO])


if __name__ == "__main__":
    workdir = r"C:\OmicData\LiGroupWork\PDBeAPI\issue"
    t0 = perf_counter()
    multi_tasks(workdir, 4).result()
    elapsed = perf_counter() - t0
    print(f'downloaded in {elapsed}s')
