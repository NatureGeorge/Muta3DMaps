# @Created Date: 2019-11-24 09:07:03 pm
# @Filename: setup.py
# @Email:  1730416009@stu.suda.edu.cn
# @Author: ZeFeng Zhu
# @Last Modified: 2019-12-23 04:27:14 pm
# @Copyright (c) 2019 MinghuiGroup, Soochow University
from setuptools import setup, find_packages, find_namespace_packages
import Muta3DMaps


with open("README.md", "rt") as f:
    readme = f.read()


setup(
      name="Muta3DMaps",
      version=Muta3DMaps.__version__,

      packages=find_packages(),
      package_data={
        '': ['*.ini'],
      },
      entry_points='''
        [console_scripts]
        Muta3DMaps=Muta3DMaps.Run:interface
      ''',
      install_requires=[
        'Click',
        'aiohttp>=3.6.2',
        'pandas>=0.25.0',
        'numpy>=1.16.0',
        'biopython==1.73',
        'wget>=3.2',
        'retrying>=1.3.0'
     ],
      license="MIT",
      author_email="minghui.li@suda.edu.cn",
      maintainer="ZeFeng Zhu",
      maintainer_email="1730416009@stu.suda.edu.cn",
      description="A Python Package that retrieves, filtering and organizes data from various databases or tools via API to process the residue-level mapping between protein sequences and protein 3D structures.",
      long_description=readme,
      long_description_content_type="text/markdown",
      url="https://github.com/NatureGeorge/Muta3DMaps",
      python_requires=">=3.5.*",
      classifiers=[
        "Operating System :: OS Independent",
        "Programming Language :: Python",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.5",
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3.7",
    ],
    # packages=find_namespace_packages(where='src'),
    # include_package_data=True,
    )

"""
Packaged By

python setup.py sdist bdist_wheel
twine upload --repository-url https://test.pypi.org/legacy/ dist/*

> https://packaging.python.org/tutorials/packaging-projects/
"""
