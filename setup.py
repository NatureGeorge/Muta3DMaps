# @Date:   2019-11-24T19:44:32+08:00
# @Email:  1730416009@stu.suda.edu.cn
# @Filename: setup.py
# @Last modified time: 2019-11-25T01:32:39+08:00
from setuptools import setup, find_packages


with open("README.md", "rt") as f:
    readme = f.read()


setup(name="Muta3DMaps",
      version="1.0",

      packages=find_packages('Muta3DMaps'),  # 'src'
      package_dir={'': 'Muta3DMaps'},
      package_data={
        '': ['*.ini'],
      },
      entry_points='''
        [console_scripts]
        Muta3DMaps=Muta3DMaps.Run:interface
      ''',
      install_requires=[
        'Click',
        'pandas',
        'numpy',
        'biopython',
        'wget',
        'retrying'
     ],
      license="MIT",
      author_email="minghui.li@suda.edu.cn",
      maintainer="ZeFeng Zhu",
      maintainer_email="1730416009@stu.suda.edu.cn",
      description="A Python Package that retrieves, filtering and organizes data from various databases or tools via API to process the residue-level mapping between protein sequences and protein 3D structures.",
      long_description=readme,
      long_description_content_type="text/markdown",
      python_requires="!=3.0.*, !=3.1.*, !=3.2.*, !=3.3.*",
      classifiers=[
        "Operating System :: OS Independent",
        "Programming Language :: Python",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.4",
        "Programming Language :: Python :: 3.5",
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3.7",
    ],
    # packages=find_namespace_packages(where='src'),
    # include_package_data=True,
    )
