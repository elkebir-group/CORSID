#! /usr/bin/env python
# -*- coding: utf-8 -*-
import setuptools

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()


setuptools.setup(
    name='corsid',
    packages=["corsid"],
    description="Core Sequence Identifier",
    long_description=long_description,
    long_description_content_type="text/markdown",
    version='0.1.2',
    url='http://github.com/elkebir-group/CORSID',
    author='Chuanyi Zhang',
    author_email='chuanyi5@illinois.edu',
    license='MIT',
    python_requires='>=3.7',
    scripts=[
        'scripts/corsid',
        'scripts/corsid_a',
    ],
    install_requires=[
        "numpy",
        "pysam",
        "pandas",
        "pytablewriter",
        "tqdm",
    ],
)