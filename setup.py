#! /usr/bin/env python
# -*- coding: utf-8 -*-
import setuptools

setuptools.setup(
    name='corsid',
    packages=["corsid"],
    description="Core Sequence Identifier",
    version='0.1.0',
    url='http://github.com/elkebir-group/CORSID ',
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
    ],
)