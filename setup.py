#!/usr/bin/env python

import setuptools
from distutils.core import setup

setup(
name='PyParanoid',
packages = ['pyparanoid'],
scripts = ['IdentifyOrthologs.py','PropagateGroups.py','BuildGroups.py','inparanoid2.pl'],
version='0.3',
description='Fast and scalable homolog detection for thousands of genomes',
long_description='PyParanoid is a set of scripts and functions for identifying homology groups in large databases of microbial genomes. See https://github.com/ryanmelnyk/PyParanoid for details.',
author='Ryan A. Melnyk',
author_email='schmelnyk@gmail.com',
license='MIT',
url='https://github.com/ryanmelnyk/PyParanoid',
install_requires=['biopython', 'numpy', 'ncbi_genome_download', 'pandas', \
	'seaborn', 'matplotlib','reportlab','ijson','ete3'],
classifiers=['Programming Language :: Python','Programming Language :: Python :: 2']
)
