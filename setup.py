#!/usr/bin/env python

import setuptools
from distutils.core import setup

setup(
name='PyParanoid',
packages = ['PyParanoid'],
scripts = ['IdentifyOrthologs.py','PropagateGroups.py','PyParanoid.py','inparanoid2.pl'],
version='0.1',
description='Fast and scalable homolog detection for thousands of genomes',
long_description=open('README.md').read(),
author='Ryan A. Melnyk',
author_email='schmelnyk@gmail.com',
license='MIT',
url='https://github.com/ryanmelnyk/PyParanoid',
install_requires=['biopython', 'numpy', 'ncbi_genome_download', 'pandas', \
	'seaborn', 'matplotlib','reportlab','ijson','ete3'],
classifiers=['Programming Language :: Python','Programming Language :: Python :: 2']
)
