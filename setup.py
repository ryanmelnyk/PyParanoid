#!/usr/bin/env python

import setuptools
from distutils.core import setup

setup(
name='PyParanoid',
packages = ['PyParanoid'],
scripts = ['AddLocalGenome.py', 'ClusterStats.py','DownloadEnsemblGenomes.py', 'DownloadNcbiGenomes.py','DumpTaxInfo.py'\
	'IdentifyOrthologs.py','PropagateGroups.py','PyParanoid.py'],
version='0.1dev',
description='Fast and scalable homolog detection for thousands of genomes',
long_description=open('README.md').read()
author='Ryan A. Melnyk',
author_email='schmelnyk@gmail.com',
license='MIT',
url='https://github.com/ryanmelnyk/PyParanoid',
install_requires=['biopython', 'numpy', 'ncbi_genome_download', 'pandas', 'seaborn', 'matplotlib','reportlab'],
)
