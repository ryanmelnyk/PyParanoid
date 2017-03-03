#Ryan A. Melnyk
#schmelnyk@gmail.com
#UBC Microbiology - Haney Lab

import argparse, os
import numpy as np
import itertools
import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt

def parse_args():
	parser = argparse.ArgumentParser(description='''
For a completed PyParanoid folder, plot the size of the groups to look at orthologs.
	''')
	parser.add_argument('outdir', type=str,help='path to PyParanoid output')
	return parser.parse_args()

def parse_matrix(strains,pypdir):
	matrixfile = open(os.path.join(pypdir,"homolog_matrix.txt"),'r')
	header = matrixfile.readline().rstrip().split("\t")
	indices = [header.index(s) for s in strains]

	lines = []
	for line in matrixfile:
		vals = line.rstrip().split("\t")
		# convert input data into boolean presence/absence
		lines.append([int(vals[i]) for i in indices])
	a = np.stack(lines)
	return a

def count_groupsizes(a):
	l = a.shape[0]
	x,y = [],[]
	for i in range(1,l):
		# for k in range(0,20):
		# 	 x.append(i)
		# 	 y.append(np.unique(np.nonzero(a[:,np.random.choice(l,size=i,replace=False)])[0]).shape[0])
		y.append(np.sum(a[i,:]))
		x.append(i)
	return x,y

def plot(x,y,a):
	x, y = pd.Series(x, name="gene families"), pd.Series(y, name="size of family")
	print x.head()
	print y.head()
	plt.axhline(y=a.shape[1])
	plt.axhline(y=91,color='r',linewidth=1)
	plt.xlim([0,len(x)])
	plt.ylim([0,y[1]*1.1])
	ax = sns.regplot(x=x,y=y,fit_reg=False)
	plt.savefig('test.png',dpi=300)
	return

def main():
	args = parse_args()
	pypdir = os.path.abspath(args.outdir)
	strains = [line.rstrip() for line in open(os.path.join(pypdir,"strainlist.txt"),'r')]

	print len(strains)
	a = parse_matrix(strains,pypdir)
	x,y = count_groupsizes(a)
	plot(x,y,a)


if __name__ == '__main__':
	main()
