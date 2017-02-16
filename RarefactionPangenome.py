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
For a completed PyParanoid folder, resample the whole genomes used to determine
how many strains it takes to cover the pangenome.
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
		lines.append([int(bool(int(vals[i]))) for i in indices])
	a = np.stack(lines)
	return a

def sample_combinations(a):
	l = a.shape[1]
	x = []
	y = []
	for i in range(1,l):
		if i % 10 == 0:
			print "on number {}".format(i)
		for k in range(0,20):
			 x.append(i)
			 y.append(np.unique(np.nonzero(a[:,np.random.choice(l,size=i,replace=False)])[0]).shape[0])

	return x,y

def plot(x,y):
	x, y = pd.Series(x, name="strains sampled"), pd.Series(y, name="genes")

	plt.xlim([0,220])
	plt.ylim([0,30000])
	ax = sns.regplot(x=x,y=y, color="g",order=0.5)
	plt.savefig('test.png',dpi=300)
	return

def main():
	args = parse_args()
	pypdir = os.path.abspath(args.outdir)
	strains = [line.rstrip() for line in open(os.path.join(pypdir,"strainlist.txt"),'r')]

	print len(strains)
	a = parse_matrix(strains,pypdir)
	x,y = sample_combinations(a)
	plot(x,y)


if __name__ == '__main__':
	main()
