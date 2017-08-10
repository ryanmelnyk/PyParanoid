#Ryan A. Melnyk
#schmelnyk@gmail.com
#UBC Microbiology - Haney Lab


import os
import argparse
import numpy as np
import seaborn as sns
import pandas as pd
import pyparanoid.pyparanoid as pp
import matplotlib.pyplot as plt

def parse_args():
	parser = argparse.ArgumentParser(description='''
Generate some statistics about the quality of the clustering.
	''')
	parser.add_argument('outdir', type=str,help='path to the PyParanoid output folder.')
	return parser.parse_args()

def parse_matrix(strains):
	matrixfile = open(os.path.join(outdir,"homolog_matrix.txt"),'r')
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
	poss_orthos = 0
	for i in range(1,l):
		y.append(np.sum(a[i,:]))
		x.append(i)
		singles = np.count_nonzero(a[i,:] == 1)
		if singles > a.shape[1]*0.95:
			poss_orthos += 1
	print poss_orthos, "groups are present in a single copy in >95% of strains."
	return x,y

def plot_group_sizes(x,y,a):
	x, y = pd.Series(x, name="gene families"), pd.Series(y, name="size of family")
	plt.axhline(y=a.shape[1])
	plt.xlim([0,len(x)])
	plt.ylim([0,y[1]*1.1])
	ax = sns.regplot(x=x,y=y,fit_reg=False)
	plt.savefig(os.path.join(outdir,"QCstats",'Group_Sizes.png'),dpi=300)
	plt.close()
	return

def sample_combinations(a):
	l = a.shape[1]
	x = []
	y = []
	for i in range(1,l):
		if i % 10 == 0:
			print "on number {}".format(i)
		for k in range(0,6):
			 x.append(i)
			 y.append(np.unique(np.nonzero(a[:,np.random.choice(l,size=i,replace=False)])[0]).shape[0])

	return x,y

def plot_rarefaction(x,y):
	x, y = pd.Series(x, name="strains sampled"), pd.Series(y, name="unique gene families")
	plt.xlim([0,max(x)*1.1])
	plt.ylim([0,max(y)*1.1])
	ax = sns.regplot(x=x,y=y,fit_reg=False)
	plt.savefig(os.path.join(outdir,"QCstats","Rarefaction.png"),dpi=300)
	plt.close()
	return

def main():
	args = parse_args()
	global outdir
	outdir = os.path.abspath(args.outdir)
	pp.createdirs(outdir,["QCstats"])

	strains = [line.rstrip() for line in open(os.path.join(outdir,"strainlist.txt"),'r')]
	a = parse_matrix(strains)
	x,y = count_groupsizes(a)
	plot_group_sizes(x,y,a)
	x,y = sample_combinations(a)
	plot_rarefaction(x,y)

if __name__ == '__main__':
	main()
