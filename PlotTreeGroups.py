#Ryan A. Melnyk
#schmelnyk@gmail.com
#UBC Microbiology - Haney Lab

import argparse, os
import seaborn as sns
import numpy as np
from ete2 import Tree
from itertools import combinations
import cPickle as pickle
import scipy.cluster.hierarchy as hclust
import pandas as pd
from Bio import SeqIO
import matplotlib.pyplot as plt

def parse_args():
	parser = argparse.ArgumentParser(description='''
Using seaborn to plot heatmap of homologs for a given tree.
	''')
	parser.add_argument('tree_loc', type=str,help='path to tree file')
	parser.add_argument('outdir',type=str,help='PyParanoid output directory')
	parser.add_argument('queryfolder',type=str,help='folder produced by QueryClusters.py')
	return parser.parse_args()

def dist_from_tree(tree_loc):
	tree = Tree(tree_loc)
	leaves = tree.get_leaf_names()
	# dist_mat = np.zeros((len(leaves),len(leaves)))
	# for a,b in combinations(leaves,2):
	# 	d = tree.get_distance(a,b)
	# 	dist_mat[leaves.index(a),leaves.index(b)] = d
	# 	dist_mat[leaves.index(b),leaves.index(a)] = d
	# print dist_mat.shape
	# lnkg = hclust.linkage(dist_mat)
	# print lnkg.shape
	# f = open(os.path.join(os.path.dirname(tree_loc),os.path.basename(tree_loc).split(".")[0]+".linkage.pkl"),'w')
	# pickle.dump([leaves,lnkg],f)
	# f.close()
	print leaves
	return leaves

def make_matrix(outdir,strains,queryfolder):
	df = {}
	desc = {}
	groups = []
	for line in open(os.path.join(queryfolder,"hit_groups.txt")):
		vals = line.rstrip().split("\t")
		present = []
		groups.append(vals[1])
		for seq in SeqIO.parse(open(os.path.join(outdir,"homolog_fasta",vals[3]+".faa"),'r'),'fasta'):
			if str(seq.id).split("|")[0] not in present:
				present.append(str(seq.id).split("|")[0])
		for seq in SeqIO.parse(open(os.path.join(outdir,"prop_homolog_faa",vals[3]+".faa"),'r'),'fasta'):
			if str(seq.id).split("|")[0] not in present:
				present.append(str(seq.id).split("|")[0])
		df[vals[1]] = pd.Series(dict(zip(strains,[present.count(s) for s in strains])))
	return pd.DataFrame(df).reindex(strains)[groups]

def plot_clustergrid(df,queryfolder):
	fig, ax = plt.subplots()
	hm = sns.heatmap(df, linewidths=.4,cbar=False,ax=ax,square=True)
	# plt.setp(cg.ax_heatmap.yaxis.get_majorticklabels(), rotation=0)

	# hm.set_yticklabels(rotation=0, ha='right')

	ax.yaxis.tick_right()
	plt.yticks(rotation=0)
	plt.xticks([])
	plt.savefig('{}.png'.format(os.path.join(queryfolder,"clusterplot")),format='png',dpi=300)
	plt.savefig('{}.svg'.format(os.path.join(queryfolder,"clusterplot")),format='svg')
	return

def main():
	args = parse_args()
	tree_loc = os.path.abspath(args.tree_loc)
	outdir = os.path.abspath(args.outdir)
	queryfolder = os.path.abspath(args.queryfolder)

	leaves = dist_from_tree(tree_loc)
	df = make_matrix(outdir,leaves,queryfolder)
	plot_clustergrid(df,queryfolder)



if __name__ == '__main__':
	main()
