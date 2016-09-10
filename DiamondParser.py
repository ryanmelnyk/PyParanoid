#Ryan A. Melnyk
#schmelnyk@gmail.com
#UBC Microbiology - Haney Lab

import os, argparse
from Bio import SeqIO

def parse_args():
	parser = argparse.ArgumentParser(description='''
Script for converting pairwise diamond hits into a format for inparanoid.
	''')
	parser.add_argument('a', type=str,help='first file')
	parser.add_argument('b', type=str,help='second file')
	return parser.parse_args()

def parse_diamond(a,b, a_genes, b_genes):
	o = open("{}.{}.out".format(a,b),'w')
	for line in open("{}.{}.m8".format(a,b),'r'):
		vals = line.rstrip().split()
		newvals = [vals[0],vals[1],vals[11],str(a_genes[vals[0]]), str(b_genes[vals[1]])]
		newvals.append(str(int(vals[7])-int(vals[6])+1))
		newvals.append(str(int(vals[9])-int(vals[8])+1))
		newvals.append(str(int(vals[7])-int(vals[6])+1))
		newvals.append(str(int(vals[9])-int(vals[8])+1))
		newvals.append("p:{}-{} h:{}-{}".format(vals[6],vals[7],vals[8],vals[9]))
		o.write("\t".join(newvals)+"\n")

	return

def get_genes(a):
	a_genes = {}
	for seq in SeqIO.parse(open(a+".faa",'r'),'fasta'):
		a_genes[seq.id] = len(str(seq.seq))
	return a_genes

def main():
	args = parse_args()
	a = args.a
	b = args.b

	a_genes = get_genes(a)
	b_genes = get_genes(b)

	parse_diamond(a, a, a_genes, a_genes)
	parse_diamond(a, b, a_genes, b_genes)
	parse_diamond(b, a, b_genes, a_genes)
	parse_diamond(b, b, b_genes, b_genes)

if __name__ == '__main__':
	main()
