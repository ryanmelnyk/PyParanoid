#Ryan A. Melnyk
#schmelnyk@gmail.com
#UBC Microbiology - Haney Lab

import argparse, os, subprocess
from Bio import SeqIO
from operator import itemgetter

def parse_args():
	parser = argparse.ArgumentParser(description='''
Given a set of fasta queries, determine if there is an associated homolog group.
	''')
	parser.add_argument('outdir', type=str,help='path to pyparanoid folder')
	parser.add_argument('query',type=str,help='path to fasta file query (protein)')
	return parser.parse_args()

def run_hmmsearch(outdir,query):
	cmds = "hmmsearch --tblout {} {} {}".format(os.path.join(outdir,os.path.basename(query).split(".")[0],"hmmhits.out"),os.path.join(outdir,"all_groups.hmm"),query)
	proc = subprocess.Popen(cmds.split())
	proc.wait()
	return

def setupdir(outdir,query):
	for f in [os.path.basename(query).split(".")[0]]:
		try:
			os.makedirs(os.path.join(os.path.join(outdir,f)))
		except OSError:
			print "Subfolder exists:", outdir
	return

def parse_hmmsearch(outdir,query):
	tophits = {}
	for line in open(os.path.join(outdir,os.path.basename(query).split(".")[0],"hmmhits.out"),'r'):
		if line.startswith("#"):
			continue
		vals = line.split()
		if vals[0] not in tophits:
			tophits[vals[0]] = (float(vals[5]),vals[2]," ".join(vals[18:]))
		elif float(vals[5]) > tophits[vals[0]][0]:
			tophits[vals[0]] = (float(vals[5]),vals[2]," ".join(vals[18:]))
		else:
			pass

	return tophits

def get_gene_lengths(query):
	gene_lengths = {}
	for seq in SeqIO.parse(open(query,"r"),'fasta'):
		gene_lengths[str(seq.id)] = len(seq.seq)
	return gene_lengths

def write_matrix(outdir,query,gene_lengths,tophits):

	o = open(os.path.join(outdir,os.path.basename(query).split(".")[0],"hit_groups.txt"),'w')
	for t in sorted(tophits.items(), key=lambda x: x[1][1]):
		#### SORT ISN'T WORKING
		o.write("{}\t{}\t{}\n".format(t[0],str(t[1][0]/float(gene_lengths[t[0]])),"\t".join([str(x) for x in t[1]])))

	return

def main():
	args = parse_args()
	outdir = os.path.abspath(args.outdir)
	query = os.path.abspath(args.query)

	setupdir(outdir,query)
	run_hmmsearch(outdir,query)
	tophits = parse_hmmsearch(outdir,query)
	gene_lengths = get_gene_lengths(query)

	write_matrix(outdir,query,gene_lengths,tophits)

if __name__ == '__main__':
	main()
