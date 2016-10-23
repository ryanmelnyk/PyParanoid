#Ryan A. Melnyk
#schmelnyk@gmail.com
#UBC Microbiology - Haney Lab

import os, argparse, subprocess, errno
from Bio import SeqIO

def parse_args():
	parser = argparse.ArgumentParser(description='''
Takes a complete PyParanoid directory (base and propagated) and generate list of orthologs.
	''')
	parser.add_argument('outdir', type=str,help='path to PyParanoid folder')
	return parser.parse_args()


def parse_matrix(outdir):
	orthos = []
	for line in open(os.path.join(outdir,"homolog_matrix.txt")):
		vals = line.rstrip().split("\t")
		if vals[0] == "":
			continue
		else:
			if set(vals[1:]) == set(["1"]):
				orthos.append(vals[0])
	print len(orthos), "orthologs found."
	return orthos

def extract_hmms(orthos,outdir):
	for o in orthos:
		cmds = "hmmfetch -o {} {} {}".format(os.path.join(outdir,"hmms",o+".hmm"),os.path.join(outdir,"all_groups.hmm"),o)
		proc = subprocess.Popen(cmds.split())
		proc.wait()
	return

def align_hmms(orthos, outdir):
	for o in orthos:
		cmds = ""
		### TODO
	return

def setupdir(outdir):
	for f in ["aligned"]:
		try:
			os.makedirs(os.path.join(outdir,f))
		except OSError as exc:
			if exc.errno == errno.EEXIST:
				print "Database folder exists:", os.path.join(outdir,f)
	return

def main():
	args = parse_args()
	outdir = os.path.abspath(args.outdir)

	setupdir(outdir)
	orthos = parse_matrix(outdir)
	# extract_hmms(orthos, outdir)
	align_hmms(orthos,outdir)

if __name__ == '__main__':
	main()
