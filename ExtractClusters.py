#Ryan A. Melnyk
#schmelnyk@gmail.com
#UBC Microbiology - Haney Lab

import argparse, os, re
from Bio import SeqIO
import subprocess

def parse_args():
	parser = argparse.ArgumentParser(description='''
For a given clustering, extracts a fasta file for each ortholog group.
	''')
	parser.add_argument('outdir', type=str,help='location of PyParanoid folder')
	return parser.parse_args()

def hash_fastas(outdir):
	seqdata = {}
	desc = {}
	for f in os.listdir(os.path.join(outdir,"faa")):
		for seq in SeqIO.parse(open(os.path.join(outdir,"faa",f),'r'),'fasta'):
			seqdata[str(seq.id)] = str(seq.seq)
			vals = str(seq.description).split()
			if vals[2].startswith("pep:"):
				match = re.search("(description:)(.*)",str(seq.description))
				desc[str(seq.id)] = match.group(2)
			else:
				desc[str(seq.id)] = vals[2]
	return seqdata, desc

def parse_groups(seqdata, desc, outdir):
	group_count = 1
	descript_out = open(os.path.join(outdir, "group_descriptions.txt"),'w')
	print "Writing fasta files and parsing descriptions..."
	###
	# by default choose clustering where inflation == 2.0
	###
	for line in open(os.path.join(outdir,"mcl","clusters.2.txt")):
		vals = line.rstrip().split()
		o = open(os.path.join(outdir,"homolog_fasta","group_{}.faa".format(str(group_count).zfill(5))),'w')
		this_descript = []
		for v in vals:
			o.write(">{}\n{}\n".format(v,seqdata[v]))
			this_descript.append(desc[v])
		descript_out.write("group_{}\t{}\n".format(str(group_count).zfill(5),"\t".join(this_descript)))
		o.close()
		group_count += 1

	return

def setupdir(outdir):
	for f in ["homolog_fasta","clustered","aligned","hmms","consensus_seqs"]:
		try:
			os.makedirs(os.path.join(os.path.join(outdir,f)))
		except OSError:
			print "Subfolder exists:", os.path.join(outdir,f)
	return

def cluster_seqs(outdir):
	print "Clustering sequences..."
	for f in os.listdir(os.path.join(outdir,"homolog_fasta")):
		cmds = "cd-hit -i {} -o {} -c 0.95".format(os.path.join(outdir,"homolog_fasta",f),os.path.join(outdir,"clustered",f))
		proc = subprocess.Popen(cmds.split())
		proc.wait()
	return

def align_groups(outdir):
	print "Aligning groups..."
	for f in os.listdir(os.path.join(outdir, "clustered")):
		if f.endswith(".clstr"):
			continue
		cmds = "muscle -in {} -out {}".format(os.path.join(outdir,"clustered",f),os.path.join(outdir,"aligned",f.split(".")[0]+'.aln'))
		proc = subprocess.Popen(cmds.split())
		proc.wait()
	return

def build_hmms(outdir):
	print "Building hmms..."
	for f in os.listdir(os.path.join(outdir, "aligned")):
		cmds = "hmmbuild {} {}".format(os.path.join(outdir,"hmms",f.split(".")[0]+".hmm"),os.path.join(outdir,"aligned",f))
		proc = subprocess.Popen(cmds.split())
		proc.wait()
	return

def emit_consensus_seqs(outdir):
	print "Emitting consensus sequences..."
	for f in os.listdir(os.path.join(outdir,"hmms")):
		cmds = "hmmemit -c -o {} {}".format(os.path.join(outdir,"consensus_seqs",f.split(".")[0]+".faa"),os.path.join(outdir,"hmms",f))
		proc = subprocess.Popen(cmds.split())
		proc.wait()
	return

def main():
	args = parse_args()
	outdir = os.path.abspath(args.outdir)

	setupdir(outdir)
	seqdata, desc = hash_fastas(outdir)
	parse_groups(seqdata, desc, outdir)
	cluster_seqs(outdir)
	align_groups(outdir)
	build_hmms(outdir)
	emit_consensus_seqs(outdir)

if __name__ == '__main__':
	main()
