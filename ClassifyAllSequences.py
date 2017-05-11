#Ryan A. Melnyk
#schmelnyk@gmail.com
#UBC Microbiology - Haney Lab

import argparse, os
from Bio import SeqIO

def parse_args():
	parser = argparse.ArgumentParser(description='''
Script for checking all fasta protein sequences for presence in either the
initial homolog groups or the propagated calls.
	''')
	parser.add_argument('outdir', type=str,help='PyParanoid directory')
	parser.add_argument('genomedb',type=str,help="location of genomes")
	return parser.parse_args()

def parse(outdir,genomedb,a,b):
	count = len(os.listdir(os.path.join(outdir,a)))
	seqdict = {}
	print "Parsing files from {}...".format(a)
	for f in os.listdir(os.path.join(outdir,a)):
		for seq in SeqIO.parse(open(os.path.join(outdir,a,f),'r'),'fasta'):
			vals = seq.id.split("|")
			if vals[0] not in seqdict:
				seqdict[vals[0]] = [vals[1]]
			else:
				seqdict[vals[0]].append(vals[1])
		count -= 1
		if count == 0:
			print "\tDone!"
		elif count % 1000 == 0:
			print "\t"+str(count), "remaining..."
		else:
			pass

	o = open(os.path.join(outdir,"unclustered_{}_seqs.faa".format(a)),'w')

	unclust = 0
	totalseqs = 0
	strains = [line.rstrip() for line in open(os.path.join(outdir,b),"r")]
	count = len(strains)
	print "Parsing files from {}...".format(b)
	for s in strains:
		for seq in SeqIO.parse(open(os.path.join(genomedb,"pep",s+".pep.fa"),'r'),'fasta'):
			if seq.id not in seqdict[s]:
				o.write(">{}|{}\n{}\n".format(s,str(seq.id),str(seq.seq)))
				unclust += 1
			else:
				pass
			totalseqs += 1
		count -= 1
		if count == 0:
			print "\tDone!"
		elif count % 100 == 0:
			print "\t"+str(count), "remaining..."
		else:
			pass

	print "\n{}\n".format("="*60)
	print "For the {} portion:".format(a)
	print totalseqs, "total sequences in the {} genomes".format(str(len(strains)))
	pct = (float(totalseqs)-float(unclust))*100/float(totalseqs)
	print "{} sequences covered by the {} homolog groups: ({}%)".format(totalseqs-unclust,str(len(os.listdir(os.path.join(outdir,a)))),str(round(pct,1)))
	return


def main():
	args = parse_args()
	outdir = os.path.abspath(args.outdir)
	genomedb = os.path.abspath(args.genomedb)

	parse(outdir,genomedb,"homolog_fasta","strainlist.txt")
	parse(outdir,genomedb,"prop_homolog_faa","prop_strainlist.txt")

if __name__ == '__main__':
	main()
