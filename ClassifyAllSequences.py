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
	return parser.parse_args()

def parse(outdir,a,b):
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
	count = len(os.listdir(os.path.join(outdir,b)))
	print "Parsing files from {}...".format(b)
	for f in os.listdir(os.path.join(outdir,b)):
		if f.startswith("CONSENSUS"):
			continue
		for seq in SeqIO.parse(open(os.path.join(outdir,b,f),'r'),'fasta'):
			vals = seq.id.split("|")
			if vals[1] not in seqdict[vals[0]]:
				o.write(">{}\n{}\n".format(str(seq.id),str(seq.seq)))
				unclust += 1
			else:
				pass
			totalseqs += 1
		count -= 1
		if count == 0:
			print "\tDone!"
		elif count % 10 == 0:
			print "\t"+str(count), "remaining..."
		else:
			pass

	print "\n{}\n".format("="*60)
	print "For the {} portion:".format(a)
	print totalseqs, "total sequences in the {} genomes".format(str(len(os.listdir(os.path.join(outdir,b)))))
	pct = (float(totalseqs)-float(unclust))*100/float(totalseqs)
	print "{} sequences covered by the {} homolog groups: ({}%)".format(totalseqs-unclust,str(len(os.listdir(os.path.join(outdir,a)))),str(round(pct,1)))
	return


def main():
	args = parse_args()
	outdir = os.path.abspath(args.outdir)

	parse(outdir,"homolog_fasta","faa")
	parse(outdir,"prop_homolog_faa","prop_faa")

if __name__ == '__main__':
	main()
