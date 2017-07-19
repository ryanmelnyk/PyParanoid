
#Ryan A. Melnyk
#schmelnyk@gmail.com
#UBC Microbiology - Haney Lab

import os
from Bio import SeqIO

def get_seqs(fastafile):
	count = 0
	for s in SeqIO.parse(open(fastafile),'r'):
		count += 1
	print count
	return


def InParanoid(a, b, outdir):
	print "Running InParanoid on {} + {}...".format(a,b)
	get_seqs(os.path.join(outdir,"faa",a+".faa"))

	return
