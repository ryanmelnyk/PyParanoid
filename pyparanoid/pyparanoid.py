
#Ryan A. Melnyk
#schmelnyk@gmail.com
#UBC Microbiology - Haney Lab

import os
# from Bio import SeqIO

### Attempting to re-implement InParanoid algorithm

# def get_seqs(fastafile):
# 	count = 0
# 	seqdict = {}
# 	seqarr = []
# 	for s in SeqIO.parse(open(fastafile,'r'),'fasta'):
# 		count += 1
# 		seqdict[s.id] = count
# 		seqarr.append(s.id)
# 	return seqdict, seqarr
#
#
# def overlap_test(vals,seq_threshold=0.5):
# 	check = True
# 	if vals[3] >= vals[4]:
# 		if vals[5] < (seq_threshold * float(vals[3])):
# 			check = False
# 	if vals[4] >= vals[3]:
# 		if vals[6] < (seq_threshold * float(vals[4])):
# 			check = False
# 	return
#
# def parse_hits(blastfile):
# 	hits = {}
# 	for line in open(blastfile,'r'):
# 		vals = line.rstrip().split("\t")
# 		if overlap_test(vals):
# 			continue
# 		if vals[0] not in hits:
# 			hits[vals[0]] = {vals[1]:vals[2]}
# 		else:
# 			hits[vals[0]][vals[1]] = vals[2]
# 	return hits
#
# def get_top_hit(a):
# 	#return the key of the top-scoring hit in an id-score dict
# 	return sorted(a, key = lambda x: float(a[x]),reverse=True)[0]
#
# def find_rbhs(a,b):
# 	orths,inparalogs = {},{}
# 	for topA in a:
# 		topB = get_top_hit(a[topA])
# 		if get_top_hit(b[topB]) == topA:
# 			orths[(topA,topB)] = max(a[topA][topB],b[topB][topA])
# 			inparalogs[(topA,topB)] =
# 	return orths
#
# def InParanoid(a, b, outdir):
# 	print "Running InParanoid on {} + {}...".format(a,b)
#
# 	adict, aarr = get_seqs(os.path.join(outdir,"faa",a+".faa"))
# 	bdict, barr = get_seqs(os.path.join(outdir,"faa",b+".faa"))
#
# 	ab_hits = parse_hits(os.path.join(outdir,"out","{}.{}.out".format(a,b)))
# 	ba_hits = parse_hits(os.path.join(outdir,"out","{}.{}.out".format(b,a)))
#
# 	orths = find_rbhs(ab_hits, ba_hits)
# 	for o in orths:
# 		print o,orths[o]
# 	return

def createdirs(outdir, folders):
	for f in folders:
		try:
			os.makedirs(os.path.join(os.path.join(outdir,f)))
		except OSError:
			print "Subfolder exists:", os.path.join(outdir,f)
	return

def cleanup(d):
	print "Cleaning up", d
	for f in os.listdir(d):
		os.remove(os.path.join(d,f))
	return
