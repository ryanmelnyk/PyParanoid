
#Ryan A. Melnyk
#schmelnyk@gmail.com
#UBC Microbiology - Haney Lab

import os
from Bio import SeqIO
import numpy as np

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
			os.makedirs(os.path.abspath(os.path.join(outdir,f)))
		except OSError:
			print "Subfolder exists:", os.path.join(outdir,f)
	return

def cleanup(d):
	print "Cleaning up", d
	for f in os.listdir(d):
		os.remove(os.path.join(d,f))
	os.rmdir(d)
	return

def dump_matrices(outdir):
	o_loc = open(os.path.join(outdir,"locustag_matrix.txt"),'w')
	o_hom = open(os.path.join(outdir,"homolog_matrix.txt"),'w')
	strains = [line.rstrip() for line in open(os.path.join(outdir,"strainlist.txt"),'r')]
	try:
		[strains.append(s) for s in [line.rstrip() for line in open(os.path.join(outdir,"prop_strainlist.txt"),'r')]]
	except IOError:
		pass

	o_loc.write("\t{}\n".format("\t".join(strains)))
	o_hom.write("\t{}\n".format("\t".join(strains)))

	groups = {}
	for seq in SeqIO.parse(open(os.path.join(outdir,"homolog.faa"),'r'),'fasta'):
		vals = str(seq.id).split("|")
		if vals[2] not in groups:
			groups[vals[2]] = {s : [] for s in strains}
		groups[vals[2]][vals[0]].append(vals[1])
	try:
		for seq in SeqIO.parse(open(os.path.join(outdir,"prop_homolog.faa"),'r'),'fasta'):
			vals = str(seq.id).split("|")
			if vals[2] not in groups:
				groups[vals[2]] = {s : [] for s in strains}
			groups[vals[2]][vals[0]].append(vals[1])
	except IOError:
		pass

	for g in sorted(groups.keys()):
		thisline = []
		counts = []
		for s in strains:
			if len(groups[g][s]) == 0:
				thisline.append("None")
				counts.append("0")
			elif len(groups[g][s]) == 1:
				thisline.append(groups[g][s][0])
				counts.append("1")
			else:
				thisline.append(";".join(groups[g][s]))
				counts.append(str(len(groups[g][s])))
		o_loc.write("{}\t{}\n".format(g,"\t".join(thisline)))
		o_hom.write("{}\t{}\n".format(g,"\t".join(counts)))
	o_loc.close()
	o_hom.close()

	return

def parse_matrix(outdir):
	strains = [line.rstrip() for line in open(os.path.join(outdir,"strainlist.txt"),'r')]
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

def get_groupsizes(relpath):
	outdir = os.path.abspath(relpath)
	a = parse_matrix(outdir)
	l = a.shape[0]
	x,y = [],[]
	poss_orthos = 0
	for i in range(0,l):
		y.append(np.sum(a[i,:]))
		x.append(i)
		singles = np.count_nonzero(a[i,:] == 1)
		if singles > a.shape[1]*0.95:
			poss_orthos += 1
	print l, "total groups in", os.path.basename(outdir)
	print poss_orthos, "groups are present in a single copy in >95% of strains."
	return (x,y,poss_orthos)

def get_rarefaction(relpath):
	outdir = os.path.abspath(relpath)
	a = parse_matrix(outdir)
	l = a.shape[1]
	x = []
	y = []
	for i in range(1,l):
		if i % 100 == 0:
			print "on genome number {} of {}".format(i,l)
		for k in range(0,3):
			 x.append(i)
			 y.append(np.unique(np.nonzero(a[:,np.random.choice(l,size=i,replace=False)])[0]).shape[0])
	print "Done!"
	return (x,y)
