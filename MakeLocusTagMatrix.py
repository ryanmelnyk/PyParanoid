#Ryan A. Melnyk
#schmelnyk@gmail.com
#UBC Microbiology - Haney Lab


import os, argparse
from Bio import SeqIO

def parse_args():
	parser = argparse.ArgumentParser(description='''
Makes a matrix that contains locus tag information.
	''')
	parser.add_argument('outdir', type=str,help='path to PyParanoid directory')
	return parser.parse_args()

def dump_matrix(outdir):
	o = open(os.path.join(outdir,"locustag_matrix.txt"),'w')
	strains = [line.rstrip() for line in open(os.path.join(outdir,"strainlist.txt"),'r')]
	[strains.append(s) for s in [line.rstrip() for line in open(os.path.join(outdir,"prop_strainlist.txt"),'r')]]

	o.write("\t{}\n".format("\t".join(strains)))

	for f in os.listdir(os.path.join(outdir,"homolog_fasta")):
		hits = {s : [] for s in strains}
		for seq in SeqIO.parse(open(os.path.join(outdir,"homolog_fasta",f),'r'),'fasta'):
			hits[seq.id.split("|")[0]].append(seq.id.split("|")[1])
		try:
			for seq in SeqIO.parse(open(os.path.join(outdir,"prop_homolog_faa",f),'r'),'fasta'):
				hits[seq.id.split("|")[0]].append(seq.id.split("|")[1])
		except IOError:
			pass
		thisline = []
		for s in strains:
			if len(hits[s]) == 0:
				thisline.append("None")
			elif len(hits[s]) == 1:
				thisline.append(hits[s][0])
			else:
				thisline.append(";".join(hits[s]))
		o.write("{}\t{}\n".format(f.split(".")[0],"\t".join(thisline)))
	return

def main():
	args = parse_args()
	dump_matrix(os.path.abspath(args.outdir))

if __name__ == '__main__':
	main()
