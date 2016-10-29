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
	print "Parsing matrix to identify orthologs..."
	for line in open(os.path.join(outdir,"homolog_matrix.txt")):
		vals = line.rstrip().split("\t")
		if vals[0] == "":
			continue
		else:
			if set(vals[1:]) == set(["1"]):
				orthos.append(vals[0])
	print len(orthos), "orthologs found."
	return orthos

def concat_orthos(orthos, outdir):
	count = len(orthos)
	print "Concatenating {} ortholog files...".format(str(count))
	for o in orthos:
		out = open(os.path.join(outdir,"concat",o+".faa"),'w')
		for seq in SeqIO.parse(open(os.path.join(outdir,"homolog_fasta",o+".faa"),'r'),'fasta'):
			out.write(">{}\n{}\n".format(str(seq.id).split("|")[0],str(seq.seq)))
		for seq in SeqIO.parse(open(os.path.join(outdir,"prop_homolog_faa",o+".faa"),'r'),'fasta'):
			out.write(">{}\n{}\n".format(str(seq.id).split("|")[0],str(seq.seq)))
		count -= 1
		if count == 0:
			print "\tDone!"
		elif count % 100 == 0:
			print "\t"+str(count), "remaining..."
		else:
			pass
	return

def align_orthos(orthos,outdir):
	count = len(orthos)
	print "Aligning {} ortholog files...".format(str(count))
	FNULL = open(os.devnull, 'w')
	for o in orthos:
		# cmds = "kalign {} {}".format(os.path.join(outdir,"concat",o+".faa"),os.path.join(outdir,"ortho_align",o+".fna"))
		cmds = "hmmalign -o {} {} {}".format(os.path.join(outdir,"ortho_align",o+".sto"),os.path.join(outdir,"hmms",o+".hmm"),os.path.join(outdir,"concat",o+".faa"))
		proc = subprocess.Popen(cmds.split(),stdout=FNULL,stderr=FNULL)
		proc.wait()
		count -= 1
		if count == 0:
			print "\tDone!"
		elif count % 100 == 0:
			print "\t"+str(count), "remaining..."
		else:
			pass
	FNULL.close()
	return

def setupdir(outdir):
	for f in ["ortho_align","concat"]:
		try:
			os.makedirs(os.path.join(outdir,f))
		except OSError as exc:
			if exc.errno == errno.EEXIST:
				print "Database folder exists:", os.path.join(outdir,f)
	return

def extract_hmms(orthos,outdir):
	count = len(orthos)
	print "Extracting {} HMM files...".format(str(count))
	FNULL = open(os.devnull, 'w')
	for o in orthos:
		# cmds = "kalign {} {}".format(os.path.join(outdir,"concat",o+".faa"),os.path.join(outdir,"ortho_align",o+".fna"))
		cmds = "hmmfetch -o {} {} {}".format(os.path.join(outdir,"hmms",o+".hmm"),os.path.join(outdir,"all_groups.hmm"),o)
		proc = subprocess.Popen(cmds.split(),stdout=FNULL,stderr=FNULL)
		proc.wait()
		count -= 1
		if count == 0:
			print "\tDone!"
		elif count % 100 == 0:
			print "\t"+str(count), "remaining..."
		else:
			pass
	FNULL.close()
	return

def create_master_alignment(orthos,outdir):
	strains = [line.rstrip() for line in open(os.path.join(outdir,"strainlist.txt"))]
	[strains.append(s) for s in [line.rstrip() for line in open(os.path.join(outdir,"prop_strainlist.txt"))]]
	align_data = {k : [] for k in strains}
	print align_data

	for o in orthos:
		for line in open(os.path.join(outdir,"ortho_align",o+".sto")):
			if line.startswith("#") or line.startswith("//"):
				continue
			else:
				vals = line.rstrip().split()
				if len(vals) < 1:
					continue
				elif vals[0] in align_data:
					align_data[vals[0]] += vals[1]
				else:
					align_data[vals[0]] = vals[1]

	o = open(os.path.join(outdir,"master_alignment.faa"),'w')
	for a in align_data:
		print a, len(align_data[a])
		o.write(">{}\n{}\n".format(a,"".join(align_data[a]).upper().replace(".","-")))
	o.close()

	return align_data

def main():
	args = parse_args()
	outdir = os.path.abspath(args.outdir)

	setupdir(outdir)
	orthos = parse_matrix(outdir)
	extract_hmms(orthos,outdir)
	concat_orthos(orthos,outdir)
	align_orthos(orthos,outdir)
	create_master_alignment(orthos,outdir)

if __name__ == '__main__':
	main()
