#Ryan A. Melnyk
#schmelnyk@gmail.com
#UBC Microbiology - Haney Lab

import os, argparse, subprocess, errno, sys, string
from Bio import SeqIO

def parse_args():
	parser = argparse.ArgumentParser(description='''
Takes a complete PyParanoid directory (base and propagated) and generate list of orthologs. Using the 'threshold'
argument relaxes the cutoff and includes homologs that occur exactly once in some fraction of all strains (e.g. 90%).
	''')
	parser.add_argument('outdir', type=str,help='path to PyParanoid folder')
	parser.add_argument('--threshold',type=float,help='proportion of strains for a homolog')
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

def parse_threshold_matrix(outdir, t):
	orthos = []
	print "Parsing matrix to identify orthologs..."
	for line in open(os.path.join(outdir,"homolog_matrix.txt")):
		vals = line.rstrip().split("\t")
		if vals[0] == "":
			continue
		else:
			if float(vals.count("1"))/float(len(vals)-1) > t:
				orthos.append(vals[0])
	print len(orthos), "orthologs found."
	return orthos

def concat_orthos(orthos,outdir,strains):
	count = len(orthos)
	print "Concatenating {} ortholog files...".format(str(count))
	for o in orthos:
		selected = []
		out = open(os.path.join(outdir,"concat",o+".faa"),'w')
		try:
			for seq in SeqIO.parse(open(os.path.join(outdir,"homolog_fasta",o+".faa"),'r'),'fasta'):
				strain = str(seq.id).split("|")[0]
				if strain not in selected:
					out.write(">{}\n{}\n".format(strain,str(seq.seq)))
					selected.append(strain)
		except IOError:
			pass
		try:
			for seq in SeqIO.parse(open(os.path.join(outdir,"prop_homolog_faa",o+".faa"),'r'),'fasta'):
				strain = str(seq.id).split("|")[0]
				if strain not in selected:
					out.write(">{}\n{}\n".format(strain,str(seq.seq)))
					selected.append(strain)
		except IOError:
			pass
		out.close()
		count -= 1
		if count == 0:
			print "\tDone!"
		elif count % 10 == 0:
			print "\t"+str(count), "remaining..."
		else:
			pass
	return

def align_orthos(orthos,outdir):
	count = len(orthos)
	print "Aligning {} ortholog files...".format(str(count))
	for o in orthos:
		# cmds = "kalign {} {}".format(os.path.join(outdir,"concat",o+".faa"),os.path.join(outdir,"ortho_align",o+".fna"))
		cmds = "hmmalign -o {} {} {}".format(os.path.join(outdir,"ortho_align",o+".sto"),os.path.join(outdir,"hmms",o+".hmm"),os.path.join(outdir,"concat",o+".faa"))
		proc = subprocess.Popen(cmds.split())
		proc.wait()
		count -= 1
		if count == 0:
			print "\tDone!"
		elif count % 10 == 0:
			print "\t"+str(count), "remaining..."
		else:
			pass
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
	present = [f.split(".")[0] for f in os.listdir(os.path.join(outdir,"hmms"))]
	print "Extracting {} HMM files...{} already found.".format(str(count),str(len(present)))
	FNULL = open(os.devnull, 'w')
	for o in orthos:
		count -= 1
		if o in present:
			pass
		else:
			cmds = "hmmfetch -o {} {} {}".format(os.path.join(outdir,"hmms",o+".hmm"),os.path.join(outdir,"all_groups.hmm"),o)
			proc = subprocess.Popen(cmds.split(),stdout=FNULL,stderr=FNULL)
			proc.wait()
		if count % 10 == 0:
			print "\t"+str(count), "remaining..."
		else:
			pass
	if count == 0:
		print "\tDone!"
	FNULL.close()
	return

def create_master_alignment(orthos,outdir,strains):

	align_data = {k : [] for k in strains}
	count = len(orthos)
	total_leng = 0 ###DEBUG
	print "Creating master alignment...Parsing {} homologs...".format(str(count))
	for o in orthos:
		count -= 1
		present = []
		for line in open(os.path.join(outdir,"hmms",o+".hmm")):
			if line.startswith("LENG"):
				length = int(line.rstrip().split()[1])
				total_leng += length ###DEBUG
				break
		for line in open(os.path.join(outdir,"ortho_align",o+".sto")):
			if line.startswith("#") or line.startswith("//"):
				continue
			else:
				vals = line.rstrip().split()
				if len(vals) < 1:
					continue
				elif vals[0] in align_data:
					align_data[vals[0]].append(vals[1].translate(None,string.ascii_lowercase).replace(".",""))
					if vals[0] not in present:
						present.append(vals[0])
				else:
					align_data[vals[0]] = [vals[1].translate(None,string.ascii_lowercase).replace(".","")]
					if vals[0] not in present:
						present.append(vals[0])
		for s in strains:
			if s not in present:
				align_data[s].append("-"*length)
			if len("".join(align_data[s])) != total_leng:
				print s, "is short!"
				print total_leng, len("".join(align_data[s]))
				print align_data[s]
				print o
				sys.exit()
		if count % 10 == 0:
			print "\t"+str(count), "remaining..."
		else:
			pass
	print "Done!"
	print "Writing alignment..."
	o = open(os.path.join(outdir,"master_alignment.faa"),'w')
	for a in align_data:
		o.write(">{}\n{}\n".format(a,"".join(align_data[a]).upper().replace(".","-")))
	o.close()

	return align_data

def cleanup(orthos,outdir):
	for o in orthos:
		os.remove(os.path.join(outdir,"concat",o+".faa"))
		os.remove(os.path.join(outdir,"ortho_align",o+".sto"))
	return

def get_strains(outdir):
	strains = [line.rstrip() for line in open(os.path.join(outdir,"strainlist.txt"))]
	[strains.append(s) for s in [line.rstrip() for line in open(os.path.join(outdir,"prop_strainlist.txt"))]]
	return strains

def main():
	args = parse_args()
	outdir = os.path.abspath(args.outdir)

	setupdir(outdir)
	if args.threshold:
		orthos = parse_threshold_matrix(outdir, args.threshold)
	else:
		orthos = parse_matrix(outdir)

	extract_hmms(orthos,outdir)
	strains = get_strains(outdir)
	concat_orthos(orthos,outdir,strains)
	# align_orthos(orthos,outdir)
	create_master_alignment(orthos,outdir,strains)
	# cleanup(orthos,outdir)

if __name__ == '__main__':
	main()
