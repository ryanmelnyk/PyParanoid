#Ryan A. Melnyk
#schmelnyk@gmail.com
#UBC Microbiology - Haney Lab


import os, argparse, subprocess, itertools
import pandas as pd

from Bio import SeqIO

def parse_args():
	parser = argparse.ArgumentParser(description='''
Script for formatting the PyParanoid output for a format compatible with MCL.
Runs the clustering process as well, providing stats on the output.
	''')
	parser.add_argument('outdir', type=str,help='location of PyParanoid output')
	parser.add_argument('--graphical',action='store_true',help='use if you want graphical output - don\'t use on bugaboo')
	return parser.parse_args()

def create_abc_file(outdir):
	o = open(os.path.join(outdir,"mcl","input.abc"),'w')
	count = len(os.listdir(os.path.join(outdir,"paranoid_output")))
	print "Parsing", count, "output files."
	for f in os.listdir(os.path.join(outdir,"paranoid_output")):
		for line in open(os.path.join(outdir,"paranoid_output",f),'r'):
			if line.startswith("Orto"):
				continue
			else:
				genes = {}
				vals = line.rstrip().split()
				index = 2
				while index < len(vals):
					genes[vals[index]] = float(vals[index+1])
					index += 2
				pairs = list(itertools.combinations(genes.keys(),2))
				for p in list(pairs):
					o.write("{}\t{}\t{}\n".format(p[0],p[1],min(genes[p[0]],genes[p[1]])))
		count -= 1
		if count == 0:
			print "\tDone!"
		elif count % 100 == 0:
			print "\t"+str(count), "remaining..."
		else:
			pass
	o.close()
	return

def clean_up(outdir):
	for f in os.listdir(os.path.join(outdir,"paranoid_output")):
		os.remove(os.path.join(outdir,"paranoid_output",f))
	return

def run_mcxload(outdir):
	cmds = ["mcxload","--stream-mirror","-abc",os.path.join(outdir,"mcl","input.abc"),"-o",os.path.join(outdir,"mcl","data.mci"), "-write-tab",os.path.join(outdir,"mcl","data.tab")]
	proc = subprocess.Popen(cmds)
	proc.wait()
	return

def cluster(outdir):
	for i in ["1.4","2","4","6"]:
		cmds = ["mcl",os.path.join(outdir,"mcl","data.mci"),"-t","8","-I",i,"-o", os.path.join(outdir,"mcl","mcl.{}.out".format(i.replace(".","")))]
		proc = subprocess.Popen(cmds)
		proc.wait()
	return

def dump(outdir):
	for i in ["1.4","2","4","6"]:
		cmds = ["mcxdump", "-imx",os.path.join(outdir,"mcl","data.mci"),"-tabr",os.path.join(outdir,"mcl","data.tab"),"-icl", os.path.join(outdir,"mcl","mcl.{}.out".format(i.replace(".",""))),"-o", os.path.join(outdir,"mcl","clusters.{}.txt".format(i.replace(".","")))]
		proc = subprocess.Popen(cmds)
		proc.wait()
	return

def parse_clusters(outdir, strains, seq_number):
	d = {}

	for i in ["1.4","2","4","6"]:
		o = open(os.path.join(outdir,"mcl","clusterstats.{}.out".format(i.replace(".",""))),'w')
		orthologs = 0
		paralogs = 0
		nghs = 0
		singletons = 0
		lengths = []
		coverage = 0
		for line in open(os.path.join(outdir,"mcl","clusters.{}.txt".format(i.replace(".",""))),"r"):
			vals = line.rstrip().split()
			coverage += len(vals)
			if set([v.split("|")[0] for v in vals]) == set(strains):
				if len(vals) > len(strains):
					paralogs += 1
				elif len(vals) == len(strains):
					orthologs += 1
			else:
				if len(vals) == 1:
					singletons += 1
				else:
					nghs += 1
			lengths.append(len(vals))
		d[i] = pd.Series(sorted(lengths,reverse=True))
		o.write("{} orthologs.\n".format(orthologs))
		o.write("{} paralogs.\n".format(paralogs))
		o.write("{} non-global homologs\n".format(nghs))
		o.write("{} singletons.\n".format(singletons))
		o.write("="*60+"\n")
		o.write("Largest cluster: {} sequences.\n".format(d[i][0]))
		o.write("Coverage: {}% of sequence DB\n".format(round(float(coverage)*100/float(seq_number),1)))
		o.write("="*60+"\n")
		o.write("{} total homolog groups.\n".format(len(lengths)))
		o.close()
	return pd.DataFrame(d)

def get_number_of_seqs(outdir, strains):
	count = 0
	for s in strains:
		for seq in SeqIO.parse(open(os.path.join(outdir,"faa",s+".faa"),'r'),'fasta'):
			count += 1
	print count, "total sequences in database."
	return count

def main():
	args = parse_args()
	outdir = os.path.abspath(args.outdir)
	strains = [line.rstrip() for line in open(os.path.join(outdir, "strainlist.txt"))]

	create_abc_file(outdir)
	clean_up(outdir)
	run_mcxload(outdir)
	cluster(outdir)
	dump(outdir)

	seq_number = get_number_of_seqs(outdir, strains)
	df = parse_clusters(outdir, strains, seq_number)

	if args.graphical:
		import matplotlib.pyplot as plt
		fig = plt.figure()
		df.plot()
		plt.savefig(os.path.join(outdir, "mcl",'ortholog_groups.png'),dpi=300)


if __name__ == '__main__':
	main()
