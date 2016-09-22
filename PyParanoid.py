#Ryan A. Melnyk
#schmelnyk@gmail.com
#UBC Microbiology - Haney Lab

import argparse, os, errno, shutil, subprocess
import itertools
from Bio import SeqIO

def parse_args():
	parser = argparse.ArgumentParser(description='''
Master script for running PyParanoid process.

Note that default diamond uses 6 threads when requesting resources.
	''')
	parser.add_argument('genomedb',type=str,help='relative path to genomedb for raw data')
	parser.add_argument('strainlist',type=str,help='text file, one strain per line. names should be "species" field from genome_metadata table')
	parser.add_argument('outdir',type=str,help='folder to store output')
	return parser.parse_args()

def setupdir(outdir,strains,genomedb):
	try:
		os.makedirs(outdir)
	except OSError as exc:
		if exc.errno == errno.EEXIST:
			print "Database folder exists:", outdir

	for f in ["faa","dmnd","m8","out","paranoid_output"]:
		try:
			os.makedirs(os.path.join(os.path.join(outdir,f)))
		except OSError:
			print "Subfolder exists:", outdir

	for s in strains:
		try:
			shutil.copy(os.path.join(genomedb,"pep",s+".pep.fa"),os.path.join(outdir,"faa",s+".faa"))
		except IOError as exc:
			if exc.errno == 2:
				print s, 'not found in database...'
	return

def make_diamond_databases(strains,outdir):
	count = len(strains)
	print "Making diamond databases for", count, "strains..."
	for s in strains:
		cmds = "diamond makedb --in {}/faa/{}.faa -d {}/dmnd/{}.dmnd --quiet --threads 6".format(outdir,s,outdir,s)
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

def run_diamond(strains,outdir):
	print "Running reflexive diamond on all", len(strains), "strains..."
	count = len(strains)
	for s in strains:
		cmds = "diamond blastp --query {}/faa/{}.faa -d {}/dmnd/{}.dmnd -o {}/m8/{}.{}.m8 -f tab --min-score 50 --quiet --threads 6".format(outdir,s,outdir,s,outdir,s,s)
		proc = subprocess.Popen(cmds.split())
		proc.wait()
		count -= 1
		if count == 0:
			print "\tDone!"
		elif (count % 10 == 0):
			print "\t"+str(count), "remaining..."
		else:
			pass

	pairs = list(itertools.permutations(strains,2))
	count = len(pairs)
	print "Running pairwise diamond on all", count, "pairwise permutations..."
	for p in list(pairs):
		cmds = "diamond blastp --query {}/faa/{}.faa -d {}/dmnd/{}.dmnd -o {}/m8/{}.{}.m8 -f tab --min-score 50 --quiet --threads 6".format(outdir,p[0],outdir,p[1],outdir,p[0],p[1])
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

def parse_diamond(genes,outdir):
	print "Parsing diamond results..."
	files = os.listdir(os.path.join(outdir,"m8"))
	count = len(files)
	for f in files:
		vals = f.split(".")
		a = vals[0]
		b = vals[1]
		o = open("{}/out/{}.{}.out".format(outdir,a,b),'w')
		for line in open(os.path.join(outdir,"m8",f),'r'):
			vals = line.rstrip().split()
			newvals = [vals[0],vals[1],vals[11],str(genes[a][vals[0]]), str(genes[b][vals[1]])]
			newvals.append(str(int(vals[7])-int(vals[6])+1))
			newvals.append(str(int(vals[9])-int(vals[8])+1))
			newvals.append(str(int(vals[7])-int(vals[6])+1))
			newvals.append(str(int(vals[9])-int(vals[8])+1))
			newvals.append("p:{}-{} h:{}-{}".format(vals[6],vals[7],vals[8],vals[9]))
			o.write("\t".join(newvals)+"\n")
			count -= 1
		if count == 0:
			print "\tDone!"
		elif count % 10 == 0:
			print "\t"+str(count), "remaining..."
		else:
			pass
	return

def get_genes(strains,outdir):
	print "Getting gene lengths..."
	genes = {}
	for s in strains:
		genes[s] = {}
		for seq in SeqIO.parse(open("{}/faa/{}.faa".format(outdir,s),'r'),'fasta'):
			genes[s][seq.id] = len(str(seq.seq))
	return genes

def run_inparanoid(strains,outdir):
	pairs = list(itertools.permutations(strains,2))
	count = len(pairs)
	print "Running InParanoid to identify orthologs for",count,"permutations..."
	for p in pairs:
		cmds = "perl inparanoid2.pl {} {} {}".format(p[0],p[1],outdir,outdir)
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

def clean_up(outdir):
	for f in os.listdir(os.path.join(outdir,"m8")):
		os.remove(os.path.join(outdir,"m8",f))
	for f in os.listdir(os.path.join(outdir,"out")):
		os.remove(os.path.join(outdir,"out",f))
	return

def main():
	args = parse_args()
	genomedb = os.path.abspath(args.genomedb)
	strains = [x.rstrip() for x in open(os.path.abspath(args.strainlist),'r')]
	outdir = os.path.abspath(args.outdir)

	setupdir(outdir,strains,genomedb)
	shutil.copy(os.path.abspath(args.strainlist),os.path.join(outdir,"strainlist.txt"))
	make_diamond_databases(strains,outdir)
	run_diamond(strains,outdir)
	genes = get_genes(strains,outdir)
	parse_diamond(genes,outdir)
	run_inparanoid(strains, outdir)
	clean_up(outdir)

if __name__ == '__main__':
	main()
