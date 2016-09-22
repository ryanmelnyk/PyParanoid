#Ryan A. Melnyk
#schmelnyk@gmail.com
#UBC Microbiology - Haney Lab

import argparse, os, errno, shutil, subprocess
import itertools

def parse_args():
	parser = argparse.ArgumentParser(description='''
Master script for running PyParanoid process.
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

	for f in ["faa","dmnd","m8"]:
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
	print "Making diamond databases..."
	for s in strains:
		cmds = "diamond makedb --in {}/faa/{}.faa -d {}/dmnd/{}.dmnd --quiet".format(outdir,s,outdir,s)
		proc = subprocess.Popen(cmds.split())
		proc.wait()
	return

def run_diamond(strains,outdir):
	print "Running reflexive diamond on all", len(strains), "strains..."
	count = len(strains)
	for s in strains:
		cmds = "diamond blastp --query {}/faa/{}.faa -d {}/dmnd/{}.dmnd -o {}/m8/{}.{}.m8 -f tab --min-score 50 --quiet".format(outdir,s,outdir,s,outdir,s,s)
		proc = subprocess.Popen(cmds.split())
		proc.wait()
		count -= 1
		print count, "remaining..."

	pairs = list(itertools.permutations(strains,2))
	count = len(pairs)
	print "Running reflexive diamond on all", count, "pairwise permutations..."
	for p in list(pairs):
		print p[0],p[1]
		cmds = "diamond blastp --query {}/faa/{}.faa -d {}/dmnd/{}.dmnd -o {}/m8/{}.{}.m8 -f tab --min-score 50 --quiet".format(outdir,p[0],outdir,p[1],outdir,p[0],p[1])
		print cmds
		proc = subprocess.Popen(cmds.split())
		proc.wait()
		count -= 1
		print count, "remaining..."

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



if __name__ == '__main__':
	main()
