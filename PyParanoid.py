#Ryan A. Melnyk
#schmelnyk@gmail.com
#UBC Microbiology - Haney Lab

import argparse, os, errno, shutil, subprocess
import itertools, sys
from Bio import SeqIO

def parse_args():
	parser = argparse.ArgumentParser(description='''
Master script for running PyParanoid process.

Note that default diamond uses 6 threads when requesting resources.
	''')
	parser.add_argument('genomedb',type=str,help='relative path to genomedb for raw data')
	parser.add_argument('strainlist',type=str,help='text file, one strain per line. names should be "species" field from genome_metadata table')
	parser.add_argument('outdir',type=str,help='folder to store output')
	parser.add_argument('--add',action='store_true',help='use if you intend to add to an existing PyParanoid output folder')
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
				print s, 'not found in database...check your strainlist.'
				sys.exit()
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

def run_diamond_on_old_strains(strains,outdir,old_strains):

	print "Running reflexive diamond on all", len(old_strains), "strains..."
	count = len(old_strains)
	for s in old_strains:
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

	count = 2 * len(strains) * len(old_strains)
	print "Running pairwise diamond against all old strains:", count, "pairwise permutations ..."
	for s in strains:
		for t in old_strains:
			cmds = "diamond blastp --query {}/faa/{}.faa -d {}/dmnd/{}.dmnd -o {}/m8/{}.{}.m8 -f tab --min-score 50 --quiet --threads 6".format(outdir,s,outdir,t,outdir,s,t)
			proc = subprocess.Popen(cmds.split())
			proc.wait()
			cmds = "diamond blastp --query {}/faa/{}.faa -d {}/dmnd/{}.dmnd -o {}/m8/{}.{}.m8 -f tab --min-score 50 --quiet --threads 6".format(outdir,t,outdir,s,outdir,t,s)
			proc = subprocess.Popen(cmds.split())
			proc.wait()
			count -= 2
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

	pairs = list(itertools.combinations(strains,2))
	count = len(pairs)
	print "Running inparanoid on", count, "pairs of strains..."
	for p in pairs:
		a = "{}.{}.txt".format(p[0],p[1])
		b = "{}.{}.txt".format(p[1],p[0])
		if p[0] == p[1]:
			pass
		elif a in os.listdir(os.path.join(outdir,"paranoid_output")):
			pass
		elif b in os.listdir(os.path.join(outdir,"paranoid_output")):
			pass
		else:
			cmds = "perl inparanoid2.pl {} {} {}".format(p[0],p[1],outdir)
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
		vals = f.split(".")
		if vals[0] == vals[1]:
			pass
		else:
			os.remove(os.path.join(outdir,"out",f))
	return

def append_strains(outdir,strainlist):
	o = open(os.path.join(outdir,"strainlist.txt"),'a')
	for line in open(strainlist,'r'):
		o.write(line)
	return

def check_unique(old_strains,strains):
	for s in strains:
		if s in old_strains:
			print "New strain already in database. Check strainlist and try again. Exiting..."
		else:
			pass
	return

def main():
	args = parse_args()
	genomedb = os.path.abspath(args.genomedb)
	strains = [x.rstrip() for x in open(os.path.abspath(args.strainlist),'r')]
	outdir = os.path.abspath(args.outdir)

	if args.add:
		print "adding"
		setupdir(outdir,strains,genomedb)
		old_strains = [x.rstrip() for x in open(os.path.join(outdir,"strainlist.txt"),'r')]
		check_unique(old_strains,strains)
		append_strains(outdir,os.path.abspath(args.strainlist))
		make_diamond_databases(strains,outdir)
		run_diamond(strains,outdir)
		run_diamond_on_old_strains(strains,outdir,old_strains)
		[strains.append(s) for s in old_strains]
		genes = get_genes(strains,outdir)
		parse_diamond(genes,outdir)
		run_inparanoid(strains, outdir)
		clean_up(outdir)
	else:
		print "not adding"
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
