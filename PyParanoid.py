#Ryan A. Melnyk
#schmelnyk@gmail.com
#UBC Microbiology - Haney Lab

import argparse, os, errno, shutil, subprocess
import multiprocessing as mp
import itertools, sys
from Bio import SeqIO

def parse_args():
	parser = argparse.ArgumentParser(description='''
Master script for running PyParanoid process.

modes: multi_setup, multi_parse
	''')
	parser.add_argument('genomedb',type=str,help='relative path to genomedb for raw data')
	parser.add_argument('strainlist',type=str,help='text file, one strain per line. names should be "species" field from genome_metadata table')
	parser.add_argument('outdir',type=str,help='folder to store output')
	parser.add_argument('--cpus',type=int,help='number of CPUs to use for tasks. Defaults to # of cores available.')
	parser.add_argument('--mode',type=str,help='mode of PyParanoid to run')
	parser.add_argument('--debug',action="store_true",help="don't clean up intermediate files")
	return parser.parse_args()

def setupdir(strains,genomedb):
	try:
		os.makedirs(outdir)
	except OSError as exc:
		if exc.errno == errno.EEXIST:
			print "Database folder exists:", outdir

	for f in ["faa","m8","out","paranoid_output","mcl","dmnd_tmp"]:
		try:
			os.makedirs(os.path.join(os.path.join(outdir,f)))
		except OSError:
			print "Subfolder exists:", outdir

	print "Formatting",len(strains), "fasta files..."
	for s in strains:
		try:
			i = open(os.path.join(genomedb,"pep",s+".pep.fa"),"r")
		except IOError as exc:
			if exc.errno == 2:
				print s, 'not found in database...check your strainlist.'
				sys.exit()
		o = open(os.path.join(outdir,"faa",s+".faa"),"w")
		for seq in SeqIO.parse(i,'fasta'):
			seq.id = s+"|"+str(seq.id)
			SeqIO.write(seq,o,'fasta')
		o.close()
	return

def make_diamond_databases(strains):
	count = len(strains)
	print "Making diamond databases for", count, "strains..."
	o = open(os.path.join(outdir,"all_strains.faa"),'w')
	for s in strains:
		o.write(open("{}/faa/{}.faa".format(outdir,s),'r').read())
	o.close()

	cmds = "diamond makedb --in {} -d {} --quiet --threads {}".format(os.path.join(outdir,"all_strains.faa"),os.path.join(outdir,"all_strains.dmnd"),cpus)
	proc = subprocess.Popen(cmds.split())
	proc.wait()
	os.remove(os.path.join(outdir,"all_strains.faa"))
	return

def run_diamond(strains,multi):
	if multi:
		jobfile = open(os.path.join(outdir,"diamond_jobs.sh"),'w')
		print "Making jobfile for DIAMOND runs..."
	else:
		print "Running DIAMOND on all", len(strains), "strains..."

	count = len(strains)
	for s in strains:
		cmds = "diamond blastp --query {}/faa/{}.faa -d {}/all_strains.dmnd -t {}/dmnd_tmp -o {}/m8/{}.m8 -f tab --max-target-seqs {} --min-score 50 --quiet --threads {}".format(outdir,s,outdir,outdir,outdir,s,len(strains)*10,cpus)
		if multi:
			jobfile.write("{}\n".format(cmds))
		else:
			proc = subprocess.Popen(cmds.split())
			proc.wait()
			count -= 1
			if count == 0:
				print "\tDone!"
			elif count % 20 == 0:
				print "\t"+str(count), "remaining..."
			else:
				pass
	return

def parse_diamond(genes,strains):

	count = len(strains)
	print "Parsing diamond results for {} strains...".format(count)
	for s in strains:
		hits = {}

		for line in open(os.path.join(outdir,"m8",s+".m8"),'r'):
			vals = line.rstrip().split()
			target_strain = vals[1].split("|")[0]
			if target_strain not in hits:
				hits[target_strain] = []
			newvals = [vals[0],vals[1],vals[11],str(genes[s][vals[0]]), str(genes[target_strain][vals[1]])]
			newvals.append(str(int(vals[7])-int(vals[6])+1))
			newvals.append(str(int(vals[9])-int(vals[8])+1))
			newvals.append(str(int(vals[7])-int(vals[6])+1))
			newvals.append(str(int(vals[9])-int(vals[8])+1))
			newvals.append("p:{}-{} h:{}-{}".format(vals[6],vals[7],vals[8],vals[9]))
			hits[target_strain].append(newvals)

		for ts in hits:
			o = open(os.path.join(outdir,"out","{}.{}.out".format(s,ts)),'w')
			for h in hits[ts]:
				o.write("{}\n".format("\t".join(h)))
			o.close()
		if not debug:
			os.remove(os.path.join(outdir,"m8",s+".m8"))
		count -= 1
		if count == 0:
			print "\tDone!"
		elif count % 20 == 0:
			print "\t"+str(count), "remaining..."
		else:
			pass
	return

def get_genes(strains):
	print "Getting gene lengths..."
	genes = {}
	for s in strains:
		genes[s] = {}
		for seq in SeqIO.parse(open("{}/faa/{}.faa".format(outdir,s),'r'),'fasta'):
			genes[s][seq.id] = len(str(seq.seq))
	return genes

def run_inparanoid(strains,multi):
	pairs = list(itertools.combinations(strains,2))

	if multi:
		jobfile = open(os.path.join(outdir,"inparanoid_jobs.sh"),'w')
		print "Prepping InParanoid batch files..."
		for p in pairs:
			cmds = "perl {}/inparanoid2.pl {} {} {}".format(pypath,p[0],p[1],outdir+"/")
			jobfile.write("{}\n".format(cmds))
		jobfile.close()
	else:
		count = len(pairs)
		print "Running InParanoid on", count, "pairs of strains..."
		pool = mp.Pool(processes=cpus)
		[pool.apply_async(IP_RUN, args=(p,pypath,outdir,)) for p in pairs]
		pool.close()
		pool.join()

	return

def IP_RUN(p):
	cmds = "perl {}/inparanoid2.pl {} {} {}".format(pypath,p[0],p[1],outdir+"/")
	proc = subprocess.Popen(cmds.split())
	proc.wait()
	return

def clean_up():
	if not debug:
		for f in os.listdir(os.path.join(outdir,"out")):
			os.remove(os.path.join(outdir,"out",f))
	return

def main():
	args = parse_args()
	genomedb = os.path.abspath(args.genomedb)
	strains = [x.rstrip() for x in open(os.path.abspath(args.strainlist),'r')]
	if len(set(strains)) != len(strains):
		print "Duplicate entry in strainlist! Exiting..."

	global outdir, pypath
	outdir = os.path.abspath(args.outdir)
	pypath = os.path.abspath(os.path.dirname(sys.argv[0]))

	if args.mode:
		if args.mode not in ["multi_setup","multi_parse"]:
			print "Unknown mode!!! Exiting..."
			sys.exit()

	global cpus
	if args.cpus:
		cpus = args.cpus
	else:
		cpus = mp.cpu_count()


	global debug
	if args.debug:
		debug = True
	else:
		debug = False


	if args.mode == "multi_setup":
		setupdir(strains,genomedb)
		shutil.copy(os.path.abspath(args.strainlist),os.path.join(outdir,"strainlist.txt"))
		make_diamond_databases(strains)
		run_diamond(strains,True)
	elif args.mode == "multi_parse":
		genes = get_genes(strains)
		parse_diamond(genes,strains)
		run_inparanoid(strains,True)
	else:
		setupdir(strains,genomedb)
		shutil.copy(os.path.abspath(args.strainlist),os.path.join(outdir,"strainlist.txt"))
		make_diamond_databases(strains)
		run_diamond(strains,False)
		genes = get_genes(strains)
		parse_diamond(genes,strains)
		run_inparanoid(strains,False)
		clean_up(outdir)

if __name__ == '__main__':
	main()
