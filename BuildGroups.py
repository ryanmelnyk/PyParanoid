#! /usr/bin/env python
#Ryan A. Melnyk
#schmelnyk@gmail.com
#UBC Microbiology - Haney Lab

import argparse
import os
import sys
import errno
import shutil
import subprocess
import itertools
import re
import pyparanoid.pyparanoid as pp
import multiprocessing as mp
from Bio import SeqIO

def parse_args():
	parser = argparse.ArgumentParser(description='''
Master script for running PyParanoid process.

modes: setup, parse, cluster, extract
	''')
	parser.add_argument('genomedb',type=str,help='relative path to genomedb for raw data')
	parser.add_argument('strainlist',type=str,help='text file, one strain per line. names should be "species" field from genome_metadata table')
	parser.add_argument('outdir',type=str,help='folder to store output')
	parser.add_argument('--cpus',type=int,help='number of CPUs to use for tasks. Defaults to # of cores available.')
	parser.add_argument('--mode',type=str,help='mode of PyParanoid to run')
	parser.add_argument('--clean',action="store_true",help="clean up intermediate files")
	parser.add_argument('--threshold',type=int,help="minimum size of groups")
	parser.add_argument('--inflate', type=float,help="inflation parameter for mcl, default = 2")
	parser.add_argument('--verbose',action="store_true",help="Print progress to STDOUT")
	parser.add_argument('--multi',action="store_true",help="use only with mode setup and parse")
	parser.add_argument('--use_MP',action="store_true",help="use the python multiprocessing module to dramatically speed up certain steps")
	return parser.parse_args()

def setupdir(strains,genomedb):
	try:
		os.makedirs(outdir)
	except OSError as exc:
		if exc.errno == errno.EEXIST:
			print "Database folder exists:", outdir

	pp.createdirs(outdir, ["faa","m8","out","paranoid_output","dmnd_tmp"])

	if not os.path.isdir(genomedb):
		print "GenomeDB folder", genomedb, "doesn't exist...exiting."
		sys.exit()

	if not os.path.isdir(os.path.join(genomedb,"pep")):
		print "GenomeDB folder is missing a 'pep' folder...exiting."
		sys.exit()

	if verbose:
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
	if verbose:
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

def run_diamond(strains):
	if multi:
		jobfile = open(os.path.join(outdir,"diamond_jobs.sh"),'w')
		if verbose:
			print "Making jobfile for DIAMOND runs..."
	else:
		if verbose:
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
				if verbose:
					print "\tDone!"
			elif count % 100 == 0:
				if verbose:
					print "\t"+str(count), "remaining..."
			else:
				pass
	return

def parse_diamond(genes,strains):

	count = len(strains)
	if verbose:
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
		if clean:
			os.remove(os.path.join(outdir,"m8",s+".m8"))
		count -= 1
		if count == 0:
			if verbose:
				print "\tDone!"
		elif count % 100 == 0:
			if verbose:
				print "\t"+str(count), "remaining..."
		else:
			pass
	return

def get_genes(strains):
	if verbose:
		print "Getting gene lengths..."
	genes = {}
	for s in strains:
		genes[s] = {}
		for seq in SeqIO.parse(open("{}/faa/{}.faa".format(outdir,s),'r'),'fasta'):
			genes[s][seq.id] = len(str(seq.seq))
	return genes

def run_inparanoid(strains):
	pairs = list(itertools.combinations(strains,2))

	if multi:
		jobfile = open(os.path.join(outdir,"inparanoid_jobs.sh"),'w')
		if verbose:
			print "Prepping InParanoid batch files..."
		for p in pairs:
			cmds = "inparanoid2.pl {} {} {}".format(p[0],p[1],outdir+"/")
			jobfile.write("{}\n".format(cmds))
		jobfile.close()
	else:
		count = len(pairs)
		if verbose:
			print "Running InParanoid on", count, "pairs of strains..."
		if use_MP:
			pool = mp.Pool(processes=cpus)
			[pool.apply_async(IP_RUN, args=(p,)) for p in pairs]
			pool.close()
			pool.join()
		else:
			if verbose:
				print "Sequential mode..."
			count = len(pairs)
			for p in pairs:
				IP_RUN(p)
				count -= 1
				if verbose:
					if count % 100 == 0:
						print "\t"+str(count), "remaining..."
	return

def IP_RUN(p):
	cmds = "inparanoid2.pl {} {} {}".format(p[0],p[1],outdir+"/")
	proc = subprocess.Popen(cmds.split())
	proc.wait()
	return

def hash_fastas():
	seqdata = {}
	desc = {}
	count = 0
	for f in os.listdir(os.path.join(outdir,"faa")):
		for seq in SeqIO.parse(open(os.path.join(outdir,"faa",f),'r'),'fasta'):
			seqdata[str(seq.id)] = str(seq.seq)
			match = re.search("(description:)(.*)",str(seq.description))
			if match == None:
				desc[str(seq.id)] = str(seq.description).replace("MULTISPECIES: ","").split(" [")[0]
			else:
				desc[str(seq.id)] = match.group(2)

			count += 1
	return seqdata, desc, count

def create_abc_file():
	o = open(os.path.join(outdir,"mcl","input.abc"),'w')
	count = len(os.listdir(os.path.join(outdir,"paranoid_output")))
	if verbose:
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
			if verbose:
				print "\tDone!"
		elif count % 1000 == 0:
			if verbose:
				print "\t"+str(count), "remaining..."
		else:
			pass
		if clean:
			os.remove(os.path.join(outdir,"paranoid_output",f))
	o.close()
	return

def mcxload():
	cmds = ["mcxload","--stream-mirror","-abc",os.path.join(outdir,"mcl","input.abc"),"-o",os.path.join(outdir,"mcl","data.mci"), "-write-tab",os.path.join(outdir,"mcl","data.tab")]
	proc = subprocess.Popen(cmds)
	proc.wait()
	return

def mcl_cluster():
	cmds = ["mcl",os.path.join(outdir,"mcl","data.mci"),"-te",str(cpus),"-I",str(inflate),"-o", os.path.join(outdir,"mcl","mcl.out")]
	proc = subprocess.Popen(cmds)
	proc.wait()
	return

def mcxdump():
	cmds = ["mcxdump", "-imx",os.path.join(outdir,"mcl","data.mci"),"-tabr",os.path.join(outdir,"mcl","data.tab"),"-icl", os.path.join(outdir,"mcl","mcl.out"),"-o", os.path.join(outdir,"mcl","clusters.txt")]
	proc = subprocess.Popen(cmds)
	proc.wait()
	return

def parse_clusters(strains, seq_number):
	o = open(os.path.join(outdir,"mcl","clusterstats.out"),'w')
	orthologs = 0
	paralogs = 0
	nghs = 0
	singletons = 0
	lengths = []
	coverage = 0
	for line in open(os.path.join(outdir,"mcl","clusters.txt"),"r"):
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
	o.write("{} orthologs.\n".format(orthologs))
	o.write("{} paralogs.\n".format(paralogs))
	o.write("{} non-global homologs\n".format(nghs))
	o.write("{} singletons.\n".format(singletons))
	o.write("="*60+"\n")
	o.write("Largest cluster: {} sequences.\n".format(sorted(lengths,reverse=True)[0]))
	o.write("Coverage: {}% of sequence DB\n".format(round(float(coverage)*100/float(seq_number),1)))
	o.write("="*60+"\n")
	o.write("{} total homolog groups.\n".format(len(lengths)))
	o.close()
	return

def parse_groups(seqdata, desc):
	group_count = 1
	descript_out = open(os.path.join(outdir, "group_descriptions.txt"),'w')
	print "Writing fasta files and parsing descriptions..."
	for line in open(os.path.join(outdir,"mcl","clusters.txt")):
		vals = line.rstrip().split()
		if len(vals) <= threshold:
			continue
		else:
			o = open(os.path.join(outdir,"homolog_faa","group_{}.faa".format(str(group_count).zfill(5))),'w')
			this_descript = []
			for v in vals:
				o.write(">{}\n{}\n".format(v,seqdata[v]))
				this_descript.append(desc[v])
			descript_out.write("group_{}\t{}\n".format(str(group_count).zfill(5),"\t".join(this_descript)))
			o.close()
			group_count += 1
	if verbose:
		print group_count, "groups equal to or larger than", threshold, "sequences."
	return

def cdhit_seqs():
	if verbose:
		print "Clustering sequences..."
	if use_MP:
		pool = mp.Pool(processes=cpus)
		[pool.apply_async(CDHIT_RUN, args=(f,)) for f in os.listdir(os.path.join(outdir,"homolog_faa"))]
		pool.close()
		pool.join()
	else:
		if verbose:
			print "Sequential mode..."
		count = len(os.listdir(os.path.join(outdir,"homolog_faa")))
		for f in os.listdir(os.path.join(outdir,"homolog_faa")):
			CDHIT_RUN(f)
			count -= 1
			if verbose:
				if count % 1000 == 0:
					print "\t"+str(count), "remaining..."
	return

def CDHIT_RUN(f):
	FNULL = open(os.devnull,'w')
	cmds = "cd-hit -i {} -o {} -c 0.95".format(os.path.join(outdir,"homolog_faa",f),os.path.join(outdir,"clustered",f))
	proc = subprocess.Popen(cmds.split(),stdout=FNULL,stderr=FNULL)
	proc.wait()
	FNULL.close()
	return

def align_groups():
	if verbose:
		print "Aligning groups..."
	files = []
	for f in os.listdir(os.path.join(outdir, "clustered")):
		if f.endswith(".clstr"):
			continue
		else:
			files.append(f)
	if use_MP:
		pool = mp.Pool(processes=cpus)
		[pool.apply_async(ALIGN_RUN, args=(f,)) for f in files]
		pool.close()
		pool.join()
	else:
		if verbose:
			print "Sequential mode..."
		count = len(files)
		for f in files:
			ALIGN_RUN(f)
			count -= 1
			if verbose:
				if count % 1000 == 0:
					print "\t"+str(count), "remaining..."
	return

def ALIGN_RUN(f):
	FNULL = open(os.devnull,'w')
	cmds = "muscle -in {} -out {}".format(os.path.join(outdir,"clustered",f),os.path.join(outdir,"aligned",f.split(".")[0]+'.aln'))
	proc = subprocess.Popen(cmds.split(),stdout=FNULL,stderr=FNULL)
	proc.wait()
	FNULL.close()
	return

def build_hmms():
	if verbose:
		print "Building hmms..."
	if use_MP:
		pool = mp.Pool(processes=cpus)
		[pool.apply_async(HMMBUILD_RUN, args=(f,)) for f in os.listdir(os.path.join(outdir, "aligned"))]
		pool.close()
		pool.join()
	else:
		if verbose:
			print "Sequential mode..."
		count = len(os.listdir(os.path.join(outdir, "aligned")))
		for f in os.listdir(os.path.join(outdir, "aligned")):
			HMMBUILD_RUN(f)
			count -= 1
			if verbose:
				if count % 1000 == 0:
					print "\t"+str(count), "remaining..."
	return

def HMMBUILD_RUN(f):
	FNULL = open(os.devnull,'w')
	cmds = "hmmbuild --cpu 1 {} {}".format(os.path.join(outdir,"hmms",f.split(".")[0]+".hmm"),os.path.join(outdir,"aligned",f))
	proc = subprocess.Popen(cmds.split(),stdout=FNULL,stderr=FNULL)
	proc.wait()
	FNULL.close()
	return

def emit_consensus_seqs():
	if verbose:
		print "Emitting consensus sequences..."
	if use_MP:
		pool = mp.Pool(processes=cpus)
		[pool.apply_async(HMMEMIT_RUN, args=(f,)) for f in os.listdir(os.path.join(outdir, "hmms"))]
		pool.close()
		pool.join()
	else:
		if verbose:
			print "Sequential mode..."
		count = len(os.listdir(os.path.join(outdir, "hmms")))
		for f in os.listdir(os.path.join(outdir, "hmms")):
			HMMEMIT_RUN(f)
			count -= 1
			if verbose:
				if count % 1000 == 0:
					print "\t"+str(count), "remaining..."

	return

def HMMEMIT_RUN(f):
	FNULL = open(os.devnull,'w')
	cmds = "hmmemit -c -o {} {}".format(os.path.join(outdir,"consensus_seqs",f.split(".")[0]+".faa"),os.path.join(outdir,"hmms",f))
	proc = subprocess.Popen(cmds.split(),stdout=FNULL,stderr=FNULL)
	proc.wait()
	FNULL.close()
	return

def combine_seqs():
	if verbose:
		print "Writing multi-hmm file..."
	o = open(os.path.join(outdir,"all_groups.hmm"),'w')
	for f in os.listdir(os.path.join(outdir,"hmms")):
		o.write(open(os.path.join(outdir,"hmms",f),'r').read())
	o.close()
	if verbose:
		print "Writing multi-fasta consensus file..."
	p = open(os.path.join(outdir,"all_groups.faa"),'w')
	for f in os.listdir(os.path.join(outdir,"consensus_seqs")):
		p.write(open(os.path.join(outdir,"consensus_seqs",f),'r').read())
	p.close()
	return

def combine_homologs():
	o = open(os.path.join(outdir,"homolog.faa"),'w')
	for f in os.listdir(os.path.join(outdir, "homolog_faa")):
		group_id = f.split(".")[0]
		for seq in SeqIO.parse(open(os.path.join(outdir,"homolog_faa",f),'r'),'fasta'):
			seq.id = seq.id+"|"+group_id
			seq.description = ""
			SeqIO.write(seq,o,'fasta')
	o.close()
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
		if args.mode not in ["multi_setup","parse","extract","cluster"]:
			print "Unknown mode!!! Exiting..."
			sys.exit()

	global cpus
	if args.cpus:
		cpus = args.cpus
	else:
		cpus = mp.cpu_count()

	global clean
	if args.clean:
		clean = True
	else:
		clean = False

	global verbose
	if args.verbose:
		verbose = True
	else:
		verbose = False

	global inflate
	if args.inflate:
		inflate = args.inflate
	else:
		inflate = 2.0

	global threshold
	if args.threshold:
		threshold = args.threshold
	else:
		threshold = 0

	global multi
	if args.multi:
		multi = True
	else:
		multi = False

	global use_MP
	if args.use_MP:
		use_MP = True
	else:
		use_MP = False

	if not args.mode or args.mode == "multi_setup":
		setupdir(strains,genomedb)
		shutil.copy(os.path.abspath(args.strainlist),os.path.join(outdir,"strainlist.txt"))
		make_diamond_databases(strains)
		run_diamond(strains)
	if not args.mode or args.mode == "parse":
		genes = get_genes(strains)
		parse_diamond(genes,strains)
		run_inparanoid(strains)
	if not args.mode or args.mode == "cluster":
		if clean:
			pp.cleanup(os.path.join(outdir,"out"))
		pp.createdirs(outdir,["mcl"])
		create_abc_file()
		mcxload()
		mcl_cluster()
		mcxdump()
	if not args.mode or args.mode == "extract":
		seqdata, desc, seq_number = hash_fastas()
		pp.createdirs(outdir,["homolog_faa","clustered","aligned","hmms","consensus_seqs"])
		parse_clusters(strains,seq_number)
		parse_groups(seqdata,desc)
		cdhit_seqs()
		align_groups()
		if clean:
			pp.cleanup(os.path.join(outdir,"clustered"))
		build_hmms()
		if clean:
			pp.cleanup(os.path.join(outdir,"aligned"))
		emit_consensus_seqs()
		combine_seqs()
		combine_homologs()
		if clean:
			pp.cleanup(os.path.join(outdir,"hmms"))
			pp.cleanup(os.path.join(outdir,"consensus_seqs"))
			pp.cleanup(os.path.join(outdir,"m8"))
			pp.cleanup(os.path.join(outdir,"paranoid_output"))
			pp.cleanup(os.path.join(outdir,"dmnd_tmp"))
			pp.cleanup(os.path.join(outdir,"faa"))
			pp.cleanup(os.path.join(outdir,"homolog_faa"))
			pp.cleanup(os.path.join(outdir,"mcl"))
			os.remove(os.path.join(outdir,"all_strains.dmnd"))
		pp.dump_matrices(outdir)


if __name__ == '__main__':
	main()
