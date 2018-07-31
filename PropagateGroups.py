#! /usr/bin/env python
#Ryan A. Melnyk
#schmelnyk@gmail.com
#UBC Microbiology - Haney Lab

import os
import argparse
import sys
import subprocess
import shutil
import pyparanoid.pyparanoid as pp
import multiprocessing as mp
from Bio import SeqIO

def parse_args():
	parser = argparse.ArgumentParser(description='''
Using consensus sequences for each homolog group to propagate to new genomes.
	''')
	parser.add_argument('genomedb', type=str,help='path to genome database')
	parser.add_argument('new_strainlist',type=str,help='path to list of new strains')
	parser.add_argument('outdir', type=str,help='path to directory containing PyParanoid output')
	parser.add_argument('--cpus',type=int,help='number of cpus to use - defaults to # available')
	parser.add_argument('--use_MP',action="store_true",help="use the python multiprocessing module to dramatically speed up certain steps.")
	return parser.parse_args()

def check_strains(new_strains,genomedb):
	old_strains = [line.rstrip() for line in open(os.path.join(outdir,"strainlist.txt"),'r')]
	try:
		for line in open(os.path.join(outdir,"prop_strainlist.txt"),'r'):
			old_strains.append(line.rstrip())
	except IOError:
		pass

	for s in new_strains:
		if s in old_strains:
			print s, "is already in strainlist...check your strainlist."
			sys.exit()

	ns = open(os.path.join(outdir,"prop_strainlist.txt"),'a')
	for s in new_strains:
		try:
			i = open(os.path.join(genomedb,"pep",s+".pep.fa"),"r")
		except IOError as exc:
			if exc.errno == 2:
				print s, 'not found in database...check your strainlist.'
				sys.exit()
		o = open(os.path.join(outdir,"prop_faa",s+".faa"),"w")
		for seq in SeqIO.parse(i,'fasta'):
			seq.id = s+"|"+str(seq.id)
			SeqIO.write(seq,o,'fasta')
		o.close()
		ns.write("{}\n".format(s))
	ns.close()
	shutil.copy(os.path.join(outdir,"all_groups.faa"),os.path.join(outdir,"prop_faa","CONSENSUS.faa"))
	return

def make_diamond_databases(strains):
	count = len(strains)
	print "Making diamond databases for", count, "strains..."
	for s in strains:
		cmds = "diamond makedb --in {}/prop_faa/{}.faa -d {}/prop_dmnd/{}.dmnd --quiet --threads {}".format(outdir,s,outdir,s,cpus)
		proc = subprocess.Popen(cmds.split())
		proc.wait()
		count -= 1
		if count == 0:
			print "\tDone!"
		elif count % 10 == 0:
			print "\t"+str(count), "remaining..."
		else:
			pass
	cmds = "diamond makedb --in {}/all_groups.faa -d {}/all_groups.dmnd --quiet --threads {}".format(outdir,outdir,cpus)
	proc = subprocess.Popen(cmds.split())
	proc.wait()
	return

def run_diamond(strains):
	print "Running diamond on all", len(strains), "strains..."
	count = len(strains)
	for s in strains:
		cmds = "diamond blastp --query {}/prop_faa/{}.faa -d {}/prop_dmnd/{}.dmnd -o {}/prop_m8/{}.{}.m8 -f tab --min-score 50 --quiet --threads {}".format(outdir,s,outdir,s,outdir,s,s,cpus)
		proc = subprocess.Popen(cmds.split())
		proc.wait()
		cmds = "diamond blastp --query {}/prop_faa/{}.faa -d {}/all_groups.dmnd -o {}/prop_m8/{}.CONSENSUS.m8 -f tab --min-score 50 --quiet --threads {}".format(outdir,s,outdir,outdir,s,cpus)
		proc = subprocess.Popen(cmds.split())
		proc.wait()
		cmds = "diamond blastp --query {}/all_groups.faa -d {}/prop_dmnd/{}.dmnd -o {}/prop_m8/CONSENSUS.{}.m8 -f tab --min-score 50 --quiet --threads {}".format(outdir,outdir,s,outdir,s,cpus)
		proc = subprocess.Popen(cmds.split())
		proc.wait()
		count -= 1
		if count == 0:
			print "\tDone!"
		elif (count % 10 == 0):
			print "\t"+str(count), "remaining..."
		else:
			pass
	cmds = "diamond blastp --query {}/all_groups.faa -d {}/all_groups.dmnd -o {}/prop_m8/CONSENSUS.CONSENSUS.m8 -f tab --min-score 50 --quiet --threads {}".format(outdir,outdir,outdir,cpus)
	proc = subprocess.Popen(cmds.split())
	proc.wait()
	return

def parse_diamond(genes):
	print "Parsing diamond results..."
	files = os.listdir(os.path.join(outdir,"prop_m8"))
	count = len(files)
	for f in files:
		vals = f.split(".")
		a = vals[0]
		b = vals[1]
		o = open("{}/prop_out/{}.{}.out".format(outdir,a,b),'w')
		for line in open(os.path.join(outdir,"prop_m8",f),'r'):
			vals = line.rstrip().split()
			newvals = [vals[0],vals[1],vals[11],str(len(genes[a][vals[0]])), str(len(genes[b][vals[1]]))]
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

def get_genes(strains):
	print "Getting gene lengths..."
	genes = {}
	for s in strains:
		genes[s] = {}
		for seq in SeqIO.parse(open("{}/prop_faa/{}.faa".format(outdir,s),'r'),'fasta'):
			genes[s][seq.id] = str(seq.seq)
	genes["CONSENSUS"] = {}
	for seq in SeqIO.parse(open("{}/all_groups.faa".format(outdir,s),'r'),'fasta'):
		genes["CONSENSUS"][seq.id] = str(seq.seq)
	return genes

def run_inparanoid(strains,pypath):
	print "Running inparanoid on", len(strains), "strains..."
	count = len(strains)
	if use_MP:
		pool = mp.Pool(processes=cpus)
		[pool.apply_async(IP_RUN, args=(s,)) for s in strains]
		pool.close()
		pool.join()
	else:
		for s in strains:
			IP_RUN(s)
	return

def IP_RUN(s):
	cmds = "inparanoid2.pl {} {} {}".format(s,"CONSENSUS",outdir+"/prop_")
	proc = subprocess.Popen(cmds.split())
	proc.wait()
	return

def parse_inparanoid(new_strains):
	group_members = {}
	for s in new_strains:
		for line in open(os.path.join(outdir,"prop_paranoid_output","{}.CONSENSUS.txt".format(s))):
			if line.startswith("Orto"):
				continue
			vals = line.rstrip().split()
			hits = []
			for v in vals:
				if v.startswith(s):
					hits.append(v)
				if v.startswith("group_"):
					group = v.split("-")[0]
					break
			if group not in group_members:
				group_members[group] = hits
			else:
				[group_members[group].append(h) for h in hits]

	return group_members

def extract_fastas(genes,group_members):
	o = open(os.path.join(outdir,"prop_homolog.faa"),'a')
	for g in group_members:
		for h in group_members[g]:
			o.write(">{}|{}\n{}\n".format(h,g.split("-")[0],genes[h.split("|")[0]][h]))
	o.close()
	return

def main():
	args = parse_args()
	global pypath
	pypath = os.path.abspath(os.path.dirname(sys.argv[0]))
	global outdir
	outdir = os.path.abspath(args.outdir)
	genomedb = os.path.abspath(args.genomedb)
	new_strains = [line.rstrip() for line in open(os.path.abspath(args.new_strainlist),'r')]

	global cpus
	if args.cpus:
		cpus = args.cpus
	else:
		cpus = mp.cpu_count()

	global use_MP
	if args.use_MP:
		use_MP = True
	else:
		use_MP = False

	pp.createdirs(outdir,["prop_faa","prop_dmnd","prop_m8","prop_out","prop_paranoid_output","prop_homolog_faa"])
	check_strains(new_strains,genomedb)
	make_diamond_databases(new_strains)
	run_diamond(new_strains)
	genes = get_genes(new_strains)
	parse_diamond(genes)
	run_inparanoid(new_strains,pypath)
	group_members = parse_inparanoid(new_strains)
	extract_fastas(genes,group_members)
	pp.dump_matrices(outdir)
	for f in ["prop_m8","prop_out","prop_dmnd","prop_paranoid_output","prop_faa","prop_homolog_faa"]:
		pp.cleanup(os.path.join(outdir,f))


if __name__ == '__main__':
	main()
