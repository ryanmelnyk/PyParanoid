#Ryan A. Melnyk
#schmelnyk@gmail.com
#UBC Microbiology - Haney Lab

import os, argparse, sys, subprocess, shutil
from Bio import SeqIO

def parse_args():
	parser = argparse.ArgumentParser(description='''
Using consensus sequences for each homolog group to propagate to new genomes.
	''')
	parser.add_argument('outdir', type=str,help='path to directory containing PyParanoid output')
	parser.add_argument('genomedb', type=str,help='path to genome database')
	parser.add_argument('new_strains',type=str,help='path to list of new strains')
	return parser.parse_args()

def check_strains(new_strains,genomedb,outdir):
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

def setupdir(outdir):
	for f in ["prop_faa","prop_dmnd","prop_m8","prop_out","prop_paranoid_output","prop_homolog_faa"]:
		try:
			os.makedirs(os.path.join(os.path.join(outdir,f)))
		except OSError:
			print "Subfolder exists:", outdir
	return

def make_diamond_databases(strains,outdir):
	count = len(strains)
	print "Making diamond databases for", count, "strains..."
	for s in strains:
		cmds = "diamond makedb --in {}/prop_faa/{}.faa -d {}/prop_dmnd/{}.dmnd --quiet --threads 6".format(outdir,s,outdir,s)
		proc = subprocess.Popen(cmds.split())
		proc.wait()
		count -= 1
		if count == 0:
			print "\tDone!"
		elif count % 10 == 0:
			print "\t"+str(count), "remaining..."
		else:
			pass
	cmds = "diamond makedb --in {}/all_groups.faa -d {}/all_groups.dmnd --quiet --threads 6".format(outdir,outdir)
	proc = subprocess.Popen(cmds.split())
	proc.wait()
	return

def run_diamond(strains,outdir):
	print "Running diamond on all", len(strains), "strains..."
	count = len(strains)
	for s in strains:
		cmds = "diamond blastp --query {}/prop_faa/{}.faa -d {}/prop_dmnd/{}.dmnd -o {}/prop_m8/{}.{}.m8 -f tab --min-score 50 --quiet --threads 6".format(outdir,s,outdir,s,outdir,s,s)
		proc = subprocess.Popen(cmds.split())
		proc.wait()
		cmds = "diamond blastp --query {}/prop_faa/{}.faa -d {}/all_groups.dmnd -o {}/prop_m8/{}.CONSENSUS.m8 -f tab --min-score 50 --quiet --threads 6".format(outdir,s,outdir,outdir,s)
		proc = subprocess.Popen(cmds.split())
		proc.wait()
		cmds = "diamond blastp --query {}/all_groups.faa -d {}/prop_dmnd/{}.dmnd -o {}/prop_m8/CONSENSUS.{}.m8 -f tab --min-score 50 --quiet --threads 6".format(outdir,outdir,s,outdir,s)
		proc = subprocess.Popen(cmds.split())
		proc.wait()
		count -= 1
		if count == 0:
			print "\tDone!"
		elif (count % 10 == 0):
			print "\t"+str(count), "remaining..."
		else:
			pass
	cmds = "diamond blastp --query {}/all_groups.faa -d {}/all_groups.dmnd -o {}/prop_m8/CONSENSUS.CONSENSUS.m8 -f tab --min-score 50 --quiet --threads 6".format(outdir,outdir,outdir)
	proc = subprocess.Popen(cmds.split())
	proc.wait()
	return

def parse_diamond(genes,outdir):
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

def get_genes(strains,outdir):
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

def run_inparanoid(strains,outdir):
	print "Running inparanoid on", len(strains), "strains..."
	count = len(strains)
	for s in strains:
		cmds = "perl inparanoid2.pl {} {} {}".format(s,"CONSENSUS",outdir+"/prop_")
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
	for f in os.listdir(os.path.join(outdir,"prop_m8")):
		os.remove(os.path.join(outdir,"prop_m8",f))
	for f in os.listdir(os.path.join(outdir,"prop_out")):
		vals = f.split(".")
		if vals[0] == vals[1]:
			pass
		else:
			os.remove(os.path.join(outdir,"prop_out",f))
	return

def parse_inparanoid(outdir,new_strains):
	group_members = {}
	print new_strains
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

	for g in sorted(group_members.keys(), key = lambda x: int(x.split("_")[1])):
		print g, group_members[g]


	return group_members

def extract_fastas(outdir,genes,group_members):
	for g in group_members:
		o = open(os.path.join(outdir,"prop_homolog_faa",g.split("-")[0]+".faa"),'w')
		for h in group_members[g]:
			o.write(">{}\n{}\n".format(h,genes[h.split("|")[0]][h]))
		o.close()
	return

def dump_matrix(outdir):
	o = open(os.path.join(outdir,"homolog_matrix.txt"),'w')
	strains = [line.rstrip() for line in open(os.path.join(outdir,"strainlist.txt"),'r')]
	[strains.append(s) for s in [line.rstrip() for line in open(os.path.join(outdir,"prop_strainlist.txt"),'r')]]
	print strains
	print len(strains)

	o.write("\t{}\n".format("\t".join(strains)))

	for f in os.listdir(os.path.join(outdir,"homolog_fasta")):
		hits = []
		for seq in SeqIO.parse(open(os.path.join(outdir,"homolog_fasta",f),'r'),'fasta'):
			hits.append(seq.id.split("|")[0])
		try:
			for seq in SeqIO.parse(open(os.path.join(outdir,"prop_homolog_faa",f),'r'),'fasta'):
				hits.append(seq.id.split("|")[0])
		except IOError:
			print "No propagated sequence found!"
		print hits
		o.write("{}\t{}\n".format(f.split(".")[0],"\t".join([str(hits.count(s)) for s in strains])))
	return

def main():
	args = parse_args()
	outdir = os.path.abspath(args.outdir)
	genomedb = os.path.abspath(args.genomedb)
	new_strains = [line.rstrip() for line in open(os.path.abspath(args.new_strains),'r')]

	setupdir(outdir)
	check_strains(new_strains,genomedb,outdir)
	make_diamond_databases(new_strains,outdir)
	run_diamond(new_strains,outdir)
	genes = get_genes(new_strains,outdir)
	parse_diamond(genes,outdir)
	run_inparanoid(new_strains, outdir)
	clean_up(outdir)
	group_members = parse_inparanoid(outdir,new_strains)
	extract_fastas(outdir,genes,group_members)
	dump_matrix(outdir)

if __name__ == '__main__':
	main()
