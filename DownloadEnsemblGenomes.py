#Ryan A. Melnyk
#schmelnyk@gmail.com
#UBC Microbiology - Haney Lab

# This script is used to download genomes from Ensembl Bacteria (http://bacteria.ensembl.org/)

import ijson
import ftplib
import datetime
import argparse, os, errno, sys, subprocess
from urllib import urlopen

def parse_args():
	parser = argparse.ArgumentParser(description='''
A script for accessing the current release of Ensembl Bacteria and downloading
complete genomes.  Will download nucleotide, amino acid, and CDS information and store
some metadata in a text file.
	''')
	parser.add_argument('genomedb', type=str,help='directory to download genomes')
	parser.add_argument('--names', type=str, help='comma-separated list of keywords to find in name of strain to download draft genomes. i.e. "pseudomonas,salmonella,syringae,K12". Only complete/exact matches will be downloaded and spelling counts!')
	parser.add_argument('--max', type=int, help="maximum number of genomes to download")
	return parser.parse_args()

def parse_json(outdir, assemblies, names, maxgen):
	ens = ftplib.FTP('ftp.ensemblgenomes.org')
	ens.login()
	ens.cwd('pub/bacteria/current')

	for j in ens.pwd().split("/"):
		if j.startswith("release"):
			print "Current release of EnsemblBacteria:", j
			o = open(os.path.join(outdir,"{}.txt".format(j)),'w')
			break

	fields = ["assembly_id",'assembly_level','base_count','name', 'strain', 'dbname','species','taxonomy_id','contigs','protein_coding_genes']
	o.write("\t".join(fields)+"\n")

	items = ijson.items(urlopen("ftp://ftp.ensemblgenomes.org/{}/species_metadata_EnsemblBacteria.json".format(ens.pwd())),'item')
	count = 0
	finished_genomes = {}

	found = 0
	for js in items:
		count += 1

		thisline = []
		thisline.append(js["assembly_id"])
		for feat in ['assembly_level','base_count','name', 'strain', 'dbname','species','taxonomy_id']:
			thisline.append(js[feat])
		thisline.append(str(len(js["sequences"])))
		thisline.append(js["annotations"]["nProteinCoding"])

		fields = js['species'].split("_")
		for n in names:
			if n in fields:
				if js["species"] not in assemblies:
					# print js["species"]
					finished_genomes[js["assembly_id"]] = {}
					for feat in ['assembly_level','base_count','name', 'strain', 'dbname','species','taxonomy_id']:
						finished_genomes[js["assembly_id"]][feat] = js[feat]
					finished_genomes[js["assembly_id"]]['contigs'] = len(js["sequences"])
					finished_genomes[js["assembly_id"]]['ngenes'] = js["annotations"]["nProteinCoding"]
					found += 1

		o.write("\t".join([str(x) for x in thisline])+"\n")

		if count % 1000 == 0:
			print count, "JSON records parsed."
		if maxgen is not None:
			if found == maxgen:
				print "{} new genomes to download found...exiting JSON parser...".format(found)
				break

	print count, "total JSON records parsed."
	print len(assemblies), "found in", os.path.basename(outdir)+"."
	print len(finished_genomes), "remaining to download."
	o.close()
	ens.close()
	return finished_genomes,j

def setupdirs(outdir):
	try:
		os.makedirs(outdir)
	except OSError as exc:
		if exc.errno == errno.EEXIST:
			print "Database folder exists:", outdir

	for f in ["pep","dna", "gff3"]:
		try:
			os.makedirs(os.path.join(os.path.join(outdir,f)))
		except OSError:
			print "Subfolder exists:", outdir

	print outdir
	return

def get_files(fg, outdir, EV):
	ens = ftplib.FTP('ftp.ensemblgenomes.org')
	ens.login()
	print "Downloading genome files..."
	count = 0
	files = os.listdir(outdir)
	if "genome_metadata.txt" not in files:
		o = open(os.path.join(outdir,"genome_metadata.txt"),'w')
		fields = ["assembly_id", "base_count", "species", "taxonomy_id", "contigs", "protein_coding_genes", "source", "date_added", "date_modified"]
		o.write("\t".join(fields)+"\n")
	else:
		o = open(os.path.join(outdir,"genome_metadata.txt"),'a')
	for f in fg:
		try:
			ens.cwd("/pub/bacteria/current/fasta/{}/{}/dna".format("_".join(fg[f]["dbname"].split("_")[0:3]),fg[f]["species"]))
			for filepath in ens.nlst():
				if filepath.endswith(".dna.toplevel.fa.gz"):
					download_and_unzip(ens,filepath,os.path.join(outdir,"dna",fg[f]["species"]+".dna.fa.gz"))
			ens.cwd("../pep")
			for filepath in ens.nlst():
				if filepath.endswith(".pep.all.fa.gz"):
					download_and_unzip(ens,filepath,os.path.join(outdir,"pep",fg[f]["species"]+".pep.fa.gz"))
			ens.cwd("/pub/bacteria/current/gff3/{}/{}".format("_".join(fg[f]["dbname"].split("_")[0:3]),fg[f]["species"]))
			for filepath in ens.nlst():
				fields = filepath.split(".")
				if ".".join(fields[3:]) == ("gff3.gz"):
					download_and_unzip(ens,filepath,os.path.join(outdir,"gff3",fg[f]["species"]+".gff3.gz"))
			vals = [f]
			for key in ["base_count", "species", "taxonomy_id", "contigs","ngenes"]:
				vals.append(fg[f][key])
			vals.append("ensembl-"+EV)
			vals.append(datetime.datetime.now())
			vals.append(datetime.datetime.now())
			o.write("\t".join([str(v) for v in vals])+"\n")
			count += 1
			print fg[f]["species"], "processed.", count, "files downloaded."
		except ftplib.error_perm:
			print "Issue with data for {}...skipping.".format(f)

	o.close()
	ens.close()
	return

def download_and_unzip(ftp,f,outfile):
	o = open(outfile,'wb')
	ftp.retrbinary("RETR " + f, o.write)
	o.close()
	cmds = ["gunzip",outfile]
	proc = subprocess.Popen(cmds)
	proc.wait()
	return

def check_db(outdir):
	assemblies = []
	if "genome_metadata.txt" not in os.listdir(outdir):
		pass
	else:
		for line in open(os.path.join(outdir,"genome_metadata.txt"),'r'):
			if line.startswith("assembly_id"):
				continue
			else:
				assemblies.append(line.rstrip().split("\t")[3])
	return assemblies

def main():
	args = parse_args()

	outdir = os.path.abspath(args.genomedb)
	if args.names == None:
		names = []
	else:
		names = args.names.split(",")
	for n in names:
		print n

	setupdirs(outdir)

	assemblies = check_db(outdir)

	finished_genomes, ENSEMBL_VERSION = parse_json(outdir,assemblies,names,args.max)
	get_files(finished_genomes, outdir, ENSEMBL_VERSION)

if __name__ == '__main__':
	main()
