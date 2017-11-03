import os
import errno
import ftplib
import ijson
import subprocess
import datetime
import shutil
import sys
from urllib import urlopen
from Bio import SeqIO

def check_db(outdir):
	assemblies = []
	if "genome_metadata.txt" not in os.listdir(outdir):
		print os.path.basename(outdir), "is empty..."
		pass
	else:
		for line in open(os.path.join(outdir,"genome_metadata.txt"),'r'):
			if line.startswith("assembly_id"):
				continue
			else:
				assemblies.append(line.rstrip().split("\t")[2])
	return assemblies

def setupdirs(outdir):
	print "Setting up", os.path.basename(outdir)
	try:
		os.makedirs(outdir)
	except OSError as exc:
		if exc.errno == errno.EEXIST:
			print "\tDatabase folder exists:", outdir

	for f in ["pep","dna","gbk"]:
		try:
			os.makedirs(os.path.join(os.path.join(outdir,f)))
		except OSError:
			print "\tSubfolder exists:", os.path.join(os.path.basename(outdir),f)

	return

def download_Ensembl_files(outdir, names=False, maxgen=100, taxids=False, complete=True):
	assemblies = check_db(outdir)
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
	genomes = {}

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

		if names:
			for n in names.split(","):
				if n in fields:
					if js["species"] not in assemblies:
						if complete:
							if js["assembly_level"] == "chromosome":
								genomes = get_data(genomes,js)
								found += 1
						else:
							genomes = get_data(genomes,js)
							found += 1
		elif taxids:
			if str(js['taxonomy_id']) in taxids:
				if js['species'] not in assemblies:
					if complete:
						if js["assembly_level"] == "chromosome":
							genomes = get_data(genomes,js)
							found += 1
					else:
						genomes = get_data(genomes,js)
						found += 1
		else:
			if js['species'] not in assemblies:
				if complete:
					if js["assembly_level"] == "chromosome":
						genomes = get_data(genomes,js)
						found += 1
				else:
					genomes = get_data(genomes,js)
					found += 1

		o.write("\t".join([str(x) for x in thisline])+"\n")

		if count % 10000 == 0:
			print "\t", count, "JSON records parsed."
		if maxgen is not None:
			if found == maxgen:
				print "{} new genomes to download found...exiting JSON parser...".format(found)
				break

	print "\t",count, "total JSON records parsed..."
	print len(assemblies), "found in", os.path.basename(outdir)+"."
	print len(genomes), "remaining to download."
	o.close()
	ens.close()
	for g in genomes:
		genomes[g]["version"] = j
	Ensembl_ftp(genomes,outdir)
	return


def get_data(genomes,js):
	if js["assembly_id"] not in genomes:
		genomes[js["assembly_id"]] = {}
		for feat in ['assembly_level','base_count','name', 'strain', 'dbname','species','taxonomy_id']:
			genomes[js["assembly_id"]][feat] = js[feat]
		genomes[js["assembly_id"]]['contigs'] = len(js["sequences"])
		genomes[js["assembly_id"]]['ngenes'] = js["annotations"]["nProteinCoding"]
	return genomes

def Ensembl_ftp(fg, outdir):
	ens = ftplib.FTP('ftp.ensemblgenomes.org')
	ens.login()
	print "Downloading {} genome files...".format(str(len(fg.keys())))
	count = 0
	files = os.listdir(outdir)
	if "genome_metadata.txt" not in files:
		print "Initializing metadata file..."
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
			vals = [f]
			for key in ["base_count", "species", "taxonomy_id", "contigs","ngenes"]:
				vals.append(fg[f][key])
			vals.append("ensembl-"+fg[f]["version"])
			vals.append(datetime.datetime.now())
			vals.append(datetime.datetime.now())
			o.write("\t".join([str(v) for v in vals])+"\n")
			count += 1
			if count % 100 == 0:
				print "\t", count, "files downloaded."
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


def check_unique(species_id,outdir):
	strains = []
	for line in open(os.path.join(outdir,"genome_metadata.txt"),'r'):
		if line.startswith("assembly_id"):
			continue
		else:
			strains.append(line.rstrip().split("\t")[2])

	if species_id not in strains:
		print "Species ID is unique! Moving on..."
		return True
	else:
		print "Species ID is not unique. Select a new ID."
		return False

def add_Prokka_genome(outdir,prokka,species_id,tax_id="2"):
	if check_unique(species_id,outdir):
		stats = {}
		print "Copying files..."
		for f in os.listdir(prokka):
			if f.endswith(".faa"):
				stats['ngenes'] = 0
				for seq in SeqIO.parse(open(os.path.join(prokka,f),'r'),'fasta'):
					stats['ngenes'] += 1
				shutil.copy(os.path.join(prokka,f),os.path.join(outdir,"pep","{}.pep.fa".format(species_id)))
			elif f.endswith(".fna"):
				basecount = 0
				contigcount = 0
				for seq in SeqIO.parse(open(os.path.join(prokka,f),'r'),'fasta'):
					basecount += len(str(seq.seq))
					contigcount += 1
				stats['basecount'] = basecount
				stats['contigcount'] = contigcount
				shutil.copy(os.path.join(prokka,f),os.path.join(outdir,"dna","{}.dna.fa".format(species_id)))
			elif f.endswith(".gbk"):
				shutil.copy(os.path.join(prokka,f),os.path.join(outdir,"gbk","{}.gbk".format(species_id)))
			else:
				pass
		o = open(os.path.join(outdir,"genome_metadata.txt"),'a')
		vals = [species_id+"_v1",stats['basecount'],species_id,tax_id,stats['contigcount'],stats['ngenes'],"prokka_in_house", datetime.datetime.now(),datetime.datetime.now()]
		o.write("\t".join([str(v) for v in vals])+"\n")
		o.close()
	else:
		pass
	return
