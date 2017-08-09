#Ryan A. Melnyk
#schmelnyk@gmail.com
#UBC Microbiology - Haney Lab

import os
import sys
import argparse
import subprocess
from Bio import SeqIO
import shutil
from datetime import datetime
import ncbi_genome_download as ngd

def parse_args():
	parser = argparse.ArgumentParser(description='''
Using Kai Blin's ncbi-genome-download to add genomes to a local genomedb.  Downloads
from RefSeq, as Genbank is more or less redundant with EnsemblGenomes.
	''')

	parser.add_argument('outdir', type=str,help='directory to download genomes')
	parser.add_argument('--name',type=str,help="strain names for ncbi-genome-download - use quotes")
	parser.add_argument('--taxid',type=str,help="ncbi ID for ncbi-genome-download")
	parser.add_argument('--cpus',type=int,help="number of cpus to use.")
	return parser.parse_args()

def check_db(outdir):
	assemblies,species_tags = [],[]
	if "genome_metadata.txt" not in os.listdir(outdir):
		pass
	else:
		for line in open(os.path.join(outdir,"genome_metadata.txt"),'r'):
			if line.startswith("assembly_id"):
				continue
			else:
				vals = line.rstrip().split("\t")
				assemblies.append(vals[0])
				species_tags.append(vals[2])


	return assemblies, species_tags

def download_files(name,taxid,outdir,cpus):
	print "Downloading files..."
	files = ["fasta","protein-fasta","gff","assembly-stats"]
	if name:
		for f in files:
			print "\tworking on {} files...".format(f)
			ngd.download(group="bacteria",genus=name,file_format=f,section="refseq",output=outdir,parallel=cpus)
	elif taxid:
		for f in files:
			print "\tworking on {} files...".format(f)
			ngd.download(group="bacteria",taxid=taxid,file_format=f,section="refseq",output=outdir,parallel=cpus)
	return

def copy_files(outdir,assemblies,species_tags):
	metadata = open(os.path.join(outdir,"genome_metadata.txt"),'a')
	if not os.path.exists(os.path.join(outdir,"refseq")):
		print "No files downloaded!"
		print "Check genus/species name or taxid"
		sys.exit()
	for f in os.listdir(os.path.join(outdir,"refseq","bacteria")):
		if len(os.listdir(os.path.join(outdir,"refseq","bacteria",f))) != 5:
			print "Files may be missing. Skipping..."
			add = False
			continue
		else:
			add = True
		assembly_id, orgname, tax_id = "","",""
		for g in os.listdir(os.path.join(outdir,"refseq","bacteria",f)):
			if g.endswith(".gz"):
				cmds = "gunzip {}".format(os.path.join(outdir,"refseq","bacteria",f,g))
				proc = subprocess.Popen(cmds.split())
				proc.wait()
			if g.endswith("_assembly_stats.txt"):
				for line in open(os.path.join(outdir,"refseq","bacteria",f,g)):
					if line.startswith("# Organism name:"):
						orgname = line.split(":")[1].split("(")[0].strip().replace(" ","_").replace("-","_").replace("/","_").replace(".","_").lower()
						print "Processing {}...".format(orgname)
					if line.startswith("# GenBank assembly accession:"):
						assembly_id = line.split(":")[1].strip()
						if assembly_id in assemblies:
							add = False
							print "\t", assembly_id, "already in", outdir
						else:
							if orgname in species_tags:
								print "\t", orgname, "isn't a unique name."
								orgname += assembly_id.replace("GCA","_gca").split(".")[0]
								print "\t", "renamed to:", orgname
					if line.startswith("# Taxid:"):
						tax_id = line.split(":")[1].strip()
		if add:
			for g in os.listdir(os.path.join(outdir,"refseq","bacteria",f)):
				if g.endswith(".fna"):
					total_length = 0
					contig_count = 0
					o = open(os.path.join(outdir,"dna",orgname+".dna.fa"),'w')
					for seq in SeqIO.parse(open(os.path.join(outdir,"refseq","bacteria",f,g),'r'),'fasta'):
						total_length += len(str(seq.seq))
						SeqIO.write(seq,o,'fasta')
						contig_count += 1
					o.close()
				if g.endswith(".faa"):
					o = open(os.path.join(outdir,"pep",orgname+".pep.fa"),'w')
					gene_count = 0
					for seq in SeqIO.parse(open(os.path.join(outdir,"refseq","bacteria",f,g),'r'),'fasta'):
						gene_count += 1
						SeqIO.write(seq,o,'fasta')
					o.close()
				if g.endswith(".gff"):
					shutil.copy(os.path.join(outdir,"refseq","bacteria",f,g),os.path.join(outdir,"gff3",orgname+".gff3"))

			meta_vals = [assembly_id,total_length,orgname,tax_id,contig_count,gene_count,"NCBI-RefSeq",datetime.now(),datetime.now()]
			metadata.write("\t".join([str(x) for x in meta_vals])+"\n")
	metadata.close()
	return

def main():
	args = parse_args()
	outdir = os.path.abspath(args.outdir)
	if args.name and args.taxid:
		print "Can't specify a name and a taxid...exiting."
		sys.exit()
	if args.name:
		name = args.name
	else:
		name = None

	if args.taxid:
		taxid = args.taxid
	else:
		taxid = None

	if args.cpus:
		cpus = args.cpus
	else:
		cpus = None

	assemblies,species_tags = check_db(outdir)
	if os.path.exists(os.path.join(outdir,"refseq")):
		print "Incomplete NCBI download in",outdir,"detected."
		print "Clean up NCBI files before proceeding."
		sys.exit()
	download_files(name,taxid,outdir,cpus)
	copy_files(outdir,assemblies,species_tags)
	shutil.rmtree(os.path.join(outdir,"refseq"))





if __name__ == '__main__':
	main()
