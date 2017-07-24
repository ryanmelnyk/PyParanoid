#Ryan A. Melnyk
#schmelnyk@gmail.com
#UBC Microbiology - Haney Lab

import argparse, os, shutil
import sys, datetime
from Bio import SeqIO, Entrez

def parse_args():
	parser = argparse.ArgumentParser(description='''
This script will add a local genome annotated with Prokka to the genomedb.""
	''')
	parser.add_argument('prokka', type=str,help='relative path to prokka directory')
	parser.add_argument('genomedb', type=str,help='directory to download genomes')
	parser.add_argument('species_id', type=str,help='species_id for strain identification (Must be unique in genomedb!)')
	parser.add_argument('tax_id', type=str, help='taxonomy id for NCBI lookup. (Pseudomonas = "286")')
	return parser.parse_args()

def check_unique(species_id,outdir):
	strains = []
	for line in open(os.path.join(outdir,"genome_metadata.txt"),'r'):
		if line.startswith("assembly_id"):
			continue
		else:
			strains.append(line.rstrip().split("\t")[3])

	if species_id not in strains:
		print "Species ID is unique! Moving on..."
	else:
		print "Species ID is not unique. Select a new ID."
		print "Exiting script..."
		sys.exit()
	return


def copy_files(outdir,prokka,species_id,tax_id):
	CDS = []
	RNA = []
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
		elif f.endswith(".gff"):
			shutil.copy(os.path.join(prokka,f),os.path.join(outdir,"gff3","{}.gff3".format(species_id)))
		else:
			pass
	o = open(os.path.join(outdir,"genome_metadata.txt"),'a')
	vals = [species_id+"_v1",stats['basecount'],species_id,tax_id,stats['contigcount'],stats['ngenes'],"prokka_in_house", datetime.datetime.now(),datetime.datetime.now()]
	o.write("\t".join([str(v) for v in vals])+"\n")
	o.close()
	return

def main():
	args = parse_args()
	prokka = os.path.abspath(args.prokka)
	outdir = os.path.abspath(args.genomedb)
	species_id = args.species_id
	tax_id = args.tax_id

	check_unique(species_id,outdir)
	stats = copy_files(outdir,prokka,species_id,tax_id)


if __name__ == '__main__':
	main()
