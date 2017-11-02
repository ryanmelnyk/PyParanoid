#Ryan A. Melnyk
#schmelnyk@gmail.com
#UBC Microbiology - Haney Lab

import argparse
import os
from Bio import Entrez
import datetime


def parse_args():
	parser = argparse.ArgumentParser(description='''
Populate the taxonomy table with taxonomy information extracted from NCBI Taxonomy.
	''')
	parser.add_argument('outdir', type=str,help='directory where genomes have been downloaded.')
	parser.add_argument('email',type=str,help="email for accessing Entrez-NCBI")
	parser.add_argument('--max',type=int,help='maximum number of taxonomy queries to perform before exiting')
	return parser.parse_args()

def fetch_tax(taxdata,outdir,email):
	Entrez.email = email
	if os.path.exists(os.path.join(outdir,"tax_info.txt")):
		o = open(os.path.join(outdir,"tax_info.txt"),'a')
	else:
		o = open(os.path.join(outdir,"tax_info.txt"),'w')
		o.write("species_id\ttaxonomy_id\tsuperkingdom\tphylum\tclass\torder\tfamily\tgenus\tspecies\tdate_added\tdate_modified\n")
	count = len(taxdata)
	print "Extracting", str(count), "taxonomy records from Entrez-NCBI..."
	for t in taxdata:
		values = [t,taxdata[t],None,None,None,None,None,None,None,datetime.datetime.now(),datetime.datetime.now()]
		records = Entrez.parse(Entrez.efetch(db="taxonomy",id=str(taxdata[t]),retmode="xml"))
		for r in records:
			for l in r["LineageEx"]:
				if l["Rank"] == "superkingdom":
					values[2] = l["ScientificName"]
				elif l["Rank"] == "phylum":
					values[3] = l["ScientificName"]
				elif l["Rank"] == "class":
					values[4] = l["ScientificName"]
				elif l["Rank"] == "order":
					values[5] = l["ScientificName"]
				elif l["Rank"] == "family":
					values[6] = l["ScientificName"]
				elif l["Rank"] == "genus":
					values[7] = l["ScientificName"]
				elif l["Rank"] == "species":
					values[8] = l["ScientificName"]
				else:
					pass
		count -= 1
		if count % 100 == 0:
			print count, "records remaining..."

		o.write("\t".join([str(x) for x in values])+"\n")
	print "Done!"

	return

def get_genomedb_data(outdir,max_queries):
	taxdata = {}
	existing = []
	if os.path.exists(os.path.join(outdir,"tax_info.txt")):
		for line in open(os.path.join(outdir,"tax_info.txt"),'r'):
			if line.startswith("species_id"):
				continue
			else:
				existing.append(line.rstrip().split("\t")[0])
	for line in open(os.path.join(outdir,"genome_metadata.txt"),'r'):
		if line.startswith("assembly_id"):
			continue
		else:
			vals = line.rstrip().split("\t")
			if max_queries:
				if len(taxdata.keys()) == max_queries:
					break
			if vals[2] not in existing:
				taxdata[vals[2]] = vals[3]
	return taxdata

def main():
	args = parse_args()
	outdir = os.path.abspath(args.outdir)
	if args.max:
		max_queries = args.max
	else:
		max_queries = None

	taxdata = get_genomedb_data(outdir,max_queries)
	fetch_tax(taxdata,outdir,email)

if __name__ == '__main__':
	main()
