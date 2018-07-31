import os
import errno
import ftplib
import ijson
import subprocess
import datetime
import shutil
import sys
from urllib import urlopen
from Bio import SeqIO, Entrez
import ncbi_genome_download as ngd

def check_db(outdir):
	assemblies,species_tags = [],[]
	if os.path.isdir(outdir):
		pass
	else:
		print os.path.abspath(outdir), "doesn't exist! Run gdb.setupdirs() first."
		sys.exit()

	for f in ["pep","dna","gbk"]:
		if os.path.isdir(os.path.join(outdir,f)):
			pass
		else:
			print os.path.abspath(os.path.join(outdir,f)), "doesn't exist! Run gdb.setupdirs() first."
			sys.exit()

	if "genome_metadata.txt" not in os.listdir(outdir):
		pass
	else:
		o = open(os.path.join(outdir,"genome_metadata.txt"),'r')
		for line in o:
			if line.startswith("assembly_id"):
				continue
			else:
				vals = line.rstrip().split("\t")
				assemblies.append(vals[0])
				species_tags.append(vals[2])
		o.close()

	return assemblies, species_tags

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

def download_Ensembl_files(outdir, names=False, maxgen=10, taxids=False, complete=True):
	assemblies,species_tags = check_db(outdir)
	ens = ftplib.FTP('ftp.ensemblgenomes.org')
	ens.login()
	ens.cwd('pub/bacteria/current')

	for j in ens.pwd().split("/"):
		if j.startswith("release"):
			print "Current release of EnsemblBacteria:", j
			break

	if j+".txt" not in os.listdir(outdir):
		o = open(os.path.join(outdir,"{}.txt".format(j)),'w')
		fields = ["assembly_id",'assembly_level','base_count','name', 'strain',\
				  'dbname','species','taxonomy_id','contigs','protein_coding_genes']
		o.write("\t".join(fields)+"\n")

		items = ijson.items(urlopen("ftp://ftp.ensemblgenomes.org/{}/species_metadata_EnsemblBacteria.json"
									.format(ens.pwd())),'item')
		count = 0
		for js in items:
			count += 1

			thisline = []
			for feat in ["assembly_accession",'assembly_level','base_count']:
				thisline.append(js['assembly'][feat])
			for feat in ['display_name', 'strain']:
				thisline.append(js['organism'][feat])
			thisline.append(js["core"]["dbname"])
			for feat in ['name', 'taxonomy_id']:
				thisline.append(js['organism'][feat])
			thisline.append(str(len(js["assembly"]["sequences"])))
			thisline.append(js["annotations"]["nProteinCoding"])

			o.write("\t".join([str(x) for x in thisline])+"\n")

			if count % 10000 == 0:
				print "\t", count, "JSON records parsed."
		print "\t",count, "total JSON records parsed..."
		o.close()
	else:
		print j+".txt", "found in", outdir+"...parsing..."
	ens.close()

	genomes = {}
	found = 0
	for line in open(os.path.join(outdir,j+".txt"),'r'):
		if line.startswith("assembly_id"):
			continue
		else:
			vals = line.rstrip().split("\t")

		fields = vals[6].split("_")

		if names:
			for n in names.split(","):
				if n in fields:
					if (vals[6] not in species_tags) and (vals[0] not in assemblies):
						if complete:
							if vals[1] == "chromosome":
								genomes = get_data(genomes,vals)
						else:
							genomes = get_data(genomes,vals)
		elif taxids:
			if str(vals[7]) in taxids:
				if (vals[6] not in species_tags) and (vals[0] not in assemblies):
					if complete:
						if vals[1] == "chromosome":
							genomes = get_data(genomes,vals)
					else:
						genomes = get_data(genomes,vals)
		elif not (names or taxids):
			if (vals[6] not in species_tags) and (vals[0] not in assemblies):
				if complete:
					if vals[1] == "chromosome":
						genomes = get_data(genomes,vals)
				else:
					genomes = get_data(genomes,vals)
		else:
			pass

		if maxgen is not None:
			if len(genomes.keys()) == maxgen:
				print "\tMaximum of {} genomes reached..."\
					.format(str(len(genomes.keys())))
				break

	print len(assemblies), "found in", os.path.basename(outdir)+"."
	print len(genomes), "genomes in queue to download..."

	for g in genomes:
		genomes[g]["version"] = j
	Ensembl_ftp(genomes,outdir)
	return

def get_data(genomes,vals):
	if vals[0] not in genomes:
		genomes[vals[0]] = {}
		genomes[vals[0]]['assembly_level'] = vals[1]
		genomes[vals[0]]['base_count'] = vals[2]
		genomes[vals[0]]['display_name'] = vals[3]
		genomes[vals[0]]['strain'] = vals[4]
		genomes[vals[0]]['dbname'] = vals[5]
		genomes[vals[0]]['name'] = vals[6]
		genomes[vals[0]]['taxonomy_id'] = vals[7]
		genomes[vals[0]]['contigs'] = vals[8]
		genomes[vals[0]]['ngenes'] = vals[9]
	return genomes

def Ensembl_ftp(fg, outdir):
	ens = ftplib.FTP('ftp.ensemblgenomes.org')
	ens.login()
	print "Downloading {} genome files...".format(str(len(fg.keys())))
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
			ens.cwd("/pub/bacteria/current/fasta/{}/{}/pep".format("_".join(fg[f]["dbname"].split("_")[0:3]),fg[f]["name"]))
			for filepath in ens.nlst():
				if filepath.endswith(".pep.all.fa.gz"):
					download_and_unzip(ens,filepath,os.path.join(outdir,"pep",fg[f]["name"]+".pep.fa.gz"))

			vals = [f]
			for key in ["base_count", "name", "taxonomy_id", "contigs","ngenes"]:
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

def add_Prokka_genome(outdir,prokka,species_id,taxid="2"):
	if check_unique(species_id,outdir):
		stats = {}
		print "Copying files for", species_id
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
				shutil.copy(os.path.join(prokka,f),os.path.join(outdir,"dna","{}.fna".format(species_id)))
			elif f.endswith(".gbk"):
				shutil.copy(os.path.join(prokka,f),os.path.join(outdir,"gbk","{}.gbk".format(species_id)))
			else:
				pass
		o = open(os.path.join(outdir,"genome_metadata.txt"),'a')
		vals = [species_id+"_v1",stats['basecount'],species_id,taxid,stats['contigcount'],stats['ngenes'],"prokka_in_house", datetime.datetime.now(),datetime.datetime.now()]
		o.write("\t".join([str(v) for v in vals])+"\n")
		o.close()
	else:
		pass
	return

def add_img_genome(outdir,img,species_id,taxid="2"):
	if check_unique(species_id,outdir):
		stats = {}
		print "Copying files for", species_id
		for f in os.listdir(img):
			vals = f.split(".")
			if f.endswith(".faa"):
				stats['ngenes'] = 0
				for seq in SeqIO.parse(open(os.path.join(img,f),'r'),'fasta'):
					stats['ngenes'] += 1
				shutil.copy(os.path.join(img,f),os.path.join(outdir,"pep","{}.pep.fa".format(species_id)))
			elif vals[1] == "fna":
				basecount = 0
				contigcount = 0
				for seq in SeqIO.parse(open(os.path.join(img,f),'r'),'fasta'):
					basecount += len(str(seq.seq))
					contigcount += 1
				stats['basecount'] = basecount
				stats['contigcount'] = contigcount
				shutil.copy(os.path.join(img,f),os.path.join(outdir,"dna","{}.fna".format(species_id)))
			else:
				pass
		o = open(os.path.join(outdir,"genome_metadata.txt"),'a')
		vals = [species_id,stats['basecount'],species_id,taxid,stats['contigcount'],stats['ngenes'],"img", datetime.datetime.now(),datetime.datetime.now()]
		o.write("\t".join([str(v) for v in vals])+"\n")
		o.close()
	else:
		pass
	return

def download_Refseq_files(outdir,cpus=1,names=False,taxids=False):
	assemblies,species_tags = check_db(outdir)

	files = ["fasta","protein-fasta","assembly-stats"]
	if not (names or taxids):
		print "Must specify a name or a taxid."
	elif os.path.exists(os.path.join(outdir,"refseq")):
		print "Refseq download already exists at", os.path.join(outdir,"refseq")
		print "Delete before proceeding."
	else:
		if names:
			for name in names.split(","):
				print "Downloading files for {}...".format(name)
				for f in files:
					print "\tworking on {} files...".format(f)
					if cpus == 1:
						ngd.download(group="bacteria",genus=name,file_format=f,section="refseq",output=outdir)
					else:
						ngd.download(group="bacteria",genus=name,file_format=f,section="refseq",output=outdir,parallel=cpus)
		if taxids:
			for taxid in taxids.split(","):
				print "Downloading files for {}...".format(str(taxid))
				for f in files:
					print "\tworking on {} files...".format(f)
					if cpus == 1:
						ngd.download(group="bacteria",taxid=taxid,file_format=f,section="refseq",output=outdir)
					else:
						ngd.download(group="bacteria",taxid=taxid,file_format=f,section="refseq",output=outdir,parallel=cpus)
		process_Refseq(outdir,assemblies,species_tags)
		if os.path.exists(os.path.join(outdir,"refseq")):
			shutil.rmtree(os.path.join(outdir,"refseq"))
	return

def process_Refseq(outdir,assemblies,species_tags):
	metadata = open(os.path.join(outdir,"genome_metadata.txt"),'a')
	if not os.path.exists(os.path.join(outdir,"refseq")):
		print "No files downloaded!"
		print "Check genus/species name or taxid"
	else:
		print len(os.listdir(os.path.join(outdir,"refseq","bacteria"))), "genomes to process."
		count = 0
		added = 0
		for f in os.listdir(os.path.join(outdir,"refseq","bacteria")):
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
						if line.startswith("# GenBank assembly accession:"):
							assembly_id = line.split(":")[1].strip()
							if assembly_id in assemblies:
								add = False
							else:
								if orgname in species_tags:
									orgname += assembly_id.replace("GCA","_gca").split(".")[0]
								add = True
						if line.startswith("# Taxid:"):
							tax_id = line.split(":")[1].strip()
			if add:
				for g in os.listdir(os.path.join(outdir,"refseq","bacteria",f)):
					if g.endswith(".fna"):
						total_length = 0
						contig_count = 0

						for seq in SeqIO.parse(open(os.path.join(outdir,"refseq","bacteria",f,g),'r'),'fasta'):
							total_length += len(str(seq.seq))
							contig_count += 1

					if g.endswith(".faa"):
						o = open(os.path.join(outdir,"pep",orgname+".pep.fa"),'w')
						gene_count = 0
						for seq in SeqIO.parse(open(os.path.join(outdir,"refseq","bacteria",f,g),'r'),'fasta'):
							gene_count += 1
							SeqIO.write(seq,o,'fasta')
						o.close()
				meta_vals = [assembly_id,total_length,orgname,tax_id,contig_count,gene_count,"NCBI-RefSeq",datetime.datetime.now(),datetime.datetime.now()]
				metadata.write("\t".join([str(x) for x in meta_vals])+"\n")
				species_tags.append(orgname)
				added += 1
			count += 1
			if count % 100 == 0:
				print count, "processed..."
		metadata.close()
		print count, "genomes processed and", added, "added to {}.".format(os.path.basename(outdir))
	return

def get_taxonomy(outdir,email="yourname@email.com",max_queries=None):
	taxdata = get_taxdata_from_genomedb(outdir,max_queries)
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

def get_taxdata_from_genomedb(outdir, max_queries):
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

def download_genbank_files(strains,genomedb):
	prokka = []
	refseq = []
	ensembl = {}
	for line in open(os.path.join(genomedb,"genome_metadata.txt")):
		vals = line.rstrip().split("\t")
		if vals[2] in strains:
			if vals[6].startswith("ensembl"):
				rel = "-".join(vals[6].split("-")[1:])
				if rel not in ensembl:
					ensembl[rel] = [vals[2]]
				else:
					ensembl[rel].append(vals[2])
			elif vals[6].startswith("NCBI"):
				refseq.append((vals[0],vals[2]))
			elif vals[6].startswith("prokka"):
				prokka.append((vals[0],vals[2]))
			else:
				pass


	p_count = 0
	for p in prokka:
		if os.path.exists(os.path.join(genomedb,"gbk",p[1]+".gbk")):
			pass
			p_count += 1
		else:
			print "Genbank from prokka assembly",p[1], "not found..."

	ens = ftplib.FTP('ftp.ensemblgenomes.org')
	ens.login()
	e_count = 0
	for e in ensembl:
		for line in open(os.path.join(genomedb,"{}.txt".format(e)),'r'):
			vals = line.rstrip().split("\t")
			if vals[6] in ensembl[e]:
				if os.path.exists(os.path.join(genomedb,"gbk",vals[6]+".gbk")):
					e_count += 1
				else:
					wd = '/pub/{}/bacteria/genbank/{}_collection/{}'.format(e,"_".join(vals[5].split("_")[0:2]),vals[6])
					ens.cwd(wd)

					for filepath in ens.nlst():
						if filepath.endswith(".dat.gz"):
							o = open(os.path.join(genomedb,"gbk",vals[6]+'.gbk.gz'),'wb')
							ens.retrbinary("RETR " + filepath, o.write)
							o.close()
							cmds = ["gunzip",os.path.join(genomedb,"gbk",vals[6]+'.gbk.gz')]
							proc = subprocess.Popen(cmds)
							proc.wait()
					e_count += 1
			else:
				pass
	ens.close()

	r_count = 0
	for r in refseq:
		if os.path.exists(os.path.join(genomedb,"gbk",r[1]+".gbk")):
			r_count += 1
		else:
			assembly = r[0].split("_")[1]
			species = r[1]
			ens = ftplib.FTP('ftp.ncbi.nlm.nih.gov')
			ens.login()
			ens.cwd("/genomes/all/GCF/{}/{}/{}".format(assembly[0:3],assembly[3:6], assembly[6:9]))
			for f in ens.nlst():
				if f.split("_")[1] == assembly:
					ens.cwd(f)
					for g in ens.nlst():
						if g.endswith("_genomic.gbff.gz"):
							filepath = os.path.join(genomedb,"gbk","{}.gbk.gz".format(species))
							gbk = open(filepath,'wb')
							ens.retrbinary("RETR " + g, gbk.write)
							gbk.close()
							cmds = ["gunzip",os.path.join(genomedb,"gbk","{}.gbk.gz".format(species))]
							proc = subprocess.Popen(cmds)
							proc.wait()
							r_count += 1
			ens.close()

	e_tot = 0
	for e in ensembl:
		e_tot += len(ensembl[e])
	print p_count,"of", len(prokka), "prokka genbank files available."
	print e_count,"of", e_tot, "Ensembl genbank files available."
	print r_count, "of", len(refseq), "refseq genbank files available."
	return

def download_dna_files(strains,genomedb):
	prokka = []
	refseq = []
	ensembl = {}
	for line in open(os.path.join(genomedb,"genome_metadata.txt")):
		vals = line.rstrip().split("\t")
		if vals[2] in strains:
			if vals[6].startswith("ensembl"):
				rel = "-".join(vals[6].split("-")[1:])
				if rel not in ensembl:
					ensembl[rel] = [vals[2]]
				else:
					ensembl[rel].append(vals[2])
			elif vals[6].startswith("NCBI"):
				refseq.append((vals[0],vals[2]))
			elif vals[6].startswith("prokka"):
				prokka.append((vals[0],vals[2]))
			else:
				pass

	p_count = 0
	for p in prokka:
		if os.path.exists(os.path.join(genomedb,"dna",p[1]+".fna")):
			pass
			p_count += 1
		else:
			print "Fasta file from prokka assembly",p[1], "not found..."

	ens = ftplib.FTP('ftp.ensemblgenomes.org')
	ens.login()
	e_count = 0
	for e in ensembl:
		if e == "release-32":
			for strain in ensembl[e]:
				print "release-32 outdated...skipping",strain
		else:
			for line in open(os.path.join(genomedb,"{}.txt".format(e)),'r'):
				vals = line.rstrip().split("\t")
				if vals[6] in ensembl[e]:
					if os.path.exists(os.path.join(genomedb,"dna",vals[6]+".fna")):
						e_count += 1
					else:
						wd = '/pub/{}/bacteria/fasta/{}_collection/{}/dna'.format(e,"_".join(vals[5].split("_")[0:2]),vals[6])
						ens.cwd(wd)

						for filepath in ens.nlst():
							if filepath.endswith(".dna.toplevel.fa.gz"):
								o = open(os.path.join(genomedb,"dna",vals[6]+'.fna.gz'),'wb')
								ens.retrbinary("RETR " + filepath, o.write)
								o.close()
								cmds = ["gunzip",os.path.join(genomedb,"dna",vals[6]+'.fna.gz')]
								proc = subprocess.Popen(cmds)
								proc.wait()
						e_count += 1
				else:
					pass
	ens.close()

	r_count = 0
	for r in refseq:
		if os.path.exists(os.path.join(genomedb,"dna",r[1]+".fna")):
			r_count += 1
		else:
			assembly = r[0].split("_")[1]
			species = r[1]
			ens = ftplib.FTP('ftp.ncbi.nlm.nih.gov')
			ens.login()
			ens.cwd("/genomes/all/GCF/{}/{}/{}".format(assembly[0:3],assembly[3:6], assembly[6:9]))
			for f in ens.nlst():
				if f.split("_")[1] == assembly:
					ens.cwd(f)
					for g in ens.nlst():
						if g.endswith("_genomic.fna.gz"):
							filepath = os.path.join(genomedb,"dna","{}.fna.gz".format(species))
							dna = open(filepath,'wb')
							ens.retrbinary("RETR " + g, dna.write)
							dna.close()
							cmds = ["gunzip",os.path.join(genomedb,"dna","{}.fna.gz".format(species))]
							proc = subprocess.Popen(cmds)
							proc.wait()
							r_count += 1
			ens.close()

	e_tot = 0
	for e in ensembl:
		e_tot += len(ensembl[e])
	print p_count,"of", len(prokka), "prokka DNA files available."
	print e_count,"of", e_tot, "Ensembl DNA files available."
	print r_count, "of", len(refseq), "refseq DNA files available."
	return
