import os

def check_db(outdir):
	assemblies = []
	if "genome_metadata.txt" not in os.listdir(outdir):
		pass
	else:
		for line in open(os.path.join(outdir,"genome_metadata.txt"),'r'):
			if line.startswith("assembly_id"):
				continue
			else:
				assemblies.append(line.rstrip().split("\t")[2])
	return assemblies

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
