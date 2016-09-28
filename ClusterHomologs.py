#Ryan A. Melnyk
#schmelnyk@gmail.com
#UBC Microbiology - Haney Lab


import os, argparse, subprocess, itertools

def parse_args():
	parser = argparse.ArgumentParser(description='''
Script for formatting the PyParanoid output for a format compatible with MCL.
Runs the clustering process as well, providing stats on the output.
	''')
	parser.add_argument('outdir', type=str,help='location of PyParanoid output')
	return parser.parse_args()

def create_abc_file(outdir):
	o = open(os.path.join(outdir,"input.abc"),'w')
	count = len(os.listdir(os.path.join(outdir,"paranoid_output")))
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
			print "\tDone!"
		elif count % 100 == 0:
			print "\t"+str(count), "remaining..."
		else:
			pass
	o.close()
	return

def run_mcxload(outdir):
	cmds = ["mcxload","--stream-mirror","-abc",os.path.join(outdir,"input.abc"),"-o",os.path.join(outdir,"data.mci"),"-write-tab",os.path.join(outdir,"data.tab")]
	proc = subprocess.Popen(cmds)
	proc.wait()
	return

def cluster(outdir):
	for i in ["1.4","2","4","6"]:
		cmds = ["mcl",os.path.join(outdir,"data.mci"),"-I",i,"-o",os.path.join(outdir,"mcl.{}.out".format(i.replace(".","")))]
		proc = subprocess.Popen(cmds)
		proc.wait()
	return

def dump(outdir):
	for i in ["1.4","2","4","6"]:
		cmds = ["mcxdump", "-imx", os.path.join(outdir,"data.mci"),"-tabr",os.path.join(outdir,"data.tab"),"-icl",os.path.join(outdir,"mcl.{}.out".format(i.replace(".",""))),"-o",os.path.join(outdir,"clusters.{}.txt".format(i.replace(".","")))]
		proc = subprocess.Popen(cmds)
		proc.wait()
	return

def parse_clusters(outdir):
	for i in ["1.4","2","4","6"]:
		for line in open(os.path.join(outdir,"clusters.{}.txt".format(i.replace(".",""))),"r"):
			vals = line.rstrip().split()
	return

def main():
	args = parse_args()
	outdir = os.path.abspath(args.outdir)

	# create_abc_file(outdir)
	# run_mcxload(outdir)
	# cluster(outdir)
	dump(outdir)
	parse_clusters(outdir)

if __name__ == '__main__':
	main()
