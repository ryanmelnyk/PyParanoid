from Bio import SeqIO
from Bio.Graphics import GenomeDiagram
from Bio import GenBank
from Bio.SeqFeature import FeatureLocation
import seaborn as sns
import random
import argparse, os, sys, math
from ete3 import Tree, TreeStyle, NodeStyle
import matplotlib.colors as colors
from reportlab.lib import colors as rcolors
import numpy as np

def match_seqs(fastafile,outdir):
	os.system("phmmer --tblout {} {} {}".format(os.path.join(os.path.abspath(outdir),".phmmer.hits"),os.path.abspath(fastafile),os.path.join(outdir,"all_groups.faa")))
	hit_scores = {}
	for line in open(os.path.join(outdir,".phmmer.hits"),'r'):
		if line.startswith("#"):
			continue
		else:
			vals = line.rstrip().split()
			hit_scores[vals[0]] = float(vals[5])
	hit_gene_lengths = {}
	for seq in SeqIO.parse(open(os.path.join(outdir,"all_groups.faa"),'r'),'fasta'):
		if seq.id in hit_scores:
			hit_gene_lengths[seq.id] = len(seq.seq)
	norm_scores = {}
	for group in hit_scores:
		norm_scores[group] = hit_scores[group]/float(hit_gene_lengths[group])
	for i in sorted(norm_scores.iteritems(), key=lambda (k,v): (v,k), reverse=True):
		if i[1] > 1.0:
			score = "Good match!"
		if i[1] > 0.5 and i[1] < 1.0:
			score = "Possible match..."
		if i[1] < 0.5:
			score = "Probably not a very good match..."
		print i[0],round(i[1],3), score
	return

def add_group_to_tree(group,treefile,outdir,to_compress=False):
	if to_compress:
		compress = to_compress.split(",")
	else:
		compress = []
	for line in open(os.path.join(outdir,"homolog_matrix.txt"),'r'):
		if line.startswith("\t"):
			header = line.rstrip().split("\t")[1:]
		if not line.startswith(group):
			continue
		else:
			vals = [int(x) for x in line.rstrip().split("\t")[1:]]
	groupdata = dict(zip(header,vals))

	ts = TreeStyle()
	tree = Tree(os.path.abspath(treefile))

	pal = sns.cubehelix_palette(rot=-.4, n_colors=13)
	for node in tree.iter_descendants("preorder"):
		this_node = []
		nstyle = NodeStyle()
		nstyle["shape"] = "circle"

		if node.is_leaf():
			try:
				if groupdata[node.name] > 0:
					nstyle["fgcolor"] = colors.rgb2hex(pal[12])
				else:
					nstyle["fgcolor"] = colors.rgb2hex(pal[0])
			except KeyError:
				nstyle["fgcolor"] = colors.rgb2hex(pal[0])
		else:
			species = {}
			for x in node.iter_descendants("preorder"):
				if x.is_leaf():
					this_node.append(x.name)
					s = x.name.split("_")[1]
					if s in species:
						species[s] += 1
					else:
						species[s] = 1
			for c in compress:
				try:
					if float(species[c])/float(len(this_node)) > 0.95:
						nstyle["draw_descendants"] = False
						node.name = "{} clade".format(c)
				except KeyError:
					pass
			count = 0
			for t in this_node:
				try:
					if groupdata[t] > 0:
						count +=1
				except KeyError:
					pass
			v = int(round(float(count)/float(len(this_node))*12))
			nstyle["fgcolor"] = colors.rgb2hex(pal[v])
			nstyle["size"] = 3*math.sqrt(len(this_node))
		node.set_style(nstyle)

	return tree


def subset_matrix(strains,outdir):
	infile = open(os.path.join(outdir,"homolog_matrix.txt"),'r')
	header_line = infile.readline().rstrip().split("\t")
	indices = [header_line.index(s) for s in strains]

	lines = []
	groups = []
	for line in infile:
		vals = line.rstrip().split("\t")
		# convert input data into boolean presence/absence
		lines.append([int(bool(int(vals[i]))) for i in indices])
		groups.append(vals[0])
	a = np.stack(lines)
	return a, groups


def dump_matrix(cc,strains,outfile):
	o = open(outfile,'w')
	o.write("\t{}\n".format("\t".join(strains)))
	for i in range(0,len(strains)):
		o.write("{}\t{}\n".format(strains[i],"\t".join([str(x) for x in cc[i]])))
	o.close()
	return

def find_unique_genes(a,strains,groups):
	unique = {}
	missing = {}
	common = []
	for s in strains:
		unique[s] = []
		missing[s] = []
	for i in range(0,a.shape[0]):
		nz = np.nonzero(a[i])
		if len(nz[0]) == 1:
			unique[strains[nz[0][0]]].append(groups[i])
		elif len(nz[0]) == len(strains)-1:
			missing[strains[np.nonzero(a[i] < 1)[0][0]]].append(groups[i])
		elif len(nz[0]) == len(strains):
			common.append(groups[i])
		else:
			pass
	for s in strains:
		print s
		print "\t", len(unique[s]), "unique"
		print "\t", len(missing[s]), "missing"
	print len(common), "common to all strains."
	return {"unique":unique,"missing":missing,"common":common}

def find_unique_loci(strain,outdir,uniq_info):
	locusfile = open(os.path.join(outdir,"locustag_matrix.txt"),'r')
	header_line = locusfile.readline().rstrip().split("\t")
	i = header_line.index(strain)

	unique_loci = []
	for line in locusfile:
		vals = line.rstrip().split("\t")
		if vals[0] in uniq_info["unique"][strain]:
			for x in vals[i].split(";"):
				unique_loci.append(x)
	return unique_loci

def plot_unique_genome_diagram(gbk, unique_loci):
	parser = GenBank.FeatureParser()
	fhandle = open(gbk, 'r')
	genbank_entry = parser.parse(fhandle)
	fhandle.close()

	gdd = GenomeDiagram.Diagram(gbk)
	gd_track_for_features = gdd.new_track(1, name="CDS",scale_smalltick_interval=100000)
	gdfs = gd_track_for_features.new_set()
	for feature in genbank_entry.features:
		if feature.type == 'CDS':
			feature.strand = 1
			if feature.qualifiers['locus_tag'][0] in unique_loci:
				gdfs.add_feature(feature, color=rcolors.HexColor("#93341F"))
			else:
				gdfs.add_feature(feature, color=rcolors.HexColor("#058F45"))
	gdd.draw(format='circular', orientation='landscape',tracklines=0, pagesize='A5', fragments=5, circular=1)
	return gdd

def synteny_check(gbk,outdir,map_to,strains,outfile):
	seqs = {}
	for seq in SeqIO.parse(open(gbk,'r'),'genbank'):
		seqs[seq.id] = []
		for feat in seq.features:
			if feat.type == "CDS":
				try:
					seqs[seq.id].append(feat.qualifiers["locus_tag"][0])
				except KeyError:
					seqs[seq.id].append(feat.qualifiers["protein_id"][0])
	groupdict = {}

	locustags = []
	for seq in seqs:
		for s in seqs[seq]:
			locustags.append(s)

	header = open(os.path.join(outdir,"locustag_matrix.txt"),'r').readline().rstrip().split("\t")
	i = header.index(map_to)
	indices = [header.index(x) for x in strains]
	for line in open(os.path.join(outdir,"locustag_matrix.txt"),'r'):
		vals = [v.split(".")[0] for v in line.rstrip().split("\t")[i].split(";")]
		for s in locustags:
			if s in vals:
				groupdict[s] = line.rstrip().split("\t")[0]

	groups = [g for g in groupdict.values()]
	datadict = {}
	for line in open(os.path.join(outdir,"locustag_matrix.txt"),'r'):
		vals = line.rstrip().split("\t")
		if vals[0] in groups:
			datadict[vals[0]] = []
			for i in indices:
				if len(vals[i].split(";")) > 1:
					datadict[vals[0]].append("Multiple")
				else:
					datadict[vals[0]].append(vals[i])

	group_annotations = {}
	for line in open(os.path.join(outdir,"group_descriptions.txt"),'r'):
		vals = line.rstrip().split("\t")
		counts = {}
		for x in set(vals[1:]):
			counts[x] = vals.count(x)
		group_annotations[vals[0]] = sorted(counts.items(), reverse=True, key = lambda x: x[1])[0][0]

	o = open(outfile,'w')
	for seq in seqs:
		o.write(">{}\n".format(seq))
		o.write("{}\tgroup\t{}\tannotation\n".format(map_to,"\t".join(strains)))
		for s in seqs[seq]:
			try:
				o.write("{}\t{}\t{}\t{}\n".format(s,groupdict[s],"\t".join(datadict[groupdict[s]]),group_annotations[groupdict[s]]))
			except KeyError:
				pass
	return

def _parse_genbank(g,genomedb):
	if g[2] == "ensembl" or g[2] == "NCBI":
		FIELD = "protein_id"
	elif g[2] == "prokka_in_house":
		FIELD = "locus_tag"
	for seq in SeqIO.parse(open(os.path.join(os.path.abspath(genomedb),"gbk",g[0]+".gbk"),'r'),"genbank"):
		for feat in seq.features:
			if feat.type == "CDS":
				try:
					if feat.qualifiers[FIELD][0] == g[1]:
						return seq, (int(feat.location.start), int(feat.location.end))
				except KeyError:
					pass
	print g
	print "locus not found. try again."
	return

def _make_tracks(seq, span, coords, g, GD, count, locus_tags, labels):

	if g[2] == "ensembl" or g[2] == "NCBI":
		FIELD = "protein_id"
	elif g[2] == "prokka_in_house":
		FIELD = "locus_tag"

	track = GD.new_track(count, height=1, name="CDS",\
		scale_ticks=False,scale_largeticks=False,scale_largetick_labels=False, scale_largetick_interval=10000,\
		scale_smallticks=False, scale_smalltick_labels=False, scale_smalltick_interval=1000,\
		greytrack=False, hide=False, scale=False
	)
	feature_set = track.new_set()


	count = 0
	for feat in seq.features:
		if feat.type == "CDS":
			if int(feat.location.start) > (coords[0]-(span/2)) and int(feat.location.end) < (coords[1]+(span/2)):
				newloc = FeatureLocation(feat.location.start-(coords[0]-(span/2)),feat.location.end-(coords[0]-(span/2)),strand=feat.strand)
				feat.location = newloc
				if g[2] == "prokka_in_house" or g[2] == "NCBI":
					feature_set.add_feature(feat, sigil="BIGARROW", arrowshaft_height=1, arrowhead_length=.4,color="#D3D3D3", \
						label=labels,name=feat.qualifiers['product'][0],label_strand=1,label_size = 8,label_position="middle", label_angle=20, \
						border=rcolors.black)
				try:
					feature_set.add_feature(feat, sigil="BIGARROW", arrowshaft_height=1, arrowhead_length=.4,color="#D3D3D3", \
						label=labels,name=feat.qualifiers[FIELD][0],label_strand=-1,label_size = 8,label_position="middle", label_angle=90, \
						border=rcolors.black)
					locus_tags[g[0].split(".")[0]].append(feat.qualifiers[FIELD][0].split(".")[0])
				except KeyError:
					pass

	return

def plot_genomic_regions(locustagfile,genomedb,pypdir,span=50000,hl_groups=[],labels=False):
	strains = []
	for line in open(os.path.abspath(locustagfile),'r'):
		vals = line.rstrip().split("\t")
		strains.append([vals[0],vals[1]])

	for line in open(os.path.join(os.path.abspath(genomedb),"genome_metadata.txt"),'r'):
		vals = line.rstrip().split("\t")
		for i in range(0,len(strains)):
			if strains[i][0] == vals[2]:
				strains[i].append(vals[6].split("-")[0])

	GD = GenomeDiagram.Diagram('gbk',"temp.pdf")
	count = 1
	locus_tags = {}
	for g in reversed(strains):
		if g[0] not in locus_tags:
			locus_tags[g[0]] = []
		contigseq, coords = _parse_genbank(g,genomedb)
		_make_tracks(contigseq, span, coords, g, GD, count, locus_tags, labels)
		count += 1

	groups = _find_homologs(GD, locus_tags,os.path.join(pypdir,"locustag_matrix.txt"),hl_groups,set([x[0] for x in strains]))
	_change_colors(GD, groups)
	return GD

def _find_homologs(GD, locus_tags, locus_mat,hl_groups,strains):
	groups = {}
	count = 0

	indices = [open(locus_mat,'r').readline().rstrip().split("\t").index(s) for s in strains]
	for line in open(locus_mat,'r'):
		vals = line.rstrip().split("\t")
		tags = set([v.split(".")[0] for v in [x.split(";") for x in [vals[i] for i in indices]] for v in v])
		found = []
		for strain in locus_tags:
			for t in locus_tags[strain]:
				if t in tags:
					found.append((strain,t))

		if vals[0] in hl_groups:
			for f in found:
				groups[f[1]] = "HIGHLIGHT"
		elif len(found) > 1:
			for f in found:
				groups[f[1]] = count
			count += 1
		else:
			pass
	return groups

def _change_colors(GD, groups):
	cl = [rcolors.HexColor(c) for c in sns.cubehelix_palette(len(set(groups.values())),dark=0.1,light=0.9,rot=2.5).as_hex()]
	random.shuffle(cl)
	for t in GD.get_tracks():
		for s in t.get_sets():
			for feat in s.get_features():
				if feat.name.split(".")[0] in groups:
					if groups[feat.name.split(".")[0]] == "HIGHLIGHT":
						feat.color = rcolors.HexColor(u"#4df92a")
					else:
						feat.color = cl[groups[feat.name.split(".")[0]]]
	return
