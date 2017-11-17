from Bio import SeqIO
from Bio.Graphics import GenomeDiagram
from Bio.SeqFeature import FeatureLocation
import seaborn as sns
import random
import argparse, os, sys, math
from reportlab.lib.units import inch
from ete3 import Tree, TreeStyle, NodeStyle
import matplotlib.colors as colors

def parse_genbank(g):
	if g[2] == "ensembl" or g[2] == "refseq":
		FIELD = "protein_id"
	elif g[2] == "prokka":
		FIELD = "locus_tag"
	for seq in SeqIO.parse(open(g[0],'r'),"genbank"):
		for feat in seq.features:
			if feat.type == "CDS":
				try:
					if feat.qualifiers[FIELD][0] == g[1]:
						print "Found", g[1], "in", seq.id
						print "loc:",feat.location
						return seq, (int(feat.location.start), int(feat.location.end))
				except KeyError:
					pass

	print "locus not found. try again."
	sys.exit()
	return

def find_homologs(GD, locus_tags, locus_mat,hl_groups):
	groups = {}
	count = 0
	print locus_tags
	for line in open(locus_mat,'r'):
		vals = [v.split(".")[0] for v in [x.split(";") for x in line.rstrip().split("\t")] for v in v]
		found = []
		for strain in locus_tags:
			for t in locus_tags[strain]:
				if t in vals:
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
	print groups
	return groups

def change_colors(GD, groups):
	cl = [colors.HexColor(c) for c in sns.cubehelix_palette(len(set(groups.values())),dark=0.1,light=0.9,rot=2.5).as_hex()]
	random.shuffle(cl)
	for t in GD.get_tracks():
		for s in t.get_sets():
			for feat in s.get_features():
				if feat.name.split(".")[0] in groups:
					if groups[feat.name.split(".")[0]] == "HIGHLIGHT":
						feat.color = colors.HexColor(u"#4df92a")
					else:
						feat.color = cl[groups[feat.name.split(".")[0]]]
	return

def make_tracks(seq, span, coords, g, GD, count, locus_tags, labels):

	if g[2] == "ensembl" or g[2] == "refseq":
		FIELD = "protein_id"
	elif g[2] == "prokka":
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
				print newloc
				if g[2] == "prokka" or g[2] == "refseq":
					feature_set.add_feature(feat, sigil="BIGARROW", arrowshaft_height=1, arrowhead_length=.4,color="#D3D3D3", \
						label=labels,name=feat.qualifiers['product'][0],label_strand=1,label_size = 8,label_position="middle", label_angle=20, \
						border=colors.black)
				try:
					feature_set.add_feature(feat, sigil="BIGARROW", arrowshaft_height=1, arrowhead_length=.4,color="#D3D3D3", \
						label=labels,name=feat.qualifiers[FIELD][0],label_strand=-1,label_size = 8,label_position="middle", label_angle=90, \
						border=colors.black)
					locus_tags[g[0].split(".")[0]].append(feat.qualifiers[FIELD][0].split(".")[0])
				except KeyError:
					pass
		elif feat.type == "gene":
			if g[2] == "ensembl" or g[2] == "refseq":
				if int(feat.location.start) > (coords[0]-(span/2)) and int(feat.location.end) < (coords[1]+(span/2)):
					newloc = FeatureLocation(feat.location.start-(coords[0]-(span/2)),feat.location.end-(coords[0]-(span/2)),strand=feat.strand)
					feat.location = newloc
					try:
						feature_set.add_feature(feat, sigil="BIGARROW", arrowshaft_height=1, arrowhead_length=.4,color="#D3D3D3", \
							label=labels,name=feat.qualifiers['note'][0],label_strand=1,label_size = 8,label_position="middle", label_angle=20, \
							border=colors.black)
					except KeyError:
						pass

	return


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
