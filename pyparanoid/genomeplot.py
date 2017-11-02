from Bio import SeqIO
from Bio.Graphics import GenomeDiagram
from Bio.SeqFeature import FeatureLocation
import seaborn as sns
import random
import argparse, os, sys
import seaborn as sns
from reportlab.lib.units import inch

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
