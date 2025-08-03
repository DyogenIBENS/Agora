#!/usr/bin/env python3

import sys
import operator
import collections

import utils.myFile
import utils.myTools
from utils.myTools import file
import utils.myGenomes
import utils.myPsOutput
import utils.myKaryoDrawer

__doc__ = """
	Compare two genomes thanks to an orthologues list (ancGenes or 2 columns orthologues gene list)
		
		-can draw dotplot or karyotype
		-can print orthologous gene list, ortologous chromosome list, gene difference, orthologs count
		
	
	Usage:	misc.compareGenomes.py genesST.Homo.sapiens.list.bz2 genesST.Mus.musculus.list.bz2 ancGenes.Euarchontoglires.list.bz2  > dotplot_Human_Mouse.ps
		misc.compareGenomes.py genesST.Homo.sapiens.list.bz2 genesST.Mus.musculus.list.bz2 ancGenes.Euarchontoglires.list.bz2 -mode=drawKaryotype -minChrSize=200 > Karyo_human_min200genes.ps
"""

modes = ["drawMatrix", "drawKaryotype", "printOrthologuesList", "printOrthologuesCount", "printGeneDiff", "printOrthologousChrom"]

arguments = utils.myTools.checkArgs(
	[("studiedGenome", file), ("referenceGenome", file), ("orthologuesList", file)],
	[("includeGaps", bool, False), ("includeScaffolds", bool, True), ("includeRandoms",bool,False), ("includeNones",bool,False),
	("reverse",bool,False),
	("mode",str,modes),
	("orthoslist:fullgenenames",bool,False),
	("orthoschr:minHomology",int,90),
	("minChrSize",int,0),
	("matrix:scaleY",bool,False), ("matrix:pointSize",float,-1), ("sortBySize",bool,False),
	("matrix:colorFile",str,""), ("matrix:defaultColor",str,"black"), ("matrix:penColor",str,"black"),
	("karyo:landscape",bool,False),
	("ps:backgroundColor",str,"")],
	__doc__
)


# Chargement des fichiers
genesAnc = utils.myGenomes.Genome(arguments["orthologuesList"])
genome1 = utils.myGenomes.Genome(arguments["studiedGenome"], ancGenes=genesAnc)
genome2 = utils.myGenomes.Genome(arguments["referenceGenome"], ancGenes=genesAnc)
if arguments["reverse"]:
	(genome1,genome2) = (genome2,genome1)

chr1 = []
chr2 = []
chr1.extend(genome1.chrList[utils.myGenomes.ContigType.Chromosome])
chr2.extend(genome2.chrList[utils.myGenomes.ContigType.Chromosome])
if arguments["includeScaffolds"]:
	chr1.extend(genome1.chrList[utils.myGenomes.ContigType.Scaffold])
	chr2.extend(genome2.chrList[utils.myGenomes.ContigType.Scaffold])
if arguments["includeRandoms"]:
	chr1.extend(genome1.chrList[utils.myGenomes.ContigType.Random])
	chr2.extend(genome2.chrList[utils.myGenomes.ContigType.Random])
if arguments["includeNones"]:
	chr1.extend(genome1.chrList[utils.myGenomes.ContigType.Unknown])
	chr2.extend(genome2.chrList[utils.myGenomes.ContigType.Unknown])
print(len(chr1), len(chr2), file=sys.stderr)

chr1 = [c for c in chr1 if len(genome1.lstGenes[c]) >= arguments["minChrSize"]]
chr2 = [c for c in chr2 if len(genome2.lstGenes[c]) >= arguments["minChrSize"]]

table12 = genome1.buildOrthosTable(chr1, genome2, chr2, arguments["includeGaps"], genesAnc)
#print >>sys.stderr, table12
table21 = genome2.buildOrthosTable(chr2, genome1, chr1, arguments["includeGaps"], genesAnc)

#
# Matrice avec les genes orthologues
######################################
def drawMatrix():

	# Matrice

	print("Affichage ", end=' ', file=sys.stderr)

	if arguments["sortBySize"]:
		chr1.sort(key=lambda c: len(genome1.lstGenes[c]), reverse=True)
		chr2.sort(key=lambda c: len(genome2.lstGenes[c]), reverse=True)

	utils.myPsOutput.printPsHeader()
	if arguments["ps:backgroundColor"] != "":
		utils.myPsOutput.drawBox(0,0, 21,29.7, arguments["ps:backgroundColor"], arguments["ps:backgroundColor"])
	sys.stderr.write('.')
	colors = utils.myGenomes.Genome(arguments["matrix:colorFile"]) if arguments["matrix:colorFile"] != "" else None

	# Initialisations
	nb = sum([len(table12[c]) for c in table12])
	scaleX = 19. / float(nb)
	scaleY = 19. / float(sum([len(table21[c]) for c in table21])) if arguments["matrix:scaleY"] else scaleX
	dp = scaleX if arguments["matrix:pointSize"] < 0 else arguments["matrix:pointSize"]
	sys.stderr.write('.')

	def prepareGenome(dicOrthos, lst, func):
		i = 0
		y = 0
		lstNum = {}
		for c in lst:
			func(c, y, len(dicOrthos[c]))
			y += len(dicOrthos[c])
			for (gene,_) in dicOrthos[c]:
				lstNum[(c,gene)] = i
				i += 1
		func(None, y, None)
		return lstNum

	dl1 = float(sum([len(table21[c]) for c in table21])) * scaleY
	def line1(c, x, l):
		utils.myPsOutput.drawLine(1 + x*scaleX, 1, 0, dl1, arguments["matrix:penColor"])
		if c:
			utils.myPsOutput.drawText(1 + (x+l/2)*scaleX, 0.7, c, arguments["matrix:penColor"])

	def line2(c, x, l):
		utils.myPsOutput.drawLine(1, 1 + x*scaleY, 19, 0, arguments["matrix:penColor"])
		if c:
			print("90 rotate")
			utils.myPsOutput.drawText(1 + (x+l/2)*scaleY, -0.9, c, arguments["matrix:penColor"])
			print("-90 rotate")

	lstNum1 = prepareGenome(table12, chr1, line1)
	sys.stderr.write('.')
	lstNum2 = prepareGenome(table21, chr2, line2)
	sys.stderr.write('.')

	print("0 setlinewidth")


	for c1 in table12:
		for (i1,t) in table12[c1]:
			xx = 1 + float(lstNum1[(c1,i1)]) * scaleX
			for (c2,i2) in t:

				coul = arguments["matrix:defaultColor"]
				if colors is not None:
					tmp = set(colors.getPosition(genome1.lstGenes[c1][i1].names + genome2.lstGenes[c2][i2].names))
					for (c,i) in genesAnc.getPosition(genome1.lstGenes[c1][i1].names + genome2.lstGenes[c2][i2].names):
						tmp.update(colors.getPosition(genesAnc.lstGenes[c][i].names))
					if len(tmp) > 0:
						coul = tmp.pop()[0]

				yy = 1 + lstNum2[(c2,i2)]*scaleY
				utils.myPsOutput.drawBox(xx, yy, dp, dp, coul, coul)

	utils.myPsOutput.drawText(4, 0.3, arguments["referenceGenome"] if arguments["reverse"] else arguments["studiedGenome"], arguments["matrix:penColor"])
	print("90 rotate")
	utils.myPsOutput.drawText(4, -0.5, arguments["studiedGenome"] if arguments["reverse"] else arguments["referenceGenome"], arguments["matrix:penColor"])
	print("-90 rotate")
	utils.myPsOutput.printPsFooter()
	print(" OK", file=sys.stderr)


#
# Dessine le caryotype de la premiere espece avec les couleurs des chromosomes de la seconde
##############################################################################################
def drawKaryotype():

	(lx,ly) = utils.myPsOutput.printPsHeader(arguments["karyo:landscape"])
	if arguments["ps:backgroundColor"] != "":
		utils.myPsOutput.drawBox(0,0, lx,ly, arguments["ps:backgroundColor"], arguments["ps:backgroundColor"])

	data = []
	for c in chr1:
		newl = []
		for (_,val) in table12.get(c, []):
			if len(val) == 0:
				newl.append( None )
			else:
				newl.append( val[0][0] )
		data.append( (c, newl) )

	color_chrom = {}
	for (i,c) in enumerate(sorted(chr2, key=lambda c: (len(genome2.lstGenes[c]), c), reverse=True)):
		color_chrom[c] = str(i+1)

	print("Affichage ...", end=' ', file=sys.stderr)
	utils.myKaryoDrawer.drawKaryo(data, arguments, x0=1, y0=1, lx=lx-2, ly=ly-2, bysize=arguments["sortBySize"], color_chrom=color_chrom)
	utils.myPsOutput.printPsFooter()
	print("OK", file=sys.stderr)


#
# Tableau avec le nombre d'orthologues pour chaque paire de chromosomes
#########################################################################
def printOrthologuesCount():

	print(utils.myFile.myTSV.printLine([""] + chr2))
	for c1 in chr1:
		count = collections.defaultdict(int)
		for (i1,t) in table12[c1]:
			for (c2,i2) in t:
				count[c2] += 1
		print(utils.myFile.myTSV.printLine([c1] + [count[c2] for c2 in chr2]))


#
# Pour chaque gene de la premiere espece, les orthologues dans l'autre espece
###############################################################################
def printOrthologuesList():
	# Fichier avec les noms des paires de genes orthologues et leurs coordonnees
	def printGene(g):
		s = list(g)
		s[-1] = "/".join(s[-1]) if arguments["orthoslist:fullgenenames"] else s[-1][0]
		return s
	for c1 in chr1:
		for (i1,t) in sorted(table12[c1]):
			g1 = genome1.lstGenes[c1][i1]
			for (c2,i2) in sorted(t):
				print(utils.myFile.myTSV.printLine( printGene(g1) + printGene(genome2.lstGenes[c2][i2]) ))

#
# Affiche les differences dans le contenu en genes
####################################################
def printGeneDiff():
	def getGeneTxt(g):
		return "/".join(g.names) + ":%s:%d-%d:%d" % g[:4]

	all = set()
	combin = utils.myTools.myCombinator()
	for c1 in table12:
		for (i1,t) in table12[c1]:
			combin.addLink( [(1,c1,i1)] + [(2,c2,i2) for (c2,i2) in t] )
	for c2 in table21:
		for (i2,t) in table21[c2]:
			combin.addLink( [(2,c2,i2)] + [(1,c1,i1) for (c1,i1) in t] )
	for g in combin:
		e1 = [getGeneTxt(genome1.lstGenes[c][i]) for (x,c,i) in g if x == 1]
		e2 = [getGeneTxt(genome2.lstGenes[c][i]) for (x,c,i) in g if x == 2]
		if len(e1) == 0:
			assert len(e2) == 1
			print("+", end=' ')
		elif len(e2) == 0:
			assert len(e1) == 1
			print("-", end=' ')
		elif (len(e1) == 1) and (len(e2) == 1):
			print("=", end=' ')
		elif (len(e1) > 1) and (len(e2) == 1):
			print("--", end=' ')
		elif (len(e1) == 1) and (len(e2) > 1):
			print("++", end=' ')
		else:
			print("**", end=' ')
		print(" ".join(e1 + e2))

#
# Affiche les rearrangements de chromosomes
#############################################
def printOrthologousChrom():

	for c1 in chr1:
		count = collections.defaultdict(int)
		for (i1,t) in table12[c1]:
			for (c2,i2) in t:
				count[c2] += 1
		res = [c1]
		t = sorted(iter(count.items()), key=operator.itemgetter(1))
		n = (sum(count.values()) * arguments["orthoschr:minHomology"]) / 100
		while n > 0:
			x = t.pop()
			res.append("%s (%d)" % x)
			n -= x[1]
		print(utils.myFile.myTSV.printLine(res))


eval(str(arguments["mode"]))()

