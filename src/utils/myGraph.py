# -*- coding: utf-8 -*-

# LibsDyogen version 1.0 (6/11/2015) -- modified for AGORA v3.0
# python v2.7 at least is needed
# Copyright © 2006-2021 IBENS/Dyogen, 2020-2021 EMBL-European Bioinformatics Institute, 2021 Genome Research Ltd : Matthieu MUFFATO, Thi Thuy Nga NGUYEN, Joseph LUCAS and Hugues ROEST CROLLIUS
# mail : agora@bio.ens.psl.eu
# Licences GLP v3 and CeCILL v2

###################################################
# Fonctions communes de traitement des diagonales #
###################################################

import collections
import itertools
import sys

from . import myFile
from . import myGenomes
from . import myMaths
from . import myTools


OrthosFilterType = myTools.Enum('NoFilter', 'InCommonAncestor', 'InBothSpecies')


def loadConservedPairsAnc(filename):
	pairwiseDiags = []
	f = myFile.openFile(filename, "r")
	for l in f:
		t = l.split("\t")
		pairwiseDiags.append(((int(t[0]), int(t[1])), (int(t[2]), int(t[3])), int(t[4])))
	f.close()
	return pairwiseDiags

#
# Chargement d'un fichier de paires conservees
################################################
def loadConservedPairs(filename, targets):
	print("Loading conserved pairs of %s ..." % filename, end=' ', file=sys.stderr)
	pairwiseDiags = collections.defaultdict(list)
	f = myFile.openFile(filename, "r")
	for l in f:
		t = l[:-1].split("\t")
		if t[0] in targets:
			pairwiseDiags[t[0]].append( ((int(t[1]),int(t[2])), (int(t[3]),int(t[4])), int(t[7])) )
	f.close()
	print("OK", file=sys.stderr)
	return pairwiseDiags


#
# Chargement de blocs integres
################################
def loadIntegr(filename):
	integr = []
	singletons = set()
	print("Loading ancestral blocks of %s ..." % filename, end=' ', file=sys.stderr)
	f = myFile.openFile(filename, "r")
	for l in f:
		t = l.split("\t")
		diagA = [int(x) for x in t[2].split()]
		diagS = [int(x) for x in t[3].split()]
		diagW = [int(x) for x in t[4].split()]
		assert len(diagA) == len(diagS)
		assert len(diagA) == (len(diagW)+1)
		if len(diagA) == 1:
			singletons.update(diagA)
		else:
			integr.append((list(zip(diagA,diagS)),diagW))
	f.close()
	print(myMaths.myStats.txtSummary([len(x[0]) for x in integr]), "+", len(singletons), "singletons OK", file=sys.stderr)
	return (integr,singletons)

#
# Extrait toutes les diagonales entre deux genomes (eventuellement des singletons)
# Pour optimiser, on demande 
#   genome1 qui est un dictionnaire qui associe a chaque chromosome 
#     la liste des numeros des genes ancestraux sur ce chromosome
#   dic2 qui associe a un numero de gene ancestral ses positions sur le genome 2
########################################################################################
def iterateDiags(genome1, dic2, sameStrand):

	l1 = []
	l2 = []
	la = []
	lastPos2 = []
	lastS1 = 0
	
	# Parcours du genome 1
	for (i1,(j1,s1)) in enumerate(genome1):
		presI2 = dic2[j1] if j1 >= 0 else []
		
		# On regarde chaque orthologue du gene
		for ((c2,i2,s2), (lastC2,lastI2,lastS2)) in itertools.product(presI2, lastPos2):
			# Chromosomes differents -> indiscutable
			if c2 != lastC2:
				continue
			# Meme brin
			if sameStrand:
				# Les brins initiaux imposent le sens de parcours (+1 ou -1)
				if i2 != lastI2 + lastS1*lastS2:
					continue
				# Le nouveau brin doit etre coherent
				if lastS1*s1 != lastS2*s2:
					continue
			else:
				# On demande juste que les deux genes soient cote a cote
				if abs(i2-lastI2) != 1:
					continue
			
			# On a passe tous les tests, c'est OK
			# On ecrit l'orthologue que l'on a choisi pour le coup d'avant (utile en cas de one2many, aucun effet si one2one)
			l2[-1] = (lastI2,lastS2)
			l2.append((i2,s2))
			lastPos2 = [(c2,i2,s2)]
			break

		# On n'a pas trouve de i2 satisfaisant, c'est la fin de la diagonale
		else:
			# On l'enregistre si elle n'est pas vide
			if len(l2) > 0:
				yield (lastPos2[0][0], l1, l2, la)
			# On recommence a zero
			lastPos2 = presI2
			l1 = []
			la = []
			# Pour que les diagonales de longueur 1 soient correctes
			l2 = [presI2[0][1:]] if len(presI2) > 0 else []
		
		l1.append((i1,s1))
		la.append(j1)
		lastS1 = s1
	
	if len(l2) > 0:
		yield (lastPos2[0][0], l1, l2, la)


#
# Proxy de generateur gerant la remise differee d'elements
############################################################
class queueWithBackup:

	def __init__(self, gen):
		self.gen = gen
		self.backup = collections.deque()
		self.todofirst = []
	
	def __iter__(self):
		return self
	
	# Le prochain element renvoye vient soit du buffer, soit du generateur principal
	def __next__(self):
		if len(self.todofirst) > 0:
			return self.todofirst.pop()
		return next(self.gen)
	
	# L'element reinsere est mis en attente
	def putBack(self, x):
		self.backup.appendleft(x)
	
	# Les elements sauvegardes sont mis dans un buffer de sortie et deviennent prioritaires
	def rewind(self):
		self.todofirst.extend(self.backup)
		self.backup = collections.deque()


#
# Lit les diagonales et les fusionne si elles sont separees par un trou pas trop grand
########################################################################################
def diagMerger(diagGen, sameStrand, largeurTrou):
	
	diagGen = queueWithBackup( (l1, l2, la, c2, l1[0][1]/l2[0][1], (min(l1)[0],max(l1)[0]), (min(l2)[0],max(l2)[0])) for (c2,l1,l2,la) in diagGen )

	# On rassemble des diagonales separees par une espace pas trop large
	for (la1,la2,laa,ca2,sa,(_,fina1),(deba2,fina2)) in diagGen:
		for curr in diagGen:
			(lb1,lb2,lba,cb2,sb,(debb1,finb1),(debb2,finb2)) = curr
			
			# Trou trop grand sur l'espece 1, aucune chance de le continuer
			if debb1 > (fina1+largeurTrou+1):
				diagGen.putBack(curr)
				break
			
			# Chromosomes differents de l'espece 2
			if ca2 != cb2:
				ok = False
			elif sameStrand:
				if sa != sb:
					ok = False
				elif sa > 0:
					ok = (fina2 < debb2 <= (fina2+largeurTrou+1))
				else:
					ok = (finb2 < deba2 <= (finb2+largeurTrou+1))
			else:
				ok = (min(abs(deba2-finb2), abs(debb2-fina2)) <= (largeurTrou+1))
			
			if ok:
				la1.extend(lb1)
				la2.extend(lb2)
				laa.extend(lba)
				fina1 = finb1
				deba2 = min(deba2,debb2)
				fina2 = max(fina2,finb2)
			else:
				diagGen.putBack(curr)

		yield (ca2,la1,la2,laa)
		diagGen.rewind()

#
# Procedure complete de calculs des diagonales a partir de 2 genomes, des orthologues et de certains parametres
################################################################################################################
def calcDiags(g1, g2, orthos, fusionThreshold=-1, sameStrand=True, orthosFilter=OrthosFilterType.NoFilter, minChromLength=0):

	# Ecrit les genomes comme suites de numeros de genes ancestraux
	def translateGenome(genome):
		newGenome = {}
		for c in genome.chrList[myGenomes.ContigType.Chromosome] + genome.chrList[myGenomes.ContigType.Scaffold]:
			tmp = [(orthos.getPositions(g.names),g.strand) for g in genome.lstGenes[c]]
			#assert set(len(x[0]) for x in tmp).issubset(set([0,1]))
			#assert set(list(x[0])[0].chromosome for x in tmp if len(x[0]) > 0).issubset([None])
			newGenome[c] = [(g.pop().index if len(g) > 0 else -1, strand) for (g,strand) in tmp]
		return newGenome
	newg1 = translateGenome(g1)
	newg2 = translateGenome(g2)

	# Marque les genes non presents dans l'intersection pour suppression future
	def markGeneIntersection():
		def usedValues(genome):
			val = set()
			for x in genome.values():
				val.update(i for (i,_) in x)
			return val
		def rewrite(genome, inters):
			for c in genome:
				back = len(genome[c])
				genome[c] = [(i,s) if i in inters else (-1,s) for (i,s) in genome[c]]
		val1 = usedValues(newg1)
		val2 = usedValues(newg2)
		inters = val1.intersection(val2)
		rewrite(newg1, inters)
		rewrite(newg2, inters)

	# Enleve les chromosomes trop petits
	def filterSize(genome):
		flag = False
		for c in list(genome.keys()):
			if len(genome[c]) < minChromLength:
				flag = True
				del genome[c]
		return flag

	# Enleve les genes sans lien ancestral
	def filterContent(genome, trans):
		for c in genome:
			tmp = [(i,x) for (i,x) in enumerate(genome[c]) if x[0] != -1]
			genome[c] = [x for (_,x) in tmp]
			trans[c] = dict((newi,trans[c].get(oldi,oldi)) for (newi,(oldi,_)) in enumerate(tmp))

	trans1 = collections.defaultdict(dict)
	trans2 = collections.defaultdict(dict)

	# Dans tous les cas, il faut filtrer sur la taille
	# On garde tous les genes
	if orthosFilter == OrthosFilterType.NoFilter:
		filterSize(newg1)
		filterSize(newg2)
	# On ne garde que les genes presents chez l'ancetre
	elif orthosFilter == OrthosFilterType.InCommonAncestor:
		filterContent(newg1, trans1)
		filterContent(newg2, trans2)
		filterSize(newg1)
		filterSize(newg2)
	# Ne conserve que les genes presents dans les deux genomes
	elif orthosFilter == OrthosFilterType.InBothSpecies:
		while True:
			markGeneIntersection()
			filterContent(newg1, trans1)
			filterContent(newg2, trans2)
			useful = filterSize(newg1) or filterSize(newg2)
			if not useful:
				break
	else:
		assert False
	
	# Pour chaque gene ancestral, ses positions dans le genome 2
	newLoc = [[] for x in range(len(orthos.lstGenes[None]))]
	for c in newg2:
		for (i,(ianc,s)) in enumerate(newg2[c]):
			if ianc != -1:
				newLoc[ianc].append( (c,i,s) )

	for c1 in newg1:
		src = iterateDiags(newg1[c1], newLoc, sameStrand)
		if (fusionThreshold > 0) or (not sameStrand):
			src = diagMerger(src, sameStrand, fusionThreshold)
		if orthosFilter != OrthosFilterType.NoFilter:
			for (c2,d1,d2,da) in src:
				yield ((c1,[(trans1[c1][i1],s1) for (i1,s1) in d1]), (c2,[(trans2[c2][i2],s2) for (i2,s2) in d2]), da)
		else:
			for (c2,d1,d2,da) in src:
				yield ((c1,d1), (c2,d2), da)

#@myTools.memoize
def revGene(gene):
	(x,sx) = gene
	# sx is an integer for the graph construction of AGORA
	# sx = +1 or -1, standard case
	# however sometimes sx is not defined, and we chose 0 and the reverse of 0 is 10
	if sx == 0:
		return (x, 10)
	elif sx == 10:
		return (x, 0)
	else:
		return (x,-sx)

#
# Un ensemble de diagonales que l'on represente comme un graphe ou les noeuds sont les genes
##############################################################################################
class WeightedDiagGraph:

	#
	# Constructeur
	################
	def __init__(self):
		# Les aretes du graphe et les orientations relatives des genes
		def newDicInt():
			return collections.defaultdict(int)
		self.aretes = collections.defaultdict(newDicInt)

	#
	# Insere un lien pondere
	##########################
	def addLink(self, xsx, ysy, weight):
		if xsx[0] == ysy[0]:
			return
		rxsx = revGene(xsx)
		rysy = revGene(ysy)
		self.aretes[xsx]
		self.aretes[rxsx]
		self.aretes[ysy]
		self.aretes[rysy]
		self.aretes[xsx][ysy] += weight
		self.aretes[rysy][rxsx] += weight
		
		
	#
	# Insere une diagonale avec poids fixe
	########################################
	def addDiag(self, diag, weight=1):
		for (xsx,ysy) in myTools.myIterator.slidingTuple(diag):
			self.addLink(xsx, ysy, weight)
		
	#
	# Insere une diagonale avec des poids variants
	################################################
	def addWeightedDiag(self, diag, weights):
		assert len(diag) == (len(weights)+1)
		for ((xsx,ysy),w) in zip(myTools.myIterator.slidingTuple(diag), weights):
			self.addLink(xsx, ysy, w)
		
	#
	# Affiche le graphe initial
	#############################
	def printIniGraph(self):
		print("INIGRAPH %d {" % len(self.aretes))
		for (xsx,l) in self.aretes.items():
			if len(l) > 0:
				print("\t", xsx, "> {%d}" % len(l))
				for (ysy,w) in l.items():
					print("\t\t", ysy, "[%s]" % w)
		print("}")


	#
	# Garde successivement les aretes de meilleur poids tant qu'elles n'introduisent pas de carrefour ou de cycle
	##############################################################################################################
	def cleanGraphTopDown(self, minimalWeight, searchLoops=True):
		allEdges = []
		res = {}
		allSucc = {}
		allPred = {}
		allNodes = set(self.aretes)
		
		if searchLoops:
			for (xsx,l) in list(self.aretes.items()):
				if len(l) != 2:
					continue
				ll = [(len(self.aretes[revGene(ysy)]),ysy) for ysy in l]
				ll.sort()
				if (ll[0][0] != 1) or (ll[1][0] != 2):
					continue
				target = ll[1][1]
				next = ll[0][1]
				if l[target] < l[next]:
					continue
				length = 0
				while next != target:
					if len(self.aretes[next]) != 1:
						break
					if len(self.aretes[revGene(next)]) != 1:
						break
					next = list(self.aretes[next])[0]
					length += 1
				else:
					print("loop", length, xsx, target)
					self.aretes[xsx].pop(target)
					self.aretes[revGene(target)].pop(revGene(xsx))

		for (xsx,l) in self.aretes.items():
			allSucc[xsx] = []
			allPred[xsx] = []
			for (ysy,c) in l.items():
				if (c >= minimalWeight) and (xsx < ysy):
					allEdges.append( (c,xsx,ysy) )
		
		def addEdge(c, xsx, ysy):
			# On ecrit l'arete
			res[xsx] = (ysy, c)
			# Tables de detection de cycles
			allSucc[xsx] = [ysy] + allSucc[ysy]
			assert len(allSucc[ysy]) == len(set(allSucc[ysy]))
			for t in allPred[xsx]:
				allSucc[t].extend(allSucc[xsx])
				assert len(allSucc[t]) == len(set(allSucc[t]))
			allPred[ysy] = [xsx] + allPred[xsx]
			assert len(allPred[ysy]) == len(set(allPred[ysy]))
			for t in allSucc[ysy]:
				allPred[t].extend(allPred[ysy])
				assert len(allPred[t]) == len(set(allPred[t]))

		allEdges.sort(reverse = True)
		for (c,xsx,ysy) in allEdges:
			rxsx = revGene(xsx)
			rysy = revGene(ysy)

			# 3 cas ambigus
			if xsx in res:
				# > 1 successeur
				print("not used /successor", xsx, ysy, c)
				continue
			if rysy in res:
				# > 1 predecesseur
				print("not used /predecessor", xsx, ysy, c)
				continue
			if xsx in allSucc[ysy]:
				# Cycle
				print("not used /cycle", xsx, ysy, c)
				continue
			assert rysy not in allSucc[rxsx]

			addEdge(c, xsx, ysy)
			addEdge(c, rysy, rxsx)
	
		allNodes.difference_update(res)
		allNodes.difference_update(ysy for (ysy,_) in res.values())

		self.aretes = res
		self.singletons = allNodes

		# Affichage du graphe
		print("GRAPH %d {" % len(self.aretes))
		for (tx,(ty,c)) in self.aretes.items():
			print("\t", tx, ">", ty, "[%s]" % c)
		print("}")



	#
	# Construit un graphe dans lequel les suites d'aretes sans carrefours ont ete reduites
	#
	def getBestDiags(self):

		# Renvoie le chemin qui part de src en parcourant les aretes
		def followSommet(xsx):
		
			res = []
			scores = []
			print("begin", xsx)
			
			# On part de src et on prend les successeurs jusqu'a la fin du chemin
			while xsx in self.aretes:
				# On marque notre passage
				res.append(xsx)
				
				# Le prochain noeud a visiter
				(ysy,c) = self.aretes.pop(xsx)
				print("pop edge", xsx, ysy, c)
				scores.append(c)

				# Traitement de l'arete inverse
				assert self.aretes.pop(revGene(ysy)) == (revGene(xsx),c)
				
				# Passage au suivant
				xsx = ysy

			res.append(xsx)
			assert revGene(xsx) not in self.aretes
			assert len(res) >= 2

			print("end", len(res), res)

			return (res,scores)

		
		# On cherche les extremites pour lancer les blocs integres
		todo = [xsx for xsx in self.aretes if revGene(xsx) not in self.aretes]
		for xsx in todo:
			if xsx in self.aretes:
				yield followSommet(xsx)
			
		assert len(self.aretes) == 0

		# On a isole quelques noeuds
		assert self.singletons == set(revGene(xsx) for xsx in self.singletons)
		for (x,sx) in self.singletons:
			if sx == 1:
				print("singleton", x)
				yield ([(x,1)],[])


