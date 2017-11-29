
#
# Fonctions communes de traitement des communautes
# Formats d'entree des donnees:
#  - Fonction de score
#  - Dictionnaire d'aretes
#  - Fichier "node1 node2 weight"
# Format de sortie des donnees:
#  - Liste pour chaque composante connexe
#    - Les noeuds de la composante
#    - Les coupures interessantes (score,alpha)
#    - Le dendogramme (liste des (alphas, fils, pere))
#    - Le dendogramme (instance de WalktrapDendogram)


import sys


# Lancer un walktrap a partir d'un graphe
# Permet de detecter les composantes connexes
###############################################
	
def doWalktrap(edges, **args):

	import utils.myTools
	import _walktrap

	print >> sys.stderr, "Computing connected components ...",
	# Les composantes connexes
	combin = utils.myTools.myCombinator()
	superdelete = []
	
	for (x,l) in edges.iteritems():
		delete = [y for y in l if l[y] <= 0.]
		for y in delete:
			del l[y]
		l.pop(x, None)
		if len(l) > 0:
			combin.addLink(l.keys() + [x])
		else:
			superdelete.append(x)
	
	res = []
	for x in superdelete:
		del edges[x]
		res.append( ([x], [], None, None) )

	print >> sys.stderr, "Launching walktrap ",
	n = len(edges)
	for nodes in combin:
		# Reindexation des noeuds
		indNodes = {}
		for (i,node) in enumerate(nodes):
			indNodes[node] = i
	
		(relevantCuts,dend) = _walktrap.doWalktrap(indNodes, edges, **args)

		# On doit revenir aux noms de noeuds originels
		def translate(x):
			if x < len(nodes):
				return nodes[x]
			else:
				return (x,)
		dend = [(cut,tuple(translate(f) for f in fils),translate(pere)) for (cut,fils,pere) in dend]
		res.append( (nodes, relevantCuts, dend, WalktrapDendogram(dend, nodes)) )
		sys.stderr.write('.')
	print >> sys.stderr, " OK"

	return res


# Lance un walktrap sur la liste de blocs l, en utilisant func pour le calcul du score
#######################################################################################
def clusterWithNb(l, funcScore, funcBestChoice, putLonelyinNone, **kwargs):

	(edges,none) = funcScore(l)
	print >> sys.stderr, "input", len(edges), len(none), len(l)

	#print >> sys.stderr, edges
	s = len(edges)
	res = doWalktrap(edges, **kwargs)

	# Chaque composante connexe
	chrOrder = []
	for (nodes,cuts,_,dend) in res:
		if len(cuts) == 0:
			chrOrder.append(nodes)
		else:
			(alpha,score,(clust,lonely)) = funcBestChoice([(alpha,score,dend.cut(alpha)) for (alpha,score) in cuts])
			assert sum(len(x) for x in clust) + len(lonely) == len(nodes)
			chrOrder.extend(clust)
			if putLonelyinNone:
				none.update(lonely)
			else:
				chrOrder.append(lonely)

	print >> sys.stderr, "output", len(edges), len(none), [len(x) for x in chrOrder]
	assert len(l) == (len(none) + sum(len(x) for x in chrOrder)), (len(l), len(none), sum(len(x) for x in chrOrder), [len(x) for x in chrOrder])
	return (chrOrder,none)


# Applique succesivement plusieurs clusterWithNb sur chaque liste de blocs dans ll
###################################################################################
def applyMultipleClust(ll, lfuncScore, lfuncBestChoice, putLonelyinNone, **kwargs):
	res = []
	none = []
	for l in ll:
		for (funcScore, funcBestChoice) in zip(lfuncScore, lfuncBestChoice):
			(l,n) = clusterWithNb(l, funcScore, funcBestChoice, putLonelyinNone, **kwargs)
			none.extend(n)
		res.extend(l)
	assert sum(len(l) for l in ll) == (len(none) + sum(len(l) for l in res))
	return (res,none)
	return res + [[x] for x in none]



###################################################
# Chargement d'un fichier de resultat de walktrap #
###################################################
def loadWalktrapOutput(f):
	
	# On charge les fusions
	allMerges = []
	lstFils = {}
	for line in f:
		if line == "\n":
			break
		
		l = line.split(':')
		scale = float(l[0])
		l = l[1].split('-->')
		
		allMerges.append( (scale,tuple([int(x) for x in l[0].split('+')]),int(l[1])) )

	allMerges.sort( reverse = True )

	lstCoup = []
	for line in f:
		try:
			# On extrait les lignes "alpha relevance"
			c = line.split()
			lstCoup.append( (float(c[0]),float(c[1])) )
		except ValueError:
			pass
	
	return (lstCoup, allMerges)


########################################################################################
# Le dendogramme resultat, que l'on peut couper a un niveau pour recuperer les classes #
########################################################################################
class WalktrapDendogram:

	def __init__(self, lstMerges, lstNodes):

		self.allMerges = lstMerges
		self.allMerges.sort(reverse = True)
		
		self.lstFils = {}
		for (_,fils,pere) in self.allMerges:
			self.lstFils[pere] = fils
		
		self.lstAll = lstNodes

	# Renvoie (la liste des clusters, les noeuds en dehors des clusters)
	def cut(self, scale):

		# On extrait les communautes correspondantes
		lstClusters = []
		fathersAlreadySeen = set()
		nodesNotSeen = set(self.lstAll)
		for (s,_,pere) in self.allMerges:
			if (s < scale) and (pere not in fathersAlreadySeen):
				cluster = []
				todo = [pere]
				while len(todo) > 0:
					father = todo.pop()
					if father in self.lstFils:
						fathersAlreadySeen.add(father)
						todo.extend(self.lstFils[father])
					else:
						nodesNotSeen.discard(father)
						cluster.append(father)
				lstClusters.append( cluster )
		return (lstClusters, list(nodesNotSeen))


####################################################
# Demande a l'utilisateur quelle partition choisir #
####################################################
def askPartitionChoice(dend, cuts):

	import utils.myMaths

	def mystr( (alpha,relevance,(clusters,lonely)) ):
		return "alpha=%f relevance=%f clusters=%d size=%d lonely=%d sizes={%s}" % \
			(alpha,relevance,len(clusters),sum([len(c) for c in clusters]),len(lonely),utils.myMaths.myStats.txtSummary([len(c) for c in clusters]))

	res = [(alpha,relevance,dend.cut(alpha)) for (alpha,relevance) in cuts]
	# Le choix par defaut
	if len(res) == 1:
		print >> sys.stderr, "1 seule possibilite"
		x = 0
	else:
		# Si on peut, on propose a l'utilisateur de choisir
		for x in res:
			print >> sys.stderr, "> " + mystr(x)
		import os
		if os.isatty(sys.stdin.fileno()):
			print >> sys.stderr, "Aucune entree utilisateur"
			x = 0
		else:
			while True:
				try:
					print >> sys.stderr, "Choix ? ",
					x = int(raw_input())
					break
				except ValueError:
					pass
				except IOError:
					x = 0
					break
				except EOFError:
					x = 0
					break
	print >> sys.stderr, "Choix de " + mystr(res[x])
	return res[x]

