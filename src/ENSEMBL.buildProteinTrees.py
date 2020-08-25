#!/usr/bin/env python2
# Copyright (c) 2010 IBENS/Dyogen : Amelie PERES, Matthieu MUFFATO, Alexandra LOUIS, Hugues ROEST CROLLIUS
# Mail : hrc@ens.fr or alouis@biologie.ens.fr
# Licences GLP v3 and CeCILL v2

__doc__ = """
	Corrige les arbres d'Ensembl en fonction du seuil minimal de duplication_score et de l'arbre des especes desire
		1: score par defaut (0 -> 1)
		2: coef multiplicateur d'un score reference 6X_species / all_species
		3: duplication_confidence_score calcule sur uniquement 6X_species
"""


import sys
import itertools
import collections

import utils.myFile
import utils.myTools
import utils.myPhylTree
import utils.myProteinTree

arguments = utils.myTools.checkArgs( \
	[("phylTree.conf",file), ("ensemblTree",file)], \
	[("cutoff",str,"-1"), ("defaultFamName",str,"GCUSGT%08d"), ("scoreMethod",int,[1,2,3]), ("newNodeID",int,100000000), ("recurs",bool,False)], \
	__doc__ \
)

phylTree = utils.myPhylTree.PhylogeneticTree(arguments["phylTree.conf"])
sys.setrecursionlimit(20000)

# Limites automatiques de score de duplication
if arguments["scoreMethod"] in [1, 3]:
	def calc(anc, val):
		return val
elif arguments["scoreMethod"] == 2:
	def calc(anc, val):
		nesp = len(phylTree.species[anc])
		n2X = len(phylTree.lstEsp2X.intersection(phylTree.species[anc]))
		# La moitie des especes non 2X a vu la duplication (au minimum 1 espece)
		return round(max(1., val*(nesp-n2X)) / nesp, 3) - 2e-3

minDuplicationScore = {}
try:
	# Une limite pour tout le monde
	val = float(arguments["cutoff"])
	for anc in phylTree.listAncestr:
		minDuplicationScore[anc] = calc(anc, val)
except ValueError:
	f = utils.myFile.openFile(arguments["cutoff"], "r")
	for l in f:
		t = l.split()
		anc = phylTree.officialName[t[0]]
		minDuplicationScore[anc] = calc(anc, float(t[1]))
	f.close()
#print >> sys.stderr, "minDuplicationScore:", minDuplicationScore

# Les scores dans l'abre pour les especes modernes valent toujours 1, on doit toujours les accepter
for esp in phylTree.listSpecies:
	minDuplicationScore[esp] = 0
#print >> sys.stderr, "minDuplicationScore:", minDuplicationScore

utils.myProteinTree.nextNodeID = arguments["newNodeID"]

@utils.myTools.memoize
def goodSpecies(anc):
	return phylTree.species[anc].difference(phylTree.lstEsp2X)

def alwaysTrue(tree, rnode):
	return True



def calculScoreConfidence(treeData, treeInfo, rnode):
       # print >> sys.stderr, 'CALCULSCORE'
	@utils.myTools.memoize
	def getSpeciesSets(node):
		if node in treeData:
                        #print "ESSAIIIIIIIIIII", treeData
			return set().union(*(getSpeciesSets(x) for (x,_) in treeData[node]))
		else:
			#print >> sys.stderr, "aaaaaaaaa", treeInfo[node]["taxon_name"]
			assert treeInfo[node]["taxon_name"] in phylTree.listSpecies
			return set([treeInfo[node]["taxon_name"]])

        #print >> sys.stderr, "AOAOA", rnode
	#if rnode not in tree.data:
		#return False

	speciessets = [getSpeciesSets(x) for (x,_) in treeData[rnode]]
        #print >> sys.stderr, "species", speciessets
	inters = set()
        #print >> sys.stderr, "in", inters


	for (s1,s2) in itertools.combinations(speciessets, 2):
                #print >> sys.stderr, "YAYAYAYAYAYAYAYAYAYAYYAYAYAYAYAYA", s1, s2
		inters.update(s1.intersection(s2))

        return



        


def hasLowScore(tree, rnode):

	@utils.myTools.memoize
	def getSpeciesSets(node):
		if node in tree.data:
                        #print "ESSAIIIIIIIIIII", tree.data
			return set().union(*(getSpeciesSets(x) for (x,_) in tree.data[node]))
		else:
			#print >> sys.stderr, "aaaaaaaaa", tree.info[node]["taxon_name"]
			assert tree.info[node]["taxon_name"] in phylTree.listSpecies
			return set([tree.info[node]["taxon_name"]])

        #print >> sys.stderr, "AOAOA", rnode
	if rnode not in tree.data:
		return False

	speciessets = [getSpeciesSets(x) for (x,_) in tree.data[rnode]]
      #  print >> sys.stderr, "species", speciessets
	inters = set()
       # print >> sys.stderr, "in", inters


	for (s1,s2) in itertools.combinations(speciessets, 2):
               # print >> sys.stderr, "spe", s1, s2
		inters.update(s1.intersection(s2))
               
               # print >> sys.stderr, "end", inters


        for caractere in inf:
            if caractere == "taxon_name":
                ancestor = inf[caractere]  
               # print >> sys.stderr, 'ancestor', ancestor
                compt = 0
                number_gene_in_each_group = {}

                for (fils,dist) in phylTree.items[ancestor]:
                    Allfils = []
                    #print >> sys.stderr, 'fils', fils
                    for desc in phylTree.allDescendants[fils]:
                        if desc in phylTree.listSpecies:
                            Allfils.append(desc)
                    
                   # print >> sys.stderr, 'allllfils', Allfils
                   
                    number_gene = 0  
                    
                   # print >> sys.stderr, 'inters', inters
                   # print >> sys.stderr, 'Allfils', Allfils
                    for bouh in inters:
                       # print >> sys.stderr, 'ininters', "/" + bouh + "/"

                        if bouh in Allfils:
                            
                            #print >> sys.stderr, 'bouh', bouh
                            number_gene +=1
                   # print >> sys.stderr, 'RSE', number_gene     
                    number_gene_in_each_group[fils] = number_gene

                #print >> sys.stderr, 'FINAL', number_gene_in_each_group
                
                global passage
                if passage == 1:
                    groupok = 0
                    for group in number_gene_in_each_group:
                        if number_gene_in_each_group[group] != 0:
                            groupok = 1
                           # print >> sys.stderr, 'GROUPPPE', group, node

                    
                    if groupok == 1:
                        for groupa in number_gene_in_each_group:
                            if number_gene_in_each_group[groupa] == 0:
                               # print  >> sys.stderr, 'NONNNNNNNNAANNNNNN'
                                global MauvaiseDup
                      
                                MauvaiseDup.append((node,groupa))
                           # print >> sys.stderr, 'mauvaisedup444444', MauvaiseDup
                   
                                   # inters = []

              
                else:
                    for group in number_gene_in_each_group:
                    #if number_gene_in_each_group[group] != 0:
                           # groupok = 1
                            

                    #print >> sys.stderr, 'GROUPPPE', group, node
                    #if groupok == 1:
                        if number_gene_in_each_group[group] == 0:
                            #print  >> sys.stderr, 'NONNNNNNNNNNNNNN'
                            global MauvaiseDup
                      
                            MauvaiseDup.append((node,group))
                                #print >> sys.stderr, 'mauvaisedup444444', MauvaiseDup
                   
                           # inters = []            







                        #for essa in Allfils:
                            #print >> sys.stderr, 'inallfils', "/" + essa + "/"
                            #if bouh in s1:
                              #  nuber_gene_in_first_group +=1


                       # for bouh in s2:
                            
                           # print >> sys.stderr, 'bouh', bouh
                            #if bouh in Allfils:
                               # nuber_gene_in_second_group +=1
                    #print >> sys.stderr, 'RES', nuber_gene_in_first_group, nuber_gene_in_second_group         
                     

                        #if bouh not in Allfils:
                          #  print >> sys.stderr, 'bouh not fils'
                        



	all = set().union(*speciessets)
	anc = tree.info[rnode]["taxon_name"]
       # print >> sys.stderr, 'anc', anc
	if arguments["scoreMethod"] == 3:
		inters.intersection_update(goodSpecies(anc))
		all.intersection_update(goodSpecies(anc))
	#print >> sys.stderr,rnode
       # print >> sys.stderr, "RESULT", (len(inters) == 0)
       # print >> sys.stderr, "score", len(inters), minDuplicationScore[anc], minDuplicationScore[anc] * len(all)
       # if (len(inters) == 0) and (minDuplicationScore[anc] == 0):
               # print >> sys.stderr, 'YESSSSSSSSSSSS'
       # if (len(inters) < (minDuplicationScore[anc] * len(all))):
               # print >> sys.stderr, 'NOOOOOOOOOOOOOOO'
	return ((len(inters) == 0) and (minDuplicationScore[anc] == 0)) or (len(inters) < (minDuplicationScore[anc] * len(all)))


def testDuplicationInLaterNode(tree, rnode):
    return








nbEdit = {"dubious": 0, "toolow": 0, "good": 0}

nbdup = 0

for (nb,tree) in enumerate(utils.myProteinTree.loadTree(arguments["ensemblTree"])):
        treename = ''
        nbEditspe = {"dubious": 0, "toolow": 0, "good": 0}
       # print >> sys.stderr, 'SPECIES', phylTree.listSpecies
        NodeID = 500000000000000
        MauvaiseDup = []
        MauvaiseDupBis = []
       # print >> sys.stderr, '\n\n\n\n\nPASSAGE\n\n\n'
       # print >> sys.stderr, nb, tree
        Taille = len(tree.info)
        #print >> sys.stderr, Taille


        

	assert max(tree.info) < arguments["newNodeID"]

        passage = 0
        nbdupspe = 0
        boucle = 0
        while boucle == 0:
                
                passage +=1
                boucle = 1
                duptocomplete = []

                Ancdata = {}
                for i in tree.data:        
                        Ancdata[i] = tree.data[i]
                NodeModif = []
               # # >> sys.stderr, '\n\n\n\nDATAAAAA', tree.data
	        # On trie les bonnes duplications des mauvaises
	        ################################################
                nodeLowscore = []
	        for (node,inf) in tree.info.iteritems():
                        if 'tree_name' in inf:
                                treename = inf['tree_name']
                        
		       # print >> sys.stderr,"hhh", node,inf

                        if inf['taxon_name'] in phylTree.listSpecies and inf['Duplication'] >= 2:
                                #on s'arrete aux feuilles
                               # print >> sys.stderr,"EROOOOOOOOR4"
                                if passage == 1:
                                        nbEdit["good"] += 1
                                        nbEditspe["good"] += 1
                                        nbdup +=1
                                        nbdupspe +=1
                                
                                
                        else:
		                if inf['Duplication'] != 0:
                                        if passage == 1:
                                               # print >> sys.stderr, 'ouii'
                                                nbdup +=1
                                                nbdupspe +=1
		
			                if 'dubious_duplication' in inf:
				                # On considere que les duplications 'dubious' ne sont pas valables pour une duplication
				                assert inf['Duplication'] == 1
				                del inf['dubious_duplication']
                                                if passage == 1:
				                        nbEdit["dubious"] += 1
                                                        nbEditspe["dubious"] += 1
                                               # print >> sys.stderr, "DUBIOUS", node,inf
                                                boucle = 0 

                  

			                if hasLowScore(tree, node):
                                                NodeModif.append(node)
				                inf['Duplication'] = 1
                                                if passage == 1:
				                        nbEdit["toolow"] += 1
                                                        nbEditspe["toolow"] += 1
                                                #print >> sys.stderr, "HasLowScore", node,inf
                                                #testDuplicationInLaterNode
                                                boucle = 0
                                                nodeLowscore.append(node)
                                                #for ba, bo in MauvaiseDup:
                                                   # if ba == node:
                                                       # print >> sys.stderr, "/n/nYESSSSSSSSSSSSSSSS"
                                                       # MauvaiseDupBis.append((ba,bo))

			                else:

				                # Il faut la passer a 2 si le score est suffisant
				                # Attention: pour les arbres d'Ensembl dont la racine est une duplication, celle-ci vaut 1 (parce qu'elle n'a pas d'outgroup)
				                if inf['Duplication'] == 1:
					                inf['Duplication'] = 3
                                                       # print >> sys.stderr,"goood", node, inf
				                else:
					                assert inf['Duplication'] in [2,3]
                                                       # print >> sys.stderr,"other", node, inf
                                                if passage == 1:
				                        nbEdit["good"] += 1
                                                        nbEditspe["good"] += 1
                
               # for (node,inf) in tree.info.iteritems():
                       # print >> sys.stderr, "End", node,inf
	        tree.flattenTree(phylTree, True)
               # for (node,inf) in tree.info.iteritems():
                       # print >> sys.stderr, "flatten", node,inf

                NodeToTest = {}
                for noeud in NodeModif:
                       # print >> sys.stderr, 'ok', noeud
                        if noeud in tree.data:
                                if noeud in nodeLowscore:
                                    NodeToTest[noeud]=[]
              #  print >> sys.stderr, 'Mauvaisedup', MauvaiseDup

                for duo in MauvaiseDup:
                         if duo[0] in NodeToTest:
                                 NodeToTest[duo[0]].append(duo[1])    
              #  print >> sys.stderr, 'Nodetotest', NodeToTest






                NewDataTreeInfo = {}
                for iii in tree.info:
                        NewDataTreeInfo[iii] = tree.info[iii]
               # print >> sys.stderr, 'TREEEEEEEEEEEEEEEE', NewDataTreeInfo
	        for (node,inf) in tree.info.iteritems():
                        if node in NodeToTest:
                                ChildToRemove = []
                               # print >> sys.stderr, 'nodeatester', node
                                

                              #  print >> sys.stderr, 'inf', inf

                              #  print >> sys.stderr, tree.data[node]

                                for aaa in tree.data[node]:
                                      #  print >> sys.stderr, 'aaaa', aaa[0], NewDataTreeInfo[aaa[0]]
                                        for informa in NewDataTreeInfo[aaa[0]]:
                                                if informa == 'taxon_name':
        #                                               # print >> sys.stderr, 'informa', tree.info[aaa[0]][informa]  
                                        #print >> sys.stderr, 'taxonname', tree.info[aaa[0][taxon_name]]
                                                        
                                                        for ancetre in NodeToTest[node]:
                                                                
                                                                #print >> sys.stderr, 'ancetre', ancetre

                                                                if phylTree.isChildOf(NewDataTreeInfo[aaa[0]][informa], ancetre):
                                                                        #print >> sys.stderr, 'ouiiiiiiiiiiiiiiiiiiii', node, ancetre
                                                                        ChildToRemove.append(aaa)



                               # if len(tree.data[node]) == len(ChildToRemove):
                                       # print >> sys.stderr, 'DONADA'
                                        #del NodeToTest[node]

 

                                #if (len(tree.data[node]) + len(ChildToRemove)) >= 2:
                                if (len(tree.data[node]) != len(ChildToRemove)):
                                        NodeID +=1
                                       # print >> sys.stderr,'iciiiiiiiiiiii', node, ChildToRemove
                                        datoo = tree.data[node]
                                       # print >> sys.stderr, 'datoo', datoo
                                       # print >> sys.stderr, 'tree.data[node]', tree.data[node]

                                       # print >> sys.stderr, 'ancdatapo', Ancdata[node]
                                        ancnodetotest = []
                                        newnodedata1 = []
                                        ChildOk = []
                                        ChildeNewnode = []
                                       # print >> sys.stderr, 'treedatapour18', tree.data[node]
                                       # print >> sys.stderr, 'Acndatapour18', Ancdata[node]
                                       # AncnodeTree = []
                                      #  AncnodeTreeNode = []
                                       # for ancnode in Ancdata[node]:
                                            #    AncnodeTree.append(ancnode)
                                            #    AncnodeTreeNode.append(ancnode[0])
                                       # for treenode in tree.data[node]:
                                             #   if treenode not in AncnodeTreeNode:
                                                  #      AncnodeTree.append(treenode)
#
                                     #   print >> sys.stderr, 'anctrenode', AncnodeTree
                                                       
                                        for ancnode in Ancdata[node]:


                                               # print >> sys.stderr, 'ancnode18', ancnode
                                                if ancnode in tree.data[node]:
                                                       # print >> sys.stderr, 'yesss', ancnode
                                                        newnodedata1.append(ancnode)
                                                        ChildOk.append(ancnode[0])
                                                else:
                                                        ancnodetotest.append(ancnode)


                                        listok = []
                                        for jhdi in tree.data[node]:
                                                listok.append(jhdi[0])


                                        if ancnodetotest != []:
                                          
                                               # print >> sys.stderr, 'node a creer', ancnodetotest
                                                for noooode in ancnodetotest:
                                                        recur = 0
                                                        NodeAA = []
                                                       # print >> sys.stderr, 'ICICICICIC', Ancdata[noooode[0]]
          
                                                        for infoancnode in Ancdata[noooode[0]]:
                                                                NodeAA.append(infoancnode)
                                                        


                                                      #  print >> sys.stderr, '1111111111111111', NodeAA
                                                        while recur == 0:
                                                       
                                                                recur = 1

                                                                #print >> sys.stderr, 'comprendre', NodeAA
                                                                Remove = []
                                                                Add = []
                                                                for nodeanc in NodeAA:
                                                                       #  print >> sys.stderr, 'NODEEE', nodeanc
                                                                         
                                                                         if nodeanc[0] in listok or nodeanc[0] not in Ancdata:
                                                                                 newnodedata1.append((nodeanc[0], nodeanc[1]))
                                                                                 
                                                                                 ChildeNewnode.append(nodeanc[0])
                                                                                 Remove.append(nodeanc)
                                                                                # print >> sys.stderr, 'Remove', Remove
                                                                                # print >> sys.stderr, '22222222222222222222', NodeAA
                                                                         else:
                                                                                 recur = 0
                                                                   #              print >> sys.stderr, 'bug', nodeanc
                                                                                 Remove.append(nodeanc)
                                                                                # print >> sys.stderr, 'Remove', Remove
                                                                                # print >> sys.stderr, 'ancdat', Ancdata[nodeanc[0]]
                                                                                 for jesaispas in Ancdata[nodeanc[0]]:
                                                                                         Add.append(jesaispas)
                                                                                        # print >> sys.stderr, 'Add', Add
                                                                                        # print >> sys.stderr, '3333333333', NodeAA
                                                                for chose in Remove:
                                                                        NodeAA.remove(chose)
                                                                for otherchose in Add:
                                                                        NodeAA.append(otherchose)
                                                               # print >> sys.stderr, 'fiiiiin', NodeAA


                                                                               
                                       # print >> sys.stderr, 'fin newnode', newnodedata1                                  
                                       # print >> sys.stderr, 'chile to remove', ChildToRemove
                                       # print >> sys.stderr, 'ChildOk', ChildOk
                                        for chilnoode in ChildToRemove:
                                                if chilnoode[0] in ChildeNewnode:
                                                        ChildeNewnode.remove(chilnoode[0])   


                                       # print >> sys.stderr, 'ChildeNewnode', ChildeNewnode


                                       # print >> sys.stderr, 'tree', tree.data[node]

                                        #print >> sys.stderr, 'datoo', datoo
                                        tree.data[node] = ChildToRemove
                                        tree.data[node].append((NodeID, 0))
                                       # print >> sys.stderr, 'DUUUP', tree.info[node]

                                        

                                       # print >> sys.stderr, len(ChildeNewnode) + len(ChildOk)
                                        if ChildeNewnode == []:
                                                
                                                tree.data[NodeID] = []
                                                for efgh in datoo:
                                                        if efgh[0] in ChildOk:
                                                                tree.data[NodeID].append(efgh)
                                                               # print >> sys.stderr, 'ESS', efgh
                                        elif (len(ChildeNewnode) + len(ChildOk)) <= 2:
                                                tree.data[NodeID] = []
                                                for efgh in datoo:
                                                        if efgh[0] in ChildOk:
                                                                tree.data[NodeID].append(efgh)
                                                               # print >> sys.stderr, 'ESS', efgh

                                                        if efgh[0] in ChildeNewnode:
                                                                tree.data[NodeID].append(efgh)
                                                              #  print >> sys.stderr, 'ESS', efgh
                                                                


                                        else:
                                                tree.data[NodeID] = []
                                                for efgh in datoo:
                                                        if efgh[0] in ChildOk:
                                                                tree.data[NodeID].append(efgh)
                                                              #  print >> sys.stderr, 'ASS', efgh
                                                tree.data[NodeID].append((NodeID +1, 0))
                                                NodeID += 1
                                                tree.data[NodeID] = []
                                                for efgh in datoo:
                                                        if efgh[0] in ChildeNewnode:
                                                                tree.data[NodeID].append(efgh)
                                                               # print >> sys.stderr, 'USS', efgh
                                        

                                      #  print >> sys.stderr, 'treelast', node, tree.data[node]
                                      #  print >> sys.stderr, 'treelast', NodeID, tree.data[NodeID]
                                       # print >> sys.stderr, tree.data









                                        speciessets2 = []

                                        if ChildeNewnode == [] or (len(ChildeNewnode) + len(ChildOk)) <= 2:
                                                NewDataTreeInfo[NodeID] = {}
                                                NewDataTreeInfo[NodeID]['taxon_name'] = phylTree.lastCommonAncestor( [NewDataTreeInfo[g]['taxon_name'] for (g,_) in tree.data[NodeID]] ) 
                                                NewDataTreeInfo[NodeID]['Duplication'] = 2
                                               # print >> sys.stderr, 'yoyoyoyoyoyoyoyoyoYOYO'
                                               # print >> sys.stderr, 'node' ,NodeID
                                              #  print >> sys.stderr, tree.data
                                                #calculScoreConfidence(tree.data, tree.info, NodeID)
	                                       # for (x,_) in tree.data[NodeID]:
                                                       # if x in tree.data: 
                                                               # print >> sys.stderr, 'ESSAIIIIIII', x
                                                               # for (y,_) in tree.data[x]:
                                                                        #speciessets = []
                                                                        #print >> sys.stderr,'ooo', y
                                                
                                                #speciessets2 = [getSpeciesSets2(x) for (x,_) in tree.data[NodeID]]
                                                #print >> sys.stderr, 'ICIIOLOLOLOL', speciessets2 
                                                duptocomplete.append(NodeID)
                                                NewDataTreeInfo[NodeID]['duplication_confidence_score'] = 0

                                        else:
                                                NewDataTreeInfo[NodeID] = {}
                                                NewDataTreeInfo[NodeID - 1] = {}
                                                NewDataTreeInfo[NodeID]['taxon_name'] = phylTree.lastCommonAncestor( [NewDataTreeInfo[g]['taxon_name'] for (g,_) in tree.data[NodeID]] ) 
                                                NewDataTreeInfo[NodeID]['Duplication'] = 0
                                                NewDataTreeInfo[NodeID - 1]['taxon_name'] = phylTree.lastCommonAncestor( [NewDataTreeInfo[g]['taxon_name'] for (g,_) in tree.data[NodeID - 1]] )  
                                                NewDataTreeInfo[NodeID - 1]['Duplication'] = 2
                                               # print >> sys.stderr, 'YOYOyoyoyoyoyoyoyo'
                                               # print >> sys.stderr, 'node' ,NodeID - 1
                                               # print >> sys.stderr, tree.data
                                               
                                                #calculScoreConfidence(tree.data, tree.info, NodeID - 1)
	                                        #for (x,_) in tree.data[NodeID]:
                                                        #if x in tree.data: 
                                                               # print >> sys.stderr, 'ESSAIIIIIII', x
                                                                #for (y,_) in tree.data[x]:
                                                                        #speciessets = []
                                                                        #print >> sys.stderr,'ooo', y
                                                #speciessets2 = [getSpeciesSets2(x) for (x,_) in tree.data[NodeID]]
                                                #print >> sys.stderr, 'ICIIOLOLOLOL', speciessets2 
                                                duptocomplete.append(NodeID - 1)
                                                NewDataTreeInfo[NodeID - 1]['duplication_confidence_score'] = 0


                                        NewDataTreeInfo[node]['taxon_name'] = phylTree.lastCommonAncestor( [NewDataTreeInfo[g]['taxon_name'] for (g,_) in tree.data[node]] ) 
                                        NewDataTreeInfo[node]['Duplication'] = 0
                                        

                     


	#@utils.myTools.memoize
	#def getSpeciesSets(node):
#		if node in treeData:
#                        print "ESSAIIIIIIIIIII", treeData#
#			return set().union(*(getSpeciesSets(x) for (x,_) in treeData[node]))
#		else:
#			print >> sys.stderr, "aaaaaaaaa", treeInfo[node]["taxon_name"]
#			assert treeInfo[node]["taxon_name"] in phylTree.listSpecies
#			return set([treeInfo[node]["taxon_name"]])
#
 #       print >> sys.stderr, "AOAOA", rnode
	#if rnode not in tree.data:
		#return False

#	speciessets = [getSpeciesSets(x) for (x,_) in treeData[rnode]]
 #       print >> sys.stderr, "species", speciessets
#	inters = set()
 #       print >> sys.stderr, "in", inters


#	for (s1,s2) in itertools.combinations(speciessets, 2):
 #               print >> sys.stderr, "YAYAYAYAYAYAYAYAYAYAYYAYAYAYAYAYA", s1, s2
#		inters.update(s1.intersection(s2))











                                #print >> sys.stderr, 'NodeToTest.remove', NodeToTest







                          
                #for i in MauvaiseDup:
                       # print >> sys.stderr, 'mauvaise dup', i
                MauvaiseDup = []
                #print >> sys.stderr, 'mauvaise dup2222', MauvaiseDup
                tree.info = NewDataTreeInfo
               # print >> sys.stderr, 'EEEENDDDYOUYOU', tree.info
               # for ijj in tree.data:
                #        print >> sys.stderr, ijj, tree.data[ijj]
               # print >> sys.stderr, Taille
               # boucle = 1
               # print >> sys.stderr, 'DUPTOCOMPLETE', duptocomplete


                CaclulScoreConfiance = {}
                @utils.myTools.memoize
                def getSpeciesSets2(bnode):
	                if bnode in tree.data:
                                #print "ESSAIIIIIIIIIII", tree.data
		                return set().union(*(getSpeciesSets2(x) for (x,_) in tree.data[bnode]))
	                else:
		              #  print >> sys.stderr, "aaaaaaaaa", tree.info[bnode]["taxon_name"]
		                assert tree.info[bnode]["taxon_name"] in phylTree.listSpecies
		                return set([tree.info[bnode]["taxon_name"]])

               # print >> sys.stderr, 'anc tree info', tree.info
                for znode in duptocomplete:
                        speciessets2 = [getSpeciesSets2(x) for (x,_) in tree.data[znode]]
                        #print >> sys.stderr, 'ICIIOLOLOLOL', znode,  speciessets2 
                        if len(speciessets2) != 2 and len(speciessets2) != 1:
                                inutil = 'a'
                        if len(speciessets2) != 1:
                                listesp1 = []
                                listesp2 = []
                                listunion = []
                                listintersection = []
                                listesp1 = list(speciessets2[0])
                                listesp2 = list(speciessets2[1])
                               # print >> sys.stderr, 'listesp1', listesp1
                               # print >> sys.stderr, 'listesp2', listesp2
                                listintersection = list(set(listesp1) & set(listesp2))
                                #print >> sys.stderr, 'listintersection', listintersection
                                listunion = list(set(listesp1) | set(listesp2))
                                #print >> sys.stderr, 'listunion', listunion
                                nb1 = float(len(listintersection))
                                nb2 = float(len(listunion))
                                CalculScoreconfidence2 = nb1/nb2
                              #  print >> sys.stderr, 'CalculScoreconfidence2', CalculScoreconfidence2, nb1, nb2
                        
                                tree.info[znode]['duplication_confidence_score'] = CalculScoreconfidence2


              #  print >> sys.stderr, 'New tree info', tree.info

#        print >> sys.stderr,'arbe', treename, 'edition', nbEditspe, 'duplication',nbdupspe

	tree.rebuildTree(phylTree, hasLowScore if arguments["recurs"] else alwaysTrue)
#        for (node,inf) in tree.info.iteritems():
#	        print >> sys.stderr, "rebuildTree", node,inf
	if "tree_name" not in tree.info[tree.root]:
		tree.info[tree.root]["tree_name"] = arguments["defaultFamName"] % nb

	tree.printTree(sys.stdout)


print >> sys.stderr, 'total edition', nbEdit
print >> sys.stderr, 'total duplication', nbdup

