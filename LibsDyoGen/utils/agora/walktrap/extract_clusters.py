import sys
import string

#######################################################################
###############################          ##############################
###############################   Init   ##############################
###############################          ##############################
#######################################################################

filename   = sys.argv[1]
scale = sys.argv[2]
index_filename = sys.argv[3]

#print "\nExtracting clusters at scale " + scale + " from file " + filename + " with index file " + index_filename + "...\n"


class Merge (object):
	scale = 0.0
	father = 0
	children = []

allMerges = []
name = {}


########################################################################
#############################               ############################
#############################   Functions   ############################
#############################               ############################
########################################################################

def findMaximumIndex ():
	max = 0
	for merge in allMerges:
		if int(merge.father) > max:
			max = int(merge.father)
	return int (max)

def printAllChildren (father):
	#print "Calling printAllChildren on father = " + father
	alreadySeen[int(father)] = 1	
	found = 0
	for merge in allMerges:
		if merge.father == father:
			found = 1
			for child in merge.children:
				printAllChildren (child)
	if found == 0:
		print name[father]

#######################################################################
####################                                ###################
####################   Parse the communities file   ###################
####################                                ###################
#######################################################################

file = open (filename, mode="r")

for line in file:
	if line == "\n":
		break
	line = line.rstrip ()
	line = line.replace (" ", "")
	list = line.split (':')
	
	merge = Merge ()
	merge.children = []
	merge.scale = list[0]

	list[1] = list[1].replace ('-->', ':')
	fatherAndChildren = list[1].split (':')
	merge.father = fatherAndChildren[1]

	children = fatherAndChildren[0].split ('+')

	for child in children:
		merge.children.append (child)

	allMerges.append (merge)

allMerges.reverse ()

#######################################################################
##############                                            #############
##############   Get the gene names from the index file   #############
##############                                            #############
#######################################################################

file = open (index_filename, mode="r")

for line in file:
	line = line.rstrip ()
	couple = line.split (' ')
	#print "Adding gene name " + couple[1] + " for index " + couple[0]
	name[couple[0]] = couple[1]

#for merge in allMerges:
#	print "scale = " + merge.scale + " -- father = " + merge.father + " -- "
#	for child in merge.children:
#		print "\tChild " + child

########################################################################
#########################                       ########################
#########################   Find the clusters   ########################
#########################                       ########################
########################################################################

cluster = 1

max = findMaximumIndex ()

alreadySeen = []

for i in xrange (max + 1):
	alreadySeen.append (0)

for merge in allMerges:
	#print "merge.father = " + merge.father
	if merge.scale < scale and alreadySeen[int(merge.father)] == 0:
		print "Cluster number " + str (cluster) + "\n================"
		cluster = cluster + 1
		printAllChildren (merge.father)
		print ""
