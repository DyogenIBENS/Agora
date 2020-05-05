#! /usr/bin/env python
# -*- coding: utf-8 -*-
# Agora v1.00
# python 2.7
# Copyright © 2016 IBENS/Dyogen : Matthieu MUFFATO, Alexandra Louis, Nga Thi tuy Nguyen, Joseph Lucas, Hugues ROEST CROLLIUS
# mail : hrc@ens.fr or alouis@biologie.ens.fr
# This is free software; you may copy, modify and/or distribute this work under the terms of the GNU General Public License, version 3 or later and the CeCiLL v2 license in France

__doc__ = """
        Convert gene trees from NHX to my protTree forest
"""

import cStringIO
import string

import utils.myFile as myFile
import utils.myPhylTree as myPhylTree
import utils.myTools as myTools

nodeid = 0
ntree = 0

def processData(data):

    tree = myPhylTree.PhylogeneticTree(cStringIO.StringIO(data))

    def printTree(indent, node):
        global nodeid, ntree
        print "%sid\t%d" % (indent, nodeid)
        nodeid += 1

        info={}


        if "B" in tree.info[node]:
            #info["Bootstrap"]=int(tree.info[node]["B"])
            info["Bootstrap"] = tree.info[node]["B"]

        if "D" in tree.info[node] and tree.info[node]["D"]=="N":
            info["Duplication"]=0
            if "E" in tree.info[node]:
                info["taxon_lost"]=tree.info[node]["E"].split("=-$")[1].split("-")
                if "S" in tree.info[node]:
                    info["taxon_name"]=string.capitalize(tree.info[node]["S"].replace("_"," ").replace("."," "))
            else:
                if "S" in tree.info[node] :
                    info["taxon_name"]=string.capitalize(tree.info[node]["S"].replace("_"," ").replace("."," "))



        elif "D" in tree.info[node] and tree.info[node]["D"]=="Y":

            if "DD" in tree.info[node] and tree.info[node]["DD"]=="Y" :
                info["Duplication"]=1
                info["dubious_duplication"]=1

                if "E" in tree.info[node] :
                    info["taxon_lost"]=tree.info[node]["E"].split("=$-")[1].split("-")

                    if "S" in tree.info[node] :
                        info["taxon_name"]=string.capitalize(tree.info[node]["S"].replace("_"," ").replace("."," "))

                        if "SIS" in tree.info[node] :
                            info["duplication_confidence_score"]=float(tree.info[node]["SIS"])/100
                else :

                    if "S" in tree.info[node]:
                        info["taxon_name"]=string.capitalize(tree.info[node]["S"].replace("_"," ").replace("."," "))

                        if "SIS" in tree.info[node]:
                            info["duplication_confidence_score"]=float(tree.info[node]["SIS"])/100
            else :
                info["Duplication"]=2
                if "E" in tree.info[node] :
                    info["taxon_lost"]=tree.info[node]["E"].split("=$-")[1].split("-")

                    if "S" in tree.info[node]:
                        info["taxon_name"]=string.capitalize(tree.info[node]["S"].replace("_"," ").replace("."," "))

                        if "SIS" in tree.info[node] :
                            info["duplication_confidence_score"]=float(tree.info[node]["SIS"])/100
                else :
                    if "S" in tree.info[node]:
                        info["taxon_name"]=string.capitalize(tree.info[node]["S"].replace("_"," ").replace("."," "))

                        if "SIS" in tree.info[node] :
                            info["duplication_confidence_score"]=float(tree.info[node]["SIS"])/100

        elif "E" in tree.info[node] :
            info["taxon_lost"]=tree.info[node]["E"].split("=$-")[1].split("-")

            if "S" in tree.info[node]:
                info["Duplication"]=0
                info["taxon_name"]=string.capitalize(tree.info[node]["S"].replace("_"," ").replace("."," "))


        elif "S" in tree.info[node]:
            info["Duplication"]=0
            info["taxon_name"]=string.capitalize(tree.info[node]["S"].replace("_"," ").replace("."," "))



        if indent == "":
            info["tree_name"] = "Fam%06d" % ntree
            ntree += 1

        if node not in tree.items:
            # modifier par alex pour garder le "_" dans le nom de gene!
            ############################################################

            #x = node.rpartition("_")[0]
            x=node

            # fin de modif Alex


            info["gene_name"] = x
        else:
            info["node_name"] = node

        print "%sinfo\t%s" % (indent, info)

        if node in tree.items:
            indent = indent + "\t"
            for (e,l) in tree.items[node]:
                print "%slen\t%g" % (indent,l)
                printTree(indent, e)

    printTree("", tree.root)


arguments = myTools.checkArgs( [("tree",file)], [], __doc__)

f = myFile.openFile(arguments["tree"], "r")
for line in f:
    if len(line.replace(" ","").replace("\n","")) == 0:
        #Do nothing : empty line
        continue
    elif line.find(";\n"):
        processData(line)
    else:
        raise NameError("The nhx tree is not formated properly. Please take care to have one tree per line. A line ends with \";\\n\"")
#       break
f.close()
