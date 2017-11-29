# -*- coding: utf-8 -*-

# LibsDyogen version 1.0 (6/11/2015)
# python v2.7 at least is needed
# Copyright © 2015 IBENS/Dyogen : Matthieu MUFFATTO, Joseph LUCAS and Hugues ROEST CROLLIUS
# mail : jlucas@ens.fr
# Licences GLP v3 and CeCILL v2

import sys
import collections

import myTools
import myFile

# Manage gene trees in the form of (data,info)

# return the ancestor names since the previousAnc
def getIntermediateAnc(phylTree, previousAnc, lastWrittenAnc, newAnc, isDuplication):

    if lastWrittenAnc:                                                      # cruise mode
        toWrite =list(phylTree.dicLinks[lastWrittenAnc][newAnc][1:])    # genes of the species between lastWritten and newAnc are recorded exepted if it is a duplication node (see after)

    elif not lastWrittenAnc:                                        # in the first species of the gene tree or at the first node outside the first species
        if previousAnc == None:                                 # at the root
            toWrite =list([newAnc])                         # the gene is recorded exepted if it is a duplication (see after)

        elif previousAnc:                                       # not at the root

            if previousAnc == newAnc:                       # still in the first species
                toWrite=list([newAnc])                  # the gene is recorded exepted if it is a duplication (see after)

            elif previousAnc != newAnc:                     # in a child species
                toWrite=list(phylTree.dicLinks[previousAnc][newAnc]) # the genes of the species between previousAnc and newAnc are recorded and newAnc is removed if the node is a duplication node (see after)

    if isDuplication:                                               # if the current node is a duplication newAnc is not recorded in the list to be written
        toWrite.remove(newAnc)

    root=False                                              # root refers to terminal genes of the first species
    if not lastWrittenAnc:                                  # in the first species or at the first node outside the first species
        if not isDuplication:                           # if it is a speciation node
            root = True
        if previousAnc != newAnc and isDuplication:     # if at the first node outside the first species and if this node is a duplication node
            root = True

    if len(toWrite) == 0:
        newLastWritten = lastWrittenAnc
    else:
        newLastWritten = toWrite[-1]

    return (toWrite,newLastWritten,root)

# class for managing a forest of gene trees
class ProteinTree:

    def __init__(self, data=None, info=None, root=None):
        self.data = {} if data is None else data
        self.info = {} if info is None else info
        self.root = root


    def doBackup(self):
        self.backRoot = self.root
        self.backData = dict((node,values[:]) for (node,values) in self.data.iteritems())
        self.backInfo = dict((node,values.copy()) for (node,values) in self.info.iteritems())



    # print the tree into the phylTree format (with tabulations)
    def printTree(self, f, node=None):

        def rec(n, node):

            indent = "\t" * n
            # id of the node
            print >> f, "%sid\t%d" % (indent, node)
            # informations
            print >> f, "%sinfo\t{%s}" % (indent, ", ".join(repr(key) + ": " + repr(value) for (key, value) in sorted(self.info[node].iteritems())))
            # children
            for (g,d) in self.data.get(node,[]):
                print >> f, "%s\tlen\t%f" % (indent, d)
                rec(n+1, g)

        rec(0, self.root if node is None else node)
        try:
            f.flush()
        except AttributeError:
            pass

    # print the tree into the Newick format (with parentheses)
    def printNewick(self, f, root=None, withDist=True, withTags=False, withAncSpeciesNames=False, withAncGenesNames=False):
        NHX = {"Duplication": "D", "Bootstrap": "B", "taxon_name": "S", "duplication_confidence_score": "SIS", "dubious_duplication": "DD"}
        def rec(node):
            if node in self.data:
                return "(" + ",".join(
                        rec(g)
                        + ((":%f" % l) if withDist else "")
                        + ("[&&NHX:" + ":".join(("%s=%s" % ((NHX[tag],self.info[g][tag]) if tag!="Duplication" else (NHX[tag],"N" if self.info[g][tag]== 0 else "Y"))).replace(" ", ".") for tag in NHX if tag in self.info[g]) + "]" if withTags else "")
                        for (g,l) in self.data[node]
                ) + ")" + (self.info[node]["taxon_name"].replace(' ', '.') if withAncSpeciesNames and ("taxon_name" in self.info[node]) else '')+(self.info[node]['family_name'].split("/")[0]if withAncGenesNames and ("taxon_name" in self.info[node]) else '')
            else:
                return self.info[node]['gene_name'].split("/")[0]

        if root is None:
            root = self.root
        print >> f, rec(root) + ("[&&NHX:" + ":".join(("%s=%s" % ((NHX[tag],self.info[root][tag]) if tag!="Duplication" else (NHX[tag],"N" if self.info[root][tag]== 0 else "Y"))).replace(" ", ".") for tag in NHX if tag in self.info[root]) + "]" if withTags else "") + ";"
        try:
            f.flush()
        except AttributeError:
            pass


    #FIXME print tree into the Newick format
    def printNewickTree(self, f, node=None):
        genes = []
        def rec(node):
            if node not in self.data:
                genes.append(self.info[node]['gene_name'])
                return self.info[node]['gene_name']
            else:
                return "(" + ",".join([rec(x) + ":" + str(l) for (x,l)  in self.data[node]]) + ") " + self.info[node]['family_name']
        tr = rec(self.root if node is None else node)
        print >> f, " ".join(genes)
        print >> f, tr, ";"


    # Compact a tree by removing intermediary nodes that have only one child
    def compactTree(self, phylTree, node=None):

        def do(node):
            # end of the process on one leaf
            if node not in self.data:
                return False

            flag = False
            # recursive calls
            for (gg,_) in self.data[node]:
                flag |= do(gg)

            if len(self.data[node]) > 1:
                return flag

            # edition of the current node
            (g,l) = self.data[node][0]
            if g in self.data:
                self.data[node] = [(gg,ll+l) for (gg,ll) in self.data[g]]
                del self.data[g]
            else:
                del self.data[node]
            self.info[node] = self.info[g]
            del self.info[g]
            return True

        return do(self.root if node is None else node)


    # rename all nodes in order to make them match with the common ancestors of their descendants
    def renameTree(self, phylTree, node=None):

        def do(node):
            # end of the process on one leaf
            if node not in self.data:
                return False

            flag = False
            # recursive calls
            for (g,_) in self.data[node]:
                flag |= do(g)

            # rename the current node
            newName = phylTree.lastCommonAncestor( [self.info[g]['taxon_name'] for (g,_) in self.data[node]] )
            flag |= (self.info[node]['taxon_name'] != newName)
            self.info[node]['taxon_name'] = newName

            return flag

        return do(self.root if node is None else node)


    # Flatten a node and direct children if they represent the same taxon and if their is no duplication
    # For instance:
    #  ((Eutheria1,Eutheria2)XA,(Eutheria3,Eutheria4)XB)XC is transformed into (Eutheria1,Eutheria2,Eutheria3,Eutheria4)
    #    only if XA, XB et XC are speciation nodes
    def flattenTree(self, phylTree, rec, node=None):

        def do(node):
            # end of the process on one leaf
            if node not in self.data:
                return False
            assert len(self.data[node]) > 0

            flag = False
            # recursive calls
            if rec:
                for (g,_) in self.data[node]:
                    flag |= do(g)

            self.info[node]['taxon_name'] = phylTree.lastCommonAncestor( [self.info[g]['taxon_name'] for (g,_) in self.data[node]] )

            # if it is a true duplication, their is nothing more to do
            if self.info[node]['Duplication'] >= 2:
                return flag

            newData = []
            taxonName = self.info[node]['taxon_name']
            for (g,d) in self.data[node]:
                inf = self.info[g]
                # 2x the same taxon and no duplication

                if (inf['taxon_name'] == taxonName) and (inf['Duplication'] < 2) and (g in self.data):

                    newData.extend([(g2,d+d2) for (g2,d2) in self.data[g]])
                    del self.data[g]

                    self.info[node].update(self.info[g])
                    del self.info[g]
                    flag = True
                else:
                    newData.append( (g,d) )

            assert len(newData) == len(set(g for (g,_) in newData)), newData
            assert len(newData) > 0, (node,self.data[node],newData)

            self.info[node]['Duplication'] = 0
            self.data[node] = newData

            if len(self.data[node]) > 0:
                assert self.info[node]['taxon_name'] == phylTree.lastCommonAncestor( [self.info[g]['taxon_name'] for (g,_) in self.data[node]] )

            return flag

        return do(self.root if node is None else node)


    # give back the expected topology to the tree (to match the species tree)
    #   gather equivalent nodes under the same child
    def rebuildTree(self, phylTree, hasLowScore, node=None):

        def do(node):

            # end of the process on one leaf
            if node not in self.data:
                return False

            flag = False

            # only if it is not a true duplication the children are changed
            if self.info[node]['Duplication'] < 2:

                # the node will be a speciation exepted special case under
                self.info[node]['Duplication'] = 0

                # the children are redefined for -our- phylogenetic tree by organising the children into packs
                children = collections.defaultdict(list)
                anc = self.info[node]['taxon_name']
                lchildren = phylTree.items.get(anc, [])
                for (g,d) in self.data[node]:
                    gname = self.info[g]['taxon_name']
                    for (a,_) in lchildren:
                        if phylTree.isChildOf(gname, a):
                            children[a].append( (g,d) )
                            break
                    else:
                        # we can be here only if g is in the same ancestor than node and if g is a duplication,
                        # this usually entails Duplication >= 2, exepted for the new created nodes
                        assert (gname == anc), "ERROR: name!=anc [%s / %s / %s]" % (node, anc, gname)
                        # false: g is not always a duplication !
                        # assert self.info[g]['Duplication'] >= 2
                        # the current node will thus be a duplication node
                        self.info[node]['Duplication'] = 3
                        children[anc].append( (g,d) )

                #print >> sys.stderr, "NEW CHILD", child
                # len(child):
                #  1 -> only anc
                #  2 or 3 among C1/C2/anc
                assert (len(children) != 1) or (anc in children), "ERROR: 1=anc [%s / %s]" % (node, children)
                assert (len(children) <= (1+len(lchildren))), "ERROR: len>(1+nbChildren) [%s / %s]" % (node, children)

                todo = []

                if len(children) > 1:
                    if anc in children:
                        lst1 = children.pop(anc)
                        lst2 = []
                        for tmp in children.itervalues():
                            lst2.extend(tmp)
                        items = [(anc,lst1), (anc,lst2)]
                    else:
                        items = children.items()

                    newData = set()
                    for (anc,lst) in items:
                        if len(lst) == 1:
                            newData.add( lst[0] )
                        elif len(lst) > 1:
                            for (g,l) in self.data[node]:
                                if (g in self.data) and (self.data[g] == lst):
                                    newData.add( (g,l) )
                                    break
                                if g in self.data:
                                    assert sorted(self.data[g]) != sorted(lst)
                            else:
                                global nextNodeID
                                nextNodeID += 1
                                length = min([d for (_,d) in lst]) / 2
                                self.data[nextNodeID] = [(g,d-length) for (g,d) in lst]
                                anc = phylTree.lastCommonAncestor([self.info[g]['taxon_name'] for (g,_) in lst])
                                self.info[nextNodeID] = {'taxon_name':anc}
                                self.info[nextNodeID]["Duplication"] = 1 if hasLowScore(self, nextNodeID) else 3
                                todo.append(nextNodeID)
                                newData.add( (nextNodeID,length) )
                                self.flattenTree(phylTree, False,  nextNodeID)
                                flag = True
                    assert len(newData) == len(set(g for (g,_) in newData)), newData
                    self.data[node] = [x for x in self.data[node] if x in newData] + list(newData.difference(self.data[node]))
                    for x in todo:
                        if hasLowScore(self, x):
                            self.info[x]["Duplication"] = 0
                            self.flattenTree(phylTree, False, x)
            # recursive calls
            for (g,_) in self.data[node]:
                flag |= do(g)
            return flag
        return do(self.root if node is None else node)

nextNodeID = -1

@myTools.deprecated
def printTree(ft, data, info, root):
    ProteinTree(data, info, root).printTree(ft)



# return the suffix associated to a duplication rank ('a' -> 'z', 'aa' -> 'az', 'ba' ...)
def getDupSuffix(n, upper=False):
    base = 64 if upper else 96
    assert 1 <= n
    s = "."
    while n > 26:
        s = s + chr(base + (n - 1) % 26)
        n = 1 + (n - 1)/26
    return s + chr(base + n)


# load the tree from a file
def loadTree(name):

    ns = myTools.Namespace()

    # read the next line of the file (and bufferise the next one)
    def nextLine():
        old = ns.curr
        try:
            l = ""
            while (l == "") or l.startswith("#"):
                # the final '\n' is removed and we cut owing to the '\t'
                l = f.next().replace('\n', '')
            l = l.split('\t')
            # the triplet (indentation,key,value) is recorded
            ns.curr = (len(l)-2, l[-2], l[-1])
        except StopIteration:
            ns.curr = None
        return old

    # the analysing process of the lines of the file
    def recLoad(tree, indent):

        # id of the point
        currID = int(nextLine()[2])
        # associated informations
        tree.info[currID] = eval(nextLine()[2])

        # children ?
        child = []
        while (ns.curr != None) and (ns.curr[0] == indent+1):
            length = float(nextLine()[2])
            child.append((recLoad(tree, indent+1), length))
        if len(child) > 0:
            tree.data[currID] = child

        return currID


    print >> sys.stderr, "Loading the forest of gene trees %s ..." % name,
    f = myFile.openFile(name, "r") if isinstance(name, str) else name
    ns.curr = None
    nextLine()
    n = (0,0,0)
    while True:
        tree = ProteinTree()
        tree.root = recLoad(tree, 0)
        yield tree
        n = (n[0]+1, n[1]+len(tree.data), n[2]+len(tree.info)-len(tree.data))
        if ns.curr == None:
            break
    print >> sys.stderr, "%d roots, %d branches, %d nodes OK" % n

    f.close()
