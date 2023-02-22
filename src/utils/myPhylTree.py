# -*- coding: utf-8 -*-

# LibsDyogen version 1.0 (6/11/2015) -- modified for AGORA v3.1
# python 3.5 at least is needed
# Copyright © 2006-2022 IBENS/Dyogen, 2020-2021 EMBL-European Bioinformatics Institute, 2021-2022 Genome Research Ltd : Matthieu MUFFATO, Joseph LUCAS and Hugues ROEST CROLLIUS
# mail : agora@bio.ens.psl.eu
# Licences GLP v3 and CeCILL v2

import os
import sys
import itertools
import collections

from . import myFile
from . import newick

SYMBOL6X = '.'
SYMBOL2X = '*'
# global counter
nodeIndex = 0

GeneSpeciesPosition = collections.namedtuple("GeneSpeciesPosition", ['species', 'chromosome', 'index'])

# class mother of all types of trees
class PhylogeneticTree:

    ParentItem = collections.namedtuple("ParentItem", ['name', 'distance'])

    def reinitTree(self, stream=open(os.devnull, 'w')):
        # object initialisations
        # self.parent = {..., son: parentItem, ...}
        self.parent = self.newCommonNamesMapperInstance()
        # self.species = {..., nodeI: [extant nodes that are descendant from nodeI], ...}
        # interior nodes are not included in the list of descendants
        self.species = self.newCommonNamesMapperInstance()
        # self.fileName = {..., node: "fileName", ...} with 'fileName' the name
        # given in genes or ancGenes files (eg Homo.sapiens for node 'Homo
        # sapiens'
        self.fileName = self.newCommonNamesMapperInstance()
        # self.dicLinks = {..., node1: {..., node2: [list of nodes on the
        # evolutive path between node1 and node 2], ...}, ...}
        # Usage:
        # self.dicLinks[node1][node2] = [list of node on the
        # evolutive path between node1 and node 2]
        self.dicLinks = self.newCommonNamesMapperInstance()
        # self.dicParents = {..., node1: {..., node2: [list of nodes that are
        # descendants of the LCA(node1, node2)], ...}, ...}
        self.dicParents = self.newCommonNamesMapperInstance()
        # self.allDescendants = {..., nodeI: [nodes that are descendant from nodeI], ...}
        # interior nodes are included in the list of descendants
        self.allDescendants = self.newCommonNamesMapperInstance()
        self.tmpS = []
        self.tmpA = []

        # analysing process of the tree
        def recInitialize(node, stream=open(os.devnull, 'w')):
            print(".", file=stream)
            self.dicLinks.setdefault(node, self.newCommonNamesMapperInstance())
            # each link between a node and itself is a list containing node
            self.dicLinks.get(node).setdefault(node, [node])
            self.dicParents.setdefault(node, self.newCommonNamesMapperInstance())
            self.dicParents.get(node).setdefault(node, node)
            family = [node]
            if node in self.items:
                # interior node (not a leaf of the tree)
                self.fileName.setdefault(node, str(node))
                s = []
                # list of subFamilies
                lsf = []
                self.tmpA.append(node)
                for (son, bLength) in self.items.get(node):
                    self.parent.setdefault(son, PhylogeneticTree.ParentItem(node, bLength))
                    subFamily = recInitialize(son, stream)
                    s.extend(self.species.get(son))
                    # we go back up
                    family.extend(subFamily)
                    lsf.append(subFamily)
                    for e in subFamily:
                        self.dicParents.get(e).setdefault(node, node)
                        self.dicParents.get(node).setdefault(e, node)
                        self.dicLinks.get(e).setdefault(node, self.dicLinks.get(e).get(son) + [node])
                        self.dicLinks.get(node).setdefault(e, [node] + self.dicLinks.get(son).get(e))
                self.species.setdefault(node, frozenset(s))
                # phylogenetic relationships
                for (sf1, sf2) in itertools.combinations(lsf, 2):
                    for e1 in sf1:
                        for e2 in sf2:
                            self.dicParents.get(e1).setdefault(e2, node)
                            self.dicParents.get(e2).setdefault(e1, node)
                            self.dicLinks.get(e1).setdefault(e2, self.dicLinks.get(e1).get(node) + self.dicLinks.get(node).get(e2)[1:])
                            self.dicLinks.get(e2).setdefault(e1, self.dicLinks.get(e2).get(node) + self.dicLinks.get(node).get(e1)[1:])
            else:
                # leaf of the tree
                if ' ' in node:
                    self.fileName.setdefault(node, str(node).replace(' ', '.'))
                self.fileName.setdefault(node, str(node))
                self.species.setdefault(node, frozenset([node]))
                self.tmpS.append(node)
            self.allDescendants.setdefault(node, frozenset(family))
            return family

        print("Data analysis ", end=' ', file=stream)
        # filling
        recInitialize(self.root, stream=stream)
        # list of extant species
        self.listSpecies = frozenset(self.species).difference(self.items)
        assert frozenset(self.tmpS) == self.listSpecies
        # list of ancestral species
        self.listAncestr = frozenset(self.items)
        assert frozenset(self.tmpA) == self.listAncestr

        # post-analysis
        self.outgroupSpecies = self.newCommonNamesMapperInstance()
        # FIXME: bad name, this is a dict that take a name as keys and return an
        # int that corresponds to the node/leaf ID
        self.indNames = self.newCommonNamesMapperInstance()
        # This is a dict that takes a name as key and return the
        # branch ID.
        # setdefault: This method returns the key value available in the
        # dictionary and if given key is not available then it will return
        # provided default value.
        self.indBranches = self.newCommonNamesMapperInstance()
        # This loop defines the order of iterations in the tree
        self.allNames = self.tmpA + self.tmpS
        for (i, e) in enumerate(self.allNames):
            self.outgroupSpecies.setdefault(e, self.listSpecies.difference(self.species.get(e)))
            self.indNames.setdefault(e, i)
            if e == self.root:
                assert i == 0
            elif e != self.root:
                # each branch correspond to the son species, the root species
                # has no associated branch
                self.indBranches.setdefault(e, i-1)
            self.officialName[i] = e
        # (e)numerateItems, all nodes excepted the root using indNames
        self.numItems = [[]] * len(self.allNames)
        for a in self.listAncestr:
            self.numItems[self.indNames.get(a)] = [(self.indNames.get(son), l) for (son, l) in self.items.get(a)]
        # (e)numerateParents, all nodes excepted the leaves using indNames
        self.numParent = [None] * len(self.allNames)
        for (e, (par, l)) in self.parent.items():
            self.numParent[self.indNames.get(e)] = (self.indNames.get(par), l)
        # self.numParent[self.indNames.get('self.root')] = None
        assert set(self.indBranches.values()) == set(range(len(self.numParent)-1))
        # official names / common names
        tmp = collections.defaultdict(list)
        for (curr, off) in self.officialName.items():
            tmp[off].append(curr)
        self.commonNames = self.newCommonNamesMapperInstance()
        self.commonNames.update(tmp)
        self.dicGenes = {}
        self.dicGenomes = self.newCommonNamesMapperInstance()

        print(" OK", file=stream)

    def __init__(self, file=None, skipInit=False, stream=open(os.devnull, 'w')):
        if isinstance(file, tuple):
            print("Creation of the phylogenetic tree ...", end=' ', file=stream)
            (self.items, self.root, self.officialName) = file
        else:
            self.officialName = {}
            self.items = self.newCommonNamesMapperInstance()
            if file is None:
                return
            print("Loading phylogenetic tree %s ..." % file, end=' ', file=stream)
            # name and instance of file
            f = myFile.openFile(file, 'r')
            try:
                self.name = f.name
            except AttributeError:
                self.name = file
            f = myFile.firstLineBuffer(f)
            if (';' in f.firstLine) or ('(' in f.firstLine):
                self.__loadFromNewick__(' '.join(f).replace('\n', '') + " ;")
            else:
                self.__loadFromMyFormat__(f)
            f.close()
            if not skipInit:
                self.reinitTree()
            else:
                print("OK", file=stream)

    # load a phylogenetic tree into the phylTree format (with tabulations)
    def __loadFromMyFormat__(self, f):

        # loading process
        def loadFile():
            lines = []
            for line in f:
                # the final "\n" is removed and we cut owing to the "\t"
                l = line.split('\t')
                # indentation level
                indent = 0
                for x in l:
                    if x != '':
                        break
                    indent += 1
                # the triplet (indentation,names,age) is recorded
                assert l[indent] == x
                names = [xx.strip() for xx in x.split('|')]
                if len(l) > (indent+1):
                    age = int(l[indent+1])
                else:
                    age = 0
                lines.append((indent, names, age))
            lines.reverse()
            return lines

        # lines analysis process
        # lines = [..., (indent, names, age), ..., (0, namesOfRoot, age)]
        def recLoad(indent, lines):
            # we continue only if the next line is indented as expected
            if len(lines) == 0 or lines[-1][0] != indent:
                return None
            # the line is loaded
            currLine = lines.pop()
            # all the children sub-trees are loaded, as long as it is possible
            children = []
            while True:
                tmp = recLoad(indent + 1, lines)
                if tmp is None:
                    break
                # record (name_of_child, evolution_time)
                children.append((tmp, currLine[2] - self.ages.get(tmp)))
            # AGORA needs a minimised tree
            if len(children) == 1:
                return children[0][0]
            # current node
            (currIndent, currNames, currAge) = currLine
            name = currNames[0]
            # No children: the node is a leaf
            if len(children) == 0:
                if name[0] == SYMBOL6X:
                    currLine[1][0] = name = currNames[0][1:]
                    self.lstEsp6X.add(name)
                elif name[0] == SYMBOL2X:
                    currLine[1][0] = name = currNames[0][1:]
                    self.lstEsp2X.add(name)
                else:
                    name = currNames[0]
                    self.lstEspFull.add(name)
            # several children, they are recorded
            else:
                self.items.setdefault(name, children)
            # standard informations
            self.ages.setdefault(name, currLine[2])
            for s in currLine[1]:
                self.officialName[s] = name
            return name

        self.ages = self.newCommonNamesMapperInstance()
        self.lstEsp2X = set()
        self.lstEsp6X = set()
        self.lstEspFull = set()
        lines = loadFile()
        self.root = recLoad(0, lines)

    # tree in the Newick format (between parentheses)
    def __loadFromNewick__(self, s):

        #FIXME read the nb next characters of the tree
        def readStr(nb):
            ret = s[self.pos:self.pos+nb]
            self.pos += nb
            return ret

        # remove forbidden characters
        def keepWhile(car):
            x = self.pos
            while s[x] in car:
                x += 1
            return readStr(x-self.pos)

        # keep the authorised characters
        def keepUntil(car):
            x = self.pos
            while s[x] not in car:
                x += 1
            return readStr(x-self.pos)

        def readTree():
            def recParseTree(n):
                nodeName = n.name
                if nodeName:
                    if nodeName[0] == SYMBOL6X:
                        nodeName = nodeName[1:]
                        self.lstEsp6X.add(nodeName)
                    elif nodeName[0] == SYMBOL2X:
                        nodeName = nodeName[1:]
                        self.lstEsp2X.add(nodeName)
                    else:
                        self.lstEspFull.add(nodeName)

                children = [recParseTree(c) for c in n.descendants]
                elt = (children, n.name)
                return (elt, n.length, n.properties)

            trees = newick.loads(s)
            if len(trees) > 1:
                raise Exception(f"Expected 1 tree in file, found {len(trees)}")
            return recParseTree(trees[0])

        def calcAges(data):
            ((children, name), _, _) = data
            if len(children) == 0:
                res = 0
            else:
                ress = []
                for child in children:
                    (eltt, length, _) = child
                    ress.append(length + calcAges((eltt, _, _)))
                #FIXME, this is false for gene trees of ensembl
                #assert all([foo == ress[0] for foo in ress]) # all the children + branch lengths to the parent should have the same age
                res = max(ress)
            self.ages[name] = res
            return res

        # fill the variables of a GenericTree
        def storeTree(data):
            ((children, name), length, info) = data
            if (name == '') or (name in self.officialName):
                global nodeIndex
                name = "NAME_%d" % nodeIndex
                nodeIndex += 1
            self.officialName[name] = name
            self.info[name] = info

            if len(children) > 0:
                items = []
                for arbre in children:
                    items.append(storeTree(arbre))
                self.items.setdefault(name, items)

            self.root = name
            return (name, length)

        #FIXME add the age of ancestors for an easier conversion into myFormat
        self.pos = 0
        self.lstEsp2X = set()
        self.lstEsp6X = set()
        self.lstEspFull = set()
        data = readTree()
        self.info = {}
        storeTree(data)
        self.ages = self.newCommonNamesMapperInstance()
        calcAges(data)

    # print into phylTree format, with tabulations
    def printPhylTree(self, fh=sys.stdout):
        def recPrint(node, indent):
            names = myFile.myTSV.printLine([self.fileName[node]] + [x for x in self.commonNames.get(node, "") if isinstance(x, str) and (x != node)], delim="|")
            if node in self.listSpecies:
                print(("\t" * indent) + str(names), file=fh)
            else:
                print(("\t" * indent) + str(names) + "\t" + str(int(self.ages[node])), file=fh)
                for (f, _) in self.items[node]:
                    recPrint(f, indent+1)
        recPrint(self.root, 0)

    def printNewick(self, fh=sys.stdout):
        def recConvert(node):
            a = self.fileName[node]
            if node in self.listSpecies:
                return a
            else:
                return "(" + ",".join([recConvert(e) + ":" + str(l) for (e,l) in self.items[node]]) + ")%s" % a
        print(recConvert(self.root) + ";", file=fh)

    # return the name of the last common ancestor of several species
    def lastCommonAncestor(self, species):
        anc = species[0]
        for e in species[1:]:
            anc = self.dicParents[anc][e]
            if anc == self.root:
                return self.root
        return anc

    # assess if 'child' is really the child of 'parent'
    def isChildOf(self, child, parent):
        return self.dicParents[child][parent] == self.officialName[parent]

    #FIXME
    def compact(self, maxlength=1e-4):
        def do(node):
            if node not in self.items:
                return
            newitems = []
            for x in self.items[node]:
                do(x[0])
                if (x[1] >= maxlength) or (x[0] not in self.items):
                    newitems.append(x)
                else:
                    print("Removing node", x[0], file=sys.stderr)
                    newitems.extend(self.items.pop(x[0]))
            self.items[node] = newitems
        do(self.root)
        self.reinitTree()

    # return the list of species to be compared
    # ancestors must be prefixed with '.' to be taken into account as genomes
    #   otherwise, that is their descendant species that are used
    def getTargets(self, s):

        # are we exclusively on the common ancestors
        allIntermediates = (s[-1] == "+")
        if allIntermediates:
            s = s[:-1]

        # list of extant species
        listSpecies = set()
        listAncestors = set()
        for x in s.split(','):
            if x[0] != '.':
                listSpecies.update(self.species[x])
            else:
                listSpecies.add(x[1:])
                listAncestors.add(x[1:])

        # list of ancestors
        for (e1, e2) in itertools.combinations(listSpecies, 2):
            if allIntermediates:
                listAncestors.update(self.dicLinks[e1][e2][1:-1])
            else:
                listAncestors.add(self.dicParents[e1][e2])
            # if e1 or e2 is already an ancestor

        return (listSpecies, listAncestors)

    # return the list of species pointed out by target
    #   anc  [anc+children]
    #   =esp [esp]
    #   /anc [outgroups]
    #   _** in order to remove instead of adding
    def getTargetsSpec(self, target):
        lesp = set()
        for x in target.split(","):
            if x.startswith("_"):
                if x.startswith("_/"):
                    lesp.difference_update(self.outgroupSpecies[x[2:]])
                elif x.startswith("_="):
                    lesp.discard(x[2:])
                else:
                    lesp.difference_update(self.species[x[1:]])
            else:
                if x.startswith("/"):
                    lesp.update(self.outgroupSpecies[x[1:]])
                elif x.startswith("="):
                    lesp.add(x[1:])
                else:
                    lesp.update(self.species[x])
        return lesp

    # return the list of ancestors pointed out by target
    #   anc  [anc+children]
    #   =anc [anc]
    #   \anc [anc+parents]
    #   /anc [anc+parents+outgroups]
    #   _** in order to remove instead of adding
    def getTargetsAnc(self, target):
        lanc = set()
        for x in target.split(","):
            if x.startswith("_"):
                if x.startswith("_/"):
                    lanc.difference_update(e for e in self.listAncestr if not self.isChildOf(e, x[2:]))
                elif x.startswith("_="):
                    lanc.discard(x[2:])
                elif x.startswith("_\\"):
                    lanc.difference_update(e for e in self.listAncestr if self.isChildOf(x[2:], e))
                else:
                    lanc.difference_update(e for e in self.listAncestr if self.isChildOf(e, x[1:]))
            else:
                if x.startswith("/"):
                    lanc.update(e for e in self.listAncestr if not self.isChildOf(e, x[1:]))
                elif x.startswith("="):
                    lanc.add(x[1:])
                elif x.startswith("\\"):
                    lanc.update(e for e in self.listAncestr if self.isChildOf(x[1:], e))
                else:
                    lanc.update(e for e in self.listAncestr if self.isChildOf(e, x))
        return lanc

    # Wrapper around getTargetsAnc and getTargetsSpec to list the species and
    # ancestors to use in pairwise comparisons
    def getTargetsForPairwise(self, target, extantSpeciesFilter):

        # All the target ancestors
        targets = self.getTargetsAnc(target)
        # print "targets", targets

        # All the allowed extant species
        if extantSpeciesFilter:
            listSpecies = set()
            for anc in self.getTargetsSpec(extantSpeciesFilter):
                listSpecies.update(self.species[anc])
        else:
            listSpecies = set(self.listSpecies)
        # print "listSpecies", listSpecies

        # All the possible ancestors given the extant species
        # and the ones required to do comparisons
        listAncestors = set()
        requiredAncestors = set()
        for (e1,e2) in itertools.combinations(listSpecies, 2):
            inter = targets.intersection(self.dicLinks[e1][e2][1:-1])
            # If the pair of species intersects a target
            if inter:
                listAncestors.update(inter)
                requiredAncestors.add(self.dicParents[e1][e2])
        # print "listAncestors", listAncestors
        # print "requiredAncestors", requiredAncestors
        accessoryAncestors = requiredAncestors.difference(listAncestors)
        # print "accessoryAncestors", accessoryAncestors

        if not listAncestors:
            print("No target ancestors. Received arguments: target=%s extantSpeciesFilter=%s" % (target, extantSpeciesFilter), file=sys.stderr)
            sys.exit(1)

        return (listSpecies, listAncestors, accessoryAncestors)

    # return the structure of the sub-tree in which are exclusively the chosen species
    def getSubTree(self, goodSpecies, rootAnc=None):
        goodAnc = set(self.dicParents[e1][e2] for (e1, e2) in itertools.combinations(goodSpecies, 2))
        newtree = collections.defaultdict(list)
        from . import myMaths

        def do(node):
            if node in goodAnc:
                for (x, _) in self.items[node]:
                    newtree[node].extend(do(x))
                return [node]
            elif node in self.items:
                return myMaths.flatten([do(x) for (x, _) in self.items[node]])
            elif node in goodSpecies:
                return [node]
            else:
                return []
        root = do(self.root if rootAnc is None else rootAnc)
        assert len(root) == 1
        return (root[0], newtree)

    # return a dict that uses internally official names of taxons
    #  but that can be used with common names
    def newCommonNamesMapperInstance(self):
        dsi = dict.__setitem__
        dgi = dict.__getitem__
        off = self.officialName

        class commonNamesMapper(dict):

            def __getitem__(self, name):
                if name in off:
                    return dgi(self, off[name])
                else:
                    return dgi(self, name)

            def __setitem__(self, name, value):
                if name in off:
                    dsi(self, off[name], value)
                else:
                    dsi(self, name, value)

            # Because it's recursive, use a non-ambiguous name for this method,
            # since 'to_dict' for instance, is already a method of
            # pandas.DataFrames.
            def common_names_mapper_2_dict(self):
                """Recursively converts the commonNamesMapper instance to a
                dictionary. (to allow pickling for example)"""
                as_dict = {}
                for key, value in self.items():
                    if hasattr(value, 'common_names_mapper_2_dict'):
                        value = value.common_names_mapper_2_dict()
                    as_dict[key] = value
                return as_dict

        return commonNamesMapper()

    # load all the species that come from an ancestor
    def loadAllSpeciesSince(self, ancestr, template, **args):
        if ancestr in self.species:
            l = self.species[ancestr]
        else:
            l = self.listSpecies
        self.loadSpeciesFromList(l, template, **args)

    # load all the outgroup species of an ancestor
    def loadAllSpeciesBefore(self, ancestr, template, **args):
        if ancestr in self.outgroupSpecies:
            l = self.outgroupSpecies[ancestr]
        else:
            l = self.listSpecies
        self.loadSpeciesFromList(l, template, **args)

    # load all the species of a list
    def loadSpeciesFromList(self, lst, template, storeGenomes=True):

        from . import myGenomes

        for esp in lst:
            esp = self.officialName[esp]
            g = myGenomes.Genome(template % self.fileName[esp])
            if storeGenomes:
                self.dicGenomes[esp] = g
            for (x, (c, i)) in g.dicGenes.items():
                #self.dicGenes[x] = (esp, c, i)
                self.dicGenes[x] = GeneSpeciesPosition(esp, c, i)

    # return the number of genes in each species for a given family
    def findFamilyComposition(self, fam):
        score = self.newCommonNamesMapperInstance()
        score.update(dict.fromkeys(self.officialName, []))
        for g in fam:
            if g in self.dicGenes:
                (e, _, _) = self.dicGenes[g]
                score[e].append(g)
        return score

    # topple over ("basculer" in french) the tree around a node
    def reroot(self, node, useOutgroups, newname=0):
        self.tmpItems = self.items.copy()
        if useOutgroups:
            anc = node
            while anc in self.parent:
                (par, l) = self.parent[anc]
                self.tmpItems[anc] = self.tmpItems[anc] + [(par, l)]
                self.tmpItems[par] = [x for x in self.tmpItems[par] if x[0] != anc]
                anc = par
        self.tmpItems[newname] = self.tmpItems[node]

    def printBranchName(self, child, stream=sys.stderr):
        # To print branch name to sys.stderr
        # Please be aware that changing this breaks the analysis of the stderr by myScore.effectiveValuesFromSimulationLogErr()
        (parent, bLength) = self.parent[child]
        foo = "# Branch %s -> %s (%s My) #" % (parent, child, bLength)
        bar = "#" * len(foo)
        print(bar, file=stream)
        print(foo, file=stream)
        print(bar, file=stream)

    # FIXME, ensure that no two species has the same acronym!!
    def speciesAcronym(self, node):
        if node in self.listSpecies:
            # First char of each word with a upper char for the first
            speciesAcronym = [word for word in node.split(' ')]
            speciesAcronym = [speciesAcronym[0][0].upper()] + [char[0].lower() for char in speciesAcronym[1:]]
        elif node in self.listAncestr:
            # The two first chars
            speciesAcronym = node[:2]
        return ''.join(speciesAcronym)
