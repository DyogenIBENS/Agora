#! /usr/bin/env python
# -*- coding: utf-8 -*-
# AGORA v3.1
# python 3.5
# Copyright © 2006-2021 IBENS/Dyogen, 2020-2021 EMBL-European Bioinformatics Institute, 2021 Genome Research Ltd : Matthieu MUFFATO, Alexandra LOUIS, Thi Thuy Nga NGUYEN, Hugues ROEST CROLLIUS
# mail : agora@bio.ens.psl.eu
# This is free software; you may copy, modify and/or distribute this work under the terms of the GNU General Public License, version 3 or later and the CeCiLL v2 license in France

__doc__ = """
    Reformat the provided orthology groups (and genes) to suit AGORA. In particular, we want to
    ensure that all gene names are unique across the entire dataset, and we do so by prefixing
    them with the species name.

    Usage:
        src/ALL.reformatGeneFamilies.py example/data/Species.nwk TODO \
                -IN.genesFiles= \
                -OUT.ancGenesFiles=example/results/ancGenes/all/ancGenes.%s.list.bz2 \
                -OUT.genesFiles=example/results/genes/genes.%s.list.bz2
"""

import collections
import os
import sys

import utils.myFile
import utils.myGenomes
import utils.myPhylTree
import utils.myTools

# Arguments
###########

arguments = utils.myTools.checkArgs(
    [("speciesTree", utils.myTools.FileArgChecker), ("orthologyGroups", utils.myTools.PatternArgChecker)],
    [("IN.genesFiles", str, ""), ("OUT.ancGenesFiles", str, ""), ("OUT.genesFiles", str, "")],
    __doc__
)

phylTree = utils.myPhylTree.PhylogeneticTree(arguments["speciesTree"])

nameMapping = collections.defaultdict(dict)
for species in sorted(phylTree.listSpecies) + sorted(phylTree.listAncestr):
    inputPath = arguments["IN.genesFiles"] % phylTree.fileName[species]
    outputPath = arguments["OUT.genesFiles"] % phylTree.fileName[species]
    fi = utils.myFile.openFile(inputPath, "r")
    fo = utils.myFile.openFile(outputPath, "w")
    for line in fi:
        t = line[:-1].split("\t")
        oldName = t[4]
        assert (oldName not in nameMapping) or (species not in nameMapping[oldName])
        newName = phylTree.fileName[species] + "." + oldName
        nameMapping[oldName][species] = newName
        print(*t[:4], newName, sep="\t", file=fo)

# ancGenes are already present, just copy them over, but add a family name for unique reference
for anc in sorted(phylTree.listAncestr):
    inputPath = arguments["orthologyGroups"] % phylTree.fileName[anc]
    descendants = set(phylTree.species[anc])
    if os.path.exists(inputPath):
        print("Copying families of ancestral genome %s ..." % anc, end=' ', file=sys.stderr)
        code = 'FAM' + anc[:4].upper()
        outputPath = arguments["OUT.ancGenesFiles"] % phylTree.fileName[anc]
        fi = utils.myFile.openFile(inputPath, "r")
        fo = utils.myFile.openFile(outputPath, "w")
        n = 0
        for l in fi:
            t = l.split()
            allNewNames = []
            for g in set(t[1:]):
                newNames = [nameMapping[g][species] for species in descendants if species in nameMapping[g]]
                #assert newNames, (anc, descendants, t)
                #assert t.count(g) == len(newNames), (anc, t, g, [species for species in descendants if species in nameMapping[g]], newNames, t.count(g), len(newNames))
                allNewNames.extend(newNames)
            n += 1
            uniqueName = phylTree.fileName[anc] + "." + t[0]
            print(uniqueName, *allNewNames, file=fo)
        fo.close()
        fi.close()
        print(n, "OK", file=sys.stderr)
    else:
        print("No file for '%s' in '%s'" % (anc, arguments["orthologyGroups"]))

