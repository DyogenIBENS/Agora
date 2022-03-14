#! /usr/bin/env python
# -*- coding: utf-8 -*-
# AGORA v3.1
# python 3.5
# Copyright © 2006-2022 IBENS/Dyogen, 2020-2021 EMBL-European Bioinformatics Institute, 2021-2022 Genome Research Ltd : Matthieu MUFFATO, Alexandra LOUIS, Thi Thuy Nga NGUYEN, Hugues ROEST CROLLIUS
# mail : agora@bio.ens.psl.eu
# This is free software; you may copy, modify and/or distribute this work under the terms of the GNU General Public License, version 3 or later and the CeCiLL v2 license in France

__doc__ = """
    Reformat the provided orthology groups (and genes) to suit AGORA. In particular, we want to
    ensure that all gene names are unique across the entire dataset, and we do so by prefixing
    them with the species name.

    Usage:
        src/ALL.reformatGeneFamilies.py example/data/Species.nwk example/data/orthologyGroups/orthologyGroups.%s.list.bz2 \
                -IN.genesFiles=example/data/genes/genes.%s.list.bz2 \
                -OUT.ancGenesFiles=example/results/ancGenes/all/ancGenes.%s.list.bz2 \
                -OUT.genesFiles=example/results/genes/genes.%s.list.bz2
"""

import collections
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

# Make sure the gene names are unique by adding the species name (both
# extant species and ancestors)
ancGeneNames = {}
nameMapping = collections.defaultdict(dict)
for species in sorted(phylTree.listSpecies) + sorted(phylTree.listAncestr):
    print("Renaming the genes of %s ..." % species, end=' ', file=sys.stderr)
    seen = set()
    if species in phylTree.listAncestr:
        ancGeneNames[species] = seen
    n = 0
    inputPath = arguments["IN.genesFiles"] % phylTree.fileName[species]
    if (not utils.myFile.hasAccess(inputPath)) and (species in phylTree.listAncestr):
        print("SKIPPING", file=sys.stderr)
        continue
    outputPath = arguments["OUT.genesFiles"] % phylTree.fileName[species]
    fi = utils.myFile.openFile(inputPath, "r")
    fo = utils.myFile.openFile(outputPath, "w")
    for line in fi:
        n += 1
        t = line[:-1].split("\t")
        assert len(t) == 5
        oldName = t[4]
        assert oldName not in seen
        seen.add(oldName)
        newName = phylTree.fileName[species] + "." + oldName
        nameMapping[oldName][species] = newName
        print(*t[:4], newName, sep="\t", file=fo)
    fi.close()
    fo.close()
    print(n, "OK", file=sys.stderr)

# Same for the ancGene: the ancGene's name itself, and its descendants
# Also restrict the list of descendants to the extant speciess
for anc in sorted(phylTree.listAncestr):
    print("Updating the ancestral families of %s ..." % anc, end=' ', file=sys.stderr)
    orthologyGroups = []
    inputPath = arguments["orthologyGroups"] % phylTree.fileName[anc]
    fi = utils.myFile.openFile(inputPath, "r")
    for l in fi:
        orthologyGroups.append(l.split())
    fi.close()
    hasAncGeneName = all(og[0] in ancGeneNames[anc] for og in orthologyGroups)
    print("with names" if hasAncGeneName else "adding names", "...", end=' ', file=sys.stderr)
    n = 0
    descendants = phylTree.species[anc]
    outputPath = arguments["OUT.ancGenesFiles"] % phylTree.fileName[anc]
    fo = utils.myFile.openFile(outputPath, "w")
    for og in orthologyGroups:
        n += 1
        ancGeneName = og.pop(0) if hasAncGeneName else str(n)
        ancGeneName = phylTree.fileName[anc] + "." + ancGeneName
        allNewNames = []
        assert len(og) == len(set(og))
        for g in og:
            newNames = [nameMapping[g][species] for species in descendants if species in nameMapping[g]]
            allNewNames.extend(newNames)
        print(ancGeneName, *allNewNames, file=fo)
    fo.close()
    print(n, "OK", file=sys.stderr)

