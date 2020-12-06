#! /usr/bin/env python
# -*- coding: utf-8 -*-
# AGORA v2.0
# python 2.7
# Copyright © 2020 IBENS/Dyogen and EMBL-European Bioinformatics Institute : Matthieu MUFFATO, Alexandra LOUIS, Thi Thuy Nga NGUYEN, Hugues ROEST CROLLIUS
# mail : agora@bio.ens.psl.eu
# This is free software; you may copy, modify and/or distribute this work under the terms of the GNU General Public License, version 3 or later and the CeCiLL v2 license in France

__doc__ = """
	Tries to insert a maximum of singletons inside blocks:
		Function 1: successive addition of the edges starting from the best (limit: ~ 30 pairs)
		Function 2: -idem- (limit: ~ 40 pairs)
		Function * 3: recursive analysis of all paths between start and end, with successive selection of the maximum (limit: ~ 100 pairs)
				31: sum of the scores
				32: length
				33: -length
				34: scores sorted in descending order
				35: scores sorted in ascending order
				36: mean of scores
		Function * 4: selection of the pair of maximum weight ensures the existence of the path and divide the problem into two sub-problems (no size limit)
				40: size criterion for the choice of resolution functions of sub-problems
				4XX: XX use for sub-problems
	-func=0,f1|size2,f2|size3,f3t|sizel (t for using a thread+timeout, sizel for a maximum size of graph)

"""

import Queue
import collections
import itertools
import multiprocessing
import sys
import threading
import time

from joblib import Parallel, delayed

import utils.myFile
import utils.myGenomes
import utils.myGraph
import utils.myMaths
import utils.myPhylTree
import utils.myTools

# Arguments
arguments = utils.myTools.checkArgs(
    [("speciesTree", file), ("target", str), ("pairwise", str)],
    [("IN.ancBlocks", str, ""), ("OUT.ancBlocks", str, ""), ("LOG.ancGraph", str, "refine_log/%s.log.bz2"),
     ("nbThreads", int, 0),
     ("minimalWeight", int, 1), ("mustExtend", bool, False), ("loop", bool, False), ("timeout", int, 150),
     ("func", str, "0,32|100,40t|10000")],
    __doc__
)


# reverse gene
###############
def rev((g, s)):
    return (g, -s)


# Decorateur pour utilisation de "with"
#########################################
class setBacktracker(list):
    # Place des elements dans des ensembles
    def __enter__(self):
        for (s, x) in self:
            assert x not in s, (x, s)
            s.add(x)

    # Et les retire
    def __exit__(self, type, value, traceback):
        for (s, x) in self:
            s.remove(x)


# Met a jour la liste des successeurs/predecesseurs
####################################################
class GraphContainer:
    def __init__(self, lstIniPairwise=[]):

        if len(lstIniPairwise) > 0:
            print "graphcall", "begin", len(lstIniPairwise)

        als = self.allSucc = collections.defaultdict(set)
        alp = self.allPred = collections.defaultdict(set)

        for (_, g1, g2) in lstIniPairwise:

            if g2 in als[g1]:
                continue

            alpg1 = alp[g1]
            als[g1].add(g1)
            alpg1.add(g1)

            alsg2 = als[g2]
            alsg2.add(g2)
            alp[g2].add(g2)

            als[g1].update(alsg2)
            for t in alpg1:
                als[t].update(alsg2)

            alp[g2].update(alpg1)
            for t in alsg2:
                alp[t].update(alpg1)

        if len(lstIniPairwise) > 0:
            print "graphcall", "end", len(lstIniPairwise)

    def findPairwise(self, gfrom, gto, lstPairwise):

        als = self.allSucc
        alp = self.allPred
        alsfrom = als[gfrom]
        for (i, (s, g1, g2)) in enumerate(lstPairwise):

            if g2 not in als[g1]:
                alpg1 = alp[g1]
                als[g1].add(g1)
                alpg1.add(g1)

                alsg2 = als[g2]
                alsg2.add(g2)
                alp[g2].add(g2)

                als[g1].update(alsg2)
                for t in alpg1:
                    als[t].update(alsg2)

                alp[g2].update(alpg1)
                for t in alsg2:
                    alp[t].update(alpg1)

            if gto in alsfrom:
                yield (i, s, g1, g2)


def bestPath1(start, end, lstPairwise):
    allSucc = GraphContainer(lstPairwise).allSucc

    forbiddeng1 = set([rev(start), rev(end), end])
    forbiddeng2 = set([rev(start), rev(end), start])
    seen = set([start, end])

    def search(i0, parts, scores):
        for (i, (s, g1, g2)) in enumerate(itertools.islice(lstPairwise, i0, None)):
            if (g1 in forbiddeng1) or (g2 in forbiddeng2):
                continue

            for (la, lb) in utils.myTools.myIterator.slidingTuple(parts):
                ga = la[-1]
                gb = lb[0]
                if (g1 in allSucc[ga]) and (gb in allSucc[g2]):

                    newparts = parts[:]
                    newscores = scores[:]
                    j = parts.index(lb)
                    g1seen = (parts[j - 1][-1] == g1)
                    g2seen = (parts[j][0] == g2)
                    if g1seen:
                        newscores[j - 1] = newscores[j - 1] + [s]
                        if g2seen:
                            newparts[j - 1] = newparts[j - 1] + newparts[j]
                            del newparts[j]
                            newscores[j - 1] = newscores[j - 1] + newscores[j]
                            del newscores[j]
                        else:
                            newparts[j - 1] = newparts[j - 1] + [g2]
                    else:
                        if g2seen:
                            newparts[j] = [g1] + newparts[j]
                            newscores[j] = [s] + newscores[j]
                        else:
                            newparts.insert(j, [g1, g2])
                            newscores.insert(j, [s])

                    # Loop detection
                    if (g1seen != (g1 in seen)) or (g2seen != (g2 in seen)):
                        continue

                    if len(newparts) == 1:
                        return (newparts[0], newscores[0])

                    todo = setBacktracker()

                    todo.append((forbiddeng1, g1))
                    if not g1seen:
                        todo.append((seen, g1))
                        todo.extend([(forbiddeng1, rev(g1)), (forbiddeng2, rev(g1))])

                    todo.append((forbiddeng2, g2))
                    if not g2seen:
                        todo.append((seen, g2))
                        todo.extend([(forbiddeng1, rev(g2)), (forbiddeng2, rev(g2))])

                    with todo:
                        r = search(i0 + i + 1, newparts, newscores)
                        if r is not None:
                            return r

    return search(0, [[start], [end]], [[], []])


def bestPath2(start, end, lstPairwise):
    weight = {}
    for (s, g1, g2) in lstPairwise:
        weight[(g1, g2)] = s
    selectedF = set()
    selectedR = set()
    links = {}
    links[None] = start
    links[end] = None

    def search(i0):
        for (i, (s, g1, g2)) in enumerate(itertools.islice(lstPairwise, i0, None)):
            if (g1 in selectedF) or (g2 in selectedR):
                continue

            path = []
            curr = g2
            while curr in links:
                path.append(curr)
                curr = links[curr]

            if (len(path) == len(links)) and (curr == g1):
                try:
                    i = path.index(None)
                    # Solution
                    path = path[i + 1:] + [g1] + path[:i]
                    scores = [weight[x] for x in utils.myTools.myIterator.slidingTuple(path)]
                    return (path, scores)
                except ValueError:
                    # Loop
                    continue

            todo = setBacktracker()
            todo.append((selectedF, g1))
            todo.append((selectedR, rev(g1)))
            todo.append((selectedR, g2))
            todo.append((selectedF, rev(g2)))

            with todo:
                links[g1] = g2
                r = search(i0 + i + 1)
                if r is not None:
                    return r
                del links[g1]

    return search(0)


def _bestPath3(start, end, lstPairwise, key=None):
    allSucc = GraphContainer(lstPairwise).allSucc

    next = collections.defaultdict(list)
    for (s, g1, g2) in lstPairwise:
        next[g1].append((g2, s))
    seen = set()

    def search(node):
        if node not in next:
            return None
        poss = []
        for (g2, s) in next[node]:
            if g2[0] in seen:
                continue
            elif g2 == end:
                poss.append(([node, end], [s]))
            elif end in allSucc[g2]:
                seen.add(g2[0])
                path = search(g2)
                if path is not None:
                    poss.append(([node] + path[0], [s] + path[1]))
                seen.remove(g2[0])
        if len(poss) > 0:
            return max(poss, key=key)

    return search(start)


def bestPath31(*args):
    return _bestPath3(*args, key=lambda x: sum(x[1]))


def bestPath32(*args):
    return _bestPath3(*args, key=lambda x: len(x[0]))


def bestPath33(*args):
    return _bestPath3(*args, key=lambda x: -len(x[0]))


def bestPath34(*args):
    return _bestPath3(*args, key=lambda x: sorted(x[1], reverse=True))


def bestPath35(*args):
    return _bestPath3(*args, key=lambda x: sorted(x[1]))


def bestPath36(*args):
    return _bestPath3(*args, key=lambda x: sum(x[1]) / float(len(x[1])))


def _bestPath4(start, end, lstPairwise, subfunc=None):
    def get(gfrom, gto, lstPairwise):

        sp = GraphContainer()
        spals = sp.allSucc
        spalsfrom = spals[gfrom]
        for (i, s, g1, g2) in sp.findPairwise(gfrom, gto, lstPairwise):

            avoid = set([g1, g2, rev(g1), rev(g2)])
            if gfrom == g1:
                r1 = ([g1], [])
            else:
                diag1 = [x for x in itertools.islice(lstPairwise, i) if
                         (x[1] in spalsfrom) and (g1 in spals[x[2]]) and (avoid.isdisjoint(x) or (x[2] == g1))]
                r1 = subfunc(gfrom, g1, diag1)
                if r1 is None:
                    continue

            if gto == g2:
                r2 = ([g2], [])
            else:
                avoid.update(r1[0])
                avoid.update(rev(x) for x in r1[0])
                diag2 = [x for x in itertools.islice(lstPairwise, i) if
                         (x[1] in spals[g2]) and (gto in spals[x[2]]) and (avoid.isdisjoint(x) or (x[1] == g2))]
                r2 = subfunc(g2, gto, diag2)
                if r2 is None:
                    continue

            return (r1[0] + r2[0], r1[1] + [s] + r2[1])

    if subfunc is None:
        subfunc = get

    return get(start, end, lstPairwise)


def bestPath40(*args):
    def reselectSize(*args):
        for (size, f, thread) in func:
            if len(args[-1]) >= size:
                break
        # print f, args[0], args[1], len(args[-1]), "subpairs", args[-1]
        print f, args[0], args[1], len(args[-1])
        return f(*args)

    return _bestPath4(*args, subfunc=reselectSize)


def bestPath41(*args):
    return _bestPath4(*args, subfunc=bestPath1)


def bestPath42(*args):
    return _bestPath4(*args, subfunc=bestPath2)


def bestPath431(*args):
    return _bestPath4(*args, subfunc=bestPath31)


def bestPath432(*args):
    return _bestPath4(*args, subfunc=bestPath32)


def bestPath433(*args):
    return _bestPath4(*args, subfunc=bestPath33)


def bestPath434(*args):
    return _bestPath4(*args, subfunc=bestPath34)


def bestPath435(*args):
    return _bestPath4(*args, subfunc=bestPath35)


def bestPath436(*args):
    return _bestPath4(*args, subfunc=bestPath36)


def bestPath44(*args):
    return _bestPath4(*args, subfunc=None)


# Rajoute les pairwise qui commencent en "start" ou finissent en "end"
########################################################################
def prepareGraph(pairwiseDiags, singletons, graphData, start, end):
    print "preparecall", start, end
    newPairwise = []
    startlink = set([start])
    endlink = set([end])

    for (g2, s) in pairwiseDiags[start].iteritems():
        if (g2[0] in singletons) or ((g2 == end) and not arguments["mustExtend"]):
            startlink.add(g2)
            startlink.update(graphData.allSucc[g2])
            newPairwise.append((s, start, g2))
            newPairwise.append((s, rev(g2), rev(start)))

    for (g1, s) in pairwiseDiags[rev(end)].iteritems():
        if g1[0] in singletons:
            endlink.add(rev(g1))
            endlink.update(graphData.allPred[rev(g1)])
            newPairwise.append((s, rev(g1), end))
            newPairwise.append((s, rev(end), g1))

    print "addpairwise", len(newPairwise)
    print "startlink", len(startlink)
    print "endlink", len(endlink)

    return (newPairwise, startlink, endlink)


def do(anc):
    # Redirect the standard output to a file
    ini_stdout = sys.stdout
    sys.stdout = utils.myFile.openFile(arguments["LOG.ancGraph"] % phylTree.fileName[anc], "w")

    pairwiseDiags = loadPairwise(arguments["pairwise"] % phylTree.fileName[anc])

    (integr, singletons) = utils.myGraph.loadIntegr(arguments["IN.ancBlocks"] % phylTree.fileName[anc])
    newintegr = integr

    while len(integr) > 0:

        print "loop"
        nsing = len(singletons)

        print >> sys.stderr, "Building graphs ...",

        # Toutes les paires entre singletons
        commonPairwise = []
        for (g1, l) in pairwiseDiags.iteritems():
            if g1[0] in singletons:
                for (g2, s) in l.iteritems():
                    if g2[0] in singletons:
                        commonPairwise.append((s, g1, g2))
        print "commonpairwise", len(commonPairwise)

        # La structure de graphe de base
        sp = GraphContainer(commonPairwise)
        # print "inisucc", len(sp.allSucc), sum(len(x) for x in sp.allSucc.itervalues())

        # Dictionnaire gene -> intervalles
        dicGenes = collections.defaultdict(set)
        # Dictionnaire intervalle -> gene
        dicInterv = collections.defaultdict(set)
        # La liste des paires utilisables pour chaque intervalle
        goodPairwise = {}
        for (i, (b, _)) in enumerate(integr):
            for (j, p) in enumerate(utils.myTools.myIterator.slidingTuple(b)):
                (newPairwise, startlink, endlink) = prepareGraph(pairwiseDiags, singletons, sp, p[0], p[1])
                if startlink.isdisjoint(endlink):
                    print "nogoodpairwise"
                else:
                    interv = (i, j)
                    goodPairwise[interv] = [x for x in commonPairwise if (x[1] in startlink) and (x[2] in endlink)]
                    goodPairwise[interv].extend(x for x in newPairwise if (x[1] in startlink) and (x[2] in endlink))
                    goodPairwise[interv].sort(reverse=True)
                    print "goodpairwise", len(goodPairwise[interv])
                    for (_, g1, g2) in goodPairwise[interv]:
                        if g1 != p[0]:
                            dicGenes[g1[0]].add(interv)
                            dicInterv[interv].add(g1[0])
                        if g2 != p[1]:
                            dicGenes[g2[0]].add(interv)
                            dicInterv[interv].add(g2[0])
        print >> sys.stderr, "OK"

        # Supprime l'intervalle de la liste de ceux a parcourir et met a jour la liste des singletons
        def applyResult(interv, filtered, todelete):
            # print "apply %d/%d" % interv
            if interv in goodPairwise:
                del goodPairwise[interv]
            if interv not in res:
                return
            t = [x[0] for x in res[interv][0][1:-1]]
            # print "included genes", t
            for x in t:
                singletons.remove(x)
                filtered.update(dicGenes[x])
            todelete.update(t)
            print len(singletons), "singletons remaining"

        def bestPathWrapper(interv):

            start = integr[interv[0]][0][interv[1]]
            end = integr[interv[0]][0][interv[1] + 1]
            lstPairwise = goodPairwise[interv]

            print "searchcall", "%d/%d" % interv, start, end, "pairs", len(lstPairwise), lstPairwise

            if len(lstPairwise) >= maxsize:
                print "toobig"

            else:
                try:
                    # Recherche du pattern de choix de fonction
                    for (size, f, thread) in func:
                        if len(lstPairwise) >= size:
                            break
                    if thread:

                        queue = Queue.Queue()
                        def bestPathInQueue(f, args):
                            queue.put(f(*args))

                        print "with thread"
                        p = threading.Thread(target=bestPathInQueue, args=(f, (start, end, lstPairwise)))
                        p.start()
                        # st = time.time()
                        try:
                            r = queue.get(True, arguments["timeout"])
                        except Queue.Empty:
                            p._Thread__stop()
                            # Au cas ou le resultat serait arrive entre temps
                            r = queue.get_nowait()
                        # et = time.time()
                        p.join()
                    else:
                        # st = time.clock()
                        r = f(start, end, lstPairwise)
                        # et = time.clock()
                    # print "duration", et - st
                    # global alltime
                    # alltime += (et - st)

                    if r is None:
                        print "nopath"
                    else:
                        print "solution", len(r[0]) - 2, sum(r[1]), r[0], r[1]
                        assert r[0][0] == start
                        assert r[0][-1] == end
                        assert len(r[0]) == len(r[1]) + 1
                        res[interv] = r

                except Queue.Empty:
                    # Graphe trop grand a examiner, il faudra revenir
                    print "queuetimeout"
                    p.join()
                    return False
            return True

        # Contient les resultats
        res = {}
        timeouts = set()

        print >> sys.stderr, "adding singletons ...",
        while len(goodPairwise) > 0:

            print "nb intervals todo", len(goodPairwise)

            # Calcul des meilleurs chemins
            paths = []
            filtered = set()
            todelete = set()
            for interv in goodPairwise.keys():
                res.pop(interv, None)
                if bestPathWrapper(interv):
                    if interv in res:
                        if len(res[interv][0]) > 2:
                            # Resultat: chemin
                            paths.append(interv)
                        else:
                            # Resultat: lien direct
                            applyResult(interv, filtered, todelete)
                    else:
                        # Resultat: pas de chemin
                        applyResult(interv, filtered, todelete)
                else:
                    # Pas de resultat
                    timeouts.add(interv)
            print len(timeouts), "timeouts"

            # Association gene -> intervalles qui l'utilisent
            counts = collections.defaultdict(list)
            for interv in paths:
                for g in res[interv][0][1:-1]:
                    counts[g[0]].append(interv)
            print len(counts), "referenced genes"

            # Regroupement des intervalles avec des resultats qui s'excluent
            comb = utils.myTools.myCombinator()
            for l in counts.itervalues():
                comb.addLink(l)

            for g in comb:
                if len(g) == 1:
                    print "autonomous interv %d/%d" % g[0]
                    applyResult(g[0], filtered, todelete)
                else:
                    print "mutually exluded interv", len(g), ["%d/%d" % x for x in g]
                    # On selectionne le resultat le mieux soutenu
                    scores = []
                    for interv in g:
                        s = 0
                        for ((g1, g2), w) in itertools.izip(utils.myTools.myIterator.slidingTuple(res[interv][0]),
                                                            res[interv][1]):
                            if (len(counts[g1[0]]) > 1) or (len(counts[g2[0]]) > 1):
                                s += w
                        scores.append((s, interv))
                    scores.sort(reverse=True)

                    # On applique son resultat
                    print "best interv", [(x[0], "%d/%d" % x[1]) for x in scores]
                    for (_, interv) in scores:
                        if todelete.isdisjoint(x[0] for x in res[interv][0][1:-1]):
                            applyResult(interv, filtered, todelete)

            # Mise a jour des intervalles en enlevant les liens pairwise qui ne sont plus disponibles
            filtered.intersection_update(goodPairwise)
            print "filtered intervals", ["%d/%d" % x for x in filtered]
            for x in filtered:
                print "filtercall", "%d/%d" % x
                goodPairwise[x] = [l for l in goodPairwise[x] if
                                   (l[1][0] not in todelete) and (l[2][0] not in todelete)]
                allSucc = GraphContainer(goodPairwise[x]).allSucc
                timeouts.discard(x)
                newstart = integr[x[0]][0][x[1]]
                newend = integr[x[0]][0][x[1] + 1]
                if (newstart not in allSucc) or (newend not in allSucc[newstart]):
                    print "nolink", "%d/%d" % x
                    del goodPairwise[x]
                    res.pop(x, None)

            for interv in timeouts:
                applyResult(interv, filtered, todelete)

        print >> sys.stderr, "OK"

        # Rassemblement des nouveaux blocs ancestraux
        print >> sys.stderr, "Final blocks of", anc,
        newintegr = []
        for (i, (b, s)) in enumerate(integr):
            newb = []
            news = []
            newintegr.append((newb, news))
            for (j, (p, w)) in enumerate(itertools.izip(utils.myTools.myIterator.slidingTuple(b), s)):
                r = res.get((i, j))
                if r is None:
                    newb.append(p[0])
                    news.append(pairwiseDiags[p[0]][p[1]] if p[1] in pairwiseDiags[p[0]] else -w)
                else:
                    newb.extend(r[0][:-1])
                    news.extend(r[1])
            newb.append(p[1])
        print >> sys.stderr, utils.myMaths.myStats.txtSummary([len(x[0]) for x in newintegr]), "+", len(
            singletons), "singletons"

        if (not arguments["loop"]) or (len(singletons) == nsing):
            break
        integr = newintegr

    # Impression des resultats finaux
    f = utils.myFile.openFile(arguments["OUT.ancBlocks"] % phylTree.fileName[anc], "w")
    for (i, (newb, news)) in enumerate(newintegr):
        print >> f, utils.myFile.myTSV.printLine(
            [anc, len(newb), utils.myFile.myTSV.printLine([x[0] for x in newb], delim=" "),
             utils.myFile.myTSV.printLine([x[1] for x in newb], delim=" "),
             utils.myFile.myTSV.printLine(news, delim=" ")])
    for x in singletons:
        print >> f, utils.myFile.myTSV.printLine([anc, 1, x, 1, ""])
    f.close()

    # Revert to the true standard output
    sys.stdout.close()
    sys.stdout = ini_stdout


start = time.time()
phylTree = utils.myPhylTree.PhylogeneticTree(arguments["speciesTree"])
targets = phylTree.getTargetsAnc(arguments["target"])

# Initialisation des fonctions de recherche
func = [x.split(",") for x in arguments["func"].split("|")]
maxsize = int(func[-1][0])
func = [(int(x[0]), eval("bestPath" + x[1].replace("t", "")), "t" in x[1]) for x in func[:-1]]
func.reverse()


def loadPairwise(file):
    pairwiseDiags = collections.defaultdict(dict)
    lstPairwise = utils.myGraph.loadConservedPairsAnc(file)
    for d in lstPairwise:
        if d[2] >= arguments["minimalWeight"]:
            pairwiseDiags[d[0]][d[1]] = d[2]
            pairwiseDiags[rev(d[1])][rev(d[0])] = d[2]
    return pairwiseDiags


print >> sys.stderr, targets

n_cpu = arguments["nbThreads"] or multiprocessing.cpu_count()
Parallel(n_jobs=n_cpu)(delayed(do)(anc) for anc in targets)

print >> sys.stderr, "total computation time", (time.time() - start)
