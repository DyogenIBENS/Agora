#! /usr/bin/env python
# -*- coding: utf-8 -*-
# Agora v1.00
# python 2.7
# Copyright © 2016 IBENS/Dyogen : Matthieu MUFFATO, Alexandra Louis, Nga Thi tuy Nguyen, Joseph Lucas, Hugues ROEST CROLLIUS
# mail : hrc@ens.fr or alouis@biologie.ens.fr
# This is free software; you may copy, modify and/or distribute this work under the terms of the GNU General Public License, version 3 or later and the CeCiLL v2 license in France

import os
import subprocess
import sys
import time

import utils.myFile
import utils.myPhylTree
import utils.myTools

__doc__ = """
    Run all the AGORA programs thanks to configuration file

    Usage:
          ../src/agora.py ../conf/agora.ini
          ../src/agora.py ../conf/agora.ini -workingDir=. -nbThreads=4
"""

arguments = utils.myTools.checkArgs(
    [("agora.conf", file)],
    [("workingDir", file, "."), ("nbThreads", int, 1),
     ("prog:synteny", str, os.path.dirname(os.path.abspath(__file__)) + "/buildSynteny.%s-%s.py"),
     ("prog:ancGenes", str, os.path.dirname(os.path.abspath(__file__)) + "/ALL.filterGeneFamilies-%s.py")
     ],
    __doc__)

# loading configuration file
################################
bysections = {}
f = utils.myFile.openFile(arguments["agora.conf"], "r")
for l in f:

    l = l.partition("#")[0].strip()
    if l.startswith(">"):
        curr = []
        bysections[l[1:].strip().lower()] = curr
    elif len(l) > 1:
        curr.append(l)
f.close()


# Cut a string for delim
##########################
def partition(s, delim):
    x = s.partition(delim)
    return (x[0].strip(), x[2].strip())


# Files and directory section
##############################
files = {}
for x in bysections["files"]:
    x = partition(x, "=")
    files[x[0].lower()] = x[1]
phylTree = utils.myPhylTree.PhylogeneticTree(files["speciestree"])


# Managing the list of programs to launch and their dependencies
#################################################################
class TaskList():
    def __init__(self):
        self.list = []
        self.dic = {}

    def __len__(self):
        return len(self.list)

    def addTask(self, name, dep, data):
        print "New task", len(self.list), name
        print dep
        print data
        print
        self.dic[name] = len(self.list)
        self.list.append((set(self.dic[x] for x in dep), data))

    def removeDep(self, i):
        for (dep, _) in self.list:
            dep.discard(i)

    def getAvailable(self):
        tmp = [i for (i, task) in enumerate(self.list) if len(task[0]) == 0]
        print "Available tasks:", tmp
        if len(tmp) > 0:
            next = tmp[0]
            self.list[next][0].add(None)
            return (next,) + self.list[next]
        else:
            return None


tasklist = TaskList()
tasklist.addTask(("ancgenes", "all"), [], (None, None, None, False))

# Ancestral genes lists Section
################################


ancGenes = {"all": "all", "0": "all"}
for x in bysections["ancgenes"]:
    x = partition(x, "=")
    (params, root) = partition(x[1], "!")
    if root == "":
        root = phylTree.root
    # index
    t = x[0].split(",")
    sizes = params.split(",")
    minsize = ""
    maxsize = ""
    for i in range(len(sizes)):
        size = sizes[i].split()
        minsize += str(size[0]) + ","
        maxsize += str(size[1]) + ","
        dirname = "size-" + str(size[0]) + "-" + str(size[1])
        ancGenes[t[i].replace("*", "")] = dirname
    minsize = minsize[:-1]
    maxsize = maxsize[:-1]

tasklist.addTask(
    ("ancgenes", "size"),
    [],
    (
        [arguments["prog:ancGenes"] % "size", files["speciestree"], root,
         files["ancgenesdata"] % {"filt": "all", "name": "%s"},
         files["ancgenesdata"] % {"filt": "size-%s-%s", "name": "%s"}] + [minsize, maxsize],
        files["ancgenesoutput"],
        files["ancgeneslog"],
        "*" not in x[0]
    )
)

# Pairwise comparison section
#############################


for x in bysections.get("pairwise", []):
    x = partition(x, "=")
    dirname = ancGenes[x[0].replace("*", "").strip()]
    (params, root) = partition(x[1], "!")
    params = params.split()
    if root == "":
        root = phylTree.root

    # Pairwise comparison tasks
    tasklist.addTask(
        ("pairwise", dirname),
        [],
        (
            [arguments["prog:synteny"] % ("pairwise", params[0]), files["speciestree"], root,
             "-ancGenesFiles=" + files["ancgenesdata"] % {"filt": dirname, "name": "%s"},
             "-genesFiles=" + files["genes"] % {"name": "%s"},
             "-OUT.pairwise=" + files["pairwiseoutput"] % {"filt": dirname, "name": "%s"}] + params[1:],
            # files["pairwiseoutput"] % {"filt": dirname},
            "abc",
            files["pairwiselog"] % {"filt": dirname},
            "*" not in x[0]
        )
    )

# Integration section
#####################
interm = {}
refMethod = {}
for x in bysections.get("integration", []):
    tolaunch = ("*" not in x)

    (params, root) = partition(x, "!")
    (params, input) = partition(params, "<")
    (params, output) = partition(params, ">")
    params = params.replace("*", "").split()

    # method's name
    currMethod = params.pop(0)[1:-1] if params[0].startswith("[") else params[0]

    # used compared pairs
    dirname = None
    if params[-1].startswith("("):
        dirname = ancGenes[params.pop()[1:-1]]
        currMethod = currMethod[:-1] if currMethod.endswith("/") else (currMethod + "-" + dirname)
    elif currMethod.endswith("/"):
        currMethod = currMethod[:-1]

    if params[0] == "denovo":
        newMethod = currMethod
        if root == "":
            root = phylTree.root
        refMethod[newMethod] = (newMethod, root)
    else:
        # input data
        prevMethod = interm[input] if input != "" else newMethod
        newMethod = currMethod[1:] if currMethod.startswith("/") else (prevMethod + "." + currMethod)
        if root == "":
            root = refMethod[prevMethod][1]
        refMethod[newMethod] = (newMethod, root) if params[0] == "refine" else refMethod[prevMethod]

    if output != "":
        interm[output] = newMethod

    # task parameters
    args = [arguments["prog:synteny"] % ("integr", params[0]), files["speciestree"], root] + params[1:] + [
        "-OUT.ancDiags=" + files["integrblocks"] % {"method": newMethod, "name": "%s"}]
    dep = []
    if dirname is not None:
        dep.append(("pairwise", dirname))
        args.append(files["pairwiseoutput"] % {"filt": dirname, "name": "%s"})


    if params[0] == "denovo":
        args.append("-ancGenesFiles=" + files["ancgenesdata"] % {"filt": "all", "name": "%s"})
    else:
        # there're input blocks only for not denovo integration ( != "denovo" )
        dep.append(("integr", prevMethod))
        args.append("-IN.ancDiags=" + files["integrblocks"] % {"method": prevMethod, "name": "%s"})

    if params[0] == "halfinsert":
        # The script needs singleton reference for "halfinsert"
        dep.append(("integr", refMethod[newMethod][0]))
        args.append("-REF.ancDiags=" + files["integrblocks"] % {"method": refMethod[newMethod][0], "name": "%s"})

    if params[0] == "groups":
        args.append("-ancGenesFiles=" + files["ancgenesdata"] % {"filt": "all", "name": "%s"})
        args.append("-genesFiles=" + files["genes"] % {"name": "%s"})

    # Integration task
    tasklist.addTask(
        ("integr", newMethod),
        dep,
        (
            args,
            files["integroutput"] % {"method": newMethod},
            files["integrlog"] % {"method": newMethod},
            tolaunch
        )
    )

# Launching tasks in multiple threads
#####################################

os.chdir(arguments["workingDir"])

import multiprocessing

manager = multiprocessing.Manager()
queue = manager.Queue()

nrun = 0
proc = {}
completed = 0


# Waiting for a task to finish
def joinnext():
    print "Waiting ...",
    sys.stdout.flush()
    i = queue.get()
    print i, "is now finished"
    tasklist.removeDep(i)
    proc.pop(i).join()
    global completed
    completed += 1
    global nrun
    nrun -= 1


# Launch  program function
def golaunch(i, args, out, log):
    stdout = utils.myFile.openFile(out, "w")
    stderr = utils.myFile.openFile(log, "w")
    p = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=stderr)
    for l in p.stdout:
        print >> stdout, l,
    p.wait()
    for l in p.stdout:
        print >> stdout, l,
    stdout.close()
    stderr.close()
    time.sleep(5)
    queue.put(i)


# Queue
while completed < len(tasklist):

    print "todo:", len(tasklist)
    print "done:", completed

    if nrun == arguments["nbThreads"]:
        joinnext()
    else:
        todo = tasklist.getAvailable()
        if todo is None:
            joinnext()
        else:
            (next, dep, (args, out, log, launch)) = todo
            print ("Launching" if launch else "Skipping"), next, args, ">", out, "2>", log
            if launch:
                proc[next] = multiprocessing.Process(target=golaunch, args=(next, args, out, log))
                proc[next].start()
                nrun += 1
            else:
                tasklist.removeDep(next)
                completed += 1

assert nrun == 0
