#! /usr/bin/env python
# -*- coding: utf-8 -*-
# AGORA v2.1
# python 2.7
# Copyright © 2021 IBENS/Dyogen and EMBL-European Bioinformatics Institute : Matthieu MUFFATO, Alexandra LOUIS, Thi Thuy Nga NGUYEN, Hugues ROEST CROLLIUS
# mail : agora@bio.ens.psl.eu
# This is free software; you may copy, modify and/or distribute this work under the terms of the GNU General Public License, version 3 or later and the CeCiLL v2 license in France

import itertools
import multiprocessing
import os
import resource
import subprocess
import sys
import time
import threading

import psutil

from . import myFile

# Managing the list of programs to launch and their dependencies
#################################################################
class TaskList():

    rusage_unit = 1 if sys.platform == "darwin" else 1024

    def __init__(self):
        self.list = []
        self.dic = {}
        self.nrun = 0
        self.proc = {}
        self.nthreads = {}
        self.completed = 0
        self.failed = 0
        manager = multiprocessing.Manager()
        self.self_pids = [os.getpid(), manager._process.pid]
        self.queue = manager.Queue()
        self.memusage = manager.dict()
        self.memlock = manager.Lock()

    def printGraphviz(self, fh):
        print >> fh, "digraph", "{"
        for (name, taskId) in self.dic.iteritems():
            print >> fh, '%d [label="%s"]' % (taskId, "/".join(name))
            for dep in self.list[taskId][0]:
                print >> fh, '%d -> %d' % (dep, taskId)
        print >> fh, "}"

    def addTask(self, name, dep, data, multithreaded=False):
        taskId = len(self.list)
        print "New task", taskId, name
        print dep
        print data
        self.list.append((set(self.dic[x] for x in dep), data, multithreaded))
        if name in self.dic:
            if name + ("1",) in self.dic:
                self.list[self.dic[name]][0].add(taskId)
                for i in itertools.count(2):
                    newName = name + (str(i),)
                    if newName not in self.dic:
                        print "! Name clash ! Renamed to", newName, "under a collector"
                        self.dic[newName] = taskId
                        break
            else:
                print "! Name clash ! Introducing a collector task"
                collectorId = self.dic[name]
                self.list.append(self.list[collectorId])
                self.list[collectorId] = (set([taskId, taskId + 1]), (None, None, None, False), False)
                self.dic[name + ("1",)] = taskId + 1
                self.dic[name + ("2",)] = taskId
        else:
            self.dic[name] = taskId
        print
        return taskId

    def removeDep(self, i):
        for (dep, _, _) in self.list:
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

    # Waiting for a task to finish
    def joinNext(self):
        print "Waiting ..."
        sys.stdout.flush()
        (i, r) = self.queue.get()
        print "task", i, "is now finished (status %d)" % r
        self.proc.pop(i).join()
        if r == 0:
            self.removeDep(i)
            self.completed += 1
        else:
            self.failed += 1
            print >> sys.stderr, ">", "Inspect", self.list[i][1][2], "for more information"
        self.nrun -= self.nthreads.pop(i)

    def getProcMemoryUsage(self, proc):
        try:
            # Slower and not available on macOS, but more accurate
            return proc.memory_full_info().uss
        except psutil.AccessDenied:
            return proc.memory_info().rss

    def getMemoryUsage(self, pid):
        try:
            proc = psutil.Process(pid)
            return self.getProcMemoryUsage(proc)
        except (psutil.NoSuchProcess, IOError):
            return None

    def getRecursiveMemoryUsage(self, pid):
        try:
            proc = psutil.Process(pid)
            total_mem = self.getProcMemoryUsage(proc)
            children = proc.children(recursive=True)
        except (psutil.NoSuchProcess, IOError):
            return None
        for subproc in children:
            try:
                mem = self.getProcMemoryUsage(subproc)
                total_mem += mem
            except (psutil.NoSuchProcess, IOError):
                pass
        return total_mem

    def updateMemoryUsage(self, pid, mem):
        # Can't do this atomically, so need to wrap it with a lock
        self.memlock.acquire()
        if (pid in self.memusage) and (mem > self.memusage[pid]):
            self.memusage[pid] = mem
        self.memlock.release()

    def memoryMonitor(self):
        while self.memusage:
            total_mem = 0
            for pid in self.memusage.keys():
                mem = self.getRecursiveMemoryUsage(pid)
                if mem:
                    total_mem += mem
                    self.updateMemoryUsage(pid, mem)
            for pid in self.self_pids:
                mem = self.getMemoryUsage(pid)
                if mem:
                    total_mem += mem
            self.updateMemoryUsage(self.self_pids[0], mem)
            time.sleep(5)

    def printCPUUsageStats(self, intro, start):
        ru = resource.getrusage(resource.RUSAGE_CHILDREN)
        elapsed = time.time() - start
        # Use the lock so that we don't remove the key in the middle of memoryMonitor using it
        self.memlock.acquire()
        mem = max(ru.ru_maxrss * self.rusage_unit, self.memusage.pop(os.getpid()))
        self.memlock.release()
        print intro, "%g sec CPU time / %g sec elapsed = %g%% CPU usage, %g MB RAM" % (ru.ru_utime + ru.ru_stime, elapsed, 100. * (ru.ru_utime + ru.ru_stime) / elapsed, mem / 1024. / 1024.)

    # Launch program function
    def goLaunch(self, i, args, out, log):
        start = time.time()
        stdout = myFile.openFile(out or os.devnull, "w")
        stderr = myFile.openFile(log or os.devnull, "w")
        # stderr must have a fileno, so must be a regular file (not a .bz2 etc)
        # stdout can be anything, incl. a .bz2
        try:
            p = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=stderr)
        except Exception as e:
            stdout.close()
            stderr.close()
            print "task %d could not start:" % i, e
            time.sleep(5)
            self.queue.put((i, -1))
            # FIXME: then it hangs in multiprocessing/managers (checking self.memusage)
            return

        # This is where the .bz2 compression would happen
        for l in p.stdout:
            print >> stdout, l,
        r = p.wait()
        for l in p.stdout:
            print >> stdout, l,
        stdout.close()
        stderr.close()
        self.printCPUUsageStats("task %d report:" % i, start)
        time.sleep(5)
        self.queue.put((i, r))

    # Launching tasks in multiple threads
    def runAll(self, nbThreads, sequential):
        start = time.time()

        self.memusage[os.getpid()] = 0
        monitorThread = threading.Thread(target=self.memoryMonitor)
        monitorThread.start()

        # Queue
        while (self.completed + self.failed) < len(self.list):

            print "Status: %d to do, %d running, %d done, %d failed -- %d total" % \
                    (len(self.list)-len(self.proc)-self.completed-self.failed, len(self.proc), self.completed, self.failed, len(self.list))

            if (self.nrun == nbThreads) or (sequential and self.nrun):
                self.joinNext()
            else:
                todo = self.getAvailable()
                if todo is None:
                    if self.nrun == 0:
                        print "Workflow stopped because of failures"
                        break
                    self.joinNext()
                else:
                    (next, dep, (args, out, log, launch), multithreaded) = todo
                    if launch:
                        print "Launching task", next, args, ">", out, "2>", log
                        if multithreaded:
                            self.nthreads[next] = nbThreads - self.nrun
                            args.append("-nbThreads=%d" % self.nthreads[next])
                            print "Using", self.nthreads[next], "threads"
                        else:
                            self.nthreads[next] = 1
                        self.proc[next] = multiprocessing.Process(target=self.goLaunch, args=(next, args, out, log))
                        self.proc[next].start()
                        self.memusage[self.proc[next].pid] = 0
                        self.nrun += self.nthreads[next]
                    else:
                        print "Skipping task", next, args
                        self.removeDep(next)
                        self.completed += 1

        assert self.nrun == 0
        if not self.failed:
            print "Workflow complete"
        self.printCPUUsageStats("Workflow report:", start)
        monitorThread.join()
        return self.failed

class AgoraWorkflow:

    # Default paths (in case not set in the configuration file)
    defaultPaths = {
        'ancGenesData': 'ancGenes/%(filt)s/ancGenes.%(name)s.list.bz2',
        'ancGenesLog': 'ancGenes/%(filt)s.log',
        'filteredBlocksData': 'filtBlocks/%(filt)s/blocks.%(name)s.list.bz2',
        'filteredBlocksLog': 'filtBlocks/%(filt)s/log',
        'geneTreesWithAncNames': 'GeneTreeForest.withAncGenes.nhx.bz2',
        'pairwiseOutput': 'pairwise/pairs-%(filt)s/%(name)s.list.bz2',
        'pairwiseLog': 'pairwise/pairs-%(filt)s/log',
        'adjacenciesOutput': 'pairwise/adjacencies-%(filt)s/%(name)s.list.bz2',
        'adjacenciesDebug': 'pairwise/adjacencies-%(filt)s/%(name)s.log.bz2',
        'adjacenciesLog': 'pairwise/adjacencies-%(filt)s/log',
        'ancBlocks': 'ancBlocks/%(method)s/blocks.%(name)s.list.bz2',
        'ancGraphs': 'ancBlocks/%(method)s/graph.%(name)s.txt.bz2',
        'ancLog': 'ancBlocks/%(method)s/log',
        'ancGenomesOutput': 'ancGenomes/%(method)s/ancGenome.%(name)s.list.bz2',
        'ancGenomesLog': 'ancGenomes/%(method)s/log',
    }
    inputParams = ["speciesTree", "geneTrees", "genes"]
    allAncGenesName = "all"


    def __init__(self, defaultRoot, defaultExtantSpeciesFilter, scriptDir, files):
        self.defaultRoot = defaultRoot
        self.defaultExtantSpeciesFilter = ["-extantSpeciesFilter=" + defaultExtantSpeciesFilter] if defaultExtantSpeciesFilter else []
        self.tasklist = TaskList()
        self.scriptDir = scriptDir
        self.files = files
        self.interm = {}
        self.refMethod = {}
        self.ancBlocksAsAncGenes = False
        self.ancGenesTaskName = "ancgenes"
        self.ancGenesFileEntryName = "ancGenesData"
        self.pairwiseFileEntryName = "pairwiseOutput"
        self.allAncGenesPath = self.files["ancGenesData"] % {"filt": self.allAncGenesName, "name": "%s"}
        self.selectionPool = []
        # With agora-*.py, people may use %s instead of %(name)s
        if '%(name)s' not in files['genes']:
            files['genes'] = files['genes'].replace('%s', '%(name)s')

    def addDummy(self, taskFullName, dependencies=[]):
        return self.tasklist.addTask(taskFullName, dependencies, (None, None, None, False))

    def addAncGenesGenerationAnalysis(self, launch=True):
        taskFullName = ("ancgenes", self.allAncGenesName)
        if launch:
            return self.tasklist.addTask(
                taskFullName,
                [],
                (
                    [
                        os.path.join(self.scriptDir, "ALL.extractGeneFamilies.py"),
                        self.files["speciesTree"],
                        self.files["geneTrees"],
                        "-OUT.ancGenesFiles=" + self.files["ancGenesData"] % {"filt": self.allAncGenesName, "name": "%s"},
                    ],
                    self.files["geneTreesWithAncNames"],
                    self.files["ancGenesLog"] % {"filt": "ancGenes"},
                    launch,
                )
            )
        else:
            return self.addDummy(taskFullName)

    # FIXME: both this and the callers implement their own naming scheme (size-0.9-1.1). Risk is that they diverge
    def addAncGenesFilterAnalysis(self, methodName, params, ancestor=None, launch=True):

        taskName = "-".join([methodName] + params)

        if self.ancBlocksAsAncGenes:
            taskName = self.blocksName + "-" + taskName
            inputName = self.blocksName + "-" + self.allAncGenesName
            scriptTemplate = "ALL.filterContigs-%s.py"
            inputPath = self.files["filteredBlocksData"] % {"filt": self.blocksName + "-" + self.allAncGenesName, "name": "%s"}
            outputPath = self.files["filteredBlocksData"] % {"filt": taskName, "name": "%s"}
            logPath = self.files["filteredBlocksLog"] % {"filt": taskName}
        else:
            scriptTemplate = "ALL.filterGeneFamilies-%s.py"
            inputName = self.allAncGenesName
            inputPath = self.files["ancGenesData"] % {"filt": self.allAncGenesName, "name": "%s"}
            outputPath = self.files["ancGenesData"] % {"filt": methodName + "-%s", "name": "%s"}
            logPath = self.files["ancGenesLog"] % {"filt": taskName}

        return self.tasklist.addTask(
            (self.ancGenesTaskName, taskName),
            [(self.ancGenesTaskName, inputName)],
            (
                [
                    os.path.join(self.scriptDir, scriptTemplate  % methodName),
                    self.files["speciesTree"],
                    ancestor or self.defaultRoot,
                    inputPath,
                    outputPath,
                ] + params,
                None,
                logPath,
                launch,
            )
        )

    def addPairwiseAnalysis(self, ancGenesName, methodName=None, params=[], ancestor=None, launch=True):

        if self.ancBlocksAsAncGenes:
            if methodName is None:
                methodName = "conservedAdjacencies"
            ancGenesName = self.blocksName + "-" + ancGenesName
            params.append("-iniAncGenesFiles=" + self.files["ancGenesData"] % {"filt": self.allAncGenesName, "name": "%s"})
            params.append("-LOG.pairwise=" + self.files["adjacenciesDebug"] % {"filt": ancGenesName, "name": "%s"})
        else:
            if methodName is None:
                methodName = "conservedPairs"

        return self.tasklist.addTask(
            ("pairwise", self.ancGenesTaskName + "-" + ancGenesName),
            [(self.ancGenesTaskName, ancGenesName)],
            (
                [
                    os.path.join(self.scriptDir, "buildSynteny.pairwise-%s.py" % methodName),
                    self.files["speciesTree"],
                    ancestor or self.defaultRoot,
                    "-ancGenesFiles=" + self.files[self.ancGenesFileEntryName] % {"filt": ancGenesName, "name": "%s"},
                    "-genesFiles=" + self.files["genes"] % {"name": "%s"},
                    "-OUT.pairwise=" + self.files[self.pairwiseFileEntryName] % {"filt": ancGenesName, "name": "%s"}
                ] + self.defaultExtantSpeciesFilter + params,
                None,
                self.files[self.pairwiseFileEntryName.replace("Output", "Log")] % {"filt": ancGenesName},
                launch,
            ),
            methodName == "conservedAdjacencies",  # Only this one is multithreaded
        )

    def addIntegrationAnalysis(self, methodName, params, pairwiseName, taskName=None, inputName=None, outputName=None, ancestor=None, launch=True):

        if taskName is None:
            taskName = methodName

        if taskName.endswith("/"):
            taskName = taskName[:-1]
        elif pairwiseName:
            taskName = taskName + "-" + pairwiseName

        if methodName == "denovo":
            newMethod = taskName
            if self.ancBlocksAsAncGenes:
                newMethod = self.blocksName + "." + newMethod
            if not ancestor:
                ancestor = self.defaultRoot
            self.refMethod[newMethod] = (newMethod, ancestor)
        else:
            if inputName:
                self.prevMethod = self.interm[inputName]
            if taskName.startswith("/"):
                newMethod = taskName[1:]
            else:
                newMethod = self.prevMethod + "." + taskName
            if not ancestor:
                ancestor = self.refMethod[self.prevMethod][1]
            if methodName == "fillin":
                self.refMethod[newMethod] = (newMethod, ancestor)
            else:
                self.refMethod[newMethod] = self.refMethod[self.prevMethod]

        if outputName:
            self.interm[outputName] = newMethod

        # task parameters
        args = [
                os.path.join(self.scriptDir, "buildSynteny.integr-%s.py" % methodName),
                self.files["speciesTree"],
                ancestor,
        ] + params

        if methodName == "publish":
            # "publish" is not an integration method
            args[0] = os.path.join(self.scriptDir, "convert.ancGenomes.blocks-to-genes.py")
            args.append("-OUT.ancGenomes=" + self.files["ancGenomesOutput"] % {"method": newMethod, "name": "%s"})
            logfile = self.files["ancGenomesLog"]
        else:
            args.append("-OUT.ancBlocks=" + self.files["ancBlocks"] % {"method": newMethod, "name": "%s"})
            logfile = self.files["ancLog"]

        dep = []
        if pairwiseName is not None:
            if self.ancBlocksAsAncGenes:
                pairwiseName = self.blocksName + "-" + pairwiseName
            dep.append(("pairwise", self.ancGenesTaskName + "-" + pairwiseName))
            args.append(self.files[self.pairwiseFileEntryName] % {"filt": pairwiseName, "name": "%s"})

        if methodName in ["denovo", "scaffolds", "publish"]:
            args.append("-ancGenesFiles=" + self.allAncGenesPath)

        # No input data to consider for the denovo method
        if methodName != "denovo":
            dep.append(("integr", self.prevMethod))
            args.append("-IN.ancBlocks=" + self.files["ancBlocks"] % {"method": self.prevMethod, "name": "%s"})

        if methodName == "insertion":
            # The script needs singleton reference for "insertion"
            dep.append(("integr", self.refMethod[newMethod][0]))
            args.append("-REF.ancBlocks=" + self.files["ancBlocks"] % {"method": self.refMethod[newMethod][0], "name": "%s"})

        if methodName == "scaffolds":
            args.append("-genesFiles=" + self.files["genes"] % {"name": "%s"})
            args.extend(self.defaultExtantSpeciesFilter)

        if methodName not in ["copy", "publish"]:
            args.append("-LOG.ancGraph=" + self.files["ancGraphs"] % {"method": newMethod, "name": "%s"})

        # Most of the methods are multithreaded
        multithreaded = methodName not in ["copy"]

        # The publish method doesn't generate ancBlocks and can't be used as an input method
        if methodName != "publish":
            self.prevMethod = newMethod

        return self.tasklist.addTask(
            ("integr" if methodName != "publish" else "publish", newMethod),
            dep,
            (
                args,
                None,
                logfile % {"method": newMethod},
                launch,
            ),
            multithreaded,
        )


    def reconstructionPassWithAncGenesFiltering(self, filteringMethod, filteringParams, ancestor=None, launch=True):
        filteringParams = map(str, filteringParams)
        filteredAncGenesDirName = filteringMethod + "-" + "-".join(filteringParams)
        self.addAncGenesFilterAnalysis(filteringMethod, filteringParams, ancestor=ancestor, launch=launch)
        # Don't run twice
        if self.ancBlocksAsAncGenes:
            pairwiseTaskName = ("pairwise", self.ancGenesTaskName + "-" + self.blocksName + "-" + self.allAncGenesName)
        else:
            pairwiseTaskName = ("pairwise", self.ancGenesTaskName + "-" + self.allAncGenesName)
        if pairwiseTaskName not in self.tasklist.dic:
            self.addPairwiseAnalysis(self.allAncGenesName, ancestor=ancestor, launch=launch)
        self.addPairwiseAnalysis(filteredAncGenesDirName, ancestor=ancestor, launch=launch)
        self.addIntegrationAnalysis("denovo", [], filteredAncGenesDirName, ancestor=ancestor, launch=launch)
        self.addIntegrationAnalysis("fillin", [], self.allAncGenesName, ancestor=ancestor, launch=launch)
        self.addIntegrationAnalysis("fusion", ["+onlySingletons"], self.allAncGenesName, ancestor=ancestor, launch=launch)
        self.addIntegrationAnalysis("insertion", [], self.allAncGenesName, ancestor=ancestor, launch=launch)


    def publishGenome(self, outputName=None, inputName=None, ancestor=None, launch=True):

        if inputName:
            self.prevMethod = self.interm[inputName]

        if outputName is None:
            outputName = self.prevMethod

        if not ancestor:
            ancestor = self.refMethod[self.prevMethod][1]

        # task parameters
        args = [
                os.path.join(self.scriptDir, "convert.ancGenomes.blocks-to-genes.py"),
                self.files["speciesTree"],
                ancestor,
                "-IN.ancBlocks=" + self.files["ancBlocks"] % {"method": self.prevMethod, "name": "%s"},
                "-ancGenesFiles=" + self.files["ancGenesData"] % {"filt": self.allAncGenesName, "name": "%s"},
                "-OUT.ancGenomes=" + self.files["ancGenomesOutput"] % {"method": outputName, "name": "%s"},
        ]

        return self.tasklist.addTask(
            ("conversion", outputName),
            [("integr", self.prevMethod)],
            (
                args,
                None,
                self.files["ancGenomesLog"] % {"method": outputName},
                launch,
            ),
            True,
        )


    def useBlocksAsAncGenes(self, ancestor=None, launch=True):
        self.ancBlocksAsAncGenes = True
        self.ancGenesTaskName = "ancblocks"
        self.ancGenesFileEntryName = "filteredBlocksData"
        self.pairwiseFileEntryName = "adjacenciesOutput"
        self.blocksName = self.prevMethod
        self.allAncGenesPath = self.files["ancBlocks"] % {"method": self.blocksName, "name": "%s"}

        return self.tasklist.addTask(
            ("ancblocks" , self.blocksName + "-" + self.allAncGenesName),
            [("integr", self.blocksName)],
            (
                [
                    os.path.join(self.scriptDir, "buildSynteny.integr-copy.py"),
                    self.files["speciesTree"],
                    ancestor or self.defaultRoot,
                    "-IN.ancBlocks=" + self.files["ancBlocks"] % {"method": self.blocksName, "name": "%s"},
                    "-OUT.ancBlocks=" + self.files["filteredBlocksData"] % {"filt": self.blocksName + "-" + self.allAncGenesName, "name": "%s"}
                ],
                None,
                self.files["filteredBlocksLog"] % {"filt": self.blocksName + "-" + self.allAncGenesName, "name": "%s"},
                launch,
            )
        )

    def convertToRealAncGenes(self, taskName=None, inputName=None, outputName=None, ancestor=None, launch=True):

        if taskName is None:
            taskName = "asAncGenes"

        if inputName:
            self.prevMethod = self.interm[inputName]

        if taskName.startswith("/"):
            newMethod = taskName[1:]
        else:
            newMethod = self.prevMethod + "." + taskName

        if not ancestor:
            ancestor = self.refMethod[self.prevMethod][1]

        self.refMethod[newMethod] = self.refMethod[self.prevMethod]

        if outputName:
            self.interm[outputName] = newMethod

        # task parameters
        deps = [("integr", self.prevMethod)]
        args = [
                os.path.join(self.scriptDir, "convert.ancGenomes.blocks-of-blocks-to-blocks-of-ancGenes.py"),
                self.files["speciesTree"],
                ancestor,
                "-IN.scaffoldsFile=" + self.files["ancBlocks"] % {"method": self.prevMethod, "name": "%s"},
                "-IN.contigsFile=" + self.files["filteredBlocksData"] % {"filt": self.blocksName + "-" + self.allAncGenesName, "name": "%s"},
                "-OUT.ancBlocksFile=" + self.files["ancBlocks"] % {"method": newMethod, "name": "%s"},
        ]

        self.prevMethod = newMethod

        return self.tasklist.addTask(
            ("integr", newMethod),
            deps,
            (
                args,
                None,
                self.files["ancLog"] % {"method": newMethod},
                launch,
            ),
            True,
        )


    def revertToRealAncGenes(self):
        self.ancBlocksAsAncGenes = False
        self.prevMethod = self.blocksName
        self.ancGenesTaskName = "ancgenes"
        self.ancGenesFileEntryName = "ancGenesData"
        self.pairwiseFileEntryName = "pairwiseOutput"
        self.allAncGenesPath = self.files["ancGenesData"] % {"filt": self.allAncGenesName, "name": "%s"}

    def markForSelection(self):
        self.selectionPool.append( ("integr", self.prevMethod) )

    def addSelectionAnalysis(self, taskName=None, outputName=None, ancestor=None, launch=True):

        if taskName:
            newMethod = taskName
        else:
            newMethod = "selection"

        if newMethod.startswith("/"):
            newMethod = newMethod[1:]
        else:
            newMethod = self.prevMethod + "." + newMethod

        if outputName:
            self.interm[outputName] = newMethod

        if not ancestor:
            ancestor = self.refMethod[self.prevMethod][1]

        self.refMethod[newMethod] = self.refMethod[self.prevMethod]

        # task parameters
        args = [
                os.path.join(self.scriptDir, "ALL.selectBestReconstruction.py"),
                self.files["speciesTree"],
                ancestor,
                "-OUT.ancBlocks=" + self.files["ancBlocks"] % {"method": newMethod, "name": "%s"},
        ] + [self.files["ancBlocks"] % {"method": method, "name": "%s"} for (_, method) in self.selectionPool]

        task = self.tasklist.addTask(
            ("integr", newMethod),
            self.selectionPool,
            (
                args,
                None,
                self.files["ancLog"] % {"method": newMethod},
                launch,
            ),
        )
        self.prevMethod = newMethod
        self.selectionPool = []

        return task
