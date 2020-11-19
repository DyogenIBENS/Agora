#! /usr/bin/env python
# -*- coding: utf-8 -*-
# Agora v1.5
# python 2.7
# Copyright © 2020 IBENS/Dyogen and EMBL-European Bioinformatics Institute : Matthieu MUFFATO, Alexandra LOUIS, Thi Thuy Nga NGUYEN, Hugues ROEST CROLLIUS
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

    def getMemoryUsage(self, pid):
        try:
            proc = psutil.Process(pid)
            return proc.memory_full_info().uss
        except (psutil.NoSuchProcess, IOError):
            return None

    def getRecursiveMemoryUsage(self, pid):
        try:
            proc = psutil.Process(pid)
            total_mem = proc.memory_full_info().uss
            children = proc.children(recursive=True)
        except (psutil.NoSuchProcess, IOError):
            return None
        for subproc in children:
            try:
                # mem = subproc.memory_info().rss
                # Slower but more accurate
                mem = subproc.memory_full_info().uss
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
        mem = max(ru.ru_maxrss * 1024, self.memusage.pop(os.getpid()))
        self.memlock.release()
        print intro, "%g sec CPU time / %g sec elapsed = %g%% CPU usage, %g MB RAM" % (ru.ru_utime + ru.ru_stime, elapsed, 100. * (ru.ru_utime + ru.ru_stime) / elapsed, mem / 1024. / 1024.)

    # Launch program function
    def goLaunch(self, i, args, out, log):
        start = time.time()
        stdout = myFile.openFile(out or os.devnull, "w")
        stderr = myFile.openFile(log or os.devnull, "w")
        # stderr must have a fileno, so must be a regular file (not a .bz2 etc)
        # stdout can be anything, incl. a .bz2
        p = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=stderr)
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
        'geneTreesWithAncNames': 'GeneTreeForest.withAncGenes.nhx.bz2',
        'pairwiseOutput': 'pairwise/pairs-%(filt)s/%(name)s.list.bz2',
        'pairwiseLog': 'pairwise/pairs-%(filt)s/log',
        'integrBlocks': 'integrDiags/%(method)s/diags.%(name)s.list.bz2',
        'integrOutput': 'integrDiags/%(method)s/graph.%(name)s.log.bz2',
        'integrLog': 'integrDiags/%(method)s/log',
        'ancGenomesOutput': 'ancGenomes/%(method)s/ancGenome.%(name)s.list.bz2',
        'ancGenomesLog': 'ancGenomes/%(method)s/log',
    }
    inputParams = ["speciesTree", "geneTrees", "genes"]


    def __init__(self, defaultRoot, defaultExtantSpeciesFilter, scriptDir, files):
        self.defaultRoot = defaultRoot
        self.defaultExtantSpeciesFilter = ["-extantSpeciesFilter=" + defaultExtantSpeciesFilter] if defaultExtantSpeciesFilter else []
        self.tasklist = TaskList()
        self.scriptDir = scriptDir
        self.files = files
        self.allAncGenesTaskName = "all"
        self.allAncGenesDirName = "all"
        self.interm = {}
        self.refMethod = {}
        # With agora-*.py, people may use %s instead of $(name)s
        if '%(name)s' not in files['genes']:
            files['genes'] = files['genes'].replace('%s', '%(name)s')

    def addDummy(self, taskFullName, dependencies=[]):
        return self.tasklist.addTask(taskFullName, dependencies, (None, None, None, False))

    def addAncGenesGenerationAnalysis(self, launch=True):
        taskFullName = ("ancgenes", self.allAncGenesTaskName)
        if launch:
            return self.tasklist.addTask(
                taskFullName,
                [],
                (
                    [
                        os.path.join(self.scriptDir, "ALL.extractGeneFamilies.py"),
                        self.files["speciesTree"],
                        self.files["geneTrees"],
                        "-OUT.ancGenesFiles=" + self.files["ancGenesData"] % {"filt": self.allAncGenesDirName, "name": "%s"},
                    ],
                    self.files["geneTreesWithAncNames"],
                    self.files["ancGenesLog"] % {"filt": "ancGenes"},
                    launch,
                )
            )
        else:
            return self.addDummy(taskFullName)

    def addAncGenesFilterAnalysis(self, taskName, methodName, params, dirnameTemplate, ancestor=None, launch=True):
        return self.tasklist.addTask(
            ("ancgenes", taskName),
            [("ancgenes", self.allAncGenesTaskName)],
            (
                [
                    os.path.join(self.scriptDir, "ALL.filterGeneFamilies-%s.py" % methodName),
                    self.files["speciesTree"],
                    ancestor or self.defaultRoot,
                    self.files["ancGenesData"] % {"filt": self.allAncGenesDirName, "name": "%s"},
                    self.files["ancGenesData"] % {"filt": dirnameTemplate, "name": "%s"}
                ] + params,
                None,
                self.files["ancGenesLog"] % {"filt": taskName},
                launch,
            )
        )

    def addPairwiseAnalysis(self, taskName, methodName="conservedPairs", params=[], ancestor=None, launch=True):
        return self.tasklist.addTask(
            ("pairwise", taskName),
            [("ancgenes",taskName)],
            (
                [
                    os.path.join(self.scriptDir, "buildSynteny.pairwise-%s.py" % methodName),
                    self.files["speciesTree"],
                    ancestor or self.defaultRoot,
                    "-ancGenesFiles=" + self.files["ancGenesData"] % {"filt": taskName, "name": "%s"},
                    "-genesFiles=" + self.files["genes"] % {"name": "%s"},
                    "-OUT.pairwise=" + self.files["pairwiseOutput"] % {"filt": taskName, "name": "%s"}
                ] + self.defaultExtantSpeciesFilter + params,
                None,
                self.files["pairwiseLog"] % {"filt": taskName},
                launch,
            )
        )

    def addIntegrationAnalysis(self, methodName, params, pairwiseName, taskName=None, inputName=None, outputName=None, ancestor=None, launch=True):

        if taskName is None:
            taskName = methodName

        if pairwiseName:
            taskName = taskName[:-1] if taskName.endswith("/") else (taskName+ "-" + pairwiseName)
        elif taskName.endswith("/"):
            taskName = taskName[:-1]

        if methodName == "denovo":
            newMethod = taskName
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
            args[0] = os.path.join(self.scriptDir, "convert.ancGenomes.diags-genes.py")
            args.append("-OUT.ancGenomes=" + self.files["ancGenomesOutput"] % {"method": newMethod, "name": "%s"})
            logfile = self.files["ancGenomesLog"]
        else:
            args.append("-OUT.ancDiags=" + self.files["integrBlocks"] % {"method": newMethod, "name": "%s"})
            logfile = self.files["integrLog"]

        dep = []
        if pairwiseName is not None:
            dep.append(("pairwise", pairwiseName))
            args.append(self.files["pairwiseOutput"] % {"filt": pairwiseName, "name": "%s"})

        if methodName in ["denovo", "groups", "publish"]:
            args.append("-ancGenesFiles=" + self.files["ancGenesData"] % {"filt": "all", "name": "%s"})

        # No input data to consider for the denovo method
        if methodName != "denovo":
            dep.append(("integr", self.prevMethod))
            args.append("-IN.ancDiags=" + self.files["integrBlocks"] % {"method": self.prevMethod, "name": "%s"})

        if methodName == "halfinsert":
            # The script needs singleton reference for "halfinsert"
            dep.append(("integr", self.refMethod[newMethod][0]))
            args.append("-REF.ancDiags=" + self.files["integrBlocks"] % {"method": self.refMethod[newMethod][0], "name": "%s"})

        if methodName == "groups":
            args.append("-genesFiles=" + self.files["genes"] % {"name": "%s"})
            args.extend(self.defaultExtantSpeciesFilter)

        if methodName not in ["copy", "publish"]:
            args.append("-LOG.ancGraph=" + self.files["integrOutput"] % {"method": newMethod, "name": "%s"})

        # Most of the methods are multithreaded
        multithreaded = methodName not in ["copy"]

        # The publish method doesn't generate integrDiags and can't be used as an input method
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
