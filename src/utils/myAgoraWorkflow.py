#! /usr/bin/env python
# -*- coding: utf-8 -*-
# AGORA v3.0
# python 2.7
# Copyright © 2006-2021 IBENS/Dyogen, 2020-2021 EMBL-European Bioinformatics Institute, 2021 Genome Research Ltd : Matthieu MUFFATO, Alexandra LOUIS, Thi Thuy Nga NGUYEN, Hugues ROEST CROLLIUS
# mail : agora@bio.ens.psl.eu
# This is free software; you may copy, modify and/or distribute this work under the terms of the GNU General Public License, version 3 or later and the CeCiLL v2 license in France

import collections
import itertools
import json
import multiprocessing
import os
import resource
import subprocess
import sys
import time
import threading

import psutil

from . import myFile

# A command that will be run. args represents the entire command-line, incl. the executable
Command = collections.namedtuple("Command", ['args', 'out', 'log'])
Task = collections.namedtuple("Task", ['dependencies', 'command', 'multithreaded'])


# Managing the list of programs to launch and their dependencies
#################################################################
class TaskList():

    rusage_unit = 1 if sys.platform == "darwin" else 1024
    status_filename = '.agora'

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
        print("digraph", "{", file=fh)
        for (name, taskId) in self.dic.items():
            print('%d [label="%s"]' % (taskId, "/".join(name)), file=fh)
            for dep in self.list[taskId].dependencies:
                print('%d -> %d' % (dep, taskId), file=fh)
        print("}", file=fh)

    def addTask(self, name, dep, command, multithreaded=False):
        taskId = len(self.list)
        print("New task", taskId, name)
        print(dep)
        print(command)
        self.list.append(Task(set(self.dic[x] for x in dep), command, multithreaded))
        if name in self.dic:
            if name + ("1",) in self.dic:
                self.list[self.dic[name]].dependencies.add(taskId)
                for i in itertools.count(2):
                    newName = name + (str(i),)
                    if newName not in self.dic:
                        print("! Name clash ! Renamed to", newName, "under a collector")
                        self.dic[newName] = taskId
                        break
            else:
                print("! Name clash ! Introducing a collector task")
                collectorId = self.dic[name]
                self.list.append(self.list[collectorId])
                self.list[collectorId] = Task(set([taskId, taskId + 1]), Command(None, None, None), False)
                self.dic[name + ("1",)] = taskId + 1
                self.dic[name + ("2",)] = taskId
        else:
            self.dic[name] = taskId
        print()
        return taskId

    def removeDep(self, i):
        for t in self.list:
            t.dependencies.discard(i)

    def getAvailable(self):
        tmp = [i for (i, task) in enumerate(self.list) if len(task.dependencies) == 0]
        print("Available tasks:", tmp)
        if len(tmp) > 0:
            taskId = tmp[0]
            # Add something (None) to the dependency list so that the task is never selected again
            # This way the task numbers remain the same
            self.list[taskId].dependencies.add(None)
            return taskId
        else:
            return None

    # Waiting for a task to finish
    def joinNext(self):
        print("Waiting ...")
        sys.stdout.flush()
        (i, r) = self.queue.get()
        print("task", i, "is now finished (status %d)" % r)
        self.proc.pop(i).join()
        if r == 0:
            self.removeDep(i)
            self.completed += 1
        else:
            self.failed += 1
            print(">", "Inspect", self.list[i].command.log, "for more information", file=sys.stderr)
        self.nrun -= self.nthreads.pop(i)

    def getJsonPath(self, i):
        command = self.list[i].command
        if command.log or command.out:
            return (command.log or command.out) + self.status_filename
        return None

    def getJsonPayload(self, i):
        command = self.list[i].command
        return {
                'command': command.args,
                'out': command.out,
                'log': command.log,
                }

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
            for pid in list(self.memusage.keys()):
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
        print(intro, "%g sec CPU time / %g sec elapsed = %g%% CPU usage, %g MB RAM" % (ru.ru_utime + ru.ru_stime, elapsed, 100. * (ru.ru_utime + ru.ru_stime) / elapsed, mem / 1024. / 1024.))

    # Launch program function
    def goLaunch(self, i, command, status_file):
        start = time.time()
        if status_file and os.path.exists(status_file):
            os.remove(status_file)
        stdout = myFile.openFile(command.out or os.devnull, "wb")
        stderr = myFile.openFile(command.log, "wb") if command.log else subprocess.DEVNULL
        # stderr must have a fileno, so must be a regular file (not a .bz2 etc)
        # stdout can be anything, incl. a .bz2
        try:
            p = subprocess.Popen(command.args, stdout=subprocess.PIPE, stderr=stderr)
        except Exception as e:
            stdout.close()
            stderr.close()
            print("task %d could not start:" % i, e)
            time.sleep(5)
            self.queue.put((i, -1))
            # FIXME: then it hangs in multiprocessing/managers (checking self.memusage)
            return

        # This is where the .bz2 compression would happen
        for l in p.stdout:
            stdout.write(l)
        r = p.wait()
        for l in p.stdout:
            stdout.write(l)
        stdout.close()
        if command.log:
            stderr.close()
        if r == 0 and status_file:
            report = self.getJsonPayload(i)
            with open(status_file, 'w') as fh:
                json.dump(report, fh)
        self.printCPUUsageStats("task %d report:" % i, start)
        time.sleep(5)
        self.queue.put((i, r))

    # Launching tasks in multiple threads
    def runAll(self, nbThreads, sequential, forceRerun):
        start = time.time()

        self.memusage[os.getpid()] = 0
        monitorThread = threading.Thread(target=self.memoryMonitor)
        monitorThread.start()

        # Queue
        while (self.completed + self.failed) < len(self.list):

            print("Status: %d to do, %d running, %d done, %d failed -- %d total" % \
                    (len(self.list)-len(self.proc)-self.completed-self.failed, len(self.proc), self.completed, self.failed, len(self.list)))

            if (self.nrun == nbThreads) or (sequential and self.nrun):
                self.joinNext()
            else:
                taskId = self.getAvailable()
                if taskId is None:
                    if self.nrun == 0:
                        print("Workflow stopped because of failures")
                        break
                    self.joinNext()
                else:
                    command = self.list[taskId].command
                    # Check if this step has already run
                    launch = True
                    if not command.args:
                        print("Dummy task")
                        launch = False
                    status_file = self.getJsonPath(taskId)
                    if forceRerun:
                        print("'forceRerun' option given, not checking the control file")
                    elif status_file:
                        print("Control file", status_file, end=' ')
                        if os.path.exists(status_file):
                            with open(status_file, 'r') as fh:
                                j = json.load(fh)
                            if j == self.getJsonPayload(taskId):
                                launch = False
                                print("present - same parameters")
                            else:
                                print("present - different parameters")
                        else:
                            print("missing")
                    else:
                        print("No control file could be identified")
                    if launch:
                        print("Launching task", taskId, command.args, ">", command.out, "2>", command.log)
                        if self.list[taskId].multithreaded:
                            self.nthreads[taskId] = nbThreads - self.nrun
                            # Creating a new list so that the original list remains available for the Json dump
                            command = Command(command.args + ["-nbThreads=%d" % self.nthreads[taskId]], command.out, command.log)
                            print("Using", self.nthreads[taskId], "threads")
                        else:
                            self.nthreads[taskId] = 1
                        self.proc[taskId] = multiprocessing.Process(target=self.goLaunch, args=(taskId, command, status_file))
                        self.proc[taskId].start()
                        self.memusage[self.proc[taskId].pid] = 0
                        self.nrun += self.nthreads[taskId]
                    else:
                        print("Skipping task", taskId)
                        self.removeDep(taskId)
                        self.completed += 1

        assert self.nrun == 0
        if not self.failed:
            print("Workflow complete")
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
        return self.tasklist.addTask(taskFullName, dependencies, Command(None, None, None), False)

    def addAncGenesGenerationAnalysis(self):
            taskFullName = ("ancgenes", self.allAncGenesName)
            return self.tasklist.addTask(
                taskFullName,
                [],
                Command(
                    [
                        os.path.join(self.scriptDir, "ALL.extractGeneFamilies.py"),
                        self.files["speciesTree"],
                        self.files["geneTrees"],
                        "-OUT.ancGenesFiles=" + self.files["ancGenesData"] % {"filt": self.allAncGenesName, "name": "%s"},
                    ],
                    self.files["geneTreesWithAncNames"],
                    self.files["ancGenesLog"] % {"filt": "ancGenes"},
                )
            )

    # FIXME: both this and the callers implement their own naming scheme (size-0.9-1.1). Risk is that they diverge
    def addAncGenesFilterAnalysis(self, methodName, params, ancestor=None):

        taskName = "-".join([methodName] + params)

        if self.ancBlocksAsAncGenes:
            taskName = self.blocksName + "-" + taskName
            inputName = self.blocksName + "-" + self.allAncGenesName
            scriptTemplate = "ALL.filterBlocks-%s.py"
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
            Command(
                [
                    os.path.join(self.scriptDir, scriptTemplate  % methodName),
                    self.files["speciesTree"],
                    ancestor or self.defaultRoot,
                    inputPath,
                    outputPath,
                ] + params,
                None,
                logPath,
            )
        )

    def addPairwiseAnalysis(self, ancGenesName, methodName=None, params=[], ancestor=None):

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
            Command(
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
            ),
            methodName == "conservedAdjacencies",  # Only this one is multithreaded
        )

    def addIntegrationAnalysis(self, methodName, params, pairwiseName, taskName=None, inputName=None, outputName=None, ancestor=None):

        # Legacy interface, still used in .ini-based workflows
        if methodName == "publish":
            return self.publishGenome(outputName=taskName, inputName=inputName, ancestor=ancestor)

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

        args.append("-OUT.ancBlocks=" + self.files["ancBlocks"] % {"method": newMethod, "name": "%s"})
        logfile = self.files["ancLog"]

        dep = []
        if pairwiseName is not None:
            if self.ancBlocksAsAncGenes:
                pairwiseName = self.blocksName + "-" + pairwiseName
            dep.append(("pairwise", self.ancGenesTaskName + "-" + pairwiseName))
            args.append(self.files[self.pairwiseFileEntryName] % {"filt": pairwiseName, "name": "%s"})

        if methodName in ["denovo", "scaffolds"]:
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

        if methodName != "copy":
            args.append("-LOG.ancGraph=" + self.files["ancGraphs"] % {"method": newMethod, "name": "%s"})

        # Most of the methods are multithreaded
        multithreaded = methodName not in ["copy"]

        self.prevMethod = newMethod

        return self.tasklist.addTask(
            ("integr", newMethod),
            dep,
            Command(
                args,
                None,
                logfile % {"method": newMethod},
            ),
            multithreaded,
        )


    def reconstructionPassWithAncGenesFiltering(self, filteringMethod, filteringParams, ancestor=None):
        filteringParams = list(map(str, filteringParams))
        filteredAncGenesDirName = filteringMethod + "-" + "-".join(filteringParams)
        self.addAncGenesFilterAnalysis(filteringMethod, filteringParams, ancestor=ancestor)
        # Don't run twice
        if self.ancBlocksAsAncGenes:
            pairwiseTaskName = ("pairwise", self.ancGenesTaskName + "-" + self.blocksName + "-" + self.allAncGenesName)
        else:
            pairwiseTaskName = ("pairwise", self.ancGenesTaskName + "-" + self.allAncGenesName)
        if pairwiseTaskName not in self.tasklist.dic:
            self.addPairwiseAnalysis(self.allAncGenesName, ancestor=ancestor)
        self.addPairwiseAnalysis(filteredAncGenesDirName, ancestor=ancestor)
        self.addIntegrationAnalysis("denovo", [], filteredAncGenesDirName, ancestor=ancestor)
        self.addIntegrationAnalysis("fillin", [], self.allAncGenesName, ancestor=ancestor)
        self.addIntegrationAnalysis("fusion", ["+onlySingletons"], self.allAncGenesName, ancestor=ancestor)
        self.addIntegrationAnalysis("insertion", [], self.allAncGenesName, ancestor=ancestor)


    def publishGenome(self, outputName=None, inputName=None, ancestor=None):

        if inputName:
            self.prevMethod = self.interm[inputName]

        if outputName is None:
            outputName = self.prevMethod
        outputName = outputName.strip("/")

        if not ancestor:
            ancestor = self.refMethod[self.prevMethod][1]

        # task parameters
        args = [
                os.path.join(self.scriptDir, "convert.ancGenomes.blocks-to-genes.py"),
                self.files["speciesTree"],
                ancestor,
                "+orderBySize",
                "-IN.ancBlocks=" + self.files["ancBlocks"] % {"method": self.prevMethod, "name": "%s"},
                "-ancGenesFiles=" + self.files["ancGenesData"] % {"filt": self.allAncGenesName, "name": "%s"},
                "-OUT.ancGenomes=" + self.files["ancGenomesOutput"] % {"method": outputName, "name": "%s"},
        ]

        return self.tasklist.addTask(
            ("conversion", outputName),
            [("integr", self.prevMethod)],
            Command(
                args,
                None,
                self.files["ancGenomesLog"] % {"method": outputName},
            ),
            True,
        )


    def useBlocksAsAncGenes(self, ancestor=None):
        self.ancBlocksAsAncGenes = True
        self.ancGenesTaskName = "ancblocks"
        self.ancGenesFileEntryName = "filteredBlocksData"
        self.pairwiseFileEntryName = "adjacenciesOutput"
        self.blocksName = self.prevMethod
        self.allAncGenesPath = self.files["ancBlocks"] % {"method": self.blocksName, "name": "%s"}

        return self.tasklist.addTask(
            ("ancblocks" , self.blocksName + "-" + self.allAncGenesName),
            [("integr", self.blocksName)],
            Command(
                [
                    os.path.join(self.scriptDir, "buildSynteny.integr-copy.py"),
                    self.files["speciesTree"],
                    ancestor or self.defaultRoot,
                    "-IN.ancBlocks=" + self.files["ancBlocks"] % {"method": self.blocksName, "name": "%s"},
                    "-OUT.ancBlocks=" + self.files["filteredBlocksData"] % {"filt": self.blocksName + "-" + self.allAncGenesName, "name": "%s"}
                ],
                None,
                self.files["filteredBlocksLog"] % {"filt": self.blocksName + "-" + self.allAncGenesName, "name": "%s"},
            )
        )

    def convertToRealAncGenes(self, taskName=None, inputName=None, outputName=None, ancestor=None):

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
                "-IN.blocksBlocksFile=" + self.files["ancBlocks"] % {"method": self.prevMethod, "name": "%s"},
                "-IN.blocksGenesFile=" + self.files["filteredBlocksData"] % {"filt": self.blocksName + "-" + self.allAncGenesName, "name": "%s"},
                "-OUT.ancBlocksFile=" + self.files["ancBlocks"] % {"method": newMethod, "name": "%s"},
        ]

        self.prevMethod = newMethod

        return self.tasklist.addTask(
            ("integr", newMethod),
            deps,
            Command(
                args,
                None,
                self.files["ancLog"] % {"method": newMethod},
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

    def addSelectionAnalysis(self, taskName=None, outputName=None, ancestor=None):

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
            Command(
                args,
                None,
                self.files["ancLog"] % {"method": newMethod},
            ),
        )
        self.prevMethod = newMethod
        self.selectionPool = []

        return task
