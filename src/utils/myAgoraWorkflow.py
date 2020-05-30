#! /usr/bin/env python
# -*- coding: utf-8 -*-
# Agora v1.00
# python 2.7
# Copyright © 2020 IBENS/Dyogen and EMBL-European Bioinformatics Institute : Matthieu MUFFATO, Alexandra LOUIS, Thi Thuy Nga NGUYEN, Hugues ROEST CROLLIUS
# mail : agora@bio.ens.psl.eu
# This is free software; you may copy, modify and/or distribute this work under the terms of the GNU General Public License, version 3 or later and the CeCiLL v2 license in France

import multiprocessing
import os
import subprocess
import sys
import time

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
        self.queue = manager.Queue()

    def addTask(self, name, dep, data, multithreaded=False):
        taskId = len(self.list)
        print "New task", taskId, name
        print dep
        print data
        print
        self.dic[name] = taskId
        self.list.append((set(self.dic[x] for x in dep), data, multithreaded))
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
        print "Waiting ...",
        sys.stdout.flush()
        (i, r) = self.queue.get()
        print "task", i, "is now finished (status", r, ")"
        if r == 0:
            self.removeDep(i)
        self.proc.pop(i).join()
        if r == 0:
            self.completed += 1
        else:
            self.failed += 1
        self.nrun -= self.nthreads.pop(i)

    # Launch program function
    def goLaunch(self, i, args, out, log):
        stdout = myFile.openFile(out, "w")
        stderr = myFile.openFile(log, "w")
        p = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=stderr)
        for l in p.stdout:
            print >> stdout, l,
        r = p.wait()
        for l in p.stdout:
            print >> stdout, l,
        stdout.close()
        stderr.close()
        time.sleep(5)
        self.queue.put((i, r))

    # Launching tasks in multiple threads
    def runAll(self, nbThreads):
        # Queue
        while (self.completed + self.failed) < len(self.list):

            print "running:", self.nrun
            print "todo:", len(self.list)
            print "done:", self.completed
            print "failed:", self.failed

            if self.nrun == nbThreads:
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
                    print ("Launching" if launch else "Skipping"), next, args, ">", out, "2>", log
                    if launch:
                        if multithreaded:
                            self.nthreads[next] = nbThreads - self.nrun
                            args.append("-nbThreads=%d" % self.nthreads[next])
                            print "Using", self.nthreads[next], "threads"
                        else:
                            self.nthreads[next] = 1
                        self.proc[next] = multiprocessing.Process(target=self.goLaunch, args=(next, args, out, log))
                        self.proc[next].start()
                        self.nrun += self.nthreads[next]
                    else:
                        self.removeDep(next)
                        self.completed += 1

        assert self.nrun == 0
        if not self.failed:
            print "Workflow complete"
        return self.failed

class AgoraWorkflow:

    def __init__(self, defaultRoot, scriptDir, files):
        self.defaultRoot = defaultRoot
        self.tasklist = TaskList()
        self.scriptDir = scriptDir
        self.files = files
        self.allAncGenesTaskName = "all"
        self.allAncGenesDirName = "all"
        self.interm = {}
        self.refMethod = {}

    def addDummy(self, taskFullName, dependencies=[]):
        return self.tasklist.addTask(taskFullName, dependencies, (None, None, None, False))

    def addAncGenesGenerationAnalysis(self, launch):
        taskFullName = ("ancgenes", self.allAncGenesTaskName)
        if launch:
            return self.tasklist.addTask(
                taskFullName,
                [],
                (
                    [
                        os.path.join(self.scriptDir, "ALL.extractGeneFamilies.py"),
                        self.files["speciestree"],
                        self.files["genetrees"],
                        "-OUT.ancGenesFiles=" + self.files["ancgenesdata"] % {"filt": self.allAncGenesDirName, "name": "%s"},
                    ],
                    self.files["genetreeswithancnames"],
                    self.files["ancgeneslog"] % {"filt": "ancGenes"},
                    launch,
                )
            )
        else:
            return self.addDummy(taskFullName)

    def addAncGenesFilterAnalysis(self, taskName, methodName, dirnameTemplate, ancestor, params, launch):
        return self.tasklist.addTask(
            ("ancgenes", taskName),
            [("ancgenes", self.allAncGenesTaskName)],
            (
                [
                    os.path.join(self.scriptDir, "ALL.filterGeneFamilies-%s.py" % methodName),
                    self.files["speciestree"],
                    ancestor or self.defaultRoot,
                    self.files["ancgenesdata"] % {"filt": self.allAncGenesDirName, "name": "%s"},
                    self.files["ancgenesdata"] % {"filt": dirnameTemplate, "name": "%s"}
                ] + params,
                os.devnull,
                self.files["ancgeneslog"] % {"filt": taskName},
                launch,
            )
        )

    def addPairwiseAnalysis(self, taskName, methodName, ancestor, params, launch):
        return self.tasklist.addTask(
            ("pairwise", taskName),
            [("ancgenes",taskName)],
            (
                [
                    os.path.join(self.scriptDir, "buildSynteny.pairwise-%s.py" % methodName),
                    self.files["speciestree"],
                    ancestor or self.defaultRoot,
                    "-ancGenesFiles=" + self.files["ancgenesdata"] % {"filt": taskName, "name": "%s"},
                    "-genesFiles=" + self.files["genes"] % {"name": "%s"},
                    "-OUT.pairwise=" + self.files["pairwiseoutput"] % {"filt": taskName, "name": "%s"}
                ] + params,
                os.devnull,
                self.files["pairwiselog"] % {"filt": taskName},
                launch,
            )
        )

    def addIntegrationAnalysis(self, taskName, methodName, inputName, outputName, pairwiseName, ancestor, params, launch):

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
            if methodName == "refine":
                self.refMethod[newMethod] = (newMethod, ancestor)
            else:
                self.refMethod[newMethod] = self.refMethod[self.prevMethod]

        if outputName:
            self.interm[outputName] = newMethod

        # task parameters
        args = [
                os.path.join(self.scriptDir, "buildSynteny.integr-%s.py" % methodName),
                self.files["speciestree"],
                ancestor,
        ] + params

        if methodName == "publish":
            # "publish" is not an integration method
            args[0] = os.path.join(self.scriptDir, "convert.ancGenomes.diags-genes.py")
            args.append("-OUT.ancGenomes=" + self.files["ancgenomesoutput"] % {"method": newMethod, "name": "%s"})
            logfile = self.files["ancgenomeslog"]
        else:
            args.append("-OUT.ancDiags=" + self.files["integrblocks"] % {"method": newMethod, "name": "%s"})
            logfile = self.files["integrlog"]

        dep = []
        if pairwiseName is not None:
            dep.append(("pairwise", pairwiseName))
            args.append(self.files["pairwiseoutput"] % {"filt": pairwiseName, "name": "%s"})

        if methodName in ["denovo", "groups", "publish"]:
            args.append("-ancGenesFiles=" + self.files["ancgenesdata"] % {"filt": "all", "name": "%s"})

        # No input data to consider for the denovo method
        if methodName != "denovo":
            dep.append(("integr", self.prevMethod))
            args.append("-IN.ancDiags=" + self.files["integrblocks"] % {"method": self.prevMethod, "name": "%s"})

        if methodName == "halfinsert":
            # The script needs singleton reference for "halfinsert"
            dep.append(("integr", self.refMethod[newMethod][0]))
            args.append("-REF.ancDiags=" + self.files["integrblocks"] % {"method": self.refMethod[newMethod][0], "name": "%s"})

        if methodName == "groups":
            args.append("-genesFiles=" + self.files["genes"] % {"name": "%s"})

        if methodName not in ["copy", "publish"]:
            args.append("-LOG.ancGraph=" + self.files["integroutput"] % {"method": newMethod, "name": "%s"})

        # Currently, all those methods are multithreaded
        multithreaded = True

        # The publish method doesn't generate integrDiags and can't be used as an input method
        if methodName != "publish":
            self.prevMethod = newMethod

        return self.tasklist.addTask(
            ("integr", newMethod),
            dep,
            (
                args,
                os.devnull,
                logfile % {"method": newMethod},
                launch,
            ),
            multithreaded,
        )
