#! /usr/bin/env python
# -*- coding: utf-8 -*-
# Agora v1.00
# python 2.7
# Copyright © 2020 IBENS/Dyogen and EMBL-European Bioinformatics Institute : Matthieu MUFFATO, Alexandra LOUIS, Thi Thuy Nga NGUYEN, Hugues ROEST CROLLIUS
# mail : agora@bio.ens.psl.eu
# This is free software; you may copy, modify and/or distribute this work under the terms of the GNU General Public License, version 3 or later and the CeCiLL v2 license in France

import multiprocessing
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
        print "New task", len(self.list), name
        print dep
        print data
        print
        self.dic[name] = len(self.list)
        self.list.append((set(self.dic[x] for x in dep), data, multithreaded))

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
