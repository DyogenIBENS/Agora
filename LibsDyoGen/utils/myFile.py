# -*- coding: utf-8 -*-

# LibsDyogen version 1.0 (6/11/2015)
# python v2.7 at least is needed
# Copyright © 2015 IBENS/Dyogen : Matthieu MUFFATTO, Joseph LUCAS and Hugues ROEST CROLLIUS
# mail : jlucas@ens.fr
# Licences GLP v3 and CeCILL v2

# file management functions

import itertools
import collections
import os
import sys

null = open(os.devnull, 'w')

# tabular file management
class myTSV:

    csvProxy = collections.namedtuple("csvProxy", ['file','csvobject'])

    # read with the csv module
    @staticmethod
    def reader(fileName, **keywords):
        import csv
        f = openFile(fileName, 'r')
        return myTSV.csvProxy(f,csv.reader(f, delimiter="\t", quoting=csv.QUOTE_NONE, lineterminator="\n", **keywords))

    # write with the csv module
    @staticmethod
    def writer(fileName):
        import csv
        f = openFile(fileName, 'w')
        return myTSV.csvProxy(f,csv.writer(f, delimiter="\t", quoting=csv.QUOTE_NONE, lineterminator="\n"))

    # return the prepared line for printing
    @staticmethod
    def printLine(line, delim = "\t", func = str):
        return delim.join(func(x) for x in line)


    # read a tabular file, convert columns separated by delim depending on the type_list
    @staticmethod
    def readTabular(filename, type_list, delim = '\t'):
        f = openFile(filename, 'r')
        # list of each column type
        new_type_list = []
        for x in type_list:
            if type(x) == type:
                new_type_list.append(x)
            else:
                new_type_list.extend([x[0]] * x[1])
        # run through the file (parcours du fichier)
        for (i,line) in enumerate(f):
            current_line = line.replace('\n','').split(delim)
            assert len(current_line) == len(new_type_list), "Error number of columns. Line:%d" % (i+1)
            yield tuple(t(x) for (x,t) in itertools.izip(current_line,new_type_list))
        f.close()

    # load MySQL dumps (join truncated lines)
    @staticmethod
    def MySQLFileLoader(f):
        tmp = ""
        for ligne in f:
            ligne = ligne.replace('\n', '')
            if ligne[-1] == '\\':
                # sign that shows that the line is not finished
                tmp = ligne[:-1]
            else:
                yield tmp + ligne
                tmp = ""
        assert (tmp == "")

    # write a mySQL dump file (\N for NULL)
    @staticmethod
    def MySQLFileWriter(data):
        return myTSV.printLine(data).replace("None", "\N")

# Read the first line of a file. Useful when you want to know the format without having to open-close it before opening it once more
class firstLineBuffer:
    def __init__(self, f):
        self.f = f
        try:
            self.firstLine = self.next()
        except StopIteration:
            self.firstLine = ""

    def __iter__(self):
        yield self.firstLine
        while True:
            yield self.next()

    def next(self):
        while True:
            l = self.f.next().replace('\n', '').replace('\r', '')
            # Suppression of the lines with comments
            if (not l.startswith("#")) and (len(l) > 0):
                return l

    def close(self):
        return self.f.close()


# existing file
def hasAccess(s):
    return os.access(os.path.expanduser(s), os.R_OK)


# open a file and decompress it if possible
# return the object 'file' and the full name of the file
def openFile(nom, mode):

    # file already open
    if type(nom) != str:
        return nom

    # file on the web
    elif nom.startswith("http://") or nom.startswith("ftp://"):
        comm = "wget %s -O -"
        # Compression bzip2
        if nom.endswith(".bz2"):
            comm += " | bunzip2"
        # Compression gzip
        elif nom.endswith(".gz"):
            comm += " | gunzip"
        # Compression lzma
        elif nom.endswith(".lzma"):
            comm += " | unlzma"
        (stdin,f,stderr) = os.popen3( comm % nom )
        stdin.close()
        stderr.close()

    # standard entry
    elif nom == "-":
        return sys.stdin

    # file on the disk
    else:
        nom = os.path.expanduser(nom)
        if ("w" in mode) or ("a" in mode):
            # create the folder for the output into files
            try:
                os.makedirs(os.path.dirname(nom))
            except OSError:
                pass
        i = nom.find(".zip/")
        if (mode == "r") and (i >= 0):
            import zipfile
            import cStringIO
            f = zipfile.ZipFile(nom[:i+4], "r")
            f = cStringIO.StringIO(f.read(nom[i+5:]))
        # Compression bzip2
        elif nom.endswith(".bz2"):
            import bz2
            f = bz2.BZ2File(nom, mode)
        # Compression gzip
        elif nom.endswith(".gz"):
            import gzip
            f = gzip.GzipFile(nom, mode)
        # Compression lzma
        elif nom.endswith(".lzma") or nom.endswith(".xz"):
            import lzma
            f = lzma.LZMAFile(nom, mode)
        else:
            f = open(nom, mode)
    return f
