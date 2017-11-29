# -*- coding: utf-8 -*-

# LibsDyogen version 1.0 (6/11/2015)
# python v2.7 at least is needed
# Copyright © 2015 IBENS/Dyogen : Matthieu MUFFATTO, Joseph LUCAS and Hugues ROEST CROLLIUS
# mail : jlucas@ens.fr
# Licences GLP v3 and CeCILL v2

import os
import re
import sys
import itertools
import time
import string
import warnings
import collections

import enum

from functools import wraps
from collections import OrderedDict, Callable

import myFile

null = open(os.devnull, 'w')

debug = null

class Namespace: pass

# http://stackoverflow.com/questions/1969005/enumerations-in-python
# Usage
#>>> x = Enum('foo', 'bar', 'baz', 'bat')
#>>> x.baz
#2
#>>> x.bat
#3
class Enum(object):
    def __init__(self, *keys):
        self.____dict__.update(zip(keys, range(len(keys))))

def applyFunctions(fun, data):
    for (f, x) in itertools.izip(fun, data):
        yield f(x)

def funcFilter(fun):
    return lambda data: (f(x) for (f, x) in itertools.izip(fun, data))

def __delitem__(self, key):
    dict.__delitem__(self, self[key])
    dict.__delitem__(self, key)

def __len__(self):
    """Returns the number of connections"""
    return dict.__len__(self) // 2

# Print a correctly spaced table, all columns are aligned
# Input table format:
# table =[[l0c0, l0c1, l0c2],
#         [l1c0, l1c1, l1c2]]
# With 'lxcy' a string that corresponds to line 'x' and column 'y'
def printTable(table, output, spaceBetweenColumns=2):
    max_lens = []
    for i in range(len(table[0])):
        max_lens.append(max([len(str(r[i])) for r in table]))
    res = "\n".join(
        ["".join([string.ljust(str(e), l + spaceBetweenColumns)
                  for e, l in zip(r, max_lens)])
         for r in table])
    print >> output, res
    return res

# FIXME: to print well in stream, the user needs to ensure that between two calls to printProgressIn, nothing had been
# written is stream
class ProgressBar:
    def __init__(self, totalLength, step=1):
        self.totalLength = totalLength
        self.listOfPercentage = range(0, 101, step)[1:]

    def printProgressIn(self, stream, currentLength, prefixLabel=None):
        progress = int(float(currentLength*100)/self.totalLength)
        if progress in self.listOfPercentage:
            if prefixLabel is None:
                stream.write("Progress: %d%%   \r" % progress)
            else:
                # No return of the cursor.
                # If they are some prints between two execution of printProgressIn, this allow to see the line in the
                # console. It avoids the removal of the message of progress by intermediary lines.
                stream.write(prefixLabel + "progress: %d%%   \n" % progress)
            stream.flush()
            self.listOfPercentage.remove(progress)


# decorator that adds a switchable verbose mode to a function
# FIXME, if the decorated function is called with more arguments that in the
# function definition, there is no error raised.
def verbose(functionToExcecute):
    @wraps(functionToExcecute) # to avoid changing the name of the function
    def modifiedFunction(*args, **kargs):
        old_sys_stderr = sys.stderr
        if 'verbose' in kargs:
            if kargs['verbose'] == True:
                res = functionToExcecute(*args, **kargs)
                # **kargs still contains verbose
            else:
                sys.stderr = open(os.devnull, 'w')
                res = functionToExcecute(*args, **kargs)
                # **kargs still contains verbose
        else:
            warnings.warn("function %s has no option verbose although it uses a verbose decorator" % functionToExcecute.__name__, category=SyntaxWarning, stacklevel=2)
            res = functionToExcecute(*args, **kargs)
        sys.stderr = old_sys_stderr
        # sys.stderr = sys.__stderr__   if you wanna come back to a verbose mode
        return res
    return  modifiedFunction

# decorator for functions that requires a minimal python version >= 2.7 for instance
# version is a tuple, for instance if the function requires python version at least 2.7, version = (2,7)
def minimalPythonVersion(version):
    def decorator(functionToExcecute):
        def modifiedFunction(*args, **kargs):
            if sys.version_info < version:
                raise Exception("Function %s needs at least python %s.%s" % (functionToExcecute.__name__,version[0],version[1]))
            else:
                return functionToExcecute(*args, **kargs)
        return modifiedFunction
    return decorator

# decorator that computes the execution time
def tictac(functionToExcecute):
    @wraps(functionToExcecute) # to avoid changing the name of the function
    def modifiedFunction(*args,**kargs):
        tic = time.time()
        res = functionToExcecute(*args,**kargs)
        tac = time.time()
        deltaTicTac = tac - tic
        print >> sys.stderr, "Function \"%s\" was executed in %s seconds" % (functionToExcecute.__name__, deltaTicTac)
        return res
    return modifiedFunction

# decorator that warns the user that the function is deprecated
def deprecated(func):
    """This is a decorator which can be used to mark functions
    as deprecated. It will result in a warning being emmitted
    when the function is used."""
    import warnings
    def newFunc(*args, **kwargs):
        warnings.warn("Call to deprecated function %s." % func.__name__, category=DeprecationWarning, stacklevel=2)
        return func(*args, **kwargs)
    newFunc.__name__ = func.__name__
    newFunc.__doc__ = func.__doc__
    newFunc.__dict__.update(func.__dict__)
    return newFunc

# record results of a function for each parameter value #
class memoize:
    """Decorator that caches a value returned by the function each time it is called.
    If called later with the same arguments, the cached value is returned, and
    not re-evaluated.
    """
    def __init__(self, func):
        self.func = func
        self.nbcall = 0
        self.cache = {}

    def __repr__(self):
        return "[%s: %d values cached, %d calls, %.2fx speedup]" %\
            (self.func.__name__, len(self.cache),
             self.nbcall,
             self.nbcall/float(len(self.cache)) if len(self.cache) > 0 else 0)

    def __call__(self, *args, **kwargs):
        self.nbcall += 1
        try:
            return self.cache[args]
        except KeyError:
            self.cache[args] = self.func(*args)
            value = self.func(*args)
            return value
        except TypeError:
            # uncachable -- for instance, passing a list as an argument.
            # Better to not cache than to blow up entirely.
            print >> sys.stderr, "Warning: %s is not cacheable (from %s/%s)" % (args, self, self.func)
            return self.func(*args)

    def reinit_stats(self):
        self.nbcall = 0
        self.cache = {}

    def __doc__(self):
        """Return the function's docstring."""
        return self.func.__doc__

class memoizeMethod(object):
    def __init__(self, function):
        self.function = function
        self.nbcall = 0
        self.cache = {}

    def __get__(self, instance, cls=None):
        self.instance = instance
        return self

    def __call__(self, *args):
        self.nbcall += 1
        if args in self.cache:
            return self.cache[args]
        else:
            res = self.cache[args] = self.function(self.instance, *args)
            return res

    def reinit_stats(self):
        self.nbcall = 0
        self.cache = {}


# FIXME, best memoizeMethod
# cf http://code.activestate.com/recipes/577452-a-memoize-decorator-for-instance-methods/
#class memoizeMethod2(object):
#    """cache the return value of a method
#
#
#    This class is meant to be used as a decorator of methods. The return value
#    from a given method invocation will be cached on the instance whose method
#    was invoked. All arguments passed to a method decorated with memoize must
#    be hashable.
#
#    If a memoized method is invoked directly on its class the result will not
#    be cached. Instead the method will be invoked like a static method:
#    class Obj(object):
#        @memoize
#        def add_to(self, arg):
#            return self + arg
#    Obj.add_to(1) # not enough arguments
#    Obj.add_to(1, 2) # returns 3, result is not cached
#    """

# Check is an excecutable is accessible
# This may be usefull to check if a plugged external programm has been added to
# the PATH environment variable.
def which(program):
    import os
    def is_exe(fpath):
        return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

    fpath, fname = os.path.split(program)
    if fpath:
        if is_exe(program):
            return program
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            path = path.strip('"')
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return exe_file
    return None


# iterator of adjacent components of a list
class myIterator:
    # sliding tuple (x1, x2, ...) of width=width
    @staticmethod
    def slidingTuple(lst, width=2):
        if len(lst) < width:
            raise StopIteration
        else:
            # idxW0, idx of the left extremity of the sliding window
            for idxW0 in xrange(len(lst) - width + 1):
                yield tuple(lst[idxW0: idxW0 + width])

# liste of partitions of size k in range(n)
@memoize
def partitions(n, k):
    if n == 1:
        if k == 1:
            return [ [[0]] ]
    if (n >= k) and (k >= 1):
        all = []
        for x in partitions(n-1, k-1):
            all.append( x + [[n-1]] )
        for x in partitions(n-1, k):
            for i in xrange(k):
                all.append( [y if i != j else y + [n-1] for (j,y) in enumerate(x)] )
        return all
    else:
        return []

# management of a parallel execution on a range of values
def getRange(s):
    if myFile.hasAccess(s):
        f = myFile.openFile(s, "r")
        lst = []
        for l in f:
            lst.extend( [int(x) for x in l.replace('\n', '').split()] )
        f.close()
        return lst
    else:
        (start,_,end) = s.partition(':')
        return range(int(start), int(end)+1)


# hashable dict class, useful to use it as a key
class hashabledict(dict):
    def __hash__(self):
        return hash(tuple(sorted(self.items())))

# hashable list class
class hashablelist(list):
    def __hash__(self):
        return hash(tuple(self))

# This class allows to group a list of elements.
# From an initial list of elements, links are added between these elements.
# The class gathers elements that are linked.
class myCombinator:

    def __init__(self, ini = []):
        self.grp = list(ini)
        self.dic = {}
        for i in xrange(len(self.grp)):
            self.grp[i] = list(set(self.grp[i]))
            for x in self.grp[i]:
                self.dic[x] = i

    # define a link between all elements of obj
    # update sets already built
    def addLink(self, obj):

        if len(obj) == 0:
            return []

        obj = set(obj)
        grp = self.grp
        dic = self.dic

        # elements of obj already present in the combinator
        d = set( dic[x] for x in obj if x in dic )

        if len(d) == 0:
            # None, the obj is added just like it is
            i = len(grp)
            grp.append(list(set(obj)))
            for x in obj:
                dic[x] = i
            return grp
        else:
            i = d.pop()
            grpiextend = grp[i].extend
            for x in d:
                grpiextend(grp[x])
                for y in grp[x]:
                    dic[y] = i
                #FIXME not del grp[x] ?
                # see reduce and "empty sets"
                grp[x] = []
            dd = [x for x in obj if x not in dic]
            for x in dd:
                dic[x] = i
            grpiextend(dd)
            return grp[i]


    # return an iterator over the data
    # empty sets are thus removed
    def __iter__(self):
        for g in self.grp:
            if len(g) > 0:
                yield g

    # remove empty sets
    def reduce(self):
        self.__init__(self)

    # reset combinator
    def reset(self):
        self.__init__()

# add more options to a specific module
__moduleoptions = []
def addModuleOptions(namespace, options):
    for (name,typ,val) in options:
        __moduleoptions.append((namespace+":"+name, typ, val))


# ask a list of file in arguments
class FileList:
    def __init__(self, value):
        self.minNbFiles = value

    def __repr__(self):
        return '<FileList(%d)>' % self.minNbFiles


# Parse arguments on the command line
#  1. requested arguments (name, builder)
#  2. options in the form of -opt=val (name, builder, default_value)
# If an error occurs, user's command line is printed as well as a short description of the bug and a brief manual of the script (info).
def checkArgs(args, options, info, showArgs=True, loadOnlyDefaultOptions=False):

    options = options + __moduleoptions
    # print error informations if wrong arguments
    def error_usage(reason):
        print >> sys.stderr, "- ERROR -", reason
        print >> sys.stderr, " Usage :", sys.argv[0]
        for (i,t) in enumerate(args):
            print >> sys.stderr, "\t", "%d:" % (i+1), t[0], t[1]
        for t in options:
            if isinstance(t[1], enum.Enum):
                print >> sys.stderr, "\t", "  -%s %s (%s)" % (t[0], t[1]._keys, t[2])
            elif t[1] == bool:
                print >> sys.stderr, "\t", "+/-%s (%s)" % (t[0],t[2])
            else:
                print >> sys.stderr, "\t", "  -%s %s (%s)" % t
        if info != "":
            print >> sys.stderr, "\n", info
        sys.exit(1)

    def putValue(typ, authorisedVals, v):
        # instantiate the value depending on the type
        if typ == bool:
            # Type booleen
            res = {"false": False, "true":True}[v.lower()]
        elif typ == file:
            # Type 'fichier': test of presence
            v = os.path.expanduser(v)
            if not myFile.hasAccess(v):
                error_usage("File '%s' innaccessible" % v)
            else:
                res = v
        elif isinstance(typ, enum.Enum):
            try:
                res = getattr(typ, v)
            except AttributeError:
                error_usage("'%s' is not among %s" % (v,typ._keys))
        else:
            # otherwise the builder is used
            res = typ(v)
            if isinstance(authorisedVals, list) and (res not in authorisedVals):
                # non authorised parameter value
                error_usage("'%s' is not among %s" % (res,myFile.myTSV.printLine(authorisedVals, '/')))
        return res

    valOpt = collections.OrderedDict()
    opt = {}
    for (name,typ,val) in options:
        opt[name] = (typ,val)
        valOpt[name] = val[0] if isinstance(val, list) else getattr(typ, val) if isinstance(typ, enum.Enum) else val
    if loadOnlyDefaultOptions:
        return valOpt

    valArg = collections.OrderedDict()
    # arguments are scanned, counted and values are extracted
    for tt in sys.argv[1:]:

        t = tt.replace('^', ' ')

        # an optional argument
        if t[0] in '-+':

            # non bool parameter
            try:
                i = t.index('=')
                s = t[1:i]

                # the parameter name must be known
                if not s in valOpt:
                    error_usage("Option '%s' unknown" % s)

                valOpt[s] = putValue(opt[s][0], opt[s][1], t[i+1:])

            # if '=' is not found, it is a bool type
            except ValueError:
                s = t[1:]
                # unexpected parameter name
                if s not in valOpt:

                    # predefined values
                    if s.startswith("psyco"):
                        if t[0] == '+':
                            try:
                                import psyco
                                psyco.full()
                            except ImportError:
                                print >> sys.stderr, "Unable to load psyco !"
                    elif s == "bz2":
                        if t[0] == '+':
                            import bz2
                            sys.stdout = bz2.BZ2File("/dev/stdout", "w")
                    elif s == "gz":
                        if t[0] == '+':
                            import gzip
                            sys.stdout = gzip.GzipFile("/dev/stdout", "w")
                    elif (s == "lzma") or (s == "xz"):
                        if t[0] == '+':
                            import lzma
                            sys.stdout = lzma.LZMAFile("/dev/stdout", "w")
                    elif s == "debug":
                        if t[0] == '+':
                            global debug
                            debug = sys.stderr
                    else:
                        error_usage("Option '%s' unknown" % s)
                elif opt[s][0] != bool:
                    error_usage("Use -%s=value" % s)
                else:
                    # Here, False is assigned
                    valOpt[s] = (t[0] == '+')
        else:
            if len(valArg) < len(args):
                (s,typ) = args[len(valArg)]
                if isinstance(typ, FileList):
                    valArg[s] = list()
                    assert len(valArg) == len(args)
                    valArg[s].append(putValue(file, None, t))
                else:
                    valArg[s] = putValue(typ, None, t)
            elif isinstance(args[-1][1], FileList):
                valArg[args[-1][0]].append(putValue(file, None, t))
            else:
                error_usage("Too many arguments on '%s'" % t)

    if len(args[-1])>0 and isinstance(args[-1][1], FileList):
        if args[-1][0] not in valArg:
            valArg[args[-1][0]] = []
        if len(valArg[args[-1][0]]) < args[-1][1].minNbFiles:
            error_usage("Not enough files for '%s'" % args[-1][0])

    # there is less than the minimal number of arguments
    #FIXME, the second part of the condition should be avoided by upstream corrections
    if len(valArg) < len(args) and not (len(args) == 1 and args[0] == ()):
        print >> sys.stderr, "valArg=", valArg
        print >> sys.stderr, "args=", args
        error_usage("Not enough arguments")

    valArg.update(valOpt)
    if showArgs:
        # print >> sys.stderr, "Arguments:", valArg
        printArguments(valArg, stream=sys.stderr)
    return valArg

def printArguments(arguments, stream=open(os.devnull, 'w')):
    longestKey = 0
    longestValue = 0
    try:
        rows, columns = getTerminalSize()
    except:
        rows = 50
        columns = 80

    for (key, value) in arguments.iteritems():
        longestKey = max(len(str(key)), longestKey)
        longestValue = max(len(str(value)), longestValue)
    longestValue = min(longestValue, rows - longestKey - 7)
    lines = []
    lines.append('| ' + 'Key'.ljust(longestKey) + ' | ' + 'Values'.ljust(longestValue) + ' |')
    for (key, value) in arguments.iteritems():
        lines.append('| ' + str(key).ljust(longestKey) + ' | ' + str(value).ljust(longestValue) + ' |')
    longestLine = min(max([len(line) for line in lines]), rows)
    print >> stream, '-' * longestLine
    print >> stream, lines[0]
    print >> stream, '-' * longestLine
    for line in lines[1:]:
        print >> stream, line
    print >> stream, '-' * longestLine

# This class is useful for recording many information (either a list of items
# or a value) for each cell of a matrix of whole genome comparisons. See
# myDiags.py for instance.
class Dict2d(collections.defaultdict):
    # Idea to make something more general
    #def multi_dimensions(n, type):
    #    """ Creates an n-dimension dictionary where the n-th dimension is of type 'type' """
    #    if n <= 1:
    #        return type()
    #    return collections.defaultdict(lambda: multi_dimensions(n-1, type))

    def __init__(self, type):
        self.type = type
        collections.defaultdict.__init__(self, lambda: collections.defaultdict(type))

    def iteritems2d(self):
        assert self.type == list or self.type == set
        for (k1, k2) in self.keys2d():
            for item in self[k1][k2]:
                yield ((k1, k2), item)

    def items2d(self):
        return list(self.iteritems2d())

    def __add__(self, other):
        assert isinstance(other, Dict2d)
        res = Dict2d(self.type)
        for (k1, k2) in self.keys2d():
            res[k1][k2] = self.type(self[k1][k2])
        for (k1, k2) in other.keys2d():
            # if 'typ' was 'list', concatenate the lists
            # if 'typ' was int, add ints in 2D cells
            res[k1][k2] += self.type(other[k1][k2])
        return res

    def keys2d(self):
        return [(k1, k2) for k1 in collections.defaultdict.__iter__(self) for k2 in collections.defaultdict.__iter__(self[k1])]

    def values2d(self):
        return [self[k1][k2] for (k1, k2) in self.keys2d()]

    def intoList(self):
        res = []
        for ((k1, k2), item) in self.iteritems2d():
            res.append(((k1, k2), item))
        return res

class OrderedDict2dOfLists(Dict2d):
    def __init__(self):
        Dict2d.__init__(self, list)
        self.id2location = {}
        self.location2id = collections.defaultdict(lambda: collections.defaultdict(list))
        self.orderedIds = []
        self.maxId = 0

    def identifyItems(self):
        # FIXME
        # reinitialise data
        self.id2location = {}
        self.location2id = collections.defaultdict(lambda: collections.defaultdict(list))
        self.orderedIds = []
        self.maxId = 0
        # fill data
        for (k1, k2) in self.keys2d():
            for (idx, item) in enumerate(self[k1][k2]):
                self.maxId += 1
                self.id2location[self.maxId] = (k1, k2, idx)
                self.location2id[k1][k2].append(self.maxId)
                self.orderedIds.append(self.maxId)
        assert len(self.location2id[k1][k2]) == len(self[k1][k2])
        return self.id2location

    def getItemById(self, id):
        (k1, k2, idx) = self.id2location[id]
        return self[k1][k2][idx]

    def getItemsAndIdsByLocation(self, (k1, k2)):
        return [(self.getItemById(id), id) for id in self.location2id[k1][k2]]

    def getItemLocationById(self, id):
        # (k1, k2, idx) = self.id2location[id]
        return self.id2location[id]

    def iterByOrderedIds(self):
        for id in self.orderedIds:
            (k1, k2, idx) = self.id2location[id]
            yield (id, (k1, k2), self.getItemById(id))

    def removeIds(self, setOfRemovedIds):
        # debug assertions
        assert len(self.id2location.keys()) == len(self.orderedIds)
        assert set(self.id2location.keys()) == set(self.orderedIds)

        copyOrderedIds = list(self.orderedIds)  # need a deep copy
        for id in copyOrderedIds:
            if id in setOfRemovedIds:
                (k1, k2, idx) = self.getItemLocationById(id)
                # remove the item of self[k1][k2] at the index 'idx'
                del self[k1][k2][idx]
                del self.id2location[id]
                del self.location2id[k1][k2][idx]
                self.orderedIds.remove(id)
                for (higherSbIdx, higherSbId) in enumerate(self.location2id[k1][k2][idx:]):
                    higherSbIdx = higherSbIdx + idx
                    self.id2location[higherSbId] = (k1, k2, higherSbIdx)
                assert len(self.location2id[k1][k2]) == len(self[k1][k2])

                # debug assertions
                assert len(self.id2location.keys()) == len(self.orderedIds)
                assert set(self.id2location.keys()) == set(self.orderedIds)

                if len(self[k1][k2]) == 0:
                    del self[k1][k2]
                    if len(self[k1]) == 0:
                        del self[k1]
        # DEBUG assertion
        #for k1 in self:
        #    for k2 in self[k1]:
        #        assert len(self[k1][k2])>0, "k1=%s, k2=%s" % (k1, k2)

    def keepIds(self, setOfKeptIds):
        allIds = set([id for (id, _, _) in self.iterByOrderedIds()])
        setOfRemovedIds = allIds - setOfKeptIds
        self.removeIds(setOfRemovedIds)

    def addToLocation(self, (k1, k2), item):
        self.maxId += 1
        id = self.maxId
        self.addToLocationWithId((k1, k2), item, id)
        return id

    def addToLocationWithId(self, (k1, k2), item, id):
        assert id not in self.orderedIds
        self.maxId = id if id > self.maxId else self.maxId
        self[k1][k2].append(item)
        self.id2location[id] = (k1, k2, len(self[k1][k2]) - 1)
        self.location2id[k1][k2].append(id)
        self.orderedIds.append(id)

    def __add__(self, other):
        assert isinstance(other, OrderedDict2dOfLists)
        res = OrderedDict2dOfLists()
        for (id, (k1, k2), item) in self.iterByOrderedIds():
            res.addToLocationWithId((k1, k2), item, id)
        minSelfId = min([sys.maxint] + self.orderedIds)
        maxSelfId = max([0] + self.orderedIds)
        for (id, (k1, k2), item) in other.iterByOrderedIds():
            assert id < minSelfId or maxSelfId < id
            res.addToLocationWithId((k1, k2), item, id)
        return res

    def intoList(self):
        res = []
        for (id, (k1, k2), item) in self.iterByOrderedIds():
            res.append(((k1, k2), item, id))
        return res

# http://stackoverflow.com/questions/9001509/how-can-i-sort-a-dictionary-by-key
# to order the keys of an OrderedDict: d = collections.OrderedDict(sorted(d.items()))

# This class is a fusion of collections.defaultdict and collections.OrderedDict
class DefaultOrderedDict(OrderedDict):
    # Source: http://stackoverflow.com/a/6190500/562769
    def __init__(self, default_factory=None, *a, **kw):
        if (default_factory is not None and
           not isinstance(default_factory, Callable)):
            raise TypeError('first argument must be callable')
        OrderedDict.__init__(self, *a, **kw)
        self.default_factory = default_factory

    def __getitem__(self, key):
        try:
            return OrderedDict.__getitem__(self, key)
        except KeyError:
            return self.__missing__(key)

    def __missing__(self, key):
        if self.default_factory is None:
            raise KeyError(key)
        self[key] = value = self.default_factory()
        return value

    def __reduce__(self):
        if self.default_factory is None:
            args = tuple()
        else:
            args = self.default_factory,
        return type(self), args, None, None, self.items()

    def copy(self):
        return self.__copy__()

    def __copy__(self):
        return type(self)(self.default_factory, self)

    def __deepcopy__(self, memo):
        import copy
        return type(self)(self.default_factory,
                          copy.deepcopy(self.items()))

    def __repr__(self):
        return 'OrderedDefaultDict(%s, %s)' % (self.default_factory,
                                               OrderedDict.__repr__(self))


#class DefaultOrderedDictOldVersion(collections.OrderedDict):
#    # Source: http://stackoverflow.com/a/6190500/562769
#    def __init__(self, type, *a, **kw):
#        self.type = type
#        if hasattr(type, '__call__'): # isinstance(type, collections.Callable): Not working
#            raise TypeError('first argument must be callable')
#        collections.OrderedDict.__init__(self, type)
#
#    def __getitem__(self, key):
#        #try:
#        return collections.OrderedDict.__getitem__(self, key)
#        #except KeyError:
#        #    return self.__missing__(key)
#
#    def __missing__(self, key):
#        self[key] = self.type()
#        return self[key]
#
#    def __reduce__(self): # optional, for pickle support
#        if self.default_factory is None:
#            args = tuple()
#        else:
#            args = self.default_factory,
#        return type(self), args, None, None, iter(self.items())
#
#    def copy(self):
#        return self.__copy__()
#
#    def __copy__(self):
#       return type(self)(self.type, self)
#
#    def __deepcopy__(self, memo):
#        import copy
#        return type(self)(self.default_factory,
#                          copy.deepcopy(self.items()))
#
#    def __repr__(self):
#        return 'OrderedDefaultDict(%s, %s)' % (self.default_factory,
#                                               collections.OrderedDict.__repr__(self))


# http://stackoverflow.com/questions/566746/how-to-get-console-window-width-in-python
# """ getTerminalSize()
#  - get width and height of console
#  - works on linux,os x,windows,cygwin(windows)
# """
def getTerminalSize():
   import platform
   current_os = platform.system()
   tuple_xy=None
   if current_os == 'Windows':
       tuple_xy = _getTerminalSize_windows()
       if tuple_xy is None:
          tuple_xy = _getTerminalSize_tput()
          # needed for window's python in cygwin's xterm!
   if current_os == 'Linux' or current_os == 'Darwin' or  current_os.startswith('CYGWIN'):
       tuple_xy = _getTerminalSize_linux()
   if tuple_xy is None:
       print "default"
       tuple_xy = (80, 25)      # default value
   return tuple_xy

def _getTerminalSize_windows():
    res=None
    try:
        from ctypes import windll, create_string_buffer

        # stdin handle is -10
        # stdout handle is -11
        # stderr handle is -12

        h = windll.kernel32.GetStdHandle(-12)
        csbi = create_string_buffer(22)
        res = windll.kernel32.GetConsoleScreenBufferInfo(h, csbi)
    except:
        return None
    if res:
        import struct
        (bufx, bufy, curx, cury, wattr,
         left, top, right, bottom, maxx, maxy) = struct.unpack("hhhhHhhhhhh", csbi.raw)
        sizex = right - left + 1
        sizey = bottom - top + 1
        return sizex, sizey
    else:
        return None

# TODO
def _getTerminalSize_tput():
    # get terminal width
    # src: http://stackoverflow.com/questions/263890/how-do-i-find-the-width-height-of-a-terminal-window
    try:
       import subprocess
       proc=subprocess.Popen(["tput", "cols"],stdin=subprocess.PIPE,stdout=subprocess.PIPE)
       output=proc.communicate(input=None)
       cols=int(output[0])
       proc=subprocess.Popen(["tput", "lines"],stdin=subprocess.PIPE,stdout=subprocess.PIPE)
       output=proc.communicate(input=None)
       rows=int(output[0])
       return (cols,rows)
    except:
       return None


def _getTerminalSize_linux():
    def ioctl_GWINSZ(fd):
        try:
            import fcntl, termios, struct, os
            cr = struct.unpack('hh', fcntl.ioctl(fd, termios.TIOCGWINSZ,'1234'))
        except:
            return None
        return cr
    cr = ioctl_GWINSZ(0) or ioctl_GWINSZ(1) or ioctl_GWINSZ(2)
    if not cr:
        try:
            fd = os.open(os.ctermid(), os.O_RDONLY)
            cr = ioctl_GWINSZ(fd)
            os.close(fd)
        except:
            pass
    return int(cr[1]), int(cr[0])

def atoi(text):
    return int(text) if text.isdigit() else text

## http://stackoverflow.com/questions/23465462/check-if-sorted-using-recursion-in-python
# def _isSorted(L, increasingOrder=True):
#     return len(L) < 2 or (L[0] <= L[1] and isSorted(L[1:]))

def isSorted(l, increasingOrder=False, stricly=False, key=lambda x: x):
    assert isinstance(l, list)
    assert isinstance(increasingOrder, bool)
    if increasingOrder:
        if stricly:
            res= all(key(l[i]) < key(l[i+1]) for i in xrange(len(l)-1))
        else:
            res = all(key(l[i]) <= key(l[i+1]) for i in xrange(len(l)-1))
    else:
        if stricly:
            res = all(key(l[i]) > key(l[i+1]) for i in xrange(len(l)-1))
        else:
            res = all(key(l[i]) >= key(l[i+1]) for i in xrange(len(l)-1))
    return res

def keyNaturalSort(chrName):
    '''
    alist.sort(key=natural_keys) sorts in human wished order
    http://nedbatchelder.com/blog/200712/human_sorting.html
    (See Toothy's implementation in the comments)
    '''
    res = [atoi(c) for c in re.split('([0-9]+)', str(chrName))]
    return res

# http://stackoverflow.com/questions/14380371/export-a-latex-table-from-pandas-dataframe
# for exemple tableIntoLatex(np.random.random((5, 5)))
# in latex insert \usepackage{booktabs} in the imports
# def tableIntoLatex(numpyArray):
#     import pandas as pd
#     df = pd.DataFrame(numpyArray)
#     return df.to_latex()

# from tabulate import tabulate
# table = [["spam",42],["eggs",451],["bacon",0]]
# headers = ["item", "qty"]
# print tabulate(table, headers, tablefmt="latex")

# Beautiful example of how to draw an empirical distribution: source :
# http://stackoverflow.com/questions/9378420/how-to-plot-cdf-in-matplotlib-in-python
# import numpy as np
# from pylab import *
#
# # Create some test data
# dx = .01
# X  = np.arange(-2,2,dx)
# Y  = exp(-X**2)
#
# # Normalize the data to a proper PDF
# LOVE THAT!
# Y /= (dx*Y).sum()
#
# (or even better: Y /= np.trapz(Y,X), if the bins are not all equal)
#
# # Compute the CDF
# CY = np.cumsum(Y*dx)
#
# # Plot both
# plot(X,Y)
# plot(X,CY,'r--')
#
# show()


