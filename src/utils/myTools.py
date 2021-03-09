# -*- coding: utf-8 -*-

# LibsDyogen version 1.0 (6/11/2015) -- modified for AGORA v3.0
# python v2.7 at least is needed
# Copyright © 2006-2021 IBENS/Dyogen, 2020-2021 EMBL-European Bioinformatics Institute, 2021 Genome Research Ltd : Matthieu MUFFATO, Joseph LUCAS and Hugues ROEST CROLLIUS
# mail : agora@bio.ens.psl.eu
# Licences GLP v3 and CeCILL v2

import os
import re
import sys
import itertools
import time
import string
import warnings
import collections

from functools import wraps

from . import myFile

class file: pass

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
        self.__dict__.update(list(zip(keys, list(range(len(keys))))))

def applyFunctions(fun, data):
    for (f, x) in zip(fun, data):
        yield f(x)

def funcFilter(fun):
    return lambda data: (f(x) for (f, x) in zip(fun, data))

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
        ["".join([str(e).ljust(l + spaceBetweenColumns)
                  for e, l in zip(r, max_lens)])
         for r in table])
    print(res, file=output)
    return res

# FIXME: to print well in stream, the user needs to ensure that between two calls to printProgressIn, nothing had been
# written is stream
class ProgressBar:
    def __init__(self, totalLength, step=1):
        self.totalLength = totalLength
        self.listOfPercentage = list(range(0, 101, step))[1:]

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

# decorator that computes the execution time
def tictac(functionToExcecute):
    @wraps(functionToExcecute) # to avoid changing the name of the function
    def modifiedFunction(*args,**kargs):
        tic = time.time()
        res = functionToExcecute(*args,**kargs)
        tac = time.time()
        deltaTicTac = tac - tic
        print("Function \"%s\" was executed in %s seconds" % (functionToExcecute.__name__, deltaTicTac), file=sys.stderr)
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
            print("Warning: %s is not cacheable (from %s/%s)" % (args, self, self.func), file=sys.stderr)
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
        if len(lst) >= width:
            # idxW0, idx of the left extremity of the sliding window
            for idxW0 in range(len(lst) - width + 1):
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
            for i in range(k):
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
        return list(range(int(start), int(end)+1))


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
        for i in range(len(self.grp)):
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


# List of command-line arguments of a certain type
class ParamList:
    def __init__(self, typ, value):
        self.typ = typ
        self.minNb = value

    def __repr__(self):
        return '<ParamList(%s, %d)>' % (self.typ, self.minNb)

# Specialised version that requires files
def FileList(value):
    return ParamList(file, value)


# Parse arguments on the command line
#  1. requested arguments (name, builder)
#  2. options in the form of -opt=val (name, builder, default_value)
# If an error occurs, user's command line is printed as well as a short description of the bug and a brief manual of the script (info).
def checkArgs(args, options, info, showArgs=True, loadOnlyDefaultOptions=False):

    options = options + __moduleoptions
    # print error informations if wrong arguments
    def error_usage(reason):
        print("- ERROR -", reason, file=sys.stderr)
        print(" Usage :", sys.argv[0], file=sys.stderr)
        for (i,t) in enumerate(args):
            print("\t", "%d:" % (i+1), t[0], list(t[1].__dict__) if isinstance(t[1], Enum) else t[1], file=sys.stderr)
        for t in options:
            if isinstance(t[1], Enum):
                print("\t", "  -%s %s (%s)" % (t[0], list(t[1].__dict__), t[2]), file=sys.stderr)
            elif t[1] == bool:
                print("\t", "+/-%s (%s)" % (t[0],t[2]), file=sys.stderr)
            else:
                print("\t", "  -%s %s (%s)" % t, file=sys.stderr)
        if info != "":
            print("\n", info, file=sys.stderr)
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
        elif isinstance(typ, Enum):
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
        valOpt[name] = val[0] if isinstance(val, list) else getattr(typ, val) if isinstance(typ, Enum) else val
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
                    if s == "bz2":
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
                if isinstance(typ, ParamList):
                    valArg[s] = list()
                    assert len(valArg) == len(args)
                    valArg[s].append(putValue(typ.typ, None, t))
                else:
                    valArg[s] = putValue(typ, None, t)
            elif isinstance(args[-1][1], ParamList):
                valArg[args[-1][0]].append(putValue(typ.typ, None, t))
            else:
                error_usage("Too many arguments on '%s'" % t)

    if len(args[-1])>0 and isinstance(args[-1][1], ParamList):
        if args[-1][0] not in valArg:
            valArg[args[-1][0]] = []
        if len(valArg[args[-1][0]]) < args[-1][1].minNb:
            error_usage("Not enough values for '%s'" % args[-1][0])

    # there is less than the minimal number of arguments
    #FIXME, the second part of the condition should be avoided by upstream corrections
    if len(valArg) < len(args) and not (len(args) == 1 and args[0] == ()):
        # print >> sys.stderr, "valArg=", valArg
        # print >> sys.stderr, "args=", args
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

    for (key, value) in arguments.items():
        longestKey = max(len(str(key)), longestKey)
        longestValue = max(len(str(value)), longestValue)
    longestValue = min(longestValue, rows - longestKey - 7)
    lines = []
    lines.append('| ' + 'Key'.ljust(longestKey) + ' | ' + 'Values'.ljust(longestValue) + ' |')
    for (key, value) in arguments.items():
        lines.append('| ' + str(key).ljust(longestKey) + ' | ' + str(value).ljust(longestValue) + ' |')
    longestLine = min(max([len(line) for line in lines]), rows)
    print('-' * longestLine, file=stream)
    print(lines[0], file=stream)
    print('-' * longestLine, file=stream)
    for line in lines[1:]:
        print(line, file=stream)
    print('-' * longestLine, file=stream)

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
       print("default")
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
            res= all(key(l[i]) < key(l[i+1]) for i in range(len(l)-1))
        else:
            res = all(key(l[i]) <= key(l[i+1]) for i in range(len(l)-1))
    else:
        if stricly:
            res = all(key(l[i]) > key(l[i+1]) for i in range(len(l)-1))
        else:
            res = all(key(l[i]) >= key(l[i+1]) for i in range(len(l)-1))
    return res

def keyNaturalSort(chrName):
    '''
    alist.sort(key=natural_keys) sorts in human wished order
    http://nedbatchelder.com/blog/200712/human_sorting.html
    (See Toothy's implementation in the comments)
    '''
    res = [atoi(c) for c in re.split('([0-9]+)', str(chrName))]
    return res

