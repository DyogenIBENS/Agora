# -*- coding: utf-8 -*-

# LibsDyogen version 1.0 (6/11/2015) -- modified for AGORA v3.0
# python v2.7 at least is needed
# Copyright © 2006-2021 IBENS/Dyogen, 2020-2021 EMBL-European Bioinformatics Institute, 2021 Genome Research Ltd : Matthieu MUFFATO, Joseph LUCAS and Hugues ROEST CROLLIUS
# mail : agora@bio.ens.psl.eu
# Licences GLP v3 and CeCILL v2

import random
import sys
import math
import operator
from . import myTools
import bisect
from functools import reduce

#############################
# Fonctions de statistiques #
#############################
class myStats:

    # Renvoie la moyenne d'une liste
    #################################
    @staticmethod
    def mean(lst):
        if len(lst) == 0:
            return None
        return float(sum(lst))/float(len(lst))

    # Renvoie l'ecart type
    #######################
    @staticmethod
    def stddev(lst, m=None, corrected=False):
        if len(lst) == 0:
            return None
        if m == None:
            m = myStats.mean(lst)
        s = sum((x-m) ** 2 for x in lst)
        return math.sqrt(float(s) / float(len(lst)-int(corrected)))

    # Renvoie la valeur comprise a x% des donnees (x=50 -> mediane)
    ################################################################
    @staticmethod
    def getValue(lst, x):
        lst.sort()
        if len(lst) == 0:
            return None
        return lst[int((x*(len(lst)-1))/100.)]

    # Renvoie la valeur telle que x% soit au dessus (x=50 -> N50)
    # Before using this function, do not forget to sort the input list
    #############################################################
    @staticmethod
    def getValueNX(lst, x):
        lst.sort()
        if len(lst) == 0:
            return None
        tmp = (sum(lst) * float(x)) / 100.
        for x in lst.__reversed__():
            tmp -= x
            if tmp <= 0:
                return x

    # Returns the weighted average of a list of lengths
    #######################################################
    @staticmethod
    def getWeightedAverage(listOfLengths):
        if len(listOfLengths) == 0:
            return None
        else:
            tmp = []
            sumLengths = sum(listOfLengths)
            for length in listOfLengths:
                weight = float(length) / sumLengths
                tmp.append(weight * length)
            res = float(sum(tmp))
            return res

    # Returns (min quart1 median quart3 N75 N50 N25 max mean stddev len)
    ####################################################################
    @staticmethod
    def valSummary(lst):
        l = list(lst)
        m = myStats.mean(l)
        return (myStats.getValue(l, 0), myStats.getValue(l, 25),myStats.getValue(l, 50),myStats.getValue(l, 75), \
                myStats.getValueNX(l, 75),myStats.getValueNX(l, 50),myStats.getValueNX(l, 25), myStats.getValue(l, 100), m, myStats.stddev(l, m), len(l))
    
    # Returns (min quart1 median quart3  N75 N50 N25 weightedAverage max mean stddev len)
    #####################################################################################
    @staticmethod
    def valSummary2(lst):
        l = list(lst)
        m = myStats.mean(l)
        return (myStats.getValue(l, 0), myStats.getValue(l, 25),myStats.getValue(l, 50),myStats.getValue(l, 75), \
                myStats.getValueNX(l, 75),myStats.getValueNX(l, 50),myStats.getValueNX(l, 25), myStats.getWeightedAverage(l), myStats.getValue(l, 100), m, myStats.stddev(l, m), len(l))

    # Returns (min N50 weightedAverage max)
    #########################################################
    @staticmethod
    def syntheticTxtSummary(lst):
        l = list(lst)
        l.sort()
        if len(l) > 0:
            res = "min=%s\tN50=%s\tWA=%.2f\tmax=%s" % (l[0], myStats.getValueNX(l, 50), myStats.getWeightedAverage(l), l[-1])
        else:
            res = "empty list, no synthetic info on the distribution is available"
        return res


    @staticmethod
    def txtSummary(lst, withN50=True):
        pattern = ["%s [%s/%s/%s] "]
        if withN50:
            pattern.append("[%s/%s/%s] ")
        pattern.append("%s [")
        pattern.append("%.2f/%.2f" if len(lst) > 0 else "%s/%s")
        pattern.append("-%d]")
        if len(lst) > 0:
            data = myStats.valSummary(lst)
            return ''.join(pattern) % (data if withN50 else data[:4] + data[-4:])
        else:
            return ''.join(pattern) % ((None,)*(10 if withN50 else 7) + (0,) )

    @staticmethod
    def distribSummary(distribution, nbFirstVals=5, nbLastVals=3):
        assert isinstance(distribution, dict)
        assert isinstance(list(distribution.values())[0], int)
        distribRepr = [" %s:%s" % (distribution[val], val) for val in sorted(distribution.keys())]
        if len(list(distribution.keys())) > nbFirstVals + nbLastVals:
            distribRepr = distribRepr[:nbFirstVals] + ["..."] + distribRepr[-nbLastVals:]
        return ','.join(distribRepr)

    # Renvoie la moyenne ponderee d'une liste
    ##############################
    @staticmethod
    def weightedMean(lst):
        sV = 0.
        sP = 0.
        for (val,poids) in lst:
            sV += val*poids
            sP += poids
        if sP == 0:
            return 0.
        return sV/sP

    # Renvoie la correlation entre les deux variables
    ##################################################
    @staticmethod
    def correlation(x, y):
        N = len(x)
        if N != len(y):
            N = min(N, len(y))
            x = x[:N]
            y = y[:N]

        sum_sq_x = 0.
        sum_sq_y = 0.
        sum_coproduct = 0.
        mean_x = x[0]
        mean_y = y[0]
        for i in range(1,N):
            sweep = i / (i + 1.0)
            delta_x = x[i] - mean_x
            delta_y = y[i] - mean_y
            sum_sq_x += delta_x * delta_x * sweep
            sum_sq_y += delta_y * delta_y * sweep
            sum_coproduct += delta_x * delta_y * sweep
            mean_x += delta_x / (i + 1.0)
            mean_y += delta_y / (i + 1.0)
        pop_sd_x = math.sqrt( sum_sq_x / N )
        pop_sd_y = math.sqrt( sum_sq_y / N )
        cov_x_y = sum_coproduct / N
        correlation = cov_x_y / (pop_sd_x * pop_sd_y)
        return correlation


    # Renvoie la quantite de chaque valeur dans la liste
    #####################################################
    @staticmethod
    def count(l):
        import collections
        d = collections.defaultdict(int)
        for x in l:
            d[x] += 1
        return d


######################################################
# Genere les permutations de la liste l              #
# Chaque permutation est referee par un index unique #
######################################################
class permutationGenerator:

    def __init__(self, l):

        self.l = tuple(l)
        self.n = len(l)
        fact = [1] * self.n
        for i in range(self.n,1,-1):
            fact[i-2] = i * fact[i-1]
        self.fact = fact
        self.a = [None] * self.n

    def getNbPerm(self, p):
        return self.fact[self.n-p-1]

    def getPermutation(self, k, p):

        assert k >= 0
        assert k < self.fact[self.n-p-1]

        for i in range(1,self.n):
            (self.a[i],k) = divmod(k, self.fact[i])

        b = self.l[:]
        for i in range(1,self.n):
            (b[i],b[self.a[i]]) = (b[self.a[i]],b[i])

        return b[-p:]


#############################
# Fonctions d'interpolation #
#############################
class myInterpolator:

    # Interpolation lineaire (1 dimension)
    ########################################
    @staticmethod
    def oneDimLinear(points, val):
        intervals = []

        for ((u,FU), (v,FV)) in myTools.myIterator.slidingTuple(list(zip(points, val))):
            A = (FU - FV) / (u - v)
            B = FU - A * u
            intervals.append((A, B))
        intervals.append(intervals[-1])

        def tmp(x):
            i = bisect.bisect_right(points, x)
            (A, B) = intervals[i-1]
            return (A * x + B)

        return tmp


    # Interpolation cubique (1 dimension)
    #######################################
    @staticmethod
    def oneDimCubic(points, val):
        intervals = []

        deriv = [0] + [(FV-FU) / (2*(v-u)) for ((u,FU), (v,FV)) in myTools.myIterator.slidingTuple(list(zip(points, val)))] + [0]
        tangeantes = [d1+d2 for (d1,d2) in myTools.myIterator.slidingTuple(deriv)]

        # Generate a cubic spline for each interpolation interval.
        for ((u,FU,DU), (v,FV,DV)) in myTools.myIterator.slidingTuple(list(zip(points,val,tangeantes))):

            denom = float((u - v)**3)

            A = ((-DV - DU) * v + (DV + DU) * u + 2 * FV - 2 * FU) / denom
            B = -((-DV - 2 * DU) * v**2  + u * ((DU - DV) * v + 3 * FV - 3 * FU) + 3 * FV * v - 3 * FU * v + (2 * DV + DU) * u**2) / denom
            C = (- DU * v**3  + u * ((- 2 * DV - DU) * v**2  + 6 * FV * v - 6 * FU * v) + (DV + 2 * DU) * u**2 * v + DV * u**3) / denom
            D = -(u *(-DU * v**3  - 3 * FU * v**2) + FU * v**3 + u**2 * ((DU - DV) * v**2 + 3 * FV * v) + u**3 * (DV * v - FV)) / denom

            intervals.append((A, B, C, D))
        intervals.append(intervals[-1])

        def tmp(x):
            i = bisect.bisect_right(points, x)
            (A, B, C, D) = intervals[i-1]
            return ((A * x + B) * x + C) * x + D

        def c_code():
            def codeChoice(intervalList):
                n = len(intervalList)
                if n < 2:
                    return ("A=%.10e;B=%.10e;C=%.10e;D=%.10e;" % intervalList[0])
                n2 = n / 2
                return ("if (x < %.10e) {%s} else {%s}" % (points[n2], codeChoice(intervalList[:n2]), codeChoice(intervalList[n2:])))
            return ("double interpolator(double x) {double A,B,C,D;" + codeChoice(intervals) + "return ((A * x + B) * x + C) * x + D;}")

        return tmp


    # Cree un interpolateur pour chaque dimension
    ###############################################
    @staticmethod
    def getMultDim(interpolator, points, val):
        linter = [interpolator(points, [x[i] for x in val]) for i in range(len(val[0]))]
        return lambda x: tuple(inter(x) for inter in linter)


####################################
# Generateur de valeurs aleatoires #
####################################
class randomValue:

    # Renvoie un indice au hasard dans l, compte tenu des valeurs de l, qui indiquent une frequence d'apparition
    #############################################################################################################
    @staticmethod
    def bisectChooser(l):
        cumulatedLs = [0] * (len(l)+1)
        for i in range(len(l)):
            cumulatedLs[i+1] = l[i] + cumulatedLs[i]
        m = cumulatedLs[-1]
        return lambda: bisect.bisect_left(cumulatedLs, random.random() * m) - 1

    # Adapted von Mises distribution in [0,1]
    ######################################################################
    @staticmethod
    def myVonMises(mean, kappa):

        # vonmisesvariate returns a value between 0 and 2*pi. We want a value
        # between -1 and 1 (mean 0)
        vmrv = random.vonmisesvariate(math.pi, kappa) / math.pi - 1

        r = vmrv * float(1 - mean)
        # r is in [- 1 + mean, + (1 - mean)]
        r = r + mean
        # r is in [-1 + 2 * mean, 1]
        if r < 0:
            # abs(r) is in [0, 1 - 2 * mean]
            return mean * abs(r) / (1 - 2 * mean)
        elif r <= 1:
            return r
        else:
            print(r, file=sys.stderr)
            raise ValueError('this case should not happen' + "mean=%s, kappa=%s and r=%s" % (mean, kappa, r))

    # Tirage aleatoire selon une densite issue d'une loi geometrique
    ##################################################################
    @staticmethod
    def geometric(p):
        return int(math.ceil(math.log(1.0 - random.random(), 1.0 - p)))


    @staticmethod
    @myTools.memoize
    def intParetoMean(alpha, precision, niter):
        s = 0.
        n = 0
        lastm = 0
        while True:
            for _ in range(niter):
                s += int(random.paretovariate(alpha))
                n += 1
            newm = s/n
            if abs(lastm-newm) < precision:
                return newm
            lastm = newm

    @staticmethod
    def paretoAlphaFromMean(m, step=0.01, precision=0.01, niter=100000):
        alpha = int((1. + 1./(m-1.)) / step) * step
        tries = {}
        while alpha not in tries:
            res = randomValue.intParetoMean(alpha, precision, niter)
            tries[alpha] = res
            alpha += -step if res < m else step
        return min((abs(y-m),x) for (x,y) in tries.items())[1]


######################
# Applatit une liste #
######################
def flatten(lst):
    res = []
    for x in lst:
        res.extend(x)
    return res


##########################################################################
# Renvoie la racine carre d'un entier, sous forme de chaine de caractere #
##########################################################################
def sqrti(n, part, nbdec):
    if (len(n) % 2) == 1:
        n = "0" + n
    backup = len(n) // 2
    iniL = len(n) + len(part)
    s = n + (part + ("00" * nbdec))[:2*nbdec]
    maxL = len(s)
    c = 0
    p = 0
    i = 0
    while ((c != 0) or (i < iniL)) and (i<maxL):
        c = 100*c + int(s[i:i+2])
        i += 2
        for x in [9,8,7,6,5,4,3,2,1,0]:
            y = (20*p+x)*x
            if y <= c:
                break
        p = 10*p+x
        c -= y
    sp = str(p)
    return (sp[:backup], sp[backup:])

################################################################
# Fast computation of combinations for probability computation #
################################################################
@myTools.memoize
def combinations(n,p):
    if n < 0:
        raise ValueError("combination(n=%s,p=%s), n cannot be < 0" % (n,p))

    if p > n:
        return 0 # there is 0 elements with p elements
    else:
        if p == 0:
            return 1 # there is only 1 combination with 0 element, the null ensembl
        elif p == n:
            return 1 # there is only 1 combination with n elements, the ensembl that conatins n elements
        else:
            try:
                # Pascal's recurrence formula, works well with the @myTools.memoize decorator
                return float(combinations(n-1,p-1) + combinations(n-1,p))
            except:
                #Pascal's recursion fails sometimes because of a "maximum recursion depth exceeded"
                num = prod([a for a in range(n-p+1,n+1)])
                den = prod([a for a in range(1,p+1)])
                try:
                    num = float(num)
                except OverflowError:
                    num = float('inf')
                try:
                    den = float(den)
                except OverflowError:
                    den = float('inf')

                if num == float('inf') and den == float('inf'):
                    return None
                elif num == float('inf') and den != float('inf'):
                    return float('inf')
                elif num != float('inf') and den == float('inf'):
                    return 0.0
                else:
                    return float(num)/den

def prod(factors):
    return reduce(operator.mul, factors, 1)

def mean(l):
    return float(sum(l)) / len(l)

def ratioAbs(ratio):
    assert ratio > 0
    return ratio if ratio >= 1 else 1.0 / ratio
