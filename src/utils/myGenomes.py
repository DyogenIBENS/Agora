# -*- coding: utf-8 -*-

# LibsDyogen version 1.0 (6/11/2015) -- modified for AGORA v3.0
# python v2.7 at least is needed
# Copyright © 2006-2021 IBENS/Dyogen, 2020-2021 EMBL-European Bioinformatics Institute, 2021 Genome Research Ltd : Matthieu MUFFATO, Joseph LUCAS and Hugues ROEST CROLLIUS
# mail : agora@bio.ens.psl.eu
# Licences GLP v3 and CeCILL v2

import sys
import itertools
import collections

from . import myFile
from . import myTools

Gene = collections.namedtuple("Gene", ['chromosome', 'beginning', 'end', 'strand', 'names'])
GenePosition = collections.namedtuple("GenePosition", ['chromosome', 'index'])

# convert chromosome names into integers whenever it is possible
def commonChrName(x):
    try:
        return int(x)
    except TypeError:
        return None
    except ValueError:
        return sys.intern(x)

ContigType = myTools.Enum('Chromosome', 'Scaffold', 'Unknown', 'Random')
def contigType(chrName):
    if chrName in [None, "Un_random", "UNKN", "Un"]:
        # chromosome named "Un" in Monodelphis.domestica corresponds to a scaffold that concatenates all unassembled
        # fragments
        return ContigType.Unknown
    try:
        x = int(chrName)
        if x < 100:
            return ContigType.Chromosome
        else:
            return ContigType.Scaffold
    except:
        chrNameLow = chrName.lower()
        if "rand" in chrNameLow:
            # FIXME : what is a random contig ? It seems that a random contig is always associated to a chromosome number.
            # Thus we know that this random chromosome is on a specific chromosome number but we do not know where.
            return ContigType.Random
        # ki in Homo.sapiens data81
        # jh in mus.musculus, Canis.lupus.familiaris
        # aaex in Canis.lupus.familiaris
        # aadn in Gallus.gallus
        keys = ["cont", "scaff", "ultra", "reftig", "_", "un", "gl", "ki", "jh", "aaex", "aadn"]
        for x in keys:
            if x in chrNameLow:
                return ContigType.Scaffold
        if (chrName in ["U", "E64", "2-micron"]) or chrName.endswith("Het"):
            return ContigType.Scaffold
        else:
            # sex chromosomes of the chicken as (Z, W) or sex chromosomes of the human (X, Y)
            return ContigType.Chromosome

class myFASTA:
    # load a FASTA file, one sequence at a time
    @staticmethod
    def iterFile(name):
        f = myFile.openFile(name, "r")
        name = None
        for ligne in f:
            ligne = ligne.replace('\n', '').strip()
            # chevrons indicate the beginning of a new sequence
            if ligne.startswith('>'):
                tmp = []
                if name != None:
                    yield (name, "".join(tmp))
                name = ligne[1:].strip()
            # lines must be concatenated
            elif name != None:
                tmp.append(ligne.upper())
        if name != None:
            yield (name, "".join(tmp))
        f.close()

    # load a whole FASTA file
    @staticmethod
    def loadFile(name):
        return dict(myFASTA.iterFile(name))

    # write a sequence into the FASTA format
    @staticmethod
    def printSeq(f, name, seq, length=60):
        print(">" + name, file=f)
        while len(seq) != 0:
            print(seq[:length], file=f)
            seq = seq[length:]
        print(file=f)

# count nucleotides
def getMonoNuc(seq, x1, x2, sel):
    n = 0
    gc = 0
    for x in range(x1, x2+1):
        if seq[x] == "N":
            continue
        n += 1
        if seq[x] in sel:
            gc += 1
    return (n,gc)

# count dinucleotides
# a "dinucleotide" is a single piece of DNA or RNA that is two nucleotides long
def getDiNuc(seq, x1, x2, sel):
    n = 0
    gc = 0
    for x in range(x1, x2):
        s = seq[x:x+2]
        if "N" in s:
            continue
        n += 1
        if s in sel:
            gc += 1
    return (n,gc)

codon2aa = {
     'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L', 'TCT': 'S',
     'TCC': 'S', 'TCA': 'S', 'TCG': 'S', 'TAT': 'Y', 'TAC': 'Y',
     'TGT': 'C', 'TGC': 'C', 'TGG': 'W', 'CTT': 'L', 'CTC': 'L',
     'CTA': 'L', 'CTG': 'L', 'CCT': 'P', 'CCC': 'P', 'CCA': 'P',
     'CCG': 'P', 'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
     'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R', 'ATT': 'I',
     'ATC': 'I', 'ATA': 'I', 'ATG': 'M', 'ACT': 'T', 'ACC': 'T',
     'ACA': 'T', 'ACG': 'T', 'AAT': 'N', 'AAC': 'N', 'AAA': 'K',
     'AAG': 'K', 'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
     'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V', 'GCT': 'A',
     'GCC': 'A', 'GCA': 'A', 'GCG': 'A', 'GAT': 'D', 'GAC': 'D',
     'GAA': 'E', 'GAG': 'E', 'GGT': 'G', 'GGC': 'G', 'GGA': 'G',
     'GGG': 'G', 'TAA': '*', 'TAG': '*', 'TGA': '*' }

aa2codon = collections.defaultdict(list)
for (c,p) in codon2aa.items():
    aa2codon[p].append(c)

# general class for genomes
class Genome:

    # builder
    def __init__(self, fichier, **kwargs):

        if isinstance(fichier, str):
            print("Loading genome of", fichier, "...", end=' ', file=sys.stderr)
            f = myFile.firstLineBuffer(myFile.openFile(fichier, 'r'))

            # list of genes per chromosome
            self.lstGenes = collections.defaultdict(list)
            info = False
            # choice of the loading function
            c = f.firstLine.split("\t")
            if f.firstLine.startswith(">") or f.firstLine.endswith("$"):
                # GRIMM-Synteny format
                ######################
                if f.firstLine.startswith(">"):
                    self.name = f.firstLine[1:].strip()
                chrom = 1
                for l in f:
                    l = l.strip()
                    if not l.endswith("$"):
                        continue
                    for (i,x) in enumerate(l.replace("$","").split()):
                        strand = -1 if x.startswith("-") else 1
                        self.addGene([x[1:] if x[0] in "-+" else x], chrom, i, i+1, strand)
                    chrom += 1
                print("(GRIMM)", end=' ', file=sys.stderr)

            elif len(c) == 1:
                # ancestral genes: "NAMES"
                ##########################
                for (i,l) in enumerate(f):
                    self.lstGenes[None].append( Gene(None, i, i+1, 0, tuple(sys.intern(x) for x in l.split())) )
                print("(ancestral genes)", end=' ', file=sys.stderr)

            elif (len(c) == 2) and not set(c[1]).issubset("01-"):
                # ancestral genome: "CHR NAMES"
                ###############################
                lastC = None
                for (i,l) in enumerate(f):
                    c = l.split("\t")
                    if lastC != c[0]:
                        lastC = c[0]
                        dec = i
                    self.addGene(c[1].split(), c[0], i-dec, i-dec+1, 0)
                print("(ancestral genome: chrom+noms)", end=' ', file=sys.stderr)

            elif (len(c) >= 5) and (" " not in c[3]) and (len(c[4]) > 0):
                # Ensembl: "CHR BEG END STRAND NAMES"
                #####################################
                for l in f:
                    c = l.replace('\n', '').split('\t')
                    if len(c) == 5:
                        (chr, beg, end, s, geneNames) = c
                    else:
                        assert len(c) == 6
                        (chr, beg, end, s, geneNames, transcriptNames) = c
                    s = int(s) if s != 'None' else None
                    (beg, end) = (int(beg), int(end))
                    # TODO ensure that the 5' extremity (transcription start site) is the beg column
                    # if not ((s in {None, +1} and beg < end) or (s == -1  and end < beg)):
                    #     pass
                    # assert (s in {None, +1} and beg < end) or (s == -1  and end < beg), 's=%s, beg=%s and end=%s' % (s, beg, end)
                    self.addGene(geneNames.split(), chr, beg, end, s)
                print("(Ensembl)", end=' ', file=sys.stderr)

            elif (len(c) == 4) and int(c[1]) < 2:
                # ancestral genome: "CHR STRAND LST-INDEX LST-STRANDS"
                ######################################################
                if 'ancGenes' in kwargs:
                    ancGenes = kwargs["ancGenes"].lstGenes[None]
                lastC = None
                for l in f:
                    c = l.split("\t")
                    if lastC != c[0]:
                        lastC = c[0]
                        pos = 0
                        currC = commonChrName(c[0])
                    data = list(zip([int(x) for x in c[2].split()], [int(x) for x in c[3].split()]))
                    if int(c[1]) < 0:
                        data = [(i,-s) for (i,s) in data.__reversed__()]
                    for (index,strand) in data:
                        if 'ancGenes' in kwargs:
                            self.lstGenes[currC].append( Gene(currC, pos, pos+1, strand, ancGenes[index].names) )
                        else:
                            self.lstGenes[currC].append( Gene(currC, pos, pos+1, strand, (index,)) )
                        pos += 1
                print("(ancestral genome: chrom+diags)", end=' ', file=sys.stderr)

            else:
                # Optional column for adjacencies support scores
                ilss = None
                if len(c) == 2:
                    (ili,ils) = (0,1)
                else:
                    assert len(c) >= 4
                    (ili,ils) = (2,3)
                    ilss = 4
                    self.ancName = c[0]

                if 'ancGenes' in kwargs:
                    ancGenes = kwargs["ancGenes"].lstGenes[None]

                if ilss:
                    self.support = {}

                # ancestral genome: "LST-INDEX LST-STRANDS"
                #############################################
                for (i,l) in enumerate(f):
                    c = l.split("\t")
                    chrom = i+1
                    lchrom = self.lstGenes[chrom]
                    if ilss:
                        self.support[chrom] = c[ilss].split()
                    for (pos,(index,strand)) in enumerate(zip(c[ili].split(), c[ils].split())):
                        if 'ancGenes' in kwargs:
                            lchrom.append( Gene(chrom, pos, pos+1, int(strand), ancGenes[int(index)].names) )
                        else:
                            lchrom.append( Gene(chrom, pos, pos+1, int(strand), (int(index),) ) )
                print("(ancestral genome: diags)", end=' ', file=sys.stderr)

            f.close()
            self.name = fichier

        else:
            genomeBase = fichier
            print("Filtering of", genomeBase.name, "...", end=' ', file=sys.stderr)
            filterIn = set(kwargs["filterIn"]) if "filterIn" in kwargs else None
            filterOut = set(kwargs["filterOut"]) if "filterOut" in kwargs else None

            def filt(gene):
                if filterIn is not None:
                    return any(s in filterIn for s in gene.names)
                if filterOut is not None:
                    return all(s not in filterOut for s in gene.names)
                return True

            self.lstGenes = {}
            for (chrom,l) in genomeBase.lstGenes.items():
                l = [gene for gene in l if filt(gene)]
                if len(l) > 0:
                    self.lstGenes[chrom] = l
            self.name = "Filter from " + genomeBase.name
            print("%d genes -> %d genes" % (sum(len(x) for x in genomeBase.lstGenes.values()), sum(len(x) for x in self.lstGenes.values())), end=' ', file=sys.stderr)

        self.init(**kwargs)
        print("OK", file=sys.stderr)
        if info:
            print('INFO: because genes may be indexed by their 5\' (which is the right choice) extremity, gene.end=3\' < gene.beg=5\' may be true', file=sys.stderr)

    # initialise the dictionary and the list of chromosomes
    def init(self, **kwargs):

        self.dicGenes = {}
        self.chrList = collections.defaultdict(list)
        self.chrSet = collections.defaultdict(set)
        withDict = kwargs.get("withDict", True)

        for chrom in self.lstGenes:

            # associate a gene name to its location in the genome
            self.lstGenes[chrom].sort()
            if withDict:
                for (i,gene) in enumerate(self.lstGenes[chrom]):
                    for s in gene.names:
                        self.dicGenes[s] = GenePosition(chrom,i)

            # classification into chromosomes/contigs/scaffolds
            t = contigType(chrom)
            self.chrList[t].append(chrom)
            self.chrSet[t].add(chrom)

        # sort chromosome names
        for t in self.chrList:
            self.chrList[t].sort(key=lambda c: (isinstance(c, str), c))


    # add a gene to the genome
    def addGene(self, names, chromosome, beg, end, strand):
        # FIXME : 0 ?
        assert strand in {-1, 0, 1, None}
        chromosome = commonChrName(chromosome)
        self.lstGenes[chromosome].append(Gene(chromosome, beg, end, strand, tuple(sys.intern(s) for s in names)) )

    # return genes in the chromosome that are between beg and end base pairs
    def getGenesAt(self, chr, beg, end, onlyInside=False):
        if chr not in self.lstGenes:
            return
        lst = self.lstGenes[chr]

        # dichotomizing search of a gene with a non-null overlap with the window
        def dichotFind(a, b):
            i = (a+b)//2
            if a == i: # si b == a+1
                return a
            if lst[i].end < beg:
                return dichotFind(i,b)
            elif lst[i].beginning > end:
                return dichotFind(a,i)
            else:
                return i
        index = dichotFind(0, len(lst)-1)

        res = []
        for i in range(index-1, -1, -1):
            g = lst[i]
            if g.end < beg:
                break
            #if g.beginning <= end: # every time true
            assert g.beginning <= end
            if (not onlyInside) or ((beg <= g.beginning) and (g.end <= end)):
                res.append(g)
        for g in reversed(res):
            yield g

        for i in range(index, len(lst)):
            g = lst[i]
            if g.beginning > end:
                break
            #if g.end >= beg: # every time true
            assert g.end >= beg
            if (not onlyInside) or ((beg <= g.beginning) and (g.end <= end)):
                yield g

    # return genes in the vicinity of another given gene (window of + or - 'l' genes around the given gene)
    def getGenesNear(self, chr, index, l):
        if chr not in self.lstGenes:
            return []

        return self.lstGenes[chr][max(0, index-l):index+l+1]

    # return all genes
    def __iter__(self):
        for t in self.lstGenes.values():
            for g in t:
                yield g

    # search gene locations given by its names
    def getPositions(self, names):
        return set((self.dicGenes[s] for s in names if s in self.dicGenes))

    # return other names of a gene (assuming that this gene is only present at one location)
    def getOtherNames(self, name):
        if name not in self.dicGenes:
            return []
        (c,i) = self.dicGenes[name]
        return [x for x in self.lstGenes[c][i].names if x != name]

    # build the list of orthologs between the two genomes
    # FIXME
    def buildOrthosTable(self, chr1, genome2, chr2, includeGaps, genesAnc):
        chr2 = set(chr2)
        # Tous les orthologues entre pour les chromosomes OK du genome 1
        res = {}
        for c1 in chr1:
            res[c1] = []
            for (i1,g1) in enumerate(self.lstGenes[c1]):
                tmp = genome2.getPositions(g1.names)
                if genesAnc != None:
                    for (c,i) in genesAnc.getPositions(g1.names):
                        tmp.update(genome2.getPositions(genesAnc.lstGenes[c][i].names))

                # On ne garde que les chromosomes OK du genome 2
                tmp = [(c,i) for (c,i) in tmp if c in chr2]
                # +/- includeGaps
                if includeGaps or (len(tmp) > 0):
                    res[c1].append( (i1,tmp) )
        return res

    # print the genome (Ensembl format)
    def printEnsembl(self, f):
        for chrom in self.lstGenes:
            for gene in self.lstGenes[chrom]:
                print(myFile.myTSV.printLine([chrom, gene.beginning, gene.end, gene.strand, " ".join(gene.names)]), file=f)

    def intoDict(self):
        # Warning : it is important to use a dict since there are sometimes a
        # jump in the numerotation of chromosomes in a genome.
        newGenome = {}
        for c in list(self.lstGenes.keys()):
            assert len(self.lstGenes[c]) >=1
            newGenome[c] = [(g.names[0],g.strand) for g in self.lstGenes[c]]
        return newGenome

    def computeMeanInterGeneLen(self, setConsideredGeneNames):
        chromosomesMeanInterConsideredGeneLen = {}
        lGeneCoord = []
        totalLenChrom = 0
        totalLenConsideredGenes = 0
        for c in self.chrList[ContigType.Chromosome]:
            lenConsideredGenes = 0
            nbConsideredGenes = 0
            for g in self.lstGenes[c]:
                assert g.beginning < g.end
                lGeneCoord.extend([g.beginning, g.end])
                assert len(g.names) == 1, g.names
                if g.names[0] in setConsideredGeneNames:
                    lenConsideredGenes += g.end - g.beginning
                    nbConsideredGenes += 1
            # convert in Mb
            minC = float(min(lGeneCoord)) / 1000000
            maxC = float(max(lGeneCoord)) / 1000000
            lenChrom = maxC - minC
            print(c, file=sys.stderr)
            chromosomesMeanInterConsideredGeneLen[str(c)] = float(lenChrom) / lenConsideredGenes if lenConsideredGenes > 0 else float(lenChrom)
            totalLenChrom += lenChrom if lenChrom > 0 else None
            totalLenConsideredGenes += lenConsideredGenes
        if totalLenChrom:
            densityInConsideredGenes = totalLenConsideredGenes / totalLenChrom
        else:
            densityInConsideredGenes = None
        return (chromosomesMeanInterConsideredGeneLen, densityInConsideredGenes)
