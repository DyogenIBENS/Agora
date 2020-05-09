#!/usr/bin/env bash
#Launch all the commands in the README file and stops on errors if any
set -e
# Any subsequent commands which fail will cause the shell script to exit
# immediately
red='\e[0;31m'
green='\e[0;32m'
NC='\e[0m' # No Color

print_and_run_commands () {
    for line in "$@"
    do
        printf "${green}${line}${NC}\n"
        eval "$line"
    done
}

# from https://stackoverflow.com/questions/5349718/how-can-i-repeat-a-character-in-bash
print_title () {
    t=$1
    l=$(eval "printf '%0.s=' {1..${#t}}")
    printf "${red}${l}${NC}\n"
    printf "${red}${t}${NC}\n"
    printf "${red}${l}${NC}\n"
}

error () {
    printf "${red}HOHOHO problem! $1${NC}\n"
    exit 1
}

#creates tmp/ directory
print_title 'creates tmp directory for testing'
print_and_run_commands "rm -rf tmp" "mkdir tmp"

#############################################
#	Check integrity of pre-processing scripts #
#############################################
print_title 'check the preprocessing scripts'

preProcessCommandLines=(
# convet a .nhx tree into a protTree (forest of gene trees)
"src/convert.geneTrees.NHX-phylTree.py example/data/GeneTreeForest.nhx.bz2 > tmp/geneTrees.protTree"
# convet a .nwk tree into a phylTree
"src/convert.speciesTree.Newick-phylTree.py example/data/Species.nwk > tmp/speciesTree.phylTree"
)
print_and_run_commands "${preProcessCommandLines[@]}"

###################################################
#	Check integrity of ALL.extractGeneFamilies.py #
###################################################
print_title 'check the ancGenes families extraction scripts'

extractGeneFamiliesCommandLines=(
"src/ALL.extractGeneFamilies.py tmp/speciesTree.phylTree tmp/geneTrees.protTree -OUT.ancGenesFiles=tmp/ancGenes/all/ancGenes.%s.list.bz2 > tmp/geneTrees.afterExtractingAncGenes.protTree"
)
print_and_run_commands "${extractGeneFamiliesCommandLines[@]}"

#########################################
#	Check integrity of agora.py		    #
#########################################
print_title 'creation of the configuration file'

sed  s,../example/data/Species.conf,speciesTree.phylTree, ./conf/agora-robust.ini > tmp/agora-robust.ini

print_title 'check the agora.py encapsulated script'
agoraCommandLines=(
# agora.py
"src/agora.py tmp/agora-robust.ini -workingDir=tmp -nbThreads=1"
)

print_and_run_commands "${agoraCommandLines[@]}"

NbAncGenomes=`ls tmp/ancGenomes/robust/ancGenome.* | wc -l`

if [ ${NbAncGenomes} != 4 ]
    then
        error 'Some ancestral genomes are missing'
fi

verifCommandLines=(
)
for i in $(seq 0 3)
do
    ANCGENOME_FILENAME="ancGenome.A${i}.list.bz2"
    verifCommandLines+=("cmp 'tmp/ancGenomes/robust/$ANCGENOME_FILENAME' 'example/results/ancGenomes/robust/$ANCGENOME_FILENAME'")
done

print_title 'output verification'

print_and_run_commands "${verifCommandLines[@]}"

echo
printf "${green} The ancestral genomes are available in tmp/ancGenomes/robust/${NC}\n"
printf "${green} Everything seems OK! Enjoy AGORA${NC}\n"

