#!/usr/bin/env bash
#Launch all the commands in the README file and stops on errors if any
set -e
# Any subsequent commands which fail will cause the shell script to exit
# immediately
red='\e[0;31m'
green='\e[0;32m'
blue='\e[0;36m'
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
    printf "${blue}${l}${NC}\n"
    printf "${blue}${t}${NC}\n"
    printf "${blue}${l}${NC}\n"
}

error () {
    printf "${red}HOHOHO problem! $1${NC}\n"
    exit 1
}

#creates tmp/ directory
print_title 'creates tmp directory for testing'
print_and_run_commands "rm -rf tmp" "mkdir tmp"


###################################################
#	Check integrity of ALL.extractGeneFamilies.py #
###################################################
print_title 'check the ancGenes families extraction scripts'

extractGeneFamiliesCommandLines=(
"src/ALL.extractGeneFamilies.py example/data/Species.nwk example/data/GeneTreeForest.nhx.bz2 -OUT.ancGenesFiles=tmp/ancGenes/all/ancGenes.%s.list.bz2 > tmp/geneTrees.afterExtractingAncGenes.protTree"
)
print_and_run_commands "${extractGeneFamiliesCommandLines[@]}"

#########################################
#	Check integrity of agora.py		    #
#########################################

print_title 'check the agora.py encapsulated script'
agoraCommandLines=(
# agora.py
"src/agora.py conf/agora-vertebrates.ini -workingDir=tmp -nbThreads=1"
)

print_and_run_commands "${agoraCommandLines[@]}"

NbAncGenomes=`ls tmp/ancGenomes/vertebrates-workflow/ancGenome.* | wc -l`

if [ ${NbAncGenomes} != 4 ]
    then
        error 'Some ancestral genomes are missing'
fi

verifCommandLines=(
)
for i in $(seq 0 3)
do
    ANCGENOME_FILENAME="ancGenome.A${i}.list.bz2"
    verifCommandLines+=("cmp 'tmp/ancGenomes/vertebrates-workflow/$ANCGENOME_FILENAME' 'example/results/ancGenomes/vertebrates-workflow/$ANCGENOME_FILENAME'")
done

print_title 'output verification'

print_and_run_commands "${verifCommandLines[@]}"

echo
printf "${green} The ancestral genomes are available in tmp/ancGenomes/vertebrates-workflow/${NC}\n"
printf "${green} Everything seems OK! Enjoy AGORA${NC}\n"

