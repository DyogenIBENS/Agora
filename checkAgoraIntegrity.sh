#!/usr/bin/env bash
#Launch all the commands in the README file and stops on errors if any
set -e
# Any subsequent commands which fail will cause the shell script to exit
# immediately
red='\e[0;31m'
green='\e[0;32m'
NC='\e[0m' # No Color

#creates tmp/ directory
printf "${red}---------------------------------${NC}\n"
printf "${red}creates tmp directory for testing${NC}\n"
printf "${red}---------------------------------${NC}\n"
if [ -d tmp ]
    then
        printf "${green}mkdir tmp${NC}\n"
        rm -r tmp
        mkdir tmp
    else
        printf "${green}mkdir tmp${NC}\n"
        mkdir tmp
fi



#############################################
#	Check integrity of pre-processing scripts #
#############################################
printf "${red}-------------------------------${NC}\n"
printf "${red}check the preprocessing scripts${NC}\n"
printf "${red}-------------------------------${NC}\n"

preProcessCommandLines=(
# convet a .nhx tree into a protTree (forest of gene trees)
"src/preprocessing/nhxGeneTrees2phylTreeGeneTrees.py example/data/GeneTreeForest.nhx.bz2 > tmp/geneTrees.protTree"
# convet a .nwk tree into a phylTree
"src/preprocessing/newickSpeciesTree2phylTreeSpeciesTree.py example/data/Species.nwk > tmp/speciesTree.phylTree"
)
for line in "${preProcessCommandLines[@]}"
	do
		printf "${green}${line}${NC}\n"
		eval ${line}
done

###################################################
#	Check integrity of ALL.extractGeneFamilies.py #
###################################################
printf "${red}----------------------------------------------${NC}\n"
printf "${red}check the ancGenes families extraction scripts${NC}\n"
printf "${red}----------------------------------------------${NC}\n"

extractGeneFamiliesCommandLines=(
"src/ALL.extractGeneFamilies.py tmp/speciesTree.phylTree tmp/geneTrees.protTree -OUT.ancGenesFiles=tmp/ancGenes/all/ancGenes.%s.list.bz2 > tmp/geneTrees.afterExtractingAncGenes.protTree"
)
for line in "${extractGeneFamiliesCommandLines[@]}"
	do
		printf "${green}${line}${NC}\n"
		eval ${line}
done
#########################################
#	Check integrity of agora.py		    #
#########################################
printf "${red}----------------------------------${NC}\n"
printf "${red}creation of the configuration file${NC}\n"
printf "${red}----------------------------------${NC}\n"


sed  s,PATH_To_SpeciesTree,${PWD}/tmp/speciesTree.phylTree, ./conf/agora-size.ini | sed s,PATH_To_GENEFILES,${PWD}/example/data/genes, > tmp/agora-size.ini



printf "${red}--------------------------------------${NC}\n"
printf "${red}check the agora.py encapsulated script${NC}\n"
printf "${red}--------------------------------------${NC}\n"
agoraCommandLines=(
# agora.py
"src/agora.py tmp/agora-size.ini -workingDir=tmp"
)

for line in "${agoraCommandLines[@]}"
	do
		printf "${green}${line}${NC}\n"
		eval ${line}
done

NbAncGenomes=`ls tmp/diags/integr/final/anc/ | wc -l`

if [ ${NbAncGenomes} == 4 ]
    then
        printf "${red} Everything seems OK! Enjoy AGORA${NC}\n"
    else
        printf "${red} HOHOHO problems! ${NC}\n"
        exit 1
fi

#########################################
#	Check integrity of postprocessing   #
#########################################

printf "${red}--------------------------------${NC}\n"
printf "${red}check the postprocessing script ${NC}\n"
printf "${red}--------------------------------${NC}\n"

mkdir tmp/ancGenomes
for i in tmp/diags/integr/final/anc/*
    do
        prefix="tmp/diags/integr/final/anc/diags."
        ancGenome=${i/.list.bz2/}; ancGenome=${ancGenome#$prefix};
        printf "${green}formating ${ancGenome}${NC}\n"
        printf "${green}src/postprocessing/misc.convertContigsToGenome.py $i tmp/ancGenes/all/ancGenes.${ancGenome}.list.bz2 > tmp/ancGenomes/ancGenome.${ancGenome}.list${NC}\n"
        src/postprocessing/misc.convertContigsToGenome.py $i tmp/ancGenes/all/ancGenes.${ancGenome}.list.bz2 > tmp/ancGenomes/ancGenome.${ancGenome}.list
    done

printf "${red} the ancestral genomes should be available in tmp/ancGenomes/${NC}\n"

