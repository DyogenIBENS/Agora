# AGORA v1.4 (2020-06-25)

## Introduction

AGORA stands for “Algorithm for Gene Order Reconstruction in Ancestors” and was developed by
Matthieu Muffato in the DYOGEN Laboratory at the École normale supérieure in Paris.

```
    // | |     //   ) )  //   ) ) //   ) )  // | |
   //__| |    //        //   / / //___/ /  //__| |
  / ___  |   //  ____  //   / / / ___ (   / ___  |
 //    | |  //    / / //   / / //   | |  //    | |
//     | | ((____/ / ((___/ / //    | | //     | |
```

AGORA has been constantly used in the group since 2010, especially to
generate ancestral genomes for the [Genomicus](https://www.genomicus.biologie.ens.fr/genomicus)
online server for comparative genomics.

## License

This code may be freely distributed and modified under the terms of the GNU General Public License version 3 (GPL v3)
and the CeCILL licence version 2 of the CNRS. These licences are contained in the files:

1. [LICENSE-GPL.txt](LICENSE-GPL.txt) (or on [www.gnu.org](https://www.gnu.org/licenses/gpl-3.0-standalone.html))
2. [LICENCE-CeCILL.txt](LICENCE-CeCILL.txt) (or on [www.cecill.info](https://cecill.info/licences/Licence_CeCILL_V2-en.html))

Copyright for this code is held jointly by the Dyogen (DYnamic and Organisation of GENomes) team
of the Institut de Biologie de l'Ecole Normale Supérieure (IBENS) 46 rue d'Ulm Paris, the European Bioinformatics Institute outstation of the European Molecular Biology Laboratory,
and the individual authors.

- Copyright © 2006-2020 IBENS/Dyogen : Alexandra LOUIS, Thi Thuy Nga NGUYEN, Matthieu MUFFATO and Hugues ROEST CROLLIUS
- Copyright © 2020 EMBL-European Bioinformatics Institute

## Contact

Email agora {at} bio {dot} ens {dot} psl {dot} eu

## Installation

To simplify deployment, AGORA already embeds a modified version of
[LibsDyogen](https://github.com/DyogenIBENS/LibsDyogen) version 1.0
(6/11/2015), a Python library
for bioinformatics and comparative genomics developed by the same group.

AGORA is written in Python 2 and is currently not compatible with Python 3.
You can install a Python 2 environment with all the dependencies with
[_conda_](https://docs.conda.io/)

```
conda create --file conda_env.yml
```

Alternatively you can add the required dependencies to an existing
environment (e.g. a Python _virtualenv_):

```
pip install -r requirements.txt
```

Once everything is installed, run this to check the installation:

```
./checkAgoraIntegrity.sh
```

It should run for a few minutes and end with this message in green:

> The ancestral genomes are available in tmp/ancGenomes/

## Usage

In a nutshell, you need to gather:

* a species tree
* the list of genes of each species
* gene trees

and then try:

```bash
src/agora1.py species-tree.nwk gene-trees.nhx genes.%s.list
```

If the ancestral genomes are too fragmented, run this otherwise:

```bash
src/agora2.py species-tree.nwk gene-trees.nhx genes.%s.list
```

Check out our [user manual](doc/HowTo.md) for more information about the
input file formats, what these two scripts do, and how to tune AGORA even
further. Also available as [docx](doc/HowTo.docx) and [pdf](doc/HowTo.pdf).
