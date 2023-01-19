# AGORA v3.1 (2022-02-05)

## Introduction

AGORA stands for “Algorithm for Gene Order Reconstruction in Ancestors” and was developed by
Matthieu Muffato in the DYOGEN Laboratory at the École normale supérieure in Paris in 2008.

```
    // | |     //   ) )  //   ) ) //   ) )  // | |
   //__| |    //        //   / / //___/ /  //__| |
  / ___  |   //  ____  //   / / / ___ (   / ___  |
 //    | |  //    / / //   / / //   | |  //    | |
//     | | ((____/ / ((___/ / //    | | //     | |
```

AGORA is used to generate ancestral genomes for the
[Genomicus](https://www.genomicus.biologie.ens.fr/genomicus) online server
for gene order comparison, and has been in constant use in the group since.

## License

This code may be freely distributed and modified under the terms of the GNU General Public License version 3 (GPL v3)
and the CeCILL licence version 2 of the CNRS. These licences are contained in the files:

1. [LICENSE-GPL.txt](LICENSE-GPL.txt) (or on [www.gnu.org](https://www.gnu.org/licenses/gpl-3.0-standalone.html))
2. [LICENCE-CeCILL.txt](LICENCE-CeCILL.txt) (or on [www.cecill.info](https://cecill.info/licences/Licence_CeCILL_V2-en.html))

Copyright for this code is held jointly by the Dyogen (DYnamic and Organisation of GENomes) team
of the Institut de Biologie de l'Ecole Normale Supérieure (IBENS) 46 rue d'Ulm Paris, the European Bioinformatics Institute outstation of the European Molecular Biology Laboratory,
Genome Research Ltd,
and the individual authors.

- Copyright © 2006-2022 IBENS/Dyogen : Alexandra LOUIS, Thi Thuy Nga NGUYEN, Matthieu MUFFATO and Hugues ROEST CROLLIUS
- Copyright © 2020-2021 EMBL-European Bioinformatics Institute
- Copyright © 2021-2022 Genome Research Ltd

## Contact

Email agora {at} bio {dot} ens {dot} psl {dot} eu

## If you use AGORA, please cite:

Matthieu Muffato, Alexandra Louis, Nga Thi Thuy Nguyen, Joseph M. Lucas, Camille Berthelot, Hugues Roest Crollius. [Reconstruction of hundreds of reference ancestral genomes across the eukaryotic kingdom.](https://doi.org/10.1038/s41559-022-01956-z) ***Nat Ecol Evol*** (Jan 2023).

## Installation

To simplify deployment, AGORA already embeds a modified version of
[LibsDyogen](https://github.com/DyogenIBENS/LibsDyogen) version 1.0
(6/11/2015), a Python library
for bioinformatics and comparative genomics developed by the same group.

AGORA is written in Python 3, which is widely available.
You can install a Python 3 environment with all the dependencies with
[_conda_](https://docs.conda.io/)

```
conda env create --file conda_env.yml
```

Alternatively you can add the required dependencies to an existing
environment (e.g. a Python _virtualenv_):

```
pip install -r requirements.txt
```

:warning: AGORA is currently not compatible with Python 3.8+ on macOS.

AGORA is compatible with [PyPy](https://www.pypy.org/) (an alternative,
faster implementation of Python) which significantly speeds up the
reconstructions, whilst using more memory.

Once everything is installed, run this to check the installation:

```
./checkAgoraIntegrity.sh
```

It should run for a few minutes and end with this message in green:

> The ancestral genomes are available in tmp/ancGenomes/

## Usage

In a nutshell, you need to gather:

* a species tree (e.g. `species-tree.nwk`)
* the list of genes of each species (e.g. matching the pattern `genes/genes.%s.list`)
* gene trees (e.g. `gene-trees.nhx`), or orthology groups for each ancestor (e.g. matching the pattern `orthologyGroups/orthologyGroups.%s.list`)

and then try one of these:

```bash
src/agora-basic.py species-tree.nwk gene-trees.nhx genes/genes.%s.list
src/agora-basic.py species-tree.nwk orthologyGroups/orthologyGroups.%s.list genes/genes.%s.list
```

If the ancestral genomes are too fragmented, run `src/agora-generic.py` instead of `src/agora-basic.py`.

Check out our [user manual](doc/HowTo.md) for more information about the
input file formats, what these two scripts do, and how to tune AGORA even
further. Also available as [docx](doc/HowTo.docx) and [pdf](doc/HowTo.pdf).
