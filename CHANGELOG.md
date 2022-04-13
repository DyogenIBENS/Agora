## 2022-03-14 - v3.2 (pending)

1. [bugfix] -- `ALL.reformatGeneFamilies.py` doesn't need the ancestral
   species to have genes
2. [bugfix] -- NHX trees can now miss internal labels, as AGORA can generate
   names by itself.
3. [bugfix] -- Species names in the gene trees are not transformed by AGORA
   anymore, and simply have to match the names in the species tree.
4. [new] -- `misc.compareGenomes.py` now defaults to colouring the
   chromosomes according to their size

## 2022-02-05 - v3.1

1. [new] -- The `agora-*.py` scripts now accept existing orthology groups
   as an alternative to gene trees.
2. [new] -- New `-compress` option to change the compression algorithm
   (or skip compression).

## 2021-03-12 - v3.0

1. [major change] -- Migrated to Python 3.
2. [change] -- Settled on using "blocks" for all reconstructions: no
   contigs / scaffolds any more.
3. [change] -- Removed the scripts used to convert to our own tree format
   since AGORA now natively supports Newick and NHX.
4. [change] -- `buildSynteny.integr-scaffolds.py` now outputs all the
   support scores, incl. the ones from the input blocks.
5. [new] -- Up-to-date help messages for all scripts.
6. [new] -- Option to sort the output genomes by decreasing block size.
7. [new] -- Workflow steps are now skipped if they have already been run
   (plus an option to force rerunning an entire workflow).
8. [new] -- Dropped the dependency on joblib (Python's `multiprocessing`
   now works well enough for our needs).
9. [bugfix] -- Translated some French comments.

## 2021-02-07 - v2.1

1. [change] -- Changed the parameters of gene-families filtering in the
   Vertebrates workflow to the values used for Genomicus production.
2. [new] -- Support for macOS (Catalina).
3. [new] -- Added the ability to do the second reconstruction pass
   (scaffolds) in a multi-integration fashion in the library
   `myAgoraWorkflow`.
4. [new] -- New preset to run reconstructions for Plants
   (`agora-plants.py`).
5. [new] -- New script to test all reconstruction parameters and
   automatically select the best. It should apply to any clade
   (`agora-generic.py`).

## 2020-12-07 - v2.0

1. [major change] -- Renamed all scripts and most options to match the
   upcoming publication.

## 2020-10-01 - v1.5

1. [bugfix] -- Added the missing `+onlySingletons` option to the
   `buildSynteny.integr-extend.py` call in `agora2.py`.
2. [change] -- New `-extantSpeciesFilter` parameter to control which extant
   species can be used to run pairwise comparisons
   (`buildSynteny.pairwise-conservedPairs.py` and
   `buildSynteny.integr-groups.py`). The option is exposed in all workflow
   scripts. In `buildSynteny.integr-groups.py` it replaces the third
   positional command-line argument.
3. [change] -- Removed the multithreading option of
   `buildSynteny.integr-copy.py` as it does not give any benefits.
4. [change] -- Run the workflow steps sequentially by default in the agora
   scripts, with an option to enable their parallelisation.
5. [new] -- Print memory usage stats when running workflows.
6. [new] -- Decreased the memory usage of `ALL.filterGeneFamilies-size.py`,
   `buildSynteny.pairwise-conservedPairs.py`, and `buildSynteny.integr-groups.py`.

## 2020-06-25 - v1.4

1. [bugfix] -- Avoid `ALL.extractGeneFamilies.py` crashing because of a
   `ValueError` in certain conditions.
2. [bugfix] -- Some mitochondrial genomes were excluded from the
   reconstructions.
3. [bugfix] -- `buildSynteny.integr-refine.py` doesn't crash with a
   `ZeroDivisionError` when there isn't a previous score to compare
   against.
4. [change] -- `buildSynteny.integr-refine.py` and
   `buildSynteny.integr-groups.py` can now run across multiple cores.
5. [change] -- The path evaluation parameter of
   `buildSynteny.integr-refine.py` now defaults to the one used in the
   paper.
6. [change] -- Removed two Python dependencies: `numpy` and `enum`.
7. [change] -- Unary nodes are now allowed in the species-tree.
8. [new] -- New library (`myAgoraWorkflow`) to simplify the creation of
   workflows around AGORA scripts.
9. [new] -- Print CPU usage stats when running workflows.
10. [new] -- New scripts (`agora1.py` and `agora2.py`) to run the default
    workflows.

## 2020-05-13 - v1.3

1. [change] -- Moved the _LibsDyoGen_ inside `scripts/utils/`
2. [change] -- `buildSynteny.integr-groups.py` now accepts `_` as a value
   for the `usedSpecies` parameter and interprets it as the same value as
   the `target` parameter.
3. [change] -- `ALL.extractGeneFamilies.py` now outputs trees in the same
   format as the input trees

## 2020-05-10 - v1.2

1. [bugfix] -- `agora.py` now makes the script it launches follow the
   `nbThreads` option.
2. [bugfix] -- Can now have multiple ancGenes _size_ instructions in the
   configuration file.
3. [bugfix] -- `agora.py` now applies correct dependencies for the pairwise
   comparisons
4. [change] -- Moved and renamed the scripts that convert file formats.
5. [new] -- All scripts now natively support Newick species tree and NHX
	 gene trees.
6. [new] -- New `publish` instruction to convert the ancestral genomes from
   the _diags_ format to _ancGenomes_.
7. [new] -- `agora.py` can now generate the first set of _ancGenes_ (_all_)

## 2020-04-27 - v1.1

1. [bugfix] -- `agora.py` was not detecting task failures.
2. [bugfix] -- `agora.py` was always creating an ancGenes-filter task even
   when it doesn't need one.
3. [change] -- Added a `-LOG.ancGraph` option to all the
   `buildSynteny.integr` scripts (instead of hardcoded paths).
   The `integrOutput` option of the configuration files now refer to this
   per-ancestor log file, rather than the standard output of the scripts
   (which is now empty).
4. [change] -- Relative paths for input genes and species-tree are now
   evaluated from the location of the configuration file.
5. [change] -- Following a rewrite of the HowTo, the structure of the output
   directories as defined in the configuration files has changed:
   * `diags/pairwise` &rarr; `pairwise`,
   * `diags/integr` &rarr; `integrDiags`.

## 2016-09-15 - v1.0

_Initial release_

