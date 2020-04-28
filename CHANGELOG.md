## 2020-04-27 - v1.1

1. [bugfix] -- `agora.py` was not detecting workflow errors
2. [change] -- Added a `-LOG.ancGraph` option to all the
   `buildSynteny.integr` scripts (instead of hardcoded paths).
   The `integrOutput` option of the configuration files now refer to this
   per-ancestor log file, rather than the standard output of the scripts
   (which is now empty).
3. [bugfix] -- `agora.py` was always creating an ancGenes-filter task even
   when it doesn't need one.
4. [change] -- Relative paths for input genes and species-tree are now
   evaluated from the location of the configuration file.
5. [change] -- Following a rewrite of the HowTo, the structure of the output
   directories as defined in the configuration files has changed:
   * `diags/pairwise` &rarr; `pairwise`
   * `diags/integr` &rarr; `integrDiags`

## 2016-09-15 - v1.0

_Initial release_

