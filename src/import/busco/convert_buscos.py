#!/usr/bin/env python3

import pathlib
import sys

import utils


# Make the gene name unique for the genes.* files
def get_unique_gene_name(species, busco_id):
    return f'{species}.{busco_id}'

if __name__ == '__main__':

    # Check command-line parameters
    if len(sys.argv) not in [3, 4]:
        print("ERROR. Usage: %s /path/to/busco_directory /path/to/species_tree.newick [/path/to/output_directory]" % sys.argv[0], file=sys.stderr)
        sys.exit(1)
    busco_dir = pathlib.Path(sys.argv[1])
    assert busco_dir.is_dir()
    species_tree_path = pathlib.Path(sys.argv[2])
    assert species_tree_path.is_file()
    if len(sys.argv) == 4:
        output_directory = sys.argv[3]
    else:
        output_directory = "."
    output_directory = pathlib.Path(output_directory)
    assert not output_directory.is_file()

    # Load the BUSCO annotations
    all_species = []
    busco_ids = {}
    all_busco_ids = set()
    genes_path = output_directory / 'genes'
    genes_path.mkdir(parents=True, exist_ok=True)
    for busco_file in busco_dir.glob('*.tsv'):
        species = busco_file.name.split('/')[-1].split('.')[0]
        all_species.append(species)
        print("Loading", species, "from", busco_file, "...", end=" ")
        busco_ids[species] = set()
        with open(busco_file) as busco_fh:
            with open(genes_path / f'genes.{species}.list', 'w') as genes_fh:
                for line in busco_fh:
                    if line.startswith("#"):
                        continue
                    t = line[:-1].split("\t")
                    if t[1] != "Complete":
                        continue
                    assert len(t) >= 6, t
                    chr_name = t[2].split(":")[0]
                    print("\t".join([chr_name, t[3], t[4], "1" if t[5] == "+" else "-1", get_unique_gene_name(species, t[0])]), file=genes_fh)
                    busco_ids[species].add(t[0])
        all_busco_ids.update(busco_ids[species])
        print(len(busco_ids[species]), "genes OK")

    # Load the species tree, and reprint it with pseudo ancestral species names
    print("Loading species tree ...", end=" ")
    ans = utils.AncestorSearch(species_tree_path)
    with open(output_directory / 'species_tree.nh', 'w') as fh:
        print(ans.as_str(), file=fh)
    if set(all_species) != set(ans.ref_species):
        print("ERROR: species mismatch", file=sys.stderr)
        print("Species missing from the species tree:", set(all_species).difference(ans.ref_species), file=sys.stderr)
        print("Species in the species tree, but without BUSCO annotation:", set(ans.ref_species).difference(all_species), file=sys.stderr)
        sys.exit(1)
    print("OK")

    # Summary of all the gene families (not used by AGORA, but convenient to have)
    print("Grouping the BUSCO genes ...", end=" ")
    anc_genes = []
    with open(output_directory / 'families.txt', 'w') as fh:
        for busco_id in sorted(all_busco_ids):
            species = [species_name for species_name in all_species if busco_id in busco_ids[species_name]]
            species.sort()
            anc_genes.append( (busco_id,species) )
            anc_name = ans.anc_level(species)
            # Print the family ID, the LCA, and the list of species that have it
            print(busco_id, anc_name, " ".join(species), sep="\t", file=fh)
    print(len(anc_genes), "families OK")

    orthology_groups_path = output_directory / 'orthologyGroups'
    orthology_groups_path.mkdir(parents=True, exist_ok=True)

    # Create the orthology groups files as required by AGORA
    # for ancestors
    for anc_name in ans.anc_desc:
        with open(orthology_groups_path / f'orthologyGroups.{anc_name}.list', 'w') as fh:
            print("Writing", fh.name, "...", end=" ")
            n = 0
            for (busco_id, species) in anc_genes:
                sub_species = ans.in_anc(species, anc_name)
                if sub_species:
                    print(f'FAMILY.{busco_id}', *[get_unique_gene_name(s, busco_id) for s in sub_species], file=fh)
                    n += 1
            print(n, "groups OK")

    # and extant species
    for species_name in all_species:
        with open(orthology_groups_path / f'orthologyGroups.{species_name}.list', 'w') as fh:
            print("Writing", fh.name, "...", end=" ")
            n = 0
            for gene_id in sorted(busco_ids[species_name]):
                print(f'FAMILY.{gene_id}', get_unique_gene_name(species_name, gene_id), file=fh)
                n += 1
            print(n, "groups OK")


