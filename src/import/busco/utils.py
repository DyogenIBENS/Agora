#!/usr/bin/env python3

import collections
import sys

import newick


# Load a genome from a BED file (5 columns)
def load_genome_5(path):
	genome = collections.defaultdict(list)
	with open(path, "r") as fh:
		for l in fh:
			t = l.split()
			genome[t[0]].append( (t[4], int(t[3])) )
	return list(genome.items())


def genome_to_genes(genome):
    genes = []
    for (c, gene_ids) in genome:
        for (i, s) in enumerate(gene_ids):
            gene_id = s[0].split('_')[-1]
            genes.append(gene_id)
    return set(genes)


# Helper class that allows searching for the LCA of any list of species
class AncestorSearch:

    def __init__(self, filename: str) -> None:
        self.ref_species = []
        self.anc_desc = {}
        with open(filename) as fh:
            ss = fh.readline().strip().rstrip(";").rstrip()
            self.tree = newick.buildTree(ss)
            self._populate_tree(self.tree)
        self.species_index = {s: i for (i, s) in enumerate(self.ref_species)}

    # Print the underlying tree as Newick, with names assigned to the ancestors
    def as_str(self):
        def _rec_as_str(d):
            if 'children' in d:
                return '(' + ','.join([_rec_as_str(x) for x in d['children']]) + ')' + d['label']
            else:
                return d['label']
        return _rec_as_str(self.tree) + ';'

    # Populate a matrix, telling for each node in which sub-tree (outgroup / child 1 / child 2) each species is
    # Also assigns a unique name to each ancestor
    def _populate_tree(self, d):
        if 'children' in d:
            if len(d['children']) == 1:
                return self._populate_tree(d['children'][0])
            # Before calling _populate_tree
            n = len(self.ref_species)
            labels_list = [self._populate_tree(c) for c in d['children']]
            # Pad the species we've already seen with zeros (all the extra species are outgroups)
            extra_zeros = [0] * sum([len(x) for x in labels_list])
            for li in self.anc_desc.values():
                li.extend(extra_zeros)
            anc_name = 'anc_%d' % (len(self.anc_desc)+1)
            #print(d['label'], file=sys.stderr)
            if d['label']:
                if '[' in d['label']:
                    # Remove NHX tags
                    d['label'] = d['label'].split('[')[0]
                # Keep the branch length
                d['label'] = anc_name + ":" + d['label'].split(':')[1]
            else:
                d['label'] = anc_name
            self.anc_desc[anc_name] = [0] * n
            ll = []
            for (i, l) in enumerate(labels_list):
                self.anc_desc[anc_name] += [i+1] * len(l)
                ll.extend(l)
            return ll
        else:
            l = d['label'].split(':')[0]
            self.ref_species.append(l)
            return [l]

    def in_anc(self, species_set, anc_name):
        ind = [self.anc_desc[anc_name][self.species_index[s]] for s in species_set]
        if len(set(ind)) >= 2:
            # print(species_set, anc_name, ind, [s for (i, s) in zip(ind, species_set) if i])
            return [s for (i, s) in zip(ind, species_set) if i]
        else:
            return []

    # Return the name of the LCA of a given species set
    def anc_level(self, species_set):
        ind = [self.species_index[s] for s in species_set]
        cand = None
        for (a, l) in self.anc_desc.items():
            s = set([l[i] for i in ind])
            # The LCA is the only ancestor that includes all species (none is an outgroup) in multiple children
            if (0 not in s) and len(s) >= 2:
                assert len(species_set) >= 2
                return a
        # Single-species sets are not captured by the above
        assert len(species_set) == 1

if __name__ == '__main__':
    ans = AncestorSearch(sys.argv[1])
    ans.ref_species = ['human', 'chimp', 'macaque', 'mouse', 'rat', 'dog', 'cow', 'monodelphis', 'chicken']
    ans.species_index = {s: i for (i, s) in enumerate(ans.ref_species)}
    ans.anc_desc = {
        'homin':   [1, 2, 0, 0, 0, 0, 0, 0, 0],
        'catar':   [1, 1, 2, 0, 0, 0, 0, 0, 0],
        'rodent':  [0, 0, 0, 1, 2, 0, 0, 0, 0],
        'euarch':  [1, 1, 1, 2, 2, 0, 0, 0, 0],
        'laura':   [0, 0, 0, 0, 0, 1, 2, 0, 0],
        'boreo':   [1, 1, 1, 1, 1, 2, 2, 0, 0],
        'theria':  [1, 1, 1, 1, 1, 1, 1, 2, 0],
        'amniota': [1, 1, 1, 1, 1, 1, 1, 1, 2],
    }
    print(ans.anc_level(['human', 'mouse']))
