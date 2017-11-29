// File: walktrap.cpp
//-----------------------------------------------------------------------------
// Walktrap v0.2 -- Finds community structure of networks using random walks
// Copyright (C) 2004-2005 Pascal Pons
//
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
//-----------------------------------------------------------------------------
// Author   : Pascal Pons 
// Email    : pons@liafa.jussieu.fr
// Web page : http://www.liafa.jussieu.fr/~pons/
// Location : Paris, France
// Time	    : June 2005
//-----------------------------------------------------------------------------
// see readme.txt for more details

#include <Python.h>
#include "graph.h"
#include "communities.h"
#include <ctime>
#include <map>
#include <set>
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <algorithm>

using namespace std;


void loadFromPython(Graph& G, PyObject* nodes, PyObject* edges) {
	
	if(G.vertices) delete[] G.vertices;

	G.nb_vertices = PyDict_Size(nodes);
	G.vertices = new Vertex[G.nb_vertices];
	G.nb_edges = 0;
	G.total_weight = 0.;

	PyObject *X, *indX;
	Py_ssize_t posX = 0;

	while (PyDict_Next(nodes, &posX, &X, &indX)) {
		long x = PyInt_AS_LONG(indX);
		PyObject* edgesFromX = PyDict_GetItem(edges, X);
		PyObject *Y, *valEdge;
		Py_ssize_t posY = 0;

		G.vertices[x].degree = PyDict_Size(edgesFromX);
		G.vertices[x].edges = new Edge[G.vertices[x].degree + 1];
		G.vertices[x].total_weight = 0.;
		G.nb_edges += G.vertices[x].degree;
		int i = 1;

		while (PyDict_Next(edgesFromX, &posY, &Y, &valEdge)) {
			PyObject *tmp = PyDict_GetItem(nodes, Y);
			long y = PyInt_AS_LONG(tmp);
			double v = PyFloat_AS_DOUBLE(valEdge);

			G.vertices[x].edges[i].neighbor = y;
			G.vertices[x].edges[i].weight = v;
			G.vertices[x].total_weight += v;
			i += 1;
		}

		G.total_weight += G.vertices[x].total_weight;
		G.vertices[x].edges[0].neighbor = x;
		G.vertices[x].edges[0].weight = G.vertices[x].total_weight/double(G.vertices[x].degree);
		G.vertices[x].total_weight += G.vertices[x].total_weight/double(G.vertices[x].degree);
		G.vertices[x].degree += 1;
		sort(G.vertices[x].edges, G.vertices[x].edges+G.vertices[x].degree);
		PyDict_DelItem(edges, X);
	}
	G.total_weight /= 2.;
	G.nb_edges /= 2;
}


void* Communities::get_hierarchy() {

	// La liste des fusions
	PyObject* lst = PyList_New(0);
	
	if (!alpha_min)
		return lst;
		
	for(int i = G->nb_vertices; sorted_alpha_min[i] != -1; i++) {
		int c = sorted_alpha_min[i];
		// Un tuple pour la liste des noeuds fusionnes
		PyObject* lstNodes = PyTuple_New(sub_communities[c][0]);
		for(int j = 1; j <= sub_communities[c][0]; j++)
			PyTuple_SET_ITEM(lstNodes, j-1, PyInt_FromLong(sub_communities[c][j]));
		
		// Un tuple pour (l'echelle alpha - les noeuds fusionnes - le nom du noeud cree)
		PyList_Append(lst, Py_BuildValue("(fOi)", alpha_min[c], lstNodes, c));
	}
	return lst;
}


const char* infoModule = 
"WalkTrap v0.3 -- Finds community structure of networks using random walks.\n \
Copyright (C) 2004 Pascal Pons & Matthieu Muffato 2007.\n\n \
";


const char* infoFunc =
"Launch walktrap clustering on graph specified with 'nodes' and 'edges'.\n \
Options:\n\
\trandomWalksLength <int>    set the length of random walks to x (default = 4)\n \
\tverboseLevel <int>    set to x the detail level of the output (0 <= x <= 4, default = 2)\n \
\tshowProgress <bool>    display the progress\n \
\tmemoryUseLimit <unsigned long>    limit the memory usage to x MB\n \
\tqualityFunction <int>    1 = sigma, 2 = modularity, 3 = perf\n \
see readme.txt for more details\n";


static PyObject* pywalktrap_doWalktrap(PyObject *self, PyObject *args, PyObject *keywds)
{
	PyObject* nodes;
	PyObject* edges;
	int length = 4;
	int detail = 0;
	int showProgress = 0;
	unsigned long max_memory = 0;
	int quality_function = 2;	// 1 = sigma, 2 = modularity, 3 = perf
	static const char* kwlist[] = {"nodes", "edges", "randomWalksLength", "verboseLevel", "showProgress", "memoryUseLimit", "qualityFunction", NULL};

	vector<float> print_partition;
	if (!PyArg_ParseTupleAndKeywords(args, keywds, "O!O!|iiiki", const_cast<char **>(kwlist), &PyDict_Type, &nodes, &PyDict_Type, &edges, &length, &detail, &showProgress, &max_memory, &quality_function))
		return NULL;
	int silent = !showProgress;

	time_t begin = clock();
	Graph* G = new Graph;
	loadFromPython(*G, nodes, edges);
	if(!silent) cerr << G->nb_vertices << " vertices and " << G->nb_edges << " edges sucessfuly loaded." << endl << endl;

	Communities C(G, length, silent, max_memory);  

	if(!silent) cerr << "merging the communities:";

	while(!C.H->is_empty()) {
		C.merge_nearest_communities();
	}

	if(!silent) cerr << endl << endl << "computing hierarchy and scale factor relevance." << endl;

	C.compute_hierarchy(quality_function);
	C.print_hierarchy(detail);

	map<float,float> M;
	C.find_best_partition(0.5, M, detail);
	PyObject* lst = PyList_New(M.size());
	int i = 0;
	for(map<float,float>::reverse_iterator it = M.rbegin(); it != M.rend(); ++it, i++)
		PyList_SET_ITEM(lst, i, Py_BuildValue("(ff)", it->second, it->first));
	
	// Le resultat final, le tuple (liste des coupures pertinentes, dendogramme)
	PyObject* ret = Py_BuildValue("(OO)", lst, C.get_hierarchy());
    
	if (!silent)
		cerr << endl << "computation successfully terminated in " << double(clock() - begin) / double(CLOCKS_PER_SEC) << "s" << endl;
	delete G;

	return ret;
}

static PyMethodDef PyWalktrapMethods[] = {
	{"doWalktrap", (PyCFunction) pywalktrap_doWalktrap, METH_VARARGS | METH_KEYWORDS, infoFunc},
	{NULL, NULL, 0, NULL}        /* Sentinel */
};


PyMODINIT_FUNC init_walktrap(void)
{
	(void) Py_InitModule3("_walktrap", PyWalktrapMethods, infoModule);
}
