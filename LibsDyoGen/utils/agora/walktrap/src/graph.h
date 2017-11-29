// File: graph.h
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


#ifndef GRAPH_H
#define GRAPH_H
#include <iostream>

using namespace std;

class Edge {			// code an edge of a given vertex
public:
  int neighbor;			// the number of the neighbor vertex
  float weight;			// the weight of the edge
};
bool operator<(const Edge& E1, const Edge& E2);


class Vertex {
public:
  Edge* edges;			// the edges of the vertex
  int degree;			// number of neighbors
  float total_weight;		// the total weight of the vertex

  Vertex();			// creates empty vertex
  ~Vertex();			// destructor
};

class Graph {
public:
  int nb_vertices;		// number of vertices
  int nb_edges;			// number of edges
  float total_weight;		// total weight of the edges
  Vertex* vertices;		// array of the vertices

  unsigned long memory();	// the total memory used in Bytes
  Graph();			// create an empty graph
  ~Graph();			// destructor
  char** index;			// to keep the real name of the vertices
  bool load_index(char* input_file);
};

istream& operator>>(istream& in, Graph& G);	// get a graph from a stream


#endif

