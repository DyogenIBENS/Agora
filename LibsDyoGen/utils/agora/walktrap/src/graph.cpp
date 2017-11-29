// File: graph.cpp
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

#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <cstring>
#include "graph.h"

using namespace std;


bool operator<(const Edge& E1, const Edge& E2) {
  return(E1.neighbor < E2.neighbor);
}


Vertex::Vertex() {
  degree = 0;
  edges = 0;
  total_weight = 0.;
}

Vertex::~Vertex() {
  if(edges) delete[] edges;
}

Graph::Graph() {
  nb_vertices = 0;
  nb_edges = 0;
  vertices = 0;
  index = 0;
  total_weight = 0.;
}

Graph::~Graph () {
  if (vertices) delete[] vertices;
}

class Edge_list {  
public:
  int* V1;
  int* V2;
  float* W;

  int size;
  int size_max;
  
  void add(int v1, int v2, float w);
  Edge_list() {
    size = 0;
    size_max = 1024;
    V1 = new int[1024];
    V2 = new int[1024];
    W = new float[1024];
  }
  ~Edge_list() {
    if(V1) delete[] V1;
    if(V2) delete[] V2;
    if(W) delete[] W;
  }
};

void Edge_list::add(int v1, int v2, float w) {
  if(size == size_max) {
    int* tmp1 = new int[2*size_max];
    int* tmp2 = new int[2*size_max];
    float* tmp3 = new float[2*size_max];
    for(int i = 0; i < size_max; i++) {
      tmp1[i] = V1[i];
      tmp2[i] = V2[i];      
      tmp3[i] = W[i];
    }
    delete[] V1;
    delete[] V2;
    delete[] W;
    V1 = tmp1;
    V2 = tmp2;
    W = tmp3;
    size_max *= 2;
  }
  V1[size] = v1;
  V2[size] = v2;
  W[size] = w;
  size++;
}




istream& operator>>(istream& stream, Graph& G) {
  if(G.vertices) delete[] G.vertices;

  int nb_line = 0;
  int max_vertex = 0;
  
  Edge_list EL;

  while (!stream.eof()) {	// loop for each line of the file
    nb_line++;
    string str;
    getline(stream, str);
    if(str[0] == '#') continue;	// a comment line
    istringstream line(str);

    int v1;
    line >> v1;
    if(line.fail()) {
      if(line.eof()) continue;
      cerr << "error : unable to read line " << nb_line << " : " << str << endl;
      exit(0);
    }
    int v2;
    line >> v2;
    if(line.fail()) {
      cerr << "error : unable to read line " << nb_line << " : " << str << endl;
      exit(0);
    }
    float w;
    line >> w;
    if(line.fail()) {
      if(line.eof()) w = 1.;
      else {
	cerr << "error : unable to read line " << nb_line << " : " << str << endl;
	exit(0);
      }
    }
    if(!line.eof()) {
      cerr << "error : line " << nb_line << " too long : " << str << endl;
      exit(0);
    }

    if(v1 > max_vertex) max_vertex = v1;
    if(v2 > max_vertex) max_vertex = v2;
    if((v1 < 0) || (v2 < 0)) {
      cerr << "error : line " << nb_line << " negative vertex number : " << str << endl;
      exit(0);
    }
    if(w < 0) {
      cerr << "error : line " << nb_line << " negative weight : " << str << endl;
      exit(0);
    }
    EL.add(v1, v2, w);
  }

  G.nb_vertices = max_vertex + 1; 
  G.vertices = new Vertex[G.nb_vertices];
  G.nb_edges = 0;
  G.total_weight = 0.;

  for(int i = 0; i < EL.size; i++) {
      G.vertices[EL.V1[i]].degree++;
      G.vertices[EL.V2[i]].degree++;
      G.vertices[EL.V1[i]].total_weight += EL.W[i];
      G.vertices[EL.V2[i]].total_weight += EL.W[i];
      G.nb_edges++;
      G.total_weight += EL.W[i];
    }

  for(int i = 0; i < G.nb_vertices; i++) {
    if(G.vertices[i].degree == 0) {
      cerr << "error : degree of vertex " << i << " is 0" << endl;
      exit(0);
    }
    G.vertices[i].edges = new Edge[G.vertices[i].degree + 1];
    G.vertices[i].edges[0].neighbor = i;
    G.vertices[i].edges[0].weight = G.vertices[i].total_weight/double(G.vertices[i].degree);
    G.vertices[i].total_weight+= G.vertices[i].total_weight/double(G.vertices[i].degree);
    G.vertices[i].degree = 1;
  }
 
  for(int i = 0; i < EL.size; i++) {
    G.vertices[EL.V1[i]].edges[G.vertices[EL.V1[i]].degree].neighbor = EL.V2[i];
    G.vertices[EL.V1[i]].edges[G.vertices[EL.V1[i]].degree].weight = EL.W[i];
    G.vertices[EL.V1[i]].degree++;
    G.vertices[EL.V2[i]].edges[G.vertices[EL.V2[i]].degree].neighbor = EL.V1[i];
    G.vertices[EL.V2[i]].edges[G.vertices[EL.V2[i]].degree].weight = EL.W[i];
    G.vertices[EL.V2[i]].degree++;
  }  
  
  for(int i = 0; i < G.nb_vertices; i++)
    sort(G.vertices[i].edges, G.vertices[i].edges+G.vertices[i].degree);

  for(int i = 0; i < G.nb_vertices; i++) {  // merge multi edges
    int a = 0;
    for(int b = 1; b < G.vertices[i].degree; b++) {
      if(G.vertices[i].edges[b].neighbor == G.vertices[i].edges[a].neighbor)
	G.vertices[i].edges[a].weight += G.vertices[i].edges[b].weight;
      else 
	G.vertices[i].edges[++a] = G.vertices[i].edges[b];
    }
    G.vertices[i].degree = a+1;
  }

  return stream;
}


unsigned long Graph::memory() {
  unsigned long m = 0;
  m += (unsigned long)(nb_vertices)*sizeof(Vertex);
  m += 2*(unsigned long)(nb_edges)*sizeof(Edge);
  m += sizeof(Graph);
  if(index != 0) {
    m += (unsigned long)(nb_vertices)*sizeof(char*);
    for(int i = 0; i < nb_vertices; i++)
      m += strlen(index[i]) + 1;
  }
  return m;
}


bool Graph::load_index(char* input_file) {
  ifstream f;
  f.open(input_file, ios::in);
  if(!f.is_open()) {
    cerr << "unable to open file " << index << " , index is ignored" << endl; 
    return false;
  }
  index = new char*[nb_vertices];
  for(int i = 0; i < nb_vertices; ++i)
    index[i] = 0;
   
  int nb_line = 0;
  while (!f.eof()) {	// loop for each line of the file
    nb_line++;
    int i;
    f >> i;
    if(f.fail()) {
      if(f.eof()) break;
      cerr << "error : unable to read line " << nb_line << " , index is ignored" << endl;
      for(int i = 0; i < nb_vertices; ++i)
	if(index[i]) delete[] index[i];
      delete[] index;
      index = 0;
      return false;
    }
    if(i < 0 || i >= nb_vertices || index[i] != 0) {
      cerr << "error : invalid vertex number at line : " << nb_line << " , index is ignored" << endl;
      for(int i = 0; i < nb_vertices; ++i)
	if(index[i]) delete[] index[i];
      delete[] index;
      index = 0;
      return false;
    }
    string str;
    getline(f, str);
    if(str.length() <= 1) {
      cerr << "error : unable to read line " << nb_line << " , index is ignored" << endl;
      for(int i = 0; i < nb_vertices; ++i)
	if(index[i]) delete[] index[i];
      delete[] index;
      index = 0;
      return false;
    }
    char* c = new char[str.size()];
    strcpy(c, str.c_str()+1);
    index[i] = c;
  }
  
  for(int i = 0; i < nb_vertices; ++i)
    if(index[i] == 0) {
      cerr << "error : vertex " << i << "not found, index is ignored" << endl;
      for(int i = 0; i < nb_vertices; ++i)
	if(index[i]) delete[] index[i];
      delete[] index;
      index = 0;
      return false;
    }
  return true;
}














