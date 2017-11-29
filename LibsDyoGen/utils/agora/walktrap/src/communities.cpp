// File: communities.cpp
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

#include "communities.h"
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <cmath>

//#include <vector>

int Probabilities::length = 0;
Communities* Probabilities::C = 0;
float* Probabilities::tmp_vector1 = 0;
float* Probabilities::tmp_vector2 = 0;
int* Probabilities::id = 0;
int* Probabilities::vertices1 = 0;
int* Probabilities::vertices2 = 0;
int Probabilities::current_id = 0;


Neighbor::Neighbor() {
  next_community1 = 0;
  previous_community1 = 0;
  next_community2 = 0;	     
  previous_community2 = 0;
  heap_index = -1;
}

Probabilities::~Probabilities() {
  C->memory_used -= memory();
  if(P) delete[] P;
  if(vertices) delete[] vertices;
}

Probabilities::Probabilities(int community) {
  Graph* G = C->G;
  int nb_vertices1 = 0;
  int nb_vertices2 = 0;

  float initial_proba = 1./float(C->communities[community].size);
  int last =  C->members[C->communities[community].last_member];  
  for(int m = C->communities[community].first_member; m != last; m = C->members[m]) {
    tmp_vector1[m] = initial_proba;
    vertices1[nb_vertices1++] = m;
  }
  
  for(int t = 0; t < length; t++) {
    current_id++;    
    if(nb_vertices1 > (G->nb_vertices/2)) {
      nb_vertices2 = G->nb_vertices;
      for(int i = 0; i < G->nb_vertices; i++)
	tmp_vector2[i] = 0.;
      if(nb_vertices1 == G->nb_vertices) {
	for(int i = 0; i < G->nb_vertices; i++) {
	  float proba = tmp_vector1[i]/G->vertices[i].total_weight;
  	  for(int j = 0; j < G->vertices[i].degree; j++)
	    tmp_vector2[G->vertices[i].edges[j].neighbor] += proba*G->vertices[i].edges[j].weight;
	}
      }
      else {
	for(int i = 0; i < nb_vertices1; i++) {
	  int v1 = vertices1[i];
	  float proba = tmp_vector1[v1]/G->vertices[v1].total_weight;
  	  for(int j = 0; j < G->vertices[v1].degree; j++)
	    tmp_vector2[G->vertices[v1].edges[j].neighbor] += proba*G->vertices[v1].edges[j].weight;
	}
      }
    }
    else {
      nb_vertices2 = 0;
      for(int i = 0; i < nb_vertices1; i++) {
	int v1 = vertices1[i];
        float proba = tmp_vector1[v1]/G->vertices[v1].total_weight;
        for(int j = 0; j < G->vertices[v1].degree; j++) {
	  int v2 = G->vertices[v1].edges[j].neighbor;
	  if(id[v2] == current_id)
  	    tmp_vector2[v2] += proba*G->vertices[v1].edges[j].weight;
  	  else {
	    tmp_vector2[v2] = proba*G->vertices[v1].edges[j].weight;
	    id[v2] = current_id;
	    vertices2[nb_vertices2++] = v2;
	  }
        }
      }
    }
    float* tmp = tmp_vector2;
    tmp_vector2 = tmp_vector1;
    tmp_vector1 = tmp;

    int* tmp2 = vertices2;
    vertices2 = vertices1;
    vertices1 = tmp2;

    nb_vertices1 = nb_vertices2;
  }

  if(nb_vertices1 > (G->nb_vertices/2)) {
    P = new float[G->nb_vertices];
    size = G->nb_vertices;
    vertices = 0;
    if(nb_vertices1 == G->nb_vertices) {
      for(int i = 0; i < G->nb_vertices; i++)
	P[i] = tmp_vector1[i]/sqrt(G->vertices[i].total_weight);
    }
    else {
      for(int i = 0; i < G->nb_vertices; i++)
	P[i] = 0.;
      for(int i = 0; i < nb_vertices1; i++)
	P[vertices1[i]] = tmp_vector1[vertices1[i]]/sqrt(G->vertices[vertices1[i]].total_weight);
    }
  }
  else {
    P = new float[nb_vertices1];
    size = nb_vertices1;
    vertices = new int[nb_vertices1];
    int j = 0;
    for(int i = 0; i < G->nb_vertices; i++) {
      if(id[i] == current_id) {
	P[j] = tmp_vector1[i]/sqrt(G->vertices[i].total_weight);
	vertices[j] = i;
	j++;
      }
    }
  }
  C->memory_used += memory();
}

Probabilities::Probabilities(int community1, int community2) {
  // The two following probability vectors must exist.
  // Do not call this function if it is not the case.
  Probabilities* P1 = C->communities[community1].P;
  Probabilities* P2 = C->communities[community2].P;
  
  float w1 = float(C->communities[community1].size)/float(C->communities[community1].size + C->communities[community2].size);
  float w2 = float(C->communities[community2].size)/float(C->communities[community1].size + C->communities[community2].size);


  if(P1->size == C->G->nb_vertices) {
    P = new float[C->G->nb_vertices];
    size = C->G->nb_vertices;
    vertices = 0;
    
    if(P2->size == C->G->nb_vertices) {	// two full vectors
      for(int i = 0; i < C->G->nb_vertices; i++)
	P[i] = P1->P[i]*w1 + P2->P[i]*w2;
    }
    else {  // P1 full vector, P2 partial vector
      int j = 0;
      for(int i = 0; i < P2->size; i++) {
	for(; j < P2->vertices[i]; j++)
	  P[j] = P1->P[j]*w1;
	P[j] = P1->P[j]*w1 + P2->P[i]*w2;
	j++;
      }
      for(; j < C->G->nb_vertices; j++)
	P[j] = P1->P[j]*w1;
    }
  }
  else {
    if(P2->size == C->G->nb_vertices) { // P1 partial vector, P2 full vector
      P = new float[C->G->nb_vertices];
      size = C->G->nb_vertices;
      vertices = 0;

      int j = 0;
      for(int i = 0; i < P1->size; i++) {
	for(; j < P1->vertices[i]; j++)
	  P[j] = P2->P[j]*w2;
	P[j] = P1->P[i]*w1 + P2->P[j]*w2;
	j++;
      }
      for(; j < C->G->nb_vertices; j++)
	P[j] = P2->P[j]*w2;
    }
    else {  // two partial vectors
      int i = 0;
      int j = 0;
      int nb_vertices1 = 0;
      while((i < P1->size) && (j < P2->size)) {
	if(P1->vertices[i] < P2->vertices[j]) {
	  tmp_vector1[P1->vertices[i]] = P1->P[i]*w1;
	  vertices1[nb_vertices1++] = P1->vertices[i];
	  i++;
	  continue;
	}
	if(P1->vertices[i] > P2->vertices[j]) {
	  tmp_vector1[P2->vertices[j]] = P2->P[j]*w2;
	  vertices1[nb_vertices1++] = P2->vertices[j];
	  j++;
	  continue;
	}
	tmp_vector1[P1->vertices[i]] = P1->P[i]*w1 + P2->P[j]*w2;
	vertices1[nb_vertices1++] = P1->vertices[i];
	i++;
	j++;
      }
      if(i == P1->size) {
	for(; j < P2->size; j++) {
	  tmp_vector1[P2->vertices[j]] = P2->P[j]*w2;
	  vertices1[nb_vertices1++] = P2->vertices[j];
	}
      }
      else {
	for(; i < P1->size; i++) {
	  tmp_vector1[P1->vertices[i]] = P1->P[i]*w1;
	  vertices1[nb_vertices1++] = P1->vertices[i];
	}
      }

      if(nb_vertices1 > (C->G->nb_vertices/2)) {
	P = new float[C->G->nb_vertices];
	size = C->G->nb_vertices;
	vertices = 0;
	for(int i = 0; i < C->G->nb_vertices; i++)
	  P[i] = 0.;
	for(int i = 0; i < nb_vertices1; i++)
	  P[vertices1[i]] = tmp_vector1[vertices1[i]];
      }
      else {
	P = new float[nb_vertices1];
	size = nb_vertices1;
	vertices = new int[nb_vertices1];
	for(int i = 0; i < nb_vertices1; i++) {
	  vertices[i] = vertices1[i];
	  P[i] = tmp_vector1[vertices1[i]];
	}
      }
    }
  }

  C->memory_used += memory();
}

double Probabilities::compute_distance(const Probabilities* P2) const {
  double r = 0.;
  if(vertices) {
    if(P2->vertices) {  // two partial vectors
      int i = 0;
      int j = 0;
      while((i < size) && (j < P2->size)) {
	if(vertices[i] < P2->vertices[j]) {
	  r += P[i]*P[i];
	  i++;
	  continue;
	}
	if(vertices[i] > P2->vertices[j]) {
	  r += P2->P[j]*P2->P[j];
	  j++;
	  continue;
	}
	r += (P[i] - P2->P[j])*(P[i] - P2->P[j]);
	i++;
	j++;
      }
      if(i == size) {
	for(; j < P2->size; j++)
	  r += P2->P[j]*P2->P[j];
      }
      else {
	for(; i < size; i++)
	  r += P[i]*P[i];
      }
    }
    else {  // P1 partial vector, P2 full vector 

      int i = 0;
      for(int j = 0; j < size; j++) {
	for(; i < vertices[j]; i++)
	  r += P2->P[i]*P2->P[i];
	r += (P[j] - P2->P[i])*(P[j] - P2->P[i]);
	i++;
      }
      for(; i < P2->size; i++)
	r += P2->P[i]*P2->P[i];      
    }
  }
  else {
    if(P2->vertices) {  // P1 full vector, P2 partial vector
      int i = 0;
      for(int j = 0; j < P2->size; j++) {
	for(; i < P2->vertices[j]; i++)
	  r += P[i]*P[i];
	r += (P[i] - P2->P[j])*(P[i] - P2->P[j]);
	i++;
      }
      for(; i < size; i++)
	r += P[i]*P[i];
    }
    else {  // two full vectors
      for(int i = 0; i < size; i++)
	r += (P[i] - P2->P[i])*(P[i] - P2->P[i]);
    }
  }
  return r;
}

unsigned long Probabilities::memory() {
  if(vertices)
    return (sizeof(Probabilities) + (unsigned long)(size)*(sizeof(float) + sizeof(int)));
  else
    return (sizeof(Probabilities) + (unsigned long)(size)*sizeof(float));
}

Community::Community() {
  P = 0;
  first_neighbor = 0;
  last_neighbor = 0;
  sub_community_of = -1;
  sub_communities[0] = -1;
  sub_communities[1] = -1;
  sigma = 0.;
  internal_weight = 0.;
  total_weight = 0.;
}

Community::~Community() {
  if(P) delete P;
}


Communities::Communities(Graph* graph, int random_walks_length, bool s, unsigned long m) {
  silent = s;
  max_memory = m;
  memory_used = 0;
  G = graph;
  alpha_min = 0;
  alpha_max = 0;
  sub_communities = 0;
  sorted_alpha_min = 0;
  
  Probabilities::C = this;
  Probabilities::length = random_walks_length;
  Probabilities::tmp_vector1 = new float[G->nb_vertices];
  Probabilities::tmp_vector2 = new float[G->nb_vertices];
  Probabilities::id = new int[G->nb_vertices];
  for(int i = 0; i < G->nb_vertices; i++) Probabilities::id[i] = 0;
  Probabilities::vertices1 = new int[G->nb_vertices];
  Probabilities::vertices2 = new int[G->nb_vertices];
  Probabilities::current_id = 0;

  
  members = new int[G->nb_vertices];  
  for(int i = 0; i < G->nb_vertices; i++)
    members[i] = -1;

  H = new Neighbor_heap(G->nb_edges);
  communities = new Community[2*G->nb_vertices];

// init the n single vertex communities

  if(max_memory > 0)
    min_delta_sigma = new Min_delta_sigma_heap(G->nb_vertices*2);
  else min_delta_sigma = 0;
  
  for(int i = 0; i < G->nb_vertices; i++) {
    communities[i].this_community = i;
    communities[i].first_member = i;
    communities[i].last_member = i;
    communities[i].size = 1;
    communities[i].sub_community_of = 0;
  }

  nb_communities = G->nb_vertices;
  nb_active_communities = G->nb_vertices;

  if(!silent) cerr << "computing random walks and the first distances:";
  for(int i = 0; i < G->nb_vertices; i++)
    for(int j = 0; j < G->vertices[i].degree; j++)
      if (i < G->vertices[i].edges[j].neighbor) {
	communities[i].total_weight += G->vertices[i].edges[j].weight/2.;
	communities[G->vertices[i].edges[j].neighbor].total_weight += G->vertices[i].edges[j].weight/2.;
	Neighbor* N = new Neighbor;
	N->community1 = i;
	N->community2 = G->vertices[i].edges[j].neighbor;
	N->delta_sigma = -1./double(min(G->vertices[i].degree,  G->vertices[G->vertices[i].edges[j].neighbor].degree));
	N->weight = G->vertices[i].edges[j].weight;
	N->exact = false;
	add_neighbor(N);
      }

  if(max_memory > 0) {
    memory_used += min_delta_sigma->memory();
    memory_used += 2*(unsigned long)(G->nb_vertices)*sizeof(Community);
    memory_used += (unsigned long)(G->nb_vertices)*(2*sizeof(float) + 3*sizeof(int)); // the static data of Probabilities class
    memory_used += H->memory() + (unsigned long)(G->nb_edges)*sizeof(Neighbor);
    memory_used += G->memory();    
  }

  int c = 0;
  Neighbor* N = H->get_first();  
  while(!N->exact) {
    update_neighbor(N, compute_delta_sigma(N->community1, N->community2));
    N->exact = true;
    N = H->get_first();
    if(max_memory > 0) manage_memory();
    if(!silent) {
      c++;
      for(int k = (500*(c-1))/G->nb_edges + 1; k <= (500*c)/G->nb_edges; k++) {
	if(k % 50 == 1) {cerr.width(2); cerr << endl << k/ 5 << "% ";}
	cerr << ".";
      }
    }
  }
  
  if(!silent) cerr << endl << endl;

}

Communities::~Communities() {
  delete[] members;
  delete[] communities;
  delete H;
  if(min_delta_sigma) delete min_delta_sigma;

  if(alpha_min) delete[] alpha_min;
  if(alpha_max) delete[] alpha_max;
  if(sub_communities) {
    for(int i = 0; i < nb_communities; i++)
      if(sub_communities[i]) delete[] sub_communities[i];
    delete[] sub_communities;
  }
  
  delete[] Probabilities::tmp_vector1;
  delete[] Probabilities::tmp_vector2;
  delete[] Probabilities::id;
  delete[] Probabilities::vertices1;
  delete[] Probabilities::vertices2;
}

float Community::min_delta_sigma() {
  float r = 1.;
  for(Neighbor* N = first_neighbor; N != 0;) {
    if(N->delta_sigma < r) r = N->delta_sigma;
    if(N->community1 == this_community)
      N = N->next_community1;
    else
      N = N->next_community2;
  }
  return r;
}


void Community::add_neighbor(Neighbor* N) { // add a new neighbor at the end of the list
  if (last_neighbor) {
    if(last_neighbor->community1 == this_community)
      last_neighbor->next_community1 = N;
    else
      last_neighbor->next_community2 = N;
    
    if(N->community1 == this_community)
      N->previous_community1 = last_neighbor;
    else
      N->previous_community2 = last_neighbor;
  }
  else {
    first_neighbor = N;
    if(N->community1 == this_community)
      N->previous_community1 = 0;
    else
      N->previous_community2 = 0;
  }
  last_neighbor = N;
}

void Community::remove_neighbor(Neighbor* N) {	// remove a neighbor from the list
  if (N->community1 == this_community) {
    if(N->next_community1) {
//      if (N->next_community1->community1 == this_community)
	N->next_community1->previous_community1 = N->previous_community1;
//      else 
//	N->next_community1->previous_community2 = N->previous_community1;
    }
    else last_neighbor = N->previous_community1;
    if(N->previous_community1) {
      if (N->previous_community1->community1 == this_community)
	N->previous_community1->next_community1 = N->next_community1;
      else 
	N->previous_community1->next_community2 = N->next_community1;
    }
    else first_neighbor = N->next_community1;
  }
  else {
    if(N->next_community2) {
      if (N->next_community2->community1 == this_community)
	N->next_community2->previous_community1 = N->previous_community2;
      else 
	N->next_community2->previous_community2 = N->previous_community2;
    }
    else last_neighbor = N->previous_community2;
    if(N->previous_community2) {
//      if (N->previous_community2->community1 == this_community)
//	N->previous_community2->next_community1 = N->next_community2;
//      else 
	N->previous_community2->next_community2 = N->next_community2;
    }
    else first_neighbor = N->next_community2;
  }
}

void Communities::remove_neighbor(Neighbor* N) {
  communities[N->community1].remove_neighbor(N);
  communities[N->community2].remove_neighbor(N);
  H->remove(N);

  if(max_memory > 0) {
    if(N->delta_sigma == min_delta_sigma->delta_sigma[N->community1]) {
      min_delta_sigma->delta_sigma[N->community1] = communities[N->community1].min_delta_sigma();
      if(communities[N->community1].P) min_delta_sigma->update(N->community1);
    }

    if(N->delta_sigma == min_delta_sigma->delta_sigma[N->community2]) {
      min_delta_sigma->delta_sigma[N->community2] = communities[N->community2].min_delta_sigma();
      if(communities[N->community2].P) min_delta_sigma->update(N->community2);
    }
  }
}

void Communities::add_neighbor(Neighbor* N) {
  communities[N->community1].add_neighbor(N);
  communities[N->community2].add_neighbor(N);
  H->add(N);

  if(max_memory > 0) {
    if(N->delta_sigma < min_delta_sigma->delta_sigma[N->community1]) {
      min_delta_sigma->delta_sigma[N->community1] = N->delta_sigma;
      if(communities[N->community1].P) min_delta_sigma->update(N->community1);
    }

    if(N->delta_sigma < min_delta_sigma->delta_sigma[N->community2]) {
      min_delta_sigma->delta_sigma[N->community2] = N->delta_sigma;
      if(communities[N->community2].P) min_delta_sigma->update(N->community2);
    }
  }
}

void Communities::update_neighbor(Neighbor* N, float new_delta_sigma) {
  if(max_memory > 0) {
    if(new_delta_sigma < min_delta_sigma->delta_sigma[N->community1]) {
      min_delta_sigma->delta_sigma[N->community1] = new_delta_sigma;
      if(communities[N->community1].P) min_delta_sigma->update(N->community1);
    }
   
    if(new_delta_sigma < min_delta_sigma->delta_sigma[N->community2]) {
      min_delta_sigma->delta_sigma[N->community2] = new_delta_sigma;
      if(communities[N->community2].P) min_delta_sigma->update(N->community2);
    }

    float old_delta_sigma = N->delta_sigma;
    N->delta_sigma = new_delta_sigma;
    H->update(N);

    if(old_delta_sigma == min_delta_sigma->delta_sigma[N->community1]) {
      min_delta_sigma->delta_sigma[N->community1] = communities[N->community1].min_delta_sigma();
      if(communities[N->community1].P) min_delta_sigma->update(N->community1);
    }

    if(old_delta_sigma == min_delta_sigma->delta_sigma[N->community2]) {
      min_delta_sigma->delta_sigma[N->community2] = communities[N->community2].min_delta_sigma();
      if(communities[N->community2].P) min_delta_sigma->update(N->community2);
    }
  }
  else {
    N->delta_sigma = new_delta_sigma;
    H->update(N);
  }
}

void Communities::manage_memory() {
  while((memory_used > max_memory) && !min_delta_sigma->is_empty()) {
    int c = min_delta_sigma->get_max_community();
    delete communities[c].P;
    communities[c].P = 0;
    min_delta_sigma->remove_community(c);
  }  
}



void Communities::merge_communities(Neighbor* merge_N) {
  int c1 = merge_N->community1;
  int c2 = merge_N->community2;
  
  communities[nb_communities].first_member = communities[c1].first_member;	// merge the 
  communities[nb_communities].last_member = communities[c2].last_member;	// two lists   
  members[communities[c1].last_member] = communities[c2].first_member;		// of members

  communities[nb_communities].size = communities[c1].size + communities[c2].size;
  communities[nb_communities].this_community = nb_communities;
  communities[nb_communities].sub_community_of = 0;
  communities[nb_communities].sub_communities[0] = c1;
  communities[nb_communities].sub_communities[1] = c2;
  communities[nb_communities].total_weight = communities[c1].total_weight + communities[c2].total_weight;
  communities[nb_communities].internal_weight = communities[c1].internal_weight + communities[c2].internal_weight + merge_N->weight;
  communities[nb_communities].sigma = communities[c1].sigma + communities[c2].sigma + merge_N->delta_sigma;
  
  communities[c1].sub_community_of = nb_communities;
  communities[c2].sub_community_of = nb_communities;

// update the new probability vector...
  
  if(communities[c1].P && communities[c2].P) communities[nb_communities].P = new Probabilities(c1, c2);

  if(communities[c1].P) {
    delete communities[c1].P; 
    communities[c1].P = 0;
    if(max_memory > 0) min_delta_sigma->remove_community(c1);
  }
  if(communities[c2].P) {
    delete communities[c2].P;
    communities[c2].P = 0;
    if(max_memory > 0) min_delta_sigma->remove_community(c2);
  }

  if(max_memory > 0) {
    min_delta_sigma->delta_sigma[c1] = -1.;		    // to avoid to update the min_delta_sigma for these communities
    min_delta_sigma->delta_sigma[c2] = -1.;		    // 
    min_delta_sigma->delta_sigma[nb_communities] = -1.;
  }
  
// update the new neighbors
// by enumerating all the neighbors of c1 and c2

  Neighbor* N1 = communities[c1].first_neighbor;
  Neighbor* N2 = communities[c2].first_neighbor;

  while(N1 && N2) { 
    int neighbor_community1;
    int neighbor_community2;
    
    if (N1->community1 == c1) neighbor_community1 = N1->community2;
    else neighbor_community1 = N1->community1;
    if (N2->community1 == c2) neighbor_community2 = N2->community2;
    else neighbor_community2 = N2->community1;

    if (neighbor_community1 < neighbor_community2) {
      Neighbor* tmp = N1;
      if (N1->community1 == c1) N1 = N1->next_community1;
      else N1 = N1->next_community2;
      remove_neighbor(tmp);
      Neighbor* N = new Neighbor;
      N->weight = tmp->weight;
      N->community1 = neighbor_community1;
      N->community2 = nb_communities;
      N->delta_sigma = (double(communities[c1].size+communities[neighbor_community1].size)*tmp->delta_sigma + double(communities[c2].size)*merge_N->delta_sigma)/(double(communities[c1].size+communities[c2].size+communities[neighbor_community1].size));//compute_delta_sigma(neighbor_community1, nb_communities);
      N->exact = false;
      delete tmp;
      add_neighbor(N); 
    }
    
    if (neighbor_community2 < neighbor_community1) {
      Neighbor* tmp = N2;
      if (N2->community1 == c2) N2 = N2->next_community1;
      else N2 = N2->next_community2;
      remove_neighbor(tmp);
      Neighbor* N = new Neighbor;
      N->weight = tmp->weight;
      N->community1 = neighbor_community2;
      N->community2 = nb_communities;
      N->delta_sigma = (double(communities[c1].size)*merge_N->delta_sigma + double(communities[c2].size+communities[neighbor_community2].size)*tmp->delta_sigma)/(double(communities[c1].size+communities[c2].size+communities[neighbor_community2].size));//compute_delta_sigma(neighbor_community2, nb_communities);
      N->exact = false;
      delete tmp;
      add_neighbor(N); 
    }
    
    if (neighbor_community1 == neighbor_community2) {
      Neighbor* tmp1 = N1;
      Neighbor* tmp2 = N2;
      bool exact = N1->exact && N2->exact;
      if (N1->community1 == c1) N1 = N1->next_community1;
      else N1 = N1->next_community2;
      if (N2->community1 == c2) N2 = N2->next_community1;
      else N2 = N2->next_community2;
      remove_neighbor(tmp1);
      remove_neighbor(tmp2);
      Neighbor* N = new Neighbor;
      N->weight = tmp1->weight + tmp2->weight;
      N->community1 = neighbor_community1;
      N->community2 = nb_communities;
      N->delta_sigma = (double(communities[c1].size+communities[neighbor_community1].size)*tmp1->delta_sigma + double(communities[c2].size+communities[neighbor_community1].size)*tmp2->delta_sigma - double(communities[neighbor_community1].size)*merge_N->delta_sigma)/(double(communities[c1].size+communities[c2].size+communities[neighbor_community1].size));
      N->exact = exact;
      delete tmp1;
      delete tmp2;
      add_neighbor(N);
    }
  }

  
  if(!N1) {
    while(N2) {
//      double delta_sigma2 = N2->delta_sigma;
      int neighbor_community;
      if (N2->community1 == c2) neighbor_community = N2->community2;
      else neighbor_community = N2->community1;
      Neighbor* tmp = N2;
      if (N2->community1 == c2) N2 = N2->next_community1;
      else N2 = N2->next_community2;
      remove_neighbor(tmp);
      Neighbor* N = new Neighbor;
      N->weight = tmp->weight;
      N->community1 = neighbor_community;
      N->community2 = nb_communities;
      N->delta_sigma = (double(communities[c1].size)*merge_N->delta_sigma + double(communities[c2].size+communities[neighbor_community].size)*tmp->delta_sigma)/(double(communities[c1].size+communities[c2].size+communities[neighbor_community].size));//compute_delta_sigma(neighbor_community, nb_communities);
      N->exact = false;
      delete tmp;
      add_neighbor(N);
    }
  }
  if(!N2) {
    while(N1) {
//      double delta_sigma1 = N1->delta_sigma;
      int neighbor_community;
      if (N1->community1 == c1) neighbor_community = N1->community2;
      else neighbor_community = N1->community1;
      Neighbor* tmp = N1;
      if (N1->community1 == c1) N1 = N1->next_community1;
      else N1 = N1->next_community2;
      remove_neighbor(tmp);
      Neighbor* N = new Neighbor;
      N->weight = tmp->weight;
      N->community1 = neighbor_community;
      N->community2 = nb_communities;
      N->delta_sigma = (double(communities[c1].size+communities[neighbor_community].size)*tmp->delta_sigma + double(communities[c2].size)*merge_N->delta_sigma)/(double(communities[c1].size+communities[c2].size+communities[neighbor_community].size));//compute_delta_sigma(neighbor_community, nb_communities);
      N->exact = false;
      delete tmp;
      add_neighbor(N);
    }
  }

  if(max_memory > 0) {
    min_delta_sigma->delta_sigma[nb_communities] = communities[nb_communities].min_delta_sigma();
    min_delta_sigma->update(nb_communities);
  } 

  nb_communities++;
  nb_active_communities--;
}

double Communities::merge_nearest_communities() {
  Neighbor* N = H->get_first();  
  while(!N->exact) {
    update_neighbor(N, compute_delta_sigma(N->community1, N->community2));
    N->exact = true;
    N = H->get_first();
    if(max_memory > 0) manage_memory();
  }

  double d = N->delta_sigma;
  remove_neighbor(N);

  merge_communities(N);
  if(max_memory > 0) manage_memory();

  delete N;

  if(!silent) {
    for(int k = (500*(G->nb_vertices - nb_active_communities - 1))/(G->nb_vertices-1) + 1; k <= (500*(G->nb_vertices - nb_active_communities))/(G->nb_vertices-1); k++) {
      if(k % 50 == 1) {cerr.width(2); cerr << endl << k/ 5 << "% ";}
      cerr << ".";
    }
  }
  return d;
}

double Communities::compute_delta_sigma(int community1, int community2) {
  if(!communities[community1].P) {
    communities[community1].P = new Probabilities(community1);
    if(max_memory > 0) min_delta_sigma->update(community1);
  }
  if(!communities[community2].P) {
    communities[community2].P = new Probabilities(community2);
    if(max_memory > 0) min_delta_sigma->update(community2);
  }
  
  return communities[community1].P->compute_distance(communities[community2].P)*double(communities[community1].size)*double(communities[community2].size)/double(communities[community1].size + communities[community2].size);
}

/*
void Communities::print_community(int c) {
  cout << "community " << c << " = {";
  for(int m = communities[c].first_member; m != members[communities[c].last_member]; m = members[m]) {
    if(G->index) cout << G->index[m];
    else cout << m;
    if(members[m] != members[communities[c].last_member]) cout << ", ";
  }
  cout << "}" << endl;
}*/


void Communities::print_community(int c) {
  for(int m = communities[c].first_member; m != members[communities[c].last_member]; m = members[m]) {
    if(G->index) cout << G->index[m];
    else cout << m;
    cout << " " << c << endl;
  }
}

void Communities::compute_hierarchy(int mode) {
  Hierarchy* hierarchy = new Hierarchy[G->nb_vertices];
  compute_hierarchy(nb_communities - 1, hierarchy, mode);
  
  vector<int> list;
  bool* B = new bool[nb_communities];
  if(alpha_min) delete[] alpha_min;
  alpha_min = new float[nb_communities];
  if(alpha_max) delete[] alpha_max;
  alpha_max = new float[nb_communities];
  if(sub_communities) {
    for(int i = 0; i < nb_communities; i++)
      if(sub_communities[i]) delete[] sub_communities[i];
    delete[] sub_communities;
  }
  sub_communities = new int*[nb_communities];
  if(sorted_alpha_min) delete[] sorted_alpha_min;
  sorted_alpha_min = new int[nb_communities+1];
  
  for(int i = 0; i < G->nb_vertices; i++) {
    B[i] = true;
    alpha_min[i] = hierarchy[0].alpha;
    alpha_max[i] = 0.;
    sub_communities[i] = 0;
    sorted_alpha_min[i] = i;
  }
  for(int i = G->nb_vertices; i < nb_communities; i++) {
    B[i] = false;
    alpha_min[i] = 0.;
    alpha_max[i] = 0.;
    sub_communities[i] = 0;
    sorted_alpha_min[i] = -1;		    // -1 will be used to mark the end of the vector
  }
  sorted_alpha_min[nb_communities] = -1;    // -1 will be used to mark the end of the vector

  int compteur = G->nb_vertices;
  for(int i = 0; hierarchy[i].alpha != 1.; i++) {
    int c = hierarchy[i+1].community;
    alpha_min[c] = hierarchy[i].alpha;
    sorted_alpha_min[compteur++] = c;
    list.clear();
    find_sub_communities(c, B, list);
    sub_communities[c] = new int[list.size() + 1];
    sub_communities[c][0] = list.size();
    for(unsigned int j = 0; j < list.size(); j++) {
      sub_communities[c][j+1] = list[j];
      alpha_max[list[j]] = hierarchy[i].alpha;
    }
    B[c] = true;
  }
  alpha_max[nb_communities - 1] = alpha_min[nb_communities - 1];    // the largest community = G

  delete[] B;
  delete[] hierarchy;
}

void Communities::find_sub_communities(int community, bool* B, vector<int>& list) {
  if(B[community]) {
    list.push_back(community);
    return;
  }
  find_sub_communities(communities[community].sub_communities[0], B, list);
  find_sub_communities(communities[community].sub_communities[1], B, list);
}


void Communities::compute_hierarchy(int community, Hierarchy* hierarchy, int mode) {
  if(communities[community].size == 1) {   
    hierarchy[0].alpha = 1.;
    switch(mode) {
      case 1:	// sigma
	hierarchy[0].s = -communities[community].sigma/communities[nb_communities-1].sigma;
	hierarchy[0].g = -1./float(G->nb_vertices);
//	hierarchy[0].g = communities[community].internal_weight/G->total_weight;
	break;
      case 2:	// modularity
	hierarchy[0].s = -communities[community].total_weight*communities[community].total_weight/G->total_weight;
	hierarchy[0].g = communities[community].internal_weight;
	break;
      case 3:	// performance
	hierarchy[0].s = communities[community].size * (G->nb_vertices - communities[community].size) * G->total_weight/G->nb_edges - (communities[community].total_weight - communities[community].internal_weight);
	hierarchy[0].g = 2*communities[community].internal_weight;
	break;
    }
    hierarchy[0].community = -1;
    return;
  }

  Hierarchy* tmpH1 = new Hierarchy[communities[communities[community].sub_communities[0]].size];
  Hierarchy* tmpH2 = new Hierarchy[communities[communities[community].sub_communities[1]].size];
  
  compute_hierarchy(communities[community].sub_communities[0], tmpH1, mode);
  compute_hierarchy(communities[community].sub_communities[1], tmpH2, mode);

  int a = 0;
  int b = 0;
  int current_community = -1;

  while((tmpH1[a].alpha != 1.) || (tmpH2[b].alpha != 1.)) {
    if(tmpH1[a].alpha < tmpH2[b].alpha) {
      hierarchy[a + b].alpha = tmpH1[a].alpha;
      hierarchy[a + b].s = tmpH1[a].s + tmpH2[b].s;
      hierarchy[a + b].g = tmpH1[a].g + tmpH2[b].g;
      hierarchy[a + b].community = current_community;
      current_community = tmpH1[++a].community;
    }
    else {
      hierarchy[a + b].alpha = tmpH2[b].alpha;
      hierarchy[a + b].s = tmpH1[a].s + tmpH2[b].s;
      hierarchy[a + b].g = tmpH1[a].g + tmpH2[b].g;
      hierarchy[a + b].community = current_community;
      current_community = tmpH2[++b].community;
    }
  }
  
  hierarchy[a + b].s = tmpH1[a].s + tmpH2[b].s;
  hierarchy[a + b].g = tmpH1[a].g + tmpH2[b].g;
  hierarchy[a + b].community = current_community;

  int c = a + b;
  float ds = 0;
  float dg = 0;
  while(true) {
    switch (mode) {
      case 1:
	ds = hierarchy[c].s - (-communities[community].sigma/communities[nb_communities-1].sigma);
	dg = -1./float(G->nb_vertices) - hierarchy[c].g;
//	dg = communities[community].internal_weight/G->total_weight - hierarchy[c].g;
	break;
      case 2:
	ds = hierarchy[c].s - (-communities[community].total_weight*communities[community].total_weight/G->total_weight);
	dg = communities[community].internal_weight - hierarchy[c].g;
	break;
      case 3:
	ds = hierarchy[c].s - (communities[community].size * (G->nb_vertices - communities[community].size) * G->total_weight/G->nb_edges - (communities[community].total_weight - communities[community].internal_weight));
	dg = 2*communities[community].internal_weight - hierarchy[c].g;
	break;
    }
    if ((c == 0) || (ds/(dg+ds) > hierarchy[c-1].alpha)) break;
    c--;
  }
  
  hierarchy[c].alpha = ds/(dg+ds);
  hierarchy[c + 1].alpha = 1.;

  switch (mode) {
    case 1:
      hierarchy[c + 1].s = -communities[community].sigma/communities[nb_communities-1].sigma;
    hierarchy[c + 1].g = -1./float(G->nb_vertices);
//      hierarchy[c + 1].g = communities[community].internal_weight/G->total_weight;
      break;
    case 2:
      hierarchy[c + 1].s = -communities[community].total_weight*communities[community].total_weight/G->total_weight;
      hierarchy[c + 1].g = communities[community].internal_weight;
      break;
    case 3:
      hierarchy[c + 1].s = communities[community].size * (G->nb_vertices - communities[community].size) * G->total_weight/G->nb_edges - (communities[community].total_weight - communities[community].internal_weight);
      hierarchy[c + 1].g = 2*communities[community].internal_weight;
      break;
  }

  hierarchy[c + 1].community = community;

  delete[] tmpH1;
  delete[] tmpH2;  
}

void Communities::print_hierarchy(int detail) {
  if(!alpha_min) return;

  if(detail >= 1) {
    for(int i = G->nb_vertices; sorted_alpha_min[i] != -1; i++) {
      int c = sorted_alpha_min[i];
      cout << alpha_min[c] << " : ";
      for(int j = 1; j < sub_communities[c][0]; j++)
	cout << sub_communities[c][j] << " + ";
      cout << sub_communities[c][sub_communities[c][0]] << " --> " << c << endl;
      if(detail >= 4) print_community(c);
    }
  }
}

void Communities::print_partition(float alpha) {
  if(!alpha_min) return;
//  cout << endl << "Partition for scale factor alpha = " << alpha << " :" << endl;

  for(int i = 0; i < nb_communities; i++)
    if(alpha_min[i] <= alpha && alpha_max[i] > alpha)
      print_community(i);
//  cout << endl;
}

void Communities::find_best_partition(float ratio, map<float, float>& M, int detail) {
  if(!alpha_min) return;

  float last_alpha = 0.;
  double r = 0.;
  double d1 = 0.;
  double d2 = 0.;

  int s = 200;		// samples : 0 .. 1 by steps 1/s
  int w = int(0.05*s);	// window size for local maxima

  float* A = new float[s];
  float* Hq = new float[s];
  for(int i = 0 ; i < s; i++) {
    Hq[i] = 0.;
    A[i] = (float(i) + 0.5)/float(s);
  }

  for(int i = 0; sorted_alpha_min[i] != -1; i++) {
    int c = sorted_alpha_min[i];

    if(alpha_min[c] > last_alpha) {
      if(r + d1*last_alpha + d2 * last_alpha * last_alpha > Hq[int(last_alpha * s)]) {
	A[int(last_alpha * s)] = last_alpha;
	Hq[int(last_alpha * s)] = r + d1*last_alpha + d2 * last_alpha * last_alpha;
      }
      if(d2 != 0.) {
	double alpha_extremum = -d1/(2.*d2);
	if (last_alpha < alpha_extremum && alpha_extremum < alpha_min[c]) 
	  if(r + d1*alpha_extremum + d2 * alpha_extremum * alpha_extremum > Hq[int(alpha_extremum * s)]) {
	    A[int(alpha_extremum * s)] = alpha_extremum;
	    Hq[int(alpha_extremum * s)] = r + d1*alpha_extremum + d2 * alpha_extremum * alpha_extremum;
	}
      }

      for(int alpha_interm = int(last_alpha * s) + 1; alpha_interm < int(alpha_min[c] * s); alpha_interm++)
	  if(r + d1*(float(alpha_interm) + 0.5)/s + d2 * (float(alpha_interm) + 0.5)/s * (float(alpha_interm) + 0.5)/s > Hq[alpha_interm]) {
	    A[alpha_interm] = (float(alpha_interm) + 0.5)/s;
	    Hq[alpha_interm] = r + d1*(float(alpha_interm) + 0.5)/s + d2 * (float(alpha_interm) + 0.5)/s * (float(alpha_interm) + 0.5)/s;
	  }
      
      if(r + d1*alpha_min[c] + d2 * alpha_min[c] * alpha_min[c] > Hq[int(alpha_min[c] * s)]) {
	A[int(alpha_min[c] * s)] = alpha_min[c];
	Hq[int(alpha_min[c] * s)] = r + d1*alpha_min[c] + d2 * alpha_min[c] * alpha_min[c];
      }
      last_alpha = alpha_min[c];
    }
    
    if(alpha_min[c] < alpha_max[c]) {
      r += (alpha_max[c]-alpha_min[c]) * double(communities[c].size)/double(G->nb_vertices) * ratio;

      r -= 4.*(alpha_max[c]*alpha_min[c])/(alpha_max[c]-alpha_min[c]) * double(communities[c].size)/double(G->nb_vertices) * (1 - ratio);
      d1 += 4.*(alpha_max[c]+alpha_min[c])/(alpha_max[c]-alpha_min[c]) * double(communities[c].size)/double(G->nb_vertices)  * (1 - ratio);
      d2 -= 4./(alpha_max[c]-alpha_min[c]) * double(communities[c].size)/double(G->nb_vertices)  * (1 - ratio);
    }

    if(sub_communities[c])
      for(int j = 1; j <= sub_communities[c][0]; j++) {
	int c2 = sub_communities[c][j];
    	if(alpha_min[c2] < alpha_max[c2]) {
  	  r -= (alpha_max[c2]-alpha_min[c2]) * double(communities[c2].size)/double(G->nb_vertices) * ratio;
	  
	  r += 4.*(alpha_max[c2]*alpha_min[c2])/(alpha_max[c2]-alpha_min[c2]) * double(communities[c2].size)/double(G->nb_vertices) * (1 - ratio);
	  d1 -= 4.*(alpha_max[c2]+alpha_min[c2])/(alpha_max[c2]-alpha_min[c2]) * double(communities[c2].size)/double(G->nb_vertices)  * (1 - ratio);
	  d2 += 4./(alpha_max[c2]-alpha_min[c2]) * double(communities[c2].size)/double(G->nb_vertices)  * (1 - ratio);
	}
      }
  }

  if (detail >= 3) {
    cout << endl << "Scale Factor Relevance:" << endl; 
    cout << "Alpha\tRelevance" << endl; 
    for(int i = 0; i < s; i++)
      cout << setprecision(5) << A[i] << "\t" << setprecision(5) << Hq[i] << endl;
    cout << endl;
  }
 
  for(int i = 0 ; i < s; i++) {
    float max_hq = -1.;
    int a_max = 0;
    for(int a = min(i+w, s-1); a >= max(i-w, 0); a--)
      if(Hq[a] > max_hq) {
	max_hq = Hq[a];
	a_max = a;
      }
    if(a_max == i && Hq[i] != 0.) {
      M.insert(pair<float,float>(Hq[i], A[i]));
    }
  }
  
  delete A;
  delete Hq;
}



