// File: communities.h
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


#ifndef COMMUNITIES_H
#define COMMUNITIES_H

#include "graph.h"
#include "heap.h"
#include <vector>
#include <map>


class Communities;
class Probabilities {
public:
  static float* tmp_vector1;	// 
  static float* tmp_vector2;	// 
  static int* id;	    // 
  static int* vertices1;    //
  static int* vertices2;    //  
  static int current_id;    // 

  static Communities* C;				    // pointer to all the communities
  static int length;					    // length of the random walks

  
  int size;						    // number of probabilities stored
  int* vertices;					    // the vertices corresponding to the stored probabilities, 0 if all the probabilities are stored
  float* P;						    // the probabilities
  
  unsigned long memory();				    // the memory (in Bytes) used by the object
  double compute_distance(const Probabilities* P2) const;   // compute the squared distance r^2 between this probability vector and P2
  Probabilities(int community);				    // compute the probability vector of a community
  Probabilities(int community1, int community2);	    // merge the probability vectors of two communities in a new one
							    // the two communities must have their probability vectors stored
							    
  ~Probabilities();					    // destructor
};

class Community {
public:
  
  Neighbor* first_neighbor;	// first item of the list of adjacent communities
  Neighbor* last_neighbor;	// last item of the list of adjacent communities
  
  int this_community;		// number of this community
  int first_member;		// number of the first vertex of the community
  int last_member;		// number of the last vertex of the community
  int size;			// number of members of the community
  
  Probabilities* P;		// the probability vector, 0 if not stored.  


  float sigma;			// sigma(C) of the community
  float internal_weight;	// sum of the weight of the internal edges
  float total_weight;		// sum of the weight of all the edges of the community (an edge between two communities is a half-edge for each community)
    
  int sub_communities[2];	// the two sub communities, -1 if no sub communities;
  int sub_community_of;		// number of the community in which this community has been merged
				// 0 if the community is active
				// -1 if the community is not used
  
  void merge(Community &C1, Community &C2);	// create a new community by merging C1 an C2
  void add_neighbor(Neighbor* N);
  void remove_neighbor(Neighbor* N);
  float min_delta_sigma();			// compute the minimal delta sigma among all the neighbors of this community
  
  Community();			// create an empty community
  ~Community();			// destructor
};

class Communities {
private:
  bool silent;		// whether the progression is displayed
  int details;		// between 0 and 3, how much details are printed
  unsigned long max_memory;	// size in Byte of maximal memory usage, -1 for no limit
  
  class Hierarchy {
    public:
    float alpha;	// the alpha_max at which the correponding partition disappears (1. = last partition)
    float s;		// = l in the paper
    float g;		// = h in the paper
    int community;	// the last community added (-1 lowest partition)
  };
  
  float* alpha_min;	    // the alpha of creation of the community
  float* alpha_max;	    // the alpha of destruction of the community
  int** sub_communities;    // the sub communities in the alpha tree
  int* sorted_alpha_min;   // the list of the community sorted by alpha_min

  void find_sub_communities(int community, bool* B, vector<int>& list);	// find the sub communities of a community that have been seen.
  void compute_hierarchy(int current_community, Hierarchy* hierarchy, int mode);	// 
  
public:
  
  unsigned long memory_used;			    // in bytes
  Min_delta_sigma_heap* min_delta_sigma;    	    // the min delta_sigma of the community with a saved probability vector (for memory management)
  
  Graph* G;		    // the graph
  int* members;		    // the members of each community represented as a chained list.
			    // a community points to the first_member the array which contains 
			    // the next member (-1 = end of the community)
  Neighbor_heap* H;	    // the distances between adjacent communities.


  Community* communities;	// array of the communities
  
  int nb_communities;		// number of valid communities 
  int nb_active_communities;	// number of active communities
  
  Communities(Graph* G, int random_walks_length = 3, bool silent = false, unsigned long max_memory = 0);    // Constructor
  ~Communities();					// Destructor


  void merge_communities(Neighbor* N);			// create a community by merging two existing communities
  double merge_nearest_communities();

  
  double compute_delta_sigma(int c1, int c2);		// compute delta_sigma(c1,c2) 

  void remove_neighbor(Neighbor* N);
  void add_neighbor(Neighbor* N);
  void update_neighbor(Neighbor* N, float new_delta_sigma);

  void manage_memory();
  
  void print_community(int c);				// print a community

  void compute_hierarchy(int mode);			// compute the hiearchical structure in alpha_min alpha_max and sub_communities.
							// mode = 1 : sigma, mode = 2 : modularity
  void print_hierarchy(int detail);				// print the whole hierarchy computed
  void print_partition(float alpha);
  void find_best_partition(float ratio, map<float,float>& M, int detail);

  void* get_hierarchy();
	  
};



#endif
