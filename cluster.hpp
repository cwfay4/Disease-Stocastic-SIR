#ifndef CLUSTER_H
#define CLUSTER_H

#include "graph.hpp"

class cluster_stat{
    public:
       int num_cluster;          //num_cl
       int size_n_lc;        //size_lc
       int size_edge_lc;  //edge_lc
       int size_core;             //ncore
       int size_core_edges;       //nEcore
       int spanning_cl;           //spcluster  
       int n_I;
       bool span;
       int l_cluster_label;
       cluster_stat (int);
       int count_clstr(graph); //returns size of largest cluster
       bool span_cluster(graph); //returns true of largest cluster spans
       void HK(graph g);  //Hoshen-Kopleman works
       void HK2(graph g);//Hoshen-Kopleman not so much
       void write_clstr_to_screen(graph g);
    private:
       std::vector<int> label;
       std::vector<int> cl_size;   //number of nodes in clusters;
       std::vector<int> cl_e_size; //number of edges in clusters;
       int HK_find_simple(int);
       int HK_find(int);
       void HK_join(int, int, int, int);
       void HK_union(int, int);       
};


#endif
