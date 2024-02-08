/***********************************************************************/
/** Hoshen-Kopelman Cluster counting                                 ***/
/**                                                                  ***/
/** C. Fay Jan 2003                                                  ***/
/** Version 2.0.0.0   05.10.2004                                     ***/
/** Hoshen-Kopleman Algorithm for counting cluster sizes             ***/
/***********************************************************************/

#include <cstdio>                                     
#include <cstdlib>
#include <cmath>
#include <iostream>
#include <fstream>
#include <ostream>
#include <iomanip>
#include <vector>
#include "graph.hpp"
#include "cluster.hpp"

int HK_find(int node1, std::vector<int> label){
	int test_node=node1, temp_label;
	while (label[test_node]!=test_node){ //find cannonical for node1
		test_node=label[test_node];
	}
	//the final temp_label is cannonical
	while (label[node1]!=node1){ //relabel tree to cannonical
	   	int step_node = label[node1];
	   	label[node1]=temp_label;
	   	node1=step_node;
	}
	return test_node;	
}
void HK_join(int node1, int node2, int node1_c, int node2_c, std::vector<int> label){
	//put node 2 into cluster of node1
	label[node1]=node1_c;
	label[node2]=node1_c;
	//cluster[node2]=node1_c;	
}
int HK_union(int node1, int node2, std::vector<int> label){
	return label[HK_find(node1,label)]=HK_find(node2,label);
}

void HK(graph g, cluster_stat cluster_data){
   int size_n=g.get_n();
   std::vector<int> label(g.get_n()+1,-2), cluster(g.get_n()+1,-1);
   int n1_cannonical, n2_cannonical, n1,n2, j=0;
   
   for(unsigned int i=0; i<size_n; i++){ //label each node as its own cluster
      cluster[i]=i; // Later: I want empty clusters to be labeled -1
      g.nodes[i].cluster=i;
   }
   cluster_data.num_cluster=0;
   cluster[size_n+1]=0;
   //the vector cluster counts the canonical clusters and points the non-cannonical to the good clusters
   for(unsigned int n1=0; n1<size_n; n1++){ //for every node n1
	  if (g.nodes[n1].edges.size()<=0) {
         label[n1]=-1;             /* label nodes of size zero to cluster -1 */
	     cluster_data.n_I++;       /* n_I = number of Isolated points */ //clusters of size 1 (1 node)
	     cluster[size_n+1]++;     //in the cluster vector element (n+1) is the number of isolated clusters. 
	  }
	  else {  //node n1 has edges
		 if (label[n1]==-2) { //unlabeled
			 label[n1]=n1;
		     cluster_data.num_cluster++;
		 }
		 n1_cannonical = HK_find(n1,label);  //find cannonical label n1
            for(j=0; j<g.nodes[n1].edges.size(); j++){ /* for each node n2 connected to n1*/ 
               n2=g.nodes[n1].edges[j]; //find node n2  //we only need to check the clusters of nodes n2<n1
               if (n2<n1) {//label n2 with the label of n1
                  if (label[n2]==-2){
					 label[n2]=n1_cannonical; //n2 not in a cluster put in cluster of n1
				     }
			      n2_cannonical=HK_find(n2, label); //find n2_cannonical  
			  // }
			      if (n1_cannonical<n2_cannonical){ //put 2 into 1
				     HK_join(n1,n2,n1_cannonical, n2_cannonical, label);
				     cluster_data.num_cluster--;
			      }
			      else if (n1_cannonical>n2_cannonical) { //put 1 into 2
				      HK_join(n2,n1,n2_cannonical, n1_cannonical, label); 
				      cluster_data.num_cluster--;
			      }
			   }
		    }
      }   
   }
   return;
}
