#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include "graph.hpp"
#include "cluster.hpp"

cluster_stat::cluster_stat (int n) { //generating function
   label.resize(n,-2);
   num_cluster=0;
   size_n_lc=1;
   size_edge_lc=0;
   size_core=0;
   size_core_edges=0;
   n_I=0;
   spanning_cl=-2;
   span=0;
   cl_size.resize(n,1);
}

bool cluster_stat::span_cluster(graph g){
  int i=0, clusternum=0; 
  bool possible=false, spancluster=false;
  double dn=(double) g.get_n(), sroot=sqrt(dn);
  int root = (int) (sroot+0.5);

  if (g.grphtype==1){//if a 2-D graph
	 dn=(double) g.get_n();
	 sroot=sqrt(dn);
	 root = (int) (sroot+0.5);
  }
  if (g.grphtype==2){//if a 3-D graph volume = root X root X root,  Area of size = root X root
	 sroot = std::pow((double)g.get_n(), 1/3.);
	 int croot = (int) (sroot+0.5); 
	 root = croot*croot;  
  }
  
  
  cout<<" root "<<root<<" "<<l_cluster_label;
  if(size_n_lc < root){
    spancluster=false;
   /* cout<<"size_lc < root"<<endl;*/
    }
  else if(size_n_lc  > g.get_n()-root){
    spancluster=true;
    cout<<"size_lc > n-root"<<endl;
    }
  else{ 
     i=0;    
     while((!possible) && i<=root-1){ //see if node i is on largest cluster
      	if ((label[i]==l_cluster_label) || (HK_find_simple(i)==l_cluster_label)) possible=true; //test node on bottom row
	    i++; 
     }
     /*cout<<" possible "<<possible<<"__";*/
     i=g.get_n()-root-1;
     //cout<<" "<<i<<" "<<root;
     while(possible && !spancluster && i<g.get_n()){ //test node on top row
        cout<<"poss "<<possible<<" span "<<spancluster<<" node "<<i<<" "<<label[i]<<" "<<l_cluster_label<<endl;
	    if ((label[i]==l_cluster_label) || (HK_find_simple(i)==l_cluster_label)) spancluster=true;
	    cout<<"poss "<<possible<<" span "<<spancluster<<" node "<<i<<endl;
	    i++;   
     }   
  }	
   return (spancluster);	
}
void cluster_stat::write_clstr_to_screen(graph g){
   for (int t=0; t<g.nodes.size(); t++) {
     std::cout<<"  node [ id "<<t<<" ]"<<" cluster [ "<<label[t]<<" ]"<<endl;
   }
   std::cout<<"\n";
   return;		
}

int cluster_stat::HK_find_simple(int node1){
	int test_node=node1, temp_label;
	while ((label[test_node]!=test_node) && (label[test_node]!=-1)){ //find cannonical for node1
		test_node=label[test_node];
	}
	//the final temp_label is cannonical
	return test_node;	
}

int cluster_stat::HK_find(int node1){
	int test_node=node1, temp_label;
	while ((label[test_node]!=test_node) && (label[test_node]!=-1)){ //find cannonical for node1
		test_node=label[test_node];
	}
	//the final temp_label is cannonical
	while (label[node1]!=node1){ //relabel tree to cannonical
	   	int step_node = label[node1];
	   	label[node1]=test_node;
	   	node1=step_node;
	}
	return test_node;	
}

void cluster_stat::HK_join(int node1, int node2, int node2_c, int node1_c){
	//put node 2 into cluster of node1
	label[node1]=node1_c;
	label[node2]=node1_c;
	label[node2_c]=node1_c;
	cl_size[node1_c]=cl_size[node1_c]+cl_size[node2_c];
	cl_e_size[node1_c]=cl_e_size[node1_c]+cl_e_size[node2_c];
	if (size_n_lc<cl_size[node1_c]){ //count the nodes in the largest cluster
		size_n_lc=cl_size[node1_c];
		l_cluster_label=node1_c;
	}
	//cluster[node2]=node1_c;	
	return;
}
void cluster_stat::HK_union(int node1, int cannonical_tree){
	label[HK_find(node1)]=HK_find(cannonical_tree);
	cl_size[cannonical_tree]=cl_size[cannonical_tree]+cl_size[HK_find(node1)];
	if (size_n_lc<cl_size[cannonical_tree]){ //count the nodes in the largest cluster
		size_n_lc=cl_size[cannonical_tree];
		l_cluster_label=cannonical_tree;
	}
	return ;
}

void cluster_stat::HK(graph g){
   int n1_cannonical=0, n2_cannonical=1, n1=0, n2=1, j=0;
   
   size_n_lc=1;
   
   for(unsigned int n1=0; n1<g.nodes.size(); n1++){  //label each node as its own cluster
     label[n1]=n1; // Later: I want empty clusters to be labeled -1
     cl_size[n1]=1;
     cl_e_size.push_back(g.nodes[n1].edges.size());
   }
   num_cluster=g.nodes.size();
   write_clstr_to_screen(g);
   
   for(unsigned int n1=0; n1<g.nodes.size(); n1++){ //for every node n1
	   //cout<<"node:"<<n1<<" size "<<g.nodes[n1].edges.size();
	  if (g.nodes[n1].edges.size()<=0) {
         label[n1]=-1;             /* label nodes of size zero to cluster -1 */
	     n_I++;       /* n_I = number of Isolated points */ //clusters of size 1 (1 node)
	     num_cluster--;
	  }
	  else {  //node n1 has edges
		 if (label[n1]==-2) { //unlabeled
		 	 label[n1]=n1;
		     num_cluster++;
		 }
         for(j=0; j<g.nodes[n1].edges.size(); j++){ /* for each node n2 connected to n1*/
			 n1_cannonical = HK_find(n1);   //find cannonical label n1
             n2=g.nodes[n1].edges[j]; //find node n2  //we only need to check the clusters of nodes n2<n1
               //if (n2<n1) {//label n2 with the label of n1
              // cout<<"edge "<<n1<<" "<<n2<<" "<<label[n2];
              
             if (label[n2]==-2) {//n2 not in a cluster put in cluster of n1
				   label[n2]=n1_cannonical;
				   cl_size[n1_cannonical]++;
				   } 
			 n2_cannonical=HK_find(n2); 
			   //cout<<" "<<n1_cannonical<<" "<<n2_cannonical<<" labels "<<label[n1]<<" "<<label[n2]<<endl;
			 if ((n1_cannonical<n2_cannonical) && (n1_cannonical!=n2_cannonical)){ //put 2 into 1
			    HK_join(n1,n2,n2_cannonical, n1_cannonical);
			   // HK_union(n2,n1);
				num_cluster--;
				 // cl_size[n1]=cl_size[n1]+cl_size[n2];
				 // cout<<"n2 add "<<n2<<" in cluster "<<n2_cannonical<<" to cluster of "<<n1<<" in cluster "<<n1_cannonical<<endl;
			 }
			 else if ((n1_cannonical>n2_cannonical) && (n1_cannonical!=n2_cannonical)) { //put 1 into 2
			    HK_join(n2,n1,n1_cannonical, n2_cannonical); 
			   // HK_union(n1,n2);
				num_cluster--;
				 // cl_size[n2]=cl_size[n2]+cl_size[n1];
				 // cout<<"n1 add "<<n1<<" in cluster "<<n1_cannonical<<" to cluster of "<<n2<<" in cluster ";
				 // cout<<n2_cannonical<<" "<<label[n1]<<endl;
			 }
		   
		}
      
      }   
      
   }
   
   size_n_lc=count_clstr(g);
   return;
}


void cluster_stat::HK2(graph g){
   int n1_cannonical=0, n2_cannonical=1, n1=0, n2=1, j=0;
   
//for(unsigned int i=0; i<g.n; i++){ //label each node as its own cluster
  //   label[i]=-2; // Later: I want empty clusters to be labeled -1
  // }
  // num_cluster=0;
  // cluster[g.n+1]=0;
   //the vector cluster counts the canonical clusters and points the non-cannonical to the good clusters
   write_clstr_to_screen(g);
   
   for(unsigned int n1=0; n1<g.nodes.size(); n1++){ //for every node n1
	   //cout<<"node:"<<n1<<" size "<<g.nodes[n1].edges.size();
	  if (g.nodes[n1].edges.size()<=0) {
         label[n1]=-1;             /* label nodes of size zero to cluster -1 */
	     n_I++;       /* n_I = number of Isolated points */ //clusters of size 1 (1 node)
	     //cl_size[n1]=1;
	     //cluster[g.n+1]++;     //in the cluster vector element (n+1) is the number of isolated clusters. 
	    // cout<<" isolated node"<<n1<<" "<<n_I;
	  }
	  else {  //node n1 has edges
		 if (label[n1]==-2) { //unlabeled
		 	 label[n1]=n1;
		     num_cluster++;
		     //cl_size[n1]++;
		  //   cout<<" new cluster "<<n1<<" "<<num_cluster<<endl;;
		 }
		 n1_cannonical = HK_find(n1);  //find cannonical label n1
            for(j=0; j<g.nodes[n1].edges.size(); j++){ /* for each node n2 connected to n1*/ 
               n2=g.nodes[n1].edges[j]; //find node n2  //we only need to check the clusters of nodes n2<n1
               //if (n2<n1) {//label n2 with the label of n1
              // cout<<"edge "<<n1<<" "<<n2<<" "<<label[n2];
               if (label[n2]==-2) {//n2 not in a cluster put in cluster of n1
				   label[n2]=n1_cannonical;
				//   cout<<" "<<label[n2];
				  // cl_size[n1]++;
				   } 
			   n2_cannonical=HK_find(n2); 
			   //cout<<" "<<n1_cannonical<<" "<<n2_cannonical<<" labels "<<label[n1]<<" "<<label[n2]<<endl;
			   if (n1_cannonical<n2_cannonical){ //put 2 into 1
				  HK_join(n1,n2,n1_cannonical, n2_cannonical);
				  num_cluster--;
				 // cl_size[n1]=cl_size[n1]+cl_size[n2];
				//  cout<<" add"<<n2<<" to cluster of "<<n1<<endl;
			   }
			   else if (n1_cannonical>n2_cannonical) { //put 1 into 2
				  HK_join(n2,n1,n2_cannonical, n1_cannonical); 
				  num_cluster--;
				 // cl_size[n2]=cl_size[n2]+cl_size[n1];
				 // cout<<" add"<<n1<<" to cluster of "<<n2<<endl;
			   }
		    }
      }   
   }
   return;
}
int cluster_stat::count_clstr(graph g){
   int max_size=1, l_cluster_label=-1, max_edge;
   //std::vector<int>cl_e_size;
   for (int t=0; t<g.nodes.size(); t++) {
	   if (label[t]!=t && label[t]!=-1){
	       HK_find(t);
	   }  
	   if (max_size<cl_size[t]) {
		   max_size=cl_size[t]; 
		   max_edge=cl_e_size[t];
		   l_cluster_label=t;
	   } 
   }
   
   size_edge_lc=max_edge/2;   
   return max_size;
}
