/***********************************************************************/
/** make graph                                                       ***/
/**     a program to make graphs                                     ***/
/**                                                                  ***/
/**                                                                  ***/
/**                                                                  ***/
/**                                                                  ***/
/** C. Fay Jan 2015                                                  ***/
/** Version 0.0.0.1   05.10.2004                                     ***/
/** Run: mkgrph number_nodes connectivity seed                       ***/
/** Output: dat files that can be opened by xmgrace as a visual      ***/
/**        representation of the graph                               ***/
/***********************************************************************/

/***********************************************************************/
/** August 26 2016                                                   ***/
/**   read and write gml appears to work                             ***/
/**   g.get_nE() returns the number of edges x2                      ***/
/**   cluster counting not checked                                   ***/
/***********************************************************************/


#include <cstdio>                                     
#include <cstdlib>
#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <iomanip>
#include <vector>
#include <time.h>

//#include "Table2D.hpp"
#include "graph.hpp"
#include "cluster.hpp"
#include "random.hpp"
//#include "HK.hpp"

using namespace std;


/***********************************************************************/
/**  ran0 from numerical recipes                                     ***/
/**                                                                  ***/
/***********************************************************************/
//#define IA 16807
//#define IM 2147483647
//#define AM (1.0/IM)
//#define IQ 127773
//#define IR 2836
//#define MASK 123459876
/***********************************************************************/
//float ran0(long *idum)
//{
//     long k;
//     float ans;
//          
//     *idum ^= MASK;
//     k=(*idum)/IQ;
//     *idum=IA*(*idum-k*IQ)-IR*k;
//     if (*idum <0) *idum += IM;
//     ans=AM*(*idum);
//     *idum ^= MASK;
//     return (ans);
//}
/***********************************************************************/
/***********************************************************************/

/************* build_ptg_n2() ***************************/
/** builds planar triangular lattice with n nodes      **/
/** written by c.fay                                   **/
/**                                                    **/
/** Returns:                                           **/
/**     nothing                                        **/
/********************************************************/

void GenEdge(int n, long &idum, int &r1, int &r2){
   int l1,l2;
   double dl;
   
   dl=fmod(ran0(&idum)*n,n)+1;
   l1=(int) dl;
   dl=fmod(ran0(&idum)*n,n)+1;
   l2=(int) dl;
   while(l1==l2){
     dl=fmod(ran0(&idum)*n,n)+1;
     l2=(int)dl;
   }   
   if(l1>l2){
     r2=l1;
     r1=l2;
   }
   else{
     r2=l2;
     r1=l1;
   }
  /* cout<<r1<<" "<<r2<<endl;*/
   return;
}

int main(int argc, char *argv[]){
   int argz = 1, num_nodes=100, nE=0, num_edge=0, seed=2579831, grphtype=1; /* 0 for rndm, 1 for pt */
   char prc='c';
   double edge_prob=3.0, c=0.0, test_p;
   bool write=true, debug_test=true;  
   graph g;
   seed = (int)( time(NULL) );
   int dummy;

   std::stringstream ss;
   std::string inputstr, inputfile;
   
   g.set_debug(debug_test);
   
   int i, j=0;
   for (i = 1; i < argc; i++) {
	  ss.clear();
	  ss.str(argv[i]);
      ss>>inputstr; 
/*********************************************************************/
      if(inputstr.compare("-pt") == 0) g.grphtype=1;         //planar triangular graph
      else if(inputstr.compare("-ptr3") == 0) g.grphtype=2;   //planar triangular graph with ran3
      else if(inputstr.compare("-fcc") == 0) g.grphtype=3;   //FCC Lattice
      else if(inputstr.compare("-cubic") == 0) g.grphtype=4;   //cubic Lattice
      else if(inputstr.compare("-sq") == 0) g.grphtype=5;   //squareLattice     
      else if(inputstr.compare("-sclf") == 0) g.grphtype=6;   //scale free Lattice
      else if(inputstr.compare("-staticsclf") == 0) g.grphtype=7;   //scale free Lattice      
      else if(inputstr.compare("-bb") == 0) g.grphtype=8;   //scale free Lattice          
//      else if(inputstr.compare("-bb") == 0) {                 //barabasi graph
//     	g.grphtype=8;
//	    i++;
//		ss>>dummy;
//		dummy=atoi(argv[i]);
//	    g.set_L(dummy);
	   // g.set_L(atoi(argv[argz]));
	    //cout<<"nkernal: "<<L<<" ";
//	  }
      else if(inputstr.compare("-rndfixedc") == 0||inputstr.compare("-rf") ==0){          //random fixed connectivity
      	g.grphtype=9;
	    i++;
		ss>>dummy;
		dummy=atoi(argv[i]);
	    g.set_L(dummy);
	    cout<<"fixed maxc: "<<g.get_L()<<" ";
	  }
      else if(inputstr.compare("-rnd2") == 0) g.grphtype=10;   //random Lattice type 2	  
      else if(inputstr.compare("-rnd") == 0) g.grphtype=0;   //random graph   
      else if(inputstr.compare("-rgml") ==0||inputstr.compare("-r") ==0){  //read a file in gml format
	     i++;
	     ss.clear();
	     ss.str(argv[i]);
	     ss>>inputfile;
	     g.grphtype=-1;
         if (debug_test) cout<<"You want to read the file "<<inputfile<<std::endl;
	  }
/*********************************************************************/	  
	  else if(j==0){
		  ss>>num_nodes;
		  num_nodes=atoi(argv[i]);
		  if (g.grphtype==1 || g.grphtype==2 || g.grphtype==3 || g.grphtype==4 || g.grphtype==5){
		     g.set_L(num_nodes);
		     if(g.grphtype==1 || g.grphtype==2 || g.grphtype==5){
				 num_nodes = pow(g.get_L(),2);
				  g.set_n(num_nodes);
		     }
		     else {
                num_nodes = pow(g.get_L(),3);
				g.set_n(num_nodes); 
		     }
	      }
	      else g.set_n(num_nodes);
		  j++;
		  cout<<num_nodes<<endl;
	  }
	  else if(j==1){
		  edge_prob=atof(argv[i]); //works
		  g.p=edge_prob;
		  cout<<"p="<<ss.str()<<" "<<edge_prob<<std::endl;
		  j++;
	  }
	  else if(j==2){
		  ss>>seed;
		  seed=atoi(argv[i]);
		  cout<<"seed="<<ss.str()<<" "<<seed<<std::endl;
		  j++;
	  }	  
//	  cout<<argv[i]<<" "<<ss<<" "<<inputstring<<endl;
   }
   
//   n = atoi(argv[argz++]);
  // sscanf(argv[argz++], "%lf", &p);
 //  seed = atoi(argv[argz++]) ;   
     if (debug_test) cout<<"Make graph: grphtype "<<g.grphtype<<" "<<prc<<" "<<num_nodes<<" "<<edge_prob<<" "<<seed<<std::endl; 
  //  if (debug_test) cout<<"Make graph "<<g.grphtype<<" "<<prc<<" "<<g.get_n()<<" "<<g.p<<" "<<seed<<std::endl; 
     if (debug_test) g.set_debug(true);    
     if (g.grphtype==-1) g.read_gml(inputfile);
     else  {
		 g.seed=seed;
		 g.build_graph(write);
	 }
 
   c=g.get_c();
   if (debug_test) {
	     cout<<"made graph "<<g.get_n()<<" nE="<<g.get_nE()<<" n="<<g.nodes.size()<<" c="<<g.c;
         cout<<std::endl;
	 }
   
  /* cout<<gv.size()<<" "<<g2.size()<<" "<<gsize.size();*/
   
   
//   if (write){
        std::string outputstr="tempcre.gml";
        g.write_gml(outputstr);
//     std::ofstream output("tempcre.gml");
 //    grphfns_write_graph(n, gv, g2, output);
//     output.close();
//      }
  // std::ofstream output2("displ.csv");
//   grphfns_dsplygrph2(n, gv, g2);
   //output2.close();
   cluster_stat clstr (g.nodes.size());
   if (debug_test) {
	    cout<<"clustr stat built"<<endl;
        cout<<g.get_n()<<" "<<g.get_nE()<<" "<<g.nodes.size()<<" "<<c<<" "<<g.get_nE()<<" "<<g.c;
        cout<<std::endl;
	}
   
   // cout<<g.n<<" "<<c<<" "<<g.nE<<" ";
   //cout<<std::endl;  
   
//   clstr.HK(g);
      
   cout<<"prepare to write cluster data"<<endl;
//   clstr.write_clstr_to_screen(g);
   
//   bool span=false;
//   cout<<"fspan="<<span<<endl;
//   span=clstr.span_cluster(g);
   
//   cout<<"found span cluster "<<span<<endl;
//   cout<<"CL stat:  num clustr "<<clstr.num_cluster<<" Isolated nodes "<<clstr.n_I;
//   cout<<" size_l_cluster "<<clstr.count_clstr(g)<<" "<<clstr.size_n_lc<<" "<<clstr.size_edge_lc;
//   cout<<" span? "<<span;
   //cout<<" size_l_cluster_edges "<<clstr.size_l_cluster_edges<<std::endl;
   cout<<endl;
   
   return 0;
}
