/***********************************************************************/
/** Core Percolation                                                 ***/
/**                                                                  ***/
/** C. Fay May 2002                                                  ***/
/** Version 2.0.0.0   07.12.2004                                     ***/
/** Run: core_percolation n p seed                                   ***/
/** Output: tempcore.gml                                             ***/
/**        And Data File Of: Nodes Edges Nodes_in_core Edges_in_core ***/
/***********************************************************************/

#include <cstdio>                                     
#include <cstdlib>
#include <cmath>
#include <iostream>
#include <fstream>
#include <ostream>
#include <iomanip>

#include "grphfns.hpp"
#include "crprcfns.hpp"
#include "clstrfns.hpp"

using namespace std;

int main(int argc, char *argv[]){
   int argz = 1, n=100, ncore=1, nEcore=0;
   int nE=0, nC=0, nl=0, num_cl=0, size_lc=0;
   int edge_lc=0, nT=0, seed=257, spcluster=0;  
   double p=3.0, c=0.0;
   bool write=false, dsply=false; 
   int grphtype=1;  /* 0 for rndm, 1 for pt */
   std::vector< int > gv(1,0), g2(1,0), gsize(1,0);
   char *filename;
   char prc='c';
   
   while((argv[argz][0] == '-')&&(!isdigit(argv[argz][1]))&&(argc!=0))
   {                                              
      if(strcmp(argv[argz], "-pt") == 0)
        grphtype=1;
      else if(strcmp(argv[argz], "-fcc") == 0)
        grphtype=2;
      else if(strcmp(argv[argz], "-rnd") == 0)
        grphtype=0;
      else if(strcmp(argv[argz], "-rp") == 0)
        prc = 'r';
      else if(strcmp(argv[argz], "-c") == 0)
        prc = 'c';
      else if(strcmp(argv[argz], "-t") ==0)
        prc = 't';
      else if(strcmp(argv[argz], "-ct") ==0)
        prc = 'b';
      else if(strcmp(argv[argz], "-rgml") ==0||strcmp(argv[argz], "-r") ==0){
	  argz++;
	  filename = argv[argz];
	  grphtype=-1;
          //cout<<filename<<" ";
	  }
      else{
        exit(1);
      }
      argz++; 
   }
   
   n = atoi(argv[argz++]);
   sscanf(argv[argz++], "%lf", &p);
   seed = atoi(argv[argz++]) ;   
   cout<<n<<" "<<p<<" "<<seed<<" ";
   
   switch (grphtype){
      case -1:{
         grphfns_read_gml(filename, gsize, gv, g2, n, nE);
         break;
      }
      default:{
        grphfns_build_grph_n(n, nE, p, seed, gsize, gv, g2, grphtype, write);
        break;
      }
   }

   if (dsply){
      std::ofstream output1("dsplyg.csv");
      grphfns_dsplygrph(n, gv, g2, output1);
      output1.close();
   }
   c=2.0*(double)nE/(double)n;
   cout<<c<<" "<<nE<<" ";
   
   switch (prc){
      case 'c':{           //core precolation
         ncore=0;
         crprcfns_LRA(n, nE, nC, gsize, gv, g2, nl);
         //cout<<endl<<" LRA: ";
         cout<<nl<<" "<<nC<<" ";
         if (dsply){
            std::ofstream output2("dsplyc.csv");
            grphfns_dsplygrph(n, gv, g2, output2);
            output2.close();
         }
         break;
      }
      case 't':{ /* tri_prc w/o creprc*/
         ncore=0;
         crprcfns_triprc(n, nE, nC, gsize, gv, g2, nl, nT);
         if (dsply){
            std::ofstream output4("dsplytc.csv");
            grphfns_dsplygrph(n, gv, g2, output4);
            output4.close();
	 }
         //cout<<endl<<" TRI: ";
         cout<<nl<<" "<<nT<<" "<<nC<<" ";
         break;
      }
      case 'b':{ /* tri_prc w creprc */
         ncore=0;
         crprcfns_LRA(n, nE, nC, gsize, gv, g2, nl);
         //cout<<endl<<" LRA: ";
         cout<<nl<<" "<<nC<<" ";
         if (dsply){
            std::ofstream output3("dsplyc.csv");
            grphfns_dsplygrph(n, gv, g2, output3);
            output3.close();
         }       
         clstrfns_HK(n, nE, gsize, gv, g2, num_cl, 
                 size_lc, ncore, nEcore, spcluster, edge_lc);
         c=0;
         if (ncore >0) c=2.0*(double)nEcore/(double)ncore;
	 cout<<ncore<<" "<<nEcore<<" "<<c<<" ";
         c=0;
         if (size_lc >0) c=2.0*(double)edge_lc/(double)size_lc;
         cout<<num_cl<<" "<<size_lc<<" ";
         cout<<edge_lc<<" "<<c<<" "<<spcluster<<" ";      
      
         crprcfns_triprc(n, nE, nC, gsize, gv, g2, nl, nT);
	 if (dsply){
            std::ofstream output4("dsplytc.csv");
            grphfns_dsplygrph(n, gv, g2, output4);
            output4.close();
	 }
         // cout<<endl<<" TRI: ";
         c=0;
         if (ncore>0) c=2.0*(double)nEcore/(double)ncore;
         cout<<nl<<" "<<nT<<" "<<nC<<" ";
         break;
      }
   }
   
   if (write){
      std::ofstream output("tempcre.gml");
      grphfns_write_graph(n, gv, g2, output);
      output.close();
      }
  // cout<<endl<<" HK:  "; 
  clstrfns_HK(n, nE, gsize, gv, g2, num_cl, 
          size_lc, ncore, nEcore, spcluster, edge_lc); 
  c=0;
  if (ncore >0) c=2.0*(double)nEcore/(double)ncore;
  cout<<ncore<<" "<<nEcore<<" "<<c<<" ";
  c=0;
  if (size_lc>0) c=2.0*(double)edge_lc/(double)size_lc;
  cout<<num_cl<<" "<<size_lc<<" "<<edge_lc<<" "<<c<<" "<<spcluster<<" ";
  cout<<endl;

  return 0;
}
  
