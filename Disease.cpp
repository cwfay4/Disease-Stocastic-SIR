/***********************************************************************/
/** LoPR Algorithm Control Program          _array                   ***/
/**                                                                  ***/
/** C. Fay August 2005                                               ***/
/** Version 2.1.0.1   22.09.2005                                     ***/
/** Run: vertcvr gml_filname seed                                    ***/
/** Output:                                                          ***/
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

#include "random.hpp"
#include "Stats.hpp"
#include "graph.hpp"
#include "diseasePOP.hpp"
#include "Disease_basic_utils.hpp"

//#include "defaults.hpp"
//#include "grphfns_a.hpp"
//#include "cvrfns_a.hpp"
//#include "Table2D.hpp"
//#include "LoPR_alg_o.hpp"
//#include "LoPR_basic_utils.hpp"

using namespace std;

#define NUM_STEPS 1000
//#define grphtypedflt 0
//#define debug true


int main(int argc, char *argv[]){
   int argz = 1;
   
//   double Psum=0;
  // double p=3.0;
//   bool write=true, dsply=false, histo=false, debug_test=false, fulloutput=false, conv=true;
//   bool style=false; //is a new graph made for each pop.iteration or not
   boolcat b;
//   b.NUM_ITER=10;
//   char *filename;
//   char LoPR;
  // int grphtype=0;  /* 0 for rndm, 1 for pt */
   long lseed[3]={0},lseedstart=-1, lseedblnk=lseedstart;
   dPOP pop(0,0,0,0,0.1,0.1,0.1,0.1,10,false,false);
 //  pop.iteration=10;	pop.lifetime=10;    pop.I=false;    pop.simple=false;
//   int dummy;
   int ITERCASE=-1;
   double ITERSTEP, ITERSTART, ITERSTOP;
   std::string dname;

   graph g;
   g.grphtype=grphtypedflt;
   
   inputprocessing(argc, argv,  g, pop, b,  ITERCASE); //this reads the command line for all of the relevant controls and values
   
   cout<<"n= "<<g.get_n()<<" c= "<<g.get_c()<<"["<<g.p<<"]"<<" nE= "<<g.get_nE()<<" p_sick= "<<pop.p_sick;
   cout<<" contagin= "<<pop.contagin<<" fatality="<<pop.fatality<<endl;
   cout<<"p_Immune="<<pop.p_Immune<<" lifetime= "<<pop.lifetime<<" p_Recovery= "<<pop.p_Recovery<<endl;
   
//****** Set output file name ************************************
   //std::ofstream dout("dout.csv");
   std::string oname=outputname(g,pop);
   //std::string oname="disease-ss.csv";
   if (b.debug_test) cout<<oname<<std::endl;
//**************************************************************** 
    
   if (b.debug_test) cout<<lseedblnk<<" "<<lseedstart<<std::endl;       
   if (lseedstart==-1) lseedstart=time(0);
   lseedblnk=lseedstart;
   if (b.debug_test) cout<<lseedblnk<<" "<<lseedstart<<std::endl;   
   for (unsigned int i=1;i<=1000;i++){
        ran2(&lseedblnk);
        }
   for (unsigned int i=0;i<3;i++){  
      lseed[i]=(-1*(long) (ran2(&lseedblnk)*1e8));
   }
   if (b.debug_test) cout<<lseedblnk<<" "<<lseedstart<<" ran2 initialized "<<std::endl;
        
   cout<<g.get_n()<<" "<<g.p<<" "<<lseedstart<<" "; 
//****************************************************************    
   
   //long seed=lseed[0];  
   if (b.debug_test) cout<<" preparing to make a graph "<<endl;  
   g.set_debug(b.debug_test);
   switch (g.grphtype){
      case -1:{
		 if (b.debug_test) cout<<dname<<std::endl; 
         g.read_gml(dname);
         break;
      }
      default:{
		 g.seed=lseed[0];
//		 if (debug_test) cout<<" build graph with "<<lseed[0]<<" "<<g.grphtype<<std::endl;
		 g.build_graph(b.write);
		 if (b.debug_test) cout<<" ok a graph is built "<<lseed[0]<<std::endl;
         break;
      }
   }

   double c=g.get_c();
      if (b.debug_test) cout<<g.c<<std::endl;
   cout<<"n= "<<g.get_n()<<" c= "<<g.get_c()<<" nE= "<<g.get_nE()<<" p_sick= "<<pop.p_sick;
   cout<<" contagin= "<<pop.contagin<<" fatality="<<pop.fatality<<endl;
   cout<<"p_Immune="<<pop.p_Immune<<" p_Recovery= "<<pop.p_Recovery<<endl;

//      cout<<n<<" "<<p<<" "<<c<<" "<<nE<<" "<<lseed[0]<<" "<<lseed[1]<<" "<<lseed[2];        
//   }
//   else cout<<c<<" "<<nE<<" "<<" ";
   //vce_precon(g, lseed[1], pop);
   pop.precondition(g, lseed[1]);
   if (b.debug_test) cout<<"preconditioned finished"<<std::endl;  
    
  // disease_pop(g, pop, lseed[2]);
    switch (pop.ltype){
        case 1:{ //only a lifetime
        	pop.p_Recovery=1.1;
	   	    pop.pop_evolve_lifetime(g, lseed[2]);
		    break;
	    }
	    case 2:{ //lifetime and p_Recovery
	   	    pop.pop_evolve(g, lseed[2]);	  	
		    break;
	    }
   	    default:{ //only p_Recovery
   	    	pop.lifetime=0;
   	     	pop.pop_evolve(g, lseed[2]); 
			break;
	    }   	
   }
    
   if (b.debug_test) cout<<"disease pop finished"<<std::endl;  
    
   cout<<"Max Sick="<<pop.n_sick_max<<endl;
   cout<<"Time of Max Sick="<<pop.n_sick_mIter<<endl;
   cout<<"n= "<<g.get_n()<<" c= "<<g.get_c()<<"["<<g.p<<"]"<<" nE= "<<g.get_nE()<<" p_sick= "<<pop.p_sick<<" p_sick_max= "<<pop.n_sick_max;
   cout<<" contagin= "<<pop.contagin<<endl;
   cout<<"fatality="<<pop.fatality<<" p_Immune= "<<pop.p_Immune<<" lifetime= "<<pop.lifetime<<" p_Recovery= "<<pop.p_Recovery<<endl;
   //cout<<g.grphtype<<endl;
   return 0;
}
