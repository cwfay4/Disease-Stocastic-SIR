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

   double Psum=0;
   double p=3.0;
   bool write=true, dsply=false, histo=false, debug_test=false, fulloutput=false, conv=true;
   char *filename;
   char LoPR;
   int grphtype=0;  /* 0 for rndm, 1 for pt */
   long lseed[3]={0},lseedstart=-1, lseedblnk=lseedstart;
   dPOP pop;
   pop.lifetime=10;
   int dummy;
   pop.simple=false;

   graph g;
   g.grphtype=grphtypedflt;
   
   std::stringstream ss;
   std::string inputstr, inputfile;
   
   int i, j=0;
   for (i = 1; i < argc; i++) {
   	  ss.clear();
	  ss.str(argv[i]);
      ss>>inputstr;
      if (inputstr.compare("-debug") == 0) debug_test=true;
      else if(inputstr.compare("-simple") == 0) pop.simple=true;
      else if(inputstr.compare("-pt") == 0) g.grphtype=1;         //planar triangular graph
      else if(inputstr.compare("-ptr3") == 0) g.grphtype=2;   //planar triangular graph with ran3
      else if(inputstr.compare("-fcc") == 0) g.grphtype=3;   //FCC Lattice
      else if(inputstr.compare("-cubic") == 0) g.grphtype=4;   //cubic Lattice
      else if(inputstr.compare("-sq") == 0) g.grphtype=5;   //squareLattice     
      else if(inputstr.compare("-sclf") == 0) g.grphtype=6;   //scale free Lattice
      else if(inputstr.compare("-rnd2") == 0) g.grphtype=7;   //random Lattice type 2
      else if(inputstr.compare("-bb") == 0) g.grphtype=8;   //scale free Lattice          
//      else if(inputstr.compare("-bb") == 0) {                 //bethe lattice?
 //     	g.grphtype=8;
//	    i++;
//		ss>>dummy;
//		dummy=atoi(argv[i]);
//	    g.set_L(dummy);
	   // g.set_L(atoi(argv[argz]));
	    //cout<<"nkernal: "<<L<<" ";
//	  }
      else if(inputstr.compare("-rndfixedc") == 0||inputstr.compare("-rf") ==0){ //random fixed connectivity
      	g.grphtype=9;
	    i++;
		ss>>dummy;
		dummy=atoi(argv[i]);
	    g.set_L(dummy);
	    if (debug_test) cout<<"fixed maxc: "<<g.get_L()<<" ";
	  }
      else if(inputstr.compare("-rnd") == 0)        g.grphtype=0;   //random graph  
      else if(inputstr.compare("-I") == 0)          pop.I=true;   //Immune after sickness     
      else if(inputstr.compare("-w") == 0)       	write=true;
      else if(inputstr.compare("-dsply") == 0)     	dsply=true;
      else if(inputstr.compare("-fullout") == 0)   	fulloutput=true;
      else if(inputstr.compare("-histo") == 0)     	histo=true; 
      else if(inputstr.compare("-conv") == 0)     	conv=true; 
            else if(inputstr.compare("-h") == 0){  //help menu
		Disease_usage(false); 
		return 0;  
	  }     	
      else if(inputstr.compare("-hv") == 0){  //help menu verbose
		Disease_usage(true); 
		return 0;  
	  }	       
      else if(inputstr.compare("-rgml") ==0||inputstr.compare("-r") ==0){  //read a file in gml format
	     i++;
	     ss.clear();
	     ss.str(argv[i]);
	     ss>>inputfile;
	     g.grphtype=-1;
         if (debug_test) cout<<"You want to read the file "<<inputfile<<std::endl;
	  }
      else if(j==0){  //input # of nodes
		  ss>>dummy;
		  dummy=atoi(argv[i]);
		  if (g.grphtype==1 || g.grphtype==2 || g.grphtype==3 || g.grphtype==4 || g.grphtype==5){
		     g.set_L(dummy);
		     if(g.grphtype==1 || g.grphtype==2 || g.grphtype==5){
				 dummy = pow(g.get_L(),2);
				  g.set_n(dummy);
		     }
		     else {
                dummy = pow(g.get_L(),3);
				g.set_n(dummy); 
		     }
	      }
	      else g.set_n(dummy);
		  j++;
		  if (debug_test) cout<<"n= "<<g.get_n()<<endl;
		 // i++;
	  }
	  else if(j==1){  //input edge probability
		  g.p=atof(argv[i]); //works
		  if (true) cout<<"p="<<ss.str()<<" "<<g.p<<std::endl;
		  j++;
	  }
	  else if(j==2){ //probability of being initially sick
		  pop.p_sick=atof(argv[i]); //works
		  if (debug_test) cout<<"p_sick="<<ss.str()<<" "<<pop.p_sick<<std::endl;
		  j++;
	  }
	  else if(j==3){ //probability of getting sick
		  pop.contagin=atof(argv[i]); //works
		  if (debug_test) cout<<"p_contagin="<<ss.str()<<" "<<pop.contagin<<std::endl;
		  j++;
	  }
	  else if(j==4){ //probability of death!
		  pop.fatality=atof(argv[i]); //works
		  if (debug_test) cout<<"fatality="<<ss.str()<<" "<<pop.fatality<<std::endl;
		  j++;
	  }	  
	  else if(j==5){ //probability of initial immunity
		  pop.p_Immune=atof(argv[i]); //works
		  if (debug_test) cout<<"p_Immune="<<ss.str()<<" "<<pop.p_Immune<<std::endl;
		  j++;
	  }
	  else if(j==6){ //sickness duration or probability of recovery
		  pop.p_Recovery=atof(argv[i]); //added a new parameter, probability of recovery
		  pop.lifetime=0; //works, initializing lifetime to be zero but in a way so it can be changed later
		  if (debug_test) cout<<"p_Recovery="<<ss.str()<<" "<<pop.p_Recovery<<std::endl;
		  j++;
	  }	  
//	  if (lseedstart==-1) lseedstart=time(0);   //starting seed
//	  else {
//		  ss>>lseedstart;
//		  lseedstart=atoi(argv[i]);
//		  j++;
//	  }		  
//	  if (debug_test) cout<<"seed="<<ss.str()<<" "<<lseedstart<<std::endl;	  	  
//	  if (g.grphtype==-1) nsamp=1;  //number of samples to run
//      else nsamp = atoi(argv[argz++]);
   
//      if (g.grphtype==8) cout<<"nkernal: "<<g.get_L()<<" "; // dimension of lattice
      if (g.grphtype==3 || g.grphtype==4) { // dimension of lattice
         g.set_L(g.get_n());
         g.set_n(std::pow(g.get_L(),3));
      }   
 //     if (gn > SIZE) { depreciated
 //        cout<<"error n to large for compiled code.  n:"<<n<<" max_n:"<<SIZE<<endl;
 //        exit(1);   
 //     }    
	  
//	  cout<<argv[i]<<" "<<ss<<" "<<inputstring<<endl;
   }
   
   cout<<"n= "<<g.get_n()<<" c= "<<g.get_c()<<"["<<g.p<<"]"<<" nE= "<<g.get_nE()<<" p_sick= "<<pop.p_sick;
   cout<<" contagin= "<<pop.contagin<<" fatality="<<pop.fatality<<endl;
   cout<<"p_Immune="<<pop.p_Immune<<" lifetime= "<<pop.lifetime<<"p_Recovery= "<<pop.p_Recovery<<endl;
//*****************************************************************        
   std::string gt="p", rnt="-r2", sbc="f";
   std::ostringstream strs;
   strs << g.get_n();
   std::string nstr = strs.str();
   strs << g.p;
   std::string pstr = strs.str();  
   if (g.bc==1) sbc="p";
   else if (g.bc==2) sbc="pp";
   if (g.grphtype == 0 || g.grphtype==7) gt="r";
   else if (g.grphtype == 2) rnt="-r3";
   else if (g.grphtype == 3) gt="f";
   else if (g.grphtype == 4) gt="c";
   else if (g.grphtype == 5) gt="sq";
   else if (g.grphtype == 6) gt="sclf";
   else if (g.grphtype == 8) gt="bb";
   
  // std::ofstream dout("dout.csv");
   std::string dname="dLo-"+gt+nstr+"-"+pstr+"-"+sbc+rnt+".csv";
   if (debug_test) cout<<dname<<std::endl;
//*****************************************************************        
   if (debug_test) cout<<lseedblnk<<" "<<lseedstart<<std::endl;       
   if (lseedstart==-1) lseedstart=time(0);
   lseedblnk=lseedstart;
   if (debug_test) cout<<lseedblnk<<" "<<lseedstart<<std::endl;   
   for (unsigned int i=1;i<=1000;i++){
        ran2(&lseedblnk);
        }
   for (unsigned int i=0;i<3;i++){  
      lseed[i]=(-1*(long) (ran2(&lseedblnk)*1e8));
   }
   if (debug_test) cout<<lseedblnk<<" "<<lseedstart<<" ran2 initialized "<<std::endl;
        
   cout<<g.get_n()<<" "<<g.p<<" "<<lseedstart<<" "; 
   //long seed=lseed[0];  
   if (debug_test) cout<<" preparing to make a graph "<<endl;  
   g.set_debug(debug_test);
   switch (grphtype){
      case -1:{
		 if (debug_test) cout<<dname<<std::endl; 
         g.read_gml(inputfile);
         break;
      }
      default:{
		 g.seed=lseed[0];
//		 if (debug_test) cout<<" build graph with "<<lseed[0]<<" "<<g.grphtype<<std::endl;
		 g.build_graph(write);
		 if (debug_test) cout<<" ok a graph is built "<<lseed[0]<<std::endl;
         break;
      }
   }

   double c=g.get_c();
      if (debug_test) cout<<g.c<<std::endl;
   cout<<"n= "<<g.get_n()<<" c= "<<g.get_c()<<" nE= "<<g.get_nE()<<" p_sick= "<<pop.p_sick;
   cout<<" contagin= "<<pop.contagin<<" fatality="<<pop.fatality<<endl;
   cout<<"p_Immune="<<pop.p_Immune<<" p_Recovery= "<<pop.p_Recovery<<endl;
//   if (false) {
//      cout<<n<<" "<<p<<" "<<c<<" "<<nE<<" "<<lseed[0]<<" "<<lseed[1]<<" "<<lseed[2];        
//   }
//   else cout<<c<<" "<<nE<<" "<<" ";
   //vce_precon(g, lseed[1], pop);
   pop.precondition(g, lseed[1]);
   if (debug_test) cout<<"preconditioned finished"<<std::endl;  
    
  // disease_pop(g, pop, lseed[2]);
   pop.pop_evolve(g, lseed[2]);   
   if (debug_test) cout<<"disease pop finished"<<std::endl;   
   
   cout<<"n= "<<g.get_n()<<" c= "<<g.get_c()<<"["<<g.p<<"]"<<" nE= "<<g.get_nE();
   cout<<" p_sick= "<<pop.p_sick<<" p_sick_max= "<<pop.n_sick_max;
   cout<<" contagin= "<<pop.contagin<<endl;
   cout<<"fatality="<<pop.fatality<<" p_Immune= "<<pop.p_Immune;
   cout<<" lifetime= "<<pop.lifetime<<" p_Recovery= "<<pop.p_Recovery<<endl;

   return 0;
}
