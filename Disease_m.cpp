/***********************************************************************/
/** Disease stocastic SIR  Control Program                           ***/
/**                                                                  ***/
/** C. Fay August 2005                                               ***/
/** Version 2.1.0.1   14.05.2005                                     ***/
/** Run: Disease_m gml_filname seed                                  ***/
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
#include "stats.hpp"
#include "graph.hpp"
#include "diseasePOP.hpp"
#include "Disease_basic_utils.hpp"

using namespace std;

//#define NUM_STEPS 500 //number of time steps replaced by POP_STEPS
//#define NUM_ITER 100
//#define grphtypedflt 0  //this might also be in the graph.hpp
//#define debug true

int main(int argc, char **argv){
   //initialization of variables
   int argz = 1;
//   double Psum=0;
  // double p=3.0;
//   bool write=true, dsply=false, histo=false, debug_test=false, fulloutput=false, conv=true;
//   bool style=false; //is a new graph made for each pop.iteration or not
   boolcat b;
   b.NUM_ITER=10;
//   char *filename;
//   char LoPR;
  // int grphtype=0;  /* 0 for rndm, 1 for pt */
   long lseed[3]={0},lseedstart=-1, lseedblnk=lseedstart;
   dPOP pop(0,0,0,0,0.1,0.1,0.1,0.1,10,false,false);
 //  pop.iteration=10;	pop.lifetime=10;    pop.I=false;    pop.simple=false;
//   int dummy;
   int ITERCASE=-1;
   double ITERSTEP, ITERSTART, ITERSTOP;
   
   graph g;
   g.grphtype=grphtypedflt;
   g.set_debug(false);
   
 //  std::stringstream ss;
 //  std::string inputstr, inputfile;
   
   inputprocessing(argc, argv,  g, pop, b,  ITERCASE); //this reads the command line for all of the relevant controls and values
   
   //testing
   
   if(true){ 
      cout<<"g.grphtype="<<g.grphtype<<" lseedstart"<<b.lseedstart<<" ";
      cout<<"output control: w"<<b.write<<" disply"<<b.dsply;
	  cout<<" full"<<b.fulloutput<<" histo"<<b.histo;
	  cout<<" con"<<b.conv<<" debug"<<b.debug_test<<" "<<endl;
	  cout<<"disease control: I"<<pop.I<<" simple"<<pop.simple;
	  cout<<" ITERCASE"<<ITERCASE<<" "<<pop.I<<" "<<pop.I<<" "<<endl;
      cout<<"n= "<<g.get_n()<<" p="<<" "<<g.p<<" p_sick="<<pop.p_sick;
	  cout<<" p_contagin="<<pop.contagin;
      cout<<" fatality="<<pop.fatality<<" p_Immune="<<pop.p_Immune;
	  cout<<" lifetime="<<pop.lifetime<<endl;
   } 

//****** Set output file name ************************************
   //std::ofstream dout("dout.csv");
   std::string oname=outputname(g,pop);
   //std::string oname="disease-ss.csv";
   if (b.debug_test) cout<<oname<<std::endl;
//**************************************************************** 
//****** Random number prep **************************************      
   if (b.debug_test) cout<<" seed"<<lseedblnk<<" "<<b.lseedstart<<std::endl;       
   if (b.lseedstart==-1) b.lseedstart=time(0);
   lseedblnk=b.lseedstart;
   if (b.debug_test) cout<<" seed"<<lseedblnk<<" "<<b.lseedstart<<std::endl;   
   for (unsigned int i=1;i<=1000;i++){
        ran2(&lseedblnk);
        }
   for (unsigned int i=0;i<3;i++){  
      lseed[i]=(-1*(long) (ran2(&lseedblnk)*1e8));
   }
   if (b.debug_test) cout<<lseedblnk<<" "<<b.lseedstart<<" ran2 initialized "<<std::endl;
   if (b.debug_test) cout<<g.get_n()<<" "<<g.p<<" "<<b.lseedstart<<" "; 
//**************************************************************** 
//****** Set iteration values ************************************  
   //g.set_debug(b.debug_test);
   if (g.grphtype ==-1) g.grphtype=0; //we are not loading a graph so make sure we are using the default - this will be replaced when we load graphs
   if (b.debug_test) cout<<" intial graph parameters set: preparing to make a graph "<<endl;    
   
   switch (ITERCASE){
      case 1:{  //iterate over contagin / virality
         ITERSTOP=pop.contagin;
         pop.contagin=0.1;
         ITERSTEP=0.01;
         ITERSTART=pop.contagin;
         break;
      }
      case 2:{ //iterate over fatality
         ITERSTOP=pop.fatality;
		 pop.fatality=0.01;
         ITERSTEP=0.01;
		 ITERSTART= pop.fatality;
         break;
      }        
      case 3:{  //iterate over lifetime (make lifetime a function of the number of time steps
		 pop.lifetime=0.01*(double)pop.POP_STEPS;
		 ITERSTART=pop.lifetime;
		 ITERSTOP=0.1*(double)pop.POP_STEPS;
		 ITERSTEP=0.001*(double)pop.POP_STEPS;
         break;
      }         
      default:{ //iterate over p
		 //ITERSTOP=g.p;
		 //g.p=0.1;
		 ITERSTEP=0.1;
         ITERSTART=g.p;
         ITERSTOP=1.5;
         if (ITERSTART>ITERSTOP) ITERSTOP=ITERSTART+ITERSTEP*10;
         if (ITERSTOP>g.get_n()-1) ITERSTOP=g.get_n()-1; 
		 if (b.debug_test) cout<<"pc ITERSTART,ITERSTOP "<<ITERSTART<<","<<ITERSTOP<<endl;		       
         break;
      }
   }
   
   if (b.debug_test) cout<<"Initial iteration parameters set"<<endl;

   // In this group(s) of iterations we track stats for two collections.  Each pop.evolve is run for POP_STEPS
   // time steps.  We run pop.Evolve_STEPS number of evolutions.  We run evolutions over each iteration step, for a total of
   // macro_it_count iteration values.

   double dn=(double)pop.POP_STEPS; //dn is the double prescision number of time steps here for calculations
   iteration_stats macrodata;  //will have macro_it_count values averaged over Evolve_STEPS number of evolutions
   // macrodata stores: n, nE, c, n_sick, n_healthy, n_immune, n_dead, n_sick_max all values are final state values.
   int macro_it_count=0;
   
   if (b.debug_test) cout<<"Prepared for iteration loop"<<endl;   
   while (ITERSTART<ITERSTOP){ //scan over interested quantity
        if (b.debug_test) cout<<"iteration loop "<<macro_it_count<<" "<<ITERSTART<<" "<<ITERSTOP<<endl; 
        
        macrodata.add_elementsof(8);//adds a new element for each macro_it_count;  it adds with 8 spaces for macro data.
        if (b.debug_test) cout<<"macrodata.add_element"<<macro_it_count+1<<endl;         
        macrodata.increase_element(macro_it_count, 0, (double)g.get_n());
        
		if (b.debug_test) cout<<"macrodata.increased_element"<<macro_it_count<<endl;             	    
        if (b.debug_test) cout<<"building graph"<<endl;   	    
	    
	    if (b.style==false){ //style controls if the average for NUM_ITER is many evolutions of one graph or one evolution of many graphs
	       g.seed=lseed[0];  //false means that one graph is made and evolved several times.
	       g.build_graph(false);
		}   
       if (b.debug_test) cout<<"Iteration loop built graph"<<endl;   
       // if (pop.Evolve_STEPS>1){ //we have more than one graph for each iterator value, thus we will average over each pop.evolve step
       iteration_stats microdata(pop.POP_STEPS, 4); //this will save the values for each pop.evolve step per one iteration;
	    
	   for (int i=0; i<pop.Evolve_STEPS; i++){  //how many graphs do we average over per each iterator
	        if (b.style==true){
	           g.seed=lseed[0]; //change this
	           g.build_graph(false);
		    }   
	        if(b.debug_test)cout<<" type: "<<g.grphtype<<" ITERNUM "<<dn<<" n: "<<g.get_n()<<" nE: "<<g.get_nE()<<" nE: "<<g.calc_nE()<<" c: "<<g.get_c()<<endl;
	  		//if(true) cout<<"Number of Edges"<<g.get_nE()<<endl;
	      //vce_precon(g, lseed[1], pop);
	      //disease_pop(g, pop, lseed[2]);
	      
	        pop.precondition(g, lseed[1]);
	      
	        if(b.debug_test)cout<<" precondition finished:  n_sick "<<pop.n_sick<<" "<<pop.n_healthy<<" "<<pop.n_immune<<endl;
	        
	        if (pop.n_sick>=1) { //just incase there is no one sick after the precondition
	            	            
			    if (pop.Evolve_STEPS<=1) {
			    	if (b.debug_test)cout<<"run pop_evolve(g,lseed[2])"<<endl;
			    	pop.pop_evolve(g, lseed[2]);
				}
	            else {
	            	if (b.debug_test)cout<<"run pop_evolve(g,lseed[2],microdata)"<<endl;
	            //	pop.debug=b.debug_test;
		            pop.pop_evolve(g,lseed[2],microdata); //problem here
		        }
	    	}
	        
	        if(b.debug_test)cout<<" pop evolve finished: "<< g.get_n()<<" "<<g.get_nE()<<" "<<g.get_c()<<" n_sick "<<pop.n_sick<<" "<<pop.n_healthy<<" "<<pop.n_immune<<" "<<pop.n_dead<<" "<<pop.n_sick_max<<endl;
            double doublenum=(double)g.get_n();
			macrodata.increase_element(macro_it_count, 1, (double)g.get_nE());
	        macrodata.increase_element(macro_it_count, 2, g.get_c());
	        macrodata.increase_element(macro_it_count, 3, (double)pop.n_sick/doublenum);
	        macrodata.increase_element(macro_it_count, 4, (double)pop.n_healthy/doublenum);
	        macrodata.increase_element(macro_it_count, 5, (double)pop.n_immune/doublenum);
	        macrodata.increase_element(macro_it_count, 6, (double)pop.n_dead/doublenum);
	        macrodata.increase_element(macro_it_count, 7, (double)pop.n_sick_max/doublenum);
	        if(b.debug_test)cout<<"macro data stored"<<endl;
							         
	   } //End of Pop.iteration loop
	   
	    if (pop.Evolve_STEPS>1){
	       microdata.calc_stats(pop.Evolve_STEPS);
	       microdata.output_elements(pop.POP_STEPS, 4, "Dconv.csv"); //This will output the average values for each step of the pop.iteration
	    }
	   
	   
	   //for(int j=0; j<a+-datum.size(); j++){ //calculate the standard error
	   //	  double var = datum[j].standard_error();
	   //}
	   
 	   cout<<"ITERNUM "<<dn<<" n: " <<macrodata.get_element(macro_it_count,0)<<" ITERSTART: "<<ITERSTART;
	   cout<<" nE: "<<macrodata.get_element(macro_it_count,1)<<" c: "<<macrodata.get_element(macro_it_count,2);
	   cout<<" n_Sick: "<<macrodata.get_element(macro_it_count,3)<<" n_H: "<<macrodata.get_element(macro_it_count,4);
	   cout<<" n_Immune: "<<macrodata.get_element(macro_it_count,5)<<" n_dead: "<<macrodata.get_element(macro_it_count,6);
	   cout<<" n_sick_max: "<<macrodata.get_element(macro_it_count,7);
	   cout<<" "<<pop.contagin<<" "<<pop.fatality<<" "<<pop.lifetime<<endl;
  
	   switch (ITERCASE){  //iterate your interested quantity
            case 1:{//iterate over contagin / virality
               pop.contagin+=ITERSTEP;
               ITERSTART=pop.contagin;
               break;
            }
            case 2:{//iterate over fatality
               pop.fatality+=ITERSTEP;
               ITERSTART=pop.fatality;
               break;
            }        
            case 3:{//iterate over lifetime
               pop.lifetime+=ITERSTEP;
               ITERSTART=pop.lifetime;
               break;
            }         
            default:{//iterate over p 
               g.p+=ITERSTEP;
               ITERSTART=g.p;
              // cout<<" g.p "<<g.p<<endl;
               break;
            }
       } //end ITERCASE Switch
         macro_it_count++;   
    }//End Iterstop while loop
    macrodata.calc_stats((double)macro_it_count);
	macrodata.output_elements(macro_it_count, 8, "d-steadystate.csv"); //This will output the average values for each step of the pop.iteration 

   return 0;
}
