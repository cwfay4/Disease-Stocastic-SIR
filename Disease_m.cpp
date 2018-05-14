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
#define NUM_ITER 100
//#define grphtypedflt 0
//#define debug true



int main(int argc, char **argv){
   int argz = 1;

   double Psum=0;
   double p=3.0;
   bool write=true, dsply=false, histo=false, debug_test=true, fulloutput=false, conv=true;
   bool style=false; //is a new graph made for each pop.iteration or not
   boolcat b;
   char *filename;
   char LoPR;
   int grphtype=0;  /* 0 for rndm, 1 for pt */
   long lseed[3]={0},lseedstart=-1, lseedblnk=lseedstart;
   dPOP pop;
   pop.lifetime=10;
   pop.I=false;
   pop.simple=false;
   int dummy;
   int ITERCASE=-1;
   double ITERSTEP, ITERSTART, ITERSTOP;

   graph g;
   g.grphtype=grphtypedflt;
   
   std::stringstream ss;
   std::string inputstr, inputfile;
   
  // std::vector<std::stringstream> args;
 //  for (int i=1;i<argc;i++){
  // 	  ss.clear();
//	  ss.str(argv[i]);
 //     args.push_back(ss);
 //  }
   
 //  (argv + 1, argv + argc + !argc);
   //pass_test(args);
   //pass_test2(argc, argv);
   inputprocessing(argc, argv,  g, pop, b,  ITERCASE);
   
   if(debug_test){ 
      cout<<"g.grphtype="<<g.grphtype<<" lseedstart"<<b.lseedstart<<" ";
      cout<<"output control: w"<<b.write<<" disply"<<b.dsply<<" full"<<b.fulloutput<<" histo"<<b.histo<<" con"<<b.conv<<" debug"<<b.debug_test<<" "<<endl;
	  cout<<"disease control: I"<<pop.I<<" simple"<<pop.simple<<" ITERCASE"<<ITERCASE<<" "<<pop.I<<" "<<pop.I<<" "<<endl;
      cout<<"n= "<<g.get_n()<<" p="<<" "<<g.p<<" p_sick="<<pop.p_sick<<" p_contagin="<<pop.contagin;
      cout<<" fatality="<<pop.fatality<<" p_Immune="<<pop.p_Immune<<" lifetime="<<pop.lifetime<<endl;
   } 

//****** Set output file name ************************************
   //std::ofstream dout("dout.csv");
   std::string oname=outputname(g,pop);
   //std::string oname="disease-ss.csv";
   if (debug_test) cout<<oname<<std::endl;
//**************************************************************** 
//****** Random number prep **************************************      
   if (debug_test) cout<<" seed"<<lseedblnk<<" "<<b.lseedstart<<std::endl;       
   if (b.lseedstart==-1) b.lseedstart=time(0);
   lseedblnk=b.lseedstart;
   if (debug_test) cout<<" seed"<<lseedblnk<<" "<<b.lseedstart<<std::endl;   
   for (unsigned int i=1;i<=1000;i++){
        ran2(&lseedblnk);
        }
   for (unsigned int i=0;i<3;i++){  
      lseed[i]=(-1*(long) (ran2(&lseedblnk)*1e8));
   }
   if (debug_test) cout<<lseedblnk<<" "<<b.lseedstart<<" ran2 initialized "<<std::endl;
   if (debug_test) cout<<g.get_n()<<" "<<g.p<<" "<<b.lseedstart<<" "; 
   if (debug_test) cout<<" preparing to make a graph "<<endl; 
   
//**************************************************************** 
//****** Set iteration values ************************************  
   g.set_debug(debug_test);
   if (g.grphtype ==-1) g.grphtype=0; //we are not loading a graph so make sure we are using the default - this will be replaced when we load graphs
   switch (ITERCASE){
      case 1:{  //iterate over contagin / virality
         pop.contagin=0.1;
         ITERSTEP=0.01;
         ITERSTART=pop.contagin;
         ITERSTOP=10;
         break;
      }
      case 2:{ //iterate over fatality
		 pop.fatality=0.01;
         ITERSTEP=0.01;
		 ITERSTART= pop.fatality;
		 ITERSTOP=1.00;
         break;
      }        
      case 3:{  //iterate over lifetime
		 pop.lifetime=0.001*(double)NUM_STEPS;
		 ITERSTART=pop.lifetime;
		 ITERSTOP=0.1*(double)NUM_STEPS;
		 ITERSTEP=0.001*(double)NUM_STEPS;
         break;
      }         
      default:{ //iterate over p
		 g.p=0.1;
		 ITERSTEP=0.1;
         ITERSTART=g.p;
         ITERSTOP=10;
         if (ITERSTOP>g.get_n()-1) ITERSTOP=g.get_n()-1; 		       
         break;
      }
   }

   pop.iteration=1;	  //controls the number of graphs averaged for each trial.  

   double dn=(double)NUM_ITER;  
   while (ITERSTART<ITERSTOP){ //scan over interested quantity
	      //std::vector <double> datum(10,0); //datum is a vector containing all of the average results, dadum sq is the square of the average results
	     // std::vector <double> datumsq(10,0);//from these two we will calculate the standard error
          //datum[0]=(double)g.get_n();
         // datumsq[0]=datum[0]*datum[0];	
          
        std::vector <Stats> datum;
        for (int i=0; i<8; i++) { //this will contain the max/min etc values for each NUM_ITER value;
            Stats data(dn);
            datum.push_back(data);
        }
        datum[0].value=(double)g.get_n();
        datum[0].value=(double)g.get_n()*(double)g.get_n();
		 		  
	    if (style==false){ //style controls if the average for NUM_ITER is many evolutions of one graph or one evolution of many graphs
	       g.seed=lseed[0]; //false means that one graph is made and evolved several times.
	       g.build_graph(false);
		}   
          
       // if (pop.iteration>1){ //we have more than one graph for each iterator value, thus we will average over each pop.evolve step
        	std::vector <Stats> evolve_steps; //this will save the values for each pop.evolve step;
            for (int i=0; i<8; i++) { //this will contain the pop.evolve steps;
            	Stats evstepsdata((double) pop.iteration);
            	evolve_steps.push_back(evstepsdata);
            }   
	    //}
	    
	   for (int i=0; i<pop.iteration; i++){  //how many graphs do we average over per each iterator
	        if (style==true){
	           g.seed=lseed[0]; //change this
	           g.build_graph(false);
		    }   
	        if(false)cout<<" type: "<<g.grphtype<<" ITERNUM "<<dn<<" n: "<<g.get_n()<<" nE: "<<g.get_nE()<<" nE: "<<g.calc_nE()<<" c: "<<g.get_c()<<endl;
	  
	      //vce_precon(g, lseed[1], pop);
	      //disease_pop(g, pop, lseed[2]);
	      
	        pop.precondition(g, lseed[1]);
	      
	        if (pop.iteration<=1) pop.pop_evolve(g, lseed[2]);
	        else {
		       pop.pop_evolve(g,lseed[2],evolve_steps);
		    }
	    
	        datum[1].add_data((double)g.get_nE());
            datum[2].add_data((double)g.get_c());
            datum[3].add_data((double)pop.n_sick/(double)g.get_n());
            datum[4].add_data((double)pop.n_healthy/(double)g.get_n());
            datum[5].add_data((double)pop.n_immune/(double)g.get_n());
            datum[6].add_data(pop.n_dead/(double)g.get_n());
            datum[7].add_data(pop.n_sick_max/(double)g.get_n());          
           //output the convergence for this value of ITERSTART
            for(int j=0; j<evolve_steps.size(); j++){ //calculate the standard error
	   	        double var = evolve_steps[j].standard_error();
	        }
	        bool conv=true;
	        if (conv){ 
	        
	            std::ofstream output("Disease_Conv.csv", ios_base::app);
	            for (int k=0;k<pop.POP_STEPS; k++){
	                output<<k+1<<", "; //for every step in pop.evolve
	                for(int j=0; j<evolve_steps.size(); j++){ //output averages 
	   	                output<<evolve_steps[1].value<<", "<<evolve_steps[1].std_error<<", ";
	                } 
	                output<<endl;
	            }                
	            output.close();
            } 
			//********************************************************************   
	   }
	   
	   for(int j=0; j<datum.size(); j++){ //calculate the standard error
	   	  double var = datum[j].standard_error();
	   }
	   
	   std::ofstream output(oname.data(), ios_base::app);
	   //std::ofstream output("d-steadystate.csv", ios_base::app);
	   output<<g.nodes.size()<<", "<<ITERSTART<<", "; //graph size and iterator
	  //********* initial parameters **************************************
	   output<<pop.p_sick<<", "<<pop.contagin<<", "<<pop.fatality<<", "<<pop.p_Immune;
       output<<", "<<pop.lifetime<<", "<<pop.I<<", ";
	   //********** runtime results and std error ***************************       
       for (int k = 1; k<=7; k++){
          output<<datum[k].value<<", "<<datum[k].std_error<<", ";  //probably could calculate the standard error here.
       }
       output<<endl;
	   output.close();  
	   cout<<"ITERNUM "<<dn<<" n: " <<datum[0].value<<" ITERSTART: "<<ITERSTART<<" nE: "<<datum[1].value/dn<<" c: "<<datum[2].value/dn<<" n_Sick: "<<datum[3].value/dn;
	   cout<<" n_H: "<<datum[4].value/dn<<" n_Immune: "<<datum[5].value/dn<<" n_dead: "<<datum[6].value/dn<<" n_sick_max: "<<datum[7].value/dn<<" "<<pop.contagin<<" "<<pop.fatality<<" "<<pop.lifetime<<endl;
  
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
       }
   }

   return 0;
}
