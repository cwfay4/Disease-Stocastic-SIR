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
#define NUM_ITER 10
//#define grphtypedflt 0
//#define debug true

void add_data(std::vector<double>& data, std::vector<double>& datasq, int j, double datum){
   //cout<<j<<endl;
   data[j]=data[j]+datum; //sum of datum
   datasq[j]=data[j]+datum*datum; //sum of square of data
return;
}
double standard_error(double data, double datasq, double dsamp){
	double var=sqrt(abs(datasq-(data*data/dsamp)))/sqrt(dsamp-1);//calculate std error;
	return var;
}
std::string outputname(graph g, dPOP pop){
   std::string gt="p", rnt="r2", sbc="f", pre="p";
   std::ostringstream strs;
   strs << g.get_n();
   std::string nstr = strs.str();
   strs << g.p;
   std::string pstr = strs.str();  
   if (g.bc==1) sbc="p";
   else if (g.bc==2) sbc="pp";
   if (g.grphtype == 0 || g.grphtype==7) gt="r";
   else if (g.grphtype == 2) rnt="r3";
   else if (g.grphtype == 3) gt="f";
   else if (g.grphtype == 4) gt="c";
   else if (g.grphtype == 5) gt="sq";
   else if (g.grphtype == 6) gt="sclf";
   else if (g.grphtype == 8) gt="bb";
   //std::string gname="d-"+gt+nstr+"-"+pstr+"-"+sbc+rnt;
   std::string gname="d"+nstr+sbc+rnt+"ss.csv";
   
   //strs << pop.p_sick;
   //std::string psstr = strs.str();
   //strs << pop.p_Immune;
   //std::string pistr = strs.str(); 
   //strs << pop.contagin;
   //std::string constr = strs.str();
   //strs << pop.fatality;
   //std::string fatalstr = strs.str(); 
  //strs << pop.lifetime;
   //std::string lifelstr = strs.str();
   //if (ITERCASE==1) pre="v";
   //else if(ITERCASE==2) pre="f";
   //else if(ITERCASE==3) pre="l";   
   //if (pop.simple) pre+="s";   
   //if (pop.I) pre+="I";
   //std::string dname=pre+"-"+psstr+"-"+pistr+"-"+constr+"-"+fatalstr+"-"+lifelstr;
   
  // std::ofstream dout("dout.csv");
   //std::string oname=gname+dname+"ss.csv";
   //std::string oname="disease"+"ss.csv";
  // cout<<gname<<endl;
   //cout<<dname<<endl;
  // if (false) cout<<dname<<std::endl;

    return gname;
}
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
   pop.I=false;
   pop.simple=false;
   int dummy;
   int ITERCASE=0;
   double ITERSTEP, ITERSTART, ITERSTOP;

   graph g;
   g.grphtype=grphtypedflt;
   
   std::stringstream ss;
   std::string inputstr, inputfile;

//******** Input processing ***********************************************************
   int i, j=0;
   for (i = 1; i < argc; i++) {
   	  ss.clear();
	  ss.str(argv[i]);
      ss>>inputstr; 
      if(inputstr.compare("-pt") == 0)        g.grphtype=1;  //planar triangular graph
      else if(inputstr.compare("-ptr3") == 0) g.grphtype=2;  //planar triangular graph with ran3
      else if(inputstr.compare("-fcc") == 0)  g.grphtype=3;  //FCC Lattice
      else if(inputstr.compare("-cubic") == 0)g.grphtype=4;  //cubic Lattice
      else if(inputstr.compare("-sq") == 0)   g.grphtype=5;  //squareLattice     
      else if(inputstr.compare("-sclf") == 0) g.grphtype=6;  //scale free Lattice
      else if(inputstr.compare("-rnd2") == 0) g.grphtype=7;  //random Lattice type 2
      else if(inputstr.compare("-bb") == 0) {                //barabasi scale free
      	g.grphtype=8;
	    i++;
		ss>>dummy;
		dummy=atoi(argv[i]);
	    g.set_L(dummy);
	   // g.set_L(atoi(argv[argz]));
	    //cout<<"nkernal: "<<L<<" ";
	  }  
      else if(inputstr.compare("-rndfixedc") == 0||inputstr.compare("-rf") ==0){ //random fixed connectivity
      	g.grphtype=9;
	    i++;
		ss>>dummy;
		dummy=atoi(argv[i]);
	    g.set_L(dummy);
	    if (debug_test) cout<<"fixed maxc: "<<g.get_L()<<" ";
	  }
      else if(inputstr.compare("-rnd") == 0)        g.grphtype=0;   //random graph 
      //************************End of Graph Types******************************** 
      //************************Boolean control***********************************
      else if(inputstr.compare("-w") == 0)       	write=true;
      else if(inputstr.compare("-dsply") == 0)     	dsply=true;
      else if(inputstr.compare("-fullout") == 0)   	fulloutput=true;
      else if(inputstr.compare("-histo") == 0)     	histo=true; 
      else if(inputstr.compare("-conv") == 0)     	conv=true;  
      else if(inputstr.compare("-debug") == 0)     	debug_test=true;
      else if(inputstr.compare("-h") == 0){  //help menu
		Disease_usage(false); 
		return 0;  
	  }     	
      else if(inputstr.compare("-hv") == 0){  //help menu verbose
		Disease_usage(true); 
		return 0;  
	  }	  
      //************************Disease control***********************************      
      else if(inputstr.compare("-I") == 0)          pop.I=true;   //Immune after sickness 
      else if(inputstr.compare("-simple") == 0)     {
		   pop.simple=true;   //no random vectors, n_sick counts immunity/death/recovered   
		   if (pop.I==true) pop.I=false;
	   }    
	  else if(inputstr.compare("-simpleI") == 0)     {
		   pop.simple=true;   //no random vectors, n_sick counts immunity/death/recovered   
		   pop.I=true;
	   }      
      else if(inputstr.compare("-pc") == 0)         ITERCASE=0; //iterate over the connectivity
      else if(inputstr.compare("-virality") == 0)   ITERCASE=1; //iterate over the virality
      else if(inputstr.compare("-fatality") == 0)   ITERCASE=2; //iterate over the fatality
      else if(inputstr.compare("-lifetime") == 0)   ITERCASE=3; //iterate over the lifetime
      //************************read graph ***********************************       
      else if(inputstr.compare("-rgml") ==0||inputstr.compare("-r") ==0){  //read a file in gml format
	     i++;
	     ss.clear();
	     ss.str(argv[i]);
	     ss>>inputfile;
	     g.grphtype=-1;
         if (debug_test) cout<<"You want to read the file "<<inputfile<<std::endl;
	  }
	  //*** Read n p p_sick contagin fatality p_Immune lifetime **************
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
		  if (debug_test) cout<<"p="<<ss.str()<<" "<<g.p<<std::endl;
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
		  if (pop.p_Immune > 1.0) {
			  cout<<"100% or more immune I am bored.  Let's start with no one immune."<<endl;
			  pop.p_Immune=0.0;
		  }
		  if (debug_test) cout<<"p_Immune="<<ss.str()<<" "<<pop.p_Immune<<std::endl;
		  j++;
	  }
	  else if(j==6){ //sickness duration
		  pop.lifetime=atof(argv[i]); //works
		  if (debug_test) cout<<"lifetime="<<ss.str()<<" "<<pop.lifetime<<std::endl;
		  j++;
	  }	  
	  //**** if bb graph give size of kernal *********************************
      if (g.grphtype==8) cout<<"nkernal: "<<g.get_L()<<" "; // dimension of lattice
      if (g.grphtype==3 || g.grphtype==4) { // dimension of lattice
         g.set_L(g.get_n());
         g.set_n(std::pow(g.get_L(),3));
      }   
   }
//******* End of input processing ********************************  
   cout<<"n= "<<g.get_n()<<" c= "<<g.get_c()<<" nE= "<<g.get_nE()<<" p_sick= "<<pop.p_sick;
   cout<<" contagin= "<<pop.contagin<<" fatality="<<pop.fatality<<" p_Immune="<<pop.p_Immune;
   cout<<" lifetime= "<<pop.lifetime<<" I.bool "<<pop.I<<endl;
//*****************************************************************        
/****** Set output file name *************************************/
   //std::ofstream dout("dout.csv");
   std::string oname=outputname(g,pop);
   //std::string oname="disease-ss.csv";
   if (debug_test) cout<<oname<<std::endl;
/*****************************************************************/ 
/****** Random number prep ***************************************/      
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
   if (debug_test) cout<<g.get_n()<<" "<<g.p<<" "<<lseedstart<<" "; 
   if (debug_test) cout<<" preparing to make a graph "<<endl; 
   
/*****************************************************************/ 
/****** Set iteration values *************************************/   
   g.set_debug(debug_test);
   if (g.grphtype ==-1) g.grphtype=0; //we are not loading a graph so make sure we are using the default 
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

   int ITERNUM=0;	      

   double dn=(double)NUM_ITER;
   while (ITERSTART<ITERSTOP){ //scan over interested quantity
	      std::vector <double> datum(10,0);
	      std::vector <double> datumsq(10,0);
          datum[0]=(double)g.get_n();
          datumsq[0]=datum[0]*datum[0];	      
	   for (i=0; i<NUM_ITER; i++){ //how many graphs do we average over
		 //cout<<"Call graph clear "; 
		  //g.clear();
		 // cout<<"clear 1 "; 
		 // g.clear();
		  //cout<<"clear 2: Num Edges: "<<g.get_nE();
	     // pop.clear(); 
	      g.seed=lseed[0]; //change this
	      g.build_graph(false);
	      if(false)cout<<" type: "<<g.grphtype<<" ITERNUM "<<dn<<" n: "<<g.get_n()<<" nE: "<<g.get_nE()<<" nE: "<<g.calc_nE()<<" c: "<<g.get_c()<<endl;
	  
	      //vce_precon(g, lseed[1], pop);
	      //disease_pop(g, pop, lseed[2]);
	      
	      pop.precondition(g, lseed[1]);
	      pop.pop_evolve(g, lseed[2]);
	      
	      datum[1]+=(double)g.get_nE();
	      datumsq[1]+=(double)g.get_nE()*(double)g.get_nE();
	      datum[2]+=g.get_c();
	      datumsq[2]+=(double)g.get_c()*(double)g.get_c();
	      datum[3]+=pop.n_sick/(double)g.get_n();
	      datumsq[3]+=(pop.n_sick/(double)g.get_n())*(pop.n_sick/(double)g.get_n());
	      datum[4]+=pop.n_healthy/(double)g.get_n();
	      datumsq[4]+=(pop.n_healthy/(double)g.get_n())*(pop.n_healthy/(double)g.get_n());	      
	      datum[5]+=pop.n_immune/(double)g.get_n();
	      datumsq[5]+=(pop.n_immune/(double)g.get_n())*(pop.n_immune/(double)g.get_n());		      
	      datum[6]+=pop.n_dead/(double)g.get_n();  
	      datumsq[6]+=(pop.n_dead/(double)g.get_n())*(pop.n_dead/(double)g.get_n());	      
	      datum[7]+=pop.n_sick_max/(double)g.get_n();
	      datumsq[7]+=(pop.n_dead/(double)g.get_n())*(pop.n_dead/(double)g.get_n());	       
	   }
	   std::ofstream output(oname.data(), ios_base::app);
	   //std::ofstream output("d-steadystate.csv", ios_base::app);
	   output<<g.nodes.size()<<" "<<ITERSTART<<" "; //graph size and iterator
	   /********** initial parameters ***************************************/
	   output<<pop.p_sick<<" "<<pop.contagin<<" "<<pop.fatality<<" "<<pop.p_Immune;
       output<<" "<<pop.lifetime<<" "<<pop.I<<" ";
	   /********** runtime results and std error ****************************/       
       for (int k = 1; k<=7; k++){
		  double var=standard_error(datum[k],datumsq[k],dn);
          output<<datum[k]/dn<<" "<<var<<" ";
       }
       output<<endl;
	   output.close();  
	   cout<<"ITERNUM "<<dn<<" n: "<<datum[0]<<" ITERSTART: "<<ITERSTART<<" nE: "<<datum[1]/dn<<" c: "<<datum[2]/dn<<" n_Sick: "<<datum[3]/dn;
	   cout<<" n_H: "<<datum[4]/dn<<" n_Immune: "<<datum[5]/dn<<" n_dead: "<<datum[6]/dn<<" "<<pop.contagin<<" "<<pop.fatality<<" "<<pop.lifetime<<endl;
  
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
