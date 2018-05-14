/***********************************************************************/
/** crprcfns                                                         ***/
/**                                                                  ***/
/** C. Fay Jan 2003                                                  ***/
/** Version 2.0.0.0   05.10.2004                                     ***/
/** Run: crprcfnsgraph gml_filename                                  ***/
/** Output: dat files that can be opened by xmgrace as a visual      ***/
/**        representation of the graph                               ***/
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

#include "diseasePOP.hpp"
#include "graph.hpp"
#include "Disease_basic_utils.hpp"
//#include "tp.hpp"

//#include "LoPR_alg_o.hpp"

using namespace std;

boolcat::boolcat(){
   write=false;
   dsply=false;
   histo=false;
   debug_test=false;
   fulloutput=false;
   conv=false;	
   ITERCASE=-1;
   lseedstart=-1;
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

void Disease_usage(bool verbose) {
  cout<<"Disease : Disease population \n";
  cout<<"   Usage: [-h][-rgml <filename>][GRAPHTYPE] n c p_sick contagin fatality p_Immune lifetime"<<endl;
  cout<<endl;
  cout<<"   Output: original parameters and debug parameters "<<endl;
  cout<<"   n=[] c=[] nE=[] p_sick=[] contagin=[] fatality=[] p_Immune=[] lifetime=[] "<<endl;  
  cout<<"   [n] [c] [seed] graph written as [graph name. gml]  "<<endl;    
  cout<<"   n=[] c=[] nE=[] p_sick=[] contagin=[] fatality=[] p_Immune=[]  "<<endl;
  cout<<endl;  
  cout<<"   Output file: Disease_Conv.csv [comma delinated]"<<endl;
  cout<<"   [Iteration] [n_healthy] [n_sick] [n_immune] [n_dead]  "<<endl; 
  cout<<endl;
  cout<<"Disease_m : Disease population averaged over multiple graphs \n";
  cout<<"   Usage: [-h][-rgml <filename>][GRAPHTYPE] n c p_sick contagin fatality p_Immune lifetime"<<endl;
  cout<<endl;
//  cout<<"   Output: original parameters and debug parameters "<<endl;
//  cout<<"   n=[] c=[] nE=[] p_sick=[] contagin=[] fatality=[] p_Immune=[] lifetime=[] "<<endl;  
//  cout<<"   [n] [c] [seed] graph written as [graph name. gml]  "<<endl;    
//  cout<<"   n=[] c=[] nE=[] p_sick=[] contagin=[] fatality=[] p_Immune=[]  "<<endl;
//  cout<<endl;   
  cout<<"   Output file: d+[n]+[bc]+[rnt]+ss.csv [comma delinated]"<<endl;
  cout<<"   [n] [initial iteration value] [p_sick] [contagin] [fatality] [Immune] [lifetime] [nE] [c] [n_healthy] [n_Immune] [n_dead]"<<endl; 
  cout<<endl;   
  cout<<endl;  
  if (!verbose)
    cout<<"Type '%s -hv' for long help\n";
  else {
   // cout<<endl<<"Revision "<<tp_version<<" Charles Fay 2006"<<endl;
   // cout<<"    grphfns_version "<<gf_version<<" Charles Fay 2006"<<endl;
   // cout<<"    LoPRfns_version "<<lp_version<<" Charles Fay 2006"<<endl<<endl;
    cout<<"   -h: help \n"<<endl;
    cout<<"GRAPH OPTIONS:"<<endl;
    cout<<"   -rnd: random graph using ran2 from numerical methods \n";
    cout<<"   -pt: diluted planar triangular lattice using ran2 \n"; 
    cout<<"   -ptr3: diluted planar triangular lattice using ran3\n";    
    cout<<"   -fcc: diluted fcc lattice using ran2 \n";
    cout<<"   -cubic: diluted cubic lattice using ran2 \n";
    cout<<"   -sq: diluted square lattice using ran2 \n";    
    cout<<"   -sclf: scale free graph using ran2\n";
    cout<<"   -staticsclf: static scale free graph of the type from cond-mat/0312336 \n";     
    cout<<"   -rnd2: random graph sampling every edge using ran2 \n";
    cout<<"   -rnd3: using same procedure as the rf algorithm \n"; 
//    cout<<"   -bb [m]: barabasi network, m is integer number of nodes in kernal \n";
//    cout<<"            (c should be an integer as well)\n"; 
    cout<<"   -bb barabasi network, 5% of nodes as part of the kernal \n";
    cout<<"            (c should be an integer as well)\n";     
    cout<<"   -rf [m]: a random graph, connectivies are controlled so that c<=m \n";
    cout<<"   -rf_t2 [m]: a random graph, connectivies are controlled so that c<=m, different generator \n";
    cout<<"   -rf_64 [m]: a random graph, connectivies are controlled so that c<=m \n";
    cout<<"               optimized for a 64bit turion processor. \n";
//    cout<<"   -staticSF: static scale free graph of the type from cond-mat/0312336 \n";
//    cout<<"   -Ksat [m]: Ksat with m variables per clause,  and c clauses (not finished)\n"; 
    cout<<endl;   
    cout<<"   -bc0 : free boundary conditions (default)"<<endl;    
    cout<<"   -bc1 : periodic one one side (only option for 2-d graphs) \n";
    cout<<"   -bc2 : periodic on two sides \n"; 
//    cout<<"LoPR OPTIONS:"<<endl;   
//    cout<<"   -sb: site LoPR then bond LoPR \n";
//    cout<<"   -sb1: site LoPR then bond LoPR with one node list\n";    
//    cout<<"   -sDig: site LoPR then site based DIG \n"; 
//    cout<<"   -sDbD: site LoPR then site based DIG, bond LoPR then site based DIG\n";
//    cout<<"   -sDbD_1: site LoPR then site based DIG, bond LoPR then site based DIG\n"; 
//    cout<<"   -sDbD_m: site LoPR then site based DIG, bond LoPR then site based DIG\n";
 //   cout<<"   -greedy: greedy cover by selection of a random node \n";
//    cout<<"   -Gmaxc: greedy cover by selection of the node with largest current connectivity\n";             
//    cout<<"GENERAL OPTIONS:"<<endl;
    cout<<"   -w: generate gml files of the graph. \n";
    cout<<"   -dsply: generate csv files that can be plotted to make a \n";
    cout<<"           visual output of the graph. \n";
    cout<<"   -fullout: generates full output. \n";
    cout<<"   -histo: generates histo graph of node probabilities. \n";
    cout<<"   -conv: generates convergence data of graphs. \n";
    cout<<"   -seed [seed]: run with this seed\n";
    cout<<"               (if comparing to a negative seed make [seed] positive)\n";
    cout<<"   -rgml [filename]: read gml format graph file  \n";
    cout<<endl;
//    cout<<"Output is of the form: \n";
//    cout<<"    n c seed1 seed2 seed3 c_actual nE \n";
//    cout<<"    sELoPR_covered_Fraction FrozenFraction \n";
//    cout<<"    bELoPR_covered_Fraction FrozenFraction \n";
    cout<<endl;
    cout<<"Variables: \n";
    cout<<"    n: number of nodes \n";
    cout<<"    c: connectivity (average number of edges per node) \n";
    cout<<"    nsamples: number of graphs to sample \n";
    cout<<"    c_actual: connectivity of graph (average number of edges per node)\n";    
    cout<<"    nE: number of Edges in graph \n";
    cout<<"    p_sick: intial number of immune as a probabilty  \n";
    cout<<"    contagin: probabitly of getting sick \n";
    cout<<"    fatality: probabilty of dying after getting sick \n";       
    cout<<"    p_Immune: intial number of immune as a probabilty \n";
    cout<<"    lifetime: lifetime of disease \n";
    cout<<"    iteration: iteration step \n";
    cout<<"    n_healthy: number healthy (fraction) \n";
    cout<<"    n_sick: number sick (fraction) \n";
    cout<<"    n_Immune: number Immune (fraction) \n";
    cout<<"    n_dead: number dead (fraction) \n";
    cout<<endl;
    cout<<endl<<"DEFAULT OPTIONS: "<<endl;
    cout<<"    Graph type: ";
    if (grphtypedflt==0) cout<<"random using ran2"<<endl;
    else if (grphtypedflt==1) cout<<"triangular lattice"<<endl;
    else if (grphtypedflt==2) cout<<"triangular lattice using ran3"<<endl;
    else if (grphtypedflt==3) cout<<"fcc lattice"<<endl;
    else if (grphtypedflt==4) cout<<"cubic lattice"<<endl;
    else if (grphtypedflt==5) cout<<"square lattice"<<endl;
    else if (grphtypedflt==6) cout<<"scale-free graph"<<endl;
    else if (grphtypedflt==7) cout<<"rnd2"<<endl;
    else if (grphtypedflt==8) cout<<"barabasi network"<<endl;
}
    cout<<endl; 
    
  exit(1);
}
void pass_test(std::vector<std::string>* args){
	return;
}
void pass_test2(int& argc, char **argv){
	return;
}


void inputprocessing(int& argc, char **argv, graph& g, dPOP& pop, boolcat& b, int& ITERCASE){
   int argz = 1;
   int dummy;  
   std::stringstream ss;
   std::string inputstr, inputfile;
   b.debug_test=false;
   
   //for (int i = 1; i < argc; i++) {
   //   cout<<"i="<<i<<" "<<argv[i]<<" ";	
   //}
   
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
      else if(inputstr.compare("-sclf") == 0) g.grphtype=6;   //scale free Lattice
      else if(inputstr.compare("-staticsclf") == 0) g.grphtype=7;   //scale free Lattice      
      else if(inputstr.compare("-bb") == 0) g.grphtype=8;   //scale free Lattice    
      
 //     else if(inputstr.compare("-bb") == 0) {                 //bethe lattice?
//      	g.grphtype=8;
// 	    i++;
// 		ss>>dummy;
// 		dummy=atoi(argv[i]);
// 		int L =dummy;
// 	    g.set_L(dummy);
//	     g.set_L(atoi(argv[argz]));
//	     cout<<"nkernal: "<<g.get_L()<<" ";
// 	  }
      else if(inputstr.compare("-rndfixedc") == 0||inputstr.compare("-rf") ==0){          //random fixed connectivity
      	g.grphtype=9;
	    i++;
		ss>>dummy;
		dummy=atoi(argv[i]);
	    g.set_L(dummy);
	    cout<<"fixed maxc: "<<g.get_L()<<" ";
	  }
      else if(inputstr.compare("-rnd2") == 0) g.grphtype=10;   //random Lattice type 2	 
      else if(inputstr.compare("-CM") == 0) g.grphtype=11;   //Configuration Model      
	  else if(inputstr.compare("-HCM") == 0) g.grphtype=12;   //Configuration Model	    
      else if(inputstr.compare("-rnd") == 0) g.grphtype=0;   //random graph   
      //************************Boundary Conditions*******************************
      else if(inputstr.compare("-bc0") == 0) g.boundarycond=0;   //random graph   
      else if(inputstr.compare("-bc1") == 0) g.boundarycond=1;   //random graph 
      else if(inputstr.compare("-bc2") == 0) g.boundarycond=2;   //random graph 
      //************************Get a Seed for the Graph************************** 
      else if(inputstr.compare("-seed") == 0){          //random fixed connectivity
        i++;
		dummy=atoi(argv[i]);
	    b.lseedstart=dummy;
	    cout<<"seed: "<<b.lseedstart<<" ";
	  }
      //************************Boolean control***********************************
      else if(inputstr.compare("-w") == 0)       	b.write=true;
      else if(inputstr.compare("-dsply") == 0)     	b.dsply=true;
      else if(inputstr.compare("-fullout") == 0)   	b.fulloutput=true;
      else if(inputstr.compare("-histo") == 0)     	b.histo=true; 
      else if(inputstr.compare("-conv") == 0)     	b.conv=true;  
      else if(inputstr.compare("-debug") == 0)     	b.debug_test=true;
      else if(inputstr.compare("-h") == 0){  //help menu
		Disease_usage(false); 
		exit(1);  
	  }     	
      else if(inputstr.compare("-hv") == 0){  //help menu verbose
		Disease_usage(true); 
		exit(1);  
	  }	  
      //************************Disease control***********************************      
      else if(inputstr.compare("-I") == 0)          pop.I=true;   //Immune after sickness    
	  else if(inputstr.compare("-simple") == 0)     pop.simple=true;   //no random vectors, n_sick counts immunity/death/recovered
	  //************************Disease Iterators*********************************
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
         if (b.debug_test) cout<<"You want to read the file "<<inputfile<<std::endl;
	  }
	  else if(inputstr.compare("-seed") ==0){  //read a file in gml format
	     i++;
	     ss.clear();
	     ss>>dummy;
	     dummy=atoi(argv[i]);

         if (b.debug_test) cout<<"You want to read the file "<<inputfile<<std::endl;
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
		  if (b.debug_test) cout<<"n= "<<g.get_n()<<endl;
		 // i++;
	  }
	  else if(j==1){  //input edge probability
		  g.p=atof(argv[i]); //works
		  if (b.debug_test) cout<<"p="<<ss.str()<<" "<<g.p<<std::endl;
		  j++;
	  }
	  else if(j==2){ //probability of being initially sick
		  pop.p_sick=atof(argv[i]); //works
		  if (b.debug_test) cout<<"p_sick="<<ss.str()<<" "<<pop.p_sick<<std::endl;
		  j++;
	  }
	  else if(j==3){ //probability of getting sick
		  pop.contagin=atof(argv[i]); //works
		  if (b.debug_test) cout<<"p_contagin="<<ss.str()<<" "<<pop.contagin<<std::endl;
		  j++;
	  }
	  else if(j==4){ //probability of death!
		  pop.fatality=atof(argv[i]); //works
		  if (b.debug_test) cout<<"fatality="<<ss.str()<<" "<<pop.fatality<<std::endl;
		  j++;
	  }	  
	  else if(j==5){ //probability of initial immunity
		  pop.p_Immune=atof(argv[i]); //works
		  if (b.debug_test) cout<<"p_Immune="<<ss.str()<<" "<<pop.p_Immune<<std::endl;
		  j++;
	  }
	  else if(j==6){ //sickness duration
		  pop.lifetime=atof(argv[i]); //works
		  if (b.debug_test) cout<<"lifetime="<<ss.str()<<" "<<pop.lifetime<<std::endl;
		  j++;
	  }	  
	  //**** if bb graph give size of kernal *********************************
      if (g.grphtype==8) cout<<"nkernal: "<<g.get_L()<<" "; // dimension of lattice
      if (g.grphtype==3 || g.grphtype==4) { // dimension of lattice
         g.set_L(g.get_n());
         g.set_n(std::pow(g.get_L(),3));
      }   
   }
   
	
	return;
}
