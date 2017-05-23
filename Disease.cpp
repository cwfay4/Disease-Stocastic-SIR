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

//#include "defaults.hpp"
//#include "grphfns_a.hpp"
//#include "cvrfns_a.hpp"
//#include "Table2D.hpp"
//#include "LoPR_alg_o.hpp"
//#include "LoPR_basic_utils.hpp"

using namespace std;

#define NUM_STEPS 1000
#define grphtypedflt 0
//#define debug true

void vce_precon(graph& g, long &seed, dPOP& a){
   long idum=seed;
     
   for (unsigned int i=1;i<=1000;i++){ //initialize random number generator
      ran2(&seed);
   }   
      
   for(int i=0; i<g.get_n(); i++){
	  g.nodes[i].stateP=0;	   
      if (g.nodes[i].edges.size()==0){ //if node has no edges set to uncovered (0)
         g.nodes[i].frz=true;
	     g.nodes[i].intState=-1; //isolated node ignore it
      }
      else if (ran2(&idum)<a.p_Immune){  // healthy immune
	     g.nodes[i].intState=0;
	     g.nodes[i].frz=true;
	     a.n_immune++;
	  }
	  else if (ran2(&idum)<a.p_sick){ //sick
		 g.nodes[i].intState=1;
	     g.nodes[i].frz=false; 
	     a.n_sick++;
      }
      else {// healthy not immune
		 g.nodes[i].frz=false;
	     g.nodes[i].intState=0;
	     a.n_healthy++;
      }
      if (g.get_debug())cout<<i<<" state"<<g.nodes[i].intState<<" frz"<<g.nodes[i].frz<<endl;
   }
   if (true) cout<<" sick "<<a.n_sick<<" Healthy "<<a.n_healthy<<" Immune "<<a.n_immune<<" dead "<<a.n_dead<<endl;
//   g.calc_sum_P();
   return;   
}
//***********************************************************************
void randomize_node_list(int n, std::vector<int>& nlist, long &idum){
   std::vector<int> nlist2(nlist.size(),0);
   int j;
   double dj;
   bool debug=false;
  
   for(unsigned int i=0; i<nlist.size(); i++){
      dj=fmod(ran2(&idum)*nlist.size(),nlist.size()); 
      j=(int) dj; // chose a random node from list 2 and place it in i in the list;	   
      std::iter_swap(nlist.begin()+i,nlist.begin()+j);
   }
 //  nlist.clear();
 //  if (debug) cout<<" nlist size after clearing "<<nlist.size()<<endl;
 //  while (nlist2.size()>0){
 //     dj=fmod(ran2(&idum)*nlist2.size(),nlist2.size());
 //     j=(int) dj; // chose a random node from list 2 and place it in i in the list;
 //     nlist.push_back(nlist2[j]);
 //     nlist.pop_front();
      //std::vector<int>::iterator 
      //     p1=find(nlist.begin(), nlist.end(), nlist2[i]); 
     // auto it = std::find(nlist2.begin(), nlist2.end(), nlist2[i]);
     // if (it != nlist.end()) {
     //     using std::swap;     // swap the one to be removed with the last element
     //     swap(*it, nlist.back()); // swap the one to be removed with the last element
     //     nlist.pop_back();        // and remove the item at the end of the container
     // }                        // to prevent moving all items after 'nlist2[i]' by one
     //nlist2.erase(nlist2[j]); 
  //   nlist2.erase(nlist2.begin()+(j-1)); 
     // nlist2.erase(find(nlist2.begin(), nlist2.end(), nlist[i])); //not right synatatically
//   }
   if (debug) cout<<" nlist size after randomizing "<<nlist.size()<<endl;
   nlist2.clear();
   return;
}
void disease_pop(graph& g, dPOP a, long seed2){
	bool conv=true;
	long idum=seed2;
	
	double dn = (double) g.get_n();
	
	std::vector<int> nlist(g.get_n(),0);
	for(unsigned int i=0; i<nlist.size(); i++){
      nlist[i]=i;
    }
    
	randomize_node_list(g.get_n(), nlist, seed2);
	
	//if (g.get_debug()) {
	//	cout<<a.n_sick<<" "<<a.n_healthy<<" "<<a.n_immune;
	//    cout<<"node list randomized ";
	// 	for(unsigned int i=0; i<nlist.size(); i++){
    //      cout<<" nlist "<<nlist[i]<<" ";
    //    } 
    //    cout<<endl;  
	// }
    
    for (unsigned it=0; it<NUM_STEPS; it++){ //iterate over trials
		randomize_node_list(g.get_n(), nlist, seed2);
		for(unsigned i=0; i<nlist.size(); i++){ //for all nodes
			//cout<<i<<" state"<<g.nodes[i].intState<<" frz"<<g.nodes[i].frz<<endl;
			if (!g.nodes[i].frz && g.nodes[i].intState==1){ //if you are sick but not dead
			//	cout<<" you are sick but not dead"<<endl;
			   if (g.nodes[i].stateP > a.lifetime){
				   if(ran3(&idum)<a.fatality){ //oops you died
					  g.nodes[i].frz =true;
					  a.n_dead++;
					  a.n_sick--;
				//	  cout<<" you died!"<<endl;
				   }
				   else {  //bam! you are healed...and immune.
					  g.nodes[i].intState=0;
					  if (a.I) {
	                    g.nodes[i].frz=true; 
	                    a.n_immune++;
	                  }
	                  else a.n_healthy++;
	                  a.n_sick--;
	                  g.nodes[i].stateP=0;
	               //   cout<<" you are healed!"<<endl;
			       }
			   }
			   else g.nodes[i].stateP++; //time counter for sickness
			}			
			else if (!g.nodes[i].frz && g.nodes[i].intState==0){ //if you are healthy, you might get sick.
			  // cout<<" you are healthy"<<endl;
			   double n_I=0;
			   for (int j=0; j < g.nodes[i].edges.size(); j++){//for all edges connedtd to i
	               int n2=g.nodes[i].edges[j];
	               if (g.nodes[n2].intState==1 && !g.nodes[n2].frz) n_I++; //a node that is frozen and sick is dead!
	           }
	           double prob = 1.0 - pow((1.0-a.contagin),n_I); //probability of getting sick
	           //cout<<prob<<" ";
	           if (ran3(&idum)<prob){ //Are you sick
                  g.nodes[i].intState=1; //apparently yes.
	              g.nodes[i].frz=false;
	              g.nodes[i].stateP=0;
	              a.n_sick++;
	              a.n_healthy--;
	              //cout<<" You are now sick"<<endl;
		       }
	           else if (ran2(&idum)<a.p_sick){  //areyou randomly getting sick?
				  g.nodes[i].intState=1; //apparently yes.
	              g.nodes[i].frz=false;
	              g.nodes[i].stateP=0;
	              a.n_sick++;
	              a.n_healthy--; 
			   }
		    }
		}
	  
      if (conv){ 
	     std::ofstream output("Disease_Conv.csv", ios_base::app);
         output<<it<<" "<<(double)a.n_healthy/dn<<" "<<(double)a.n_sick/dn<<" "<<(double)a.n_immune/dn<<" "<<(double)a.n_dead/dn<<endl;
	     output.close();
      }
      if (g.get_debug()) cout<<a.n_healthy<<" "<<a.n_sick<<" "<<a.n_immune<<endl;
   }
   cout<<" sick "<<a.n_sick<<" Healthy "<<a.n_healthy<<" Immune "<<a.n_immune<<" dead "<<a.n_dead<<endl;
   cout<<" sick "<<(double)a.n_sick/dn<<" Healthy "<<(double)a.n_healthy/dn<<" Immune "<<(double)a.n_immune/dn<<" dead "<<(double)a.n_dead<<endl;
   return;	
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
   int dummy;

   graph g;
   g.grphtype=grphtypedflt;
   
   std::stringstream ss;
   std::string inputstr, inputfile;
   
   int i, j=0;
   for (i = 1; i < argc; i++) {
   	  ss.clear();
	  ss.str(argv[i]);
      ss>>inputstr; 
      if(inputstr.compare("-pt") == 0) g.grphtype=1;         //planar triangular graph
      else if(inputstr.compare("-ptr3") == 0) g.grphtype=2;   //planar triangular graph with ran3
      else if(inputstr.compare("-fcc") == 0) g.grphtype=3;   //FCC Lattice
      else if(inputstr.compare("-cubic") == 0) g.grphtype=4;   //cubic Lattice
      else if(inputstr.compare("-sq") == 0) g.grphtype=5;   //squareLattice     
      else if(inputstr.compare("-sclf") == 0) g.grphtype=6;   //scale free Lattice
      else if(inputstr.compare("-rnd2") == 0) g.grphtype=7;   //random Lattice type 2
      else if(inputstr.compare("-bb") == 0) {                 //bethe lattice?
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
      else if(inputstr.compare("-I") == 0)          pop.I=true;   //Immune after sickness     
      else if(inputstr.compare("-w") == 0)       	write=true;
      else if(inputstr.compare("-dsply") == 0)     	dsply=true;
      else if(inputstr.compare("-fullout") == 0)   	fulloutput=true;
      else if(inputstr.compare("-histo") == 0)     	histo=true; 
      else if(inputstr.compare("-conv") == 0)     	conv=true;      
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
		  if (debug_test) cout<<"p_Immune="<<ss.str()<<" "<<pop.p_Immune<<std::endl;
		  j++;
	  }
	  else if(j==6){ //sickness duration
		  pop.lifetime=atof(argv[i]); //works
		  if (debug_test) cout<<"lifetime="<<ss.str()<<" "<<pop.lifetime<<std::endl;
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
   
      if (g.grphtype==8) cout<<"nkernal: "<<g.get_L()<<" "; // dimension of lattice
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
   
   cout<<"n= "<<g.get_n()<<" c= "<<g.get_c()<<" nE= "<<g.get_nE()<<" p_sick= "<<pop.p_sick;
   cout<<" contagin= "<<pop.contagin<<" fatality="<<pop.fatality<<" p_Immune="<<pop.p_Immune;
   cout<<" lifetime= "<<pop.lifetime<<endl;
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
   cout<<" contagin= "<<pop.contagin<<" fatality="<<pop.fatality<<" p_Immune="<<pop.p_Immune<<endl;
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
   
   return 0;
}
