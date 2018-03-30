/** vcELoPR_alg                                                      ***/
/**                                                                  ***/
/** C. Fay Jan 2003                                                  ***/
/** Version 2.1.0.0   05.10.2004                                     ***/
/** Run: crprcfnsgraph gml_filename                                  ***/
/** Output: dat files that can be opened by xmgrace as a visual      ***/
/**        representation of the graph                               ***/
/***********************************************************************/

#include <cstdio>                                     
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <vector>
#include <list>
#include <string>

#include "random.hpp"
#include "graph.hpp"

#include "defaults.hpp"
//#include "grphfns_a.hpp"
//#include "cvrfns_a.hpp"
//#include "Table2D.hpp"
#include "LoPR_alg_o.hpp"
#include "LoPR_basic_utils.hpp"

using namespace std;
//**************************************************************************
// makes histogram of P[node] occurance in G
//void vce_histogram(graph g, double numbox, std::string ss){ 
//      std::vector<int> histogram((int)numbox+1,0);
//      double dj, boxsize=1.0/(double) numbox;
 //     int j;
      //std::string name (s);
 //     name=ss.data()+"P_histo.csv";
//      
//      for(int i=0; i<g.get_n();i++){
//         dj=g.node[i].stateP/boxsize;
//	     j = (int) dj;
	 //cout<<P[i]<<" "<<j<<endl;
//	     while(g.node[i].stateP>(double)j*boxsize) j++;
//	     if (j<=(int)numbox) histogram[j]++;
//	     else cout<<"error P[i] too large \n";      
 //     }
 //     std::ofstream output(name.data());
 //     for(int i=0; i<=(int) numbox; i++){  //in the output 1 = in MIS, 0 = in VC
//         double di=(double)i*boxsize;
 //        output<<1-di<<" "<<histogram[i]<<endl;
	 //cout<<1-di<<" "<<histogram[i]<<endl;
//      }
//      output.close();
//      return;
//}
//**************************************************************************
void vce_precon(graph g, long &seed){
      long idum=seed;
      int n2, nstop, nstart, t;
      bool frzlf=false;
     
       for (unsigned int i=1;i<=1000;i++){ //initialize random number generator
        ran2(&seed);
       }   
      
   for(int i=0; i<g.get_n(); i++){
      if (g.nodes[i].edges.size()==0){ //if node has no edges set to uncovered (0)
         g.nodes[i].frz=true;
         //P[i]=ran2(&idum);
	     g.nodes[i].stateP=0;
      }
      else if (g.nodes[i].edges.size()==1 && frzlf){  // if node has one edge it is a leaf
	     g.nodes[i].stateP=0;
	     int j = g.nodes[i].edges.back(); //int j = g.node[i].edges.first();
	     g.nodes[j].stateP=1;
	     g.nodes[i].frz=true;
	     g.nodes[j].frz=true;
	  }
      else {//if (!g.nodes[i].frz){
		 g.nodes[i].frz=false;
         g.nodes[i].stateP=ran2(&idum);
      }
   }
   g.calc_sum_P();
   return;   
}
//***********************************************************************
void vce_randomize_node_list(int n, std::vector<int>& nlist, long &idum){
   std::vector<int> nlist2(nlist.size(),0);
   int j;
   double dj;
   
   for(unsigned int i=0; i<nlist.size(); i++){
      nlist2[i]=i;
   }
   
   for(int i=0; i<nlist.size(); i++){
      dj=fmod(ran2(&idum)*nlist2.size(),nlist2.size());
      j=(int) dj; // chose a random node from list 2 and place it in i in the list;
      nlist[i]=nlist2[j];
      //std::vector<int>::iterator 
      //     p1=find(nlist.begin(), nlist.end(), nlist2[i]); 
      auto it = std::find(nlist2.begin(), nlist2.end(), nlist2[i]);
      if (it != nlist.end()) {
          using std::swap;     // swap the one to be removed with the last element
          swap(*it, nlist.back()); // swap the one to be removed with the last element
          nlist.pop_back();        // and remove the item at the end of the container
      }                        // to prevent moving all items after 'nlist2[i]' by one

     // nlist2.erase(find(nlist2.begin(), nlist2.end(), nlist[i])); //not right synatatically
   }
   
   nlist2.clear();
   return;
}
//***********************************************************************
//***********************************************************************
void vce_siteELoPRalg(graph g, long &idum, bool conv, bool rlists){
    double *P_old, *dP;

    std::vector< int >  nlist(g.get_n(),0); //list of nodes that will be randomized
    double P_prod=1, delta_ave_P=100, delta_P_ave=100, Pave=1, dPsum=0;
    double Pave_old=0;
    int n1, n2;

    int l=0;
    
    int n = g.get_n();
    double dn = double (n);
    vce_randomize_node_list(n, nlist, idum);
    
    while (delta_P_ave>BND || delta_ave_P>BND){
       dPsum=0;
    
       for(unsigned i=0; i<nlist.size(); i++){
            n1=nlist[i];
	    //cout<<n1<<endl;
	        P_old[n1]=g.nodes[n1].stateP;
            if(!g.nodes[n1].frz){ //if the node is not frozen find a new stateP
	           P_prod=1;
               for (int j=0; j < g.nodes[n1].edges.size(); j++){//for all edges connedtd to n1
	               n2=g.nodes[n1].edges[j];
	               if (n2!=0) P_prod=P_prod*g.nodes[n2].stateP;
	           } 
               g.nodes[n1].stateP=1-P_prod;
	        }
	        dPsum+=fabs(g.nodes[n1].stateP-P_old[n1]);
       }
       Pave_old=g.get_Psum()/dn;
       g.calc_sum_P();
       Pave=g.get_Psum()/dn;
       delta_ave_P=fabs(Pave-Pave_old);
       delta_P_ave=dPsum/dn;
       if (conv){  //print convergence as a function of loop number
         l++;
	     std::ofstream output("sConv.csv", ios_base::app);
         output<<l<<" "<<Pave<<" "<<delta_ave_P<<" "<<delta_P_ave<<endl;
	     output.close();
       }
       if (rlists) vce_randomize_node_list(n, nlist, idum); 
    }
    g.calc_sum_P_F();
        
   delete[] P_old;
    return;
}

//***********************************************************************
double vce_calc_PI_b(graph g, int n1, int n2){ 
   // given n1, n2 defining the edge (n1,n2) find X1 for that edge
   int n3=0;
   double xl=1, xk=1, S=0;
   int nstart, nstop, t;
   
   for (int j=0; j < g.nodes[n1].edges.size(); j++){//for all edges connected to n1
	  n3=g.nodes[n1].edges[j];
      if (n3!=n2 && n3!=0) xl=xl*g.nodes[n3].stateP;
   }
   
   for (int j=0; j < g.nodes[n2].edges.size(); j++){//for all edges connected to n2
	  n3=g.nodes[n2].edges[j];
      if (n3!=n1 && n3!=0) xk=xk*g.nodes[n3].stateP;
   }
   S=xl-(0.5*xk*xl);
   return S;   
}
void vce_bELoPRalg(graph g, long &idum, bool conv, bool rlists){
   double S=0, X=0, delta_ave_P=100, Pave=1;
   double Pave_old, dPsum=0, delta_P_ave;
   double dP;
   std::vector <double> P_old(g.get_n(),0);    
   for(unsigned int i=0;i<=SIZE;i++){
       P_old[i]=0;
    }
   std::vector< int >  nlist(g.get_n(),0);
   int n1, n2, nstart, nstop, t;

   int l=0;
   int n = g.get_n();
   double dn = double (n);
   if(!rlists) vce_randomize_node_list(n, nlist, idum);  
   
   while(delta_P_ave>BND || delta_ave_P>BND){
      dPsum=0;
      if(rlists) vce_randomize_node_list(n, nlist, idum);
      for(unsigned i=0; i<nlist.size(); i++){ //for all nodes
         n1=nlist[i];
	     P_old[n1]=g.nodes[n1].stateP; //store previous state
	     //if (g.nodes[n1].edges.size()==0) g.nodes[n1].stateP=0; //if the node is isolated
	     //else {
	     if(!g.nodes[n1].frz){
	        S=0;
            for (int j=0; j < g.nodes[n1].edges.size(); j++){//for all edges connedtd to n1
	           n2=g.nodes[n1].edges[j];
               X=vce_calc_PI_b(g, n1, n2);
		       S+=X;
	        }
	        g.nodes[n1].stateP=1.0-S/(double)g.nodes[n1].edges.size();
	     }	    
	     dPsum+=fabs(g.nodes[n1].stateP-P_old[n1-1]);
      }
      Pave_old=g.get_Psum()/dn; //average P_sum from previous run
      g.calc_sum_P(); 
      Pave=g.get_Psum()/dn;   //average P_sum from current run
      delta_ave_P=fabs(Pave-Pave_old);  //change in average P_sum
      delta_P_ave=dPsum/dn;  //average of site changes
      if (conv){
         l++;
         std::ofstream output("bConv.csv", ios_base::app);
         output<<l<<" "<<Pave<<" "<<delta_ave_P<<" "<<delta_P_ave<<endl;
         output.close();
      }
   }
   g.calc_sum_P_F();
   return;
}
void vce_bELoPRalg2(graph g, long &idum){
   double S=0, X=0, delta_ave_P=100, Pave=0, Pave_old;
   int n1, n2;
   int j;
   double dj;
   int n = g.get_n();
   std::vector <double> P_old(n,0);  
   double dn = double (n);
   
   std::vector< int >  nlist(n,0);
   std::vector<int> nlist2, nlist3(n,0);// nlist4(n,0);
   for(int i=0; i<n; i++){
      nlist[i]=i;
   }
//nlist.assign(nlist4.begin(),nlist4.end());
   while(delta_ave_P>BND){
      nlist2.assign(nlist.begin(),nlist.end());
      for(unsigned i=0; i<nlist.size(); i++){ //for all nodes in list nlist
	     P_old[i]=g.nodes[i].stateP;
         dj=fmod(ran2(&idum)*nlist2.size(),nlist2.size()); //choose a node 2
         j=(int) dj;
         nlist3[i]=nlist2[j]; //nlist3 is a randomized list.
         nlist2.erase(find(nlist2.begin(), nlist2.end(), nlist3[i]));
         n1=nlist[i];
         if(!g.nodes[n1].frz){
	        S=0;
            for (int j=0; j < g.nodes[n1].edges.size(); j++){//for all edges connedtd to n1
	           n2=g.nodes[n1].edges[j];
               X=vce_calc_PI_b(g, n1, n2);
		       S+=X;
	        }
	        g.nodes[n1].stateP=1.0-S/(double)g.nodes[n1].edges.size();
	     }	
	  }
      nlist.assign(nlist3.begin(),nlist3.end()); //list 1 is now the randomized list
      Pave_old=g.get_Psum()/dn;
      g.calc_sum_P();
      Pave=g.get_Psum()/dn;
      delta_ave_P=fabs(Pave-Pave_old);
      //delta_P_ave=dPsum/dn;
      //Pave=vce_sum_P(P,n)/(double) n;
      //delta_ave_P=fabs(Pave-Pave_old);
   }
   g.calc_sum_P_F();
   return;
}
//***********************************************************************
void vce_listdegen(graph g, std::list< int > &degen){
/* lists the number of degenerate nodes */

   degen.clear();
   for(int i=0; i<g.get_n(); i++){
      if (g.nodes[i].stateP !=0 && g.nodes[i].stateP != 1){
         degen.push_back(i);
      }
   }
   return;
}
//***********************************************************************
void vce_siteDIG_weighted(graph g, long &seed, bool rlists, bool weight){  //weighted
/* uncovered 1 degenerate node */
   std::list< int >  degen_nodes;
   int node1;
   bool conv=convdflt;
   
   vce_listdegen(g, degen_nodes);
   while (degen_nodes.size()>0){
	 //choose the first node with stateP<0.2 and and uncover it.
	 if (weight){  //weight switched by the variable weight
        while(g.nodes[degen_nodes.front()].stateP>0.2 && degen_nodes.size()>1){ 
           degen_nodes.pop_front();
        }
     }
     node1=degen_nodes.front(); 
     g.nodes[node1].stateP=0; //uncover node 1
     //degen_nodes.clear();
     vce_siteELoPRalg(g, seed, conv, rlists); //rerun siteEloPR
     vce_listdegen(g, degen_nodes);
     }
   g.calc_sum_P_F();
   degen_nodes.clear();
   return;
}
void vce_siteDIG(graph g, long &seed, bool rlists){
/* uncovered 1 degenerate node */
   std::list< int >  degen_nodes;
   int node1;
   bool conv=convdflt;
   
   vce_listdegen(g, degen_nodes);
   while (degen_nodes.size()>0){
     node1=degen_nodes.front();
     g.nodes[node1].stateP=0;
     // frz[node1]=true;
     // degen_nodes.clear();
     vce_siteELoPRalg(g, seed, conv, rlists); //rerun siteEloPR
     vce_listdegen(g, degen_nodes);
   }
   g.calc_sum_P_F();
   degen_nodes.clear();   
   return;
}
//***********************************************************************
//***********************************************************************
//***********************************************************************
//***********************************************************************
//***********************************************************************
//***** control subroutines *********************************************
//***********************************************************************
void vcd_random_number_prep(std::vector< long > &seed){
   //********************************** Setting up random numbers ********
   for (unsigned int z=1;z<=1000;z++){//initialize random number generator
       float fseed=ran2(&seed[0]);
       }
  // std::vector< long > seed(3,0); 
   for (unsigned int z=1;z<=3;z++){
       seed[z]=(-1*(int) (ran2(&seed[0])*1e8));
       }
   //long lseed1 = seed[0];
   //long lseed2 = seed[1];
   //long lseed3 = seed[2]; 
   
   return;	
}
//***********************************************************************
void vce_sb(graph g, bool write, bool dsply, bool fulloutput,
      std::string &dname, std::vector<double> &datum, long &lseedin, bool histo, 
      bool conv, bool rlists, char* filename){
   int nE=0;
   double Psum=0, Fsum=0; 
   std::vector< long > seed(4,0); 
   seed[0]=lseedin;
  // long lseed=lseedin; 
   //********************************** Setting up random numbers ********
   //for (unsigned int z=1;z<=1000;z++){//initialize random number generator
   //    float fseed=ran2(&lseed);
   //    }
   //std::vector< long > seed(3,0); 
   //for (unsigned int z=1;z<=3;z++){
    //   seed[z-1]=(-1*(int) (ran2(&lseed)*1e8));
    //   }
    
    
   long lseed1 = seed[1];
   long lseed2 = seed[2];
   long lseed3 = seed[3]; 
   //*********************************************************************
   g.seed=lseed1;
   g.build_graph(write);
   if (dsply) g.dsplygrph();
   
   double c=g.get_c();
   if (fulloutput) {
      std::ofstream dout(dname.data(), ios_base::app);
      dout<<g.get_n()<<" "<<g.p<<" "<<c<<" "<<g.get_nE()<<" "<<lseedin;
      dout.close();        
   }
   
   //saves data
   datum[0]=(double)g.get_n();
   datum[1]=g.p;
   datum[2]=g.c;
   datum[3]=(double)g.get_nE();
         
   //double *V; //V is the node states not needed.
   //V=new double[SIZE+1];
   //bool *frz;
   //frz=new bool[SIZE+1];   
   //for (int i=0; i<=SIZE; i++){
   //   V[i]=0
   //   frz[i]=false;
   //}
   
   vce_precon(g, lseed2);
   //vce_precon(gv, g2, gsize, n, V, frz, lseed2);
   
   vce_siteELoPRalg(g, lseed3, conv, rlists); //rerun siteEloPR
   //g.graph_calc_sum_P_F(); //now in the alg's
   Psum=g.get_Psum();
   Fsum=g.get_Fsum();   
   
   if (histo) g.StateHistogram(100, "s");
   if (fulloutput) {
      std::ofstream dout(dname.data(), ios_base::app);
      dout<<" sELoPR: "<<Psum/(double) g.get_n();
      dout.close();
   }
   datum[4]=(double)Psum/(double) g.get_n();
   datum[5]=(double)Fsum/(double) g.get_n();      
   
   vce_precon(g, lseed2);  //precondition for a independent bELoPR run
   vce_bELoPRalg(g, lseed1, conv, rlists);
   //g.graph_calc_sum_P_F(); //now in the alg's
   Psum=g.get_Psum();
   Fsum=g.get_Fsum();  
   if (histo) g.StateHistogram(100, "s");
   if (fulloutput) {
      std::ofstream dout(dname.data(), ios_base::app);
      dout<<" bELoPR: "<<Psum/dn<<endl;
      dout.close();
   }
   datum[6]=(double)Psum/(double) g.get_n();
   datum[7]=(double)Fsum/(double) g.get_n(); 
   
   return;
}

//***********************************************************************
void vce_sDig(graph g, bool write, bool dsply, bool fulloutput, bool rlists,
      std::string &dname, std::vector<double> &datum, 
      long &lseedin, bool histo, bool conv, char* filename){

   int nE=0;
   double Psum=0, Fsum=0; 
   
   std::vector< long > seed(4,0);  //seed prep
   seed[0]=lseedin;
   vcd_random_number_prep(seed);
   long lseed1 = seed[1];
   long lseed2 = seed[2];
   long lseed3 = seed[3]; 
   
   g.seed=lseed1;
   g.build_graph(write);
   if (dsply) g.dsplygrph();

   c=g.get_c(); //calculates c
   if (fulloutput) {
      std::ofstream dout(dname.data(), ios_base::app);
      dout<<g.get_n()<<" "<<g.p<<" "<<g.c<<" "<<g.get_nE()<<" "<<lseedin; 
     dout.close();       
   }
   datum[0]=(double)g.get_n();
   datum[1]=g.p;
   datum[2]=g.c;
   datum[3]=(double)g.get_nE();
          
   vce_precon(g, lseed2);
   vce_siteELoPRalg(g, lseed3, conv, rlists);
   //g.graph_calc_sum_P_F();
   Psum = g.get_Psum();
   Fsum = g.get_Fsum();        
  // if (dsply){
  //    std::ofstream output1("dsplyg.csv");
   //   grphfns_dsplygrph_wstate(n, gv, g2, V, output1, "s");
   //   output1.close();
   //}
   if (histo) g.StateHistogram(100, "s");
   if (fulloutput){
      std::ofstream dout(dname.data(), ios_base::app);
      dout<<" sELoPR: "<<Psum/n<<" FrozFrac: "<<Fsum/dn; 
      dout.close();   
   }
   datum[4]=(double)Psum/(double)g.get_n();
   datum[5]=(double)Fsum/(double)g.get_n();  
   
   vce_siteDIG(g, frz, lseed1, rlists);
//   g.graph_calc_sum_P_F();
   Psum = g.get_Psum();
   Fsum = g.get_Fsum();       
   //if (dsply){
   //   std::ofstream output1("dsplyg.csv");
   //   grphfns_dsplygrph_wstate(n, gv, g2, V, output1, "sd");
   //   output1.close();
   //}
   if (histo) g.StateHistogram(100, "s");
   if (fulloutput) {
      std::ofstream dout(dname.data(), ios_base::app);
      dout<<" sDig: "<<Psum/dn<<" FrozFrac: "<<Fsum/dn<<endl;
      dout.close();   
   }
   datum[6]=(double)Psum/(double)g.get_n();
   datum[7]=(double)Fsum/(double)g.get_n();
      

   return;     
}
//***********************************************************************
void vce_sDbD(graph g, bool write, bool dsply, bool fulloutput, bool rlists,
      std::string &dname, std::vector<double> &datum, 
      long &lseedin, bool histo, bool conv, char* filename){
   int nE=0;
   double Psum=0, Fsum=0; 
   
   std::vector< long > seed(4,0);  //seed prep
   seed[0]=lseedin;
   vcd_random_number_prep(seed);
   long lseed1 = seed[1];
   long lseed2 = seed[2];
   long lseed3 = seed[3]; 
   
   g.seed=lseed1;
   g.build_graph(write);
   if (dsply) g.dsplygrph();

   double c=g.get_c(); //calculates c
   
   if (fulloutput) {
      std::ofstream dout(dname.data(), ios_base::app);
      dout<<g.get_n()<<" "<<g.p<<" "<<c<<" "<<g_get_nE()<<" "<<lseedin;
      dout.close();        
   }
   datum[0]=(double)g.get_n();
   datum[1]=g.p;
   datum[2]=g.c;
   datum[3]=(double)g.get_nE();
          
   vce_precon(g, lseed2);
   vce_siteELoPRalg(g, lseed3, conv, rlists);
   Psum = g.get_Psum();
   Fsum = g.get_Fsum(); 
   
//   if (dsply){
//      std::ofstream output1("dsplyg.csv");
//      grphfns_dsplygrph_wstate(n, gv, g2, V, output1, "s");
//      output1.close();
//   }
   if (histo) g.StateHistogram(100, "s");
   if (fulloutput) {
      std::ofstream dout(dname.data(), ios_base::app);
      dout<<" sELoPR: "<<Psum/n<<" FrozFrac: "<<Fsum/n;
      dout.close();
   }
   datum[4]=(double)Psum/(double)n;
   datum[5]=(double)Fsum/(double)n;      
      
   vce_siteDIG(g, lseed1, rlists);
   Psum = g.get_Psum();
   Fsum = g.get_Fsum();     
 //  if (dsply){
 //     std::ofstream output1("dsplyg.csv");
 //     grphfns_dsplygrph_wstate(n, gv, g2, V, output1, "sd");
 //     output1.close();
//   }
   if (fulloutput) {
      std::ofstream dout(dname.data(), ios_base::app);
      dout<<" sDig: "<<Psum/dn<<" FrozFrac: "<<Fsum/dn;
      dout.close();   
   }
   datum[6]=(double)Psum/(double)g.get_n();
   datum[7]=(double)Fsum/(double)g.get_n();
   
   vce_precon(g, lseed2);  //precondition for a independent bELoPR run
   vce_bELoPRalg(g, lseed1, conv, rlists);
   //g.graph_calc_sum_P_F(); //now in the alg's
   Psum=g.get_Psum();
   Fsum=g.get_Fsum(); 
   
//   if (dsply){
 //     std::ofstream output1("dsplyg.csv");
 //     grphfns_dsplygrph_wstate(n, gv, g2, V, output1, "b");
 //     output1.close();
 //  }
   if (histo) g.StateHistogram(100, "s");
   if (fulloutput) {
      std::ofstream dout(dname.data(), ios_base::app);
      dout<<" bELoPR: "<<Psum/dn<<" FrozFrac: "<<Fsum/dn;
      dout.close();   
   }
   datum[8]=(double)Psum/(double)g.get_n();
   datum[9]=(double)Fsum/(double)g.get_n();      
      
   vce_siteDIG(g, lseed1, rlists);
   Psum = g.get_Psum();
   Fsum = g.get_Fsum();      
//   if (dsply){
//      std::ofstream output1("dsplyg.csv");
//      grphfns_dsplygrph_wstate(n, gv, g2, V, output1, "bd");
//      output1.close();
//   }
   if (fulloutput) {
      std::ofstream dout(dname.data(), ios_base::app);
      dout<<" b_sDig: "<<Psum/dn<<" FrozFrac: "<<Fsum/dn<<endl;
      dout.close();   
   }
   datum[10]=(double)Psum/(double)g.get_n();
   datum[11]=(double)Fsum/(double)g.get_n();
      
   seed.clear();   
   return;     
}

/*****************************************************************************/
void vce_sDbD_1list_many(graph g, bool write, bool dsply, bool fulloutput,
      std::string &dname, std::vector<double> &datum, 
      long &lseedin, bool histo, bool conv, char* filename){
   int nE=0;
   double Psum=0, Fsum=0; 
   
   std::vector< long > seed(4,0);  //seed prep
   seed[0]=lseedin;
   vcd_random_number_prep(seed);
   long lseed1 = seed[1];
   long lseed2 = seed[2];
   long lseed3 = seed[3]; 
   
   g.seed=lseed1;
   g.build_graph(write);
   if (dsply) g.displaygraph();

   c=g.get_c(); //calculates c
   
 //     if (dsply){
 //        std::ofstream output1("dsplyg.csv");
 //        grphfns_dsplygrph(n, gv, g2, output1);
 //        output1.close();
 //     }
       if (fulloutput) {
         std::ofstream dout(dname.data(), ios_base::app);
         dout<<n<<" "<<p<<" "<<c<<" "<<nE<<" "<<lseedin;
         dout.close();        
      }
   datum[0]=(double)g.get_n();
   datum[1]=g.p;
   datum[2]=g.c;
   datum[3]=(double)g.get_nE();

                      // 1 list sLoPR
          
   vce_precon(g, lseed2);
   vce_siteELoPRalg(g, lseed3, conv, false);
   Psum = g.get_Psum();
   Fsum = g.get_Fsum(); 
//      if (dsply){
 //        std::ofstream output1("dsplyg.csv");
 //        grphfns_dsplygrph_wstate(n, gv, g2, V, output1, "s");
//         output1.close();
//      }
   if (histo) g.StateHistogram(100, "s");
      if (fulloutput) {
         std::ofstream dout(dname.data(), ios_base::app);
         dout<<" sELoPR_1: "<<Psum/n<<" FrozFrac_1: "<<Fsum/n;
         dout.close();
      }
      datum[4]=(double)Psum/(double)g.get_n();
      datum[5]=(double)Fsum/(double)g.get_n();      
      
   vce_siteDIG(g, frz, lseed1, false);
   Psum = g.get_Psum();
   Fsum = g.get_Fsum();  
//      if (dsply){
//         std::ofstream output1("dsplyg.csv");
//         grphfns_dsplygrph_wstate(n, gv, g2, V, output1, "sd");
//         output1.close();
//      }
      if (fulloutput) {
         std::ofstream dout(dname.data(), ios_base::app);
         dout<<" sDig_1: "<<Psum/n<<" FrozFrac_1: "<<Fsum/n;
         dout.close();   
      }
      datum[6]=(double)Psum/(double)g.get_n();
      datum[7]=(double)Fsum/(double)g.get_n();
   
                              // 1 list bLoPR
   
   vce_precon(g, lseed2);  //precondition for a independent bELoPR run
   vce_bELoPRalg(g, lseed1, conv, false);
   //g.graph_calc_sum_P_F(); //now in the alg's
   Psum=g.get_Psum();
   Fsum=g.get_Fsum(); 
   
//   if (dsply){
 //     std::ofstream output1("dsplyg.csv");
//      grphfns_dsplygrph_wstate(n, gv, g2, V, output1, "b");
//      output1.close();
//   }
   if (histo) g.StateHistogram(100, "b");
   if (fulloutput) {
      std::ofstream dout(dname.data(), ios_base::app);
      dout<<" bELoPR_1: "<<Psum/n<<" FrozFrac_1: "<<Fsum/n;
      dout.close();   
   }
   datum[8]=(double)Psum/(double)g.get_n();
   datum[9]=(double)Fsum/(double)g.get_n();      
      
   vce_siteDIG(g, frz, lseed1, false);
   Psum = g.get_Psum();
   Fsum = g.get_Fsum();    
//  if (dsply){
//      std::ofstream output1("dsplyg.csv");
//      grphfns_dsplygrph_wstate(n, gv, g2, V, output1, "bd");
//     output1.close();
//  }
   if (fulloutput) {
      std::ofstream dout(dname.data(), ios_base::app);
      dout<<" b_sDig_1: "<<Psum/n<<" FrozFrac_1: "<<Fsum/n<<endl;
      dout.close();   
   }
   datum[10]=(double)Psum/(double)g.get_n();
   datum[11]=(double)Fsum/(double)g.get_n();
   
                             // sLoPR (standard)
   vce_precon(g, lseed2);
   vce_siteELoPRalg(g, lseed3, conv, true);
   Psum = g.get_Psum();
   Fsum = g.get_Fsum(); 
//   if (dsply){
//     std::ofstream output1("dsplyg.csv");
//     grphfns_dsplygrph_wstate(n, gv, g2, V, output1, "s");
//     output1.close();
//   }
   if (histo) g.StateHistogram(100, "s");
   if (fulloutput) {
      std::ofstream dout(dname.data(), ios_base::app);
      dout<<" sELoPR: "<<Psum/n<<" FrozFrac: "<<Fsum/n;
      dout.close();
   }
   datum[12]=(double)Psum/(double)g.get_n();
   datum[13]=(double)Fsum/(double)g.get_n();      
      
   vce_siteDIG(g, frz, lseed1, true);
   Psum = g.get_Psum();
   Fsum = g.get_Fsum();     
 //  if (dsply){
 //     std::ofstream output1("dsplyg.csv");
 //     grphfns_dsplygrph_wstate(n, gv, g2, V, output1, "sd");
 //     output1.close();
 //  }
   if (fulloutput) {
      std::ofstream dout(dname.data(), ios_base::app);
      dout<<" sDig: "<<Psum/n<<" FrozFrac: "<<Fsum/n;
      dout.close();   
   }
   datum[14]=(double)Psum/(double)g.get_n();
   datum[15]=(double)Fsum/(double)g.get_n();
   
                             // bLoPR standard
   
   vce_precon(g, lseed2);  //precondition for a independent bELoPR run
   vce_bELoPRalg(g, lseed1, conv, true);
   //g.graph_calc_sum_P_F(); //now in the alg's
   Psum=g.get_Psum();
   Fsum=g.get_Fsum(); 
   
//   if (dsply){
//     std::ofstream output1("dsplyg.csv");
//      grphfns_dsplygrph_wstate(n, gv, g2, V, output1, "b");
//      output1.close();
//   }
   if (histo) g.StateHistogram(100, "b");
   if (fulloutput) {
      std::ofstream dout(dname.data(), ios_base::app);
      dout<<" bELoPR: "<<Psum/n<<" FrozFrac: "<<Fsum/n;
      dout.close();   
   }
   datum[16]=(double)Psum/(double)g.get_n();
   datum[17]=(double)Fsum/(double)g.get_n();      
      
      
   vce_siteDIG(g, frz, lseed1, true);
   Psum = g.get_Psum();
   Fsum = g.get_Fsum();   
//   if (dsply){
//      std::ofstream output1("dsplyg.csv");
//      grphfns_dsplygrph_wstate(n, gv, g2, V, output1, "bd");
//      output1.close();
//   }
   if (fulloutput) {
      std::ofstream dout(dname.data(), ios_base::app);
      dout<<" b_sDig: "<<Psum/n<<" FrozFrac: "<<Fsum/n<<endl;
      dout.close();   
   }
   datum[18]=(double)Psum/(double)g.get_n();
   datum[19]=(double)Fsum/(double)g.get_n();
  
//   delete[] V;
   delete[] seed;   
   return;     
}
//***********************************************************************
void vce_sDig_best(graph g,bool dsply, bool write, bool fulloutput, bool rlists,
       long &lseedin, bool conv, std::string &dname, char* filename){
   unsigned NUMBER=10;
   double Psum=0, Fsum=0, Pbest=0.66, Pave=0; 
   long lseed=lseedin; 
   int nE=0;
   
   std::vector< long > seed(4,0);  //seed prep
   seed[0]=lseedin;
   vcd_random_number_prep(seed);
   long lseed1 = seed[1];
   long lseed2 = seed[2];
   long lseed3 = seed[3]; 
   
   g.seed=lseed1;
   g.build_graph(write);
   c=g.get_c(); //calculates 
   if (dsply) g.displaygraph();

   if (fulloutput) {
      std::ofstream dout(dname.data(), ios_base::app);
      dout<<n<<" "<<p<<" "<<c<<" "<<nE<<" "<<lseedin;
      dout.close();        
   }
   
   std::vector< double > V2(g.get_n(),0);
   
   for (unsigned int i=0; i<=NUMBER; i++){
      //cout<<i<<endl;
      vce_precon(g, lseed2);
      vce_siteELoPRalg(g, lseed3, conv, true);
      Psum = g.get_Psum();
      Fsum = g.get_Fsum(); 
      if (i==0){
          cout<<"Vo "<<Psum/(double)g.get_n()<<" Fo "<<Fsum/(double)g.get_n()<<" ";
      }
      for(int j=0;j<=n;j++){//we want to average P overall the runs
         V2[j]=V2[j]+V[j];
      }
   } 
   for(int j=0;j<=g.get_n();j++){ //sum and average P over all the runs
         g.nodes[j].statP=V2[j]/(double)NUMBER; //This is the average P per site
	 if(V[j]==1 || V[j]==0) g.nodes[j].frz=true;
	 else g.nodes[j].frz=false;
	 //V2[j]=g.nodes[j].statP; //copy average P per site into V2
   }
   g.graph_calc_sum_P_F();
   Pbest=g.get_Psum()/(double)g.get_n();
   cout<<"V_1 "<<Psum/(double)g.get_n()<<" F_1 "<<Fsum/(double)g.get_n()<<" ";

   vce_siteDIG(g, lseed1, rlists); //run DIG on average P per site
   Psum = g.get_Psum();
   Fsum = g.get_Fsum();    
   cout<<"V_dig "<<Psum/(double)g.get_n()<<" F_dig "<<Fsum/(double)g.get_n()<<" ";
   
   for (unsigned int i=0; i<=NUMBER; i++){  //keep running s EloPR to see if we get better than the previous average
  //    cout<<i<<endl;
      vce_precon(g, lseed2);
      vce_siteELoPRalg(g, lseed3, conv, true);
      Psum = g.get_Psum();
      //Fsum = g.get_Fsum(); 
      Pave=Psum/(double)g.get_n();
      if (Pave < Pbest) { //if you get a state better than Pbest make it the new best state
         Pbest=Pave;
	     for(int j=0;j<=g.get_n();j++){
            V2[j]=g.nodes[j].statP;//save the best state.
         }
      }
   }
   
   for(int j=0;j<=g.get_n();j++){
      g.nodes[j].statP=V2[j];//save the best state.
   }
   g.graph_calc_sum_P_F();
   Psum = g.get_Psum();
   Fsum = g.get_Fsum(); 
   cout<<"V_f "<<Psum/(double)g.get_n()<<" F_f "<<Fsum/(double)g.get_n()<<" ";
     
   vce_siteDIG(g, frz, lseed1, true);
   Psum = g.get_Psum();
   Fsum = g.get_Fsum();      
   cout<<"V_dig "<<Psum/(double)g.get_n()<<" F_dig "<<Fsum/(double)g.get_n()<<endl;
      
   delete[] V2;
   delete[] seed;   
   return;     
}
//***********************************************************************
//void vce_Flip(graph g){ //put in graph.cpp?
 //  double dummy;
//   for(int j=0;j<g.get_n();j++){
 //     dummy=g.nodes[j].stateP;
//      g.nodes[j].stateP=1-dummy;
//   } 
//   return;	
//}
//***********************************************************************
void vce_sFlip(graph g,bool dsply, bool write, bool fulloutput, bool rlists,
       long &lseedin, bool conv, std::string &dname, char* filename){
   double Psum=0, Fsum=0, Pbest=0.66, Pave=0; 
   int l=0;
   int nE=0;
   
   std::vector< long > seed(4,0);  //seed prep
   seed[0]=lseedin;
   vcd_random_number_prep(seed);
   long lseed1 = seed[1];
   long lseed2 = seed[2];
   long lseed3 = seed[3];   
   
   g.seed=lseed1;
   g.build_graph(write);
   double c=g.get_c(); //calculates 
   if (dsply) g.displaygrph();
   
   if (fulloutput) {
      std::ofstream dout(dname.data(), ios_base::app);
      dout<<g.get_n()<<" "<<g.p<<" "<<g.c<<" "<<g_get_nE()<<" "<<lseedin;
      dout.close();        
   }
   
   vce_precon(g, lseed2);
   vce_siteELoPRalg(g, lseed3, conv, rlists);
   Psum = g.get_Psum();
   Fsum = g.get_Fsum(); 
   Pave=Psum/(double)g.get_n();
   cout<<"<V> "<<Psum/(double)g.get_n()<<" <F> "<<Fsum/(double)g.get_n()<<endl;
   do{    // Keep flipping and running EloPR until you get the best run. 2times min.
      l++;
      Pbest=Pave;
      g.flip();
      vce_siteELoPRalg(g, lseed3, conv, rlists);
      Psum = g.get_Psum();

      Pave=Psum/(double)g.get_n();
     // cout<<Pave<<" "<<Pbest<<" "<<l<<endl;
   }while(Pave<Pbest);
   
   Fsum = g.get_Fsum(); 
   cout<<"<V> "<<Psum/(double)g.get_n()<<" <F> "<<Fsum/(double)g.get_n()<<" "<<l<<endl;  
     
   vce_siteDIG(g, lseed1, true);
   Psum = g.get_Psum();
   Fsum = g.get_Fsum();     
   cout<<"<V_dig> "<<Psum/(double)g.get_n()<<" <F_dig> "<<Fsum/(double)g.get_n()<<endl;
   seed.clear();
   
   return;
}
//***********************************************************************
void vce_sDb_flip(graph g,bool dsply, bool write, bool fulloutput, bool rlists,
       long &lseedin, bool conv, std::string &dname, char* filename){
   unsigned NUMBER=100;
   double Psum=0, Fsum=0, Pbest=100, Pave2=0, Pave=0; 
   int nE=0;
   
   std::vector< long > seed(4,0);  //seed prep
   seed[0]=lseedin;
   vcd_random_number_prep(seed);
   long lseed1 = seed[1];
   long lseed2 = seed[2];
   long lseed3 = seed[3];   

   std::vector< double > V2(g.get_n(),0);

   g.seed=lseed1;
   g.build_graph(write);
   c=g.get_c(); //calculates 
   if (dsply) g.displaygraph();
   
   if (fulloutput) {
      std::ofstream dout(dname.data(), ios_base::app);
      dout<<n<<" "<<p<<" "<<c<<" "<<nE<<" "<<lseedin;
      dout.close();        
   }
   //setting Pbest from an average over a set of runs
   for (unsigned int i=0; i<=NUMBER; i++){//Run ELoPR a NUMBER of times
      //cout<<i<<endl;
      vce_precon(g, lseed2);
      vce_siteELoPRalg(g, lseed3, conv, rlists);
      for(int j=0;j<=n;j++){
         V2[j]=V2[j]+g.nodes[j].stateP;
      }
   } 
   for(int j=0;j<g.get_n();j++){ //set P to the average of the state P
         g.nodes[j].stateP=V2[j]/(double)NUMBER;
	 if(g.nodes[j].stateP==1 || g.nodes[j].stateP==0) g.nodes[j].frz=true;
	 else g.nodes[j].frz=false;
	// V2[j]=V[j];
   }
   Psum = g.get_Psum();
   Fsum = g.get_Fsum();   
   Pbest=Psum/(double)g.get_n();
   cout<<"<V> "<<Pbest<<" <F> "<<Fsum/(double)g.get_n()<<endl;

   vce_siteDIG(g, lseed1, true);
   Psum = g.get_Psum();
   Fsum = g.get_Fsum();     
   cout<<"<V_dig> "<<Psum/(double)g.get_n()<<" <F_dig> "<<Fsum/(double)g.get_n()<<endl;
   int l=0;
   
   for(int j=0;j<=g.get_n();j++){ //save Pbest state
      V2[j]=g.nodes[j].stateP;
   }
   
   do{
      l++;
      g.flip();
      vce_bELoPRalg(g, lseed3, conv, rlists);
      Psum = g.get_Psum();
      Fsum = g.get_Fsum();  
      Pave=Psum/(double)g.get_n();
     // cout<<Pave<<" "<<Pbest<<" "<<l<<endl;
   }while(Pave>Pbest && l<NUMBER); //we should have an out incase we never get better than Pbest
   if (Pave>Pbest){ //if we didn't find better than Pbest restore Pbest
	   
	  for(int j=0;j<=g.get_n();j++){ //save Pbest state
        g.nodes[j].stateP=V2[j];
      }   
   }
   else Pbest=Pave;
   //repeat
   for (unsigned int i=0; i<=NUMBER; i++){
  //    cout<<i<<endl;
      vce_precon(g, lseed2);
      vce_siteELoPRalg(g, lseed3, conv, rlists);
      Psum = g.get_Psum();
      Fsum = g.get_Fsum();  
      Pave=Psum/(double)g.get_n();
      if (Pave < Pbest) {//save the best state
         Pave2=Pbest;
         Pbest=Pave;
         for(int j=0;j<=g.get_n();j++){ //save Pbest state
            V2[j]=g.nodes[j].stateP;
          }
      }
   }
   for(int j=0;j<=g.get_n();j++){ //restore the best state
      g.nodes[j].statP=V2[j];//save the best state.
   }
   g.graph_calc_sum_P_F();
   Psum = g.get_Psum();
   Fsum = g.get_Fsum(); 
   cout<<"<V> "<<Psum/(double)g.get_n()<<" <F> "<<Fsum/(double)g.get_n()<<endl;
     
   vce_siteDIG(g, frz, lseed1, true);
   Psum = g.get_Psum();
   Fsum = g.get_Fsum();    
   cout<<"<V_dig> "<<Psum/(double)g.get_n()<<" <F_dig> "<<Fsum/(double)g.get_n()<<endl;
      
   delete[] V2;
   delete[] seed;   
   return;     
}



/***********************************************************************/
/** Find a greedy cover                                                **/
/***********************************************************************/

/************* crprcfns_ers_lf() ************************/
/** chooses node with max connectivity                 **/
/** written by c.fay                                   **/
/**                                                    **/
/** Returns:  nothing                                  **/
/********************************************************/
int cvrfns_choose_maxc(graph g, std::vector<int> nlist, int &maxc){
   int n1, n1s, l1=nlist.size();
   bool found=false;
   if (maxc==0){
	 n1=-1; 
     return n1;
   }
   do{   
       n1=nodes[l1];
       n1s=g.nodes[n1].edges.size();
       l1--;
   }while(n1s<maxc && l1>=0)||(!found))
   if (n1s >= maxc) found=true;
   if (!found && maxc >0) {
      maxc--;
      n1 = cvrfns_choose_maxc(g, nlist, maxc);
      }
   if (
   //maxc=g.nodes[n1].edges.size(); I don't think I need this line
  // cout<<maxc<<" "<<endl;
   return (n1);
}
/************* cvrfns_greedy_maxc()**********************/
/** Greedy cover by choosing node of maximum c         **/
/** written by c.fay                                   **/
/**                                                    **/
/** Returns:  size of cover (int)                      **/
/********************************************************/
int cvrfns_greedy_maxc(graph g, nlist){ //put in graph.cpp i think
   int l=0,n1,n2,n2size,n3,n3size,lsize,n4;
   int t, n1start, n1stop, n3start, n3stop;
   int nCtest=0, nL,E=nE, 
   int maxc=g.get_n()-1; //max possible n
   
  // for (int i=0; i<g.get_n(); i++){
  //   nlist[i]=i;
  //   }

   E=g.calc_nE(); /*number of edges and max c*/
  // cout<<maxc<<endl;
   
   while(E>0 || nlist.size()>0){ // while there are Edges 
      n1 = cvrfns_choose_maxc2(g, nlist, maxc);   //select a node = n1 with max connectivity
    //  cout<<n1<<" "<<maxc<<" "<<E<<endl;
      g.nodes[n1].stateP=1; //cover node 1
      g.nodes[n1].frz=true; //set as frz 
      for (int j=0; j < g.nodes[n1].edges.size(); j++){//for all edges conneded to n1
	     n2=g.nodes[n1].edges[j];
         g.eraseEdge(n1, n2);  //Erase all edges of n1
         E-- 
	  }  
      nCtest++;
  }
   int E2=g.calc_nE(); 
   if (E2>0) nCtest=nCtest+cvrfns_greedy_maxc(g, nlist);

   return (nCtest);
}
/************* cvrfns_greedy() **************************/
/** Leaf Removal Algorithm                             **/
/** written by c.fay                                   **/
/**                                                    **/
/** Returns:  nothing                                  **/
/********************************************************/
int cvrfns_greedy(graph g, std::vector<int> nlist){
   int l=0,n1,n2,n2size,n3,n3size,lsize,n4;
   int t, n1start, n1stop, n3start, n3stop;
   int nCtest=0, nL, maxc=0;
   int *nodes=new int[SIZE+1];
   int E=g.get_nE(); /*number of edges*/
   
   l=g.get_n();
   while(E>0 || nlist.size()>0){ /* while there are Edges */
      n1=nlist.end();
	  nlist.pop_back();  //select a node from nlist
 //     n1 = cvrfns_fnd_node(l, size, nodes);   /* select a node = n1*/
      g.nodes[n1].stateP=1; //cover node 1
      g.nodes[n1].frz=true; //set as frz    
      for (int j=0; j < g.nodes[n1].edges.size(); j++){ /* look at neighbors of n1 */
      	 n3=g.nodes[n1].edges[j];
      	// g.nodes[n3].stateP=0; //uncover n3
       //  g.nodes[n3].frz=true; //set as frz 
      	 g.eraseEdge(n1, n3);//erase edge n1,n3
	     E--;
	    // if (size[n3-1]==0) cvrfns_ers_node(n3, l, nodes);/*remove n3 from list */
         }	    
      }
      nCtest++;  
   }
   int E2=g.calc_nE(); 
   if (E2>0) nCtest=nCtest+cvrfns_greedy(g, nlist);

//   delete[] nodes;
   return (nCtest);
}
/************************************************************************************/
void cvf_Greedy(graph g, bool write, bool dsply, bool fulloutput,
      std::string &dname, std::vector<double> &datum, long &lseedin, bool histo, 
      bool conv, bool rlists, char* filename){  //rlists will be the trigger between Greedy and Greedy_maxc

   int nE=0;
   double Psum=0; 
   long lseed=lseedin; 
   
   std::vector< long > seed(4,0);  //seed prep
   seed[0]=lseedin;
   vcd_random_number_prep(seed);
   long lseed1 = seed[1];
   long lseed2 = seed[2];
   long lseed3 = seed[3];  

   g.seed=lseed1;
   g.build_graph(write);
   c=g.get_c(); //calculates 
   if (dsply) g.displaygraph();
   
   if (fulloutput) {
      std::ofstream dout(dname.data(), ios_base::app);
      dout<<g.get_n()<<" "<<g.p<<" "<<g.c<<" "<<g.get_nE()<<" "<<lseedin; 
     dout.close();       
   }
   datum[0]=(double)g.get_n();
   datum[1]=g.p;
   datum[2]=g.c;
   datum[3]=(double)g.get_nE();
   
   std::vector< int >  nlist(g.get_n(),0);
   vce_randomize_node_list(n, nlist, idum); //make a list of nodes and randomize

   if (rlists) datum[4]=(double) cvrfns_greedy_maxc(g, nlist)/(double)n;
   else  datum[4]=(double) cvrfns_greedy(g, nlist)/(double) n;
   
   if (dsply){
      std::ofstream output1("dsplygreedy.csv");
      grphfns_dsplygrph_wstate(n, gv, g2, V, output1, "G");
      output1.close();
   }
   if (fulloutput) {
      std::ofstream dout(dname.data(), ios_base::app);
      dout<<" Psum_greedy "<<Psum/n;
      dout.close();   
   }  
       
   return;     
}

