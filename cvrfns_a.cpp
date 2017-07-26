/***********************************************************************/
/** cvrfns  _array                                                   ***/
/**                                                                  ***/
/** C. Fay Jan 2003                                                  ***/
/** Version 2.1.2.0   04.27.2007                                     ***/
/***********************************************************************/

#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <string>
#include "tp.hpp"
#include "grphfns_a.hpp"
#include "Table2D.hpp"

using namespace std;
/***********************************************************************/
/** Find a greed cover                                                **/
/***********************************************************************/
/************* crprcfns_cnt_lv() ************************/
/** counts number of leaves                            **/
/** written by c.fay                                   **/
/**                                                    **/
/** Returns:  nothing                                  **/
/********************************************************/
int cvrfns_cnt_edges(int n, int size[], int nodes[], int &maxc){
   int n1s=0, i=0, l=0;
   maxc=0;
   double ns=0;
   for (i=0; i<n; i++){ /*count # of egdes and make a list of nodes*/
      n1s=n1s+size[i];
      nodes[l]=i+1;
      l++;
      if (size[i]>maxc) maxc=size[i];
   }
   ns=(double) n1s/ 2.0; 
   n1s=(int) ns;
   return (n1s);
}
/************* crprcfns_cnt_lv() ************************/
/** counts number of leaves                            **/
/** written by c.fay                                   **/
/**                                                    **/
/** Returns:  nothing                                  **/
/********************************************************/
int cvrfns_find_maxc(int n, int size[]){
   int maxc=0;
   for (int i=0; i<n; i++){ /*count # of egdes and make a list of nodes*/
      if (size[i]>maxc) maxc=size[i];
   }
   return (maxc);
}
/************* crprcfns_fnd_lv() ************************/
/** finds a leaf                                       **/
/** written by c.fay                                   **/
/**                                                    **/
/** Returns:  nothing                                  **/
/********************************************************/
int cvrfns_fnd_node(int &l, int size[], int nodes[]){
   int n1, n1s;
   n1=nodes[l-1];
   n1s=size[n1-1];
   nodes[l-1]=0;
   l--;
   if (n1s==0) n1 = cvrfns_fnd_node(l, size, nodes); 
   return (n1);
}
/************* crprcfns_ers_lf() ************************/
/** erases a leaf                                      **/
/** written by c.fay                                   **/
/**                                                    **/
/** Returns:  nothing                                  **/
/********************************************************/
void cvrfns_ers_node(int n2, int &l, int nodes[]){
   int i, lsize,test;
   lsize=l;	       
   for (i=0; i<=l-1; i++){
      test=nodes[i];
      if (test==n2){
         nodes[i]=nodes[l-1];
         nodes[l-1]=0;
         lsize--;
       }
   }
   /*cout<<"erase: "<<n2<<" ";*/
   l=lsize;
   return;  
}
/************* crprcfns_ers_lf() ************************/
/** chooses node with max connectivity                 **/
/** written by c.fay                                   **/
/**                                                    **/
/** Returns:  nothing                                  **/
/********************************************************/
int cvrfns_choose_maxc(int &l, int size[], int nodes[], int &maxc){
   int n1, n1s, l1=l;
   //int maxc0=maxc
   n1=nodes[l-1];
   n1s=size[n1-1];
   bool found=false;
   while (n1s<maxc && l1>0){
       l1--;
       n1=nodes[l1-1];
       n1s=size[n1-1];
   }
   if (n1s >= maxc) found=true;
   if (!found && maxc >0) {
      maxc--;
      n1 = cvrfns_choose_maxc(l, size, nodes, maxc);
      }
   cvrfns_ers_node(n1, l, nodes);
   maxc=size[n1-1];
  // cout<<maxc<<" "<<endl;
   return (n1);
}
int cvrfns_choose_maxc2(int n, int size[], int &maxc){
   int n1, n1s=0, i=0;
  // cout<<maxc<<" ";
   while (n1s<maxc && i<n){
      i++;
      n1s=size[i-1];
    //  cout<<i<<" "<<n1s<<endl;
   }
   n1=i;
  // cout<<n1<<" "<<n1s<<" ";
   if (n1s < maxc) {
      if (maxc>0) maxc--;
      else maxc=cvrfns_find_maxc(n, size);
      n1 = cvrfns_choose_maxc2(n, size, maxc);
      }
   maxc=size[n1-1]; 
  // cout<<n1<<" "<<n1s<<" "<<maxc<<endl;
   return (n1);
}   
/************* crprcfns_ers_lf() ************************/
/** chooses node with max connectivity                 **/
/** written by c.fay                                   **/
/**                                                    **/
/** Returns:  nothing                                  **/
/********************************************************/
void vce_randomize_node_list(int n, int nodes[], long &idum){
   std::vector<int> nlist2(n,0);
   int j;
   double dj;
   
   for(unsigned int i=0; i<n; i++){
      nlist2[i]=i+1;
   }
   
   for(unsigned int i=0; i<n; i++){
      dj=fmod(ran2(&idum)*nlist2.size(),nlist2.size());
      j=(int) dj;
      nodes[i]=nlist2[j];
      //std::vector<int>::iterator 
      //     p1=find(nlist.begin(), nlist.end(), nlist2[i]);      
      nlist2.erase(find(nlist2.begin(), nlist2.end(), nodes[i])); //not right synatatically
   }
   return;
}
/************* cvrfns_greedy() **************************/
/** Leaf Removal Algorithm                             **/
/** written by c.fay                                   **/
/**                                                    **/
/** Returns:  nothing                                  **/
/********************************************************/
int cvrfns_greedy(int n, int nE, int size[], int gv[], int g2[], double V[]){
   int l=0,n1,n2,n2size,n3,n3size,lsize,n4;
   int t, n1start, n1stop, n3start, n3stop;
   int nCtest=0, nL, maxc=0;
   int *nodes=new int[SIZE+1];
   int E=nE;
   
   for (int i=0; i<=SIZE; i++){
     nodes[i]=0;
     }

   E=cvrfns_cnt_edges(n, size, nodes, maxc); /*number of initial leaves*/
   
   l=n;
   while(E>0){ /* while there are Edges */
      n1 = cvrfns_fnd_node(l, size, nodes);   /* select a node = n1*/
      grphfns_find_start_stop(n1-1,g2,n1start,n1stop,t);
      for (int k1=n1start; k1<n1stop; k1++){ /* look at neighbors of n1 */
      	 n3=gv[k1];
         if (n3 != 0){
	    gv[k1]=0; /*remove out edge n1-n3*/
	    size[n1-1]--;
	    E--;
            grphfns_find_start_stop(n3-1,g2,n3start,n3stop,t);
	    for(int k3=n3start; k3< n3stop; k3++){  /*look at neighbors of n3 */
               n2=gv[k3];
	         if (n2==n1){       /*remove in edge n3-n1(n2)*/
                   gv[k3]=0;
                   size[n3-1]--;
                 }
	    }
	    if (size[n3-1]==0) cvrfns_ers_node(n3, l, nodes);/*remove n3 from list */
         }	    
      }
      nCtest++;
      size[n1-1]=0;
      V[n1-1]=1.0;     
      /* if (size[n1-1]==0)crprcfns_ers_node(n1, l, node);*/
   }

   delete[] nodes;
   return (nCtest);
}
/************* cvrfns_greedy_maxc()**********************/
/** Greedy cover by choosing node of maximum c         **/
/** written by c.fay                                   **/
/**                                                    **/
/** Returns:  size of cover (int)                      **/
/********************************************************/
int cvrfns_greedy_maxc(int n, int nE, int size[], int gv[], int g2[], double V[]){
   int l=0,n1,n2,n2size,n3,n3size,lsize,n4;
   int t, n1start, n1stop, n3start, n3stop;
   int nCtest=0, nL,E=nE, maxc=0;
   
   int *nodes=new int[SIZE+1];
   for (int i=0; i<=SIZE; i++){
     nodes[i]=0;
     }

   E=cvrfns_cnt_edges(n, size, nodes, maxc); /*number of initial leaves*/
  // cout<<maxc<<endl;
   
   while(E>0){ // while there are Edges 
      n1 = cvrfns_choose_maxc2(n, size, maxc);   //select a node = n1
    //  cout<<n1<<" "<<maxc<<" "<<E<<endl;
      grphfns_find_start_stop(n1-1,g2,n1start,n1stop,t);
      for (int k1=n1start; k1<n1stop; k1++){ // look at neighbors of n1 
      	 n3=gv[k1];
         if (n3 != 0){
	    gv[k1]=0; //remove out edge n1-n3
	    //size[n1-1]--;
	    //E--;
            grphfns_find_start_stop(n3-1,g2,n3start,n3stop,t);
	    for(int k3=n3start; k3< n3stop; k3++){  //*look at neighbors of n3 
               n2=gv[k3];
	         if (n2==n1){       //*remove in edge n3-n1(n2)
                   gv[k3]=0;
                   size[n3-1]--;
                 }
	    }
	    if (size[n3-1]==0) cvrfns_ers_node(n3, l, nodes);//*remove n3 from list 
         }	    
      }
      nCtest++;
      E=E-size[n1-1];
      size[n1-1]=0;
      V[n1-1]=1.0;     
   }

   delete[] nodes;
   return (nCtest);
}
/*********************************************************************************/
void cvf_Greedy(int gv[], int g2[], int gsize[], int n, int L,
      int grphtype, bool write, bool dsply, bool fulloutput,
      std::string &dname, std::vector<double> &datum, double p, 
      long &lseedin, bool histo, 
      bool conv, int bc, char* filename){

   int nE=0;
   double Psum=0; 
   long lseed=lseedin; 
   
   for (unsigned int z=1;z<=1000;z++){
         float fseed=ran2(&lseed);
         }
      std::vector< long > seed(3,0); 
      for (unsigned int z=1;z<=3;z++){
         seed[z-1]=(-1*(int) (ran2(&lseed)*1e8));
         }
   
      long lseed1 = seed[0];
      long lseed2 = seed[1];
      long lseed3 = seed[2]; 
   
   grphfns_build_grph_n(n, nE, p, lseed1, gsize, gv, g2, grphtype, write,L, bc, filename);
   if (dsply){
      std::ofstream output1("dsplyg.csv");\
      grphfns_dsplygrph(n, gv, g2, output1);
      output1.close();
   }
   
   double c=2.0*(double)nE/(double)n;
   if (fulloutput) {
      std::ofstream dout(dname.data(), ios_base::app);
      dout<<n<<" "<<p<<" "<<c<<" "<<nE<<" "<<lseedin; 
     dout.close();       
   }
   datum[0]=(double)n;
   datum[1]=p;
   datum[2]=c;
   datum[3]=(double)nE;
   
   double *V;
      V=new double[SIZE+1];
   for (int i=0; i<=SIZE; i++){
      V[i]=0;
   }
   
   datum[4]=(double) cvrfns_greedy(n, nE, gsize, gv, g2, V)/(double) n;
   
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
void cvf_Greedy_maxc(int gv[], int g2[], int gsize[], int n, int L,
      int grphtype, bool write, bool dsply, bool fulloutput,
      std::string &dname, std::vector<double> &datum, double p, 
      long &lseedin, bool histo, 
      bool conv, int bc, char* filename){

   int nE=0;
   double Psum=0; 
   long lseed=lseedin; 
   
   for (unsigned int z=1;z<=1000;z++){
         float fseed=ran2(&lseed);
         }
      std::vector< long > seed(3,0); 
      for (unsigned int z=1;z<=3;z++){
         seed[z-1]=(-1*(int) (ran2(&lseed)*1e8));
         }
   
      long lseed1 = seed[0];
      long lseed2 = seed[1];
      long lseed3 = seed[2]; 
   
   grphfns_build_grph_n(n, nE, p, lseed1, gsize, gv, g2, grphtype, write,L, bc, filename);
   if (dsply){
      std::ofstream output1("dsplyg.csv");\
      grphfns_dsplygrph(n, gv, g2, output1);
      output1.close();
   }
   
   double c=2.0*(double)nE/(double)n;
   if (fulloutput) {
      std::ofstream dout(dname.data(), ios_base::app);
      dout<<n<<" "<<p<<" "<<c<<" "<<nE<<" "<<lseedin; 
     dout.close();       
   }
   datum[0]=(double)n;
   datum[1]=p;
   datum[2]=c;
   datum[3]=(double)nE;
   
   double *V;
      V=new double[SIZE+1];
   for (int i=0; i<=SIZE; i++){
      V[i]=0;
   }
   
   datum[4]=(double) cvrfns_greedy_maxc(n, nE, gsize, gv, g2, V)/(double)n;
   
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
