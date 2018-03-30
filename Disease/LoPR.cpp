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
#include <ostream>
#include <iomanip>

#include "tp.hpp"
#include "grphfns_a.hpp"
#include "LoPR_alg_a.hpp"
#include "LoPR_basic_utils.hpp"

using namespace std;

int main(int argc, char *argv[]){
   int argz = 1, n, nE, L;
   //int do_DIG=0, do_DIG_best=0, do_bE=0, do_sb=0, do_sbD=0;
   //int do_bDIG=0, do_HELP=0, do_TRI=0, do_sFbF=0;
   //int numfrz=0, numfrz2=0, numfrz3=0, numfrz4=0;
   //long idum, lseed, lseed2;
   int *gsize, *gv, *g2;
   gsize=new int[SIZE+1];
   g2=new int[SIZE+1];
   
   double Psum=0;
   double p=3.0;
   bool write=false, dsply=false, histo=false;
   char *filename;
   char LoPR;
   int grphtype=1;  /* 0 for rndm, 1 for pt */
   long lseed[3]={0},lseedstart=-1;
   
   while((argv[argz][0] == '-')&&(!isdigit(argv[argz][1]))&&(argc!=0))
   {                                              
      if(strcmp(argv[argz], "-rnd") == 0)
        grphtype=0;
      else if(strcmp(argv[argz], "-pt") == 0)
        grphtype=1;      
      else if(strcmp(argv[argz], "-ptr3") == 0)
      	grphtype=2;
      else if(strcmp(argv[argz], "-fcc") == 0)
        grphtype=3;
      else if(strcmp(argv[argz], "-cubic") == 0)
      	grphtype=4;
      else if(strcmp(argv[argz], "-sq") == 0)
      	grphtype=5;
      else if(strcmp(argv[argz], "-sclf") == 0)
      	grphtype=6;
      else if(strcmp(argv[argz], "-rndr3") == 0)
      	grphtype=7;
      else if(strcmp(argv[argz], "-rndg3") == 0)
      	grphtype=8;      
      else if(strcmp(argv[argz], "-sb") == 0)
     	LoPR='1';
      else if(strcmp(argv[argz], "-sDig") == 0)
      	LoPR='D';
      else if(strcmp(argv[argz], "-sDbest") == 0)
      	LoPR='T';
      else if(strcmp(argv[argz], "-sFlip") == 0)
      	LoPR='F';
      else if(strcmp(argv[argz], "-w") == 0)
      	write=true;
      else if(strcmp(argv[argz], "-histo") == 0){
      	if (!histo) { 
	   histo=true;
	   cout<<"<P> histogram on ";
	}
	else {
	   histo=false;
	   cout<<"<P> histogram off ";
	}
	}
      else if(strcmp(argv[argz], "-seed") == 0){
        argz++;
        lseedstart=-1*atoi(argv[argz]);
        }
      else if(strcmp(argv[argz], "-dsply") == 0)
      	if (dsply) dsply=false;
	else dsply=true;
      else if(strcmp(argv[argz], "-rgml") ==0||strcmp(argv[argz], "-r") ==0){
	  argz++;
	  filename = argv[argz];
	  cout<<"read "<< filename<<" ";
	  grphtype=-1;
	  }      
      else if(strcmp(argv[argz], "-h") ==0)
         LoPR_usage(1);	
      else{
        exit(1);
      }
      argz++; 
   }

   n = atoi(argv[argz++]);
   sscanf(argv[argz++], "%lf", &p);
   L=n;
    
   if (grphtype==3 || grphtype==4) {
      n=L*L*L;
     if (grphtype==3) gv=new int[12*SIZE+1];
   }
   else gv=new int[6*SIZE+1];
   for (int i=0; i<=SIZE; i++){
     g2[i]=0;
     gsize[i]=0;
     gv[i]=0;
     gv[i+SIZE]=0;
     gv[i+2*SIZE]=0;
     gv[i+3*SIZE]=0;    
     gv[i+4*SIZE]=0;
     gv[i+5*SIZE]=0;
     if (grphtype==3){
        gv[i+6*SIZE]=0;
        gv[i+7*SIZE]=0;
        gv[i+8*SIZE]=0;
        gv[i+9*SIZE]=0;
        gv[i+10*SIZE]=0;
	gv[i+11*SIZE]=0;
        }
     }
      
   if (lseedstart==-1) lseedstart=time(0);
   long lseedblnk=lseedstart;
   for (unsigned int i=1;i<=1000;i++){
        ran2(&lseedblnk);
        }
   for (unsigned int i=0;i<3;i++){  
      lseed[i]=(-1*(long) (ran2(&lseedblnk)*1e8));
   }
   
   cout<<n<<" "<<p<<" "<<lseedstart<<" "; 
   long seed=lseed[0];  
   switch (grphtype){
      case -1:{
         grphfns_read_gml(filename, gsize, gv, g2, n, nE);
         break;
      }
      default:{
         grphfns_build_grph_n(n, nE, p, seed, gsize, gv, g2, grphtype, write, L);
         break;
      }
   }
   double c=2.0*(double)nE/(double)n;
   if (false) {
      cout<<n<<" "<<p<<" "<<c<<" "<<nE<<" "<<lseed[0]<<" "<<lseed[1]<<" "<<lseed[2];        
   }
   else cout<<c<<" "<<nE<<" "<<" ";
   
   switch  (LoPR){
      case 'T':{
      //cout<<"sDig_best \n";
      vce_sDig_best(gv, g2, gsize, n, write, true, p, lseed[0], lseed[1], lseed[2]);
      break;
      } 
      case 'F':{
      cout<<"sFlip \n";
      vce_sFlip(gv, g2, gsize, n, write, true, p, lseed[0], lseed[1], lseed[2]);
      break; 
      }
      default: {
      cout<<"sDbD \n";
      double *V;
      V=new double[SIZE+1];
      bool *frz;
      frz=new bool[SIZE+1];
      for (int i=0; i<=SIZE; i++){
         V[i]=0;
         frz[i]=false;
      }
      std::vector< int >  nlist(n,0);

      seed=lseed[1];
      vce_precon(gv, g2, gsize, n, V, frz, seed);
   //vce_randomize_node_list(n, nlist, lseed);
   
   //for(unsigned int i=0; i<nlist.size(); i++){
   //   cout<<nlist[i]<<endl;  
  // }
   
      double Fsum=0;
      seed=lseed[2];
      vce_siteELoPRalg(gv, g2, gsize, n, V, frz, nlist, seed);
      vce_sum_P_F(V, Psum, Fsum, n);
      if (dsply){
         std::ofstream output1("dsplyg.csv");
         grphfns_dsplygrph_wstate(n, gv, g2, V, output1, "s");
         output1.close();
      }
      if (histo) vce_histogram(gv, g2, V, 100, n, "s");
      cout<<"<V_s> "<<Psum/(double)n<<" <F_s> "<<Fsum/(double)n;
   
      vce_siteDIG(gv, g2, gsize, n, V, frz, nlist, seed);
      vce_sum_P_F(V, Psum, Fsum, n);   
      if (dsply){
         std::ofstream output1("dsplyg.csv");
         grphfns_dsplygrph_wstate(n, gv, g2, V, output1, "sd");
         output1.close();
      }
      if (histo) vce_histogram(gv, g2, V, 100, n, "s");
      cout<<" <V_dig> "<<Psum/(double)n<<" <F_dig> "<<Fsum/(double)n<<endl;
   
      Fsum=0;
      seed=lseed[1];
      vce_precon(gv, g2, gsize, n, V, frz, seed);
      seed=lseed[0];
      vce_bELoPRalg(gv, g2, gsize, n, V, frz, nlist, seed);
      vce_sum_P_F(V, Psum, Fsum, n);
      if (dsply){
         std::ofstream output1("dsplyg.csv");
         grphfns_dsplygrph_wstate(n, gv, g2, V, output1, "b");
         output1.close();
      }
      if (histo) vce_histogram(gv, g2, V, 100, n, "b"); 
      cout<<"<V_b> "<<Psum/n<<" <F_b> "<<Fsum/n;
   
      vce_siteDIG(gv, g2, gsize, n, V, frz, nlist, seed);
      vce_sum_P_F(V, Psum, Fsum, n);   
      if (dsply){
         std::ofstream output1("dsplyg.csv");
         grphfns_dsplygrph_wstate(n, gv, g2, V, output1, "bsd");
         output1.close();
      }
      if (histo) vce_histogram(gv, g2, V, 100, n, "s");
      cout<<" <V_dig> "<<Psum/n<<" <F_dig> "<<Fsum/n<<endl;
      break;
   }
   }

   return 0;
}
