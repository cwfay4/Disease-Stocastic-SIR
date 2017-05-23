/***********************************************************************/
/** Vertex Cover triangle Algorithm                                  ***/
/**                                                                  ***/
/** C. Fay January 2016                                               ***/
/** Version 3.0.0.0   22.08.2003                                     ***/
/** Run: vertcvr gml_filname seed                                    ***/
/** Output:                                                          ***/
/***********************************************************************/
#include <cstdio>                                     
#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>

#include "random.hpp"
#include "graph.hpp"

#include "defaults.hpp"
//#include "grphfns_a.hpp"
//#include "cvrfns_a.hpp"
//#include "Table2D.hpp"
#include "LoPR_alg_o.hpp"
#include "LoPR_basic_utils.hpp"

using namespace std;

void add_data(std::vector<double>& data, std::vector<double>& datasq, int j, double datum){
   //cout<<j<<endl;
   data[j]=data[j]+datum; //sum of datum
   datasq[j]=data[j]+datum*datum; //sum of square of data
return;
}
double standard_error(double data, double datasq, double dsamp){
	double var=sqrt(datasq-(data*data/dsamp))/sqrt(dsamp-1);//calculate std error;
	return var;
}
int main(int argc, char *argv[]){
   int argz = 1, nsamp; 
   //int do_DIG=0, do_DIG_best=0, do_bE=0, do_sb=0, do_sbD=0;
   //int do_bDIG=0, do_HELP=0, do_TRI=0, do_sFbF=0;
   //int numfrz=0, numfrz2=0, numfrz3=0, numfrz4=0;
   //long idum, lseed, lseed2;
   double Psum, Fsum;
   double p=3.0;
   bool write=writedflt, dsply=dsplydflt, fulloutput=fulloutdflt, histo=clstrdflt;
   bool conv=convdflt;
   bool debug_test=false;
  // int grphtype=grphtypedflt;  /* 0 for rndm, 1 for pt */
   int bc=bcdflt;
   //Table2D data(20,2);            
   char LoPR=LoPRdflt;
   char *filename;
   long lseedstart=-1;
   
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
      else if(inputstr.compare("-bb") == 0) {
      	g.grphtype=8;
	    argz++;
	    g.set_L(atoi(argv[argz]));
	    //cout<<"nkernal: "<<L<<" ";
	  }
      else if(inputstr.compare("-rndfixedc") == 0){
      	g.grphtype=9;
	    argz++;
	    g.set_L(atoi(argv[argz]));
	    //cout<<"fixed maxc: "<<L<<" ";
	  }
      else if(inputstr.compare("-rnd") == 0) g.grphtype=0;   //random graph      
     // else if(inputstr.compare("-rp") == 0) prc = 'r';     //bond percolation
     // else if(inputstr.compare("-c") == 0)  prc = 'c';     //core percolation by leaf removal
     // else if(inputstr.compare("-t") ==0)    prc = 't';    //core percolation by triangle removal
     // else if(inputstr.compare("-ct") ==0)   prc = 'b';    //core percolation by leaf and triangle removal
      else if(inputstr.compare("-bc0") == 0)       	bc=0;
      else if(inputstr.compare("-bc1") == 0)       	bc=1;	
      else if(inputstr.compare("-bc2") == 0)       	bc=2;
      else if(inputstr.compare("-sb") == 0)      	LoPR='1';
      else if(inputstr.compare("-sb1") == 0)      	LoPR='L';
      else if(inputstr.compare("-sDig") == 0)       LoPR='D';
      else if(inputstr.compare("-sDbD") == 0)      	LoPR='2';
      else if(inputstr.compare("-sDbD1") == 0)      LoPR='3';
      else if(inputstr.compare("-sDbD1_m") == 0)   	LoPR='Z';
      else if(inputstr.compare("-greedy") == 0)     LoPR='G';
      else if(inputstr.compare("-Gmaxc") == 0)      LoPR='M';
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
     else if(j==0){
		  //ss>>n;
		  g.set_n(atoi(argv[i]));
		  j++;
		  cout<<g.get_n()<<endl;
	  }
	  else if(j==1){
		  g.p=atof(argv[i]); //works
		 // cout<<"p="<<ss.str()<<" "<<g.p<<std::endl;
		  j++;
	  }   
	  if (g.grphtype==-1) nsamp=1;
      else nsamp = atoi(argv[argz++]);
   
      if (g.grphtype==8) cout<<"nkernal: "<<g.get_L()<<" ";
      if (g.grphtype==3 || g.grphtype==4) {
         g.set_L(g.get_n());
         g.set_n(std::pow(g.get_L(),3));
      }   
 //     if (gn > SIZE) { depreciated
 //        cout<<"error n to large for compiled code.  n:"<<n<<" max_n:"<<SIZE<<endl;
 //        exit(1);   
 //     }    
      if (lseedstart==-1) lseedstart=time(0);
	  else {
		  ss>>lseedstart;
		  lseedstart=atoi(argv[i]);
		  cout<<"seed="<<ss.str()<<" "<<lseedstart<<std::endl;
		  j++;
	  }	  
//	  cout<<argv[i]<<" "<<ss<<" "<<inputstring<<endl;

   }
   
   long lseedblnk=lseedstart;
   //std::vector< long > seed(3*nsamp,0);
   std::vector< long > lseed(nsamp,0);
   if (nsamp!=1){
   for (unsigned int z=1;z<=1000;z++){
      float fseed=ran2(&lseedblnk);
      }
  // for (int i=1; i<=3*nsamp;i++){ 
  //    seed[i-1]=(-1*(int) (ran2(&lseedblnk)*1e8));
   //   }
   for (int i=1; i<=nsamp;i++){
      lseed[i-1]=(-1*(int) (ran2(&lseedblnk)*1e8));
      }
   }
   else lseed[0]=lseedblnk;

   std::string gt="p", rnt="-r2", sbc="f";
   std::ostringstream strs;
   strs << g.get_n();
   std::string nstr = strs.str();
   strs << g.p;
   std::string pstr = strs.str();  
   if (bc==1) sbc="p";
   else if (bc==2) sbc="pp";
   if (g.grphtype == 0 || g.grphtype==7) gt="r";
   else if (g.grphtype == 2) rnt="-r3";
   else if (g.grphtype == 3) gt="f";
   else if (g.grphtype == 4) gt="c";
   else if (g.grphtype == 5) gt="sq";
   else if (g.grphtype == 6) gt="sclf";
   else if (g.grphtype == 8) gt="bb";
   
  // std::ofstream dout("dout.csv");
   std::string dname="dLo-"+gt+nstr+"-"+pstr+"-"+sbc+rnt+".csv";
   //**********************************************LoPR run****************
   std::vector<double> data(5,0);
   std::vector<double> datasq(5,0);
   
   for (int i=0; i<nsamp; i++){
      long nseed=lseed[i];

      switch  (LoPR){
         case 'D':{ //s Dig
	     std::vector<double> datum(8,0);
	     data.resize(8,0);
	     datasq.resize(8,0);
            vce_sDig(g, write, dsply, fulloutput, true, dname, datum, nseed, histo, conv, filename);
            for (int j=0; j<8;j++){
               add_data(data,datasq,j,datum[j]);
            }
	    break;
	 }
	 case '2':{ //s Dig b Dig
	    std::vector<double> datum(12,0);
	     data.resize(12,0);
	     datasq.resize(12,0);
	        bool rlist = false;
            vce_sDbD(g, write, dsply, fulloutput, rlist, dname, datum, nseed, histo, conv, filename);
            for (int j=0; j<12;j++){
               add_data(data,datasq,j,datum[j]);
            }
	    break;
	 }
	 case '3':{ //s Dig b Dig 1list
	    std::vector<double> datum(12,0);
	     data.resize(12,0);
	     datasq.resize(12,0);
            vce_sDbD(g, write, dsply, fulloutput, false, dname, datum, nseed, histo, conv, filename);
            for (int j=0; j<12;j++){
               add_data(data,datasq,j,datum[j]);
            }
	    break;
	 }
	 case 'L':{ //sb 1list  outputs to default
            std::vector<double> datum(8,0);
	     data.resize(8,0);
	     datasq.resize(8,0);
            vce_sb(g, write, dsply, fulloutput, dname, datum, nseed, histo, conv, false, filename);
            for (int j=0; j<8;j++){
               add_data(data,datasq,j,datum[j]);
            }
	    break;
	 } 
	 case 'Z':{ //sb 1list  outputs to default
            std::vector<double> datum(21,0);
	     data.resize(21,0);
	     datasq.resize(21,0);
            vce_sDbD_1list_many(g, write, dsply, fulloutput, dname, datum, nseed, histo, conv, filename);
            for (int j=0; j<21;j++){
               add_data(data,datasq,j,datum[j]);
            }
	    break;
	 } 
	 case 'G':{ //s Dig
	    std::vector<double> datum(5,0);
            cvf_Greedy(g, write, dsply, fulloutput, dname, datum, nseed, histo, conv, false, filename);
            for (int j=0; j<5;j++){
               add_data(data,datasq,j,datum[j]);
            }
	    break;
	 }
	 case 'M':{ //s Dig
	    std::vector<double> datum(5,0);;
            cvf_Greedy(g, write, dsply, fulloutput, dname, datum, nseed, histo, conv, true, filename);
            for (int j=0; j<5;j++){
               add_data(data,datasq,j,datum[j]);
            }
	    break;
	 }	 
	 default:{ //sb
        std::vector<double> datum(8,0);
	     data.resize(8,0);
	     datasq.resize(8,0);
            vce_sb(g, write, dsply, fulloutput, dname, datum, nseed, histo, conv, true, filename);
            for (int j=0; j<8;j++){
               add_data(data,datasq,j,datum[j]);
            }
	    break;
	 }
     }
     //need to clear g here.
   }
   //*******************************************LoPR output********* 
   double dum, dumsq, dsamp = (double) nsamp;
   switch  (LoPR){
      case 'D':{
		 std::string label[8][10]={{"N"},{"p"},{"p_a"},{"NE"},
                 {"sLoPR_Psum"},{"Fsum"},{"sDig_Psum"},{"Fsum"}}; 
         for(int k=0; k<8; k++){
            double var=standard_error(data[k],datasq[k],dsamp);//calculate std error;
            if (k==0 || k==1) cout<<label[k][0]<<" "<<data[k]/dsamp<<" "; 
            else cout<<label[k][0]<<" "<<data[k]/dsamp<<" "<<var<<" ";
         }
	  break;
      }
      case '2':{
		 std::string label[12][10]={{"N"},{"p"},{"p_a"},{"NE"},
                 {"sLoPR_Psum"},{"Fsum"},{"sDig_Psum"},{"Fsum"},
		 {"bLoPR_Psum"},{"Fsum"},{"b_sDig_Psum"},{"Fsum"}}; 
         for(int k=0; k<12; k++){
            double var=standard_error(data[k],datasq[k],dsamp);//calculate std error;
            if (k==0 || k==1) cout<<label[k][0]<<" "<<data[k]/dsamp<<" "; 
            else cout<<label[k][0]<<" "<<data[k]/dsamp<<" "<<var<<" ";
         }
	 break;
      }
      case '3':{
		 std::string label[12][10]={{"N"},{"p"},{"p_a"},{"NE"},
                 {"sLoPR1_Psum"},{"Fsum"},{"sDig1_Psum"},{"Fsum"},
		 {"bLoPR1_Psum"},{"Fsum"},{"b_sDig1_Psum"},{"Fsum"}};
         for(int k=0; k<12; k++){
            double var=standard_error(data[k],datasq[k],dsamp);//calculate std error;
            if (k==0 || k==1) cout<<label[k][0]<<" "<<data[k]/dsamp<<" "; 
            else cout<<label[k][0]<<" "<<data[k]/dsamp<<" "<<var<<" ";
         }
	 break;
      } 
      case 'Z':{
		 std::string label[20][10]={{"N"},{"p"},{"p_a"},{"NE"},
                 {"sLoPR1_Psum"},{"Fsum"},{"sDig1_Psum"},{"Fsum"},
		 {"bLoPR1_Psum"},{"Fsum"},{"b_sDig1_Psum"},{"Fsum"},
		 {"sLoPR_Psum"},{"Fsum"},{"sDig_Psum"},{"Fsum"},
		 {"bLoPR_Psum"},{"Fsum"},{"b_sDig_Psum"},{"Fsum"}};
         for(int k=0; k<20; k++){
            double var=standard_error(data[k],datasq[k],dsamp);//calculate std error;
            if (k==0 || k==1) cout<<label[k][0]<<" "<<data[k]/dsamp<<" "; 
            else cout<<label[k][0]<<" "<<data[k]/dsamp<<" "<<var<<" ";
         }
	 break;
      }
      case 'G':{
		 std::string label[5][10]={{"N"},{"p"},{"p_a"},{"NE"},
                 {"Psum_greedy"}};
         for(int k=0; k<5; k++){
            double var=standard_error(data[k],datasq[k],dsamp);//calculate std error;
            if (k==0 || k==1) cout<<label[k][0]<<" "<<data[k]/dsamp<<" "; 
            else cout<<label[k][0]<<" "<<data[k]/dsamp<<" "<<var<<" ";
         }
	 break;
      }
     case 'M':{
		 std::string label[5][10]={{"N"},{"p"},{"p_a"},{"NE"},
                 {"Psum_Gmaxc"}};
         for(int k=0; k<5; k++){
            double var=standard_error(data[k],datasq[k],dsamp);//calculate std error;
            if (k==0 || k==1) cout<<label[k][0]<<" "<<data[k]/dsamp<<" "; 
            else cout<<label[k][0]<<" "<<data[k]/dsamp<<" "<<var<<" ";
         }
	 break;
      }
      default:{
		 std::string label[8][10]={{"N"},{"p"},{"p_a"},{"NE"},
                 {"sLoPR_Psum"},{"Fsum"},{"bLoPR_Psum"},{"Fsum"}};
         for(int k=0; k<8; k++){
            double var=standard_error(data[k],datasq[k],dsamp);//calculate std error;
            if (k==0 || k==1) cout<<label[k][0]<<" "<<data[k]/dsamp<<" "; 
            else cout<<label[k][0]<<" "<<data[k]/dsamp<<" "<<var<<" ";
         }
	 break;
      }
   }
   cout<<"nsamples: "<<nsamp<<endl;
   
   return 0;
}

