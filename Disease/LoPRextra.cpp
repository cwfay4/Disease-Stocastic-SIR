//***********************************************************************
void vce_siteELoPRalg_1list(graph g, long &idum, bool conv){
    double *P_old, *dP;
    P_old=new double[SIZE+1];
    for(unsigned int i=0;i<=SIZE;i++){ //make a new vector(array) of the state probability from the previous run
       P_old[i]=0;
    }
    std::vector< int >  nlist(n,0); //list of nodes that will be randomized
    double P_prod=1, delta_ave_P=100, delta_P_ave=100, Pave=1, dPsum=0;
    double Pave_old=0;
    int n1;
    int n2, nstart, nstop, t;

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
    }
    return;
}

void vce_sb(graph g, bool write, bool dsply, bool fulloutput,
      std::string &dname, std::vector<double> &datum, long &lseedin, bool histo, 
      bool conv, char* filename){
   int nE=0;
   double Psum=0, Fsum=0; 
   
   long lseed=lseedin; 
   //********************************** Setting up random numbers ********
   for (unsigned int z=1;z<=1000;z++){//initialize random number generator
       float fseed=ran2(&lseed);
       }
   std::vector< long > seed(3,0); 
   for (unsigned int z=1;z<=3;z++){
       seed[z-1]=(-1*(int) (ran2(&lseed)*1e8));
       }
   long lseed1 = seed[0];
   long lseed2 = seed[1];
   long lseed3 = seed[2]; 
   //*********************************************************************
   g.build_graph(write);
   if (dsply){
      std::ofstream output1("dsplyg.csv");\
      grphfns_dsplygrph(g, output1);
      output1.close();
   }
   c=g.get_c();
   if (fulloutput) {
      std::ofstream dout(dname.data(), ios_base::app);
      dout<<n<<" "<<p<<" "<<c<<" "<<nE<<" "<<lseedin;
      dout.close();        
   }
   
   //saves data
   datum[0]=(double)n;
   datum[1]=p;
   datum[2]=c;
   datum[3]=(double)nE;
         
   //double *V; //V is the node states not needed.
   //V=new double[SIZE+1];
   //bool *frz;
   //frz=new bool[SIZE+1];   
   //for (int i=0; i<=SIZE; i++){
   //   V[i]=0;
   //   frz[i]=false;
   //}
   
   vce_precon(g, &seed)
   //vce_precon(gv, g2, gsize, n, V, frz, lseed2);
   
   vce_siteELoPRalg(g, seed, conv); //rerun siteEloPR
   g.graph_calc_sum_P_F();
   Psum=g.get_Psum();
   Fsum=g.get_Fsum();   
   
   if (histo) vce_histogram(g, 100, "s");
   if (fulloutput) {
      std::ofstream dout(dname.data(), ios_base::app);
      dout<<" sELoPR: "<<Psum/(double) g.get_n();
      dout.close();
   }
   datum[4]=(double)Psum/(double) g.get_n();
   datum[5]=(double)Fsum/(double) g.get_n();      
   
   vce_precon(g, lseed2);
   vce_bELoPRalg(g, lseed1, conv);
   g.graph_calc_sum_P_F();
   Psum=g.get_Psum();
   Fsum=g.get_Fsum();  
   if (histo) vce_histogram(g, 100, "s");
   if (fulloutput) {
      std::ofstream dout(dname.data(), ios_base::app);
      dout<<" bELoPR: "<<Psum/n<<endl;
      dout.close();
   }
   datum[6]=(double)Psum/(double) g.get_n();
   datum[7]=(double)Fsum/(double) g.get_n(); 
   
   return;
}
}
