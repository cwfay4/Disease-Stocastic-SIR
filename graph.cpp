#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include "random.hpp"
#include "graph.hpp"

using namespace std;


node::node() {
    size=0; //local connectivity
    id=0;
    cluster =0;
    frz=false;
    
   // edges.reserve(5);	
}
node::node(int a, int b) {  //I don't know if I am using this
    size=a;
    id=b;
    cluster =0;
    frz=false;
  
   // edges.reserve(5);	
}
graph::graph(){
	//n=0;
	num_n=0;
	nE=0;
	grphtype=0;
	boundarycond=0;
	p=0;
	c=0;
	nL=0;
	nT=0;
	//nodes=NULL;
}
graph::graph(int a){
	num_n=a;
	nE=0;
	//grphtype=0;
	//boundarycond=0;
	p=0;
	c=0;
	nL=0;
	nT=0;
	for (int i=0; i<num_n; i++){ //put num_n nodes into the graph
		addNode();
	}
}
/************************************************************************************/
int graph::get_n(){
   return num_n;	
}
int graph::get_nE(){
   return nE;	
}
int graph::calc_nE(){//set number of edges in the graph
   int sum=0;
   for(int i=0; i<num_n; i++){
       sum+=nodes[i].edges.size();   
   }
   nE=sum/2;
   return nE;
}
int graph::get_L(){
   return L;	
}
double graph::get_c(){
	c=2.0*(double)nE/(double)nodes.size();
	return c;
}
bool graph::get_debug(){
   return debug;	
}
/************************************************************************************/
void graph::set_n(int a){ //set number of nodes in the graph
   num_n=a;
   return;
}
void graph::set_nE(int a){//set number of edges in the graph
   nE=a;
   return;
}
void graph::set_L(int a){//I don't know why I wrote this
   L=a;
   return;
}
void graph::set_bc(int a){//set the boundary conditions
   bc=a;
   return;
}
void graph::set_debug(bool a){//set the boundary conditions
   debug=a;
   return;
}
/************************************************************************************/
/**** subroutines for finding the total node state Psum and number of frozen ********/
/**** sites                                                                  ********/
/************************************************************************************/
double graph::sum_del_P(double P_old[]){
   double Psum=0, Psum_old=0, sum_del_P=0;
   
   for(int i=0; i<num_n; i++){ 
      sum_del_P+=fabs(P_old[i]-nodes[i].stateP);
      Psum+=nodes[i].stateP;
      Psum_old+=P_old[i];
   }
   return sum_del_P;	
}
double graph::calc_sum_P(){
   Psum=0;
   for(int i=0; i<num_n; i++){ 
      Psum+=nodes[i].stateP;
   }
   return Psum;
}
void graph::calc_sum_P_F(){
   Psum=0;
   Fsum=0;
   for(int i=0; i<num_n; i++){ 
      Psum+=nodes[i].stateP;
      if (nodes[i].frz) Fsum+=1;
   }
   return;
}
double graph::get_Psum(){
   return Psum;	
}
double graph::get_Fsum(){
   return Fsum;	
}
void graph::flip(){ //put in graph.cpp?
   double dummy;
   for(int j=0;j<num_n;j++){
      dummy=nodes[j].stateP;
      nodes[j].stateP=1-dummy;
   } 
   return;	
}
void graph::StateHistogram(double numbox, std::string ss){
      std::vector<int> histogram((int)numbox+1,0);
      double dj, boxsize=1.0/(double) numbox;
      int j;
      std::string name, type ("P_histo.csv");
      name=ss+type;
      
      for(int i=0; i<num_n;i++){
         dj=nodes[i].stateP/boxsize;
	     j = (int) dj;
	 //cout<<P[i]<<" "<<j<<endl;
	     while(nodes[i].stateP>(double)j*boxsize) j++;
	     if (j<=(int)numbox) histogram[j]++;
	     else cout<<"error P[i] too large \n";      
      }
      std::ofstream output(name.data());
      for(int i=0; i<=(int) numbox; i++){  //in the output 1 = in MIS, 0 = in VC
         double di=(double)i*boxsize;
         output<<1-di<<" "<<histogram[i]<<endl;
	 //cout<<1-di<<" "<<histogram[i]<<endl;
      }
      output.close();
      return;	
}
/************************************************************************************/
void graph::addNode(){
    node node1;
	nodes.push_back(node1);
	
	return;	
}

void graph::populate_graph(){ //adding twice?
	//cout<<endl<<" populate "<<num_n<<" "<<nodes.size()<<endl;
	while(num_n>nodes.size()){ //add nodes until num_n = nodes.size()
		addNode();
		
	}
	//cout<<endl<<" populate "<<num_n<<" "<<nodes.size()<<endl;	
	while(num_n<nodes.size()){ //remove nodes unitl num_n = nodes.size()
		nodes.pop_back();
	}
    //cout<<endl<<" populate "<<num_n<<" "<<nodes.size()<<endl;
	return;
}
void graph::populate_graph(int nkernal){
	//cout<<endl<<" populate "<<num_n<<" "<<nodes.size()<<endl;
	while(num_n>nkernal){
		addNode();
	}
	//cout<<endl<<" populate "<<num_n<<" "<<nodes.size()<<endl;	
	while(num_n<nkernal){
		nodes.pop_back();
	}
    //cout<<endl<<" populate "<<num_n<<" "<<nodes.size()<<endl;
	return;
}

void graph::addEdge(int head, int tail){
//	if(debug) std::cout<<"edge adding:"<<head<<" "<<tail<<" "<<nE<<endl;
    nodes[head].edges.push_back(tail); //add tail to list of nodes connected to head
    nodes[tail].edges.push_back(head);
    nE++;
    if(debug) std::cout<<"edge added:"<<head<<" "<<tail<<" "<<nE<<endl;
   
   return;
}
void graph::eraseEdge(int head, int tail)
{
  int j = 0;
  for (int i = 0; i < nodes[head].edges.size(); ++i) {
    if (nodes[head].edges[i] != tail) nodes[head].edges[j++] = nodes[head].edges[i];
  }
  nodes[head].edges.resize(j);
  j = 0;
  for (int i = 0; i < nodes[tail].edges.size(); ++i) {
    if (nodes[tail].edges[i] != tail) nodes[tail].edges[j++] = nodes[tail].edges[i];
  }
  nodes[tail].edges.resize(j);
  
  return;
}

/************************************************************************************/
/**** write gml output **************************************************************/
/************************************************************************************/
void graph::write_gml(std::string outfile){
   std::ofstream output(outfile.c_str()); 

   output<<"Creator \"LEDA write_gml\" \n"<<"\n";
   output<<"graph [\n"<<"\n";
   output<<"  directed 0 \n"<<"\n";

   for (int t=0; t<nodes.size(); t++) {
      output<<"  node [ id "<<t<<" ]\n";
   }
   output<<"\n";
	
	for (unsigned int i=0; i<nodes.size(); i++){ //for every node
	    for(unsigned int j=0; j<nodes[i].edges.size(); j++){
			//cout<<"connectivity "<<nodelist[i].edge_list.size()<<" "<<i<<" "<<nodelist[i].edge_list[j]<<endl;
		   if (i<nodes[i].edges[j]){
				//cout<<i<<" "<<nodes[i].edges[j]<<endl;
				output<<"  edge [ source "<<i<<" target "<<nodes[i].edges[j]<<" ]\n";
			}	
		}	
    }
	output<<"\n]";
	output.close();
	cout<<"graph written as "<<outfile<<endl;
   return;	
}
/***********************************************************************************/
void graph::write_gml_to_screen(){
   //std::ofstream output(outfile.c_str()); 

   std::cout<<"Creator \"LEDA write_gml\" \n"<<"\n";
   std::cout<<"graph [\n"<<"\n";
   std::cout<<"  directed 0 \n"<<"\n";

   for (int t=0; t<nodes.size(); t++) {
      std::cout<<"  node [ id "<<t<<" ]\n";
   }
   std::cout<<"\n";
	
	for (unsigned int i=0; i<nodes.size(); i++){ //for every node
	    for(unsigned int j=0; j<nodes[i].edges.size(); j++){
			//cout<<"connectivity "<<nodelist[i].edge_list.size()<<" "<<i<<" "<<nodelist[i].edge_list[j]<<endl;
		   if (i<nodes[i].edges[j]){
				//cout<<i<<" "<<nodes[i].edges[j]<<endl;
				std::cout<<"  edge [ source "<<i<<"  target "<<nodes[i].edges[j]<<" ]\n";
			}	
		}	
    }
	std::cout<<"\n]";
   return;	
	
}
/***********************************************************************************/
void graph::read_gml(string filename){
   std::string temp;
   int n;
   int c, node_l, edge_l, edge_l2, l=0, max_n;
   bool node_check, edge_check;
   std::vector< int > head, tail;
  // graph g;
   
   max_n=0;
   nE=0;
   
   std::ifstream gml(filename.c_str());
   do{
      node_check=false; //are we looking at the list of nodes
      node_l=-1;
      edge_l=-1;
      edge_l2=30;
      edge_check=false; //are we looking at the list of nodes
      l=0;
      do{
      l++;
      if ((!node_check||l!=node_l)&&(!edge_check||(l!=edge_l&&l!=edge_l2))){
        gml>>temp;
      } else if (node_check && l==node_l) {
        gml >> n;
	    n++; //n counts the number of nodes
	   // cout<<"prepare to insert node"<<endl;
	   // node node1;
	   // g.nodes.push_back(node1);
	    addNode();
	   // cout<<"insert node"<<endl;
	    if (n>max_n) max_n=n;
      } else if (edge_check && l==edge_l){
        gml >> n;
	    head.push_back(n);
	    if (n>max_n) max_n=n;
      } else if (edge_check && l==edge_l2){
        gml >> n;
	    tail.push_back(n);
	    if (n>max_n) max_n=n;
      }
     
      if (temp == "node"){
        node_l=l+3;
	    node_check=true;
      } else if (temp == "edge"){
        edge_l=l+3;
	    edge_l2=l+5;
	    edge_check=true;
      }
      
      }while (temp!="]"); 
   } while((c=gml.get()) != EOF);
   
   if(max_n!=nodes.size()) cout<<"error in number of nodes \n";
   if(tail.size()!=head.size()) return; 
   nE=tail.size();
   n=max_n;
   
   //cout<<max_n;
   cout<<" graph read:: number of nodes:"<<nodes.size()<<" number of edges:"<<nE<<endl;
  
    int n1, n2, nstart, nstop, t, nE_total=nE;
    nE=0;
    for (unsigned int i=0; i<head.size(); i++){ //for all edges insert edges into nodelist.edge_list
	   n1=head[i];
	   n2=tail[i];
	   addEdge(n1, n2);
    }
    
    c=2.0*(double)nE/(double)n;
    num_n=n;
    
    if (nE != nE_total) cout<<"Error in calculating number of edges"<<nE<<" "<<nE_total<<endl;
    bool test = false; //debugging option
    if (test) {
       write_gml_to_screen();
       string outfile="temp.gml";
       write_gml(outfile);
       cout<<"read successful"<<endl;
       c=2.0*(double)nE/(double)n;
       cout<<n<<" "<<nodes.size()<<" "<<c<<" "<<nE<<" ";
       cout<<std::endl;
    }

   return;
}
/***********************************************************************************/
/******* Subroutines for Random graphs *********************************************/
/***********************************************************************************/
/***********************************************************************************/
void graph::get_ij(int n1,int root, int& i,int& j){ //get (x,y) coordinate for a node 
	double dj;
	  dj=((double)(n1+1.0)/(double)root+0.9);
      j=(int)dj;
      i=n1+1-(j-1)*root;
      if (i==0) i=1;
   return;
}
void graph::dsplygrph(){ //output a graph into two csv files so they can be plotted
   int root, n2=0, n;
   int i, j, l, m;
   double sroot, dn;
   
   n = get_n(); 
   dn = (double) n;                                     
   sroot = sqrt(dn);                  
   root = (int) (sroot+0.5);
   
   std::ofstream output1("displ1.csv");  //location of all nodes
   for (int n1=0; n1<n-1; n1++) { //for all nodes (i,j) for each node
      get_ij(n1,root,i,j);
      output1<<i<<" "<<j<<"\n";
   }   
   output1<<endl;
   output1.close(); 
   std::ofstream output2("displ2.csv");  //location of all edges
   for (int n1=0; n1<=n-1; n1++) { //for all edges
	  for (unsigned int k=0; k<nodes[i].edges.size(); k++){
          n2=nodes[n1].edges[j];  
          if(n1<n2 && n2!=0){
	           get_ij(n1,root, i,j);
	           get_ij(n2,root, l,m);
	           output2<<i<<" "<<j<<"\n";
	           output2<<l<<" "<<m<<endl<<endl;
	      }
      }
   } 
   output2.close();
   return;
}
/***********************************************************************************/
/**** Generate a head and tail set *************************************************/
void graph::GenEdge(long &idum, int &r1, int &r2){
   int l1,l2;
   double dl;
   
   dl=fmod(ran0(&idum)*num_n,num_n);
   l1=(int) dl;
   dl=fmod(ran0(&idum)*num_n,num_n);
   l2=(int) dl;
   while(l1==l2){
     dl=fmod(ran0(&idum)*num_n,num_n);
     l2=(int)dl;
   }   
   if(l1>l2){
     r2=l1;
     r1=l2;
   }
   else{
     r2=l2;
     r1=l1;
   }
  /* cout<<r1<<" "<<r2<<endl;*/
   return;
}
/***********************************************************************************/
/***********************************************************************************/
void graph::GenEdge3(long &idum, int &r1, int &r2){
   int l1,l2;
   double dl;
  
   dl=fmod(ran3(&idum)*num_n,num_n);
   l1=(int) dl;
   dl=fmod(ran3(&idum)*num_n,num_n);
   l2=(int) dl;
   while(l1==l2){
     dl=fmod(ran3(&idum)*num_n,num_n);
     l2=(int)dl;
   }   
   if(l1>l2){
     r2=l1;
     r1=l2;
   }
   else{
     r2=l2;
     r1=l1;
   }
  /* cout<<r1<<" "<<r2<<endl;*/
   return;
}
/***********************************************************************************/
/******** Generate m random edges in a graph ***************************************/
void graph::build_random(){
	long idum=seed;
	double dm=(double)num_n*p/2.0;
	double dm_max=((double)(num_n*(num_n-1)))/2.0;
    int m=int(dm+0.5), mmax=int(dm_max+0.5);
    bool no_edge=true;
    int n1, n2;
	
	populate_graph();	
	
	if (debug) cout<<" N:"<<num_n<<" no. Edges:"<<m<<" no. E max:"<<mmax<<" ";
	if (m>mmax) m=mmax;
	while (nE<m) {
		GenEdge(idum, n1, n2);
		//if (debug) cout<<idum<<" "<<n1<<" "<<n2<<" ";
		no_edge=true;
		for(int j=0; j<nodes[n1].edges.size(); j++){
			if (nodes[n1].edges[j]==n2) no_edge=false;
			//if (debug) cout<<" no edge "<<no_edge<<" "<<false<<" ";
		}
		if (no_edge) {
			addEdge(n1,n2); //adds edge incriments nE++;
		    //if (debug) cout<<"Add edge "<<n1<<" "<<n2<<" "<<nE<<" "<<m<<endl;
		}
	}
    cout<<"End of build random _ no. Edges:"<<nE<<" ";
    return;	
}
void graph::build_random(int nkernal){ //builds a random kernal
	long idum=seed;
	double dm=(double)nkernal*p/2.0;
    int m=int(dm+0.5);
    bool no_edge=true;
    int n1, n2;
	
	populate_graph(nkernal);	
	
	while (nE<m) {
		GenEdge(idum, n1, n2);
		if (debug) cout<<idum<<" "<<n1<<" "<<n2<<" ";
		no_edge=true;
		for(int j=0; j<nodes[n1].edges.size(); j++){
			if (nodes[n1].edges[j]==n2) no_edge=false;
			//if (debug) cout<<" no edge "<<no_edge<<" "<<false<<" ";
		}
		if (no_edge) {
			addEdge(n1,n2); //adds edge incriments nE++;
		    //if (debug) cout<<"Add edge "<<n1<<" "<<n2<<" "<<nE<<" "<<m<<endl;
		}
	}

    return;	//need to modify so that gen edge works with only the nkernal # of edges
}


/***********************************************************************************/
/******* Build a random graph testing each possible edge ***************************/
void graph::build_random_naive(){	
	long idum=seed;
	double Rand;
	int n1, n2;
	populate_graph();
	
	double p_E=p/((double)num_n-1.0); //edge probabiility
	
   // double dm=1.0*(double)n*p/2.0;
   // int m = (int)(dm+0.5);
	
	for (n1=0; n1<num_n-1; n1++){
		for (n2=n1+1; n2<num_n; n2++){
			Rand=ran0(&idum);
			if (Rand < p_E){ //put in edge n1,n2
				addEdge(n1,n2);
				if(debug) cout<<"add edge "<<n1<<" "<<n2<<endl;
			}
		}
	}
		
    return;	
}
void graph::build_random_naive(int nkernal){	
	long idum=seed;
	double Rand;
	int n1, n2;
	populate_graph();
	
	double p_E=p/((double)nkernal-1.0); //edge probabiility - is this right?
	
   // double dm=1.0*(double)n*p/2.0;
   // int m = (int)(dm+0.5);
	
	for (n1=0; n1<nkernal-1; n1++){
		for (n2=n1+1; n2<nkernal; n2++){
			Rand=ran0(&idum);
			if (Rand < p_E){ //put in edge n1,n2
				addEdge(n1,n2);
				if(debug) cout<<"add edge "<<n1<<" "<<n2<<endl;
			}
		}
	}
		
    return;	
}

/***********************************************************************************/
/***** Build a random graph capping the maximum connectivity on each node **********/
void graph::build_random_fixed(){
	long idum=seed;
	populate_graph();
	bool no_edge=true, size_ok=false;
	int n1, n2;
		
	unsigned long long l=0, m, mE, maxc=L;
	double dm, dl;
	dm = double(num_n)*c/2.0;
	m =(unsigned long long)(dm+0.5);
	
//	mE = num_n*(num_n-1)/2;
	
	for (int i=0; i<1000; i++) float blank=ran3(&idum);  //initialize random
	
	nE =0;
	
	while (nE<m) {
		GenEdge(idum,n1,n2);
		no_edge=true;
		size_ok=false;
		for(int j=0; j<nodes[n1].edges.size(); j++){ //is edge already in the graph
			if (nodes[n1].edges[j]==n2) no_edge=false;
		}
		if (nodes[n1].edges.size()<maxc && nodes[n2].edges.size()<maxc) size_ok=true;
		if (!no_edge && size_ok) addEdge(n1,n2); //adds edge incriments nE++;
	}
	
	return;	
}
/***********************************************************************************/
/*** for a random graph there are n(n-1)/2 possible edges, here we choose one ******/
/*** and make a  head to tail (n1,n2) pair *****************************************/
unsigned long long graph::gen_edge_number(unsigned long long maxE, long &idum){
    double dl;
    
    dl = fmod(ran3(&idum)*(double(maxE)),(double)maxE); //find an edge between 1 and maxE	
	unsigned long long l=(unsigned long long) dl;
	
	return l;		
}
void graph::get_n1_n2_from_l(unsigned long long l, int &n1, int &n2){
    int i=0, n1carry, n2carry;
    unsigned long long t = (unsigned long long)num_n-1;
    while ((t<l) && (i<num_n)){
		i++;
		t=t+(num_n-i-1);
	}
	n1=i;
	t=t-(num_n-i-1)+1;
	n2=n1+1;
	while (t<l){
		t++;
		n2++;		
    }
	if (l==0) {
	    n1=0;
	    n2=1;	
	}
	return;	
}
/***********************************************************************************/
/***** Build a random graph capping the maximum connectivity on each node **********/
/***** using the gen_edge_number routine *******************************************/
void graph::build_random_fixed_2(){
	long idum=seed;
	populate_graph();
	bool no_edge=true, size_ok=false;
	int n1, n2;
	int maxc=L;

	unsigned long long l=0, m, mE;
//	double dm, dl, ;
	double dm = double(num_n)*c/2.0;
	m =(unsigned long long)(dm+0.5);
	mE = num_n*(num_n-1)/2;
	if(debug) cout<<" initialize ran3 ";
	for (int i=0; i<1000; i++) float blank=ran3(&idum);  //initialize random
	
	nE =0;
	if(debug) cout<<" prepare to add edges ";
	while (nE<m) {
		l=gen_edge_number(mE, idum);
		get_n1_n2_from_l(l, n1, n2);
		no_edge=true;
		size_ok=false;
		for(int j=0; j<nodes[n1].edges.size(); j++){ //is edge already in the graph
			if (nodes[n1].edges[j]==n2) no_edge=false;
		}
		if (nodes[n1].edges.size()<maxc && nodes[n2].edges.size()<maxc) size_ok=true;
		if (!no_edge && size_ok) addEdge(n1,n2); //adds edge incriments nE++;
	}
	
	return;	
}
/***********************************************************************************/
/****** Scale-free random graphs       *********************************************/
/***********************************************************************************/
/***********************************************************************************/
/********************************************************************/
/** builds scale free random graph                                 **/
/**      static model  cond-mat/0312336                            **/
/********************************************************************/
double graph::prob_staticSF(int n1, int n2, double norm, double alpha) {
   double p_E;

   p_E= (pow((double)n1,-alpha)/norm)*(pow((double)n2,-alpha)/norm);
   //p_E= (pow((double)n1+1.0,-p)*pow((double)n2+1.0,-p))/(norm*norm);
      //p_E= (pow((double)n1+1.0,p)/norm)*(pow((double)n2+1.0,p)/norm);
   return p_E;
}
/********************************************************************/
/** builds scale free random graph                                 **/
/**      static model  cond-mat/0312336                            **/
/**      p is a control parameter [0,1) if p=0 graph is ER         **/
/********************************************************************/
void graph::build_static_sclfr(double alpha){
   double p_E=p/(double)num_n;
   float Rand;
   long idum=seed;     
   nE=0;
   
   double dm=(double)num_n*p/2.0;
   int m=int(dm+0.5);
   
   bool no_edge=true;
   int n1, n2;
   
   populate_graph();
   
   if(debug) cout<<" static sclfree populated "<<num_n<<" "<<m<<endl;	
   
   // initialize norm
   while (p>=1) p=p/num_n;
   double norm=0;
   for (n1=1; n1<=num_n; n1++){
	   norm +=pow((double)n1, -alpha);
	  // norm +=pow((double)n1+1.0, p);
   }
   if(debug) cout<<"norm: "<<norm<<" "<<p_E<<" ";
   if (debug) p_E=prob_staticSF(100,1,norm, alpha);
   if(debug) cout<<p_E<<" ";
  
	if(debug) cout<<" prepare to add edges ";	
	
   while (nE<m) {
		GenEdge(idum, n1, n2);
		if(nodes[n2].edges.size()<nodes[n1].edges.size()){
			int j=n2;
			n2=n1;
			n1=j;
		}
		//if (debug) cout<<idum<<" "<<n1<<" "<<n2<<" ";
		no_edge=true;
		int j=0;
		while (no_edge && j<nodes[n1].edges.size()){ //is this already an edge?
			if (n2 == nodes[n1].edges[j]) no_edge=false;
		//	if (debug) cout<<" no edge "<<no_edge<<" "<<false<<" "<<j<<" ";
			j++;
		}
		if (no_edge) {
			p_E=prob_staticSF(n1,n2,norm,alpha);
		//	if (debug) cout<<" "<<n1<<" "<<n2<<" "<<p_E<<" "<<endl;
			if (Rand < p_E) {
				addEdge(n1,n2); //adds edge incriments nE++;				
			 //   if(debug) cout<<"add edge "<<n1<<" "<<n2<<endl;
			    if (debug) cout<<"Add edge "<<n1<<" "<<n2<<" "<<nE<<" "<<m<<" "<<Rand<<" "<<p_E<<endl;
			}
		    //if (debug) cout<<"Add edge "<<n1<<" "<<n2<<" "<<nE<<" "<<m<<endl;
		}
	}

   return;
}
/********************************************************************/
/** builds power law network                                       **/
/**                                                                **/
/********************************************************************/
double graph::prob_PI(int n2) {

   double p_E;

   p_E = p*double(nodes[n2].edges.size())/(2.0 * (double) nE)/2.0;
   return p_E;
}
void graph::build_sclfr_grph(){
   double p_E=c/(double)num_n;
   float Rand;
   long idum=seed;     
 //  nE=0;
    populate_graph();
   if(debug) cout<<" populated ";	
 	if(debug) cout<<" prepare to add edges ";	  	
	for (int n1=1; n1<num_n; n1++){ //for all nodes
		for (int n2=0; n2<n1; n2++){ //for all nodes below n1 
	        Rand=ran2(&idum);  
	        if (nE ==0) p_E=1;
	        else p_E=prob_PI(n2);
	        if (Rand < p_E && nodes[n1].edges.size()<100 && nodes[n1].edges.size()<100){ //put in edge n1,n2
				addEdge(n1,n2);
				if(debug) cout<<"add edge "<<n1<<" "<<n2<<" "<<nE<<endl;
			}
	    } 			
    } 

   return;
}
/********************************************************************/
/** builds barabosi network                                        **/
/**                                                                **/
/** in a barabosi style network there are nkernal nodes placed w/  **/
/** bonds as the kernal of the graph, then nodes are added to the  **/
/** graph always with c # of edges, these are connected randomly to**/
/** previous nodes based on the probability for each previous node **/
/** which = (connectivity of node i)/total number of edges         **/
/**                                                                **/
/** This probability is generated by placing the head and tail node**/
/** of each edge into the vector prob and chosing randomly a       **/
/** reference in the vector prob as the node that is now connected **/
/** to the new node.                                               **/
/********************************************************************/
void graph::build_barabosi_network(int nkernal){ //kernal is the initial size of a random core.
	int c_curr=0, n2=0;
	int num_n_total=num_n; //the is a dummy to hold the total number of nodes
	std::vector< int > prob;  //a probabilty for each node store the 
	long idum=seed;
	bool s=true;

	//build a random core of nkernal nodes
	num_n=nkernal;
	build_random(); //this builds a random graph of nkernal number of nodes
//	build_random_naive();
//	build_random_naive(nkernal);		
	//cur_c=get_c();
	// would like to change this so that one node in the kernal is highly connected 
	//(all connected to node 0 with some probability)
//	num_n=num_n_total;
	//if (debug)cout<<"built kernal.  no. node:"<<num_n<<" no. Edges:"<<nE<<" "<<endl;


	for (int i=0; i<nodes.size(); i++){ //initialize prob vector
		for(int j=0; j<nodes[i].edges.size(); j++){  //every head and tail goes into prob vector
			n2=nodes[i].edges[j];
			if (i<n2) {
			   prob.push_back(i);
			   prob.push_back(n2);
		    }
		}
	}
	if (debug)cout<<"built kernal.  no. node:"<<num_n<<" no. Edges:"<<nE<<" "<<prob.size()<<" ";

	num_n=num_n_total;
	populate_graph();  //add the rest of the nodes to the graph
	if (debug)cout<<"populated.  no. node:"<<num_n<<" in graph"<<nodes.size()<<" c:"<<p<<endl;

    int c_stop=p;

	for (int n1=nkernal-1; n1<num_n; n1++){ //insert an edge with c connections
		c_curr=0;
		while (c_curr<c_stop){
     	  //	if(debug)cout<<" start while"<<" ";
		  // if(debug)cout<<" s:"<<s<<" fm "<<(int)fmod(ran2(&idum)*prob.size(),prob.size())<<" "<<prob.size()<<" "<<prob[(int)fmod(ran2(&idum)*prob.size(),prob.size())];
		   do {
		     int j=(int)fmod(ran2(&idum)*prob.size(),prob.size());
		     n2 = prob[j];
	       }  while (n1==n2);
		   if(debug)cout<<" n2:"<<n1<<" n2:"<<n2<<" ";
		   int j=0;
		   s=true;
		   while (s && j<nodes[n2].edges.size()){
			   if (n1==nodes[n2].edges[j]) {
				   s=false; //edge n1 n2 exists
				   if(debug)cout<<" n2:"<<n1<<" n2:"<<n2<<" exists ";
			   }
			   j++;
		   }
		   if (s) {
			   addEdge(n1,n2);
			   c_curr++;
			   if(debug)cout<<"edge added: c_curr:"<<c_curr;
		       prob.push_back(n1);
			   prob.push_back(n2);
			   if(debug)cout<<"edges added added to prob vector "<<prob.size();
		   }
		   if(debug)cout<<" return"<<endl;
 		}
	}				
   if (debug)cout<<"built kernal.  no. node:"<<num_n<<" no. Edges:"<<nE<<" "<<endl;
		
   return;
}
/***********************************************************************************/
/****** Routines for building lattices *********************************************/
/***********************************************************************************/
/***********************************************************************************/
void graph::build_triangular_lattice(){
	int n1, n2, n3, n4;
	long idum=seed;
	double Rand, p_E=p;
	
	num_n=L*L;
	if (debug) cout<<"pt lattice num_n="<<get_n()<<" "<<num_n<<" "<<nodes.size()<<endl;
	populate_graph(); //add num_n nodes to graph.
	
	for (int row=0; row<L; row++){
		for (int col=0; col<L; col++){
		   n1=row*L+col;
		   n2=n1+1;
		   n3=n1+L;
		   n4=n1+L+1;
		   if (col!=(L-1)){
		      Rand=ran3(&idum);
	   	      if (Rand < p_E){
				 addEdge(n1,n2); 
			  }
		   }	
		   if (row!=(L-1)){
		      Rand=ran3(&idum);
	   	      if (Rand < p_E){
				 addEdge(n1,n3); 
			  }			
		   }	
		   if (col!=(L-1) && ! row!=(L-1)){
		      Rand=ran3(&idum);
	   	      if (Rand < p_E){
				 addEdge(n1,n4); 
			  }		
		   }
		}
	}
	return;	
}
/***********************************************************************************/
/***********************************************************************************/
void graph::build_fcc(){
   	int num_n=L*L*L; //make an (L x L) lattice
	int n1, n2, n3, n4, n5, n6, n7;
    double Rand;
	double p_E=p;
	long idum=seed;
		
	populate_graph(); //add num_n nodes to graph.
	
	for (int layer=0; layer<L; layer++){
		for (int row=0; row<L; row++){
			for (int col=0; col<L; col++){
		        n1=layer*L*L+row*L+col;
		        n2=n1+1;
		        n3=n1+L;
		        n4=n1+L+1;
		        n5=n1+L*L; //above n1
		        n6=n1+L*L+1; // next to n5 (above n2)
		        n7=n1+L*L+L+1; //above n4
		        if (col!=(L-1)){ //not last col
		           Rand=ran3(&idum);
	   	           if (Rand < p_E){
				      addEdge(n1,n2); 
			       }
		         }	
		        if (row!=(L-1)){//not last row
		           Rand=ran3(&idum);
		           if (Rand < p_E){
                      addEdge(n1,n3); 
		           }			
		        }
			    if (col!=(L-1) && ! row!=(L-1)){//not last row or last col
		           Rand=ran3(&idum);
	   	           if (Rand < p_E){
				      addEdge(n1,n4); 
			       }		
		        }
		        if (layer!=(L-1)){ //not last layer
		           Rand=ran3(&idum);
	   	           if (Rand < p_E){
				      addEdge(n1,n5); 
			       }
			       if (col!=(L-1)){//not last col (and not last layer)
		               Rand=ran3(&idum);
		               if (Rand < p_E){
                          addEdge(n1,n6); 
		               }			
		           }
			       if (col!=(L-1) && row!=(L-1)){//not last row and last col(and not last layer)
		               Rand=ran3(&idum);
		               if (Rand < p_E){
                          addEdge(n1,n7); 
		               }			
		           }			       
			       
		         }	
		      }  	        	
		}
	}
   return;
}
/***********************************************************************************/
/***********************************************************************************/
void graph::build_square_lattice(){
	//same as build triangular but w/o the diagonals.
	int num_n=L*L; //make an (L x L) lattice
	double Rand;
	double p_E=p;
	long idum=seed;
	int n1, n2, n3;
	
	populate_graph(); //add num_n nodes to graph.
	
	for (int row=0; row<L; row++){
		
		for (int col=0; col<L; col++){
		   n1=row*L+col;
		   n2=n1+1;
		   n3=n1+L;
		   if (col!=(L-1)){
		      Rand=ran3(&idum);
	   	      if (Rand < p_E){
				 addEdge(n1,n2); 
			  }
		   }	
		   if (row!=(L-1)){
		      Rand=ran3(&idum);
	   	      if (Rand < p_E){
				 addEdge(n1,n3); 
			  }			
		   }	
		}
	}
	return;
}
/***********************************************************************************/
/***********************************************************************************/
void graph::build_cubic(){
   	int num_n=L*L*L; //make an (L x L) lattice
	int n1, n2, n3, n4, n5, n6, n7;
	long idum=seed;
	double Rand;
	double p_E=p;
		
	populate_graph(); //add num_n nodes to graph.
	
	for (int layer=0; layer<L; layer++){
		for (int row=0; row<L; row++){
			for (int col=0; col<L; col++){
		        n1=layer*L*L+row*L+col;
		        n2=n1+1;
		        n3=n1+L;
		        n4=n1+L+1;
		        n5=n1+L*L; //above n1
		        if (col!=(L-1)){ //not last col
		           Rand=ran3(&idum);
	   	           if (Rand < p_E){
				      addEdge(n1,n2); 
			       }
		         }	
		        if (row!=(L-1)){//not last row
		           Rand=ran3(&idum);
		           if (Rand < p_E){
                      addEdge(n1,n3); 
		           }			
		        }
			    if (col!=(L-1) && ! row!=(L-1)){//not last row or last col
		           Rand=ran3(&idum);
	   	           if (Rand < p_E){
				      addEdge(n1,n4); 
			       }		
		        }
		        if (layer!=(L-1)){ //not last layer
		           Rand=ran3(&idum);
	   	           if (Rand < p_E){
				      addEdge(n1,n5); 
			       }       
		        }	
		    }  	        	
		}
	}
	
	return;	
}
/***********************************************************************************/
/***********************************************************************************/
void graph::clear(){
	//cout<<nodes.size()<<" "<<get_nE();
	int i=0;
	//populate_graph();
	while (!nodes.empty()){
		while (!nodes[i].edges.empty()){
		   nodes[i].edges.pop_back(); 
		 //  cout<<i<<" removing edges"<<endl;
		// nodes[i].edges.resize();
	    }
		nodes.pop_back();
	   // while (!nodes[i].edges.empty()){
		// nodes[i].edges.pop_back(); 
		// nodes[i].edges.resize();
	 // }	
	 // cout<<i<<" not empty"<<endl;	
	  i++;
	}

	nodes.clear();
	std::vector<node>(nodes).swap(nodes);
	//nodes.swap(std::vector<node>());
	set_nE(0);
	nL=0;
	nT=0;
	Psum=0;
	Fsum=0;
	//		  cout<<i<<" clear nodes"<<endl;
	return;
}
/***********************************************************************************/
/***********************************************************************************/
/***********************************************************************************/
void graph::build_graph(bool write){    

   clear();//this should make sure than when you construct a 
             //graph it is empty. It will still contain
			 //the number of nodes
   if (debug) cout<<" build graph subroutine  g.type:"<<grphtype<<"  N:"<<num_n<<" p:"<<p<<" c:"<<c<<" seed:"<<seed<<std::endl; 
   switch(grphtype){
	  case 0:{
         if (p/((double)num_n-1)>.75) build_random_naive();
         else build_random();
         if(debug)cout<<"built random"<<" n: "<<get_n()<<" nE: "<<get_nE()<<endl;
         break;
      }
      case 1:{ //-pt
         build_triangular_lattice();
         if(debug)cout<<"built pt"<<endl;
         break;		  
	  }
	  case 2:{ //ptr3
		  build_triangular_lattice();
         if(debug)cout<<"built pt r3"<<endl;		  	
         break;
         	  }
      case 3:{ //-fcc
		  build_fcc();
         if(debug)cout<<"built fcc"<<endl;
         break;
         	  }
      case 4:{ //-cubic
		  build_cubic();
         if(debug)cout<<"built cubic"<<endl;
         break;		  
	  } 
	  case 5:{ //-sq
		  build_square_lattice();
         if(debug)cout<<"built sq"<<endl;		  
         break;		  
	  }      
	  case 6:{ //-sclf
	      if(debug)cout<<"building sclf ";	  
		  build_sclfr_grph();
          if(debug)cout<<"built sclf"<<endl;
        break;
	  }      
	  case 7:{ //-staticsclf
		  double alpha=0.5;
		  build_static_sclfr(alpha);
          if(debug)cout<<"built sclf"<<endl;		  
          break;
	  }      
	  case 8:{ //-bb
		  int nkernal=0;
		  double percent=0.00, step=1/(double)num_n;

		  while (nkernal<=10*p){ //kernal must be bigger than the number of edges added per node

			  percent=percent+step;
			  nkernal =(int) (percent*(double)num_n);
			 // if(debug)cout<<"nkernal"<<nkernal<<endl; 
		  }
		  if(debug)cout<<"nkernal"<<nkernal<<endl;
		  build_barabosi_network(nkernal);
          if(debug)cout<<"built bb"<<endl;
          break;        
	  } 	  	  	  
      case 9:{
		  build_random_fixed();
          if(debug)cout<<"built random fixed"<<endl;
		  
      break;		  
	  }      
      default:{
		  build_triangular_lattice();
         if(debug)cout<<"built pt"<<endl; 
      break;
      }
   }
  // nE=get_nE();
  // get_c();
   
   c=2.0*(double)nE/(double)num_n;
   if(debug)cout<<"built graph"<<" n: "<<get_n()<<" nE: "<<get_nE()<<endl;
   if (write) {
      std::string output="tempgrph.gml";
      write_gml(output);
   }
   return;
}


