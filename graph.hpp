#ifndef GRAPH_H
#define GRAPH_H

#include <cmath>

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include "graph.hpp"
#include "random.hpp"

using namespace std;

/***********************************************************************/
/**  Classes                                                         ***/
/**                                                                  ***/
/***********************************************************************/

class node{ //node class containing cluster number and a vector of edges
	public:
	    int size;   //int connectivity;
        int id;  
        int cluster;
        //int state; //integer state covered not covered etc
        bool frz; //boolean frozen or not
        double stateP;  //site state probability
        int intState;
	    std::vector<int> edges;
	   // node ();
	    node(int, int);
	    node(); //default constructor
};
/***********************************************************************/
/**    Graph class                                                   ***/
/**       contains # of edges, # of nodes, add a node or add an edge ***/
/**       and functions for creating different graph types           ***/
/**                                                                  ***/
/***********************************************************************/
class graph { 
	 private:
	    int num_n;
	    int nE;
	    int L; //L is the number of nodes in one dimension of a lattice - functions as maxc for randomfixed graph
	    //double c; //average connectivity of the graph
//	    int maxc;
	    double Psum;
	    double Fsum;
	    bool debug;
	    void addNode();  //adds a new node
	    void addEdge(int head, int tail); //adds and edge between head and tail
	    void get_ij(int n1,int root, int& i,int& j);  //get a head and tail from an edge number (only for random graphs)
	    void GenEdge(long &idum, int &r1, int &r2);  //create a random edge using ran0
	    void GenEdge3(long &idum, int &r1, int &r2);//create a random edge using ran3
	    void populate_graph(); 
	    void populate_graph(int); 	    
	    void build_random();               //build a random graph
	    void build_random(int);
	    void build_random_fixed();         //build a random graph with a fixed connectivity
	    void build_random_naive();         //build a random graph testing all edges
	    void build_random_naive(int);
	    unsigned long long gen_edge_number(unsigned long long maxE, long &idum);
	    void get_n1_n2_from_l(unsigned long long l, int &n1, int &n2);
	    void build_random_fixed_2();       //build a random graph with a fixed connectivity 
        double prob_staticSF(int n1, int n2, double norm); //node connection probability
	    void build_static_sclfr();  //static model  cond-mat/0312336
	    double prob_PI(int n2);     //probability for scale free power law network
	    void build_sclfr_grph();    //scale free power law network
	    void build_barabosi_network(int nkernal); //builds barabosi network
	    void build_triangular_lattice();   //build a triangular lattice
	    void build_fcc();                  //build a FCC lattice
	    void build_square_lattice();       //build a squarte lattice
	    void build_cubic();                //build a cubic lattice
	public:
	   // int n; //number of nodes in graph
	    long seed;
	    int grphtype; //integer used a switch to produce a graph of type
	    double p; //generator edge_prob
	    double c;  //actual connectivity of graph	    
	    int nL; //number of leaves
	    int nT; //number of triangles (incident on boundary of graph)
	    int bc;
	    std::vector<node> nodes;	    
	    graph();
	    graph(int); //initial based off the number of nodes
	    int get_n(); //get the number of nodes
	    int get_nE();//get the number of edges
	    int calc_nE();
	    int get_L();
	    double get_c();
	    bool get_debug();
        double sum_del_P(double P_old[]); 
        double calc_sum_P(); 
        void calc_sum_P_F();
	    double get_Psum(); //sum of P's from nodes
	    double get_Fsum(); //sum of P's from nodes
	    void flip();
	    void StateHistogram(double numbox, std::string ss);
	    void set_n(int);
	    void set_nE(int);
	    void set_L(int);
	    void set_bc(int);
	    void set_debug(bool);
	    void eraseEdge(int head, int tail);  //remove an edge from the graph
	    void write_gml(string);//write a graph saved in gml format
	    void write_gml_to_screen();//write a graph saved in gml format to the screen
	    void read_gml(string); //read a graph saved in gml format
	    void dsplygrph();
	    void clear(); 
	    void build_graph(bool);
};

#endif
