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
#include <iostream>
#include <fstream>
#include "LoPR_basic_utils.hpp"
//#include "tp.hpp"
#include "graph.hpp"
#include "LoPR_alg_o.hpp"

using namespace std;

void LoPR_usage(int verbose) {
  cout<<"LoPR_a_ma : LoPR array averaged over multiple samples \n";
  cout<<"Usage: [-h][-rgml <filename>][GRAPHTYPE] n c nsamples"<<endl;
  if (!verbose)
    cout<<"Type '%s -h' for long help\n";
  else
    cout<<endl<<"Revision "<<tp_version<<" Charles Fay 2006"<<endl;
    cout<<"    grphfns_version "<<gf_version<<" Charles Fay 2006"<<endl;
    cout<<"    LoPRfns_version "<<lp_version<<" Charles Fay 2006"<<endl<<endl;
    cout<<"   -h: help \n"<<endl;
    cout<<"GRAPH OPTIONS:"<<endl;
    cout<<"   -rnd: random graph using ran2 from numerical methods \n";
    cout<<"   -pt: diluted planar triangular lattice using ran2 \n"; 
    cout<<"   -ptr3: diluted planar triangular lattice using ran3\n";    
    cout<<"   -fcc: diluted fcc lattice using ran2 \n";
    cout<<"   -cubic: diluted cubic lattice using ran2 \n";
    cout<<"   -sq: diluted square lattice using ran2 \n";    
    cout<<"   -sclf: scale free graph using ran2 (barabasi network?) \n";
    cout<<"   -rnd2: random graph sampling every edge using ran2 \n";
    cout<<"   -rnd3: using same procedure as the rf algorithm \n"; 
    cout<<"   -bb [m]: barabasi network, m is integer number of nodes in kernal \n";
    cout<<"            (c should be an integer as well)\n"; 
    cout<<"   -rf [m]: a random graph, connectivies are controlled so that c<=m \n";
    cout<<"   -rf_t2 [m]: a random graph, connectivies are controlled so that c<=m, different generator \n";
    cout<<"   -rf_64 [m]: a random graph, connectivies are controlled so that c<=m \n";
    cout<<"               optimized for a 64bit turion processor. \n";
    cout<<"   -staticSF: static scale free graph of the type from cond-mat/0312336 \n";
    cout<<"   -Ksat [m]: Ksat with m variables per clause,  and c clauses (not finished)\n"; 
    cout<<endl;   
    cout<<"   -bc0 : free boundary conditions (default)"<<endl;    
    cout<<"   -bc1 : periodic one one side (only option for 2-d graphs) \n";
    cout<<"   -bc2 : periodic on two sides \n"; 
    cout<<"LoPR OPTIONS:"<<endl;   
    cout<<"   -sb: site LoPR then bond LoPR \n";
    cout<<"   -sb1: site LoPR then bond LoPR with one node list\n";    
    cout<<"   -sDig: site LoPR then site based DIG \n"; 
    cout<<"   -sDbD: site LoPR then site based DIG, bond LoPR then site based DIG\n";
    cout<<"   -sDbD_1: site LoPR then site based DIG, bond LoPR then site based DIG\n"; 
    cout<<"   -sDbD_m: site LoPR then site based DIG, bond LoPR then site based DIG\n";
    cout<<"   -greedy: greedy cover by selection of a random node \n";
    cout<<"   -Gmaxc: greedy cover by selection of the node with largest current connectivity\n";             
    cout<<"GENERAL OPTIONS:"<<endl;
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
    cout<<"Output is of the form: \n";
    cout<<"    n c seed1 seed2 seed3 c_actual nE \n";
    cout<<"    sELoPR_covered_Fraction FrozenFraction \n";
    cout<<"    bELoPR_covered_Fraction FrozenFraction \n";
    cout<<endl;
    cout<<"Variables: \n";
    cout<<"    n: number of nodes \n";
    cout<<"    c: connectivity \n";
    cout<<"    nsamples: number of graphs to sample \n";
    cout<<"    c_actual: connectivity of graph \n";    
    cout<<"    nE: number of Edges in graph \n";
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
    cout<<"    Boundary conditions: ";
    if (bcdflt==0) cout<<"free boundaries"<<endl;
    else if (bcdflt==1) cout<<"periodic on one boundary"<<endl;
    else if (bcdflt==2) cout<<"periodic on two boundaries"<<endl;
    cout<<"    Full output: ";
    if (fulloutdflt) cout<<"TRUE"<<endl;
    else cout<<"FALSE"<<endl;
    cout<<"    Write: ";
    if (writedflt) cout<<"TRUE"<<endl;
    else cout<<"FALSE"<<endl;
    cout<<"    Histo: ";
    if (clstrdflt) cout<<"TRUE"<<endl;
    else cout<<"FALSE"<<endl;
    cout<<"    Display: ";
    if (dsplydflt) cout<<"TRUE"<<endl;
    else cout<<"FALSE"<<endl;
    cout<<"    Convergence output: ";
    if (convdflt) cout<<"TRUE"<<endl;
    else cout<<"FALSE"<<endl;
    cout<<"    Convergence Bound: "<<BND;    
    cout<<"    Max size array: "<<SIZE<<endl;
    cout<<endl; 
    
  exit(1);
}
