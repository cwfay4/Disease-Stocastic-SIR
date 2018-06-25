/*** include file for clstrfns.cpp ***/

#ifndef _Disease_basic_utils_HPP
#define _Disease_basic_utils_HPP

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


//#define NUM_STEPS 25
//#define MASK 123459876
#define grphtypedflt 000

/** macros **/


/** types and classes**/

class boolcat{
   public:
      bool write;
      bool dsply;
      bool histo;
      bool debug_test;
      bool fulloutput;
      bool conv;
      bool style;  //this controls one or many graphs averaged over
	  int ITERCASE;	
	  int NUM_ITER;
//	  int NUM_GRAPH;
	  long lseedstart;
	  boolcat();
};

/** function prototypes **/
extern std::string outputname(graph g, dPOP pop);
extern void Disease_usage(bool verbose);
extern void pass_test(std::vector<std::string>* );
extern void pass_test(int& argc, char **argv);
extern void inputprocessing(int& , char** , graph& , dPOP& , boolcat& , int& ); 
     
#endif
