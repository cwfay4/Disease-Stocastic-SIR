/*** include file for clstrfns.cpp ***/

#ifndef _Disease_basic_utils_H_
#define _Disease_basic_utils_H_

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
	  boolcat();
};

/** function prototypes **/
extern void Disease_usage(bool verbose);
extern void Disease_input(int argc, char *argv[], graph g, dPOP pop, boolcat b); 
     
#endif
