/*** include file for build.c ***/

#ifndef _crprcfnsgraph_H_
#define _crprcfnsgraph_H_

#include <cstdio>                                     
#include <cstdlib>
#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>

#include "tp.hpp"
#include "grphfns_a.hpp"
#include "crprcfns_a.hpp"

#define cf_version "2.0.0.0"

/** macros **/


/** types and classes**/


/** function prototypes **/
extern void crprcfns_cnt_cre_dta(int n, int &ncore, int &nEcore, int gv[], 
   int g2[]);
extern void crprcfns_LRA(int n, int nE, int &nC, int size[], int gv[], 
   int g2[], int &num_leaves);
extern void crprcfns_LRA2(int n, int nE, int &nC, int size[], int gv[], 
   int g2[], int &num_leaves);  
extern void crprcfns_triprc(int n, int nE, int &nC, int size[], int gv[], 
   int g2[], int &num_leaves, int &num_tri);  
#endif
