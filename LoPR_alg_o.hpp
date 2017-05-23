/*** include file for build.c ***/

#ifndef _vcELoPR_alg_H_
#define _vcELoPR_alg_H_

#include <cstdio>                                     
#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <list>

#include "defaults.hpp"
#include "random.hpp"
#include "graph.hpp"
//#include "crprcfns_a.hpp"

#define lp_version "2.1.1.1"
/** macros **/


/** types and classes**/


/** function prototypes **/
//extern void vce_histogram(int gv[], int g2[],
//      double P[], double numbox, int n, char* s);
extern void vce_precon(graph g, long &seed);
extern void vce_randomize_node_list(int n, std::vector<int>& nlist, long &idum);
extern void vce_siteELoPRalg(graph g, long &idum, bool conv, bool rlists);
extern double vce_calc_PI_b(graph g, int n1, int n2);
extern void vce_bELoPRalg(graph g, long &idum, bool conv, bool rlists);   
extern void vce_bELoPRalg2(graph g, long &idum); 
extern void vce_listdegen(graph g, std::list< int > &degen);
extern void vce_siteDIG_weighted(graph g, long &seed, bool rlists, bool weight);
extern void vce_siteDIG(graph g, long &seed, bool rlists);
extern void vcd_random_number_prep(std::vector< long > &seed);
//control functions      
extern void vce_sb(graph g, bool write, bool dsply, bool fulloutput,
      std::string &dname, std::vector<double> &datum, long &lseedin, bool histo, 
      bool conv, bool rlists, char* filename);
extern void vce_sDig(graph g, bool write, bool dsply, bool fulloutput, bool rlists,
      std::string &dname, std::vector<double> &datum, 
      long &lseedin, bool histo, bool conv, char* filename); 
extern void vce_sDbD(graph g, bool write, bool dsply, bool fulloutput, bool rlists,
      std::string &dname, std::vector<double> &datum, 
      long &lseedin, bool histo, bool conv, char* filename);
extern void vce_sDbD_1list_many(graph g, bool write, bool dsply, bool fulloutput,
      std::string &dname, std::vector<double> &datum, 
      long &lseedin, bool histo, bool conv, char* filename);         
extern void vce_sDig_best(graph g,bool dsply, bool write, bool fulloutput, bool rlists,
       long &lseedin, bool conv, std::string &dname, char* filename);
extern void vce_sFlip(graph g,bool dsply, bool write, bool fulloutput, bool rlists,
       long &lseedin, bool conv, std::string &dname, char* filename); 
extern void vce_sDb_flip(graph g,bool dsply, bool write, bool fulloutput, bool rlists,
       long &lseedin, bool conv, std::string &dname, char* filename);
//cover greedy
extern int cvrfns_choose_maxc(graph g, std::vector<int>& nlist, int &maxc);
extern int cvrfns_greedy_maxc(graph g, std::vector<int>& nlist);
extern int cvrfns_greedy(graph g, std::vector<int>& nlist);
extern void cvf_Greedy(graph g, bool write, bool dsply, bool fulloutput,
      std::string &dname, std::vector<double> &datum, long &lseedin, bool histo, 
      bool conv, bool rlists, char* filename);
#endif
