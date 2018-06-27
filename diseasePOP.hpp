#ifndef DiseasePOP_HPP
#define DiseasePOP_HPP

//#include <cmath>

//#include <iostream>
//#include <fstream>
//#include <sstream>
//#include <string>
//#include <vector>
#include "stats.hpp"
#include "graph.hpp"
//#include "random.hpp"
//#include "diseasePOP.hpp"

#define P_STEPS 500

class dPOP{
	private:
	   void randomize_node_list(int, std::vector<int>& , long&);
	   
	public:
	   int n_immune;
	   int n_sick;
	   int n_sick_max;
	   int n_sick_mIter;
	   int n_healthy;
	   int n_dead;
	   int POP_STEPS;  //number of time steps
	   int Evolve_STEPS; //number of similar evolutions
	   double lifetime;
	   double p_Immune;
	   double p_sick;
	   double contagin;
	   double fatality;
	   double p_Recovery;	
	   bool I;  //switch for immunity and death
	  // int iteration;
	   bool randomvectors;
	   bool simple; //switch for immunity and additional random new illness;
	   bool debug;
	   dPOP (int, int, int, int, double, double, double, double, int, bool, bool);
	   dPOP ();
	   void clear ();
	   void precondition(graph& , long &);
	   void pop_evolve_lifetime(graph&, long int);
	   void pop_evolve(graph&, long int);
	   void pop_evolve(graph&, long int, iteration_stats&);

};


#endif
