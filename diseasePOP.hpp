#ifndef DiseasePOP_H
#define DiseasePOP_H

#include <cmath>

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include "graph.hpp"
#include "random.hpp"
#include "diseasePOP.hpp"

#define POP_STEPS 500

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
	   double lifetime;
	   double p_Immune;
	   double p_sick;
	   double contagin;
	   double fatality;	
	   bool I;  //switch for immunity and death
	   bool randomvectors;
	   bool simple; //switch for immunity and additional random new illness;
	   dPOP (int, int, int, int, double, double, double, double);
	   dPOP ();
	   void clear ();
	   void precondition(graph& , long &);
	   void pop_evolve(graph&, long int);
};


#endif
