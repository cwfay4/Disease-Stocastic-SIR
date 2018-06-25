#ifndef stats_HPP
#define stats_HPP

#include <cmath>
#include <vector>
#include <iostream>
#include <fstream>
#include <string>

// stats collects a variable and its square from the variable, its square and the number of samples the variable is averaged 
// over you can calculate the standard error.
class Stats{
   private:
   	    double value;  //running summation of the value
        double value_sq;  //running summation of the value squared
        double std_error;  //std Error = sqrt [(<x^2>-<x>^2)/N]
        double dsamp;  //N the number of samples
        double average_value;
   public:
   	    void set_value(double);
   	    void set_average_value(double);
   	    void set_value_sq(double);
   	    void set_std_error(double);
   	    void set_NumberSample(double);
        void add_data(double);
        double get_average();
        double get_average(double);
        int get_NumberSample();
        double standard_error();
        double standard_error(double);      
        Stats();
        Stats(double);
        Stats(double, double, double, double, double);
};
// a class creating a vector where each element in the vector is a Stats class item
class group_stats{
	public:	
		std::vector<Stats> stats_group;
		void addcategory();
		void addcategory(double, double, double, double, double);
   	    void set_value(int, double);
   	    void set_average_value(int, double);
   	    void set_value_sq(int, double);
   	    void set_std_error(int, double);
   	    void set_NumberSample(int, double);				
		void increase_element(int, double);
		double get_value(int);
		double get_value(int, double);
		int get_NumberSample(int);		
		double get_std_error(int);
		double get_std_error(int, double);		
		void output_elements(int, double, std::string);
		void calc_stats();		
		void calc_stats(int);
		group_stats();
		group_stats(int);
};


class iteration_stats{
    public:
    	std::vector<group_stats> it_stats_vect;
    	iteration_stats();
    	iteration_stats(int, int);
    	void add_element(int);
    	void set_value(int, double);
   	    void set_average_value(int, double);
   	    void set_value_sq(int, double);
   	    void set_std_error(int, double);
   	    void set_NumberSample(int, double);	
    	void increase_element(int, int, double);
    	double get_element();
    	double get_element(int);
		double get_element(int, int);
		int get_NumberSample(int);		
		double get_std_error(int);
		double get_std_error(int, double);
    	void output_elements(int, int, std::string);
    	void calc_stats(int);
};


#endif
