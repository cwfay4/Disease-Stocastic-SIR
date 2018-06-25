// Stats Class

#include <cmath>
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include "stats.hpp"

using namespace std;

Stats::Stats (){
   value=0;
   value_sq=0;
   average_value=0;
   std_error=0;
   dsamp=0;
   return;
}
Stats::Stats (double dn){
   value=0;
   value_sq=0;
   average_value=0;
   std_error=0;
   dsamp=dn;
   return;	
}
Stats::Stats (double a, double b, double c, double d, double e){
   value=a;
   value_sq=b;
   average_value=d;
   std_error=d;
   dsamp=e;
   return;	
}
void Stats::set_value(double a){
	value = a;
	return;
}
void Stats::set_average_value(double a){
	average_value=a;
	return;
}
void Stats::set_value_sq(double a){
	value_sq=a;
	return;
}
void Stats::set_std_error(double a){
	std_error=a;
	return;
}
void Stats::set_NumberSample(double a){
	dsamp=a;
	return;
}
void Stats::add_data(double datum){ //dsamp is the number of graphs being run over
   //cout<<j<<endl;
   value=value+datum; //sum of datum : should be an average
   value_sq=value_sq+value*value; //sum of square of data
   dsamp=dsamp+1; //running count of the number of data points
   average_value=value/dsamp;
   return;
}
double Stats::get_average(){
	average_value = value/dsamp;
	return average_value;
}
double Stats::get_average(double a){
	dsamp=a;
	average_value = value/dsamp;
	return average_value;	
}
int Stats::get_NumberSample(){
	int Nsample=(int)dsamp;
	return Nsample;
}
double Stats::standard_error(){
	double av_sq=value_sq/dsamp;
	double var=sqrt(abs(value_sq-(average_value*average_value))/(dsamp-1));//calculate std error;
	std_error=var;
	return var;
}
double Stats::standard_error(double dsamp1){
	dsamp=dsamp1;
	double var=sqrt(abs(value_sq-(value*value)))/sqrt(dsamp-1);//calculate std error;
	std_error=var;
	return var;
}
/**********************************************************************************/
group_stats::group_stats(){
	return;
}
group_stats::group_stats(int a){
	for (int i=0; i<a; i++){ //put num_n nodes into the graph
		addcategory();
	}
	return;
}
void group_stats::addcategory(){
	//add a new elment to stats_list
    Stats Stat1;
	stats_group.push_back(Stat1);
	
	return;	
}
void group_stats::addcategory(double a, double b, double c, double d, double e){
	//add a new elment to stats_group
    Stats Stat1(a, b, c, d, e);
	stats_group.push_back(Stat1);
	
	return;	
}
void group_stats::set_value(int a, double db){
	stats_group[a].set_value(db);
	return;
}
void group_stats::set_average_value(int a, double db){
	stats_group[a].set_value(db);
	return;
}
void group_stats::set_value_sq(int a, double db){
	stats_group[a].set_value_sq(db);
	return;
}
void group_stats::set_std_error(int a, double db){
	stats_group[a].set_std_error(db);
	return;
}
void group_stats::set_NumberSample(int a, double db){
	stats_group[a].set_NumberSample(db);
	return;
}
void group_stats::increase_element(int a, double value){
	//add value to the running sum of element [a] 
	stats_group[a].add_data(value);
    return;
}
double group_stats::get_value(int a){
	return stats_group[a].get_average();
}
double group_stats::get_value(int a, double db){
	return stats_group[a].get_average(db);
}
int group_stats::get_NumberSample(int a){
	return stats_group[a].get_NumberSample();
}
double group_stats::get_std_error(int a){
	return stats_group[a].standard_error();
}
double group_stats::get_std_error(int a, double db){
	return stats_group[a].standard_error(db);
}
//double 
void group_stats::output_elements(int stepvalue, double dsamp1, std::string oname){
	std::ofstream output(oname.c_str(), std::ofstream::app);
	   
	for (int k=0;k<stats_group.size();k++){
	   	output<<k<<" ";
       	output<<stats_group[k].get_average(dsamp1)<<" "<<stats_group[k].standard_error(dsamp1)<<" ";
	}
	output<<endl;
	output.close();
	return;
}
void group_stats::calc_stats(){ //calculates the standard error for every value in the it_stats_vector
	       for (int j=0;j<stats_group.size();j++){
	       	   stats_group[j].standard_error();
		   }  	
	return;
}
void group_stats::calc_stats(int sample_number){ //calculates the standard error for every value in the it_stats_vector
	       for (int j=0;j<stats_group.size();j++){
	       	   stats_group[j].standard_error((double)sample_number);
		   }  
	return;
}
/**********************************************************************************/
iteration_stats::iteration_stats(){
	return;
}
// a is the number of iterations b is the number of variables to keep stats on at each iteration
iteration_stats::iteration_stats(int a, int b){
//	group_stats gpstats;
  	for (int i=0; i<a; i++){ //put num_n nodes into the graph
		group_stats gpstats(b);
		it_stats_vect.push_back(gpstats);
	}
	return;   
}
void iteration_stats::add_element(int a){
//	group_stats gpstats;
  	for (int i=0; i<a; i++){ //put num_n nodes into the graph
		group_stats gpstats;
		it_stats_vect.push_back(gpstats);
	}  	
	return;
}
    	void set_value(int, double);
   	    void set_average_value(int, double);
   	    void set_value_sq(int, double);
   	    void set_std_error(int, double);
   	    void set_NumberSample(int, double);
void iteration_stats::increase_element(int stepvalue, int a, double value){
	it_stats_vect[stepvalue].stats_group[a].add_data(value);
    return;
}

    	double get_value(int);
		double get_value(int, double);
		int get_NumberSample(int);		
		double get_std_error(int);
		double get_std_error(int, double);
double iteration_stats::get_element(int stepvalue, int a){
	return it_stats_vect[stepvalue].stats_group[a].get_average();
}
void iteration_stats::output_elements(int stepvalue, int a, std::string oname){
	   std::ofstream output(oname.c_str(), std::ofstream::app);
	   
	   for (int k=0;k<it_stats_vect.size();k++){
	   	   output<<k<<" ";
	       for (int j=0;j<it_stats_vect[k].stats_group.size();j++){
	       	   output<<it_stats_vect[k].stats_group[j].get_average()<<" "<<it_stats_vect[k].stats_group[j].standard_error()<<" ";
		   }
	   	   output<<endl;
	   	   	   	
	   }
	output.close();
	return;
}
void iteration_stats::calc_stats(int sample_number){ //calculates the standard error for every value in the it_stats_vector

	   for (int k=0;k<it_stats_vect.size();k++){
	       for (int j=0;j<it_stats_vect[k].stats_group.size();j++){
	       	   it_stats_vect[k].stats_group[j].standard_error((double)sample_number);
		   }  	
	   }
	return;
}


