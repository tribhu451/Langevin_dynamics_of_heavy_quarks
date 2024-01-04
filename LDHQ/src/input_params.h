#pragma once
#include<iostream>
#include<fstream>
#include<sstream>
#include<string>
#include<chrono>
#include<stdio.h>
#include<string.h>

using std::cout;
using std::ofstream;
using std::endl;
using namespace std::chrono;
using std::string;
using std::cin;
using std::fstream;
using std::ios;
using std::istringstream;

typedef struct{

  int nquarks ; 
  std::string bulk_file ; 
  std::string ncoll_prof_file ; 
  std::string glauber_input_file ; 
  std::string pythia_file ; 
  double g0 ; 
  double t0 ;
  double p0 ; 
  double dt ;  
  int Do_set_etas_assuming_boost_invariance ; 
  int Do_evolve_primordial_ccbar_pair ;
  double delta_etas_gap ;  
} input_paramters ; 



class read_parameters{

 public :
  
  read_parameters();
  ~read_parameters();
  void read_parameters_from_file(input_paramters &iparam, string input_file_name);
  
};
