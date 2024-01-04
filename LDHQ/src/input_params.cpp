#include "input_params.h"

read_parameters::read_parameters(){
}

read_parameters::~read_parameters(){
}

// This functions reads the input file and sets the input parameters in a structure for global use in the code. 
void read_parameters::read_parameters_from_file(input_paramters &iparam, std::string input_file_name)
{
  std::cout << "reading input paramters ... " << std::endl ; 
  string param_name; char param_value[50];
  
  istringstream* iss;
  char   buff[1000];
  
  fstream File0;
  File0.open(input_file_name.c_str(),ios::in);
  if(!File0){
     cout<<"No input parameter file found.\nexit(1)";
     exit(1);
  }
  
  int number = 0;
  while(!File0.eof())
    {
      File0.getline(buff,850);
      
      if (!(*buff) || (*buff == '#')){
	number ++; 
	continue;
      }
      
      iss = new istringstream(buff);
      *iss >> param_name >> param_value ;
      delete iss ; 
      
      if(param_name == "nquarks" ){
	iparam.nquarks = atof(param_value);
      }

      if(param_name == "g0" ){
	iparam.g0 = atof(param_value);
      }

      if(param_name == "t0" ){
	iparam.t0 = atof(param_value);
      }

      if(param_name == "p0" ){
	iparam.p0 = atof(param_value);
      }

      if(param_name == "dt" ){
	iparam.dt = atof(param_value);
      }

      if(param_name == "delta_etas_gap" ){
	iparam.delta_etas_gap = atof(param_value);
      }
      
      if(param_name == "bulk_file"  ){
	iparam.bulk_file = param_value;
      }

      if(param_name == "ncoll_prof_file"  ){
	iparam.ncoll_prof_file = param_value;
      }

      if(param_name == "glauber_input_file"  ){
	iparam.glauber_input_file = param_value;
      }

      if(param_name == "pythia_file"  ){
	iparam.pythia_file = param_value;
      }

      if(param_name == "Do_set_etas_assuming_boost_invariance" ){
	iparam.Do_set_etas_assuming_boost_invariance = atof(param_value);
      }

      if(param_name == "Do_evolve_primordial_ccbar_pair" ){
	iparam.Do_evolve_primordial_ccbar_pair = atof(param_value);
      }

      if(param_name == "end"  ){
         break ; 
      }

    number++; 
    
    } 
  File0.close();
}




