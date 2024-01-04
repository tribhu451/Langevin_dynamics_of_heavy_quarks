#pragma once
#include<iostream>
#include<fstream>
#include<sstream>
#include<string>
#include<cmath>
#include "random.h"
#include "input_params.h"

using std::istringstream ;

class initialisation{


  private :
    std::string ncoll_profile_filename ; 
    std::string glauber_input_filename ; 
    random_gen* rand ; 
    random_gen* rand2 ; 
    input_paramters &iparam ; 


    // these variables are required to read the ncoll profile
    // which will be helpful to sample the position of heavy quarks.
    int nx        ; 
    int ny        ;
    double xmin   ; 
    double ymin   ;
    double dx     ;
    double dy     ;
    
    double matter_eta_plateau ; 
    double matter_eta_fall ; 

    istringstream* iss ;
    char   buff[500]   ;

    double **ncollxy ; 

    inline int floor_ix(double xx){
      return (int)((xx - xmin) / dx) ; 
    }
    inline int floor_iy(double yy){
      return (int)((yy - ymin) / dy) ; 
    }

   int theta(double xx);
   double parameterized_pt_dist_200GeV(double pt);

  public :
    initialisation(input_paramters &iparam_);
    ~initialisation();
    void read_and_store_ncoll_profile();
    double interpolate_and_find_ncoll(double xx, double yy);
    void sample_and_assign_initial_positions_to_heavy_quarks(int nq, double* xc, double* yc, double* etasc);
    void sample_and_assign_initial_momenta_to_heavy_quarks(int nq, double* pxx, double* pyy, double* pzz, int* pid, int* status,
         double* xx, double* yy, double* etass);
    int get_no_of_heavy_quarks_to_be_evolved();

};








