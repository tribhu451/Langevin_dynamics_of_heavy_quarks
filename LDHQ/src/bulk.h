#pragma once
#include<iostream>
#include<fstream>
#include<sstream>
#include<string>
#include "input_params.h"

using std::istringstream ;

class bulk{

  private :

    input_paramters &iparam ; 
    std::string bulk_file_name ; 
    int ntau      ;
    int nx        ; 
    int ny        ;
    int neta      ;
    double taumin ;
    double xmin   ; 
    double ymin   ;
    double etamin ;
    double dtau   ;
    double dx     ;
    double dy     ;
    double deta   ;

    double ****temperature ; 
    double ****vx          ; 
    double ****vy          ; 
    double ****vz          ; 

    inline int floor_itau(double xx){
      return (int)((xx - taumin) / dtau) ; 
    }

    inline int floor_ix(double xx){
      return (int)((xx - xmin) / dx)     ; 
    }

    inline int floor_iy(double xx){
      return (int)((xx - ymin) / dy)     ; 
    }

    inline int floor_ieta(double xx){
      return (int)((xx - etamin) / deta) ; 
    }


   istringstream* iss ;
   char   buff[500]   ;

   // four dimensional interpolation routine
   double four_dimension_linear_interpolation(
            double* lattice_spacing, double fraction[2][4], double**** cube) ;
 

  public :

   bulk(input_paramters &iparam) ;
   ~bulk() ;
   void read_bulk_info() ;
   double get_temperature(double mtau, double mx, double my, double meta) ;
   double get_vx(double mtau, double mx, double my, double meta) ;
   double get_vy(double mtau, double mx, double my, double meta) ;
   double get_vz(double mtau, double mx, double my, double meta) ;
   double get_hydro_bulk_tau0();

};






