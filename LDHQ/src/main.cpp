#include<iostream>
#include<fstream>
#include "bulk.h"
#include "langevin.h"
#include "random.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "initialisation.h"
#include "hadronise.h"
#include "input_params.h"

/*
void sampling_check(){

  std::ofstream file ; 
  const int nq = 100000 ; 
  double xc[nq] = {0.};
  double yc[nq] = {0.};
  double etasc[nq] = {0.};
  initialisation* ic = new initialisation();
  ic->read_and_store_ncoll_profile();
  ic->sample_and_assign_initial_positions_to_heavy_quarks(nq, xc, yc, etasc);
  file.open("sampled_positions.dat");
  for(int ii=0; ii<nq; ii++){
    file << xc[ii] << "  " << yc[ii] << "  " << etasc[ii] << std::endl ; 
  }
  file.close();

  double x = 2.3 ; 
  double y = -1.1 ;

  file.open("interp_point_check.dat");
  file << x << "  " << y << "  " << ic->interpolate_and_find_ncoll(x,y) << std::endl ;
  file.close();
}
*/

int main(int argc, char **argv){

  if(argc != 1){
    std::cout << "input file required ..." << std::endl ;
    exit(1); 
  }


  std::cout << "#===================================#" << std::endl;
  std::cout << "# Langevin dynamics of heavy quarks #" << std::endl;
  std::cout << "#===================================#" << std::endl;

  input_paramters iparams ;
  read_parameters* r = new read_parameters();
  r->read_parameters_from_file(iparams, "input_parameters"); 

  // block for realistic evolution.
  initialisation* ic = new initialisation(iparams);
  ic->read_and_store_ncoll_profile();

  bulk* blk = new bulk(iparams) ;
  blk->read_bulk_info();
  
  hadronise* hadron = new hadronise();
  langevin* lgv = new langevin(blk,ic,hadron,iparams);
  lgv->evolve_heavy_quark_over_hydro_fluid() ; 
  

  // block for random number generation check //
  /*
  std::ofstream mfile;
  mfile.open("outputs/random_number_distribution.dat");
  random_gen* rnd = new random_gen();
  for(int ii=0; ii<100000; ii++){
   mfile << rnd->random_gaussian(-0.5, 3) << std::endl ;
  } 
  mfile.close();
  */


  // block for checks in simple background //
  /*
  langevin* lgv = new langevin();
  lgv->evolve_on_a_static_background_1D() ;
  lgv->evolve_on_a_static_background_3D() ;
  lgv->evolve_on_a_background_with_constant_1D_flow() ; 
  */ 

  //sampling_check();

  //delete blk ; 
  return 0 ; 
}




