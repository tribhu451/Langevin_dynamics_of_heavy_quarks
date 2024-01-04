#pragma once
#include<iostream>
#include<cmath>
#include <iomanip> 
#include "bulk.h"
#include "random.h"
#include "TProfile.h"
#include "initialisation.h"
#include "hadronise.h"
#include "input_params.h"

class langevin{

  private :
    random_gen* rand ; 
    bulk* blk ; 
    initialisation* ic ;
    hadronise* hadron; 
    input_paramters &iparam ; 
    double G0, P0, T0, dt, mass, tau0, taumax, quark_mass ;
    void setup_params();
    double get_drag_coefficient(double px, double py, double pz, double temperature); 
    void lorentz_transformation_x( double t, double x, double y, double z, 
                              double vx, double vy, double vz,
                              double &tprime, double &xprime, double &yprime, double &zprime );
    void lorentz_transformation_p( double E, double px, double py, double pz, 
                                   double vx, double vy, double vz,
                                   double &Eprime, double &pxprime, double &pyprime, double &pzprime );
    void update_position( double oldt, double oldx, double oldy, double oldz, double E, double px,
                          double py, double pz, double dt, double &nwt, double &nwx, double &nwy, double &nwz );

    void update_momentum( double oldpx, double oldpy, double oldpz,
                          double gamma, double diffusion_coeff , double rhox,
                          double rhoy, double rhoz, double dt, double &nwpx, double &nwpy, double &nwpz );

  public :
   langevin(bulk* , initialisation*, hadronise* , input_paramters &iparam_);
   langevin();
   ~langevin();
   void evolve();
   void evolve_on_a_static_background_1D();
   void evolve_on_a_static_background_3D();
   void evolve_on_a_background_with_constant_1D_flow();
   void evolve_heavy_quark_over_hydro_fluid();



};





