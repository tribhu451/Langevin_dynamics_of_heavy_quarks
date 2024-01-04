#include "langevin.h"

langevin::langevin(bulk* blk_, initialisation* ic_, hadronise* hadron_, input_paramters &iparam_ ):iparam(iparam_){
   
  ic = ic_ ; 
  blk  = blk_ ; 
  hadron = hadron_;
  setup_params();
  tau0 = blk->get_hydro_bulk_tau0() ; // in fm 
}

/*
langevin::langevin(input_paramters &iparam_){
  tau0 = 1.0 ; // in fm 
  setup_params();
}
*/

void langevin::setup_params(){
  rand = new random_gen();
  G0 = iparam.g0 ;
  P0 = iparam.p0 ;
  T0 = iparam.t0 ;
  dt = iparam.dt ;    // in fm 
  quark_mass = 1.5 ;  // in GeV
  taumax = 30.0 ;     // in fm
} 


langevin::~langevin(){
}


// general lorentz transformation when velocity direction is arbitrary. 
void langevin::lorentz_transformation_x( double t, double x, double y, double z, 
                                         double vx, double vy, double vz,
                                         double &tprime, double &xprime,
                                         double &yprime, double &zprime ){
  
  double v2    = vx * vx + vy * vy + vz * vz ; 
  double gamma = 1.0 / sqrt(1.0 - v2 ) ; 
  double vdotr = vx * x + vy * y + vz * z ; 
  
  tprime = gamma * t - gamma * vdotr ; 
  xprime = x - gamma * vx * t  + ( gamma - 1 ) * vx / v2 * vdotr ;  
  yprime = y - gamma * vy * t  + ( gamma - 1 ) * vy / v2 * vdotr ;  
  zprime = z - gamma * vz * t  + ( gamma - 1 ) * vz / v2 * vdotr ;  
}


void langevin::lorentz_transformation_p( double E, double px, double py, double pz, 
                                         double vx, double vy, double vz,
                                         double &Eprime, double &pxprime,
                                         double &pyprime, double &pzprime ){
  
  double v2    = vx * vx + vy * vy + vz * vz ; 
  double gamma = 1.0 / sqrt(1.0 - v2 ) ; 
  double vdotp = vx * px + vy * py + vz * pz ; 
  
  pxprime = px - gamma * vx * E  + ( gamma - 1 ) * vx / v2 * vdotp ;  
  pyprime = py - gamma * vy * E  + ( gamma - 1 ) * vy / v2 * vdotp ;  
  pzprime = pz - gamma * vz * E  + ( gamma - 1 ) * vz / v2 * vdotp ;  
  Eprime = gamma * E - gamma * vdotp ; 
}


void langevin::update_position( double oldt, double oldx, double oldy, double oldz,
                                double E, double px, double py, double pz, double dt,
                                double &nwt, double &nwx, double &nwy, double &nwz ){
  
  nwt = oldt + dt ; 
  nwx = oldx + ( px / E ) * dt ; 
  nwy = oldy + ( py / E ) * dt ; 
  nwz = oldz + ( pz / E ) * dt ; 
}


void langevin::update_momentum( double oldpx, double oldpy, double oldpz,
				double gamma, double diffusion_coeff , double rhox,
				double rhoy, double rhoz, double dt,
				double &nwpx, double &nwpy, double &nwpz ){
  
  nwpx = oldpx - gamma * oldpx * dt + rhox * sqrt( 2. * diffusion_coeff * dt ) ;
  nwpy = oldpy - gamma * oldpy * dt + rhoy * sqrt( 2. * diffusion_coeff * dt ) ;
  nwpz = oldpz - gamma * oldpz * dt + rhoz * sqrt( 2. * diffusion_coeff * dt ) ;
}


double langevin::get_drag_coefficient(double px, double py, double pz, double temperature){
  return G0 * pow(sqrt(px * px + py * py + pz * pz), P0) * temperature * pow(temperature/quark_mass, T0) ;
} 

/*
  void langevin::evolve(){
  // In this function, every local rest frame variable is denoted with "star"
  // wheareas all lab frame variables are with as usual notation.
  
  // All the calculations and variables of heavy quarks are in cartesian coordinate. 
  
  double tau, etas ; // heavy quark position in Lab frame(milne) 
  double t, x, y, z; // heavy quark position in Lab frame(cartesian)

  tau  = tau0 ; 
  etas = 0. ; 
  x = 1.3 ; 
  y = 0.3 ;
 
  // Heavy quark momentum in Lab frame
  double E, px, py, pz ; 
  px = 0 ; 
  py = 0 ; 
  pz = 0 ;

  for(int ii=0; ii<10; ii++){

    // Background information. Velocities are in cartesian coordinate. 
    double T, vx, vy, vz ; 
    T  = blk->get_temperature(tau, x, y, etas);
    vx = blk->get_vx(tau, x, y, etas);
    vy = blk->get_vy(tau, x, y, etas);
    vz = blk->get_vz(tau, x, y, etas);
    double lorentz_gamma = 1. / sqrt( 1. - vx * vx - vy * vy - vz * vz ) ; 
    std::cout << "🌟 T = " << T << ", vx = " << vx << ", vy = " << vy << ", vz = " << vz << std::endl ; 

    if(T < 0.150) break ; 


    t = tau * cosh(etas) ; 
    z = tau * sinh(etas) ; 
    E  = sqrt( px * px + py * py + pz * pz + quark_mass * quark_mass ) ; 

    // Transormation from Lab frame to local rest frame of the fluid.
    double tstar=0., xstar=0., ystar=0., zstar=0., Estar=0., pxstar=0., pystar=0., pzstar=0. ; 
    lorentz_transformation_x( t, x, y, z, vx, vy, vz, tstar, xstar, ystar, zstar );
    lorentz_transformation_p( E, px, py, pz, vx, vy, vz, Estar, pxstar, pystar, pzstar );

    // evolve the heavy quark position and momentum in local rest frame.
    double dtstar = 1.0 / lorentz_gamma * dt ;
    update_position( tstar, xstar, ystar, zstar, Estar, pxstar,
                        pystar, pzstar, dtstar, tstar, xstar, ystar, zstar );
    double drag_coeff      = get_drag_coefficient( pxstar,  pystar, pzstar, T);
    double diffusion_coeff = 2 * drag_coeff * Estar * T ; 
    double rhox       = rand->random_gaussian(0.,1.);
    double rhoy       = rand->random_gaussian(0.,1.);
    double rhoz       = rand->random_gaussian(0.,1.);
    update_momentum( pxstar, pystar, pzstar, drag_coeff, diffusion_coeff , rhox,
                        rhoy, rhoz, dtstar, pxstar, pystar, pzstar );


    // Transormation from local rest frame to Lab frame of the fluid.
    lorentz_transformation_x( tstar, xstar, ystar, zstar, -vx, -vy, -vz, t, x, y, z );
    lorentz_transformation_p( Estar, pxstar, pystar, pzstar, -vx, -vy, -vz, E, px, py, pz );

    tau  = sqrt(t*t - z*z);
    etas = atanh(z/t);


    std::cout << "tau  =  "  << sqrt(t*t - z*z) << ", x = " << x 
              << ", y = " << y << ", etas = " << atanh(z/t) << std::endl ; 
    std::cout << "E  =  "  << E << ", px = " << px 
              << ", py = " << py << ", pz = " << pz << " ... " << "\n" << std::endl ;
 
  } // for loop of ii
}

*/


void langevin::evolve_on_a_static_background_1D(){
  // All the calculations and variables of heavy quarks are in cartesian coordinate. 
  double bin_start  = -4   ; 
  double bin_end    = 4  ;
  int    no_of_bins = 80 ; // always even
  double bin_width  = 
    ( bin_end - bin_start ) / no_of_bins ; 
  
  double hist_dist[no_of_bins] ; 
  for(int ii=0; ii<no_of_bins; ii++)
    hist_dist[ii] = 0. ; 
  
  
  const int Nheavyquarks = 50000  ; 
  int Ntimestep    = 1000 ; 
  
  double *tau  = new double[Nheavyquarks];
  double *t    = new double[Nheavyquarks];
  double *etas = new double[Nheavyquarks];
  double *z    = new double[Nheavyquarks];
  double *x    = new double[Nheavyquarks];
  double *y    = new double[Nheavyquarks];
  double *px   = new double[Nheavyquarks];
  double *py   = new double[Nheavyquarks];
  double *pz   = new double[Nheavyquarks];
  
  for(int ii=0; ii<Nheavyquarks; ii++){
    tau[ii]  = tau0 ; 
    etas[ii] = 0. ; 
    x[ii]    = 0. ; 
    y[ii]    = 0. ; 
    t[ii]    = tau[ii] * cosh(etas[ii]) ; 
    z[ii]    = tau[ii] * sinh(etas[ii]) ; 
    px[ii]   = 0.0 ; 
    py[ii]   = 0.0 ; 
    pz[ii]   = 0.0 ; 
  }
  
  // Background information. Velocities are in cartesian coordinate. 
  double T ; 
  T  = 0.240 ; // in GeV
  
#pragma omp parallel
  {  
#pragma omp for schedule(guided)
    for(int ii=0; ii<Nheavyquarks; ii++){ // loop over quarks
      for(int jj=0; jj<Ntimestep; jj++){ // loop over time step
	double E  = sqrt( px[ii] * px[ii] + py[ii] * py[ii] + pz[ii] * pz[ii] + quark_mass * quark_mass ) ;
	
	// evolve the heavy quark position in local rest frame.
	update_position( t[ii], x[ii], y[ii], z[ii], E, px[ii],
			 py[ii], pz[ii], dt, t[ii], x[ii], y[ii], z[ii] );
	
	// the white noise and drag coefficient
	double rhox       = rand->random_gaussian(0.,1.);
	double rhoy       = 0.;
	double rhoz       = 0.;
	double drag_coeff  = get_drag_coefficient(px[ii],py[ii], pz[ii],T) ; 
        
	// calculate the diffusion coefficient D(p+dp) at next timestep //
	double temp_diffusion_coeff  =  drag_coeff * E * T ;
	double temp_px, temp_py, temp_pz, temp_e ;
	update_momentum( px[ii], py[ii], pz[ii], drag_coeff, temp_diffusion_coeff,
			 rhox, rhoy, rhoz, dt, temp_px, temp_py, temp_pz );
	temp_e = sqrt(temp_px*temp_px + temp_py*temp_py + temp_pz*temp_pz + quark_mass * quark_mass);
	double temp_drag_coeff = get_drag_coefficient(temp_px,temp_py,temp_pz,T) ;
	double diffusion_coeff  =  temp_drag_coeff * temp_e * T ; // GeV^2 * fm^{-1} 
	
	// update the momentum with post point realization //
	update_momentum( px[ii], py[ii], pz[ii], drag_coeff, diffusion_coeff,
			 rhox, rhoy, rhoz, dt, px[ii], py[ii], pz[ii] );
	
	tau[ii]  = sqrt(t[ii]*t[ii] - z[ii]*z[ii]);
	etas[ii] = atanh(z[ii]/t[ii]);
      } // loop over timestep
    } // loop over quarks
  } // omp block
  
  
  for(int ii=0; ii<Nheavyquarks; ii++ ){
    double Q  = px[ii]  ;
    int fill_bin_no = ( Q - bin_start ) / bin_width ;
    if( fill_bin_no < no_of_bins && fill_bin_no > 0 ){
      hist_dist[fill_bin_no] += 1 ; 
    }
  }
  
  std::ofstream out_file;
  std::stringstream output_filename;
  output_filename.str("");
  output_filename << "outputs/quark_distribution_in_1D_static_background_of_temperature_" ;
  output_filename << T ;
  output_filename << ".dat";
  
  out_file.open(output_filename.str().c_str(), std::ios::out );
  for(int ii=0; ii<no_of_bins; ii++){
    out_file <<  ii  << "  " << bin_start + (ii) * bin_width + (bin_width/2.) << "  " 
             << "  " << hist_dist[ii] / Nheavyquarks
             << "  " << ( hist_dist[ii] / Nheavyquarks ) / bin_width  << std::endl ;  
  }
  out_file.close();
  
  delete tau ;
  delete etas;
  delete t ; 
  delete z ;
  delete x ;
  delete y ; 
  delete px ; 
  delete py ; 
  delete pz ; 
}


void langevin::evolve_on_a_static_background_3D(){
  // All the calculations and variables of heavy quarks are in cartesian coordinate. 
  double bin_start  = -4   ; 
  double bin_end    = 4  ;
  int    no_of_bins = 80 ; // always even
  double bin_width  = 
    ( bin_end - bin_start ) / no_of_bins ; 
  
  double hist_dist[no_of_bins] ; 
  for(int ii=0; ii<no_of_bins; ii++)
    hist_dist[ii] = 0. ; 
  
  
  const int Nheavyquarks = 50000  ; 
  int Ntimestep    = 1000 ; 
  
  double *tau  = new double[Nheavyquarks];
  double *t    = new double[Nheavyquarks];
  double *etas = new double[Nheavyquarks];
  double *z    = new double[Nheavyquarks];
  double *x    = new double[Nheavyquarks];
  double *y    = new double[Nheavyquarks];
  double *px   = new double[Nheavyquarks];
  double *py   = new double[Nheavyquarks];
  double *pz   = new double[Nheavyquarks];
 
  for(int ii=0; ii<Nheavyquarks; ii++){
    tau[ii]  = tau0 ; 
    etas[ii] = 0. ; 
    x[ii]    = 0. ; 
    y[ii]    = 0. ; 
    t[ii]    = tau[ii] * cosh(etas[ii]) ; 
    z[ii]    = tau[ii] * sinh(etas[ii]) ; 
    px[ii]   = 0.0 ; 
    py[ii]   = 0.0 ; 
    pz[ii]   = 0.0 ; 
  }

  // Background information. Velocities are in cartesian coordinate. 
  double T ; 
  T  = 0.240 ; // in GeV
  
#pragma omp parallel
  {  
#pragma omp for schedule(guided)
    for(int ii=0; ii<Nheavyquarks; ii++){ // loop over quarks
      for(int jj=0; jj<Ntimestep; jj++){ // loop over time step
	double E  = sqrt( px[ii] * px[ii] + py[ii] * py[ii] + pz[ii] * pz[ii] + quark_mass * quark_mass ) ;
	
	// evolve the heavy quark position in local rest frame.
	update_position( t[ii], x[ii], y[ii], z[ii], E, px[ii],
			 py[ii], pz[ii], dt, t[ii], x[ii], y[ii], z[ii] );
	
	// the white noise and drag coefficient
	double rhox       = rand->random_gaussian(0.,1.);
	double rhoy       = rand->random_gaussian(0.,1.);
	double rhoz       = rand->random_gaussian(0.,1.);
	double drag_coeff  = get_drag_coefficient(px[ii],py[ii], pz[ii],T) ; 
        
	// calculate the diffusion coefficient D(p+dp) at next timestep //
	double temp_diffusion_coeff  =  drag_coeff * E * T ;
	double temp_px, temp_py, temp_pz, temp_e ;
	update_momentum( px[ii], py[ii], pz[ii], drag_coeff, temp_diffusion_coeff,
			 rhox, rhoy, rhoz, dt, temp_px, temp_py, temp_pz );
	temp_e = sqrt(temp_px*temp_px + temp_py*temp_py + temp_pz*temp_pz + quark_mass * quark_mass);
	double temp_drag_coeff = get_drag_coefficient(temp_px,temp_py,temp_pz,T) ;
	double diffusion_coeff  =  temp_drag_coeff * temp_e * T ; // GeV^2 * fm^{-1} 
	
	// update the momentum with post point realization //
	update_momentum( px[ii], py[ii], pz[ii], drag_coeff, diffusion_coeff,
			 rhox, rhoy, rhoz, dt, px[ii], py[ii], pz[ii] );
	
	tau[ii]  = sqrt(t[ii]*t[ii] - z[ii]*z[ii]);
	etas[ii] = atanh(z[ii]/t[ii]);
      } // loop over timestep
    } // loop over quarks
  } // omp block
  
  
  for(int ii=0; ii<Nheavyquarks; ii++ ){
    double Q  = sqrt( px[ii] * px[ii] + py[ii] * py[ii] + pz[ii] * pz[ii] ) ;
    int fill_bin_no = ( Q - bin_start ) / bin_width ;
    if( fill_bin_no < no_of_bins && fill_bin_no > 0 ){
      hist_dist[fill_bin_no] += 1 ; 
    }
  }
  
  std::ofstream out_file;
  std::stringstream output_filename;
  output_filename.str("");
  output_filename << "outputs/quark_distribution_in_3D_static_background_of_temperature_" ;
  output_filename << T ;
  output_filename << ".dat";

  out_file.open(output_filename.str().c_str(), std::ios::out );
  for(int ii=0; ii<no_of_bins; ii++){
    out_file <<  ii  << "  " << bin_start + (ii) * bin_width + (bin_width/2.) << "  " 
             << "  " << hist_dist[ii] / Nheavyquarks
             << "  " << ( hist_dist[ii] / Nheavyquarks ) / bin_width  << std::endl ;  
  }
  out_file.close();

   
  delete tau ;
  delete etas;
  delete t ; 
  delete z ;
  delete x ;
  delete y ; 
  delete px ; 
  delete py ; 
  delete pz ; 
}


void langevin::evolve_on_a_background_with_constant_1D_flow(){
  // All the calculations and variables of heavy quarks are in cartesian coordinate. 
  double bin_start  = 0   ; 
  double bin_end    = 10  ;
  int    no_of_bins = 40 ; // always even
  double bin_width  = 
    ( bin_end - bin_start ) / no_of_bins ; 
  
  double hist_dist[no_of_bins] ; 
  for(int ii=0; ii<no_of_bins; ii++)
    hist_dist[ii] = 0. ; 


  const int Nheavyquarks = iparam.nquarks  ; 
  int Ntimestep    = 2000 ; 
  std::cout << "Heavy quarks to evolve = " << Nheavyquarks << std::endl ; 

  double *tau  = new double[Nheavyquarks];
  double *t    = new double[Nheavyquarks];
  double *etas = new double[Nheavyquarks];
  double *z    = new double[Nheavyquarks];
  double *x    = new double[Nheavyquarks];
  double *y    = new double[Nheavyquarks];
  double *px   = new double[Nheavyquarks];
  double *py   = new double[Nheavyquarks];
  double *pz   = new double[Nheavyquarks];
  double *E    = new double[Nheavyquarks];

 
  for(int ii=0; ii<Nheavyquarks; ii++){
    tau[ii]  = 0 ; 
    etas[ii] = 0. ; 
    x[ii]    = 0. ; 
    y[ii]    = 0. ; 
    t[ii]    = 0 ; 
    z[ii]    = 0 ; 
    px[ii]   = 0.0 ; 
    py[ii]   = 0.0 ; 
    pz[ii]   = 0.0 ; 
    E[ii] = sqrt( px[ii] * px[ii] + py[ii] * py[ii] 
      + pz[ii] * pz[ii] + quark_mass * quark_mass ) ;
  }

  // Background information. Velocities are in cartesian coordinate. 
  double T = 0.180 ; 
  double vx = 0.9 ;
  double vy = 0. ; 
  double vz = 0. ; 

#pragma omp parallel
  {  
#pragma omp for schedule(guided)
    for(int ii=0; ii<Nheavyquarks; ii++){ // loop over quarks
      for(int jj=0; jj<Ntimestep; jj++){ // loop over time step
	
	// Transormation from Lab frame to local rest frame of the fluid.
	double tstar=0., xstar=0., ystar=0., zstar=0., Estar=0., pxstar=0., pystar=0., pzstar=0. ; 
	lorentz_transformation_x( t[ii], x[ii], y[ii], z[ii], vx, vy, vz, tstar, xstar, ystar, zstar );
	lorentz_transformation_p( E[ii], px[ii], py[ii], pz[ii], vx, vy, vz, Estar, pxstar, pystar, pzstar );
	
	//evolve the heavy quark position and momentum in local rest frame.
	//double lorentz_gamma = 1. / sqrt( 1. - vx * vx - vy * vy - vz * vz ) ; 
	//double dtstar = (1 / lorentz_gamma) * dt ;
	double dtstar = (Estar / E[ii]) * dt ;
	update_position( tstar, xstar, ystar, zstar, Estar, pxstar,
			 pystar, pzstar, dtstar, tstar, xstar, ystar, zstar );
	double rhox = rand->random_gaussian(0.,1.);
	double rhoy = 0.;
	double rhoz = 0.;
	double drag_coeff = get_drag_coefficient(pxstar, pystar, pzstar, T);
	
	double temp_diffusion_coeff = drag_coeff * Estar * T ; 
	double temp_pxstar, temp_pystar, temp_pzstar, temp_estar ; 
	update_momentum(pxstar, pystar, pzstar, drag_coeff, temp_diffusion_coeff , rhox,
                        rhoy, rhoz, dtstar, temp_pxstar, temp_pystar, temp_pzstar );
	temp_estar = sqrt(temp_pxstar*temp_pxstar + temp_pystar*temp_pystar + 
			  temp_pzstar*temp_pzstar + quark_mass*quark_mass);
	
	double temp_drag_coeff = get_drag_coefficient(temp_pxstar,temp_pystar,temp_pzstar,T) ;
	double diffusion_coeff  =  temp_drag_coeff * temp_estar * T ; 
	
	// update the momentum with post point realization //
	update_momentum(pxstar, pystar, pzstar, drag_coeff, diffusion_coeff,
			rhox, rhoy, rhoz, dtstar, pxstar, pystar, pzstar);
	Estar = sqrt( pxstar*pxstar + pystar*pystar + 
		      pzstar*pzstar + quark_mass*quark_mass);
	
	
	// Transormation from local rest frame to Lab frame of the fluid.
         lorentz_transformation_x( tstar, xstar, ystar, zstar, -vx, -vy, -vz, t[ii], x[ii], y[ii], z[ii] );
         lorentz_transformation_p( Estar, pxstar, pystar, pzstar, -vx, -vy, -vz, E[ii], px[ii], py[ii], pz[ii] );
	 
         tau[ii]  = sqrt(t[ii]*t[ii] - z[ii]*z[ii]);
         etas[ii] = atanh(z[ii]/t[ii]);
	 
      } // loop over timestep
        // lorentz_transformation_p( E[ii], px[ii], py[ii], pz[ii], vx, vy, vz, E[ii], px[ii], py[ii], pz[ii] );  
    } // loop over quarks
  } // omp block
  
  
  
  for(int ii=0; ii<Nheavyquarks; ii++ ){
    double Q = px[ii] ;
    int fill_bin_no = ( Q - bin_start ) / bin_width ;
    if( fill_bin_no < no_of_bins && fill_bin_no > 0 ){
      hist_dist[fill_bin_no] += 1 ; 
    }
  }
 
  std::ofstream out_file;
  std::stringstream output_filename;
  output_filename.str("");
  output_filename << "outputs/quark_distribution_on_1D_evolving_background_of_temperature_" ;
  output_filename << T ;
  output_filename << ".dat";

  out_file.open(output_filename.str().c_str(), std::ios::out );
  for(int ii=0; ii<no_of_bins; ii++){
    out_file <<  ii  << "  " << bin_start + (ii) * bin_width + (bin_width/2.) << "  " 
             << "  " << hist_dist[ii] / Nheavyquarks
             << "  " << ( hist_dist[ii] / Nheavyquarks ) / bin_width  << std::endl ;  
  }
  out_file.close();
   
  delete tau ;
  delete etas;
  delete t ; 
  delete z ;
  delete x ;
  delete y ; 
  delete px ; 
  delete py ; 
  delete pz ; 
  delete E ; 
}



void langevin::evolve_heavy_quark_over_hydro_fluid(){
  // In this function, every local rest frame variable is denoted with "star"
  // wheareas all lab frame variables are with as usual notation.
  // All the calculations and variables of heavy quarks are in cartesian coordinate. 
  
  std::cout << "Langevin dynamics of heavy" 
            << " quarks on 3+1D hydro background ... " 
            << std::endl ; 
  double T, vx, vy, vz ; 
  int Nheavyquarks = iparam.nquarks  ; 
  std::cout << "No. of Charm quarks to evolve = " << Nheavyquarks << std::endl ; 
  std::cout << "g0 = " << G0 << std::endl ; 
  std::cout << "t0 = " << T0 << std::endl ; 
  std::cout << "p0 = " << P0 << std::endl ; 
  std::cout << "dt = " << dt << std::endl ; 

  // when taking input from pythia
  Nheavyquarks = ic->get_no_of_heavy_quarks_to_be_evolved() ; 
  std::cout << "No. of Charm quarks to evolve = " << Nheavyquarks << std::endl ; 
  
  double *tau  = new double[Nheavyquarks];
  double *t    = new double[Nheavyquarks];
  double *etas = new double[Nheavyquarks];
  double *z    = new double[Nheavyquarks];
  double *x    = new double[Nheavyquarks];
  double *y    = new double[Nheavyquarks];
  double *px   = new double[Nheavyquarks];
  double *py   = new double[Nheavyquarks];
  double *pz   = new double[Nheavyquarks];
  double *E    = new double[Nheavyquarks];
  int    *PID  = new int[Nheavyquarks];
  int    *status = new int[Nheavyquarks];
 
  // quantities to be asigned
  for(int ii=0; ii<Nheavyquarks; ii++){
    PID[ii]  = 0 ; 
    etas[ii] = 0.0 ; 
    x[ii]    = 0.0 ; 
    y[ii]    = 0.0 ; 
    px[ii]   = 0.0 ; 
    py[ii]   = 0.0 ; 
    pz[ii]   = 0.0 ; 
  }
  
  ic->sample_and_assign_initial_positions_to_heavy_quarks(Nheavyquarks, x, y, etas);
  ic->sample_and_assign_initial_momenta_to_heavy_quarks(Nheavyquarks, px, py, pz, PID, status, x, y, etas);

  // derived quantities
  for(int ii=0; ii<Nheavyquarks; ii++){
    tau[ii] = tau0 ; 
    t[ii] = tau0*cosh(etas[ii]) ; 
    z[ii] = tau0*sinh(etas[ii]) ; 
    E[ii] = sqrt( px[ii] * px[ii] + py[ii] * py[ii] 
            + pz[ii] * pz[ii] + quark_mass * quark_mass ) ;
  }


  // write output to binary file //
  std::ofstream outbin1;
  std::stringstream output_filename11;
  output_filename11.str("");
  output_filename11 <<"init_quark_list.bin";
  outbin1.open(output_filename11.str().c_str(), std::ios::out | std::ios::binary);
  for(int ii=0; ii<Nheavyquarks; ii++){
    float particle_array1[] = { static_cast<float>(PID[ii]),
                                static_cast<float>(status[ii]),
                                static_cast<float>(t[ii]),
                                static_cast<float>(x[ii]),
                                static_cast<float>(y[ii]),
                                static_cast<float>(z[ii]),
                                static_cast<float>(E[ii]),
                                static_cast<float>(px[ii]),
                                static_cast<float>(py[ii]),
                                static_cast<float>(pz[ii])
                             };
    for(int jj=0; jj<10;jj++){
        outbin1.write( reinterpret_cast<char *>(&(particle_array1[jj])),
                      sizeof(float) );
    }
  }
  outbin1.close() ; 


//#pragma omp parallel
  {  
//#pragma omp for schedule(guided)
    for(int ii=0; ii<Nheavyquarks; ii++){ // loop over quarks
      //for(int jj=0; jj<Ntimestep; jj++){ // loop over time step
       if(ii%1000==0){
         std::cout << "Working on Quark No. " << ii << std::endl ; 
       }
      T = 1.;

      // Things added on 14th July 2023//
      T  = blk->get_temperature(tau[ii], x[ii], y[ii], etas[ii]);
      if( fabs(etas[ii]) > 7.95 ){
        T  = 0.01;
      }

      while(T>0.15){

	// get the local field values //
	T  = blk->get_temperature(tau[ii], x[ii], y[ii], etas[ii]);
	vx = blk->get_vx(tau[ii], x[ii], y[ii], etas[ii]);
	vy = blk->get_vy(tau[ii], x[ii], y[ii], etas[ii]);
	vz = blk->get_vz(tau[ii], x[ii], y[ii], etas[ii]);

	// vz regulation //
        // if( fabs(vz) > 0.999000 ){
        //   if(vz < 0 ){
        //     vz = -0.999000 ;
	//   }
	//   else{
        //     vz =  0.999000 ;
	//   }
	// }

        
        //std::cout << std::fixed << std::setprecision(3) << "tau=" << tau[ii] << " fm,   T= " << T << " GeV,  " 
        //         << "x= " << x[ii] << " fm,  y= " << y[ii] <<   " fm,  etas= " << etas[ii] << ",  " 
        //         << "px= " << px[ii] << " GeV,  py= " << py[ii] <<   " GeV,  pz= " << pz[ii] << " GeV ... "
        //         << std::endl ; 
        

	// Transormation from Lab frame to local rest frame of the fluid.
	double tstar=0., xstar=0., ystar=0., zstar=0., Estar=0., pxstar=0., pystar=0., pzstar=0. ; 
	lorentz_transformation_x( t[ii], x[ii], y[ii], z[ii], vx, vy, vz, tstar, xstar, ystar, zstar );
	lorentz_transformation_p( E[ii], px[ii], py[ii], pz[ii], vx, vy, vz, Estar, pxstar, pystar, pzstar );
	
	//evolve the heavy quark position and momentum in local rest frame.
	double dtstar = (Estar / E[ii]) * dt ;
	update_position( tstar, xstar, ystar, zstar, Estar, pxstar,
			 pystar, pzstar, dtstar, tstar, xstar, ystar, zstar );
	double rhox = rand->random_gaussian(0.,1.);
	double rhoy = rand->random_gaussian(0.,1.);
	double rhoz = rand->random_gaussian(0.,1.);
	double drag_coeff = get_drag_coefficient(pxstar, pystar, pzstar, T);

        //if(dtstar>(1./drag_coeff)){
        //  std::cout << "dt > Relaxation_time(=1/drag_coeff)" << std::endl ; 
        //  exit(1);
        //}
	
	double temp_diffusion_coeff = drag_coeff * Estar * T ; 
	double temp_pxstar, temp_pystar, temp_pzstar, temp_estar ; 
	update_momentum(pxstar, pystar, pzstar, drag_coeff, temp_diffusion_coeff , rhox,
                        rhoy, rhoz, dtstar, temp_pxstar, temp_pystar, temp_pzstar );
	temp_estar = sqrt(temp_pxstar*temp_pxstar + temp_pystar*temp_pystar + 
			  temp_pzstar*temp_pzstar + quark_mass*quark_mass);
	
	double temp_drag_coeff = get_drag_coefficient(temp_pxstar,temp_pystar,temp_pzstar,T) ;
	double diffusion_coeff  =  temp_drag_coeff * temp_estar * T ; 
	
	// update the momentum with post point realization //
	update_momentum(pxstar, pystar, pzstar, drag_coeff, diffusion_coeff,
			rhox, rhoy, rhoz, dtstar, pxstar, pystar, pzstar);
	Estar = sqrt( pxstar*pxstar + pystar*pystar + 
		      pzstar*pzstar + quark_mass*quark_mass);

	// Transormation from local rest frame to Lab frame of the fluid.
	lorentz_transformation_x( tstar, xstar, ystar, zstar, -vx, -vy, -vz, t[ii], x[ii], y[ii], z[ii] );
	lorentz_transformation_p( Estar, pxstar, pystar, pzstar, -vx, -vy, -vz, E[ii], px[ii], py[ii], pz[ii] );
	
        // findout tau and etas 
	tau[ii]  = sqrt(t[ii]*t[ii] - z[ii]*z[ii]);
	etas[ii] = atanh(z[ii]/t[ii]);
        if( std::isinf(tau[ii]) || std::isnan(tau[ii]) ){
          std::cout << "tau = inf/nan" << std::endl ; 
          exit(1);
        }
        if( std::isinf(etas[ii]) || std::isnan(etas[ii]) ){
          std::cout << "etas = inf/nan" << std::endl ; 
          exit(1);
        }	

        if(tau[ii] > (tau0+15) ){
          std::cout << "evolving for a very long time" 
                    << " but not freezeing out ..." 
                    << std::endl ; 
          std::cout << "Current T = " << T << " GeV." << std::endl ; 
          std::cout << "Current tau = " << tau[ii] << " fm, tau0 = " << tau0 << " fm" << std::endl ; 
          exit(1);
        }

      } // loop over timestep
    } // loop over quarks
  } // omp block
  

  for(int ii=0; ii<Nheavyquarks; ii++){
   E[ii] = sqrt( px[ii] * px[ii] + py[ii] * py[ii] 
            + pz[ii] * pz[ii] + quark_mass * quark_mass ) ;
  }


  // write output to binary file //
  std::ofstream outbin;
  std::stringstream output_filename1;
  output_filename1.str("");
  output_filename1 <<"heavy_quark_list.bin";
  outbin.open(output_filename1.str().c_str(), std::ios::out | std::ios::binary);
  for(int ii=0; ii<Nheavyquarks; ii++){
    float particle_array[] = { static_cast<float>(PID[ii]),
                               static_cast<float>(status[ii]),
                               static_cast<float>(t[ii]),
                               static_cast<float>(x[ii]),
                               static_cast<float>(y[ii]),
                               static_cast<float>(z[ii]),
                               static_cast<float>(E[ii]),
                               static_cast<float>(px[ii]),
                               static_cast<float>(py[ii]),
                               static_cast<float>(pz[ii])
                             };
    for(int jj=0; jj<10;jj++){
        outbin.write( reinterpret_cast<char *>(&(particle_array[jj])),
                      sizeof(float) );
    }
  }
  outbin.close() ; 

  // petersen fragmentation
  hadron->fragment(Nheavyquarks, px, py);
  

  for(int ii=0; ii<Nheavyquarks; ii++){
   E[ii] = sqrt( px[ii] * px[ii] + py[ii] * py[ii] 
            + pz[ii] * pz[ii] + quark_mass * quark_mass ) ;
  }


  output_filename1.str("");
  output_filename1 <<"heavy_hadron_list.bin";
  outbin.open(output_filename1.str().c_str(), std::ios::out | std::ios::binary);
  for(int ii=0; ii<Nheavyquarks; ii++){
    float particle_array[] = { static_cast<float>(PID[ii]),
                               static_cast<float>(status[ii]),
                               static_cast<float>(t[ii]),
                               static_cast<float>(x[ii]),
                               static_cast<float>(y[ii]),
                               static_cast<float>(z[ii]),
                               static_cast<float>(E[ii]),
                               static_cast<float>(px[ii]),
                               static_cast<float>(py[ii]),
                               static_cast<float>(pz[ii])
                             };
    for(int jj=0; jj<10;jj++){
        outbin.write( reinterpret_cast<char *>(&(particle_array[jj])),
                      sizeof(float) );
    }
  }
  outbin.close() ; 




  /*
  std::ofstream outdat;
  outdat.open("heavy_hadron_list.dat");
  for(int ii=0; ii<Nheavyquarks; ii++){
    outdat << ii << "  " << t[ii] << "  " << x[ii] <<
                 << "  " << y[ii] << "  " << z[ii] << 
                 << "  " << E[ii] << "  " << px[ii] << 
                 << "  " << py[ii] << "  " << pz[ii] <<
                 std::endl ;  
  }
  outdat.close();
  */
  

  delete tau ;
  delete etas;
  delete t ; 
  delete z ;
  delete x ;
  delete y ; 
  delete px ; 
  delete py ; 
  delete pz ; 
  delete E ; 
  delete PID ;
  delete status ;  
}










