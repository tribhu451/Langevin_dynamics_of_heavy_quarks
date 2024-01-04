#include "bulk.h"


bulk::bulk(input_paramters &iparam_):iparam(iparam_),ntau(2),nx(2),ny(2),neta(2){
  bulk_file_name = iparam.bulk_file.c_str() ; 
  std::string dummy ;
  double  Ataumin;
  double  Adtau;
  int     Antau;
  double  Axmin;
  double  Adx;
  int     Anx; 
  double  Aymin;
  double  Ady;
  int     Any;
  double  Aetamin;
  double  Adeta;
  int     Aneta;
 
  std::ifstream ifile;
  ifile.open(bulk_file_name.c_str(), std::ios::in );
  if(!ifile){
    std::cout << "bulk information input" <<
    " file not found ... " << std::endl ; 
    exit(1);
  }
  else{
   std::cout << "bulk file : " << bulk_file_name.c_str() << std::endl ; 
  }
  ifile.getline(buff,350);
  iss = new istringstream(buff);
  *iss >> dummy >> Ataumin >> dummy 
       >> Adtau >> dummy >> Antau 
       >> dummy >> Axmin   >> dummy 
       >> Adx   >> dummy >> Anx 
       >> dummy >> Aymin   >> dummy 
       >> Ady   >> dummy >> Any 
       >> dummy >> Aetamin >> dummy 
       >> Adeta >> dummy >> Aneta ;
  delete iss;
  ifile.close();

  ntau   =  Antau   ;
  nx     =  Anx     ; 
  ny     =  Any     ;
  neta   =  Aneta   ;
  taumin =  Ataumin ;
  xmin   =  Axmin   ; 
  ymin   =  Aymin   ;
  etamin =  Aetamin ;
  dtau   =  Adtau   ;
  dx     =  Adx     ;
  dy     =  Ady     ;
  deta   =  Adeta   ;
 
  std::cout << "==============" << std::endl ; 
  std::cout << "taumin = " << taumin << std::endl ;  
  std::cout << "xmin   = " << xmin << std::endl ;  
  std::cout << "ymin   = " << ymin << std::endl ;  
  std::cout << "etamin = " << etamin << std::endl ;  
  std::cout << "ntau = " << ntau << std::endl ; 
  std::cout << "nx   = " << nx << std::endl ; 
  std::cout << "ny   = " << ny << std::endl ; 
  std::cout << "neta = " << neta << std::endl ; 
  std::cout << "dtau = " << dtau << std::endl ; 
  std::cout << "dx   = " << dx << std::endl ; 
  std::cout << "dy   = " << dy << std::endl ; 
  std::cout << "deta = " << deta << std::endl ; 
  std::cout << "==============" << std::endl ; 

  // initialize the hyper-grid to store bulk information.
  // temperature
  temperature = new double*** [ntau];
  for (int itau = 0; itau < ntau; itau++) {
    temperature[itau] = new double** [nx];
    for (int ix = 0; ix < nx; ix++) {
      temperature[itau][ix] = new double* [ny];
      for (int iy = 0; iy < ny; iy++) {
        temperature[itau][ix][iy] = new double[neta];
        for (int ieta = 0; ieta < neta; ieta++){
          temperature[itau][ix][iy][ieta] = 0.0;
        }  
      } 
    } 
  }  


 //  velocity-x 
  vx = new double*** [ntau];
  for (int itau = 0; itau < ntau; itau++) {
    vx[itau] = new double** [nx];
    for (int ix = 0; ix < nx; ix++) {
      vx[itau][ix] = new double* [ny];
      for (int iy = 0; iy < ny; iy++) {
        vx[itau][ix][iy] = new double[neta];
        for (int ieta = 0; ieta < neta; ieta++){
          vx[itau][ix][iy][ieta] = 0.0;
        }  
      } 
    } 
  } 


 //  velocity-y 
  vy = new double*** [ntau];
  for (int itau = 0; itau < ntau; itau++) {
    vy[itau] = new double** [nx];
    for (int ix = 0; ix < nx; ix++) {
      vy[itau][ix] = new double* [ny];
      for (int iy = 0; iy < ny; iy++) {
        vy[itau][ix][iy] = new double[neta];
        for (int ieta = 0; ieta < neta; ieta++){
          vy[itau][ix][iy][ieta] = 0.0;
        }  
      } 
    } 
  } 


 //  velocity-z
  vz = new double*** [ntau];
  for (int itau = 0; itau < ntau; itau++) {
    vz[itau] = new double** [nx];
    for (int ix = 0; ix < nx; ix++) {
      vz[itau][ix] = new double* [ny];
      for (int iy = 0; iy < ny; iy++) {
        vz[itau][ix][iy] = new double[neta];
        for (int ieta = 0; ieta < neta; ieta++){
          vz[itau][ix][iy][ieta] = 0.0;
        }  // ieta loop
      } // iy loop
    } // ix loop
  } // itau loop


}



bulk::~bulk(){

  // clean up
  for (int i = 0; i < ntau; i++){
    for (int j = 0; j < nx; j++){
       for (int k = 0; k < ny; k++){
         delete [] temperature[i][j][k];
         delete [] vx[i][j][k];
         delete [] vy[i][j][k];
         delete [] vz[i][j][k];
       }
      delete [] temperature[i][j];
      delete [] vx[i][j];
      delete [] vy[i][j];
      delete [] vz[i][j];
    }
    delete [] temperature[i];
    delete [] vx[i];
    delete [] vy[i];
    delete [] vz[i];
  }
  delete [] temperature;
  delete [] vx;
  delete [] vy;
  delete [] vz;

}


double bulk::get_hydro_bulk_tau0(){
  return taumin ; 
}


void bulk::read_bulk_info(){

  double tau_temp, x_temp, y_temp, eta_temp ; 
  double epsilon_temp, rhob_temp ;
  double temperature_temp, mub_temp ;
  double vx_temp, vy_temp, vz_temp ; 

  std::ifstream ifile;
  ifile.open(bulk_file_name.c_str(), std::ios::in );

  ifile.getline(buff,300) ;
  std::cout << "reading bulk information from input file ... " << std::endl ; 
  for(int itau=0; itau<ntau; itau++){
    //if( itau%10 == 0){
     //std::cout << "🌏 itau(ntau)= " << itau << "(" << ntau << ")" << std::endl ;
    //} 
    for(int ieta=0; ieta<neta; ieta++){
      for(int iy=0; iy<ny; iy++){
        for(int ix=0; ix<nx; ix++){

          ifile.getline(buff,200);
          iss = new istringstream(buff);

          *iss >> tau_temp >> x_temp >> y_temp >> eta_temp 
                >> epsilon_temp >> rhob_temp >> temperature_temp 
                >> mub_temp >> vx_temp >> vy_temp >> vz_temp ;

           delete iss ; 

           //if( itau%10 == 0 && ix == 0 && iy == 0 && ieta == 0 ){
              //std::cout << "🌏tau = " << tau_temp << " fm." << std::endl ; 
           //}

           if( itau%10 == 0 && ix == nx/2 && iy == ny/2 && ieta == neta/2 ){
             std::cout << "🌏 tau= " << tau_temp << " fm, x= " << x_temp
                           << " fm, y= " << y_temp << " fm,  etas= " << eta_temp
                             << ", T= " << temperature_temp << " GeV  ... " << std::endl ; 
           }
 
           if(temperature_temp < 0){
             std::cout << "negetive temperature in the bulk." << std::endl ; 
             std::cout << "May be file format is not compatible !!!" << std::endl ; 
             exit(1);
           }


                // velocities are in cartesian coordinate 
          temperature  [itau][ix][iy][ieta]  =  temperature_temp ; 
          vx           [itau][ix][iy][ieta]  =           vx_temp ; 
          vy           [itau][ix][iy][ieta]  =           vy_temp ; 
          vz           [itau][ix][iy][ieta]  =           vz_temp ; 


  
        }  // ix loop
      }  // iy loop
    }  // ieta loop
  }  // itau loop

  ifile.close(); 
  std::cout << "reading bulk information completed ... 👍 " << std::endl ; 
}


double bulk::get_temperature(double mtau, double mx, double my, double meta){

  double tau_0000 = taumin + floor_itau(mtau) * dtau ;  
  double x_0000   = xmin   + floor_ix(mx)     * dx   ;  
  double y_0000   = ymin   + floor_iy(my)     * dy   ;  
  double eta_0000 = etamin + floor_ieta(meta) * deta ;  

  //std::cout << "tau_floor = " << tau_0000 << ",  x_floor = " 
      //<< x_0000 << ",  y_floor = " << y_0000 << ",  eta_floor = "
         //<< eta_0000 << std::endl ; 

  double lattice_spacing[4] = {dtau, dx, dy, deta};
  double x_fraction[2][4];
  x_fraction[1][0] = mtau - tau_0000 ; 
  x_fraction[1][1] = mx   - x_0000   ; 
  x_fraction[1][2] = my   - y_0000   ; 
  x_fraction[1][3] = meta - eta_0000 ;

  for (int ii = 0; ii < 4; ii++) {
    x_fraction[0][ii] = lattice_spacing[ii] - x_fraction[1][ii];
  }

   // initialize the hyper-cube for Cornelius
   double ****cube = new double*** [2];
    for (int i = 0; i < 2; i++) {
      cube[i] = new double** [2];
      for (int j = 0; j < 2; j++) {
        cube[i][j] = new double* [2];
        for (int k = 0; k < 2; k++) {
           cube[i][j][k] = new double[2];
           for (int l = 0; l < 2; l++)
             cube[i][j][k][l] = 0.0;
       }
     }
   }

  int jtau = floor_itau(mtau) ;
  int jx   = floor_ix(mx) ; 
  int jy   = floor_iy(my) ; 
  int jeta = floor_ieta(meta) ;

  //std::cout << "floor_tau_index = " << jtau << ",  floor_x_index = " 
      //<< jx << ",  floor_y_index = " << jy << ",  floor_eta_index = "
         //<< jeta << std::endl ;   

  cube[0][0][0][0] = temperature[jtau][jx][jy][jeta];
  cube[0][0][1][0] = temperature[jtau][jx][jy+1][jeta];
  cube[0][1][0][0] = temperature[jtau][jx+1][jy][jeta];
  cube[0][1][1][0] = temperature[jtau][jx+1][jy+1][jeta];
  cube[1][0][0][0] = temperature[jtau+1][jx][jy][jeta];
  cube[1][0][1][0] = temperature[jtau+1][jx][jy+1][jeta];
  cube[1][1][0][0] = temperature[jtau+1][jx+1][jy][jeta];
  cube[1][1][1][0] = temperature[jtau+1][jx+1][jy+1][jeta];
  cube[0][0][0][1] = temperature[jtau][jx][jy][jeta+1];
  cube[0][0][1][1] = temperature[jtau][jx][jy+1][jeta+1];
  cube[0][1][0][1] = temperature[jtau][jx+1][jy][jeta+1];
  cube[0][1][1][1] = temperature[jtau][jx+1][jy+1][jeta+1];
  cube[1][0][0][1] = temperature[jtau+1][jx][jy][jeta+1];
  cube[1][0][1][1] = temperature[jtau+1][jx][jy+1][jeta+1];
  cube[1][1][0][1] = temperature[jtau+1][jx+1][jy][jeta+1];
  cube[1][1][1][1] = temperature[jtau+1][jx+1][jy+1][jeta+1];

  double interpolated_value =  
      four_dimension_linear_interpolation(lattice_spacing, x_fraction, cube);

    // delete cube 
    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 2; j++) {
            for (int k = 0; k < 2; k++)
                delete [] cube[i][j][k];
            delete [] cube[i][j];
        }
        delete [] cube[i];
    }
    delete [] cube;

  return interpolated_value ; 
}



double bulk::get_vx(double mtau, double mx, double my, double meta){

  double tau_0000 = taumin + floor_itau(mtau) * dtau ;  
  double x_0000   = xmin   + floor_ix(mx)     * dx   ;  
  double y_0000   = ymin   + floor_iy(my)     * dy   ;  
  double eta_0000 = etamin + floor_ieta(meta) * deta ;  

  double lattice_spacing[4] = {dtau, dx, dy, deta};
  double x_fraction[2][4];
  x_fraction[1][0] = mtau - tau_0000 ; 
  x_fraction[1][1] = mx   - x_0000   ; 
  x_fraction[1][2] = my   - y_0000   ; 
  x_fraction[1][3] = meta - eta_0000 ;

  for (int ii = 0; ii < 4; ii++) {
    x_fraction[0][ii] = lattice_spacing[ii] - x_fraction[1][ii];
  }

   // initialize the hyper-cube for Cornelius
   double ****cube = new double*** [2];
    for (int i = 0; i < 2; i++) {
      cube[i] = new double** [2];
      for (int j = 0; j < 2; j++) {
        cube[i][j] = new double* [2];
        for (int k = 0; k < 2; k++) {
           cube[i][j][k] = new double[2];
           for (int l = 0; l < 2; l++)
             cube[i][j][k][l] = 0.0;
       }
     }
   }

  int jtau = floor_itau(mtau) ;
  int jx   = floor_ix(mx) ; 
  int jy   = floor_iy(my) ; 
  int jeta = floor_ieta(meta) ;

  cube[0][0][0][0] = vx[jtau][jx][jy][jeta];
  cube[0][0][1][0] = vx[jtau][jx][jy+1][jeta];
  cube[0][1][0][0] = vx[jtau][jx+1][jy][jeta];
  cube[0][1][1][0] = vx[jtau][jx+1][jy+1][jeta];
  cube[1][0][0][0] = vx[jtau+1][jx][jy][jeta];
  cube[1][0][1][0] = vx[jtau+1][jx][jy+1][jeta];
  cube[1][1][0][0] = vx[jtau+1][jx+1][jy][jeta];
  cube[1][1][1][0] = vx[jtau+1][jx+1][jy+1][jeta];
  cube[0][0][0][1] = vx[jtau][jx][jy][jeta+1];
  cube[0][0][1][1] = vx[jtau][jx][jy+1][jeta+1];
  cube[0][1][0][1] = vx[jtau][jx+1][jy][jeta+1];
  cube[0][1][1][1] = vx[jtau][jx+1][jy+1][jeta+1];
  cube[1][0][0][1] = vx[jtau+1][jx][jy][jeta+1];
  cube[1][0][1][1] = vx[jtau+1][jx][jy+1][jeta+1];
  cube[1][1][0][1] = vx[jtau+1][jx+1][jy][jeta+1];
  cube[1][1][1][1] = vx[jtau+1][jx+1][jy+1][jeta+1];

  double interpolated_value =  
      four_dimension_linear_interpolation(lattice_spacing, x_fraction, cube);

    // delete cube 
    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 2; j++) {
            for (int k = 0; k < 2; k++)
                delete [] cube[i][j][k];
            delete [] cube[i][j];
        }
        delete [] cube[i];
    }
    delete [] cube;

  return interpolated_value ; 
}



double bulk::get_vy(double mtau, double mx, double my, double meta){

  double tau_0000 = taumin + floor_itau(mtau) * dtau ;  
  double x_0000   = xmin   + floor_ix(mx)     * dx   ;  
  double y_0000   = ymin   + floor_iy(my)     * dy   ;  
  double eta_0000 = etamin + floor_ieta(meta) * deta ;  

  double lattice_spacing[4] = {dtau, dx, dy, deta};
  double x_fraction[2][4];
  x_fraction[1][0] = mtau - tau_0000 ; 
  x_fraction[1][1] = mx   - x_0000   ; 
  x_fraction[1][2] = my   - y_0000   ; 
  x_fraction[1][3] = meta - eta_0000 ;

  for (int ii = 0; ii < 4; ii++) {
    x_fraction[0][ii] = lattice_spacing[ii] - x_fraction[1][ii];
  }

   // initialize the hyper-cube for Cornelius
   double ****cube = new double*** [2];
    for (int i = 0; i < 2; i++) {
      cube[i] = new double** [2];
      for (int j = 0; j < 2; j++) {
        cube[i][j] = new double* [2];
        for (int k = 0; k < 2; k++) {
           cube[i][j][k] = new double[2];
           for (int l = 0; l < 2; l++)
             cube[i][j][k][l] = 0.0;
       }
     }
   }

  int jtau = floor_itau(mtau) ;
  int jx   = floor_ix(mx) ; 
  int jy   = floor_iy(my) ; 
  int jeta = floor_ieta(meta) ;

  cube[0][0][0][0] = vy[jtau][jx][jy][jeta];
  cube[0][0][1][0] = vy[jtau][jx][jy+1][jeta];
  cube[0][1][0][0] = vy[jtau][jx+1][jy][jeta];
  cube[0][1][1][0] = vy[jtau][jx+1][jy+1][jeta];
  cube[1][0][0][0] = vy[jtau+1][jx][jy][jeta];
  cube[1][0][1][0] = vy[jtau+1][jx][jy+1][jeta];
  cube[1][1][0][0] = vy[jtau+1][jx+1][jy][jeta];
  cube[1][1][1][0] = vy[jtau+1][jx+1][jy+1][jeta];
  cube[0][0][0][1] = vy[jtau][jx][jy][jeta+1];
  cube[0][0][1][1] = vy[jtau][jx][jy+1][jeta+1];
  cube[0][1][0][1] = vy[jtau][jx+1][jy][jeta+1];
  cube[0][1][1][1] = vy[jtau][jx+1][jy+1][jeta+1];
  cube[1][0][0][1] = vy[jtau+1][jx][jy][jeta+1];
  cube[1][0][1][1] = vy[jtau+1][jx][jy+1][jeta+1];
  cube[1][1][0][1] = vy[jtau+1][jx+1][jy][jeta+1];
  cube[1][1][1][1] = vy[jtau+1][jx+1][jy+1][jeta+1];

  double interpolated_value =  
      four_dimension_linear_interpolation(lattice_spacing, x_fraction, cube);

    // delete cube 
    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 2; j++) {
            for (int k = 0; k < 2; k++)
                delete [] cube[i][j][k];
            delete [] cube[i][j];
        }
        delete [] cube[i];
    }
    delete [] cube;

  return interpolated_value ; 
}



double bulk::get_vz(double mtau, double mx, double my, double meta){

  double tau_0000 = taumin + floor_itau(mtau) * dtau ; // nearest lower grid point in tau   
  double x_0000   = xmin   + floor_ix(mx)     * dx   ;  
  double y_0000   = ymin   + floor_iy(my)     * dy   ;  
  double eta_0000 = etamin + floor_ieta(meta) * deta ;  

  double lattice_spacing[4] = {dtau, dx, dy, deta};
  double x_fraction[2][4];
  x_fraction[1][0] = mtau - tau_0000 ; 
  x_fraction[1][1] = mx   - x_0000   ; 
  x_fraction[1][2] = my   - y_0000   ; 
  x_fraction[1][3] = meta - eta_0000 ;

  for (int ii = 0; ii < 4; ii++) {
    x_fraction[0][ii] = lattice_spacing[ii] - x_fraction[1][ii];
  }

   // initialize the hyper-cube for Cornelius
   double ****cube = new double*** [2];
    for (int i = 0; i < 2; i++) {
      cube[i] = new double** [2];
      for (int j = 0; j < 2; j++) {
        cube[i][j] = new double* [2];
        for (int k = 0; k < 2; k++) {
           cube[i][j][k] = new double[2];
           for (int l = 0; l < 2; l++)
             cube[i][j][k][l] = 0.0;
       }
     }
   }

  int jtau = floor_itau(mtau) ;
  int jx   = floor_ix(mx) ; 
  int jy   = floor_iy(my) ; 
  int jeta = floor_ieta(meta) ;

  cube[0][0][0][0] = vz[jtau][jx][jy][jeta];
  cube[0][0][1][0] = vz[jtau][jx][jy+1][jeta];
  cube[0][1][0][0] = vz[jtau][jx+1][jy][jeta];
  cube[0][1][1][0] = vz[jtau][jx+1][jy+1][jeta];
  cube[1][0][0][0] = vz[jtau+1][jx][jy][jeta];
  cube[1][0][1][0] = vz[jtau+1][jx][jy+1][jeta];
  cube[1][1][0][0] = vz[jtau+1][jx+1][jy][jeta];
  cube[1][1][1][0] = vz[jtau+1][jx+1][jy+1][jeta];
  cube[0][0][0][1] = vz[jtau][jx][jy][jeta+1];
  cube[0][0][1][1] = vz[jtau][jx][jy+1][jeta+1];
  cube[0][1][0][1] = vz[jtau][jx+1][jy][jeta+1];
  cube[0][1][1][1] = vz[jtau][jx+1][jy+1][jeta+1];
  cube[1][0][0][1] = vz[jtau+1][jx][jy][jeta+1];
  cube[1][0][1][1] = vz[jtau+1][jx][jy+1][jeta+1];
  cube[1][1][0][1] = vz[jtau+1][jx+1][jy][jeta+1];
  cube[1][1][1][1] = vz[jtau+1][jx+1][jy+1][jeta+1];

  double interpolated_value =  
      four_dimension_linear_interpolation(lattice_spacing, x_fraction, cube);

    // delete cube 
    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 2; j++) {
            for (int k = 0; k < 2; k++)
                delete [] cube[i][j][k];
            delete [] cube[i][j];
        }
        delete [] cube[i];
    }
    delete [] cube;

  return interpolated_value ; 
}



// four dimensional interpolation routine
double bulk::four_dimension_linear_interpolation(
        double* lattice_spacing, double fraction[2][4], double**** cube) {
   double denorm = 1.0;
   double results = 0.0;
   for (int i = 0; i < 4; i++) {
       denorm *= lattice_spacing[i];
   }
   for (int i = 0; i < 2; i++) {
       for (int j = 0; j < 2; j++) {
           for (int k = 0; k < 2; k++) {
               for (int l = 0; l < 2; l++) {
                   results += (cube[i][j][k][l]*fraction[i][0]*fraction[j][1]
                               *fraction[k][2]*fraction[l][3]);
               }
           }
       }
   }
   results = results/denorm;
   return (results);
}

















