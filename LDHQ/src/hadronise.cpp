#include "hadronise.h"

hadronise::hadronise(){
  rand = new random_gen();
}

hadronise::~hadronise(){
}

double hadronise::zprob(double z){
  return 1. / (6.3 * z * (1 - 1 / z - 0.04 / (1 - z)) * (1 - 1 / z - 0.04 / (1 - z)));
}

double hadronise::samplez(){
  double z, frac, prob ;
  int d=0;
  int counter = 0 ; 
  while(d==0){
    z = rand->rand_uniform();
    frac = zprob(z);
    prob = rand->rand_uniform();
    if(prob < frac){
      d=1;
    }
    else{ 
     counter += 1 ;
     if(counter > 10000){
       std::cout << "so many trials but could not "  
                 << "sample the z fraction from petersen fragmentation ... " 
                 << std::endl ; 
       exit(1);
     }
     continue ; 
    }
  };
  return z;
}



void hadronise::fragment(int nq, double* pxx, double* pyy){
  double z ; 
  for(int ii=0; ii<nq; ii++){
    z = samplez();
    pxx[ii] = z * pxx[ii];
    pyy[ii] = z * pyy[ii];
  }

}




