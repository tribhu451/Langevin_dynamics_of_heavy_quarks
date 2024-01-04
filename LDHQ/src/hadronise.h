#pragma once
#include<iostream>
#include<fstream>
#include<sstream>
#include<string>
#include<cmath>
#include "random.h"

using std::istringstream ;

class hadronise{


  private :

    random_gen* rand ;
    double zprob(double z);
    double samplez();

  public :
    hadronise();
    ~hadronise();
    void fragment(int nq, double* pxx, double* pyy);



};

