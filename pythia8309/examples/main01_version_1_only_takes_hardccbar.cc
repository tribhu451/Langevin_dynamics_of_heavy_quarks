// main01.cc is a part of the PYTHIA event generator.
// Copyright (C) 2023 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Keywords: basic usage; charged multiplicity;

// This is a simple test program. It fits on one slide in a talk.
// It studies the charged multiplicity distribution at the LHC.

/*
#include "Pythia8/Pythia.h"
using namespace Pythia8;
int main() {
  // Generator. Process selection. LHC initialization. Histogram.
  Pythia pythia;
  pythia.readString("Beams:eCM = 8000.");
  pythia.readString("HardQCD:all = on");
  pythia.readString("PhaseSpace:pTHatMin = 20.");
  pythia.init();
  Hist mult("charged multiplicity", 100, -0.5, 799.5);
  // Begin event loop. Generate event. Skip if error. List first one.
  for (int iEvent = 0; iEvent < 100; ++iEvent) {
    if (!pythia.next()) continue;
    // Find number of all final charged particles and fill histogram.
    int nCharged = 0;
    for (int i = 0; i < pythia.event.size(); ++i)
      if (pythia.event[i].isFinal() && pythia.event[i].isCharged())
        ++nCharged;
    mult.fill( nCharged );
  // End of event loop. Statistics. Histogram. Done.
  }
  pythia.stat();
  cout << mult;
  return 0;

}
*/





#include "Pythia8/Pythia.h"
#include <vector>

using namespace Pythia8;
int main() {
  // Generator. Process selection. LHC initialization. Histogram.
  Pythia pythia;
  pythia.readString("Beams:idA = 2212");
  pythia.readString("Beams:idB = 2212");
  pythia.readString("Beams:eCM = 2760");
  //pythia.readString("HardQCD:all = on");
  //pythia.readString("SoftQCD:all = on");
  //pythia.readString("SoftQCD:inelastic = on");
  
  
  //pythia.readString("HardQCD:gg2ccbar = on");
  //pythia.readString("HardQCD:qqbar2ccbar = on");
  pythia.readString("HardQCD:hardccbar = on");
  
  
  
  // Pick new random number seed for each run, based on clock.
  pythia.readString("Random:setSeed = on");
  pythia.readString("Random:seed = 0");
  pythia.init();
  
  int Nevents = 10000 ; // No of PYTHIA events to be generated.

  std::ofstream file;
  file.open("init_ccbar_info.dat");
  file << Nevents << std::endl ;
  file << "#EventID  MCID  status SrNoPrimaryC SrNoPrimaryMother1 SrNoPrimaryMother2 Prim_e Prim_px Prim_py Prim_Pz SrNoFinalProduct " << 
     "SrNoFinalProductsMother1 SrNoFinalProductMother2  e  px  py  pz" << std::endl ; 

  std::ofstream ofile;
  ofile.open("delta_phi_and_delta_eta.dat");
  
  int produced_charm_quark_number = 0 ; 
  int produced_charm_quark_pid[200] ;
  int produced_charm_quark_eventID[200] ;
  int produced_charm_quark_status[200] ;
  int produced_charm_quark_index[200] ;
  double produced_charm_quark_mass[200] ;
  int c_mother1[200];  
  int c_mother2[200];  
  double c_e[200];  
  double c_px[200];  
  double c_py[200];  
  double c_pz[200];  
 
  int counter_primary_cs = 0 ; 
  int index_of_primary_cs[20] ;
  int mother1_of_primary_cs[20] ;
  int mother2_of_primary_cs[20] ; 
  int index_of_final_product_cs[20] ; 
  int pairing_flag_of_cs[20] ; 
  
  // Begin event loop. Generate event. Skip if error. List first one.
  for (int iEvent = 0; iEvent < Nevents; ++iEvent) {
    
    counter_primary_cs = 0 ; 
    //file << "  " << std::endl ;
    
    if (!pythia.next()) continue;
    produced_charm_quark_number = 0 ;  
    for (int i = 0; i < pythia.event.size(); ++i){
      
      int mother1 =  pythia.event[i].mother1() ; 
      int mother2 = pythia.event[i].mother2() ;
      int mother1_pid = pythia.event[mother1].id() ;  
      int mother2_pid = pythia.event[mother2].id() ;  
      int pid = pythia.event[i].id() ;   
      int status = pythia.event[i].status() ;
      double e = pythia.event[i].e() ;    
      double px = pythia.event[i].px() ;    
      double py = pythia.event[i].py() ;    
      double pz = pythia.event[i].pz() ;    
      //double mass = sqrt(e*e - px * px - py * py - pz * pz) ;    
      double mass = pythia.event[i].m() ;    
      
      // if charm quarks
      if(abs(pid) == 4){ 
	int inclus_stat = 0 ;  
	
	// produced only through the outgoing process
	if(fabs(status) == 23){
	  inclus_stat += 1 ; 
	}
	
	// produced from gluon shower
	if( mother2 == 0 && mother1_pid == 21  ){
          inclus_stat += 1 ;
	}
	
	// storing the primordially produced charm quarks
	if(inclus_stat > 0 ){
	  index_of_primary_cs[counter_primary_cs] = i ; 
	  mother1_of_primary_cs[counter_primary_cs]=mother1 ;
	  mother2_of_primary_cs[counter_primary_cs]=mother2 ;
	  index_of_final_product_cs[counter_primary_cs] = i ;
	  counter_primary_cs += 1 ; 
	}
        
	produced_charm_quark_pid[produced_charm_quark_number] = pid ;
	produced_charm_quark_mass[produced_charm_quark_number] = mass ;
	produced_charm_quark_index[produced_charm_quark_number] = i ;
	produced_charm_quark_eventID[produced_charm_quark_number] = iEvent ;
	produced_charm_quark_status[produced_charm_quark_number] = status ;
        c_mother1[produced_charm_quark_number] = mother1 ; 	
        c_mother2[produced_charm_quark_number] = mother2 ; 	
	c_e[produced_charm_quark_number] = e ; 
	c_px[produced_charm_quark_number] = px ; 
	c_py[produced_charm_quark_number] = py ; 
	c_pz[produced_charm_quark_number] = pz ; 
	produced_charm_quark_number += 1 ;
	
	// write down the entire history of all the produced charm quarks  
	// file << iEvent << "  " << i << "  "
	//      << pid << "  " << status << "  " 
	//      << mother1 << "  " << mother2 << "  "
	//      << e << "  " <<  px << "  " 
	//      << py << "  " << pz << "  " 
	//      << sqrt(fabs( e*e - px*px - py*py - pz*pz)) 
	//      << std::endl ;
        
      } // if charm quark
    } // particle loop
    
    
    
    
    
    
    // finout the final product
    file << counter_primary_cs << std::endl ; 
    for(int ii=0; ii<counter_primary_cs; ii++){
      pairing_flag_of_cs[ii] = 0 ; 
      int counter = 0 ; 
      while( counter < produced_charm_quark_number ){
        if( c_mother1[counter] == index_of_final_product_cs[ii] || 
            c_mother2[counter] == index_of_final_product_cs[ii] ){
	  index_of_final_product_cs[ii] = produced_charm_quark_index[counter] ; 
	  counter = 0 ; 
	}
	counter++ ; 
      } // while loop
    }
    
    
    
    
    
    // findout the pairs and write them in file
    for(int ii=0; ii<counter_primary_cs; ii++){
      if(pairing_flag_of_cs[ii]>0){
	continue ; 
      }
      for(int jj=0; jj<counter_primary_cs; jj++){
	if(ii==jj){
	  continue ; 
	}
	if( mother1_of_primary_cs[ii] == mother1_of_primary_cs[jj] &&
	    mother2_of_primary_cs[ii] == mother2_of_primary_cs[jj] ){ // pair found
	  
	  // first find out the actual array position of the final charm pairs in current event
	  int f1, f2 ;
	  f1 = 100000 ; 
	  f2 = 100000 ;  
	  for(int ll=0; ll<produced_charm_quark_number; ll++){
	    if( produced_charm_quark_index[ll] == index_of_final_product_cs[ii] ){
	      f1 = ll ; 
	    }
	  }
	  
	  for(int ll=0; ll<produced_charm_quark_number; ll++){
	    if( produced_charm_quark_index[ll] == index_of_final_product_cs[jj] ){
	      f2 = ll ; 
	    }
	  }
	  
	  if( f1 > 9999 || f2 > 9999){
	    std::cout << "unable to get actual index of charm quark ..." << std::endl ; 
	    exit(1);
	  }
	  
	  int i1, i2 ; 
	  i1 = 100000 ; 
	  i2 = 100000 ; 
	  
	  for(int ll=0; ll<produced_charm_quark_number; ll++){
	    if( produced_charm_quark_index[ll] == index_of_primary_cs[ii] ){
	      i1 = ll ; 
	    }
	  }
	  
	  for(int ll=0; ll<produced_charm_quark_number; ll++){
	    if( produced_charm_quark_index[ll] == index_of_primary_cs[jj] ){
	      i2 = ll ; 
	    }
	  }
	  
	  if(produced_charm_quark_pid[f1] == 4){   
	    file << produced_charm_quark_eventID[f1] << "  " 
		 << produced_charm_quark_pid[f1]<< "  "
		 << produced_charm_quark_status[i1] << "  "
		 << index_of_primary_cs[ii] << "  " 
		 << c_mother1[i1] << "  " 
		 << c_mother2[i1] << "  " 
		 << c_e[i1] << "  " 
		 << c_px[i1] << "  " << c_py[i1] << "  " << c_pz[i1] << "  " 
		 << index_of_final_product_cs[ii] << "  "
		 << c_mother1[f1] << "  "	
		 << c_mother2[f1] << "  "	
		 << c_e[f1] << "  "
		 << c_px[f1] << "  " << c_py[f1] << "  " << c_pz[f1] << std::endl ;
	  
	    file << produced_charm_quark_eventID[f2] << "  " 
		 << produced_charm_quark_pid[f2]<< "  "
		 << produced_charm_quark_status[i2] << "  "
		 << index_of_primary_cs[jj] << "  " 
		 << c_mother1[i2] << "  " 
		 << c_mother2[i2] << "  " 
		 << c_e[i2] << "  " 
		 << c_px[i2] << "  " << c_py[i2] << "  " << c_pz[i2] << "  " 
		 << index_of_final_product_cs[jj] << "  "
		 << c_mother1[f2] << "  "	
		 << c_mother2[f2] << "  "	
		 << c_e[f2] << "  "
		 << c_px[f2] << "  " << c_py[f2] << "  " << c_pz[f2] << std::endl ;
	  
	  }
	  else{
	    
	   
	    file << produced_charm_quark_eventID[f2] << "  " 
		 << produced_charm_quark_pid[f2]<< "  "
		 << produced_charm_quark_status[i2] << "  "
		 << index_of_primary_cs[jj] << "  " 
		 << c_mother1[i2] << "  " 
		 << c_mother2[i2] << "  " 
		 << c_e[i2] << "  " 
		 << c_px[i2] << "  " << c_py[i2] << "  " << c_pz[i2] << "  " 
		 << index_of_final_product_cs[jj] << "  "
		 << c_mother1[f2] << "  "	
		 << c_mother2[f2] << "  "	
		 << c_e[f2] << "  "
		 << c_px[f2] << "  " << c_py[f2] << "  " << c_pz[f2] << std::endl ;


	    file << produced_charm_quark_eventID[f1] << "  " 
		 << produced_charm_quark_pid[f1]<< "  "
		 << produced_charm_quark_status[i1] << "  "
		 << index_of_primary_cs[ii] << "  " 
		 << c_mother1[i1] << "  " 
		 << c_mother2[i1] << "  " 
		 << c_e[i1] << "  " 
		 << c_px[i1] << "  " << c_py[i1] << "  " << c_pz[i1] << "  " 
		 << index_of_final_product_cs[ii] << "  "
		 << c_mother1[f1] << "  "	
		 << c_mother2[f1] << "  "	
		 << c_e[f1] << "  "
		 << c_px[f1] << "  " << c_py[f1] << "  " << c_pz[f1] << std::endl ;
	    
	  }
	  pairing_flag_of_cs[ii] = 1 ; 
	  pairing_flag_of_cs[jj] = 1 ;
	  
	} // if pair found
      } // jj loop
    } // ii loop
    
    
    
    
    
    
  } // End of event loop.
  
  file.close() ; 
  ofile.close() ; 
  pythia.stat();
  return 0;
  
}
















