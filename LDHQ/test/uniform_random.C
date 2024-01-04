
//double rand_uniform(){
//  return(rand_uniform_dist(*ran_generator));
//}


//int get_seed() const {
//  return(seed);
//}





void uniform_random(){

  TH1D* h1 = new TH1D("h1", "h1", 10, 0., 1. );

  long seed;
  std::random_device ran_dev;
  std::unique_ptr<std::mt19937> ran_generator;
  std::uniform_real_distribution<double> rand_uniform_dist(0,1);

  seed = ran_dev();
  //seed = 1234 ; 
  ran_generator = std::unique_ptr<std::mt19937>(new std::mt19937(seed));

  for(int ii=0; ii<100; ii++){
    //std::cout << rand_uniform_dist(*ran_generator) << std::endl ;
   h1->Fill(rand_uniform_dist(*ran_generator));
  } 

  h1->Draw() ; 

}
