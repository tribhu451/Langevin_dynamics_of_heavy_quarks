void read_and_plot_random_dist(){

  double num;

  TH1D* h1 = new TH1D("h1", "h1", 100, -10, 10 );

  ifstream mFile;
  mFile.open("../random_number_distribution.dat");
  if(!mFile){
    cout << "Input file doesn't exit ... " << endl;
    exit(1); 
  }

  int icount=0;
  while(!mFile.eof()){
    mFile >> num ;
    h1->Fill(num);
    icount++ ; 
  }

  h1->Draw() ; 

}
