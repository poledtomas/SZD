void draw(TH1D* hist1,TH1D* hist2,TH1D* hist3, int k){

  TCanvas* histogramy = new TCanvas("histogramy","histogramy",1800,800);

  if (k ==1) {
    histogramy->Divide(3,1);

    histogramy->cd(1);
    hist1->Draw();

    histogramy->cd(2);
    hist2->Draw();

    histogramy->cd(3);
    hist3->Draw();
    }

   else {
    histogramy->Divide(1,1);
    hist3->Draw();
  }

}
void contour(TH2D* hist) {
  TCanvas* contour = new TCanvas("contour","contour",800,800);
  gStyle->SetPalette(kRainBow);
  hist->SetStats(0);
  hist->SetContour(100);
  hist->Draw("COLZ");

}
float Maximum(TH2D* hist_contour, int k){

  double bin_max = hist_contour->GetMaximumBin();
  int A_max,B_max,z_max;
  hist_contour->GetBinXYZ(bin_max,A_max,B_max,z_max);

  if (k==0) {
    return A_max;
  } else {
    return B_max;
  }

}
void MLM() {

  int Nbins_sbg = 20000;
  int nBins = 100;
  double mean = 4.0, sigma = 2.0, xmin = 0, xmax = 15.0;
  int count = 0;

  double Amin = 0;
  double Bmin = 0;
  double diffA = 1;
  double diffB = 1;

  int NA = 1000;
  int NB = 1000;
  int Amax = Amin+(NA-1)*diffA;
  int Bmax = Amin+(NB-1)*diffA;

  TH1D *hist_signal = new TH1D("Signal","Signal",nBins,xmin,xmax);
  TH1D *hist_background = new TH1D("Background","Background",nBins,xmin,xmax);
  TH1D *hist_data = new TH1D("Data = Signal + Background","Data",nBins,xmin,xmax);
  TH2D *hist_contour = new TH2D("Contour","Contour",NA,Amin,Amax,NB,Bmin,Bmax);

  //Systematics
  if(gRandom) delete gRandom;
  gRandom = new TRandom3(0);

  while(count<Nbins_sbg){
    double maf=gRandom->Uniform(xmin,xmax);
    double pom=gRandom->Gaus(mean,sigma);
    if(maf<xmax && maf>xmin && pom<xmax && pom>xmin){
      hist_background->Fill(maf);
      hist_signal->Fill(pom);
      hist_data->Fill(maf);
      hist_data->Fill(pom);
      count=count+1;
      }
	}

  double A[NA];
  double B[NB];
  double sum;
  double sigma_sigma = 2*sigma*sigma;

  for(int i = 0; i < NA; i++){
  	A[i] = Amin + i*diffA;

  	for(int j = 0; j < NB; j++){
  		B[j] = Bmin + j*diffB;
  		sum = 0;
  		for(int k = 1; k <= nBins; k++){

          double x_k = hist_data->GetBinCenter(k); // hodota v kte binu
          double N_k =hist_data->GetBinContent(k); // pocet eventu v binu
          double D_k = A[i]*TMath::Exp(-(x_k-mean)*(x_k-mean)/(sigma_sigma))+B[j];
          sum = sum + N_k*TMath::Log(D_k)-D_k;
  			}

            hist_contour->SetBinContent(hist_contour->FindBin(A[i], B[j]), sum);
            }
  	}

    draw(hist_signal,hist_background,hist_data,1);
    contour(hist_contour);

    float A_max = Maximum(hist_contour,0);
    float B_max = Maximum(hist_contour,1);

    cout << "Maximum A je: "<<A_max*diffA<< " Maximum B je: " << B_max*diffB <<endl;

    float Bevents = B_max*(xmax-xmin);
    float Dataevents = hist_data->Integral(xmin,xmax);
    float Aevents = Dataevents-Bevents;

    cout << "#Bevents je: "<<Bevents<< " #Aevents je: " << Aevents <<" #Dataevents: "<<Dataevents<<" Ratio B/S : "<<Bevents/Aevents<<endl;

  }
