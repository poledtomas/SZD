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
void contour(TH2D* hist1,TH2D* hist2, float k) {

  if (k==0) {
    TCanvas* contour = new TCanvas("contour","contour",800,800);
    gStyle->SetPalette(kRainBow);
    hist1->SetStats(0);
    hist1->SetContour(100);
    hist1->Draw("COLZ");
  } else {
    TCanvas* contour1 = new TCanvas("Metoda MLM","Metoda MLm",800,800);
    gStyle->SetPalette(kRainBow);
    hist1->SetStats(0);
    hist1->SetContour(100);
    hist1->Draw("COLZ");

    TCanvas* contour2 = new TCanvas("Metoda CHI2","Metoda CHI2",800,800);
    gStyle->SetPalette(kRainBow);
    hist2->SetStats(0);
    hist2->SetContour(100);
    hist2->Draw("COLZ");
  }

}
void chi_fit(){

  double mean_teo = 50, sigma_teo = 20, xmin = 0,xmax = 100;
  int nBins = 100, Nbins_sbg = 100000, pocet_mean = 1000, pocet_sigma = 1000, count = 0;

  TH1D *hist_data = new TH1D("Gauss","Data+Background",nBins,xmin,xmax);
  TH1D *hist_gauss = new TH1D("Gauss","Signal",nBins,xmin,xmax);
  TH1D *hist1 = new TH1D("Hist1","Uniform1",nBins,xmin,xmax);
  TH1D *hist2 = new TH1D("Hist2","Uniform2",nBins,xmin,xmax);

  if(gRandom) delete gRandom;
  gRandom = new TRandom3(0);

  while(count<Nbins_sbg){
    double maf=gRandom->Uniform(xmin,xmax);
    double pom=gRandom->Uniform(xmin,xmax);
    double top=gRandom->Gaus(mean_teo,sigma_teo);
    if(maf<xmax && maf>xmin && pom<xmax && pom>xmin&& top<xmax && top>xmin){
      hist1->Fill(maf);
      hist2->Fill(pom);
      hist_gauss->Fill(top);
      hist_data->Fill(maf);
      hist_data->Fill(top);
      hist_data->Fill(pom,-1);
      count=count+1;
    }
	}

  count = 1;
  while(count<=nBins){
    if(hist_data->GetBinContent(count)<0) hist_data->SetBinContent(count,0);
    count=count+1;
  }

  double vek_mean[pocet_mean];
  double vek_sigma[pocet_sigma];

  double sigma_min =15;
  double mean_min = 40;
  double diff= 0.01;

  double mean_max = mean_min+(pocet_mean-1)*diff;
  double sigma_max = sigma_min+(pocet_sigma-1)*diff;

  Double_t sum_chi,sum, max_sum, max_sum_chi;
  max_sum = 10e6;
  max_sum_chi = 10e6;
  double mean_mlm, sigma_mlm,mean_chi, sigma_chi;

  TH2D *contour_mlm = new TH2D("contour_mlm","contour_mlm",pocet_mean,mean_min,mean_max,pocet_sigma,sigma_min,sigma_max);
  TH2D *contour_chi = new TH2D("contour_chi","contour_chi",pocet_mean,mean_min,mean_max,pocet_sigma,sigma_min,sigma_max);

  for (int i = 0; i < pocet_mean; i++){
    vek_mean[i] = mean_min + i*diff;

    for (int j = 0; j < pocet_sigma; j++){

      vek_sigma[j] = sigma_min+j*diff;
      sum = 0;
      sum_chi = 0;

      for (int k = 1; k < nBins; k++){

        double x_k = hist_data->GetBinCenter(k);
        double A = hist_data->GetMaximum();
        double N_k = hist_data->GetBinContent(k);
        double sigma_error_k = hist_data->GetBinError(k);
        double D = A*TMath::Exp(-(x_k-vek_mean[i])*(x_k-vek_mean[i])/(2*vek_sigma[j]*vek_sigma[j]));
        sum =sum+ N_k*(TMath::Log(2.5*vek_sigma[j])+((x_k-vek_mean[i])/vek_sigma[j])*((x_k-vek_mean[i])/vek_sigma[j]));
        sum_chi = sum_chi+ ((N_k-D)/sigma_error_k)*((N_k-D)/sigma_error_k);
      }

      contour_mlm->SetBinContent(contour_mlm->FindBin(vek_mean[i],vek_sigma[j]), sum);
      contour_chi->SetBinContent(contour_mlm->FindBin(vek_mean[i],vek_sigma[j]), sum_chi);

      if(sum<max_sum){
        max_sum=sum;
        mean_mlm=vek_mean[i];
        sigma_mlm=vek_sigma[j];
      }

      if(sum_chi<max_sum_chi){
        max_sum_chi=sum_chi;
        mean_chi=vek_mean[i];
        sigma_chi=vek_sigma[j];
      }

    }
  }

  cout << "Vstupni paramtery mean a sigma: (" <<  mean_teo << "," <<  sigma_teo <<")"<<"\n"<<endl;
  cout << "Mean a sigma z metody MLM: (" << mean_mlm << "," <<  sigma_mlm <<")"<<"\n"<< "Pomer (ratio) mean_mlm/mean_teo: "<<mean_mlm/mean_teo<<" Pomer(ratio) sigma_mlm/sigma_teo: "<<sigma_mlm/sigma_teo<<"\n" <<endl;
  cout << "Mean a sigma z metody CHI2: (" << mean_chi << "," <<  sigma_chi <<")"<<"\n"<<"Pomer (ratio) mean_chi/mean_teo: "<<mean_chi/mean_teo<<" Pomer(ratio) sigma_chi/sigma_teo: "<<sigma_chi/sigma_teo<< endl;
  contour(contour_mlm,contour_chi,1);
  draw(hist1,hist_gauss,hist_data,1);
}
