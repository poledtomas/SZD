#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <stdio.h>
#include <stdlib.h>

void draw(TH1D* histuhel,TH1D* histx){

  TCanvas* plot = new TCanvas("plot","Hisogramy",1000,500);
  plot->Divide(2,1);

  plot->cd(1);
  histuhel->Draw();
  plot->cd(2);
  histx->Draw();
}

void Maximum(TH2D* hist_contour, float diffx,float diffy){

  double bin_max = hist_contour->GetMaximumBin();
  int x_max,y_max,z_max;
  hist_contour->GetBinXYZ(bin_max,x_max,y_max,z_max);
  cout << "Max hodnota je: "<<bin_max<<" V bode: ("<<x_max*diffx<<","<<y_max*diffy<<")"<<endl;
}

void contour(TH2D* hist) {
  TCanvas* contour = new TCanvas("contour","contour",600,600);
  gStyle->SetPalette(kRainBow);
  hist->SetStats(0);
  hist->SetContour(100);
  hist->Draw("COLZ");

}
void zaric (){

  Double_t y_0 = 1.0;
  Double_t x_0 = 3.0;
  Double_t xmax = 5.0;
  Double_t xmin = 0;
  Double_t ymax = 5.0;
  Double_t ymin = 0;

  Int_t Nx = 80;
  Int_t Ny = 80;
  Double_t diffx = (xmax-xmin)/Nx;
  Double_t diffy = (ymax-ymin)/Ny;

  TH1D* histuhel=new TH1D("histuhel","histogram uhel",100,-3.,3.);
  TH1D* histx = new TH1D("histx","histogram x",100,xmin,xmax);
  TH2D *histxy = new TH2D("histxy","Contour plot",Nx,xmin,ymax,Ny,ymin,ymax);

  UInt_t nevents =  10000;
  Int_t count = 0;

  vector<double> xpozice;

  while (count < nevents) {

    Double_t uhel = gRandom->Uniform(0,TMath::Pi());
    uhel =uhel - (TMath::Pi())/2;
    histuhel->Fill(uhel);

    Double_t x = y_0*TMath::Tan(uhel)+x_0;
    if (xmin < x && xmax >x) {
      xpozice.push_back(x);
      count=count+1;
      histx->Fill(x);
    }

  }

  Double_t x[Nx];
  Double_t y[Ny];

  for (int i  = 0;  i< Nx; i++) {
    x[i]=xmin+i*diffx;

    for (int j = 0; j < Ny; j++) {
      y[j]=ymin+j*diffy;

      Double_t sum = 0;

      for(Int_t l = 1; l <= xpozice.size(); l++){
        sum=sum+TMath::Log((y[j])/((xpozice[l]-x[i])*(xpozice[l]-x[i])+y[j]*y[j]));
      }
      histxy->SetBinContent(histxy->FindBin(x[i], y[j]), sum);
                 }
    }

    draw(histuhel,histx);
    contour(histxy);
    Maximum(histxy,diffx,diffy);

}
