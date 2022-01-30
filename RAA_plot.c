void RAA_plot(){

  Double_t xAxis1[31] = {-100, -80, -60, -40, -20, -10, -5, -3, -1, 1, 3, 4, 5, 6, 7, 8, 9, 10, 12, 14, 16, 18, 20, 25, 30, 40, 50, 60, 70, 80, 100};
  TH1D *hpp = new TH1D("hpp","measured data pp",30, xAxis1);
  hpp->SetBinContent(9,1.118516);
  hpp->SetBinContent(10,0.1920894);
  hpp->SetBinContent(11,0.007120608);
  hpp->SetBinContent(12,0.001577695);
  hpp->SetBinContent(13,0.0004146155);
  hpp->SetBinContent(14,0.0001296542);
  hpp->SetBinContent(15,4.644327e-05);
  hpp->SetBinContent(16,2.057404e-05);
  hpp->SetBinContent(17,9.958874e-06);
  hpp->SetBinContent(18,4.334099e-06);
  hpp->SetBinContent(19,1.442999e-06);
  hpp->SetBinContent(20,5.654266e-07);
  hpp->SetBinContent(21,2.407638e-07);
  hpp->SetBinContent(22,1.094223e-07);
  hpp->SetBinContent(23,3.354351e-08);
  hpp->SetBinContent(24,6.04691e-09);
  hpp->SetBinContent(25,7.423793e-10);
  hpp->SetBinContent(26,3.350337e-11);
  hpp->SetBinError(9,0.0006107765);
  hpp->SetBinError(10,0.0003700558);
  hpp->SetBinError(11,4.49924e-05);
  hpp->SetBinError(12,1.821935e-05);
  hpp->SetBinError(13,7.447584e-06);
  hpp->SetBinError(14,2.939223e-06);
  hpp->SetBinError(15,1.155389e-06);
  hpp->SetBinError(16,5.875219e-07);
  hpp->SetBinError(17,1.855936e-07);
  hpp->SetBinError(18,1.141886e-07);
  hpp->SetBinError(19,3.0114e-08);
  hpp->SetBinError(20,5.526975e-09);
  hpp->SetBinError(21,2.701332e-09);
  hpp->SetBinError(22,1.525211e-09);
  hpp->SetBinError(23,6.343152e-10);
  hpp->SetBinError(24,9.460454e-11);
  hpp->SetBinError(25,2.709424e-11);
  hpp->SetBinError(26,4.271672e-13);
  hpp->SetEntries(1000);


  TH1D *hAA = new TH1D("hAA","measured data AA",30, xAxis1);
  hAA->SetBinContent(12,0.0004135252);
  hAA->SetBinContent(13,0.0001128615);
  hAA->SetBinContent(14,5.329878e-05);
  hAA->SetBinContent(15,2.740152e-05);
  hAA->SetBinContent(16,1.417387e-05);
  hAA->SetBinContent(17,8.002939e-06);
  hAA->SetBinContent(18,3.765388e-06);
  hAA->SetBinContent(19,1.418623e-06);
  hAA->SetBinContent(20,5.724047e-07);
  hAA->SetBinContent(21,2.704585e-07);
  hAA->SetBinContent(22,1.308473e-07);
  hAA->SetBinContent(23,3.482406e-08);
  hAA->SetBinContent(24,5.969517e-09);
  hAA->SetBinContent(25,8.159869e-10);
  hAA->SetBinContent(26,1.976192e-11);
  hAA->SetBinError(12,1.06535e-06);
  hAA->SetBinError(13,4.957864e-07);
  hAA->SetBinError(14,3.151375e-07);
  hAA->SetBinError(15,2.136749e-07);
  hAA->SetBinError(16,1.451485e-07);
  hAA->SetBinError(17,1.005136e-07);
  hAA->SetBinError(18,5.988379e-08);
  hAA->SetBinError(19,3.659272e-08);
  hAA->SetBinError(20,2.201622e-08);
  hAA->SetBinError(21,1.464843e-08);
  hAA->SetBinError(22,9.694462e-09);
  hAA->SetBinError(23,3.96629e-09);
  hAA->SetBinError(24,1.470905e-09);
  hAA->SetBinError(25,4.142606e-10);
  hAA->SetBinError(26,1.499208e-11);
  hAA->SetEntries(60);




  // Extracting spectra
  // vybereme si 3 biny v rozmezi 20 az 40 abychom kodily v loopech a bylo to jasne
  // jde o to ze budeme vykreslovat 3 obrazky pro kazdy bin a pro 30 puvodnich binu by to byl mazec ...

  // we choose 3 bins for pt from 20 to 40 so we have to code it in loop (like a real case)
  // but since we want to visualize all plots, 3 are enough (out of 30 above)

  const Int_t iNumBins = 3;

  Double_t dAxis[iNumBins+1] = {20,25,30,40};

  TH1D* h_data_pp = new TH1D("h_data_pp","pp",iNumBins,dAxis);
  TH1D* h_data_AA = new TH1D("h_data_AA","AA",iNumBins,dAxis);
  TH1D* h_data_RAA = new TH1D("h_data_RAA","RAA",iNumBins,dAxis);

  // pro kazdy bin nactu hodnoty a chyby v danem binu a ulozim to do histogramu,
  // spoctu pomer hodnot a ulozim do histogramu - ale co s chybama v binech?

  // load the values to new histograms for the 3 bins
  // we have ratio = the means, but how to calculate the sigma?

  for(Int_t i = 1; i <= iNumBins; i++)
  {
    Double_t dValuePP = hpp->GetBinContent(22+i);
    Double_t dSigmaPP = hpp->GetBinError(22+i);
    h_data_pp->SetBinContent(i,dValuePP);
    h_data_pp->SetBinError(i,dSigmaPP);

    Double_t dValueAA = hAA->GetBinContent(22+i);
    Double_t dSigmaAA = hAA->GetBinError(22+i);
    h_data_AA->SetBinContent(i,dValueAA);
    h_data_AA->SetBinError(i,dSigmaAA);

    Double_t dValueRAA = dValueAA/dValuePP;

    h_data_RAA->SetBinContent(i,dValueRAA);
  }


  //Systematics
   if(gRandom) delete gRandom;
  gRandom = new TRandom3(0);

  TCanvas* cAsym = new TCanvas("cAsym","Systematics",1300,600);
  cAsym->Divide(3,3);		// 3 obrazky pro 3 biny

  // a trick - for meaningful values - scaling of the values to be roughly 1.
  // can do this as this will cancel in the ratio
  Double_t Factor[iNumBins] = {1e8,1e9,1e10};
  // then we work with multiplied values: h_data_pp->GetBinContent(k+1)*Factor[k] , h_data_pp->GetBinError(k+1)*Factor[k] ,  kde k = 0,1,2


  // 3 histograms for 3 bins - for random values and their ratio = Err; THIS IS WHAT I WANT TO PLOT and RAA histogram fitted
  TH1D* h_toy_pp[iNumBins]; //pp
  TH1D* h_toy_AA[iNumBins]; //AA
  TH1D* h_toy_RAA[iNumBins]; //ratio for fit

  // histogram binning and ranges
  Double_t MinBinVal = -5;
  Double_t MaxBinVal = 20;
  Double_t MinBinErr = -1;
  Double_t MaxBinErr = 3.5;
  Int_t    nBinsVal = 1001;
  Int_t    nBinsErr = 1001;

  	//  define funcion of Gauss at the given range MinBinVal-MaxBinVal - for value generation
	// we need these functions to generate toy data - mean is the ration, sigma of gauss = bin error of data
// ** CODE HERE
TF1* gauss_pp = new TF1("gauss_pp","(0.39894228040/[1])*exp(-0.5*(x-[0])*(x-[0])/([1]*[1]))",MinBinVal,MaxBinVal);
TF1* gauss_AA = new TF1("gauss_AA","(0.39894228040/[1])*exp(-0.5*(x-[0])*(x-[0])/([1]*[1]))",MinBinVal,MaxBinVal);

      //  define 3 functions leftGauss & rightGauss for fitting both sides of ratio histogram
      //  3x because you have 3 bins to fit
// ** CODE HERE

TF1* leftGauss = new TF1("leftGauss","((0.39894228040*[2])/[1])*exp(-0.5*(x-[0])*(x-[0])/([1]*[1]))",MinBinErr,MaxBinErr);
TF1* rightGauss = new TF1("rightGauss","((0.39894228040*[2])/[1])*exp(-0.5*(x-[0])*(x-[0])/([1]*[1]))",MinBinErr,MaxBinErr);


  	// here store the results
  	Double_t RAAval[iNumBins];	// RAA from data - the mean
        RAAval[0] = h_data_RAA->GetBinContent(1);
        RAAval[1] = h_data_RAA->GetBinContent(2);
        RAAval[2] = h_data_RAA->GetBinContent(3);
	//MAIN TASK - FILL THESE VARIABLES
  	Double_t sigma_left[iNumBins];
  	Double_t sigma_right[iNumBins];


  // loop over bins:
  for(Int_t k = 0; k < iNumBins; k++)
  {

      h_toy_pp[k] = new TH1D("h_toy_pp","h_toy_pp",nBinsVal,MinBinVal,MaxBinVal);
      h_toy_AA[k] = new TH1D("h_toy_AA","h_toy_AA",nBinsVal,MinBinVal,MaxBinVal);
      h_toy_RAA[k] = new TH1D("h_toy_RAA","h_toy_RAA",nBinsErr,MinBinErr,MaxBinErr);

      // since we loop from 0, the first bin has index 1 => k+1 bin
      // here is the trick to set every value around 1 using the binContent*Factor[k] and binErr*Factor[k]Factor
      // this is the crucial step - to set parameters of the gauss according to the data


// ** CODE HERE
Double_t hodnota_pp = (h_data_pp->GetBinContent(k+1))*Factor[k];
Double_t hodnota_AA = (h_data_AA->GetBinContent(k+1))*Factor[k];

Double_t binErr_pp = (h_data_pp->GetBinError(k+1))*Factor[k];
Double_t binErr_AA = (h_data_AA->GetBinError(k+1))*Factor[k];

gauss_pp->SetParameters(hodnota_pp,binErr_pp);
gauss_AA->SetParameters(hodnota_AA,binErr_AA);

      for(Int_t p = 0; p < 100000+1; p++)
      {

        	// fill toy histograms with randomly sampled numbers

// ** CODE HERE
h_toy_pp[k]->Fill(gauss_pp->GetRandom(MinBinVal,MaxBinVal));
h_toy_AA[k]->Fill(gauss_AA->GetRandom(MinBinVal,MaxBinVal));
h_toy_RAA[k]->Fill((gauss_AA->GetRandom(MinBinVal,MaxBinVal))/(gauss_pp->GetRandom(MinBinVal,MaxBinVal)));

      }

  // so we have pseudodata for the ratio - let's fit
      // plotting toy data:
      Int_t N = k*3;
      cAsym->cd(N+1);
      h_toy_pp[k]->Draw();
      cAsym->cd(N+2);
      h_toy_AA[k]->Draw();
      cAsym->cd(N+3);
      h_toy_RAA[k]->Draw();

	// set parameters for 3x gauss left, 3x gause right (i.e. for 3 bins)
	// FixParameter()   for mean
	// dont forget to set the normalization!!! -- it depends on tthe amount of toy data
	// fit the toy RAA
	// fill the pre-prepared variables for sigmas

// ** CODE HERE

  leftGauss->SetLineColor(kRed);
  rightGauss->SetLineColor(kGreen);
	leftGauss->SetParameters(h_toy_RAA[k]->GetMean(),h_toy_RAA[k]->GetRMS(),h_toy_RAA[k]->Integral("width"));
	rightGauss->SetParameters(h_toy_RAA[k]->GetMean(),h_toy_RAA[k]->GetRMS(),h_toy_RAA[k]->Integral("width"));

  leftGauss->FixParameter(0,h_toy_RAA[k]->GetMean());
	rightGauss->FixParameter(0,h_toy_RAA[k]->GetMean());

	h_toy_RAA[k]->Fit("leftGauss","leftGauss","",MinBinErr,h_toy_RAA[k]->GetMean());
	h_toy_RAA[k]->Fit("rightGauss","rightGauss+","",h_toy_RAA[k]->GetMean(),MaxBinErr);

	sigma_left[k]=leftGauss->GetParameter(1);
	sigma_right[k]=rightGauss->GetParameter(1);


  }   // end of loop over bins


  std::cout << "### Results Review ###" << std::endl;
  for(Int_t m = 0; m < iNumBins; m++)
  {
      std::cout << "RAA:" << RAAval[m] <<" AsymLeft:" << sigma_left[m] << "  AsymRight:" << sigma_right[m] << endl;
  }
  std::cout << "=======================" << std::endl;

  //Konverzion of histogram to points in graph:
  Double_t dX[3] = {22.5,27.5,35.};
  Double_t dXEL[3] = {2.5,2.5,5.};
  Double_t dXEH[3] = {2.5,2.5,5.};
  Double_t dYEH[3];
  Double_t dYEL[3];

  for(Int_t i = 0; i < 4; i++)
  {
    dYEL[i] = TMath::Abs(sigma_left[i]);
    dYEH[i] = TMath::Abs(sigma_right[i]);
  }


  TGraphAsymmErrors* geRAA = new TGraphAsymmErrors(3,dX,RAAval,dXEH,dXEL,dYEH,dYEL);
  geRAA->SetTitle("RAA");


  //Plotting
  TCanvas* cBin = new TCanvas("cBin","Bins",1300,600);

  cBin->Divide(2,1);

  cBin->cd(1);
  cBin->GetPad(1)->SetLogy();

  h_data_pp->SetMarkerColor(3);
  h_data_pp->SetLineColor(3);
  h_data_pp->Draw("e");

  h_data_AA->SetMarkerColor(2);
  h_data_pp->SetLineColor(2);
  h_data_AA->Draw("same");


  cBin->cd(2);
  //h_data_RAA->Draw("E");
  geRAA->Draw("ap");














}
