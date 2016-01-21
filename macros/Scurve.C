#define Scurve_cxx
#include "Scurve.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

void Scurve::Loop()
{
    gSystem->Load("../lib/libSCT_plot.dylib");
    
    gSystem->Load("/Users/lucreziastellabruni/Dropbox/SCt_correlations/landau_gauss/lib/libLandauGaussFit.dylib");
    Double_t x[14], y[14], dy[14], dx[14], d[14];
    int n= 14;
    
    TH2D *h2= new TH2D("h2","",50,0, 300,50,0,1);
    if (fChain == 0) return;
    
    Long64_t nentries = fChain->GetEntriesFast();
    
    Long64_t nbytes = 0, nb = 0;
    for (Long64_t jentry=0; jentry<nentries;jentry++) {
        Long64_t ientry = LoadTree(jentry);
        if (ientry < 0) break;
        nb = fChain->GetEntry(jentry);   nbytes += nb;
        
        x[jentry] = Threshold;
        y[jentry] = total_efficiency;
        std::cout<<"thr:  "<<x[jentry]<<" eff:  "<<y[jentry]<<std::endl;
 
        h2->Fill(x[jentry],y[jentry]);
        /* dx[jentry]= x[jentry+1]-x[jentry];
         dy[jentry]= abs(y[jentry+1]-y[jentry]);
         d[jentry]=(dy[jentry]/dx[jentry]);
         
         std::cout<<"d= "<<d[jentry]<<std::endl;
         h1->Fill(d[jentry]);*/
    }
    
    
    // then you have to define a fit object
    landgausFit fit;
    
    // and set the start parameter
    // by default the amplitude is set to 1
    fit.setStartAmplitude(1);
    fit.setLimits_Amplitude(0.7, 1.0);
    
    fit.setStartGaussSigma(30);
    fit.setStartLandauMean(110);
    fit.setStartLandauSigma(2.3);
    TCanvas *c4 = new TCanvas("c4","",1000, 400);
    
    
    TGraph *g = new TGraph(n,x,y);
    
    
    // then you can fit you TGraph object
    //fit(h2);  // fit the data set
    
    //when you have fitted the data set you can draw it.
    //fit.DrawfitFunction(); //see attachment
    
    //or print out the parameters
    //fit.printResults(cout);
    //TF1 *f1;
    //f1 = fit.getLandauGauss();
    //f1->SetParameter("Amplitude",-1);
    
    //TGraph *g1 = new TGraph(f1,"");
    h2->SetMarkerStyle(7);
    h2->SetMarkerColor(kRed);
    h2->SetMarkerSize(5);
    h2->Draw();
    
    
    TCanvas *c3 = new TCanvas("c3","",1000, 400);
    //TGraph *g1 = (TGraph*)f1->DrawDerivative();
    /*TH1D* h4 = new TH1D("h4", "h4", 20, 0, 0);
     TH1D* h5 = new TH1D("h5", "h5", 20, 0, 0);
     
     c3->Divide(3,1);
     c3->cd(1);
     pl->Draw(strip_nearest_fit_cut, S_DrawOption().draw_y().output_object(h4));
     c3->cd(2);
     pl->Draw(dut_fit_cut(), S_DrawOption().draw_y().output_object(h5));
     c3->cd(3);
     TH1D *h6=(TH1D*)h4->Clone();
     h6->Divide(h5);
     h6->Draw();
     
     
     TGraph *gr = new TGraph(n,x,y);
     TGraph *gr2 = new TGraph(n,x,d);
     gr2->Draw("AC*");
     gr->Draw("same AC*");
     */
    
}
