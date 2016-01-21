
void run(){
    gSystem->Load("../lib/libSCT_plot.dylib");
    gSystem->Load("/Users/lucreziastellabruni/Dropbox/SCt_correlations/landau_gauss/lib/libLandauGaussFit.dylib");
    const char* inFileName="/Users/lucreziastellabruni/fitter694.root";
    TFile* _file0 = new TFile(inFileName);
    TFile* _file1 = new TFile("output694.root","RECREATE");
    
    
    const char* name4 = "plot4";
    const char* strip = "strip";

    s_plane_collection strip_nearest = pl->addPlot(sct_plot::find_nearest_strip(y_axis_def, 1000,strip), sct_coll::DUT_hit(), sct_coll::DUT_fitted());
    //s_plane_collection hitmap__= pl->addPlot(sct_plot::cut_x_y(S_XCut(-2.8,5.2)+S_YCut(-3.6,4.2), ""),strip_nearest.get(2).getX_def(),strip_nearest.get(0).getY_def()  );
    //s_plane_collection hitmap__= pl->addPlot(sct_plot::cut_x_y(S_XCut(-6,6)+S_YCut(-6,6),""),strip_nearest.get(2) );
    
    s_plane_collection hitmap__= pl->addPlot(sct_plot::cut_x_y( S_XCut(-1,3)+S_YCut(-3.6,4.2),""),strip_nearest.get(2) );
   // s_plane_collection rot__= pl->addPlot(sct_plot::cut_x_y(S_XCut(-2.8,5.2)+S_YCut(-3.6,4.2), ""),strip_nearest.get(2).getX_def(),strip_nearest.get(0).getY_def()  );
    
    auto apix_true_hits = pl->addPlot(sct_plot::find_nearest( 0.1, 0.2), sct_coll::apix_fitted(), sct_coll::apix_hit()).get(1);
    
    auto dut_fitted_tracks = pl->addPlot(sct_plot::find_nearest( 1, 1), sct_coll::DUT_fitted(), apix_true_hits).get(1);
    
    s_plane_collection dut_fit_cut = pl->addPlot(sct_plot::cut_x_y(S_XCut(-1, 3) + S_YCut(-3.6, 4.2)), dut_fitted_tracks);
    
    S_plane_def strip_nearest_fit_cut = pl->addPlot(sct_plot::find_nearest_strip(y_axis_def, 1000), sct_coll::DUT_hit(), dut_fit_cut()).get(2);
    
    
    pl->loop();
    

    
    TCanvas *c1 = new TCanvas("c1","c1",600,900);
    c1->Divide(1,2);
    
    TH1D* h1 = new TH1D("h1", "h1", 100, -0.2, 0.2);
    TH2D* h2 = new TH2D("h2", "h2", 100,-5, 6, 100, 0, 0);
    
    c1->cd(2);
    
    pl->Draw(hitmap__(), S_DrawOption("colz", "", "y:x>>h2"));
    c1->cd(1);
    auto h3 = SCT_helpers::HistogrammSilhouette(h2, y_axis_def);
    h3->Fit("gaus");
    TF1 *fit1 = (TF1*)h3->GetFunction("gaus");
    
    Double_t max, Z_cut;
    max = fit1->GetParameter(0) ;
    Z_cut=max*0.2;
    std::cout<<" max: "<<max<<"  Z_cut: "<<Z_cut<<std::endl;
    SCT_helpers::CutTH2(h2, S_Cut_BinContent(Z_cut));
    
    h3->Write();
    h2->Write();
    
    TCanvas *c2 = new TCanvas("c2","show profile",600,900);
    c2->Divide(1,2);
    TProfile *prof1 = h2->ProfileX("_pfx1", 1, 100,"");
    prof1->SetLineColor(1);
    prof1->SetLineWidth(4);
    prof1->Fit("pol1","0");
    c2->cd(1);
    prof1->Draw();
    
    TF1 *f1 =(TF1*)prof1->GetFunction("pol1");
    c2->cd(2);

    auto f=sct_analyis::Draw_New_aligment(*_file0, f1);
    
    
    TCanvas *c3 = new TCanvas("c3","",1000, 400);
    TH2D* h4 = new TH2D("h4", "h4", 25, -3, 4 ,25, -1, 3 );
    TH2D* h5 = new TH2D("h5", "h5", 25, -3, 4,25,-1, 3);
    Double_t n1, n2;
    n1 = pl->Draw(strip_nearest_fit_cut, S_DrawOption().draw_x_VS_y().output_object(h4));
    c3->Divide(3,2);
    c3->cd(1);
    pl->Draw(strip_nearest_fit_cut, S_DrawOption().draw_x_VS_y().output_object(h4));
     std::cout << "h4 = " << n1 << std::endl;
    c3->cd(2);
    pl->Draw(dut_fit_cut(), S_DrawOption().draw_x_VS_y().output_object(h5));
    n2 =  pl->Draw(dut_fit_cut(), S_DrawOption().draw_x_VS_y().output_object(h5));
    std::cout << "h5 = " <<n2 << std::endl;
    std::cout<<"eff:  "<<n1/n2<<std::endl;
    h4->Sumw2();
    h5->Sumw2();
    c3->cd(3);
    TH2D *h6=(TH2D*)h4->Clone();
    h6->Divide(h5);
    h6->Draw("TEXT");//"colz");
    c3->cd(4);
    TH1D *px = h6->ProjectionX("px", 0, 25);
    px->Scale(0.04);
    px->Draw();
    c3->cd(5);
    TH1D *py = h6->ProjectionY("py", 0, 25);
    py->Scale(0.04);
    py->Draw();

    
    
    
    // then you have to define a fit object
    landgausFit fit;
    
    // and set the start parameter
    // by default the amplitude is set to 1
    fit.setStartAmplitude(100);
    fit.setLimits_Amplitude(95, 100);
    
    
    TTree *T = new TTree("ntuple","data from ascii file");
    T->ReadFile("test.dat","x:y:z");
    
    // after you read in the ascii file you have to convert it to a tgraph
    Int_t n= T->Draw("z:y","","*");
    
    TGraph *g = new TGraph(n,T->GetV2(),T->GetV1());
    
    // then you can fit you TGraph object
    fit(g);  // fit the data set
  
    //when you have fitted the data set you can draw it.
    fit.DrawfitFunction(); //see attachment
    
    //or print out the parameters
    fit.printResults(cout);
    TCanvas *c4 = new TCanvas("c4","",1000, 400);

    g->Draw("");
    
   /* TCanvas *c3 = new TCanvas("c3","",1000, 400);
    TH1D* h4 = new TH1D("h4", "h4", 20, 0, 0);
    TH1D* h5 = new TH1D("h5", "h5", 20, 0, 0);
    
    c3->Divide(3,1);
    c3->cd(1);
    pl->Draw(strip_nearest_fit_cut, S_DrawOption().draw_y().output_object(h4));
    c3->cd(2);
    pl->Draw(dut_fit_cut(), S_DrawOption().draw_y().output_object(h5));
    c3->cd(3);
    TH1D *h6=(TH1D*)h4->Clone();
    h6->Divide(h5);
    h6->Draw(); */
    
}







