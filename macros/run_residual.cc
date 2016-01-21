//To Run it first upload the lib /gSystem->Load("libSCT_plot.dylib");
s_process_files* run_residual(){
    
    const char* inFileName ="/Users/lucreziastellabruni/fitter694.root";

    TFile* _file1 = new TFile("Residual_run703.root", "RECREATE");
    
    s_process_files* p=new s_process_files();
    p->setOutputName("Residual_run703.root");
    p->SetNumberOfBins(400);
    p->AddCut(S_YCut(-4, 2));
    p->setActiveArea(280, 360);
    p->SetRotation(0.0171112);
    p->SetPosition(3.16097e+002, 0);
    p->push_files(inFileName, 50, 703, 350);
    p->process();
    //p->DrawResidual(-2, 2);
   // p->DrawResidualVsMissingCordinate(-2,2);
    p->Draw_Efficinecy_map();
    return p;
}







