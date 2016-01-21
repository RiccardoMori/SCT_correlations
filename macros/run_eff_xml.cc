//To Run it first upload the lib /gSystem->Load("libSCT_plot.dylib");
s_process_files* run_eff_xml(){
    
    
    s_process_files* p=new s_process_files();
    xmlImputFiles::ADDRun(p,"SCT_correlations/data/dataFileList.xml","/Users/lucreziastellabruni/");
        return p;
}







