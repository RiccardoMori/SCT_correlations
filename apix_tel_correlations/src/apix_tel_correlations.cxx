
#include "internal/inStripEfficiency.hh"
#include <iostream>
#include <sstream>
#include <thread>

#include "processorBase.h"
#include "TError.h"




#include "tclap/CmdLine.h"
#include "xml_helpers/xml_fileList.hh"
#include "TApplication.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TProfile.h"
#include "TF1.h"
#include "TBrowser.h"
#include "TH2.h"
#include "geometry/setup_description.hh"
#include "s_file_base.h"
#include "sct_processors.h"
#include "SCT_helpers.h"
#include "internal/hit_efficiency.hh"
#include "internal/inStripClusterSize.hh"

#include "internal/residual_efficienct.hh"
#include "processors/find_second_nearest_strip.hh"
#include "processors/find_nearest.hh"
#include "sct_processors.h"
#include "processors/modulo.hh"
#include "processors/find_nearest_strip.hh"
#include "s_plot_collection.h"
#include "plane_def.h"


using namespace xml_util;

using namespace rapidxml;

std::vector<TCanvas*> gCanvas;
TBrowser* gBrowser = NULL;


using namespace std;
using namespace TCLAP;

void remove_root_printouts() {
  gErrorIgnoreLevel = kError;  // ignoring root printouts (replace of histograms) 
}

struct  inParam {
  int argc;
  char **argv;
};


void add2File(TFile* file_, TH2* h, const char* name) {
  h->SetTitle(name);
  h->SetName(name);
  file_->Add(h);
}
using namespace sct_type;
int asyncMain(void *arg) {
  int argc_asdasddsa = 1;
  char **argv_sadsda = new char*[1];
  argv_sadsda[0] = "dummy";
#ifdef _DEBUG
  TApplication theApp("App", &argc_asdasddsa, argv_sadsda);
#endif // _DEBUG

  inParam* para = static_cast<inParam *>(arg);
  int argc = para->argc;
  char **argv = para->argv;


  CmdLine cmd("ProcessFile", ' ', "0.1");

  ValueArg<std::string> FileNameArg("i", "inFile", "xml filename", true, "", "string");
  cmd.add(FileNameArg);
  ValueArg<std::string>  inPath("p", "inPath", "path to the root files", true, ".", "string");
  cmd.add(inPath);
  ValueArg<std::string>  output_path("o", "outPath", "output path", false, "output.root", "string");
  ValueArg<int> element("e", "element", "element of interest  in the XML file", true, 1, "int");
  cmd.add(element);
  cmd.parse(argc, argv);  //terminates on error


  auto m_input_files_xml = std::make_shared<xmlImputFiles::XML_imput_file>(FileNameArg.getValue().c_str());

  TFile * out_file = new TFile(output_path.getValue().c_str(), "recreate");

  auto file_path = (inPath.getValue() + "/" + m_input_files_xml->fileList().at(element.getValue()).name());
  std::cout << "opening file path: " << file_path << std::endl;
    sct_files::fitter_file fitter_file_(file_path.c_str(), m_input_files_xml->globalConfig().gearFile().c_str());

  auto pl1 = fitter_file_.get_collection();
  pl1->setOutputFile(out_file);

  auto apix_local = sct_corr::processor::convert_zs_data_to_hits_GBL(
    fitter_file_.apix_zs_data(),
    *fitter_file_.get_gear()->detector.layer_by_ID(20),
    s_plot_prob().doNotSaveToDisk()
    );
  auto apix_global = sct_corr::processor::convert_local_to_global(
    apix_local,
    *fitter_file_.get_gear()->detector.layer_by_ID(20),
    s_plot_prob().doNotSaveToDisk()
    );
  auto apix_local_8 = sct_corr::processor::convert_global_to_local(
    apix_global,
    *fitter_file_.get_gear()->detector.layer_by_ID(8),
    s_plot_prob().doNotSaveToDisk()
    );
  auto corr_xx = sct_corr::processor::correlation(apix_local_8.getX_def(), fitter_file_.DUT_fitted_local_GBL().getX_def(), s_plot_prob().SaveToDisk());
  auto corr_yy = sct_corr::processor::correlation(apix_local_8.getY_def(), fitter_file_.DUT_fitted_local_GBL().getY_def(), s_plot_prob().SaveToDisk());

  pl1->loop();
  auto h =SCT_helpers::Draw<TH2>(corr_xx,S_DrawOption());
  auto h_local = (TH2*)h->Clone("corr_xx");
  auto p = SCT_helpers::LinearFit_Of_Profile(h_local, sct_type::procent_t(20));
  h_local->Draw("colz");
  p.Draw("same");
  p.Print();


  new TCanvas();


  auto h1 = SCT_helpers::Draw<TH2>(corr_yy, S_DrawOption());
  auto h1_local = (TH2*)h1->Clone("corr_yy");
  auto p1 = SCT_helpers::LinearFit_Of_Profile(h1_local, sct_type::procent_t(20));
  h1_local->Draw("colz");
  p1.Draw("same");
  p1.Print();
  
  out_file->Add(h_local);
  out_file->Add(h1_local);
  out_file->Write(); 
#ifdef _DEBUG
  theApp.Run();
#endif // _DEBUG
  return 0;
}
int main(int argc, char **argv) {

  remove_root_printouts();


  inParam para;
  para.argc = argc;
  para.argv = argv;


  asyncMain(&para);




  return 0;

}

