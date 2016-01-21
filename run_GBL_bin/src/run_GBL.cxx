
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
#include "internal/exceptions.hh"


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



using namespace sct_type;
int asyncMain(void *arg) {


  inParam* para = static_cast<inParam *>(arg);
  int argc = para->argc;
  char **argv = para->argv;
#ifdef _DEBUG
  TApplication theApp("App", &argc, argv);
#endif // _DEBUG

  std::string path_ = "D:/GBL/DEVICE_1_ASIC_on_Position_7_Jim_350V/";
  std::string name_ = "run000688_";
  std::string name_suffix = "fitter";
  std::string extension_ = "root";
  std::string fullName = path_ + name_ + name_suffix + "." + extension_;

  //TFile * file_ = new TFile(fullName.c_str());
  //   rapidxml::file<> m_file("D:/GBL/DEVICE_1_ASIC_on_Position_7_Jim_350V/alignedGear-check-iter2-run000703_with_plane20.xml");
  //   rapidxml::xml_document<> m_doc;
  //   m_doc.parse<0>(m_file.data());
  // 
  //   auto gear = sct_corr::Xgear(m_doc.first_node("gear"));
  //  auto = sct_corr::load_gear("D:/GBL/DEVICE_1_ASIC_on_Position_7_Jim_350V/alignedGear-check-iter2-run000703_with_plane20.xml");

  // r_plot_collection pl(file_);
  //  pl.setOutputFile(out_file);

#if 1
  try {
  auto m_input_files_xml = std::make_shared<xmlImputFiles::XML_imput_file>("D:/data_2015_10/root/GBL_LS_Position_3_SHORTSTRIPS_395V.xml");
  TFile * out_file = new TFile("output.root", "recreate");
  sct_files::fitter_file fitter_file_("d:/data_2015_10/root/run000624_fitter.root", "D:/data_2015_10/root/alignedGear-check-iter1-run000624.xml");
  auto pl1 = fitter_file_.get_collection();
  pl1->setOutputFile(out_file);

  auto gbl_collection1 = fitter_file_.get_correlations_channel(
    S_YCut(-42010, 10000),
    residualCut_t(10000000),
    rot_angle_t(0),
    move_t(161),
    s_plot_prob("GBL").SaveToDisk()
    );
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
  auto apix_local_8= sct_corr::processor::convert_global_to_local(
    apix_global,
    *fitter_file_.get_gear()->detector.layer_by_ID(8),
    s_plot_prob().doNotSaveToDisk()
    );
  auto trueHits = sct_processor::find_nearest(
    apix_local_8,
   fitter_file_.DUT_fitted_local_GBL(),
    1, // residual cut x
    1,  // residual cut y
    s_plot_prob().SaveToDisk()
    ).getHitOnPlaneB();

  auto trueHitsInStrips = sct_corr::processor::convert_hits_to_zs_data_GBL(
    trueHits,
    *fitter_file_.get_gear()->detector.layer_by_ID(8),
    s_plot_prob().doNotSaveToDisk()
    );
  auto moved_true = sct_corr::processor::coordinate_transform_move(
    trueHitsInStrips,
    161,
    0,
    s_plot_prob().doNotSaveToDisk()
    );

    auto corr = sct_corr::processor::correlation(moved_true.getX_def(), fitter_file_.DUT_zs_data().getX_def(), s_plot_prob().SaveToDisk());
  auto corr1 = sct_corr::processor::residual(moved_true.getX_def(), fitter_file_.DUT_zs_data().getX_def(), s_plot_prob().SaveToDisk());
  pl1->loop(10000);
  new TCanvas();
  TH2D h1("h1", "DUT x Vs Tel y hist", 100, 0, 0, 100, 0, 0);
  SCT_helpers::Draw<TH2>(corr, S_DrawOption().opt_colz().draw_y_VS_x());
  new TCanvas();
 // SCT_helpers::Draw<TH1>(corr1, S_DrawOption().opt_colz().draw_x().cut_x(-5,5));
 SCT_helpers::Draw<TH2>(gbl_collection1.getResidualVSmissing(), S_DrawOption().opt_colz().draw_x_VS_y().cut_x(-5,5));
  gBrowser = new TBrowser();
#ifdef _DEBUG
  theApp.Run();
#endif // _DEBUG

  } catch (...) {

  return sct_corr::handleExceptions();
}

#endif
#if 0
  try {
  TFile * file_ = new TFile("d:/data_2015_10/root/run000624_fitter.root");

  TFile * file_1 = new TFile("Dummy.root", "recreate");
  auto pl1 = sct_corr::create_plot_collection();


  pl1->addFile(file_);
  pl1->setOutputFile(file_1);
  sct_corr::plane_def sz_dut(sct_type::collectionName_t("local_hit"), sct_type::ID_t(8), nullptr);
  sz_dut.set_plot_collection(pl1);
  sct_corr::plane_def sz_tel(sct_type::collectionName_t("zsdata_apix"), sct_type::ID_t(20), nullptr);
  sz_tel.set_plot_collection(pl1);


  // auto cuted = sct_corr::processor::cut_x_y(sz_tel, S_XCut(200) + S_YCut(200));
  // auto cuted_dut = sct_corr::processor::cut_x_y(sz_dut, S_XCut(1060, 1180));
  auto corr = sct_corr::processor::correlation(sz_dut.getX_def(), sz_tel.getY_def(), s_plot_prob().SaveToDisk());
  auto corryx = sct_corr::processor::correlation(sz_dut.getY_def(), sz_tel.getX_def(), s_plot_prob().SaveToDisk());
  auto corrxx = sct_corr::processor::correlation(sz_dut.getX_def(), sz_tel.getX_def(), s_plot_prob().SaveToDisk());
  auto corryy = sct_corr::processor::correlation(sz_dut.getY_def(), sz_tel.getY_def(), s_plot_prob().SaveToDisk());

  auto sz_dut_normalized = sct_corr::processor::coordinate_transform(sz_dut, 1, 0, 1, 0, s_plot_prob().doNotSaveToDisk());
  auto sz_tel_normalized = sct_corr::processor::coordinate_transform(sz_tel, 1, 0, 1, 0, s_plot_prob().doNotSaveToDisk());
  auto resyy = sct_corr::processor::residual(sz_dut_normalized.getY_def(), sz_tel_normalized.getY_def());

  auto resxy = sct_corr::processor::residual(sz_dut_normalized.getX_def(), sz_tel_normalized.getY_def());

  auto resyx = sct_corr::processor::residual(sz_dut_normalized.getY_def(), sz_tel_normalized.getX_def());

  auto resxx = sct_corr::processor::residual(sz_dut_normalized.getX_def(), sz_tel_normalized.getX_def());

  pl1->loop(100000);
  new TCanvas();
  TH2D h1("h1", "DUT x Vs Tel y hist", 100, 0, 0, 100, 0, 0);
  SCT_helpers::Draw<TH2>(corr, S_DrawOption().opt_colz())->SetTitle("DUT x VS Tel y");
  //  auto p = SCT_helpers::LinearFit_Of_Profile(&h1, sct_type::procent_t(70));
  // h1.Draw("colz");
  // p.Draw("same");
  // p.Print();
  new TCanvas();

  TH2D h2("h2", "DUT x Vs Tel y hist", 1000, 0, 0, 1000, 0, 0);

  SCT_helpers::Draw<TH2>(corryx, S_DrawOption().opt_colz())->SetTitle("DUT y VS Tel x");
  new TCanvas();
  SCT_helpers::Draw<TH2>(corrxx, S_DrawOption().opt_colz())->SetTitle("DUT x VS Tel x");
  new TCanvas();
  SCT_helpers::Draw<TH2>(corryy, S_DrawOption().opt_colz())->SetTitle("DUT y VS Tel Y");
  new TCanvas();

  TH2D h3("h3", "DUT x Vs Tel y hist", 100, 0, 0, 100, 0, 0);

  SCT_helpers::Draw<TH2>(resyy, S_DrawOption().opt_colz().draw_x_VS_y())->SetTitle("DUT y - TEL y");


  new TCanvas();
  SCT_helpers::Draw<TH2>(resyx, S_DrawOption().opt_colz().draw_x_VS_y())->SetTitle("DUT y - TEL x");
  new TCanvas();
  SCT_helpers::Draw<TH2>(resxy, S_DrawOption().opt_colz().draw_x_VS_y())->SetTitle("DUT x - TEL y");
  new TCanvas();
  SCT_helpers::Draw<TH2>(resxx, S_DrawOption().opt_colz().draw_x_VS_y())->SetTitle("DUT x - TEL x");


  gBrowser = new TBrowser();
    } catch (...) {

    return sct_corr::handleExceptions();
  }
#ifdef _DEBUG
    theApp.Run();
#endif // _DEBUG



#endif
#if 0

  try {
    //TFile * file_ = new TFile("C:/Users/Argg/Documents/Neuer Ordner/gitHub/eudaq/bin/test000580.root");
    TFile * file_ = new TFile("d:/data_2015_10/root/run000634_fitter.root");
    TFile * file_1 = new TFile("Dummy.root", "recreate");
    auto pl1 = sct_corr::create_plot_collection();


    pl1->addFile(file_);
    pl1->setOutputFile(file_1);
    sct_corr::plane_def sz_dut(sct_type::collectionName_t("szData"), sct_type::ID_t(8), nullptr);
    sz_dut.set_plot_collection(pl1);
    sct_corr::plane_def sz_tel(sct_type::collectionName_t("szData"), sct_type::ID_t(7), nullptr);
    sz_tel.set_plot_collection(pl1);


    auto cuted = sct_corr::processor::cut_x_y(sz_tel, S_XCut(0));
    auto cuted_dut = sct_corr::processor::cut_x_y(sz_dut, S_XCut(1060, 1180));
    auto corr = sct_corr::processor::correlation(cuted_dut.getX_def(), cuted.getY_def(), s_plot_prob().SaveToDisk());
    auto corryx = sct_corr::processor::correlation(cuted_dut.getY_def(), cuted.getX_def(), s_plot_prob().SaveToDisk());
    auto corryy = sct_corr::processor::correlation(cuted_dut.getY_def(), cuted.getY_def(), s_plot_prob().SaveToDisk());

    auto sz_dut_normalized = sct_corr::processor::coordinate_transform(sz_dut, 1, 0, 1, 0, s_plot_prob().doNotSaveToDisk());
    auto sz_tel_normalized = sct_corr::processor::coordinate_transform(cuted, 3.19042, -978.631, 1, 0, s_plot_prob().doNotSaveToDisk());
    auto res = sct_corr::processor::residual(sz_dut_normalized.getX_def(), sz_tel_normalized.getX_def());

    auto corrxx = sct_corr::processor::correlation(sz_dut_normalized.getX_def(), sz_tel_normalized.getX_def(), s_plot_prob().SaveToDisk());
    pl1->loop(100000);
    new TCanvas();
    SCT_helpers::Draw<TH2>(corr, S_DrawOption().opt_colz());

    new TCanvas();
    SCT_helpers::Draw<TH2>(corryx, S_DrawOption().opt_colz())->SetTitle("DUT y VS Tel x");
    new TCanvas();
    TH2D h1("h1", "DUT x Vs Tel y hist", 100, 0, 0, 100, 0, 0);
    SCT_helpers::Draw<TH2>(corrxx, S_DrawOption().opt_colz().output_object(&h1))->SetTitle("DUT x VS Tel x");
    auto p = SCT_helpers::LinearFit_Of_Profile(&h1, sct_type::procent_t(40));
    h1.Draw("colz");
    p.Draw("same");
    p.Print();


    new TCanvas();
    SCT_helpers::Draw<TH2>(corryy, S_DrawOption().opt_colz())->SetTitle("DUT y VS Tel Y");
    new TCanvas();
    SCT_helpers::Draw(res, S_DrawOption().opt_colz().draw_x_VS_y());

    gBrowser = new TBrowser();
  } catch (...) {

    return sct_corr::handleExceptions();
  }
#ifdef _DEBUG
  theApp.Run();
#endif // _DEBUG
#endif // 0
#if 0


  sct_files::fitter_file file___(
    fullName.c_str(),
    "D:/GBL/DEVICE_1_ASIC_on_Position_7_Jim_350V/alignedGear-check-iter2-run000703_with_plane20.xml"
    );


  TFile * out_file = new TFile("output.root", "recreate");
  auto pl = file___.get_collection();
  pl->setOutputFile(out_file);



  auto gbl_collection = file___.get_correlations_channel(
    S_YCut(-42, -36),
    residualCut_t(10),
    rot_angle_t(0),
    move_t(-6.36702e-001),
    s_plot_prob("GBL").SaveToDisk()
    );



  auto trueHit_with_tdc_cut = sct_corr::processor::if_a_get_b(
    file___.DUT_TTC_data(),
    gbl_collection.getTotalTrueHits(),
    S_XCut(0, 100),
    s_plot_prob("withTDC_CUT")
    );


  auto  second_nearest = sct_processor::find_nearest_strip(
    gbl_collection.getTotalTrueHits(),
    file___.DUT_zs_data(),
    x_axis_def,
    3,
    "first"
    );

  auto mod = sct_processor::modulo(
    second_nearest.getHitOnPlaneA(),
    modulo_t(3),
    x_axis_def
    );


  auto rmap = sct_corr::processor::correlation(
    mod.getModulo().getX_def(),
    second_nearest.getResidual().getX_def()
    );


  sct_corr::inStripEfficiency instrip(
    gbl_collection.getTotalTrueHits(),
    second_nearest.getHitOnPlaneA(),
    S_XCut(280, 360),
    x_axis_def,
    modulo_t(3),
    s_plot_prob("inStripEffi")
    );


  sct_corr::hit_efficiency eff(
    gbl_collection.getTotalTrueHits(),
    second_nearest.getHitOnPlaneA(),
    s_plot_prob("effi")
    );
  sct_corr::inStripClusterSize cl_instrip(
    gbl_collection.getTrueHitsWithDUT(),
    file___.DUT_zs_data(), 10,
    x_axis_def,
    modulo_t(3),
    s_plot_prob("cluster_size_instrip").SaveToDisk()
    );

  sct_corr::residual_efficienct res_eff(
    gbl_collection.getTotalTrueHits(),
    second_nearest.getHitOnPlaneB(),
    S_XCut(280, 360),
    stripNr_t(400),
    x_axis_def,
    s_plot_prob("Res_efficiency")
    );

  pl->loop();

  //pl->loop();

  gCanvas.push_back(new TCanvas());
  SCT_helpers::Draw(rmap, S_DrawOption().draw_y_VS_x());

  gCanvas.push_back(new TCanvas());

  res_eff.Draw();
  SCT_helpers::saveTH1_as_txt(*res_eff.getEfficiency_map(), (path_ + name_ + "residual_efficiency" + "." + "txt").c_str());

  gCanvas.push_back(new TCanvas());
  instrip.Draw(S_DrawOption());
  std::cout << "Efficiency " << instrip.get_total_efficiency() << std::endl;
  SCT_helpers::saveTH1_as_txt(*instrip.getEfficiency_map(), (path_ + name_ + "instripEffi" + "." + "txt").c_str());

  gCanvas.push_back(new TCanvas());

  eff.Draw();

  gCanvas.push_back(new TCanvas());

  cl_instrip.Draw();
  SCT_helpers::saveTH1_as_txt(*cl_instrip.getProfile(), (path_ + name_ + "clusterSize" + "." + "txt").c_str());
  gBrowser = new TBrowser();
#ifdef _DEBUG
  theApp.Run();
#endif // _DEBUG


#endif
  return 0;
}
int main(int argc, char **argv) {

  remove_root_printouts();






  inParam para;
  para.argc = argc;
  para.argv = argv;

 
 asyncMain( &para);


  return 0;

}

