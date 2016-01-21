#include "processorBase.h"
#include "TFile.h"

#include "plane_def.h"
#include <string>
#include <iostream>
#include "treeCollection.h"
#include "TF1.h"
#include "TCanvas.h"
#include "geometry/setup_description.hh"
#include "xml_helpers/xml_fileList.hh"
#include "SCT_helpers.h"
#include "s_file_base.h"
#include "sct_processors.h"
#include "internal/hit_efficiency.hh"
#include "internal/inStripClusterSize.hh"
#include "internal/residual_efficienct.hh"
#include "internal/inStripEfficiency.hh"
#include "sct_types.h"
#include "internal/exceptions.hh"
#include "processors/find_nearest_strip.hh"
#define  Class_factory_Utilities_THROW(x) SCT_THROW(x)

std::string to_string(const std::string& par) {
  return par;
}
#include "factoryDef.hh"

bool gDo_print = false;
int gPos = 0;
using namespace sct_type;
std::ostream* m_out = &std::cout;

void drawResidual(s_process_collection_standard& p, const xmlImputFiles::MinMaxRange<double> * range_ = nullptr) {
#ifdef _DEBUG
  new TCanvas();
#endif // _DEBUG

  if (range_) {
    p.DrawResidual(range_->getMin(), range_->getMax());

  } else {
    p.DrawResidual();

  }

}

void Draw_Track_hits(s_process_collection_standard& p) {
#ifdef _DEBUG
  new TCanvas();
#endif // _DEBUG

  p.Draw_Hit_map();
}
void Draw_DUT_hits(s_process_collection_standard& p) {
#ifdef _DEBUG
  new TCanvas();
#endif // _DEBUG

  p.Draw_DUT_Hits_map();
}
void draw_efficiency_map(s_process_collection_standard& p) {
#ifdef _DEBUG
  new TCanvas();
#endif // _DEBUG

  p.Draw_Efficinecy_map();
}

void Draw_missing_coordinate(s_process_collection_standard& p, const xmlImputFiles::MinMaxRange<double> * range_ = nullptr) {
#ifdef _DEBUG
  new TCanvas();
#endif // _DEBUG

  if (range_) {
    p.DrawResidualVsMissingCordinate(range_->getMin(), range_->getMax());
  } else {
    p.DrawResidualVsMissingCordinate();
  }
}
void Draw_Residual_VS_N(s_process_collection_standard& p, const xmlImputFiles::MinMaxRange<double> * range_ = nullptr) {
#ifdef _DEBUG
  new TCanvas();
#endif // _DEBUG

  if (range_) {
    p.DrawResidualVsEvent(range_->getMin(), range_->getMax());
  } else {
    p.DrawResidualVsEvent();
  }
}

class xml_print {

public:
  static void setOutPutStream(std::ostream& out) {
    m_out = &out;
  }
  xml_print(const char* name) :m_name(name), m_pos(gPos) {
    gPos += 2;
    if (!gDo_print) {
      return;
    }
    *m_out << std::string(m_pos, ' ') << "<" << m_name << ">" << std::endl;
  }

  template<typename T>
  xml_print(const char* name, const T& value) :xml_print(name) {
    print(value);
    close();
  }
  template<typename T>
  void print(const T& value) {
    if (!gDo_print) {
      return;
    }
    if (first_element) {
      auto s = std::string(m_pos + 2, ' ');
      *m_out << s << value << std::endl;
      first_element = false;
    } else {

      *m_out << std::string(m_pos + 2, ' ') << ", " << value << std::endl;
    }


  }
  virtual ~xml_print() {
    if (!end_printed) {
      close();
    }


  }
  void close() {
    gPos -= 2;
    if (!gDo_print) {
      return;
    }

    *m_out << std::string(m_pos, ' ') << "</" << m_name << ">" << std::endl;
    end_printed = true;
  }
private:

  std::string m_name;
  bool end_printed = false, first_element = true;
  int m_pos;
};

double BinNomialSigma(double totalHits, double DUTHits) {
  return sqrt((DUTHits / totalHits)*(1 - (DUTHits / totalHits))*(1 / totalHits));
}

void pushChannel(sct_corr::rootEventRunOutput& outEvent, Double_t channel_x,
                 Double_t channel_y,
                 Double_t Effi,
                 Double_t NumberOfEvents,
                 Double_t Error_Effi,
                 const sct_type::ID_t& ID
                 ) {
  outEvent.getData(x_axis_def)->push_back(channel_x);
  outEvent.getData(y_axis_def)->push_back(channel_y);
  outEvent.getData(Occupancy_axis_def)->push_back(Effi);
  outEvent.getData(Occupancy_error_axis_def)->push_back(Error_Effi);
  outEvent.getData(NumOfEvents_axis_def)->push_back(NumberOfEvents);
  outEvent.getData(getIDString())->push_back(necessary_CONVERSION(ID));


}


void push2outputEvent(sct_corr::rootEventRunOutput& outEvent, const TH1D& ouantity, const TH1D& numOfEvents, const sct_type::ID_t& ID) {
  for (Int_t i = 0; i < ouantity.GetNbinsX(); ++i) {
    pushChannel(outEvent,
                ouantity.GetBinCenter(i),        //xPosition
                1,                               //yPosition   
                ouantity.GetBinContent(i),       //Efficiency
                numOfEvents.GetBinContent(i),    //Total True Hits

                BinNomialSigma(
                numOfEvents.GetBinContent(i),
                ouantity.GetBinContent(i)
                ),
                ID
                );
  }

}



namespace sct_corr {





Processor::~Processor() {

}

void Processor::push_processor_batch(std::weak_ptr<processorBatch> batch) {
  m_batch = batch;
}



std::weak_ptr<processorBatch> Processor::get_batch() {
  return m_batch;
}


















processorBase::processorBase() : m_outputl(sct_type::collectionName_t("out")) {





}


processorBase::~processorBase() {


}

void processorBase::setOutputName(const char* name) {
  m_outname = name;
}

void processorBase::push_files(TFile* _file, double Threshold, double runNumber) {
  FileProberties p;
  p.setTFile(_file);
  p.m_Threshold = Threshold;
  p.m_runNumber = runNumber;
  m_files.push_back(p);
}

void processorBase::push_files(const char* _fileName, double Threshold, double runNumber, double HV) {
  FileProberties p;
  p.setTFile(std::shared_ptr<TFile>(new TFile(_fileName)));
  if (!p.getTfile()->IsOpen()) {

    return;
  }

  p.m_Threshold = Threshold;
  p.m_runNumber = runNumber;
  p.m_HV = HV;
  m_files.push_back(p);

}



int processorBase::Add_XML_RunList(const std::string& xmlInputFileName, std::string path__, std::string outputPath /*= "."*/, int element /* = -1*/) {
  path__ += "/";
  m_input_files_xml = std::make_shared<xmlImputFiles::XML_imput_file>(xmlInputFileName.c_str());

  if (m_input_files_xml->fileList().empty()) {
    return -1;
  }
  auto collname = m_input_files_xml->globalConfig().CollectionName();
  if (element != -1) {
    outputPath += "/" + collname + "_" + get_suffix() + "_" + std::to_string(element) + ".root";
  } else {
    outputPath += "/" + collname + "_" + get_suffix() + ".root";
  }


  setOutputName(outputPath.c_str());


  setGearFile(m_input_files_xml->globalConfig().gearFile().c_str());


  if (element != -1) {
    if (element >= (int)m_input_files_xml->fileList().size()) {
      SCT_THROW("out of boundary. Selected element number larger then file list size");
    }
    auto& e = m_input_files_xml->fileList()[element];
    push_files((path__ + std::string(e.name())).c_str(), e.threshold(), e.runNumber(), e.HV());

  } else {
    for (auto& e : m_input_files_xml->fileList()) {
      push_files((path__ + std::string(e.name())).c_str(), e.threshold(), e.runNumber(), e.HV());
    }

  }
  if (m_files.empty()) {
    SCT_THROW("not input file found");
  }

  return 0;
}



void processorBase::setGearFile(const char* name) {

  rapidxml::file<> m_file(name);
  rapidxml::xml_document<> m_doc;
  m_doc.parse<0>(m_file.data());

  m_gear = std::make_shared<sct_corr::Xgear>(m_doc.first_node("gear"));
}

void processorBase::setPrintout(bool print) {
  gDo_print = print;
}


const xmlImputFiles::XML_imput_file* processorBase::get_xml_input() const {
  return m_input_files_xml.get();
}

const sct_corr::Xgear* processorBase::get_gear() const {
  return m_gear.get();
}

void processorBase::process_set_run_prob(const FileProberties& fileP) {
  xml_print("fileName", fileP.getTfile()->GetName());
  m_outputl.reset();

  xml_print("m_runNumber", fileP.m_runNumber);
  m_outputl.set_RunNumber(fileP.m_runNumber);


  xml_print("threshold", fileP.m_Threshold);
  m_outputl.set_Threshold(fileP.m_Threshold);

  xml_print("HV", fileP.m_HV);
  m_outputl.set_HV(fileP.m_HV);
}

void processorBase::start_collection(output_TFile_ptr file__) {
  m_outputTree = std::make_shared<sct_corr::treeCollection_ouput>(
    m_outputl,
    &m_buffer,
    true
    );

  m_outputTree->getTTree()->SetDirectory(necessary_CONVERSION(file__)->GetDirectory("/"));
}

bool processorBase::process() {

  TCanvas c;

  auto files = xml_print("files1");


  auto _file1 = sct_corr::output_TFile_ptr(new TFile(
    m_outname.c_str(),
    "recreate"
    ));
  start_collection(_file1);



  for (auto &e : m_files) {

    process_set_run_prob(e);
    process_file(&e);
    m_outputTree->fill();
  }
  end_collection();

  necessary_CONVERSION(_file1)->Write();
  //delete necessary_CONVERSION(_file1);
  return true;
}


Processor::ProcessState processorBatch::ProceessEvent() {
  for (auto& e : m_processors) {
    auto ret = e->ProceessEvent();
    if (ret != ok) {
      return ret;
    }
  }
  return ok;
}

void processorBatch::push_processor(Processor* processor_) {
  if (processor_) {
    m_processors.push_back(processor_);
    processor_->push_processor_batch(get_batch());
  }
}

void processorBatch::loop() {
  while (ProceessEvent() == ok) {

  }
}

std::shared_ptr<processorBatch> create_batch() {
  auto ret = std::make_shared<processorBatch>();

  ret->push_processor_batch(ret);

  return ret;
}

ProcessorXML_loader::ProcessorXML_loader(const std::string& xmlInputFileName, std::string path__, std::string outputPath /*= "."*/) :m_path(path__ + "/") {
  m_input_files_xml = std::make_shared<xmlImputFiles::XML_imput_file>(xmlInputFileName.c_str());


  auto collname = m_input_files_xml->globalConfig().CollectionName();


  m_outname = outputPath += "/" + collname + ".root";






  rapidxml::file<> m_file(m_input_files_xml->globalConfig().gearFile().c_str());
  rapidxml::xml_document<> m_doc;
  m_doc.parse<0>(m_file.data());

  m_gear = std::make_shared<sct_corr::Xgear>(m_doc.first_node("gear"));


  m_files = std::make_shared<FileProberties>();

  m_owndBatch = create_batch();
  m_files->m_batch = m_owndBatch;
  push_processor_batch(m_owndBatch);
  m_owndBatch->push_processor(this);
}

Processor::ProcessState ProcessorXML_loader::ProceessEvent() {

  while (m_iterrator <= m_input_files_xml->fileList().size()) {
    auto& e = m_input_files_xml->fileList()[m_iterrator++];
    m_files->m_Threshold = e.threshold();
    m_files->m_HV = e.HV();
    m_files->m_runNumber = e.runNumber();
    m_files->setTFile(std::shared_ptr<TFile>(new TFile((m_path + std::string(e.name())).c_str())));
    if (m_files->getTfile()->IsOpen()) {
      return ok;
    }
  }

  return done;
}

const FileProberties* ProcessorXML_loader::getData() const {
  return m_files.get();
}

processorEfficiency::processorEfficiency(FileProberties* fileProb) :m_file(fileProb) {

}

Processor::ProcessState processorEfficiency::ProceessEvent() {
  process_reset();
  auto file_PRINTOUT = xml_print("file");

  process_set_run_prob();


  m_file_fitter.reset();





  //   m_plotCollection = sct_corr::create_plot_collection();
  //   m_plotCollection->addFile(m_file->getTfile());
  //   m_plotCollection->setOutputFile(m_dummy);
  // 
  //   m_file_fitter = std::make_shared<sct_files::fitter_file>(m_plotCollection, get_gear());

  //   m_output_planes = m_file_fitter->get_correlations_channel(
  //     get_xml_input()->globalConfig().cut(),
  //     residualCut_t(get_xml_input()->globalConfig().residual_cut()),
  //     rot_angle_t(get_xml_input()->globalConfig().Rotation()),
  //     move_t(get_xml_input()->globalConfig().Position_value()),
  //     s_plot_prob("Hits_in_channels")
  //     .SaveToDisk()
  //     );



  auto res = sct_corr::processor::residual(
    m_file_fitter->DUT_fitted_local_GBL().getX_def(),
    m_file_fitter->DUT_hit_local().getX_def(),
    s_plot_prob("residualVSEvent").SaveToDisk()
    );

  m_res_VS_event.push_back(res);

#ifdef _DEBUG
  m_plotCollection->loop(40000);
#else
  m_plotCollection->loop();
#endif // _DEBUG



  //   Draw_Efficinecy_map();
  // 
  //   extract_hitMap();
  //   extract_efficiency();
  //   extract_residual();
  //   extract_rotation();

  m_outputTree->fill();
  return ok;
}

void processorEfficiency::process_reset() {
  m_plotCollection.reset();
  m_res_VS_event.clear();
  //  m_outputl.reset();
}

void processorEfficiency::process_set_run_prob() {
  xml_print("fileName", m_file->getTfile()->GetName());


  xml_print("m_runNumber", m_file->m_runNumber);
  // m_outputl.set_RunNumber(m_file->m_runNumber);


  xml_print("threshold", m_file->m_Threshold);
  //  m_outputl.set_Threshold(m_file->m_Threshold);

  xml_print("HV", m_file->m_HV);
  //  m_outputl.set_HV(m_file->m_HV);
}

}
void display_registered_processors(std::ostream& out) {
  auto ret = Class_factory_Utilities::Factory<sct_corr::processorBase>::GetTypes();
  out << "registered Processors: \n";
  for (auto e : ret) {
    out << e << std::endl;

  }

}
std::unique_ptr<sct_corr::processorBase> create_processor(const std::string& processorName) {
  std::string dummy;
  return Class_factory_Utilities::Factory<sct_corr::processorBase>::Create(processorName, dummy);

}

void s_process_collection_standard::extract_efficiency() {

  double totalHits = (double)m_plotCollection->Draw(m_output_planes.getTotalTrueHits(),
                                                    S_DrawOption()
                                                    .cut_x(get_xml_input()->globalConfig().AvtiveStrips().getMin(),
                                                    get_xml_input()->globalConfig().AvtiveStrips().getMax())
                                                    );

  xml_print("TotalNumOfEvents", totalHits);
  m_outputl.set_TotalNumOfEvents(totalHits);

  double DUTHits = (double)m_plotCollection->Draw(m_output_planes.getTrueHitsWithDUT(),
                                                  S_DrawOption()
                                                  .cut_x(get_xml_input()->globalConfig().AvtiveStrips().getMin(),
                                                  get_xml_input()->globalConfig().AvtiveStrips().getMax())
                                                  );
  xml_print("DUTHits", DUTHits);

  xml_print("Efficiency", DUTHits / totalHits);
  m_outputl.set_Total_efficiency(DUTHits / totalHits);


  auto Error_efficiency = BinNomialSigma(totalHits,
                                         DUTHits);

  xml_print("Error_efficiency", Error_efficiency);
  m_outputl.set_Error_efficiency(Error_efficiency);


}

void s_process_collection_standard::extract_hitMap() {
  push2outputEvent(m_outputl, *m_Efficieny_map, *m_Efficieny_trueHits, sct_type::ID_t(0));
}

void s_process_collection_standard::extract_residual() {
  DrawResidual(-3, 3);


  TF1 f("f1", "gaus");

  m_Residual->Fit(&f, "Q");
  {
    auto residual_sigma = f.GetParameter("Sigma");
    xml_print("residual_sigma", residual_sigma);
    m_outputl.set_residual(residual_sigma);
  }
  {
    auto residual_mean = f.GetParameter("Mean");
    xml_print("residual_mean", residual_mean);
    m_outputl.set_offset(residual_mean);
  }
}

void s_process_collection_standard::extract_rotation() {
  DrawResidualVsMissingCordinate(-10, 10);
  auto h = getResidualVsMissingCordinate();
  auto f1 = SCT_helpers::LinearFit_Of_Profile(h, sct_type::procent_t(1));
  auto rot = TMath::ATan(f1.GetParameter("p1"));
  m_outputl.set_rotation(rot);
  xml_print("rotation", rot);
}

void s_process_collection_standard::process_reset() {
  m_plotCollection.reset();
  m_res_VS_event.clear();
}


bool s_process_collection_standard::process_file(FileProberties* fileP) {
  process_reset();
  auto file_PRINTOUT = xml_print("file");



  m_file_fitter.reset();





  m_plotCollection = sct_corr::create_plot_collection();
  m_plotCollection->addFile(fileP->getTfile());
  m_plotCollection->setOutputFile(m_dummy);

  m_file_fitter = std::make_shared<sct_files::fitter_file>(m_plotCollection, get_gear());

  m_output_planes = m_file_fitter->get_correlations_channel(
    *get_xml_input(),
    s_plot_prob("Hits_in_channels")
    .SaveToDisk()
    );



  auto res = sct_corr::processor::residual(
    m_file_fitter->DUT_fitted_local_GBL().getX_def(),
    m_file_fitter->DUT_hit_local().getX_def(),
    s_plot_prob("residualVSEvent").SaveToDisk()
    );

  m_res_VS_event.push_back(res);

#ifdef _DEBUG
  m_plotCollection->loop(20000);
#else
  m_plotCollection->loop();
#endif // _DEBUG



  Draw_Efficinecy_map();

  extract_hitMap();
  extract_efficiency();
  extract_residual();
  extract_rotation();

  m_outputTree->fill();
  return true;
}



s_process_collection_standard::s_process_collection_standard(Parameter_ref) {
  m_dummy = new TFile("dummy1.root", "recreate");
  if (!m_dummy->IsOpen()) {
    SCT_THROW("unable to open file: dummy1.root");
  }
}

s_process_collection_standard::~s_process_collection_standard() {
  if (m_dummy) {
    m_dummy->Write();
  }

}

Long64_t s_process_collection_standard::DrawResidual(Double_t min_X, Double_t max_X) {
  m_Residual = std::make_shared<TH1D>(
    "residual",
    "residual",
    100,
    min_X,
    max_X
    );


  return m_plotCollection->Draw(
    m_output_planes.getResidual(),
    S_DrawOption()
    .draw_x()
    .cut_x(min_X, max_X)
    .output_object(m_Residual.get())
    );
}



Long64_t s_process_collection_standard::DrawResidual() {

  auto ret = m_plotCollection->Draw(
    m_output_planes.getResidual(),
    S_DrawOption().draw_x()
    );

  TH1* h = dynamic_cast<TH1*>(gPad->GetPrimitive("htemp"));
  h->SetTitle("residual");
  return ret;
}

Long64_t s_process_collection_standard::DrawResidualVsEvent(Double_t min_X, Double_t max_X) {
  m_ResidualVsEvent = std::make_shared<TH2D>(
    "ResidualVsEvent",
    "Residual Vs Event",
    100, 0, 0,
    100, min_X, max_X
    );

  return m_plotCollection->Draw(
    m_res_VS_event(),
    S_DrawOption()
    .draw_x_VS_y()
    .cut_x(min_X, max_X)
    .output_object(m_ResidualVsEvent.get())
    .opt_colz()
    );
}

Long64_t s_process_collection_standard::DrawResidualVsEvent() {

  auto ret = m_plotCollection->Draw(
    m_res_VS_event(),
    S_DrawOption().draw_x_VS_y()
    );

  TH2* h = dynamic_cast<TH2*>(gPad->GetPrimitive("htemp"));
  h->SetTitle("Residual Vs Event");

  return ret;

}

Long64_t s_process_collection_standard::DrawResidualVsMissingCordinate(Double_t min_X, Double_t max_X) {
  m_resVSMissing = std::make_shared<TH2D>(
    "ResidualVsMissingCordinate",
    "Residual Vs Missing Coordinate",
    100, 0, 0,
    100, min_X, max_X
    );




  auto ret = m_plotCollection->Draw(
    m_output_planes.getResidualVSmissing(),
    S_DrawOption()
    .draw_x_VS_y()
    .cut_x(min_X, max_X)
    .output_object(m_resVSMissing.get())
    .opt_colz()
    );
  auto f = new TF1(SCT_helpers::LinearFit_Of_Profile(m_resVSMissing.get(), sct_type::procent_t(0)));
  //   std::cout << f->GetParameter("p1") << std::endl;
  //   std::cout << f->GetParameter("p0") << std::endl;
  m_plotCollection->Draw(
    m_output_planes.getResidualVSmissing(),
    S_DrawOption()
    .draw_x_VS_y()
    .cut_x(min_X, max_X)
    .output_object(m_resVSMissing.get())
    .opt_colz()
    );

  f->Draw("same");

  return ret;
}

Long64_t s_process_collection_standard::DrawResidualVsMissingCordinate() {

  auto ret = m_plotCollection->Draw(
    m_output_planes.getResidualVSmissing(),
    S_DrawOption()
    .draw_x_VS_y()
    .opt_colz()
    );


  //TH2* h = dynamic_cast<TH2*>(gPad->GetPrimitive("htemp"));
  //auto f = new TF1(SCT_helpers::LinearFit_Of_Profile(h, 0.2));
  //f->Draw("same");
  // h->SetTitle("Residual Vs Missing Coordinate");

  return ret;
}

Long64_t s_process_collection_standard::Draw_Efficinecy_map() {

  m_Efficieny_trueHits = std::make_shared<TH1D>(
    "total",
    "total",
    get_xml_input()->globalConfig().NumberOfBins(), -0.5, get_xml_input()->globalConfig().NumberOfStrips() - 0.5
    );

  m_plotCollection->Draw(
    m_output_planes.getTotalTrueHits(),
    S_DrawOption()
    .draw_x()
    .output_object(m_Efficieny_trueHits.get())
    );

  m_Efficieny_map = std::make_shared<TH1D>(
    "Efficiency",
    "Efficiency",
    get_xml_input()->globalConfig().NumberOfBins(), -0.5, get_xml_input()->globalConfig().NumberOfStrips() - 0.5
    );

  auto n = m_plotCollection->Draw(
    m_output_planes.getTrueHitsWithDUT(),
    S_DrawOption()
    .draw_x()
    .output_object(m_Efficieny_map.get())
    );
  auto e = SCT_helpers::calc_efficiency(m_Efficieny_trueHits.get(), m_Efficieny_map.get());
  auto eth2d = dynamic_cast<TH1D*>(e);

  m_Efficieny_map = std::shared_ptr<TH1D>(eth2d);
  // m_Efficieny_map->Divide(m_Efficieny_trueHits.get());

  m_Efficieny_map->Draw();
  return n;
}

Long64_t s_process_collection_standard::Draw_Hit_map() {
  m_Hits_total = std::make_shared<TH1D>(
    "total",
    "total",
    get_xml_input()->globalConfig().NumberOfBins(), -0.5, get_xml_input()->globalConfig().NumberOfStrips() - 0.5
    );

  return  m_plotCollection->Draw(
    m_output_planes.getTotalTrueHits(),
    S_DrawOption()
    .draw_x()
    .output_object(m_Hits_total.get())
    );
}

Long64_t s_process_collection_standard::Draw_DUT_Hits_map() {
  m_Hits_with_DUT_Hits = std::make_shared<TH1D>(
    "DUT",
    "DUT",
    get_xml_input()->globalConfig().NumberOfBins(), -0.5, get_xml_input()->globalConfig().NumberOfStrips() - 0.5
    );

  return m_plotCollection->Draw(
    m_output_planes.getTrueHitsWithDUT(),
    S_DrawOption()
    .draw_x()
    .output_object(m_Hits_with_DUT_Hits.get())
    );
}

TH2D* s_process_collection_standard::getResidualVsMissingCordinate() {
  return m_resVSMissing.get();
}

void s_process_collection_standard::saveHistograms(TFile* outPutFile /*= nullptr*/, xmlImputFiles::MinMaxRange<double>* residual_cut /*= nullptr*/) {
  drawResidual(*this, residual_cut);


  Draw_DUT_hits(*this);

  Draw_Track_hits(*this);

  draw_efficiency_map(*this);
  Draw_Residual_VS_N(*this, residual_cut);
  Draw_missing_coordinate(*this, residual_cut);
  if (outPutFile) {
    outPutFile->Add(m_Efficieny_map.get());
    outPutFile->Add(m_Efficieny_trueHits.get());
    outPutFile->Add(m_Hits_total.get());
    outPutFile->Add(m_Hits_with_DUT_Hits.get());
    outPutFile->Add(m_ResidualVsEvent.get());
    outPutFile->Add(m_resVSMissing.get());
    outPutFile->Add(m_Residual.get());
  }

}

std::string s_process_collection_standard::get_suffix() const {
  return "standard";
}

TFile* FileProberties::getTfile() const {
  if (m_fileOwnd) {
    return m_fileOwnd.get();
  }

  if (m_file) {
    return m_file;
  }
  return nullptr;
}

void FileProberties::setTFile(std::shared_ptr<TFile> file) {
  m_fileOwnd = file;
}

void FileProberties::setTFile(TFile* file) {
  m_file = file;
}

s_process_collection_modulo::s_process_collection_modulo(Parameter_ref) {
  m_dummy = new TFile("dummy1.root", "recreate");
}

s_process_collection_modulo::~s_process_collection_modulo() {

}


bool s_process_collection_modulo::process_file(FileProberties* fileP) {


  m_plotCollection = sct_corr::create_plot_collection();
  m_plotCollection->addFile(fileP->getTfile());
  m_plotCollection->setOutputFile(m_dummy);
  m_file_fitter = std::make_shared<sct_files::fitter_file>(m_plotCollection, get_gear());

  auto pl = m_file_fitter->get_collection();



  m_gbl_collection = m_file_fitter->get_correlations_channel(
    *get_xml_input(),
    s_plot_prob("Hits_in_channels")
    .SaveToDisk()
    );
  auto ActiveStrips = get_xml_input()->globalConfig().AvtiveStrips();



  m_instripEfficiency = std::make_shared<sct_corr::inStripEfficiency>(
    m_gbl_collection.getTotalTrueHits(),
    m_gbl_collection.getTrueHitsWithDUT(),
    S_XCut(ActiveStrips.getMin(), ActiveStrips.getMax()),
    x_axis_def,
    sct_type::modulo_t(3),
    s_plot_prob("inStripEffi")
    );

  m_instripClusterSize = std::make_shared<sct_corr::inStripClusterSize>(
    m_gbl_collection.getTrueHitsWithDUT(),
    m_file_fitter->DUT_zs_data(),
    10,
    x_axis_def,
    sct_type::modulo_t(3),
    s_plot_prob("cluster_size_instrip").SaveToDisk()
    );


  // auto second_hit = sct_corr::processor::remove_closest(m_file_fitter->DUT_zs_data(), m_gbl_collection.getTotalTrueHits(), x_axis_def, s_plot_prob().doNotSaveToDisk());

  m_residualEffieciency = std::make_shared<
    sct_corr::residual_efficienct>(
    m_gbl_collection.getTotalTrueHits(),
    m_file_fitter->DUT_zs_data(),
    S_XCut(ActiveStrips.getMin(), ActiveStrips.getMax()),
    sct_type::stripNr_t(ActiveStrips.getMax() + 20),
    x_axis_def,
    s_plot_prob("Res_efficiency")
    );

#ifdef _DEBUG
  pl->loop(20000);
#else
  pl->loop();
#endif // _DEBUG


  m_residualEffieciency->Draw();
  m_instripEfficiency->Draw();
  push2outputEvent(m_outputl, *m_residualEffieciency->getEfficiency_map(), *m_residualEffieciency->get_total(), ID_t(0));
  push2outputEvent(m_outputl, *m_instripEfficiency->getEfficiency_map(), *m_instripEfficiency->getHits(), ID_t(1));
  return true;
}

void s_process_collection_modulo::saveHistograms(TFile* outPutFile /*= nullptr */, xmlImputFiles::MinMaxRange<double>* residual_cut /*= nullptr */) {

#ifdef _DEBUG
  new TCanvas();
#endif
  m_instripClusterSize->Draw();
#ifdef _DEBUG
  new TCanvas();
#endif // _DEBUG
  m_instripEfficiency->Draw();
#ifdef _DEBUG
  new TCanvas();
#endif // _DEBUG
  m_residualEffieciency->Draw();

  if (outPutFile) {
    outPutFile->Add(m_instripClusterSize->getHistogram());
    outPutFile->Add(m_instripEfficiency->getEfficiency_map());
    outPutFile->Add(m_residualEffieciency->getEfficiency_map());
    outPutFile->Add(m_residualEffieciency->get_total());
  }
}

std::string s_process_collection_modulo::get_suffix() const {
  return "modulo";
}

s_process_collection_modulo_second::s_process_collection_modulo_second(Parameter_ref) {
  m_dummy = new TFile("dummy1.root", "recreate");
  if (!m_dummy->IsOpen()) {
    SCT_THROW("unable to open file: dummy1.root");
  }
}

s_process_collection_modulo_second::~s_process_collection_modulo_second() {
  if (m_dummy) {
    m_dummy->Write();
  }
}

std::string s_process_collection_modulo_second::get_suffix() const {
  return "modulo_second";
}

bool s_process_collection_modulo_second::process_file(FileProberties* fileP) {
  m_plotCollection = sct_corr::create_plot_collection();
  m_plotCollection->addFile(fileP->getTfile());
  m_plotCollection->setOutputFile(m_dummy);
  m_file_fitter = std::make_shared<sct_files::fitter_file>(m_plotCollection, get_gear());

  auto pl = m_file_fitter->get_collection();



  m_gbl_collection = m_file_fitter->get_correlations_channel(
    *get_xml_input(),
    s_plot_prob("Hits_in_channels")
    .SaveToDisk()
    );
  auto ActiveStrips = get_xml_input()->globalConfig().AvtiveStrips();

  auto second_hit1 = sct_corr::processor::remove_closest(m_file_fitter->DUT_zs_data(), m_gbl_collection.getTotalTrueHits(), x_axis_def, s_plot_prob().doNotSaveToDisk());


  auto find_closest = sct_processor::find_nearest_strip(
    m_gbl_collection.getTotalTrueHits(),
    second_hit1,
    x_axis_def,
    get_xml_input()->globalConfig().residual_cut(),
    s_plot_prob("second_highest")
    );


  m_instripEfficiency = std::make_shared<sct_corr::inStripEfficiency>(
    m_gbl_collection.getTotalTrueHits(),
    find_closest.getHitOnPlaneA(),
    S_XCut(ActiveStrips.getMin(), ActiveStrips.getMax()),
    x_axis_def,
    sct_type::modulo_t(3),
    s_plot_prob("inStripEffi")
    );






#ifdef _DEBUG
  pl->loop(20000);
#else
  pl->loop();
#endif // _DEBUG


  m_instripEfficiency->Draw();

  push2outputEvent(m_outputl, *m_instripEfficiency->getEfficiency_map(), *m_instripEfficiency->getHits(), ID_t(0));
  return true;
}

void s_process_collection_modulo_second::saveHistograms(TFile* outPutFile /*= nullptr*/, xmlImputFiles::MinMaxRange<double>* residual_cut /*= nullptr*/) {
#ifdef _DEBUG
  new TCanvas();
#endif // _DEBUG
  m_instripEfficiency->Draw();

  if (outPutFile) {
    outPutFile->Add(m_instripEfficiency->getEfficiency_map());
    outPutFile->Add(m_instripEfficiency->getHits());
  }
}

void s_process_collection_standard_second::saveHistograms(TFile* outPutFile /*= nullptr*/, xmlImputFiles::MinMaxRange<double>* residual_cut /*= nullptr*/) {
  DrawResidual(residual_cut->getMin(), residual_cut->getMax());

  Draw_DUT_Hits_map();

  Draw_Hit_map();

  Draw_Efficinecy_map();


  if (outPutFile) {
    outPutFile->Add(m_Efficieny_map.get());
    outPutFile->Add(m_Efficieny_trueHits.get());
    outPutFile->Add(m_Hits_total.get());
    outPutFile->Add(m_Hits_with_DUT_Hits.get());
    outPutFile->Add(m_Residual.get());
  }

}

s_process_collection_standard_second::s_process_collection_standard_second(Parameter_ref) {
  m_dummy = new TFile("dummy1.root", "recreate");
  if (!m_dummy->IsOpen()) {
    SCT_THROW("unable to open file: dummy1.root");
  }
}

s_process_collection_standard_second::~s_process_collection_standard_second() {
  if (m_dummy) {
    m_dummy->Write();
  }
}

bool s_process_collection_standard_second::process_file(FileProberties* fileP) {
  process_reset();
  auto file_PRINTOUT = xml_print("file");



  m_file_fitter.reset();





  m_plotCollection = sct_corr::create_plot_collection();
  m_plotCollection->addFile(fileP->getTfile());
  m_plotCollection->setOutputFile(m_dummy);

  m_file_fitter = std::make_shared<sct_files::fitter_file>(m_plotCollection, get_gear());

  m_gbl_collection = m_file_fitter->get_correlations_channel(
    *get_xml_input(),
    s_plot_prob("Hits_in_channels")
    .SaveToDisk()
    );

  auto second_hit1 = sct_corr::processor::remove_closest(m_file_fitter->DUT_zs_data(), m_gbl_collection.getTotalTrueHits(), x_axis_def, s_plot_prob().doNotSaveToDisk());


  auto find_closest = sct_processor::find_nearest_strip(
    m_gbl_collection.getTotalTrueHits(),
    second_hit1,
    x_axis_def,
    get_xml_input()->globalConfig().residual_cut(),
    s_plot_prob("second_highest")
    );


  m_collection.push_back("TrueHits_with_DUT", find_closest.getHitOnPlaneA());
  m_collection.push_back("residual", find_closest.getResidual());



#ifdef _DEBUG
  m_plotCollection->loop(20000);
#else
  m_plotCollection->loop();
#endif // _DEBUG



  Draw_Efficinecy_map();

  extract_hitMap();
  extract_efficiency();
  extract_residual();


  m_outputTree->fill();
  return true;
}

Long64_t s_process_collection_standard_second::DrawResidual(Double_t min_X, Double_t max_X) {
  m_Residual = std::make_shared<TH1D>(
    "residual",
    "residual",
    100,
    min_X,
    max_X
    );


  return m_plotCollection->Draw(
    m_collection.get("residual")(),
    S_DrawOption()
    .draw_x()
    .cut_x(min_X, max_X)
    .output_object(m_Residual.get())
    );
}

Long64_t s_process_collection_standard_second::Draw_Efficinecy_map() {

  m_Efficieny_trueHits = std::make_shared<TH1D>(
    "total",
    "total",
    get_xml_input()->globalConfig().NumberOfBins(), -0.5, get_xml_input()->globalConfig().NumberOfStrips() - 0.5
    );

  m_plotCollection->Draw(
    m_gbl_collection.getTotalTrueHits(),
    S_DrawOption()
    .draw_x()
    .output_object(m_Efficieny_trueHits.get())
    );

  m_Efficieny_map = std::make_shared<TH1D>(
    "Efficiency",
    "Efficiency",
    get_xml_input()->globalConfig().NumberOfBins(), -0.5, get_xml_input()->globalConfig().NumberOfStrips() - 0.5
    );

  auto n = m_plotCollection->Draw(
    m_collection.get("TrueHits_with_DUT")(),
    S_DrawOption()
    .draw_x()
    .output_object(m_Efficieny_map.get())
    );
  auto e = SCT_helpers::calc_efficiency(m_Efficieny_trueHits.get(), m_Efficieny_map.get());
  auto eth2d = dynamic_cast<TH1D*>(e);

  m_Efficieny_map = std::shared_ptr<TH1D>(eth2d);
  // m_Efficieny_map->Divide(m_Efficieny_trueHits.get());

  m_Efficieny_map->Draw();
  return n;
}

std::string s_process_collection_standard_second::get_suffix() const {
  return "standard_second";
}

void s_process_collection_standard_second::process_reset() {
  m_plotCollection.reset();
  m_collection.clear();
}

Long64_t s_process_collection_standard_second::Draw_Hit_map() {
  m_Hits_total = std::make_shared<TH1D>(
    "total",
    "total",
    get_xml_input()->globalConfig().NumberOfBins(), -0.5, get_xml_input()->globalConfig().NumberOfStrips() - 0.5
    );

  return  m_plotCollection->Draw(
    m_gbl_collection.getTotalTrueHits(),
    S_DrawOption()
    .draw_x()
    .output_object(m_Hits_total.get())
    );
}

Long64_t s_process_collection_standard_second::Draw_DUT_Hits_map() {
  m_Hits_with_DUT_Hits = std::make_shared<TH1D>(
    "DUT",
    "DUT",
    get_xml_input()->globalConfig().NumberOfBins(), -0.5, get_xml_input()->globalConfig().NumberOfStrips() - 0.5
    );

  return m_plotCollection->Draw(
    m_collection.get("TrueHits_with_DUT")(),
    S_DrawOption()
    .draw_x()
    .output_object(m_Hits_with_DUT_Hits.get())
    );
}

void s_process_collection_standard_second::extract_efficiency() {
  double totalHits = (double)m_plotCollection->Draw(m_gbl_collection.getTotalTrueHits(),
                                                    S_DrawOption()
                                                    .cut_x(get_xml_input()->globalConfig().AvtiveStrips().getMin(),
                                                    get_xml_input()->globalConfig().AvtiveStrips().getMax())
                                                    );

  xml_print("TotalNumOfEvents", totalHits);
  m_outputl.set_TotalNumOfEvents(totalHits);

  double DUTHits = (double)m_plotCollection->Draw(m_collection.get("TrueHits_with_DUT")(),
                                                  S_DrawOption()
                                                  .cut_x(get_xml_input()->globalConfig().AvtiveStrips().getMin(),
                                                  get_xml_input()->globalConfig().AvtiveStrips().getMax())
                                                  );
  xml_print("DUTHits", DUTHits);

  xml_print("Efficiency", DUTHits / totalHits);
  m_outputl.set_Total_efficiency(DUTHits / totalHits);


  auto Error_efficiency = BinNomialSigma(totalHits,
                                         DUTHits);

  xml_print("Error_efficiency", Error_efficiency);
  m_outputl.set_Error_efficiency(Error_efficiency);
}

void s_process_collection_standard_second::extract_hitMap() {
  push2outputEvent(m_outputl, *m_Efficieny_map, *m_Efficieny_trueHits, sct_type::ID_t(0));
}

void s_process_collection_standard_second::extract_residual() {
  DrawResidual(-3, 3);


  TF1 f("f1", "gaus");

  m_Residual->Fit(&f, "Q");
  {
    auto residual_sigma = f.GetParameter("Sigma");
    xml_print("residual_sigma", residual_sigma);
    m_outputl.set_residual(residual_sigma);
  }
  {
    auto residual_mean = f.GetParameter("Mean");
    xml_print("residual_mean", residual_mean);
    m_outputl.set_offset(residual_mean);
  }
}

void s_process_collection_modulo_ex::saveHistograms(TFile* outPutFile /*= nullptr*/, xmlImputFiles::MinMaxRange<double>* residual_cut /*= nullptr*/) {
#ifdef _DEBUG
  new TCanvas();
#endif // _DEBUG
  m_instripEfficiency->Draw();
  if (outPutFile) {
    outPutFile->Add(m_instripEfficiency->getEfficiency_map());
    outPutFile->Add(m_instripEfficiency->getHits());
  }
}

std::string s_process_collection_modulo_ex::get_suffix() const {
  return "modulo_ex";
}

bool s_process_collection_modulo_ex::process_file(FileProberties* fileP) {
  m_plotCollection = sct_corr::create_plot_collection();
  m_plotCollection->addFile(fileP->getTfile());
  m_plotCollection->setOutputFile(m_dummy);
  m_file_fitter = std::make_shared<sct_files::fitter_file>(m_plotCollection, get_gear());

  auto pl = m_file_fitter->get_collection();



  m_gbl_collection = m_file_fitter->get_correlations_channel(
    *get_xml_input(),
    s_plot_prob("Hits_in_channels")
    .SaveToDisk()
    );
  auto ActiveStrips = get_xml_input()->globalConfig().AvtiveStrips();



  m_instripEfficiency = std::make_shared<sct_corr::inStripEfficiency>(
    m_gbl_collection.getTotalTrueHits(),
    m_gbl_collection.getTrueHitsWithDUT(),
    S_XCut(ActiveStrips.getMin(), ActiveStrips.getMax()),
    x_axis_def,
    sct_type::modulo_t(3),
    s_plot_prob("inStripEffi")
    );




#ifdef _DEBUG
  pl->loop(20000);
#else
  pl->loop();
#endif // _DEBUG



  m_instripEfficiency->Draw();
  push2outputEvent(m_outputl, *m_instripEfficiency->getEfficiency_map(), *m_instripEfficiency->getHits(), ID_t(0));
  return true;
}

s_process_collection_modulo_ex::~s_process_collection_modulo_ex() {

}

s_process_collection_modulo_ex::s_process_collection_modulo_ex(Parameter_ref) {
  m_dummy = new TFile("dummy1.root", "recreate");
}

s_process_collection_residual::s_process_collection_residual(Parameter_ref) {
  m_dummy = new TFile("dummy1.root", "recreate");
}

s_process_collection_residual::~s_process_collection_residual() {

}

void s_process_collection_residual::saveHistograms(TFile* outPutFile /*= nullptr*/, xmlImputFiles::MinMaxRange<double>* residual_cut /*= nullptr*/) {

#ifdef _DEBUG
  new TCanvas();
#endif // _DEBUG
  m_residualEffieciency->Draw();

  if (outPutFile) {
    outPutFile->Add(m_residualEffieciency->getEfficiency_map());
    outPutFile->Add(m_residualEffieciency->get_total());
  }
}

std::string s_process_collection_residual::get_suffix() const {
  return "residual";
}

bool s_process_collection_residual::process_file(FileProberties* fileP) {

  m_plotCollection = sct_corr::create_plot_collection();
  m_plotCollection->addFile(fileP->getTfile());
  m_plotCollection->setOutputFile(m_dummy);
  m_file_fitter = std::make_shared<sct_files::fitter_file>(m_plotCollection, get_gear());

  auto pl = m_file_fitter->get_collection();



  m_gbl_collection = m_file_fitter->get_correlations_channel(
    *get_xml_input(),
    s_plot_prob("Hits_in_channels")
    .SaveToDisk()
    );
  auto ActiveStrips = get_xml_input()->globalConfig().AvtiveStrips();


  m_residualEffieciency = std::make_shared<
    sct_corr::residual_efficienct>(
    m_gbl_collection.getTotalTrueHits(),
    m_file_fitter->DUT_zs_data(),
    S_XCut(ActiveStrips.getMin(), ActiveStrips.getMax()),
    sct_type::stripNr_t(ActiveStrips.getMax() + 20),
    x_axis_def,
    s_plot_prob("Res_efficiency")
    );

#ifdef _DEBUG
  pl->loop(20000);
#else
  pl->loop();
#endif // _DEBUG


  m_residualEffieciency->Draw();
  push2outputEvent(m_outputl, *m_residualEffieciency->getEfficiency_map(), *m_residualEffieciency->get_total(), ID_t(0));
  return true;
}

using namespace sct_corr;
registerBaseClassDef(processorBase);

registerProcessor("Modulo", s_process_collection_modulo);
registerProcessor("Modulo_ex", s_process_collection_modulo_ex);
registerProcessor("Standard", s_process_collection_standard);
registerProcessor("Standard_second", s_process_collection_standard_second);
registerProcessor("Modulo_second", s_process_collection_modulo_second);
registerProcessor("Residual", s_process_collection_residual);
registerProcessor("Residual_second", s_process_collection_residual_second);

s_process_collection_residual_second::s_process_collection_residual_second(Parameter_ref par) {

  m_dummy = new TFile("dummy1.root", "recreate");
  if (!m_dummy->IsOpen()) {
    SCT_THROW("unable to open file: dummy1.root");
  }

}

s_process_collection_residual_second::~s_process_collection_residual_second() {
  if (m_dummy) {
    m_dummy->Write();
  }

}

void s_process_collection_residual_second::saveHistograms(TFile* outPutFile /*= nullptr*/, xmlImputFiles::MinMaxRange<double>* residual_cut /*= nullptr*/) {
#ifdef _DEBUG
  new TCanvas();
#endif // _DEBUG
  m_residualEffieciency->Draw();
  if (outPutFile) {
    outPutFile->Add(m_residualEffieciency->getEfficiency_map());
    outPutFile->Add(m_residualEffieciency->get_total());
  }
}

std::string s_process_collection_residual_second::get_suffix() const {
  return "residual_second";
}

bool s_process_collection_residual_second::process_file(FileProberties* fileP) {
  m_plotCollection = sct_corr::create_plot_collection();
  m_plotCollection->addFile(fileP->getTfile());
  m_plotCollection->setOutputFile(m_dummy);
  m_file_fitter = std::make_shared<sct_files::fitter_file>(m_plotCollection, get_gear());

  auto pl = m_file_fitter->get_collection();



  m_gbl_collection = m_file_fitter->get_correlations_channel(
    *get_xml_input(),
    s_plot_prob("Hits_in_channels")
    .SaveToDisk()
    );
  auto ActiveStrips = get_xml_input()->globalConfig().AvtiveStrips();

  auto second_hit1 = sct_corr::processor::remove_closest(m_file_fitter->DUT_zs_data(), m_gbl_collection.getTotalTrueHits(), x_axis_def, s_plot_prob().doNotSaveToDisk());


  auto find_closest = sct_processor::find_nearest_strip(
    m_gbl_collection.getTotalTrueHits(),
    second_hit1,
    x_axis_def,
    get_xml_input()->globalConfig().residual_cut(),
    s_plot_prob("second_highest")
    );


  m_residualEffieciency = std::make_shared<sct_corr::residual_efficienct>(
    m_gbl_collection.getTotalTrueHits(),
    find_closest.getHitOnPlaneA(),
    S_XCut(ActiveStrips.getMin(), ActiveStrips.getMax()),
    sct_type::stripNr_t(ActiveStrips.getMax() + 20),
    x_axis_def,
    s_plot_prob("Res_efficiency")
    );






#ifdef _DEBUG
  pl->loop(20000);
#else
  pl->loop();
#endif // _DEBUG


  m_residualEffieciency->Draw();

  push2outputEvent(m_outputl, *m_residualEffieciency->getEfficiency_map(), *m_residualEffieciency->get_total(), ID_t(0));
  return true;
}
