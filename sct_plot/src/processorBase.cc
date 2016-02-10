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

#define _DEBUG 0

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
  //printf("%s\n", S_DrawOption().cut_x(get_xml_input()->globalConfig().AvtiveStrips().getMin(), get_xml_input()->globalConfig().AvtiveStrips().getMax()).getCut().GetTitle());

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
 
	if (m_dummy){
		delete m_dummy;
	}
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
  m_fit = std::make_shared<TF1>(SCT_helpers::LinearFit_Of_Profile(m_resVSMissing.get(), sct_type::procent_t(0)));
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

  m_fit->Draw("same");

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

////////////////////////////////////////////////////////////////////////////////////////////Riccardo////////////////////////////////////////////////////////////////////////////
///////////////////////////// Centering
s_process_collection_centering::s_process_collection_centering(Parameter_ref par) {
	m_dummy = new TFile("dummy1.root", "recreate");

	TTree *treeout = new TTree("centering", "Tree of centered strip");
	treeout->Branch("strip", &stripref, "strip/I");
	treeout->Branch("x_seedcenter", &x_seedcenter, "x_seedcenter/D");
	treeout->Write("centering",2);

	hcentering = new TH2D("histo_centering", "histo for centering", (int)((xmax - xmin) / Dx), xmin, xmax, (int)((ymax - ymin) / Dy), ymin, ymax);
	stripGauss = new TF1("fitStripGauss", "gaus", xmin, xmax);

	Mask.push_back(785);	// TO DO: read from gear
	Mask.push_back(784);	// TO DO: read from gear
	Mask.push_back(820);	// TO DO: read from gear
	Mask.push_back(821);	// TO DO: read from gear
	Mask.push_back(822);	// TO DO: read from gear
	Mask.push_back(823);	// TO DO: read from gear
	Mask.push_back(824);	// TO DO: read from gear
	Mask.push_back(863);	// TO DO: read from gear
	Mask.push_back(864);	// TO DO: read from gear
	Mask.push_back(865);	// TO DO: read from gear
	Mask.push_back(866);	// TO DO: read from gear

}

s_process_collection_centering::~s_process_collection_centering() {

}

bool s_process_collection_centering::process_file(FileProberties* fileP) {

	m_plotCollection = sct_corr::create_plot_collection();
	m_plotCollection->addFile(fileP->getTfile());
	m_plotCollection->setOutputFile(m_dummy);
	m_file_fitter = std::make_shared<sct_files::fitter_file>(m_plotCollection, get_gear());

	auto pl = m_file_fitter->get_collection();

	printf("- Filename: %s.\n", fileP->getTfile()->GetName());

	// Load of trees and initializations
	TTree *treeout = (TTree*)m_dummy->Get("centering");
	treeout->Branch("strip", &stripref, "strip/I");
	treeout->Branch("x_seedcenter", &x_seedcenter, "x_seedcenter/D");


	TTree *tree = (TTree*)fileP->getTfile()->Get("zsdata_strip");
	TTree *tree2 = (TTree*)fileP->getTfile()->Get("GBL_tracks");
	std::vector<double > *ID = 0;
	std::vector<double > *x = 0;
	std::vector<double > *y = 0;
	Int_t event_nr = 0;

	std::vector<double > *ID2 = 0;
	std::vector<double > *x2 = 0;
	std::vector<double > *y2 = 0;
	Int_t event_nr2 = 0;
	std::vector<double > *phi = 0;
	std::vector<double > *chi2 = 0;
	std::vector<double > *ndf = 0;
	std::vector<double > *theta = 0;

	tree->SetBranchAddress("ID", &ID);
	tree->SetBranchAddress("x", &x);
	tree->SetBranchAddress("y", &y);
	tree->SetBranchAddress("event_nr", &event_nr);

	tree2->SetBranchAddress("ID", &ID2);
	tree2->SetBranchAddress("x", &x2);
	tree2->SetBranchAddress("y", &y2);
	tree2->SetBranchAddress("event_nr", &event_nr2);
	tree2->SetBranchAddress("phi", &phi);
	tree2->SetBranchAddress("chi2", &chi2);
	tree2->SetBranchAddress("ndf", &ndf);
	tree2->SetBranchAddress("theta", &theta);
	int n_entries = tree->GetEntries();
	int n_entries2 = tree2->GetEntries();

	auto ActiveStrips = get_xml_input()->globalConfig().AvtiveStrips();

	int plane;
	int id_inclu;
	double tiltedpos_x, tiltedpos_y;
	int closest_strip;
	double d;

	// Loop searching the offset
	stripref = ActiveStrips.getMin() + 1;
	for (int i = 0; i < n_entries; i++) {
		tree->GetEntry(i);
		plane = 0;
		if (x->size() == 0)
			continue;
		tree2->GetEntry(i);
		if (x2->size() == 0)
			continue;
		while (ID->at(0) != ID2->at(plane))
			plane++;
		tiltedpos_x = (TMath::Cos(alpha) * x2->at(plane)) - (TMath::Sin(alpha) * y2->at(plane));
		tiltedpos_y = (TMath::Sin(alpha) * x2->at(plane)) + (TMath::Cos(alpha) * y2->at(plane));
		for (int id_hit = 0; id_hit < x->size(); id_hit++)
			if (x->at(id_hit) == stripref)
				hcentering->Fill(tiltedpos_x, tiltedpos_y);
	}
	stripGauss->SetRange(hcentering->ProjectionX("p_x")->GetBinCenter(hcentering->ProjectionX("p_x")->GetMaximumBin()) - 2 * p, hcentering->ProjectionX("p_x")->GetBinCenter(hcentering->ProjectionX("p_x")->GetMaximumBin()) + 2 * p);
	hcentering->ProjectionX("p_x")->Fit(stripGauss, "R");
	x_seedcenter = stripGauss->GetParameter(1);//hcentering->ProjectionX("p_x")->GetBinCenter(hcentering->ProjectionX("p_x")->GetMaximumBin());	// TO DO: fit with a constrained gaussian????

#ifdef _DEBUG
	printf("- Active strips: %d - %d.\n", ActiveStrips.getMin(), ActiveStrips.getMax());
	std::cout << "Taking the channel " << stripref << " as a reference" << std::endl;
	std::cout << "Center of first strip found at " << x_seedcenter << std::endl;
	new TCanvas();
	hcentering->Draw();
	hcentering->Write("hcentering", 2);
#endif

#ifdef _DEBUG
	new TCanvas();
	hcentering->Draw();
#endif
	hcentering->Write("hcentering", 2);

	treeout->Fill();
	printf("- Filling tree: \n\t- (strip,center)=%d, %f.\n", stripref, x_seedcenter);
	x_seedcenter = x_seedcenter - ((double)stripref)*p;
	stripref = 0;
	treeout->Fill();
	printf("\t- (strip,center)=%d, %f.\n", stripref, x_seedcenter);
	TTree *newtreeout = treeout->CloneTree();
	newtreeout->Print();
	newtreeout->Write("centering", TObject::kOverwrite);
	delete newtreeout;

	/*
	m_residualEffieciency->Draw();
	m_instripEfficiency->Draw();
	push2outputEvent(m_outputl, *m_residualEffieciency->getEfficiency_map(), *m_residualEffieciency->get_total(), ID_t(0));
	push2outputEvent(m_outputl, *m_instripEfficiency->getEfficiency_map(), *m_instripEfficiency->getHits(), ID_t(1));
	*/
	return true;
}
bool s_process_collection_centering::isMasked(int channel) {
	for (int i = 0;i < Mask.size();i++) {
		if (channel == Mask.at(i))
			return true;
	}
	return false;
}


void s_process_collection_centering::saveHistograms(TFile* outPutFile /*= nullptr */, xmlImputFiles::MinMaxRange<double>* residual_cut /*= nullptr */) {

#ifdef _DEBUG
	new TCanvas();
#endif
	//	m_instripClusterSize->Draw();
#ifdef _DEBUG
	new TCanvas();
#endif // _DEBUG
	//	m_instripEfficiency->Draw();
#ifdef _DEBUG
	new TCanvas();
#endif // _DEBUG
	//	m_residualEffieciency->Draw();

	/*	if (outPutFile) {
	outPutFile->Add(m_instripClusterSize->getHistogram());
	outPutFile->Add(m_instripEfficiency->getEfficiency_map());
	outPutFile->Add(m_residualEffieciency->getEfficiency_map());
	outPutFile->Add(m_residualEffieciency->get_total());
	}*/
}

std::string s_process_collection_centering::get_suffix() const {
	return "Centering";
}
//////////////////////////////////////////////
//////////////////////////// Gain vs. position
s_process_collection_modulo_forPositions::s_process_collection_modulo_forPositions(Parameter_ref par) {
	// Retriving alignment
	TFile* m_centering = new TFile("endcap_450V_centering.root", "read");	// TO DO: take name of collection
	TTree *t_centering = new TTree();
	t_centering = (TTree*)m_centering->Get("centering");
	int stripref = 0;
	double x_seedcenter = 0;
	centering Centering;
	t_centering->SetBranchAddress("strip", &stripref);
	t_centering->SetBranchAddress("x_seedcenter", &x_seedcenter);
	t_centering->GetEntry(t_centering->GetEntries()-1);
	Centering.stripref = stripref;
	Centering.x_seedcenter = x_seedcenter;
	Dxcentering = Centering.x_seedcenter - Centering.stripref*p;
	printf("- Centering: strip %d at %f, traslation of %f.\n", Centering.stripref, Centering.x_seedcenter, Dxcentering);
	// Retriving max corr events
	TFile* m_correvents = new TFile("endcap_450V_Correlation_check.root", "read");	// TO DO: take name of collection
	TTree *t_correlated_events = new TTree();
	t_correlated_events = (TTree*)m_correvents->Get("correlated_events");
	int run = 0;
	int max_correvent = 0;
	corr x;
	t_correlated_events->SetBranchAddress("run", &run);
	t_correlated_events->SetBranchAddress("max_correvent", &max_correvent);
	for (int i = 0;i < t_correlated_events->GetEntries();i++) {
		t_correlated_events->GetEntry(i);
		x.n_correntries = max_correvent;
		x.corrrun = run;
		Corr.push_back(x);
	}

	m_dummy = new TFile("dummy1.root", "recreate");
	t_correlated_events->Print();
	t_correlated_events->Write();

	h_neighbleft = new TH2D("Reb_hits_in_x_neighbleft_occupancy", "Summed Hitmap of every 3*74.5=223.5um in x; x in mm (Resolution 5um), left neighbour, occupancy", (int) (6*p/Dx)+1, -3*p, +3*p,N_thresholdbins,min_threshold,max_threshold);
	h_neighbright = new TH2D("Reb_hits_in_x_neighbright_occupancy", "Summed Hitmap of every 3*74.5=223.5um in x; x in mm (Resolution 5um), right neighbour, occupancy", (int)(6 * p / Dx) + 1, -3*p, +3*p, N_thresholdbins, min_threshold, max_threshold);
	h_seed = new TH2D("Reb_hits_in_x_seed_occupancy", "Summed Hitmap of every 3*74.5=223.5um in x; x in mm (Resolution 5um), seed, occupancy", (int)(6* p / Dx) + 1, -3*p, +3*p, N_thresholdbins, min_threshold, max_threshold);
	h_neighbleft_eff = new TH2D("Reb_hits_in_x_neighbleft_efficiency", "Summed Hitmap of every 3*74.5=223.5um in x; x in mm (Resolution 5um), left neighbour, efficiency", (int)(6 * p / Dx) + 1, -3 * p, +3 * p, N_thresholdbins, min_threshold, max_threshold);
	h_neighbright_eff = new TH2D("Reb_hits_in_x_neighbright_efficiency", "Summed Hitmap of every 3*74.5=223.5um in x; x in mm (Resolution 5um), right neighbout, efficiency", (int)(6 * p / Dx) + 1, -3 * p, +3 * p, N_thresholdbins, min_threshold, max_threshold);
	h_seed_eff = new TH2D("Reb_hits_in_x_seed_efficiency", "Summed Hitmap of every 3*74.5=223.5um in x; x in mm (Resolution 5um), seed, efficiency", (int)(6 * p / Dx) + 1, -3 * p, +3 * p, N_thresholdbins, min_threshold, max_threshold);
	h_hitx = new TH2D("Hits_in_x", "Summed Hitmap of every 3*74.5=223.5um in x; x in mm (Resolution 5um), total", (int)(6 * p / Dx) + 1, -3 * p, +3 * p, N_thresholdbins, min_threshold, max_threshold);
	h_hitmap_DUT = new TH2D("h_hitmap_DUT", "Hitmap on DUT", (int)((xmax-xmin)/ Dx) + 1, xmin, xmax, (int)((ymax - ymin) / Dx) + 1, ymin, ymax);
	h_hitmap_DUT_multiple = new TH2D("h_hitmap_DUT_multiple", "Hitmap on DUT with multiplicity", (int)((xmax - xmin) / Dx) + 1, xmin, xmax, (int)((ymax - ymin) / Dx) + 1, ymin, ymax);
	h_hitmap_m26 = new TH2D("h_hitmap_m26", "hitmap on mimosa", (int)((xmax - xmin) / Dx) + 1, xmin, xmax, (int)((ymax - ymin) / Dx) + 1, ymin, ymax);
	h_hitmap_losthits = new TH2D("h_hitmap_losthits", "Hitmap of lost hits", (int)((xmax - xmin) / Dx) + 1, xmin, xmax, (int)((ymax - ymin) / Dx) + 1, ymin, ymax);

	Mask.push_back(785);	// TO DO: read from gear
	Mask.push_back(784);	// TO DO: read from gear
	Mask.push_back(820);	// TO DO: read from gear
	Mask.push_back(821);	// TO DO: read from gear
	Mask.push_back(822);	// TO DO: read from gear
	Mask.push_back(823);	// TO DO: read from gear
	Mask.push_back(824);	// TO DO: read from gear
	Mask.push_back(863);	// TO DO: read from gear
	Mask.push_back(864);	// TO DO: read from gear
	Mask.push_back(865);	// TO DO: read from gear
	Mask.push_back(866);	// TO DO: read from gear

}

s_process_collection_modulo_forPositions::~s_process_collection_modulo_forPositions() {

}

bool s_process_collection_modulo_forPositions::process_file(FileProberties* fileP) {

	m_plotCollection = sct_corr::create_plot_collection();
	m_plotCollection->addFile(fileP->getTfile());
	m_plotCollection->setOutputFile(m_dummy);
	m_file_fitter = std::make_shared<sct_files::fitter_file>(m_plotCollection, get_gear());

	auto pl = m_file_fitter->get_collection();

	printf("%s\n", TFile::CurrentFile()->GetName());
	printf("- Filename: %s.\n", fileP->getTfile()->GetName());

	// Load of trees and initializations
	h_hitmap_DUT->Reset();
	h_hitmap_DUT_multiple->Reset();
	h_hitmap_m26->Reset();
	h_hitmap_losthits->Reset();

	FILE *filetxt;
	char temp[128];
	sprintf(temp,"tracks_%d.txt",(int)fileP->m_runNumber);
	filetxt=fopen(temp,"w");
	TTree *tree = (TTree*)fileP->getTfile()->Get("zsdata_strip");
	TTree *tree2 = (TTree*)fileP->getTfile()->Get("GBL_tracks");
	std::vector<double > *ID = 0;
	std::vector<double > *x = 0;
	std::vector<double > *y = 0;
	Int_t event_nr = 0;

	std::vector<double > *ID2 = 0;
	std::vector<double > *x2 = 0;
	std::vector<double > *y2 = 0;
	Int_t event_nr2 = 0;
	std::vector<double > *phi = 0;
	std::vector<double > *chi2 = 0;
	std::vector<double > *ndf = 0;
	std::vector<double > *theta = 0;

	tree->SetBranchAddress("ID", &ID);
	tree->SetBranchAddress("x", &x);
	tree->SetBranchAddress("y", &y);
	tree->SetBranchAddress("event_nr", &event_nr);

	tree2->SetBranchAddress("ID", &ID2);
	tree2->SetBranchAddress("x", &x2);
	tree2->SetBranchAddress("y", &y2);
	tree2->SetBranchAddress("event_nr", &event_nr2);
	tree2->SetBranchAddress("phi", &phi);
	tree2->SetBranchAddress("chi2", &chi2);
	tree2->SetBranchAddress("ndf", &ndf);
	tree2->SetBranchAddress("theta", &theta);
	int n_entries = tree->GetEntries();
	int n_entries2 = tree2->GetEntries();
/*
	TTree *t_correlated_events = new TTree();
	t_correlated_events = (TTree*)m_dummy->Get("correlated_events");
	int run= 0;
	int max_correvent = 0;
	t_correlated_events->SetBranchAddress("run", &run);
	t_correlated_events->SetBranchAddress("max_correvent", &max_correvent);
	for (int i = 0;i < t_correlated_events->GetEntries();i++) {
		t_correlated_events->GetEntry(i);
		if (run == (int)fileP->m_runNumber) {
			n_entries = max_correvent;
			break;
		}
	}
	*/
	for (int i = 0;i < Corr.size();i++) {
		if (Corr.at(i).corrrun== (int)fileP->m_runNumber) {
			n_entries = Corr.at(i).n_correntries;
			break;
		}
	}

printf("- Max correlated events for run %d = %d.\n", (int)fileP->m_runNumber, n_entries);

	auto ActiveStrips = get_xml_input()->globalConfig().AvtiveStrips();

	int plane;
	int id_inclu;
	double tiltedpos_x, tiltedpos_y;
	int closest_strip;
	double d;

	// Modulo 3
	for (int i = 0; i < n_entries; i++) {
		tree->GetEntry(i);
		plane = 0;
		if (x->size() == 0)
			continue;
		tree2->GetEntry(i);
		if (x2->size() == 0)
			continue;
		while (ID->at(0) != ID2->at(plane))
			plane++;
		tiltedpos_x = (TMath::Cos(alpha) * x2->at(plane)) - (TMath::Sin(alpha) * y2->at(plane));
		tiltedpos_y = (TMath::Sin(alpha) * x2->at(plane)) + (TMath::Cos(alpha) * y2->at(plane));
		if (tiltedpos_x >= XCutmax || tiltedpos_x <= XCutmin || tiltedpos_y >= YCutmax || tiltedpos_y <= YCutmin) {	// excluding hit out of cut and on masked channels
//			printf("Hit at %f,%f, excluded.\n", tiltedpos_x, tiltedpos_y);
			continue;
		}
		// Fill hit in x histo
		closest_strip = (int)ActiveStrips.getMin();
		d = tiltedpos_x - (Dxcentering + (closest_strip)*p);
		for (int strip = (int)ActiveStrips.getMin()+1;strip < (int)ActiveStrips.getMax();strip++)
			if (abs(tiltedpos_x - (Dxcentering + (strip)*p)) < abs(d)) {
				closest_strip = strip;
				d = tiltedpos_x - (Dxcentering + (strip)*p);
			}
		if (isMasked(closest_strip) || closest_strip >= (ActiveStrips.getMin() + ((int)(ActiveStrips.getMax() - ActiveStrips.getMin()) / 3) * 3) || (closest_strip > ActiveStrips.getMax()) || (closest_strip < ActiveStrips.getMin())) {
//			printf("- Hit of track %d at %f around strip %d removed.\n", i,tiltedpos_x,closest_strip);
		}
		else {
			id_inclu = ((int)((closest_strip - ActiveStrips.getMin()) - ((int)((closest_strip - ActiveStrips.getMin()) / 3)) * 3))-1;
			h_hitmap_m26->Fill(tiltedpos_x, tiltedpos_y);
			if (x->size() == 0)
				h_hitmap_losthits->Fill(tiltedpos_x, tiltedpos_y);
			else {
				for (int id_hit = 0; id_hit < x->size(); id_hit++)
					if ((isMasked((int)x->at(id_hit)) || x->at(id_hit) >= (ActiveStrips.getMin() + ((int)(ActiveStrips.getMax() - ActiveStrips.getMin()) / 3) * 3) || (x->at(id_hit) > ActiveStrips.getMax()) || (x->at(id_hit) < ActiveStrips.getMin())) == 0) {	// excluding hit out of cut and on masked channels
						h_hitmap_losthits->Fill(tiltedpos_x, tiltedpos_y);
						printf("-(id,strip)= %d, %f\n", id_hit, x->at(id_hit));
					}
			}
			h_hitx->Fill(d + (id_inclu)*p, fileP->m_Threshold);
			fprintf(filetxt,"- Track: (xposition,yposition)=%f , %f \t (closest strip, stripinclu, distance)=%d , %d, %f)\n", tiltedpos_x, tiltedpos_y, closest_strip, id_inclu, d);
		}

		int flag = 0;	// for counting hits withour considering multiplicity on h_hitmap_DUT
		for (int id_hit = 0; id_hit < x->size(); id_hit++) {
			if (isMasked((int)x->at(id_hit)) || x->at(id_hit) >= (ActiveStrips.getMin() + ((int)(ActiveStrips.getMax() - ActiveStrips.getMin()) / 3) * 3) || (x->at(id_hit) > ActiveStrips.getMax()) || (x->at(id_hit) < ActiveStrips.getMin()))	// excluding hit out of cut and on masked channels
				continue;
			h_hitmap_DUT_multiple->Fill(tiltedpos_x, tiltedpos_y);
			if (flag == 0) {
				h_hitmap_DUT->Fill(tiltedpos_x, tiltedpos_y);
				flag++;
			}
			id_inclu = ((int)((x->at(id_hit) - ActiveStrips.getMin()) - ((int)((x->at(id_hit) - ActiveStrips.getMin()) / 3)) * 3))-1;
			switch (id_inclu) {
			case -1:
				h_neighbleft->Fill(tiltedpos_x - (Dxcentering + (x->at(id_hit))*p) + (id_inclu)*p, fileP->m_Threshold);
				break;
			case 0:
				h_seed->Fill(tiltedpos_x - (Dxcentering + (x->at(id_hit))*p) + (id_inclu)*p, fileP->m_Threshold);
				break;
			case +1:
				h_neighbright->Fill(tiltedpos_x - (Dxcentering + (x->at(id_hit))*p) + (id_inclu)*p, fileP->m_Threshold);
				break;
			default:
				perror("process_file(): id_inclu out of range.");
				break;
			}
			fprintf(filetxt, "\t- Hit: (xposition,yposition)=%f , %f \t (strip, stripinclu, multiplicity)=%f , %d , %d)\n", tiltedpos_x, tiltedpos_y, x->at(id_hit), id_inclu, id_hit);
//			printf("- (x, id, tilte pos x, id_inclu,closest_strip,d) = (%f, %f, %f, %d,%d,%f)\n", x2->at(plane), x->at(id_hit), tiltedpos_x, id_inclu,closest_strip,d);
//			printf("- tiltedpos_x - (x_seedcenter + (x->at(id_hit) - ActiveStrips.getMin())*p + 0)  =  %f - (%f + (%f - %d)*%f + 0) = %f\n", tiltedpos_x, x_seedcenter, x->at(id_hit), ActiveStrips.getMin(),p, tiltedpos_x - (x_seedcenter + (x->at(id_hit) - ActiveStrips.getMin())*p + 0));
		}
	}
	// Efficiency computation
	h_neighbleft_eff->Divide(h_neighbleft,h_hitx);
	h_neighbright_eff->Divide(h_neighbright,h_hitx);
	h_seed_eff->Divide(h_seed,h_hitx);

#ifdef _DEBUG
	new TCanvas();
	h_neighbleft->Draw("colz");
	new TCanvas();
	h_seed->Draw("colz");
	new TCanvas();
	h_neighbright->Draw("colz");
	new TCanvas();
	h_hitx->Draw("colz");
	new TCanvas();
	h_neighbleft_eff->Draw("colz");
	new TCanvas();
	h_seed_eff->Draw("colz");
	new TCanvas();
	h_neighbright_eff->Draw("colz");

	h_neighbleft->Write("h_neighbleft", 2);
	h_neighbright->Write("h_neighbright", 2);
	h_seed->Write("h_seed",2);
	h_hitx->Write("h_hitx", 2);
	h_neighbleft_eff->Write("h_neighbleft_eff", 2);
	h_neighbright_eff->Write("h_neighbright_eff", 2);
	h_seed_eff->Write("h_seed_eff",2);
	h_hitmap_m26->Write("h_hitmap_m26");
	h_hitmap_DUT->Write("h_hitmap_DUT");
	h_hitmap_DUT_multiple->Write("h_hitmap_DUT_multiple");
	h_hitmap_losthits->Write("h_hitmap_losthits");
#endif
	fclose(filetxt);
	
/*
	m_residualEffieciency->Draw();
	m_instripEfficiency->Draw();
	push2outputEvent(m_outputl, *m_residualEffieciency->getEfficiency_map(), *m_residualEffieciency->get_total(), ID_t(0));
	push2outputEvent(m_outputl, *m_instripEfficiency->getEfficiency_map(), *m_instripEfficiency->getHits(), ID_t(1));
	*/
	return true;
}
bool s_process_collection_modulo_forPositions::isMasked(int channel) {
	for (int i = 0;i < Mask.size();i++) {
		if (channel == Mask.at(i))
			return true;
	}
	return false;
}


void s_process_collection_modulo_forPositions::saveHistograms(TFile* outPutFile /*= nullptr */, xmlImputFiles::MinMaxRange<double>* residual_cut /*= nullptr */) {

#ifdef _DEBUG
	new TCanvas();
#endif
//	m_instripClusterSize->Draw();
#ifdef _DEBUG
	new TCanvas();
#endif // _DEBUG
//	m_instripEfficiency->Draw();
#ifdef _DEBUG
	new TCanvas();
#endif // _DEBUG
//	m_residualEffieciency->Draw();

/*	if (outPutFile) {
		outPutFile->Add(m_instripClusterSize->getHistogram());
		outPutFile->Add(m_instripEfficiency->getEfficiency_map());
		outPutFile->Add(m_residualEffieciency->getEfficiency_map());
		outPutFile->Add(m_residualEffieciency->get_total());
	}*/
}

std::string s_process_collection_modulo_forPositions::get_suffix() const {
	return "Modulo_forPositions";
}
///////////////////////////// Correlation check

s_process_collection_correlation_check::s_process_collection_correlation_check(Parameter_ref par) {

	// Retriving alignment
	TFile* m_centering = new TFile("endcap_450V_centering.root", "read");	// TO DO: take name of collection
	TTree *t_centering = new TTree();
	t_centering = (TTree*)m_centering->Get("centering");
	int stripref = 0;
	double x_seedcenter = 0;
	centering Centering;
	t_centering->SetBranchAddress("strip", &stripref);
	t_centering->SetBranchAddress("x_seedcenter", &x_seedcenter);
	t_centering->GetEntry(t_centering->GetEntries()-1);
	Centering.stripref = stripref;
	Centering.x_seedcenter = x_seedcenter;
	Dxcentering = Centering.x_seedcenter - Centering.stripref*p;
	printf("- Centering: strip %d at %f, traslation of %f.\n", Centering.stripref, Centering.x_seedcenter, Dxcentering);

	m_dummy = new TFile("dummy1.root", "recreate");
	
	corrGauss = new TF1("fitCorrGauss", "gaus", xmin_corr, xmax_corr);
	treeout = new TTree("correlated_events","Correlated event limit per run");
	treeout->Branch("run", &run, "run/I");
	treeout->Branch("max_correvent", &max_correvent, "max_correvent/I");
	treeout->Write("correlated_events");
	delete treeout;

	Mask.push_back(785);	// TO DO: read from gear
	Mask.push_back(784);	// TO DO: read from gear
	Mask.push_back(820);	// TO DO: read from gear
	Mask.push_back(821);	// TO DO: read from gear
	Mask.push_back(822);	// TO DO: read from gear
	Mask.push_back(823);	// TO DO: read from gear
	Mask.push_back(824);	// TO DO: read from gear
	Mask.push_back(863);	// TO DO: read from gear
	Mask.push_back(864);	// TO DO: read from gear
	Mask.push_back(865);	// TO DO: read from gear
	Mask.push_back(866);	// TO DO: read from gear

}

s_process_collection_correlation_check::~s_process_collection_correlation_check() {

}

bool s_process_collection_correlation_check::process_file(FileProberties* fileP) {

	printf("- Current directory %s...\n", gDirectory->CurrentDirectory()->GetName());
	m_plotCollection = sct_corr::create_plot_collection();
	m_plotCollection->addFile(fileP->getTfile());
	m_plotCollection->setOutputFile(m_dummy);
	m_file_fitter = std::make_shared<sct_files::fitter_file>(m_plotCollection, get_gear());

	auto pl = m_file_fitter->get_collection();

	printf("- Filename: %s.\n", fileP->getTfile()->GetName());
	// Load of trees and initializations
	treeout = (TTree*)m_dummy->Get("correlated_events");
	treeout->SetBranchAddress("run", &run);
	treeout->SetBranchAddress("max_correvent", &max_correvent);

	TTree *tree = (TTree*)fileP->getTfile()->Get("zsdata_strip");
	TTree *tree2 = (TTree*)fileP->getTfile()->Get("GBL_tracks");
	std::vector<double > *ID = 0;
	std::vector<double > *x = 0;
	std::vector<double > *y = 0;
	Int_t event_nr = 0;

	std::vector<double > *ID2 = 0;
	std::vector<double > *x2 = 0;
	std::vector<double > *y2 = 0;
	Int_t event_nr2 = 0;
	std::vector<double > *phi = 0;
	std::vector<double > *chi2 = 0;
	std::vector<double > *ndf = 0;
	std::vector<double > *theta = 0;

	tree->SetBranchAddress("ID", &ID);
	tree->SetBranchAddress("x", &x);
	tree->SetBranchAddress("y", &y);
	tree->SetBranchAddress("event_nr", &event_nr);

	tree2->SetBranchAddress("ID", &ID2);
	tree2->SetBranchAddress("x", &x2);
	tree2->SetBranchAddress("y", &y2);
	tree2->SetBranchAddress("event_nr", &event_nr2);
	tree2->SetBranchAddress("phi", &phi);
	tree2->SetBranchAddress("chi2", &chi2);
	tree2->SetBranchAddress("ndf", &ndf);
	tree2->SetBranchAddress("theta", &theta);
	int n_entries = tree->GetEntries();
	int n_entries2 = tree2->GetEntries();

	auto ActiveStrips = get_xml_input()->globalConfig().AvtiveStrips();

	int plane;
	int id_inclu;
	double tiltedpos_x, tiltedpos_y;
	int closest_strip;
	double d;

	h_corr_DUT = new TH2D("Hits_in_x_events_DUT", "Summed Hitmap for orrelation check, DUT", (int)((xmax - xmin) / p), xmin, xmax, (int)(n_entries/Nslot),0, (int)(n_entries / Nslot));
	h_corr_m26 = new TH2D("Hits_in_x_events_m26", "Summed Hitmap for orrelation check, m26", (int)((xmax - xmin) / p), xmin, xmax, (int)(n_entries / Nslot), 0, (int)(n_entries / Nslot));
	h_corr = new TH2D("correlation_in_x_events", "Correlation between hitmaps for correlation check", (int)((xmax_corr - xmin_corr) / p), xmin_corr, xmax_corr, (int)(n_entries / Nslot), 0, (int)(n_entries / Nslot));
	h_corr2 = new TH2D("correlation_in_x_events_2", "Correlation between for correlation check_2", (int)((xmax_corr - xmin_corr) / p), xmin_corr, xmax_corr, (int)(n_entries / Nslot), 0, (int)(n_entries / Nslot));

	// Correlation check
	for (int j = 0;j < (int)(n_entries / Nslot);j++) {
//		h_corr_DUT->Reset();
//		h_corr_m26->Reset();
//		h_corr->Reset();

		for (int i = 0; i < Nslot; i++) {
			if ((j*Nslot + i) >= n_entries)
				break;
			tree->GetEntry((j*Nslot + i));
			if (x->size() == 0)
				continue;
			plane = 0;
			tree2->GetEntry((j*Nslot + i));
			if (x2->size() == 0)
				continue;
			while (ID->at(0) != ID2->at(plane))
				plane++;
			tiltedpos_x = (TMath::Cos(alpha) * x2->at(plane)) - (TMath::Sin(alpha) * y2->at(plane));
			tiltedpos_y = (TMath::Sin(alpha) * x2->at(plane)) + (TMath::Cos(alpha) * y2->at(plane));
			if (tiltedpos_x >= XCutmax || tiltedpos_x <= XCutmin || tiltedpos_y >= YCutmax || tiltedpos_y <= YCutmin)	// excluding hit out of cut and on masked channels
				continue;
			h_corr_m26->Fill(x2->at(plane),j);
			for (int id_hit = 0; id_hit < x->size(); id_hit++) {
				if (isMasked((int)x->at(id_hit)) || (x->at(id_hit) > ActiveStrips.getMax()) || (x->at(id_hit) < ActiveStrips.getMin()))	// excluding hit out of cut and on masked channels
					continue;
				h_corr_DUT->Fill(Dxcentering + ((x->at(id_hit))*p),j);
				h_corr2->Fill(Dxcentering + ((x->at(id_hit))*p) - x2->at(plane), j);
			}
		}
		// Cross-correlation
		double corr;
		for (int k = 0;k < h_corr->GetXaxis()->GetNbins();k++) {
			corr = 0;
			for (int l = 0;l < h_corr_m26->GetXaxis()->GetNbins();l++) {
				if ((h_corr_DUT->GetXaxis()->FindBin(h_corr->GetXaxis()->GetBinCenter(k) - h_corr_m26->GetXaxis()->GetBinCenter(l)))>=0 && (h_corr_DUT->GetXaxis()->FindBin(h_corr->GetXaxis()->GetBinCenter(k) - h_corr_m26->GetXaxis()->GetBinCenter(l)))<h_corr_DUT->GetXaxis()->GetNbins())
					corr += (h_corr_m26->GetBinContent(l,j)*h_corr_DUT->GetBinContent(h_corr_DUT->GetXaxis()->FindBin(h_corr->GetXaxis()->GetBinCenter(k)-h_corr_m26->GetXaxis()->GetBinCenter(l)),j));
			}
			h_corr->SetBinContent(k,j,corr);
		}
	}
	// Fit and check
	run = (int)fileP->m_runNumber;
	max_correvent = n_entries;
	int x_max,z_max,j_max;
	int max = h_corr2->GetMaximum();
	h_corr2->GetBinXYZ(h_corr2->GetMaximumBin(), x_max, j_max, z_max);
	h_corr2->ProjectionX("p_x_corr", j_max, j_max)->Fit(corrGauss, "R");
	corrGauss->SetRange(Dxcentering - 10*p, Dxcentering + 10 * p);
	int bin_min, bin_max;
	bin_min = h_corr2->GetXaxis()->FindBin( - 0.5);
	bin_max = h_corr2->GetXaxis()->FindBin( + 0.5);

//	corrGauss->SetParLimits(1, corrGauss->GetParameter(1) - 5 * p, corrGauss->GetParameter(1) + 5 * p);
	for (int j = 1;j < (int)(n_entries / Nslot)+1;j++) {
		h_corr2->ProjectionX("p_x_corr", j, j)->Fit(corrGauss, "R");
/*
#ifdef _DEBUG
		new TCanvas();
		h_corr2->ProjectionX("p_x_corr", j, j)->Draw();
		h_corr2->ProjectionX("p_x_corr", j, j)->Write("h_corr2_p_x");
#endif
*/
#ifdef _DEBUG
		printf("- bins: %d, %d\n", bin_min, bin_max);
		printf("- Maximum correlated event (slot), sigma, max sigma, max of histo at = %d (%d), %f, %f, %d at %d\n", max_correvent, j - 1, corrGauss->GetParameter(2), sigma_max,max,j_max);
		printf("- Range, Sigma limits= %f - %f, %f - %f\n", corrGauss->GetParameter(1) - 10 * p, corrGauss->GetParameter(1) + 10 * p, corrGauss->GetParameter(1) - 5 * p, corrGauss->GetParameter(1) + 5 * p);
#endif
		if (h_corr2->ProjectionX("p_x_corr", j, j)->GetMaximum()!=0 && ((h_corr2->ProjectionX("p_x_corr", j, j)->GetMaximum()<0.1*max) || (h_corr2->ProjectionX("p_x_corr", j, j)->GetMaximumBin()<bin_min || h_corr2->ProjectionX("p_x_corr", j, j)->GetMaximumBin()>bin_max))) {//corrGauss->GetParameter(2) > sigma_max && 
			max_correvent = (j - 1)*Nslot;
#ifdef _DEBUG
			printf("- Maximum correlated event %d for run number %d.\n", max_correvent, fileP->m_runNumber);
#endif
			break;
		}
	}

#ifdef _DEBUG
//	new TCanvas();
//	h_corr_m26->Draw();
//	h_corr_m26->Write("h_corr_m26");
//	new TCanvas();
//	h_corr_DUT->Draw();
//	h_corr_DUT->Write("h_corr_DUT");
//	new TCanvas();
//	h_corr->Draw();
//	h_corr->Write("h_corr");
	TCanvas *c=new TCanvas();
	h_corr2->Draw();
	TLine *l = new TLine(h_corr2->GetXaxis()->GetXmin(), (max_correvent/Nslot), h_corr2->GetXaxis()->GetXmax(), (max_correvent / Nslot));
	printf("%f, %f, %f, %f");
	l->SetLineWidth(8);
	l->SetLineColor(kRed);
	l->Draw();
	c->Write("corr_show");
	h_corr2->Write("h_corr2");
	h_corr2->Reset();
#endif
	treeout->Fill();
	TTree *newtreeout = treeout->CloneTree();
	newtreeout->Print();
	newtreeout->Write("correlated_events",TObject::kOverwrite);
	delete newtreeout;

	delete h_corr_m26;
	delete h_corr_DUT;
	delete h_corr;
	delete h_corr2;
printf("- Run number %d: maximum event=%d.\n",run,max_correvent);
/*
char *outfile = (char *)gDirectory->GetFile()->GetName();
gDirectory->GetFile()->Close();
	// Save only correlated event files
	TTree **Tree = new TTree*[fileP->getTfile()->GetListOfKeys()->GetSize()];
	for (int i = 0;i < fileP->getTfile()->GetListOfKeys()->GetSize();i++)
		Tree[i] = (TTree*)fileP->getTfile()->Get(fileP->getTfile()->GetListOfKeys()->At(i)->GetName());
	char filename[128];
	char temp[128];
	char dir[128];
	strcpy(temp, fileP->getTfile()->GetName());
	const char *filepart=strpbrk(temp,"\/");
	strcpy(filename,"C:\\data\\2015_10_endcap");
	strcat(filename, "\/part_");
	strcat(filename, filepart+1);

	TFile *f_copy = new TFile(filename,"recreate");
	printf("- Current directory %s...\n", gDirectory->CurrentDirectory()->GetName());
	TTree **Tree_copy = new TTree*[fileP->getTfile()->GetListOfKeys()->GetSize()];
	for (int i = 0;i < fileP->getTfile()->GetListOfKeys()->GetSize();i++) {
		Tree_copy[i] = Tree[i]->CloneTree(max_correvent);
#ifdef _DEBUG
		Tree_copy[i]->Print();
#endif
		Tree_copy[i]->AutoSave();
	}
	printf("- Saving %s...\n", filename);
	f_copy->Close();
	delete f_copy;
	// Go back
	TFile *f=new TFile("C:\\Users\\Riccardo\\Documents\\Git\\SCT_correlations\\bin\\endcap_450V_Correlation_check.root","recreate");
	printf("- Current directory %s...\n", gDirectory->CurrentDirectory()->GetName());
*/
	return true;
}
bool s_process_collection_correlation_check::isMasked(int channel) {
	for (int i = 0;i < Mask.size();i++) {
		if (channel == Mask.at(i))
			return true;
	}
	return false;
}

void s_process_collection_correlation_check::saveHistograms(TFile* outPutFile /*= nullptr */, xmlImputFiles::MinMaxRange<double>* residual_cut /*= nullptr */) {

#ifdef _DEBUG
	new TCanvas();
#endif
	//	m_instripClusterSize->Draw();
#ifdef _DEBUG
	new TCanvas();
#endif // _DEBUG
	//	m_instripEfficiency->Draw();
#ifdef _DEBUG
	new TCanvas();
#endif // _DEBUG
	//	m_residualEffieciency->Draw();

	/*	if (outPutFile) {
	outPutFile->Add(m_instripClusterSize->getHistogram());
	outPutFile->Add(m_instripEfficiency->getEfficiency_map());
	outPutFile->Add(m_residualEffieciency->getEfficiency_map());
	outPutFile->Add(m_residualEffieciency->get_total());
	}*/
}

std::string s_process_collection_correlation_check::get_suffix() const {
	return "Correlation_check";
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


using namespace sct_corr;
registerBaseClassDef(processorBase);

registerProcessor("Modulo", s_process_collection_modulo);
registerProcessor("Modulo_ex", s_process_collection_modulo_ex);
registerProcessor("Standard", s_process_collection_standard);
registerProcessor("Standard_second", s_process_collection_standard_second);
registerProcessor("Modulo_second", s_process_collection_modulo_second);
registerProcessor("Residual", s_process_collection_residual);
registerProcessor("Residual_second", s_process_collection_residual_second);
registerProcessor("Modulo_forPositions", s_process_collection_modulo_forPositions);
registerProcessor("Correlation_check", s_process_collection_correlation_check);
registerProcessor("Centering", s_process_collection_centering);

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
