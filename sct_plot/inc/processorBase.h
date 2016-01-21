#ifndef s_process_files_h__
#define s_process_files_h__
#include "Rtypes.h"
#include "internal/platform.hh"
#include "TH1.h"


#include <memory>
#include "s_plot_collection.h"
#include "s_cuts.h"
#include "sct_events/rootEventRunOutput.hh"
#include "sct_event_buffer.h"
#include "plane_def.h"
#include "xml_helpers/xml_fileList.hh"

#include "factory.hh"
#define  registerProcessor(processorName,processorType) registerClass(processorBase,processorType,processorName)
namespace sct_corr {
class treeCollection_ouput;
struct Xgear;
class inStripEfficiency;

class inStripClusterSize;
class residual_efficienct;
class processorBatch;
}
namespace sct_files {
class fitter_file;
}

namespace xmlImputFiles {
class XML_imput_file;
}
class FileProberties {
public:
  TFile* getTfile() const;
  void setTFile(TFile* file);
  void setTFile(std::shared_ptr<TFile> file);
  double m_Threshold = 0;
  double m_runNumber = 0;
  double m_HV = 0;
  std::weak_ptr<sct_corr::processorBatch> m_batch;
private:
  std::shared_ptr<TFile> m_fileOwnd;
  TFile* m_file = nullptr;
};

#include "TTree.h"
#include "TH2.h"

namespace sct_corr {



class DllExport Processor {
public:
  enum ProcessState {
    ok,
    done
  };
  virtual ~Processor();

  virtual ProcessState ProceessEvent() = 0;

  void push_processor_batch(std::weak_ptr<processorBatch> batch);
  std::weak_ptr<processorBatch> get_batch();
private:

  std::weak_ptr<processorBatch> m_batch;
};


class DllExport processorBatch :public Processor {
public:
  processorBatch() {}
  virtual ~processorBatch() {}

  virtual ProcessState ProceessEvent() override;
  void push_processor(Processor* processor_);
  void loop();
private:
  std::vector<Processor*> m_processors;
};
std::shared_ptr<processorBatch> create_batch();

class DllExport ProcessorXML_loader:public Processor {
public:
  ProcessorXML_loader(const std::string& xmlInputFileName, std::string path__, std::string outputPath = ".");
  virtual ~ProcessorXML_loader() {}


  virtual ProcessState ProceessEvent() override;

  const FileProberties* getData() const;
  
private:

  std::shared_ptr<xmlImputFiles::XML_imput_file> m_input_files_xml;



  std::shared_ptr<sct_corr::Xgear> m_gear;
  std::shared_ptr<FileProberties> m_files;


  std::shared_ptr<processorBatch> m_owndBatch;
  TFile* m_outpuFile = nullptr;
  std::string m_outname,m_path;
  size_t m_iterrator = 0;
};


class DllExport processorEfficiency :public Processor {
public:
  processorEfficiency(FileProberties* fileProb);
  virtual ProcessState ProceessEvent() override;

private:

  void process_reset();
  void process_set_run_prob();
  FileProberties* m_file;
  //sct_corr::rootEventRunOutput m_outputl;
  std::shared_ptr<sct_corr::treeCollection_ouput> m_outputTree;
  sct_corr::sct_event_buffer m_buffer;



  std::shared_ptr<TH1D> m_Residual;
  std::shared_ptr<TH1D> m_Hits_total;
  std::shared_ptr<TH1D> m_Hits_with_DUT_Hits;
  std::shared_ptr<TH1D> m_Efficieny_map;
  std::shared_ptr<TH1D> m_Efficieny_trueHits;
  std::shared_ptr<TH2D> m_resVSMissing;
  std::shared_ptr<TH2D> m_ResidualVsEvent;
  s_plane_collection m_res_VS_event;
  s_plane_collection_correlations m_output_planes;

  std::shared_ptr<sct_corr::plot_collection> m_plotCollection;
  std::shared_ptr <sct_files::fitter_file> m_file_fitter;

};

TYPE_CLASS_PTR(output_TFile_ptr, TFile*);
class DllExport processorBase {
public:
  using MainType = std::string;
  using Parameter_t = std::string;
  using Parameter_ref = const std::string&;
  processorBase();
  virtual ~processorBase();
  void setOutputName(const char* name);
  void push_files(TFile* _file, double Threshold, double runNumber);
  void push_files(const char* _fileName, double Threshold, double runNumber, double HV);

  int Add_XML_RunList(const std::string& xmlInputFileName, std::string path__, std::string outputPath = ".", int element = -1);

  void setGearFile(const char* name);
  void setPrintout(bool print);
  bool process();
  virtual void saveHistograms(TFile* outPutFile = nullptr, xmlImputFiles::MinMaxRange<double>* residual_cut = nullptr) = 0;

  const xmlImputFiles::XML_imput_file* get_xml_input() const;
  const sct_corr::Xgear* get_gear() const;
protected:
  sct_corr::rootEventRunOutput m_outputl;
  std::shared_ptr<sct_corr::treeCollection_ouput> m_outputTree;
  sct_corr::sct_event_buffer m_buffer;
private:
  virtual std::string get_suffix() const = 0;

  std::shared_ptr<xmlImputFiles::XML_imput_file> m_input_files_xml;

  void process_set_run_prob(const FileProberties& fileP);


  std::shared_ptr<sct_corr::Xgear> m_gear;


  void start_collection(output_TFile_ptr file__) ;
  virtual  bool process_file(FileProberties* fileP) = 0;
  virtual void end_collection() {}


  std::vector<FileProberties> m_files;



  std::string m_outname;




};
}
class DllExport s_process_collection_standard :public sct_corr::processorBase {
public:
  s_process_collection_standard(Parameter_ref par);
  virtual ~s_process_collection_standard();
  Long64_t DrawResidual(Double_t min_X, Double_t max_X);
  Long64_t DrawResidual();
  Long64_t DrawResidualVsEvent(Double_t min_X, Double_t max_X);
  Long64_t DrawResidualVsEvent();
  Long64_t DrawResidualVsMissingCordinate(Double_t min_X, Double_t max_X);
  Long64_t DrawResidualVsMissingCordinate();
  Long64_t Draw_Efficinecy_map();
  Long64_t Draw_Hit_map();
  Long64_t Draw_DUT_Hits_map();
  TH2D* getResidualVsMissingCordinate();
  virtual void saveHistograms(TFile* outPutFile = nullptr, xmlImputFiles::MinMaxRange<double>* residual_cut = nullptr);

private:

  virtual  bool process_file(FileProberties* fileP) override;
  void extract_efficiency();
  void extract_hitMap();
  void extract_residual();
  void extract_rotation();
  void process_reset();



  virtual std::string get_suffix() const override;


  std::shared_ptr<TH1D> m_Residual;
  std::shared_ptr<TH1D> m_Hits_total;
  std::shared_ptr<TH1D> m_Hits_with_DUT_Hits;
  std::shared_ptr<TH1D> m_Efficieny_map;
  std::shared_ptr<TH1D> m_Efficieny_trueHits;
  std::shared_ptr<TH2D> m_resVSMissing;
  std::shared_ptr<TH2D> m_ResidualVsEvent;

  s_plane_collection m_res_VS_event;
  s_plane_collection_correlations m_output_planes;

  std::shared_ptr<sct_corr::plot_collection> m_plotCollection;
  std::shared_ptr <sct_files::fitter_file> m_file_fitter;
  TFile* m_dummy = nullptr;


};







class DllExport s_process_collection_modulo : public sct_corr::processorBase {
public:
  s_process_collection_modulo(Parameter_ref par);
  virtual ~s_process_collection_modulo();

  virtual void saveHistograms(TFile* outPutFile  = nullptr , xmlImputFiles::MinMaxRange<double>* residual_cut  = nullptr ) override;
private:
  
  virtual std::string get_suffix() const override;

  virtual  bool process_file(FileProberties* fileP) override;



  TFile* m_dummy = nullptr;

  s_plane_collection_correlations m_gbl_collection;

  std::shared_ptr<sct_corr::plot_collection> m_plotCollection;
  std::shared_ptr<sct_files::fitter_file> m_file_fitter;
  std::shared_ptr<sct_corr::inStripEfficiency> m_instripEfficiency;
  std::shared_ptr<sct_corr::inStripClusterSize> m_instripClusterSize;
  std::shared_ptr<sct_corr::residual_efficienct> m_residualEffieciency;
};

class DllExport s_process_collection_modulo_ex : public sct_corr::processorBase {
public:
  s_process_collection_modulo_ex(Parameter_ref);
  virtual ~s_process_collection_modulo_ex();

  virtual void saveHistograms(TFile* outPutFile = nullptr, xmlImputFiles::MinMaxRange<double>* residual_cut = nullptr) override;
private:

  virtual std::string get_suffix() const override;

  virtual  bool process_file(FileProberties* fileP) override;



  TFile* m_dummy = nullptr;

  s_plane_collection_correlations m_gbl_collection;

  std::shared_ptr<sct_corr::plot_collection> m_plotCollection;
  std::shared_ptr<sct_files::fitter_file> m_file_fitter;
  std::shared_ptr<sct_corr::inStripEfficiency> m_instripEfficiency;
};

class DllExport s_process_collection_residual : public sct_corr::processorBase {
public:
  s_process_collection_residual(Parameter_ref par);
  virtual ~s_process_collection_residual();

  virtual void saveHistograms(TFile* outPutFile = nullptr, xmlImputFiles::MinMaxRange<double>* residual_cut = nullptr) override;
private:

  virtual std::string get_suffix() const override;

  virtual  bool process_file(FileProberties* fileP) override;



  TFile* m_dummy = nullptr;

  s_plane_collection_correlations m_gbl_collection;

  std::shared_ptr<sct_corr::plot_collection> m_plotCollection;
  std::shared_ptr<sct_files::fitter_file> m_file_fitter;
  std::shared_ptr<sct_corr::residual_efficienct> m_residualEffieciency;
};

class DllExport s_process_collection_residual_second : public sct_corr::processorBase {
public:
  s_process_collection_residual_second(Parameter_ref par);
  virtual ~s_process_collection_residual_second();

  virtual void saveHistograms(TFile* outPutFile = nullptr, xmlImputFiles::MinMaxRange<double>* residual_cut = nullptr) override;
private:

  virtual std::string get_suffix() const override;

  virtual  bool process_file(FileProberties* fileP) override;



  TFile* m_dummy = nullptr;

  s_plane_collection_correlations m_gbl_collection;

  std::shared_ptr<sct_corr::plot_collection> m_plotCollection;
  std::shared_ptr<sct_files::fitter_file> m_file_fitter;
  std::shared_ptr<sct_corr::residual_efficienct> m_residualEffieciency;
};

class DllExport s_process_collection_modulo_second : public sct_corr::processorBase {
public:
  s_process_collection_modulo_second(Parameter_ref par);
  virtual ~s_process_collection_modulo_second();

  virtual void saveHistograms(TFile* outPutFile = nullptr, xmlImputFiles::MinMaxRange<double>* residual_cut = nullptr) override;
private:

  virtual std::string get_suffix() const override;

  virtual  bool process_file(FileProberties* fileP) override;



  TFile* m_dummy = nullptr;

  s_plane_collection_correlations m_gbl_collection;

  std::shared_ptr<sct_corr::plot_collection> m_plotCollection;
  std::shared_ptr<sct_files::fitter_file> m_file_fitter;
  std::shared_ptr<sct_corr::inStripEfficiency> m_instripEfficiency;
};


class DllExport s_process_collection_standard_second :public sct_corr::processorBase {
public:
  s_process_collection_standard_second(Parameter_ref par);
  virtual ~s_process_collection_standard_second();
  Long64_t DrawResidual(Double_t min_X, Double_t max_X);


  Long64_t Draw_Efficinecy_map();
  Long64_t Draw_Hit_map();
  Long64_t Draw_DUT_Hits_map();

  virtual void saveHistograms(TFile* outPutFile = nullptr, xmlImputFiles::MinMaxRange<double>* residual_cut = nullptr);

private:

  virtual  bool process_file(FileProberties* fileP) override;
  void extract_efficiency();
  void extract_hitMap();
  void extract_residual();
  
  void process_reset();



  virtual std::string get_suffix() const override;


  std::shared_ptr<TH1D> m_Residual;
  std::shared_ptr<TH1D> m_Hits_total;
  std::shared_ptr<TH1D> m_Hits_with_DUT_Hits;
  std::shared_ptr<TH1D> m_Efficieny_map;
  std::shared_ptr<TH1D> m_Efficieny_trueHits;


  s_plane_collection_correlations m_gbl_collection;

  s_plane_collection m_collection;
  std::shared_ptr<sct_corr::plot_collection> m_plotCollection;
  std::shared_ptr <sct_files::fitter_file> m_file_fitter;
  TFile* m_dummy = nullptr;


};

DllExport void display_registered_processors(std::ostream&);
DllExport std::unique_ptr<sct_corr::processorBase> create_processor(const std::string& processorName);


#endif // s_process_files_h__
