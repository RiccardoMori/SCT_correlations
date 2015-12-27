#ifndef scurve_fit_collection_h__
#define scurve_fit_collection_h__



#include "TGraphErrors.h"
#include "TCanvas.h"
#include "landgausFit.h"
#include <fstream>
#include <memory>
#include "S_Axis.h"
#include "getScurvesFromTTree.hh"

class TTree;

struct output_data {
  output_data() {}
  double mpv_fit = 0, gauss_sigma = 0, landau_sigma = 0, max_efficiency = 0, mpv_estimated = 0, sigma_estimated = 0, x_pos = 0, chiSquare = 0, mpv_fit_error = 0;

  void print(std::ostream& out) {
    out << x_pos << ";  " << max_efficiency << ";  " << mpv_fit << ";  " << mpv_fit_error << ";   " << landau_sigma << ";   " << gauss_sigma << ";  " << chiSquare << ";  " << mpv_estimated << ";  " << sigma_estimated << std::endl;
  }
  void printHeeader(std::ostream& out) {
    out << "x_pos" << ";  " << "max_efficiency" << ";  " << "mpv_fit" << ";  " << "mpv_fit_error" << ";   " << "landau_sigma" << ";   " << "gauss_sigma" << ";  " << "chiSquare" << ";  " << "mpv_estimated" << ";  " << "sigma_estimated" << std::endl;
  }
};
class scurve_fit_collection {
public:
  scurve_fit_collection(const char* outputFileName);
  void set_mpv_start_low_cut(double mpv_start_low_cut);
  void set_number_of_min_events(int N_min);
  void setTree(TTree* tree);
  void setBeamRuns(TTree* noise);
  void setNoiseRun(TTree* hits);
  void setStartValues();
  void check();


  void setStartMPV(double mpv);

  ////////////////PROCESS SCURVE EFFICIENCY VS THRESHOLD
  void processTotal(const char * total_name);

  /////////PROCESS SINGLE STRIP EFF VS THRESHOLD
  void processStrip(const char* stripName, const std::vector<unsigned>& ignoreIndex);

  TTree *m_tree = nullptr;
  double m_start_amp = 1,
    m_start_mean = 8.79884,
    m_start_gaus_sigma = 30,
    m_start_landau_sigma = 0.386907,
    m_gaus_sigma_lower_boundary = 10,
    m_x_pos = 0,
    m_mpv_start_low_cut = 280;

   int  m_Number_of_min_events= 1000;
  output_data m_out_data;
  bool m_use_total_efficiency_fit_as_start = true;

  TGraphErrors g;
  TGraphErrors pluto;
  TCanvas m_c;
  landgausFit f;
  std::ofstream m_out;
  std::string m_threshold_Cut;
  std::vector<double> m_axis;
  std::shared_ptr<getScurvesFromTTree> m_tree_sc_curve;
};


#endif // scurve_fit_collection_h__