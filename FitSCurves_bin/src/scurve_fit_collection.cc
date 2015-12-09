#include "scurve_fit_collection.hh"
#include "SCT_helpers.h"
#include "s_DrawOption.h"
#include <algorithm>
#include "Math/Interpolator.h"

double sum(const std::vector<double> & vec) {
  double N_event = 0;
  for (const auto&e : vec) {
    N_event += e;
  }
  return N_event;
}

double interpolate(double xa, double xb, double ya, double yb, double x) {
  
  return ya + (yb - ya) / (xb - xa)*(x - xa);

}

std::pair< std::vector<double>, std::vector<double> > remove_outlier(const std::vector<double> & x, const std::vector<double> & y) {
 std::pair< std::vector<double>, std::vector<double> > ret;
  if (x.size()!=y.size()|| x.size() < 3) {
    std::cout << "remove_outlier" << std::endl;
    return ret;
  }
  double x_ = x[0];
  double y_ = y[0];
  ret.first.push_back(x[0]);
  ret.second.push_back(y[0]);
  for (size_t i = 1; i < x.size()-1;++i)  {
    auto y_in = interpolate(x_, x[i + 1], y_, y[i + 1], x[i]);
    std::cout << TMath::Abs(y_in - y[i]) << std::endl;
    if (TMath::Abs( y_in - y[i] ) < 0.1 ) {
      ret.first.push_back(x[i]);
      ret.second.push_back(y[i]);
      y_ = y[i];
      x_ = x[i];
    }
    ret.first.push_back(x.back());
    ret.second.push_back(y.back());
  }
  return ret;
}

scurve_fit_collection::scurve_fit_collection(const char* outputFileName) : m_out(outputFileName) {
  setStartValues();
  m_out_data.printHeeader(m_out);
}

void scurve_fit_collection::set_mpv_start_low_cut(double mpv_start_low_cut) {
  m_mpv_start_low_cut = mpv_start_low_cut;
}

void scurve_fit_collection::setTree(TTree* tree) {
  m_tree = tree;
  m_tree_sc_curve = std::make_shared<getScurvesFromTTree>(tree);
}

void scurve_fit_collection::setBeamRuns(TTree* hits) {
  m_use_total_efficiency_fit_as_start = true;
  m_start_amp = 1;
  if (!m_use_total_efficiency_fit_as_start) {
    m_start_mean = 110;
  }
  m_start_gaus_sigma = 52;
  m_start_landau_sigma = 10;
  m_gaus_sigma_lower_boundary = 20;
  m_mpv_start_low_cut = 280;
  m_threshold_Cut = "Threshold>=25";
  setThreshold_linear_trans(1, 0);
  setTree(hits);
  setStartValues();
}

void scurve_fit_collection::setNoiseRun(TTree* hits) {
  m_start_amp = 1;
  m_start_mean = 8;
  m_start_gaus_sigma = 12;
  m_start_landau_sigma = 1.2;
  m_gaus_sigma_lower_boundary = 1;
  m_use_total_efficiency_fit_as_start = false;
  m_mpv_start_low_cut = 10;
  m_threshold_Cut = "Threshold>=0";
  setThreshold_linear_trans(2.772, 12.39);
  setTree(hits);
  setStartValues();
}

void scurve_fit_collection::setStartValues() {
  f.setFitRange(0, 600);
  f.setLimits_Amplitude(0.7, 1.0);
  f.setStartAmplitude(m_start_amp);
  f.setStartGaussSigma(m_start_gaus_sigma);
  f.setStartLandauMean(m_start_mean);
  f.setStartLandauSigma(m_start_landau_sigma);
  f.setLimits_GaussSigma(m_gaus_sigma_lower_boundary, 1000);
}

void scurve_fit_collection::check() {
  if (f.getLandauSigma() == m_gaus_sigma_lower_boundary) {
    std::cout << "fit hit boundary \n";
  }
}


void scurve_fit_collection::setStartMPV(double mpv) {
  m_start_mean = mpv;
  m_use_total_efficiency_fit_as_start = false;
}

void scurve_fit_collection::processTotal(const char * total_name) {
  
  TCanvas canv("canv", "", 500, 500);
  setStartValues();
  SCT_helpers::DrawTTree(
    m_tree,
    S_DrawOption()
    .output_object(&g)
    .draw_axis(total_name)
    .opt_bar()
    .cut(m_threshold_Cut.c_str())
    );

  f(&g);
  f.DrawfitFunction();
  //push_to_file();
  if (m_use_total_efficiency_fit_as_start) {
    m_start_amp = f.getAmplitude();

    m_start_mean = f.getLandauMostProbable();
  }
  m_out_data = output_data();

  m_out_data.chiSquare = f.getChiSqare();
  m_out_data.gauss_sigma = f.getGaussSigma();
  m_out_data.landau_sigma = f.getLandauSigma();
  m_out_data.max_efficiency = f.getAmplitude();
  m_out_data.mpv_fit = f.getLandauMostProbable();
  m_out_data.print(m_out);
  m_start_gaus_sigma = f.getGaussSigma();
  m_start_landau_sigma = f.getLandauSigma();
  canv.SaveAs("total.png");
  canv.Close();
}

double get_startPoint(const std::vector<double> & thr, const std::vector<double> & occ) {
  const double occ_begin = occ[0];
  for (size_t i = 0; i < occ.size(); ++i) {
    if (occ.at(i) < 0.5 * occ_begin) {
      return interpolate(occ.at(i - 1), occ.at(i), thr.at(i - 1), thr.at(i), 0.5 * occ_begin);
      
    }

  }
  return -1;
}
double get_startSigma(const std::vector<double> & thr, const std::vector<double> & occ) {

  const double occ_begin = occ[0];
  double max_thr = 0, min_thr = 0;
  for (size_t i = 0; i < occ.size(); ++i) {
    if (occ.at(i) <(1 - 0.158650) * occ_begin) {
      min_thr = interpolate(occ.at(i - 1), occ.at(i), thr.at(i - 1), thr.at(i), (1 - 0.158650) * occ_begin);
      break;
    }

  }
  for (size_t i = 0; i < occ.size(); ++i) {
    if (occ.at(i) < 0.158650 * occ_begin) {
      
      max_thr = interpolate(occ.at(i - 1), occ.at(i), thr.at(i - 1), thr.at(i), 0.158650 * occ_begin);
      
      break;
    }

  }
  return max_thr - min_thr;
  
  

}

class scurve_buffer{
public:
  scurve_buffer(output_data& outbuffer_, const std::vector<double> &x_, const std::vector<double> &y_, const std::vector<double> &N_) :outbuffer(outbuffer_), x(x_), y(y_),N(N_) {}
  output_data& outbuffer;
  const std::vector<double> &x;
  const std::vector<double> &y;
  const std::vector<double> &N;
};

template <typename T, typename... Args> 
T __check(Args&&... args);


template <typename T>
auto setBuffer(scurve_buffer&  inBuffer, T& task_) -> decltype(__check<T&>(task_.buffer = &inBuffer)) {
  task_.buffer = &inBuffer;
  return task_;
}
class loop_over_scurves {
public:
  loop_over_scurves(getScurvesFromTTree& scurves) :m_scurves(scurves) {}
  getScurvesFromTTree& m_scurves;
};

TASK_DEFINITION(const loop_over_scurves& task_) {
  const auto max_index = task_.m_scurves.size();

  output_data outbuffer;
  for (size_t i = 0; i < max_index; ++i) {
    outbuffer = output_data();
    task_.m_scurves.set_pos(i);
    outbuffer.x_pos = task_.m_scurves.getX();
    scurve_buffer buf(outbuffer, *task_.m_scurves.get_Threshold(), *task_.m_scurves.get_Ocuupancy(), *task_.m_scurves.get_total_hits());
    LOOP_TASK_NEXT(buf);
    std::cout << outbuffer.x_pos << std::endl;
  }
  RETURN_OK;
}


class print__ {
public:
  print__(std::ostream& out):m_out(out) {}
  std::ostream& m_out;
  scurve_buffer* buffer=nullptr;
};

TASK_DEFINITION(print__& task_) {
  auto ret = RUN_TASK(*task_.buffer);
  task_.buffer->outbuffer.print(task_.m_out);
  return ret;
}


class filter_more_than{
public:
  filter_more_than(int minNumber) :m_minElements(minNumber) {}
  scurve_buffer* buffer = nullptr;
  int m_minElements = 0;
};


TASK_DEFINITION(const filter_more_than& task_) {
  double N_event = sum(task_.buffer->N);
  if (N_event<task_.m_minElements) {
    return returnTypes::return_skip;
  }
  RETURN_NEXT_TASK(*task_.buffer);

}

class filter_max_occupancy_higher_than{
public:
  filter_max_occupancy_higher_than(double max_occ) :m_max_occupancy(max_occ){}
  scurve_buffer* buffer = nullptr;
  double m_max_occupancy;
};

TASK_DEFINITION(const filter_max_occupancy_higher_than& task_) {
  
  auto& occ = task_.buffer->y;
  if (occ.empty()) {
    return return_skip;
  }
  auto d = std::max_element(occ.begin(), occ.end());
  if (*d < task_.m_max_occupancy) {
    return return_skip;
  }
  RETURN_NEXT_TASK(*task_.buffer);
}

class startValue {
public:
  startValue(double mpv_, double sigma_):mpv(mpv_),sigma(sigma_) {

  }
  scurve_buffer* buffer = nullptr;
  double mpv = 0,sigma = 0;
};

TASK_DEFINITION(const startValue&  task_) {

  task_.buffer->outbuffer.sigma_estimated = task_.sigma;
  task_.buffer->outbuffer.mpv_estimated = task_.mpv;
  RETURN_NEXT_TASK(*task_.buffer);
}


class estimateStartValue {
public:
  estimateStartValue(double lower_cut_):lower_cut(lower_cut_){

  }
  scurve_buffer* buffer = nullptr;
  double lower_cut = 0;
};

TASK_DEFINITION(const estimateStartValue& task_) {
  
  const auto mpv_start = get_startPoint(task_.buffer->x, task_.buffer->y);
  if (mpv_start < task_.lower_cut) {
    return return_skip;
  }

  auto sigma = get_startSigma(task_.buffer->x, task_.buffer->y);
  task_.buffer->outbuffer.sigma_estimated = sigma;
  task_.buffer->outbuffer.mpv_estimated = mpv_start;
  RETURN_NEXT_TASK(*task_.buffer);
}

class fit_scurve {
public:
  fit_scurve(landgausFit& f_):f(f_) {

  }
  scurve_buffer* buffer = nullptr;
  landgausFit& f;
};

TASK_DEFINITION(const fit_scurve& task_) {

  TGraph g(task_.buffer->x.size(), task_.buffer->x.data(), task_.buffer->y.data());
  auto & f = task_.f;
  f(&g);

  auto& m_out_data = task_.buffer->outbuffer;
  m_out_data.mpv_fit = f.getLandauMostProbable();
  m_out_data.mpv_fit_error = f.getErrorOfLandauMP();
  m_out_data.max_efficiency = f.getAmplitude();
  m_out_data.gauss_sigma = f.getGaussSigma();
  m_out_data.landau_sigma = f.getLandauSigma();
  m_out_data.chiSquare = f.getChiSqare();
  f.DrawfitFunction();
  RETURN_NEXT_TASK(*task_.buffer);
}
class save_fit {
public:
  save_fit(TCanvas& c_) :c(c_) {

  }
  scurve_buffer* buffer = nullptr;
  TCanvas& c;
  int i = 0;
};

TASK_DEFINITION(save_fit&  task_) {
  std::string cname = std::string("scurve") + std::to_string(task_.i++) + std::string(".png");
  task_.c.SaveAs(cname.c_str());


  RETURN_NEXT_TASK(*task_.buffer);
}

class removeOutlier {
public:
  removeOutlier(double max_outlier):max_out(max_outlier) {}
  std::vector<double> x, y, N;
  scurve_buffer* buffer = nullptr;
  double max_out;
 
};


TASK_DEFINITION(removeOutlier&  task_) {
  task_.x.clear();
  task_.y.clear();
  task_.N.clear();

  auto & x = task_.buffer->x;
  auto & y = task_.buffer->y;
  auto & N = task_.buffer->N;
  if (x.size() != y.size() || x.size() < 3) {
    std::cout << "remove_outlier" << std::endl;
    RETURN_SKIP;
  }
  double x_ = x[0];
  double y_ = y[0];
  task_.x.push_back(x[0]);
  task_.y.push_back(y[0]);
  task_.N.push_back(N[0]);
  for (size_t i = 1; i < x.size() - 1; ++i) {
    auto y_in = interpolate(x_, x[i + 1], y_, y[i + 1], x[i]);
  //  std::cout << TMath::Abs(y_in - y[i]) << std::endl;
    if (TMath::Abs(y_in - y[i]) <task_.max_out) {
      task_.x.push_back(x[i]);
      task_.y.push_back(y[i]);
      task_.N.push_back(N[i]);
      y_ = y[i];
      x_ = x[i];
    }

  }
  task_.x.push_back(x.back());
  task_.y.push_back(y.back());
  task_.N.push_back(N.back());

  scurve_buffer buffer(task_.buffer->outbuffer, task_.x, task_.y, task_.N);
  RETURN_NEXT_TASK(buffer);
}

class removeElements {
public:
  removeElements(std::vector<unsigned> index_) :index(index_) {

  }
  scurve_buffer* buffer = nullptr;
  std::vector<double> x, y, N;
  std::vector<unsigned> index;
};

TASK_DEFINITION(removeElements& task_) {
  task_.x.clear();
  task_.y.clear();
  task_.N.clear();

  auto & x = task_.buffer->x;
  auto & y = task_.buffer->y;
  auto & N = task_.buffer->N;
  for (size_t i = 0; i < x.size() ; ++i) {
   
    if (std::find(task_.index.begin(), task_.index.end(), i) == task_.index.end()) {
      task_.x.push_back(x[i]);
      task_.y.push_back(y[i]);
      task_.N.push_back(N[i]);
    }

  }

  scurve_buffer buffer(task_.buffer->outbuffer, task_.x, task_.y, task_.N);
  RETURN_NEXT_TASK(buffer);
}

class stop_ {
public:
  stop_() {}
  
};
 template <typename T>
 stop_& setBuffer(T&&, stop_& task_) {
   return task_;
 }
returnTypes runTask(const stop_& task_) {
  RETURN_OK;
}
void scurve_fit_collection::processStrip(const char* stripName, const std::vector<unsigned>& ignoreIndex) {
  
 
  runTask(
    loop_over_scurves(*m_tree_sc_curve),
    print__(m_out),
    filter_more_than(1000),
    filter_max_occupancy_higher_than(0.5),
    estimateStartValue(150),
    removeElements(ignoreIndex),
   // removeOutlier(0.1),
    fit_scurve(f),
    save_fit(m_c),
    stop_()
    );
 

}
