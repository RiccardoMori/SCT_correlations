#ifndef residual_efficienct_h__
#define residual_efficienct_h__

#include "internal/platform.hh"
#include "s_plot_prob.h"
#include <memory>
#include "S_Axis.h"
#include "sct_types.h"


class S_Cut;
class s_plot_prob;
class TH1D;
class S_DrawOption;
namespace sct_corr {
class plane_def;



class DllExport residual_efficienct {
public:
  residual_efficienct(
    const sct_corr::plane_def& trueHits,
    const sct_corr::plane_def& sz_data,
    const S_Cut& cut_,
    const sct_type::stripNr_t& strips,
    axis_def search_axis,
    const s_plot_prob& plot_prob
    );
  residual_efficienct(
    const sct_corr::plane_def& trueHits,
    const sct_corr::plane_def& sz_data,
    const sct_type::stripNr_t&  strips,
    axis_def search_axis,
    const s_plot_prob& plot_prob
    );
  Long64_t Draw_true_hits(const S_DrawOption& d_option);

  Long64_t Draw_DUT_hits(const S_DrawOption& d_option);
  Long64_t Draw();
  Long64_t Draw(const S_DrawOption& d_option);
  TH1D* getEfficiency_map() const;
  TH1D* get_total() const;
  double get_efficiency() const;
private:
  std::shared_ptr<sct_corr::plane_def> m_trueHits;
  std::shared_ptr<sct_corr::plane_def> m_DUT_hits;
  std::shared_ptr<TH1D> m_efficiency,m_total_hits;
  s_plot_prob m_plot_prob;
  Long64_t m_total = -1;
  Long64_t m_dut_count = 1;


};

}


#endif // residual_efficienct_h__
