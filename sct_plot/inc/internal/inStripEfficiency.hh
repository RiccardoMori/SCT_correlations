#ifndef inStripEfficiency_h__
#define inStripEfficiency_h__
#include <memory>

#include "S_Axis.h"
#include "Rtypes.h"
#include "internal/platform.hh"
#include "s_plot_prob.h"
#include "sct_types.h"

class S_Cut;

class TH1D;
class S_DrawOption;

namespace sct_corr {
class plane_def;

class DllExport inStripEfficiency {
public:
  inStripEfficiency(
    const sct_corr::plane_def& trueHits,
    const sct_corr::plane_def& trueHits_with_dut,
    const S_Cut& cut_,
    axis_def search_axis,
    sct_type::modulo_t mod_,
    const s_plot_prob& plot_prob
    );
  inStripEfficiency(
    const sct_corr::plane_def& trueHits,
    const sct_corr::plane_def& trueHits_with_dut,
    axis_def search_axis,
    sct_type::modulo_t mod_,
    const s_plot_prob& plot_prob
    );
  Long64_t Draw_true_hits(const S_DrawOption& d_option);

  Long64_t Draw_DUT_hits(const S_DrawOption& d_option);
  Long64_t Draw();
  Long64_t Draw(const S_DrawOption& d_option);

  TH1D* getEfficiency_map() const;
  TH1D* getHits() const;
  double get_total_efficiency() const;
private:
  std::shared_ptr<sct_corr::plane_def> m_trueHits;
  std::shared_ptr<sct_corr::plane_def> m_DUT_hits;
  std::shared_ptr<TH1D> m_efficiency,m_total_hits;
  sct_type::modulo_t m_mod;
  axis_def m_search_axis;
  s_plot_prob m_plot_prob;
  double m_efficiency_total;
};
}
#endif // inStripEfficiency_h__
