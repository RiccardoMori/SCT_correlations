#ifndef hit_efficiency_h__
#define hit_efficiency_h__
#include "Rtypes.h"
#include <memory>
#include "internal/platform.hh"
#include "s_plot_prob.h"



class S_Cut;
class s_plot_prob;
class TH1D;
class S_DrawOption;
namespace sct_corr {
class plane_def;

class DllExport hit_efficiency {
  public:
    hit_efficiency(
      const plane_def& trueHits,
      const plane_def& trueHits_with_dut,
      const s_plot_prob& plot_prob
      );
    Long64_t Draw_true_hits(const S_DrawOption& d_option);

    Long64_t Draw_DUT_hits(const S_DrawOption& d_option);
    Long64_t Draw();
    Long64_t Draw(const S_DrawOption& d_option);
    TH1D* getEfficiency_map() const;
    double get_efficiency() const;
  private:
    std::shared_ptr<plane_def> m_trueHits;
    std::shared_ptr<plane_def> m_DUT_hits;
    std::shared_ptr<TH1D> m_efficiency;
    s_plot_prob m_plot_prob;
    Long64_t m_total = -1;
    Long64_t m_dut_count=1;
};
}

#endif // hit_efficiency_h__
