#ifndef inStripClusterSize_h__
#define inStripClusterSize_h__
#include <memory>
#include "S_Axis.h"
#include "s_plot_prob.h"
#include "internal/platform.hh"
#include "sct_types.h"
class TH2D;
class TH2;
class TProfile;
class S_Cut;
class s_plot_prob;
class S_DrawOption;
namespace sct_corr {
class plane_def;
class DllExport inStripClusterSize {
public:
  inStripClusterSize(
    const sct_corr::plane_def& trueHits_with_dut,
    const sct_corr::plane_def& sz_data,
    int max_cluster_size,
    axis_def search_axis,
    const sct_type::modulo_t &,
    const s_plot_prob& plot_prob = ""
    );

  Long64_t Draw();
  Long64_t Draw(const S_DrawOption& d_option);

  TProfile* getProfile() const;
  TH2* getHistogram() const; 

private:
  std::shared_ptr<sct_corr::plane_def> m_cluster_size_vs_rel_pos;
  std::shared_ptr<TH2D> m_cluster_size_plot;
  std::shared_ptr<TProfile> m_profile;
  sct_type::modulo_t m_mod;
  axis_def m_search_axis;
  s_plot_prob m_plot_prob;
  int m_max_cluster;
};
}

#endif // inStripClusterSize_h__
