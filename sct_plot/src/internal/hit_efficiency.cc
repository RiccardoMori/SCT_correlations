#include "internal/hit_efficiency.hh"
#include "plane_def.h"
#include "s_DrawOption.h"
#include "s_plot_prob.h"
#include "TH1.h"
#include "SCT_helpers.h"

sct_corr::hit_efficiency::hit_efficiency(
  const sct_corr::plane_def& trueHits, 
  const sct_corr::plane_def& trueHits_with_dut, 
  const s_plot_prob& plot_prob
  ):m_plot_prob(plot_prob) {
  m_DUT_hits = trueHits_with_dut.copy();
  m_trueHits = trueHits.copy();
}

Long64_t sct_corr::hit_efficiency::Draw_true_hits(
  const S_DrawOption& d_option
  ) {
  return m_trueHits->Draw(d_option);

}

Long64_t sct_corr::hit_efficiency::Draw_DUT_hits(const S_DrawOption& d_option) {
  return m_DUT_hits->Draw(d_option);
}

Long64_t sct_corr::hit_efficiency::Draw() {
 return Draw(S_DrawOption().draw_x());
}

TH1D* sct_corr::hit_efficiency::getEfficiency_map() const {
  
  if (!m_efficiency)
  {
    return nullptr;
  }
  return m_efficiency.get();
}

double sct_corr::hit_efficiency::get_efficiency() const {
  return static_cast<double>(m_dut_count) / m_total;
}

Long64_t sct_corr::hit_efficiency::Draw(
  const S_DrawOption& d_option
  ) {
  
  TH1D trueHitsMod("mod_true", "true hits", 400, 0, 400);
  TH1D DUtHitsMod("mod_dut", "DUT hits", 400, 0, 400);
  S_DrawOption local_(d_option);
  m_total = Draw_true_hits(local_.output_object(&trueHitsMod));
  auto m_dut_count = Draw_DUT_hits(local_.output_object(&DUtHitsMod));


  m_efficiency = std::shared_ptr<TH1D>(dynamic_cast<TH1D*>(SCT_helpers::calc_efficiency(&trueHitsMod, &DUtHitsMod)));
  if (!m_efficiency) {
    return -1;
  }
  m_efficiency->SetTitle(necessary_CONVERSION(m_plot_prob.getName()).c_str());
  m_efficiency->Draw(d_option.getOptions());

  return m_dut_count;
}
