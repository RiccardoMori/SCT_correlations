#include "internal/residual_efficienct.hh"
#include "plane_def.h"
#include "internal/PLaneVsPlane_plots.hh"
#include "TMath.h"
#include "s_plot.h"
#include "SCT_helpers.h"
#include "sct_processors.h"
#include "s_cuts.h"
#include "sct_types.h"
#include "internal/exceptions.hh"






namespace sct_corr {
class residual_efficiency_processor :public plotPlaneVsPlane {
public:
  residual_efficiency_processor(axis_def search_axis, int num_of_channel, const s_plot_prob& = "");

  virtual void processHit(const plane_hit&  p1, const plane_hit&  p2);


  virtual void processHit_1(const plane_hit&  p1) override;


  double  make_residual(double true_hit, double dut_hit);



  int m_strips;
  const axis_def m_axis;
  virtual s_plane_collection getOutputcollection();
  virtual const char* getType() const override;
};

residual_efficiency_processor::residual_efficiency_processor(
  axis_def search_axis,
  int num_of_channel,
  const s_plot_prob&  plot_prob/*= ""*/
  )
  :plotPlaneVsPlane(plot_prob),
  m_axis(search_axis),
  m_strips(num_of_channel) {

}




void residual_efficiency_processor::processHit_1(const plane_hit&  p1) {
  for (int i = 0; i < m_strips; ++i) {

    if (m_axis == x_axis_def) {

      auto r = make_residual(p1.x, i);

      pushHit(r, p1.y, 0);

    } else if (m_axis == y_axis_def) {
      auto r = make_residual(p1.y, i);
      pushHit(r, p1.x, 0);

    }
  }
}
void residual_efficiency_processor::processHit(const plane_hit& p1, const plane_hit& p2) {

  if (m_axis == x_axis_def) {

    auto r = make_residual(p1.x, p2.x);

    pushHit(r, p1.y, 1);

  } else if (m_axis == y_axis_def) {
    auto r = make_residual(p1.y, p2.y);
    pushHit(r, p1.x, 1);

  }



}



double  residual_efficiency_processor::make_residual(double true_hit, double dut_hit) {
  return TMath::Floor(dut_hit) - true_hit;
}




s_plane_collection residual_efficiency_processor::getOutputcollection() {
  s_plane_collection ret;
  ret.m_planes.push_back(std::make_pair(std::string("all_residuals"), plane_def(getOutputName(), sct_type::ID_t(0))));
  ret.m_planes.push_back(std::make_pair(std::string("dut_hits"), plane_def(getOutputName(), sct_type::ID_t(1))));
  return ret;
}

const char* residual_efficiency_processor::getType() const {
  return "residual_efficiency_processor";
}



s_plane_collection make_residual_efficiency_processor(
  const plane_def& planeA,
  const plane_def& planeB,
  axis_def search_axis,
  const sct_type::stripNr_t& strips,
  const s_plot_prob& plot_prob_ /*= "" */
  ) {
  auto plA = planeA.get_plot();
  auto plB = planeB.get_plot();

  if (plA.get() != plB.get()) {
    SCT_THROW("referencing to different plot collection");
  }



  auto collection_ = plA->addPlot(
    S_plot(new residual_efficiency_processor(
    search_axis,
    Un_necessary_CONVERSION(strips),
    plot_prob_)),
    planeA,
    planeB
    );

 
  return collection_;
}


residual_efficienct::residual_efficienct(
  const plane_def& trueHits,
  const plane_def& sz_data,
  const sct_type::stripNr_t&  strips,
  axis_def search_axis,
  const s_plot_prob& plot_prob
  ) :m_plot_prob(plot_prob) {


  auto r = make_residual_efficiency_processor(trueHits, sz_data, search_axis, strips, m_plot_prob);
  m_DUT_hits = r.get("dut_hits")().copy();
  m_trueHits = r.get("all_residuals")().copy();




}

residual_efficienct::residual_efficienct(
  const plane_def& trueHits, 
  const plane_def& sz_data, 
  const S_Cut& cut_, 
  const sct_type::stripNr_t&  strips,
  axis_def search_axis, 
  const s_plot_prob& plot_prob
  ) :residual_efficienct(
  processor::cut_x_y(trueHits, cut_, s_plot_prob().doNotSaveToDisk()),
  sz_data,
  strips,
  search_axis,
  plot_prob 
  )
{

}

Long64_t residual_efficienct::Draw_true_hits(const S_DrawOption& d_option) {
  return m_trueHits->Draw(d_option);
}

Long64_t residual_efficienct::Draw_DUT_hits(const S_DrawOption& d_option) {
  return m_DUT_hits->Draw(d_option);
}


Long64_t residual_efficienct::Draw() {
  return Draw(S_DrawOption().draw_x());
}

Long64_t residual_efficienct::Draw(const S_DrawOption& d_option) {
  m_total_hits = std::make_shared<TH1D>("res_total_true_hits", "total true hits", 100, -2, 2);
  TH1D DUtHitsMod("res_DUT_hits", "DUT hits", 100, -2, 2);
  S_DrawOption local_(d_option);
  m_total = Draw_true_hits(local_.output_object(m_total_hits.get()));
  auto m_dut_count = Draw_DUT_hits(local_.output_object(&DUtHitsMod));


  m_efficiency = std::shared_ptr<TH1D>(dynamic_cast<TH1D*>(SCT_helpers::calc_efficiency(m_total_hits.get(), &DUtHitsMod)));
  m_efficiency->SetName("res_effieciency");
  if (!m_efficiency) {
    return -1;
  }
  m_efficiency->SetTitle(necessary_CONVERSION(m_plot_prob.getName()).c_str());
  m_efficiency->Draw(d_option.getOptions());
  m_efficiency->GetXaxis()->SetTitle("residual");
  m_efficiency->GetYaxis()->SetTitle("efficiency");
  
  return m_dut_count;
}

TH1D* residual_efficienct::getEfficiency_map() const {
  return m_efficiency.get();
}

TH1D* residual_efficienct::get_total() const {
  return m_total_hits.get();
}

double residual_efficienct::get_efficiency() const {
  return (double)m_dut_count / m_total;
}

}
