#include "internal/PLaneVsPlane_plots.hh"
#include "sct_plots.h"

namespace sct_corr{
  class plot_a_if_b_has_a_hit :public plotPlaneVsPlane{
  public:
    plot_a_if_b_has_a_hit(const s_plot_prob& = "");
    virtual void processHit(const plane_hit&  p1, const plane_hit&  p2) override;
    virtual s_plane_collection getOutputcollection();
    virtual const char* getType() const override;

  };


  plot_a_if_b_has_a_hit::plot_a_if_b_has_a_hit(const s_plot_prob& plot_prob) : plotPlaneVsPlane(plot_prob)
  {

  }

  void plot_a_if_b_has_a_hit::processHit(const plane_hit& p1, const plane_hit& p2)
  {
    pushHit(p1.x, p1.y, 0);
    pushHit(p2.x, p2.y, 1);
  }

  s_plane_collection plot_a_if_b_has_a_hit::getOutputcollection()
  {
    s_plane_collection ret;
    ret.m_planes.push_back(std::make_pair(std::string("A_if_B_has_a_hit_get_A"), plane_def(getOutputName(), sct_type::ID_t(0))));
    ret.m_planes.push_back(std::make_pair(std::string("A_if_B_has_a_hit_get_B"), plane_def(getOutputName(), sct_type::ID_t(1))));
    return ret;
  }


  const char* plot_a_if_b_has_a_hit::getType() const
  {
    return "plot_a_if_b_has_a_hit__";
  }

}
sct_corr::S_plot sct_plot::A_if_B(const s_plot_prob& plot_prob)
  {
    return sct_corr::S_plot(new sct_corr::plot_a_if_b_has_a_hit(plot_prob));
  }