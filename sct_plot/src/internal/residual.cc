#include "internal/plot_corr2d.hh"
#include "treeCollection.h"

namespace sct_corr{
  class residual : public plot_corr2d {
  public:
    residual(const s_plot_prob& = "");
    virtual void processHit(double x, double y) override;


    virtual s_plane_collection getOutputcollection();

    virtual const char* getType() const override;
  private:
    int m_index = 0;
  };


  residual::residual(const s_plot_prob& plot_prob) : plot_corr2d(plot_prob)
  {

  }

  void residual::processHit(double x, double y)
  {
    pushHit(x - y, get_current_event());
  }


  s_plane_collection residual::getOutputcollection()
  {
    s_plane_collection ret;
    ret.m_planes.push_back(std::make_pair(std::string("residual"), plane_def(getOutputName(), sct_type::ID_t(0))));
    return ret;
  }

  const char* residual::getType() const
  {
    return "residual___";
  }
}
sct_corr::S_plot sct_plot::residual(const s_plot_prob& plot_prob)
{
  return sct_corr::S_plot(new sct_corr::residual(plot_prob));
}
