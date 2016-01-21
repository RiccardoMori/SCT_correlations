#include "internal/plot_hit2d.hh"
#include "internal/plane.hh"
#include <iostream>
#include "s_cuts.h"
namespace sct_corr{
  class cut_x_y :public plot_hit2d{
  public:
    cut_x_y(const S_Cut& cut_, const s_plot_prob& = "");
    virtual void processHit(double x, double y) override;
    virtual s_plane_collection getOutputcollection();

    std::shared_ptr<S_Cut> m_cut;


    virtual const char* getType() const override;
  };



  cut_x_y::cut_x_y(const S_Cut& cut_, const s_plot_prob& plot_prob) :plot_hit2d(plot_prob), m_cut(cut_.copy())
  {

  }

  void cut_x_y::processHit(double x, double y)
  {

    if (m_cut->isOutOfRange(1, x, y))
    {
      return;
    }

   
    pushHit(x, y);
  }

  s_plane_collection cut_x_y::getOutputcollection()
  {
    s_plane_collection ret;
    ret.m_planes.push_back(std::make_pair(std::string("cut_x_y"), plane_def(getOutputName(), sct_type::ID_t(0))));
    return ret;
  }



  const char* cut_x_y::getType() const
  {
    return "cut_x_y__";
  }
}
sct_corr::S_plot sct_plot::cut_x_y(const S_Cut& cut_, const s_plot_prob& plot_prob)
{

  return sct_corr::S_plot(new sct_corr::cut_x_y(cut_, plot_prob));
}


