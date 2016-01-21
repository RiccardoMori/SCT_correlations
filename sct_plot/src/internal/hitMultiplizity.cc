#include "internal/plot_hit2d.hh"
#include "sct_plots.h"
namespace sct_corr{

  class hitMultiplizity :public plot_hit2d{
  public:
    hitMultiplizity(const s_plot_prob& = "");

    virtual void processEventStart();
    virtual void processHit(double x, double y);;
    virtual void processEventEnd();

    Int_t m_counter = 0;
    virtual s_plane_collection getOutputcollection();

    virtual const char* getType() const override;
  };




  hitMultiplizity::hitMultiplizity(const s_plot_prob& plot_prob) : plot_hit2d(plot_prob)
  {

  }

  void hitMultiplizity::processEventStart()
  {
    m_counter = 0;
  }

  void hitMultiplizity::processHit(double x, double y)
  {
    ++m_counter;
  }

  void hitMultiplizity::processEventEnd()
  {
    pushHit(m_counter, 0);
  }

  s_plane_collection hitMultiplizity::getOutputcollection()
  {
    s_plane_collection ret;
    ret.m_planes.push_back(std::make_pair(std::string("HitMultiplicity"), plane_def(getOutputName(), sct_type::ID_t(0))));
    return ret;
  }

  const char* hitMultiplizity::getType() const
  {
    return "hitMultiplizity___";
  }
}
sct_corr::S_plot sct_plot::hitMultiplizity(const s_plot_prob& plot_prob)
{
  return sct_corr::S_plot(new sct_corr::hitMultiplizity(plot_prob));
}
