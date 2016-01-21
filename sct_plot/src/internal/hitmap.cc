#include "internal/plot_hit2d.hh"
#include "internal/plane.hh"

namespace sct_corr{
  class hitmap :public plot_hit2d{
  public:
    hitmap(const s_plot_prob& = "");
    virtual void processHit(double x, double y) override;
    virtual s_plane_collection getOutputcollection();
    virtual const char* getType() const override;
  };

  hitmap::hitmap(const s_plot_prob& plot_prob) :plot_hit2d(plot_prob)
  {

  }

  void hitmap::processHit(double x, double y)
  {
    pushHit(x, y);
  }

  s_plane_collection hitmap::getOutputcollection()
  {
    s_plane_collection ret;
    ret.m_planes.push_back(std::make_pair(std::string("hitmap"), plane_def(getOutputName(), sct_type::ID_t(0))));
    return ret;
  }

  const char* hitmap::getType() const
  {
    return "hitmap__";
  }
}
sct_corr::S_plot sct_plot::hitmap(const s_plot_prob& plot_prob)
{
  return sct_corr::S_plot(new sct_corr::hitmap(plot_prob));
}
