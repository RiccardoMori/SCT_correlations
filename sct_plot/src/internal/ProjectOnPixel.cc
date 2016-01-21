#include "internal/plot_hit2d.hh"
#include "internal/plane.hh"

namespace sct_corr{
  class ProjectOnPixel :public plot_hit2d{
  public:
    ProjectOnPixel(const s_plot_prob& = "");

    virtual void processHit(double x, double y) override;

    double m_x_pixelsize = 0.074, m_ypixelsize = 100000;
    virtual s_plane_collection getOutputcollection() override;
    virtual const char* getType() const override;
  };


  ProjectOnPixel::ProjectOnPixel(const s_plot_prob& plot_prob) : plot_hit2d(plot_prob)
  {

  }

  void ProjectOnPixel::processHit(double x, double y)
  {
    pushHit((int)(m_x->get() / m_x_pixelsize * 100) % 100, (int)(m_y->get() / m_ypixelsize * 100) % 100);
  }

  s_plane_collection ProjectOnPixel::getOutputcollection()
  {
    s_plane_collection ret;
    ret.m_planes.push_back(std::make_pair(std::string("ProjectOnPixel"), plane_def(getOutputName(), sct_type::ID_t(0))));
    return ret;
  }

  const char* ProjectOnPixel::getType() const
  {
    return "ProjectOnPixel___";
  }
}
sct_corr::S_plot sct_plot::projectOnPixel(const s_plot_prob& plot_prob)
{
  return sct_corr::S_plot(new sct_corr::ProjectOnPixel(plot_prob));
}