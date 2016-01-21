#include "internal/plot_hit2d.hh"
#include "internal/plane.hh"
#include <math.h>

namespace sct_corr{
  class moduloHitMap :public plot_hit2d{
  public:
    moduloHitMap(double modulo_param_x,double modulo_y,const s_plot_prob& = "");
    virtual void processHit(double x, double y) override;
    virtual s_plane_collection getOutputcollection();
    virtual const char* getType() const override;
    double m_modulo_x , m_mod_y;
  };



  moduloHitMap::moduloHitMap(double modulo_param_x, double modulo_y, const s_plot_prob& plot_prob/*= ""*/) : plot_hit2d(plot_prob), m_modulo_x(modulo_param_x), m_mod_y(modulo_y)
  {

  }

  void moduloHitMap::processHit(double x, double y)
  {
    auto newx = fmod(x, m_modulo_x);
    auto newy = fmod(y, m_mod_y);
    pushHit(newx,newy);
  }

  s_plane_collection moduloHitMap::getOutputcollection()
  {
    s_plane_collection ret;
    ret.m_planes.push_back(std::make_pair(std::string("hitmap"), plane_def(getOutputName(), sct_type::ID_t(0))));
    return ret;
  }

  const char* moduloHitMap::getType() const
  {
    return "moduloHitMap__";
  }
}
sct_corr::S_plot sct_plot::moduloHitMap(double modulo_parameter_x, double m_modulo_y, const s_plot_prob& plot_prob/*= ""*/) {
  return sct_corr::S_plot(new sct_corr::moduloHitMap(modulo_parameter_x, m_modulo_y, plot_prob));
}
