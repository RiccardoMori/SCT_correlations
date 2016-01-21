#include "internal/PLaneVsPlane_plots.hh"
#include <iostream>
#include "TMath.h"
#include "plane_def.h"
#include <math.h>
#include "sct_plots.h"



namespace sct_corr{

  class modulo_find_nearest_strip :public plotPlaneVsPlane{
  public:
    modulo_find_nearest_strip(axis_def search_axis, Double_t cutOfff, Double_t modulo_param, const s_plot_prob& = "");
    virtual void processEventStart();
    virtual void processHit(const plane_hit&  p1, const plane_hit&  p2);

    virtual void processEventEnd();
    double m_cutOff,m_mod;
    plane_hit dist = plane_hit(0, 0), h1 = plane_hit(0, 0), h2 = plane_hit(0, 0);
    double r = 10000000;
    axis_def m_axis;
    virtual s_plane_collection getOutputcollection();
    virtual const char* getType() const override;
    plane_hit modulo(const plane_hit&  p1);
  };




  modulo_find_nearest_strip::modulo_find_nearest_strip(axis_def search_axis, Double_t cutOfff, Double_t modulo_param, const s_plot_prob& plot_prob/*= ""*/)
    :plotPlaneVsPlane(plot_prob), 
    m_axis(search_axis), 
    m_cutOff(cutOfff),
    m_mod(modulo_param)
  {

  }

  void modulo_find_nearest_strip::processEventStart()
  {

    dist = plane_hit(0, 0);
    h1 = plane_hit(0, 0);
    h2 = plane_hit(0, 0);
    r = m_cutOff;
  }

  void modulo_find_nearest_strip::processHit(const plane_hit& p1, const plane_hit& p2)
  {
    plane_hit e(p1.x - p2.x, p1.y - p2.y);
    double r1 = 0;
    if (m_axis == x_axis_def)
    {
      r1 = TMath::Abs(e.x);
    }
    else if (m_axis == y_axis_def) {
      r1 = TMath::Abs(e.y);
    }
    else{
      std::cout << "unknown axis" << std::endl;
    }
    if (r1 > 0 && r1 < r)
    {
      r = r1;
      dist = e;
      h1 = modulo(p1);
      h2 = modulo(p2);
    }
    if (r1 == 0)
    {
      std::cout << "[find_nearest_strip] error r1== 0" << std::endl;
    }
  }

  void modulo_find_nearest_strip::processEventEnd()
  {
    if (r < m_cutOff)
    {

      pushHit(dist.x, dist.y, 0);
      pushHit(h1.x, h1.y, 1);
      pushHit(h2.x, h2.y, 2);
    }
  }

  s_plane_collection modulo_find_nearest_strip::getOutputcollection()
  {
    s_plane_collection ret;
    ret.m_planes.push_back(std::make_pair(std::string("nearest_strip_distance"), plane_def(getOutputName(), sct_type::ID_t(0))));
    ret.m_planes.push_back(std::make_pair(std::string("nearest_strip_plane1"), plane_def(getOutputName(), sct_type::ID_t(1))));
    ret.m_planes.push_back(std::make_pair(std::string("nearest_strip_plane2"), plane_def(getOutputName(), sct_type::ID_t(2))));
    return ret;
  }


  const char* modulo_find_nearest_strip::getType() const
  {
    return "modulo_find_nearest_strip___";
  }

  sct_corr::plane_hit modulo_find_nearest_strip::modulo(const plane_hit& p1) {
    
    if (m_axis= x_axis_def)
    {
      auto x_mod = fmod(p1.x,m_mod);
      return sct_corr::plane_hit ( x_mod,p1.y);
    }

    if (m_axis==y_axis_def)
    {
      auto y_mod = fmod(p1.y,m_mod);
      return sct_corr::plane_hit(p1.x,y_mod);
    }
  
    return  sct_corr::plane_hit(0, 0);
  }

}

sct_corr::S_plot sct_plot::modulo_find_nearest_strip(axis_def search_axis, double modulo_param, Double_t cutOfff /*= 100000*/, const s_plot_prob& plot_prob/*= ""*/) {

  return sct_corr::S_plot(new sct_corr::modulo_find_nearest_strip(search_axis, cutOfff, modulo_param, plot_prob));
}
