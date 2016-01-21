#include "internal/PLaneVsPlane_plots.hh"
#include <iostream>
#include "TMath.h"
#include "plane_def.h"
#include "internal/find_nearest_strip_processor.hh"
#include "sct_plots.h"


namespace sct_corr{


  find_nearest_strip_processor::find_nearest_strip_processor(axis_def search_axis, Double_t cutOfff, const s_plot_prob& plot_prob) :plotPlaneVsPlane(plot_prob), m_axis(search_axis), m_cutOff(cutOfff)
  {

  }

  void find_nearest_strip_processor::processEventStart()
  {

    dist = plane_hit(0, 0);
    h1 = plane_hit(0, 0);
    h2 = plane_hit(0, 0);
    r = m_cutOff;
  }

  void find_nearest_strip_processor::processHit(const plane_hit& p1, const plane_hit& p2)
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
      h1 = p1;
      h2 = p2;
    }
    if (r1 == 0)
    {
      std::cout << "[find_nearest_strip] error r1== 0" << std::endl;
    }
  }

  void find_nearest_strip_processor::processEventEnd()
  {
    if (r < m_cutOff)
    {

      pushHit(dist.x, dist.y, 0);
      pushHit(h1.x, h1.y, 1);
      pushHit(h2.x, h2.y, 2);
    }
  }

  s_plane_collection find_nearest_strip_processor::getOutputcollection()
  {
    s_plane_collection ret;
    ret.m_planes.push_back(std::make_pair(std::string(nearest_strip_distance_name), plane_def(getOutputName(), sct_type::ID_t(0))));
    ret.m_planes.push_back(std::make_pair(std::string(nearest_strip_plane1_name), plane_def(getOutputName(), sct_type::ID_t(1))));
    ret.m_planes.push_back(std::make_pair(std::string(nearest_strip_plane2_name), plane_def(getOutputName(), sct_type::ID_t(2))));
    return ret;
  }


  const char* find_nearest_strip_processor::getType() const
  {
    return "find_nearest_strip___";
  }
}
sct_corr::S_plot sct_plot::find_nearest_strip(axis_def search_axis, Double_t cutOfff /*=100000*/, const s_plot_prob& plot_prob)
{
  return sct_corr::S_plot(new sct_corr::find_nearest_strip_processor(search_axis, cutOfff, plot_prob));
}
