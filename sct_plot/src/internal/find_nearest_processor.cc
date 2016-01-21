#include "internal/PLaneVsPlane_plots.hh"
#include <iostream>
#include "TMath.h"
#include "internal/find_nearest_processor.hh"
#include "sct_plots.h"

namespace sct_corr{


#define nearest_distance_pos 0
#define nearest_plane1_pos 1
#define nearest_plane2_pos 2

find_nearest_processor::find_nearest_processor(double x_cutOff, double y_cutOff, const s_plot_prob& plot_prob) :plotPlaneVsPlane(plot_prob), m_x_cutOff(x_cutOff), m_y_cutOff(y_cutOff)
  {

  }

  void find_nearest_processor::processEventStart()
  {
    dist = plane_hit(0, 0);
    h1 = plane_hit(0, 0);
    h2 = plane_hit(0, 0);
    r = m_c_noHit;
  }

  void find_nearest_processor::processHit(const plane_hit& p1, const plane_hit& p2)
  {


    auto e = plane_hit((p1.x - p2.x), (p1.y - p2.y));
    auto r1 = TMath::Sqrt((e.x)*(e.x) + (e.y)*(e.y));

    if (
      r1 < r
      &&
      TMath::Abs(e.x) < m_x_cutOff
      &&
      TMath::Abs(e.y) < m_y_cutOff)
    {
      r = r1;
      dist = e;
      h1 = p1;
      h2 = p2;
    }

  }

  void find_nearest_processor::processEventEnd()
  {
    if (r < m_c_noHit)
    {
      pushHit(dist.x, dist.y, nearest_distance_pos);
      pushHit(h1.x, h1.y, nearest_plane1_pos);
      pushHit(h2.x, h1.y, nearest_plane2_pos);
    }
    else{
      auto i = "noHIt";
    }
  }

  s_plane_collection find_nearest_processor::getOutputcollection()
  {
    s_plane_collection ret;
    ret.m_planes.push_back(std::make_pair(std::string("nearest_distance"), plane_def(getOutputName(), sct_type::ID_t(nearest_distance_pos))));
    ret.m_planes.push_back(std::make_pair(std::string("nearest_plane1"), plane_def(getOutputName(), sct_type::ID_t(nearest_plane1_pos))));
    ret.m_planes.push_back(std::make_pair(std::string("nearest_plane2"), plane_def(getOutputName(), sct_type::ID_t(nearest_plane2_pos))));
    return ret;
  }




  const char* find_nearest_processor::getType() const
  {
    return "find_nearest__";
  }
}
sct_corr::S_plot sct_plot::find_nearest(Double_t x_cutoff, Double_t y_cutoff, const s_plot_prob& plot_prob)
{
  return sct_corr::S_plot(new sct_corr::find_nearest_processor(x_cutoff, y_cutoff, plot_prob));
}
