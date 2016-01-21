#include "internal/find_coressponding.hh"
#include "sct_plots.h"


namespace sct_corr{
  plot_find_correspondingX::plot_find_correspondingX(const s_plot_prob& plot_prob) : plotPlaneVsPlane(plot_prob)
  {

  }

  void plot_find_correspondingX::processHit(const plane_hit& p1, const plane_hit& p2)
  {
    if (TMath::Abs(p1.x - p2.x) < 0.1)
    {
      pushHit(p1.x, p1.y);
    }
  }

  s_plane_collection plot_find_correspondingX::getOutputcollection()
  {
    s_plane_collection ret;
    ret.m_planes.push_back(std::make_pair(std::string("plot_find_correspondingX"), plane_def(getOutputName(), sct_type::ID_t(0))));
    return ret;
  }




  const char* plot_find_correspondingX::getType() const
  {
    return  "plot_find_correspondingX___";
  }

  plot_find_correspondingXY::plot_find_correspondingXY(const s_plot_prob& plot_prob) : plotPlaneVsPlane(plot_prob)
  {

  }




  void plot_find_correspondingXY::processHit(const plane_hit&  p1, const plane_hit&  p2) {

    if (TMath::Abs(p1.x - p2.x) < 0.1 && TMath::Abs(p1.y - p2.y) < 0.1)
    {
      pushHit(p1.x, p1.y);
    }

  }

  s_plane_collection plot_find_correspondingXY::getOutputcollection() {

    s_plane_collection ret;
    ret.m_planes.push_back(std::make_pair(std::string("correspondingXY"), plane_def(getOutputName(), sct_type::ID_t(0))));
    return ret;
  }



  const char* plot_find_correspondingXY::getType() const
  {
    return "plot_find_correspondingXY__";
  }
}
sct_corr::S_plot sct_plot::find_correspondingX(const s_plot_prob& plot_prob)
 {
   return sct_corr::S_plot(new sct_corr::plot_find_correspondingX(plot_prob));
 }

sct_corr::S_plot sct_plot::find_correspondingXY(const s_plot_prob& plot_prob)
 {  
   return sct_corr::S_plot(new sct_corr::plot_find_correspondingXY(plot_prob));
 }
