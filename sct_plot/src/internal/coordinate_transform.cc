#include "internal/plot_hit2d.hh"
#include <iostream>
#include "sct_plots.h"

namespace sct_corr{
  class coordinate_transform : public plot_hit2d {
  public:
    coordinate_transform(Double_t x_slope, Double_t x_offset, Double_t y_slope, Double_t y_offset, const s_plot_prob& = "");
    virtual void processHit(double x, double y) override;
    double m_x_slope = 0;
    double m_y_slope = 0;
    double m_x_offset = 0;
    double m_y_offset = 0;

    virtual s_plane_collection getOutputcollection();

    virtual const char* getType() const override;
  };





  coordinate_transform::coordinate_transform(Double_t x_slope, Double_t x_offset, Double_t y_slope, Double_t y_offset, const s_plot_prob& plot_prob)
    :plot_hit2d(plot_prob),
    m_x_slope(x_slope),
    m_y_slope(y_slope),
    m_x_offset(x_offset),
    m_y_offset(y_offset)

  {

  }

  void coordinate_transform::processHit(double x, double y)
  {
    pushHit(x*m_x_slope + m_x_offset, y*m_y_slope + m_y_offset);
  }

  s_plane_collection coordinate_transform::getOutputcollection()
  {
    s_plane_collection ret;
    ret.m_planes.push_back(std::make_pair(std::string("coordinate_transform"), plane_def(getOutputName(), sct_type::ID_t(0))));
    return ret;
  }

  const char* coordinate_transform::getType() const
  {
    return "coordinate_transform____";
  }
}
sct_corr::S_plot sct_plot::coordinate_transform(Double_t x_slope, Double_t x_offset, Double_t y_slope, Double_t y_offset, const s_plot_prob& plot_prob)
{

  return sct_corr::S_plot(new sct_corr::coordinate_transform(x_slope, x_offset, y_slope, y_offset, plot_prob));
}

sct_corr::S_plot sct_plot::coordinate_transform_move(Double_t x_offset, Double_t y_offset, const s_plot_prob& plot_prob)
{

  return   sct_plot::coordinate_transform(1, x_offset, 1, y_offset, plot_prob);
}