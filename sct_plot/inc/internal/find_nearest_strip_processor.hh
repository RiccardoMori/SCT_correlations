#ifndef find_nearest_strip_processor_h__
#define find_nearest_strip_processor_h__
#include "internal/PLaneVsPlane_plots.hh"
#define nearest_strip_distance_name "nearest_strip_distance"
#define nearest_strip_plane1_name "nearest_strip_plane1"
#define nearest_strip_plane2_name "nearest_strip_plane2"
namespace sct_corr {

class find_nearest_strip_processor :public plotPlaneVsPlane {
public:
  find_nearest_strip_processor(axis_def search_axis, Double_t cutOfff, const s_plot_prob& = "");
  virtual void processEventStart();
  virtual void processHit(const plane_hit&  p1, const plane_hit&  p2);

  virtual void processEventEnd();
  double m_cutOff;
  plane_hit dist = plane_hit(0, 0), h1 = plane_hit(0, 0), h2 = plane_hit(0, 0);
  double r = 10000000;
  axis_def m_axis;
  virtual s_plane_collection getOutputcollection();
  virtual const char* getType() const override;
};

}
#endif // find_nearest_strip_processor_h__
