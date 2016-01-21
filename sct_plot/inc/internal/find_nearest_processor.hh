#ifndef find_nearest_processor_h__
#define find_nearest_processor_h__
#include "internal/PLaneVsPlane_plots.hh"
namespace sct_corr {
class find_nearest_processor :public plotPlaneVsPlane {
public:


  find_nearest_processor(double x_cutOff, double y_cutOff, const s_plot_prob& = "");
  virtual void processEventStart();
  virtual void processHit(const plane_hit&  p1, const plane_hit&  p2);

  virtual void processEventEnd();
  const double m_c_noHit = 100000000000;
  double m_x_cutOff, m_y_cutOff;
  double r;
  plane_hit  dist = plane_hit(0, 0), h1 = plane_hit(0, 0), h2 = plane_hit(0, 0);
  virtual s_plane_collection getOutputcollection();

  virtual const char* getType() const override;
};
}
#endif // find_nearest_processor_h__
