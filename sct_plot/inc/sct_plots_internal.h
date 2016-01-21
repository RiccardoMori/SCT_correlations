#ifndef sct_plots_internal_h__
#define sct_plots_internal_h__
#include "sct_plots.h"
#include "TMath.h"
#include "internal/plane_hit.hh"

namespace sct_corr{
  plane_hit rotate(const  plane_hit& h, double Angle);





  plot* create_plot(const char* type, const char* name, axis_ref* x, axis_ref* y);

  plot* create_plot(const char* type, const char* name, S_plane* x, S_plane* y);
  plot* create_plot(const S_plot_def& plot_def, S_plane* x, S_plane* y);
  plot* create_plot(const S_plot_def& plot_def, axis_ref* x, axis_ref* y);
}
#endif // sct_plots_internal_h__
