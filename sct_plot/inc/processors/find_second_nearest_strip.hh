#ifndef find_secound_nearest_strip_h__
#define find_secound_nearest_strip_h__
#include "plane_def.h"
#include "s_plot_prob.h"
#include "Rtypes.h"
#include "S_Axis.h"


#include "internal/platform.hh"

#include "internal/PLaneVsPlane_plots.hh"

namespace sct_processor {


class DllExport find_second_nearest_strip {
public:
  find_second_nearest_strip(const sct_corr::plane_def& planeA,
                     const sct_corr::plane_def& planeB,
                     axis_def search_axis,
                     Double_t cutOfff,
                     const s_plot_prob& plot_prob_ = "");
  sct_corr::plane_def getResidual() const;
  sct_corr::plane_def getHitOnPlaneA() const;
  sct_corr::plane_def getHitOnPlaneB() const;

  sct_corr::plane_def getSecondResidual() const;
  sct_corr::plane_def getSecondHitOnPlaneA() const;
  sct_corr::plane_def getsecondHitOnPlaneB() const;
private:
  s_plane_collection m_planes;
};
}

#endif // find_secound_nearest_strip_h__
