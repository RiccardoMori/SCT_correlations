#ifndef find_nearest_strip_h__
#define find_nearest_strip_h__

#include "plane_def.h"
#include "s_plot_prob.h"
#include "S_Axis.h"
#include "Rtypes.h"

namespace sct_processor {

class DllExport find_nearest_strip {
public:
  find_nearest_strip(const sct_corr::plane_def& planeA,
                     const sct_corr::plane_def& planeB,
                     axis_def search_axis,
                     Double_t cutOfff,
                     const s_plot_prob& plot_prob_ = "" );
  sct_corr::plane_def getResidual() const;
  sct_corr::plane_def getHitOnPlaneA() const;
  sct_corr::plane_def getHitOnPlaneB() const;

private:
  s_plane_collection_find_closest m_planes;
};


}
#endif // find_nearest_strip_h__
