#ifndef modulo_h__
#define modulo_h__
#include "plane_def.h"
#include "sct_types.h"
#include "s_plot_prob.h"
#include "S_Axis.h"
#include "internal/platform.hh"

namespace sct_processor {

class DllExport modulo {
public:
  modulo(const sct_corr::plane_def& planeA,
         const sct_type::modulo_t& modulo_value,
         axis_def mod_axis,
         const s_plot_prob& plot_prob_ = ""
         );
  sct_corr::plane_def getModulo() const;
private:
  s_plane_collection m_planes;
};



}

#endif // modulo_h__
