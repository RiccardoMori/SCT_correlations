#ifndef s_plane_h__
#define s_plane_h__

#include "internal/platform.hh"
#include "internal/plane_hit.hh"
#include "Rtypes.h"
#include <memory>
#include "S_Axis.h"

class S_Axis;
class S_treeCollection;


namespace sct_corr{
  class treeCollection;
  class axis_ref;
  class plane;
  class plane_def;
  class rootEventBase;
}



class DllExport S_plane{
public:
  S_plane();


  S_plane(const sct_corr::plane_def& plane_def, sct_corr::treeCollection* hits);
  S_plane(const sct_corr::plane_def& plane_def, sct_corr::rootEventBase* hits);

  sct_corr::plane* getPlane() const;





  const sct_type::collectionName_t& getName() const;
  const sct_type::ID_t& getID() const;
  bool next();



  const sct_corr::axis_ref* getX() const;
  const sct_corr::axis_ref* getY() const;
  const sct_corr::axis_ref* getAxis(axis_def ax) const;

  std::shared_ptr<sct_corr::plane_def> m_plane_def;
  std::shared_ptr<sct_corr::plane> m_plane;


};






#endif // s_plane_h__
