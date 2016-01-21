#include "sct_events/rootEvent_X_Y_hits.hh"
#include "S_Axis.h"
#include "internal/planeX_Y_hits.hh"
#include "internal/plane_hit.hh"


namespace sct_corr{
  std::vector<sct_type::AxesName_t> get_rootEvent_X_Y_hitsList(){
    return {
      getIDString(),
      axis2String(x_axis_def),
      axis2String(y_axis_def) };
   
  }
  rootEvent_X_Y_hits::rootEvent_X_Y_hits(const sct_type::collectionName_t& name) :rootEventBase(name, get_rootEvent_X_Y_hitsList()),
    m_x(getData(x_axis_def)),
    m_y(getData(y_axis_def)), 
    m_id(getData(getIDString()))
  {

  }



  rootEvent_X_Y_hits::rootEvent_X_Y_hits()
  {

  }

  rootEvent_X_Y_hits::~rootEvent_X_Y_hits()
  {

  }

  std::shared_ptr<plane> rootEvent_X_Y_hits::createPlane(const sct_type::ID_t& ID)
  {
    return std::dynamic_pointer_cast<plane>(std::make_shared<planeX_Y>( plane_struct(ID, this)));
  }



  void rootEvent_X_Y_hits::push_Hit(const plane_hit& h, double ID/*=0*/)
  {
    m_x->push_back(h.x);
    m_y->push_back(h.y);
    m_id->push_back(ID);
    get_event_nr()++;
  }

}