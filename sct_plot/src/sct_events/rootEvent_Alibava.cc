#include "sct_events/rootEvent_Alibava.hh"
#include "internal/planeX_Y_hits.hh"

namespace sct_corr {
std::vector<std::string> get_rootEvent_Alibava() {
  std::vector<std::string> ret;
  ret.emplace_back(getIDString());
  ret.emplace_back(axis2String(x_axis_def));
  ret.emplace_back(axis2String(y_axis_def));
  ret.emplace_back(axis2String(charge_axis_def));
  return ret;
}

rootEvent_Alibava::rootEvent_Alibava(const char* name):rootEventBase(name,get_rootEvent_Alibava()),
m_x(getData(x_axis_def)),
m_y(getData(y_axis_def)),
m_id(getData(getIDString())),
m_charge(getData(charge_axis_def))
{

}

rootEvent_Alibava::rootEvent_Alibava() {

}

rootEvent_Alibava::~rootEvent_Alibava() {

}

void rootEvent_Alibava::push_Hit(double x, double y) {
  push_Hit(
    Alibava_hit(
    x,
    y,
    1
    )
    );
}

void rootEvent_Alibava::push_Hit(double x, double y, double ID) {
  push_Hit(
    Alibava_hit(
    x,
    y,
    1
    ),
    ID
    );
}

void rootEvent_Alibava::push_Hit(const plane_hit& h, double ID/*=0*/) {
  push_Hit(
    Alibava_hit(
    h.x,
    h.y,
    1
    ),
    ID
  );

}

void rootEvent_Alibava::push_Hit(const Alibava_hit& h, double ID /*= 0*/) {
  m_x->push_back(h.x);
  m_y->push_back(h.y);
  m_id->push_back(ID);
  m_charge->push_back(h.charge);
  get_event_nr()++;
}

std::shared_ptr<plane> sct_corr::rootEvent_Alibava::createPlane(double ID) {
  return std::dynamic_pointer_cast<plane>(std::make_shared<planeX_Y>(plane_struct(ID, this)));
}

}