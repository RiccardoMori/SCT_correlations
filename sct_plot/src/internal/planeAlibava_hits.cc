#include "internal/planeAlibava_hit.hh"
#include <iostream>
#include "factory.hh"
namespace sct_corr{

plane::MainType getPlane_Alibava_names() {
  plane::MainType ret;
  ret.m_names.emplace_back(getIDString());
  ret.m_names.emplace_back(axis2String(x_axis_def));
  ret.m_names.emplace_back(axis2String(y_axis_def));
  ret.m_names.emplace_back(axis2String(charge_axis_def));

  return ret;
}

  planeAlibava_hit::planeAlibava_hit(Parameter_ref pl) :plane(pl) {
    m_x = getAxis(x_axis_def);
    m_y = getAxis(y_axis_def);
    m_charge = getAxis(charge_axis_def);
    m_hit = std::make_shared<Alibava_hit>(0, 0, 0);

  }


  std::string planeAlibava_hit::getType() const {
    return "planeAlibava_hits";
  }

  bool planeAlibava_hit::next() const {
    bool ret = m_axis.next();
    if (ret) {
      m_hit->x = m_x->get();
      m_hit->y = m_y->get();
      m_hit->charge = m_charge->get();
    }
    return ret;
  }

  const plane_hit * planeAlibava_hit::getHit() const {
    return dynamic_cast<const plane_hit *>(m_hit.get());
  }

  registerPlane(planeAlibava_hit, getPlane_Alibava_names());
}