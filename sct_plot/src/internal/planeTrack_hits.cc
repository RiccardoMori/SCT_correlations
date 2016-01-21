#include "internal/planeTrack_hits.hh"
#include <iostream>
#include "factory.hh"
namespace sct_corr{

  plane::MainType getPlane_track_names() {
    plane::MainType ret;

    ret.m_names.push_back(sct_type::AxesName_t("x"));
    ret.m_names.push_back(sct_type::AxesName_t("y"));
    ret.m_names.push_back(sct_type::AxesName_t("ID"));
    ret.m_names.push_back(sct_type::AxesName_t("chi2"));
    ret.m_names.push_back(sct_type::AxesName_t("ndf"));
    ret.m_names.push_back(sct_type::AxesName_t("phi"));
    ret.m_names.push_back(sct_type::AxesName_t("theta"));
    return ret;
  }
  
  planeTrack_hits::planeTrack_hits(Parameter_ref pl) :plane(pl)
  {
    m_x     = getAxis(x_axis_def);
    m_y     = getAxis(y_axis_def);
    m_chi2  = getAxis(chi2_axis_def);
    m_phi   = getAxis(phi_axis_def);
    m_theta = getAxis(theta_axis_def);
    m_ndf   = getAxis(Ndf_axis_def);
    m_hit = std::make_shared<track_hits>(0, 0, 0, 0, 0,0);
  }




  std::string planeTrack_hits::getType() const {
    return "planeTrack_hits";
  }

  bool planeTrack_hits::next() const
  {
    bool ret = m_axis.next();
    if (ret)
    {
      m_hit->x = m_x->get();
      m_hit->y = m_y->get();
      m_hit->chi2 = m_chi2->get();
      m_hit->phi = m_phi->get();
      m_hit->theta = m_theta->get();
      m_hit->ndf = m_ndf->get();
    }
    return ret;
  }

  const plane_hit * planeTrack_hits::getHit() const
  {
    return dynamic_cast<const plane_hit *>(m_hit.get());
  }
  registerPlane(planeTrack_hits, getPlane_track_names());
}