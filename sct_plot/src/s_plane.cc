#include <exception>
#include <iostream>

#include "sct_plots_internal.h"
#include "internal/plane.hh"
#include "S_Axis.h"

#include "treeCollection.h"
#include "internal/plane_hit.hh"
#include "s_DrawOption.h"
#include "s_plane.h"
#include "plane_def.h"
#include "geometry/setup_description.hh"
#include "internal/exceptions.hh"


S_plane::S_plane() : m_plane(nullptr) {
  SCT_THROW("unsupported default constructor do not use");
}




S_plane::S_plane(const sct_corr::plane_def& plane_def, sct_corr::treeCollection* hits) : S_plane(plane_def,&hits->m_rootBuffer) {
  
}





S_plane::S_plane(const sct_corr::plane_def& plane_def, sct_corr::rootEventBase* hits) : m_plane_def(plane_def.copy()) {
  m_plane = hits->createPlane(m_plane_def->getID());
}

sct_corr::plane* S_plane::getPlane() const {
  return m_plane.get();
  
}

const sct_type::collectionName_t& S_plane::getName() const {
  return m_plane_def->getName();
}
const sct_type::ID_t& S_plane::getID() const {
  return m_plane_def->getID();
}

bool S_plane::next() {

  return m_plane->next();
}




const sct_corr::axis_ref* S_plane::getX() const {
  if (!m_plane) {
    SCT_THROW("plane not set ");
  }
  return m_plane->getAxis(x_axis_def);
}

const sct_corr::axis_ref* S_plane::getY() const {
  if (!m_plane) {
    SCT_THROW("plane not set ");

  }
  return m_plane->getAxis(y_axis_def);
}

const sct_corr::axis_ref* S_plane::getAxis(axis_def ax) const {
  if (!m_plane) {
    SCT_THROW("plane not set ");
  }
  return m_plane->getAxis(ax);
}

namespace sct_corr {
const sct_corr::Xlayer* plane_def::getLayer() const {
  if (m_layer) {
    return m_layer.get();
  }

  return nullptr;
}
std::shared_ptr<sct_corr::plot_collection> plane_def::get_plot() const {
  std::shared_ptr<sct_corr::plot_collection> ret = m_plot.lock();
  if (!ret) {
    SCT_THROW("plot is expired");
  } 
  return ret;
}

std::shared_ptr<plane_def> plane_def::copy() const {
  return std::shared_ptr<plane_def>(new plane_def(*this));
}
void plane_def::set_plot_collection(std::weak_ptr<sct_corr::plot_collection> plot_collection___) {
#ifdef _DEBUG
  if (m_plot.lock()) {
    SCT_THROW("m_plot was already defined ");
  }

#endif // _DEBUG
  m_plot = plot_collection___;
  for (auto&e : m_axis) {
    e.set_plot_collection(plot_collection___);

  }
}
Long64_t plane_def::Draw(const S_DrawOption& opt) const {
  return get_plot()->Draw(*this, opt);
}

plane_def::plane_def(const sct_type::collectionName_t& name, const sct_type::ID_t& ID, const sct_corr::Xlayer* layer_) :m_name(name), m_ID(ID) {
  if (layer_) {
    m_layer = std::unique_ptr<sct_corr::Xlayer>(new sct_corr::Xlayer(*layer_));
  }
  m_axis.emplace_back(name, ID, x_axis_def);
  m_axis.emplace_back(name, ID, y_axis_def);
}
const sct_type::ID_t& plane_def::getID() const {
  return m_ID;
}

const sct_type::collectionName_t& plane_def::getName() const {
  return m_name;
}


S_Axis plane_def::getX_def() const {
  return get_Axis(x_axis_def);
}
S_Axis plane_def::getY_def() const {
  return get_Axis(y_axis_def);
}

S_Axis plane_def::get_Axis(axis_def ax) const {
  for (const auto&e : m_axis) {
    if (e.m_axis == ax) {
      return e;
    }
  }
  SCT_THROW("unknown type");
  return S_Axis();
}

const std::vector<S_Axis>& plane_def::get_axis_defs() const {
  return m_axis;
}

}
std::shared_ptr<sct_corr::plane_def> sct_corr::plane_def_GBL::copy() const {
  return std::shared_ptr<sct_corr::plane_def>(new sct_corr::plane_def_GBL(*this));
}
sct_corr::plane_def_GBL::plane_def_GBL(const sct_type::collectionName_t& name, const sct_type::ID_t& ID, const sct_corr::Xlayer* layer_) :sct_corr::plane_def(name, ID, layer_) {
  m_axis.emplace_back(name, ID, chi2_axis_def);
  m_axis.emplace_back(name, ID, Ndf_axis_def);
  m_axis.emplace_back(name, ID, phi_axis_def);
  m_axis.emplace_back(name, ID, theta_axis_def);
}

S_Axis sct_corr::plane_def_GBL::getChi2_def() const {
  return get_Axis(chi2_axis_def);
}

S_Axis sct_corr::plane_def_GBL::getNdf_def() const {
  return get_Axis(Ndf_axis_def);
}

S_Axis sct_corr::plane_def_GBL::getPhi_def() const {
  return get_Axis(phi_axis_def);
}

S_Axis sct_corr::plane_def_GBL::getTheta_def() const {
  return get_Axis(theta_axis_def);
}


sct_corr::plane_def_Alibava::plane_def_Alibava(const sct_type::collectionName_t& name, const sct_type::ID_t& ID, const sct_corr::Xlayer* layer_) :sct_corr::plane_def(name, ID, layer_) {
  m_axis.emplace_back(name, ID, charge_axis_def);
}

S_Axis sct_corr::plane_def_Alibava::getCharge_def() const {
  return get_Axis(charge_axis_def);
}


std::shared_ptr<sct_corr::plane_def> sct_corr::plane_def_Alibava::copy() const {
  return std::shared_ptr<sct_corr::plane_def>(new sct_corr::plane_def_Alibava(*this));
}



s_plane_collection_correlations::s_plane_collection_correlations() :s_plane_collection() {

}

void s_plane_collection_correlations::setResidual(const sct_corr::plane_def& plane_) {
  push_back("Residual", plane_);
}

void s_plane_collection_correlations::setResidualVsMissing(const sct_corr::plane_def& plane_) {
  push_back("ResidualVsMissing", plane_);
}

void s_plane_collection_correlations::setTotalTrueHits(const sct_corr::plane_def& plane_) {
  push_back("TotalTrueHits", plane_);
}

void s_plane_collection_correlations::setTrueHitsWithDUT(const sct_corr::plane_def& plane_) {
  push_back("TrueHitsWithDUT", plane_);
}

void s_plane_collection_correlations::setDUT_Hits(const sct_corr::plane_def& plane_) {
  push_back("DUT_Hits", plane_);
}

sct_corr::plane_def s_plane_collection_correlations::getResidual() const {
  return getByType("Residual")();
}

sct_corr::plane_def s_plane_collection_correlations::getResidualVSmissing() const {
  return getByType("ResidualVsMissing")();
}

sct_corr::plane_def s_plane_collection_correlations::getTotalTrueHits() const {
  return getByType("TotalTrueHits")();
}

sct_corr::plane_def s_plane_collection_correlations::getTrueHitsWithDUT() const {
  return getByType("TrueHitsWithDUT")();
}

sct_corr::plane_def s_plane_collection_correlations::getDUT_Hits() const {
  return getByType("DUT_Hits")();
}

s_plane_collection_find_closest::s_plane_collection_find_closest() :s_plane_collection() {

}

void s_plane_collection_find_closest::setResidual(const sct_corr::plane_def& plane_) {
  push_back("Residual", plane_);
}

void s_plane_collection_find_closest::setHitOnPlaneA(const sct_corr::plane_def& plane_) {
  push_back("HitOnPlaneA", plane_);
}

void s_plane_collection_find_closest::setHitOnPlaneB(const sct_corr::plane_def& plane_) {
  push_back("HitOnPlaneB", plane_);
}

sct_corr::plane_def s_plane_collection_find_closest::getResidual() const {
  return getByType("Residual")();
}

sct_corr::plane_def s_plane_collection_find_closest::getHitOnPlaneA() const {
  return getByType("HitOnPlaneA")();
}

sct_corr::plane_def s_plane_collection_find_closest::getHitOnPlaneB() const {
  return getByType("HitOnPlaneB")();
}

void s_plane_collection::clear() {
  m_planes.clear();
}
sct_corr::plane_def sct_corr::error_plane_def() {
  SCT_THROW("use of error plane");
  return sct_corr::plane_def(sct_type::collectionName_t("error"), sct_type::ID_t(0));
}
