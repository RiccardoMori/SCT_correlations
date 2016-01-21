#include "processors/find_second_nearest_strip.hh"





#include "internal/PLaneVsPlane_plots.hh"
#include "S_Axis.h"
#include <iostream>
#define nearest_strip_distance_name "nearest_strip_distance"
#define nearest_strip_plane1_name "nearest_strip_plane1"
#define nearest_strip_plane2_name "nearest_strip_plane2"

#define second_nearest_strip_distance_name "second_nearest_strip_distance"
#define second_nearest_strip_plane1_name "second_nearest_strip_plane1"
#define second_nearest_strip_plane2_name "second_nearest_strip_plane2"

#define MAX_R_VALUE  100000000000

namespace sct_corr {

class find_second_nearest_strip_processor :public plotPlaneVsPlane {
public:
  find_second_nearest_strip_processor(axis_def search_axis, Double_t cutOfff, const s_plot_prob& = "");
  virtual void processEventStart();
  virtual void processHit(const plane_hit&  p1, const plane_hit&  p2);

  virtual void processEventEnd();
  double m_cutOff;
  class closest {
  public:
  plane_hit dist = plane_hit(0, 0), h1 = plane_hit(0, 0), h2 = plane_hit(0, 0);
  double r = 10000000;
  } m_closest,m_second;
  axis_def m_axis;
  double m_min = 10000;
  virtual s_plane_collection getOutputcollection();
  virtual const char* getType() const override;
};

find_second_nearest_strip_processor::find_second_nearest_strip_processor(axis_def search_axis, Double_t cutOfff, const s_plot_prob& plot_prob) :plotPlaneVsPlane(plot_prob), m_axis(search_axis), m_cutOff(cutOfff) {

}

void find_second_nearest_strip_processor::processEventStart() {

  m_closest = closest();
  m_second = closest();
}

void find_second_nearest_strip_processor::processHit(const plane_hit& p1, const plane_hit& p2) {
  plane_hit e(p1.x - p2.x, p1.y - p2.y);
  double r1 = MAX_R_VALUE;
  if (m_axis == x_axis_def) {
    r1 = TMath::Abs(e.x);
    m_min = TMath::Min(m_min, e.x);
  } else if (m_axis == y_axis_def) {
    r1 = TMath::Abs(e.y);
    m_min = TMath::Min(m_min, e.y);
  } else {
    std::cout << "unknown axis" << std::endl;
  }
  if ( r1 <= m_closest.r) {
    if (m_closest.r<=m_second.r)
    {
      m_second = m_closest;
    }
    m_closest.r = r1;
    m_closest.dist = e;
    m_closest.h1 = p1;
    m_closest.h2 = p2;
  } else if ( r1 <= m_second.r)
  {
    m_second.r = r1;
    m_second.dist = e;
    m_second.h1 = p1;
    m_second.h2 = p2;
  }
  if (r1 == MAX_R_VALUE) {
    std::cout << "[find_nearest_strip] error r1== MAX_R_VALUE" << std::endl;
  }
}

void find_second_nearest_strip_processor::processEventEnd() {
  if (m_second.r < m_cutOff) {
    
    if (m_closest.h1.x!=m_second.h1.x)
    {
      std::cout << "more than one hit" << std::endl;
    }
      

    pushHit(m_closest.dist, 0);
    pushHit(m_closest.h1, 1);
    pushHit(m_closest.h2, 2);



    pushHit(m_second.dist, 10);
    pushHit(m_second.h1, 11);
    pushHit(m_second.h2, 12);
  }
}

s_plane_collection find_second_nearest_strip_processor::getOutputcollection() {
  s_plane_collection ret;
  ret.m_planes.push_back(std::make_pair(std::string(nearest_strip_distance_name), plane_def(getOutputName(), sct_type::ID_t(0))));
  ret.m_planes.push_back(std::make_pair(std::string(nearest_strip_plane1_name), plane_def(getOutputName(), sct_type::ID_t(1))));
  ret.m_planes.push_back(std::make_pair(std::string(nearest_strip_plane2_name), plane_def(getOutputName(), sct_type::ID_t(2))));


  ret.m_planes.push_back(std::make_pair(std::string(second_nearest_strip_distance_name), plane_def(getOutputName(), sct_type::ID_t(10))));
  ret.m_planes.push_back(std::make_pair(std::string(second_nearest_strip_plane1_name), plane_def(getOutputName(), sct_type::ID_t(11))));
  ret.m_planes.push_back(std::make_pair(std::string(second_nearest_strip_plane2_name), plane_def(getOutputName(), sct_type::ID_t(12))));
  return ret;
}


const char* find_second_nearest_strip_processor::getType() const {
  return "find_second_nearest_strip___";
}
}





namespace sct_processor {

find_second_nearest_strip::find_second_nearest_strip(
  const sct_corr::plane_def& planeA,
  const sct_corr::plane_def& planeB, 
  axis_def search_axis, 
  Double_t cutOfff, 
  const s_plot_prob& plot_prob_ /*= ""*/) {
  auto pl = getPlotCollectionIfUnique(planeA, planeB);

  if (!pl) {
    std::cout << "[s_plane_collection_find_closest processor::find_nearest] referencing to different plot collection\n";

  }
  m_planes = pl->addPlot(sct_corr::S_plot(new sct_corr::find_second_nearest_strip_processor(search_axis, cutOfff, plot_prob_)), planeA, planeB);




}

sct_corr::plane_def find_second_nearest_strip::getResidual() const {
  return m_planes.getByType(nearest_strip_distance_name)();
}

sct_corr::plane_def find_second_nearest_strip::getHitOnPlaneA() const {
  return m_planes.getByType(nearest_strip_plane1_name)();
}

sct_corr::plane_def find_second_nearest_strip::getHitOnPlaneB() const {
  return m_planes.getByType(nearest_strip_plane2_name)();
}

sct_corr::plane_def find_second_nearest_strip::getSecondResidual() const {
  return m_planes.getByType(second_nearest_strip_distance_name)();
}

sct_corr::plane_def find_second_nearest_strip::getSecondHitOnPlaneA() const {
  return m_planes.getByType(second_nearest_strip_plane1_name)();
}

sct_corr::plane_def find_second_nearest_strip::getsecondHitOnPlaneB() const {
  return m_planes.getByType(second_nearest_strip_plane2_name)();
}

}