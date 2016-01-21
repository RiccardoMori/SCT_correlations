#include "internal/PLaneVsPlane_plots.hh"
#include <iostream>
#include "TMath.h"
#include "plane_def.h"
#include "internal/find_nearest_strip_processor.hh"
#include "sct_plots.h"
#include "sct_processors.h"
#include "internal/exceptions.hh"




namespace sct_corr {

class remove_closest :public plotPlaneVsPlane {
public:
  remove_closest(axis_def search_axis):m_search_axis(search_axis) {
    m_dist = 10000;
  }

  virtual void processEventStart() override;
  virtual void processHit(const plane_hit&  p1, const plane_hit&  p2)override;

  virtual void processEventEnd() override {}
  virtual void processHit_1(const plane_hit&  p1) override;
  virtual const char* getType() const override {
    return "remove_closest";
  }
  virtual s_plane_collection getOutputcollection() override;
private:
  double m_x, m_y, m_dist;
  axis_def m_search_axis;
};

void remove_closest::processEventStart() {
  m_dist = 10000;
}

void remove_closest::processHit(const plane_hit& p1, const plane_hit& p2) {

  if (m_search_axis == x_axis_def) {
    if (TMath::Abs(p1.x-p2.x) <m_dist )   {
      m_dist = TMath::Abs(p1.x - p2.x);
      m_x = p1.x;
      m_y = p1.y;
    }


  } else  if (m_search_axis==y_axis_def){
    if (TMath::Abs(p1.y - p2.y) < m_dist) {
      m_dist = TMath::Abs(p1.y - p2.y);
      m_x = p1.x;
      m_y = p1.y;
    }
  }
}



void remove_closest::processHit_1(const plane_hit& p1) {
  if (p1.x == m_x && p1.y== m_y)  {
    return; //skip closest
  }
  pushHit(p1,0);

}

s_plane_collection remove_closest::getOutputcollection() {
  s_plane_collection ret;
  ret.m_planes.push_back(std::make_pair(std::string("remove_closest"), plane_def(getOutputName(), sct_type::ID_t(0))));
  return ret;
}

}
sct_corr::plane_def sct_corr::processor::remove_closest(const plane_def& planeA, const plane_def& planeB, const axis_def search_axis, const s_plot_prob& /*= "" */) {

  auto pl = getPlotCollectionIfUnique(planeA, planeB);

  if (!pl) {
    SCT_THROW("referencing to different plot collection");
  }
  
  auto ret = pl->addPlot(S_plot(new sct_corr::remove_closest(search_axis)), planeA, planeB)();

  return ret;
}
