#include "processors/find_nearest_strip.hh"
#include "plane_def.h"
#include <iostream>
#include "s_plot_collection.h"
#include "s_plot_prob.h"
#include "s_plot.h"
#include "internal/find_nearest_strip_processor.hh"


namespace sct_processor {
find_nearest_strip::find_nearest_strip(
  const sct_corr::plane_def& planeA,
  const sct_corr::plane_def& planeB,
  axis_def search_axis,
  Double_t cutOfff,
  const s_plot_prob& plot_prob_ /*= "" */
  ) {
  auto pl = getPlotCollectionIfUnique(planeA, planeB);

  if (!pl) {
    std::cout << "[s_plane_collection_find_closest processor::find_nearest] referencing to different plot collection\n";

  }
  auto collection_ = pl->addPlot(sct_corr::S_plot(new sct_corr::find_nearest_strip_processor(search_axis, cutOfff, plot_prob_)), planeA, planeB);
  
  
  m_planes.setResidual(collection_.getByType(nearest_strip_distance_name)());
  m_planes.setHitOnPlaneA(collection_.getByType(nearest_strip_plane1_name)());
  m_planes.setHitOnPlaneB(collection_.getByType(nearest_strip_plane2_name)());

}


sct_corr::plane_def find_nearest_strip::getResidual() const {
  return m_planes.getResidual();
}

sct_corr::plane_def find_nearest_strip::getHitOnPlaneA() const {
  return m_planes.getHitOnPlaneA();
}

sct_corr::plane_def find_nearest_strip::getHitOnPlaneB() const {
  return m_planes.getHitOnPlaneB();
}
}