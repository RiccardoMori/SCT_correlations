#include "processors/find_nearest.hh"
#include "plane_def.h"
#include <iostream>
#include "s_plot_collection.h"
#include "internal/find_nearest_processor.hh"
#include "s_plot_prob.h"
#include "s_plot.h"


namespace sct_processor {
find_nearest::find_nearest(
  const sct_corr::plane_def& planeA,
  const sct_corr::plane_def& planeB,
  Double_t x_cutoff,
  Double_t y_cutoff,
  const s_plot_prob& plot_prob__/*= "" */
  ) {
  auto pl = getPlotCollectionIfUnique(planeA, planeB);

  if (!pl) {
    std::cout << "[s_plane_collection_find_closest processor::find_nearest] referencing to different plot collection\n";

  }

  auto collection_ = pl->addPlot(sct_corr::S_plot(new sct_corr::find_nearest_processor(x_cutoff, y_cutoff, plot_prob__)), planeA, planeB);


  m_planes.setResidual(collection_.getByType("nearest_distance")());
  m_planes.setHitOnPlaneA(collection_.getByType("nearest_plane1")());
  m_planes.setHitOnPlaneB(collection_.getByType("nearest_plane2")());

}


sct_corr::plane_def find_nearest::getResidual() const {
  return m_planes.getResidual();
}

sct_corr::plane_def find_nearest::getHitOnPlaneA() const {
  return m_planes.getHitOnPlaneA();
}

sct_corr::plane_def find_nearest::getHitOnPlaneB() const {
  return m_planes.getHitOnPlaneB();
}
}