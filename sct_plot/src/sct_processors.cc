#include "sct_processors.h"
#include "sct_plots.h"
#include "TMath.h"
#include "geometry/setup_description.hh"
#include "internal/exceptions.hh"


namespace sct_corr {


plane_def processor::cut_x_y(
  const plane_def& planeA,
  const S_Cut& cut_,
  const s_plot_prob& plot_prob/*= ""*/
  ) {
  auto pl = planeA.get_plot();

  plane_def ret = pl->addPlot(sct_plot::cut_x_y(cut_, plot_prob), planeA)();

  return ret;
}

sct_corr::plane_def processor::if_a_get_b(
  const plane_def& planeA, 
  const plane_def& planeB, 
  const S_Cut& cut_, 
  const s_plot_prob& plot_prob_ /*= "" */
  ) {
  auto pl = getPlotCollectionIfUnique(planeA, planeB);

  if (!pl) {
    SCT_THROW("referencing to different plot collection");
  }

  plane_def a_cutted = processor::cut_x_y(planeA, cut_, s_plot_prob().doNotSaveToDisk());
  plane_def ret = pl->addPlot(sct_plot::A_if_B(plot_prob_), planeB, a_cutted)();
  return ret;
}

plane_def processor::cluster_strip(
  const plane_def& hits,
  axis_def search_axis,
  double cluster_distance,
  const s_plot_prob& plot_prob
  ) {
  auto pl = hits.get_plot();
  plane_def ret = pl->addPlot(sct_plot::clustering_strip(search_axis, cluster_distance, plot_prob), hits)();


  return ret;

}



plane_def processor::coordinate_transform(
  const plane_def& planeA,
  Double_t x_slope, Double_t x_offset,
  Double_t y_slope, Double_t y_offset,
  const s_plot_prob& plot_prob
  ) {
  auto pl = planeA.get_plot();
  plane_def ret = pl->addPlot(
    sct_plot::coordinate_transform(
    x_slope, x_offset,
    y_slope, y_offset,
    plot_prob), planeA
    )();


  return ret;
}

plane_def processor::coordinate_transform_move(
  const plane_def& planeA,
  Double_t x_offset, Double_t y_offset,
  const s_plot_prob& plot_prob__
  ) {
  auto pl = planeA.get_plot();
  plane_def ret = pl->addPlot(
    sct_plot::coordinate_transform_move(
    x_offset, y_offset,
    plot_prob__), planeA
    )();


  return ret;
}





s_plane_collection_find_closest processor::modulo_find_nearest_strip(
  const plane_def& planeA,
  const plane_def& planeB,
  axis_def search_axis,
  Double_t cutOfff,
  Double_t modulo_parameter,
  const s_plot_prob& plot_prob_ /*= "" */
  ) {

  auto pl = getPlotCollectionIfUnique(planeA, planeB);

  if (!pl) {
    SCT_THROW("referencing to different plot collection");
  }

  auto collection_ = pl->addPlot(
    sct_plot::modulo_find_nearest_strip(search_axis, modulo_parameter, cutOfff, plot_prob_),
    planeA,
    planeB
    );

  s_plane_collection_find_closest ret;
  ret.setResidual(collection_.getByType("nearest_strip_distance")());
  ret.setHitOnPlaneA(collection_.getByType("nearest_strip_plane1")());
  ret.setHitOnPlaneB(collection_.getByType("nearest_strip_plane2")());

  return ret;
}


plane_def processor::correlation(
  const S_Axis& axisA,
  const S_Axis& axisB,
  const s_plot_prob& plot_prob__/*= "" */
  ) {
  auto pl = getPlotCollectionIfUnique(axisA, axisB);

  if (!pl) {
    SCT_THROW("referencing to different plot collection");
  }

  auto ret = pl->addPlot(sct_plot::correlation(plot_prob__), axisA, axisB)();

  return ret;
}


plane_def processor::correlation(
  const plane_def& planeA,
  const s_plot_prob& plot_prob__/*= "" */
  ) {
  auto plA = planeA.get_plot();
  auto ret = plA->addPlot(sct_plot::correlation(plot_prob__), planeA)();

  return ret;
}

plane_def processor::rotate(
  const plane_def& planeA,
  Double_t angle,
  const s_plot_prob& plot_prob__/*= "" */
  ) {

  auto plA = planeA.get_plot();
  auto ret = plA->addPlot(sct_plot::rotated(angle, plot_prob__), planeA)();


  return ret;
}

plane_def processor::linear_trans(
  const plane_def& planeA,
  Double_t a11, Double_t a21,
  Double_t a12, Double_t a22,
  const s_plot_prob& plot_prob__
  ) {
  auto plA = planeA.get_plot();
  auto ret = plA->addPlot(sct_plot::linear_trans(
    a11, a21,
    a12, a22,
    plot_prob__), planeA
    )();


  return ret;
}

plane_def processor::hitmap(
  const plane_def& hits,
  const s_plot_prob& plot_prob__
  ) {
  auto plA = hits.get_plot();
  auto ret = plA->addPlot(sct_plot::hitmap(plot_prob__), hits)();

  return ret;

}

plane_def processor::hitmap(
  const S_Axis& axisA,
  const S_Axis& axisB,
  const s_plot_prob& plot_prob__
  ) {
  auto pl = getPlotCollectionIfUnique(axisA, axisB);

  if (!pl) {
    SCT_THROW("referencing to different plot collection");
  }

  auto ret = pl->addPlot(sct_plot::hitmap(plot_prob__), axisA, axisB)();

  return ret;
}

plane_def processor::moduloHitMap(
  const S_Axis& axisA,
  const S_Axis& axisB,
  double mod_x, double mod_y,
  const s_plot_prob& plot_prob
  ) {
  auto pl = getPlotCollectionIfUnique(axisA, axisB);

  if (!pl) {
    SCT_THROW("referencing to different plot collection");
  }
  auto ret = pl->addPlot(sct_plot::moduloHitMap(mod_x, mod_y, plot_prob), axisA, axisB)();

  return ret;
}

plane_def processor::moduloHitMap(
  const plane_def& hits,
  double mod_x,
  double mod_y,
  const s_plot_prob& plot_prob
  ) {
  auto plA = hits.get_plot();

  auto ret = plA->addPlot(sct_plot::moduloHitMap(mod_x, mod_y, plot_prob), hits)();

  return ret;
}

plane_def processor::residual(
  const S_Axis& axisA,
  const S_Axis& axisB,
  const s_plot_prob& plot_prob
  ) {
  auto plA = axisA.get_plot();
  auto ret = plA->addPlot(sct_plot::residual(plot_prob), axisA, axisB)();

  return ret;
}



plane_def processor::convert_local_to_global(
  const plane_def& local_hits,
  const sct_corr::Xlayer& layer,
  const s_plot_prob& plot_prob_
  ) {

  auto rot = rotate(
    local_hits,
    layer.ladder.rotationXY / 360 * TMath::Pi() * 2,
    s_plot_prob().doNotSaveToDisk()
    );
  auto a = linear_trans(
    rot,
    layer.sensitive.rotation1,
    layer.sensitive.rotation2,
    layer.sensitive.rotation3,
    layer.sensitive.rotation4,
    s_plot_prob().doNotSaveToDisk()
    );


  auto b = coordinate_transform_move(
    a,
    layer.ladder.positionX,
    layer.ladder.positionY,
    plot_prob_
    );

  return b;
}


plane_def processor::convert_global_to_local(
  const plane_def& global_hits,
  const sct_corr::Xlayer& layer,
  const s_plot_prob& plot_prob_/*= ""*/
  ) {


  auto b = coordinate_transform_move(
    global_hits,
    -layer.ladder.positionX,
    -layer.ladder.positionY,
    s_plot_prob().doNotSaveToDisk()
    );

  auto det_A = layer.sensitive.rotation1*layer.sensitive.rotation4
    - layer.sensitive.rotation2*layer.sensitive.rotation3;

  if (det_A == 0) {
    SCT_THROW("unable to invert matrix : Det_A ==0 ");

  }


  auto a = linear_trans(
    b,
    layer.sensitive.rotation4 / det_A,
    -layer.sensitive.rotation2 / det_A,
    -layer.sensitive.rotation3 / det_A,
    layer.sensitive.rotation1 / det_A,
    s_plot_prob().doNotSaveToDisk()
    );

  auto local_hits = rotate(
    a,
    -layer.ladder.rotationXY / 360 * TMath::Pi() * 2,
    plot_prob_
    );

  return local_hits;
}

plane_def processor::convert_hits_to_zs_data_GBL(
  const plane_def& hits,
  const sct_corr::Xlayer& layer,
  const s_plot_prob& plot_prob_ /*= ""*/
  ) {

  return coordinate_transform(
    hits,
    1 / layer.sensitive.pitchX,                                                         // x_slope
    -(layer.sensitive.positionX - layer.sensitive.sizeX / 2) / layer.sensitive.pitchX,  //x_offset
    1 / layer.sensitive.pitchY,                                                         //y_slope
    -(layer.sensitive.positionY - layer.sensitive.sizeY / 2) / layer.sensitive.pitchY,  //y_offset
    plot_prob_
    );
}

plane_def processor::convert_zs_data_to_hits_GBL(const plane_def& sz_data, const sct_corr::Xlayer& layer, const s_plot_prob& plot_prob_/*= ""*/) {
  return coordinate_transform(
    sz_data,
    layer.sensitive.pitchX,                                     // x_slope
    layer.sensitive.positionX - layer.sensitive.sizeX / 2,      //x_offset
    layer.sensitive.pitchY,                                     //y_slope
    layer.sensitive.positionY - layer.sensitive.sizeY / 2,      //y_offset
    plot_prob_
    );
}

plane_def processor::convert_zs_data_to_hits_DAF(
  const plane_def& sz_data,
  const sct_corr::Xlayer& layer,
  const s_plot_prob& plot_prob_/*= ""*/
  ) {
  return coordinate_transform(
    sz_data,
    layer.sensitive.pitchX,          // x_slope
    -layer.sensitive.sizeX / 2,      //x_offset
    layer.sensitive.pitchY,          //y_slope
    -layer.sensitive.sizeY / 2,      //y_offset
    plot_prob_
    );
}

plane_def processor::convert_hits_to_zs_data_DAF(
  const plane_def& hits,
  const sct_corr::Xlayer& layer,
  const s_plot_prob& plot_prob_
  ) {
  return coordinate_transform(
    hits,
    1 / layer.sensitive.pitchX,                            // x_slope
    (layer.sensitive.sizeX / 2) / layer.sensitive.pitchX,  //x_offset
    1 / layer.sensitive.pitchY,                             //y_slope
    (layer.sensitive.sizeY / 2) / layer.sensitive.pitchY,  //y_offset
    plot_prob_
    );
}
}