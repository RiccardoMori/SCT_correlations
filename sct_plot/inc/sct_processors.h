#ifndef sct_processors_h__
#define sct_processors_h__
#include "plane_def.h"
#include "internal/platform.hh"
#include "s_cuts.h"
#include "s_plot_prob.h"



namespace sct_corr {
class DllExport processor {
public:
  static  plane_def cut_x_y(
    const plane_def& planeA,
    const  S_Cut& cut_,
    const s_plot_prob& = ""
    );

  static  plane_def remove_closest(
    const plane_def& planeA,
    const plane_def& planeB,
    const  axis_def search_axis,
    const s_plot_prob& = ""
    );

  static plane_def if_a_get_b(
    const plane_def& planeA,
    const plane_def& planeB,
    const S_Cut& cut_,
    const s_plot_prob& plot_prob_ = ""
    );
  static  plane_def cluster_strip(
    const plane_def& hits,
    axis_def search_axis,
    double cluster_distance,
    const s_plot_prob& = ""
    );


  static  plane_def coordinate_transform(
    const plane_def& planeA,
    Double_t x_slope,
    Double_t x_offset,
    Double_t y_slope,
    Double_t y_offset,
    const s_plot_prob& = ""
    );


  static  plane_def coordinate_transform_move(
    const plane_def& planeA,
    Double_t x_offset,
    Double_t y_offset,
    const s_plot_prob& = ""
    );





  static s_plane_collection_find_closest modulo_find_nearest_strip(
    const plane_def& planeA,
    const plane_def& planeB,
    axis_def search_axis,
    Double_t cutOfff,
    Double_t modulo_parameter,
    const s_plot_prob& plot_prob_ = ""
    );




  static  plane_def correlation(
    const S_Axis& axisA,
    const S_Axis& axisB,
    const s_plot_prob& = ""
    );


  static  plane_def correlation(
    const plane_def& planeA,
    const s_plot_prob& = ""
    );

  static  plane_def rotate(
    const plane_def& planeA,
    Double_t angle,
    const s_plot_prob& = ""
    );

  static  plane_def linear_trans(
    const plane_def& planeA,
    Double_t a11, Double_t a21,
    Double_t a12, Double_t a22,
    const s_plot_prob& = ""
    );

  static plane_def hitmap(
    const plane_def& hits,
    const s_plot_prob& = ""
    );


  static  plane_def hitmap(
    const S_Axis& axisA,
    const S_Axis& axisB,
    const s_plot_prob& = ""
    );

  static  plane_def moduloHitMap(
    const S_Axis& axisA,
    const S_Axis& axisB,
    double mod_x, double mod_y,
    const s_plot_prob& = ""
    );


  static  plane_def moduloHitMap(
    const plane_def& hits,
    double mod_x, double mod_y,
    const s_plot_prob& = ""
    );




  static  plane_def residual(
    const S_Axis& axisA,
    const S_Axis& axisB,
    const s_plot_prob& = ""
    );

  static  sct_corr::plane_def_Alibava residual_with_charge(
    const sct_corr::plane_def_Alibava& hits_A,
    const plane_def& hits_B,
    bool full_cluster = true,
    const s_plot_prob& = ""
    );
  static plane_def convert_local_to_global(
    const plane_def& local_hits,
    const sct_corr::Xlayer&,
    const s_plot_prob& = ""
    );
  static plane_def convert_global_to_local(
    const plane_def& local_hits,
    const sct_corr::Xlayer&,
    const s_plot_prob& = ""
    );
  static plane_def convert_hits_to_zs_data_GBL(
    const plane_def& hits,
    const sct_corr::Xlayer& layer,
    const s_plot_prob& plot_prob_ /*= ""*/
    );


  static  plane_def convert_zs_data_to_hits_GBL(
    const plane_def& sz_data,
    const sct_corr::Xlayer& layer,
    const s_plot_prob& plot_prob_/*= ""*/
    );

  static plane_def convert_zs_data_to_hits_DAF(
    const plane_def& sz_data,
    const sct_corr::Xlayer& layer,
    const s_plot_prob& plot_prob_/*= ""*/
    );

  static  plane_def convert_hits_to_zs_data_DAF(
    const plane_def& hits,
    const sct_corr::Xlayer& layer,
    const s_plot_prob& plot_prob_
    );
};

}
#endif // sct_processors_h__
