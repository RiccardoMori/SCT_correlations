#ifndef sct_plots_h__
#define sct_plots_h__

#ifndef __CINT__
#include <memory>
#include <vector>
#endif


#include "TH2.h"
#include "TTree.h"
#include "TFile.h"
#include "TH1.h"
#include "TString.h"
#include "TF1.h"
#include "Rtypes.h"
#include "TCut.h"



#include "internal/platform.hh"
#include "s_plot.h"
#include "S_Axis.h"
#include "s_plot_prob.h"
#include "s_plot_collection.h"
#include "s_DrawOption.h"
#include "s_plane_def.h"



class TH2D;

class S_plane;
class S_DrawOption;
class S_Axis;
class S_treeCollection;
class S_Cut;
class S_plot_collection;
class s_plane_collection;
class s_process_collection;
namespace sct_corr {
  class treeCollection;

  class axis_ref;
  class plot;
  class plane;
  class sct_event_buffer;
  struct Xgear;
  struct Xlayer;
}

class TH2;
class DllExport SCT_helpers {
public:
  static void CutTH2(TH2* h2, const S_Cut& cut_);
  static void CutTH1(TH1* h1, const S_Cut& cut_);
  static TH1* HistogrammSilhouette(TH2* h2, axis_def ax);
  static Long64_t DrawTTree(TTree * tree, const S_DrawOption& opt);
  static TF1 LinearFit_Of_Profile(TH2* h2,Double_t cut_prozent = 0);
  static TH1* calc_efficiency(TH1* trueHits, TH1* dutHits);
};
class S_plane_def;
class DllExport sct_coll {
public:
  static  S_plane_def     apix_hit_local();
  static  S_plane_def     apix_hit();
  static  S_plane_def     apix_zs_data();
  static  S_plane_def     apix_fitted();
  static  S_plane_def     apix_fitted_local();
  static  S_plane_def     DUT_hit_local();
  static  S_plane_def     DUT_hit();
  static  S_plane_def     DUT_zs_data();
  static  S_plane_def     DUT_fitted();
  static  S_plane_def     DUT_fitted_local();
  static  S_plane_def_GBL DUT_fitted_local_GBL();
  static  S_plane_def     DUT_TTC_data();
  static  S_plane_def     DUT_Timestamp();
  static  S_plane_def     DUT_Timestamp_L0ID();
  static  S_plane_def     DUT_TDC_L0ID();
  static  S_plane_def     DUT_TLU_TLUID();
  static  S_plane_def     tel_hit_local(double ID);
  static  S_plane_def     tel_hit(double ID);
  static  S_plane_def     tel_zs_data(double ID);
  static  S_plane_def     tel_fitted(double ID);
  static  S_plane_def     tel_fitted_local(double ID);
  static  S_plane_def_GBL tel_fitted_local_GBL(double ID);
};


class DllExport sct {
public:
  static const char* col_hit();
  static const char* col_zsdata_strip();
  static const char* col_zsdata_m26();
  static const char* col_zsdata_apix();
  static const char* col_fitpoints();
  static const char* col_fitpoints_local();
  static const char* col_local_hit();
  static const char* col_DUT_TTC();
  static const char* col_GBL_fitted_points();

  static const char* plot_hitmap();
  static const char* plot_modulo();
  static const char* plot_correlation();
  static const char* plot_residual();
  static const char* plot_clusterSize();
  static const char* plot_clusterSize_strip();
  static const char* plot_projectOnPixel();
  static const char* plot_find_correspondingX();
  static const char* plot_find_correspondingXY();
  static const char* plot_Event_size();
  static const char* plot_find_nearest();
  static const char* plot_plane_distance();
  static const char* plot_find_nearest_strip();
  static const char* plot_A_if_B();
  static const char* plot_rotated();
  static const char* plot_linear_trans();
  static const char* plot_coordinate_transform();
  static const char* plot_hitMultiplizity();
  static const char* plot_efficieny_map();
  static const char* plot_cut_x_y();
};


class S_plot_def;
class S_plot;
class DllExport sct_plot {
public:
  static  S_plot hitmap(const s_plot_prob& = "");
  static  S_plot moduloHitMap(double modulo_parameter_x, double modulo_parameter_y, const s_plot_prob& = "");
  static  S_plot correlation(const s_plot_prob& = "");
  static  S_plot residual(const s_plot_prob& = "");
  static  S_plot clustering(Double_t Pixel_distance = 2, const s_plot_prob& = "");
  static  S_plot clustering_strip(axis_def searchAxis, Double_t Pixel_distance = 2, const s_plot_prob& = "");
  static  S_plot projectOnPixel(const s_plot_prob& = "");
  static  S_plot find_correspondingX(const s_plot_prob& = "");
  static  S_plot find_correspondingXY(const s_plot_prob& = "");
  static  S_plot Event_size(const s_plot_prob& = "");
  static  S_plot find_nearest(Double_t x_cutoff, Double_t y_cutoff, const s_plot_prob& = "");
  static  S_plot find_nearest_strip(axis_def search_axis, Double_t cutOfff = 100000, const s_plot_prob& = "");
  static  S_plot modulo_find_nearest_strip(axis_def search_axis, double modulo_param, Double_t cutOfff = 100000, const s_plot_prob& = "");

  static  S_plot plane_distance(const s_plot_prob& = "");
  static  S_plot A_if_B(const s_plot_prob& = "");
  static  S_plot rotated(Double_t angle, const s_plot_prob& = "");
  static  S_plot coordinate_transform(Double_t x_slope, Double_t x_offset, Double_t y_slope, Double_t y_offset, const s_plot_prob& = "");
  static  S_plot coordinate_transform_move(Double_t x_offset, Double_t y_offset, const s_plot_prob& = "");

  static  S_plot efficiency_map(Double_t x_bin, Double_t y_bin, Double_t x_cut, Double_t y_cut, const s_plot_prob& = "");
  static  S_plot efficiency_map(Double_t x_bin, Double_t y_bin, const s_plot_prob& = "");
  static  S_plot hitMultiplizity(const s_plot_prob& = "");

  static  S_plot cut_x_y(const  S_Cut& cut_, const s_plot_prob& = "");


  static  S_plot save2LCIO(const char* filename, unsigned runnum, const s_plot_prob& = "");

  static  S_plot linear_trans(Double_t a11, Double_t a21 ,Double_t a12,Double_t a22, const s_plot_prob& = "");

public:
  static  S_plane_def coordinate_transform(S_plot_collection& pl, Double_t x_slope, Double_t x_offset, Double_t y_slope, Double_t y_offset, const S_plane_def& planeA, const s_plot_prob& = "");
  static  S_plane_def coordinate_transform_move(S_plot_collection& pl, Double_t x_offset, Double_t y_offset, const S_plane_def& planeA, const s_plot_prob& = "");

  static s_plane_collection_find_closest find_nearest(S_plot_collection& pl, Double_t x_cutoff, Double_t y_cutoff, const S_plane_def& planeA, const S_plane_def& planeB, const s_plot_prob& = "");
  static s_plane_collection_find_closest find_nearest_strip(S_plot_collection& pl, axis_def search_axis, Double_t cutOfff, const S_plane_def& planeA, const S_plane_def& planeB, const s_plot_prob& plot_prob_= "");
  
  
  static s_plane_collection_find_closest modulo_find_nearest_strip(
    S_plot_collection& pl, 
    axis_def search_axis, 
    Double_t cutOfff, 
    Double_t modulo_parameter,
    const S_plane_def& planeA, 
    const S_plane_def& planeB, 
    const s_plot_prob& plot_prob_ = ""
    );

  static s_plane_collection misalignment_strip(S_plot_collection& pl, S_plane_def fitted_plane, S_plane_def plane_hit_, axis_def Unknown_axis, const s_plot_prob& = "");
  static s_plane_collection misalignment_pixel(S_plot_collection& pl, S_plane_def fitted_plane, S_plane_def plane_hit_, const  s_plot_prob& = "");
  static s_plane_collection plane_merger(S_plot_collection& pl, s_plane_collection planes_, const  s_plot_prob& = "");
  static S_plane_def Crate_True_Fitted_DUT_Hits(S_plot_collection& pl, const  s_plot_prob& = "");
  static S_plane_def Create_True_Fitted_DUT_Hits_in_channels(S_plot_collection& pl, double pitchSize, double rotate, double move_x, double move_y, const s_plot_prob& = "");
  static s_plane_collection_correlations Create_Correlations_of_true_Fitted_hits_with_DUT_Hits_in_channels(
                                                                                                            S_plot_collection& pl,
                                                                                                            const S_Cut& fiducial_cut_,
                                                                                                            double residualCut,
                                                                                                            const  sct_corr::Xgear& gear,
                                                                                                            double rotate_angle,
                                                                                                            double move_x,
                                                                                                            const s_plot_prob& plot_prob_/*= ""*/
                                                                                                            );

  static s_plane_collection_correlations DAF_Create_Correlations_of_true_Fitted_hits_with_DUT_Hits_in_channels(
                                                                                                                S_plot_collection& pl,
                                                                                                                const S_Cut& fiducial_cut_,
                                                                                                                double residualCut,
                                                                                                                const  sct_corr::Xgear& gear,
                                                                                                                double rotate_angle,
                                                                                                                double move_x,
                                                                                                                const s_plot_prob& plot_prob_/*= ""*/
                                                                                                                );
  
  static s_plane_collection_correlations GBL_Create_Correlations_of_true_Fitted_hits_with_DUT_Hits_in_channels(S_plot_collection& pl, 
                                                                                                               const S_Cut& fiducial_cut_, 
                                                                                                               double residualCut, 
                                                                                                               const  sct_corr::Xgear& gear, 
                                                                                                               double rotate_angle,
                                                                                                               double move_x,
                                                                                                               const s_plot_prob& plot_prob_ = ""
                                                                                                               );
  
  static s_plane_collection_correlations GBL_Create_Correlations_of_true_Fitted_hits_with_DUT_Hits(S_plot_collection& pl, 
                                                                                                   const S_Cut& fiducial_cut_, 
                                                                                                   double residualCut, 
                                                                                                   const  sct_corr::Xgear&, 
                                                                                                   double rotate_angle,
                                                                                                   double move_x, 
                                                                                                   const s_plot_prob& = ""
                                                                                                   );

  static S_plane_def convert_zs_data_to_hits_GBL(
    S_plot_collection& pl, 
    const sct_corr::Xlayer&, 
    const S_plane_def& sz_data, 
    const s_plot_prob& = ""
    );


  static S_plane_def convert_hits_to_zs_data_GBL(
    S_plot_collection& pl, 
    const sct_corr::Xlayer&, 
    const S_plane_def& hits, 
    const s_plot_prob& = ""
    );

  static S_plane_def convert_zs_data_to_hits_DAF(
    S_plot_collection& pl, 
    const sct_corr::Xlayer&, 
    const S_plane_def& sz_data, 
    const s_plot_prob& = ""
    );

  static S_plane_def convert_hits_to_zs_data_DAF(S_plot_collection& pl, const sct_corr::Xlayer&, const S_plane_def& hits, const s_plot_prob& = "");

  static S_plane_def convert_local_to_global(S_plot_collection& pl, const sct_corr::Xlayer&, const S_plane_def& local_hits, const s_plot_prob& = "");
  static S_plane_def convert_global_to_local(S_plot_collection& pl, const sct_corr::Xlayer&, const S_plane_def& local_hits, const s_plot_prob& = "");
  static  S_plane_def correlation(S_plot_collection& pl, const S_Axis& axisA, const S_Axis& axisB, const s_plot_prob& = "");
  static  S_plane_def correlation(S_plot_collection& pl, const S_plane_def& planeA,const s_plot_prob& = "");

  static  S_plane_def rotate(S_plot_collection& pl, Double_t angle, const S_plane_def& planeA, const s_plot_prob& = "");

  static  S_plane_def linear_trans(S_plot_collection& pl, Double_t a11, Double_t a21, Double_t a12, Double_t a22, const S_plane_def& planeA, const s_plot_prob& = "");
  static  S_plane_def hitmap(S_plot_collection& pl, const S_plane_def& hits, const s_plot_prob& = "");
  static  S_plane_def hitmap(S_plot_collection& pl, const S_Axis& axisA, const S_Axis& axisB, const s_plot_prob& = "");
  
  static  S_plane_def moduloHitMap(S_plot_collection& pl, double mod_x, double mod_y , const S_Axis& axisA, const S_Axis& axisB, const s_plot_prob& = "");
  static  S_plane_def moduloHitMap(S_plot_collection& pl, double mod_x, double mod_y, const S_plane_def& hits, const s_plot_prob& = "");
 
  static  S_plane_def cluster_strip(S_plot_collection& pl, axis_def search_axis,double cluster_distance, const S_plane_def& hits, const s_plot_prob& = "");


  static  S_plane_def cut_x_y(S_plot_collection& pl, const  S_Cut& cut_, const S_plane_def& planeA, const s_plot_prob& = "");





  static  S_plane_def residual(S_plot_collection& pl, const S_Axis& axisA, const S_Axis& axisB, const s_plot_prob& = "");
};


class DllExport sct_analyis { // all this function have memory leaks. they are not designed to be called in a loop 

public:
  static TH2* misalignment(S_plot_collection& treeColl);
  static TF1* Draw_misalignment(TFile& __file0);
  static TF1* Draw_New_aligment(TFile& __file0, TF1 *);
  static S_Cut* Get_fiducial_area(TFile& __file0);
};





#ifdef __CINT__

#pragma link C++ class sct;
#pragma link C++ class SCT_helpers;

#pragma link C++ class sct_analyis;
#pragma link C++ class sct_plot;
#pragma link C++ class sct_coll;


#endif







#endif // sct_plots_h__
