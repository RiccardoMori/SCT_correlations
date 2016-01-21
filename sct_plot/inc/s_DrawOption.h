#ifndef s_DrawOption_h__
#define s_DrawOption_h__

#include "internal/platform.hh"
#include "Rtypes.h"
#include "TCut.h"
#include <string>
#include "S_Axis.h"

class TTree;
class TObject;
class DllExport S_DrawOption{
public:

  S_DrawOption(const char* options, const char* cuts = "", const char * axis = "y:x");
  S_DrawOption(const char* options, TCut cuts, const char * axis = "y:x");
  S_DrawOption();
  S_DrawOption& options(const char* option);
  S_DrawOption& opt_colz();
  S_DrawOption& opt_star();
  S_DrawOption& opt_bar();
  S_DrawOption& opt_same();
  S_DrawOption& opt_color(Color_t);
  S_DrawOption& color_red();


  S_DrawOption& color_white();
  S_DrawOption& color_black();
  S_DrawOption& color_gray();
  S_DrawOption& color_green();
  S_DrawOption& color_blue();
  S_DrawOption& color_yellow();
  S_DrawOption& color_magenta();
  S_DrawOption& color_cyan();
  S_DrawOption& color_orange();
  S_DrawOption& color_spring();
  S_DrawOption& color_teal();
  S_DrawOption& color_azure();
  S_DrawOption& color_violet();
  S_DrawOption& color_pink();


  S_DrawOption& cut(const char* cut_);
  S_DrawOption& cut_min(axis_def ax,Double_t min_);
  S_DrawOption& cut_max(axis_def ax, Double_t max_);
  S_DrawOption& cut(axis_def ax, Double_t min_,Double_t max_);
  S_DrawOption& cut_add(const TCut& cut_);
  S_DrawOption& cut_x(Double_t min_, Double_t max_);
  S_DrawOption& cut_x_min(Double_t min_);
  S_DrawOption& cut_x_max(Double_t max_);
  S_DrawOption& cut_y(Double_t min_, Double_t max_);
  S_DrawOption& cut_y_min(Double_t min_);
  S_DrawOption& cut_y_max(Double_t max_);
  S_DrawOption& draw_axis(const char* axis_);
  S_DrawOption& draw_axis(axis_def);
  S_DrawOption& draw_axis(axis_def,axis_def);
  S_DrawOption& draw_axis(axis_def,axis_def,axis_def);
  S_DrawOption& draw_x();
  S_DrawOption& draw_y();
  S_DrawOption& draw_x_VS_y();
  S_DrawOption& draw_y_VS_x();
  S_DrawOption& output_object(TObject* out_);
  TObject* get_output_object() const;
  Long64_t Draw(TTree * tree) const;
  const char* getOptions() const;
  const char* getAxis() const;
  TCut getCut() const;
private:
#ifndef __CINT__
  void push_output_to_TGraph(Long64_t numberOfPoints,TTree* output_tree)const;
     void push_output_to_TGraphErrors(Long64_t numberOfPoints,TTree* output_tree)const;
  std::string m_options = "colz", m_axis = "x:y";
  mutable std::string  m_out_option;
  TCut m_cut;
  TObject* m_output_object = nullptr;
  mutable std::string m_axis_dummy;
  mutable int m_numOfAxis = 2;
  mutable Color_t m_color = kBlack;
  bool m_same = false;
#endif // !__CINT__


};

#endif // s_DrawOption_h__
