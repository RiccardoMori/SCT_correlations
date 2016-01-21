#ifndef SCT_helpers_h__
#define SCT_helpers_h__
#include "internal/platform.hh"
#include "Rtypes.h"
#include "S_Axis.h"
#include "TF1.h"
#include "plane_def.h"
#include "TProfile.h"
#include "sct_types.h"
#include "TH2.h"
#include "TH1.h"
class TH2D;
class TH2;
class TH1;
class TTree;
class S_Cut;
class S_DrawOption;


namespace  SCT_helpers {

  DllExport void CutTH2(TH2* h2, const S_Cut& cut_);
  DllExport void CutTH1(TH1* h1, const S_Cut& cut_);
  DllExport TH1* HistogrammSilhouette(TH2* h2, axis_def ax);
  DllExport Long64_t DrawTTree(TTree * tree, const S_DrawOption& opt);
  //static TF1 LinearFit_Of_Profile(TH2* h2, Double_t cut_prozent = 0);
  DllExport TF1 LinearFit_Of_Profile(TH2* h2, const sct_type::procent_t& cut_prozent);
  DllExport TH1* calc_efficiency(TH1* trueHits, TH1* dutHits);
  DllExport void saveTH1_as_txt(const TH1& h1, const char* nameTXT);


  DllExport TObject* Draw_obj(const sct_corr::plane_def& plane_, const S_DrawOption& opt);



  template <typename returnType>
  inline returnType* Draw(const sct_corr::plane_def& plane_, const S_DrawOption& opt) {
    return dynamic_cast<returnType*>(Draw_obj(plane_, opt));
  }

  template <>
  inline TH1* Draw(const sct_corr::plane_def& plane_, const S_DrawOption& opt){
    auto ret = dynamic_cast<TH1*>(Draw_obj(plane_, opt));
  if (ret) {
    ret->SetTitle(necessary_CONVERSION(plane_.getName()).c_str());
  }
  return ret;
  }
  template <>
  inline TH2* Draw(const sct_corr::plane_def& plane_, const S_DrawOption& opt){
    auto ret = dynamic_cast<TH2*>(Draw_obj(plane_, opt));
  if (ret) {
    ret->SetTitle(necessary_CONVERSION(plane_.getName()).c_str());
  }
  return ret;
  }

  DllExport void saveTH1_as_txt(const TProfile& h1, const char* nameTXT);
  
};



#endif // SCT_helpers_h__
