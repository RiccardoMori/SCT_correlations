#include "sct_plots.h"

#include "treeCollection.h"
#include <iostream>
#include "cluster.h"
#include "TCanvas.h"
#include "TVector.h"
#include "TMath.h"
#include "TH2.h"
#include "TH2D.h"
#include "TAxis.h"
#include "TH1.h"
#include "TTree.h"
#include "s_cuts.h"
#include "s_plot_collection.h"
#include "plane_def.h"
#include "s_DrawOption.h"


#include "processorBase.h"
#include "sct_plots_internal.h"

#include "TF1.h"
#include "TProfile.h"
#include "geometry/setup_description.hh"
#include "plane_def.h"
#include "SCT_helpers.h"









namespace SCT_helpers{

void CutTH2(TH2* h, const S_Cut& cut_)
{

  if (!h)
  {
    std::cout << "Histogram is empty\n";
    return;
  }
  Double_t x = 0, y = 0, binContent = 0;

  for (Int_t y_bin = 0; y_bin <= h->GetNbinsY(); ++y_bin)
  {
    for (Int_t x_bin = 0; x_bin <= h->GetNbinsX(); ++x_bin)
    {
      x = h->GetXaxis()->GetBinCenter(x_bin);
      auto bin = h->GetBin(x_bin, y_bin);
      y = h->GetYaxis()->GetBinCenter(y_bin);
      binContent = h->GetBinContent(bin);


      if (cut_.isOutOfRange(binContent, x, y))
      {
        h->SetBinContent(bin, 0);
      }
    }

  }
}

void CutTH1(TH1* h1, const S_Cut& cut_)
{
  Double_t x = 0, binContent = 0;


  for (Int_t x_bin = 0; x_bin <= h1->GetNbinsX(); ++x_bin)
  {
    x = h1->GetXaxis()->GetBinCenter(x_bin);
    auto bin = h1->GetBin(x_bin);

    binContent = h1->GetBinContent(bin);


    if (cut_.isOutOfRange(binContent, x))
    {
      h1->SetBinContent(bin, 0);
    }
  }


}

TH1* HistogrammSilhouette(TH2* h2, axis_def ax)
{

  if (ax == x_axis_def)
  {

    std::string title = std::string(h2->GetTitle()) + std::string("_sil_x");
    std::string name = std::string(h2->GetName()) + std::string("_sil_x");
    TH1* ret = new TH1D(name.c_str(), title.c_str(), h2->GetNbinsX(), h2->GetXaxis()->GetBinCenter(0), h2->GetXaxis()->GetBinCenter(h2->GetNbinsX()));


    for (Int_t x_bin = 0; x_bin <= h2->GetNbinsX(); ++x_bin)
    {
      Double_t max_binContent = 0;
      for (Int_t y_bin = 0; y_bin <= h2->GetNbinsY(); ++y_bin)
      {


        auto bin = h2->GetBin(x_bin, y_bin);
        max_binContent = std::max(h2->GetBinContent(bin), max_binContent);

      }

      ret->SetBinContent(x_bin, max_binContent);
    }

    return ret;
  }
  else if (ax == y_axis_def){

    std::string title = std::string(h2->GetTitle()) + std::string("_sil_y");
    std::string name = std::string(h2->GetName()) + std::string("_sil_y");
    TH1* ret = new TH1D(name.c_str(), title.c_str(), h2->GetNbinsY(), h2->GetYaxis()->GetBinCenter(0), h2->GetYaxis()->GetBinCenter(h2->GetNbinsY()));

    for (Int_t y_bin = 0; y_bin <= h2->GetNbinsY(); ++y_bin)
    {
      Double_t max_binContent = 0;
      for (Int_t x_bin = 0; x_bin <= h2->GetNbinsX(); ++x_bin)
      {




        auto bin = h2->GetBin(x_bin, y_bin);
        max_binContent = std::max(h2->GetBinContent(bin), max_binContent);

      }

      ret->SetBinContent(y_bin, max_binContent);
    }

    return ret;
  }
  std::cout << "[SCT_helpers::HistogrammSilhouette] Unknown axis \n";
  return 0;
}



Long64_t DrawTTree(TTree * tree, const S_DrawOption& opt)
{
  return opt.Draw(tree);
}



TF1 LinearFit_Of_Profile(TH2* h2, const sct_type::procent_t& cut_prozent) {


  auto max_ = h2->GetMaximum();
  SCT_helpers::CutTH2(h2, S_Cut_BinContent(max_* necessary_CONVERSION( cut_prozent)/100));
  TProfile* p1 = h2->ProfileX();
  TF1 f("f1", "pol1", h2->GetXaxis()->GetBinCenter(0), h2->GetXaxis()->GetBinCenter(h2->GetNbinsX()));

  p1->Fit(&f, "Q");
  return f;
}





TH1* calc_efficiency(TH1* trueHits, TH1* dutHits) {

  TH1* effi = (TH1*) dutHits->Clone("efficiency");
  effi->Divide(trueHits);
  return effi;

}
void saveTH1_as_txt(const TH1& h1, const char* nameTXT) {
  std::ofstream out(nameTXT);
  out << h1.GetTitle() << std::endl;
  out << h1.GetXaxis()->GetTitle() << " ; " << h1.GetYaxis()->GetTitle()<<std::endl;

  for (int i = 0; i < h1.GetNbinsX(); ++i) {

    out << h1.GetBinCenter(i) << " ; " << h1.GetBinContent(i) << std::endl;
  }

}


void saveTH1_as_txt(const TProfile& h1, const char* nameTXT) {
  std::ofstream out(nameTXT);
  out << h1.GetTitle() << std::endl;
  out << h1.GetXaxis()->GetTitle() << " ; " << h1.GetYaxis()->GetTitle() << std::endl;

  for (int i = 0; i < h1.GetNbinsX(); ++i) {

    out << h1.GetBinCenter(i) << " ; " << h1.GetBinContent(i) << std::endl;
  }

}

TObject* Draw_obj(const sct_corr::plane_def& plane_, const S_DrawOption& opt) {
  auto pl = plane_.get_plot();
  pl->Draw(plane_, opt);
  
  if (opt.get_output_object())
  {
    return opt.get_output_object();
  }
 return gPad->GetPrimitive("htemp");
}

}	