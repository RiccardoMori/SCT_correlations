#include "internal/PLaneVsPlane_plots.hh"
#include <iostream>
#include "TMath.h"
#include "TH2.h"
#include "TCanvas.h"
#include "internal/plane_hit.hh"
#include "sct_plots.h"

namespace sct_corr{
  class efficieny_map :public plotPlaneVsPlane{
  public:
    efficieny_map(Double_t x_bin_, Double_t y_bin_, Double_t x_cut, Double_t y_cut, const s_plot_prob& = "");
    virtual void processHit(const plane_hit&  p1, const plane_hit&  p2) override;
    virtual s_plane_collection getOutputcollection() override;
    virtual Long64_t Draw(const char* options, const char* cuts /* = "" */, const char* axis /* = "y:x" */) override;


    Double_t m_xbin_size = 1, m_ybin_size = 1;
    Double_t m_xCut = 1, m_yCut = 1;
    Double_t m_xMin = 10000, m_xMax = -10000;
    Double_t m_yMin = 10000, m_yMax = -10000;

    plane_hit make_bin_hit(const plane_hit& h) const;
    void set_min_max(const plane_hit& h);

    virtual const char* getType() const override;
  };


  efficieny_map::efficieny_map(Double_t x_bin_, Double_t y_bin_, Double_t x_cut, Double_t y_cut, const s_plot_prob& plot_prob) :
    plotPlaneVsPlane(plot_prob), m_xbin_size(x_bin_), m_ybin_size(y_bin_), m_xCut(x_cut), m_yCut(y_cut)
  {

  }

  void efficieny_map::processHit(const plane_hit& p1, const plane_hit& p2)
  {
    auto binp1 = make_bin_hit(p1);
    auto binp2 = make_bin_hit(p2);
    auto x = binp2.x;
    auto y = binp2.y;

    if (TMath::Abs(p1.x - p2.x) < m_xCut&& TMath::Abs(p1.y - p2.y) < m_yCut)
    {

      pushHit(binp1.x, binp1.y, 0);
      pushHit(binp2.x, binp2.y, 1);

      set_min_max(binp1);
      set_min_max(binp2);
    }


    pushHit(binp1.x, binp1.y, 2);
    pushHit(binp2.x, binp2.y, 3);
  }

  s_plane_collection efficieny_map::getOutputcollection()
  {
    s_plane_collection ret;
    ret.m_planes.push_back(std::make_pair(std::string("A_and_B_with_A_binning"), plane_def(getOutputName(), sct_type::ID_t(0))));
    ret.m_planes.push_back(std::make_pair(std::string("A_and_B_with_B_binning"), plane_def(getOutputName(), sct_type::ID_t(1))));
    ret.m_planes.push_back(std::make_pair(std::string("A_binned"), plane_def(getOutputName(), sct_type::ID_t(2))));
    ret.m_planes.push_back(std::make_pair(std::string("B_binned"), plane_def(getOutputName(), sct_type::ID_t(3))));

    return ret;
  }

  Long64_t efficieny_map::Draw(const char* options, const char* cuts /* = "" */, const char* axis /* = "y:x" */)
  {

    //   auto x_bins = (m_xMax - m_xMin) / m_xbin_size;
    //   auto y_bins = (m_yMax - m_yMin) / m_ybin_size;
    //   TH2D* h1 = new TH2D("A_and_B", "A_and_B", x_bins, m_xMin, m_xMax, y_bins, m_yMin, m_yMax);
    //  // new TCanvas();
    // auto elements=  plotPlaneVsPlane::Draw(options, "ID==1", "y:x>>A_and_B");
    // 
    // 
    //     TH2D* h2 = new TH2D("B", "B", x_bins, m_xMin, m_xMax, y_bins, m_yMin, m_yMax);
    // //     new TCanvas();
    //    plotPlaneVsPlane::Draw(options, "ID==3", "y:x>>B");
    // 
    //    h1->Divide(h2);
    // //   new TCanvas();
    //    h1->Draw(options);
    //   return elements;
    //todo
    return plotPlaneVsPlane::Draw(options, "ID==3", "y:x");
  }


  plane_hit efficieny_map::make_bin_hit(const plane_hit& h) const
  {
    return plane_hit(static_cast<int>((h.x + 10000) / m_xbin_size)*m_xbin_size - 10000, static_cast<int>((h.y + 10000) / m_ybin_size)*m_ybin_size - 10000);
  }

  void efficieny_map::set_min_max(const plane_hit& h)
  {
    m_xMin = TMath::Min(m_xMin, h.x);
    m_xMax = TMath::Max(m_xMax, h.x);

    m_yMin = TMath::Min(m_yMin, h.y);
    m_yMax = TMath::Max(m_yMax, h.y);


  }



  const char* efficieny_map::getType() const
  {
    return "efficieny_map__";
  }
}
sct_corr::S_plot sct_plot::efficiency_map(Double_t x_bin_, Double_t y_bin_, Double_t x_cut, Double_t y_cut, const s_plot_prob& plot_prob)
{

  return sct_corr::S_plot(new sct_corr::efficieny_map(x_bin_, y_bin_, x_cut, y_cut, plot_prob));
}

sct_corr::S_plot sct_plot::efficiency_map(Double_t x_bin_, Double_t y_bin_, const s_plot_prob& plot_prob)
{

  return efficiency_map( x_bin_, y_bin_, x_bin_, y_bin_, plot_prob);
}



