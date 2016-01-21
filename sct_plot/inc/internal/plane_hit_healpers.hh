#ifndef plane_hit_healpers_h__
#define plane_hit_healpers_h__
#include "sct_plots.h"
#include "TMath.h"
#include "plane_hit.hh"
namespace sct_corr{
  inline plane_hit rotate(const  plane_hit& h, double Angle){

    double x = h.x*TMath::Cos(Angle) - h.y*TMath::Sin(Angle);
    double y = h.x*TMath::Sin(Angle) + h.y*TMath::Cos(Angle);
    return plane_hit(x, y);
  }
}

#endif // plane_hit_healpers_h__
