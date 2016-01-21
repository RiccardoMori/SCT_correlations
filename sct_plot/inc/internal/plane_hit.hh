#ifndef plane_hit_h__
#define plane_hit_h__
#include "Rtypes.h"


namespace sct_corr{
  class plane_hit
  {
  public:
    plane_hit(Double_t x_, Double_t y_) :x(x_), y(y_){}
    virtual ~plane_hit() {}
    Double_t x, y;
  };

  class track_hits :public plane_hit{
  public:
    track_hits(Double_t x_, Double_t y_, Double_t chi2_, Double_t ndf_, Double_t phi_, Double_t theta_) :plane_hit(x_, y_), chi2(chi2_), phi(phi_), theta(theta_), ndf(ndf_) {}
    virtual ~track_hits() {}

    Double_t chi2, phi,theta,ndf;
  };

  class Alibava_hit :public plane_hit {
  public:
    Alibava_hit() :Alibava_hit(0, 0, 0) {}
    Alibava_hit(Double_t x_, Double_t y_, Double_t charge__) :plane_hit(x_, y_), charge(charge__) {}

    virtual ~Alibava_hit() {}
    Double_t charge;
  };
}


#endif // plane_hit_h__
