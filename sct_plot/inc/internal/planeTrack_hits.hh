#ifndef planeTrack_hits_h__
#define planeTrack_hits_h__

#include "internal/plane.hh"

namespace sct_corr{
  class planeTrack_hits :public plane{
  public:
    planeTrack_hits(Parameter_ref pl);

    virtual std::string getType() const override;
    virtual bool next() const override;

    virtual const plane_hit * getHit() const override;


  private:
    const axis_ref* m_x = nullptr;
    const axis_ref* m_y = nullptr;
    const axis_ref* m_chi2 = nullptr;
    const axis_ref* m_phi = nullptr;
    const axis_ref* m_theta = nullptr;
    const axis_ref* m_ndf = nullptr;
    std::shared_ptr<track_hits> m_hit;
    
  };

}


#endif // planeTrack_hits_h__
