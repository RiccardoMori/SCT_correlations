#ifndef planeAlibava_hit_h__
#define planeAlibava_hit_h__


#include "internal/plane.hh"

namespace sct_corr{
  class planeAlibava_hit :public plane{
  public:
    planeAlibava_hit(Parameter_ref pl);

    virtual std::string getType() const override;
    virtual bool next() const override;

    virtual const plane_hit * getHit() const override;


  private:
    const axis_ref* m_x = nullptr;
    const axis_ref* m_y = nullptr;
    const axis_ref* m_charge = nullptr;
    std::shared_ptr<Alibava_hit> m_hit;
    
  };

}

#endif // planeAlibava_hit_h__


