#ifndef planeX_Y_hits_h__
#define planeX_Y_hits_h__
#include "internal/plane.hh"


namespace sct_corr{

  class planeX_Y :public plane{
  public:
    planeX_Y(Parameter_ref pl);
   const axis_ref* getX() const;
   const axis_ref* getY() const;
   virtual std::string getType() const override;
   virtual bool next() const override;

    virtual const plane_hit * getHit() const override;


  private:
    const axis_ref* m_x = nullptr;
    const axis_ref* m_y = nullptr;
    std::shared_ptr<plane_hit> m_hit = std::make_shared<plane_hit>(0,0);

  };

}

#endif // planeX_Y_hits_h__
