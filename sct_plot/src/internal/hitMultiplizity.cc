#include "internal/plot_hit2d.hh"
namespace sct_corr{

  class hitMultiplizity :public plot_hit2d{
  public:
    hitMultiplizity(const s_plot_prob& = "");

    virtual void processEventStart();
    virtual void processHit(double x, double y);;
    virtual void processEventEnd();

    Int_t m_counter = 0;
    virtual s_plane_collection getOutputcollection();

    virtual const char* getType() const override;
  };




  hitMultiplizity::hitMultiplizity(const s_plot_prob& plot_prob) : plot_hit2d(plot_prob)
  {

  }

  void hitMultiplizity::processEventStart()
  {
    m_counter = 0;
  }

  void hitMultiplizity::processHit(double x, double y)
  {
    ++m_counter;
  }

  void hitMultiplizity::processEventEnd()
  {
    pushHit(m_counter, 0);
  }

  s_plane_collection hitMultiplizity::getOutputcollection()
  {
    s_plane_collection ret;
    ret.m_planes.push_back(std::make_pair(std::string("HitMultiplicity"), S_plane_def(getOutputName(), 0)));
    return ret;
  }

  const char* hitMultiplizity::getType() const
  {
    return "hitMultiplizity___";
  }
}
S_plot sct_plot::hitMultiplizity(const s_plot_prob& plot_prob)
{
  return S_plot(new sct_corr::hitMultiplizity(plot_prob));
}
