#ifndef plot_corr2d_h__
#define plot_corr2d_h__
#include "internal/plot2d.hh"
namespace sct_corr{
  class plot_corr2d :public plot2d{
  public:
    plot_corr2d(const s_plot_prob& = "");
    virtual void processEventStart();
    virtual void processHit(double x, double y) = 0;
    virtual void processEventEnd();
    virtual void ProcessEvent() override final;
  };

}
#endif // plot_corr2d_h__
