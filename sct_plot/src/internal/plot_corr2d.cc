#include "internal/plot_corr2d.hh"
#include "internal/plane.hh"

namespace sct_corr{
  plot_corr2d::plot_corr2d(const s_plot_prob& prob_plot) : plot2d(prob_plot)
  {

  }

  void plot_corr2d::processEventStart()
  {

  }

  void plot_corr2d::processEventEnd()
  {

  }

  void plot_corr2d::ProcessEvent()
  {
    processEventStart();
    while (m_y->next()){
      while (m_x->next()){
        //      std::cout << m_x->get()<< "  "<< m_y->get() << std::endl;
        processHit(m_x->get(), m_y->get());
      }
    }
    processEventEnd();
  }
}