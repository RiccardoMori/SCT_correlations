#include "processors/modulo.hh"

#include "internal/plot_hit2d.hh"




namespace sct_corr {
class modulo_processor :public plot_hit2d {

public:
  modulo_processor(const sct_type::modulo_t& modulo_value, axis_def mod_axis, const s_plot_prob& plot_prob_) 
    :plot_hit2d(plot_prob_),
    m_ax(mod_axis),
    m_mod(modulo_value) 
  {
  
  
  }
  virtual void processHit(double x, double y) override;

  virtual s_plane_collection getOutputcollection() override;

  virtual const char* getType() const override;
  const axis_def m_ax;
  const sct_type::modulo_t m_mod;

};

void modulo_processor::processHit(double x, double y) {
  if (m_ax==x_axis_def)
  {
    x = fmod(x, necessary_CONVERSION(m_mod));
  }else if (m_ax==y_axis_def)
  {
    y = fmod(y, necessary_CONVERSION(m_mod));
  }
  pushHit(x, y);
}

s_plane_collection modulo_processor::getOutputcollection() {
  s_plane_collection ret;
  ret.m_planes.push_back(std::make_pair(std::string("modulo"), plane_def(getOutputName(), sct_type::ID_t(0))));
  return ret;
}

const char* modulo_processor::getType() const {
  return "modulo_____";
}

}




namespace sct_processor {

modulo::modulo(const sct_corr::plane_def& planeA, const sct_type::modulo_t& modulo_value, axis_def mod_axis, const s_plot_prob& plot_prob_ /*= "" */) {
  auto pl=planeA.get_plot();
  m_planes = pl->addPlot(sct_corr::S_plot(new sct_corr::modulo_processor(modulo_value, mod_axis, plot_prob_)), planeA);
}

sct_corr::plane_def modulo::getModulo() const {
  return m_planes.getByType("modulo")();
}
}