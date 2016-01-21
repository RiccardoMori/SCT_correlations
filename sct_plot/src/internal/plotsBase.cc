#include "internal/plotsBase.hh"
#include <string>

namespace sct_corr{
  int g_plot_count_intern = 0;
  plot::plot(const s_plot_prob& plot_proberties) :m_prob(plot_proberties), m_name(plot_proberties.getName())
  {
   
  }

  void plot::setParameter(const char* tag, const char * value)
  {

  }

  s_plane_collection plot::getOutputcollection()
  {
    s_plane_collection ret;
    return ret;
  }

   sct_type::collectionName_t plot::getName() const
  {
    if (necessary_CONVERSION(m_name).empty())
    {
      necessary_CONVERSION(m_name) = getType() + std::to_string(g_plot_count_intern++);
    }
    return m_name;
  }

  bool plot::getSave2disk() const
  {
    return m_prob.getPlotSaveOption()==save_to_disk;
  }

}