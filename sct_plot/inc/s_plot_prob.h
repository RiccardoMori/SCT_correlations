#ifndef s_plot_prob_h__
#define s_plot_prob_h__
#include <string>
#include "Rtypes.h"

#include "internal/platform.hh"
#include "plot_save_option_def.h"
#include "sct_types.h"
class DllExport s_plot_prob{
public:
  s_plot_prob();
  s_plot_prob(const char* name);
  s_plot_prob& setName(const char* name);
  s_plot_prob& setSaveOptione(plot_save_option_def save_);
  s_plot_prob& doNotSaveToDisk();
  s_plot_prob& SaveToDisk();
  sct_type::collectionName_t getName() const;
  size_t getNameSize() const;
  plot_save_option_def getPlotSaveOption() const;

private:
  sct_type::collectionName_t m_name = sct_type::collectionName_t("");
  plot_save_option_def m_save = save_to_disk;


};


#endif // s_plot_prob_h__
