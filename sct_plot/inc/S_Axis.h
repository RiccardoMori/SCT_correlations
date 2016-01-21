#ifndef s_axis_h__
#define s_axis_h__
#include "internal/platform.hh"
#include <memory>
#include <string>
#include "sct_types.h"

enum  axis_def  :int{
  x_axis_def,
  y_axis_def,
  chi2_axis_def,
  Ndf_axis_def,
  phi_axis_def,
  theta_axis_def,
  Occupancy_axis_def,
  Occupancy_error_axis_def,
  NumOfEvents_axis_def,
  charge_axis_def
};
sct_type::AxesName_t axis2String(axis_def ax);
sct_type::AxesName_t getIDString();
namespace sct_corr {
  class plot_collection;
}
class DllExport S_Axis{
public:
  S_Axis(sct_type::collectionName_t collctionName, sct_type::ID_t planeID, axis_def axis);

  S_Axis(sct_type::collectionName_t collctionName, sct_type::ID_t planeID, axis_def axis, std::weak_ptr<sct_corr::plot_collection> plot_collection_);
  std::shared_ptr<sct_corr::plot_collection> get_plot() const;
  void set_plot_collection(std::weak_ptr<sct_corr::plot_collection> plot_collection_);


  S_Axis();
  bool isValid() const;
  sct_type::collectionName_t  m_collectionName;
  sct_type::ID_t m_planeID;
  axis_def m_axis;
private:

  std::weak_ptr<sct_corr::plot_collection> m_plot_collection;
  bool m_isValid;

 
};


#endif // s_axis_h__
