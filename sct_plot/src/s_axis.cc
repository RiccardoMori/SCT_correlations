#include "S_Axis.h"
#include <iostream>
#include "s_plot_collection.h"


S_Axis::S_Axis(sct_type::collectionName_t collctionName, sct_type::ID_t planeID, axis_def axis) :m_collectionName(collctionName), m_planeID(planeID), m_axis(axis), m_isValid(true)
{
    
}

S_Axis::S_Axis():m_isValid(false),m_collectionName(""),m_planeID(0) {
  std::cout << "invalid axis" << std::endl;
}

S_Axis::S_Axis(sct_type::collectionName_t collctionName, sct_type::ID_t planeID, axis_def axis, std::weak_ptr<sct_corr::plot_collection> plot_collection_)
  : S_Axis(collctionName, planeID, axis) 
{
  m_plot_collection = plot_collection_;
}





std::shared_ptr<sct_corr::plot_collection> S_Axis::get_plot() const {
  if (std::shared_ptr<sct_corr::plot_collection> ret = m_plot_collection.lock()) {

    return ret;
  } else {
    std::cout << "plot is expired\n";
  }
  return nullptr;
}

void S_Axis::set_plot_collection(std::weak_ptr<sct_corr::plot_collection> plot_collection_) {
#ifdef _DEBUG
  if (m_plot_collection.lock())
  {
    std::cout << "m_plot_collection was already defined " << std::endl;
  }

#endif // _DEBUG

  m_plot_collection = plot_collection_;
}

bool S_Axis::isValid() const {
  return m_isValid;
}

sct_type::AxesName_t axis2String(axis_def ax)
{
    if (ax==x_axis_def)
    {
      return sct_type::AxesName_t("x");
    }
    
    if (ax==y_axis_def)
    {
      return sct_type::AxesName_t("y");
    }
    if (ax==Ndf_axis_def)
    {
      return sct_type::AxesName_t("ndf");
    }
    if (ax == chi2_axis_def )
    {
      return sct_type::AxesName_t("chi2");
    }
    if (ax == phi_axis_def)
    {
      return sct_type::AxesName_t("phi");
    }
    if (ax == theta_axis_def)
    {
      return sct_type::AxesName_t("theta");
    }
    if (ax==Occupancy_axis_def)
    {
      return sct_type::AxesName_t("Occupancy");
    }
    if (ax==Occupancy_error_axis_def)//lb
    {
      return sct_type::AxesName_t("Occupancy_error");
    }
    
    if (ax==NumOfEvents_axis_def)
    {
      return sct_type::AxesName_t("NumOfEvents");
    }
    if (ax == charge_axis_def) {
      return sct_type::AxesName_t("Charge");
    }
    return sct_type::AxesName_t("");
}

sct_type::AxesName_t getIDString()
{
  return sct_type::AxesName_t("ID");
}
