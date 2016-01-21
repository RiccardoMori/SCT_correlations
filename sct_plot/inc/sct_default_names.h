#ifndef sct_default_names_h__
#define sct_default_names_h__
#include "internal/platform.hh"
#include "sct_types.h"

class DllExport sct {
public:
  static sct_type::collectionName_t col_hit();
  static sct_type::collectionName_t col_zsdata_strip();
  static sct_type::collectionName_t col_zsdata_m26();
  static sct_type::collectionName_t col_zsdata_apix();
  static sct_type::collectionName_t col_fitpoints();
  static sct_type::collectionName_t col_fitpoints_local();
  static sct_type::collectionName_t col_local_hit();
  static sct_type::collectionName_t col_DUT_TTC();
  static sct_type::collectionName_t col_GBL_fitted_points();
};

#ifdef __CINT__
  #pragma link C++ class sct;
#endif

#endif // sct_default_names_h__
