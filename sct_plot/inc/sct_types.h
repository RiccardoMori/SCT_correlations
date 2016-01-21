#ifndef sct_types_h__
#define sct_types_h__
#include <string>
#include "internal/platform.hh"
#ifdef _DEBUG


#define necessary_CONVERSION(x) x.value
#define Un_necessary_CONVERSION(x) x.value

#define TYPE_CLASS(name,type) \
class DllExport name { \
public: \
  explicit name(const type& param_) :value(param_) {}\
  type value; \
}
#define TYPE_CLASS_PTR(name,type) \
class DllExport name { \
public: \
  explicit name(type param_) :value(param_) {}\
  type value; \
}

#else 
#define necessary_CONVERSION(x) x
#define Un_necessary_CONVERSION(x) x

#define TYPE_CLASS(name,type) \
using name = type;
#define TYPE_CLASS_PTR(name,type) \
using name = type;

#endif // _DEBUG
namespace sct_type {

  TYPE_CLASS(modulo_x, double);
  TYPE_CLASS(modulo_y, double);
  TYPE_CLASS(modulo_t, double);


  TYPE_CLASS(ID_t, double);
  TYPE_CLASS(stripNr_t, int);
  TYPE_CLASS(residualCut_t, double);
  TYPE_CLASS(rot_angle_t, double);

  TYPE_CLASS(move_t, double);
  TYPE_CLASS(collectionName_t, std::string);
  TYPE_CLASS(AxesName_t, std::string);
  TYPE_CLASS(procent_t, double);

}

#endif // sct_types_h__
