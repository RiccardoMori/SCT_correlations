#ifndef plane_h__
#define plane_h__

#include <vector>

#include "internal/axis_ref.hh"
#include "sct_plots.h"
#include "plane_hit.hh"
#include "sct_events/rootEventBase.hh"
#include "factory.hh"

#define registerPlane(DerivedFileWriter,ID)  registerClass(plane,DerivedFileWriter,ID)

namespace sct_corr{

  struct plane_struct {
  public:
    plane_struct(const sct_type::ID_t& ID_, rootEventBase* buffer_) :ID(ID_), buffer(buffer_) {}
    const  sct_type::ID_t ID;
    const rootEventBase* buffer;

  };
  class axis_names {
  public:
    axis_names() {}
    axis_names(const std::vector<sct_type::AxesName_t>& names) :m_names(names) {}
    std::vector<sct_type::AxesName_t> m_names;

  };
  std::string to_string(const axis_names& ax);
  bool operator==(const axis_names& a1, const axis_names& a2);
  bool operator==(const axis_names& a1, const std::vector<std::string>& a2);
  bool operator==(const std::vector<sct_type::AxesName_t>& a1, const axis_names& a2);
  bool operator<(const axis_names&  a1, const axis_names& a2);


  std::string operator+(const axis_names& a1, const sct_type::AxesName_t& str);
  std::string operator+(const sct_type::AxesName_t& str, const axis_names& a1);
  std::string operator+(const  std::string& str, const axis_names& a1);
  std::string operator+(const axis_names& a1, const  std::string& str);
  class plane{
  public:


    using MainType = axis_names;
    using Parameter_t = plane_struct;
    using Parameter_ref = const Parameter_t&;


    plane(Parameter_ref pl);
    static std::unique_ptr<plane> create(const sct_type::ID_t& ID, rootEventBase* buffer);
    virtual bool next() const =0;
    
    virtual const plane_hit * getHit() const =0;
    virtual std::string getType() const = 0;

    

   const axis_ref* getAxis(const sct_type::AxesName_t& name) const;
   const axis_ref* getAxis(axis_def ax) const;
  protected:
    rootEventBaseAxisCollection m_axis;
  private:

    

    const sct_type::ID_t m_ID = sct_type::ID_t(0);


  };

}
#endif // plane_h__
