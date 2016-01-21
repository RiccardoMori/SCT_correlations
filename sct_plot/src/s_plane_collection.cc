
#include <iostream>
#include "plane_def.h"
#include "s_DrawOption.h"
#include "plane_def.h"
#include "internal/exceptions.hh"


sct_corr::plane_def s_plane_collection::get(Int_t i) const
{
  if (i< static_cast<Int_t>(m_planes.size()))
  {
    return m_planes[i].second;
  }
  SCT_THROW("plane not found");
  std::cout << "[s_plane_collection] out of range" << std::endl;
  return sct_corr::error_plane_def();
}

s_plane_collection s_plane_collection::getByType(const char* type) const
{


  s_plane_collection ret;
  for (auto& e : m_planes)
  {
    if (std::string(e.first) == std::string(type)) {
      ret.push_back(e.first.c_str(), e.second);
    }

  }
//   if (ret.m_planes.empty()){
//     std::cout << "[s_plane_collection] unknown type = \"" << type << "\"" << std::endl;
// 
//   }
  return ret;




}
s_plane_collection s_plane_collection::getByName(const sct_type::collectionName_t& name) const
{
  s_plane_collection ret;
  for (auto& e : m_planes)
  {
    if (Un_necessary_CONVERSION(e.second.getName())== Un_necessary_CONVERSION(name)) {
      ret.push_back(e.first.c_str(), e.second);
    }

  }
//   if (ret.m_planes.empty()){
//     std::cout << "[s_plane_collection] unknown name = \"" << name << "\"" << std::endl;
// 
//   }
  return ret;
}
s_plane_collection s_plane_collection::get(const char* nameOrType) const
{
  auto pl = getByName(sct_type::collectionName_t(nameOrType));
  auto pl_type = getByType(nameOrType);
 
  return pl + pl_type;

}

sct_corr::plane_def s_plane_collection::get(const sct_type::collectionName_t&  name, const char* type) const
{
  auto pl=getByName(name).getByType(type);
  if (pl.m_planes.empty())
  {
      std::cout << "[s_plane_collection] unknown name = \"" << necessary_CONVERSION(name) << "\"" << std::endl;
      return sct_corr::error_plane_def();
  }
  if (pl.m_planes.size()>1)
  {
    std::cout << "[s_plane_collection] multiple combinations of type = " << type << " and name = " << necessary_CONVERSION(name) << std::endl;
  }
  return pl.get(0);
}
sct_corr::plane_def s_plane_collection::operator()() const
{
  
  return get(0);
}
s_plane_collection::s_plane_collection(const sct_corr::plane_def& plane_)
{
  push_back(plane_);
}
void s_plane_collection::set_plot_collection(std::weak_ptr<sct_corr::plot_collection> plot_collection___) {
  for (auto & e : m_planes) {
    e.second.set_plot_collection(plot_collection___);
  }
}
const char* s_plane_collection::getName(Int_t i) const
{
  if (i< static_cast<Int_t>(m_planes.size()))
  {
    return m_planes.at(i).first.c_str();
  }
  return "collection not found";
}


void s_plane_collection::showNames() const
{
  for (auto& e: m_planes)
  {
    std::cout << e.first<< std::endl;
  }
}
Int_t s_plane_collection::size() const
{
  return m_planes.size();
}

void s_plane_collection::push_back(const char* name, const sct_corr::plane_def& pl)
{
  m_planes.push_back(std::make_pair(std::string(name), pl));
}

void s_plane_collection::push_back(const sct_corr::plane_def& pl)
{
  push_back(Un_necessary_CONVERSION(pl.getName()).c_str(), pl);
}






 s_plane_collection operator+(s_plane_collection pl1, const s_plane_collection& pl2)
{
  pl1.m_planes.insert(pl1.m_planes.end(), pl2.m_planes.begin(), pl2.m_planes.end());


  return pl1;
}

 s_plane_collection operator+(s_plane_collection pl1, const sct_corr::plane_def& pl2)
 {
   pl1.push_back(pl2);
   return pl1;
 }






 s_plane_collection operator+(const sct_corr::plane_def& pl1, const sct_corr::plane_def& pl2)
  {
    s_plane_collection ret;
    ret.push_back(pl1);
    ret.push_back(pl2);
    return ret;
  }

