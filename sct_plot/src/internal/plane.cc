#include "internal/plane.hh"

#include "treeCollection.h"
#include <iostream>
#include "internal/plane_hit.hh"
#include "internal/exceptions.hh"
#include "factoryDef.hh"
#include <algorithm>
namespace sct_corr{


  plane::plane(Parameter_ref pl) : m_axis(*pl.buffer, pl.ID), m_ID(pl.ID)
  {

  }






  std::unique_ptr<plane> plane::create(const sct_type::ID_t& ID, rootEventBase* buffer) {
    const plane_struct pl(ID, buffer);
    plane::MainType type = buffer->getDataNames();
    return Class_factory_Utilities::Factory<plane>::Create( type,pl);
  }



 const axis_ref* plane::getAxis(axis_def ax) const
  {
    return m_axis.getAxis(ax);
  }

 const axis_ref* plane::getAxis(const sct_type::AxesName_t& name) const {
   return m_axis.getAxis(name);
 }


 bool isPermutation(const std::vector<sct_type::AxesName_t>& a1, const std::vector<sct_type::AxesName_t>& a2) {
 
   
   if (a1.size() != a2.size())
 {
   return false;
 }
   std::vector<std::string> m_dumm1, m_dumm2;
   for (auto&e : a1) {
     m_dumm1.push_back(Un_necessary_CONVERSION(e));
   }
   for (auto&e : a2) {
     m_dumm2.push_back(Un_necessary_CONVERSION(e));
   }

   return std::is_permutation(m_dumm1.begin(), m_dumm1.end(), m_dumm2.begin());
 }
 bool operator==(const axis_names& a1, const axis_names& a2) {
   return isPermutation(a1.m_names, a2.m_names);
 }

 bool operator==(const axis_names& a1, const std::vector<sct_type::AxesName_t>& a2) {
   return isPermutation(a1.m_names, a2);
 }

 bool operator==(const std::vector<sct_type::AxesName_t>& a1, const axis_names& a2) {
   return isPermutation(a1, a2.m_names);
 }


 std::string to_string(const axis_names& ax) {

   std::string ret;

   for (auto &e:ax.m_names )   {

     ret += " " + necessary_CONVERSION(e);
   }
   return ret;
 }

 bool operator<(const axis_names&  a1, const axis_names& a2) {
   std::vector<std::string> m_dumm1, m_dumm2;
   for (auto&e:a1.m_names)
   {
     m_dumm1.push_back(Un_necessary_CONVERSION(e));
   }
   for (auto&e : a2.m_names) {
     m_dumm2.push_back(Un_necessary_CONVERSION(e));
   }

   return m_dumm1 < m_dumm2;
 }
#ifdef _DEBUG


 std::string operator+(const axis_names& a1, const sct_type::AxesName_t& str) {
   std::string ret;

   for (const auto& e :a1.m_names)
   {
     ret += " " + necessary_CONVERSION(e)+ " ";
   }

   ret += necessary_CONVERSION(str);
   return ret;
 }

 std::string operator+(const sct_type::AxesName_t& str, const axis_names& a1) {
   std::string ret=necessary_CONVERSION(str);

   for (const auto& e : a1.m_names) {
     ret += " " + necessary_CONVERSION(e)+ " ";
   }

   return ret;
 }

#endif // _DEBUG

 std::string operator+(const axis_names& a1, const  std::string& str) {
   std::string ret;

   for (const auto& e : a1.m_names) {
     ret += " " + necessary_CONVERSION(e)+ " ";
   }

   ret += str;
   return ret;
 }

 std::string operator+(const  std::string& str, const axis_names& a1) {
   std::string ret = str;

   for (const auto& e : a1.m_names) {
     ret += " " + necessary_CONVERSION(e)+ " ";
   }

   return ret;
 }

 registerBaseClassDef(plane);
}

