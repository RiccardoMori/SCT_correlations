#include "sct_plots.h"
#include <iostream>
S_plot_def::S_plot_def(const char* type, const char* name) :m_type(type), m_name(name)
{

}

S_plot_def::S_plot_def()
{
  std::cout << "[S_plot_def] unsupported default constructor do not call \n";
}

void S_plot_def::setParameter(const char* tag, const char* value)
{
  m_tags[tag] = value;
}

void S_plot_def::setParameter(const std::string & tag, const std::string& value)
{
  m_tags[tag] = value;
}

const char * S_plot_def::getParameter(const char* tag, const char* default_value)
{
  auto it = m_tags.find(tag);
  if (it == m_tags.end())
  {
    return default_value;
  }

  return it->second.c_str();
}

std::string S_plot_def::getParameter(const std::string & tag, const std::string & default_value)
{
  auto it = m_tags.find(tag);
  if (it == m_tags.end())
  {
    return default_value;
  }

  return it->second;
}
