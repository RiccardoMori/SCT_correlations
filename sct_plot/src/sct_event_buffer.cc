#include "sct_event_buffer.h"
#include "treeCollection.h"
#include "internal/exceptions.hh"



namespace sct_corr{


void sct_event_buffer::set(const sct_type::collectionName_t& name, rootEventBase* ev)
  {
    if (IsCollection(name))
    {
      SCT_THROW(std::string("name already used. name =  ") + necessary_CONVERSION(name));
    }
    m_events.push_back(map_t(name, *ev));
  }

bool sct_event_buffer::get(const sct_type::collectionName_t& name, rootEventBase* ev)
  {
    if (!IsCollection(name))
    {
      SCT_THROW(std::string("Collection not found in buffer. Collection name ") + necessary_CONVERSION(name));
    }
    for (auto&e:m_events)
    {
      if (Un_necessary_CONVERSION(e.first) == Un_necessary_CONVERSION(name))
      {
        *ev = e.second;
      }
    }
    

    return true;
  }

rootEventBase* sct_event_buffer::get(const sct_type::collectionName_t& name) {
  if (!IsCollection(name)) {
    SCT_THROW(std::string("Collection not found in buffer. Collection name ") + necessary_CONVERSION(name));
  }
  for (auto&e : m_events) {
    if (Un_necessary_CONVERSION(e.first) == Un_necessary_CONVERSION(name)) {
     return &e.second;
    }
  }

}

TFile* sct_event_buffer::getOutputFile()
  {
    return m_outputFile;
  }

  void sct_event_buffer::setOutputFile(TFile* file)
  {
    m_outputFile = file;
  }

  bool sct_event_buffer::IsCollection(const sct_type::collectionName_t& name)
  {

    for (auto&e :m_events)
    {
      if (Un_necessary_CONVERSION(e.first) ==Un_necessary_CONVERSION(name)) {
        return true;
      }
    }
    return false;
  }

  void sct_event_buffer::reset()
  {
    m_events.clear();
  }



}