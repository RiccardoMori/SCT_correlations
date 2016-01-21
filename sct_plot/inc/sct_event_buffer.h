#ifndef sct_event_buffer_h__
#define sct_event_buffer_h__

#include <vector>

#include "Rtypes.h"
#include <memory>
#include "internal/axis_ref.hh"
#include "sct_plots.h"

#include "sct_events/rootEventBase.hh"
#include "sct_types.h"
#include <utility>


class TFile;
namespace sct_corr{



  class sct_event_buffer{
  public:
    void set(const sct_type::collectionName_t& name, rootEventBase* ev);
    bool get(const sct_type::collectionName_t& name, rootEventBase* ev);
    rootEventBase* get(const sct_type::collectionName_t& name);
    TFile* getOutputFile();
    void setOutputFile(TFile* file);
    bool IsCollection(const sct_type::collectionName_t& name);
    void reset();


    using map_t = std::pair<sct_type::collectionName_t, rootEventBase>;
    std::vector<map_t> m_events;
    TFile* m_outputFile = nullptr;
  };
}
#endif // sct_event_buffer_h__
