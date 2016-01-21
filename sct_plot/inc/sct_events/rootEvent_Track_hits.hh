#ifndef rootEvent_Track_hits_h__
#define rootEvent_Track_hits_h__
#include <memory>

#include "rootEventBase.hh"
#include "internal/plane_hit.hh"
#include "internal/platform.hh"

namespace sct_corr{
  class DllExport rootEvent_Track_hits :public rootEventBase{
  public:
    rootEvent_Track_hits(const sct_type::collectionName_t& name);
    rootEvent_Track_hits();

    ~rootEvent_Track_hits();
    virtual std::shared_ptr<plane> createPlane(const sct_type::ID_t& ID) override;
    void push_Hit(const plane_hit& h, double ID=0);

    void push_Hit(const track_hits& h, double ID=0);
  private:
    std::vector<double> *m_x     = nullptr, // these are just shortcuts 
                        *m_y     = nullptr,
                        *m_id    = nullptr,
                        *m_chi2  = nullptr,
                        *m_phi   = nullptr,
                        *m_theta = nullptr,
                        *m_ndf = nullptr;
  }; 
}
#endif // rootEvent_Track_hits_h__
