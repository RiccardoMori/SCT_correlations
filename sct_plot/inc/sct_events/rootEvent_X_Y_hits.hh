#ifndef rootEvent_X_Y_hits_h__
#define rootEvent_X_Y_hits_h__
#include "sct_events/rootEventBase.hh"
#include <vector>
#include "internal/plane_hit.hh"
namespace sct_corr{
  class DllExport rootEvent_X_Y_hits :public rootEventBase{
  public:
    rootEvent_X_Y_hits(const sct_type::collectionName_t& name);
    rootEvent_X_Y_hits();

    ~rootEvent_X_Y_hits();
    virtual std::shared_ptr<plane> createPlane(const sct_type::ID_t& ID);

    void push_Hit(const plane_hit& h, double ID=0);
  private:
    std::vector<double> *m_x=nullptr, // these are just shortcuts 
                        *m_y =nullptr,
                        *m_id =nullptr;
  };

}
#endif // rootEvent_X_Y_hits_h__