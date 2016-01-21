#ifndef rootEventRunOutput_h__
#define rootEventRunOutput_h__
#include "sct_events/rootEventBase.hh"
#include <memory>
#include "internal/platform.hh"

class TTree;
namespace sct_corr{
  class DllExport rootEventRunOutput :public rootEventBase{
  public:
    rootEventRunOutput(TTree* tree);
    rootEventRunOutput(const sct_type::collectionName_t& collectionName);
    rootEventRunOutput();
    virtual ~rootEventRunOutput(){}
    virtual std::shared_ptr<plane> createPlane(const sct_type::ID_t& ID) override;
    virtual void Save2Tree(TTree* outputTree) override;
    void set_TotalNumOfEvents(double numOfEvents);
    void set_Total_efficiency(double effi);
    void set_Error_efficiency(double error_effi); 
    void set_Threshold(double thr);
    void set_rotation(double rotation_);
    void set_HV(double HV__);
    void set_residual(double res);
    void set_offset(double offset);
    void set_RunNumber(double RunNum);
    virtual void reset();
  private:
    std::shared_ptr<double> m_rotation;
    std::shared_ptr<double> m_totalNumOfEvents ;
    std::shared_ptr<double> m_total_efficiency ;
    std::shared_ptr<double> m_error_efficiency ; 
    std::shared_ptr<double> m_residual ;
    std::shared_ptr<double> m_offset ;
    std::shared_ptr<double> m_Threshold ;
    std::shared_ptr<double> m_RunNumber ;
    std::shared_ptr<double> m_HV ;
    std::vector<double> *m_x = nullptr, // these are just shortcuts 
      *m_y = nullptr,
      *m_Occupancy=nullptr,
      *m_NumOfEvents=nullptr,
      *m_id = nullptr;

  };



}


#endif // rootEventRunOutput_h__
