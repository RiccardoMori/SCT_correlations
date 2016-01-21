#ifndef rootEventBase_h__
#define rootEventBase_h__
#include "internal/platform.hh"
#include <vector>
#include <string>
#include "internal/axis_ref.hh"

class TTree;

#ifdef _DEBUG
class ReleaseVectorDef;
#endif // _DEBUG


#include "S_Axis.h"
#include <memory>
#include <vector>
#include "sct_types.h"

namespace sct_corr{
  class plane;

  class DllExport  TTreeVectorExtractor{
  public:
    TTreeVectorExtractor(const sct_type::AxesName_t& name, TTree* tree);
    TTreeVectorExtractor(const sct_type::AxesName_t& name);
    TTreeVectorExtractor();
    bool push2TTree(TTree* tree);
    ~TTreeVectorExtractor();
    std::vector<double>* getVec() const;
    sct_type::AxesName_t getName() const;
    void reset();
#ifdef _DEBUG
    void loadFromVector();
    void PushToVector() const;
    std::shared_ptr<ReleaseVectorDef> m_vecRel;
#endif // _DEBUG

  private:
    sct_type::AxesName_t m_name = sct_type::AxesName_t("");
    std::shared_ptr<std::vector<double>> m_owend_vector;
    std::vector<double>* m_vec=NULL;
    
  };


  class DllExport rootEventBase{
  public:
    rootEventBase(TTree* tree);
    rootEventBase(const sct_type::collectionName_t&  collectionName, const  std::vector<sct_type::AxesName_t>& axis_list);

    rootEventBase();
    virtual ~rootEventBase(){}
    virtual std::shared_ptr<plane> createPlane(const sct_type::ID_t& ID);
    virtual void Save2Tree(TTree* outputTree);
    std::vector<double>* getData(const sct_type::AxesName_t& name) const;

    std::vector<double>* getData(axis_def ax) const;
    std::vector<sct_type::AxesName_t> getDataNames() const;
    const sct_type::collectionName_t& getName() const;
    virtual void reset();
#ifdef _DEBUG
    void loadFromVector();
    void PushToVector() const;
#endif
    int& get_event_nr();
  private:
    sct_type::collectionName_t m_name =  sct_type::collectionName_t("");
    std::vector<TTreeVectorExtractor> m_data;
    std::shared_ptr<int> m_event_nr_ownd;
    int* m_event_nr = NULL;
  };

  class  rootEventBaseAxis :public axis_ref{
  public:
    rootEventBaseAxis(const rootEventBase& ev, axis_def ax,double id_);
    rootEventBaseAxis(const rootEventBase& ev, const sct_type::AxesName_t& ax_name, const sct_type::ID_t& id_);
    virtual bool next() const override;
    virtual double get() const override;
    const sct_type::AxesName_t&  getName() const;
    std::vector<double>* m_axis;
    std::vector<double>* m_ID;
    const sct_type::ID_t m_planeID = sct_type::ID_t(0);
    mutable int m_curr = -1;
    sct_type::AxesName_t m_name;
  };


  class  rootEventBaseAxisCollection{
  public:

    rootEventBaseAxisCollection(const rootEventBase& ev, const sct_type::ID_t& id_);
    const axis_ref* getAxis(axis_def ax) const;
    const axis_ref* getAxis(const sct_type::AxesName_t& axisName) const;
    bool next() const;
  private:
    std::vector<rootEventBaseAxis> m_axis;
  };
}

#endif // rootEventBase_h__
