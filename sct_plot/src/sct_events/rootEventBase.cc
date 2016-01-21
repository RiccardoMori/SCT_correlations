#include "sct_events/rootEventBase.hh"
#include "TTree.h"
#include <iostream>
#include "S_Axis.h"
#include "internal/plane.hh"
#include "internal/planeX_Y_hits.hh"
#ifdef  _DEBUG
#include "ReleaseVectorDef.h"
#endif //  _DEBUG


namespace sct_corr{




#ifdef _DEBUG   
TTreeVectorExtractor::TTreeVectorExtractor(const sct_type::AxesName_t& name, TTree* tree) :
    m_vecRel(std::make_shared<ReleaseVectorDef>(necessary_CONVERSION(name).c_str())),
    m_name(name)
   
  {


    m_vecRel->load_from_TTree(tree);
    m_owend_vector= std::make_shared<std::vector<double> >();
    m_vec = m_owend_vector.get();
  }


TTreeVectorExtractor::TTreeVectorExtractor(const sct_type::AxesName_t& name) :m_vecRel(std::make_shared<ReleaseVectorDef>(necessary_CONVERSION(name).c_str())), m_name(name)
  {
    m_owend_vector = std::make_shared<std::vector<double> >();
    m_vec = m_owend_vector.get();
  }

  TTreeVectorExtractor::TTreeVectorExtractor() 
  {

  }

  TTreeVectorExtractor::~TTreeVectorExtractor()
  {

    
  }
  bool TTreeVectorExtractor::push2TTree(TTree* tree)
  {
    return m_vecRel->push2TTree(tree);
  }


  void TTreeVectorExtractor::loadFromVector()
  {


    m_vec->clear();
    for (size_t i = 0; i < m_vecRel->size();++i)
    {
      m_vec->push_back(m_vecRel->at(i));
    }
  }


  void TTreeVectorExtractor::PushToVector() const
  {
    m_vecRel->clear();
    for (auto&e:*m_vec)
    {
      m_vecRel->push_back(e);
    }
  }
  std::vector<double>* TTreeVectorExtractor::getVec() const
  {
    return m_vec;
  }



#else  //release

TTreeVectorExtractor::TTreeVectorExtractor(const sct_type::AxesName_t&  name, TTree* tree) :
    m_name(name)
   {

    tree->SetBranchAddress(necessary_CONVERSION(m_name).c_str(), &m_vec);


  }


TTreeVectorExtractor::TTreeVectorExtractor(const sct_type::AxesName_t&  name) : m_name(name)
  {
    m_owend_vector = std::make_shared<std::vector<double>>();
    m_vec = m_owend_vector.get();
  }
  TTreeVectorExtractor::~TTreeVectorExtractor()
  {


  }
  bool TTreeVectorExtractor::push2TTree(TTree* tree)
  {

    if (!tree)
    {
      std::cout << "unable to push \"" << necessary_CONVERSION(m_name) << "\" to TTree. Tree is empty " << std::endl;
      return false;
    }

    tree->Branch(necessary_CONVERSION(m_name).c_str(), &m_vec);

    return true;
  }



  std::vector<double>* TTreeVectorExtractor::getVec() const
  {
    return m_vec;
  }



#endif


  sct_type::AxesName_t TTreeVectorExtractor::getName() const
  {
    return m_name;
  }

  void TTreeVectorExtractor::reset()
  {
    if (!m_owend_vector)
    {
      std::cout << "unable to clear not owned vector \n";
    }
    m_vec->clear();
    
  }

  rootEventBase::rootEventBase(TTree* tree) :m_name(tree->GetName())
  {
    
    auto entries =  tree->GetListOfBranches()->GetEntries();

    for (int i = 0; i < entries; i++)
    {
      auto br = dynamic_cast<TBranch*>( tree->GetListOfBranches()->At(i));
      if (!br)
      {
        continue;
      }
      auto className = br->GetClassName();
      if (strcmp("vector<double>" ,className)!=0)
      {
        continue;
      }
      
      m_data.push_back(TTreeVectorExtractor(sct_type::AxesName_t(br->GetName()), tree));
    }
    m_event_nr_ownd = std::make_shared<int>(0);
   
    tree->SetBranchAddress("event_nr", m_event_nr_ownd.get());

  }

  rootEventBase::rootEventBase(const sct_type::collectionName_t& collectionName,const  std::vector<sct_type::AxesName_t>& axis_list) :m_name(collectionName)
  {
    for (auto&e : axis_list)
    {
      m_data.emplace_back(e);
    }
    m_event_nr_ownd = std::make_shared<int>(0);
    m_event_nr = m_event_nr_ownd.get();
  }


  rootEventBase::rootEventBase()
  {

  }

  std::shared_ptr<plane> rootEventBase::createPlane(const sct_type::ID_t& ID)
  {
    return  std::shared_ptr<plane>(std::move(plane::create(ID, this)));
  }

  void rootEventBase::Save2Tree(TTree* outputTree)
  {
    for (auto& e:m_data)
    {
      e.push2TTree(outputTree);
    }
    outputTree->Branch("event_nr", m_event_nr_ownd.get());
  }





  std::vector<double>* rootEventBase::getData(axis_def ax) const
  {
    return getData(sct_type::AxesName_t(axis2String(ax)));
  }

  std::vector<double>* rootEventBase::getData(const sct_type::AxesName_t& name) const {
    for (auto&e : m_data) {
      if (Un_necessary_CONVERSION(name) == Un_necessary_CONVERSION(e.getName())) { // the opperator== is missing for this class 
        return e.getVec();
      }

    }
    return nullptr;
  }

  std::vector<sct_type::AxesName_t> rootEventBase::getDataNames() const
  {
    std::vector<sct_type::AxesName_t> ret;
    for (auto&e : m_data)
    {
      ret.push_back(e.getName());
    }
    return ret;
  }

  const sct_type::collectionName_t& rootEventBase::getName() const
  {
    return m_name;
  }

  void rootEventBase::reset()
  {
    for (auto&e : m_data)
    {
      e.reset();
    }
  }

#ifdef _DEBUG
  void rootEventBase::loadFromVector()
  {
    for (auto& e:m_data)
    {
      e.loadFromVector();
    }
  }

  void rootEventBase::PushToVector() const
  {
    for (auto& e:m_data)
    {
      e.PushToVector();
    }
  }



#endif // _DEBUG

  int& rootEventBase::get_event_nr()
  {
    return *m_event_nr;
  }

  rootEventBaseAxis::rootEventBaseAxis(const rootEventBase& ev, axis_def ax, double id_) :m_planeID(id_), m_name(axis2String(ax))
  {
    m_axis = ev.getData(m_name);
    m_ID = ev.getData(sct_type::AxesName_t("ID"));
  }

  rootEventBaseAxis::rootEventBaseAxis(const rootEventBase& ev, const sct_type::AxesName_t& ax_name, const sct_type::ID_t& id_) :m_planeID(id_), m_name(ax_name)
  {
    m_axis = ev.getData(m_name);
    m_ID = ev.getData(sct_type::AxesName_t("ID"));
  }



  bool rootEventBaseAxis::next() const
  {
    do{
      if (++m_curr >= static_cast<int>(m_ID->size()))
      {
        m_curr = -1;
        return false;
      }

    } while ((m_ID->at(m_curr) != necessary_CONVERSION(m_planeID)));

    return true;
  }

  double rootEventBaseAxis::get() const
  {
    return m_axis->at(m_curr);
  }

  const sct_type::AxesName_t&  rootEventBaseAxis::getName() const
  {
    return m_name;
  }

  rootEventBaseAxisCollection::rootEventBaseAxisCollection(const rootEventBase& ev, const sct_type::ID_t& id_)
  {
    auto names = ev.getDataNames();
    for (auto& e : names){
      m_axis.emplace_back(ev, e, id_);
    }
  }

  const axis_ref* rootEventBaseAxisCollection::getAxis(axis_def ax) const
  {
    return getAxis(sct_type::AxesName_t(axis2String(ax)));
  }

  const axis_ref* rootEventBaseAxisCollection::getAxis(const sct_type::AxesName_t&  axisName) const
  {
    for (auto& e:m_axis)
    {
      if (Un_necessary_CONVERSION(axisName) == Un_necessary_CONVERSION(e.getName()))
      {
        return dynamic_cast<const axis_ref*>(&e);
      }
    }

    return nullptr;
  }

  bool rootEventBaseAxisCollection::next() const
  {
    bool first = true;
    int num = 0;
    for (auto&e:m_axis)
    {
      if (e.next())
      {
        ++num;
      }
    }

    if (num==0|| num == m_axis.size())
    {
      return num != 0;
    }

    std::cout << "vector have different length" << std::endl;
    return false;
  }

}