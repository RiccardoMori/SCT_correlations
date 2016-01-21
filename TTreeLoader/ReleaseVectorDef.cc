#include "ReleaseVectorDef.h"
#include <string>
#include <vector>
#include "TTree.h"
#include <iostream>


class ReleaseVectorDef::vectorImpl{

public:
  vectorImpl(const char* name):m_name(name){}
  bool load_from_TTree(TTree* tree){
    if (!tree)
    {
      std::cout << "unable to load \"" << m_name << "\" from TTree. Tree is empty " << std::endl;
      return false;
    }
    if (m_vec)
    {
      std::cout << "unable to load \"" << m_name << "\" from TTree \"" << tree->GetName() << "\" vector is already set" << std::endl;
      return false;
    }
    tree->SetBranchAddress(m_name.c_str(), &m_vec);
    if (!m_vec)
    {
      std::cout << "unable to load \"" << m_name << "\" from TTree \"" << tree->GetName() << "\"" << std::endl;
      return false;
    }

    return true;
  }

  bool push2TTree(TTree* tree){
    if (!tree)
    {
      std::cout << "unable to push \"" << m_name << "\" to TTree. Tree is empty " << std::endl;
      return false;
    }


    m_vec = new std::vector<double>();
    tree->Branch(m_name.c_str(), &m_vec);

    return true;
  }
  ~vectorImpl(){
    if (!m_readFromTTree)
    {
      delete m_vec;
    }
  }
  std::vector<double> *m_vec = NULL;
  std::string m_name;
  bool m_readFromTTree=true;
};


ReleaseVectorDef::ReleaseVectorDef(const char* name) :m_vecImpl(std::make_shared<vectorImpl>(name))
{

}

bool ReleaseVectorDef::load_from_TTree(TTree* tree)
{
  return m_vecImpl->load_from_TTree(tree);
}

bool ReleaseVectorDef::push2TTree(TTree* tree)
{
  return m_vecImpl->push2TTree(tree);
}

double ReleaseVectorDef::at(size_t i) const
{
  if (!m_vecImpl->m_vec)
  {
    std::cout << "vector not set! name: " << m_vecImpl->m_name << std::endl;
    return 0;
  }

  return m_vecImpl->m_vec->at(i);
}

size_t ReleaseVectorDef::size() const
{
  if (!m_vecImpl->m_vec)
  {
    std::cout << "vector not set! name: " << m_vecImpl->m_name << std::endl;
    return 0;
  }
  return m_vecImpl->m_vec->size();
}

void ReleaseVectorDef::push_back(double val)
{
  if (!m_vecImpl->m_vec)
  {
    std::cout << "vector not set! name: " << m_vecImpl->m_name << std::endl;
    return ;
  }

  m_vecImpl->m_vec->push_back(val);

}

void ReleaseVectorDef::clear()
{
  if (!m_vecImpl->m_vec)
  {
    std::cout << "vector not set! name: " << m_vecImpl->m_name << std::endl;
    return;
  }
  m_vecImpl->m_vec->clear();
}
