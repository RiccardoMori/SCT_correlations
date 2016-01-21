
#include "TTreeLoader.h"
#include "fitpoints.h"
#include <iostream>

Hit_extractor::Hit_extractor(TTree* tree) :m_hits(new fitpoints(tree))
{

}

Hit_extractor::~Hit_extractor()
{
  delete m_hits;
}

double Hit_extractor::getX()
{
  return m_hits->x->at(m_entry);
}

double Hit_extractor::getY()
{
  return m_hits->y->at(m_entry);
}

Int_t Hit_extractor::getEventNr()
{
  return m_hits->event_nr;
}

double Hit_extractor::getID()
{
  return m_hits->ID->at(m_entry);
}

void Hit_extractor::getEntry(Long64_t entry)
{
  m_entry = entry;
}

Int_t Hit_extractor::GetNumberOfEntries()
{
  Int_t ret = m_hits->x->size();
  
  if (ret != m_hits->y->size())
  {
    std::cout << "[Hit_extractor] vectors dont have the same size " << std::endl;
    return -1;
  }

  if (ret!= m_hits->ID->size())
  {
    std::cout << "[Hit_extractor] vectors dont have the same size " << std::endl;
    return -1;
  }
  return ret;
}

Int_t Hit_extractor::GetEvent(Long64_t eventNr)
{
  
  auto ret=m_hits->GetEntry(eventNr);
  getEntry(0);
  return ret;
}

Int_t Hit_extractor::GetNumberOfEvents()
{

  return m_hits->GetEntries();
}



const char* Hit_extractor::getName() const
{
  return m_hits->getName();
}

Hit_output::Hit_output(const char* name) :m_hits(new Hit_output_impl(name))
{

}

Hit_output::~Hit_output()
{
  delete m_hits;
}

void Hit_output::set(Double_t x, Double_t y, Double_t ID)
{
  m_hits->set(x, y, ID);
}

void Hit_output::setEventNR(Int_t eventNr)
{
  m_hits->setEventNr(eventNr);
}

void Hit_output::fill()
{
  m_hits->fill();
}

Int_t Hit_output::Draw(const char* axis, const char* cuts, const char * options)
{
 return  m_hits->Draw(axis, cuts, options);
}

TTree* Hit_output::getTTree()
{
  return m_hits->fChain;
}
