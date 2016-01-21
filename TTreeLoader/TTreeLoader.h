#ifndef TTreeLoader_h__
#define TTreeLoader_h__

#define DllImport   __declspec( dllimport )
#define DllExport   __declspec( dllexport )

#include "TTree.h"
class fitpoints;
class Hit_output_impl;
class DllExport Hit_extractor{
public:
  

  Hit_extractor(TTree* tree);
  virtual ~Hit_extractor();


  double getX();
  double getY();
  Int_t getEventNr();
  double getID();

  void getEntry(Long64_t entry);
  Int_t GetNumberOfEntries();
  Int_t GetEvent(Long64_t eventNr);
  Int_t GetNumberOfEvents();
  const char* getName() const;
private:
  fitpoints *m_hits;
  Int_t m_entry=0;
};


class DllExport Hit_output{
public:


  Hit_output(const char* name);
  virtual ~Hit_output();


  void set(Double_t x,Double_t y,Double_t ID);
  void setEventNR(Int_t eventNr);

  void fill();
  Int_t Draw(const char* axis, const char* cuts, const char * options);
  TTree* getTTree();
private:
  Hit_output_impl *m_hits;
};
#endif // TTreeLoader_h__
