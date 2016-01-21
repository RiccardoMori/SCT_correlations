//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed May 20 16:34:27 2015 by ROOT version 5.34/18
// from TTree fitpoints/fitpoints
// found on file: C:/Users/Argg/OneDrive/alibava/fitter.root
//////////////////////////////////////////////////////////

#ifndef fitpoints_h
#define fitpoints_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include <vector>

// Fixed size dimensions of array or collections stored in the TTree if any.

class fitpoints {
public:
  TTree          *fChain;   //!pointer to the analyzed TTree or TChain
  Int_t           fCurrent; //!current Tree number in a TChain

  // Declaration of leaf types
  std::vector<double>  *ID;
  std::vector<double>  *x;
  std::vector<double>  *y;
  Int_t           event_nr;

  // List of branches
  TBranch        *b_ID;   //!
  TBranch        *b_x;   //!
  TBranch        *b_y;   //!
  TBranch        *b_event_nr;   //!

  fitpoints(TTree *tree = 0);
  virtual ~fitpoints();
  virtual Int_t    Cut(Long64_t entry);
  virtual Int_t    GetEntry(Long64_t entry);
  virtual Int_t    GetEntries() const;
  virtual Long64_t LoadTree(Long64_t entry);
  virtual void     Init(TTree *tree);
  virtual void     Loop();
  virtual Bool_t   Notify();
  virtual void     Show(Long64_t entry = -1);
  const char* getName() const;
};


class Hit_output_impl{
public:
  TTree          *fChain;   //!pointer to the analyzed TTree or TChain
 

  // Declaration of leaf types
  std::vector<double>  m_ID;
  std::vector<double>  m_x;
  std::vector<double>  m_y;
  Int_t           event_nr;
  void set(double x, double y, double id);
  void setEventNr(Int_t eventNR);
  void fill();
  Int_t Draw(const char* axis, const char* cuts, const char * options);
  Hit_output_impl(const char* name);
  ~Hit_output_impl();
};
#endif



