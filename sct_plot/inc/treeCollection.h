#ifndef treeCollection_h__
#define treeCollection_h__



#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include <vector>
#include "Rtypes.h"
#include "sct_event_buffer.h"
#include "sct_events/rootEventBase.hh"
#include "sct_types.h"
#ifdef _DEBUG
  class Hit_extractor;
  class Hit_output;
#endif // _DEBUG

namespace sct_corr{




  class DllExport treeCollection {
  public:


    rootEventBase m_rootBuffer;
    // Declaration of leaf types
   

    treeCollection(TTree *tree);
    treeCollection(const sct_type::collectionName_t& name, sct_event_buffer* outputBuffer);
    virtual ~treeCollection();
    virtual Int_t    GetEntry(Long64_t entry);
    virtual Long64_t    GetEntries() const;
    sct_type::collectionName_t  getName() const;
  private:

    TTree* fChain = NULL;


    sct_type::collectionName_t m_name = sct_type::collectionName_t("");
  };



  class treeCollection_ouput {
  public:

    treeCollection_ouput(rootEventBase& ev, sct_event_buffer* outputBuffer, bool save2disk);
    virtual ~treeCollection_ouput();
    void fill();
    Long64_t Draw(const char* axis, const char* cuts, const char * options);
    TTree * getTTree();

    rootEventBase m_rootBuffer;



    TTree          *m_tree=nullptr;   //!pointer to the analyzed TTree or TChain



    sct_type::collectionName_t  m_name = sct_type::collectionName_t("");
  };
}
#endif// treeCollection_h__



