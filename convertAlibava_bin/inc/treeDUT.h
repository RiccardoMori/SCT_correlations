//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed Sep 09 10:44:19 2015 by ROOT version 5.34/18
// from TTree treeDUT/tree with DUT hit and cluster information
// found on file: D:/alibava/dutTree_350V.root
//////////////////////////////////////////////////////////

#ifndef treeDUT_h
#define treeDUT_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
namespace sct_corr {
class rootEvent_Track_hits;
class rootEvent_Alibava;
}
// Header file for the classes stored in the TTree if any.

// Fixed size dimensions of array or collections stored in the TTree if any.

class treeDUT {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Int_t           runnum;
   Int_t           evtnum;
   Int_t           hitnum;
   Int_t           dutID;
   Float_t         locxpos;
   Float_t         locypos;
   Float_t         locxslope;
   Float_t         locyslope;
   Int_t           clustersize;
   Int_t           clusterxchan[9];   //[clustersize]
   Int_t           clusterychan[9];   //[clustersize]
   Float_t         clustercharge[9];   //[clustersize]

   // List of branches
   TBranch        *b_runum;   //!
   TBranch        *b_evtnum;   //!
   TBranch        *b_hitnum;   //!
   TBranch        *b_dutID;   //!
   TBranch        *b_locxpos;   //!
   TBranch        *b_locypos;   //!
   TBranch        *b_locxslope;   //!
   TBranch        *b_locyslope;   //!
   TBranch        *b_clustersize;   //!
   TBranch        *b_clusterxchan;   //!
   TBranch        *b_clusterychan;   //!
   TBranch        *b_clustercharge;   //!
   

   //////////////////////////////////////////////////////////////////////////
   // output 
   //////////////////////////////////////////////////////////////////////////


   TFile          *outFile=nullptr;
   TTree          *telescope = nullptr;
   TTree          *alibava = nullptr;

   sct_corr::rootEvent_Track_hits *m_telescope=nullptr;
   sct_corr::rootEvent_Alibava    *m_alibava=nullptr;

   treeDUT(TTree *tree=0);
   treeDUT(TTree *tree,const char* outputName);
   virtual ~treeDUT();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t GetEntries();
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif
