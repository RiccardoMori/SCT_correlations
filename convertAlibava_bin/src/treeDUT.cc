#include "../inc/treeDUT.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include "TTree.h"
#include "sct_events/rootEvent_Track_hits.hh"
#include "sct_events/rootEvent_Alibava.hh"
#include "internal/plane_hit.hh"





treeDUT::treeDUT(TTree *tree) : fChain(0) {
  // if parameter tree is not specified (or zero), connect the file
  // used to generate this class and read the Tree.
  if (tree == 0) {
    TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("D:/alibava/dutTree_350V.root");
    if (!f || !f->IsOpen()) {
      f = new TFile("D:/alibava/dutTree_350V.root");
    }
    f->GetObject("treeDUT", tree);

  }
  Init(tree);
}

treeDUT::treeDUT(TTree *tree, const char* outputName) {
  // if parameter tree is not specified (or zero), connect the file
  // used to generate this class and read the Tree.
  if (tree == 0) {
    TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("D:/alibava/dutTree_350V.root");
    if (!f || !f->IsOpen()) {
      f = new TFile("D:/alibava/dutTree_350V.root");
    }
    f->GetObject("treeDUT", tree);

  }
  Init(tree);
  outFile = new TFile(outputName, "recreate");

  telescope = new TTree("telescope", "telescope");

  m_telescope = new sct_corr::rootEvent_Track_hits(sct_type::collectionName_t("telescope"));
  m_telescope->Save2Tree(telescope);

  alibava = new TTree("alibava", "alibava");
  m_alibava = new sct_corr::rootEvent_Alibava(sct_type::collectionName_t("alibava"));
  m_alibava->Save2Tree(alibava);
}

treeDUT::~treeDUT() {
  if (!fChain) return;

  outFile->Write();
 // delete fChain->GetCurrentFile();
}

Int_t treeDUT::GetEntry(Long64_t entry) {
  // Read contents of entry.
  if (!fChain) return 0;
  return fChain->GetEntry(entry);
}

Long64_t treeDUT::GetEntries() {
  return fChain->GetEntries();
}

Long64_t treeDUT::LoadTree(Long64_t entry) {
  // Set the environment to read one entry
  if (!fChain) return -5;
  Long64_t centry = fChain->LoadTree(entry);
  if (centry < 0) return centry;
  if (fChain->GetTreeNumber() != fCurrent) {
    fCurrent = fChain->GetTreeNumber();
    Notify();
  }
  return centry;
}

void treeDUT::Init(TTree *tree) {
  // The Init() function is called when the selector needs to initialize
  // a new tree or chain. Typically here the branch addresses and branch
  // pointers of the tree will be set.
  // It is normally not necessary to make changes to the generated
  // code, but the routine can be extended by the user if needed.
  // Init() will be called many times when running on PROOF
  // (once per file to be processed).

  // Set branch addresses and branch pointers
  if (!tree) return;
  fChain = tree;
  fCurrent = -1;
  fChain->SetMakeClass(1);

  fChain->SetBranchAddress("runnum", &runnum, &b_runum);
  fChain->SetBranchAddress("evtnum", &evtnum, &b_evtnum);
  fChain->SetBranchAddress("hitnum", &hitnum, &b_hitnum);
  fChain->SetBranchAddress("dutID", &dutID, &b_dutID);
  fChain->SetBranchAddress("locxpos", &locxpos, &b_locxpos);
  fChain->SetBranchAddress("locypos", &locypos, &b_locypos);
  fChain->SetBranchAddress("locxslope", &locxslope, &b_locxslope);
  fChain->SetBranchAddress("locyslope", &locyslope, &b_locyslope);
  fChain->SetBranchAddress("clustersize", &clustersize, &b_clustersize);
  fChain->SetBranchAddress("clusterxchan", clusterxchan, &b_clusterxchan);
  fChain->SetBranchAddress("clusterychan", clusterychan, &b_clusterychan);
  fChain->SetBranchAddress("clustercharge", clustercharge, &b_clustercharge);
  Notify();
}

Bool_t treeDUT::Notify() {
  // The Notify() function is called when a new file is opened. This
  // can be either for a new TTree in a TChain or when when a new TTree
  // is started when using PROOF. It is normally not necessary to make changes
  // to the generated code, but the routine can be extended by the
  // user if needed. The return value is currently not used.

  return kTRUE;
}

void treeDUT::Show(Long64_t entry) {
  // Print contents of entry.
  // If entry is not specified, print current entry
  if (!fChain) return;
  fChain->Show(entry);
}
Int_t treeDUT::Cut(Long64_t entry) {
  // This function may be called from Loop.
  // returns  1 if entry is accepted.
  // returns -1 otherwise.
  return 1;
}



void treeDUT::Loop()
{
//   In a ROOT session, you can do:
//      Root > .L treeDUT.C
//      Root > treeDUT t
//      Root > t.GetEntry(12); // Fill t data members with entry number 12
//      Root > t.Show();       // Show values of entry 12
//      Root > t.Show(16);     // Read and show values of entry 16
//      Root > t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   sct_corr::track_hits tel(0,0,0,0,0,0);
   sct_corr::Alibava_hit alib(0, 0, 0);
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
     m_alibava->reset();
     m_telescope->reset();
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;


      tel.x = locxpos;
      tel.y = locypos;
      tel.phi = locxslope; // not sure if correct
      tel.theta = locyslope; // not sure if correct
      tel.chi2 = 1;
      tel.ndf = 1;

      m_telescope->push_Hit(tel,0);

      for (int i = 0; i < clustersize;++i)
      {
        alib.x = clusterxchan[i];
        alib.y = clusterychan[i];
        alib.charge = clustercharge[i];
        m_alibava->push_Hit(alib,dutID);
      }
#ifdef _DEBUG
      m_telescope->PushToVector();
      m_alibava->PushToVector();

#endif // _DEBUG

      alibava->Fill();
      telescope->Fill();
   }
}
