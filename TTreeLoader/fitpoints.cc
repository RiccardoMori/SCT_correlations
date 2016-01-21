#define fitpoints_cxx
#include "fitpoints.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>



fitpoints::fitpoints(TTree *tree) : fChain(0)
{
  // if parameter tree is not specified (or zero), connect the file
  // used to generate this class and read the Tree.
  if (tree == 0) {
    TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("C:/Users/Argg/OneDrive/alibava/fitter.root");
    if (!f || !f->IsOpen()) {
      f = new TFile("C:/Users/Argg/OneDrive/alibava/fitter.root");
    }
    f->GetObject("fitpoints", tree);

  }
  Init(tree);
}

fitpoints::~fitpoints()
{
  if (!fChain) return;
  delete fChain->GetCurrentFile();
}

Int_t fitpoints::GetEntry(Long64_t entry)
{
  // Read contents of entry.
  if (!fChain) return 0;
  return fChain->GetEntry(entry);
}

Int_t fitpoints::GetEntries() const
{
  return fChain->GetEntries();
}

Long64_t fitpoints::LoadTree(Long64_t entry)
{
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

void fitpoints::Init(TTree *tree)
{
  // The Init() function is called when the selector needs to initialize
  // a new tree or chain. Typically here the branch addresses and branch
  // pointers of the tree will be set.
  // It is normally not necessary to make changes to the generated
  // code, but the routine can be extended by the user if needed.
  // Init() will be called many times when running on PROOF
  // (once per file to be processed).

  // Set object pointer
  ID = 0;
  x = 0;
  y = 0;
  // Set branch addresses and branch pointers
  if (!tree) return;
  fChain = tree;
  fCurrent = -1;
  fChain->SetMakeClass(1);

  fChain->SetBranchAddress("ID", &ID, &b_ID);
  fChain->SetBranchAddress("x", &x, &b_x);
  fChain->SetBranchAddress("y", &y, &b_y);
  fChain->SetBranchAddress("event_nr", &event_nr, &b_event_nr);
  Notify();
}

Bool_t fitpoints::Notify()
{
  // The Notify() function is called when a new file is opened. This
  // can be either for a new TTree in a TChain or when when a new TTree
  // is started when using PROOF. It is normally not necessary to make changes
  // to the generated code, but the routine can be extended by the
  // user if needed. The return value is currently not used.

  return kTRUE;
}

void fitpoints::Show(Long64_t entry)
{
  // Print contents of entry.
  // If entry is not specified, print current entry
  if (!fChain) return;
  fChain->Show(entry);
}

const char* fitpoints::getName() const
{
 return fChain->GetName();
}

Int_t fitpoints::Cut(Long64_t entry)
{
  // This function may be called from Loop.
  // returns  1 if entry is accepted.
  // returns -1 otherwise.
  return 1;
}


void fitpoints::Loop()
{
  //   In a ROOT session, you can do:
  //      Root > .L fitpoints.C
  //      Root > fitpoints t
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
  for (Long64_t jentry = 0; jentry < nentries; jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    // if (Cut(ientry) < 0) continue;
    std::cout << x->size() << std::endl;
    for (auto&e : *x)
    {
      std::cout << e << " ";
    }
    std::cout << std::endl;
  }
}

Hit_output_impl::Hit_output_impl(const char* name)
{
  fChain = new TTree(name, name);
  fChain->Branch("x", &m_x);
  fChain->Branch("y", &m_y);
  fChain->Branch("ID", &m_ID);

  fChain->Branch("event_nr", &event_nr);
}

void Hit_output_impl::set(double x, double y, double id)
{
  m_x.push_back(x);
  m_y.push_back(y);
  m_ID.push_back(id);
}

void Hit_output_impl::setEventNr(Int_t eventNR)
{
  event_nr = eventNR;
}

void Hit_output_impl::fill()
{
  fChain->Fill();
  event_nr = 0;
  m_ID.clear();
  m_x.clear();
  m_y.clear();
}

Int_t Hit_output_impl::Draw(const char* axis, const char* cuts, const char * options)
{
  return  fChain->Draw(axis, cuts, options);
}

Hit_output_impl::~Hit_output_impl()
{
  delete fChain;
}
