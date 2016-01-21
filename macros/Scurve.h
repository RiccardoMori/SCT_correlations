//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Thu Jul  9 12:17:18 2015 by ROOT version 6.02/05
// from TTree out/out
// found on file: DEVICE_1_ASIC_on_Position_7_Jim_350V.root
//////////////////////////////////////////////////////////

#ifndef Scurve_h
#define Scurve_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "vector"

class Scurve {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   vector<double>  *ID;
   vector<double>  *x;
   vector<double>  *y;
   vector<double>  *Occupancy;
   vector<double>  *NumOfEvents;
   Int_t           event_nr;
   Double_t        totalNumOfEvents;
   Double_t        total_efficiency;
   Double_t        residual;
   Double_t        offset;
   Double_t        Threshold;
   Double_t        RunNumber;
   Double_t        HV;

   // List of branches
   TBranch        *b_ID;   //!
   TBranch        *b_x;   //!
   TBranch        *b_y;   //!
   TBranch        *b_Occupancy;   //!
   TBranch        *b_NumOfEvents;   //!
   TBranch        *b_event_nr;   //!
   TBranch        *b_totalNumOfEvents;   //!
   TBranch        *b_total_efficiency;   //!
   TBranch        *b_residual;   //!
   TBranch        *b_offset;   //!
   TBranch        *b_Threshold;   //!
   TBranch        *b_RunNumber;   //!
   TBranch        *b_HV;   //!

   Scurve(TTree *tree=0);
   virtual ~Scurve();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef Scurve_cxx
Scurve::Scurve(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("/Users/lucreziastellabruni/Desktop/data_350/DEVICE_1_ASIC_on_Position_7_Jim_350V.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("/Users/lucreziastellabruni/Desktop/data_350/DEVICE_1_ASIC_on_Position_7_Jim_350V.root");
      }
      f->GetObject("out",tree);

   }
   Init(tree);
}

Scurve::~Scurve()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t Scurve::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t Scurve::LoadTree(Long64_t entry)
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

void Scurve::Init(TTree *tree)
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
   Occupancy = 0;
   NumOfEvents = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("ID", &ID, &b_ID);
   fChain->SetBranchAddress("x", &x, &b_x);
   fChain->SetBranchAddress("y", &y, &b_y);
   fChain->SetBranchAddress("Occupancy", &Occupancy, &b_Occupancy);
   fChain->SetBranchAddress("NumOfEvents", &NumOfEvents, &b_NumOfEvents);
   fChain->SetBranchAddress("event_nr", &event_nr, &b_event_nr);
   fChain->SetBranchAddress("totalNumOfEvents", &totalNumOfEvents, &b_totalNumOfEvents);
   fChain->SetBranchAddress("total_efficiency", &total_efficiency, &b_total_efficiency);
   fChain->SetBranchAddress("residual", &residual, &b_residual);
   fChain->SetBranchAddress("offset", &offset, &b_offset);
   fChain->SetBranchAddress("Threshold", &Threshold, &b_Threshold);
   fChain->SetBranchAddress("RunNumber", &RunNumber, &b_RunNumber);
   fChain->SetBranchAddress("HV", &HV, &b_HV);
   Notify();
}

Bool_t Scurve::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void Scurve::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t Scurve::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef Scurve_cxx
