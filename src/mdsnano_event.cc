#include "mdsnano_event.h"

#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

mdsnano_event::mdsnano_event(TTree* tree) : fChain(0) {
  if (tree == 0) return;
  Init(tree);
}

mdsnano_event::~mdsnano_event() {
  if (!fChain) return;
}

Int_t mdsnano_event::GetEntry(Long64_t entry) {
  if (!fChain) return 0;
  return fChain->GetEntry(entry);
}

Long64_t mdsnano_event::LoadTree(Long64_t entry) {
  if (!fChain) return -5;
  Long64_t centry = fChain->LoadTree(entry);
  if (centry < 0) return centry;
  if (fChain->GetTreeNumber() != fCurrent) {
    fCurrent = fChain->GetTreeNumber();
    Notify();
  }
  return centry;
}

void mdsnano_event::Init(TTree* tree) {
  if (!tree) return;
  fChain = tree;
  fCurrent = -1;
  fChain->SetMakeClass(1);

  #define MDS_BRANCH_SCALAR(TYPE, NAME) fChain->SetBranchAddress(#NAME, &NAME, &b_##NAME);
  #define MDS_BRANCH_SCALAR_INIT(TYPE, NAME, INIT) fChain->SetBranchAddress(#NAME, &NAME, &b_##NAME);
  #define MDS_BRANCH_ARRAY(TYPE, NAME, SIZE) fChain->SetBranchAddress(#NAME, NAME, &b_##NAME);
  #include "mdsnano_branches.def"
  #undef MDS_BRANCH_ARRAY
  #undef MDS_BRANCH_SCALAR_INIT
  #undef MDS_BRANCH_SCALAR
  Notify();
}

Bool_t mdsnano_event::Notify() {
  return kTRUE;
}

void mdsnano_event::Show(Long64_t entry) {
  if (!fChain) return;
  fChain->Show(entry);
}

Int_t mdsnano_event::Cut(Long64_t entry) {
  (void)entry;
  return 1;
}

void mdsnano_event::Loop() {
  if (fChain == 0) return;
  Long64_t nentries = fChain->GetEntriesFast();
  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry = 0; jentry < nentries; jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry); nbytes += nb;
    if (Cut(ientry) < 0) continue;
  }
  (void)nbytes;
}
