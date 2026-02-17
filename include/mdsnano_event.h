#ifndef mdsnano_event_h
#define mdsnano_event_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

class mdsnano_event {
 public:
  TTree* fChain; //! pointer to the analyzed TTree or TChain
  Int_t fCurrent; //! current Tree number in a TChain

  // Branch payload used by llp_MuonSystem_CA_mdsnano
  #define MDS_BRANCH_SCALAR(TYPE, NAME) TYPE NAME;
  #define MDS_BRANCH_SCALAR_INIT(TYPE, NAME, INIT) TYPE NAME = INIT;
  #define MDS_BRANCH_ARRAY(TYPE, NAME, SIZE) TYPE NAME[SIZE];
  #include "mdsnano_branches.def"
  #undef MDS_BRANCH_ARRAY
  #undef MDS_BRANCH_SCALAR_INIT
  #undef MDS_BRANCH_SCALAR

  // List of branches
  #define MDS_BRANCH_SCALAR(TYPE, NAME) TBranch* b_##NAME;
  #define MDS_BRANCH_SCALAR_INIT(TYPE, NAME, INIT) TBranch* b_##NAME;
  #define MDS_BRANCH_ARRAY(TYPE, NAME, SIZE) TBranch* b_##NAME;
  #include "mdsnano_branches.def"
  #undef MDS_BRANCH_ARRAY
  #undef MDS_BRANCH_SCALAR_INIT
  #undef MDS_BRANCH_SCALAR

  mdsnano_event(TTree* tree = 0);
  virtual ~mdsnano_event();
  virtual Int_t Cut(Long64_t entry);
  virtual Int_t GetEntry(Long64_t entry);
  virtual Long64_t LoadTree(Long64_t entry);
  virtual void Init(TTree* tree);
  virtual void Loop();
  virtual Bool_t Notify();
  virtual void Show(Long64_t entry = -1);
};

#endif
