#include "RazorHelper.h"
#include "TreeMuonSystemCombination.h"
#include "assert.h"
#include "TTree.h"

namespace {
constexpr int DEFAULT_INT = -999;
constexpr float DEFAULT_FLOAT = -999.0f;
constexpr int COUNTER_DEFAULT = 0;
constexpr int INDEX_DEFAULT = -1;
}

// Constructor
TreeMuonSystemCombination::TreeMuonSystemCombination() {
  InitVariables();
};
TreeMuonSystemCombination::~TreeMuonSystemCombination() {
  if (f_)
    f_->Close();
};
void TreeMuonSystemCombination::InitVariables() {
  #define NTUPLE_BRANCH_SCALAR(TYPE, NAME, LEAF, DEFAULT) \
    NAME = static_cast<TYPE>(DEFAULT);
  #define NTUPLE_BRANCH_ARRAY(TYPE, NAME, CAPACITY, SIZE_EXPR, LEAF, DEFAULT) \
    for (int i = 0; i < (CAPACITY); ++i) NAME[i] = static_cast<TYPE>(DEFAULT);
  #include "ntuple_branches.def"
  #undef NTUPLE_BRANCH_ARRAY
  #undef NTUPLE_BRANCH_SCALAR
};

void TreeMuonSystemCombination::InitTree() {
  assert(tree_);
  InitVariables();

  #define NTUPLE_BRANCH_SCALAR(TYPE, NAME, LEAF, DEFAULT) tree_->SetBranchAddress(#NAME, &NAME);
  #define NTUPLE_BRANCH_ARRAY(TYPE, NAME, CAPACITY, SIZE_EXPR, LEAF, DEFAULT) tree_->SetBranchAddress(#NAME, NAME);
  #include "ntuple_branches.def"
  #undef NTUPLE_BRANCH_ARRAY
  #undef NTUPLE_BRANCH_SCALAR

};

void TreeMuonSystemCombination::LoadTree(const char* file) {
  f_ = TFile::Open(file);
  assert(f_);
  tree_ = dynamic_cast<TTree*>(f_->Get("MuonSystem"));
  InitTree();
  assert(tree_);
};

void TreeMuonSystemCombination::CreateTree() {
  tree_ = new TTree("MuonSystem", "MuonSystem");
  f_ = 0;

  #define NTUPLE_BRANCH_SCALAR(TYPE, NAME, LEAF, DEFAULT) tree_->Branch(#NAME, &NAME, #NAME "/" #LEAF);
  #define NTUPLE_BRANCH_ARRAY(TYPE, NAME, CAPACITY, SIZE_EXPR, LEAF, DEFAULT) tree_->Branch(#NAME, NAME, #NAME "[" #SIZE_EXPR "]/" #LEAF);
  #include "ntuple_branches.def"
  #undef NTUPLE_BRANCH_ARRAY
  #undef NTUPLE_BRANCH_SCALAR

};
