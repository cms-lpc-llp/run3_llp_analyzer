// Class to manage files for b-tag scale factors, lepton scale factors, pileup weights, and other information

#ifndef TreeMuonSystemCombination_H
#define TreeMuonSystemCombination_H

#define N_MAX_LEPTONS 100
#define N_MAX_JETS 100
#define N_MAX_CSC 200
#define N_MAX_CSCRECHITS 5000
#define N_MAX_DTRECHITS 20000
#define N_MAX_RPCRECHITS 2000
#define NTriggersMAX 1201 // Number of trigger in the .dat file
#define N_CSC_CUT 20
#define JET_PT_CUT 10
#define MUON_PT_CUT 20
#define N_MAX_GPARTICLES 500
#define N_MAX_LLP 200
#define N_MAX_GTAU 100

#include <iostream>
#include <string>
#include <sys/stat.h>
#include "assert.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2D.h"
#include "TRandom.h"
#include "TTree.h"

#include "RazorAnalyzer.h"

#include "RazorHelper.h"

class TreeMuonSystemCombination {
 public:
  TreeMuonSystemCombination();
  ~TreeMuonSystemCombination();
  // TreeMuonSystemCombination::TreeMuonSystemCombination()
  // {
  //   InitVariables();
  // };
  // TreeMuonSystemCombination::~TreeMuonSystemCombination()
  // {
  //   if (f_) f_->Close();
  // };
  TTree* tree_;
  TFile* f_;

  // Ntuple payload and output branch schema (single source of truth).
  #define NTUPLE_BRANCH_SCALAR(TYPE, NAME, LEAF, DEFAULT) TYPE NAME;
  #define NTUPLE_BRANCH_ARRAY(TYPE, NAME, CAPACITY, SIZE_EXPR, LEAF, DEFAULT) TYPE NAME[CAPACITY];
  #include "ntuple_branches.def"
  #undef NTUPLE_BRANCH_ARRAY
  #undef NTUPLE_BRANCH_SCALAR

  void InitVariables();
  void InitTree();
  void LoadTree(const char* file);
  void CreateTree();
};
#endif
