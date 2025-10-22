// Class to manage files for b-tag scale factors, lepton scale factors, pileup weights, and other information

#ifndef TreeMuonSystemCombination_TnP_H
#define TreeMuonSystemCombination_TnP_H

#define N_MAX_LEPTONS 100
#define N_MAX_JETS 100
#define N_MAX_CSC 200
#define N_MAX_CSCRECHITS 5000
#define N_MAX_DTRECHITS 20000
#define NTriggersMAX 1201 // Number of trigger in the .dat file
#define N_CSC_CUT 20
#define JET_PT_CUT 10
#define MUON_PT_CUT 20
#define N_MAX_GPARTICLES 5000
#define N_MAX_LLP 200

#include <iostream>
#include <string>
#include <sys/stat.h>
#include "assert.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2D.h"
#include "TRandom.h"
#include "TTree.h"
#include "DBSCAN.h"

#include "RazorAnalyzer.h"

#include "RazorHelper.h"

class TreeMuonSystemCombination_TnP {
 public:
  TreeMuonSystemCombination_TnP();
  ~TreeMuonSystemCombination_TnP();
  // TreeMuonSystemCombination_TnP::TreeMuonSystemCombination_TnP()
  // {
  //   InitVariables();
  // };
  // TreeMuonSystemCombination_TnP::~TreeMuonSystemCombination_TnP()
  // {
  //   if (f_) f_->Close();
  // };
  TTree* tree_;
  TFile* f_;

  UInt_t runNum, lumiSec, evtNum, MC_condition;
  UInt_t npv, npu;
  float rho, weight;
  float pileupWeight;

  float met, metPhi;
  bool Flag_HBHENoiseFilter, Flag_HBHEIsoNoiseFilter, Flag_BadPFMuonFilter, Flag_globalSuperTightHalo2016Filter,
      Flag_CSCTightHaloFilter, Flag_BadChargedCandidateFilter, Flag_eeBadScFilter, Flag_goodVertices, Flag_ecalBadCalibFilter, Flag_all;
  int mH, mX, ctau;

  bool Flag2_HBHENoiseFilter, Flag2_HBHEIsoNoiseFilter, Flag2_BadPFMuonFilter, Flag2_globalSuperTightHalo2016Filter,
      Flag2_globalTightHalo2016Filter, Flag2_BadChargedCandidateFilter, Flag2_EcalDeadCellTriggerPrimitiveFilter,
      Flag2_ecalBadCalibFilter, Flag2_eeBadScFilter, Flag2_all;

  //csc
  int nCscRechits;
  int nCscRings;
  int nDTRechits;
  int nDtRings;

  int nDtRechitClusters;
  bool dtRechitClusterOverlap[N_MAX_CSC];
  int dtRechitClusterNSegStation1[N_MAX_CSC];
  int dtRechitClusterNSegStation2[N_MAX_CSC];
  int dtRechitClusterNSegStation3[N_MAX_CSC];
  int dtRechitClusterNSegStation4[N_MAX_CSC];
  int dtRechitClusterNOppositeSegStation1[N_MAX_CSC];
  int dtRechitClusterNOppositeSegStation2[N_MAX_CSC];
  int dtRechitClusterNOppositeSegStation3[N_MAX_CSC];
  int dtRechitClusterNOppositeSegStation4[N_MAX_CSC];
  int dtRechitClusterNHitStation1[N_MAX_CSC];
  int dtRechitClusterNHitStation2[N_MAX_CSC];
  int dtRechitClusterNHitStation3[N_MAX_CSC];
  int dtRechitClusterNHitStation4[N_MAX_CSC];

  int dtRechitCluster_match_MB1hits_0p4[N_MAX_CSC];
  int dtRechitCluster_match_MB1hits_0p5[N_MAX_CSC];
  int dtRechitCluster_match_MB1hits_cosmics_plus[N_MAX_CSC];
  int dtRechitCluster_match_MB1hits_cosmics_minus[N_MAX_CSC];
  int dtRechitCluster_match_MB1Seg_0p4[N_MAX_CSC];
  int dtRechitCluster_match_MB1Seg_0p5[N_MAX_CSC];
  int dtRechitCluster_match_RPCBx_dPhi0p5[N_MAX_CSC];
  int dtRechitCluster_match_RPChits_dPhi0p5[N_MAX_CSC];

  int dtRechitCluster_match_RB1_0p4[N_MAX_CSC];
  int dtRechitCluster_match_RB1_dPhi0p5[N_MAX_CSC];

  bool dtRechitCluster_match_gLLP[N_MAX_CSC];
  int dtRechitCluster_match_gLLP_index[N_MAX_CSC];
  float dtRechitCluster_match_gLLP_minDeltaR[N_MAX_CSC];
  float dtRechitCluster_match_gLLP_eta[N_MAX_CSC];
  float dtRechitCluster_match_gLLP_phi[N_MAX_CSC];
  float dtRechitCluster_match_gLLP_decay_r[N_MAX_CSC];
  float dtRechitCluster_match_gLLP_decay_z[N_MAX_CSC];
  bool dtRechitCluster_match_gLLP_csc[N_MAX_CSC];
  bool dtRechitCluster_match_gLLP_dt[N_MAX_CSC];
  float dtRechitCluster_match_gLLP_e[N_MAX_CSC];

  float dtRechitClusterX[N_MAX_CSC]; //[nCsc]
  float dtRechitClusterY[N_MAX_CSC]; //[nCsc]
  float dtRechitClusterZ[N_MAX_CSC]; //[nCsc]

  int dtRechitClusterWheel[N_MAX_CSC];

  float dtRechitClusterEta[N_MAX_CSC]; //[nCsc]
  float dtRechitClusterPhi[N_MAX_CSC]; //[nCsc]
  int dtRechitClusterSize[N_MAX_CSC];
  int dtRechitClusterNoiseHit[N_MAX_CSC];
  int dtRechitClusterNoiseHitStation1[N_MAX_CSC];
  int dtRechitClusterNoiseHitStation2[N_MAX_CSC];
  int dtRechitClusterNoiseHitStation3[N_MAX_CSC];
  int dtRechitClusterNoiseHitStation4[N_MAX_CSC];

  float dtRechitClusterMaxStationRatio[N_MAX_CSC]; //[nCsc]
  int dtRechitClusterMaxStation[N_MAX_CSC]; //[nCsc]
  int dtRechitClusterNStation10[N_MAX_CSC];
  float dtRechitClusterAvgStation10[N_MAX_CSC];
  float dtRechitClusterMaxChamberRatio[N_MAX_CSC]; //[nCsc]
  int dtRechitClusterMaxChamber[N_MAX_CSC]; //[nCsc]
  int dtRechitClusterNChamber[N_MAX_CSC];

  float dtRechitClusterJetVetoPt[N_MAX_CSC];
  float dtRechitClusterJetVetoE[N_MAX_CSC];
  bool dtRechitClusterJetVetoLooseId[N_MAX_CSC];
  bool dtRechitClusterJetVetoTightId[N_MAX_CSC];

  float dtRechitClusterMuonVetoPt[N_MAX_CSC];
  float dtRechitClusterMuonVetoE[N_MAX_CSC];

  bool dtRechitClusterMuonVetoTightId[N_MAX_CSC];
  bool dtRechitClusterMuonVetoLooseId[N_MAX_CSC];
  bool dtRechitClusterMuonVetoGlobal[N_MAX_CSC];

  float dtRechitClusterMet_dPhi[N_MAX_CSC];

  int nCscRechitClusters;

  bool cscRechitCluster_match_gLLP[N_MAX_CSC];
  int cscRechitCluster_match_gLLP_index[N_MAX_CSC];
  float cscRechitCluster_match_gLLP_minDeltaR[N_MAX_CSC];
  float cscRechitCluster_match_gLLP_eta[N_MAX_CSC];
  float cscRechitCluster_match_gLLP_phi[N_MAX_CSC];
  float cscRechitCluster_match_gLLP_decay_r[N_MAX_CSC];
  float cscRechitCluster_match_gLLP_decay_z[N_MAX_CSC];
  bool cscRechitCluster_match_gLLP_csc[N_MAX_CSC];
  bool cscRechitCluster_match_gLLP_dt[N_MAX_CSC];
  float cscRechitCluster_match_gLLP_e[N_MAX_CSC];

  float cscRechitClusterX[N_MAX_CSC]; //[nCsc]
  float cscRechitClusterY[N_MAX_CSC]; //[nCsc]
  float cscRechitClusterZ[N_MAX_CSC]; //[nCsc]
  float cscRechitClusterTimeWeighted[N_MAX_CSC];
  float cscRechitClusterTimeSpreadWeightedAll[N_MAX_CSC];
  float cscRechitClusterTime[N_MAX_CSC];
  float cscRechitClusterTimeSpread[N_MAX_CSC];
  float cscRechitClusterEta[N_MAX_CSC]; //[nCsc]
  float cscRechitClusterPhi[N_MAX_CSC]; //[nCsc]
  int cscRechitClusterSize[N_MAX_CSC];

  float cscRechitClusterMaxStationRatio[N_MAX_CSC]; //[nCsc]
  int cscRechitClusterMaxStation[N_MAX_CSC]; //[nCsc]
  int cscRechitClusterNStation10[N_MAX_CSC];
  float cscRechitClusterAvgStation10[N_MAX_CSC];
  float cscRechitClusterMaxChamberRatio[N_MAX_CSC]; //[nCsc]
  int cscRechitClusterMaxChamber[N_MAX_CSC]; //[nCsc]
  int cscRechitClusterNChamber[N_MAX_CSC];

  float cscRechitClusterJetVetoPt[N_MAX_CSC];
  bool cscRechitClusterJetVetoLooseId[N_MAX_CSC];
  bool cscRechitClusterJetVetoTightId[N_MAX_CSC];
  float cscRechitClusterJetVetoE[N_MAX_CSC];

  float cscRechitClusterMuonVetoPt[N_MAX_CSC];
  float cscRechitClusterMuonVetoE[N_MAX_CSC];
  bool cscRechitClusterMuonVetoLooseId[N_MAX_CSC];
  bool cscRechitClusterMuonVetoGlobal[N_MAX_CSC];

  int cscRechitCluster_match_dtSeg_0p4[N_MAX_CSC];
  int cscRechitCluster_match_MB1Seg_0p4[N_MAX_CSC];
  int cscRechitCluster_match_RE12_0p4[N_MAX_CSC];
  int cscRechitCluster_match_RB1_0p4[N_MAX_CSC];

  int cscRechitClusterNRechitChamberPlus11[N_MAX_CSC];
  int cscRechitClusterNRechitChamberPlus12[N_MAX_CSC];
  int cscRechitClusterNRechitChamberPlus13[N_MAX_CSC];
  int cscRechitClusterNRechitChamberPlus21[N_MAX_CSC];
  int cscRechitClusterNRechitChamberPlus22[N_MAX_CSC];
  int cscRechitClusterNRechitChamberPlus31[N_MAX_CSC];
  int cscRechitClusterNRechitChamberPlus32[N_MAX_CSC];
  int cscRechitClusterNRechitChamberPlus41[N_MAX_CSC];
  int cscRechitClusterNRechitChamberPlus42[N_MAX_CSC];

  int cscRechitClusterNRechitChamberMinus11[N_MAX_CSC];
  int cscRechitClusterNRechitChamberMinus12[N_MAX_CSC];
  int cscRechitClusterNRechitChamberMinus13[N_MAX_CSC];
  int cscRechitClusterNRechitChamberMinus21[N_MAX_CSC];
  int cscRechitClusterNRechitChamberMinus22[N_MAX_CSC];
  int cscRechitClusterNRechitChamberMinus31[N_MAX_CSC];
  int cscRechitClusterNRechitChamberMinus32[N_MAX_CSC];
  int cscRechitClusterNRechitChamberMinus41[N_MAX_CSC];
  int cscRechitClusterNRechitChamberMinus42[N_MAX_CSC];
  float cscRechitClusterMet_dPhi[N_MAX_CSC];

  //gLLP
  int nGLLP;
  float gLLP_eta[N_MAX_LLP];
  float gLLP_phi[N_MAX_LLP];
  float gLLP_csc[N_MAX_LLP];
  float gLLP_dt[N_MAX_LLP];
  float gLLP_beta[N_MAX_LLP];
  float gLLP_e[N_MAX_LLP];
  float gLLP_pt[N_MAX_LLP];
  float gLLP_ctau[N_MAX_LLP];
  float gLLP_decay_vertex_r[N_MAX_LLP];
  float gLLP_decay_vertex_x[N_MAX_LLP];
  float gLLP_decay_vertex_y[N_MAX_LLP];
  float gLLP_decay_vertex_z[N_MAX_LLP];

  float gHiggsPt;
  float gHiggsEta;
  float gHiggsPhi;
  float gHiggsE;

  //leptons

  int nLeptons;
  float lepE[N_MAX_LEPTONS];
  float lepPt[N_MAX_LEPTONS];
  float lepEta[N_MAX_LEPTONS];
  float lepPhi[N_MAX_LEPTONS];
  int lepPdgId[N_MAX_LEPTONS];
  float lepDZ[N_MAX_LEPTONS];

  bool lepTightId[N_MAX_LEPTONS];
  bool lepPassLooseIso[N_MAX_LEPTONS];
  bool lepPassTightIso[N_MAX_LEPTONS];
  bool lepPassVTightIso[N_MAX_LEPTONS];
  bool lepPassVVTightIso[N_MAX_LEPTONS];

  //jets
  int nJets;
  float jetE[N_MAX_JETS];
  float jetPt[N_MAX_JETS];
  float jetEta[N_MAX_JETS];
  float jetPhi[N_MAX_JETS];

  bool jetTightPassId[N_MAX_JETS];
  bool HLTDecision[NTriggersMAX];

  void InitVariables();
  void InitTree();
  void LoadTree(const char* file);
  void CreateTree();
};
#endif
