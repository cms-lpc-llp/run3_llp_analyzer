// Class to manage files for b-tag scale factors, lepton scale factors, pileup weights, and other information

#ifndef TreeMuonSystemCombination_H
#define TreeMuonSystemCombination_H

#define N_MAX_LEPTONS 100
#define N_MAX_JETS 100
#define N_MAX_CSC 200
#define N_MAX_CSCRECHITS 5000
#define N_MAX_DTRECHITS 20000
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
#include "DBSCAN.h"

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

  UInt_t runNum, lumiSec, evtNum, MC_condition;
  UInt_t npv, npu;
  float rho, weight;
  float pileupWeight, pileupWeightUp, pileupWeightDown;
  bool HLT_CSCCSC, HLT_CSCDT, jetVeto;
  float met, metPhi, Puppimet, PuppimetPhi;
  float PuppimetJESUp, PuppimetJESDown, PuppimetPhiJESDown, PuppimetPhiJESUp;
  float metJESUp, metJESDown, metPhiJESDown, metPhiJESUp;
  bool Flag_goodVertices, Flag_EcalDeadCellTriggerPrimitiveFilter, Flag_BadPFMuonFilter, Flag_BadPFMuonDzFilter, Flag_globalSuperTightHalo2016Filter,
      Flag_hfNoisyHitsFilter, Flag_eeBadScFilter, Flag_ecalBadCalibFilter, Flag_all;
  int mH, mX, ctau;
  float LHEScaleWeight[9];

  //csc
  int nCscRechits;
  int nCscRings;
  int nDtRechits;
  int nDtRings;
  bool HLT_CscCluster100_PNetTauhPFJet10_Loose;
  bool HLT_CscCluster100_Ele5;
  bool HLT_CscCluster100_Mu5;
  bool HLT_CscCluster50_Photon30Unseeded;
  bool HLT_CscCluster50_Photon20Unseeded;
  bool HLT_PFMET120_PFMHT120_IDTight;

  float MetTriggerEff;
  float MetTriggerEffUp;
  float MetTriggerEffDown;

  int nDtRechitClusters;
  int nDtRechitClusters_nocut;

  float DtRechitsEta[N_MAX_DTRECHITS];
  float DtRechitsPhi[N_MAX_DTRECHITS];
  float DtRechitsX[N_MAX_DTRECHITS];
  float DtRechitsY[N_MAX_DTRECHITS];
  float DtRechitsZ[N_MAX_DTRECHITS];
  int DtRechitsClusterId[N_MAX_DTRECHITS];

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
  int dtRechitClusterNHitWheel0[N_MAX_CSC];
  int dtRechitClusterNHitWheel1[N_MAX_CSC];
  int dtRechitClusterNHitWheel2[N_MAX_CSC];

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
  float dtRechitClusterJetVetoPtJESDown[N_MAX_CSC];
  float dtRechitClusterJetVetoPtJESUp[N_MAX_CSC];
  float dtRechitClusterJetVetoE[N_MAX_CSC];
  bool dtRechitClusterJetVetoLooseId[N_MAX_CSC];
  bool dtRechitClusterJetVetoTightId[N_MAX_CSC];

  float dtRechitClusterMuonVetoPt[N_MAX_CSC];
  float dtRechitClusterMuonVetoE[N_MAX_CSC];

  bool dtRechitClusterMuonVetoTightId[N_MAX_CSC];
  bool dtRechitClusterMuonVetoLooseId[N_MAX_CSC];
  bool dtRechitClusterMuonVetoGlobal[N_MAX_CSC];

  float dtRechitClusterMet_dPhi[N_MAX_CSC];
  float dtRechitClusterMetJESUp_dPhi[N_MAX_CSC];
  float dtRechitClusterMetJESDown_dPhi[N_MAX_CSC];
  float dtRechitClusterPuppiMet_dPhi[N_MAX_CSC];

  int dtRechitClusternXY[N_MAX_CSC];
  int dtRechitClusternZ[N_MAX_CSC];
  float dtRechitClusterXSpread[N_MAX_CSC];
  float dtRechitClusterYSpread[N_MAX_CSC];
  float dtRechitClusterZSpread[N_MAX_CSC];
  float dtRechitClusterXYSpread[N_MAX_CSC];
  float dtRechitClusterRSpread[N_MAX_CSC];
  float dtRechitClusterEtaPhiSpread[N_MAX_CSC];
  float dtRechitClusterEtaSpread[N_MAX_CSC];
  float dtRechitClusterPhiSpread[N_MAX_CSC];
  float dtRechitClusterDeltaRSpread[N_MAX_CSC];
  float dtRechitClusterMajorAxis[N_MAX_CSC];
  float dtRechitClusterMinorAxis[N_MAX_CSC];
  float dtRechitClusterSkewX[N_MAX_CSC];
  float dtRechitClusterSkewY[N_MAX_CSC];
  float dtRechitClusterSkewZ[N_MAX_CSC];
  float dtRechitClusterKurtX[N_MAX_CSC];
  float dtRechitClusterKurtY[N_MAX_CSC];
  float dtRechitClusterKurtZ[N_MAX_CSC];

  int nCscRechitClusters;
  int nCscRechitClusters_nocut;

  float CscRechitsEta[N_MAX_CSCRECHITS];
  float CscRechitsPhi[N_MAX_CSCRECHITS];
  float CscRechitsX[N_MAX_CSCRECHITS];
  float CscRechitsY[N_MAX_CSCRECHITS];
  float CscRechitsZ[N_MAX_CSCRECHITS];
  int CscRechitsClusterId[N_MAX_CSCRECHITS];

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
  int cscRechitClusternXY[N_MAX_CSC];
  int cscRechitClusternZ[N_MAX_CSC];
  float cscRechitClusterXSpread[N_MAX_CSC];
  float cscRechitClusterYSpread[N_MAX_CSC];
  float cscRechitClusterZSpread[N_MAX_CSC];
  float cscRechitClusterXYSpread[N_MAX_CSC];
  float cscRechitClusterRSpread[N_MAX_CSC];
  float cscRechitClusterEtaPhiSpread[N_MAX_CSC];
  float cscRechitClusterEtaSpread[N_MAX_CSC];
  float cscRechitClusterPhiSpread[N_MAX_CSC];
  float cscRechitClusterDeltaRSpread[N_MAX_CSC];
  float cscRechitClusterMajorAxis[N_MAX_CSC];
  float cscRechitClusterMinorAxis[N_MAX_CSC];
  float cscRechitClusterSkewX[N_MAX_CSC];
  float cscRechitClusterSkewY[N_MAX_CSC];
  float cscRechitClusterSkewZ[N_MAX_CSC];
  float cscRechitClusterKurtX[N_MAX_CSC];
  float cscRechitClusterKurtY[N_MAX_CSC];
  float cscRechitClusterKurtZ[N_MAX_CSC];

  float cscRechitClusterMaxStationRatio[N_MAX_CSC]; //[nCsc]
  int cscRechitClusterMaxStation[N_MAX_CSC]; //[nCsc]
  int cscRechitClusterNStation10[N_MAX_CSC];
  float cscRechitClusterAvgStation10[N_MAX_CSC];
  float cscRechitClusterMaxChamberRatio[N_MAX_CSC]; //[nCsc]
  int cscRechitClusterMaxChamber[N_MAX_CSC]; //[nCsc]
  int cscRechitClusterNChamber[N_MAX_CSC];

  float cscRechitClusterJetVetoPt[N_MAX_CSC];
  float cscRechitClusterJetVetoPtJESDown[N_MAX_CSC];
  float cscRechitClusterJetVetoPtJESUp[N_MAX_CSC];
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
  float cscRechitClusterHMTEfficiency[N_MAX_CSC];

  // int           cscRechitClusterHMTEffPlus11[N_MAX_CSC];
  // int           cscRechitClusterHMTEffPlus12[N_MAX_CSC];
  // int           cscRechitClusterHMTEffPlus13[N_MAX_CSC];
  // int           cscRechitClusterHMTEffPlus21[N_MAX_CSC];
  // int           cscRechitClusterHMTEffPlus22[N_MAX_CSC];
  // int           cscRechitClusterHMTEffPlus31[N_MAX_CSC];
  // int           cscRechitClusterHMTEffPlus32[N_MAX_CSC];
  // int           cscRechitClusterHMTEffPlus41[N_MAX_CSC];
  // int           cscRechitClusterHMTEffPlus42[N_MAX_CSC];
  // int           cscRechitClusterHMTEffMinus11[N_MAX_CSC];
  // int           cscRechitClusterHMTEffMinus12[N_MAX_CSC];
  // int           cscRechitClusterHMTEffMinus13[N_MAX_CSC];
  // int           cscRechitClusterHMTEffMinus21[N_MAX_CSC];
  // int           cscRechitClusterHMTEffMinus22[N_MAX_CSC];
  // int           cscRechitClusterHMTEffMinus31[N_MAX_CSC];
  // int           cscRechitClusterHMTEffMinus32[N_MAX_CSC];
  // int           cscRechitClusterHMTEffMinus41[N_MAX_CSC];
  // int           cscRechitClusterHMTEffMinus42[N_MAX_CSC];

  float cscRechitClusterMet_dPhi[N_MAX_CSC];
  float cscRechitClusterMetJESDown_dPhi[N_MAX_CSC];
  float cscRechitClusterMetJESUp_dPhi[N_MAX_CSC];
  float cscRechitClusterPuppiMet_dPhi[N_MAX_CSC];

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

  int nGenParticles;
  float gParticleEta[N_MAX_GPARTICLES];
  float gParticlePhi[N_MAX_GPARTICLES];
  float gParticlePt[N_MAX_GPARTICLES];
  int gParticleId[N_MAX_GPARTICLES];
  float gParticle_decay_vertex_r[N_MAX_GPARTICLES];
  float gParticle_decay_vertex_x[N_MAX_GPARTICLES];
  float gParticle_decay_vertex_y[N_MAX_GPARTICLES];
  float gParticle_decay_vertex_z[N_MAX_GPARTICLES];

  float gParticle_ctau[N_MAX_GPARTICLES];
  float gParticle_beta[N_MAX_GPARTICLES];

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
  float jetPtJESUp[N_MAX_JETS];
  float jetPtJESDown[N_MAX_JETS];
  float jetEta[N_MAX_JETS];
  float jetPhi[N_MAX_JETS];

  bool jetTightPassId[N_MAX_JETS];
  // bool HLTDecision[NTriggersMAX];

  int nTaus; //other variable for ID taus
  float tauM[N_MAX_JETS];
  float tauPt[N_MAX_JETS];
  float tauEta[N_MAX_JETS];
  float tauPhi[N_MAX_JETS];
  float tauE[N_MAX_JETS];
  float tauDeltaR[N_MAX_JETS];
  int tauDecayMode[N_MAX_JETS];
  float tauDz[N_MAX_JETS];
  bool tauIsVVVLoose[N_MAX_JETS];
  bool tauIsVVLoose[N_MAX_JETS];
  bool tauIsVLoose[N_MAX_JETS];
  bool tauIsLoose[N_MAX_JETS];
  bool tauIsMedium[N_MAX_JETS];
  bool tauIsTight[N_MAX_JETS];
  bool tauIsVTight[N_MAX_JETS];
  bool tauIsVVTight[N_MAX_JETS];
  int tauGenPartFlav[N_MAX_JETS];
  float deltaR_GenTauRecoTau[N_MAX_JETS];
  int tauIdDeepTau2018v2p5VSjet[N_MAX_JETS];
  int tauIdDeepTau2018v2p5VSe[N_MAX_JETS];
  int tauIdDeepTau2018v2p5VSmu[N_MAX_JETS];
  float tauFractionOfGenVisEnergy[N_MAX_JETS];
  float tauFractionOfGenVisPt[N_MAX_JETS];

  //gen Tau branches for HNL+Tau
  int nGenTaus;
  int nGenVisTau;
  int gTauPdgId[N_MAX_GTAU];
  float gTauEta[N_MAX_GTAU];
  float gTauPhi[N_MAX_GTAU];
  float gTauPt[N_MAX_GTAU];
  float gTauE[N_MAX_GTAU];
  int gVisTauDecayMode[N_MAX_GTAU];
  bool gTauHadronicDecay[N_MAX_GTAU];
  bool gTauEDecay[N_MAX_GTAU];
  bool gTauMuDecay[N_MAX_GTAU];
  float gVisTauEta[N_MAX_GTAU];
  float gVisTauPhi[N_MAX_GTAU];
  float gVisTauPt[N_MAX_GTAU];
  float gVisTauE[N_MAX_GTAU];
  float gVisTauFractionOfTotalPt[N_MAX_GTAU];
  float gVisTauFractionOfTotalEnergy[N_MAX_GTAU];

  void InitVariables();
  void InitTree();
  void LoadTree(const char* file);
  void CreateTree();
};
#endif
