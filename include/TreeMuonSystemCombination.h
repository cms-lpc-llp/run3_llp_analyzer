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

class TreeMuonSystemCombination
{

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
  TTree *tree_;
  TFile *f_;

  UInt_t runNum, lumiSec, evtNum, MC_condition;
  UInt_t npv, npu;
  float rho, weight;
  float pileupWeight;

  float met, metPhi;

  // rechit info
  // Int_t nCscRechits;
  // Int_t nDtRechits;
  Int_t cscRechitsClusterId[20000];
  Int_t dtRechitsClusterId[20000];
  Float_t cscRechitsX[20000];
  Float_t cscRechitsY[20000];
  Float_t cscRechitsZ[20000];
  Float_t dtRechitsX[20000];
  Float_t dtRechitsY[20000];
  Float_t dtRechitsZ[20000];
  Float_t cscRechitsTime[20000];
  Float_t dtRechitsTime[20000];
  Float_t dtRechitsTimeW[20000];
  Float_t cscRechitsTimeW[20000];

  bool Flag_goodVertices;
  bool Flag_globalSuperTightHalo2016Filter;
  bool Flag_EcalDeadCellTriggerPrimitiveFilter;
  bool Flag_BadPFMuonFilter;
  bool Flag_BadPFMuonDzFilter;
  bool Flag_hfNoisyHitsFilter;
  bool Flag_eeBadScFilter;
  bool Flag_ecalBadCalibFilter;

  bool Flag_BadChargedCandidateFilter;
  bool Flag_HBHENoiseFilter;
  bool Flag_HBHEIsoNoiseFilter;
  bool Flag_CSCTightHaloFilter;
  bool jetVeto;

  bool Flag_all;

  int mH, mX, ctau;

  bool Flag2_HBHENoiseFilter, Flag2_HBHEIsoNoiseFilter, Flag2_BadPFMuonFilter, Flag2_globalSuperTightHalo2016Filter,
      Flag2_globalTightHalo2016Filter, Flag2_BadChargedCandidateFilter, Flag2_EcalDeadCellTriggerPrimitiveFilter,
      Flag2_ecalBadCalibFilter, Flag2_eeBadScFilter, Flag2_all;

  // csc
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
  int dtRechitClusterNHitWheel0[N_MAX_CSC];
  int dtRechitClusterNHitWheel1[N_MAX_CSC];
  int dtRechitClusterNHitWheel2[N_MAX_CSC];
  int nDTRechitsSector[4][5][12]; //[NStattion][NWheels][NSector]

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
  int dtRechitClusterMaxStation[N_MAX_CSC];        //[nCsc]
  int dtRechitClusterNStation10[N_MAX_CSC];
  float dtRechitClusterAvgStation10[N_MAX_CSC];
  float dtRechitClusterMaxChamberRatio[N_MAX_CSC]; //[nCsc]
  int dtRechitClusterMaxChamber[N_MAX_CSC];        //[nCsc]
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
  int cscRechitClusterMaxStation[N_MAX_CSC];        //[nCsc]
  int cscRechitClusterNStation10[N_MAX_CSC];
  float cscRechitClusterAvgStation10[N_MAX_CSC];
  float cscRechitClusterMaxChamberRatio[N_MAX_CSC]; //[nCsc]
  int cscRechitClusterMaxChamber[N_MAX_CSC];        //[nCsc]
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

  // gLLP
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

  // leptons

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

  // jets
  int nJets;
  float jetE[N_MAX_JETS];
  float jetPt[N_MAX_JETS];
  float jetEta[N_MAX_JETS];
  float jetPhi[N_MAX_JETS];

  bool jetTightPassId[N_MAX_JETS];
  bool HLTDecision[NTriggersMAX];

  Int_t nGenJetAK8 = 0;
  Float_t GenJetAK8_eta[6 * 10] = {0};  //[nGenJetAK8]
  Float_t GenJetAK8_mass[6 * 10] = {0}; //[nGenJetAK8]
  Float_t GenJetAK8_phi[6 * 10] = {0};  //[nGenJetAK8]
  Float_t GenJetAK8_pt[6 * 10] = {0};   //[nGenJetAK8]
  Int_t nGenJet = 0;
  Float_t GenJet_eta[19 * 10] = {0};  //[nGenJet]
  Float_t GenJet_mass[19 * 10] = {0}; //[nGenJet]
  Float_t GenJet_phi[19 * 10] = {0};  //[nGenJet]
  Float_t GenJet_pt[19 * 10] = {0};   //[nGenJet]
  Int_t nGenPart = 0;
  Short_t GenPart_genPartIdxMother[125 * 10] = {0}; //[nGenPart]
  UShort_t GenPart_statusFlags[125 * 10] = {0};     //[nGenPart]
  Int_t GenPart_pdgId[125 * 10] = {0};              //[nGenPart]
  Int_t GenPart_status[125 * 10] = {0};             //[nGenPart]
  Float_t GenPart_eta[125 * 10] = {0};              //[nGenPart]
  Float_t GenPart_mass[125 * 10] = {0};             //[nGenPart]
  Float_t GenPart_phi[125 * 10] = {0};              //[nGenPart]
  Float_t GenPart_pt[125 * 10] = {0};               //[nGenPart]
  Int_t nGenProton = 0;
  Bool_t GenProton_isPU[11 * 10] = {0}; //[nGenProton]
  Float_t GenProton_px[11 * 10] = {0};  //[nGenProton]
  Float_t GenProton_py[11 * 10] = {0};  //[nGenProton]
  Float_t GenProton_pz[11 * 10] = {0};  //[nGenProton]
  Float_t GenProton_vz[11 * 10] = {0};  //[nGenProton]
  Int_t nSubGenJetAK8 = 0;
  Float_t SubGenJetAK8_eta[12 * 10] = {0};  //[nSubGenJetAK8]
  Float_t SubGenJetAK8_mass[12 * 10] = {0}; //[nSubGenJetAK8]
  Float_t SubGenJetAK8_phi[12 * 10] = {0};  //[nSubGenJetAK8]
  Float_t SubGenJetAK8_pt[12 * 10] = {0};   //[nSubGenJetAK8]
  Int_t Generator_id1 = 0;
  Int_t Generator_id2 = 0;
  Float_t Generator_binvar = 0;
  Float_t Generator_scalePDF = 0;
  Float_t Generator_weight = 0;
  Float_t Generator_x1 = 0;
  Float_t Generator_x2 = 0;
  Float_t Generator_xpdf1 = 0;
  Float_t Generator_xpdf2 = 0;
  Float_t GenVtx_x = 0;
  Float_t GenVtx_y = 0;
  Float_t GenVtx_z = 0;
  Int_t nGenVisTau = 0;
  UChar_t GenVisTau_status[2 * 10] = {0};           //[nGenVisTau]
  Short_t GenVisTau_charge[2 * 10] = {0};           //[nGenVisTau]
  Short_t GenVisTau_genPartIdxMother[2 * 10] = {0}; //[nGenVisTau]
  Float_t GenVisTau_eta[2 * 10] = {0};              //[nGenVisTau]
  Float_t GenVisTau_mass[2 * 10] = {0};             //[nGenVisTau]
  Float_t GenVisTau_phi[2 * 10] = {0};              //[nGenVisTau]
  Float_t GenVisTau_pt[2 * 10] = {0};               //[nGenVisTau]
  Float_t genWeight = 0;
  Float_t GenMET_phi = 0;
  Float_t GenMET_pt = 0;
  Int_t nGenDressedLepton = 0;
  Bool_t GenDressedLepton_hasTauAnc[1 * 10] = {0}; //[nGenDressedLepton]
  Int_t GenDressedLepton_pdgId[1 * 10] = {0};      //[nGenDressedLepton]
  Float_t GenDressedLepton_eta[1 * 10] = {0};      //[nGenDressedLepton]
  Float_t GenDressedLepton_mass[1 * 10] = {0};     //[nGenDressedLepton]
  Float_t GenDressedLepton_phi[1 * 10] = {0};      //[nGenDressedLepton]
  Float_t GenDressedLepton_pt[1 * 10] = {0};       //[nGenDressedLepton]
  Int_t nGenIsolatedPhoton = 0;
  Float_t GenIsolatedPhoton_eta[2 * 10] = {0};  //[nGenIsolatedPhoton]
  Float_t GenIsolatedPhoton_mass[2 * 10] = {0}; //[nGenIsolatedPhoton]
  Float_t GenIsolatedPhoton_phi[2 * 10] = {0};  //[nGenIsolatedPhoton]
  Float_t GenIsolatedPhoton_pt[2 * 10] = {0};   //[nGenIsolatedPhoton]
  Int_t genTtbarId = 0;
  Int_t nboostedTau = 0;
  UChar_t boostedTau_genPartFlav[2 * 10] = {0};  //[nboostedTau]
  Short_t boostedTau_genPartIdx[2 * 10] = {0};   //[nboostedTau]
  Int_t nsoftActivityVH = 0;
  Int_t nElectron = 0;
  UChar_t Electron_genPartFlav[5 * 10] = {0};    //[nElectron]
  Short_t Electron_genPartIdx[5 * 10] = {0};     //[nElectron]
  Int_t nFatJet = 0;
  Short_t FatJet_genJetAK8Idx[4 * 10] = {0};     //[nFatJet]
  UChar_t GenJetAK8_hadronFlavour[6 * 10] = {0}; //[nGenJetAK8]
  Short_t GenJetAK8_partonFlavour[6 * 10] = {0}; //[nGenJetAK8]
  UChar_t GenJet_hadronFlavour[19 * 10] = {0};   //[nGenJet]
  Short_t GenJet_partonFlavour[19 * 10] = {0};   //[nGenJet]
  Float_t GenVtx_t0 = 0;
  Short_t Jet_genJetIdx[16 * 10] = {0};            //[nJet]
  Int_t nLowPtElectron = 0;
  UChar_t LowPtElectron_genPartFlav[6 * 10] = {0}; //[nLowPtElectron]
  Short_t LowPtElectron_genPartIdx[6 * 10] = {0};  //[nLowPtElectron]
  Int_t nMuon = 0;
  UChar_t Muon_genPartFlav[8 * 10] = {0};          //[nMuon]
  Short_t Muon_genPartIdx[8 * 10] = {0};           //[nMuon]
  Int_t nPhoton = 0;
  UChar_t Photon_genPartFlav[6 * 10] = {0};        //[nPhoton]
  Short_t Photon_genPartIdx[6 * 10] = {0};         //[nPhoton]
  Float_t MET_fiducialGenPhi = 0;
  Float_t MET_fiducialGenPt = 0;
  Int_t nTau = 0;
  UChar_t Tau_genPartFlav[5 * 10] = {0}; //[nTau]
  Short_t Tau_genPartIdx[5 * 10] = {0};  //[nTau]

  void InitVariables();
  void InitTree();
  void LoadTree(const char *file);
  void CreateTree();
};
#endif
