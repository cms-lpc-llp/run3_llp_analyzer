/* #region: includes and usings */
#include "llp_MuonSystem_CA_mdsnano.h"
#include "RazorHelper.h"
#include "TreeMuonSystemCombination.h"

#include "CACluster.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include "TMath.h"
#include "Math/Vector3D.h"
#include <array>
#include <algorithm>
#include <cctype>
#include <cmath>
#include <iostream>
#include <limits>
#include <memory>
#include <random>

//C++ includes
#include "assert.h"
//ROOT includes
#include "TH1F.h"

//using namespace fastjet;
using namespace std::chrono;
using namespace std;
using namespace ROOT::Math;
/* #endregion */

const float ELE_MASS = 0.000511;
const float MU_MASS = 0.105658;
const int INDEX_DEFAULT = -1;

/* #region: Usefull templates/functions/structs  */
// Defining a helper for making single binned histograms.
std::unique_ptr<TH1F> makeCounterHist(const std::string& name) {
  auto h = std::make_unique<TH1F>(name.c_str(), name.c_str(), 1, 1, 2);
  h->SetDirectory(nullptr); // important: avoid ROOT owning/deleting it
  return h;
}

// Defining leptons and jets as structures
struct leptons {
  TLorentzVector lepton;
  #define NTUPLE_RECO_LEPTON_TO_STRUCT(TYPE, BRANCH, STRUCT_FIELD, LEAF, DEFAULT) \
    TYPE STRUCT_FIELD{};
  NTUPLE_RECO_LEPTON_FIELD_TABLE(NTUPLE_RECO_LEPTON_TO_STRUCT)
  #undef NTUPLE_RECO_LEPTON_TO_STRUCT
  // bool passLooseId;
  // bool passMediumId;
  bool passVetoId = false;
};
struct jets {
  TLorentzVector jet;
  float time;
  bool passId;
  // bool passLooseId;
  // bool passMediumId;
  // bool passTightId;
  bool isCSVL;
  int ecalNRechits;
  float ecalRechitE;
  float jetChargedEMEnergyFraction;
  float jetNeutralEMEnergyFraction;
  float jetChargedHadronEnergyFraction;
  float jetNeutralHadronEnergyFraction;
  bool jetPassMuFrac;
  float jetPtJESUp;
  float jetPtJESDown;
  float jetEJESUp;
  float jetEJESDown;
  float JecUnc;

  float electronEnergyFraction;
  float neutralEmEnergyFraction;
  float chargedHadronEnergyFraction;
  float neutralHadronEnergyFraction;
  float muonEnergyFraction;
};

//lepton highest pt comparator
struct largest_pt {
  inline bool operator()(const leptons& p1, const leptons& p2) {
    return p1.lepton.Pt() > p2.lepton.Pt();
  }
} my_largest_pt;

//jet highest pt comparator
struct largest_pt_jet {
  inline bool operator()(const jets& p1, const jets& p2) {
    return p1.jet.Pt() > p2.jet.Pt();
  }
} my_largest_pt_jet;
/* #endregion */

constexpr int kUnclassifiedClusterId = -1; //? is this redundant with ntuple_branches.def?
using SignalPointKey = std::pair<int, std::string>;

namespace {

  struct EventSynthesis {
    std::array<bool, 10> elePassCutBasedIDTight{};
    std::array<bool, 10> elePassCutBasedIDVeto{};
    std::array<float, 150> jetE{};
    std::array<bool, 150> jetPassIDTightLepVeto{};
    std::array<bool, 150> jetPassIDTight{};
    std::array<float, 50> muonE{};
  };

  struct EventCutState {
    bool flagGoodVertices = false;
    bool flagGlobalSuperTightHalo2016Filter = false;
    bool flagEcalDeadCellTriggerPrimitiveFilter = false;
    bool flagBadPFMuonFilter = false;
    bool flagBadPFMuonDzFilter = false;
    bool flagHfNoisyHitsFilter = false;
    bool flagEeBadScFilter = false;
    bool flagAll = false;
    bool flagEcalBadCalibFilter = false;
    bool jetVeto = true;
    bool hltCscClusterLoose = false;
    bool hltL1CSCShowerDTCluster50 = false;
  };

  struct JetStageResult {
    std::vector<jets> selectedJets;
    float metXJesDown = 0.0f;
    float metYJesDown = 0.0f;
    float metXJesUp = 0.0f;
    float metYJesUp = 0.0f;
  };

  struct GLLPMatchResult {
    float minDeltaR = 15.0f;
    int index = INDEX_DEFAULT;
  };

  EventSynthesis buildEventSynthesis(
      RazorHelper* helper,
      const std::string& analysisTag,
      int nElectron,
      const UChar_t* electronCutBased,
      int nJet,
      const float* jetEta,
      const float* jetPt,
      const float* jetMass,
      const float* jetNeHEF,
      const float* jetNeEmEF,
      const float* jetChEmEF,
      const float* jetMuEF,
      const float* jetChHEF,
      const UChar_t* jetChMultiplicity,
      const UChar_t* jetNeMultiplicity,
      const UChar_t* jetId,
      int nMuon,
      const float* muonEta,
      const float* muonPt) {
    EventSynthesis synth;

    for (int i = 0; i < nElectron; i++) {
      synth.elePassCutBasedIDTight[i] = electronCutBased[i] >= 4;
      synth.elePassCutBasedIDVeto[i] = electronCutBased[i] >= 1;
    }

    for (int i = 0; i < nJet; ++i) {
      auto eta = jetEta[i];
      auto pt = jetPt[i];
      auto pz = pt * TMath::SinH(eta);
      auto mass = jetMass[i];
      synth.jetE[i] = TMath::Sqrt(mass * mass + pt * pt + pz * pz);
    }

    for (int i = 0; i < nJet; ++i) {
      synth.jetPassIDTight[i] = helper->jetTightLepVeto(analysisTag, false, jetNeHEF[i], jetNeEmEF[i], jetChEmEF[i], jetMuEF[i], jetChHEF[i], jetChMultiplicity[i], jetNeMultiplicity[i], jetEta[i], jetId[i]);
      synth.jetPassIDTightLepVeto[i] = helper->jetTightLepVeto(analysisTag, true, jetNeHEF[i], jetNeEmEF[i], jetChEmEF[i], jetMuEF[i], jetChHEF[i], jetChMultiplicity[i], jetNeMultiplicity[i], jetEta[i], jetId[i]);
    }

    for (int i = 0; i < nMuon; ++i) {
      auto eta = muonEta[i];
      auto pt = muonPt[i];
      auto pz = pt * TMath::SinH(eta);
      auto mass = MU_MASS;
      synth.muonE[i] = TMath::Sqrt(mass * mass + pt * pt + pz * pz);
    }
    return synth;
  }

  void fillLeptonBranches(
      TreeMuonSystemCombination* muonSystem,
      const std::vector<leptons>& leptons) {
    for (const auto& tmp : leptons) {
      muonSystem->lepE[muonSystem->nLeptons] = tmp.lepton.E();
      muonSystem->lepPt[muonSystem->nLeptons] = tmp.lepton.Pt();
      muonSystem->lepEta[muonSystem->nLeptons] = tmp.lepton.Eta();
      muonSystem->lepPhi[muonSystem->nLeptons] = tmp.lepton.Phi();
      #define NTUPLE_RECO_LEPTON_TO_FILL(TYPE, BRANCH, STRUCT_FIELD, LEAF, DEFAULT) \
        muonSystem->BRANCH[muonSystem->nLeptons] = tmp.STRUCT_FIELD;
      NTUPLE_RECO_LEPTON_FIELD_TABLE(NTUPLE_RECO_LEPTON_TO_FILL)
      #undef NTUPLE_RECO_LEPTON_TO_FILL
      muonSystem->nLeptons++;
    }
  }

  std::vector<leptons> buildSelectedLeptons(
      RazorAnalyzerMerged& analyzer,
      const EventSynthesis& synth,
      int nMuon,
      const Bool_t* muonLooseId,
      const float* muonPt,
      const float* muonEta,
      const float* muonPhi,
      const Int_t* muonCharge,
      const float* muonDz,
      const Bool_t* muonTightId,
      const float* muonPfRelIso04All,
      int nElectron,
      const float* electronPt,
      const float* electronEta,
      const float* electronPhi,
      const Int_t* electronCharge,
      const float* electronDz) {
    std::vector<leptons> leptonsOut;

    //-------------------------------
    // Muons
    //-------------------------------
    for (int i = 0; i < nMuon; i++) {
      if (!muonLooseId[i])
        continue;
      if (muonPt[i] < 25)
        continue;
      if (fabs(muonEta[i]) > 2.4)
        continue;

      // remove overlaps
      bool overlap = false;
      for (const auto& lep : leptonsOut) {
        if (analyzer.deltaR(muonEta[i], muonPhi[i], lep.lepton.Eta(), lep.lepton.Phi()) < 0.3)
          overlap = true;
      }
      if (overlap)
        continue;

      leptons tmpMuon;
      tmpMuon.lepton.SetPtEtaPhiM(muonPt[i], muonEta[i], muonPhi[i], MU_MASS);
      tmpMuon.pdgId = 13 * -1 * muonCharge[i];
      tmpMuon.dZ = muonDz[i];
      tmpMuon.passId = muonTightId[i];

      tmpMuon.passLooseIso = muonPfRelIso04All[i] < 0.25;
      tmpMuon.passTightIso = muonPfRelIso04All[i] < 0.15;
      tmpMuon.passVTightIso = muonPfRelIso04All[i] < 0.10;
      tmpMuon.passVVTightIso = muonPfRelIso04All[i] < 0.05;

      tmpMuon.passVetoId = false;
      leptonsOut.push_back(tmpMuon);
    }

    //-------------------------------
    // Electrons
    //-------------------------------
    for (int i = 0; i < nElectron; i++) {
      if (!synth.elePassCutBasedIDVeto[i])
        continue;
      if (electronPt[i] < 35)
        continue;
      if (fabs(electronEta[i]) > 2.5)
        continue;

      // remove overlaps
      bool overlap = false;
      for (const auto& lep : leptonsOut) {
        if (analyzer.deltaR(electronEta[i], electronPhi[i], lep.lepton.Eta(), lep.lepton.Phi()) < 0.3)
          overlap = true;
      }
      if (overlap)
        continue;

      leptons tmpElectron;
      tmpElectron.lepton.SetPtEtaPhiM(electronPt[i], electronEta[i], electronPhi[i], ELE_MASS);
      tmpElectron.pdgId = 11 * -1 * electronCharge[i];
      tmpElectron.dZ = electronDz[i];
      tmpElectron.passId = synth.elePassCutBasedIDTight[i];
      leptonsOut.push_back(tmpElectron);
    }

    sort(leptonsOut.begin(), leptonsOut.end(), my_largest_pt);
    return leptonsOut;
  }

  void fillJetBranches(
      TreeMuonSystemCombination* muonSystem,
      const std::vector<jets>& jets) {
    for (const auto& tmp : jets) {
      if (tmp.jet.Pt() < 30)
        continue;

      muonSystem->jetE[muonSystem->nJets] = tmp.jet.E();
      muonSystem->jetPt[muonSystem->nJets] = tmp.jet.Pt();
      muonSystem->jetPtJESUp[muonSystem->nJets] = tmp.jetPtJESUp;
      muonSystem->jetPtJESDown[muonSystem->nJets] = tmp.jetPtJESDown;
      muonSystem->jetEta[muonSystem->nJets] = tmp.jet.Eta();
      muonSystem->jetPhi[muonSystem->nJets] = tmp.jet.Phi();
      muonSystem->jetTightPassId[muonSystem->nJets] = tmp.passId;

      muonSystem->nJets++;
    }
  }

  JetStageResult buildJetStageResult(
      RazorAnalyzerMerged& analyzer,
      RazorHelper* helper,
      int runNumber,
      int nJet,
      const float* jetEta,
      const float* jetPhi,
      const float* jetPt,
      const EventSynthesis& synth,
      const std::vector<leptons>& leptons) {
    JetStageResult result;

    for (int i = 0; i < nJet; ++i) {
      if (std::fabs(jetEta[i]) >= 3.0f)
        continue;
      if (jetPt[i] < 20)
        continue;
      if (!synth.jetPassIDTight[i] && !synth.jetPassIDTightLepVeto[i])
        continue;

      // Exclude selected leptons from the jet collection.
      double minDeltaR = -1.0;
      for (const auto& lep : leptons) {
        const double thisDeltaR = analyzer.deltaR(
            jetEta[i], jetPhi[i], lep.lepton.Eta(), lep.lepton.Phi());
        if (minDeltaR < 0.0 || thisDeltaR < minDeltaR)
          minDeltaR = thisDeltaR;
      }
      if (minDeltaR > 0.0 && minDeltaR < 0.4)
        continue;

      const TLorentzVector nominalJet =
          analyzer.makeTLorentzVector(jetPt[i], jetEta[i], jetPhi[i], synth.jetE[i]);

      jets outJet;
      outJet.jet = nominalJet;
      outJet.passId = synth.jetPassIDTightLepVeto[i];

      const double unc = helper->getJecUnc(jetPt[i], jetEta[i], runNumber);
      outJet.jetPtJESUp = jetPt[i] * (1 + unc);
      outJet.jetPtJESDown = jetPt[i] * (1 - unc);
      outJet.jetEJESUp = synth.jetE[i] * (1 + unc);
      outJet.jetEJESDown = synth.jetE[i] * (1 - unc);
      outJet.JecUnc = unc;

      const TLorentzVector jetJesUp =
          analyzer.makeTLorentzVector(outJet.jetPtJESUp, jetEta[i], jetPhi[i], outJet.jetEJESUp);
      const TLorentzVector jetJesDown =
          analyzer.makeTLorentzVector(outJet.jetPtJESDown, jetEta[i], jetPhi[i], outJet.jetEJESDown);

      result.metXJesUp += -1.0f * static_cast<float>(jetJesUp.Px() - nominalJet.Px());
      result.metYJesUp += -1.0f * static_cast<float>(jetJesUp.Py() - nominalJet.Py());
      result.metXJesDown += -1.0f * static_cast<float>(jetJesDown.Px() - nominalJet.Px());
      result.metYJesDown += -1.0f * static_cast<float>(jetJesDown.Py() - nominalJet.Py());

      result.selectedJets.push_back(outJet);
    }

    sort(result.selectedJets.begin(), result.selectedJets.end(), my_largest_pt_jet);
    return result;
  }

  void fillPuppiMetJesFromShift(
      RazorAnalyzerMerged& analyzer,
      TreeMuonSystemCombination* muonSystem,
      const JetStageResult& jetStage) {
    const TLorentzVector puppiMetVec =
        analyzer.makeTLorentzVectorPtEtaPhiM(muonSystem->PuppiMET_pt, 0, muonSystem->PuppiMET_phi, 0);

    const float metXJesUp = static_cast<float>(puppiMetVec.Px()) + jetStage.metXJesUp;
    const float metYJesUp = static_cast<float>(puppiMetVec.Py()) + jetStage.metYJesUp;
    muonSystem->PuppimetJESUp = std::sqrt(std::pow(metXJesUp, 2) + std::pow(metYJesUp, 2));
    muonSystem->PuppimetPhiJESUp = std::atan(metYJesUp / metXJesUp);
    if (metXJesUp < 0.0f)
      muonSystem->PuppimetPhiJESUp =
          analyzer.deltaPhi(TMath::Pi() + muonSystem->PuppimetPhiJESUp, 0.0);

    const float metXJesDown = static_cast<float>(puppiMetVec.Px()) + jetStage.metXJesDown;
    const float metYJesDown = static_cast<float>(puppiMetVec.Py()) + jetStage.metYJesDown;
    muonSystem->PuppimetJESDown = std::sqrt(std::pow(metXJesDown, 2) + std::pow(metYJesDown, 2));
    muonSystem->PuppimetPhiJESDown = std::atan(metYJesDown / metXJesDown);
    if (metXJesDown < 0.0f)
      muonSystem->PuppimetPhiJESDown =
          analyzer.deltaPhi(TMath::Pi() + muonSystem->PuppimetPhiJESDown, 0.0);
  }

  void fillRawRechits(
      TreeMuonSystemCombination* muonSystem,
      int ncscRechits,
      int ndtRecHits,
      int nrpcRecHits,
      const int* cscRechitsQuality,
      const int* cscRechitsChamber,
      const int* cscRechitsStation,
      const float* cscRechitsEta,
      const float* cscRechitsPhi,
      const float* cscRechitsX,
      const float* cscRechitsY,
      const float* cscRechitsZ,
      const float* cscRechitsTpeak,
      const float* cscRechitsTwire,
      const int* cscRechitsIChamber,
      const int* cscRechitsNStrips,
      const int* cscRechitsWGroupsBX,
      const int* cscRechitsHitWire,
      const int* cscRechitsNWireGroups,
      const float* cscRechitsE,
      const int* dtRecHitsLayer,
      const int* dtRecHitsSuperLayer,
      const int* dtRecHitsStation,
      const int* dtRecHitsWheel,
      const float* dtRecHitsEta,
      const float* dtRecHitsPhi,
      const float* dtRecHitsX,
      const float* dtRecHitsY,
      const float* dtRecHitsZ,
      const int* dtRecHitsSector,
      const int* rpcRecHitsBx,
      const int* rpcRecHitsRegion,
      const int* rpcRecHitsRing,
      const int* rpcRecHitsLayer,
      const int* rpcRecHitsStation,
      const int* rpcRecHitsSector,
      const float* rpcRecHitsX,
      const float* rpcRecHitsY,
      const float* rpcRecHitsZ,
      const float* rpcRecHitsPhi,
      const float* rpcRecHitsEta,
      const float* rpcRecHitsTime,
      const float* rpcRecHitsTimeError) {
    muonSystem->ndtRecHits = ndtRecHits;
    muonSystem->nrpcRecHits = nrpcRecHits;
    muonSystem->ncscRechits = ncscRechits;
    if (muonSystem->ncscRechits > N_MAX_CSCRECHITS)
      muonSystem->ncscRechits = N_MAX_CSCRECHITS;
    if (muonSystem->ndtRecHits > N_MAX_DTRECHITS)
      muonSystem->ndtRecHits = N_MAX_DTRECHITS;
    if (muonSystem->nrpcRecHits > N_MAX_RPCRECHITS)
      muonSystem->nrpcRecHits = N_MAX_RPCRECHITS;

    for (int i = 0; i < muonSystem->ncscRechits; i++) {
      muonSystem->cscRechits_ClusterId[i] = INDEX_DEFAULT;
      muonSystem->cscRechits_Quality[i] = cscRechitsQuality[i];
      muonSystem->cscRechits_Chamber[i] = cscRechitsChamber[i];
      muonSystem->cscRechits_Station[i] = cscRechitsStation[i];
      muonSystem->cscRechits_Eta[i] = cscRechitsEta[i];
      muonSystem->cscRechits_Phi[i] = cscRechitsPhi[i];
      muonSystem->cscRechits_X[i] = cscRechitsX[i];
      muonSystem->cscRechits_Y[i] = cscRechitsY[i];
      muonSystem->cscRechits_Z[i] = cscRechitsZ[i];
      muonSystem->cscRechits_Tpeak[i] = cscRechitsTpeak[i];
      muonSystem->cscRechits_Twire[i] = cscRechitsTwire[i];

      muonSystem->cscRechits_IChamber[i] = cscRechitsIChamber[i];
      muonSystem->cscRechits_NStrips[i] = cscRechitsNStrips[i];
      muonSystem->cscRechits_WGroupsBX[i] = cscRechitsWGroupsBX[i];
      muonSystem->cscRechits_HitWire[i] = cscRechitsHitWire[i];
      muonSystem->cscRechits_NWireGroups[i] = cscRechitsNWireGroups[i];
      muonSystem->cscRechits_E[i] = cscRechitsE[i];
    }
    for (int i = 0; i < muonSystem->ndtRecHits; i++) {
      muonSystem->dtRecHits_ClusterId[i] = INDEX_DEFAULT;
      muonSystem->dtRecHits_Layer[i] = dtRecHitsLayer[i];
      muonSystem->dtRecHits_SuperLayer[i] = dtRecHitsSuperLayer[i];
      muonSystem->dtRecHits_Station[i] = dtRecHitsStation[i];
      muonSystem->dtRecHits_Wheel[i] = dtRecHitsWheel[i];
      muonSystem->dtRecHits_Eta[i] = dtRecHitsEta[i];
      muonSystem->dtRecHits_Phi[i] = dtRecHitsPhi[i];
      muonSystem->dtRecHits_X[i] = dtRecHitsX[i];
      muonSystem->dtRecHits_Y[i] = dtRecHitsY[i];
      muonSystem->dtRecHits_Z[i] = dtRecHitsZ[i];
      muonSystem->dtRecHits_Sector[i] = dtRecHitsSector[i];
    }
    for (int i = 0; i < muonSystem->nrpcRecHits; i++) {
      muonSystem->rpcRecHits_ClusterId[i] = INDEX_DEFAULT;
      muonSystem->rpcRecHits_Bx[i] = rpcRecHitsBx[i];
      muonSystem->rpcRecHits_Region[i] = rpcRecHitsRegion[i];
      muonSystem->rpcRecHits_Ring[i] = rpcRecHitsRing[i];
      muonSystem->rpcRecHits_Layer[i] = rpcRecHitsLayer[i];
      muonSystem->rpcRecHits_Station[i] = rpcRecHitsStation[i];
      muonSystem->rpcRecHits_Sector[i] = rpcRecHitsSector[i];
      muonSystem->rpcRecHits_X[i] = rpcRecHitsX[i];
      muonSystem->rpcRecHits_Y[i] = rpcRecHitsY[i];
      muonSystem->rpcRecHits_Z[i] = rpcRecHitsZ[i];
      muonSystem->rpcRecHits_Phi[i] = rpcRecHitsPhi[i];
      muonSystem->rpcRecHits_Eta[i] = rpcRecHitsEta[i];
      muonSystem->rpcRecHits_Time[i] = rpcRecHitsTime[i];
      muonSystem->rpcRecHits_TimeError[i] = rpcRecHitsTimeError[i];
    }
  }

  int countDtRingsFromRecHits(
      int ndtRecHits,
      const int* dtRecHitsStation,
      const int* dtRecHitsWheel,
      int threshold = 50) {
    int nDTRechitsChamberMinus12 = 0;
    int nDTRechitsChamberMinus11 = 0;
    int nDTRechitsChamber10 = 0;
    int nDTRechitsChamberPlus11 = 0;
    int nDTRechitsChamberPlus12 = 0;
    int nDTRechitsChamberMinus22 = 0;
    int nDTRechitsChamberMinus21 = 0;
    int nDTRechitsChamber20 = 0;
    int nDTRechitsChamberPlus21 = 0;
    int nDTRechitsChamberPlus22 = 0;
    int nDTRechitsChamberMinus32 = 0;
    int nDTRechitsChamberMinus31 = 0;
    int nDTRechitsChamber30 = 0;
    int nDTRechitsChamberPlus31 = 0;
    int nDTRechitsChamberPlus32 = 0;
    int nDTRechitsChamberMinus42 = 0;
    int nDTRechitsChamberMinus41 = 0;
    int nDTRechitsChamber40 = 0;
    int nDTRechitsChamberPlus41 = 0;
    int nDTRechitsChamberPlus42 = 0;

    for (int i = 0; i < ndtRecHits; i++) {
      if (dtRecHitsStation[i] == 1 && dtRecHitsWheel[i] == -2)
        nDTRechitsChamberMinus12++;
      if (dtRecHitsStation[i] == 1 && dtRecHitsWheel[i] == -1)
        nDTRechitsChamberMinus11++;
      if (dtRecHitsStation[i] == 1 && dtRecHitsWheel[i] == 0)
        nDTRechitsChamber10++;
      if (dtRecHitsStation[i] == 1 && dtRecHitsWheel[i] == 1)
        nDTRechitsChamberPlus11++;
      if (dtRecHitsStation[i] == 1 && dtRecHitsWheel[i] == 2)
        nDTRechitsChamberPlus12++;
      if (dtRecHitsStation[i] == 2 && dtRecHitsWheel[i] == -2)
        nDTRechitsChamberMinus22++;
      if (dtRecHitsStation[i] == 2 && dtRecHitsWheel[i] == -1)
        nDTRechitsChamberMinus21++;
      if (dtRecHitsStation[i] == 2 && dtRecHitsWheel[i] == 0)
        nDTRechitsChamber20++;
      if (dtRecHitsStation[i] == 2 && dtRecHitsWheel[i] == 1)
        nDTRechitsChamberPlus21++;
      if (dtRecHitsStation[i] == 2 && dtRecHitsWheel[i] == 2)
        nDTRechitsChamberPlus22++;
      if (dtRecHitsStation[i] == 3 && dtRecHitsWheel[i] == -2)
        nDTRechitsChamberMinus32++;
      if (dtRecHitsStation[i] == 3 && dtRecHitsWheel[i] == -1)
        nDTRechitsChamberMinus31++;
      if (dtRecHitsStation[i] == 3 && dtRecHitsWheel[i] == 0)
        nDTRechitsChamber30++;
      if (dtRecHitsStation[i] == 3 && dtRecHitsWheel[i] == 1)
        nDTRechitsChamberPlus31++;
      if (dtRecHitsStation[i] == 3 && dtRecHitsWheel[i] == 2)
        nDTRechitsChamberPlus32++;
      if (dtRecHitsStation[i] == 4 && dtRecHitsWheel[i] == -2)
        nDTRechitsChamberMinus42++;
      if (dtRecHitsStation[i] == 4 && dtRecHitsWheel[i] == -1)
        nDTRechitsChamberMinus41++;
      if (dtRecHitsStation[i] == 4 && dtRecHitsWheel[i] == 0)
        nDTRechitsChamber40++;
      if (dtRecHitsStation[i] == 4 && dtRecHitsWheel[i] == 1)
        nDTRechitsChamberPlus41++;
      if (dtRecHitsStation[i] == 4 && dtRecHitsWheel[i] == 2)
        nDTRechitsChamberPlus42++;
    }

    int nDtRings = 0;
    if (nDTRechitsChamberMinus12 > threshold)
      nDtRings++;
    if (nDTRechitsChamberMinus11 > threshold)
      nDtRings++;
    if (nDTRechitsChamber10 > threshold)
      nDtRings++;
    if (nDTRechitsChamberPlus11 > threshold)
      nDtRings++;
    if (nDTRechitsChamberPlus12 > threshold)
      nDtRings++;
    if (nDTRechitsChamberMinus22 > threshold)
      nDtRings++;
    if (nDTRechitsChamberMinus21 > threshold)
      nDtRings++;
    if (nDTRechitsChamber20 > threshold)
      nDtRings++;
    if (nDTRechitsChamberPlus21 > threshold)
      nDtRings++;
    if (nDTRechitsChamberPlus22 > threshold)
      nDtRings++;
    if (nDTRechitsChamberMinus32 > threshold)
      nDtRings++;
    if (nDTRechitsChamberMinus31 > threshold)
      nDtRings++;
    if (nDTRechitsChamber30 > threshold)
      nDtRings++;
    if (nDTRechitsChamberPlus31 > threshold)
      nDtRings++;
    if (nDTRechitsChamberPlus32 > threshold)
      nDtRings++;
    if (nDTRechitsChamberMinus42 > threshold)
      nDtRings++;
    if (nDTRechitsChamberMinus41 > threshold)
      nDtRings++;
    if (nDTRechitsChamber40 > threshold)
      nDtRings++;
    if (nDTRechitsChamberPlus41 > threshold)
      nDtRings++;
    if (nDTRechitsChamberPlus42 > threshold)
      nDtRings++;

    return nDtRings;
  }

  int countCscRingsFromRecHits(
      int ncscRechits,
      const int* cscRechitsChamber,
      int threshold = 50) {
    int nCscRechitsChamberPlus11 = 0;
    int nCscRechitsChamberPlus12 = 0;
    int nCscRechitsChamberPlus13 = 0;
    int nCscRechitsChamberPlus21 = 0;
    int nCscRechitsChamberPlus22 = 0;
    int nCscRechitsChamberPlus31 = 0;
    int nCscRechitsChamberPlus32 = 0;
    int nCscRechitsChamberPlus41 = 0;
    int nCscRechitsChamberPlus42 = 0;

    int nCscRechitsChamberMinus11 = 0;
    int nCscRechitsChamberMinus12 = 0;
    int nCscRechitsChamberMinus13 = 0;
    int nCscRechitsChamberMinus21 = 0;
    int nCscRechitsChamberMinus22 = 0;
    int nCscRechitsChamberMinus31 = 0;
    int nCscRechitsChamberMinus32 = 0;
    int nCscRechitsChamberMinus41 = 0;
    int nCscRechitsChamberMinus42 = 0;

    for (int i = 0; i < ncscRechits; i++) {
      if (cscRechitsChamber[i] == 11)
        nCscRechitsChamberPlus11++;
      if (cscRechitsChamber[i] == 12)
        nCscRechitsChamberPlus12++;
      if (cscRechitsChamber[i] == 13)
        nCscRechitsChamberPlus13++;
      if (cscRechitsChamber[i] == 21)
        nCscRechitsChamberPlus21++;
      if (cscRechitsChamber[i] == 22)
        nCscRechitsChamberPlus22++;
      if (cscRechitsChamber[i] == 31)
        nCscRechitsChamberPlus31++;
      if (cscRechitsChamber[i] == 32)
        nCscRechitsChamberPlus32++;
      if (cscRechitsChamber[i] == 41)
        nCscRechitsChamberPlus41++;
      if (cscRechitsChamber[i] == 42)
        nCscRechitsChamberPlus42++;
      if (cscRechitsChamber[i] == -11)
        nCscRechitsChamberMinus11++;
      if (cscRechitsChamber[i] == -12)
        nCscRechitsChamberMinus12++;
      if (cscRechitsChamber[i] == -13)
        nCscRechitsChamberMinus13++;
      if (cscRechitsChamber[i] == -21)
        nCscRechitsChamberMinus21++;
      if (cscRechitsChamber[i] == -22)
        nCscRechitsChamberMinus22++;
      if (cscRechitsChamber[i] == -31)
        nCscRechitsChamberMinus31++;
      if (cscRechitsChamber[i] == -32)
        nCscRechitsChamberMinus32++;
      if (cscRechitsChamber[i] == -41)
        nCscRechitsChamberMinus41++;
      if (cscRechitsChamber[i] == -42)
        nCscRechitsChamberMinus42++;
    }

    int nCscRings = 0;
    if (nCscRechitsChamberPlus11 > threshold)
      nCscRings++;
    if (nCscRechitsChamberPlus12 > threshold)
      nCscRings++;
    if (nCscRechitsChamberPlus13 > threshold)
      nCscRings++;
    if (nCscRechitsChamberPlus21 > threshold)
      nCscRings++;
    if (nCscRechitsChamberPlus22 > threshold)
      nCscRings++;
    if (nCscRechitsChamberPlus31 > threshold)
      nCscRings++;
    if (nCscRechitsChamberPlus32 > threshold)
      nCscRings++;
    if (nCscRechitsChamberPlus41 > threshold)
      nCscRings++;
    if (nCscRechitsChamberPlus42 > threshold)
      nCscRings++;
    if (nCscRechitsChamberMinus11 > threshold)
      nCscRings++;
    if (nCscRechitsChamberMinus12 > threshold)
      nCscRings++;
    if (nCscRechitsChamberMinus13 > threshold)
      nCscRings++;
    if (nCscRechitsChamberMinus21 > threshold)
      nCscRings++;
    if (nCscRechitsChamberMinus22 > threshold)
      nCscRings++;
    if (nCscRechitsChamberMinus31 > threshold)
      nCscRings++;
    if (nCscRechitsChamberMinus32 > threshold)
      nCscRings++;
    if (nCscRechitsChamberMinus41 > threshold)
      nCscRings++;
    if (nCscRechitsChamberMinus42 > threshold)
      nCscRings++;

    return nCscRings;
  }

  std::vector<Rechits> buildCscRechitPoints(
      int ncscRechits,
      const float* cscRechitsPhi,
      const float* cscRechitsEta,
      const float* cscRechitsX,
      const float* cscRechitsY,
      const float* cscRechitsZ,
      const float* cscRechitsTpeak,
      const float* cscRechitsTwire,
      const int* cscRechitsStation,
      const int* cscRechitsChamber) {
    std::vector<Rechits> points;
    points.reserve(ncscRechits);

    for (int i = 0; i < ncscRechits; i++) {
      int layer = 0;
      Rechits p;
      p.phi = cscRechitsPhi[i];
      p.eta = cscRechitsEta[i];
      p.x = cscRechitsX[i];
      p.y = cscRechitsY[i];
      p.z = cscRechitsZ[i];
      p.t = cscRechitsTpeak[i];
      p.twire = cscRechitsTwire[i];
      p.station = cscRechitsStation[i];
      p.chamber = cscRechitsChamber[i];
      p.layer = layer;
      p.superlayer = 0;
      p.wheel = 0;
      p.clusterID = kUnclassifiedClusterId;
      points.push_back(p);
    }

    return points;
  }

  std::vector<Rechits> buildDtRecHitPoints(
      int ndtRecHits,
      const float* dtRecHitsPhi,
      const float* dtRecHitsEta,
      const float* dtRecHitsX,
      const float* dtRecHitsY,
      const float* dtRecHitsZ,
      const int* dtRecHitsStation,
      const int* dtRecHitsWheel,
      const int* dtRecHitsSuperLayer) {
    std::vector<Rechits> points;
    points.reserve(ndtRecHits);

    for (int i = 0; i < ndtRecHits; i++) {
      Rechits p;
      p.phi = dtRecHitsPhi[i];
      p.eta = dtRecHitsEta[i];
      p.x = dtRecHitsX[i];
      p.y = dtRecHitsY[i];
      p.z = dtRecHitsZ[i];
      p.t = -999.;
      p.twire = -999.;
      p.station = dtRecHitsStation[i];
      p.chamber = dtRecHitsWheel[i];
      p.superlayer = dtRecHitsSuperLayer[i];
      p.wheel = dtRecHitsWheel[i];
      p.clusterID = kUnclassifiedClusterId;
      points.push_back(p);
    }

    return points;
  }

  std::vector<Rechits> buildRpcRecHitPoints(
      int nrpcRecHits,
      const float* rpcRecHitsPhi,
      const float* rpcRecHitsEta,
      const float* rpcRecHitsX,
      const float* rpcRecHitsY,
      const float* rpcRecHitsZ,
      const float* rpcRecHitsTime,
      const int* rpcRecHitsStation,
      const int* rpcRecHitsSector,
      const int* rpcRecHitsLayer,
      const int* rpcRecHitsRing) {
    std::vector<Rechits> points;
    points.reserve(nrpcRecHits);

    for (int i = 0; i < nrpcRecHits; i++) {
      Rechits p;
      p.phi = rpcRecHitsPhi[i];
      p.eta = rpcRecHitsEta[i];
      p.x = rpcRecHitsX[i];
      p.y = rpcRecHitsY[i];
      p.z = rpcRecHitsZ[i];
      p.t = rpcRecHitsTime[i];
      p.twire = rpcRecHitsTime[i];
      p.station = rpcRecHitsStation[i];
      p.chamber = rpcRecHitsSector[i];
      p.layer = rpcRecHitsLayer[i];
      p.superlayer = 0;
      p.wheel = rpcRecHitsRing[i];
      p.clusterID = kUnclassifiedClusterId;
      points.push_back(p);
    }

    return points;
  }

  void runCAClustering(CACluster& clusterer) {
    clusterer.run();
    clusterer.clusterProperties();
    clusterer.sort_clusters();
  }

  void fillClusteredRecHitIds(
      int nStoredRecHits,
      const std::vector<Rechits>& clusteredPoints,
      int* outClusterIds) {
    int nPointsOut = nStoredRecHits;
    if (nPointsOut > static_cast<int>(clusteredPoints.size()))
      nPointsOut = clusteredPoints.size();
    for (int i = 0; i < nPointsOut; i++) {
      outClusterIds[i] = clusteredPoints[i].clusterID;
    }
  }

  EventCutState buildEventCutState(
      RazorAnalyzerMerged& analyzer,
      RazorHelper* helper,
      const std::string& analysisTag,
      bool isData,
      int runNumber,
      float puppiMetPt,
      float puppiMetPhi,
      int nJet,
      const float* jetPt,
      const float* jetEta,
      const float* jetPhi,
      const float* jetNeEmEF,
      const float* jetChEmEF,
      int nMuon,
      const Bool_t* muonIsPFcand,
      const float* muonEta,
      const float* muonPhi,
      const EventSynthesis& synth,
      Bool_t flagGoodVertices,
      Bool_t flagGlobalSuperTightHalo2016Filter,
      Bool_t flagEcalDeadCellTriggerPrimitiveFilter,
      Bool_t flagBadPFMuonFilter,
      Bool_t flagBadPFMuonDzFilter,
      Bool_t flagHfNoisyHitsFilter,
      Bool_t flagEeBadScFilter,
      Bool_t flagEcalBadCalibFilter,
      Bool_t hltCscClusterLoose,
      Bool_t hltL1CSCShowerDTCluster50) {
    EventCutState cuts;

    // noise filters
    cuts.flagGoodVertices = flagGoodVertices;
    cuts.flagGlobalSuperTightHalo2016Filter = flagGlobalSuperTightHalo2016Filter;
    cuts.flagEcalDeadCellTriggerPrimitiveFilter = flagEcalDeadCellTriggerPrimitiveFilter;
    cuts.flagBadPFMuonFilter = flagBadPFMuonFilter;
    cuts.flagBadPFMuonDzFilter = flagBadPFMuonDzFilter;
    cuts.flagHfNoisyHitsFilter = flagHfNoisyHitsFilter;
    cuts.flagEeBadScFilter = flagEeBadScFilter;
    cuts.flagAll = flagEeBadScFilter && flagHfNoisyHitsFilter &&
                   flagBadPFMuonDzFilter && flagBadPFMuonFilter &&
                   flagEcalDeadCellTriggerPrimitiveFilter &&
                   flagGlobalSuperTightHalo2016Filter && flagGoodVertices;
    if (analysisTag == "Summer24")
      cuts.flagEcalBadCalibFilter = flagEcalBadCalibFilter;

    // Flag_ecalBadCalibFilter for nanoAOD:
    // https://twiki.cern.ch/twiki/bin/view/CMS/MissingETOptionalFiltersRun2#ECal_BadCalibration_Filter_Flag
    if (analysisTag == "Summer24")
      cuts.flagEcalBadCalibFilter = flagEcalBadCalibFilter;
    else {
      cuts.flagEcalBadCalibFilter = true;
      if (isData && runNumber >= 362433 && runNumber <= 367144) {
        if (puppiMetPt > 100) {
          for (int i = 0; i < nJet; i++) {
            if (jetPt[i] < 50)
              continue;
            if (!(jetEta[i] <= -0.1f && jetEta[i] >= -0.5f && jetPhi[i] < -1.8f && jetPhi[i] > -2.1f))
              continue;
            if (!(jetNeEmEF[i] > 0.9f || jetChEmEF[i] > 0.9f))
              continue;
            if (analyzer.deltaPhi(puppiMetPhi, jetPhi[i]) < 2.9)
              continue;
            cuts.flagEcalBadCalibFilter = false;
          }
        }
      }
    }

    // jet veto map, following selections here:
    // https://cms-jerc.web.cern.ch/Recommendations/#jet-veto-maps
    cuts.jetVeto = true;
    for (int i = 0; i < nJet; i++) {
      if (jetPt[i] <= 15)
        continue;
      if (jetNeEmEF[i] + jetChEmEF[i] >= 0.9f)
        continue;
      if (!synth.jetPassIDTight[i])
        continue;
      //remove overlaps
      bool overlap = false;
      for (int j = 0; j < nMuon; j++) {
        if (!muonIsPFcand[j])
          continue;
        if (analyzer.deltaR(jetEta[i], jetPhi[i], muonEta[j], muonPhi[j]) < 0.2)
          overlap = true;
      }
      if (overlap)
        continue;
      helper->getJetVetoMap(0, 1);
      if (helper->getJetVetoMap(jetEta[i], jetPhi[i]) > 0.0)
        cuts.jetVeto = false;
      if (analysisTag == "Summer24" && helper->getJetVetoFpixMap(jetEta[i], jetPhi[i]) > 0.0)
        cuts.jetVeto = false;
    }

    cuts.hltCscClusterLoose = hltCscClusterLoose;
    cuts.hltL1CSCShowerDTCluster50 = hltL1CSCShowerDTCluster50;
    return cuts;
  }

  void writeEventCutState(TreeMuonSystemCombination* muonSystem, const EventCutState& cuts) {
    muonSystem->Flag_goodVertices = cuts.flagGoodVertices;
    muonSystem->Flag_globalSuperTightHalo2016Filter = cuts.flagGlobalSuperTightHalo2016Filter;
    muonSystem->Flag_EcalDeadCellTriggerPrimitiveFilter = cuts.flagEcalDeadCellTriggerPrimitiveFilter;
    muonSystem->Flag_BadPFMuonFilter = cuts.flagBadPFMuonFilter;
    muonSystem->Flag_BadPFMuonDzFilter = cuts.flagBadPFMuonDzFilter;
    muonSystem->Flag_hfNoisyHitsFilter = cuts.flagHfNoisyHitsFilter;
    muonSystem->Flag_eeBadScFilter = cuts.flagEeBadScFilter;
    muonSystem->Flag_all = cuts.flagAll;
    muonSystem->Flag_ecalBadCalibFilter = cuts.flagEcalBadCalibFilter;
    muonSystem->jetVeto = cuts.jetVeto;
    muonSystem->HLT_CscCluster_Loose = cuts.hltCscClusterLoose;
    muonSystem->HLT_L1CSCShower_DTCluster50 = cuts.hltL1CSCShowerDTCluster50;
  }

  void fillClusterJetVeto(
      RazorAnalyzerMerged& analyzer,
      RazorHelper* helper,
      int runNumber,
      int nJet,
      const float* jetEta,
      const float* jetPhi,
      const float* jetPt,
      const float* jetE,
      const bool* jetPassIDTight,
      float clusterEta,
      float clusterPhi,
      float& outJetVetoPt,
      float& outJetVetoE,
      bool& outJetVetoTightId,
      float& outJetVetoPtJESUp,
      float& outJetVetoPtJESDown) {
    outJetVetoPt = 0.0f;
    outJetVetoE = 0.0f;
    outJetVetoTightId = false;
    outJetVetoPtJESUp = 0.0f;
    outJetVetoPtJESDown = 0.0f;

    for (int i = 0; i < nJet; ++i) {
      if (std::fabs(jetEta[i]) > 3.0f)
        continue;
      if (analyzer.deltaR(jetEta[i], jetPhi[i], clusterEta, clusterPhi) < 0.4 &&
          jetPt[i] > outJetVetoPt) {
        outJetVetoPt = jetPt[i];
        outJetVetoE = jetE[i];
        outJetVetoTightId = jetPassIDTight[i];

        const double unc = helper->getJecUnc(jetPt[i], jetEta[i], runNumber);
        outJetVetoPtJESUp = jetPt[i] * (1 + unc);
        outJetVetoPtJESDown = jetPt[i] * (1 - unc);
      }
    }
  }

  void fillMatchedGLLPFields(
      TreeMuonSystemCombination* muonSystem,
      int clusterIndex,
      int gllpIndex,
      float* outEta,
      float* outPhi,
      float* outDecayR,
      float* outDecayZ,
      bool* outCsc,
      bool* outDt,
      float* outE) {
    outEta[clusterIndex] = muonSystem->gLLP_eta[gllpIndex];
    outPhi[clusterIndex] = muonSystem->gLLP_phi[gllpIndex];
    outDecayR[clusterIndex] = muonSystem->gLLP_decay_vertex_r[gllpIndex];
    outDecayZ[clusterIndex] = muonSystem->gLLP_decay_vertex_z[gllpIndex];
    outCsc[clusterIndex] = muonSystem->gLLP_csc[gllpIndex];
    outDt[clusterIndex] = muonSystem->gLLP_dt[gllpIndex];
    outE[clusterIndex] = muonSystem->gLLP_e[gllpIndex];
  }

  GLLPMatchResult findNearestGLLPMatch(
      RazorAnalyzerMerged& analyzer,
      float clusterEta,
      float clusterPhi,
      int nGLLP,
      const float* gLLPEta,
      const float* gLLPPhi) {
    GLLPMatchResult result;
    for (int j = 0; j < nGLLP; ++j) {
      const float currentDeltaR =
          static_cast<float>(analyzer.deltaR(clusterEta, clusterPhi, gLLPEta[j], gLLPPhi[j]));
      if (currentDeltaR < result.minDeltaR) {
        result.minDeltaR = currentDeltaR;
        result.index = j;
      }
    }
    return result;
  }

} // namespace

bool parseSignalPointFromInputPath(const std::string& inputPath, int& mh, int& mx, float& ctau, std::string& ctauTokenTag) {
  auto extractTokenAfterKey = [&inputPath](const std::string& key, const std::string& allowedChars, std::string& token) -> bool {
    const size_t keyPos = inputPath.find(key);
    if (keyPos == std::string::npos)
      return false;
    const size_t valueStart = keyPos + key.size();
    if (valueStart >= inputPath.size())
      return false;
    size_t valueEnd = valueStart;
    while (valueEnd < inputPath.size() && allowedChars.find(inputPath[valueEnd]) != std::string::npos) {
      ++valueEnd;
    }
    if (valueEnd == valueStart)
      return false;
    token = inputPath.substr(valueStart, valueEnd - valueStart);
    return true;
  };

  auto parseNonNegativeIntToken = [](const std::string& token, int& out) -> bool {
    if (token.empty())
      return false;
    for (const char c : token) {
      if (!std::isdigit(static_cast<unsigned char>(c)))
        return false;
    }
    try {
      const long long parsed = std::stoll(token);
      if (parsed < 0 || parsed > std::numeric_limits<int>::max())
        return false;
      out = static_cast<int>(parsed);
    } catch (...) {
      return false;
    }
    return true;
  };

  auto parseCtauTokenToFloat = [](const std::string& token, float& out) -> bool {
    if (token.empty())
      return false;
    std::string normalized = token;
    std::replace(normalized.begin(), normalized.end(), 'p', '.');
    std::replace(normalized.begin(), normalized.end(), 'P', '.');
    try {
      size_t parsedChars = 0;
      const float parsedValue = std::stof(normalized, &parsedChars);
      if (parsedChars != normalized.size() || parsedValue < 0.0f ||
          parsedValue > std::numeric_limits<float>::max()) {
        return false;
      }
      out = parsedValue;
    } catch (...) {
      return false;
    }
    return true;
  };

  std::string mhToken;
  std::string mxToken;
  std::string ctauToken;
  if (!extractTokenAfterKey("MH-", "0123456789", mhToken))
    return false;
  if (!extractTokenAfterKey("MS-", "0123456789", mxToken))
    return false;
  if (!extractTokenAfterKey("ctauS-", "0123456789pP.", ctauToken))
    return false;

  if (!parseNonNegativeIntToken(mhToken, mh))
    return false;
  if (!parseNonNegativeIntToken(mxToken, mx))
    return false;

  ctauTokenTag = ctauToken;
  std::replace(ctauTokenTag.begin(), ctauTokenTag.end(), 'P', 'p');
  std::replace(ctauTokenTag.begin(), ctauTokenTag.end(), '.', 'p');
  if (!parseCtauTokenToFloat(ctauToken, ctau))
    return false;
  return true;
}

void llp_MuonSystem_CA_mdsnano::Analyze(bool isData, int options, string outputfilename, string analysisTag) {

  /* #region: fChain check and nEntries declaration.*/
  if (fChain == nullptr) {
    cout << "fChain is a null pointer." << endl;
    return;
  };
  const Long64_t nEntries = fChain->GetEntries();
  cout << "Initializing..." << endl;
  cout << "IsData = " << isData << "\n";
  cout << "options = " << options << "\n";
  cout << "Total Events: " << nEntries << "\n";
  /* #endregion */
  
  /* #region: options format: MH/MX/ctau/condor: 1000/300/0/1 */
  // mh can be 3-4 digits, mx is always 3 digits, ctau is one digit(number of zeros), last digit is condor option
  // mh can be 3-4 digits, mx is always 3 digits, ctau is 2 digit(number of zeros), last digit is condor option

  // int mx = int(options/1000)%1000;
  // int mh = options/1000000;
  // int ctau = pow(10, int(options/10)%10) * int(int(options/100)%10);

  // cout<<"mh "<<mh<<", mx "<<mx<<", ctau "<<ctau<<endl;
  /* #endregion */

  /* #region: Startup decode/defaults for this analyzer invocation.*/
  bool signalScan = int(options / 10) == 1;
  int option = options % 10;
  if (analysisTag.empty()) {
    analysisTag = "Summer24";
    std::cout << "[WARNING]: empty analysisTag received, defaulting to " << analysisTag << '\n';
  }
  /* #endregion */

  /* #region: declares lookup tables for signalScan*/
  map<SignalPointKey, std::unique_ptr<TFile>> Files2D;
  map<SignalPointKey, TTree*> Trees2D;
  map<SignalPointKey, TH1F*> NEvents2D;
  map<SignalPointKey, TH1F*> accep2D;
  map<SignalPointKey, TH1F*> accep_met2D;
  map<SignalPointKey, TH1F*> Total2D;
  /* #endregion */

  /* #region: run-mode logging*/
  std::cout << "[INFO]: running on " << (isData ? "data" : "MC")
            << " with option: " << option << '\n';
  std::cout << "[INFO]: running "
            << (signalScan ? "with Signal scan" : "without Signal scan ")
            << (signalScan ? "" : std::to_string(option))
            << '\n';
  /* #endregion */

  /* #region: setting up output file/tree */
  std::unique_ptr<TFile> outFile;
  if (isData || !signalScan) {
    string outfilename = outputfilename.empty() ? "MuonSystem_Tree.root" : outputfilename;
    outFile = std::make_unique<TFile>(outfilename.c_str(), "RECREATE");
    if (!outFile || outFile->IsZombie()) {
      cerr << "[ERROR]: failed to create output file: " << outfilename << endl;
      return;
    }
    outFile->cd();
  }

  TreeMuonSystemCombination* MuonSystem = new TreeMuonSystemCombination;
  MuonSystem->CreateTree();
  MuonSystem->tree_->SetAutoFlush(0);
  MuonSystem->InitTree();
  /* #endregion */

  /* #region: histogram containing total number of processed events (for normalization)*/
  auto NEvents = makeCounterHist("NEvents");
  auto Total = makeCounterHist("Total");

  auto accep = makeCounterHist("accep");
  auto accep_csccsc = makeCounterHist("accep_csccsc");
  auto accep_cscdt = makeCounterHist("accep_cscdt");
  auto accep_met = makeCounterHist("accep_met");

  auto Nmet200 = makeCounterHist("Nmet200");
  auto NmetFilter = makeCounterHist("NmetFilter");
  auto Nlep0 = makeCounterHist("Nlep0");
  auto Njet1 = makeCounterHist("Njet1");
  auto NcosmicVeto = makeCounterHist("NcosmicVeto");
  /* #endregion */

  // initialize helper in memory safe way
  auto helper = std::make_unique<RazorHelper>(analysisTag, isData);
  // [1] Event loop
  /* #region */
  cout << "[INFO]: Loop Starting" << endl;
  auto lastReport = steady_clock::now(); // mr. timer
  std::string cachedSignalSourcePath;
  int cachedSignalMh = 0;
  int cachedSignalMx = 0;
  float cachedSignalCtau = 0.0f;
  std::string cachedSignalCtauTag;
  bool hasCachedSignalPoint = false;
  for (Long64_t jentry = 0; jentry < nEntries; jentry++) {

    // progress logging
    if (jentry % 1000 == 0) {
      const auto now = steady_clock::now();
      const double dt = duration<double>(now - lastReport).count();
      cout << "Processing entry " << jentry << " | dt=" << dt << " s\n";
      lastReport = now;
    }

    // [2] Read event entry and reset output event variables
    /* #region */
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0)
      break;
    fChain->GetEntry(jentry);

    //fill normalization histogram
    MuonSystem->InitVariables();
    /* #endregion */

    // [3] Build local aliases and per-event helper inputs
    /* #region */
    const EventSynthesis synth = buildEventSynthesis(
        helper.get(),
        analysisTag,
        nElectron,
        Electron_cutBased,
        nJet,
        Jet_eta,
        Jet_pt,
        Jet_mass,
        Jet_neHEF,
        Jet_neEmEF,
        Jet_chEmEF,
        Jet_muEF,
        Jet_chHEF,
        Jet_chMultiplicity,
        Jet_neMultiplicity,
        Jet_jetId,
        nMuon,
        Muon_eta,
        Muon_pt);

    /* #endregion */

    // [4] Signal-scan routing (split output by model point)
    /* #region */
    SignalPointKey eventSignalKey = SignalPointKey(0, std::string());
    bool hasEventSignalKey = false;
    if (!isData && signalScan) {
      std::string inputPath;
      if (fChain != nullptr && fChain->GetCurrentFile() != nullptr) {
        inputPath = fChain->GetCurrentFile()->GetName();
      }
      if (!hasCachedSignalPoint || inputPath != cachedSignalSourcePath) {
        std::string ctauTokenTag;
        if (!parseSignalPointFromInputPath(inputPath, cachedSignalMh, cachedSignalMx, cachedSignalCtau,
                                           ctauTokenTag)) {
          std::cerr << "[ERROR]: signalScan requested but failed to parse MH/MS/ctauS from input path: '"
                    << inputPath << "'. Exiting to avoid mislabeled output." << std::endl;
          return;
        }
        cachedSignalCtauTag = ctauTokenTag;
        cachedSignalSourcePath = inputPath;
        hasCachedSignalPoint = true;
      }
      MuonSystem->mH = cachedSignalMh;
      MuonSystem->mX = cachedSignalMx;
      MuonSystem->ctau = cachedSignalCtau;

      eventSignalKey = make_pair(cachedSignalMx, cachedSignalCtauTag);
      hasEventSignalKey = true;

      if (Files2D.count(eventSignalKey) == 0) { //create file and tree
        //format file name
        string thisFileName = outputfilename.empty() ? "MuonSystem_Tree.root" : outputfilename;
        thisFileName.erase(thisFileName.end() - 5, thisFileName.end());
        thisFileName += "_MH-" + to_string(cachedSignalMh)
                     + "_MS-" + to_string(cachedSignalMx)
                     + "_ctauS-" + cachedSignalCtauTag + ".root";

        Files2D[eventSignalKey] = std::make_unique<TFile>(thisFileName.c_str(), "recreate");
        Files2D[eventSignalKey]->cd();
        Trees2D[eventSignalKey] = MuonSystem->tree_->CloneTree(0);
        const std::string signalSuffix = std::to_string(cachedSignalMx) + "_" + cachedSignalCtauTag;
        NEvents2D[eventSignalKey] = new TH1F(("NEvents_" + signalSuffix).c_str(), "NEvents", 1, 0.5, 1.5);
        Total2D[eventSignalKey] = new TH1F(("Total_" + signalSuffix).c_str(), "Total", 1, 0.5, 1.5);
        accep2D[eventSignalKey] = new TH1F(("accep2D_" + signalSuffix).c_str(), "accep", 1, 0.5, 1.5);
        accep_met2D[eventSignalKey] = new TH1F(("accep_met2D_" + signalSuffix).c_str(), "accep_met", 1, 0.5, 1.5);

        cout << "Created new output file " << thisFileName << endl;
      }
      //Fill NEvents hist
      NEvents2D[eventSignalKey]->Fill(1.0, Generator_weight);
    }
    /* #endregion */

    // [5] Event metadata and per-event nominal weight
    /* #region */
    //event info
    if (isData) { //? can this be better organized?
      NEvents->Fill(1);
      MuonSystem->weight = 1;
    } else {
      MuonSystem->weight = Generator_weight;
      NEvents->Fill(1, Generator_weight);
    }
    MuonSystem->run = run; //? is there any particular reason these need to be here?
    MuonSystem->luminosityBlock = luminosityBlock;
    MuonSystem->event = event;

    if (isData && run < 360019) continue; //? should this be moved somewhere higher to save compute?
    /* #endregion */

    // [6] MC-only truth matching inputs and MC weights
    /* #region */
    if (!isData) {
      //for DS model
      /* #region  */
      MuonSystem->nGenParticles = 0;
      for (int i = 0; i < nGenPart; i++) {
        if (abs(GenPart_pdgId[i]) >= 999999) {
          MuonSystem->gParticleEta[MuonSystem->nGenParticles] = GenPart_eta[i];
          MuonSystem->gParticlePhi[MuonSystem->nGenParticles] = GenPart_phi[i];
          MuonSystem->gParticlePt[MuonSystem->nGenParticles] = GenPart_pt[i];
          MuonSystem->gParticleId[MuonSystem->nGenParticles] = GenPart_pdgId[i];
          float decay_x;
          float decay_y;
          float decay_z;
          bool foundDaughter = false;
          for (int j = 0; j < nGenPart; j++) {
            if (GenPart_genPartIdxMother[j] == i) {
              decay_x = GenPart_vx[j];
              decay_y = GenPart_vy[j];
              decay_z = GenPart_vz[j];
              foundDaughter = true;
              break;
            }
          }
          if(foundDaughter){
            MuonSystem->gParticle_decay_vertex_r[MuonSystem->nGenParticles] = sqrt(decay_x * decay_x + decay_y * decay_y);
            MuonSystem->gParticle_decay_vertex_x[MuonSystem->nGenParticles] = decay_x;
            MuonSystem->gParticle_decay_vertex_y[MuonSystem->nGenParticles] = decay_y;
            MuonSystem->gParticle_decay_vertex_z[MuonSystem->nGenParticles] = decay_z;
          } else {
            cout << "didn't find HV LLP daughter " << event << endl;
          }

          TLorentzVector particle = makeTLorentzVectorPtEtaPhiM(GenPart_pt[i], GenPart_eta[i], GenPart_phi[i], GenPart_mass[i]);

          float gParticle_decay_vertex = sqrt(pow(MuonSystem->gParticle_decay_vertex_r[MuonSystem->nGenParticles], 2) + pow(MuonSystem->gParticle_decay_vertex_z[MuonSystem->nGenParticles], 2));

          MuonSystem->gParticle_ctau[MuonSystem->nGenParticles] = gParticle_decay_vertex / (particle.Beta() * particle.Gamma());
          MuonSystem->gParticle_beta[MuonSystem->nGenParticles] = particle.Beta();
          MuonSystem->nGenParticles++;
        }
      }
      /* #endregion */

      //for Twin higgs model
      /* #region  */
      MuonSystem->nGLLP = 0;
      for (int i = 0; i < nGenPart; i++) {
        if (abs(GenPart_pdgId[i]) == 9000006) {
          MuonSystem->gLLP_eta[MuonSystem->nGLLP] = GenPart_eta[i];
          MuonSystem->gLLP_phi[MuonSystem->nGLLP] = GenPart_phi[i];
          MuonSystem->gLLP_pt[MuonSystem->nGLLP] = GenPart_pt[i];
          float decay_x;
          float decay_y;
          float decay_z;
          bool foundDaughter = false;
          for (int j = 0; j < nGenPart; j++) {
            if (GenPart_genPartIdxMother[j] == i) {
              decay_x = GenPart_vx[j];
              decay_y = GenPart_vy[j];
              decay_z = GenPart_vz[j];
              foundDaughter = true;
              break;
            }
          }
          if (foundDaughter){
            MuonSystem->gLLP_decay_vertex_r[MuonSystem->nGLLP] = sqrt(decay_x * decay_x + decay_y * decay_y);
            MuonSystem->gLLP_decay_vertex_x[MuonSystem->nGLLP] = decay_x;
            MuonSystem->gLLP_decay_vertex_y[MuonSystem->nGLLP] = decay_y;
            MuonSystem->gLLP_decay_vertex_z[MuonSystem->nGLLP] = decay_z;
          } else {
            cout << "didn't find LLP daughter " << event << endl;
          }

          TLorentzVector LLP = makeTLorentzVectorPtEtaPhiM(GenPart_pt[i], GenPart_eta[i], GenPart_phi[i], GenPart_mass[i]);
          MuonSystem->gLLP_e[MuonSystem->nGLLP] = LLP.E();

          float gLLP_decay_vertex = sqrt(pow(MuonSystem->gLLP_decay_vertex_r[MuonSystem->nGLLP], 2) + pow(MuonSystem->gLLP_decay_vertex_z[MuonSystem->nGLLP], 2));

          MuonSystem->gLLP_ctau[MuonSystem->nGLLP] = gLLP_decay_vertex / (LLP.Beta() * LLP.Gamma());
          MuonSystem->gLLP_beta[MuonSystem->nGLLP] = LLP.Beta();
          MuonSystem->gLLP_csc[MuonSystem->nGLLP] =
              abs(MuonSystem->gLLP_eta[MuonSystem->nGLLP]) < 2.4 &&
              abs(MuonSystem->gLLP_decay_vertex_z[MuonSystem->nGLLP]) < 1100 &&
              abs(MuonSystem->gLLP_decay_vertex_z[MuonSystem->nGLLP]) > 400 &&
              MuonSystem->gLLP_decay_vertex_r[MuonSystem->nGLLP] < 695.5;
          MuonSystem->gLLP_dt[MuonSystem->nGLLP] =
              abs(MuonSystem->gLLP_decay_vertex_z[MuonSystem->nGLLP]) < 661.0 &&
              MuonSystem->gLLP_decay_vertex_r[MuonSystem->nGLLP] < 800 &&
              MuonSystem->gLLP_decay_vertex_r[MuonSystem->nGLLP] > 200.0;
          MuonSystem->nGLLP++;
        }
      }
      /* #endregion */

      MuonSystem->Pileup_nTrueInt = Pileup_nTrueInt;
      MuonSystem->pileupWeight = helper->getPileupWeight(Pileup_nTrueInt);
      MuonSystem->pileupWeightUp = helper->getPileupWeightUp(Pileup_nTrueInt) / MuonSystem->pileupWeight;
      MuonSystem->pileupWeightDown = helper->getPileupWeightDown(Pileup_nTrueInt) / MuonSystem->pileupWeight;

      for (unsigned int i = 0; i < 9; i++) {
        MuonSystem->LHEScaleWeight[i] = LHEScaleWeight[i];
      }

    }
    /* #endregion */

    // [7] Event-level observables and acceptance bookkeeping
    /* #region */
    //get NPU
    MuonSystem->PV_npvs = PV_npvs;
    MuonSystem->Rho_fixedGridRhoFastjetAll = Rho_fixedGridRhoFastjetAll;
    MuonSystem->PFMET_pt = PFMET_pt;
    MuonSystem->PFMET_phi = PFMET_phi;

    MuonSystem->PuppiMET_pt = PuppiMET_pt;
    MuonSystem->PuppiMET_phi = PuppiMET_phi;

    if (signalScan && !isData && hasEventSignalKey)
      Total2D[eventSignalKey]->Fill(1.0, Generator_weight * MuonSystem->pileupWeight);
    if (signalScan && !isData && hasEventSignalKey) {
      accep2D[eventSignalKey]->Fill(1.0, Generator_weight * MuonSystem->pileupWeight);
    } else if (!isData) {
      if (MuonSystem->gLLP_csc[0] && MuonSystem->gLLP_csc[1])
        accep_csccsc->Fill(1.0, Generator_weight * MuonSystem->pileupWeight);
      if ((MuonSystem->gLLP_dt[0] && MuonSystem->gLLP_csc[1]) || (MuonSystem->gLLP_dt[1] && MuonSystem->gLLP_csc[0]))
        accep_cscdt->Fill(1.0, Generator_weight * MuonSystem->pileupWeight);
    }
    /* #endregion */

    // [8] Event filters, trigger bits, and jet veto maps
    /* #region */
    const EventCutState cuts = buildEventCutState(
        *this,
        helper.get(),
        analysisTag,
        isData,
        run,
        PuppiMET_pt,
        PuppiMET_phi,
        nJet,
        Jet_pt,
        Jet_eta,
        Jet_phi,
        Jet_neEmEF,
        Jet_chEmEF,
        nMuon,
        Muon_isPFcand,
        Muon_eta,
        Muon_phi,
        synth,
        Flag_goodVertices,
        Flag_globalSuperTightHalo2016Filter,
        Flag_EcalDeadCellTriggerPrimitiveFilter,
        Flag_BadPFMuonFilter,
        Flag_BadPFMuonDzFilter,
        Flag_hfNoisyHitsFilter,
        Flag_eeBadScFilter,
        Flag_ecalBadCalibFilter,
        HLT_CscCluster_Loose,
        HLT_L1CSCShower_DTCluster50);
    writeEventCutState(MuonSystem, cuts);
    /* #endregion */

    // [9] Lepton object selection and lepton branch fill
    /* #region */
    //*************************************************************************
    //Start Object Selection
    //*************************************************************************
    std::vector<leptons> Leptons = buildSelectedLeptons(
        *this,
        synth,
        nMuon,
        Muon_looseId,
        Muon_pt,
        Muon_eta,
        Muon_phi,
        Muon_charge,
        Muon_dz,
        Muon_tightId,
        Muon_pfRelIso04_all,
        nElectron,
        Electron_pt,
        Electron_eta,
        Electron_phi,
        Electron_charge,
        Electron_dz);
    fillLeptonBranches(MuonSystem, Leptons);
    /* #endregion */

    // [10] Jet selection, JES propagation, and jet branch fill
    /* #region */
    const JetStageResult jetStage = buildJetStageResult(
        *this,
        helper.get(),
        run,
        nJet,
        Jet_eta,
        Jet_phi,
        Jet_pt,
        synth,
        Leptons);

    fillJetBranches(MuonSystem, jetStage.selectedJets);
    fillPuppiMetJesFromShift(*this, MuonSystem, jetStage);
    /* #endregion */

    // [11] Rechit-level pass-through fill (CSC/DT/RPC + optional aliases)
    /* #region */
    fillRawRechits(
        MuonSystem,
        ncscRechits,
        ndtRecHits,
        nrpcRecHits,
        cscRechits_Quality,
        cscRechits_Chamber,
        cscRechits_Station,
        cscRechits_Eta,
        cscRechits_Phi,
        cscRechits_X,
        cscRechits_Y,
        cscRechits_Z,
        cscRechits_Tpeak,
        cscRechits_Twire,
        cscRechits_IChamber,
        cscRechits_NStrips,
        cscRechits_WGroupsBX,
        cscRechits_HitWire,
        cscRechits_NWireGroups,
        cscRechits_E,
        dtRecHits_Layer,
        dtRecHits_SuperLayer,
        dtRecHits_Station,
        dtRecHits_Wheel,
        dtRecHits_Eta,
        dtRecHits_Phi,
        dtRecHits_X,
        dtRecHits_Y,
        dtRecHits_Z,
        dtRecHits_Sector,
        rpcRecHits_Bx,
        rpcRecHits_Region,
        rpcRecHits_Ring,
        rpcRecHits_Layer,
        rpcRecHits_Station,
        rpcRecHits_Sector,
        rpcRecHits_X,
        rpcRecHits_Y,
        rpcRecHits_Z,
        rpcRecHits_Phi,
        rpcRecHits_Eta,
        rpcRecHits_Time,
        rpcRecHits_TimeError);
    /* #endregion */

    // [12] Ring occupancy summaries from raw rechits
    /* #region */
    MuonSystem->nDtRings = countDtRingsFromRecHits(
        ndtRecHits,
        dtRecHits_Station,
        dtRecHits_Wheel);
    /* #endregion */

    // [13] CSC rechit clustering, feature fill, and matching/veto variables
    /* #region */
    vector<Rechits> points = buildCscRechitPoints(
        ncscRechits,
        cscRechits_Phi,
        cscRechits_Eta,
        cscRechits_X,
        cscRechits_Y,
        cscRechits_Z,
        cscRechits_Tpeak,
        cscRechits_Twire,
        cscRechits_Station,
        cscRechits_Chamber);
    vector<int> cscRechitsClusterId;
    cscRechitsClusterId.reserve(ncscRechits);
    for (int i = 0; i < ncscRechits; i++)
      cscRechitsClusterId.push_back(-1);
    MuonSystem->nCscRings = countCscRingsFromRecHits(
        ncscRechits,
        cscRechits_Chamber);
    // Do CA clustering

    int min_point = 50; //minimum number of Rechitss to call it a cluster
    float epsilon = 0.4; //cluster radius parameter
    CACluster ds(min_point, epsilon, points);
    runCAClustering(ds);
    const auto& cscPointsClustered = ds.points();
    fillClusteredRecHitIds(
        MuonSystem->ncscRechits,
        cscPointsClustered,
        MuonSystem->cscRechits_ClusterId);

    MuonSystem->nCscRechitClusters = 0;
    MuonSystem->nCscRechitClusters_nocut = 0;
    for (auto& tmp : ds.clusters) {
      MuonSystem->nCscRechitClusters_nocut++;
      MuonSystem->cscRechitClusterX[MuonSystem->nCscRechitClusters] = tmp.x;
      MuonSystem->cscRechitClusterY[MuonSystem->nCscRechitClusters] = tmp.y;
      MuonSystem->cscRechitClusterZ[MuonSystem->nCscRechitClusters] = tmp.z;
      MuonSystem->cscRechitClusterTimeWeighted[MuonSystem->nCscRechitClusters] = tmp.tWeighted;
      MuonSystem->cscRechitClusterTimeSpreadWeightedAll[MuonSystem->nCscRechitClusters] = tmp.TSpreadWeightedAll;

      MuonSystem->cscRechitClusterTime[MuonSystem->nCscRechitClusters] = tmp.tTotal;
      MuonSystem->cscRechitClusterTimeSpread[MuonSystem->nCscRechitClusters] = tmp.TSpread;

      MuonSystem->cscRechitClusterEta[MuonSystem->nCscRechitClusters] = tmp.eta;
      MuonSystem->cscRechitClusterPhi[MuonSystem->nCscRechitClusters] = tmp.phi;

      MuonSystem->cscRechitClusterSize[MuonSystem->nCscRechitClusters] = tmp.nhits;
      MuonSystem->cscRechitClusternXY[MuonSystem->nCscRechitClusters] = tmp.nXY;
      MuonSystem->cscRechitClusternZ[MuonSystem->nCscRechitClusters] = tmp.nZ;
      MuonSystem->cscRechitClusterXSpread[MuonSystem->nCscRechitClusters] = tmp.XSpread;
      MuonSystem->cscRechitClusterYSpread[MuonSystem->nCscRechitClusters] = tmp.YSpread;
      MuonSystem->cscRechitClusterZSpread[MuonSystem->nCscRechitClusters] = tmp.ZSpread;
      MuonSystem->cscRechitClusterXYSpread[MuonSystem->nCscRechitClusters] = tmp.XYSpread;
      MuonSystem->cscRechitClusterRSpread[MuonSystem->nCscRechitClusters] = tmp.RSpread;
      MuonSystem->cscRechitClusterEtaPhiSpread[MuonSystem->nCscRechitClusters] = tmp.EtaPhiSpread;
      MuonSystem->cscRechitClusterEtaSpread[MuonSystem->nCscRechitClusters] = tmp.EtaSpread;
      MuonSystem->cscRechitClusterPhiSpread[MuonSystem->nCscRechitClusters] = tmp.PhiSpread;
      MuonSystem->cscRechitClusterDeltaRSpread[MuonSystem->nCscRechitClusters] = tmp.DeltaRSpread;
      MuonSystem->cscRechitClusterMajorAxis[MuonSystem->nCscRechitClusters] = tmp.MajorAxis;
      MuonSystem->cscRechitClusterMinorAxis[MuonSystem->nCscRechitClusters] = tmp.MinorAxis;
      MuonSystem->cscRechitClusterSkewX[MuonSystem->nCscRechitClusters] = tmp.SkewX;
      MuonSystem->cscRechitClusterSkewY[MuonSystem->nCscRechitClusters] = tmp.SkewY;
      MuonSystem->cscRechitClusterSkewZ[MuonSystem->nCscRechitClusters] = tmp.SkewZ;
      MuonSystem->cscRechitClusterKurtX[MuonSystem->nCscRechitClusters] = tmp.KurtX;
      MuonSystem->cscRechitClusterKurtY[MuonSystem->nCscRechitClusters] = tmp.KurtY;
      MuonSystem->cscRechitClusterKurtZ[MuonSystem->nCscRechitClusters] = tmp.KurtZ;

      MuonSystem->cscRechitClusterNRechitChamberPlus11[MuonSystem->nCscRechitClusters] = tmp.nCscRechitsChamberPlus11;
      MuonSystem->cscRechitClusterNRechitChamberPlus12[MuonSystem->nCscRechitClusters] = tmp.nCscRechitsChamberPlus12;
      MuonSystem->cscRechitClusterNRechitChamberPlus13[MuonSystem->nCscRechitClusters] = tmp.nCscRechitsChamberPlus13;
      MuonSystem->cscRechitClusterNRechitChamberPlus21[MuonSystem->nCscRechitClusters] = tmp.nCscRechitsChamberPlus21;
      MuonSystem->cscRechitClusterNRechitChamberPlus22[MuonSystem->nCscRechitClusters] = tmp.nCscRechitsChamberPlus22;
      MuonSystem->cscRechitClusterNRechitChamberPlus31[MuonSystem->nCscRechitClusters] = tmp.nCscRechitsChamberPlus31;
      MuonSystem->cscRechitClusterNRechitChamberPlus32[MuonSystem->nCscRechitClusters] = tmp.nCscRechitsChamberPlus32;
      MuonSystem->cscRechitClusterNRechitChamberPlus41[MuonSystem->nCscRechitClusters] = tmp.nCscRechitsChamberPlus41;
      MuonSystem->cscRechitClusterNRechitChamberPlus42[MuonSystem->nCscRechitClusters] = tmp.nCscRechitsChamberPlus42;
      MuonSystem->cscRechitClusterNRechitChamberMinus11[MuonSystem->nCscRechitClusters] = tmp.nCscRechitsChamberMinus11;
      MuonSystem->cscRechitClusterNRechitChamberMinus12[MuonSystem->nCscRechitClusters] = tmp.nCscRechitsChamberMinus12;
      MuonSystem->cscRechitClusterNRechitChamberMinus13[MuonSystem->nCscRechitClusters] = tmp.nCscRechitsChamberMinus13;
      MuonSystem->cscRechitClusterNRechitChamberMinus21[MuonSystem->nCscRechitClusters] = tmp.nCscRechitsChamberMinus21;
      MuonSystem->cscRechitClusterNRechitChamberMinus22[MuonSystem->nCscRechitClusters] = tmp.nCscRechitsChamberMinus22;
      MuonSystem->cscRechitClusterNRechitChamberMinus31[MuonSystem->nCscRechitClusters] = tmp.nCscRechitsChamberMinus31;
      MuonSystem->cscRechitClusterNRechitChamberMinus32[MuonSystem->nCscRechitClusters] = tmp.nCscRechitsChamberMinus32;
      MuonSystem->cscRechitClusterNRechitChamberMinus41[MuonSystem->nCscRechitClusters] = tmp.nCscRechitsChamberMinus41;
      MuonSystem->cscRechitClusterNRechitChamberMinus42[MuonSystem->nCscRechitClusters] = tmp.nCscRechitsChamberMinus42;
      MuonSystem->cscRechitClusterNRechitME1112[MuonSystem->nCscRechitClusters] = tmp.nCscRechitsChamberMinus11 + tmp.nCscRechitsChamberMinus12 + tmp.nCscRechitsChamberPlus11 + tmp.nCscRechitsChamberPlus12;
      MuonSystem->cscRechitClusterNRechitStation12[MuonSystem->nCscRechitClusters] = MuonSystem->cscRechitClusterNRechitME1112[MuonSystem->nCscRechitClusters] + tmp.nCscRechitsChamberMinus13 + tmp.nCscRechitsChamberPlus13 + tmp.nCscRechitsChamberMinus21 + tmp.nCscRechitsChamberPlus21 + tmp.nCscRechitsChamberMinus22 + tmp.nCscRechitsChamberPlus22;

      float efficiencies[] = {
          helper->getHMTTriggerEff(13, tmp.nCscRechitsChamberPlus13),
          helper->getHMTTriggerEff(21, tmp.nCscRechitsChamberPlus21),
          helper->getHMTTriggerEff(22, tmp.nCscRechitsChamberPlus22),
          helper->getHMTTriggerEff(31, tmp.nCscRechitsChamberPlus31),
          helper->getHMTTriggerEff(32, tmp.nCscRechitsChamberPlus32),
          helper->getHMTTriggerEff(41, tmp.nCscRechitsChamberPlus41),
          helper->getHMTTriggerEff(42, tmp.nCscRechitsChamberPlus42),
          helper->getHMTTriggerEff(13, tmp.nCscRechitsChamberMinus13),
          helper->getHMTTriggerEff(21, tmp.nCscRechitsChamberMinus21),
          helper->getHMTTriggerEff(22, tmp.nCscRechitsChamberMinus22),
          helper->getHMTTriggerEff(31, tmp.nCscRechitsChamberMinus31),
          helper->getHMTTriggerEff(32, tmp.nCscRechitsChamberMinus32),
          helper->getHMTTriggerEff(41, tmp.nCscRechitsChamberMinus41),
          helper->getHMTTriggerEff(42, tmp.nCscRechitsChamberMinus42),
      };

      MuonSystem->cscRechitClusterHMTEfficiency[MuonSystem->nCscRechitClusters] = 1.0;
      for (const float& eff : efficiencies) {
        MuonSystem->cscRechitClusterHMTEfficiency[MuonSystem->nCscRechitClusters] *= (1 - eff);
      }
      MuonSystem->cscRechitClusterHMTEfficiency[MuonSystem->nCscRechitClusters] = 1 - MuonSystem->cscRechitClusterHMTEfficiency[MuonSystem->nCscRechitClusters];

      MuonSystem->cscRechitClusterMaxChamber[MuonSystem->nCscRechitClusters] = tmp.maxChamber;
      MuonSystem->cscRechitClusterMaxChamberRatio[MuonSystem->nCscRechitClusters] = 1.0 * tmp.maxChamberRechits / tmp.nhits;
      MuonSystem->cscRechitClusterNChamber[MuonSystem->nCscRechitClusters] = tmp.nChamber;
      MuonSystem->cscRechitClusterMaxStation[MuonSystem->nCscRechitClusters] = tmp.maxStation;
      MuonSystem->cscRechitClusterMaxStationRatio[MuonSystem->nCscRechitClusters] = 1.0 * tmp.maxStationRechits / tmp.nhits;

      MuonSystem->cscRechitClusterNStation10[MuonSystem->nCscRechitClusters] = tmp.nStation10;
      MuonSystem->cscRechitClusterAvgStation10[MuonSystem->nCscRechitClusters] = tmp.avgStation10;

      //Jet veto/ muon veto
      MuonSystem->cscRechitClusterMuonVetoPt[MuonSystem->nCscRechitClusters] = 0.0;
      MuonSystem->cscRechitClusterMuonVetoE[MuonSystem->nCscRechitClusters] = 0.0;

      fillClusterJetVeto(
          *this,
          helper.get(),
          run,
          nJet,
          Jet_eta,
          Jet_phi,
          Jet_pt,
          synth.jetE.data(),
          synth.jetPassIDTight.data(),
          MuonSystem->cscRechitClusterEta[MuonSystem->nCscRechitClusters],
          MuonSystem->cscRechitClusterPhi[MuonSystem->nCscRechitClusters],
          MuonSystem->cscRechitClusterJetVetoPt[MuonSystem->nCscRechitClusters],
          MuonSystem->cscRechitClusterJetVetoE[MuonSystem->nCscRechitClusters],
          MuonSystem->cscRechitClusterJetVetoTightId[MuonSystem->nCscRechitClusters],
          MuonSystem->cscRechitClusterJetVetoPtJESUp[MuonSystem->nCscRechitClusters],
          MuonSystem->cscRechitClusterJetVetoPtJESDown[MuonSystem->nCscRechitClusters]);
      for (int i = 0; i < nMuon; i++) {
        if (fabs(Muon_eta[i]) > 3.0)
          continue;
        if (RazorAnalyzerMerged::deltaR(Muon_eta[i], Muon_phi[i], MuonSystem->cscRechitClusterEta[MuonSystem->nCscRechitClusters], MuonSystem->cscRechitClusterPhi[MuonSystem->nCscRechitClusters]) < 0.4 && Muon_pt[i] > MuonSystem->cscRechitClusterMuonVetoPt[MuonSystem->nCscRechitClusters]) {
          MuonSystem->cscRechitClusterMuonVetoPt[MuonSystem->nCscRechitClusters] = Muon_pt[i];
          MuonSystem->cscRechitClusterMuonVetoE[MuonSystem->nCscRechitClusters] = synth.muonE[i];
          MuonSystem->cscRechitClusterMuonVetoGlobal[MuonSystem->nCscRechitClusters] = Muon_isGlobal[i];
          MuonSystem->cscRechitClusterMuonVetoLooseId[MuonSystem->nCscRechitClusters] = Muon_looseId[i];
        }
        if (RazorAnalyzerMerged::deltaR(Muon_eta[i], Muon_phi[i], MuonSystem->cscRechitClusterEta[MuonSystem->nCscRechitClusters], MuonSystem->cscRechitClusterPhi[MuonSystem->nCscRechitClusters]) < 0.8 && Muon_pt[i] > MuonSystem->cscRechitClusterMuonVetoPt0p8Thresh[MuonSystem->nCscRechitClusters]) {
          MuonSystem->cscRechitClusterMuonVetoPt0p8Thresh[MuonSystem->nCscRechitClusters] = Muon_pt[i];
          MuonSystem->cscRechitClusterMuonVetoE0p8Thresh[MuonSystem->nCscRechitClusters] = synth.muonE[i];
          MuonSystem->cscRechitClusterMuonVetoGlobal0p8Thresh[MuonSystem->nCscRechitClusters] = Muon_isGlobal[i];
          MuonSystem->cscRechitClusterMuonVetoLooseId0p8Thresh[MuonSystem->nCscRechitClusters] = Muon_looseId[i];
        }
      }
      if (!isData) {
        // match to gen level LLP
        const GLLPMatchResult match = findNearestGLLPMatch(
            *this,
            MuonSystem->cscRechitClusterEta[MuonSystem->nCscRechitClusters],
            MuonSystem->cscRechitClusterPhi[MuonSystem->nCscRechitClusters],
            MuonSystem->nGLLP,
            MuonSystem->gLLP_eta,
            MuonSystem->gLLP_phi);

        MuonSystem->cscRechitCluster_match_gLLP[MuonSystem->nCscRechitClusters] = (match.minDeltaR < 0.4f);
        MuonSystem->cscRechitCluster_match_gLLP_minDeltaR[MuonSystem->nCscRechitClusters] = match.minDeltaR;
        MuonSystem->cscRechitCluster_match_gLLP_index[MuonSystem->nCscRechitClusters] = match.index;
        if (match.index >= 0 && match.index < MuonSystem->nGLLP) {
          fillMatchedGLLPFields(
              MuonSystem,
              MuonSystem->nCscRechitClusters,
              match.index,
              MuonSystem->cscRechitCluster_match_gLLP_eta,
              MuonSystem->cscRechitCluster_match_gLLP_phi,
              MuonSystem->cscRechitCluster_match_gLLP_decay_r,
              MuonSystem->cscRechitCluster_match_gLLP_decay_z,
              MuonSystem->cscRechitCluster_match_gLLP_csc,
              MuonSystem->cscRechitCluster_match_gLLP_dt,
              MuonSystem->cscRechitCluster_match_gLLP_e);
        }
      }

      //match to MB1 DT segments
      for (int i = 0; i < ndtSegments; i++) {
        if (RazorAnalyzerMerged::deltaR(dtSegments_Eta[i], dtSegments_Phi[i], MuonSystem->cscRechitClusterEta[MuonSystem->nCscRechitClusters], MuonSystem->cscRechitClusterPhi[MuonSystem->nCscRechitClusters]) < 0.4) {
          MuonSystem->cscRechitCluster_match_dtSeg_0p4[MuonSystem->nCscRechitClusters]++;
          if (dtSegments_Station[i] == 1)
            MuonSystem->cscRechitCluster_match_MB1Seg_0p4[MuonSystem->nCscRechitClusters]++;
        }
      }

      //match to RPC hits in RE1/2
      for (int i = 0; i < nrpcRecHits; i++) {
        float rpcR = sqrt(rpcRecHits_X[i] * rpcRecHits_X[i] + rpcRecHits_Y[i] * rpcRecHits_Y[i]);
        if (RazorAnalyzerMerged::deltaR(rpcRecHits_Eta[i], rpcRecHits_Phi[i], MuonSystem->cscRechitClusterEta[MuonSystem->nCscRechitClusters], MuonSystem->cscRechitClusterPhi[MuonSystem->nCscRechitClusters]) < 0.4) {
          if (rpcR < 461.0 && rpcR > 275 && abs(rpcRecHits_Z[i]) > 663 && abs(rpcRecHits_Z[i]) < 730)
            MuonSystem->cscRechitCluster_match_RE12_0p4[MuonSystem->nCscRechitClusters]++;
          if (rpcR < 470 && rpcR > 380 && abs(rpcRecHits_Z[i]) < 661)
            MuonSystem->cscRechitCluster_match_RB1_0p4[MuonSystem->nCscRechitClusters]++;
        }
      }

      MuonSystem->cscRechitClusterMet_dPhi[MuonSystem->nCscRechitClusters] = RazorAnalyzerMerged::deltaPhi(MuonSystem->cscRechitClusterPhi[MuonSystem->nCscRechitClusters], MuonSystem->PFMET_phi);
      MuonSystem->cscRechitClusterPuppiMet_dPhi[MuonSystem->nCscRechitClusters] = RazorAnalyzerMerged::deltaPhi(MuonSystem->cscRechitClusterPhi[MuonSystem->nCscRechitClusters], MuonSystem->PuppiMET_phi);
      if (MuonSystem->nTaus == 1) {
        MuonSystem->cscRechitClusterPromptTauDeltaEta[MuonSystem->nCscRechitClusters] = MuonSystem->tauEta[0] - MuonSystem->cscRechitClusterEta[MuonSystem->nCscRechitClusters];
        MuonSystem->cscRechitClusterPromptTauDeltaPhi[MuonSystem->nCscRechitClusters] = RazorAnalyzerMerged::deltaPhi(MuonSystem->tauPhi[0], MuonSystem->cscRechitClusterPhi[MuonSystem->nCscRechitClusters]);
        MuonSystem->cscRechitClusterPromptTauDeltaR[MuonSystem->nCscRechitClusters] = RazorAnalyzerMerged::deltaR(MuonSystem->tauEta[0], MuonSystem->tauPhi[0], MuonSystem->cscRechitClusterEta[MuonSystem->nCscRechitClusters], MuonSystem->cscRechitClusterPhi[MuonSystem->nCscRechitClusters]);
      }
      if (MuonSystem->nLeptons > 0) {
        if (abs(MuonSystem->lepPdgId[0]) == 11) {
          MuonSystem->cscRechitClusterPromptEleDeltaEta[MuonSystem->nCscRechitClusters] = MuonSystem->lepEta[0] - MuonSystem->cscRechitClusterEta[MuonSystem->nCscRechitClusters];
          MuonSystem->cscRechitClusterPromptEleDeltaPhi[MuonSystem->nCscRechitClusters] = RazorAnalyzerMerged::deltaPhi(MuonSystem->lepPhi[0], MuonSystem->cscRechitClusterPhi[MuonSystem->nCscRechitClusters]);
          MuonSystem->cscRechitClusterPromptEleDeltaR[MuonSystem->nCscRechitClusters] = RazorAnalyzerMerged::deltaR(MuonSystem->lepEta[0], MuonSystem->lepPhi[0], MuonSystem->cscRechitClusterEta[MuonSystem->nCscRechitClusters], MuonSystem->cscRechitClusterPhi[MuonSystem->nCscRechitClusters]);
        }
        if (abs(MuonSystem->lepPdgId[0]) == 13) {
          MuonSystem->cscRechitClusterPromptMuDeltaEta[MuonSystem->nCscRechitClusters] = MuonSystem->lepEta[0] - MuonSystem->cscRechitClusterEta[MuonSystem->nCscRechitClusters];
          MuonSystem->cscRechitClusterPromptMuDeltaPhi[MuonSystem->nCscRechitClusters] = RazorAnalyzerMerged::deltaPhi(MuonSystem->lepPhi[0], MuonSystem->cscRechitClusterPhi[MuonSystem->nCscRechitClusters]);
          MuonSystem->cscRechitClusterPromptMuDeltaR[MuonSystem->nCscRechitClusters] = RazorAnalyzerMerged::deltaR(MuonSystem->lepEta[0], MuonSystem->lepPhi[0], MuonSystem->cscRechitClusterEta[MuonSystem->nCscRechitClusters], MuonSystem->cscRechitClusterPhi[MuonSystem->nCscRechitClusters]);
        }
      }

      MuonSystem->nCscRechitClusters++;
    }
    /* #endregion */

    // [14] DT/RPC clustering, DT feature fill, and matching/veto variables
    /* #region */
    // DT cluster

    points = buildDtRecHitPoints(
        ndtRecHits,
        dtRecHits_Phi,
        dtRecHits_Eta,
        dtRecHits_X,
        dtRecHits_Y,
        dtRecHits_Z,
        dtRecHits_Station,
        dtRecHits_Wheel,
        dtRecHits_SuperLayer);

    // Do CA clustering
    int min_point_dt = 50; //minimum number of segments to call it a cluster
    float epsilon_dt = 0.2; //cluster radius parameter
    CACluster ds_dtRechit(min_point_dt, epsilon_dt, points);
    runCAClustering(ds_dtRechit);
    const auto& dtPointsClustered = ds_dtRechit.points();
    fillClusteredRecHitIds(
        MuonSystem->ndtRecHits,
        dtPointsClustered,
        MuonSystem->dtRecHits_ClusterId);

    // RPC cluster IDs
    points = buildRpcRecHitPoints(
        MuonSystem->nrpcRecHits,
        rpcRecHits_Phi,
        rpcRecHits_Eta,
        rpcRecHits_X,
        rpcRecHits_Y,
        rpcRecHits_Z,
        rpcRecHits_Time,
        rpcRecHits_Station,
        rpcRecHits_Sector,
        rpcRecHits_Layer,
        rpcRecHits_Ring);

    int min_point_rpc = 5;
    float epsilon_rpc = 0.2;
    CACluster ds_rpcRechit(min_point_rpc, epsilon_rpc, points);
    runCAClustering(ds_rpcRechit);
    const auto& rpcPointsClustered = ds_rpcRechit.points();
    fillClusteredRecHitIds(
        MuonSystem->nrpcRecHits,
        rpcPointsClustered,
        MuonSystem->rpcRecHits_ClusterId);

    MuonSystem->nDtRechitClusters = 0;
    MuonSystem->nDtRechitClusters_nocut = 0;
    for (auto& tmp : ds_dtRechit.clusters) {
      //remove overlaps
      bool overlap = false;
      for (int i = 0; i < MuonSystem->nCscRechitClusters; i++) {
        if (RazorAnalyzerMerged::deltaR(MuonSystem->cscRechitClusterEta[i], MuonSystem->cscRechitClusterPhi[i], tmp.eta, tmp.phi) < 0.4)
          overlap = true;
      }
      if (overlap)
        continue;
      MuonSystem->nDtRechitClusters_nocut++;
      MuonSystem->dtRechitClusterX[MuonSystem->nDtRechitClusters] = tmp.x;
      MuonSystem->dtRechitClusterY[MuonSystem->nDtRechitClusters] = tmp.y;
      MuonSystem->dtRechitClusterZ[MuonSystem->nDtRechitClusters] = tmp.z;
      if (abs(tmp.z) < 126.8)
        MuonSystem->dtRechitClusterWheel[MuonSystem->nDtRechitClusters] = 0;
      else if (tmp.z > 126.8 && tmp.z < 395.4)
        MuonSystem->dtRechitClusterWheel[MuonSystem->nDtRechitClusters] = 1;
      else if (tmp.z < -126.8 && tmp.z > -395.4)
        MuonSystem->dtRechitClusterWheel[MuonSystem->nDtRechitClusters] = -1;
      else if (tmp.z < 0)
        MuonSystem->dtRechitClusterWheel[MuonSystem->nDtRechitClusters] = -2;
      else
        MuonSystem->dtRechitClusterWheel[MuonSystem->nDtRechitClusters] = 2;
      MuonSystem->dtRechitClusterEta[MuonSystem->nDtRechitClusters] = tmp.eta;
      MuonSystem->dtRechitClusterPhi[MuonSystem->nDtRechitClusters] = tmp.phi;
      MuonSystem->dtRechitClusterSize[MuonSystem->nDtRechitClusters] = tmp.nhits;
      MuonSystem->dtRechitClusternXY[MuonSystem->nDtRechitClusters] = tmp.nXY;
      MuonSystem->dtRechitClusternZ[MuonSystem->nDtRechitClusters] = tmp.nZ;
      MuonSystem->dtRechitClusterXSpread[MuonSystem->nDtRechitClusters] = tmp.XSpread;
      MuonSystem->dtRechitClusterYSpread[MuonSystem->nDtRechitClusters] = tmp.YSpread;
      MuonSystem->dtRechitClusterZSpread[MuonSystem->nDtRechitClusters] = tmp.ZSpread;
      MuonSystem->dtRechitClusterXYSpread[MuonSystem->nDtRechitClusters] = tmp.XYSpread;
      MuonSystem->dtRechitClusterRSpread[MuonSystem->nDtRechitClusters] = tmp.RSpread;
      MuonSystem->dtRechitClusterEtaPhiSpread[MuonSystem->nDtRechitClusters] = tmp.EtaPhiSpread;
      MuonSystem->dtRechitClusterEtaSpread[MuonSystem->nDtRechitClusters] = tmp.EtaSpread;
      MuonSystem->dtRechitClusterPhiSpread[MuonSystem->nDtRechitClusters] = tmp.PhiSpread;
      MuonSystem->dtRechitClusterDeltaRSpread[MuonSystem->nDtRechitClusters] = tmp.DeltaRSpread;
      MuonSystem->dtRechitClusterMajorAxis[MuonSystem->nDtRechitClusters] = tmp.MajorAxis;
      MuonSystem->dtRechitClusterMinorAxis[MuonSystem->nDtRechitClusters] = tmp.MinorAxis;
      MuonSystem->dtRechitClusterSkewX[MuonSystem->nDtRechitClusters] = tmp.SkewX;
      MuonSystem->dtRechitClusterSkewY[MuonSystem->nDtRechitClusters] = tmp.SkewY;
      MuonSystem->dtRechitClusterSkewZ[MuonSystem->nDtRechitClusters] = tmp.SkewZ;
      MuonSystem->dtRechitClusterKurtX[MuonSystem->nDtRechitClusters] = tmp.KurtX;
      MuonSystem->dtRechitClusterKurtY[MuonSystem->nDtRechitClusters] = tmp.KurtY;
      MuonSystem->dtRechitClusterKurtZ[MuonSystem->nDtRechitClusters] = tmp.KurtZ;

      unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
      default_random_engine generator(seed);

      // default_random_engine generator;
      uniform_real_distribution<double> distribution(0.0, 1.0);
      float prob = 0.03;
      for (int i = 0; i < 12; ++i) {
        if (distribution(generator) < prob)
          MuonSystem->dtRechitClusterNoiseHitStation1[MuonSystem->nDtRechitClusters]++;
      }
      for (int i = 0; i < 12; ++i) {
        if (distribution(generator) < prob)
          MuonSystem->dtRechitClusterNoiseHitStation2[MuonSystem->nDtRechitClusters]++;
      }
      for (int i = 0; i < 12; ++i) {
        if (distribution(generator) < prob)
          MuonSystem->dtRechitClusterNoiseHitStation3[MuonSystem->nDtRechitClusters]++;
      }
      for (int i = 0; i < 8; ++i) {
        if (distribution(generator) < prob)
          MuonSystem->dtRechitClusterNoiseHitStation4[MuonSystem->nDtRechitClusters]++;
      }

      MuonSystem->dtRechitClusterNoiseHit[MuonSystem->nDtRechitClusters] = MuonSystem->dtRechitClusterNoiseHitStation1[MuonSystem->nDtRechitClusters] +
                                                                           MuonSystem->dtRechitClusterNoiseHitStation2[MuonSystem->nDtRechitClusters] +
                                                                           MuonSystem->dtRechitClusterNoiseHitStation3[MuonSystem->nDtRechitClusters] +
                                                                           MuonSystem->dtRechitClusterNoiseHitStation4[MuonSystem->nDtRechitClusters];

      MuonSystem->dtRechitClusterNHitStation1[MuonSystem->nDtRechitClusters] = tmp.nDtRechitsStation1;
      MuonSystem->dtRechitClusterNHitStation2[MuonSystem->nDtRechitClusters] = tmp.nDtRechitsStation2;
      MuonSystem->dtRechitClusterNHitStation3[MuonSystem->nDtRechitClusters] = tmp.nDtRechitsStation3;
      MuonSystem->dtRechitClusterNHitStation4[MuonSystem->nDtRechitClusters] = tmp.nDtRechitsStation4;

      MuonSystem->dtRechitClusterNHitWheel0[MuonSystem->nDtRechitClusters] = tmp.nDtRechitsWheel0;
      MuonSystem->dtRechitClusterNHitWheel1[MuonSystem->nDtRechitClusters] = tmp.nDtRechitsWheel1;
      MuonSystem->dtRechitClusterNHitWheel2[MuonSystem->nDtRechitClusters] = tmp.nDtRechitsWheel2;

      MuonSystem->dtRechitClusterMaxChamber[MuonSystem->nDtRechitClusters] = tmp.maxChamber;
      MuonSystem->dtRechitClusterMaxChamberRatio[MuonSystem->nDtRechitClusters] = 1.0 * tmp.maxChamberRechits / tmp.nhits;
      MuonSystem->dtRechitClusterNChamber[MuonSystem->nDtRechitClusters] = tmp.nChamber;
      MuonSystem->dtRechitClusterMaxStation[MuonSystem->nDtRechitClusters] = tmp.maxStation;
      MuonSystem->dtRechitClusterMaxStationRatio[MuonSystem->nDtRechitClusters] = 1.0 * tmp.maxStationRechits / tmp.nhits;
      MuonSystem->dtRechitClusterNStation10[MuonSystem->nDtRechitClusters] = tmp.nStation10;
      MuonSystem->dtRechitClusterAvgStation10[MuonSystem->nDtRechitClusters] = tmp.avgStation10;

      //Jet veto/ muon veto
      MuonSystem->dtRechitClusterMuonVetoPt[MuonSystem->nDtRechitClusters] = 0.0;
      MuonSystem->dtRechitClusterMuonVetoE[MuonSystem->nDtRechitClusters] = 0.0;

      fillClusterJetVeto(
          *this,
          helper.get(),
          run,
          nJet,
          Jet_eta,
          Jet_phi,
          Jet_pt,
          synth.jetE.data(),
          synth.jetPassIDTight.data(),
          MuonSystem->dtRechitClusterEta[MuonSystem->nDtRechitClusters],
          MuonSystem->dtRechitClusterPhi[MuonSystem->nDtRechitClusters],
          MuonSystem->dtRechitClusterJetVetoPt[MuonSystem->nDtRechitClusters],
          MuonSystem->dtRechitClusterJetVetoE[MuonSystem->nDtRechitClusters],
          MuonSystem->dtRechitClusterJetVetoTightId[MuonSystem->nDtRechitClusters],
          MuonSystem->dtRechitClusterJetVetoPtJESUp[MuonSystem->nDtRechitClusters],
          MuonSystem->dtRechitClusterJetVetoPtJESDown[MuonSystem->nDtRechitClusters]);

      for (int i = 0; i < nMuon; i++) {
        if (fabs(Muon_eta[i]) > 3.0)
          continue;
        if (RazorAnalyzerMerged::deltaR(Muon_eta[i], Muon_phi[i], MuonSystem->dtRechitClusterEta[MuonSystem->nDtRechitClusters], MuonSystem->dtRechitClusterPhi[MuonSystem->nDtRechitClusters]) < 0.4 && Muon_pt[i] > MuonSystem->dtRechitClusterMuonVetoPt[MuonSystem->nDtRechitClusters]) {
          MuonSystem->dtRechitClusterMuonVetoPt[MuonSystem->nDtRechitClusters] = Muon_pt[i];
          MuonSystem->dtRechitClusterMuonVetoE[MuonSystem->nDtRechitClusters] = synth.muonE[i];
          MuonSystem->dtRechitClusterMuonVetoGlobal[MuonSystem->nDtRechitClusters] = Muon_isGlobal[i];
          MuonSystem->dtRechitClusterMuonVetoLooseId[MuonSystem->nDtRechitClusters] = Muon_looseId[i];
        }
      }

      for (int i = 0; i < ndtSegments; i++) {
        if (RazorAnalyzerMerged::deltaR(dtSegments_Eta[i], dtSegments_Phi[i], MuonSystem->dtRechitClusterEta[MuonSystem->nDtRechitClusters], MuonSystem->dtRechitClusterPhi[MuonSystem->nDtRechitClusters]) < 0.4) {
          if (dtSegments_Station[i] == 1)
            MuonSystem->dtRechitClusterNSegStation1[MuonSystem->nDtRechitClusters] += 1;
          if (dtSegments_Station[i] == 2)
            MuonSystem->dtRechitClusterNSegStation2[MuonSystem->nDtRechitClusters] += 1;
          if (dtSegments_Station[i] == 3)
            MuonSystem->dtRechitClusterNSegStation3[MuonSystem->nDtRechitClusters] += 1;
          if (dtSegments_Station[i] == 4)
            MuonSystem->dtRechitClusterNSegStation4[MuonSystem->nDtRechitClusters] += 1;
        }
        if (abs(RazorAnalyzerMerged::deltaPhi(dtSegments_Phi[i], MuonSystem->dtRechitClusterPhi[MuonSystem->nDtRechitClusters])) > 2) {
          if (dtSegments_Station[i] == 1)
            MuonSystem->dtRechitClusterNOppositeSegStation1[MuonSystem->nDtRechitClusters] += 1;
          if (dtSegments_Station[i] == 2)
            MuonSystem->dtRechitClusterNOppositeSegStation2[MuonSystem->nDtRechitClusters] += 1;
          if (dtSegments_Station[i] == 3)
            MuonSystem->dtRechitClusterNOppositeSegStation3[MuonSystem->nDtRechitClusters] += 1;
          if (dtSegments_Station[i] == 4)
            MuonSystem->dtRechitClusterNOppositeSegStation4[MuonSystem->nDtRechitClusters] += 1;
        }
      }

      // match to gen-level LLP
      if (!isData) {
        const GLLPMatchResult match = findNearestGLLPMatch(
            *this,
            tmp.eta,
            tmp.phi,
            MuonSystem->nGLLP,
            MuonSystem->gLLP_eta,
            MuonSystem->gLLP_phi);

        MuonSystem->dtRechitCluster_match_gLLP[MuonSystem->nDtRechitClusters] = (match.minDeltaR < 0.4f);
        MuonSystem->dtRechitCluster_match_gLLP_minDeltaR[MuonSystem->nDtRechitClusters] = match.minDeltaR;
        MuonSystem->dtRechitCluster_match_gLLP_index[MuonSystem->nDtRechitClusters] = match.index;
        if (match.index >= 0 && match.index < MuonSystem->nGLLP) {
          fillMatchedGLLPFields(
              MuonSystem,
              MuonSystem->nDtRechitClusters,
              match.index,
              MuonSystem->dtRechitCluster_match_gLLP_eta,
              MuonSystem->dtRechitCluster_match_gLLP_phi,
              MuonSystem->dtRechitCluster_match_gLLP_decay_r,
              MuonSystem->dtRechitCluster_match_gLLP_decay_z,
              MuonSystem->dtRechitCluster_match_gLLP_csc,
              MuonSystem->dtRechitCluster_match_gLLP_dt,
              MuonSystem->dtRechitCluster_match_gLLP_e);
        }
      }

      //match to MB1 DT segments
      for (int i = 0; i < ndtRecHits; i++) {
        if (RazorAnalyzerMerged::deltaR(dtRecHits_Eta[i], dtRecHits_Phi[i], MuonSystem->dtRechitClusterEta[MuonSystem->nDtRechitClusters], MuonSystem->dtRechitClusterPhi[MuonSystem->nDtRechitClusters]) < 0.5) {
          if (dtRecHits_Station[i] == 1)
            MuonSystem->dtRechitCluster_match_MB1hits_0p5[MuonSystem->nDtRechitClusters]++;
        }
        if (RazorAnalyzerMerged::deltaR(dtRecHits_Eta[i], dtRecHits_Phi[i], MuonSystem->dtRechitClusterEta[MuonSystem->nDtRechitClusters], MuonSystem->dtRechitClusterPhi[MuonSystem->nDtRechitClusters]) < 0.4) {
          if (dtRecHits_Station[i] == 1)
            MuonSystem->dtRechitCluster_match_MB1hits_0p4[MuonSystem->nDtRechitClusters]++;
        }
        if (abs(dtRecHits_Wheel[i] - MuonSystem->dtRechitClusterWheel[MuonSystem->nDtRechitClusters]) == 1 && dtRecHits_Station[i] == 1) {
          if (abs(RazorAnalyzerMerged::deltaPhi(dtRecHits_Phi[i], MuonSystem->dtRechitClusterPhi[MuonSystem->nDtRechitClusters])) < TMath::Pi() / 4.0) {
            if (dtRecHits_Wheel[i] - MuonSystem->dtRechitClusterWheel[MuonSystem->nDtRechitClusters] == 1)
              MuonSystem->dtRechitCluster_match_MB1hits_cosmics_plus[MuonSystem->nDtRechitClusters]++;
            else
              MuonSystem->dtRechitCluster_match_MB1hits_cosmics_minus[MuonSystem->nDtRechitClusters]++;
          }
        }
      }

      std::vector<int> dtRechitCluster_match_rpcBx;

      //match to RPC hits with dPhi<0.5 and same wheel in DT
      for (int i = 0; i < nrpcRecHits; i++) {
        float rpcR = sqrt(rpcRecHits_X[i] * rpcRecHits_X[i] + rpcRecHits_Y[i] * rpcRecHits_Y[i]);
        if (rpcRecHits_Region[i] != 0)
          continue;
        if (abs(RazorAnalyzerMerged::deltaPhi(rpcRecHits_Phi[i], MuonSystem->dtRechitClusterPhi[MuonSystem->nDtRechitClusters])) < 0.5) {
          if (rpcRecHits_Ring[i] == MuonSystem->dtRechitClusterWheel[MuonSystem->nDtRechitClusters]) {
            dtRechitCluster_match_rpcBx.push_back(rpcRecHits_Bx[i]);
            MuonSystem->dtRechitCluster_match_RPChits_dPhi0p5[MuonSystem->nDtRechitClusters]++;
            if (rpcR < 470 && rpcR > 380 && abs(rpcRecHits_Z[i]) < 661)
              MuonSystem->dtRechitCluster_match_RB1_dPhi0p5[MuonSystem->nDtRechitClusters]++;
          }
        }
        if (RazorAnalyzerMerged::deltaR(rpcRecHits_Eta[i], rpcRecHits_Phi[i], MuonSystem->dtRechitClusterEta[MuonSystem->nDtRechitClusters], MuonSystem->dtRechitClusterPhi[MuonSystem->nDtRechitClusters]) < 0.4) {
          if (rpcR < 470 && rpcR > 380 && abs(rpcRecHits_Z[i]) < 661)
            MuonSystem->dtRechitCluster_match_RB1_0p4[MuonSystem->nDtRechitClusters]++;
        }
      }
      int max_occurence = 0;
      int max_bx = -999;
      for (unsigned int l = 0; l < dtRechitCluster_match_rpcBx.size(); l++) {
        int counter = 0;
        for (unsigned int j = 0; j < dtRechitCluster_match_rpcBx.size(); j++) {
          if (dtRechitCluster_match_rpcBx[j] == dtRechitCluster_match_rpcBx[l])
            counter++;
        }
        if (counter > max_occurence) {
          max_occurence = counter;
          max_bx = dtRechitCluster_match_rpcBx[l];
        }
      }
      MuonSystem->dtRechitCluster_match_RPCBx_dPhi0p5[MuonSystem->nDtRechitClusters] = max_bx;

      MuonSystem->dtRechitClusterMet_dPhi[MuonSystem->nDtRechitClusters] = RazorAnalyzerMerged::deltaPhi(MuonSystem->dtRechitClusterPhi[MuonSystem->nDtRechitClusters], MuonSystem->PFMET_phi);
      MuonSystem->dtRechitClusterPuppiMet_dPhi[MuonSystem->nDtRechitClusters] = RazorAnalyzerMerged::deltaPhi(MuonSystem->dtRechitClusterPhi[MuonSystem->nDtRechitClusters], MuonSystem->PuppiMET_phi);

      MuonSystem->nDtRechitClusters++;
    }
    /* #endregion */

    // [15] Final per-event tree fill
    /* #region */
    if (!isData && signalScan) {
      if (!hasEventSignalKey) {
        std::cerr << "[ERROR]: signalScan event is missing a parsed signal-point key. Exiting." << std::endl;
        return;
      }
      Trees2D[eventSignalKey]->Fill();
    } else {
      MuonSystem->tree_->Fill();
    }
    /* #endregion */
  }
  /* #endregion */

  // [16] End-of-job writeout for trees and histograms
  /* #region */
  if (!isData && signalScan) {
    for (auto& filePtr : Files2D) {
      cout << "Writing output tree (" << filePtr.second->GetName() << ")" << endl;
      filePtr.second->cd();
      Trees2D[filePtr.first]->Write();
      NEvents2D[filePtr.first]->Write("NEvents");
      Total2D[filePtr.first]->Write("Total");
      accep2D[filePtr.first]->Write("accep");
      accep_met2D[filePtr.first]->Write("accep_met");
      filePtr.second->Close();
    }
  } else if (!isData) {
    cout << "Filled Total of " << NEvents->GetBinContent(1) << " Events\n";
    cout << "Writing output trees..." << endl;
    outFile->cd();
    MuonSystem->tree_->Write();
    NEvents->Write();
    accep->Write("accep");
    accep_csccsc->Write("accep_csccsc");
    accep_cscdt->Write("accep_cscdt");
    accep_met->Write("accep_met");
    outFile->Close();
  } else {
    cout << "Filled Total of " << NEvents->GetBinContent(1) << " Events\n";
    cout << "Writing output trees..." << endl;
    outFile->cd();
    MuonSystem->tree_->Write();
    Nmet200->Write();
    NmetFilter->Write();
    Nlep0->Write();
    Njet1->Write();
    NcosmicVeto->Write();
    NEvents->Write();
    // outFile->Write();
    outFile->Close();
  }
  /* #endregion */
}
          
