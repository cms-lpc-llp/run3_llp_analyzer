#ifndef MUONSYSTEM_PHASE_TYPES_H
#define MUONSYSTEM_PHASE_TYPES_H

/* #region: includes */
#include <array>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "TH1F.h"
#include "TLorentzVector.h"
#include "TreeMuonSystemCombination.h"
/* #endregion */

/* #region: shared constants */
constexpr float ELE_MASS = 0.000511f;
constexpr float MU_MASS = 0.105658f;
constexpr int INDEX_DEFAULT = -1;
constexpr int kUnclassifiedClusterId = -1;
/* #endregion */

/* #region: phase data types */
struct LeptonCandidate {
  TLorentzVector lepton;
  #define NTUPLE_RECO_LEPTON_TO_STRUCT(TYPE, BRANCH, STRUCT_FIELD, LEAF, DEFAULT) \
    TYPE STRUCT_FIELD{};
  NTUPLE_RECO_LEPTON_FIELD_TABLE(NTUPLE_RECO_LEPTON_TO_STRUCT)
  #undef NTUPLE_RECO_LEPTON_TO_STRUCT
  bool passVetoId = false;
};

struct JetCandidate {
  TLorentzVector jet;
  float time = 0.0f;
  bool passId = false;
  bool isCSVL = false;
  int ecalNRechits = 0;
  float ecalRechitE = 0.0f;
  float jetChargedEMEnergyFraction = 0.0f;
  float jetNeutralEMEnergyFraction = 0.0f;
  float jetChargedHadronEnergyFraction = 0.0f;
  float jetNeutralHadronEnergyFraction = 0.0f;
  bool jetPassMuFrac = false;
  float jetPtJESUp = 0.0f;
  float jetPtJESDown = 0.0f;
  float jetEJESUp = 0.0f;
  float jetEJESDown = 0.0f;
  float JecUnc = 0.0f;

  float electronEnergyFraction = 0.0f;
  float neutralEmEnergyFraction = 0.0f;
  float chargedHadronEnergyFraction = 0.0f;
  float neutralHadronEnergyFraction = 0.0f;
  float muonEnergyFraction = 0.0f;
};

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
  std::vector<JetCandidate> selectedJets;
  float metXJesDown = 0.0f;
  float metYJesDown = 0.0f;
  float metXJesUp = 0.0f;
  float metYJesUp = 0.0f;
};

struct GLLPMatchResult {
  float minDeltaR = 15.0f;
  int index = INDEX_DEFAULT;
};

using SignalPointKey = std::pair<int, std::string>;
/* #endregion */

/* #region: utility declarations */
std::unique_ptr<TH1F> makeCounterHist(const std::string& name);
/* #endregion */

#endif
