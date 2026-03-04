/* #region: includes */
#include "MuonSystemCuttingPhase.h"

#include <cmath>
/* #endregion */

/* #region: event cutting builders */
EventCutState buildEventCutState(
    RazorAnalyzerMerged& analyzer,
    RazorHelper* helper,
    const std::string& analysisTag,
    bool isData,
    int runNumber,
    const EventSynthesis& synth) {
  return buildEventCutState(
      analyzer,
      helper,
      analysisTag,
      isData,
      runNumber,
      analyzer.PuppiMET_pt,
      analyzer.PuppiMET_phi,
      analyzer.nJet,
      analyzer.Jet_pt,
      analyzer.Jet_eta,
      analyzer.Jet_phi,
      analyzer.Jet_neEmEF,
      analyzer.Jet_chEmEF,
      analyzer.nMuon,
      analyzer.Muon_isPFcand,
      analyzer.Muon_eta,
      analyzer.Muon_phi,
      synth,
      analyzer.Flag_goodVertices,
      analyzer.Flag_globalSuperTightHalo2016Filter,
      analyzer.Flag_EcalDeadCellTriggerPrimitiveFilter,
      analyzer.Flag_BadPFMuonFilter,
      analyzer.Flag_BadPFMuonDzFilter,
      analyzer.Flag_hfNoisyHitsFilter,
      analyzer.Flag_eeBadScFilter,
      analyzer.Flag_ecalBadCalibFilter,
      analyzer.HLT_CscCluster_Loose,
      analyzer.HLT_L1CSCShower_DTCluster50);
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
/* #endregion */

/* #region: cluster-level helper utilities */
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

int dtWheelFromClusterZ(float z) {
  if (std::abs(z) < 126.8f)
    return 0;
  if (z > 126.8f && z < 395.4f)
    return 1;
  if (z < -126.8f && z > -395.4f)
    return -1;
  if (z < 0.0f)
    return -2;
  return 2;
}
/* #endregion */
