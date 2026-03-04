#ifndef MUONSYSTEM_CUTTING_PHASE_H
#define MUONSYSTEM_CUTTING_PHASE_H

/* #region: includes */
#include "MuonSystemPhaseTypes.h"
#include "RazorAnalyzer.h"
#include "RazorHelper.h"
/* #endregion */

/* #region: cutting-phase API declarations */
EventCutState buildEventCutState(
    RazorAnalyzerMerged& analyzer,
    RazorHelper* helper,
    const std::string& analysisTag,
    bool isData,
    int runNumber,
    const EventSynthesis& synth);

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
    Bool_t hltL1CSCShowerDTCluster50);

void writeEventCutState(TreeMuonSystemCombination* muonSystem, const EventCutState& cuts);

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
    float& outJetVetoPtJESDown);

GLLPMatchResult findNearestGLLPMatch(
    RazorAnalyzerMerged& analyzer,
    float clusterEta,
    float clusterPhi,
    int nGLLP,
    const float* gLLPEta,
    const float* gLLPPhi);

int dtWheelFromClusterZ(float z);
/* #endregion */

#endif
