#ifndef MUONSYSTEM_FILLING_PHASE_H
#define MUONSYSTEM_FILLING_PHASE_H

/* #region: includes */
#include <vector>

#include "CACluster.h"
#include "MuonSystemPhaseTypes.h"
#include "MuonSystemSignalScanManager.h"
#include "RazorAnalyzer.h"
#include "RazorHelper.h"
/* #endregion */

/* #region: filling-phase API declarations */
void fillLeptonBranches(
    TreeMuonSystemCombination* muonSystem,
    const std::vector<LeptonCandidate>& leptons);

void fillJetBranches(
    TreeMuonSystemCombination* muonSystem,
    const std::vector<JetCandidate>& jets);

void fillPuppiMetJesFromShift(
    RazorAnalyzerMerged& analyzer,
    TreeMuonSystemCombination* muonSystem,
    const JetStageResult& jetStage);

void fillEventIdentityAndNominalWeight(
    TreeMuonSystemCombination* muonSystem,
    TH1F* nEvents,
    bool isData,
    float generatorWeight,
    unsigned int runNumber,
    unsigned int luminosityBlock,
    unsigned long long eventNumber);

bool shouldSkipDataRun(
    bool isData,
    unsigned int runNumber);

void fillMcTruthAndPileupWeights(
    RazorAnalyzerMerged& analyzer,
    RazorHelper* helper,
    TreeMuonSystemCombination* muonSystem);

void fillEventObservables(
    TreeMuonSystemCombination* muonSystem,
    int pvNpvs,
    float rhoFixedGrid,
    float pfMetPt,
    float pfMetPhi,
    float puppiMetPt,
    float puppiMetPhi);

void fillAcceptanceCounters(
    bool isData,
    bool signalScan,
    const SignalEventState& signalState,
    MuonSystemSignalScanManager& signalScanManager,
    TreeMuonSystemCombination* muonSystem,
    TH1F* accepCsccsc,
    TH1F* accepCscdt,
    float generatorWeight);

void fillEventTreeForCurrentMode(
    bool isData,
    bool signalScan,
    TreeMuonSystemCombination* muonSystem,
    MuonSystemSignalScanManager& signalScanManager,
    const SignalEventState& signalState);

void runCAClustering(CACluster& clusterer);

void fillClusteredRecHitIds(
    int nStoredRecHits,
    const std::vector<Rechits>& clusteredPoints,
    int* outClusterIds);

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
    float* outE);

void processCscClusterStage(
    RazorAnalyzerMerged& analyzer,
    RazorHelper* helper,
    TreeMuonSystemCombination* muonSystem,
    const EventSynthesis& synth,
    bool isData,
    int runNumber);

void processDtRpcClusterStage(
    RazorAnalyzerMerged& analyzer,
    RazorHelper* helper,
    TreeMuonSystemCombination* muonSystem,
    const EventSynthesis& synth,
    bool isData,
    int runNumber);
/* #endregion */

#endif
