#ifndef MUONSYSTEM_SYNTHESIS_PHASE_H
#define MUONSYSTEM_SYNTHESIS_PHASE_H

/* #region: includes */
#include <string>
#include <vector>

#include "CACluster.h"
#include "MuonSystemPhaseTypes.h"
#include "RazorAnalyzer.h"
#include "RazorHelper.h"
/* #endregion */

/* #region: synthesis-phase API declarations */
EventSynthesis buildEventSynthesis(
    RazorAnalyzerMerged& analyzer,
    RazorHelper* helper,
    const std::string& analysisTag);

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
    const float* muonPt);

std::vector<LeptonCandidate> buildSelectedLeptons(
    RazorAnalyzerMerged& analyzer,
    const EventSynthesis& synth);

std::vector<LeptonCandidate> buildSelectedLeptons(
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
    const float* electronDz);

JetStageResult buildJetStageResult(
    RazorAnalyzerMerged& analyzer,
    RazorHelper* helper,
    int runNumber,
    const EventSynthesis& synth,
    const std::vector<LeptonCandidate>& leptons);

JetStageResult buildJetStageResult(
    RazorAnalyzerMerged& analyzer,
    RazorHelper* helper,
    int runNumber,
    int nJet,
    const float* jetEta,
    const float* jetPhi,
    const float* jetPt,
    const EventSynthesis& synth,
    const std::vector<LeptonCandidate>& leptons);

void fillRawRechits(
    TreeMuonSystemCombination* muonSystem,
    RazorAnalyzerMerged& analyzer);

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
    const float* rpcRecHitsTimeError);

int countDtRingsFromRecHits(
    int ndtRecHits,
    const int* dtRecHitsStation,
    const int* dtRecHitsWheel,
    int threshold = 50);

int countCscRingsFromRecHits(
    int ncscRechits,
    const int* cscRechitsChamber,
    int threshold = 50);

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
    const int* cscRechitsChamber);

std::vector<Rechits> buildDtRecHitPoints(
    int ndtRecHits,
    const float* dtRecHitsPhi,
    const float* dtRecHitsEta,
    const float* dtRecHitsX,
    const float* dtRecHitsY,
    const float* dtRecHitsZ,
    const int* dtRecHitsStation,
    const int* dtRecHitsWheel,
    const int* dtRecHitsSuperLayer);

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
    const int* rpcRecHitsRing);

bool parseSignalPointFromInputPath(
    const std::string& inputPath,
    int& mh,
    int& mx,
    float& ctau,
    std::string& ctauTokenTag);
/* #endregion */

#endif
