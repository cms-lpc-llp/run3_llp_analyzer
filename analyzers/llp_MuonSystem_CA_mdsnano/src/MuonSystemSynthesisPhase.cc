/* #region: includes */
#include "MuonSystemSynthesisPhase.h"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <limits>

#include "TMath.h"
/* #endregion */

/* #region: local comparators */
namespace {

struct largest_pt {
  inline bool operator()(const LeptonCandidate& p1, const LeptonCandidate& p2) {
    return p1.lepton.Pt() > p2.lepton.Pt();
  }
} my_largest_pt;

struct largest_pt_jet {
  inline bool operator()(const JetCandidate& p1, const JetCandidate& p2) {
    return p1.jet.Pt() > p2.jet.Pt();
  }
} my_largest_pt_jet;

} // namespace
/* #endregion */

/* #region: event-object synthesis */
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
    const float* electronDz) {
  std::vector<LeptonCandidate> leptonsOut;

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

    LeptonCandidate tmpMuon;
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

    LeptonCandidate tmpElectron;
    tmpElectron.lepton.SetPtEtaPhiM(electronPt[i], electronEta[i], electronPhi[i], ELE_MASS);
    tmpElectron.pdgId = 11 * -1 * electronCharge[i];
    tmpElectron.dZ = electronDz[i];
    tmpElectron.passId = synth.elePassCutBasedIDTight[i];
    leptonsOut.push_back(tmpElectron);
  }

  sort(leptonsOut.begin(), leptonsOut.end(), my_largest_pt);
  return leptonsOut;
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
    const std::vector<LeptonCandidate>& leptons) {
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

    JetCandidate outJet;
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
/* #endregion */

/* #region: raw rechit copy and ring summaries */
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
    int threshold) {
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
    int threshold) {
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
/* #endregion */

/* #region: rechit point builders */
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
/* #endregion */

/* #region: signal-point parsing */
bool parseSignalPointFromInputPath(
    const std::string& inputPath,
    int& mh,
    int& mx,
    float& ctau,
    std::string& ctauTokenTag) {
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
/* #endregion */
