/* #region: includes */
#include "MuonSystemFillingPhase.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <memory>
#include <random>
#include <vector>

#include "TMath.h"

#include "MuonSystemCuttingPhase.h"
#include "MuonSystemSynthesisPhase.h"
/* #endregion */

/* #region: generic branch-fill utilities */
std::unique_ptr<TH1F> makeCounterHist(const std::string& name) {
  auto h = std::make_unique<TH1F>(name.c_str(), name.c_str(), 1, 1, 2);
  h->SetDirectory(nullptr);
  return h;
}

void fillLeptonBranches(
    TreeMuonSystemCombination* muonSystem,
    const std::vector<LeptonCandidate>& leptons) {
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

void fillJetBranches(
    TreeMuonSystemCombination* muonSystem,
    const std::vector<JetCandidate>& jets) {
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
/* #endregion */

/* #region: CSC clustering and fill stage */
void processCscClusterStage(
    RazorAnalyzerMerged& analyzer,
    RazorHelper* helper,
    TreeMuonSystemCombination* muonSystem,
    const EventSynthesis& synth,
    bool isData,
    int runNumber) {
  TreeMuonSystemCombination* MuonSystem = muonSystem;

  std::vector<Rechits> points = buildCscRechitPoints(
      analyzer.ncscRechits,
      analyzer.cscRechits_Phi,
      analyzer.cscRechits_Eta,
      analyzer.cscRechits_X,
      analyzer.cscRechits_Y,
      analyzer.cscRechits_Z,
      analyzer.cscRechits_Tpeak,
      analyzer.cscRechits_Twire,
      analyzer.cscRechits_Station,
      analyzer.cscRechits_Chamber);
  MuonSystem->nCscRings = countCscRingsFromRecHits(
      analyzer.ncscRechits,
      analyzer.cscRechits_Chamber);
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
        analyzer,
        helper,
        runNumber,
        analyzer.nJet,
        analyzer.Jet_eta,
        analyzer.Jet_phi,
        analyzer.Jet_pt,
        synth.jetE.data(),
        synth.jetPassIDTight.data(),
        MuonSystem->cscRechitClusterEta[MuonSystem->nCscRechitClusters],
        MuonSystem->cscRechitClusterPhi[MuonSystem->nCscRechitClusters],
        MuonSystem->cscRechitClusterJetVetoPt[MuonSystem->nCscRechitClusters],
        MuonSystem->cscRechitClusterJetVetoE[MuonSystem->nCscRechitClusters],
        MuonSystem->cscRechitClusterJetVetoTightId[MuonSystem->nCscRechitClusters],
        MuonSystem->cscRechitClusterJetVetoPtJESUp[MuonSystem->nCscRechitClusters],
        MuonSystem->cscRechitClusterJetVetoPtJESDown[MuonSystem->nCscRechitClusters]);
    for (int i = 0; i < analyzer.nMuon; i++) {
      if (fabs(analyzer.Muon_eta[i]) > 3.0)
        continue;
      if (analyzer.deltaR(analyzer.Muon_eta[i], analyzer.Muon_phi[i], MuonSystem->cscRechitClusterEta[MuonSystem->nCscRechitClusters], MuonSystem->cscRechitClusterPhi[MuonSystem->nCscRechitClusters]) < 0.4 && analyzer.Muon_pt[i] > MuonSystem->cscRechitClusterMuonVetoPt[MuonSystem->nCscRechitClusters]) {
        MuonSystem->cscRechitClusterMuonVetoPt[MuonSystem->nCscRechitClusters] = analyzer.Muon_pt[i];
        MuonSystem->cscRechitClusterMuonVetoE[MuonSystem->nCscRechitClusters] = synth.muonE[i];
        MuonSystem->cscRechitClusterMuonVetoGlobal[MuonSystem->nCscRechitClusters] = analyzer.Muon_isGlobal[i];
        MuonSystem->cscRechitClusterMuonVetoLooseId[MuonSystem->nCscRechitClusters] = analyzer.Muon_looseId[i];
      }
      if (analyzer.deltaR(analyzer.Muon_eta[i], analyzer.Muon_phi[i], MuonSystem->cscRechitClusterEta[MuonSystem->nCscRechitClusters], MuonSystem->cscRechitClusterPhi[MuonSystem->nCscRechitClusters]) < 0.8 && analyzer.Muon_pt[i] > MuonSystem->cscRechitClusterMuonVetoPt0p8Thresh[MuonSystem->nCscRechitClusters]) {
        MuonSystem->cscRechitClusterMuonVetoPt0p8Thresh[MuonSystem->nCscRechitClusters] = analyzer.Muon_pt[i];
        MuonSystem->cscRechitClusterMuonVetoE0p8Thresh[MuonSystem->nCscRechitClusters] = synth.muonE[i];
        MuonSystem->cscRechitClusterMuonVetoGlobal0p8Thresh[MuonSystem->nCscRechitClusters] = analyzer.Muon_isGlobal[i];
        MuonSystem->cscRechitClusterMuonVetoLooseId0p8Thresh[MuonSystem->nCscRechitClusters] = analyzer.Muon_looseId[i];
      }
    }
    if (!isData) {
      // match to gen level LLP
      const GLLPMatchResult match = findNearestGLLPMatch(
          analyzer,
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
    for (int i = 0; i < analyzer.ndtSegments; i++) {
      if (analyzer.deltaR(analyzer.dtSegments_Eta[i], analyzer.dtSegments_Phi[i], MuonSystem->cscRechitClusterEta[MuonSystem->nCscRechitClusters], MuonSystem->cscRechitClusterPhi[MuonSystem->nCscRechitClusters]) < 0.4) {
        MuonSystem->cscRechitCluster_match_dtSeg_0p4[MuonSystem->nCscRechitClusters]++;
        if (analyzer.dtSegments_Station[i] == 1)
          MuonSystem->cscRechitCluster_match_MB1Seg_0p4[MuonSystem->nCscRechitClusters]++;
      }
    }

    //match to RPC hits in RE1/2
    for (int i = 0; i < analyzer.nrpcRecHits; i++) {
      float rpcR = sqrt(analyzer.rpcRecHits_X[i] * analyzer.rpcRecHits_X[i] + analyzer.rpcRecHits_Y[i] * analyzer.rpcRecHits_Y[i]);
      if (analyzer.deltaR(analyzer.rpcRecHits_Eta[i], analyzer.rpcRecHits_Phi[i], MuonSystem->cscRechitClusterEta[MuonSystem->nCscRechitClusters], MuonSystem->cscRechitClusterPhi[MuonSystem->nCscRechitClusters]) < 0.4) {
        if (rpcR < 461.0 && rpcR > 275 && abs(analyzer.rpcRecHits_Z[i]) > 663 && abs(analyzer.rpcRecHits_Z[i]) < 730)
          MuonSystem->cscRechitCluster_match_RE12_0p4[MuonSystem->nCscRechitClusters]++;
        if (rpcR < 470 && rpcR > 380 && abs(analyzer.rpcRecHits_Z[i]) < 661)
          MuonSystem->cscRechitCluster_match_RB1_0p4[MuonSystem->nCscRechitClusters]++;
      }
    }

    MuonSystem->cscRechitClusterMet_dPhi[MuonSystem->nCscRechitClusters] = analyzer.deltaPhi(MuonSystem->cscRechitClusterPhi[MuonSystem->nCscRechitClusters], MuonSystem->PFMET_phi);
    MuonSystem->cscRechitClusterPuppiMet_dPhi[MuonSystem->nCscRechitClusters] = analyzer.deltaPhi(MuonSystem->cscRechitClusterPhi[MuonSystem->nCscRechitClusters], MuonSystem->PuppiMET_phi);
    if (MuonSystem->nTaus == 1) {
      MuonSystem->cscRechitClusterPromptTauDeltaEta[MuonSystem->nCscRechitClusters] = MuonSystem->tauEta[0] - MuonSystem->cscRechitClusterEta[MuonSystem->nCscRechitClusters];
      MuonSystem->cscRechitClusterPromptTauDeltaPhi[MuonSystem->nCscRechitClusters] = analyzer.deltaPhi(MuonSystem->tauPhi[0], MuonSystem->cscRechitClusterPhi[MuonSystem->nCscRechitClusters]);
      MuonSystem->cscRechitClusterPromptTauDeltaR[MuonSystem->nCscRechitClusters] = analyzer.deltaR(MuonSystem->tauEta[0], MuonSystem->tauPhi[0], MuonSystem->cscRechitClusterEta[MuonSystem->nCscRechitClusters], MuonSystem->cscRechitClusterPhi[MuonSystem->nCscRechitClusters]);
    }
    if (MuonSystem->nLeptons > 0) {
      if (abs(MuonSystem->lepPdgId[0]) == 11) {
        MuonSystem->cscRechitClusterPromptEleDeltaEta[MuonSystem->nCscRechitClusters] = MuonSystem->lepEta[0] - MuonSystem->cscRechitClusterEta[MuonSystem->nCscRechitClusters];
        MuonSystem->cscRechitClusterPromptEleDeltaPhi[MuonSystem->nCscRechitClusters] = analyzer.deltaPhi(MuonSystem->lepPhi[0], MuonSystem->cscRechitClusterPhi[MuonSystem->nCscRechitClusters]);
        MuonSystem->cscRechitClusterPromptEleDeltaR[MuonSystem->nCscRechitClusters] = analyzer.deltaR(MuonSystem->lepEta[0], MuonSystem->lepPhi[0], MuonSystem->cscRechitClusterEta[MuonSystem->nCscRechitClusters], MuonSystem->cscRechitClusterPhi[MuonSystem->nCscRechitClusters]);
      }
      if (abs(MuonSystem->lepPdgId[0]) == 13) {
        MuonSystem->cscRechitClusterPromptMuDeltaEta[MuonSystem->nCscRechitClusters] = MuonSystem->lepEta[0] - MuonSystem->cscRechitClusterEta[MuonSystem->nCscRechitClusters];
        MuonSystem->cscRechitClusterPromptMuDeltaPhi[MuonSystem->nCscRechitClusters] = analyzer.deltaPhi(MuonSystem->lepPhi[0], MuonSystem->cscRechitClusterPhi[MuonSystem->nCscRechitClusters]);
        MuonSystem->cscRechitClusterPromptMuDeltaR[MuonSystem->nCscRechitClusters] = analyzer.deltaR(MuonSystem->lepEta[0], MuonSystem->lepPhi[0], MuonSystem->cscRechitClusterEta[MuonSystem->nCscRechitClusters], MuonSystem->cscRechitClusterPhi[MuonSystem->nCscRechitClusters]);
      }
    }

    MuonSystem->nCscRechitClusters++;
  }
}
/* #endregion */

/* #region: DT/RPC clustering and fill stage */
void processDtRpcClusterStage(
    RazorAnalyzerMerged& analyzer,
    RazorHelper* helper,
    TreeMuonSystemCombination* muonSystem,
    const EventSynthesis& synth,
    bool isData,
    int runNumber) {
  TreeMuonSystemCombination* MuonSystem = muonSystem;

  // DT cluster
  std::vector<Rechits> points = buildDtRecHitPoints(
      analyzer.ndtRecHits,
      analyzer.dtRecHits_Phi,
      analyzer.dtRecHits_Eta,
      analyzer.dtRecHits_X,
      analyzer.dtRecHits_Y,
      analyzer.dtRecHits_Z,
      analyzer.dtRecHits_Station,
      analyzer.dtRecHits_Wheel,
      analyzer.dtRecHits_SuperLayer);

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
      analyzer.rpcRecHits_Phi,
      analyzer.rpcRecHits_Eta,
      analyzer.rpcRecHits_X,
      analyzer.rpcRecHits_Y,
      analyzer.rpcRecHits_Z,
      analyzer.rpcRecHits_Time,
      analyzer.rpcRecHits_Station,
      analyzer.rpcRecHits_Sector,
      analyzer.rpcRecHits_Layer,
      analyzer.rpcRecHits_Ring);

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
      if (analyzer.deltaR(MuonSystem->cscRechitClusterEta[i], MuonSystem->cscRechitClusterPhi[i], tmp.eta, tmp.phi) < 0.4)
        overlap = true;
    }
    if (overlap)
      continue;
    MuonSystem->nDtRechitClusters_nocut++;
    MuonSystem->dtRechitClusterX[MuonSystem->nDtRechitClusters] = tmp.x;
    MuonSystem->dtRechitClusterY[MuonSystem->nDtRechitClusters] = tmp.y;
    MuonSystem->dtRechitClusterZ[MuonSystem->nDtRechitClusters] = tmp.z;
    MuonSystem->dtRechitClusterWheel[MuonSystem->nDtRechitClusters] = dtWheelFromClusterZ(tmp.z);
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
    std::default_random_engine generator(seed);

    std::uniform_real_distribution<double> distribution(0.0, 1.0);
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
        analyzer,
        helper,
        runNumber,
        analyzer.nJet,
        analyzer.Jet_eta,
        analyzer.Jet_phi,
        analyzer.Jet_pt,
        synth.jetE.data(),
        synth.jetPassIDTight.data(),
        MuonSystem->dtRechitClusterEta[MuonSystem->nDtRechitClusters],
        MuonSystem->dtRechitClusterPhi[MuonSystem->nDtRechitClusters],
        MuonSystem->dtRechitClusterJetVetoPt[MuonSystem->nDtRechitClusters],
        MuonSystem->dtRechitClusterJetVetoE[MuonSystem->nDtRechitClusters],
        MuonSystem->dtRechitClusterJetVetoTightId[MuonSystem->nDtRechitClusters],
        MuonSystem->dtRechitClusterJetVetoPtJESUp[MuonSystem->nDtRechitClusters],
        MuonSystem->dtRechitClusterJetVetoPtJESDown[MuonSystem->nDtRechitClusters]);

    for (int i = 0; i < analyzer.nMuon; i++) {
      if (fabs(analyzer.Muon_eta[i]) > 3.0)
        continue;
      if (analyzer.deltaR(analyzer.Muon_eta[i], analyzer.Muon_phi[i], MuonSystem->dtRechitClusterEta[MuonSystem->nDtRechitClusters], MuonSystem->dtRechitClusterPhi[MuonSystem->nDtRechitClusters]) < 0.4 && analyzer.Muon_pt[i] > MuonSystem->dtRechitClusterMuonVetoPt[MuonSystem->nDtRechitClusters]) {
        MuonSystem->dtRechitClusterMuonVetoPt[MuonSystem->nDtRechitClusters] = analyzer.Muon_pt[i];
        MuonSystem->dtRechitClusterMuonVetoE[MuonSystem->nDtRechitClusters] = synth.muonE[i];
        MuonSystem->dtRechitClusterMuonVetoGlobal[MuonSystem->nDtRechitClusters] = analyzer.Muon_isGlobal[i];
        MuonSystem->dtRechitClusterMuonVetoLooseId[MuonSystem->nDtRechitClusters] = analyzer.Muon_looseId[i];
      }
    }

    for (int i = 0; i < analyzer.ndtSegments; i++) {
      if (analyzer.deltaR(analyzer.dtSegments_Eta[i], analyzer.dtSegments_Phi[i], MuonSystem->dtRechitClusterEta[MuonSystem->nDtRechitClusters], MuonSystem->dtRechitClusterPhi[MuonSystem->nDtRechitClusters]) < 0.4) {
        if (analyzer.dtSegments_Station[i] == 1)
          MuonSystem->dtRechitClusterNSegStation1[MuonSystem->nDtRechitClusters] += 1;
        if (analyzer.dtSegments_Station[i] == 2)
          MuonSystem->dtRechitClusterNSegStation2[MuonSystem->nDtRechitClusters] += 1;
        if (analyzer.dtSegments_Station[i] == 3)
          MuonSystem->dtRechitClusterNSegStation3[MuonSystem->nDtRechitClusters] += 1;
        if (analyzer.dtSegments_Station[i] == 4)
          MuonSystem->dtRechitClusterNSegStation4[MuonSystem->nDtRechitClusters] += 1;
      }
      if (abs(analyzer.deltaPhi(analyzer.dtSegments_Phi[i], MuonSystem->dtRechitClusterPhi[MuonSystem->nDtRechitClusters])) > 2) {
        if (analyzer.dtSegments_Station[i] == 1)
          MuonSystem->dtRechitClusterNOppositeSegStation1[MuonSystem->nDtRechitClusters] += 1;
        if (analyzer.dtSegments_Station[i] == 2)
          MuonSystem->dtRechitClusterNOppositeSegStation2[MuonSystem->nDtRechitClusters] += 1;
        if (analyzer.dtSegments_Station[i] == 3)
          MuonSystem->dtRechitClusterNOppositeSegStation3[MuonSystem->nDtRechitClusters] += 1;
        if (analyzer.dtSegments_Station[i] == 4)
          MuonSystem->dtRechitClusterNOppositeSegStation4[MuonSystem->nDtRechitClusters] += 1;
      }
    }

    // match to gen-level LLP
    if (!isData) {
      const GLLPMatchResult match = findNearestGLLPMatch(
          analyzer,
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
    for (int i = 0; i < analyzer.ndtRecHits; i++) {
      if (analyzer.deltaR(analyzer.dtRecHits_Eta[i], analyzer.dtRecHits_Phi[i], MuonSystem->dtRechitClusterEta[MuonSystem->nDtRechitClusters], MuonSystem->dtRechitClusterPhi[MuonSystem->nDtRechitClusters]) < 0.5) {
        if (analyzer.dtRecHits_Station[i] == 1)
          MuonSystem->dtRechitCluster_match_MB1hits_0p5[MuonSystem->nDtRechitClusters]++;
      }
      if (analyzer.deltaR(analyzer.dtRecHits_Eta[i], analyzer.dtRecHits_Phi[i], MuonSystem->dtRechitClusterEta[MuonSystem->nDtRechitClusters], MuonSystem->dtRechitClusterPhi[MuonSystem->nDtRechitClusters]) < 0.4) {
        if (analyzer.dtRecHits_Station[i] == 1)
          MuonSystem->dtRechitCluster_match_MB1hits_0p4[MuonSystem->nDtRechitClusters]++;
      }
      if (abs(analyzer.dtRecHits_Wheel[i] - MuonSystem->dtRechitClusterWheel[MuonSystem->nDtRechitClusters]) == 1 && analyzer.dtRecHits_Station[i] == 1) {
        if (abs(analyzer.deltaPhi(analyzer.dtRecHits_Phi[i], MuonSystem->dtRechitClusterPhi[MuonSystem->nDtRechitClusters])) < TMath::Pi() / 4.0) {
          if (analyzer.dtRecHits_Wheel[i] - MuonSystem->dtRechitClusterWheel[MuonSystem->nDtRechitClusters] == 1)
            MuonSystem->dtRechitCluster_match_MB1hits_cosmics_plus[MuonSystem->nDtRechitClusters]++;
          else
            MuonSystem->dtRechitCluster_match_MB1hits_cosmics_minus[MuonSystem->nDtRechitClusters]++;
        }
      }
    }

    std::vector<int> dtRechitCluster_match_rpcBx;

    //match to RPC hits with dPhi<0.5 and same wheel in DT
    for (int i = 0; i < analyzer.nrpcRecHits; i++) {
      float rpcR = sqrt(analyzer.rpcRecHits_X[i] * analyzer.rpcRecHits_X[i] + analyzer.rpcRecHits_Y[i] * analyzer.rpcRecHits_Y[i]);
      if (analyzer.rpcRecHits_Region[i] != 0)
        continue;
      if (abs(analyzer.deltaPhi(analyzer.rpcRecHits_Phi[i], MuonSystem->dtRechitClusterPhi[MuonSystem->nDtRechitClusters])) < 0.5) {
        if (analyzer.rpcRecHits_Ring[i] == MuonSystem->dtRechitClusterWheel[MuonSystem->nDtRechitClusters]) {
          dtRechitCluster_match_rpcBx.push_back(analyzer.rpcRecHits_Bx[i]);
          MuonSystem->dtRechitCluster_match_RPChits_dPhi0p5[MuonSystem->nDtRechitClusters]++;
          if (rpcR < 470 && rpcR > 380 && abs(analyzer.rpcRecHits_Z[i]) < 661)
            MuonSystem->dtRechitCluster_match_RB1_dPhi0p5[MuonSystem->nDtRechitClusters]++;
        }
      }
      if (analyzer.deltaR(analyzer.rpcRecHits_Eta[i], analyzer.rpcRecHits_Phi[i], MuonSystem->dtRechitClusterEta[MuonSystem->nDtRechitClusters], MuonSystem->dtRechitClusterPhi[MuonSystem->nDtRechitClusters]) < 0.4) {
        if (rpcR < 470 && rpcR > 380 && abs(analyzer.rpcRecHits_Z[i]) < 661)
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

    MuonSystem->dtRechitClusterMet_dPhi[MuonSystem->nDtRechitClusters] = analyzer.deltaPhi(MuonSystem->dtRechitClusterPhi[MuonSystem->nDtRechitClusters], MuonSystem->PFMET_phi);
    MuonSystem->dtRechitClusterPuppiMet_dPhi[MuonSystem->nDtRechitClusters] = analyzer.deltaPhi(MuonSystem->dtRechitClusterPhi[MuonSystem->nDtRechitClusters], MuonSystem->PuppiMET_phi);

    MuonSystem->nDtRechitClusters++;
  }
}
/* #endregion */
