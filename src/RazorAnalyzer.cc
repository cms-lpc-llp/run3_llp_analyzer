#include "RazorAnalyzer.h"
#include "TLorentzVector.h"

using namespace std;

// ======================================================================================================
// This part is for the "analyzer" for the NANO AOD merged ntuples
// Dirty implementation, but not more dirty than the rest of the codebase
// ======================================================================================================

RazorAnalyzerMerged::RazorAnalyzerMerged(TTree* tree)
    : merged_event(tree) {
  // turn off all branches
  //  fChain->SetBranchStatus("*", 1);
  fChain->SetBranchStatus("*", 0);
}

RazorAnalyzerMerged::~RazorAnalyzerMerged() {
}

void RazorAnalyzerMerged::EnableAll() {
  fChain->SetBranchStatus("*", 1);
}

void RazorAnalyzerMerged::Analyze(bool isData, int option, string outputFileName, string label) {
  std::cerr << "Virtual RazorAnalyzerMerged::Analyze" << std::endl;
  assert(0);
}

// As there are already 10^9 copies of the same definitions in the codebase
// I hope no one will be mad if I just add yet another meaningless copy of the same code :)

double RazorAnalyzerMerged::deltaPhi(double phi1, double phi2) {
  double dphi = phi1 - phi2;
  while (dphi > TMath::Pi())
    dphi -= TMath::TwoPi();
  while (dphi <= -TMath::Pi())
    dphi += TMath::TwoPi();
  return dphi;
}

double RazorAnalyzerMerged::deltaR(double eta1, double phi1, double eta2, double phi2) {
  double dphi = deltaPhi(phi1, phi2);
  double deta = eta1 - eta2;
  return sqrt(dphi * dphi + deta * deta);
}

TLorentzVector RazorAnalyzerMerged::makeTLorentzVector(double pt, double eta, double phi, double energy) {
  TLorentzVector vec;
  vec.SetPtEtaPhiE(pt, eta, phi, energy);
  return vec;
}

TLorentzVector RazorAnalyzerMerged::makeTLorentzVectorPtEtaPhiM(double pt, double eta, double phi, double mass) {
  TLorentzVector vec;
  vec.SetPtEtaPhiM(pt, eta, phi, mass);
  return vec;
}

// ======================================================================================================
// The original analyzer for the plain ntuples
// ======================================================================================================

RazorAnalyzer::RazorAnalyzer(TTree* tree)
    : llp_event(tree) {
  // turn off all branches
  //  fChain->SetBranchStatus("*", 1);
  fChain->SetBranchStatus("*", 0);
}

RazorAnalyzer::~RazorAnalyzer() {
}

void RazorAnalyzer::Analyze(bool isData, int option, string outputFileName, string label) {
  cout << "Analyze method called on base RazorAnalyzer instance.  Parameters were: " << isData << " " << option << " " << outputFileName << " " << label << endl;
}

// NOTE: the functions below need to be maintained by hand.  If variables are added or removed from the ntuple, these functions need to be updated to reflect the changes.

void RazorAnalyzer::EnableAll() {
  // fChain->SetBranchStatus("*", 1);
  EnableEventInfo();
  EnableMuons();
  EnableElectrons();
  EnableCSC();
  EnableDT();
  EnableJets();
  EnableMet();
  EnablePileup();
  EnableGenParticles();
  EnableLLP();
  EnableMC();
}

void RazorAnalyzer::EnableEventInfo() {
  fChain->SetBranchStatus("nPV", 1);
  fChain->SetBranchStatus("pvX", 1);
  fChain->SetBranchStatus("pvY", 1);
  fChain->SetBranchStatus("pvZ", 1);
  fChain->SetBranchStatus("isData", 1);
  fChain->SetBranchStatus("runNum", 1);
  fChain->SetBranchStatus("lumiNum", 1);
  fChain->SetBranchStatus("eventNum", 1);
  fChain->SetBranchStatus("eventTime", 1);
  fChain->SetBranchStatus("fixedGridRhoAll", 1);
  fChain->SetBranchStatus("fixedGridRhoFastjetAll", 1);
  fChain->SetBranchStatus("fixedGridRhoFastjetAllCalo", 1);
  fChain->SetBranchStatus("fixedGridRhoFastjetCentralCalo", 1);
  fChain->SetBranchStatus("fixedGridRhoFastjetCentralChargedPileUp", 1);
  fChain->SetBranchStatus("fixedGridRhoFastjetCentralNeutral", 1);
  fChain->SetBranchStatus("HLTDecision", 1);
  fChain->SetBranchStatus("HLTPrescale", 1);
  fChain->SetBranchStatus("lheComments", 1);
}

void RazorAnalyzer::EnablePVAll() {
  fChain->SetBranchStatus("nPVAll", 1);
  fChain->SetBranchStatus("pvAllX", 1);
  fChain->SetBranchStatus("pvAllY", 1);
  fChain->SetBranchStatus("pvAllZ", 1);
}

void RazorAnalyzer::EnablePileup() {
  fChain->SetBranchStatus("nBunchXing", 1);
  fChain->SetBranchStatus("BunchXing", 1);
  fChain->SetBranchStatus("nPU", 1);
  fChain->SetBranchStatus("nPUmean", 1);
}

void RazorAnalyzer::EnableMuons() {
  fChain->SetBranchStatus("nMuons", 1);
  fChain->SetBranchStatus("muonE", 1);
  fChain->SetBranchStatus("muonPt", 1);
  fChain->SetBranchStatus("muonEta", 1);
  fChain->SetBranchStatus("muonPhi", 1);
  fChain->SetBranchStatus("muonCharge", 1);
  fChain->SetBranchStatus("muonIsLoose", 1);
  fChain->SetBranchStatus("muonIsMedium", 1);
  fChain->SetBranchStatus("muonIsTight", 1);
  fChain->SetBranchStatus("muon_d0", 1);
  fChain->SetBranchStatus("muon_dZ", 1);
  fChain->SetBranchStatus("muon_ip3d", 1);
  fChain->SetBranchStatus("muon_ip3dSignificance", 1);
  fChain->SetBranchStatus("muonType", 1);
  fChain->SetBranchStatus("muonQuality", 1);
  fChain->SetBranchStatus("muon_pileupIso", 1);
  fChain->SetBranchStatus("muon_chargedIso", 1);
  fChain->SetBranchStatus("muon_photonIso", 1);
  fChain->SetBranchStatus("muon_neutralHadIso", 1);
  fChain->SetBranchStatus("muon_ptrel", 1);
  fChain->SetBranchStatus("muon_chargedMiniIso", 1);
  fChain->SetBranchStatus("muon_photonAndNeutralHadronMiniIso", 1);
  fChain->SetBranchStatus("muon_chargedPileupMiniIso", 1);
  fChain->SetBranchStatus("muon_activityMiniIsoAnnulus", 1);
  fChain->SetBranchStatus("muon_validFractionTrackerHits", 1);
  fChain->SetBranchStatus("muon_isGlobal", 1);
  fChain->SetBranchStatus("muon_normChi2", 1);
  fChain->SetBranchStatus("muon_chi2LocalPosition", 1);
  fChain->SetBranchStatus("muon_kinkFinder", 1);
  fChain->SetBranchStatus("muon_segmentCompatability", 1);
  fChain->SetBranchStatus("muonIsICHEPMedium", 1);
  fChain->SetBranchStatus("muon_passSingleMuTagFilter", 1);
  fChain->SetBranchStatus("muon_passHLTFilter", 1);
}

void RazorAnalyzer::EnableElectrons() {
  fChain->SetBranchStatus("nElectrons", 1);
  fChain->SetBranchStatus("eleE", 1);
  fChain->SetBranchStatus("elePt", 1);
  fChain->SetBranchStatus("eleEta", 1);
  fChain->SetBranchStatus("elePhi", 1);
  fChain->SetBranchStatus("eleCharge", 1);
  fChain->SetBranchStatus("eleEta_SC", 1);

  fChain->SetBranchStatus("eleSigmaIetaIeta", 1);
  fChain->SetBranchStatus("eleFull5x5SigmaIetaIeta", 1);
  fChain->SetBranchStatus("eleR9", 1);
  fChain->SetBranchStatus("ele_dEta", 1);
  fChain->SetBranchStatus("ele_dPhi", 1);
  fChain->SetBranchStatus("ele_HoverE", 1);
  fChain->SetBranchStatus("ele_d0", 1);
  fChain->SetBranchStatus("ele_dZ", 1);
  fChain->SetBranchStatus("ele_passCutBasedIDLoose", 1);
  fChain->SetBranchStatus("ele_passCutBasedIDTight", 1);

  fChain->SetBranchStatus("ele_ip3d", 1);
  fChain->SetBranchStatus("ele_ip3dSignificance", 1);
  fChain->SetBranchStatus("ele_pileupIso", 1);
  fChain->SetBranchStatus("ele_chargedIso", 1);
  fChain->SetBranchStatus("ele_photonIso", 1);
  fChain->SetBranchStatus("ele_neutralHadIso", 1);
  fChain->SetBranchStatus("ele_MissHits", 1);
  fChain->SetBranchStatus("ele_PassConvVeto", 1);
  fChain->SetBranchStatus("ele_OneOverEminusOneOverP", 1);
  fChain->SetBranchStatus("ele_IDMVAHZZ", 1);
  fChain->SetBranchStatus("ele_IDMVAGeneralPurpose", 1);
  fChain->SetBranchStatus("ele_RegressionE", 1);
  fChain->SetBranchStatus("ele_CombineP4", 1);
  fChain->SetBranchStatus("ele_ptrel", 1);
  fChain->SetBranchStatus("ele_chargedMiniIso", 1);
  fChain->SetBranchStatus("ele_photonAndNeutralHadronMiniIso", 1);
  fChain->SetBranchStatus("ele_chargedPileupMiniIso", 1);
  fChain->SetBranchStatus("ele_activityMiniIsoAnnulus", 1);
  fChain->SetBranchStatus("ele_passSingleEleTagFilter", 1);
  fChain->SetBranchStatus("ele_passTPOneTagFilter", 1);
  fChain->SetBranchStatus("ele_passTPTwoTagFilter", 1);
  fChain->SetBranchStatus("ele_passTPOneProbeFilter", 1);
  fChain->SetBranchStatus("ele_passTPTwoProbeFilter", 1);
  fChain->SetBranchStatus("ele_passHLTFilter", 1);
}

void RazorAnalyzer::EnableTaus() {
  fChain->SetBranchStatus("nTaus", 1);
  fChain->SetBranchStatus("tauE", 1);
  fChain->SetBranchStatus("tauPt", 1);
  fChain->SetBranchStatus("tauEta", 1);
  fChain->SetBranchStatus("tauPhi", 1);
  fChain->SetBranchStatus("tau_IsLoose", 1);
  fChain->SetBranchStatus("tau_IsMedium", 1);
  fChain->SetBranchStatus("tau_IsTight", 1);
  fChain->SetBranchStatus("tau_passEleVetoLoose", 1);
  fChain->SetBranchStatus("tau_passEleVetoMedium", 1);
  fChain->SetBranchStatus("tau_passEleVetoTight", 1);
  fChain->SetBranchStatus("tau_passMuVetoLoose", 1);
  fChain->SetBranchStatus("tau_passMuVetoMedium", 1);
  fChain->SetBranchStatus("tau_passMuVetoTight", 1);
  fChain->SetBranchStatus("tau_ID", 1);
  fChain->SetBranchStatus("tau_combinedIsoDeltaBetaCorr3Hits", 1);
  fChain->SetBranchStatus("tau_chargedIsoPtSum", 1);
  fChain->SetBranchStatus("tau_neutralIsoPtSum", 1);
  fChain->SetBranchStatus("tau_puCorrPtSum", 1);
  fChain->SetBranchStatus("tau_eleVetoMVA", 1);
  fChain->SetBranchStatus("tau_eleVetoCategory", 1);
  fChain->SetBranchStatus("tau_muonVetoMVA", 1);
  fChain->SetBranchStatus("tau_isoMVAnewDMwLT", 1);
  fChain->SetBranchStatus("tau_isoMVAnewDMwoLT", 1);
  fChain->SetBranchStatus("tau_leadCandPt", 1);
  fChain->SetBranchStatus("tau_leadCandID", 1);
  fChain->SetBranchStatus("tau_leadChargedHadrCandPt", 1);
  fChain->SetBranchStatus("tau_leadChargedHadrCandID", 1);
}

void RazorAnalyzer::EnableIsoPFCandidates() {
  // fChain->SetBranchStatus("nIsoPFCandidates", 1);
  // fChain->SetBranchStatus("isoPFCandidatePt", 1);
  // fChain->SetBranchStatus("isoPFCandidateEta", 1);
  // fChain->SetBranchStatus("isoPFCandidatePhi", 1);
  // fChain->SetBranchStatus("isoPFCandidateIso04", 1);
  // fChain->SetBranchStatus("isoPFCandidateD0", 1);
  // fChain->SetBranchStatus("isoPFCandidatePdgId", 1);
}

void RazorAnalyzer::EnablePhotons() {
  fChain->SetBranchStatus("nPhotons", 1);
  fChain->SetBranchStatus("phoE", 1);
  fChain->SetBranchStatus("phoPt", 1);
  fChain->SetBranchStatus("phoEta", 1);
  fChain->SetBranchStatus("phoPhi", 1);
  fChain->SetBranchStatus("phoSigmaIetaIeta", 1);
  fChain->SetBranchStatus("phoFull5x5SigmaIetaIeta", 1);
  fChain->SetBranchStatus("phoR9", 1);
  fChain->SetBranchStatus("pho_HoverE", 1);
  fChain->SetBranchStatus("pho_sumChargedHadronPt", 1);
  fChain->SetBranchStatus("pho_sumNeutralHadronEt", 1);
  fChain->SetBranchStatus("pho_sumPhotonEt", 1);
  fChain->SetBranchStatus("pho_ecalPFClusterIso", 1);
  fChain->SetBranchStatus("pho_hcalPFClusterIso", 1);
  fChain->SetBranchStatus("pho_trkSumPtHollowConeDR03", 1);
  fChain->SetBranchStatus("pho_sumWorstVertexChargedHadronPt", 1);
  fChain->SetBranchStatus("pho_pfIsoChargedHadronIso", 1);
  fChain->SetBranchStatus("pho_pfIsoChargedHadronIsoWrongVtx", 1);
  fChain->SetBranchStatus("pho_pfIsoNeutralHadronIso", 1);
  fChain->SetBranchStatus("pho_pfIsoPhotonIso", 1);
  fChain->SetBranchStatus("pho_pfIsoModFrixione", 1);
  fChain->SetBranchStatus("pho_pfIsoSumPUPt", 1);
  fChain->SetBranchStatus("pho_isConversion", 1);
  fChain->SetBranchStatus("pho_passEleVeto", 1);
  fChain->SetBranchStatus("pho_RegressionE", 1);
  fChain->SetBranchStatus("pho_RegressionEUncertainty", 1);
  fChain->SetBranchStatus("pho_IDMVA", 1);
  fChain->SetBranchStatus("pho_superClusterEnergy", 1);
  fChain->SetBranchStatus("pho_superClusterRawEnergy", 1);
  fChain->SetBranchStatus("pho_superClusterEta", 1);
  fChain->SetBranchStatus("pho_superClusterPhi", 1);
  fChain->SetBranchStatus("pho_superClusterX", 1);
  fChain->SetBranchStatus("pho_superClusterY", 1);
  fChain->SetBranchStatus("pho_superClusterZ", 1);
  fChain->SetBranchStatus("pho_hasPixelSeed", 1);
  fChain->SetBranchStatus("pho_isStandardPhoton", 1);
  fChain->SetBranchStatus("pho_passHLTFilter", 1);
  fChain->SetBranchStatus("pho_convType", 1);
  fChain->SetBranchStatus("pho_convTrkZ", 1);
  fChain->SetBranchStatus("pho_convTrkClusZ", 1);
  fChain->SetBranchStatus("pho_vtxSumPx", 1);
  fChain->SetBranchStatus("pho_vtxSumPy", 1);
  fChain->SetBranchStatus("pho_seedRecHitSwitchToGain6", 1);
  fChain->SetBranchStatus("pho_seedRecHitSwitchToGain1", 1);
  fChain->SetBranchStatus("pho_anyRecHitSwitchToGain6", 1);
  fChain->SetBranchStatus("pho_anyRecHitSwitchToGain1", 1);
};

void RazorAnalyzer::EnableCaloJets() {
  fChain->SetBranchStatus("nCaloJets", 1);
  fChain->SetBranchStatus("calojetE", 1);
  fChain->SetBranchStatus("calojetPt", 1);
  fChain->SetBranchStatus("calojetEta", 1);
  fChain->SetBranchStatus("calojetPhi", 1);
};

void RazorAnalyzer::EnableJets() {
  fChain->SetBranchStatus("nJets", 1);
  fChain->SetBranchStatus("jetE", 1);
  fChain->SetBranchStatus("jetPt", 1);
  fChain->SetBranchStatus("jetEta", 1);
  fChain->SetBranchStatus("jetPhi", 1);
  fChain->SetBranchStatus("jetCSV", 1);
  fChain->SetBranchStatus("jetCISV", 1);
  fChain->SetBranchStatus("jetMass", 1);
  fChain->SetBranchStatus("jetJetArea", 1);
  fChain->SetBranchStatus("jetPileupE", 1);
  fChain->SetBranchStatus("jetPileupId", 1);
  fChain->SetBranchStatus("jetPileupIdFlag", 1);
  fChain->SetBranchStatus("jetPassIDLoose", 1);
  fChain->SetBranchStatus("jetPassIDTight", 1);
  fChain->SetBranchStatus("jetPassMuFrac", 1);
  fChain->SetBranchStatus("jetPassEleFrac", 1);
  fChain->SetBranchStatus("jetPartonFlavor", 1);
  fChain->SetBranchStatus("jetHadronFlavor", 1);
  fChain->SetBranchStatus("jetChargedHadronEnergyFraction", 1);
  fChain->SetBranchStatus("jetNeutralHadronEnergyFraction", 1);
  fChain->SetBranchStatus("jetMuonEnergyFraction", 1);
  fChain->SetBranchStatus("jetHOEnergyFraction", 1);
  fChain->SetBranchStatus("jetHFHadronEnergyFraction", 1);
  fChain->SetBranchStatus("jetHFEMEnergyFraction", 1);
  fChain->SetBranchStatus("jetElectronEnergyFraction", 1);
  fChain->SetBranchStatus("jetPhotonEnergyFraction", 1);
  fChain->SetBranchStatus("jetChargedHadronMultiplicity", 1);
  fChain->SetBranchStatus("jetNeutralHadronMultiplicity", 1);
  fChain->SetBranchStatus("jetPhotonMultiplicity", 1);
  fChain->SetBranchStatus("jetElectronMultiplicity", 1);
  fChain->SetBranchStatus("jetMuonMultiplicity", 1);
  fChain->SetBranchStatus("jetAllMuonPt", 1);
  fChain->SetBranchStatus("jetAllMuonEta", 1);
  fChain->SetBranchStatus("jetAllMuonPhi", 1);
  fChain->SetBranchStatus("jetAllMuonM", 1);
  fChain->SetBranchStatus("jetPtWeightedDZ", 1);
  fChain->SetBranchStatus("jetNRechits", 1);
  fChain->SetBranchStatus("jetRechitE", 1);
  fChain->SetBranchStatus("jetRechitT", 1);
};

void RazorAnalyzer::EnableFatJets() {
  fChain->SetBranchStatus("nFatJets", 1);
  fChain->SetBranchStatus("fatJetE", 1);
  fChain->SetBranchStatus("fatJetPt", 1);
  fChain->SetBranchStatus("fatJetEta", 1);
  fChain->SetBranchStatus("fatJetPhi", 1);
  fChain->SetBranchStatus("fatJetCorrectedPt", 1);
  fChain->SetBranchStatus("fatJetPrunedM", 1);
  fChain->SetBranchStatus("fatJetTrimmedM", 1);
  fChain->SetBranchStatus("fatJetFilteredM", 1);
  fChain->SetBranchStatus("fatJetSoftDropM", 1);
  fChain->SetBranchStatus("fatJetCorrectedSoftDropM", 1);
  fChain->SetBranchStatus("fatJetUncorrectedSoftDropM", 1);
  fChain->SetBranchStatus("fatJetTau1", 1);
  fChain->SetBranchStatus("fatJetTau2", 1);
  fChain->SetBranchStatus("fatJetTau3", 1);
  fChain->SetBranchStatus("fatJetMaxSubjetCSV", 1);
  fChain->SetBranchStatus("fatJetPassIDLoose", 1);
  fChain->SetBranchStatus("fatJetPassIDTight", 1);
}

void RazorAnalyzer::EnableMet() {
  fChain->SetBranchStatus("metPt", 1);
  fChain->SetBranchStatus("metPhi", 1);
  fChain->SetBranchStatus("metType1Pt", 1);
  fChain->SetBranchStatus("metType1Phi", 1);

  fChain->SetBranchStatus("metPuppiPt", 1);
  fChain->SetBranchStatus("metPuppiPhi", 1);
  fChain->SetBranchStatus("metCaloPt", 1);
  fChain->SetBranchStatus("metCaloPhi", 1);
  fChain->SetBranchStatus("sumMET", 1);
  fChain->SetBranchStatus("Flag_HBHENoiseFilter", 1);
  fChain->SetBranchStatus("Flag_HBHETightNoiseFilter", 1);
  fChain->SetBranchStatus("Flag_HBHEIsoNoiseFilter", 1);
  fChain->SetBranchStatus("Flag_badChargedCandidateFilter", 1);
  fChain->SetBranchStatus("Flag_BadChargedCandidateFilter", 1);
  fChain->SetBranchStatus("Flag_badMuonFilter", 1);
  fChain->SetBranchStatus("Flag_BadPFMuonFilter", 1);
  fChain->SetBranchStatus("Flag_badGlobalMuonFilter", 1);
  fChain->SetBranchStatus("Flag_duplicateMuonFilter", 1);
  fChain->SetBranchStatus("Flag_globalSuperTightHalo2016Filter", 1);
  fChain->SetBranchStatus("Flag_hcalLaserEventFilter", 1);
  fChain->SetBranchStatus("Flag_EcalDeadCellTriggerPrimitiveFilter", 1);
  fChain->SetBranchStatus("Flag_EcalDeadCellBoundaryEnergyFilter", 1);
  fChain->SetBranchStatus("Flag_goodVertices", 1);
  fChain->SetBranchStatus("Flag_trackingFailureFilter", 1);
  fChain->SetBranchStatus("Flag_eeBadScFilter", 1);
  fChain->SetBranchStatus("Flag_ecalLaserCorrFilter", 1);
  fChain->SetBranchStatus("Flag_trkPOGFilters", 1);
  fChain->SetBranchStatus("Flag_trkPOG_manystripclus53X", 1);
  fChain->SetBranchStatus("Flag_trkPOG_toomanystripclus53X", 1);
  fChain->SetBranchStatus("Flag_trkPOG_logErrorTooManyClusters", 1);
  fChain->SetBranchStatus("Flag_METFilters", 1);
  fChain->SetBranchStatus("Flag_ecalBadCalibFilter", 1);

  fChain->SetBranchStatus("Flag2_HBHENoiseFilter", 1);
  fChain->SetBranchStatus("Flag2_HBHEIsoNoiseFilter", 1);
  fChain->SetBranchStatus("Flag2_BadPFMuonFilter", 1);
  fChain->SetBranchStatus("Flag2_globalSuperTightHalo2016Filter", 1);
  fChain->SetBranchStatus("Flag2_globalTightHalo2016Filter", 1);
  fChain->SetBranchStatus("Flag2_BadChargedCandidateFilter", 1);
  fChain->SetBranchStatus("Flag2_EcalDeadCellTriggerPrimitiveFilter", 1);
  fChain->SetBranchStatus("Flag2_ecalBadCalibFilter", 1);
  fChain->SetBranchStatus("Flag2_eeBadScFilter", 1);
}

void RazorAnalyzer::EnableRazor() {

};
void RazorAnalyzer::EnableDT() {
  fChain->SetBranchStatus("nDtRechits", 1);
  fChain->SetBranchStatus("dtRechitStation", 1);
  fChain->SetBranchStatus("dtRechitWheel", 1);
  fChain->SetBranchStatus("dtRechitCorrectX", 1);
  fChain->SetBranchStatus("dtRechitCorrectY", 1);
  fChain->SetBranchStatus("dtRechitCorrectZ", 1);
  fChain->SetBranchStatus("dtRechitCorrectEta", 1);
  fChain->SetBranchStatus("dtRechitCorrectPhi", 1);
  fChain->SetBranchStatus("dtRechitTime", 1);
  fChain->SetBranchStatus("dtRechitSuperLayer", 1);
  fChain->SetBranchStatus("dtRechitLayer", 1);

  fChain->SetBranchStatus("nDtSeg", 1);
  fChain->SetBranchStatus("dtSegPhi", 1);
  fChain->SetBranchStatus("dtSegEta", 1);
  fChain->SetBranchStatus("dtSegX", 1);
  fChain->SetBranchStatus("dtSegY", 1);
  fChain->SetBranchStatus("dtSegZ", 1);
  fChain->SetBranchStatus("dtSegStation", 1);
  fChain->SetBranchStatus("dtSegWheel", 1);
  fChain->SetBranchStatus("dtSegTime", 1);

  fChain->SetBranchStatus("nRpc", 1);
  fChain->SetBranchStatus("rpcPhi", 1);
  fChain->SetBranchStatus("rpcEta", 1);
  fChain->SetBranchStatus("rpcX", 1);
  fChain->SetBranchStatus("rpcY", 1);
  fChain->SetBranchStatus("rpcZ", 1);
  fChain->SetBranchStatus("rpcT", 1);
  fChain->SetBranchStatus("rpcBx", 1);
  fChain->SetBranchStatus("rpcTError", 1);
  fChain->SetBranchStatus("rpcRegion", 1);
  fChain->SetBranchStatus("rpcRing", 1);
  fChain->SetBranchStatus("rpcStation", 1);
  fChain->SetBranchStatus("rpcSector", 1);
  fChain->SetBranchStatus("rpcLayer", 1);
  fChain->SetBranchStatus("rpcTError", 1);
};
void RazorAnalyzer::EnableCSC() {
  fChain->SetBranchStatus("ncscRechits", 1);
  fChain->SetBranchStatus("cscRechitsPhi", 1);
  fChain->SetBranchStatus("cscRechitsEta", 1);
  fChain->SetBranchStatus("cscRechitsQuality", 1);

  fChain->SetBranchStatus("cscRechitsX", 1);
  fChain->SetBranchStatus("cscRechitsY", 1);
  fChain->SetBranchStatus("cscRechitsZ", 1);

  fChain->SetBranchStatus("cscRechitsStation", 1);
  fChain->SetBranchStatus("cscRechitsChamber", 1);
  fChain->SetBranchStatus("cscRechitsTwire", 1);
  fChain->SetBranchStatus("cscRechitsTpeak", 1);
  fChain->SetBranchStatus("cscRechitsDetId", 1);

  fChain->SetBranchStatus("nCscSeg", 1);
  fChain->SetBranchStatus("cscSegPhi", 1);
  fChain->SetBranchStatus("cscSegEta", 1);
  fChain->SetBranchStatus("cscSegX", 1);
  fChain->SetBranchStatus("cscSegY", 1);
  fChain->SetBranchStatus("cscSegZ", 1);
  fChain->SetBranchStatus("cscSegT", 1);
  fChain->SetBranchStatus("cscSegChi2", 1);
  fChain->SetBranchStatus("cscSegChamber", 1);
  fChain->SetBranchStatus("cscSegStation", 1);
  fChain->SetBranchStatus("cscSegNRecHits", 1);
};

void RazorAnalyzer::EnableMC() {
  fChain->SetBranchStatus("nGenJets", 1);
  fChain->SetBranchStatus("genJetE", 1);
  fChain->SetBranchStatus("genJetPt", 1);
  fChain->SetBranchStatus("genJetEta", 1);
  fChain->SetBranchStatus("genJetPhi", 1);
  fChain->SetBranchStatus("genJetMET", 1);

  fChain->SetBranchStatus("genMetPtCalo", 1);
  fChain->SetBranchStatus("genMetPhiCalo", 1);
  fChain->SetBranchStatus("genMetPtTrue", 1);
  fChain->SetBranchStatus("genMetPhiTrue", 1);
  fChain->SetBranchStatus("genVertexX", 1);
  fChain->SetBranchStatus("genVertexY", 1);
  fChain->SetBranchStatus("genVertexZ", 1);
  fChain->SetBranchStatus("genVertexT", 1);
  fChain->SetBranchStatus("genWeight", 1);
  fChain->SetBranchStatus("genSignalProcessID", 1);
  fChain->SetBranchStatus("genQScale", 1);
  fChain->SetBranchStatus("genAlphaQCD", 1);
  fChain->SetBranchStatus("genAlphaQED", 1);
  // fChain->SetBranchStatus("lheComments", 1);
  fChain->SetBranchStatus("scaleWeights", 1);
  fChain->SetBranchStatus("pdfWeights", 1);
  fChain->SetBranchStatus("alphasWeights", 1);
}

void RazorAnalyzer::EnableGenParticles() {
  fChain->SetBranchStatus("nGenParticle", 1);
  fChain->SetBranchStatus("gParticleMotherId", 1);
  fChain->SetBranchStatus("gParticleMotherIndex", 1);
  fChain->SetBranchStatus("gParticleId", 1);
  fChain->SetBranchStatus("gParticleStatus", 1);
  fChain->SetBranchStatus("gParticleE", 1);
  fChain->SetBranchStatus("gParticlePt", 1);
  fChain->SetBranchStatus("gParticleEta", 1);
  fChain->SetBranchStatus("gParticlePhi", 1);
  fChain->SetBranchStatus("gParticleProdVertexX", 1);
  fChain->SetBranchStatus("gParticleProdVertexY", 1);
  fChain->SetBranchStatus("gParticleProdVertexZ", 1);
}
void RazorAnalyzer::EnableLLP() {
  fChain->SetBranchStatus("gLLP_eta", 1);
  fChain->SetBranchStatus("gLLP_phi", 1);
  fChain->SetBranchStatus("gLLP_beta", 1);
  fChain->SetBranchStatus("gLLP_e", 1);
  fChain->SetBranchStatus("gLLP_pt", 1);

  fChain->SetBranchStatus("gLLP_decay_vertex_x", 1);
  fChain->SetBranchStatus("gLLP_decay_vertex_y", 1);
  fChain->SetBranchStatus("gLLP_decay_vertex_z", 1);
  fChain->SetBranchStatus("gLLP_daughter_eta", 1);
  fChain->SetBranchStatus("gLLP_daughter_phi", 1);
  fChain->SetBranchStatus("gLLP_daughter_e", 1);
  fChain->SetBranchStatus("gLLP_daughter_pt", 1);
  fChain->SetBranchStatus("gLLP_daughter_mass", 1);
  fChain->SetBranchStatus("gLLP_daughter_id", 1);
}

void RazorAnalyzer::EnableEcalRechits() {
  fChain->SetBranchStatus("ecalRechit_Eta", 1);
  fChain->SetBranchStatus("ecalRechit_Phi", 1);
  fChain->SetBranchStatus("ecalRechit_X", 1);
  fChain->SetBranchStatus("ecalRechit_Y", 1);
  fChain->SetBranchStatus("ecalRechit_Z", 1);
  fChain->SetBranchStatus("ecalRechit_E", 1);
  fChain->SetBranchStatus("ecalRechit_T", 1);
  fChain->SetBranchStatus("ecalRechit_ID", 1);
  fChain->SetBranchStatus("ecalRechit_FlagOOT", 1);
  fChain->SetBranchStatus("ecalRechit_GainSwitch1", 1);
  fChain->SetBranchStatus("ecalRechit_GainSwitch6", 1);
  fChain->SetBranchStatus("ecalRechit_transpCorr", 1);
}

void RazorAnalyzer::EnableTracks() {
  fChain->SetBranchStatus("track_Pt", 1);
  fChain->SetBranchStatus("track_Eta", 1);
  fChain->SetBranchStatus("track_Phi", 1);
  fChain->SetBranchStatus("track_charge", 1);
  fChain->SetBranchStatus("track_bestVertexIndex", 1);
  fChain->SetBranchStatus("track_nMissingInnerHits", 1);
  fChain->SetBranchStatus("track_nMissingOuterHits", 1);
  fChain->SetBranchStatus("track_nPixelHits", 1);
  fChain->SetBranchStatus("track_nHits", 1);
  fChain->SetBranchStatus("track_dxyToBS", 1);
  fChain->SetBranchStatus("track_dxyErr", 1);
  fChain->SetBranchStatus("track_dzToPV", 1);
  fChain->SetBranchStatus("track_dzErr", 1);
  fChain->SetBranchStatus("track_chi2", 1);
  fChain->SetBranchStatus("track_ndof", 1);
}

double RazorAnalyzer::deltaPhi(double phi1, double phi2) {
  double dphi = phi1 - phi2;
  while (dphi > TMath::Pi())
    dphi -= TMath::TwoPi();
  while (dphi <= -TMath::Pi())
    dphi += TMath::TwoPi();
  return dphi;
}

double RazorAnalyzer::deltaR(double eta1, double phi1, double eta2, double phi2) {
  double dphi = deltaPhi(phi1, phi2);
  double deta = eta1 - eta2;
  return sqrt(dphi * dphi + deta * deta);
}

TLorentzVector RazorAnalyzer::makeTLorentzVector(double pt, double eta, double phi, double energy) {
  TLorentzVector vec;
  vec.SetPtEtaPhiE(pt, eta, phi, energy);
  return vec;
}

TLorentzVector RazorAnalyzer::makeTLorentzVectorPtEtaPhiM(double pt, double eta, double phi, double mass) {
  TLorentzVector vec;
  vec.SetPtEtaPhiM(pt, eta, phi, mass);
  return vec;
}
