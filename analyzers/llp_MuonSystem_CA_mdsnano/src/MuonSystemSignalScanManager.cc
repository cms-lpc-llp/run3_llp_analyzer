/* #region: includes */
#include "MuonSystemSignalScanManager.h"

#include <iostream>
#include <string>

#include "MuonSystemSynthesisPhase.h"
/* #endregion */

/* #region: lifecycle and routing */
MuonSystemSignalScanManager::MuonSystemSignalScanManager(
    bool enabled,
    const std::string& outputFileName,
    TreeMuonSystemCombination* muonSystem)
    : enabled_(enabled),
      outputFileName_(outputFileName),
      muonSystem_(muonSystem) {}

void MuonSystemSignalScanManager::ensureSignalOutput(
    const SignalPointKey& key,
    int mh,
    int mx,
    const std::string& ctauTokenTag) {
  if (files2D_.count(key) != 0)
    return;

  std::string thisFileName = outputFileName_.empty() ? "MuonSystem_Tree.root" : outputFileName_;
  thisFileName.erase(thisFileName.end() - 5, thisFileName.end());
  thisFileName += "_MH-" + std::to_string(mh)
               + "_MS-" + std::to_string(mx)
               + "_ctauS-" + ctauTokenTag + ".root";

  files2D_[key] = std::make_unique<TFile>(thisFileName.c_str(), "recreate");
  files2D_[key]->cd();
  trees2D_[key] = muonSystem_->tree_->CloneTree(0);
  const std::string signalSuffix = std::to_string(mx) + "_" + ctauTokenTag;
  nEvents2D_[key] = new TH1F(("NEvents_" + signalSuffix).c_str(), "NEvents", 1, 0.5, 1.5);
  total2D_[key] = new TH1F(("Total_" + signalSuffix).c_str(), "Total", 1, 0.5, 1.5);
  accep2D_[key] = new TH1F(("accep2D_" + signalSuffix).c_str(), "accep", 1, 0.5, 1.5);
  accepMet2D_[key] = new TH1F(("accep_met2D_" + signalSuffix).c_str(), "accep_met", 1, 0.5, 1.5);

  std::cout << "Created new output file " << thisFileName << std::endl;
}

SignalEventState MuonSystemSignalScanManager::beginEvent(RazorAnalyzerMerged& analyzer, float generatorWeight) {
  SignalEventState state;
  if (!enabled_)
    return state;

  std::string inputPath;
  if (analyzer.fChain != nullptr && analyzer.fChain->GetCurrentFile() != nullptr) {
    inputPath = analyzer.fChain->GetCurrentFile()->GetName();
  }
  if (!hasCachedSignalPoint_ || inputPath != cachedSignalSourcePath_) {
    std::string ctauTokenTag;
    if (!parseSignalPointFromInputPath(inputPath, cachedSignalMh_, cachedSignalMx_, cachedSignalCtau_,
                                       ctauTokenTag)) {
      std::cerr << "[ERROR]: signalScan requested but failed to parse MH/MS/ctauS from input path: '"
                << inputPath << "'. Exiting to avoid mislabeled output." << std::endl;
      return state;
    }
    cachedSignalCtauTag_ = ctauTokenTag;
    cachedSignalSourcePath_ = inputPath;
    hasCachedSignalPoint_ = true;
  }

  muonSystem_->mH = cachedSignalMh_;
  muonSystem_->mX = cachedSignalMx_;
  muonSystem_->ctau = cachedSignalCtau_;

  state.key = std::make_pair(cachedSignalMx_, cachedSignalCtauTag_);
  state.hasKey = true;

  ensureSignalOutput(state.key, cachedSignalMh_, cachedSignalMx_, cachedSignalCtauTag_);
  nEvents2D_[state.key]->Fill(1.0, generatorWeight);
  return state;
}
/* #endregion */

/* #region: per-event fills */
void MuonSystemSignalScanManager::fillTotal(const SignalEventState& state, float weight) {
  if (!enabled_ || !state.hasKey)
    return;
  total2D_[state.key]->Fill(1.0, weight);
}

void MuonSystemSignalScanManager::fillAccep(const SignalEventState& state, float weight) {
  if (!enabled_ || !state.hasKey)
    return;
  accep2D_[state.key]->Fill(1.0, weight);
}

void MuonSystemSignalScanManager::fillEventTree(const SignalEventState& state) {
  if (!enabled_ || !state.hasKey)
    return;
  trees2D_[state.key]->Fill();
}
/* #endregion */

/* #region: finalize outputs */
void MuonSystemSignalScanManager::finalize() {
  if (!enabled_)
    return;

  for (auto& filePtr : files2D_) {
    std::cout << "Writing output tree (" << filePtr.second->GetName() << ")" << std::endl;
    filePtr.second->cd();
    trees2D_[filePtr.first]->Write();
    nEvents2D_[filePtr.first]->Write("NEvents");
    total2D_[filePtr.first]->Write("Total");
    accep2D_[filePtr.first]->Write("accep");
    accepMet2D_[filePtr.first]->Write("accep_met");
    filePtr.second->Close();
  }
}
/* #endregion */
