#ifndef MUONSYSTEM_SIGNAL_SCAN_MANAGER_H
#define MUONSYSTEM_SIGNAL_SCAN_MANAGER_H

/* #region: includes */
#include <map>
#include <memory>
#include <string>

#include "TFile.h"
#include "TH1F.h"
#include "TTree.h"

#include "MuonSystemPhaseTypes.h"
#include "RazorAnalyzer.h"
#include "TreeMuonSystemCombination.h"
/* #endregion */

/* #region: signal event state */
struct SignalEventState {
  SignalPointKey key = SignalPointKey(0, std::string());
  bool hasKey = false;
};
/* #endregion */

/* #region: signal scan manager declaration */
class MuonSystemSignalScanManager {
 public:
  MuonSystemSignalScanManager(
      bool enabled,
      const std::string& outputFileName,
      TreeMuonSystemCombination* muonSystem);

  SignalEventState beginEvent(RazorAnalyzerMerged& analyzer, float generatorWeight);
  void fillTotal(const SignalEventState& state, float weight);
  void fillAccep(const SignalEventState& state, float weight);
  void fillEventTree(const SignalEventState& state);
  void finalize();
  bool enabled() const { return enabled_; }

 private:
  void ensureSignalOutput(
      const SignalPointKey& key,
      int mh,
      int mx,
      const std::string& ctauTokenTag);

  bool enabled_ = false;
  std::string outputFileName_;
  TreeMuonSystemCombination* muonSystem_ = nullptr;

  std::string cachedSignalSourcePath_;
  int cachedSignalMh_ = 0;
  int cachedSignalMx_ = 0;
  float cachedSignalCtau_ = 0.0f;
  std::string cachedSignalCtauTag_;
  bool hasCachedSignalPoint_ = false;

  std::map<SignalPointKey, std::unique_ptr<TFile>> files2D_;
  std::map<SignalPointKey, TTree*> trees2D_;
  std::map<SignalPointKey, TH1F*> nEvents2D_;
  std::map<SignalPointKey, TH1F*> accep2D_;
  std::map<SignalPointKey, TH1F*> accepMet2D_;
  std::map<SignalPointKey, TH1F*> total2D_;
};
/* #endregion */

#endif
