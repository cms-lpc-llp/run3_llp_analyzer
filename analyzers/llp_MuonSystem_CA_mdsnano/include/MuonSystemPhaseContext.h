#ifndef MUONSYSTEM_PHASE_CONTEXT_H
#define MUONSYSTEM_PHASE_CONTEXT_H

/* #region: includes */
#include <string>
#include <vector>

#include "MuonSystemPhaseTypes.h"
#include "RazorAnalyzer.h"
#include "RazorHelper.h"
/* #endregion */

/* #region: runtime context types */
struct PhaseRuntimeContext {
  RazorAnalyzerMerged& analyzer;
  RazorHelper& helper;
  TreeMuonSystemCombination& out;
  bool isData;
  int run;
  const std::string& analysisTag;
};

struct EventTransientState {
  EventSynthesis synth;
  std::vector<LeptonCandidate> leptons;
  JetStageResult jetStage;
};
/* #endregion */

#endif
