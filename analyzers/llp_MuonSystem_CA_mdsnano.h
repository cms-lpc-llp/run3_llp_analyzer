#ifndef DEF_llp_MuonSystem_CA_mdsnano
#define DEF_llp_MuonSystem_CA_mdsnano

#include "RazorAnalyzer.h"

class llp_MuonSystem_CA_mdsnano : public RazorAnalyzerMerged {
 public:
  llp_MuonSystem_CA_mdsnano(TTree* tree = 0)
      : RazorAnalyzerMerged(tree) {}
  void Analyze(bool isData, int option, string outputFileName, string label);
};

#endif
