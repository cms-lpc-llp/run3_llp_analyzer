#ifndef DEF_llp_MuonSystem_CA_TnP
#define DEF_llp_MuonSystem_CA_TnP

#include "RazorAnalyzer.h"

class llp_MuonSystem_CA_TnP : public RazorAnalyzer {
 public:
  llp_MuonSystem_CA_TnP(TTree* tree = 0)
      : RazorAnalyzer(tree) {}
  void Analyze(bool isData, int option, string outputFileName, string label);
};

#endif
