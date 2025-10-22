#ifndef DEF_llp_MuonSystem_CA
#define DEF_llp_MuonSystem_CA

#include "RazorAnalyzer.h"

class llp_MuonSystem_CA : public RazorAnalyzer {
 public:
  llp_MuonSystem_CA(TTree* tree = 0)
      : RazorAnalyzer(tree) {}
  void Analyze(bool isData, int option, string outputFileName, string label);
};

#endif
