#ifndef DEF_llp_MuonSystem_SUEP
#define DEF_llp_MuonSystem_SUEP

#include "RazorAnalyzerSUEP.h"

class llp_MuonSystem_SUEP: public RazorAnalyzerSUEP {
    public: 
        llp_MuonSystem_SUEP(TTree *tree=0): RazorAnalyzerSUEP(tree) { }
  void Analyze(bool isData, int option, string outputFileName, string label, int mh, int mx, int ctau, int T);
};

#endif
