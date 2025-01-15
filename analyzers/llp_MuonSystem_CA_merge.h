#ifndef DEF_llp_MuonSystem_CA_merge
#define DEF_llp_MuonSystem_CA_merge

#include "RazorAnalyzer.h"

class llp_MuonSystem_CA_merge: public RazorAnalyzer {
    public: 
        llp_MuonSystem_CA_merge(TTree *tree=0): RazorAnalyzer(tree) { }
        void Analyze(bool isData, int option, string outputFileName, string label);
};

#endif
