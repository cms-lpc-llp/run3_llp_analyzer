#ifndef DEF_llp_MuonSystem_CA_rechits
#define DEF_llp_MuonSystem_CA_rechits

#include "RazorAnalyzer.h"

class llp_MuonSystem_CA_rechits: public RazorAnalyzer {
    public: 
        llp_MuonSystem_CA_rechits(TTree *tree=0): RazorAnalyzer(tree) { }
        void Analyze(bool isData, int option, string outputFileName, string label);
};

#endif
