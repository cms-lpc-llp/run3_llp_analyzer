#ifndef DEF_llp_MuonSystem_CA_mdsnano_hnl
#define DEF_llp_MuonSystem_CA_mdsnano_hnl

#include "RazorAnalyzer.h"

class llp_MuonSystem_CA_mdsnano_hnl: public RazorAnalyzerMerged {
    public: 
        llp_MuonSystem_CA_mdsnano_hnl(TTree *tree=0): RazorAnalyzerMerged(tree) { }
        void Analyze(bool isData, int option, string outputFileName, string label);
};

#endif
