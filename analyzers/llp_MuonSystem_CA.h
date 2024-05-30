#ifndef DEF_llp_MuonSystem_CA
#define DEF_llp_MuonSystem_CA

#include "RazorAnalyzer.h"
#include "TreeMuonSystemCombination.h"

class llp_MuonSystem_CA : public RazorAnalyzer
{
public:
    llp_MuonSystem_CA(TTree *tree = 0) : RazorAnalyzer(tree) {}
    void Analyze(bool isData, int option, string outputFileName, string label);
};

class llp_MuonSystem_CAM : public RazorAnalyzerMerged
{
public:
    llp_MuonSystem_CAM(TTree *tree = 0) : RazorAnalyzerMerged(tree) {}
    void Analyze(bool isData, int option, string outputFileName, string label);

private:
    void fillHLT(TreeMuonSystemCombination *MuonSystem);
    void fillMetFilter(TreeMuonSystemCombination *MuonSystem);
};

#endif
