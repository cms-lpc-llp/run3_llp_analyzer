// Class to manage files for b-tag scale factors, lepton scale factors, pileup weights, and other information

#ifndef RazorHelper_H
#define RazorHelper_H

#include <iostream>
#include <string>
#include <sys/stat.h>

#include "TFile.h"
#include "TH1F.h"
#include "TH2D.h"
#include "TRandom.h"

#include "RazorAnalyzer.h"

class RazorHelper
{

public:
    // constructor takes a string specifying which set of files to load.
    RazorHelper(std::string tag_, bool isData_);
    virtual ~RazorHelper();

    std::pair<double, double> METXYCorr_Met_MetPhi(double uncormet, double uncormet_phi, int runnb, int year, bool isMC, int npv);

    double getMetTriggerSF(float met);
    // retrieve pileup weights (nominal, up, and down versions)
    double getPileupWeight(int NPU);
    double getPileupWeightUp(int NPU);
    double getPileupWeightDown(int NPU);

private:
    // member functions

    void loadTag_Razor2018_17SeptEarlyReReco();
    void loadTag_Null(); // Default when tag is not provided
    void loadCMSSWPath();

    // for Razor Razor2018
    void loadPileup_Razor2018_17SeptEarlyReReco();
    void loadTrigger_Razor2018_17SeptEarlyReReco();
    void loadJECs_Razor2018_17SeptEarlyReReco();

    // member data
    std::string tag;
    bool isData;
    std::string cmsswPath;

    TFile *metTriggerSFFile;
    TH1F *metTriggerSFHist;

    // for pileup reweighting
    TFile *pileupWeightFile;
    TH1F *pileupWeightHist;
    TH1F *pileupWeightSysUpHist;
    TH1F *pileupWeightSysDownHist;
};

#endif
