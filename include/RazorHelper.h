// Class to manage files for b-tag scale factors, lepton scale factors, pileup weights, and other information

#ifndef RazorHelper_H
#define RazorHelper_H

#include <iostream>
#include <string>
#include <sys/stat.h>

#include "TFile.h"
#include "TH1F.h"
#include "TEfficiency.h"
#include "TH2D.h"
#include "TRandom.h"


#include "JetCorrectionUncertainty.h"



#include "RazorAnalyzer.h"

class RazorHelper {

    public:
        // constructor takes a string specifying which set of files to load.
        RazorHelper(std::string tag_, bool isData_);
        virtual ~RazorHelper();

        std::pair<double,double> METXYCorr_Met_MetPhi(double uncormet, double uncormet_phi, int runnb, int year, bool isMC, int npv);



        double getHMTTriggerEff(int chamber, int nhits);
        // retrieve pileup weights (nominal, up, and down versions)
        double getPileupWeight(int NPU);
        double getPileupWeightUp(int NPU);
        double getPileupWeightDown(int NPU);
        double getMetTriggerEff(float met);
        double getMetTriggerEffUp(float met);
        double getMetTriggerEffDown(float met);
        double getJetVetoMap(float eta, float phi);
        double getJetVetoFpixMap(float eta, float phi);
        bool jetTightLepVeto(std::string tag, bool tightVeto, float Jet_neHEF, float Jet_neEmEF, float Jet_chEmEF, float Jet_muEF, float Jet_chHEF, UChar_t Jet_chMultiplicity,UChar_t Jet_neMultiplicity, float Jet_eta, bool Jet_jetId);
        double getJecUnc( float pt, float eta , int run);


    private:
        // member functions

        void loadTag_Summer22();
        void loadTag_Summer22EE();
        void loadTag_Summer23();
        void loadTag_Summer23BPix();
        void loadTag_Summer24();

        void loadTag_Null(); // Default when tag is not provided
        void loadCMSSWPath();


        void loadPileup_Summer22();
        void loadPileup_Summer22EE();
        void loadPileup_Summer23();
        void loadPileup_Summer23BPix();
        void loadPileup_Summer24();
        void loadHMTEfficiency2223();
        void loadHMTEfficiency24();

        void loadJetVeto_Summer22();
        void loadJetVeto_Summer22EE();
        void loadJetVeto_Summer23();
        void loadJetVeto_Summer23BPix();
        void loadJetVeto_Summer24();

        void loadMetTrigger_Summer24();
        void loadMetTrigger_Summer22();
        void loadMetTrigger_Summer22EE();
        void loadMetTrigger_Summer23();
        void loadMetTrigger_Summer23BPix();
        
        void loadJECs();

        //for Razor Razor2018
        void loadPileup_Razor2018_17SeptEarlyReReco();
        void loadTrigger_Razor2018_17SeptEarlyReReco();
        void loadJECs_Razor2018_17SeptEarlyReReco();

        // member data
        std::string tag;
        bool isData;
        std::string cmsswPath;



        // for pileup reweighting
        TFile *pileupWeightFile;
        TH1F *pileupWeightHist;
        TH1F *pileupWeightSysUpHist;
        TH1F *pileupWeightSysDownHist;

        TFile *JetVetoFile;
        TH2D *JetVetoHist;
        TH2D *JetVetoFpixHist;


        TFile *MetTriggerFile;
        TH1F *MetTriggerHist;
        TH1F *MetTriggerSysUpHist;
        TH1F *MetTriggerSysDownHist;



        std::vector<JetCorrectionUncertainty*> jecUnc;
	    std::vector<std::pair<int,int> > JetCorrectionsIOV;


        TFile *HMTEffFile;

        map<int, TEfficiency*> HMTEffHist;


};

#endif
