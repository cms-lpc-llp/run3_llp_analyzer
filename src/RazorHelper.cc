#include "RazorHelper.h"

// Constructor
RazorHelper::RazorHelper(std::string tag_, bool isData_):
        tag(tag_), isData(isData_){
    std::cout << "RazorHelper initializing with tag " << tag << std::endl;


    // check that CMSSW is set up
    loadCMSSWPath();
    if (cmsswPath == "") {
        loadTag_Null();
        return;
    }


    // tag for Run 3
    if (tag == "Summer22")loadTag_Summer22();
    else if (tag == "Summer22EE") loadTag_Summer22EE();
    else if (tag == "Summer23")loadTag_Summer23();
    else if (tag == "Summer23BPix")loadTag_Summer23BPix();
   // tag not found
    else {
        std::cout << "Error in RazorHelper::RazorHelper : specified tag " << tag << " is not supported!" << std::endl;
        loadTag_Null();
        return;
    }

}

// Destructor -- close all of the open TFile objects
RazorHelper::~RazorHelper() {
    if (pileupWeightFile) {
        pileupWeightFile->Close();
        delete pileupWeightFile;
    }

}

// Retrieves CMSSW_BASE and stores in variable cmsswPath
void RazorHelper::loadCMSSWPath() {
    char* cmsswPathChar = getenv("CMSSW_BASE");
    if (cmsswPathChar == NULL) {
        std::cout << "Warning in RazorHelper::loadCMSSWPath : CMSSW_BASE not detected." << std::endl;
        cmsswPath = "";
    }
    cmsswPath = std::string(cmsswPathChar);
}

void RazorHelper::loadTag_Null() {
    std::cout << "Warning: initializing all RazorHelper files and histograms to 0" << std::endl;

    // pileup weights
    pileupWeightFile = 0;
    pileupWeightHist = 0;
    pileupWeightSysUpHist = 0;
    pileupWeightSysDownHist = 0;

}


////////////////////////////////////////////////
//  Summer 22
////////////////////////////////////////////////
void RazorHelper::loadTag_Summer22() {
  loadPileup_Summer22();
}

void RazorHelper::loadPileup_Summer22() {
    // pileup weights
    std::cout << "RazorHelper: loading pileup weight histograms" << std::endl;
    pileupWeightFile = TFile::Open("PileupReweight_Summer22.root");
    pileupWeightHist = (TH1F*)pileupWeightFile->Get("npu_nominal");
    pileupWeightSysUpHist = (TH1F*)pileupWeightFile->Get("npu_up");
    pileupWeightSysDownHist = (TH1F*)pileupWeightFile->Get("npu_down");
}

////////////////////////////////////////////////
//  Summer 22EE
////////////////////////////////////////////////
void RazorHelper::loadTag_Summer22EE() {
  loadPileup_Summer22EE();
}

void RazorHelper::loadPileup_Summer22EE() {
    // pileup weights
    std::cout << "RazorHelper: loading pileup weight histograms" << std::endl;
    pileupWeightFile = TFile::Open("PileupReweight_Summer22EE.root");
    pileupWeightHist = (TH1F*)pileupWeightFile->Get("npu_nominal");
    pileupWeightSysUpHist = (TH1F*)pileupWeightFile->Get("npu_up");
    pileupWeightSysDownHist = (TH1F*)pileupWeightFile->Get("npu_down");
}
////////////////////////////////////////////////
//  Summer 23
////////////////////////////////////////////////
void RazorHelper::loadTag_Summer23() {
  loadPileup_Summer23();
}

void RazorHelper::loadPileup_Summer23() {
    // pileup weights
    std::cout << "RazorHelper: loading pileup weight histograms" << std::endl;
    pileupWeightFile = TFile::Open("PileupReweight_Summer23.root");
    pileupWeightHist = (TH1F*)pileupWeightFile->Get("npu_nominal");
    pileupWeightSysUpHist = (TH1F*)pileupWeightFile->Get("npu_up");
    pileupWeightSysDownHist = (TH1F*)pileupWeightFile->Get("npu_down");
}

////////////////////////////////////////////////
//  Summer 23BPix
////////////////////////////////////////////////
void RazorHelper::loadTag_Summer23BPix() {
  loadPileup_Summer23BPix();
}

void RazorHelper::loadPileup_Summer23BPix() {
    // pileup weights
    std::cout << "RazorHelper: loading pileup weight histograms" << std::endl;
    pileupWeightFile = TFile::Open("PileupReweight_Summer23BPix.root");
    pileupWeightHist = (TH1F*)pileupWeightFile->Get("npu_nominal");
    pileupWeightSysUpHist = (TH1F*)pileupWeightFile->Get("npu_up");
    pileupWeightSysDownHist = (TH1F*)pileupWeightFile->Get("npu_down");
}


////////////////////////////////////////////////
//  Utilities
////////////////////////////////////////////////




double RazorHelper::getPileupWeight(int NPU) {
    if (pileupWeightHist) {
        return pileupWeightHist->GetBinContent(pileupWeightHist->GetXaxis()->FindFixBin(NPU));
    }
    else {
        std::cout << "RazorHelper error: pileup weight requested, but no histogram available!" << std::endl;
        return 0;
    }
}

double RazorHelper::getPileupWeightUp(int NPU) {
    if (pileupWeightSysUpHist) {
        return pileupWeightSysUpHist->GetBinContent(pileupWeightSysUpHist->GetXaxis()->FindFixBin(NPU));
    }
    else {
        std::cout << "RazorHelper error: 'up' pileup weight requested, but no histogram available!" << std::endl;
        return 0;
    }
}

double RazorHelper::getPileupWeightDown(int NPU) {
    if (pileupWeightSysDownHist) {
        return pileupWeightSysDownHist->GetBinContent(pileupWeightSysDownHist->GetXaxis()->FindFixBin(NPU));
    }
    else {
        std::cout << "RazorHelper error: 'down' pileup weight requested, but no histogram available!" << std::endl;
        return 0;
    }
}
