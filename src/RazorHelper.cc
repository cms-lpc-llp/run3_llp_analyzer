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
    else if (tag == "Summer24")loadTag_Summer24();
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


void RazorHelper::loadHMTEfficiency() {
    // pileup weights
    std::cout << "RazorHelper: loading HMT L1 efficiency histograms" << std::endl;
    HMTEffFile = TFile::Open("L1_efficiencies_2022_2023_082324-TEff.root");
    
    HMTEffHist[11] = (TEfficiency*)HMTEffFile->Get("ME11");
    HMTEffHist[12] = (TEfficiency*)HMTEffFile->Get("ME12");
    HMTEffHist[13] = (TEfficiency*)HMTEffFile->Get("ME13");
    HMTEffHist[21] = (TEfficiency*)HMTEffFile->Get("ME21");
    HMTEffHist[22] = (TEfficiency*)HMTEffFile->Get("ME22");
    HMTEffHist[31] = (TEfficiency*)HMTEffFile->Get("ME31");
    HMTEffHist[32] = (TEfficiency*)HMTEffFile->Get("ME32");
    HMTEffHist[41] = (TEfficiency*)HMTEffFile->Get("ME41");
    HMTEffHist[42] = (TEfficiency*)HMTEffFile->Get("ME42");

}


////////////////////////////////////////////////
//  Summer 22
////////////////////////////////////////////////
void RazorHelper::loadTag_Summer22() {
  loadPileup_Summer22();
  loadHMTEfficiency();
  loadJetVeto_Summer22();
}

void RazorHelper::loadPileup_Summer22() {
    // pileup weights
    std::cout << "RazorHelper: loading pileup weight histograms" << std::endl;
    pileupWeightFile = TFile::Open("PileupReweight_Summer22.root");
    pileupWeightHist = (TH1F*)pileupWeightFile->Get("npu_nominal");
    pileupWeightSysUpHist = (TH1F*)pileupWeightFile->Get("npu_up");
    pileupWeightSysDownHist = (TH1F*)pileupWeightFile->Get("npu_down");
}

void RazorHelper::loadJetVeto_Summer22() {
    // pileup weights
    std::cout << "RazorHelper: loading jet veto map histograms" << std::endl;
    JetVetoFile = TFile::Open("Summer22_23Sep2023_RunCD_v1.root");
    if (!JetVetoFile) cout<<"Jet Veto File Not Found"<<endl;
    JetVetoHist = (TH2F*)JetVetoFile->Get("jetvetomap");
}
////////////////////////////////////////////////
//  Summer 22EE
////////////////////////////////////////////////
void RazorHelper::loadTag_Summer22EE() {
  loadPileup_Summer22EE();
   loadHMTEfficiency();
    loadJetVeto_Summer22EE();
}

void RazorHelper::loadPileup_Summer22EE() {
    // pileup weights
    std::cout << "RazorHelper: loading pileup weight histograms" << std::endl;
    pileupWeightFile = TFile::Open("PileupReweight_Summer22EE.root");
    pileupWeightHist = (TH1F*)pileupWeightFile->Get("npu_nominal");
    pileupWeightSysUpHist = (TH1F*)pileupWeightFile->Get("npu_up");
    pileupWeightSysDownHist = (TH1F*)pileupWeightFile->Get("npu_down");
}

void RazorHelper::loadJetVeto_Summer22EE() {
    // pileup weights
    std::cout << "RazorHelper: loading pileup weight histograms" << std::endl;
    JetVetoFile = TFile::Open("Summer22EE_23Sep2023_RunEFG_v1.root");
    JetVetoHist = (TH2F*)JetVetoFile->Get("jetvetomap");
}
////////////////////////////////////////////////
//  Summer 23
////////////////////////////////////////////////
void RazorHelper::loadTag_Summer23() {
  loadPileup_Summer23();
    loadHMTEfficiency();
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
void RazorHelper::loadJetVeto_Summer23() {
    // pileup weights
    std::cout << "RazorHelper: loading pileup weight histograms" << std::endl;
    JetVetoFile = TFile::Open("Summer23Prompt23_RunC_v1.root");
    JetVetoHist = (TH2F*)JetVetoFile->Get("jetvetomap");
}
////////////////////////////////////////////////
//  Summer 23BPix
////////////////////////////////////////////////
void RazorHelper::loadTag_Summer23BPix() {
  loadPileup_Summer23BPix();
    loadHMTEfficiency();
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
void RazorHelper::loadJetVeto_Summer23BPix() {
    // pileup weights
    std::cout << "RazorHelper: loading pileup weight histograms" << std::endl;
    JetVetoFile = TFile::Open("Summer23BPixPrompt23_RunD_v1.root");
    JetVetoHist = (TH2F*)JetVetoFile->Get("jetvetomap");
}
////////////////////////////////////////////////
//  Summer 24
////////////////////////////////////////////////
void RazorHelper::loadTag_Summer24() {
  loadPileup_Summer24();
  loadHMTEfficiency();
  loadJetVeto_Summer24();
}

void RazorHelper::loadPileup_Summer24() {
    // pileup weights
    std::cout << "RazorHelper: loading pileup weight histograms" << std::endl;
    pileupWeightFile = TFile::Open("PileupReweight_Summer23BPix.root");
    pileupWeightHist = (TH1F*)pileupWeightFile->Get("npu_nominal");
    pileupWeightSysUpHist = (TH1F*)pileupWeightFile->Get("npu_up");
    pileupWeightSysDownHist = (TH1F*)pileupWeightFile->Get("npu_down");
}

void RazorHelper::loadJetVeto_Summer24() {
    // pileup weights
    std::cout << "RazorHelper: loading pileup weight histograms" << std::endl;
    JetVetoFile = TFile::Open("Winter24Prompt24_2024BCDEFGHI.root");
    cout<<"here"<<endl;
    if (!JetVetoFile) cout<<"Jet Veto File Not Found"<<endl;
    JetVetoHist = (TH2F*)JetVetoFile->Get("jetvetomap");
    JetVetoFpixHist = (TH2F*)JetVetoFile->Get("jetvetomap_fpix");
}

////////////////////////////////////////////////
//  Utilities
////////////////////////////////////////////////



double RazorHelper::getHMTTriggerEff(int chamber, int nhits){

    map<int, int> hist_cutoff;
    hist_cutoff[11] = 0;
    hist_cutoff[12] = 0;
    hist_cutoff[13] = 600;
    hist_cutoff[21] = 900;
    hist_cutoff[22] = 800;
    hist_cutoff[31] = 900;
    hist_cutoff[32] = 500;
    hist_cutoff[41] = 900;
    hist_cutoff[42] = 500;

    if (HMTEffHist[chamber]){
        return HMTEffHist[chamber]->GetEfficiency(HMTEffHist[chamber]->GetTotalHistogram()->GetXaxis()->FindFixBin(min(nhits,hist_cutoff[chamber])));
    }
    std::cout << "RazorHelper error: HMT efficiency requested, but no histogram available!" << std::endl;
    return 0;
}



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


double RazorHelper::getJetVetoMap(float eta, float phi) {
    if (JetVetoHist) {
        return JetVetoHist->GetBinContent(JetVetoHist->GetXaxis()->FindFixBin(eta), JetVetoHist->GetYaxis()->FindFixBin(phi));
    }
    else {
        std::cout << "RazorHelper error: jet veto map requested, but no histogram available!" << std::endl;
        return 0;
    }
}


double RazorHelper::getJetVetoFpixMap(float eta, float phi) {
    if (JetVetoHist) {
        return JetVetoFpixHist->GetBinContent(JetVetoFpixHist->GetXaxis()->FindFixBin(eta), JetVetoFpixHist->GetYaxis()->FindFixBin(phi));
    }
    else {
        std::cout << "RazorHelper error: jet veto map requested, but no histogram available!" << std::endl;
        return 0;
    }
}