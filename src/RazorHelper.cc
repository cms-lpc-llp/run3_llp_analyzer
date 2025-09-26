#include "RazorHelper.h"

// Constructor
RazorHelper::RazorHelper(std::string tag_, bool isData_)
    : tag(tag_), isData(isData_) {
  std::cout << "RazorHelper initializing with tag " << tag << std::endl;

  // check that CMSSW is set up
  loadCMSSWPath();
  if (cmsswPath == "") {
    loadTag_Null();
    return;
  }

  // tag for Run 3
  if (tag == "Summer22")
    loadTag_Summer22();
  else if (tag == "Summer22EE")
    loadTag_Summer22EE();
  else if (tag == "Summer23")
    loadTag_Summer23();
  else if (tag == "Summer23BPix")
    loadTag_Summer23BPix();
  else if (tag == "Summer24")
    loadTag_Summer24();
  // tag not found
  else {
    std::cout << "Error in RazorHelper::RazorHelper : specified tag " << tag
              << " is not supported!" << std::endl;
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
  char *cmsswPathChar = getenv("CMSSW_BASE");
  if (cmsswPathChar == NULL) {
    std::cout
        << "Warning in RazorHelper::loadCMSSWPath : CMSSW_BASE not detected."
        << std::endl;
    cmsswPath = "";
  }
  cmsswPath = std::string(cmsswPathChar);
}

void RazorHelper::loadTag_Null() {
  std::cout << "Warning: initializing all RazorHelper files and histograms to 0"
            << std::endl;

  // pileup weights
  pileupWeightFile = 0;
  pileupWeightHist = 0;
  pileupWeightSysUpHist = 0;
  pileupWeightSysDownHist = 0;
}

void RazorHelper::loadHMTEfficiency2223() {
  // pileup weights
  std::cout << "RazorHelper: loading HMT L1 efficiency histograms" << std::endl;
  HMTEffFile = TFile::Open("L1_efficiencies_2022_2023_032625-Hists-TEff.root");

  HMTEffHist[11] = (TEfficiency *)HMTEffFile->Get("ME11");
  HMTEffHist[12] = (TEfficiency *)HMTEffFile->Get("ME12");
  HMTEffHist[13] = (TEfficiency *)HMTEffFile->Get("ME13");
  HMTEffHist[21] = (TEfficiency *)HMTEffFile->Get("ME21");
  HMTEffHist[22] = (TEfficiency *)HMTEffFile->Get("ME22");
  HMTEffHist[31] = (TEfficiency *)HMTEffFile->Get("ME31");
  HMTEffHist[32] = (TEfficiency *)HMTEffFile->Get("ME32");
  HMTEffHist[41] = (TEfficiency *)HMTEffFile->Get("ME41");
  HMTEffHist[42] = (TEfficiency *)HMTEffFile->Get("ME42");
}
void RazorHelper::loadHMTEfficiency24() {
  // pileup weights
  std::cout << "RazorHelper: loading HMT L1 efficiency histograms" << std::endl;
  HMTEffFile = TFile::Open("HMT_Efficiencies_2024.root");

  HMTEffHist[11] = (TEfficiency *)HMTEffFile->Get("ME11");
  HMTEffHist[12] = (TEfficiency *)HMTEffFile->Get("ME12");
  HMTEffHist[13] = (TEfficiency *)HMTEffFile->Get("ME13");
  HMTEffHist[21] = (TEfficiency *)HMTEffFile->Get("ME21");
  HMTEffHist[22] = (TEfficiency *)HMTEffFile->Get("ME22");
  HMTEffHist[31] = (TEfficiency *)HMTEffFile->Get("ME31");
  HMTEffHist[32] = (TEfficiency *)HMTEffFile->Get("ME32");
  HMTEffHist[41] = (TEfficiency *)HMTEffFile->Get("ME41");
  HMTEffHist[42] = (TEfficiency *)HMTEffFile->Get("ME42");
}

////////////////////////////////////////////////
//  Summer 22
////////////////////////////////////////////////
void RazorHelper::loadTag_Summer22() {
  loadPileup_Summer22();
  loadHMTEfficiency2223();
  loadJetVeto_Summer22();
  loadJECs();
  loadMetTrigger_Summer22();
}

void RazorHelper::loadMetTrigger_Summer22() {
  std::cout << "RazorHelper: loading met trigger histograms" << std::endl;
  MetTriggerFile = TFile::Open("METTriggerEff_Summer22.root");
  MetTriggerHist = (TH1F *)MetTriggerFile->Get("eff_exp");
  MetTriggerSysUpHist = (TH1F *)MetTriggerFile->Get("eff_high");
  MetTriggerSysDownHist = (TH1F *)MetTriggerFile->Get("eff_low");
}

void RazorHelper::loadPileup_Summer22() {
  // pileup weights
  std::cout << "RazorHelper: loading pileup weight histograms" << std::endl;
  pileupWeightFile = TFile::Open("PileupReweight_Summer22.root");
  pileupWeightHist = (TH1F *)pileupWeightFile->Get("npu_nominal");
  pileupWeightSysUpHist = (TH1F *)pileupWeightFile->Get("npu_up");
  pileupWeightSysDownHist = (TH1F *)pileupWeightFile->Get("npu_down");
}

void RazorHelper::loadJetVeto_Summer22() {
  // pileup weights
  std::cout << "RazorHelper: loading jet veto map histograms" << std::endl;
  JetVetoFile = TFile::Open("Summer22_23Sep2023_RunCD_v1.root");
  if (!JetVetoFile)
    cout << "Jet Veto File Not Found" << endl;
  JetVetoHist = (TH2D *)JetVetoFile->Get("jetvetomap");
}
////////////////////////////////////////////////
//  Summer 22EE
////////////////////////////////////////////////
void RazorHelper::loadTag_Summer22EE() {
  loadPileup_Summer22EE();
  loadHMTEfficiency2223();
  loadJetVeto_Summer22EE();
  loadJECs();
  loadMetTrigger_Summer22EE();
}

void RazorHelper::loadMetTrigger_Summer22EE() {
  std::cout << "RazorHelper: loading met trigger histograms" << std::endl;
  MetTriggerFile = TFile::Open("METTriggerEff_Summer22EE.root");
  MetTriggerHist = (TH1F *)MetTriggerFile->Get("eff_exp");
  MetTriggerSysUpHist = (TH1F *)MetTriggerFile->Get("eff_high");
  MetTriggerSysDownHist = (TH1F *)MetTriggerFile->Get("eff_low");
}
void RazorHelper::loadPileup_Summer22EE() {
  // pileup weights
  std::cout << "RazorHelper: loading pileup weight histograms" << std::endl;
  pileupWeightFile = TFile::Open("PileupReweight_Summer22EE.root");
  pileupWeightHist = (TH1F *)pileupWeightFile->Get("npu_nominal");
  pileupWeightSysUpHist = (TH1F *)pileupWeightFile->Get("npu_up");
  pileupWeightSysDownHist = (TH1F *)pileupWeightFile->Get("npu_down");
}

void RazorHelper::loadJetVeto_Summer22EE() {
  // pileup weights
  std::cout << "RazorHelper: loading pileup weight histograms" << std::endl;
  JetVetoFile = TFile::Open("Summer22EE_23Sep2023_RunEFG_v1.root");
  JetVetoHist = (TH2D *)JetVetoFile->Get("jetvetomap");
}
////////////////////////////////////////////////
//  Summer 23
////////////////////////////////////////////////
void RazorHelper::loadTag_Summer23() {
  loadPileup_Summer23();
  loadHMTEfficiency2223();
  loadJetVeto_Summer23();
  loadJECs();
  loadMetTrigger_Summer23();
}

void RazorHelper::loadMetTrigger_Summer23() {
  std::cout << "RazorHelper: loading met trigger histograms" << std::endl;
  MetTriggerFile = TFile::Open("METTriggerEff_Summer23.root");
  MetTriggerHist = (TH1F *)MetTriggerFile->Get("eff_exp");
  MetTriggerSysUpHist = (TH1F *)MetTriggerFile->Get("eff_high");
  MetTriggerSysDownHist = (TH1F *)MetTriggerFile->Get("eff_low");
}
void RazorHelper::loadPileup_Summer23() {
  // pileup weights
  std::cout << "RazorHelper: loading pileup weight histograms" << std::endl;
  pileupWeightFile = TFile::Open("PileupReweight_Summer23.root");
  pileupWeightHist = (TH1F *)pileupWeightFile->Get("npu_nominal");
  pileupWeightSysUpHist = (TH1F *)pileupWeightFile->Get("npu_up");
  pileupWeightSysDownHist = (TH1F *)pileupWeightFile->Get("npu_down");
}
void RazorHelper::loadJetVeto_Summer23() {
  // pileup weights
  std::cout << "RazorHelper: loading pileup weight histograms" << std::endl;
  JetVetoFile = TFile::Open("Summer23Prompt23_RunC_v1.root");
  JetVetoHist = (TH2D *)JetVetoFile->Get("jetvetomap");
}
////////////////////////////////////////////////
//  Summer 23BPix
////////////////////////////////////////////////
void RazorHelper::loadTag_Summer23BPix() {
  loadPileup_Summer23BPix();
  loadHMTEfficiency2223();
  loadJetVeto_Summer23BPix();
  loadJECs();
  loadMetTrigger_Summer23BPix();
}

void RazorHelper::loadMetTrigger_Summer23BPix() {
  std::cout << "RazorHelper: loading met trigger histograms" << std::endl;
  MetTriggerFile = TFile::Open("METTriggerEff_Summer23BPix.root");
  MetTriggerHist = (TH1F *)MetTriggerFile->Get("eff_exp");
  MetTriggerSysUpHist = (TH1F *)MetTriggerFile->Get("eff_high");
  MetTriggerSysDownHist = (TH1F *)MetTriggerFile->Get("eff_low");
}
void RazorHelper::loadPileup_Summer23BPix() {
  // pileup weights
  std::cout << "RazorHelper: loading pileup weight histograms" << std::endl;
  pileupWeightFile = TFile::Open("PileupReweight_Summer23BPix.root");
  pileupWeightHist = (TH1F *)pileupWeightFile->Get("npu_nominal");
  pileupWeightSysUpHist = (TH1F *)pileupWeightFile->Get("npu_up");
  pileupWeightSysDownHist = (TH1F *)pileupWeightFile->Get("npu_down");
}
void RazorHelper::loadJetVeto_Summer23BPix() {
  // pileup weights
  std::cout << "RazorHelper: loading pileup weight histograms" << std::endl;
  JetVetoFile = TFile::Open("Summer23BPixPrompt23_RunD_v1.root");
  JetVetoHist = (TH2D *)JetVetoFile->Get("jetvetomap");
}
////////////////////////////////////////////////
//  Summer 24
////////////////////////////////////////////////
void RazorHelper::loadTag_Summer24() {
  loadPileup_Summer24();
  loadHMTEfficiency24();
  loadJetVeto_Summer24();
  loadJECs();
  loadMetTrigger_Summer24();
}

void RazorHelper::loadMetTrigger_Summer24() {
  std::cout << "RazorHelper: loading met trigger histograms" << std::endl;
  MetTriggerFile = TFile::Open("METTriggerEff_Summer24.root");
  MetTriggerHist = (TH1F *)MetTriggerFile->Get("eff_exp");
  MetTriggerSysUpHist = (TH1F *)MetTriggerFile->Get("eff_high");
  MetTriggerSysDownHist = (TH1F *)MetTriggerFile->Get("eff_low");
}
void RazorHelper::loadPileup_Summer24() {
  // pileup weights
  std::cout << "RazorHelper: loading pileup weight histograms" << std::endl;
  pileupWeightFile = TFile::Open("PileupReweight_Summer24.root");
  pileupWeightHist = (TH1F *)pileupWeightFile->Get("npu_nominal");
  pileupWeightSysUpHist = (TH1F *)pileupWeightFile->Get("npu_up");
  pileupWeightSysDownHist = (TH1F *)pileupWeightFile->Get("npu_down");
}

void RazorHelper::loadJetVeto_Summer24() {
  // pileup weights
  std::cout << "RazorHelper: loading pileup weight histograms" << std::endl;
  JetVetoFile = TFile::Open("Winter24Prompt24_2024BCDEFGHI.root");
  if (!JetVetoFile)
    cout << "Jet Veto File Not Found" << endl;
  JetVetoHist = (TH2D *)JetVetoFile->Get("jetvetomap");
  JetVetoFpixHist = (TH2D *)JetVetoFile->Get("jetvetomap_fpix");
}

////////////////////////////////////////////////
//  Utilities
////////////////////////////////////////////////

double RazorHelper::getMetTriggerEff(float met) {
  if (MetTriggerHist) {
    return MetTriggerHist->GetBinContent(
        MetTriggerHist->GetXaxis()->FindFixBin(met));
  } else {
    std::cout << "RazorHelper error: MET trigger eff requested, but no "
                 "histogram available!"
              << std::endl;
    return 0;
  }
}

double RazorHelper::getMetTriggerEffUp(float met) {
  if (MetTriggerSysUpHist) {
    return MetTriggerSysUpHist->GetBinContent(
        MetTriggerSysUpHist->GetXaxis()->FindFixBin(met));
  } else {
    std::cout << "RazorHelper error: MET trigger eff requested, but no "
                 "histogram available!"
              << std::endl;
    return 0;
  }
}

double RazorHelper::getMetTriggerEffDown(float met) {
  if (MetTriggerSysDownHist) {
    return MetTriggerSysDownHist->GetBinContent(
        MetTriggerSysDownHist->GetXaxis()->FindFixBin(met));
  } else {
    std::cout << "RazorHelper error: MET trigger eff requested, but no "
                 "histogram available!"
              << std::endl;
    return 0;
  }
}

//// GET JET UNC for 22-24
double RazorHelper::getJecUnc(float pt, float eta, int run) {

  int foundIndex = -1;
  for (uint i = 0; i < JetCorrectionsIOV.size(); i++) {
    if (run >= JetCorrectionsIOV[i].first &&
        run <= JetCorrectionsIOV[i].second) {
      foundIndex = i;
    }
  }
  if (foundIndex == -1) {
    std::cout << "Warning: run = " << run
              << " was not found in any valid IOV range. use default index = 0 "
                 "for Jet energy corrections. \n";
    foundIndex = 0;
  }

  jecUnc[foundIndex]->setJetPt(pt);
  jecUnc[foundIndex]->setJetEta(eta);
  return jecUnc[foundIndex]->getUncertainty(true);
}

void RazorHelper::loadJECs() {
  std::cout << "RazorHelper: loading jet energy correction constants, using "
               "2018_17SeptEarlyReReco"
            << std::endl;
  // initialize
  std::string jecPathname = "JEC/";

  jecUnc = std::vector<JetCorrectionUncertainty *>();
  JetCorrectionsIOV = std::vector<std::pair<int, int>>();
  if (isData) {
    std::string jecUncPathA =
        jecPathname +
        "/Summer22_22Sep2023_V2_MC_Uncertainty_AK4PFPuppi"; // place holder for
                                                            // now, since we
                                                            // dont evaluate on
                                                            // data
    JetCorrectionUncertainty *jecUncA =
        new JetCorrectionUncertainty(jecUncPathA);
    jecUnc.push_back(jecUncA);
    JetCorrectionsIOV.push_back(std::pair<int, int>(352319, 387121));
  } else {
    std::cout << "Loading Jet Energy Corrections: 22-24 \n";
    std::string jecUncPath =
        jecPathname +
        "/Summer22_22Sep2023_V2_MC_Uncertainty_AK4PFPuppi.txt"; // same file
                                                                // works for
                                                                // 22-24
                                                                // 2025/04/25
    JetCorrectionUncertainty *jecUncMC =
        new JetCorrectionUncertainty(jecUncPath);
    jecUnc.push_back(jecUncMC);
    JetCorrectionsIOV.push_back(std::pair<int, int>(-1, 99999999));
  }
}

// implemented based on these:
// https://gitlab.cern.ch/cms-jetmet/coordination/coordination/-/issues/117
bool RazorHelper::jetTightLepVeto(std::string tag, bool tightVeto,
                                  float Jet_neHEF, float Jet_neEmEF,
                                  float Jet_chEmEF, float Jet_muEF,
                                  float Jet_chHEF, UChar_t Jet_chMultiplicity,
                                  UChar_t Jet_neMultiplicity, float Jet_eta,
                                  bool Jet_jetId) {
  bool Jet_passJetIdTightLepVeto = false;
  bool Jet_passJetIdTight = false;

  if (tag == "Summer24") {
    if (abs(Jet_eta) <= 2.6)
      Jet_passJetIdTight = (Jet_neHEF < 0.99) && (Jet_neEmEF < 0.9) &&
                           (Jet_chMultiplicity + Jet_neMultiplicity > 1) &&
                           (Jet_chHEF > 0.01) && (Jet_chMultiplicity > 0);
    else if (abs(Jet_eta) > 2.6 && abs(Jet_eta) <= 2.7)
      Jet_passJetIdTight = (Jet_neHEF < 0.90) && (Jet_neEmEF < 0.99);
    else if (abs(Jet_eta) > 2.7 && abs(Jet_eta) <= 3.0)
      Jet_passJetIdTight = (Jet_neHEF < 0.99);
    else if (abs(Jet_eta) > 3.0)
      Jet_passJetIdTight = (Jet_neMultiplicity >= 2) && (Jet_neEmEF < 0.4);

    if (abs(Jet_eta) <= 2.7)
      Jet_passJetIdTightLepVeto =
          Jet_passJetIdTight && (Jet_muEF < 0.8) && (Jet_chEmEF < 0.8);
    else
      Jet_passJetIdTightLepVeto = Jet_passJetIdTight;
  } else {
    if (abs(Jet_eta) <= 2.7)
      Jet_passJetIdTight = Jet_jetId & (1 << 1);
    else if (abs(Jet_eta) > 2.7 && abs(Jet_eta) <= 3.0)
      Jet_passJetIdTight = (Jet_jetId & (1 << 1)) && (Jet_neHEF < 0.99);
    else if (abs(Jet_eta) > 3.0)
      Jet_passJetIdTight = (Jet_jetId & (1 << 1)) && (Jet_neEmEF < 0.4);

    if (abs(Jet_eta) <= 2.7)
      Jet_passJetIdTightLepVeto =
          Jet_passJetIdTight && (Jet_muEF < 0.8) && (Jet_chEmEF < 0.8);
    else
      Jet_passJetIdTightLepVeto = Jet_passJetIdTight;
  }
  if (tightVeto)
    return Jet_passJetIdTightLepVeto;
  else
    return Jet_passJetIdTight;
}

double RazorHelper::getHMTTriggerEff(int chamber, int nhits) {

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

  if (HMTEffHist[chamber]) {
    return HMTEffHist[chamber]->GetEfficiency(
        HMTEffHist[chamber]->GetTotalHistogram()->GetXaxis()->FindFixBin(
            min(nhits, hist_cutoff[chamber])));
  }
  std::cout << "RazorHelper error: HMT efficiency requested, but no histogram "
               "available!"
            << std::endl;
  return 0;
}

double RazorHelper::getPileupWeight(int NPU) {
  if (pileupWeightHist) {
    return pileupWeightHist->GetBinContent(
        pileupWeightHist->GetXaxis()->FindFixBin(NPU));
  } else {
    std::cout << "RazorHelper error: pileup weight requested, but no histogram "
                 "available!"
              << std::endl;
    return 0;
  }
}

double RazorHelper::getPileupWeightUp(int NPU) {
  if (pileupWeightSysUpHist) {
    return pileupWeightSysUpHist->GetBinContent(
        pileupWeightSysUpHist->GetXaxis()->FindFixBin(NPU));
  } else {
    std::cout << "RazorHelper error: 'up' pileup weight requested, but no "
                 "histogram available!"
              << std::endl;
    return 0;
  }
}

double RazorHelper::getPileupWeightDown(int NPU) {
  if (pileupWeightSysDownHist) {
    return pileupWeightSysDownHist->GetBinContent(
        pileupWeightSysDownHist->GetXaxis()->FindFixBin(NPU));
  } else {
    std::cout << "RazorHelper error: 'down' pileup weight requested, but no "
                 "histogram available!"
              << std::endl;
    return 0;
  }
}

double RazorHelper::getJetVetoMap(float eta, float phi) {
  if (JetVetoHist) {
    return JetVetoHist->GetBinContent(JetVetoHist->GetXaxis()->FindFixBin(eta),
                                      JetVetoHist->GetYaxis()->FindFixBin(phi));
  } else {
    std::cout << "RazorHelper error: jet veto map requested, but no histogram "
                 "available!"
              << std::endl;
    return 0;
  }
}

double RazorHelper::getJetVetoFpixMap(float eta, float phi) {
  if (JetVetoFpixHist) {
    return JetVetoFpixHist->GetBinContent(
        JetVetoFpixHist->GetXaxis()->FindFixBin(eta),
        JetVetoFpixHist->GetYaxis()->FindFixBin(phi));
  } else {
    std::cout << "RazorHelper error: jet veto map requested, but no histogram "
                 "available!"
              << std::endl;
    return 0;
  }
}
