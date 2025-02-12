#include <fstream>
#include <sstream>
#include <iterator>
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TKey.h"
#include <assert.h>
#include <TRandom3.h>
#include "TTreeFormula.h"
#include <iostream>
#include "RazorAnalyzer.h"
#include "HNLMuonSystemTree.h"
#include "RazorHelper.h"
#include <vector>
#include "TString.h"
using namespace std;

std::string ParseCommandLine(int argc, char *argv[], std::string opt)
{
  for (int i = 1; i < argc; i++)
  {
    std::string tmp(argv[i]);
    if (tmp.find(opt) != std::string::npos)
    {
      if (tmp.find("=") != std::string::npos)
        return tmp.substr(tmp.find_last_of("=") + 1);
      if (tmp.find("--") != std::string::npos)
        return "yes";
    }
  }

  return "";
};

struct leptons
{
  TLorentzVector lepton;
  int pdgId;
  float dZ;
  // bool passLooseId;
  // bool passMediumId;
  bool passId;
  bool passVetoId;
  bool passLooseIso;
  bool passTightIso;
  bool passVTightIso;
  bool passVVTightIso;
};

struct jets
{
  TLorentzVector jet;
  float time;
  bool passId;
  // bool passLooseId;
  // bool passMediumId;
  // bool passTightId;
  bool isCSVL;
  int ecalNRechits;
  float jetCISV;
  float jetCMVA;
  float ecalRechitE;
  float jetChargedEMEnergyFraction;
  float jetNeutralEMEnergyFraction;
  float jetChargedHadronEnergyFraction;
  float jetNeutralHadronEnergyFraction;
  bool jetPassMuFrac;
  float jetPtJESUp;
  float jetPtJESDown;
  float jetEJESUp;
  float jetEJESDown;
  float JecUnc;

  float electronEnergyFraction;
  float neutralEmEnergyFraction;
  float chargedHadronEnergyFraction;
  float neutralHadronEnergyFraction;
  float muonEnergyFraction;
};

// lepton highest pt comparator
struct largest_pt
{
  inline bool operator()(const leptons &p1, const leptons &p2) { return p1.lepton.Pt() > p2.lepton.Pt(); }
} my_largest_pt;

// jet highest pt comparator
struct largest_pt_jet
{
  inline bool operator()(const jets &p1, const jets &p2) { return p1.jet.Pt() > p2.jet.Pt(); }
} my_largest_pt_jet;

void setBranchValues(int numMuons, Float_t razorVar[], std::vector<Float_t> inputVar)
{
  inputVar.clear();
  for (int i = 0; i < numMuons; i++)
  {
    inputVar.push_back(razorVar[i]);
    // cout << razorVar[i] << "\n";
  }
}
// get list of files to open, add normalization branch to the tree in each file
int main(int argc, char *argv[])
{
  // parse input list to get names of ROOT files
  if (argc < 4)
  {
    cerr << "usage MergeNtuple [inputfile1] [inputfile2] [outputfile] [isData]" << endl;
    return -1;
  }
  string filenameNTuplers(argv[1]);
  string filenameNanoAODList(argv[2]);
  string outputfilename(argv[3]);
  string analysisTag(argv[4]); // analysisTag

  const float ELE_MASS = 0.000511;
  const float MU_MASS = 0.105658;
  const float Z_MASS = 91.2;

  if (analysisTag == "")
  {
    analysisTag = "Razor2016_80X";
  }

  // create output file
  TFile *outputFile = new TFile(outputfilename.c_str(), "RECREATE");

  // open files
  // TFile *inputFile1 = TFile::Open(filenameNTuplers.c_str(), "READ");
  HNLMuonSystemTree *MuonSystem = new HNLMuonSystemTree;
  // llp_event *MuonSystem = new llp_event;
  MuonSystem->LoadTree(filenameNTuplers.c_str());
  if (!MuonSystem) {
    cout << "Input Tree not found in file " << filenameNTuplers << "\n";
  }
  assert(MuonSystem);
  // exit()

  // MuonSystem->tree_->Print();
  // TFile *inputFile2 = TFile::Open(filenameNanoAODList.c_str(), "READ");
  // TTree *inputTree2 = (TTree*)inputFile2->Get("tree");
  // if (!inputTree2) cout << "Input Tree not found in file " << filenameNanoAODList << "\n";
  // assert(inputTree2);

  // build the TChain
  // tree name is set give the structure in the first root file, see while loop below
  TChain *theChain = new TChain();
  string curFileName;
  cout << "filenameNanoAODList.c_str() = " << filenameNanoAODList.c_str() << endl;
  ifstream inputFile(filenameNanoAODList.c_str());
  int NFilesLoaded = 0;
  if (!inputFile)
  {
    cerr << "Error: input file not found!" << endl;
    return -1;
  }

  while (getline(inputFile, curFileName))
  {
    cout << "curFileName = " << curFileName << endl;
    if (NFilesLoaded == 0)
    {
      /*
        checks root file structure and add first file
      */
      std::cout << "[INFO]: loading file: " << curFileName.c_str() << std::endl;
      TFile *f_0 = TFile::Open(curFileName.c_str());
      f_0->ls();
      theChain->SetName("Events");
      std::cout << "[INFO]: setting name of tchain = Events" << std::endl;
      theChain->Add(curFileName.c_str());
      std::cout << "Loaded  " << theChain->GetEntries() << " events\n";
      delete f_0;
    }
    else
    {
      // Addind remaining files after file structure is decided
      theChain->Add(curFileName.c_str());
    }
    NFilesLoaded++;
  }
  std::cout << "Loaded Total of " << NFilesLoaded << " files\n";
  std::cout << "Loaded Total of " << theChain->GetEntries() << " events\n";
  if (theChain == NULL) {
    return -1;
  }

  uint NEventsTree1 = MuonSystem->tree_->GetEntries();
  cout << "Num nTuple Events:" << NEventsTree1 << endl;
  uint NEventsTree2 = theChain->GetEntries();
  cout << "Num AOD Events:" << NEventsTree2 << endl;
  //*****************************************************************************************
  // Make map of event number in tree 1 to event number in tree 2
  //*****************************************************************************************
  int numberOfNonMatchedEvents = 0;
  std::map<uint, uint> EventIndexToEventIndexMap;

  // Read AOD input with analyzer
  RazorAnalyzer ra = RazorAnalyzer(theChain);
  ra.EnableAll();

  ra.fChain->SetBranchStatus("event", 1);
  ra.fChain->SetBranchAddress("event", &ra.evtNumLong);

  ra.fChain->SetBranchStatus("run", 1);
  ra.fChain->SetBranchAddress("run", &ra.runNum);

  // cout << ra.muonPt->
  // loop over tree2
  std::vector<std::pair<uint, uint>> eventList2;
  std::vector<bool> matchedevent;

  cout << "building tree 2 index map (razor) \n";
  for (UInt_t m = 0; m < NEventsTree2; m++)
  {
    if (m % 10000 == 0)
      cout << "Event " << m << "\n";
    ra.fChain->GetEntry(m);
    std::pair<UInt_t, ULong64_t> p(ra.runNum, ra.evtNumLong);
    // cout << " runNum = " <<p.first<<" LS ="<<ra.lumiNum<< " evtNum = "<< p.second << "  \n";
    eventList2.push_back(p);
  }

  cout << "Looping over tree 1 (MuonSystem):  \n";
  cout << "Total Input Entries: " << NEventsTree1 << "\n";
  // loop over tree1
  for (uint n = 0; n < NEventsTree1; n++)
  {
    MuonSystem->tree_->GetEntry(n);
    if (n % 10000 == 0) {
      cout << "Event " << n << "\n";
      // cout << "MuonSystem[" << n << "].runNum = " << MuonSystem->runNum << " MuonSystem[" << n << "].evtNumLong = " << MuonSystem->evtNumLong << " MuonSystem[" << n << "].lumiNum = " << MuonSystem->lumiNum << "\n";
    }
    bool matchFound = false;
    for (uint m = 0; m < eventList2.size(); m++)
    {
      // if (n % 10000 == 0 && m % 2000 == 0) {
      //   // cout << "n = " << n << " m = " << m << "\n";
      //   cout << "eventList2[" << m << "].first = " << eventList2[m].first << " eventList2[" << m << "].second = " << eventList2[m].second << "\n";
      // }
      if (eventList2[m].first != MuonSystem->runNum) {
        continue;
      }
      if (eventList2[m].second != MuonSystem->evtNumLong) {
        continue;
      }
      // cout << "MuonSystem[" << n << "].runNum = " << MuonSystem->runNum << " MuonSystem[" << n << "].evtNumLong = " << MuonSystem->evtNumLong << "\n";
      // cout << "eventList2[" << m << "].first = " << eventList2[m].first << " eventList2[" << m << "].second = " << eventList2[m].second << "\n";
      // if matched, then add entry to the map
      EventIndexToEventIndexMap[n] = m;
      // cout << "Match  at Index: " << n << " --> " << m << "\n";
      matchFound = true;
      break;
    }

    if (!matchFound)
    {
      numberOfNonMatchedEvents++;
      matchedevent.push_back(false);
    }
    else
    {
      matchedevent.push_back(true);
    }
  }

  cout << "Matched events    = " << NEventsTree1 - numberOfNonMatchedEvents << " / " << NEventsTree1 << " \n";
  cout << "Un-Matched events = " << numberOfNonMatchedEvents << " / " << NEventsTree1 << " \n";

  //*****************************************************************************************
  // Produce Output Tree
  //*****************************************************************************************
  cout << "Merging tree " << MuonSystem->tree_->GetName() << endl;

  cout << "Events in the ntuple: " << MuonSystem->tree_->GetEntries() << endl;
  // create output tree
  outputFile->cd();

  std::vector<std::string> removeBranches{
      "nMuons"
      "muonE",
      "muonPt",
      "muonEta",
      "muonPhi",
      "muonCharge",
      "muonIsLoose",
      "muonIsMedium",
      "muonIsTight",
      "muon_d0",
      "muon_dZ",
      "muon_ip3d",
      "muon_ip3dSignificance",
      "muonType",
      "muonQuality",
      "muon_pileupIso",
      "muon_chargedIso",
      "muon_photonIso",
      "muon_neutralHadIso",
      "muon_ptrel",
      "muon_chargedMiniIso",
      "muon_photonAndNeutralHadronMiniIso",
      "muon_chargedPileupMiniIso",
      "muon_activityMiniIsoAnnulus",
      "muon_passSingleMuTagFilter",
      "muon_passHLTFilter",
      "muon_validFractionTrackerHits",
      "muon_normChi2",
      "muon_chi2LocalPosition",
      "muon_kinkFinder",
      "muon_segmentCompatability",
      "muonIsICHEPMedium",
  };

  for (std::string branch : removeBranches)
  {
    MuonSystem->tree_->SetBranchStatus(branch.c_str(), 0);
  }
  // clone tree with 0 entries, copy all the branch addresses only
  TTree *outputTree = MuonSystem->tree_->CloneTree(0);
  RazorHelper *helper = 0;
  bool isData = true;

  // helper = new RazorHelper(analysisTag, isData, false);

  // loop over Tree1 and add all the branches from tree2
  int numFloatBranches = 25;
  int numIntBranches = 7;
  int numBoolBranches = 13;
  int numCharBranches = 10;

  // nMuons branch
  Int_t fillnMuons;
  ra.fChain->SetBranchStatus("nMuon", 1);
  ra.fChain->SetBranchAddress("nMuon", &ra.nMuons);
  cout << "Adding Branch: "
        << "nMuon"
        << "\n";
  outputTree->Branch("nMuon", &fillnMuons, "nMuon/I");

  // bool branches
  Bool_t fill_HLT_CscCluster_Loose;
  ra.fChain->SetBranchStatus("HLT_CscCluster_Loose", 1);
  ra.fChain->SetBranchAddress("HLT_CscCluster_Loose", &ra.HLT_CscCluster_Loose);
  cout << "Adding Branch: "
        << "HLT_CscCluster_Loose"
        << "\n";
  outputTree->Branch("HLT_CscCluster_Loose", &fill_HLT_CscCluster_Loose, "HLT_CscCluster_Loose/O");

  Bool_t fill_HLT_CscCluster_Medium;
  ra.fChain->SetBranchStatus("HLT_CscCluster_Medium", 1);
  ra.fChain->SetBranchAddress("HLT_CscCluster_Medium", &ra.HLT_CscCluster_Medium);
  cout << "Adding Branch: "
        << "HLT_CscCluster_Medium"
        << "\n";
  outputTree->Branch("HLT_CscCluster_Medium", &fill_HLT_CscCluster_Medium, "HLT_CscCluster_Medium/O");

  Bool_t fill_HLT_CscCluster_Tight;
  ra.fChain->SetBranchStatus("HLT_CscCluster_Tight", 1);
  ra.fChain->SetBranchAddress("HLT_CscCluster_Tight", &ra.HLT_CscCluster_Tight);
  cout << "Adding Branch: "
        << "HLT_CscCluster_Tight"
        << "\n";
  outputTree->Branch("HLT_CscCluster_Tight", &fill_HLT_CscCluster_Tight, "HLT_CscCluster_Tight/O");

  const char *addBranchNamesFloat[numFloatBranches]{
      "Muon_dxy", "Muon_dxyErr", "Muon_dxybs", "Muon_dz", "Muon_dzErr", "Muon_eta",
      "Muon_ip3d", "Muon_jetPtRelv2", "Muon_jetRelIso", "Muon_mass", "Muon_miniPFRelIso_all",
      "Muon_miniPFRelIso_chg", "Muon_mvaLowPt", "Muon_mvaTTH", "Muon_pfRelIso03_all",
      "Muon_pfRelIso03_chg", "Muon_pfRelIso04_all", "Muon_phi", "Muon_pt", "Muon_ptErr", "Muon_segmentComp",
      "Muon_sip3d", "Muon_softMva", "Muon_tkRelIso", "Muon_tunepRelPt"};
  int nMuon = 200;
  Float_t addBranchesInputVarFloat[numFloatBranches][nMuon] = {0};

  Float_t *addBranchesRazorVarFloat[numFloatBranches]{
      ra.muon_dxy, ra.muon_dxyErr, ra.muon_dxybs, ra.muon_dz, ra.muon_dzErr, ra.muon_eta,
      ra.muon_ip3d, ra.muon_jetPtRelv2, ra.muon_jetRelIso, ra.muon_mass, ra.muon_miniPFRelIso_all,
      ra.muon_miniPFRelIso_chg, ra.muon_mvaLowPt, ra.muon_mvaTTH, ra.muon_pfRelIso03_all,
      ra.muon_pfRelIso03_chg, ra.muon_pfRelIso04_all, ra.muon_phi,
      ra.muon_pt, ra.muon_ptErr, ra.muon_segmentComp, ra.muon_sip3d, ra.muon_softMva,
      ra.muon_tkRelIso, ra.muon_tunepRelPt};

  for (int i = 0; i < numFloatBranches; i++)
  {
    cout << "Adding Branch: " << addBranchNamesFloat[i] << "\n";
    // std::string branchName = addBranchNamesFloat[i];
    outputTree->Branch(addBranchNamesFloat[i], addBranchesInputVarFloat[i], TString::Format("%s[nMuon]/F", addBranchNamesFloat[i]));
    ra.fChain->SetBranchStatus(addBranchNamesFloat[i], 1);
    ra.fChain->SetBranchAddress(addBranchNamesFloat[i], addBranchesRazorVarFloat[i]);
  }

  // Same for ints!
  const char *addBranchNamesInt[numIntBranches]{
      "Muon_charge", "Muon_fsrPhotonIdx", "Muon_jetIdx", "Muon_nStations",
      "Muon_nTrackerLayers", "Muon_pdgId", "Muon_tightCharge"};

  Int_t addBranchesInputVarInt[numIntBranches][nMuon] = {0};

  Int_t *addBranchesRazorVarInt[numIntBranches]{
      ra.muon_charge, ra.muon_fsrPhotonIdx, ra.muon_jetIdx, ra.muon_nStations,
      ra.muon_nTrackerLayers, ra.muon_pdgId, ra.muon_tightCharge};

  for (int i = 0; i < numIntBranches; i++)
  {
    cout << "Adding Branch: " << addBranchNamesInt[i] << "\n";
    outputTree->Branch(addBranchNamesInt[i], addBranchesInputVarInt[i], TString::Format("%s[nMuon]/I", addBranchNamesInt[i]));
    ra.fChain->SetBranchStatus(addBranchNamesInt[i], 1);
    ra.fChain->SetBranchAddress(addBranchNamesInt[i], addBranchesRazorVarInt[i]);
  };
  // bools
  const char *addBranchNamesBool[numBoolBranches]{
      "Muon_highPurity", "Muon_inTimeMuon", "Muon_isGlobal", "Muon_isPFcand",
      "Muon_isStandalone", "Muon_isTracker", "Muon_looseId", "Muon_mediumId",
      "Muon_mediumPromptId", "Muon_softId", "Muon_softMvaId", "Muon_tightId",
      "Muon_triggerIdLoose"};

  Bool_t addBranchesInputVarBool[numBoolBranches][nMuon] = {0};

  Bool_t *addBranchesRazorVarBool[numBoolBranches]{
      ra.muon_highPurity, ra.muon_inTimeMuon, ra.muon_isGlobal, ra.muon_isPFcand,
      ra.muon_isStandalone, ra.muon_isTracker, ra.muon_looseId, ra.muon_mediumId,
      ra.muon_mediumPromptId, ra.muon_softId, ra.muon_softMvaId, ra.muon_tightId,
      ra.muon_triggerIdLoose};

  for (int i = 0; i < numBoolBranches; i++)
  {
    cout << "Adding Branch: " << addBranchNamesBool[i] << "\n";
    outputTree->Branch(addBranchNamesBool[i], addBranchesInputVarBool[i], TString::Format("%s[nMuon]/O", addBranchNamesBool[i]));
    ra.fChain->SetBranchStatus(addBranchNamesBool[i], 1);
    ra.fChain->SetBranchAddress(addBranchNamesBool[i], addBranchesRazorVarBool[i]);
  };

  // chars
  const char *addBranchNamesChar[numCharBranches]{
      "Muon_cleanmask", "Muon_highPtId", "Muon_jetNDauCharged", "Muon_miniIsoId",
      "Muon_multiIsoId", "Muon_mvaId", "Muon_mvaLowPtId", "Muon_pfIsoId",
      "Muon_puppiIsoId", "Muon_tkIsoId"};

  UChar_t addBranchesInputVarChar[numCharBranches][nMuon] = {0};

  UChar_t *addBranchesRazorVarChar[numCharBranches]{
      ra.muon_cleanmask, ra.muon_highPtId, ra.muon_jetNDauCharged, ra.muon_miniIsoId,
      ra.muon_multiIsoId, ra.muon_mvaId, ra.muon_mvaLowPtId, ra.muon_pfIsoId,
      ra.muon_puppiIsoId, ra.muon_tkIsoId};

  for (int i = 0; i < numCharBranches; i++)
  {
    cout << "Adding Branch: " << addBranchNamesChar[i] << "\n";
    outputTree->Branch(addBranchNamesChar[i], addBranchesInputVarChar[i], TString::Format("%s[nMuon]/b", addBranchNamesChar[i]));
    ra.fChain->SetBranchStatus(addBranchNamesChar[i], 1);
    ra.fChain->SetBranchAddress(addBranchNamesChar[i], addBranchesRazorVarChar[i]);
  };

  int numEventsProcess = NEventsTree1;
  for (uint n = 0; n < numEventsProcess; n++)
  {
    if (n % 10000 == 0)
      cout << "Processed Event " << n << "\n";

    // Check if found a match
    if (matchedevent[n] == false)
      continue; // #UNCOMMENT#######

    // Get entries
    MuonSystem->tree_->GetEntry(n);
    ra.fChain->GetEntry(EventIndexToEventIndexMap[n]); // #UNCOMMENT#######
    // ra.fChain->GetEntry(n); //#COMMENT#######

    fillnMuons = ra.nMuons;
    fill_HLT_CscCluster_Loose = ra.HLT_CscCluster_Loose;
    fill_HLT_CscCluster_Medium = ra.HLT_CscCluster_Medium;
    fill_HLT_CscCluster_Tight = ra.HLT_CscCluster_Tight;

    for (int i = 0; i < numFloatBranches; i++)
    {
      for (int j = 0; j < ra.nMuons; j++)
      {
        addBranchesInputVarFloat[i][j] = addBranchesRazorVarFloat[i][j];
      }
    }

    for (int i = 0; i < numIntBranches; i++)
    {
      for (int j = 0; j < ra.nMuons; j++)
      {
        addBranchesInputVarInt[i][j] = addBranchesRazorVarInt[i][j];
      }
    }

    for (int i = 0; i < numBoolBranches; i++)
    {
      for (int j = 0; j < ra.nMuons; j++)
      {
        addBranchesInputVarBool[i][j] = addBranchesRazorVarBool[i][j];
      }
    }

    for (int i = 0; i < numCharBranches; i++)
    {
      // addBranchesInputVarChar[i].clear();
      for (int j = 0; j < ra.nMuons; j++)
      {
        addBranchesInputVarChar[i][j] = addBranchesRazorVarChar[i][j];
      }
    }

    // Fill out tree if conditions are met
    outputTree->Fill();
    addBranchesInputVarFloat[numFloatBranches][nMuon] = {0};
    addBranchesInputVarInt[numFloatBranches][nMuon] = {0};
    addBranchesInputVarBool[numFloatBranches][nMuon] = {0};
    addBranchesInputVarChar[numFloatBranches][nMuon] = {0};
  }
  // save information
  outputFile->cd();
  outputTree->Write();
  outputFile->Close();
  delete outputFile;
}