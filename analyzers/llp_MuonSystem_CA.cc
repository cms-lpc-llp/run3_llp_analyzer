// #include "fastjet/ClusterSequence.hh"
// #include "fastjet/ClusterSequence.hh"
// #include "fastjet/Selector.hh"
#include "llp_MuonSystem_CA.h"
#include "RazorHelper.h"
#include "TreeMuonSystemCombination.h"

#include "CACluster.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include "TMath.h"
#include "Math/Vector3D.h"
#include <iostream>
#include <random>

// C++ includes
#include "assert.h"

// ROOT includes
#include "TH1F.h"
#include "TTree.h"
#include "TFile.h"

// using namespace fastjet;
using namespace std::chrono;
using namespace std;
using namespace ROOT::Math;

struct greater_than_pt
{
  inline bool operator()(const TLorentzVector &p1, const TLorentzVector &p2) { return p1.Pt() > p2.Pt(); }
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

void llp_MuonSystem_CA::Analyze(bool isData, int options, string outputfilename, string analysisTag)
{
  // initialization: create one TTree for each analysis box
  cout << "Initializing..." << endl;
  cout << "IsData = " << isData << "\n";
  cout << "options = " << options << "\n";

  //---------------------------
  // options format: MH/MX/ctau/condor: 1000/300/0/1
  // mh can be 3-4 digits, mx is always 3 digits, ctau is one digit(number of zeros), last digit is condor option
  // mh can be 3-4 digits, mx is always 3 digits, ctau is 2 digit(number of zeros), last digit is condor option
  //
  //
  // int mx = int(options/1000)%1000;
  // int mh = options/1000000;
  // int ctau = pow(10, int(options/10)%10) * int(int(options/100)%10);
  //
  // cout<<"mh "<<mh<<", mx "<<mx<<", ctau "<<ctau<<endl;

  bool signalScan = int(options / 10) == 1;
  int option = options % 10;
  // if (options % 1){
  //   option = 1; // used when running condor
  // }
  // else{
  //   option = 0;// used when running locally
  // }

  if (isData)
  {
    std::cout << "[INFO]: running on data with option: " << option << std::endl;
  }
  else
  {
    std::cout << "[INFO]: running on MC with option: " << option << std::endl;
  }
  if (signalScan)
  {
    std::cout << "[INFO]: running with Signal scan" << std::endl;
  }
  else
  {
    std::cout << "[INFO]: running without Signal scan " << option << std::endl;
  }

  const float ELE_MASS = 0.000511;
  const float MU_MASS = 0.105658;
  const float Z_MASS = 91.2;

  if (analysisTag == "")
  {
    analysisTag = "Razor2016_80X";
  }

  //-----------------------------------------------
  // Set up Output File
  //-----------------------------------------------
  string outfilename = outputfilename;
  if (outfilename == "")
    outfilename = "MuonSystem_Tree.root";
  TFile *outFile;
  if (isData || !signalScan)
    outFile = new TFile(outfilename.c_str(), "RECREATE");

  TreeMuonSystemCombination *MuonSystem = new TreeMuonSystemCombination;
  MuonSystem->CreateTree();
  MuonSystem->tree_->SetAutoFlush(0);
  MuonSystem->InitTree();

  // for signals, need one output file for each signal point
  map<pair<int, int>, TFile *> Files2D;
  map<pair<int, int>, TTree *> Trees2D;
  map<pair<int, int>, TH1F *> NEvents2D;
  map<pair<int, int>, TH1F *> accep2D;
  map<pair<int, int>, TH1F *> accep_met2D;
  map<pair<int, int>, TH1F *> Total2D;

  // histogram containing total number of processed events (for normalization)
  TH1F *NEvents = new TH1F("NEvents", "NEvents", 1, 1, 2);
  TH1F *Total = new TH1F("Total", "Total", 1, 1, 2);

  TH1F *accep = new TH1F("accep", "acceptance", 1, 1, 2);
  TH1F *accep_met = new TH1F("accep_met", "acceptance_met", 1, 1, 2);

  TH1F *Nmet200 = new TH1F("Nmet200", "Nmet200", 1, 1, 2);
  TH1F *NmetFilter = new TH1F("NmetFilter", "NmetFilter", 1, 1, 2);
  TH1F *Nlep0 = new TH1F("Nlep0", "Nlep0", 1, 1, 2);
  TH1F *Njet1 = new TH1F("Njet1", "Njet1", 1, 1, 2);
  TH1F *NcosmicVeto = new TH1F("NcosmicVeto", "NcosmicVeto", 1, 1, 2);

  // JetDefinition jet_def( antikt_algorithm, .4 );
  // fastjet::JetDefinition jet_def(fastjet::cambridge_algorithm, 0.4);

  // vector<fastjet::PseudoJet> input_particles;

  char *cmsswPath;
  cmsswPath = getenv("CMSSW_BASE");
  string pathname;
  if (cmsswPath != NULL)
    pathname = string(cmsswPath) + "/src/llp_analyzer/data/JEC/";
  if (cmsswPath != NULL and option == 1)
    pathname = "JEC/"; // run on condor if option == 1

  //--------------------------------
  // Initialize helper
  //--------------------------------
  RazorHelper *helper = 0;
  helper = new RazorHelper(analysisTag, isData);

  //*************************************************************************
  // Look over Input File Events
  //*************************************************************************
  if (fChain == 0)
    return;
  cout << "Total Events: " << fChain->GetEntries() << "\n";
  Long64_t nbytes = 0, nb = 0;
  clock_t start, end;
  start = clock();
  for (Long64_t jentry = 0; jentry < fChain->GetEntriesFast(); jentry++)
  {
    // begin event
    if (jentry % 1000 == 0)
    {
      end = clock();
      double time_taken = double(end - start) / double(CLOCKS_PER_SEC);
      cout << "Processing entry " << jentry << endl;
      cout << "Time taken by program is : " << time_taken << endl;
      start = clock();
    }
    // if (jentry<4000)continue;
    // cout<<jentry<<endl;
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0)
      break;
    // GetEntry(ientry);
    nb = fChain->GetEntry(jentry);
    nbytes += nb;

    // fill normalization histogram
    MuonSystem->InitVariables();
    // std::cout << "deb1 " << jentry << std::endl;

    if (!isData && signalScan)
    {

      string mh_substring = lheComments->substr(lheComments->find("MH-") + 3);
      int mh = stoi(mh_substring.substr(0, mh_substring.find('_')));
      string mx_substring = lheComments->substr(lheComments->find("MS-") + 3);
      int mx = stoi(mx_substring.substr(0, mx_substring.find('_')));
      string ctau_substring = lheComments->substr(lheComments->find("ctauS-") + 6);
      int ctau = stoi(ctau_substring.substr(0, ctau_substring.find('_')));
      MuonSystem->mH = mh;
      MuonSystem->mX = mx;
      MuonSystem->ctau = ctau;

      // if (mh2 != mh || mx2!=mx || ctau2!=ctau) continue;
      // cout<<*lheComments<<endl;

      pair<int, int> signalPair = make_pair(mx, ctau);

      if (Files2D.count(signalPair) == 0)
      { // create file and tree
        // format file name
        string thisFileName = outfilename;
        thisFileName.erase(thisFileName.end() - 5, thisFileName.end());
        thisFileName += "_" + to_string(mx) + "_" + to_string(ctau) + ".root";

        Files2D[signalPair] = new TFile(thisFileName.c_str(), "recreate");
        Trees2D[signalPair] = MuonSystem->tree_->CloneTree(0);
        NEvents2D[signalPair] = new TH1F(Form("NEvents%d%d", mx, ctau), "NEvents", 1, 0.5, 1.5);
        Total2D[signalPair] = new TH1F(Form("Total%d%d", mx, ctau), "Total", 1, 0.5, 1.5);
        accep2D[signalPair] = new TH1F(Form("accep2D%d%d", mx, ctau), "acceptance", 1, 0.5, 1.5);
        accep_met2D[signalPair] = new TH1F(Form("accep_met2D%d%d", mx, ctau), "acceptance_met", 1, 0.5, 1.5);

        cout << "Created new output file " << thisFileName << endl;
      }
      // Fill NEvents hist
      NEvents2D[signalPair]->Fill(1.0, genWeight);
    }
    // event info
    if (isData)
    {
      NEvents->Fill(1);
      MuonSystem->weight = 1;
    }
    else
    {
      MuonSystem->weight = genWeight;
      NEvents->Fill(1, genWeight);
    }
    MuonSystem->runNum = runNum;
    MuonSystem->lumiSec = lumiNum;
    MuonSystem->evtNum = eventNum;

    if (!isData)
    {
      for (int i = 0; i < 2; i++)
      {
        MuonSystem->gLLP_eta[MuonSystem->nGLLP] = gLLP_eta[i];
        MuonSystem->gLLP_phi[MuonSystem->nGLLP] = gLLP_phi[i];
        MuonSystem->gLLP_e[MuonSystem->nGLLP] = gLLP_e[i];
        MuonSystem->gLLP_pt[MuonSystem->nGLLP] = gLLP_pt[i];

        MuonSystem->gLLP_decay_vertex_r[MuonSystem->nGLLP] = sqrt(gLLP_decay_vertex_x[i] * gLLP_decay_vertex_x[i] + gLLP_decay_vertex_y[i] * gLLP_decay_vertex_y[i]);
        MuonSystem->gLLP_decay_vertex_x[MuonSystem->nGLLP] = gLLP_decay_vertex_x[i];
        MuonSystem->gLLP_decay_vertex_y[MuonSystem->nGLLP] = gLLP_decay_vertex_y[i];
        MuonSystem->gLLP_decay_vertex_z[MuonSystem->nGLLP] = gLLP_decay_vertex_z[i];
        float beta = gLLP_beta[i];
        float gLLP_decay_vertex = sqrt(pow(MuonSystem->gLLP_decay_vertex_r[i], 2) + pow(MuonSystem->gLLP_decay_vertex_z[i], 2));
        float gamma = 1.0 / sqrt(1 - beta * beta);
        MuonSystem->gLLP_ctau[MuonSystem->nGLLP] = gLLP_decay_vertex / (beta * gamma);
        MuonSystem->gLLP_beta[MuonSystem->nGLLP] = gLLP_beta[i];

        if (abs(MuonSystem->gLLP_eta[i]) < 2.4 && abs(MuonSystem->gLLP_decay_vertex_z[i]) < 1100 && abs(MuonSystem->gLLP_decay_vertex_z[i]) > 400 && MuonSystem->gLLP_decay_vertex_r[i] < 695.5)
          MuonSystem->gLLP_csc[MuonSystem->nGLLP] = true;
        if (abs(MuonSystem->gLLP_decay_vertex_z[i]) < 661.0 && MuonSystem->gLLP_decay_vertex_r[i] < 800 && MuonSystem->gLLP_decay_vertex_r[i] > 200.0)
          MuonSystem->gLLP_dt[MuonSystem->nGLLP] = true;

        MuonSystem->nGLLP++;
      }

    } // end of isData

    // get NPU
    MuonSystem->npv = nPV;
    MuonSystem->rho = fixedGridRhoFastjetAll;
    MuonSystem->met = metType1Pt;
    MuonSystem->metPhi = metType1Phi;

    if (signalScan && !isData)
      Total2D[make_pair(MuonSystem->mX, MuonSystem->ctau)]->Fill(1.0, genWeight * MuonSystem->pileupWeight);
    if (signalScan && !isData)
    {
      accep2D[make_pair(MuonSystem->mX, MuonSystem->ctau)]->Fill(1.0, genWeight * MuonSystem->pileupWeight);
    }
    else if (!isData)
    {
      accep->Fill(1.0, genWeight * MuonSystem->pileupWeight);
    }

    // Triggers
    for (int i = 0; i < NTriggersMAX; i++)
    {
      MuonSystem->HLTDecision[i] = HLTDecision[i];
    }

    //*************************************************************************
    // Start Object Selection
    //*************************************************************************

    std::vector<leptons> Leptons;
    //-------------------------------
    // Muons
    //-------------------------------
    for (int i = 0; i < nMuons; i++)
    {

      if (!muonIsLoose[i])
        continue;
      if (muonPt[i] < 25)
        continue;
      if (fabs(muonEta[i]) > 2.4)
        continue;

      // remove overlaps
      bool overlap = false;
      for (auto &lep : Leptons)
      {
        if (RazorAnalyzer::deltaR(muonEta[i], muonPhi[i], lep.lepton.Eta(), lep.lepton.Phi()) < 0.3)
          overlap = true;
      }
      if (overlap)
        continue;

      leptons tmpMuon;
      tmpMuon.lepton.SetPtEtaPhiM(muonPt[i], muonEta[i], muonPhi[i], MU_MASS);
      tmpMuon.pdgId = 13 * -1 * muonCharge[i];
      tmpMuon.dZ = muon_dZ[i];
      tmpMuon.passId = muonIsTight[i];
      float muonIso = (muon_chargedIso[i] + fmax(0.0, muon_photonIso[i] + muon_neutralHadIso[i] - 0.5 * muon_pileupIso[i])) / muonPt[i];

      tmpMuon.passLooseIso = muonIso < 0.25;
      tmpMuon.passTightIso = muonIso < 0.15;
      tmpMuon.passVTightIso = muonIso < 0.10;
      tmpMuon.passVVTightIso = muonIso < 0.05;

      tmpMuon.passVetoId = false;
      Leptons.push_back(tmpMuon);
    }

    //-------------------------------
    // Electrons
    //-------------------------------
    for (int i = 0; i < nElectrons; i++)
    {
      if (!ele_passCutBasedIDVeto[i])
        continue;
      if (elePt[i] < 35)
        continue;
      if (fabs(eleEta[i]) > 2.5)
        continue;

      // remove overlaps
      bool overlap = false;
      for (auto &lep : Leptons)
      {
        if (RazorAnalyzer::deltaR(eleEta[i], elePhi[i], lep.lepton.Eta(), lep.lepton.Phi()) < 0.3)
          overlap = true;
      }
      if (overlap)
        continue;
      leptons tmpElectron;
      tmpElectron.lepton.SetPtEtaPhiM(elePt[i], eleEta[i], elePhi[i], ELE_MASS);
      tmpElectron.pdgId = 11 * -1 * eleCharge[i];
      tmpElectron.dZ = ele_dZ[i];
      tmpElectron.passId = ele_passCutBasedIDTight[i];
      Leptons.push_back(tmpElectron);
    }

    sort(Leptons.begin(), Leptons.end(), my_largest_pt);

    for (auto &tmp : Leptons)
    {
      MuonSystem->lepE[MuonSystem->nLeptons] = tmp.lepton.E();
      MuonSystem->lepPt[MuonSystem->nLeptons] = tmp.lepton.Pt();
      MuonSystem->lepEta[MuonSystem->nLeptons] = tmp.lepton.Eta();
      MuonSystem->lepPhi[MuonSystem->nLeptons] = tmp.lepton.Phi();
      MuonSystem->lepPdgId[MuonSystem->nLeptons] = tmp.pdgId;
      MuonSystem->lepDZ[MuonSystem->nLeptons] = tmp.dZ;
      MuonSystem->lepTightId[MuonSystem->nLeptons] = tmp.passId;
      MuonSystem->lepPassLooseIso[MuonSystem->nLeptons] = tmp.passLooseIso;
      MuonSystem->lepPassTightIso[MuonSystem->nLeptons] = tmp.passTightIso;
      MuonSystem->lepPassVTightIso[MuonSystem->nLeptons] = tmp.passVTightIso;
      MuonSystem->lepPassVVTightIso[MuonSystem->nLeptons] = tmp.passVVTightIso;
      MuonSystem->nLeptons++;
    }

    //-----------------------------------------------
    // Select Jets
    //-----------------------------------------------

    std::vector<jets> Jets;

    for (int i = 0; i < nJets; i++)
    {
      if (fabs(jetEta[i]) >= 3.0)
        continue;
      if (jetPt[i] < 20)
        continue;
      if (!jetPassIDLoose[i])
        continue;
      //------------------------------------------------------------
      // exclude selected muons and electrons from the jet collection
      //------------------------------------------------------------
      double deltaR = -1;
      for (auto &lep : Leptons)
      {
        double thisDR = RazorAnalyzer::deltaR(jetEta[i], jetPhi[i], lep.lepton.Eta(), lep.lepton.Phi());
        if (deltaR < 0 || thisDR < deltaR)
          deltaR = thisDR;
      }
      if (deltaR > 0 && deltaR < 0.4)
        continue; // jet matches a selected lepton

      TLorentzVector thisJet = makeTLorentzVector(jetPt[i], jetEta[i], jetPhi[i], jetE[i]);

      jets tmpJet;
      tmpJet.jet = thisJet;
      tmpJet.passId = jetPassIDTight[i];

      Jets.push_back(tmpJet);
    }

    sort(Jets.begin(), Jets.end(), my_largest_pt_jet);

    for (auto &tmp : Jets)
    {
      if (tmp.jet.Pt() < 30)
        continue;

      MuonSystem->jetE[MuonSystem->nJets] = tmp.jet.E();
      MuonSystem->jetPt[MuonSystem->nJets] = tmp.jet.Pt();
      MuonSystem->jetEta[MuonSystem->nJets] = tmp.jet.Eta();
      MuonSystem->jetPhi[MuonSystem->nJets] = tmp.jet.Phi();
      MuonSystem->jetTightPassId[MuonSystem->nJets] = tmp.passId;

      MuonSystem->nJets++;
    }

    MuonSystem->nDTRechits = 0;

    int nDTRechitsChamberMinus12 = 0;
    int nDTRechitsChamberMinus11 = 0;
    int nDTRechitsChamber10 = 0;
    int nDTRechitsChamberPlus11 = 0;
    int nDTRechitsChamberPlus12 = 0;
    int nDTRechitsChamberMinus22 = 0;
    int nDTRechitsChamberMinus21 = 0;
    int nDTRechitsChamber20 = 0;
    int nDTRechitsChamberPlus21 = 0;
    int nDTRechitsChamberPlus22 = 0;
    int nDTRechitsChamberMinus32 = 0;
    int nDTRechitsChamberMinus31 = 0;
    int nDTRechitsChamber30 = 0;
    int nDTRechitsChamberPlus31 = 0;
    int nDTRechitsChamberPlus32 = 0;
    int nDTRechitsChamberMinus42 = 0;
    int nDTRechitsChamberMinus41 = 0;
    int nDTRechitsChamber40 = 0;
    int nDTRechitsChamberPlus41 = 0;
    int nDTRechitsChamberPlus42 = 0;

    for (int i = 0; i < min(N_MAX_DTRECHITS, nDtRechits); i++)
    {

      if (dtRechitStation[i] == 1 && dtRechitWheel[i] == -2)
        nDTRechitsChamberMinus12++;
      if (dtRechitStation[i] == 1 && dtRechitWheel[i] == -1)
        nDTRechitsChamberMinus11++;
      if (dtRechitStation[i] == 1 && dtRechitWheel[i] == 0)
        nDTRechitsChamber10++;
      if (dtRechitStation[i] == 1 && dtRechitWheel[i] == 1)
        nDTRechitsChamberPlus11++;
      if (dtRechitStation[i] == 1 && dtRechitWheel[i] == 2)
        nDTRechitsChamberPlus12++;
      if (dtRechitStation[i] == 2 && dtRechitWheel[i] == -2)
        nDTRechitsChamberMinus22++;
      if (dtRechitStation[i] == 2 && dtRechitWheel[i] == -1)
        nDTRechitsChamberMinus21++;
      if (dtRechitStation[i] == 2 && dtRechitWheel[i] == 0)
        nDTRechitsChamber20++;
      if (dtRechitStation[i] == 2 && dtRechitWheel[i] == 1)
        nDTRechitsChamberPlus21++;
      if (dtRechitStation[i] == 2 && dtRechitWheel[i] == 2)
        nDTRechitsChamberPlus22++;
      if (dtRechitStation[i] == 3 && dtRechitWheel[i] == -2)
        nDTRechitsChamberMinus32++;
      if (dtRechitStation[i] == 3 && dtRechitWheel[i] == -1)
        nDTRechitsChamberMinus31++;
      if (dtRechitStation[i] == 3 && dtRechitWheel[i] == 0)
        nDTRechitsChamber30++;
      if (dtRechitStation[i] == 3 && dtRechitWheel[i] == 1)
        nDTRechitsChamberPlus31++;
      if (dtRechitStation[i] == 3 && dtRechitWheel[i] == 2)
        nDTRechitsChamberPlus32++;
      if (dtRechitStation[i] == 4 && dtRechitWheel[i] == -2)
        nDTRechitsChamberMinus42++;
      if (dtRechitStation[i] == 4 && dtRechitWheel[i] == -1)
        nDTRechitsChamberMinus41++;
      if (dtRechitStation[i] == 4 && dtRechitWheel[i] == 0)
        nDTRechitsChamber40++;
      if (dtRechitStation[i] == 4 && dtRechitWheel[i] == 1)
        nDTRechitsChamberPlus41++;
      if (dtRechitStation[i] == 4 && dtRechitWheel[i] == 2)
        nDTRechitsChamberPlus42++;

      double half_sector = TMath::Pi() / 12.0; // 12 sector of DT in 360 degree

      if (dtRechitCorrectPhi[i] < 1 * half_sector && dtRechitCorrectPhi[i] >= 1 * half_sector)
        MuonSystem->nDTRechitsSector[dtRechitStation[i] - 1][dtRechitWheel[i] + 2][0]++;
      if (dtRechitCorrectPhi[i] < 3 * half_sector && dtRechitCorrectPhi[i] >= 1 * half_sector)
        MuonSystem->nDTRechitsSector[dtRechitStation[i] - 1][dtRechitWheel[i] + 2][1]++;
      if (dtRechitCorrectPhi[i] < 5 * half_sector && dtRechitCorrectPhi[i] >= 3 * half_sector)
        MuonSystem->nDTRechitsSector[dtRechitStation[i] - 1][dtRechitWheel[i] + 2][2]++;
      if (dtRechitCorrectPhi[i] < 7 * half_sector && dtRechitCorrectPhi[i] >= 5 * half_sector)
        MuonSystem->nDTRechitsSector[dtRechitStation[i] - 1][dtRechitWheel[i] + 2][3]++;
      if (dtRechitCorrectPhi[i] < 9 * half_sector && dtRechitCorrectPhi[i] >= 7 * half_sector)
        MuonSystem->nDTRechitsSector[dtRechitStation[i] - 1][dtRechitWheel[i] + 2][4]++;
      if (dtRechitCorrectPhi[i] < 11 * half_sector && dtRechitCorrectPhi[i] >= 9 * half_sector)
        MuonSystem->nDTRechitsSector[dtRechitStation[i] - 1][dtRechitWheel[i] + 2][5]++;
      if (dtRechitCorrectPhi[i] < -11 * half_sector && dtRechitCorrectPhi[i] >= 11 * half_sector)
        MuonSystem->nDTRechitsSector[dtRechitStation[i] - 1][dtRechitWheel[i] + 2][6]++;
      if (dtRechitCorrectPhi[i] < -9 * half_sector && dtRechitCorrectPhi[i] >= -11 * half_sector)
        MuonSystem->nDTRechitsSector[dtRechitStation[i] - 1][dtRechitWheel[i] + 2][7]++;
      if (dtRechitCorrectPhi[i] < -7 * half_sector && dtRechitCorrectPhi[i] >= -9 * half_sector)
        MuonSystem->nDTRechitsSector[dtRechitStation[i] - 1][dtRechitWheel[i] + 2][8]++;
      if (dtRechitCorrectPhi[i] < -5 * half_sector && dtRechitCorrectPhi[i] >= -7 * half_sector)
        MuonSystem->nDTRechitsSector[dtRechitStation[i] - 1][dtRechitWheel[i] + 2][9]++;
      if (dtRechitCorrectPhi[i] < -3 * half_sector && dtRechitCorrectPhi[i] >= -5 * half_sector)
        MuonSystem->nDTRechitsSector[dtRechitStation[i] - 1][dtRechitWheel[i] + 2][10]++;
      if (dtRechitCorrectPhi[i] < -1 * half_sector && dtRechitCorrectPhi[i] >= -3 * half_sector)
        MuonSystem->nDTRechitsSector[dtRechitStation[i] - 1][dtRechitWheel[i] + 2][11]++;
    }

    if (nDTRechitsChamberMinus12 > 50)
      MuonSystem->nDtRings++;
    if (nDTRechitsChamberMinus11 > 50)
      MuonSystem->nDtRings++;
    if (nDTRechitsChamber10 > 50)
      MuonSystem->nDtRings++;
    if (nDTRechitsChamberPlus11 > 50)
      MuonSystem->nDtRings++;
    if (nDTRechitsChamberPlus12 > 50)
      MuonSystem->nDtRings++;
    if (nDTRechitsChamberMinus22 > 50)
      MuonSystem->nDtRings++;
    if (nDTRechitsChamberMinus21 > 50)
      MuonSystem->nDtRings++;
    if (nDTRechitsChamber20 > 50)
      MuonSystem->nDtRings++;
    if (nDTRechitsChamberPlus21 > 50)
      MuonSystem->nDtRings++;
    if (nDTRechitsChamberPlus22 > 50)
      MuonSystem->nDtRings++;
    if (nDTRechitsChamberMinus32 > 50)
      MuonSystem->nDtRings++;
    if (nDTRechitsChamberMinus31 > 50)
      MuonSystem->nDtRings++;
    if (nDTRechitsChamber30 > 50)
      MuonSystem->nDtRings++;
    if (nDTRechitsChamberPlus31 > 50)
      MuonSystem->nDtRings++;
    if (nDTRechitsChamberPlus32 > 50)
      MuonSystem->nDtRings++;
    if (nDTRechitsChamberMinus42 > 50)
      MuonSystem->nDtRings++;
    if (nDTRechitsChamberMinus41 > 50)
      MuonSystem->nDtRings++;
    if (nDTRechitsChamber40 > 50)
      MuonSystem->nDtRings++;
    if (nDTRechitsChamberPlus41 > 50)
      MuonSystem->nDtRings++;
    if (nDTRechitsChamberPlus42 > 50)
      MuonSystem->nDtRings++;

    vector<Rechits> points;
    vector<int> cscRechitsClusterId;
    points.clear();
    int nCscRechitsChamberPlus11 = 0;
    int nCscRechitsChamberPlus12 = 0;
    int nCscRechitsChamberPlus13 = 0;
    int nCscRechitsChamberPlus21 = 0;
    int nCscRechitsChamberPlus22 = 0;
    int nCscRechitsChamberPlus31 = 0;
    int nCscRechitsChamberPlus32 = 0;
    int nCscRechitsChamberPlus41 = 0;
    int nCscRechitsChamberPlus42 = 0;

    int nCscRechitsChamberMinus11 = 0;
    int nCscRechitsChamberMinus12 = 0;
    int nCscRechitsChamberMinus13 = 0;
    int nCscRechitsChamberMinus21 = 0;
    int nCscRechitsChamberMinus22 = 0;
    int nCscRechitsChamberMinus31 = 0;
    int nCscRechitsChamberMinus32 = 0;
    int nCscRechitsChamberMinus41 = 0;
    int nCscRechitsChamberMinus42 = 0;

    for (int i = 0; i < min(N_MAX_CSCRECHITS, ncscRechits); i++)
    {

      // pick out the right bits for chamber
      int chamber = ((cscRechitsDetId[i] >> 3) & 077); // https://github.com/cms-sw/cmssw/blob/master/DataFormats/MuonDetId/interface/CSCDetId.h#L147

      int layer = (cscRechitsDetId[i] & 07);
      Rechits p;
      p.phi = cscRechitsPhi[i];
      p.eta = cscRechitsEta[i];
      p.x = cscRechitsX[i];
      p.y = cscRechitsY[i];
      p.z = cscRechitsZ[i];
      p.t = cscRechitsTpeak[i];
      p.twire = cscRechitsTwire[i];
      p.station = cscRechitsStation[i];
      p.chamber = cscRechitsChamber[i];
      p.layer = layer;
      p.superlayer = 0;
      p.wheel = 0;
      p.clusterID = UNCLASSIFIED;
      points.push_back(p);
      cscRechitsClusterId.push_back(-1);

      if (cscRechitsChamber[i] == 11)
        nCscRechitsChamberPlus11++;
      if (cscRechitsChamber[i] == 12)
        nCscRechitsChamberPlus12++;
      if (cscRechitsChamber[i] == 13)
        nCscRechitsChamberPlus13++;
      if (cscRechitsChamber[i] == 21)
        nCscRechitsChamberPlus21++;
      if (cscRechitsChamber[i] == 22)
        nCscRechitsChamberPlus22++;
      if (cscRechitsChamber[i] == 31)
        nCscRechitsChamberPlus31++;
      if (cscRechitsChamber[i] == 32)
        nCscRechitsChamberPlus32++;
      if (cscRechitsChamber[i] == 41)
        nCscRechitsChamberPlus41++;
      if (cscRechitsChamber[i] == 42)
        nCscRechitsChamberPlus42++;
      if (cscRechitsChamber[i] == -11)
        nCscRechitsChamberMinus11++;
      if (cscRechitsChamber[i] == -12)
        nCscRechitsChamberMinus12++;
      if (cscRechitsChamber[i] == -13)
        nCscRechitsChamberMinus13++;
      if (cscRechitsChamber[i] == -21)
        nCscRechitsChamberMinus21++;
      if (cscRechitsChamber[i] == -22)
        nCscRechitsChamberMinus22++;
      if (cscRechitsChamber[i] == -31)
        nCscRechitsChamberMinus31++;
      if (cscRechitsChamber[i] == -32)
        nCscRechitsChamberMinus32++;
      if (cscRechitsChamber[i] == -41)
        nCscRechitsChamberMinus41++;
      if (cscRechitsChamber[i] == -42)
        nCscRechitsChamberMinus42++;
    }
    MuonSystem->nCscRings = 0;
    if (nCscRechitsChamberPlus11 > 50)
      MuonSystem->nCscRings++;
    if (nCscRechitsChamberPlus12 > 50)
      MuonSystem->nCscRings++;
    if (nCscRechitsChamberPlus13 > 50)
      MuonSystem->nCscRings++;
    if (nCscRechitsChamberPlus21 > 50)
      MuonSystem->nCscRings++;
    if (nCscRechitsChamberPlus22 > 50)
      MuonSystem->nCscRings++;
    if (nCscRechitsChamberPlus31 > 50)
      MuonSystem->nCscRings++;
    if (nCscRechitsChamberPlus32 > 50)
      MuonSystem->nCscRings++;
    if (nCscRechitsChamberPlus41 > 50)
      MuonSystem->nCscRings++;
    if (nCscRechitsChamberPlus42 > 50)
      MuonSystem->nCscRings++;
    if (nCscRechitsChamberMinus11 > 50)
      MuonSystem->nCscRings++;
    if (nCscRechitsChamberMinus12 > 50)
      MuonSystem->nCscRings++;
    if (nCscRechitsChamberMinus13 > 50)
      MuonSystem->nCscRings++;
    if (nCscRechitsChamberMinus21 > 50)
      MuonSystem->nCscRings++;
    if (nCscRechitsChamberMinus22 > 50)
      MuonSystem->nCscRings++;
    if (nCscRechitsChamberMinus31 > 50)
      MuonSystem->nCscRings++;
    if (nCscRechitsChamberMinus32 > 50)
      MuonSystem->nCscRings++;
    if (nCscRechitsChamberMinus41 > 50)
      MuonSystem->nCscRings++;
    if (nCscRechitsChamberMinus42 > 50)
      MuonSystem->nCscRings++;
    // Do DBSCAN Clustering

    int min_point = 50;  // minimum number of Rechitss to call it a cluster
    float epsilon = 0.4; // cluster radius parameter
    CACluster ds(min_point, epsilon, points);
    ds.run();
    ds.clusterProperties();
    // ds.merge_clusters();
    // ds.clusterProperties();
    ds.sort_clusters();

    MuonSystem->nCscRechitClusters = 0;
    for (auto &tmp : ds.clusters)
    {
      MuonSystem->cscRechitClusterX[MuonSystem->nCscRechitClusters] = tmp.x;
      MuonSystem->cscRechitClusterY[MuonSystem->nCscRechitClusters] = tmp.y;
      MuonSystem->cscRechitClusterZ[MuonSystem->nCscRechitClusters] = tmp.z;
      MuonSystem->cscRechitClusterTimeWeighted[MuonSystem->nCscRechitClusters] = tmp.tWeighted;
      MuonSystem->cscRechitClusterTimeSpreadWeightedAll[MuonSystem->nCscRechitClusters] = tmp.TSpreadWeightedAll;

      MuonSystem->cscRechitClusterTime[MuonSystem->nCscRechitClusters] = tmp.tTotal;
      MuonSystem->cscRechitClusterTimeSpread[MuonSystem->nCscRechitClusters] = tmp.TSpread;

      MuonSystem->cscRechitClusterEta[MuonSystem->nCscRechitClusters] = tmp.eta;
      MuonSystem->cscRechitClusterPhi[MuonSystem->nCscRechitClusters] = tmp.phi;

      MuonSystem->cscRechitClusterSize[MuonSystem->nCscRechitClusters] = tmp.nhits;
      MuonSystem->cscRechitClusternXY[MuonSystem->nCscRechitClusters] = tmp.nXY;
      MuonSystem->cscRechitClusternZ[MuonSystem->nCscRechitClusters] = tmp.nZ;
      MuonSystem->cscRechitClusterXSpread[MuonSystem->nCscRechitClusters] = tmp.XSpread;
      MuonSystem->cscRechitClusterYSpread[MuonSystem->nCscRechitClusters] = tmp.YSpread;
      MuonSystem->cscRechitClusterZSpread[MuonSystem->nCscRechitClusters] = tmp.ZSpread;
      MuonSystem->cscRechitClusterXYSpread[MuonSystem->nCscRechitClusters] = tmp.XYSpread;
      MuonSystem->cscRechitClusterRSpread[MuonSystem->nCscRechitClusters] = tmp.RSpread;
      MuonSystem->cscRechitClusterEtaPhiSpread[MuonSystem->nCscRechitClusters] = tmp.EtaPhiSpread;
      MuonSystem->cscRechitClusterEtaSpread[MuonSystem->nCscRechitClusters] = tmp.EtaSpread;
      MuonSystem->cscRechitClusterPhiSpread[MuonSystem->nCscRechitClusters] = tmp.PhiSpread;
      MuonSystem->cscRechitClusterDeltaRSpread[MuonSystem->nCscRechitClusters] = tmp.DeltaRSpread;
      MuonSystem->cscRechitClusterMajorAxis[MuonSystem->nCscRechitClusters] = tmp.MajorAxis;
      MuonSystem->cscRechitClusterMinorAxis[MuonSystem->nCscRechitClusters] = tmp.MinorAxis;
      MuonSystem->cscRechitClusterSkewX[MuonSystem->nCscRechitClusters] = tmp.SkewX;
      MuonSystem->cscRechitClusterSkewY[MuonSystem->nCscRechitClusters] = tmp.SkewY;
      MuonSystem->cscRechitClusterSkewZ[MuonSystem->nCscRechitClusters] = tmp.SkewZ;
      MuonSystem->cscRechitClusterKurtX[MuonSystem->nCscRechitClusters] = tmp.KurtX;
      MuonSystem->cscRechitClusterKurtY[MuonSystem->nCscRechitClusters] = tmp.KurtY;
      MuonSystem->cscRechitClusterKurtZ[MuonSystem->nCscRechitClusters] = tmp.KurtZ;

      MuonSystem->cscRechitClusterNRechitChamberPlus11[MuonSystem->nCscRechitClusters] = tmp.nCscRechitsChamberPlus11;
      MuonSystem->cscRechitClusterNRechitChamberPlus12[MuonSystem->nCscRechitClusters] = tmp.nCscRechitsChamberPlus12;
      MuonSystem->cscRechitClusterNRechitChamberPlus13[MuonSystem->nCscRechitClusters] = tmp.nCscRechitsChamberPlus13;
      MuonSystem->cscRechitClusterNRechitChamberPlus21[MuonSystem->nCscRechitClusters] = tmp.nCscRechitsChamberPlus21;
      MuonSystem->cscRechitClusterNRechitChamberPlus22[MuonSystem->nCscRechitClusters] = tmp.nCscRechitsChamberPlus22;
      MuonSystem->cscRechitClusterNRechitChamberPlus31[MuonSystem->nCscRechitClusters] = tmp.nCscRechitsChamberPlus31;
      MuonSystem->cscRechitClusterNRechitChamberPlus32[MuonSystem->nCscRechitClusters] = tmp.nCscRechitsChamberPlus32;
      MuonSystem->cscRechitClusterNRechitChamberPlus41[MuonSystem->nCscRechitClusters] = tmp.nCscRechitsChamberPlus41;
      MuonSystem->cscRechitClusterNRechitChamberPlus42[MuonSystem->nCscRechitClusters] = tmp.nCscRechitsChamberPlus42;
      MuonSystem->cscRechitClusterNRechitChamberMinus11[MuonSystem->nCscRechitClusters] = tmp.nCscRechitsChamberMinus11;
      MuonSystem->cscRechitClusterNRechitChamberMinus12[MuonSystem->nCscRechitClusters] = tmp.nCscRechitsChamberMinus12;
      MuonSystem->cscRechitClusterNRechitChamberMinus13[MuonSystem->nCscRechitClusters] = tmp.nCscRechitsChamberMinus13;
      MuonSystem->cscRechitClusterNRechitChamberMinus21[MuonSystem->nCscRechitClusters] = tmp.nCscRechitsChamberMinus21;
      MuonSystem->cscRechitClusterNRechitChamberMinus22[MuonSystem->nCscRechitClusters] = tmp.nCscRechitsChamberMinus22;
      MuonSystem->cscRechitClusterNRechitChamberMinus31[MuonSystem->nCscRechitClusters] = tmp.nCscRechitsChamberMinus31;
      MuonSystem->cscRechitClusterNRechitChamberMinus32[MuonSystem->nCscRechitClusters] = tmp.nCscRechitsChamberMinus32;
      MuonSystem->cscRechitClusterNRechitChamberMinus41[MuonSystem->nCscRechitClusters] = tmp.nCscRechitsChamberMinus41;
      MuonSystem->cscRechitClusterNRechitChamberMinus42[MuonSystem->nCscRechitClusters] = tmp.nCscRechitsChamberMinus42;
      MuonSystem->cscRechitClusterMaxChamber[MuonSystem->nCscRechitClusters] = tmp.maxChamber;
      MuonSystem->cscRechitClusterMaxChamberRatio[MuonSystem->nCscRechitClusters] = 1.0 * tmp.maxChamberRechits / tmp.nhits;
      MuonSystem->cscRechitClusterNChamber[MuonSystem->nCscRechitClusters] = tmp.nChamber;
      MuonSystem->cscRechitClusterMaxStation[MuonSystem->nCscRechitClusters] = tmp.maxStation;
      MuonSystem->cscRechitClusterMaxStationRatio[MuonSystem->nCscRechitClusters] = 1.0 * tmp.maxStationRechits / tmp.nhits;

      MuonSystem->cscRechitClusterNStation10[MuonSystem->nCscRechitClusters] = tmp.nStation10;
      MuonSystem->cscRechitClusterAvgStation10[MuonSystem->nCscRechitClusters] = tmp.avgStation10;

      // Jet veto/ muon veto
      MuonSystem->cscRechitClusterJetVetoPt[MuonSystem->nCscRechitClusters] = 0.0;
      MuonSystem->cscRechitClusterJetVetoE[MuonSystem->nCscRechitClusters] = 0.0;
      MuonSystem->cscRechitClusterMuonVetoPt[MuonSystem->nCscRechitClusters] = 0.0;
      MuonSystem->cscRechitClusterMuonVetoE[MuonSystem->nCscRechitClusters] = 0.0;

      // jet veto
      for (int i = 0; i < nJets; i++)
      {
        if (fabs(jetEta[i] > 3.0))
          continue;
        if (RazorAnalyzer::deltaR(jetEta[i], jetPhi[i], MuonSystem->cscRechitClusterEta[MuonSystem->nCscRechitClusters], MuonSystem->cscRechitClusterPhi[MuonSystem->nCscRechitClusters]) < 0.4 && jetPt[i] > MuonSystem->cscRechitClusterJetVetoPt[MuonSystem->nCscRechitClusters])
        {
          MuonSystem->cscRechitClusterJetVetoPt[MuonSystem->nCscRechitClusters] = jetPt[i];
          MuonSystem->cscRechitClusterJetVetoE[MuonSystem->nCscRechitClusters] = jetE[i];
          MuonSystem->cscRechitClusterJetVetoTightId[MuonSystem->nCscRechitClusters] = jetPassIDTight[i];
          MuonSystem->cscRechitClusterJetVetoLooseId[MuonSystem->nCscRechitClusters] = jetPassIDLoose[i];
        }
      }
      float min_deltaR = 15.;
      int index = 999;

      for (int i = 0; i < nMuons; i++)
      {
        if (fabs(muonEta[i] > 3.0))
          continue;
        float muonIso = (muon_chargedIso[i] + fmax(0.0, muon_photonIso[i] + muon_neutralHadIso[i] - 0.5 * muon_pileupIso[i])) / muonPt[i];
        if (RazorAnalyzer::deltaR(muonEta[i], muonPhi[i], MuonSystem->cscRechitClusterEta[MuonSystem->nCscRechitClusters], MuonSystem->cscRechitClusterPhi[MuonSystem->nCscRechitClusters]) < 0.4 && muonPt[i] > MuonSystem->cscRechitClusterMuonVetoPt[MuonSystem->nCscRechitClusters])
        {
          MuonSystem->cscRechitClusterMuonVetoPt[MuonSystem->nCscRechitClusters] = muonPt[i];
          MuonSystem->cscRechitClusterMuonVetoE[MuonSystem->nCscRechitClusters] = muonE[i];
          MuonSystem->cscRechitClusterMuonVetoGlobal[MuonSystem->nCscRechitClusters] = muon_isGlobal[i];
          MuonSystem->cscRechitClusterMuonVetoLooseId[MuonSystem->nCscRechitClusters] = muonIsLoose[i];
        }
      }
      if (!isData)
      {
        // match to gen level LLP
        min_deltaR = 15.;
        index = 999;
        for (int j = 0; j < MuonSystem->nGLLP; j++)
        {

          double current_delta_r = RazorAnalyzer::deltaR(MuonSystem->cscRechitClusterEta[MuonSystem->nCscRechitClusters], MuonSystem->cscRechitClusterPhi[MuonSystem->nCscRechitClusters], MuonSystem->gLLP_eta[j], MuonSystem->gLLP_phi[j]);
          if (current_delta_r < min_deltaR)
          {
            min_deltaR = current_delta_r;
            index = j;
          }
        }
        if (min_deltaR < 0.4)
          MuonSystem->cscRechitCluster_match_gLLP[MuonSystem->nCscRechitClusters] = true;
        else
          MuonSystem->cscRechitCluster_match_gLLP[MuonSystem->nCscRechitClusters] = false;

        MuonSystem->cscRechitCluster_match_gLLP_minDeltaR[MuonSystem->nCscRechitClusters] = min_deltaR;
        MuonSystem->cscRechitCluster_match_gLLP_index[MuonSystem->nCscRechitClusters] = index;
        MuonSystem->cscRechitCluster_match_gLLP_eta[MuonSystem->nCscRechitClusters] = MuonSystem->gLLP_eta[index];
        MuonSystem->cscRechitCluster_match_gLLP_phi[MuonSystem->nCscRechitClusters] = MuonSystem->gLLP_phi[index];
        MuonSystem->cscRechitCluster_match_gLLP_decay_r[MuonSystem->nCscRechitClusters] = MuonSystem->gLLP_decay_vertex_r[index];
        MuonSystem->cscRechitCluster_match_gLLP_decay_z[MuonSystem->nCscRechitClusters] = MuonSystem->gLLP_decay_vertex_z[index];
        MuonSystem->cscRechitCluster_match_gLLP_csc[MuonSystem->nCscRechitClusters] = MuonSystem->gLLP_csc[index];
        MuonSystem->cscRechitCluster_match_gLLP_dt[MuonSystem->nCscRechitClusters] = MuonSystem->gLLP_dt[index];
        MuonSystem->cscRechitCluster_match_gLLP_e[MuonSystem->nCscRechitClusters] = MuonSystem->gLLP_e[index];
      }

      // match to MB1 DT segments
      for (int i = 0; i < nDtSeg; i++)
      {
        if (RazorAnalyzer::deltaR(dtSegEta[i], dtSegPhi[i], MuonSystem->cscRechitClusterEta[MuonSystem->nCscRechitClusters], MuonSystem->cscRechitClusterPhi[MuonSystem->nCscRechitClusters]) < 0.4)
        {
          MuonSystem->cscRechitCluster_match_dtSeg_0p4[MuonSystem->nCscRechitClusters]++;
          if (dtSegStation[i] == 1)
            MuonSystem->cscRechitCluster_match_MB1Seg_0p4[MuonSystem->nCscRechitClusters]++;
        }
      }

      // match to RPC hits in RE1/2
      for (int i = 0; i < nRpc; i++)
      {
        float rpcR = sqrt(rpcX[i] * rpcX[i] + rpcY[i] * rpcY[i]);
        if (RazorAnalyzer::deltaR(rpcEta[i], rpcPhi[i], MuonSystem->cscRechitClusterEta[MuonSystem->nCscRechitClusters], MuonSystem->cscRechitClusterPhi[MuonSystem->nCscRechitClusters]) < 0.4)
        {
          if (rpcR < 461.0 && rpcR > 275 && abs(rpcZ[i]) > 663 && abs(rpcZ[i]) < 730)
            MuonSystem->cscRechitCluster_match_RE12_0p4[MuonSystem->nCscRechitClusters]++;
          if (rpcR < 470 && rpcR > 380 && abs(rpcZ[i]) < 661)
            MuonSystem->cscRechitCluster_match_RB1_0p4[MuonSystem->nCscRechitClusters]++;
        }
      }

      MuonSystem->cscRechitClusterMet_dPhi[MuonSystem->nCscRechitClusters] = RazorAnalyzer::deltaPhi(MuonSystem->cscRechitClusterPhi[MuonSystem->nCscRechitClusters], MuonSystem->metPhi);
      MuonSystem->nCscRechitClusters++;
    }

    // DT cluster

    points.clear();
    // cout<<"here"<<endl;

    for (int i = 0; i < min(N_MAX_DTRECHITS, nDtRechits); i++)
    {
      Rechits p;

      p.phi = dtRechitCorrectPhi[i];
      p.eta = dtRechitCorrectEta[i];
      p.x = dtRechitCorrectX[i];
      p.y = dtRechitCorrectY[i];
      p.z = dtRechitCorrectZ[i];
      p.t = dtRechitTime[i];
      p.twire = dtRechitTime[i];
      p.station = dtRechitStation[i];
      p.chamber = dtRechitWheel[i];
      p.superlayer = dtRechitSuperLayer[i];
      p.wheel = dtRechitWheel[i];
      p.clusterID = UNCLASSIFIED;
      points.push_back(p);
    }
    // cout<<"here"<<endl;

    // Do DBSCAN Clustering
    int min_point_dt = 50;  // minimum number of segments to call it a cluster
    float epsilon_dt = 0.2; // cluster radius parameter
    CACluster ds_dtRechit(min_point_dt, epsilon_dt, points);
    ds_dtRechit.run();
    ds_dtRechit.clusterProperties();
    // ds_dtRechit.merge_clusters();
    // ds_dtRechit.clusterProperties();

    ds_dtRechit.sort_clusters();

    MuonSystem->nDtRechitClusters = 0;

    for (auto &tmp : ds_dtRechit.clusters)
    {

      // remove overlaps
      bool overlap = false;
      for (int i = 0; i < MuonSystem->nCscRechitClusters; i++)
      {
        if (RazorAnalyzer::deltaR(MuonSystem->cscRechitClusterEta[i], MuonSystem->cscRechitClusterPhi[i], tmp.eta, tmp.phi) < 0.4)
          overlap = true;
      }
      if (overlap)
        continue;

      MuonSystem->dtRechitClusterX[MuonSystem->nDtRechitClusters] = tmp.x;
      MuonSystem->dtRechitClusterY[MuonSystem->nDtRechitClusters] = tmp.y;
      MuonSystem->dtRechitClusterZ[MuonSystem->nDtRechitClusters] = tmp.z;
      if (abs(tmp.z) < 126.8)
        MuonSystem->dtRechitClusterWheel[MuonSystem->nDtRechitClusters] = 0;
      else if (tmp.z > 126.8 && tmp.z < 395.4)
        MuonSystem->dtRechitClusterWheel[MuonSystem->nDtRechitClusters] = 1;
      else if (tmp.z < -126.8 && tmp.z > -395.4)
        MuonSystem->dtRechitClusterWheel[MuonSystem->nDtRechitClusters] = -1;
      else if (tmp.z < 0)
        MuonSystem->dtRechitClusterWheel[MuonSystem->nDtRechitClusters] = -2;
      else
        MuonSystem->dtRechitClusterWheel[MuonSystem->nDtRechitClusters] = 2;
      MuonSystem->dtRechitClusterEta[MuonSystem->nDtRechitClusters] = tmp.eta;
      MuonSystem->dtRechitClusterPhi[MuonSystem->nDtRechitClusters] = tmp.phi;
      MuonSystem->dtRechitClusterSize[MuonSystem->nDtRechitClusters] = tmp.nhits;
      MuonSystem->dtRechitClusternXY[MuonSystem->nDtRechitClusters] = tmp.nXY;
      MuonSystem->dtRechitClusternZ[MuonSystem->nDtRechitClusters] = tmp.nZ;
      MuonSystem->dtRechitClusterXSpread[MuonSystem->nDtRechitClusters] = tmp.XSpread;
      MuonSystem->dtRechitClusterYSpread[MuonSystem->nDtRechitClusters] = tmp.YSpread;
      MuonSystem->dtRechitClusterZSpread[MuonSystem->nDtRechitClusters] = tmp.ZSpread;
      MuonSystem->dtRechitClusterXYSpread[MuonSystem->nDtRechitClusters] = tmp.XYSpread;
      MuonSystem->dtRechitClusterRSpread[MuonSystem->nDtRechitClusters] = tmp.RSpread;
      MuonSystem->dtRechitClusterEtaPhiSpread[MuonSystem->nDtRechitClusters] = tmp.EtaPhiSpread;
      MuonSystem->dtRechitClusterEtaSpread[MuonSystem->nDtRechitClusters] = tmp.EtaSpread;
      MuonSystem->dtRechitClusterPhiSpread[MuonSystem->nDtRechitClusters] = tmp.PhiSpread;
      MuonSystem->dtRechitClusterDeltaRSpread[MuonSystem->nDtRechitClusters] = tmp.DeltaRSpread;
      MuonSystem->dtRechitClusterMajorAxis[MuonSystem->nDtRechitClusters] = tmp.MajorAxis;
      MuonSystem->dtRechitClusterMinorAxis[MuonSystem->nDtRechitClusters] = tmp.MinorAxis;
      MuonSystem->dtRechitClusterSkewX[MuonSystem->nDtRechitClusters] = tmp.SkewX;
      MuonSystem->dtRechitClusterSkewY[MuonSystem->nDtRechitClusters] = tmp.SkewY;
      MuonSystem->dtRechitClusterSkewZ[MuonSystem->nDtRechitClusters] = tmp.SkewZ;
      MuonSystem->dtRechitClusterKurtX[MuonSystem->nDtRechitClusters] = tmp.KurtX;
      MuonSystem->dtRechitClusterKurtY[MuonSystem->nDtRechitClusters] = tmp.KurtY;
      MuonSystem->dtRechitClusterKurtZ[MuonSystem->nDtRechitClusters] = tmp.KurtZ;

      unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
      default_random_engine generator(seed);

      // default_random_engine generator;
      uniform_real_distribution<double> distribution(0.0, 1.0);
      float prob = 0.03;
      for (int i = 0; i < 12; ++i)
      {
        if (distribution(generator) < prob)
          MuonSystem->dtRechitClusterNoiseHitStation1[MuonSystem->nDtRechitClusters]++;
      }
      for (int i = 0; i < 12; ++i)
      {
        if (distribution(generator) < prob)
          MuonSystem->dtRechitClusterNoiseHitStation2[MuonSystem->nDtRechitClusters]++;
      }
      for (int i = 0; i < 12; ++i)
      {
        if (distribution(generator) < prob)
          MuonSystem->dtRechitClusterNoiseHitStation3[MuonSystem->nDtRechitClusters]++;
      }
      for (int i = 0; i < 8; ++i)
      {
        if (distribution(generator) < prob)
          MuonSystem->dtRechitClusterNoiseHitStation4[MuonSystem->nDtRechitClusters]++;
      }

      MuonSystem->dtRechitClusterNoiseHit[MuonSystem->nDtRechitClusters] = MuonSystem->dtRechitClusterNoiseHitStation1[MuonSystem->nDtRechitClusters] +
                                                                           MuonSystem->dtRechitClusterNoiseHitStation2[MuonSystem->nDtRechitClusters] +
                                                                           MuonSystem->dtRechitClusterNoiseHitStation3[MuonSystem->nDtRechitClusters] +
                                                                           MuonSystem->dtRechitClusterNoiseHitStation4[MuonSystem->nDtRechitClusters];

      MuonSystem->dtRechitClusterNHitStation1[MuonSystem->nDtRechitClusters] = tmp.nDtRechitsStation1;
      MuonSystem->dtRechitClusterNHitStation2[MuonSystem->nDtRechitClusters] = tmp.nDtRechitsStation2;
      MuonSystem->dtRechitClusterNHitStation3[MuonSystem->nDtRechitClusters] = tmp.nDtRechitsStation3;
      MuonSystem->dtRechitClusterNHitStation4[MuonSystem->nDtRechitClusters] = tmp.nDtRechitsStation4;

      MuonSystem->dtRechitClusterNHitWheel0[MuonSystem->nDtRechitClusters] = tmp.nDtRechitsWheel0;
      MuonSystem->dtRechitClusterNHitWheel1[MuonSystem->nDtRechitClusters] = tmp.nDtRechitsWheel1;
      MuonSystem->dtRechitClusterNHitWheel2[MuonSystem->nDtRechitClusters] = tmp.nDtRechitsWheel2;

      MuonSystem->dtRechitClusterMaxChamber[MuonSystem->nDtRechitClusters] = tmp.maxChamber;
      MuonSystem->dtRechitClusterMaxChamberRatio[MuonSystem->nDtRechitClusters] = 1.0 * tmp.maxChamberRechits / tmp.nhits;
      MuonSystem->dtRechitClusterNChamber[MuonSystem->nDtRechitClusters] = tmp.nChamber;
      MuonSystem->dtRechitClusterMaxStation[MuonSystem->nDtRechitClusters] = tmp.maxStation;
      MuonSystem->dtRechitClusterMaxStationRatio[MuonSystem->nDtRechitClusters] = 1.0 * tmp.maxStationRechits / tmp.nhits;
      MuonSystem->dtRechitClusterNStation10[MuonSystem->nDtRechitClusters] = tmp.nStation10;
      MuonSystem->dtRechitClusterAvgStation10[MuonSystem->nDtRechitClusters] = tmp.avgStation10;

      // Jet veto/ muon veto
      MuonSystem->dtRechitClusterJetVetoPt[MuonSystem->nDtRechitClusters] = 0.0;
      MuonSystem->dtRechitClusterJetVetoE[MuonSystem->nDtRechitClusters] = 0.0;
      MuonSystem->dtRechitClusterMuonVetoPt[MuonSystem->nDtRechitClusters] = 0.0;
      MuonSystem->dtRechitClusterMuonVetoE[MuonSystem->nDtRechitClusters] = 0.0;

      // jet veto
      for (int i = 0; i < nJets; i++)
      {
        if (fabs(jetEta[i] > 3.0))
          continue;
        if (RazorAnalyzer::deltaR(jetEta[i], jetPhi[i], MuonSystem->dtRechitClusterEta[MuonSystem->nDtRechitClusters], MuonSystem->dtRechitClusterPhi[MuonSystem->nDtRechitClusters]) < 0.4 && jetPt[i] > MuonSystem->dtRechitClusterJetVetoPt[MuonSystem->nDtRechitClusters])
        {
          MuonSystem->dtRechitClusterJetVetoPt[MuonSystem->nDtRechitClusters] = jetPt[i];
          MuonSystem->dtRechitClusterJetVetoE[MuonSystem->nDtRechitClusters] = jetE[i];
          MuonSystem->dtRechitClusterJetVetoLooseId[MuonSystem->nDtRechitClusters] = jetPassIDLoose[i];
          MuonSystem->dtRechitClusterJetVetoTightId[MuonSystem->nDtRechitClusters] = jetPassIDTight[i];
        }
      }

      for (int i = 0; i < nMuons; i++)
      {
        if (fabs(muonEta[i] > 3.0))
          continue;
        float muonIso = (muon_chargedIso[i] + fmax(0.0, muon_photonIso[i] + muon_neutralHadIso[i] - 0.5 * muon_pileupIso[i])) / muonPt[i];
        if (RazorAnalyzer::deltaR(muonEta[i], muonPhi[i], MuonSystem->dtRechitClusterEta[MuonSystem->nDtRechitClusters], MuonSystem->dtRechitClusterPhi[MuonSystem->nDtRechitClusters]) < 0.4 && muonPt[i] > MuonSystem->dtRechitClusterMuonVetoPt[MuonSystem->nDtRechitClusters])
        {
          MuonSystem->dtRechitClusterMuonVetoPt[MuonSystem->nDtRechitClusters] = muonPt[i];
          MuonSystem->dtRechitClusterMuonVetoE[MuonSystem->nDtRechitClusters] = muonE[i];
          MuonSystem->dtRechitClusterMuonVetoGlobal[MuonSystem->nDtRechitClusters] = muon_isGlobal[i];
          MuonSystem->dtRechitClusterMuonVetoLooseId[MuonSystem->nDtRechitClusters] = muonIsLoose[i];
        }
      }

      for (int i = 0; i < nDtSeg; i++)
      {
        if (RazorAnalyzer::deltaR(dtSegEta[i], dtSegPhi[i], MuonSystem->dtRechitClusterEta[MuonSystem->nDtRechitClusters], MuonSystem->dtRechitClusterPhi[MuonSystem->nDtRechitClusters]) < 0.4)
        {
          if (dtSegStation[i] == 1)
            MuonSystem->dtRechitClusterNSegStation1[MuonSystem->nDtRechitClusters] += 1;
          if (dtSegStation[i] == 2)
            MuonSystem->dtRechitClusterNSegStation2[MuonSystem->nDtRechitClusters] += 1;
          if (dtSegStation[i] == 3)
            MuonSystem->dtRechitClusterNSegStation3[MuonSystem->nDtRechitClusters] += 1;
          if (dtSegStation[i] == 4)
            MuonSystem->dtRechitClusterNSegStation4[MuonSystem->nDtRechitClusters] += 1;
        }
        if (abs(RazorAnalyzer::deltaPhi(dtSegPhi[i], MuonSystem->dtRechitClusterPhi[MuonSystem->nDtRechitClusters])) > 2)
        {
          if (dtSegStation[i] == 1)
            MuonSystem->dtRechitClusterNOppositeSegStation1[MuonSystem->nDtRechitClusters] += 1;
          if (dtSegStation[i] == 2)
            MuonSystem->dtRechitClusterNOppositeSegStation2[MuonSystem->nDtRechitClusters] += 1;
          if (dtSegStation[i] == 3)
            MuonSystem->dtRechitClusterNOppositeSegStation3[MuonSystem->nDtRechitClusters] += 1;
          if (dtSegStation[i] == 4)
            MuonSystem->dtRechitClusterNOppositeSegStation4[MuonSystem->nDtRechitClusters] += 1;
        }
      }

      // match to gen-level LLP
      float min_deltaR = 15.;
      int index = 999;
      if (!isData)
      {
        for (int j = 0; j < MuonSystem->nGLLP; j++)
        {
          double current_delta_r = RazorAnalyzer::deltaR(tmp.eta, tmp.phi, MuonSystem->gLLP_eta[j], MuonSystem->gLLP_phi[j]);
          if (current_delta_r < min_deltaR)
          {
            min_deltaR = current_delta_r;
            index = j;
          }
        }
        if (min_deltaR < 0.4)
          MuonSystem->dtRechitCluster_match_gLLP[MuonSystem->nDtRechitClusters] = true;
        else
          MuonSystem->dtRechitCluster_match_gLLP[MuonSystem->nDtRechitClusters] = false;

        MuonSystem->dtRechitCluster_match_gLLP_minDeltaR[MuonSystem->nDtRechitClusters] = min_deltaR;
        MuonSystem->dtRechitCluster_match_gLLP_index[MuonSystem->nDtRechitClusters] = index;
        MuonSystem->dtRechitCluster_match_gLLP_eta[MuonSystem->nDtRechitClusters] = MuonSystem->gLLP_eta[index];
        MuonSystem->dtRechitCluster_match_gLLP_phi[MuonSystem->nDtRechitClusters] = MuonSystem->gLLP_phi[index];
        MuonSystem->dtRechitCluster_match_gLLP_decay_r[MuonSystem->nDtRechitClusters] = MuonSystem->gLLP_decay_vertex_r[index];
        MuonSystem->dtRechitCluster_match_gLLP_decay_z[MuonSystem->nDtRechitClusters] = MuonSystem->gLLP_decay_vertex_z[index];
        MuonSystem->dtRechitCluster_match_gLLP_csc[MuonSystem->nDtRechitClusters] = MuonSystem->gLLP_csc[index];
        MuonSystem->dtRechitCluster_match_gLLP_dt[MuonSystem->nDtRechitClusters] = MuonSystem->gLLP_dt[index];
        MuonSystem->dtRechitCluster_match_gLLP_e[MuonSystem->nDtRechitClusters] = MuonSystem->gLLP_e[index];
      }

      // match to MB1 DT segments
      MuonSystem->nCscRechits = ncscRechits;

      for (int i = 0; i < min(N_MAX_DTRECHITS, nDtRechits); i++)
      {
        if (RazorAnalyzer::deltaR(dtRechitCorrectEta[i], dtRechitCorrectPhi[i], MuonSystem->dtRechitClusterEta[MuonSystem->nDtRechitClusters], MuonSystem->dtRechitClusterPhi[MuonSystem->nDtRechitClusters]) < 0.5)
        {
          if (dtRechitStation[i] == 1)
            MuonSystem->dtRechitCluster_match_MB1hits_0p5[MuonSystem->nDtRechitClusters]++;
        }
        if (RazorAnalyzer::deltaR(dtRechitCorrectEta[i], dtRechitCorrectPhi[i], MuonSystem->dtRechitClusterEta[MuonSystem->nDtRechitClusters], MuonSystem->dtRechitClusterPhi[MuonSystem->nDtRechitClusters]) < 0.4)
        {
          if (dtRechitStation[i] == 1)
            MuonSystem->dtRechitCluster_match_MB1hits_0p4[MuonSystem->nDtRechitClusters]++;
        }
        if (abs(dtRechitWheel[i] - MuonSystem->dtRechitClusterWheel[MuonSystem->nDtRechitClusters]) == 1 && dtRechitStation[i] == 1)
        {
          if (abs(RazorAnalyzer::deltaPhi(dtRechitCorrectPhi[i], MuonSystem->dtRechitClusterPhi[MuonSystem->nDtRechitClusters])) < TMath::Pi() / 4.0)
          {
            if (dtRechitWheel[i] - MuonSystem->dtRechitClusterWheel[MuonSystem->nDtRechitClusters] == 1)
              MuonSystem->dtRechitCluster_match_MB1hits_cosmics_plus[MuonSystem->nDtRechitClusters]++;
            else
              MuonSystem->dtRechitCluster_match_MB1hits_cosmics_minus[MuonSystem->nDtRechitClusters]++;
          }
        }
      }

      std::vector<int> dtRechitCluster_match_rpcBx;

      // match to RPC hits with dPhi<0.5 and same wheel in DT
      for (int i = 0; i < nRpc; i++)
      {
        float rpcR = sqrt(rpcX[i] * rpcX[i] + rpcY[i] * rpcY[i]);
        if (rpcRegion[i] != 0)
          continue;
        if (abs(RazorAnalyzer::deltaPhi(rpcPhi[i], MuonSystem->dtRechitClusterPhi[MuonSystem->nDtRechitClusters])) < 0.5)
        {
          if (rpcRing[i] == MuonSystem->dtRechitClusterWheel[MuonSystem->nDtRechitClusters])
          {
            dtRechitCluster_match_rpcBx.push_back(rpcBx[i]);
            MuonSystem->dtRechitCluster_match_RPChits_dPhi0p5[MuonSystem->nDtRechitClusters]++;
            if (rpcR < 470 && rpcR > 380 && abs(rpcZ[i]) < 661)
              MuonSystem->dtRechitCluster_match_RB1_dPhi0p5[MuonSystem->nDtRechitClusters]++;
          }
        }
        if (RazorAnalyzer::deltaR(rpcEta[i], rpcPhi[i], MuonSystem->dtRechitClusterEta[MuonSystem->nDtRechitClusters], MuonSystem->dtRechitClusterPhi[MuonSystem->nDtRechitClusters]) < 0.4)
        {
          if (rpcR < 470 && rpcR > 380 && abs(rpcZ[i]) < 661)
            MuonSystem->dtRechitCluster_match_RB1_0p4[MuonSystem->nDtRechitClusters]++;
        }
      }
      int max_occurence = 0;
      int max_bx = -999;
      for (unsigned int l = 0; l < dtRechitCluster_match_rpcBx.size(); l++)
      {
        int counter = 0;
        for (unsigned int j = 0; j < dtRechitCluster_match_rpcBx.size(); j++)
        {
          if (dtRechitCluster_match_rpcBx[j] == dtRechitCluster_match_rpcBx[l])
            counter++;
        }
        if (counter > max_occurence)
        {
          max_occurence = counter;
          max_bx = dtRechitCluster_match_rpcBx[l];
        }
      }
      MuonSystem->dtRechitCluster_match_RPCBx_dPhi0p5[MuonSystem->nDtRechitClusters] = max_bx;

      MuonSystem->dtRechitClusterMet_dPhi[MuonSystem->nDtRechitClusters] = RazorAnalyzer::deltaPhi(MuonSystem->dtRechitClusterPhi[MuonSystem->nDtRechitClusters], MuonSystem->metPhi);

      MuonSystem->nDtRechitClusters++;
    }

    // if (isData && MuonSystem->nDtRechitClusters + MuonSystem->nCscRechitClusters < 1) continue;

    if (!isData && signalScan)
    {
      pair<int, int> smsPair = make_pair(MuonSystem->mX, MuonSystem->ctau);
      Trees2D[smsPair]->Fill();
    }
    else
    {
      MuonSystem->tree_->Fill();
    }
  }
  if (!isData && signalScan)
  {
    for (auto &filePtr : Files2D)
    {
      cout << "Writing output tree (" << filePtr.second->GetName() << ")" << endl;
      filePtr.second->cd();
      Trees2D[filePtr.first]->Write();
      NEvents2D[filePtr.first]->Write("NEvents");
      Total2D[filePtr.first]->Write("Total");
      accep2D[filePtr.first]->Write("acceptance");
      accep_met2D[filePtr.first]->Write("acceptance_met");
      filePtr.second->Close();
    }
  }
  else if (!isData)
  {
    cout << "Filled Total of " << NEvents->GetBinContent(1) << " Events\n";
    cout << "Writing output trees..." << endl;
    outFile->cd();
    MuonSystem->tree_->Write();
    NEvents->Write();
    accep->Write("acceptance");
    accep_met->Write("acceptance_met");
    outFile->Close();
  }

  else
  {
    cout << "Filled Total of " << NEvents->GetBinContent(1) << " Events\n";
    cout << "Writing output trees..." << endl;
    outFile->cd();
    MuonSystem->tree_->Write();
    Nmet200->Write();
    NmetFilter->Write();
    Nlep0->Write();
    Njet1->Write();
    NcosmicVeto->Write();
    NEvents->Write();
    // outFile->Write();
    outFile->Close();
  }
}

void llp_MuonSystem_CAM::Analyze(bool isData, int options, string outputfilename, string analysisTag)
{
  // initialization: create one TTree for each analysis box
  cout << "Initializing..." << endl;
  cout << "IsData = " << isData << "\n";
  cout << "options = " << options << "\n";

  //---------------------------
  // options format: MH/MX/ctau/condor: 1000/300/0/1
  // mh can be 3-4 digits, mx is always 3 digits, ctau is one digit(number of zeros), last digit is condor option
  // mh can be 3-4 digits, mx is always 3 digits, ctau is 2 digit(number of zeros), last digit is condor option
  //
  //
  // int mx = int(options/1000)%1000;
  // int mh = options/1000000;
  // int ctau = pow(10, int(options/10)%10) * int(int(options/100)%10);
  //
  // cout<<"mh "<<mh<<", mx "<<mx<<", ctau "<<ctau<<endl;

  bool signalScan = int(options / 10) == 1;
  int option = options % 10;
  // if (options % 1){
  //   option = 1; // used when running condor
  // }
  // else{
  //   option = 0;// used when running locally
  // }

  if (isData)
  {
    std::cout << "[INFO]: running on data with option: " << option << std::endl;
  }
  else
  {
    std::cout << "[INFO]: running on MC with option: " << option << std::endl;
  }
  if (signalScan)
  {
    std::cout << "[INFO]: running with Signal scan" << std::endl;
  }
  else
  {
    std::cout << "[INFO]: running without Signal scan " << option << std::endl;
  }

  const float ELE_MASS = 0.000511;
  const float MU_MASS = 0.105658;
  const float Z_MASS = 91.2;

  if (analysisTag == "")
  {
    analysisTag = "Razor2016_80X";
  }

  //-----------------------------------------------
  // Set up Output File
  //-----------------------------------------------
  string outfilename = outputfilename;
  if (outfilename == "")
    outfilename = "MuonSystem_Tree.root";
  TFile *outFile;
  if (isData || !signalScan)
    outFile = new TFile(outfilename.c_str(), "RECREATE");

  TreeMuonSystemCombination *MuonSystem = new TreeMuonSystemCombination;
  MuonSystem->CreateTree();
  MuonSystem->tree_->SetAutoFlush(0);
  MuonSystem->InitTree();

  // for signals, need one output file for each signal point
  map<pair<int, int>, TFile *> Files2D;
  map<pair<int, int>, TTree *> Trees2D;
  map<pair<int, int>, TH1F *> NEvents2D;
  map<pair<int, int>, TH1F *> accep2D;
  map<pair<int, int>, TH1F *> accep_met2D;
  map<pair<int, int>, TH1F *> Total2D;

  // histogram containing total number of processed events (for normalization)
  TH1F *NEvents = new TH1F("NEvents", "NEvents", 1, 1, 2);
  TH1F *Total = new TH1F("Total", "Total", 1, 1, 2);

  TH1F *accep = new TH1F("accep", "acceptance", 1, 1, 2);
  TH1F *accep_met = new TH1F("accep_met", "acceptance_met", 1, 1, 2);

  TH1F *Nmet200 = new TH1F("Nmet200", "Nmet200", 1, 1, 2);
  TH1F *NmetFilter = new TH1F("NmetFilter", "NmetFilter", 1, 1, 2);
  TH1F *Nlep0 = new TH1F("Nlep0", "Nlep0", 1, 1, 2);
  TH1F *Njet1 = new TH1F("Njet1", "Njet1", 1, 1, 2);
  TH1F *NcosmicVeto = new TH1F("NcosmicVeto", "NcosmicVeto", 1, 1, 2);

  // JetDefinition jet_def( antikt_algorithm, .4 );
  // fastjet::JetDefinition jet_def(fastjet::cambridge_algorithm, 0.4);

  // vector<fastjet::PseudoJet> input_particles;

  char *cmsswPath;
  cmsswPath = getenv("CMSSW_BASE");
  string pathname;
  if (cmsswPath != NULL)
    pathname = string(cmsswPath) + "/src/llp_analyzer/data/JEC/";
  if (cmsswPath != NULL and option == 1)
    pathname = "JEC/"; // run on condor if option == 1

  //--------------------------------
  // Initialize helper
  //--------------------------------
  RazorHelper *helper = 0;
  helper = new RazorHelper(analysisTag, isData);

  //*************************************************************************
  // Look over Input File Events
  //*************************************************************************
  if (fChain == 0)
    return;
  cout << "Total Events: " << fChain->GetEntries() << "\n";
  Long64_t nbytes = 0, nb = 0;
  clock_t start, end;
  start = clock();
  for (Long64_t jentry = 0; jentry < fChain->GetEntries(); jentry++)
  {
    // begin event
    if (jentry % 1000 == 0)
    {
      end = clock();
      double time_taken = double(end - start) / double(CLOCKS_PER_SEC);
      cout << "Processing entry " << jentry << endl;
      cout << "Time taken by program is : " << time_taken << endl;
      start = clock();
    }
    // if (jentry<4000)continue;
    // cout<<jentry<<endl;
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0)
      break;
    // GetEntry(ientry);
    nb = fChain->GetEntry(jentry);
    nbytes += nb;

    // fill normalization histogram
    MuonSystem->InitVariables();
    // std::cout << "deb1 " << jentry << std::endl;

    // auto dtRechitTime = dtrechit;                             // dtRechitTime
    auto eleCharge = Electron_charge;       // eleCharge
    auto eleEta = Electron_eta;             // eleEta
    auto elePhi = Electron_phi;             // elePhi
    auto elePt = Electron_pt;               // elePt
    auto ele_dZ = Electron_dz;              // ele_dZ
    bool ele_passCutBasedIDTight[10] = {0}; // ele_passCutBasedIDTight
    bool ele_passCutBasedIDVeto[10] = {0};  // ele_passCutBasedIDVeto

    for (int i = 0; i < nElectron; i++)
    {
      ele_passCutBasedIDTight[i] = Electron_cutBased[i] >= 4;
      ele_passCutBasedIDVeto[i] = Electron_cutBased[i] >= 1;
    }
    auto eventNum = event;                                    // eventNum
    auto fixedGridRhoFastjetAll = Rho_fixedGridRhoFastjetAll; // fixedGridRhoFastjetAll

    Float_t jetE[21] = {0}; // jetE
    for (int i = 0; i < nJet; ++i)
    {
      auto eta = Jet_eta[i];
      auto pt = Jet_pt[i];
      auto pz = pt * TMath::SinH(eta);
      auto mass = Jet_mass[i];
      jetE[i] = TMath::Sqrt(mass * mass + pt * pt + pz * pz);
    }
    auto jetEta = Jet_eta;         // jetEta
    bool jetPassIDLoose[21] = {0}; // jetPassIDLoose
    bool jetPassIDTight[21] = {0}; // jetPassIDTight
    for (int i = 0; i < nJet; ++i)
    {
      jetPassIDLoose[i] = Jet_jetId[i] | 1 != 0;
      jetPassIDTight[i] = Jet_jetId[i] | 2 != 0;
    }
    auto jetPhi = Jet_phi;          // jetPhi
    auto jetPt = Jet_pt;            // jetPt
    auto lumiNum = luminosityBlock; // lumiNum
    auto metType1Phi = MET_phi;     // metType1Phi
    auto metType1Pt = MET_pt;       // metType1Pt
    auto muonCharge = Muon_charge;  // muonCharge
    Float_t muonE[15] = {0};        // muonE
    for (int i = 0; i < nMuon; ++i)
    {
      auto eta = Muon_eta[i];
      auto pt = Muon_pt[i];
      auto pz = pt * TMath::SinH(eta);
      auto mass = MU_MASS;
      muonE[i] = TMath::Sqrt(mass * mass + pt * pt + pz * pz);
    }
    auto muonEta = Muon_eta;                    // muonEta
    auto muonIsLoose = Muon_looseId;            // muonIsLoose
    auto muonIsTight = Muon_tightId;            // muonIsTight
    auto muonPhi = Muon_phi;                    // muonPhi
    auto muonPt = Muon_pt;                      // muonPt
    auto muon_chargedIso = Muon_pfRelIso04_all; // muon_chargedIso
    auto muon_dZ = Muon_dz;                     // muon_dZ
    auto muon_isGlobal = Muon_isGlobal;         // muon_isGlobal
    auto nElectrons = nElectron;                // nElectrons
    auto nJets = nJet;                          // nJets
    auto nMuons = nMuon;                        // nMuons
    int nPV = PV_npvs;                          // nPV
    auto ncscRechits = nCscRechits;             // ncscRechits
    auto runNum = run;                          // runNum

    auto *lheComments = (string *)"123";
    int gLLP_beta[1] = {-1};           // gLLP_beta
    int gLLP_decay_vertex_x[1] = {-1}; // gLLP_decay_vertex_x
    int gLLP_decay_vertex_y[1] = {-1}; // gLLP_decay_vertex_y
    int gLLP_decay_vertex_z[1] = {-1}; // gLLP_decay_vertex_z
    int gLLP_e[1] = {-1};              // gLLP_e
    int gLLP_eta[1] = {-1};            // gLLP_eta
    int gLLP_phi[1] = {-1};            // gLLP_phi
    int gLLP_pt[1] = {-1};             // gLLP_pt
    auto genWeight = -1;               // genWeigh

    // event info

    if (!isData && signalScan)
    {

      string mh_substring = lheComments->substr(lheComments->find("MH-") + 3);
      int mh = stoi(mh_substring.substr(0, mh_substring.find('_')));
      string mx_substring = lheComments->substr(lheComments->find("MS-") + 3);
      int mx = stoi(mx_substring.substr(0, mx_substring.find('_')));
      string ctau_substring = lheComments->substr(lheComments->find("ctauS-") + 6);
      int ctau = stoi(ctau_substring.substr(0, ctau_substring.find('_')));
      MuonSystem->mH = mh;
      MuonSystem->mX = mx;
      MuonSystem->ctau = ctau;

      // if (mh2 != mh || mx2!=mx || ctau2!=ctau) continue;
      // cout<<*lheComments<<endl;

      pair<int, int> signalPair = make_pair(mx, ctau);

      if (Files2D.count(signalPair) == 0)
      { // create file and tree
        // format file name
        string thisFileName = outfilename;
        thisFileName.erase(thisFileName.end() - 5, thisFileName.end());
        thisFileName += "_" + to_string(mx) + "_" + to_string(ctau) + ".root";

        Files2D[signalPair] = new TFile(thisFileName.c_str(), "recreate");
        Trees2D[signalPair] = MuonSystem->tree_->CloneTree(0);
        NEvents2D[signalPair] = new TH1F(Form("NEvents%d%d", mx, ctau), "NEvents", 1, 0.5, 1.5);
        Total2D[signalPair] = new TH1F(Form("Total%d%d", mx, ctau), "Total", 1, 0.5, 1.5);
        accep2D[signalPair] = new TH1F(Form("accep2D%d%d", mx, ctau), "acceptance", 1, 0.5, 1.5);
        accep_met2D[signalPair] = new TH1F(Form("accep_met2D%d%d", mx, ctau), "acceptance_met", 1, 0.5, 1.5);

        cout << "Created new output file " << thisFileName << endl;
      }
      // Fill NEvents hist

      NEvents2D[signalPair]->Fill(1.0, genWeight);
    }
    // event info
    if (isData)
    {
      NEvents->Fill(1);
      MuonSystem->weight = 1;
    }
    else
    {
      MuonSystem->weight = genWeight;
      NEvents->Fill(1, genWeight);
    }

    MuonSystem->runNum = runNum;
    MuonSystem->lumiSec = lumiNum;
    MuonSystem->evtNum = eventNum;

    if (!isData)
    {
      for (int i = 0; i < 2; i++)
      {
        MuonSystem->gLLP_eta[MuonSystem->nGLLP] = gLLP_eta[i];
        MuonSystem->gLLP_phi[MuonSystem->nGLLP] = gLLP_phi[i];
        MuonSystem->gLLP_e[MuonSystem->nGLLP] = gLLP_e[i];
        MuonSystem->gLLP_pt[MuonSystem->nGLLP] = gLLP_pt[i];

        MuonSystem->gLLP_decay_vertex_r[MuonSystem->nGLLP] = sqrt(gLLP_decay_vertex_x[i] * gLLP_decay_vertex_x[i] + gLLP_decay_vertex_y[i] * gLLP_decay_vertex_y[i]);
        MuonSystem->gLLP_decay_vertex_x[MuonSystem->nGLLP] = gLLP_decay_vertex_x[i];
        MuonSystem->gLLP_decay_vertex_y[MuonSystem->nGLLP] = gLLP_decay_vertex_y[i];
        MuonSystem->gLLP_decay_vertex_z[MuonSystem->nGLLP] = gLLP_decay_vertex_z[i];
        float beta = gLLP_beta[i];
        float gLLP_decay_vertex = sqrt(pow(MuonSystem->gLLP_decay_vertex_r[i], 2) + pow(MuonSystem->gLLP_decay_vertex_z[i], 2));
        float gamma = 1.0 / sqrt(1 - beta * beta);
        MuonSystem->gLLP_ctau[MuonSystem->nGLLP] = gLLP_decay_vertex / (beta * gamma);
        MuonSystem->gLLP_beta[MuonSystem->nGLLP] = gLLP_beta[i];

        if (abs(MuonSystem->gLLP_eta[i]) < 2.4 && abs(MuonSystem->gLLP_decay_vertex_z[i]) < 1100 && abs(MuonSystem->gLLP_decay_vertex_z[i]) > 400 && MuonSystem->gLLP_decay_vertex_r[i] < 695.5)
          MuonSystem->gLLP_csc[MuonSystem->nGLLP] = true;
        if (abs(MuonSystem->gLLP_decay_vertex_z[i]) < 661.0 && MuonSystem->gLLP_decay_vertex_r[i] < 800 && MuonSystem->gLLP_decay_vertex_r[i] > 200.0)
          MuonSystem->gLLP_dt[MuonSystem->nGLLP] = true;

        MuonSystem->nGLLP++;
      }

    } // end of isData

    // get NPU
    MuonSystem->npv = nPV;
    MuonSystem->rho = fixedGridRhoFastjetAll;
    MuonSystem->met = metType1Pt;
    MuonSystem->metPhi = metType1Phi;

    if (signalScan && !isData)
      Total2D[make_pair(MuonSystem->mX, MuonSystem->ctau)]->Fill(1.0, genWeight * MuonSystem->pileupWeight);
    if (signalScan && !isData)
    {
      accep2D[make_pair(MuonSystem->mX, MuonSystem->ctau)]->Fill(1.0, genWeight * MuonSystem->pileupWeight);
    }
    else if (!isData)
    {
      accep->Fill(1.0, genWeight * MuonSystem->pileupWeight);
    }

    fillHLT(MuonSystem);
    fillMetFilter(MuonSystem);
    //*************************************************************************
    // Start Object Selection
    //*************************************************************************

    std::vector<leptons> Leptons;
    //-------------------------------
    // Muons
    //-------------------------------
    for (int i = 0; i < nMuons; i++)
    {

      if (!muonIsLoose[i])
        continue;
      if (muonPt[i] < 25)
        continue;
      if (fabs(muonEta[i]) > 2.4)
        continue;

      // remove overlaps
      bool overlap = false;
      for (auto &lep : Leptons)
      {
        if (RazorAnalyzerMerged::deltaR(muonEta[i], muonPhi[i], lep.lepton.Eta(), lep.lepton.Phi()) < 0.3)
          overlap = true;
      }
      if (overlap)
        continue;

      leptons tmpMuon;
      tmpMuon.lepton.SetPtEtaPhiM(muonPt[i], muonEta[i], muonPhi[i], MU_MASS);
      tmpMuon.pdgId = 13 * -1 * muonCharge[i];
      tmpMuon.dZ = muon_dZ[i];
      tmpMuon.passId = muonIsTight[i];
      float muonIso = Muon_pfRelIso04_all[i];

      tmpMuon.passLooseIso = muonIso < 0.25;
      tmpMuon.passTightIso = muonIso < 0.15;
      tmpMuon.passVTightIso = muonIso < 0.10;
      tmpMuon.passVVTightIso = muonIso < 0.05;

      tmpMuon.passVetoId = false;
      Leptons.push_back(tmpMuon);
    }

    //-------------------------------
    // Electrons
    //-------------------------------
    for (int i = 0; i < nElectrons; i++)
    {
      if (!ele_passCutBasedIDVeto[i])
        continue;
      if (elePt[i] < 35)
        continue;
      if (fabs(eleEta[i]) > 2.5)
        continue;

      // remove overlaps
      bool overlap = false;
      for (auto &lep : Leptons)
      {
        if (RazorAnalyzerMerged::deltaR(eleEta[i], elePhi[i], lep.lepton.Eta(), lep.lepton.Phi()) < 0.3)
          overlap = true;
      }
      if (overlap)
        continue;
      leptons tmpElectron;
      tmpElectron.lepton.SetPtEtaPhiM(elePt[i], eleEta[i], elePhi[i], ELE_MASS);
      tmpElectron.pdgId = 11 * -1 * eleCharge[i];
      tmpElectron.dZ = ele_dZ[i];
      tmpElectron.passId = ele_passCutBasedIDTight[i];
      Leptons.push_back(tmpElectron);
    }

    sort(Leptons.begin(), Leptons.end(), my_largest_pt);

    for (auto &tmp : Leptons)
    {
      MuonSystem->lepE[MuonSystem->nLeptons] = tmp.lepton.E();
      MuonSystem->lepPt[MuonSystem->nLeptons] = tmp.lepton.Pt();
      MuonSystem->lepEta[MuonSystem->nLeptons] = tmp.lepton.Eta();
      MuonSystem->lepPhi[MuonSystem->nLeptons] = tmp.lepton.Phi();
      MuonSystem->lepPdgId[MuonSystem->nLeptons] = tmp.pdgId;
      MuonSystem->lepDZ[MuonSystem->nLeptons] = tmp.dZ;
      MuonSystem->lepTightId[MuonSystem->nLeptons] = tmp.passId;
      MuonSystem->lepPassLooseIso[MuonSystem->nLeptons] = tmp.passLooseIso;
      MuonSystem->lepPassTightIso[MuonSystem->nLeptons] = tmp.passTightIso;
      MuonSystem->lepPassVTightIso[MuonSystem->nLeptons] = tmp.passVTightIso;
      MuonSystem->lepPassVVTightIso[MuonSystem->nLeptons] = tmp.passVVTightIso;
      MuonSystem->nLeptons++;
    }

    //-----------------------------------------------
    // Select Jets
    //-----------------------------------------------

    std::vector<jets> Jets;

    for (int i = 0; i < nJets; i++)
    {
      if (fabs(jetEta[i]) >= 3.0)
        continue;
      if (jetPt[i] < 20)
        continue;
      if (!jetPassIDLoose[i])
        continue;
      //------------------------------------------------------------
      // exclude selected muons and electrons from the jet collection
      //------------------------------------------------------------
      double deltaR = -1;
      for (auto &lep : Leptons)
      {
        double thisDR = RazorAnalyzerMerged::deltaR(jetEta[i], jetPhi[i], lep.lepton.Eta(), lep.lepton.Phi());
        if (deltaR < 0 || thisDR < deltaR)
          deltaR = thisDR;
      }
      if (deltaR > 0 && deltaR < 0.4)
        continue; // jet matches a selected lepton

      TLorentzVector thisJet = makeTLorentzVector(jetPt[i], jetEta[i], jetPhi[i], jetE[i]);

      jets tmpJet;
      tmpJet.jet = thisJet;
      tmpJet.passId = jetPassIDTight[i];

      Jets.push_back(tmpJet);
    }

    sort(Jets.begin(), Jets.end(), my_largest_pt_jet);

    for (auto &tmp : Jets)
    {
      if (tmp.jet.Pt() < 30)
        continue;

      MuonSystem->jetE[MuonSystem->nJets] = tmp.jet.E();
      MuonSystem->jetPt[MuonSystem->nJets] = tmp.jet.Pt();
      MuonSystem->jetEta[MuonSystem->nJets] = tmp.jet.Eta();
      MuonSystem->jetPhi[MuonSystem->nJets] = tmp.jet.Phi();
      MuonSystem->jetTightPassId[MuonSystem->nJets] = tmp.passId;

      MuonSystem->nJets++;
    }

    MuonSystem->nDTRechits = 0;

    int nDTRechitsChamberMinus12 = 0;
    int nDTRechitsChamberMinus11 = 0;
    int nDTRechitsChamber10 = 0;
    int nDTRechitsChamberPlus11 = 0;
    int nDTRechitsChamberPlus12 = 0;
    int nDTRechitsChamberMinus22 = 0;
    int nDTRechitsChamberMinus21 = 0;
    int nDTRechitsChamber20 = 0;
    int nDTRechitsChamberPlus21 = 0;
    int nDTRechitsChamberPlus22 = 0;
    int nDTRechitsChamberMinus32 = 0;
    int nDTRechitsChamberMinus31 = 0;
    int nDTRechitsChamber30 = 0;
    int nDTRechitsChamberPlus31 = 0;
    int nDTRechitsChamberPlus32 = 0;
    int nDTRechitsChamberMinus42 = 0;
    int nDTRechitsChamberMinus41 = 0;
    int nDTRechitsChamber40 = 0;
    int nDTRechitsChamberPlus41 = 0;
    int nDTRechitsChamberPlus42 = 0;

    for (int i = 0; i < min(N_MAX_DTRECHITS, nDtRechits); i++)
    {

      if (dtRechitStation[i] == 1 && dtRechitWheel[i] == -2)
        nDTRechitsChamberMinus12++;
      if (dtRechitStation[i] == 1 && dtRechitWheel[i] == -1)
        nDTRechitsChamberMinus11++;
      if (dtRechitStation[i] == 1 && dtRechitWheel[i] == 0)
        nDTRechitsChamber10++;
      if (dtRechitStation[i] == 1 && dtRechitWheel[i] == 1)
        nDTRechitsChamberPlus11++;
      if (dtRechitStation[i] == 1 && dtRechitWheel[i] == 2)
        nDTRechitsChamberPlus12++;
      if (dtRechitStation[i] == 2 && dtRechitWheel[i] == -2)
        nDTRechitsChamberMinus22++;
      if (dtRechitStation[i] == 2 && dtRechitWheel[i] == -1)
        nDTRechitsChamberMinus21++;
      if (dtRechitStation[i] == 2 && dtRechitWheel[i] == 0)
        nDTRechitsChamber20++;
      if (dtRechitStation[i] == 2 && dtRechitWheel[i] == 1)
        nDTRechitsChamberPlus21++;
      if (dtRechitStation[i] == 2 && dtRechitWheel[i] == 2)
        nDTRechitsChamberPlus22++;
      if (dtRechitStation[i] == 3 && dtRechitWheel[i] == -2)
        nDTRechitsChamberMinus32++;
      if (dtRechitStation[i] == 3 && dtRechitWheel[i] == -1)
        nDTRechitsChamberMinus31++;
      if (dtRechitStation[i] == 3 && dtRechitWheel[i] == 0)
        nDTRechitsChamber30++;
      if (dtRechitStation[i] == 3 && dtRechitWheel[i] == 1)
        nDTRechitsChamberPlus31++;
      if (dtRechitStation[i] == 3 && dtRechitWheel[i] == 2)
        nDTRechitsChamberPlus32++;
      if (dtRechitStation[i] == 4 && dtRechitWheel[i] == -2)
        nDTRechitsChamberMinus42++;
      if (dtRechitStation[i] == 4 && dtRechitWheel[i] == -1)
        nDTRechitsChamberMinus41++;
      if (dtRechitStation[i] == 4 && dtRechitWheel[i] == 0)
        nDTRechitsChamber40++;
      if (dtRechitStation[i] == 4 && dtRechitWheel[i] == 1)
        nDTRechitsChamberPlus41++;
      if (dtRechitStation[i] == 4 && dtRechitWheel[i] == 2)
        nDTRechitsChamberPlus42++;

      double half_sector = TMath::Pi() / 12.0; // 12 sector of DT in 360 degree

      if (dtRechitCorrectPhi[i] < 1 * half_sector && dtRechitCorrectPhi[i] >= 1 * half_sector)
        MuonSystem->nDTRechitsSector[dtRechitStation[i] - 1][dtRechitWheel[i] + 2][0]++;
      if (dtRechitCorrectPhi[i] < 3 * half_sector && dtRechitCorrectPhi[i] >= 1 * half_sector)
        MuonSystem->nDTRechitsSector[dtRechitStation[i] - 1][dtRechitWheel[i] + 2][1]++;
      if (dtRechitCorrectPhi[i] < 5 * half_sector && dtRechitCorrectPhi[i] >= 3 * half_sector)
        MuonSystem->nDTRechitsSector[dtRechitStation[i] - 1][dtRechitWheel[i] + 2][2]++;
      if (dtRechitCorrectPhi[i] < 7 * half_sector && dtRechitCorrectPhi[i] >= 5 * half_sector)
        MuonSystem->nDTRechitsSector[dtRechitStation[i] - 1][dtRechitWheel[i] + 2][3]++;
      if (dtRechitCorrectPhi[i] < 9 * half_sector && dtRechitCorrectPhi[i] >= 7 * half_sector)
        MuonSystem->nDTRechitsSector[dtRechitStation[i] - 1][dtRechitWheel[i] + 2][4]++;
      if (dtRechitCorrectPhi[i] < 11 * half_sector && dtRechitCorrectPhi[i] >= 9 * half_sector)
        MuonSystem->nDTRechitsSector[dtRechitStation[i] - 1][dtRechitWheel[i] + 2][5]++;
      if (dtRechitCorrectPhi[i] < -11 * half_sector && dtRechitCorrectPhi[i] >= 11 * half_sector)
        MuonSystem->nDTRechitsSector[dtRechitStation[i] - 1][dtRechitWheel[i] + 2][6]++;
      if (dtRechitCorrectPhi[i] < -9 * half_sector && dtRechitCorrectPhi[i] >= -11 * half_sector)
        MuonSystem->nDTRechitsSector[dtRechitStation[i] - 1][dtRechitWheel[i] + 2][7]++;
      if (dtRechitCorrectPhi[i] < -7 * half_sector && dtRechitCorrectPhi[i] >= -9 * half_sector)
        MuonSystem->nDTRechitsSector[dtRechitStation[i] - 1][dtRechitWheel[i] + 2][8]++;
      if (dtRechitCorrectPhi[i] < -5 * half_sector && dtRechitCorrectPhi[i] >= -7 * half_sector)
        MuonSystem->nDTRechitsSector[dtRechitStation[i] - 1][dtRechitWheel[i] + 2][9]++;
      if (dtRechitCorrectPhi[i] < -3 * half_sector && dtRechitCorrectPhi[i] >= -5 * half_sector)
        MuonSystem->nDTRechitsSector[dtRechitStation[i] - 1][dtRechitWheel[i] + 2][10]++;
      if (dtRechitCorrectPhi[i] < -1 * half_sector && dtRechitCorrectPhi[i] >= -3 * half_sector)
        MuonSystem->nDTRechitsSector[dtRechitStation[i] - 1][dtRechitWheel[i] + 2][11]++;
    }

    if (nDTRechitsChamberMinus12 > 50)
      MuonSystem->nDtRings++;
    if (nDTRechitsChamberMinus11 > 50)
      MuonSystem->nDtRings++;
    if (nDTRechitsChamber10 > 50)
      MuonSystem->nDtRings++;
    if (nDTRechitsChamberPlus11 > 50)
      MuonSystem->nDtRings++;
    if (nDTRechitsChamberPlus12 > 50)
      MuonSystem->nDtRings++;
    if (nDTRechitsChamberMinus22 > 50)
      MuonSystem->nDtRings++;
    if (nDTRechitsChamberMinus21 > 50)
      MuonSystem->nDtRings++;
    if (nDTRechitsChamber20 > 50)
      MuonSystem->nDtRings++;
    if (nDTRechitsChamberPlus21 > 50)
      MuonSystem->nDtRings++;
    if (nDTRechitsChamberPlus22 > 50)
      MuonSystem->nDtRings++;
    if (nDTRechitsChamberMinus32 > 50)
      MuonSystem->nDtRings++;
    if (nDTRechitsChamberMinus31 > 50)
      MuonSystem->nDtRings++;
    if (nDTRechitsChamber30 > 50)
      MuonSystem->nDtRings++;
    if (nDTRechitsChamberPlus31 > 50)
      MuonSystem->nDtRings++;
    if (nDTRechitsChamberPlus32 > 50)
      MuonSystem->nDtRings++;
    if (nDTRechitsChamberMinus42 > 50)
      MuonSystem->nDtRings++;
    if (nDTRechitsChamberMinus41 > 50)
      MuonSystem->nDtRings++;
    if (nDTRechitsChamber40 > 50)
      MuonSystem->nDtRings++;
    if (nDTRechitsChamberPlus41 > 50)
      MuonSystem->nDtRings++;
    if (nDTRechitsChamberPlus42 > 50)
      MuonSystem->nDtRings++;

    vector<Rechits> points;
    vector<int> cscRechitsClusterId;
    points.clear();
    int nCscRechitsChamberPlus11 = 0;
    int nCscRechitsChamberPlus12 = 0;
    int nCscRechitsChamberPlus13 = 0;
    int nCscRechitsChamberPlus21 = 0;
    int nCscRechitsChamberPlus22 = 0;
    int nCscRechitsChamberPlus31 = 0;
    int nCscRechitsChamberPlus32 = 0;
    int nCscRechitsChamberPlus41 = 0;
    int nCscRechitsChamberPlus42 = 0;

    int nCscRechitsChamberMinus11 = 0;
    int nCscRechitsChamberMinus12 = 0;
    int nCscRechitsChamberMinus13 = 0;
    int nCscRechitsChamberMinus21 = 0;
    int nCscRechitsChamberMinus22 = 0;
    int nCscRechitsChamberMinus31 = 0;
    int nCscRechitsChamberMinus32 = 0;
    int nCscRechitsChamberMinus41 = 0;
    int nCscRechitsChamberMinus42 = 0;

    for (int i = 0; i < min(N_MAX_CSCRECHITS, ncscRechits); i++)
    {

      // pick out the right bits for chamber
      int chamber = ((cscRechitsDetId[i] >> 3) & 077); // https://github.com/cms-sw/cmssw/blob/master/DataFormats/MuonDetId/interface/CSCDetId.h#L147

      int layer = (cscRechitsDetId[i] & 07);
      Rechits p;
      p.phi = cscRechitsPhi[i];
      p.eta = cscRechitsEta[i];
      p.x = cscRechitsX[i];
      p.y = cscRechitsY[i];
      p.z = cscRechitsZ[i];
      p.t = cscRechitsTpeak[i];
      p.twire = cscRechitsTwire[i];
      p.station = cscRechitsStation[i];
      p.chamber = cscRechitsChamber[i];
      p.layer = layer;
      p.superlayer = 0;
      p.wheel = 0;
      p.clusterID = UNCLASSIFIED;
      points.push_back(p);
      cscRechitsClusterId.push_back(-1);

      if (cscRechitsChamber[i] == 11)
        nCscRechitsChamberPlus11++;
      if (cscRechitsChamber[i] == 12)
        nCscRechitsChamberPlus12++;
      if (cscRechitsChamber[i] == 13)
        nCscRechitsChamberPlus13++;
      if (cscRechitsChamber[i] == 21)
        nCscRechitsChamberPlus21++;
      if (cscRechitsChamber[i] == 22)
        nCscRechitsChamberPlus22++;
      if (cscRechitsChamber[i] == 31)
        nCscRechitsChamberPlus31++;
      if (cscRechitsChamber[i] == 32)
        nCscRechitsChamberPlus32++;
      if (cscRechitsChamber[i] == 41)
        nCscRechitsChamberPlus41++;
      if (cscRechitsChamber[i] == 42)
        nCscRechitsChamberPlus42++;
      if (cscRechitsChamber[i] == -11)
        nCscRechitsChamberMinus11++;
      if (cscRechitsChamber[i] == -12)
        nCscRechitsChamberMinus12++;
      if (cscRechitsChamber[i] == -13)
        nCscRechitsChamberMinus13++;
      if (cscRechitsChamber[i] == -21)
        nCscRechitsChamberMinus21++;
      if (cscRechitsChamber[i] == -22)
        nCscRechitsChamberMinus22++;
      if (cscRechitsChamber[i] == -31)
        nCscRechitsChamberMinus31++;
      if (cscRechitsChamber[i] == -32)
        nCscRechitsChamberMinus32++;
      if (cscRechitsChamber[i] == -41)
        nCscRechitsChamberMinus41++;
      if (cscRechitsChamber[i] == -42)
        nCscRechitsChamberMinus42++;
    }
    MuonSystem->nCscRings = 0;
    if (nCscRechitsChamberPlus11 > 50)
      MuonSystem->nCscRings++;
    if (nCscRechitsChamberPlus12 > 50)
      MuonSystem->nCscRings++;
    if (nCscRechitsChamberPlus13 > 50)
      MuonSystem->nCscRings++;
    if (nCscRechitsChamberPlus21 > 50)
      MuonSystem->nCscRings++;
    if (nCscRechitsChamberPlus22 > 50)
      MuonSystem->nCscRings++;
    if (nCscRechitsChamberPlus31 > 50)
      MuonSystem->nCscRings++;
    if (nCscRechitsChamberPlus32 > 50)
      MuonSystem->nCscRings++;
    if (nCscRechitsChamberPlus41 > 50)
      MuonSystem->nCscRings++;
    if (nCscRechitsChamberPlus42 > 50)
      MuonSystem->nCscRings++;
    if (nCscRechitsChamberMinus11 > 50)
      MuonSystem->nCscRings++;
    if (nCscRechitsChamberMinus12 > 50)
      MuonSystem->nCscRings++;
    if (nCscRechitsChamberMinus13 > 50)
      MuonSystem->nCscRings++;
    if (nCscRechitsChamberMinus21 > 50)
      MuonSystem->nCscRings++;
    if (nCscRechitsChamberMinus22 > 50)
      MuonSystem->nCscRings++;
    if (nCscRechitsChamberMinus31 > 50)
      MuonSystem->nCscRings++;
    if (nCscRechitsChamberMinus32 > 50)
      MuonSystem->nCscRings++;
    if (nCscRechitsChamberMinus41 > 50)
      MuonSystem->nCscRings++;
    if (nCscRechitsChamberMinus42 > 50)
      MuonSystem->nCscRings++;
    // Do DBSCAN Clustering

    int min_point = 50;  // minimum number of Rechitss to call it a cluster
    float epsilon = 0.4; // cluster radius parameter
    CACluster ds(min_point, epsilon, points);
    ds.run();
    ds.clusterProperties();
    // ds.merge_clusters();
    // ds.clusterProperties();
    ds.sort_clusters();

    MuonSystem->nCscRechitClusters = 0;
    for (auto &tmp : ds.clusters)
    {
      MuonSystem->cscRechitClusterX[MuonSystem->nCscRechitClusters] = tmp.x;
      MuonSystem->cscRechitClusterY[MuonSystem->nCscRechitClusters] = tmp.y;
      MuonSystem->cscRechitClusterZ[MuonSystem->nCscRechitClusters] = tmp.z;
      MuonSystem->cscRechitClusterTimeWeighted[MuonSystem->nCscRechitClusters] = tmp.tWeighted;
      MuonSystem->cscRechitClusterTimeSpreadWeightedAll[MuonSystem->nCscRechitClusters] = tmp.TSpreadWeightedAll;

      MuonSystem->cscRechitClusterTime[MuonSystem->nCscRechitClusters] = tmp.tTotal;
      MuonSystem->cscRechitClusterTimeSpread[MuonSystem->nCscRechitClusters] = tmp.TSpread;

      MuonSystem->cscRechitClusterEta[MuonSystem->nCscRechitClusters] = tmp.eta;
      MuonSystem->cscRechitClusterPhi[MuonSystem->nCscRechitClusters] = tmp.phi;

      MuonSystem->cscRechitClusterSize[MuonSystem->nCscRechitClusters] = tmp.nhits;
      MuonSystem->cscRechitClusternXY[MuonSystem->nCscRechitClusters] = tmp.nXY;
      MuonSystem->cscRechitClusternZ[MuonSystem->nCscRechitClusters] = tmp.nZ;
      MuonSystem->cscRechitClusterXSpread[MuonSystem->nCscRechitClusters] = tmp.XSpread;
      MuonSystem->cscRechitClusterYSpread[MuonSystem->nCscRechitClusters] = tmp.YSpread;
      MuonSystem->cscRechitClusterZSpread[MuonSystem->nCscRechitClusters] = tmp.ZSpread;
      MuonSystem->cscRechitClusterXYSpread[MuonSystem->nCscRechitClusters] = tmp.XYSpread;
      MuonSystem->cscRechitClusterRSpread[MuonSystem->nCscRechitClusters] = tmp.RSpread;
      MuonSystem->cscRechitClusterEtaPhiSpread[MuonSystem->nCscRechitClusters] = tmp.EtaPhiSpread;
      MuonSystem->cscRechitClusterEtaSpread[MuonSystem->nCscRechitClusters] = tmp.EtaSpread;
      MuonSystem->cscRechitClusterPhiSpread[MuonSystem->nCscRechitClusters] = tmp.PhiSpread;
      MuonSystem->cscRechitClusterDeltaRSpread[MuonSystem->nCscRechitClusters] = tmp.DeltaRSpread;
      MuonSystem->cscRechitClusterMajorAxis[MuonSystem->nCscRechitClusters] = tmp.MajorAxis;
      MuonSystem->cscRechitClusterMinorAxis[MuonSystem->nCscRechitClusters] = tmp.MinorAxis;
      MuonSystem->cscRechitClusterSkewX[MuonSystem->nCscRechitClusters] = tmp.SkewX;
      MuonSystem->cscRechitClusterSkewY[MuonSystem->nCscRechitClusters] = tmp.SkewY;
      MuonSystem->cscRechitClusterSkewZ[MuonSystem->nCscRechitClusters] = tmp.SkewZ;
      MuonSystem->cscRechitClusterKurtX[MuonSystem->nCscRechitClusters] = tmp.KurtX;
      MuonSystem->cscRechitClusterKurtY[MuonSystem->nCscRechitClusters] = tmp.KurtY;
      MuonSystem->cscRechitClusterKurtZ[MuonSystem->nCscRechitClusters] = tmp.KurtZ;

      MuonSystem->cscRechitClusterNRechitChamberPlus11[MuonSystem->nCscRechitClusters] = tmp.nCscRechitsChamberPlus11;
      MuonSystem->cscRechitClusterNRechitChamberPlus12[MuonSystem->nCscRechitClusters] = tmp.nCscRechitsChamberPlus12;
      MuonSystem->cscRechitClusterNRechitChamberPlus13[MuonSystem->nCscRechitClusters] = tmp.nCscRechitsChamberPlus13;
      MuonSystem->cscRechitClusterNRechitChamberPlus21[MuonSystem->nCscRechitClusters] = tmp.nCscRechitsChamberPlus21;
      MuonSystem->cscRechitClusterNRechitChamberPlus22[MuonSystem->nCscRechitClusters] = tmp.nCscRechitsChamberPlus22;
      MuonSystem->cscRechitClusterNRechitChamberPlus31[MuonSystem->nCscRechitClusters] = tmp.nCscRechitsChamberPlus31;
      MuonSystem->cscRechitClusterNRechitChamberPlus32[MuonSystem->nCscRechitClusters] = tmp.nCscRechitsChamberPlus32;
      MuonSystem->cscRechitClusterNRechitChamberPlus41[MuonSystem->nCscRechitClusters] = tmp.nCscRechitsChamberPlus41;
      MuonSystem->cscRechitClusterNRechitChamberPlus42[MuonSystem->nCscRechitClusters] = tmp.nCscRechitsChamberPlus42;
      MuonSystem->cscRechitClusterNRechitChamberMinus11[MuonSystem->nCscRechitClusters] = tmp.nCscRechitsChamberMinus11;
      MuonSystem->cscRechitClusterNRechitChamberMinus12[MuonSystem->nCscRechitClusters] = tmp.nCscRechitsChamberMinus12;
      MuonSystem->cscRechitClusterNRechitChamberMinus13[MuonSystem->nCscRechitClusters] = tmp.nCscRechitsChamberMinus13;
      MuonSystem->cscRechitClusterNRechitChamberMinus21[MuonSystem->nCscRechitClusters] = tmp.nCscRechitsChamberMinus21;
      MuonSystem->cscRechitClusterNRechitChamberMinus22[MuonSystem->nCscRechitClusters] = tmp.nCscRechitsChamberMinus22;
      MuonSystem->cscRechitClusterNRechitChamberMinus31[MuonSystem->nCscRechitClusters] = tmp.nCscRechitsChamberMinus31;
      MuonSystem->cscRechitClusterNRechitChamberMinus32[MuonSystem->nCscRechitClusters] = tmp.nCscRechitsChamberMinus32;
      MuonSystem->cscRechitClusterNRechitChamberMinus41[MuonSystem->nCscRechitClusters] = tmp.nCscRechitsChamberMinus41;
      MuonSystem->cscRechitClusterNRechitChamberMinus42[MuonSystem->nCscRechitClusters] = tmp.nCscRechitsChamberMinus42;
      MuonSystem->cscRechitClusterMaxChamber[MuonSystem->nCscRechitClusters] = tmp.maxChamber;
      MuonSystem->cscRechitClusterMaxChamberRatio[MuonSystem->nCscRechitClusters] = 1.0 * tmp.maxChamberRechits / tmp.nhits;
      MuonSystem->cscRechitClusterNChamber[MuonSystem->nCscRechitClusters] = tmp.nChamber;
      MuonSystem->cscRechitClusterMaxStation[MuonSystem->nCscRechitClusters] = tmp.maxStation;
      MuonSystem->cscRechitClusterMaxStationRatio[MuonSystem->nCscRechitClusters] = 1.0 * tmp.maxStationRechits / tmp.nhits;

      MuonSystem->cscRechitClusterNStation10[MuonSystem->nCscRechitClusters] = tmp.nStation10;
      MuonSystem->cscRechitClusterAvgStation10[MuonSystem->nCscRechitClusters] = tmp.avgStation10;

      // Jet veto/ muon veto
      MuonSystem->cscRechitClusterJetVetoPt[MuonSystem->nCscRechitClusters] = 0.0;
      MuonSystem->cscRechitClusterJetVetoE[MuonSystem->nCscRechitClusters] = 0.0;
      MuonSystem->cscRechitClusterMuonVetoPt[MuonSystem->nCscRechitClusters] = 0.0;
      MuonSystem->cscRechitClusterMuonVetoE[MuonSystem->nCscRechitClusters] = 0.0;

      // jet veto
      for (int i = 0; i < nJets; i++)
      {
        if (fabs(jetEta[i] > 3.0))
          continue;
        if (RazorAnalyzerMerged::deltaR(jetEta[i], jetPhi[i], MuonSystem->cscRechitClusterEta[MuonSystem->nCscRechitClusters], MuonSystem->cscRechitClusterPhi[MuonSystem->nCscRechitClusters]) < 0.4 && jetPt[i] > MuonSystem->cscRechitClusterJetVetoPt[MuonSystem->nCscRechitClusters])
        {
          MuonSystem->cscRechitClusterJetVetoPt[MuonSystem->nCscRechitClusters] = jetPt[i];
          MuonSystem->cscRechitClusterJetVetoE[MuonSystem->nCscRechitClusters] = jetE[i];
          MuonSystem->cscRechitClusterJetVetoTightId[MuonSystem->nCscRechitClusters] = jetPassIDTight[i];
          MuonSystem->cscRechitClusterJetVetoLooseId[MuonSystem->nCscRechitClusters] = jetPassIDLoose[i];
        }
      }
      float min_deltaR = 15.;
      int index = 999;

      for (int i = 0; i < nMuons; i++)
      {
        if (fabs(muonEta[i] > 3.0))
          continue;
        float muonIso = Muon_pfRelIso04_all[i];
        if (RazorAnalyzerMerged::deltaR(muonEta[i], muonPhi[i], MuonSystem->cscRechitClusterEta[MuonSystem->nCscRechitClusters], MuonSystem->cscRechitClusterPhi[MuonSystem->nCscRechitClusters]) < 0.4 && muonPt[i] > MuonSystem->cscRechitClusterMuonVetoPt[MuonSystem->nCscRechitClusters])
        {
          MuonSystem->cscRechitClusterMuonVetoPt[MuonSystem->nCscRechitClusters] = muonPt[i];
          MuonSystem->cscRechitClusterMuonVetoE[MuonSystem->nCscRechitClusters] = muonE[i];
          MuonSystem->cscRechitClusterMuonVetoGlobal[MuonSystem->nCscRechitClusters] = muon_isGlobal[i];
          MuonSystem->cscRechitClusterMuonVetoLooseId[MuonSystem->nCscRechitClusters] = muonIsLoose[i];
        }
      }
      if (!isData)
      {
        // match to gen level LLP
        min_deltaR = 15.;
        index = 999;
        for (int j = 0; j < MuonSystem->nGLLP; j++)
        {

          double current_delta_r = RazorAnalyzerMerged::deltaR(MuonSystem->cscRechitClusterEta[MuonSystem->nCscRechitClusters], MuonSystem->cscRechitClusterPhi[MuonSystem->nCscRechitClusters], MuonSystem->gLLP_eta[j], MuonSystem->gLLP_phi[j]);
          if (current_delta_r < min_deltaR)
          {
            min_deltaR = current_delta_r;
            index = j;
          }
        }
        if (min_deltaR < 0.4)
          MuonSystem->cscRechitCluster_match_gLLP[MuonSystem->nCscRechitClusters] = true;
        else
          MuonSystem->cscRechitCluster_match_gLLP[MuonSystem->nCscRechitClusters] = false;

        MuonSystem->cscRechitCluster_match_gLLP_minDeltaR[MuonSystem->nCscRechitClusters] = min_deltaR;
        MuonSystem->cscRechitCluster_match_gLLP_index[MuonSystem->nCscRechitClusters] = index;
        MuonSystem->cscRechitCluster_match_gLLP_eta[MuonSystem->nCscRechitClusters] = MuonSystem->gLLP_eta[index];
        MuonSystem->cscRechitCluster_match_gLLP_phi[MuonSystem->nCscRechitClusters] = MuonSystem->gLLP_phi[index];
        MuonSystem->cscRechitCluster_match_gLLP_decay_r[MuonSystem->nCscRechitClusters] = MuonSystem->gLLP_decay_vertex_r[index];
        MuonSystem->cscRechitCluster_match_gLLP_decay_z[MuonSystem->nCscRechitClusters] = MuonSystem->gLLP_decay_vertex_z[index];
        MuonSystem->cscRechitCluster_match_gLLP_csc[MuonSystem->nCscRechitClusters] = MuonSystem->gLLP_csc[index];
        MuonSystem->cscRechitCluster_match_gLLP_dt[MuonSystem->nCscRechitClusters] = MuonSystem->gLLP_dt[index];
        MuonSystem->cscRechitCluster_match_gLLP_e[MuonSystem->nCscRechitClusters] = MuonSystem->gLLP_e[index];
      }

      // match to MB1 DT segments
      for (int i = 0; i < nDtSeg; i++)
      {
        if (RazorAnalyzerMerged::deltaR(dtSegEta[i], dtSegPhi[i], MuonSystem->cscRechitClusterEta[MuonSystem->nCscRechitClusters], MuonSystem->cscRechitClusterPhi[MuonSystem->nCscRechitClusters]) < 0.4)
        {
          MuonSystem->cscRechitCluster_match_dtSeg_0p4[MuonSystem->nCscRechitClusters]++;
          if (dtSegStation[i] == 1)
            MuonSystem->cscRechitCluster_match_MB1Seg_0p4[MuonSystem->nCscRechitClusters]++;
        }
      }

      // match to RPC hits in RE1/2
      for (int i = 0; i < nRpc; i++)
      {
        float rpcR = sqrt(rpcX[i] * rpcX[i] + rpcY[i] * rpcY[i]);
        if (RazorAnalyzerMerged::deltaR(rpcEta[i], rpcPhi[i], MuonSystem->cscRechitClusterEta[MuonSystem->nCscRechitClusters], MuonSystem->cscRechitClusterPhi[MuonSystem->nCscRechitClusters]) < 0.4)
        {
          if (rpcR < 461.0 && rpcR > 275 && abs(rpcZ[i]) > 663 && abs(rpcZ[i]) < 730)
            MuonSystem->cscRechitCluster_match_RE12_0p4[MuonSystem->nCscRechitClusters]++;
          if (rpcR < 470 && rpcR > 380 && abs(rpcZ[i]) < 661)
            MuonSystem->cscRechitCluster_match_RB1_0p4[MuonSystem->nCscRechitClusters]++;
        }
      }

      MuonSystem->cscRechitClusterMet_dPhi[MuonSystem->nCscRechitClusters] = RazorAnalyzerMerged::deltaPhi(MuonSystem->cscRechitClusterPhi[MuonSystem->nCscRechitClusters], MuonSystem->metPhi);
      MuonSystem->nCscRechitClusters++;
    }

    // DT cluster

    points.clear();
    // cout<<"here"<<endl;

    for (int i = 0; i < min(N_MAX_DTRECHITS, nDtRechits); i++)
    {
      Rechits p;

      p.phi = dtRechitCorrectPhi[i];
      p.eta = dtRechitCorrectEta[i];
      p.x = dtRechitCorrectX[i];
      p.y = dtRechitCorrectY[i];
      p.z = dtRechitCorrectZ[i];
      p.t = -1;
      p.twire = -1;
      p.station = dtRechitStation[i];
      p.chamber = dtRechitWheel[i];
      p.superlayer = dtRechitSuperLayer[i];
      p.wheel = dtRechitWheel[i];
      p.clusterID = UNCLASSIFIED;
      points.push_back(p);
    }
    // cout<<"here"<<endl;

    // Do DBSCAN Clustering
    int min_point_dt = 50;  // minimum number of segments to call it a cluster
    float epsilon_dt = 0.2; // cluster radius parameter
    CACluster ds_dtRechit(min_point_dt, epsilon_dt, points);
    ds_dtRechit.run();
    ds_dtRechit.clusterProperties();
    // ds_dtRechit.merge_clusters();
    // ds_dtRechit.clusterProperties();

    ds_dtRechit.sort_clusters();

    MuonSystem->nDtRechitClusters = 0;

    for (auto &tmp : ds_dtRechit.clusters)
    {

      // remove overlaps
      bool overlap = false;
      for (int i = 0; i < MuonSystem->nCscRechitClusters; i++)
      {
        if (RazorAnalyzerMerged::deltaR(MuonSystem->cscRechitClusterEta[i], MuonSystem->cscRechitClusterPhi[i], tmp.eta, tmp.phi) < 0.4)
          overlap = true;
      }
      if (overlap)
        continue;

      MuonSystem->dtRechitClusterX[MuonSystem->nDtRechitClusters] = tmp.x;
      MuonSystem->dtRechitClusterY[MuonSystem->nDtRechitClusters] = tmp.y;
      MuonSystem->dtRechitClusterZ[MuonSystem->nDtRechitClusters] = tmp.z;
      if (abs(tmp.z) < 126.8)
        MuonSystem->dtRechitClusterWheel[MuonSystem->nDtRechitClusters] = 0;
      else if (tmp.z > 126.8 && tmp.z < 395.4)
        MuonSystem->dtRechitClusterWheel[MuonSystem->nDtRechitClusters] = 1;
      else if (tmp.z < -126.8 && tmp.z > -395.4)
        MuonSystem->dtRechitClusterWheel[MuonSystem->nDtRechitClusters] = -1;
      else if (tmp.z < 0)
        MuonSystem->dtRechitClusterWheel[MuonSystem->nDtRechitClusters] = -2;
      else
        MuonSystem->dtRechitClusterWheel[MuonSystem->nDtRechitClusters] = 2;
      MuonSystem->dtRechitClusterEta[MuonSystem->nDtRechitClusters] = tmp.eta;
      MuonSystem->dtRechitClusterPhi[MuonSystem->nDtRechitClusters] = tmp.phi;
      MuonSystem->dtRechitClusterSize[MuonSystem->nDtRechitClusters] = tmp.nhits;
      MuonSystem->dtRechitClusternXY[MuonSystem->nDtRechitClusters] = tmp.nXY;
      MuonSystem->dtRechitClusternZ[MuonSystem->nDtRechitClusters] = tmp.nZ;
      MuonSystem->dtRechitClusterXSpread[MuonSystem->nDtRechitClusters] = tmp.XSpread;
      MuonSystem->dtRechitClusterYSpread[MuonSystem->nDtRechitClusters] = tmp.YSpread;
      MuonSystem->dtRechitClusterZSpread[MuonSystem->nDtRechitClusters] = tmp.ZSpread;
      MuonSystem->dtRechitClusterXYSpread[MuonSystem->nDtRechitClusters] = tmp.XYSpread;
      MuonSystem->dtRechitClusterRSpread[MuonSystem->nDtRechitClusters] = tmp.RSpread;
      MuonSystem->dtRechitClusterEtaPhiSpread[MuonSystem->nDtRechitClusters] = tmp.EtaPhiSpread;
      MuonSystem->dtRechitClusterEtaSpread[MuonSystem->nDtRechitClusters] = tmp.EtaSpread;
      MuonSystem->dtRechitClusterPhiSpread[MuonSystem->nDtRechitClusters] = tmp.PhiSpread;
      MuonSystem->dtRechitClusterDeltaRSpread[MuonSystem->nDtRechitClusters] = tmp.DeltaRSpread;
      MuonSystem->dtRechitClusterMajorAxis[MuonSystem->nDtRechitClusters] = tmp.MajorAxis;
      MuonSystem->dtRechitClusterMinorAxis[MuonSystem->nDtRechitClusters] = tmp.MinorAxis;
      MuonSystem->dtRechitClusterSkewX[MuonSystem->nDtRechitClusters] = tmp.SkewX;
      MuonSystem->dtRechitClusterSkewY[MuonSystem->nDtRechitClusters] = tmp.SkewY;
      MuonSystem->dtRechitClusterSkewZ[MuonSystem->nDtRechitClusters] = tmp.SkewZ;
      MuonSystem->dtRechitClusterKurtX[MuonSystem->nDtRechitClusters] = tmp.KurtX;
      MuonSystem->dtRechitClusterKurtY[MuonSystem->nDtRechitClusters] = tmp.KurtY;
      MuonSystem->dtRechitClusterKurtZ[MuonSystem->nDtRechitClusters] = tmp.KurtZ;

      unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
      default_random_engine generator(seed);

      // default_random_engine generator;
      uniform_real_distribution<double> distribution(0.0, 1.0);
      float prob = 0.03;
      for (int i = 0; i < 12; ++i)
      {
        if (distribution(generator) < prob)
          MuonSystem->dtRechitClusterNoiseHitStation1[MuonSystem->nDtRechitClusters]++;
      }
      for (int i = 0; i < 12; ++i)
      {
        if (distribution(generator) < prob)
          MuonSystem->dtRechitClusterNoiseHitStation2[MuonSystem->nDtRechitClusters]++;
      }
      for (int i = 0; i < 12; ++i)
      {
        if (distribution(generator) < prob)
          MuonSystem->dtRechitClusterNoiseHitStation3[MuonSystem->nDtRechitClusters]++;
      }
      for (int i = 0; i < 8; ++i)
      {
        if (distribution(generator) < prob)
          MuonSystem->dtRechitClusterNoiseHitStation4[MuonSystem->nDtRechitClusters]++;
      }

      MuonSystem->dtRechitClusterNoiseHit[MuonSystem->nDtRechitClusters] = MuonSystem->dtRechitClusterNoiseHitStation1[MuonSystem->nDtRechitClusters] +
                                                                           MuonSystem->dtRechitClusterNoiseHitStation2[MuonSystem->nDtRechitClusters] +
                                                                           MuonSystem->dtRechitClusterNoiseHitStation3[MuonSystem->nDtRechitClusters] +
                                                                           MuonSystem->dtRechitClusterNoiseHitStation4[MuonSystem->nDtRechitClusters];

      MuonSystem->dtRechitClusterNHitStation1[MuonSystem->nDtRechitClusters] = tmp.nDtRechitsStation1;
      MuonSystem->dtRechitClusterNHitStation2[MuonSystem->nDtRechitClusters] = tmp.nDtRechitsStation2;
      MuonSystem->dtRechitClusterNHitStation3[MuonSystem->nDtRechitClusters] = tmp.nDtRechitsStation3;
      MuonSystem->dtRechitClusterNHitStation4[MuonSystem->nDtRechitClusters] = tmp.nDtRechitsStation4;

      MuonSystem->dtRechitClusterNHitWheel0[MuonSystem->nDtRechitClusters] = tmp.nDtRechitsWheel0;
      MuonSystem->dtRechitClusterNHitWheel1[MuonSystem->nDtRechitClusters] = tmp.nDtRechitsWheel1;
      MuonSystem->dtRechitClusterNHitWheel2[MuonSystem->nDtRechitClusters] = tmp.nDtRechitsWheel2;

      MuonSystem->dtRechitClusterMaxChamber[MuonSystem->nDtRechitClusters] = tmp.maxChamber;
      MuonSystem->dtRechitClusterMaxChamberRatio[MuonSystem->nDtRechitClusters] = 1.0 * tmp.maxChamberRechits / tmp.nhits;
      MuonSystem->dtRechitClusterNChamber[MuonSystem->nDtRechitClusters] = tmp.nChamber;
      MuonSystem->dtRechitClusterMaxStation[MuonSystem->nDtRechitClusters] = tmp.maxStation;
      MuonSystem->dtRechitClusterMaxStationRatio[MuonSystem->nDtRechitClusters] = 1.0 * tmp.maxStationRechits / tmp.nhits;
      MuonSystem->dtRechitClusterNStation10[MuonSystem->nDtRechitClusters] = tmp.nStation10;
      MuonSystem->dtRechitClusterAvgStation10[MuonSystem->nDtRechitClusters] = tmp.avgStation10;

      // Jet veto/ muon veto
      MuonSystem->dtRechitClusterJetVetoPt[MuonSystem->nDtRechitClusters] = 0.0;
      MuonSystem->dtRechitClusterJetVetoE[MuonSystem->nDtRechitClusters] = 0.0;
      MuonSystem->dtRechitClusterMuonVetoPt[MuonSystem->nDtRechitClusters] = 0.0;
      MuonSystem->dtRechitClusterMuonVetoE[MuonSystem->nDtRechitClusters] = 0.0;

      // jet veto
      for (int i = 0; i < nJets; i++)
      {
        if (fabs(jetEta[i] > 3.0))
          continue;
        if (RazorAnalyzerMerged::deltaR(jetEta[i], jetPhi[i], MuonSystem->dtRechitClusterEta[MuonSystem->nDtRechitClusters], MuonSystem->dtRechitClusterPhi[MuonSystem->nDtRechitClusters]) < 0.4 && jetPt[i] > MuonSystem->dtRechitClusterJetVetoPt[MuonSystem->nDtRechitClusters])
        {
          MuonSystem->dtRechitClusterJetVetoPt[MuonSystem->nDtRechitClusters] = jetPt[i];
          MuonSystem->dtRechitClusterJetVetoE[MuonSystem->nDtRechitClusters] = jetE[i];
          MuonSystem->dtRechitClusterJetVetoLooseId[MuonSystem->nDtRechitClusters] = jetPassIDLoose[i];
          MuonSystem->dtRechitClusterJetVetoTightId[MuonSystem->nDtRechitClusters] = jetPassIDTight[i];
        }
      }

      for (int i = 0; i < nMuons; i++)
      {
        if (fabs(muonEta[i] > 3.0))
          continue;
        float muonIso = Muon_pfRelIso04_all[i];
        if (RazorAnalyzerMerged::deltaR(muonEta[i], muonPhi[i], MuonSystem->dtRechitClusterEta[MuonSystem->nDtRechitClusters], MuonSystem->dtRechitClusterPhi[MuonSystem->nDtRechitClusters]) < 0.4 && muonPt[i] > MuonSystem->dtRechitClusterMuonVetoPt[MuonSystem->nDtRechitClusters])
        {
          MuonSystem->dtRechitClusterMuonVetoPt[MuonSystem->nDtRechitClusters] = muonPt[i];
          MuonSystem->dtRechitClusterMuonVetoE[MuonSystem->nDtRechitClusters] = muonE[i];
          MuonSystem->dtRechitClusterMuonVetoGlobal[MuonSystem->nDtRechitClusters] = muon_isGlobal[i];
          MuonSystem->dtRechitClusterMuonVetoLooseId[MuonSystem->nDtRechitClusters] = muonIsLoose[i];
        }
      }

      for (int i = 0; i < nDtSeg; i++)
      {
        if (RazorAnalyzerMerged::deltaR(dtSegEta[i], dtSegPhi[i], MuonSystem->dtRechitClusterEta[MuonSystem->nDtRechitClusters], MuonSystem->dtRechitClusterPhi[MuonSystem->nDtRechitClusters]) < 0.4)
        {
          if (dtSegStation[i] == 1)
            MuonSystem->dtRechitClusterNSegStation1[MuonSystem->nDtRechitClusters] += 1;
          if (dtSegStation[i] == 2)
            MuonSystem->dtRechitClusterNSegStation2[MuonSystem->nDtRechitClusters] += 1;
          if (dtSegStation[i] == 3)
            MuonSystem->dtRechitClusterNSegStation3[MuonSystem->nDtRechitClusters] += 1;
          if (dtSegStation[i] == 4)
            MuonSystem->dtRechitClusterNSegStation4[MuonSystem->nDtRechitClusters] += 1;
        }
        if (abs(RazorAnalyzerMerged::deltaPhi(dtSegPhi[i], MuonSystem->dtRechitClusterPhi[MuonSystem->nDtRechitClusters])) > 2)
        {
          if (dtSegStation[i] == 1)
            MuonSystem->dtRechitClusterNOppositeSegStation1[MuonSystem->nDtRechitClusters] += 1;
          if (dtSegStation[i] == 2)
            MuonSystem->dtRechitClusterNOppositeSegStation2[MuonSystem->nDtRechitClusters] += 1;
          if (dtSegStation[i] == 3)
            MuonSystem->dtRechitClusterNOppositeSegStation3[MuonSystem->nDtRechitClusters] += 1;
          if (dtSegStation[i] == 4)
            MuonSystem->dtRechitClusterNOppositeSegStation4[MuonSystem->nDtRechitClusters] += 1;
        }
      }

      // match to gen-level LLP
      float min_deltaR = 15.;
      int index = 999;
      if (!isData)
      {
        for (int j = 0; j < MuonSystem->nGLLP; j++)
        {
          double current_delta_r = RazorAnalyzerMerged::deltaR(tmp.eta, tmp.phi, MuonSystem->gLLP_eta[j], MuonSystem->gLLP_phi[j]);
          if (current_delta_r < min_deltaR)
          {
            min_deltaR = current_delta_r;
            index = j;
          }
        }
        if (min_deltaR < 0.4)
          MuonSystem->dtRechitCluster_match_gLLP[MuonSystem->nDtRechitClusters] = true;
        else
          MuonSystem->dtRechitCluster_match_gLLP[MuonSystem->nDtRechitClusters] = false;

        MuonSystem->dtRechitCluster_match_gLLP_minDeltaR[MuonSystem->nDtRechitClusters] = min_deltaR;
        MuonSystem->dtRechitCluster_match_gLLP_index[MuonSystem->nDtRechitClusters] = index;
        MuonSystem->dtRechitCluster_match_gLLP_eta[MuonSystem->nDtRechitClusters] = MuonSystem->gLLP_eta[index];
        MuonSystem->dtRechitCluster_match_gLLP_phi[MuonSystem->nDtRechitClusters] = MuonSystem->gLLP_phi[index];
        MuonSystem->dtRechitCluster_match_gLLP_decay_r[MuonSystem->nDtRechitClusters] = MuonSystem->gLLP_decay_vertex_r[index];
        MuonSystem->dtRechitCluster_match_gLLP_decay_z[MuonSystem->nDtRechitClusters] = MuonSystem->gLLP_decay_vertex_z[index];
        MuonSystem->dtRechitCluster_match_gLLP_csc[MuonSystem->nDtRechitClusters] = MuonSystem->gLLP_csc[index];
        MuonSystem->dtRechitCluster_match_gLLP_dt[MuonSystem->nDtRechitClusters] = MuonSystem->gLLP_dt[index];
        MuonSystem->dtRechitCluster_match_gLLP_e[MuonSystem->nDtRechitClusters] = MuonSystem->gLLP_e[index];
      }

      // match to MB1 DT segments
      MuonSystem->nCscRechits = ncscRechits;

      for (int i = 0; i < min(N_MAX_DTRECHITS, nDtRechits); i++)
      {
        if (RazorAnalyzerMerged::deltaR(dtRechitCorrectEta[i], dtRechitCorrectPhi[i], MuonSystem->dtRechitClusterEta[MuonSystem->nDtRechitClusters], MuonSystem->dtRechitClusterPhi[MuonSystem->nDtRechitClusters]) < 0.5)
        {
          if (dtRechitStation[i] == 1)
            MuonSystem->dtRechitCluster_match_MB1hits_0p5[MuonSystem->nDtRechitClusters]++;
        }
        if (RazorAnalyzerMerged::deltaR(dtRechitCorrectEta[i], dtRechitCorrectPhi[i], MuonSystem->dtRechitClusterEta[MuonSystem->nDtRechitClusters], MuonSystem->dtRechitClusterPhi[MuonSystem->nDtRechitClusters]) < 0.4)
        {
          if (dtRechitStation[i] == 1)
            MuonSystem->dtRechitCluster_match_MB1hits_0p4[MuonSystem->nDtRechitClusters]++;
        }
        if (abs(dtRechitWheel[i] - MuonSystem->dtRechitClusterWheel[MuonSystem->nDtRechitClusters]) == 1 && dtRechitStation[i] == 1)
        {
          if (abs(RazorAnalyzerMerged::deltaPhi(dtRechitCorrectPhi[i], MuonSystem->dtRechitClusterPhi[MuonSystem->nDtRechitClusters])) < TMath::Pi() / 4.0)
          {
            if (dtRechitWheel[i] - MuonSystem->dtRechitClusterWheel[MuonSystem->nDtRechitClusters] == 1)
              MuonSystem->dtRechitCluster_match_MB1hits_cosmics_plus[MuonSystem->nDtRechitClusters]++;
            else
              MuonSystem->dtRechitCluster_match_MB1hits_cosmics_minus[MuonSystem->nDtRechitClusters]++;
          }
        }
      }

      std::vector<int> dtRechitCluster_match_rpcBx;

      // match to RPC hits with dPhi<0.5 and same wheel in DT
      for (int i = 0; i < nRpc; i++)
      {
        float rpcR = sqrt(rpcX[i] * rpcX[i] + rpcY[i] * rpcY[i]);
        if (rpcRegion[i] != 0)
          continue;
        if (abs(RazorAnalyzerMerged::deltaPhi(rpcPhi[i], MuonSystem->dtRechitClusterPhi[MuonSystem->nDtRechitClusters])) < 0.5)
        {
          if (rpcRing[i] == MuonSystem->dtRechitClusterWheel[MuonSystem->nDtRechitClusters])
          {
            dtRechitCluster_match_rpcBx.push_back(rpcBx[i]);
            MuonSystem->dtRechitCluster_match_RPChits_dPhi0p5[MuonSystem->nDtRechitClusters]++;
            if (rpcR < 470 && rpcR > 380 && abs(rpcZ[i]) < 661)
              MuonSystem->dtRechitCluster_match_RB1_dPhi0p5[MuonSystem->nDtRechitClusters]++;
          }
        }
        if (RazorAnalyzerMerged::deltaR(rpcEta[i], rpcPhi[i], MuonSystem->dtRechitClusterEta[MuonSystem->nDtRechitClusters], MuonSystem->dtRechitClusterPhi[MuonSystem->nDtRechitClusters]) < 0.4)
        {
          if (rpcR < 470 && rpcR > 380 && abs(rpcZ[i]) < 661)
            MuonSystem->dtRechitCluster_match_RB1_0p4[MuonSystem->nDtRechitClusters]++;
        }
      }
      int max_occurence = 0;
      int max_bx = -999;
      for (unsigned int l = 0; l < dtRechitCluster_match_rpcBx.size(); l++)
      {
        int counter = 0;
        for (unsigned int j = 0; j < dtRechitCluster_match_rpcBx.size(); j++)
        {
          if (dtRechitCluster_match_rpcBx[j] == dtRechitCluster_match_rpcBx[l])
            counter++;
        }
        if (counter > max_occurence)
        {
          max_occurence = counter;
          max_bx = dtRechitCluster_match_rpcBx[l];
        }
      }
      MuonSystem->dtRechitCluster_match_RPCBx_dPhi0p5[MuonSystem->nDtRechitClusters] = max_bx;

      MuonSystem->dtRechitClusterMet_dPhi[MuonSystem->nDtRechitClusters] = RazorAnalyzerMerged::deltaPhi(MuonSystem->dtRechitClusterPhi[MuonSystem->nDtRechitClusters], MuonSystem->metPhi);

      MuonSystem->nDtRechitClusters++;
    }

    if (isData && MuonSystem->nDtRechitClusters + MuonSystem->nCscRechitClusters < 1)
      continue;

    if (!isData && signalScan)
    {
      pair<int, int> smsPair = make_pair(MuonSystem->mX, MuonSystem->ctau);
      Trees2D[smsPair]->Fill();
    }
    else
    {
      MuonSystem->tree_->Fill();
    }
  }
  if (!isData && signalScan)
  {
    for (auto &filePtr : Files2D)
    {
      cout << "Writing output tree (" << filePtr.second->GetName() << ")" << endl;
      filePtr.second->cd();
      Trees2D[filePtr.first]->Write();
      NEvents2D[filePtr.first]->Write("NEvents");
      Total2D[filePtr.first]->Write("Total");
      accep2D[filePtr.first]->Write("acceptance");
      accep_met2D[filePtr.first]->Write("acceptance_met");
      filePtr.second->Close();
    }
  }
  else if (!isData)
  {
    cout << "Filled Total of " << NEvents->GetBinContent(1) << " Events\n";
    cout << "Writing output trees..." << endl;
    outFile->cd();
    MuonSystem->tree_->Write();
    NEvents->Write();
    accep->Write("acceptance");
    accep_met->Write("acceptance_met");
    outFile->Close();
  }

  else
  {
    cout << "Filled Total of " << NEvents->GetBinContent(1) << " Events\n";
    cout << "Writing output trees..." << endl;
    outFile->cd();
    MuonSystem->tree_->Write();
    Nmet200->Write();
    NmetFilter->Write();
    Nlep0->Write();
    Njet1->Write();
    NcosmicVeto->Write();
    NEvents->Write();
    // outFile->Write();
    outFile->Close();
  }
}

void llp_MuonSystem_CAM::fillMetFilter(TreeMuonSystemCombination *MuonSystem)
{
  MuonSystem->Flag_goodVertices = Flag_goodVertices;
  MuonSystem->Flag_globalSuperTightHalo2016Filter = Flag_globalSuperTightHalo2016Filter;
  MuonSystem->Flag_EcalDeadCellTriggerPrimitiveFilter = Flag_EcalDeadCellTriggerPrimitiveFilter;
  MuonSystem->Flag_BadPFMuonFilter = Flag_BadPFMuonFilter;
  MuonSystem->Flag_BadPFMuonDzFilter = Flag_BadPFMuonDzFilter;
  MuonSystem->Flag_hfNoisyHitsFilter = Flag_hfNoisyHitsFilter;
  MuonSystem->Flag_eeBadScFilter = Flag_eeBadScFilter;

  bool flag_all = Flag_globalTightHalo2016Filter &&
                  Flag_HBHENoiseFilter &&
                  Flag_HBHENoiseIsoFilter &&
                  Flag_EcalDeadCellTriggerPrimitiveFilter &&
                  Flag_BadPFMuonFilter &&
                  Flag_BadChargedCandidateFilter &&
                  Flag_eeBadScFilter &&
                  Flag_ecalBadCalibFilter &&
                  Flag_hfNoisyHitsFilter;
  MuonSystem->Flag_all = flag_all;

  MuonSystem->Flag_ecalBadCalibFilter = Flag_ecalBadCalibFilter;
  MuonSystem->Flag_HBHENoiseFilter = Flag_HBHENoiseFilter;
  MuonSystem->Flag_HBHEIsoNoiseFilter = false; // Non-exist artifact
  MuonSystem->Flag_CSCTightHaloFilter = Flag_CSCTightHaloFilter;
  MuonSystem->Flag_BadChargedCandidateFilter = Flag_BadChargedCandidateFilter;
}

void llp_MuonSystem_CAM::fillHLT(TreeMuonSystemCombination *MuonSystem)
{
  // Triggers
  auto arr = MuonSystem->HLTDecision;
  arr[0] = HLT_AK8PFJet360_TrimMass30;
  arr[1] = HLT_AK8PFJet380_TrimMass30;
  arr[2] = HLT_AK8PFJet400_TrimMass30;
  arr[3] = HLT_AK8PFJet420_TrimMass30;
  arr[4] = HLT_AK8PFJet400_MassSD30;
  arr[5] = HLT_AK8PFJet420_MassSD30;
  arr[6] = HLT_AK8PFJet450_MassSD30;
  arr[7] = HLT_AK8DiPFJet250_250_MassSD30;
  arr[8] = HLT_AK8DiPFJet250_250_MassSD50;
  arr[9] = HLT_AK8DiPFJet260_260_MassSD30;
  arr[10] = HLT_AK8DiPFJet270_270_MassSD30;
  arr[11] = HLT_AK8PFHT750_TrimMass50;
  arr[12] = HLT_AK8PFHT800_TrimMass50;
  arr[13] = HLT_AK8PFHT850_TrimMass50;
  arr[14] = HLT_AK8PFHT900_TrimMass50;
  arr[15] = HLT_CaloJet500_NoJetID;
  arr[16] = HLT_CaloJet550_NoJetID;
  arr[17] = HLT_DoubleMu5_Upsilon_DoubleEle3_CaloIdL_TrackIdL;
  arr[18] = HLT_DoubleMu3_DoubleEle7p5_CaloIdL_TrackIdL_Upsilon;
  arr[19] = HLT_Trimuon5_3p5_2_Upsilon_Muon;
  arr[20] = HLT_TrimuonOpen_5_3p5_2_Upsilon_Muon;
  arr[21] = HLT_DoubleEle25_CaloIdL_MW;
  arr[22] = HLT_DoubleEle27_CaloIdL_MW;
  arr[23] = HLT_DoubleEle33_CaloIdL_MW;
  arr[24] = HLT_DoubleEle24_eta2p1_WPTight_Gsf;
  arr[25] = HLT_DoubleEle8_CaloIdM_TrackIdM_Mass8_DZ_PFHT350;
  arr[26] = HLT_DoubleEle8_CaloIdM_TrackIdM_Mass8_PFHT350;
  arr[27] = HLT_Mu27_Ele37_CaloIdL_MW;
  arr[28] = HLT_Mu37_Ele27_CaloIdL_MW;
  arr[29] = HLT_Mu37_TkMu27;
  arr[30] = HLT_DoubleMu4_3_Bs;
  arr[31] = HLT_DoubleMu4_3_Jpsi;
  arr[32] = HLT_DoubleMu4_3_LowMass;
  arr[33] = HLT_DoubleMu4_LowMass_Displaced;
  arr[34] = HLT_Mu0_L1DoubleMu;
  arr[35] = HLT_Mu4_L1DoubleMu;
  arr[36] = HLT_DoubleMu4_3_Photon4_BsToMMG;
  arr[37] = HLT_DoubleMu4_3_Displaced_Photon4_BsToMMG;
  arr[38] = HLT_DoubleMu3_Trk_Tau3mu;
  arr[39] = HLT_DoubleMu3_TkMu_DsTau3Mu;
  arr[40] = HLT_DoubleMu4_Mass3p8_DZ_PFHT350;
  arr[41] = HLT_DoubleMu4_MuMuTrk_Displaced;
  arr[42] = HLT_Mu3_PFJet40;
  arr[43] = HLT_Mu7p5_L2Mu2_Jpsi;
  arr[44] = HLT_Mu7p5_L2Mu2_Upsilon;
  arr[45] = HLT_Mu3_L1SingleMu5orSingleMu7;
  arr[46] = HLT_DoublePhoton33_CaloIdL;
  arr[47] = HLT_DoublePhoton70;
  arr[48] = HLT_DoublePhoton85;
  arr[49] = HLT_Ele15_WPLoose_Gsf;
  arr[50] = HLT_Ele20_WPLoose_Gsf;
  arr[51] = HLT_DiEle27_WPTightCaloOnly_L1DoubleEG;
  arr[52] = HLT_Ele27_WPTight_Gsf;
  arr[53] = HLT_Ele28_WPTight_Gsf;
  arr[54] = HLT_Ele30_WPTight_Gsf;
  arr[55] = HLT_Ele32_WPTight_Gsf;
  arr[56] = HLT_Ele35_WPTight_Gsf;
  arr[57] = HLT_Ele35_WPTight_Gsf_L1EGMT;
  arr[58] = HLT_Ele38_WPTight_Gsf;
  arr[59] = HLT_Ele40_WPTight_Gsf;
  arr[60] = HLT_Ele32_WPTight_Gsf_L1DoubleEG;
  arr[61] = HLT_HT300_Beamspot;
  arr[62] = HLT_ZeroBias_Beamspot;
  arr[63] = HLT_IsoMu20_eta2p1_TightChargedIsoPFTauHPS27_eta2p1_TightID_CrossL1;
  arr[64] = HLT_IsoMu27_MediumDeepTauPFTauHPS20_eta2p1_SingleL1;
  arr[65] = HLT_IsoMu20;
  arr[66] = HLT_IsoMu24;
  arr[67] = HLT_IsoMu24_eta2p1;
  arr[68] = HLT_IsoMu27;
  arr[69] = HLT_UncorrectedJetE30_NoBPTX;
  arr[70] = HLT_UncorrectedJetE30_NoBPTX3BX;
  arr[71] = HLT_UncorrectedJetE60_NoBPTX3BX;
  arr[72] = HLT_UncorrectedJetE70_NoBPTX3BX;
  arr[73] = HLT_L1SingleMu18;
  arr[74] = HLT_L1SingleMu25;
  arr[75] = HLT_L1SingleMuCosmics;
  arr[76] = HLT_L2Mu10_NoVertex_NoBPTX3BX;
  arr[77] = HLT_L2Mu10_NoVertex_NoBPTX;
  arr[78] = HLT_L2Mu45_NoVertex_3Sta_NoBPTX3BX;
  arr[79] = HLT_L2Mu40_NoVertex_3Sta_NoBPTX3BX;
  arr[80] = HLT_L2Mu23NoVtx_2Cha;
  arr[81] = HLT_L2Mu23NoVtx_2Cha_CosmicSeed;
  arr[82] = HLT_DoubleL2Mu30NoVtx_2Cha_CosmicSeed_Eta2p4;
  arr[83] = HLT_DoubleL2Mu30NoVtx_2Cha_Eta2p4;
  arr[84] = HLT_DoubleL2Mu50;
  arr[85] = HLT_DoubleL2Mu23NoVtx_2Cha_CosmicSeed;
  arr[86] = HLT_DoubleL2Mu25NoVtx_2Cha_CosmicSeed;
  arr[87] = HLT_DoubleL2Mu25NoVtx_2Cha_CosmicSeed_Eta2p4;
  arr[88] = HLT_DoubleL2Mu23NoVtx_2Cha;
  arr[89] = HLT_DoubleL2Mu25NoVtx_2Cha;
  arr[90] = HLT_DoubleL2Mu25NoVtx_2Cha_Eta2p4;
  arr[91] = HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL;
  arr[92] = HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL;
  arr[93] = HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ;
  arr[94] = HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ;
  arr[95] = HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8;
  arr[96] = HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_Mass8;
  arr[97] = HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8;
  arr[98] = HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_Mass3p8;
  arr[99] = HLT_Mu25_TkMu0_Onia;
  arr[100] = HLT_Mu30_TkMu0_Psi;
  arr[101] = HLT_Mu30_TkMu0_Upsilon;
  arr[102] = HLT_Mu20_TkMu0_Phi;
  arr[103] = HLT_Mu25_TkMu0_Phi;
  arr[104] = HLT_Mu15;
  arr[105] = HLT_Mu20;
  arr[106] = HLT_Mu27;
  arr[107] = HLT_Mu50;
  arr[108] = HLT_Mu55;
  arr[109] = HLT_CascadeMu100;
  arr[110] = HLT_HighPtTkMu100;
  arr[111] = HLT_DiPFJetAve40;
  arr[112] = HLT_DiPFJetAve60;
  arr[113] = HLT_DiPFJetAve80;
  arr[114] = HLT_DiPFJetAve140;
  arr[115] = HLT_DiPFJetAve200;
  arr[116] = HLT_DiPFJetAve260;
  arr[117] = HLT_DiPFJetAve320;
  arr[118] = HLT_DiPFJetAve400;
  arr[119] = HLT_DiPFJetAve500;
  arr[120] = HLT_DiPFJetAve60_HFJEC;
  arr[121] = HLT_DiPFJetAve80_HFJEC;
  arr[122] = HLT_DiPFJetAve100_HFJEC;
  arr[123] = HLT_DiPFJetAve160_HFJEC;
  arr[124] = HLT_DiPFJetAve220_HFJEC;
  arr[125] = HLT_DiPFJetAve300_HFJEC;
  arr[126] = HLT_AK8PFJet40;
  arr[127] = HLT_AK8PFJet60;
  arr[128] = HLT_AK8PFJet80;
  arr[129] = HLT_AK8PFJet140;
  arr[130] = HLT_AK8PFJet200;
  arr[131] = HLT_AK8PFJet260;
  arr[132] = HLT_AK8PFJet320;
  arr[133] = HLT_AK8PFJet400;
  arr[134] = HLT_AK8PFJet450;
  arr[135] = HLT_AK8PFJet500;
  arr[136] = HLT_AK8PFJet550;
  arr[137] = HLT_PFJet40;
  arr[138] = HLT_PFJet60;
  arr[139] = HLT_PFJet80;
  arr[140] = HLT_PFJet110;
  arr[141] = HLT_PFJet140;
  arr[142] = HLT_PFJet200;
  arr[143] = HLT_PFJet260;
  arr[144] = HLT_PFJet320;
  arr[145] = HLT_PFJet400;
  arr[146] = HLT_PFJet450;
  arr[147] = HLT_PFJet500;
  arr[148] = HLT_PFJet550;
  arr[149] = HLT_PFJetFwd15;
  arr[150] = HLT_PFJetFwd25;
  arr[151] = HLT_PFJetFwd40;
  arr[152] = HLT_PFJetFwd60;
  arr[153] = HLT_PFJetFwd80;
  arr[154] = HLT_PFJetFwd140;
  arr[155] = HLT_PFJetFwd200;
  arr[156] = HLT_PFJetFwd260;
  arr[157] = HLT_PFJetFwd320;
  arr[158] = HLT_PFJetFwd400;
  arr[159] = HLT_PFJetFwd450;
  arr[160] = HLT_PFJetFwd500;
  arr[161] = HLT_AK8PFJetFwd15;
  arr[162] = HLT_AK8PFJetFwd25;
  arr[163] = HLT_AK8PFJetFwd40;
  arr[164] = HLT_AK8PFJetFwd60;
  arr[165] = HLT_AK8PFJetFwd80;
  arr[166] = HLT_AK8PFJetFwd140;
  arr[167] = HLT_AK8PFJetFwd200;
  arr[168] = HLT_AK8PFJetFwd260;
  arr[169] = HLT_AK8PFJetFwd320;
  arr[170] = HLT_AK8PFJetFwd400;
  arr[171] = HLT_AK8PFJetFwd450;
  arr[172] = HLT_AK8PFJetFwd500;
  arr[173] = HLT_PFHT180;
  arr[174] = HLT_PFHT250;
  arr[175] = HLT_PFHT370;
  arr[176] = HLT_PFHT430;
  arr[177] = HLT_PFHT510;
  arr[178] = HLT_PFHT590;
  arr[179] = HLT_PFHT680;
  arr[180] = HLT_PFHT780;
  arr[181] = HLT_PFHT890;
  arr[182] = HLT_PFHT1050;
  arr[183] = HLT_PFHT500_PFMET100_PFMHT100_IDTight;
  arr[184] = HLT_PFHT500_PFMET110_PFMHT110_IDTight;
  arr[185] = HLT_PFHT700_PFMET85_PFMHT85_IDTight;
  arr[186] = HLT_PFHT800_PFMET75_PFMHT75_IDTight;
  arr[187] = HLT_PFMET110_PFMHT110_IDTight;
  arr[188] = HLT_PFMET120_PFMHT120_IDTight;
  arr[189] = HLT_PFMET130_PFMHT130_IDTight;
  arr[190] = HLT_PFMET140_PFMHT140_IDTight;
  arr[191] = HLT_PFMET120_PFMHT120_IDTight_PFHT60;
  arr[192] = HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60;
  arr[193] = HLT_PFMETTypeOne120_PFMHT120_IDTight_PFHT60;
  arr[194] = HLT_PFMETTypeOne110_PFMHT110_IDTight;
  arr[195] = HLT_PFMETTypeOne120_PFMHT120_IDTight;
  arr[196] = HLT_PFMETTypeOne130_PFMHT130_IDTight;
  arr[197] = HLT_PFMETTypeOne140_PFMHT140_IDTight;
  arr[198] = HLT_PFMETNoMu110_PFMHTNoMu110_IDTight;
  arr[199] = HLT_PFMETNoMu120_PFMHTNoMu120_IDTight;
  arr[200] = HLT_PFMETNoMu130_PFMHTNoMu130_IDTight;
  arr[201] = HLT_PFMETNoMu140_PFMHTNoMu140_IDTight;
  arr[202] = HLT_PFMETNoMu110_PFMHTNoMu110_IDTight_FilterHF;
  arr[203] = HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_FilterHF;
  arr[204] = HLT_PFMETNoMu130_PFMHTNoMu130_IDTight_FilterHF;
  arr[205] = HLT_PFMETNoMu140_PFMHTNoMu140_IDTight_FilterHF;
  arr[206] = HLT_L1ETMHadSeeds;
  arr[207] = HLT_CaloMHT90;
  arr[208] = HLT_CaloMET90_NotCleaned;
  arr[209] = HLT_CaloMET350_NotCleaned;
  arr[210] = HLT_PFMET200_NotCleaned;
  arr[211] = HLT_PFMET250_NotCleaned;
  arr[212] = HLT_PFMET300_NotCleaned;
  arr[213] = HLT_PFMET200_BeamHaloCleaned;
  arr[214] = HLT_PFMETTypeOne200_BeamHaloCleaned;
  arr[215] = HLT_MET105_IsoTrk50;
  arr[216] = HLT_MET120_IsoTrk50;
  arr[217] = HLT_SingleJet30_Mu12_SinglePFJet40;
  arr[218] = HLT_Mu12eta2p3;
  arr[219] = HLT_Mu12eta2p3_PFJet40;
  arr[220] = HLT_Mu12_DoublePFJets40_PFBTagDeepCSV_p71;
  arr[221] = HLT_Mu12_DoublePFJets100_PFBTagDeepCSV_p71;
  arr[222] = HLT_Mu12_DoublePFJets200_PFBTagDeepCSV_p71;
  arr[223] = HLT_Mu12_DoublePFJets350_PFBTagDeepCSV_p71;
  arr[224] = HLT_Mu12_DoublePFJets40MaxDeta1p6_DoublePFBTagDeepCSV_p71;
  arr[225] = HLT_Mu12_DoublePFJets54MaxDeta1p6_DoublePFBTagDeepCSV_p71;
  arr[226] = HLT_DoublePFJets40_PFBTagDeepCSV_p71;
  arr[227] = HLT_DoublePFJets100_PFBTagDeepCSV_p71;
  arr[228] = HLT_DoublePFJets200_PFBTagDeepCSV_p71;
  arr[229] = HLT_DoublePFJets350_PFBTagDeepCSV_p71;
  arr[230] = HLT_DoublePFJets116MaxDeta1p6_DoublePFBTagDeepCSV_p71;
  arr[231] = HLT_DoublePFJets128MaxDeta1p6_DoublePFBTagDeepCSV_p71;
  arr[232] = HLT_Mu12_DoublePFJets40_PFBTagDeepJet_p71;
  arr[233] = HLT_Mu12_DoublePFJets100_PFBTagDeepJet_p71;
  arr[234] = HLT_Mu12_DoublePFJets200_PFBTagDeepJet_p71;
  arr[235] = HLT_Mu12_DoublePFJets350_PFBTagDeepJet_p71;
  arr[236] = HLT_Mu12_DoublePFJets40MaxDeta1p6_DoublePFBTagDeepJet_p71;
  arr[237] = HLT_Mu12_DoublePFJets54MaxDeta1p6_DoublePFBTagDeepJet_p71;
  arr[238] = HLT_DoublePFJets40_PFBTagDeepJet_p71;
  arr[239] = HLT_DoublePFJets100_PFBTagDeepJet_p71;
  arr[240] = HLT_DoublePFJets200_PFBTagDeepJet_p71;
  arr[241] = HLT_DoublePFJets350_PFBTagDeepJet_p71;
  arr[242] = HLT_DoublePFJets116MaxDeta1p6_DoublePFBTagDeepJet_p71;
  arr[243] = HLT_DoublePFJets128MaxDeta1p6_DoublePFBTagDeepJet_p71;
  arr[244] = HLT_Photon300_NoHE;
  arr[245] = HLT_Mu8_TrkIsoVVL;
  arr[246] = HLT_Mu8_DiEle12_CaloIdL_TrackIdL_DZ;
  arr[247] = HLT_Mu8_DiEle12_CaloIdL_TrackIdL;
  arr[248] = HLT_Mu8_Ele8_CaloIdM_TrackIdM_Mass8_PFHT350_DZ;
  arr[249] = HLT_Mu8_Ele8_CaloIdM_TrackIdM_Mass8_PFHT350;
  arr[250] = HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ;
  arr[251] = HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_PFDiJet30;
  arr[252] = HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_CaloDiJet30;
  arr[253] = HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_PFDiJet30_PFBtagDeepCSV_1p5;
  arr[254] = HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_CaloDiJet30_CaloBtagDeepCSV_1p5;
  arr[255] = HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL;
  arr[256] = HLT_Mu17_TrkIsoVVL;
  arr[257] = HLT_Mu19_TrkIsoVVL;
  arr[258] = HLT_BTagMu_AK4DiJet20_Mu5;
  arr[259] = HLT_BTagMu_AK4DiJet40_Mu5;
  arr[260] = HLT_BTagMu_AK4DiJet70_Mu5;
  arr[261] = HLT_BTagMu_AK4DiJet110_Mu5;
  arr[262] = HLT_BTagMu_AK4DiJet170_Mu5;
  arr[263] = HLT_BTagMu_AK4Jet300_Mu5;
  arr[264] = HLT_BTagMu_AK8DiJet170_Mu5;
  arr[265] = HLT_BTagMu_AK8Jet170_DoubleMu5;
  arr[266] = HLT_BTagMu_AK8Jet300_Mu5;
  arr[267] = HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ;
  arr[268] = HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL;
  arr[269] = HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ;
  arr[270] = HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL;
  arr[271] = HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL;
  arr[272] = HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ;
  arr[273] = HLT_Photon20;
  arr[274] = HLT_Photon33;
  arr[275] = HLT_Photon50;
  arr[276] = HLT_Photon75;
  arr[277] = HLT_Photon90;
  arr[278] = HLT_Photon120;
  arr[279] = HLT_Photon150;
  arr[280] = HLT_Photon175;
  arr[281] = HLT_Photon200;
  arr[282] = HLT_Photon30EB_TightID_TightIso;
  arr[283] = HLT_Photon110EB_TightID_TightIso;
  arr[284] = HLT_Photon100EBHE10;
  arr[285] = HLT_Photon50_R9Id90_HE10_IsoM;
  arr[286] = HLT_Photon75_R9Id90_HE10_IsoM;
  arr[287] = HLT_Photon90_R9Id90_HE10_IsoM;
  arr[288] = HLT_Photon120_R9Id90_HE10_IsoM;
  arr[289] = HLT_Photon165_R9Id90_HE10_IsoM;
  arr[290] = HLT_Diphoton30_22_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass90;
  arr[291] = HLT_Diphoton30_22_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass95;
  arr[292] = HLT_Photon35_TwoProngs35;
  arr[293] = HLT_IsoMu24_TwoProngs35;
  arr[294] = HLT_Dimuon0_Jpsi_L1_NoOS;
  arr[295] = HLT_Dimuon0_Jpsi_NoVertexing_NoOS;
  arr[296] = HLT_Dimuon0_Jpsi;
  arr[297] = HLT_Dimuon0_Jpsi_NoVertexing;
  arr[298] = HLT_Dimuon0_Jpsi_L1_4R_0er1p5R;
  arr[299] = HLT_Dimuon0_Jpsi_NoVertexing_L1_4R_0er1p5R;
  arr[300] = HLT_Dimuon0_Jpsi3p5_Muon2;
  arr[301] = HLT_Dimuon0_Upsilon_L1_4p5;
  arr[302] = HLT_Dimuon0_Upsilon_L1_5;
  arr[303] = HLT_Dimuon0_Upsilon_L1_4p5NoOS;
  arr[304] = HLT_Dimuon0_Upsilon_L1_4p5er2p0;
  arr[305] = HLT_Dimuon0_Upsilon_L1_4p5er2p0M;
  arr[306] = HLT_Dimuon0_Upsilon_NoVertexing;
  arr[307] = HLT_Dimuon0_Upsilon_L1_5M;
  arr[308] = HLT_Dimuon0_LowMass_L1_0er1p5R;
  arr[309] = HLT_Dimuon0_LowMass_L1_0er1p5;
  arr[310] = HLT_Dimuon0_LowMass;
  arr[311] = HLT_Dimuon0_LowMass_L1_4;
  arr[312] = HLT_Dimuon0_LowMass_L1_4R;
  arr[313] = HLT_Dimuon0_LowMass_L1_TM530;
  arr[314] = HLT_Dimuon0_Upsilon_Muon_L1_TM0;
  arr[315] = HLT_Dimuon0_Upsilon_Muon_NoL1Mass;
  arr[316] = HLT_TripleMu_5_3_3_Mass3p8_DZ;
  arr[317] = HLT_TripleMu_10_5_5_DZ;
  arr[318] = HLT_TripleMu_12_10_5;
  arr[319] = HLT_Tau3Mu_Mu7_Mu1_TkMu1_Tau15;
  arr[320] = HLT_Tau3Mu_Mu7_Mu1_TkMu1_Tau15_Charge1;
  arr[321] = HLT_Tau3Mu_Mu7_Mu1_TkMu1_IsoTau15;
  arr[322] = HLT_Tau3Mu_Mu7_Mu1_TkMu1_IsoTau15_Charge1;
  arr[323] = HLT_DoubleMu3_DZ_PFMET50_PFMHT60;
  arr[324] = HLT_DoubleMu3_DZ_PFMET70_PFMHT70;
  arr[325] = HLT_DoubleMu3_DZ_PFMET90_PFMHT90;
  arr[326] = HLT_DoubleMu3_Trk_Tau3mu_NoL1Mass;
  arr[327] = HLT_DoubleMu4_Jpsi_Displaced;
  arr[328] = HLT_DoubleMu4_Jpsi_NoVertexing;
  arr[329] = HLT_DoubleMu4_JpsiTrkTrk_Displaced;
  arr[330] = HLT_DoubleMu4_JpsiTrk_Bc;
  arr[331] = HLT_DoubleMu43NoFiltersNoVtx;
  arr[332] = HLT_DoubleMu48NoFiltersNoVtx;
  arr[333] = HLT_Mu43NoFiltersNoVtx_Photon43_CaloIdL;
  arr[334] = HLT_Mu48NoFiltersNoVtx_Photon48_CaloIdL;
  arr[335] = HLT_Mu38NoFiltersNoVtxDisplaced_Photon38_CaloIdL;
  arr[336] = HLT_Mu43NoFiltersNoVtxDisplaced_Photon43_CaloIdL;
  arr[337] = HLT_HT425;
  arr[338] = HLT_HT430_DisplacedDijet40_DisplacedTrack;
  arr[339] = HLT_HT500_DisplacedDijet40_DisplacedTrack;
  arr[340] = HLT_HT430_DisplacedDijet60_DisplacedTrack;
  arr[341] = HLT_HT400_DisplacedDijet40_DisplacedTrack;
  arr[342] = HLT_HT650_DisplacedDijet60_Inclusive;
  arr[343] = HLT_HT550_DisplacedDijet60_Inclusive;
  arr[344] = HLT_DiJet110_35_Mjj650_PFMET110;
  arr[345] = HLT_DiJet110_35_Mjj650_PFMET120;
  arr[346] = HLT_DiJet110_35_Mjj650_PFMET130;
  arr[347] = HLT_TripleJet110_35_35_Mjj650_PFMET110;
  arr[348] = HLT_TripleJet110_35_35_Mjj650_PFMET120;
  arr[349] = HLT_TripleJet110_35_35_Mjj650_PFMET130;
  arr[350] = HLT_Ele30_eta2p1_WPTight_Gsf_CentralPFJet35_EleCleaned;
  arr[351] = HLT_Ele28_eta2p1_WPTight_Gsf_HT150;
  arr[352] = HLT_Ele28_HighEta_SC20_Mass55;
  arr[353] = HLT_Ele15_IsoVVVL_PFHT450_CaloBTagDeepCSV_4p5;
  arr[354] = HLT_Ele15_IsoVVVL_PFHT450_PFMET50;
  arr[355] = HLT_Ele15_IsoVVVL_PFHT450;
  arr[356] = HLT_Ele50_IsoVVVL_PFHT450;
  arr[357] = HLT_Ele15_IsoVVVL_PFHT600;
  arr[358] = HLT_Mu4_TrkIsoVVL_DiPFJet90_40_DEta3p5_MJJ750_HTT300_PFMETNoMu60;
  arr[359] = HLT_Mu8_TrkIsoVVL_DiPFJet40_DEta3p5_MJJ750_HTT300_PFMETNoMu60;
  arr[360] = HLT_Mu10_TrkIsoVVL_DiPFJet40_DEta3p5_MJJ750_HTT350_PFMETNoMu60;
  arr[361] = HLT_Mu15_IsoVVVL_PFHT450_CaloBTagDeepCSV_4p5;
  arr[362] = HLT_Mu15_IsoVVVL_PFHT450_PFMET50;
  arr[363] = HLT_Mu15_IsoVVVL_PFHT450;
  arr[364] = HLT_Mu50_IsoVVVL_PFHT450;
  arr[365] = HLT_Mu15_IsoVVVL_PFHT600;
  arr[366] = HLT_Mu3er1p5_PFJet100er2p5_PFMET70_PFMHT70_IDTight;
  arr[367] = HLT_Mu3er1p5_PFJet100er2p5_PFMET80_PFMHT80_IDTight;
  arr[368] = HLT_Mu3er1p5_PFJet100er2p5_PFMET90_PFMHT90_IDTight;
  arr[369] = HLT_Mu3er1p5_PFJet100er2p5_PFMET100_PFMHT100_IDTight;
  arr[370] = HLT_Mu3er1p5_PFJet100er2p5_PFMETNoMu70_PFMHTNoMu70_IDTight;
  arr[371] = HLT_Mu3er1p5_PFJet100er2p5_PFMETNoMu80_PFMHTNoMu80_IDTight;
  arr[372] = HLT_Mu3er1p5_PFJet100er2p5_PFMETNoMu90_PFMHTNoMu90_IDTight;
  arr[373] = HLT_Mu3er1p5_PFJet100er2p5_PFMETNoMu100_PFMHTNoMu100_IDTight;
  arr[374] = HLT_Dimuon10_PsiPrime_Barrel_Seagulls;
  arr[375] = HLT_Dimuon20_Jpsi_Barrel_Seagulls;
  arr[376] = HLT_Dimuon10_Upsilon_y1p4;
  arr[377] = HLT_Dimuon12_Upsilon_y1p4;
  arr[378] = HLT_Dimuon14_Phi_Barrel_Seagulls;
  arr[379] = HLT_Dimuon25_Jpsi;
  arr[380] = HLT_Dimuon14_PsiPrime;
  arr[381] = HLT_Dimuon14_PsiPrime_noCorrL1;
  arr[382] = HLT_Dimuon18_PsiPrime;
  arr[383] = HLT_Dimuon18_PsiPrime_noCorrL1;
  arr[384] = HLT_Dimuon24_Upsilon_noCorrL1;
  arr[385] = HLT_Dimuon24_Phi_noCorrL1;
  arr[386] = HLT_Dimuon25_Jpsi_noCorrL1;
  arr[387] = HLT_DiMu4_Ele9_CaloIdL_TrackIdL_DZ_Mass3p8;
  arr[388] = HLT_DiMu9_Ele9_CaloIdL_TrackIdL_DZ;
  arr[389] = HLT_DiMu9_Ele9_CaloIdL_TrackIdL;
  arr[390] = HLT_DoubleIsoMu20_eta2p1;
  arr[391] = HLT_TrkMu12_DoubleTrkMu5NoFiltersNoVtx;
  arr[392] = HLT_Mu8;
  arr[393] = HLT_Mu17;
  arr[394] = HLT_Mu19;
  arr[395] = HLT_Mu17_Photon30_IsoCaloId;
  arr[396] = HLT_Ele8_CaloIdL_TrackIdL_IsoVL_PFJet30;
  arr[397] = HLT_Ele12_CaloIdL_TrackIdL_IsoVL_PFJet30;
  arr[398] = HLT_Ele23_CaloIdL_TrackIdL_IsoVL_PFJet30;
  arr[399] = HLT_Ele8_CaloIdM_TrackIdM_PFJet30;
  arr[400] = HLT_Ele17_CaloIdM_TrackIdM_PFJet30;
  arr[401] = HLT_Ele23_CaloIdM_TrackIdM_PFJet30;
  arr[402] = HLT_Ele50_CaloIdVT_GsfTrkIdT_PFJet165;
  arr[403] = HLT_Ele115_CaloIdVT_GsfTrkIdT;
  arr[404] = HLT_Ele135_CaloIdVT_GsfTrkIdT;
  arr[405] = HLT_PFHT330PT30_QuadPFJet_75_60_45_40_TriplePFBTagDeepCSV_4p5;
  arr[406] = HLT_PFHT330PT30_QuadPFJet_75_60_45_40;
  arr[407] = HLT_PFHT400_SixPFJet32_DoublePFBTagDeepCSV_2p94;
  arr[408] = HLT_PFHT400_SixPFJet32;
  arr[409] = HLT_PFHT450_SixPFJet36_PFBTagDeepCSV_1p59;
  arr[410] = HLT_PFHT450_SixPFJet36;
  arr[411] = HLT_PFHT400_FivePFJet_100_100_60_30_30;
  arr[412] = HLT_PFHT400_FivePFJet_100_100_60_30_30_DoublePFBTagDeepCSV_4p5;
  arr[413] = HLT_PFHT400_FivePFJet_120_120_60_30_30_DoublePFBTagDeepCSV_4p5;
  arr[414] = HLT_PFHT350;
  arr[415] = HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_PFHT350MinPFJet15;
  arr[416] = HLT_ECALHT800;
  arr[417] = HLT_DiSC30_18_EIso_AND_HE_Mass70;
  arr[418] = HLT_Physics;
  arr[419] = HLT_Random;
  arr[420] = HLT_ZeroBias;
  arr[421] = HLT_ZeroBias_Alignment;
  arr[422] = HLT_Photon20_HoverELoose;
  arr[423] = HLT_Photon30_HoverELoose;
  arr[424] = HLT_EcalCalibration;
  arr[425] = HLT_HcalCalibration;
  arr[426] = HLT_L1UnpairedBunchBptxMinus;
  arr[427] = HLT_L1UnpairedBunchBptxPlus;
  arr[428] = HLT_L1NotBptxOR;
  arr[429] = HLT_L1_CDC_SingleMu_3_er1p2_TOP120_DPHI2p618_3p142;
  arr[430] = HLT_CDC_L2cosmic_10_er1p0;
  arr[431] = HLT_CDC_L2cosmic_5p5_er1p0;
  arr[432] = HLT_HcalNZS;
  arr[433] = HLT_HcalPhiSym;
  arr[434] = HLT_HcalIsolatedbunch;
  arr[435] = HLT_IsoTrackHB;
  arr[436] = HLT_IsoTrackHE;
  arr[437] = HLT_ZeroBias_FirstCollisionAfterAbortGap;
  arr[438] = HLT_ZeroBias_IsolatedBunches;
  arr[439] = HLT_ZeroBias_FirstCollisionInTrain;
  arr[440] = HLT_ZeroBias_LastCollisionInTrain;
  arr[441] = HLT_ZeroBias_FirstBXAfterTrain;
  arr[442] = HLT_Ele16_Ele12_Ele8_CaloIdL_TrackIdL;
  arr[443] = HLT_VBF_DoubleMediumDeepTauPFTauHPS20_eta2p1;
  arr[444] = HLT_Photon75_R9Id90_HE10_IsoM_EBOnly_PFJetsMJJ300DEta3;
  arr[445] = HLT_PFMET100_PFMHT100_IDTight_PFHT60;
  arr[446] = HLT_PFMETNoMu100_PFMHTNoMu100_IDTight_PFHT60;
  arr[447] = HLT_PFMETTypeOne100_PFMHT100_IDTight_PFHT60;
  arr[448] = HLT_Mu18_Mu9_SameSign;
  arr[449] = HLT_DoubleMu2_Jpsi_DoubleTrk1_Phi1p05;
  arr[450] = HLT_DoubleMu3_DCA_PFMET50_PFMHT60;
  arr[451] = HLT_TripleMu_5_3_3_Mass3p8_DCA;
  arr[452] = HLT_QuadPFJet103_88_75_15_DoublePFBTagDeepCSV_1p3_7p7_VBF1;
  arr[453] = HLT_QuadPFJet111_90_80_15_DoublePFBTagDeepCSV_1p3_7p7_VBF1;
  arr[454] = HLT_QuadPFJet103_88_75_15_PFBTagDeepCSV_1p3_VBF2;
  arr[455] = HLT_QuadPFJet105_88_76_15_PFBTagDeepCSV_1p3_VBF2;
  arr[456] = HLT_QuadPFJet111_90_80_15_PFBTagDeepCSV_1p3_VBF2;
  arr[457] = HLT_QuadPFJet103_88_75_15;
  arr[458] = HLT_QuadPFJet105_88_76_15;
  arr[459] = HLT_QuadPFJet111_90_80_15;
  arr[460] = HLT_AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_p02;
  arr[461] = HLT_AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_np2;
  arr[462] = HLT_AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_np4;
  arr[463] = HLT_Diphoton30_18_R9IdL_AND_HE_AND_IsoCaloId_Mass55;
  arr[464] = HLT_Diphoton30_18_R9IdL_AND_HE_AND_IsoCaloId;
  arr[465] = HLT_Mu12_IP6;
  arr[466] = HLT_QuadPFJet105_88_76_15_DoublePFBTagDeepCSV_1p3_7p7_VBF1;
  arr[467] = HLT_DoubleMediumDeepTauPFTauHPS35_L2NN_eta2p1;
  arr[468] = HLT_Ele24_eta2p1_WPTight_Gsf_LooseDeepTauPFTauHPS30_eta2p1_CrossL1;
  arr[469] = HLT_IsoMu20_eta2p1_LooseDeepTauPFTauHPS27_eta2p1_CrossL1;
  arr[470] = HLT_IsoMu24_eta2p1_LooseDeepTauPFTauHPS180_eta2p1;
  arr[471] = HLT_IsoMu24_eta2p1_LooseDeepTauPFTauHPS30_eta2p1_CrossL1;
  arr[472] = HLT_IsoMu24_eta2p1_MediumDeepTauPFTauHPS35_L2NN_eta2p1_CrossL1;
  arr[473] = HLT_LooseDeepTauPFTauHPS180_L2NN_eta2p1;
  arr[474] = HLT_PFHT330PT30_QuadPFJet_75_60_45_40_TriplePFBTagDeepJet_4p5;
  arr[475] = HLT_PFHT400_FivePFJet_100_100_60_30_30_DoublePFBTagDeepJet_4p5;
  arr[476] = HLT_PFHT400_FivePFJet_120_120_60_30_30_DoublePFBTagDeepJet_4p5;
  arr[477] = HLT_PFHT400_SixPFJet32_DoublePFBTagDeepJet_2p94;
  arr[478] = HLT_PFHT450_SixPFJet36_PFBTagDeepJet_1p59;
  arr[479] = HLT_QuadPFJet103_88_75_15_DoublePFBTagDeepJet_1p3_7p7_VBF1;
  arr[480] = HLT_QuadPFJet103_88_75_15_PFBTagDeepJet_1p3_VBF2;
  arr[481] = HLT_QuadPFJet105_88_76_15_DoublePFBTagDeepJet_1p3_7p7_VBF1;
  arr[482] = HLT_QuadPFJet105_88_76_15_PFBTagDeepJet_1p3_VBF2;
  arr[483] = HLT_QuadPFJet111_90_80_15_DoublePFBTagDeepJet_1p3_7p7_VBF1;
  arr[484] = HLT_QuadPFJet111_90_80_15_PFBTagDeepJet_1p3_VBF2;
  arr[485] = HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_PFDiJet30_PFBtagDeepJet_1p5;
  arr[486] = HLT_QuadPFJet70_50_40_30;
  arr[487] = HLT_QuadPFJet70_50_40_30_PFBTagParticleNet_2BTagSum0p65;
  arr[488] = HLT_QuadPFJet70_50_40_35_PFBTagParticleNet_2BTagSum0p65;
  arr[489] = HLT_QuadPFJet70_50_45_35_PFBTagParticleNet_2BTagSum0p65;
  arr[490] = HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_PFDiJet30_PFBTagParticleNet_2BTagSum0p65;
  arr[491] = HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_QuadPFJet70_50_40_30;
  arr[492] = HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_QuadPFJet70_50_40_30_PFBTagParticleNet_2BTagSum0p65;
  arr[493] = HLT_AK8PFJet230_SoftDropMass40;
  arr[494] = HLT_AK8PFJet230_SoftDropMass40_PFAK8ParticleNetBB0p35;
  arr[495] = HLT_AK8PFJet250_SoftDropMass40_PFAK8ParticleNetBB0p35;
  arr[496] = HLT_AK8PFJet275_SoftDropMass40_PFAK8ParticleNetBB0p35;
  arr[497] = HLT_AK8PFJet400_SoftDropMass40;
  arr[498] = HLT_AK8PFJet425_SoftDropMass40;
  arr[499] = HLT_AK8PFJet450_SoftDropMass40;
  arr[500] = HLT_AK8PFJet230_SoftDropMass40_PFAK8ParticleNetTauTau0p30;
  arr[501] = HLT_AK8PFJet250_SoftDropMass40_PFAK8ParticleNetTauTau0p30;
  arr[502] = HLT_AK8PFJet275_SoftDropMass40_PFAK8ParticleNetTauTau0p30;
  arr[503] = HLT_IsoMu50_AK8PFJet230_SoftDropMass40;
  arr[504] = HLT_IsoMu50_AK8PFJet230_SoftDropMass40_PFAK8ParticleNetBB0p35;
  arr[505] = HLT_Ele50_CaloIdVT_GsfTrkIdT_AK8PFJet230_SoftDropMass40;
  arr[506] = HLT_Ele50_CaloIdVT_GsfTrkIdT_AK8PFJet230_SoftDropMass40_PFAK8ParticleNetBB0p35;
  arr[507] = HLT_DoubleMediumDeepTauPFTauHPS30_L2NN_eta2p1_PFJet60;
  arr[508] = HLT_DoubleMediumDeepTauPFTauHPS30_L2NN_eta2p1_PFJet75;
  arr[509] = HLT_IsoMu24_eta2p1_MediumDeepTauPFTauHPS30_L2NN_eta2p1_CrossL1;
  arr[510] = HLT_IsoMu24_eta2p1_MediumDeepTauPFTauHPS30_L2NN_eta2p1_PFJet60_CrossL1;
  arr[511] = HLT_IsoMu24_eta2p1_MediumDeepTauPFTauHPS30_L2NN_eta2p1_PFJet75_CrossL1;
  arr[512] = HLT_DoublePFJets40_Mass500_MediumDeepTauPFTauHPS45_L2NN_MediumDeepTauPFTauHPS20_eta2p1;
  arr[513] = HLT_IsoMu24_eta2p1_MediumDeepTauPFTauHPS20_eta2p1_SingleL1;
  arr[514] = HLT_IsoMu24_eta2p1_MediumDeepTauPFTauHPS45_L2NN_eta2p1_CrossL1;
  arr[515] = HLT_DoubleL2Mu10NoVtx_2Cha_VetoL3Mu0DxyMax1cm;
  arr[516] = HLT_DoubleL2Mu12NoVtx_2Cha_VetoL3Mu0DxyMax1cm;
  arr[517] = HLT_DoubleL2Mu14NoVtx_2Cha_VetoL3Mu0DxyMax1cm;
  arr[518] = HLT_DoubleL3Mu16_10NoVtx_DxyMin0p01cm;
  arr[519] = HLT_DoubleL3Mu18_10NoVtx_DxyMin0p01cm;
  arr[520] = HLT_DoubleL3Mu20_10NoVtx_DxyMin0p01cm;
  arr[521] = HLT_L2Mu10NoVtx_2Cha;
  arr[522] = HLT_L2Mu10NoVtx_2Cha_VetoL3Mu0DxyMax1cm;
  arr[523] = HLT_L3Mu10NoVtx;
  arr[524] = HLT_L3Mu10NoVtx_DxyMin0p01cm;
  arr[525] = HLT_DoubleL2Mu_L3Mu16NoVtx_VetoL3Mu0DxyMax0p1cm;
  arr[526] = HLT_DoubleL2Mu_L3Mu18NoVtx_VetoL3Mu0DxyMax0p1cm;
  arr[527] = HLT_DoubleL2Mu10NoVtx_2Cha_CosmicSeed_VetoL3Mu0DxyMax1cm;
  arr[528] = HLT_DoubleL2Mu12NoVtx_2Cha_CosmicSeed_VetoL3Mu0DxyMax1cm;
  arr[529] = HLT_L2Mu10NoVtx_2Cha_CosmicSeed;
  arr[530] = HLT_L2Mu10NoVtx_2Cha_CosmicSeed_VetoL3Mu0DxyMax1cm;
  arr[531] = HLT_DoubleL3dTksMu16_10NoVtx_DxyMin0p01cm;
  arr[532] = HLT_L3dTksMu10_NoVtx_DxyMin0p01cm;
  arr[533] = HLT_Mu20NoFiltersNoVtxDisplaced_Photon20_CaloCustomId;
  arr[534] = HLT_DoubleMediumChargedIsoDisplacedPFTauHPS32_Trk1_eta2p1;
  arr[535] = HLT_HT430_DelayedJet40_DoubleDelay0p5nsTrackless;
  arr[536] = HLT_HT430_DelayedJet40_DoubleDelay1nsInclusive;
  arr[537] = HLT_HT430_DelayedJet40_SingleDelay1nsTrackless;
  arr[538] = HLT_HT430_DelayedJet40_SingleDelay2nsInclusive;
  arr[539] = HLT_L1Mu6HT240;
  arr[540] = HLT_Mu6HT240_DisplacedDijet30_Inclusive0PtrkShortSig5;
  arr[541] = HLT_Mu6HT240_DisplacedDijet30_Inclusive1PtrkShortSig5_DisplacedLoose;
  arr[542] = HLT_Mu6HT240_DisplacedDijet35_Inclusive0PtrkShortSig5;
  arr[543] = HLT_Mu6HT240_DisplacedDijet35_Inclusive1PtrkShortSig5_DisplacedLoose;
  arr[544] = HLT_Mu6HT240_DisplacedDijet40_Inclusive0PtrkShortSig5;
  arr[545] = HLT_Mu6HT240_DisplacedDijet40_Inclusive1PtrkShortSig5_DisplacedLoose;
  arr[546] = HLT_HT430_DisplacedDijet30_Inclusive1PtrkShortSig5;
  arr[547] = HLT_HT430_DisplacedDijet35_Inclusive1PtrkShortSig5;
  arr[548] = HLT_HT430_DisplacedDijet40_Inclusive1PtrkShortSig5;
  arr[549] = HLT_CaloMET60_DTCluster50;
  arr[550] = HLT_CaloMET60_DTClusterNoMB1S50;
  arr[551] = HLT_L1MET_DTCluster50;
  arr[552] = HLT_L1MET_DTClusterNoMB1S50;
  arr[553] = HLT_CscCluster_Loose;
  arr[554] = HLT_CscCluster_Medium;
  arr[555] = HLT_CscCluster_Tight;
  arr[556] = HLT_L1CSCShower_DTCluster50;
  arr[557] = HLT_L1CSCShower_DTCluster75;
  arr[558] = HLT_PFMET105_IsoTrk50;
  arr[559] = HLT_PFMET110_PFJet100;
  arr[560] = HLT_HT170_L1SingleLLPJet_DisplacedDijet40_DisplacedTrack;
  arr[561] = HLT_HT200_L1SingleLLPJet_DisplacedDijet40_DisplacedTrack;
  arr[562] = HLT_HT200_L1SingleLLPJet_DisplacedDijet60_DisplacedTrack;
  arr[563] = HLT_HT270_L1SingleLLPJet_DisplacedDijet40_DisplacedTrack;
  arr[564] = HLT_HT320_L1SingleLLPJet_DisplacedDijet60_Inclusive;
  arr[565] = HLT_HT420_L1SingleLLPJet_DisplacedDijet60_Inclusive;
  arr[566] = HLT_HT200_L1SingleLLPJet_DelayedJet40_SingleDelay1nsTrackless;
  arr[567] = HLT_HT200_L1SingleLLPJet_DelayedJet40_SingleDelay2nsInclusive;
  arr[568] = HLT_HT200_L1SingleLLPJet_DelayedJet40_DoubleDelay0p5nsTrackless;
  arr[569] = HLT_HT200_L1SingleLLPJet_DelayedJet40_DoubleDelay1nsInclusive;
  arr[570] = HLT_HT200_L1SingleLLPJet_DisplacedDijet30_Inclusive1PtrkShortSig5;
  arr[571] = HLT_HT200_L1SingleLLPJet_DisplacedDijet35_Inclusive1PtrkShortSig5;
  arr[572] = HLT_HT200_L1SingleLLPJet_DisplacedDijet40_Inclusive1PtrkShortSig5;
  arr[573] = HLT_DiPhoton10Time1p4ns;
  arr[574] = HLT_DiPhoton10Time1p6ns;
  arr[575] = HLT_DiPhoton10Time1p8ns;
  arr[576] = HLT_DiPhoton10Time2ns;
  arr[577] = HLT_DiPhoton10sminlt0p1;
  arr[578] = HLT_DiPhoton10sminlt0p12;
  arr[579] = HLT_DiPhoton10_CaloIdL;
  arr[580] = HLT_DoubleEle4_eta1p22_mMax6;
  arr[581] = HLT_DoubleEle4p5_eta1p22_mMax6;
  arr[582] = HLT_DoubleEle5_eta1p22_mMax6;
  arr[583] = HLT_DoubleEle5p5_eta1p22_mMax6;
  arr[584] = HLT_DoubleEle6_eta1p22_mMax6;
  arr[585] = HLT_DoubleEle6p5_eta1p22_mMax6;
  arr[586] = HLT_DoubleEle7_eta1p22_mMax6;
  arr[587] = HLT_DoubleEle7p5_eta1p22_mMax6;
  arr[588] = HLT_DoubleEle8_eta1p22_mMax6;
  arr[589] = HLT_DoubleEle8p5_eta1p22_mMax6;
  arr[590] = HLT_DoubleEle9_eta1p22_mMax6;
  arr[591] = HLT_DoubleEle9p5_eta1p22_mMax6;
  arr[592] = HLT_DoubleEle10_eta1p22_mMax6;
  arr[593] = HLT_Diphoton20_14_eta1p5_R9IdL_AND_HE_AND_IsoTCaloIdT;
  arr[594] = HLT_Diphoton20_14_eta1p5_R9IdL_AND_HET_AND_IsoTCaloIdT;
  arr[595] = HLT_Diphoton22_14_eta1p5_R9IdL_AND_HE_AND_IsoTCaloIdT;
  arr[596] = HLT_Diphoton22_14_eta1p5_R9IdL_AND_HET_AND_IsoTCaloIdT;
  arr[597] = HLT_ExpressMuons;
  arr[598] = HLT_OnlineMonitorGroup;
  arr[599] = HLT_PPSMaxTracksPerArm1;
  arr[600] = HLT_PPSMaxTracksPerRP4;
  arr[601] = HLT_EphemeralPhysics;
  arr[602] = HLT_EphemeralZeroBias;
}
