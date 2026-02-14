/* #region: includes and usings */
#include "llp_MuonSystem_CA_mdsnano.h"
#include "RazorHelper.h"
#include "TreeMuonSystemCombination.h"

#include "CACluster.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include "TMath.h"
#include "Math/Vector3D.h"
#include <initializer_list>
#include <iostream>
#include <memory>
#include <random>

//C++ includes
#include "assert.h"
//ROOT includes
#include "TH1F.h"

//using namespace fastjet;
using namespace std::chrono;
using namespace std;
using namespace ROOT::Math;
/* #endregion */

/* #region: Usefull templates/functions/structs  */
// Defining a helper for making single binned histograms.
std::unique_ptr<TH1F> makeCounterHist(const std::string& name) {
  auto h = std::make_unique<TH1F>(name.c_str(), name.c_str(), 1, 1, 2);
  h->SetDirectory(nullptr); // important: avoid ROOT owning/deleting it
  return h;
}

// Some rechit fields are optional and/or renamed across ntuple versions.
// They are read into temporary buffers here (instead of relying only on the
// auto-generated merged_event bindings) and copied later only if present.
// Defining a helper that binds branches not included in merged_event
// Presence flags prevent copying garbage when a branch is missing; output
// values then remain at the TreeMuonSystemCombination InitTree defaults.
template <typename T>
bool bindOptionalBranch(TTree* chain, T* address, std::initializer_list<const char*> aliases) {
  for (const char* branchName : aliases) {
    if (chain->GetBranch(branchName)) {
      chain->SetBranchAddress(branchName, address);
      return true;
    }
  }
  return false;
}

// Defining leptons and jets as structures
struct leptons {
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
struct jets {
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

//lepton highest pt comparator
struct largest_pt {
  inline bool operator()(const leptons& p1, const leptons& p2) {
    return p1.lepton.Pt() > p2.lepton.Pt();
  }
} my_largest_pt;

//jet highest pt comparator
struct largest_pt_jet {
  inline bool operator()(const jets& p1, const jets& p2) {
    return p1.jet.Pt() > p2.jet.Pt();
  }
} my_largest_pt_jet;
/* #endregion */

void llp_MuonSystem_CA_mdsnano::Analyze(bool isData, int options, string outputfilename, string analysisTag) {

  /* #region: fChain check and nEntries declaration.*/
  if (fChain == nullptr) {
    cout << "fChain is a null pointer." << endl;
    return;
  };
  const Long64_t nEntries = fChain->GetEntries();
  cout << "Initializing..." << endl;
  cout << "IsData = " << isData << "\n";
  cout << "options = " << options << "\n";
  cout << "Total Events: " << nEntries << "\n";
  /* #endregion */
  
  /* #region: options format: MH/MX/ctau/condor: 1000/300/0/1 */
  // mh can be 3-4 digits, mx is always 3 digits, ctau is one digit(number of zeros), last digit is condor option
  // mh can be 3-4 digits, mx is always 3 digits, ctau is 2 digit(number of zeros), last digit is condor option

  // int mx = int(options/1000)%1000;
  // int mh = options/1000000;
  // int ctau = pow(10, int(options/10)%10) * int(int(options/100)%10);

  // cout<<"mh "<<mh<<", mx "<<mx<<", ctau "<<ctau<<endl;
  /* #endregion */

  /* #region: Startup decode/defaults for this analyzer invocation.*/
  const float ELE_MASS = 0.000511;
  const float MU_MASS = 0.105658;
  const float Z_MASS = 91.2;
  bool signalScan = int(options / 10) == 1;
  int option = options % 10;
  if (analysisTag.empty()) analysisTag = "Razor2016_80X";
  /* #endregion */

  /* #region: declares lookup tables for signalScan*/
  map<pair<int, int>, std::unique_ptr<TFile>> Files2D;
  map<pair<int, int>, TTree*> Trees2D;
  map<pair<int, int>, TH1F*> NEvents2D;
  map<pair<int, int>, TH1F*> accep2D;
  map<pair<int, int>, TH1F*> accep_met2D;
  map<pair<int, int>, TH1F*> Total2D;
  /* #endregion */

  /* #region: run-mode logging*/
  std::cout << "[INFO]: running on " << (isData ? "data" : "MC")
            << " with option: " << option << '\n';
  std::cout << "[INFO]: running "
            << (signalScan ? "with Signal scan" : "without Signal scan ")
            << (signalScan ? "" : std::to_string(option))
            << '\n';
  /* #endregion */

  /* #region: setting up output file */
  string outfilename = outputfilename;
  if (outfilename == "") outfilename = "MuonSystem_Tree.root";
  
  TFile* outFile;
  if (isData || !signalScan) outFile = new TFile(outfilename.c_str(), "RECREATE");
  
  TreeMuonSystemCombination* MuonSystem = new TreeMuonSystemCombination;
  MuonSystem->CreateTree();
  MuonSystem->tree_->SetAutoFlush(0);
  MuonSystem->InitTree();
  /* #endregion */

  /* #region: histogram containing total number of processed events (for normalization)*/
  auto NEvents = makeCounterHist("NEvents");
  auto Total = makeCounterHist("Total");

  auto accep = makeCounterHist("accep");
  auto accep_csccsc = makeCounterHist("accep_csccsc");
  auto accep_cscdt = makeCounterHist("accep_cscdt");
  auto accep_met = makeCounterHist("accep_met");

  auto Nmet200 = makeCounterHist("Nmet200");
  auto NmetFilter = makeCounterHist("NmetFilter");
  auto Nlep0 = makeCounterHist("Nlep0");
  auto Njet1 = makeCounterHist("Njet1");
  auto NcosmicVeto = makeCounterHist("NcosmicVeto");
  /* #endregion */

  // initialize helper in memory safe way
  auto helper = std::make_unique<RazorHelper>(analysisTag, isData);

  /* #region: optional branches, Accept both snake_case and camelCase branch names for compatibility with  */
  Int_t cscRechitsIChamber_opt[20000] = {0};
  Int_t cscRechitsNStrips_opt[20000] = {0};
  Int_t cscRechitsWGroupsBX_opt[20000] = {0};
  Int_t cscRechitsHitWire_opt[20000] = {0};
  Int_t cscRechitsNWireGroups_opt[20000] = {0};
  Float_t cscRechitsE_opt[20000] = {0};
  Int_t dtRechitsSector_opt[20000] = {0};
  const auto hasCscRechitsIChamber = bindOptionalBranch(fChain, cscRechitsIChamber_opt, {"cscRechits_IChamber", "cscRechitsIChamber"});
  const auto hasCscRechitsNStrips = bindOptionalBranch(fChain, cscRechitsNStrips_opt, {"cscRechits_NStrips", "cscRechitsNStrips"});
  const auto hasCscRechitsWGroupsBX = bindOptionalBranch(fChain, cscRechitsWGroupsBX_opt, {"cscRechits_WGroupsBX", "cscRechitsWGroupsBX"});
  const auto hasCscRechitsHitWire = bindOptionalBranch(fChain, cscRechitsHitWire_opt, {"cscRechits_HitWire", "cscRechitsHitWire"});
  const auto hasCscRechitsNWireGroups = bindOptionalBranch(fChain, cscRechitsNWireGroups_opt, {"cscRechits_NWireGroups", "cscRechitsNWireGroups"});
  const auto hasCscRechitsE = bindOptionalBranch(fChain, cscRechitsE_opt, {"cscRechits_E", "cscRechitsE"});
  const auto hasDtRechitsSector = bindOptionalBranch(fChain, dtRechitsSector_opt, {"dtRecHits_Sector", "dtRecHitsSector", "dtRechitSector"});
  /* #endregion */
  
  // [1] Event loop
  /* #region */
  auto lastReport = steady_clock::now(); // mr. timer
  for (Long64_t jentry = 0; jentry < nEntries; jentry++) {

    // progress logging
    if (jentry % 1000 == 0) {
      const auto now = steady_clock::now();
      const double dt = duration<double>(now - lastReport).count();
      cout << "Processing entry " << jentry << " | dt=" << dt << " s\n";
      lastReport = now;
    }

    // [2] Read event entry and reset output event variables
    /* #region */
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0)
      break;
    fChain->GetEntry(jentry);

    //fill normalization histogram
    MuonSystem->InitVariables();
    /* #endregion */

    // [3] Build local aliases and per-event helper inputs
    /* #region */
    auto nCscSeg = ncscSegments;
    auto nDtSeg = ndtSegments;
    auto nDtRechits = ndtRecHits;
    auto nRpc = nrpcRecHits;
    auto cscRechitsQuality = cscRechits_Quality;
    auto cscRechitsPhi = cscRechits_Phi;
    auto cscRechitsEta = cscRechits_Eta;
    auto cscRechitsX = cscRechits_X;
    auto cscRechitsY = cscRechits_Y;
    auto cscRechitsZ = cscRechits_Z;
    auto cscRechitsTpeak = cscRechits_Tpeak;
    auto cscRechitsTwire = cscRechits_Twire;
    auto cscRechitsChamber = cscRechits_Chamber;
    auto cscRechitsStation = cscRechits_Station;

    auto dtRechitCorrectX = dtRecHits_X;
    auto dtRechitCorrectY = dtRecHits_Y;
    auto dtRechitCorrectZ = dtRecHits_Z;
    auto dtRechitCorrectEta = dtRecHits_Eta;
    auto dtRechitCorrectPhi = dtRecHits_Phi;
    auto dtRechitLayer = dtRecHits_Layer;
    auto dtRechitStation = dtRecHits_Station;
    auto dtRechitWheel = dtRecHits_Wheel;
    auto dtRechitSuperLayer = dtRecHits_SuperLayer;

    auto dtSegStation = dtSegments_Station;
    auto dtSegPhi = dtSegments_Phi;
    auto dtSegEta = dtSegments_Eta;
    auto rpcPhi = rpcRecHits_Phi;
    auto rpcEta = rpcRecHits_Eta;
    auto rpcX = rpcRecHits_X;
    auto rpcY = rpcRecHits_Y;
    auto rpcZ = rpcRecHits_Z;

    auto rpcBx = rpcRecHits_Bx;
    auto rpcRegion = rpcRecHits_Region;
    auto rpcRing = rpcRecHits_Ring;
    auto rpcLayer = rpcRecHits_Layer;
    auto rpcStation = rpcRecHits_Station;
    auto rpcSector = rpcRecHits_Sector;
    auto rpcTime = rpcRecHits_Time;
    auto rpcTimeError = rpcRecHits_TimeError;

    auto eleCharge = Electron_charge; // eleCharge
    auto eleEta = Electron_eta; // eleEta
    auto elePhi = Electron_phi; // elePhi
    auto elePt = Electron_pt; // elePt
    auto ele_dZ = Electron_dz; // ele_dZ

    
    bool ele_passCutBasedIDTight[10] = {0}; // ele_passCutBasedIDTight
    bool ele_passCutBasedIDVeto[10] = {0}; // ele_passCutBasedIDVeto

    for (int i = 0; i < nElectron; i++) {
      ele_passCutBasedIDTight[i] = Electron_cutBased[i] >= 4;
      ele_passCutBasedIDVeto[i] = Electron_cutBased[i] >= 1;
    }
    auto eventNum = event; // eventNum
    auto fixedGridRhoFastjetAll = Rho_fixedGridRhoFastjetAll; // fixedGridRhoFastjetAll

    Float_t jetE[150] = {0}; // jetE
    for (int i = 0; i < nJet; ++i) {
      auto eta = Jet_eta[i];
      auto pt = Jet_pt[i];
      auto pz = pt * TMath::SinH(eta);
      auto mass = Jet_mass[i];
      jetE[i] = TMath::Sqrt(mass * mass + pt * pt + pz * pz);
    }
    auto jetEta = Jet_eta; // jetEta
    bool jetPassIDTightLepVeto[150] = {0}; // jetPassIDLoose
    bool jetPassIDTight[150] = {0}; // jetPassIDTight
    for (int i = 0; i < nJet; ++i) {
      jetPassIDTight[i] = helper->jetTightLepVeto(analysisTag, false, Jet_neHEF[i], Jet_neEmEF[i], Jet_chEmEF[i], Jet_muEF[i], Jet_chHEF[i], Jet_chMultiplicity[i], Jet_neMultiplicity[i], Jet_eta[i], Jet_jetId[i]);
      jetPassIDTightLepVeto[i] = helper->jetTightLepVeto(analysisTag, true, Jet_neHEF[i], Jet_neEmEF[i], Jet_chEmEF[i], Jet_muEF[i], Jet_chHEF[i], Jet_chMultiplicity[i], Jet_neMultiplicity[i], Jet_eta[i], Jet_jetId[i]);
      // cout<<jetPassIDTight[i]<<","<<jetPassIDTightLepVeto[i]<<", "<< Jet_neHEF[i]<<", "<<Jet_neEmEF[i]<<", "<<Jet_chEmEF[i]<<", "<<Jet_muEF[i]<<", "<<Jet_chHEF[i]<<", "<<Jet_chMultiplicity[i]<<", "<<Jet_neMultiplicity[i]<<", "<< Jet_eta[i]<<", "<< Jet_jetId[i]<<endl;
    }
    auto jetPhi = Jet_phi; // jetPhi
    auto jetPt = Jet_pt; // jetPt
    auto lumiNum = luminosityBlock; // lumiNum
    auto metType1Phi = MET_phi; // metType1Phi
    auto metType1Pt = MET_pt; // metType1Pt
    if (analysisTag == "Summer24") {
      metType1Phi = PFMET_phi; // metType1Phi
      metType1Pt = PFMET_pt; // metType1Pt
    }

    auto muonCharge = Muon_charge; // muonCharge
    Float_t muonE[50] = {0}; // muonE
    for (int i = 0; i < nMuon; ++i) {
      auto eta = Muon_eta[i];
      auto pt = Muon_pt[i];
      auto pz = pt * TMath::SinH(eta);
      auto mass = MU_MASS;
      muonE[i] = TMath::Sqrt(mass * mass + pt * pt + pz * pz);
    }
    auto muonEta = Muon_eta; // muonEta
    auto muonIsLoose = Muon_looseId; // muonIsLoose
    auto muonIsTight = Muon_tightId; // muonIsTight
    auto muonPhi = Muon_phi; // muonPhi
    auto muonPt = Muon_pt; // muonPt
    auto muon_pfRelIso04_all = Muon_pfRelIso04_all; // muon_chargedIso
    auto muon_dZ = Muon_dz; // muon_dZ
    auto muon_isGlobal = Muon_isGlobal; // muon_isGlobal
    auto nElectrons = nElectron; // nElectrons
    auto nJets = nJet; // nJets
    auto nMuons = nMuon; // nMuons
    int nPV = PV_npvs; // nPV
    auto runNum = run; // runNum

    auto* lheComments = (string*)"123";
    auto genWeight = Generator_weight; // genWeight
    /* #endregion */

    // [4] Signal-scan routing (split output by model point)
    /* #region */
    if (!isData && signalScan) {
      string mh_substring = lheComments->substr(lheComments->find("MH-") + 3);
      int mh = stoi(mh_substring.substr(0, mh_substring.find('_')));
      string mx_substring = lheComments->substr(lheComments->find("MS-") + 3);
      int mx = stoi(mx_substring.substr(0, mx_substring.find('_')));
      string ctau_substring = lheComments->substr(lheComments->find("ctauS-") + 6);
      int ctau = stoi(ctau_substring.substr(0, ctau_substring.find('_')));
      MuonSystem->mH = mh;
      MuonSystem->mX = mx;
      MuonSystem->ctau = ctau;

      pair<int, int> signalPair = make_pair(mx, ctau);

      if (Files2D.count(signalPair) == 0) { //create file and tree
        //format file name
        string thisFileName = outfilename;
        thisFileName.erase(thisFileName.end() - 5, thisFileName.end());
        thisFileName += "_" + to_string(mx) + "_" + to_string(ctau) + ".root";

        Files2D[signalPair] = std::make_unique<TFile>(thisFileName.c_str(), "recreate");
        Files2D[signalPair]->cd();
        Trees2D[signalPair] = MuonSystem->tree_->CloneTree(0);
        NEvents2D[signalPair] = new TH1F(Form("NEvents%d%d", mx, ctau), "NEvents", 1, 0.5, 1.5);
        Total2D[signalPair] = new TH1F(Form("Total%d%d", mx, ctau), "Total", 1, 0.5, 1.5);
        accep2D[signalPair] = new TH1F(Form("accep2D%d%d", mx, ctau), "acceptance", 1, 0.5, 1.5);
        accep_met2D[signalPair] = new TH1F(Form("accep_met2D%d%d", mx, ctau), "acceptance_met", 1, 0.5, 1.5);

        cout << "Created new output file " << thisFileName << endl;
      }
      //Fill NEvents hist
      NEvents2D[signalPair]->Fill(1.0, genWeight);
    }
    /* #endregion */

    // [5] Event metadata and per-event nominal weight
    /* #region */
    //event info
    if (isData) { //? can this be better organized?
      NEvents->Fill(1);
      MuonSystem->weight = 1;
    } else {
      MuonSystem->weight = genWeight;
      NEvents->Fill(1, genWeight);
    }
    MuonSystem->runNum = runNum; //? is there any particular reason these need to be here?
    MuonSystem->lumiSec = lumiNum;
    MuonSystem->evtNum = eventNum;

    if (isData && runNum < 360019) continue; //? should this be moved somewhere higher to save compute?
    /* #endregion */

    // [6] MC-only truth matching inputs and MC weights
    /* #region */
    if (!isData) {
      //for DS model
      MuonSystem->nGenParticles = 0;
      for (int i = 0; i < nGenPart; i++) {
        if (abs(GenPart_pdgId[i]) >= 999999) {
          // cout<<eventNum<<","<<gParticleId[i]<<","<<MuonSystem->nGLLP<<endl;
          MuonSystem->gParticleEta[MuonSystem->nGenParticles] = GenPart_eta[i];
          MuonSystem->gParticlePhi[MuonSystem->nGenParticles] = GenPart_phi[i];
          MuonSystem->gParticlePt[MuonSystem->nGenParticles] = GenPart_pt[i];
          MuonSystem->gParticleId[MuonSystem->nGenParticles] = GenPart_pdgId[i];
          float decay_x;
          float decay_y;
          float decay_z;
          bool foundDaughter = false;
          for (int j = 0; j < nGenPart; j++) {
            if (GenPart_genPartIdxMother[j] == i) {
              decay_x = GenPart_vx[j];
              decay_y = GenPart_vy[j];
              decay_z = GenPart_vz[j];
              foundDaughter = true;
              break;
            }
          }
          if (!foundDaughter)
            cout << "didn't find HV LLP daughter " << eventNum << endl;
          MuonSystem->gParticle_decay_vertex_r[MuonSystem->nGenParticles] = sqrt(decay_x * decay_x + decay_y * decay_y);
          MuonSystem->gParticle_decay_vertex_x[MuonSystem->nGenParticles] = decay_x;
          MuonSystem->gParticle_decay_vertex_y[MuonSystem->nGenParticles] = decay_y;
          MuonSystem->gParticle_decay_vertex_z[MuonSystem->nGenParticles] = decay_z;

          TLorentzVector particle = makeTLorentzVectorPtEtaPhiM(GenPart_pt[i], GenPart_eta[i], GenPart_phi[i], GenPart_mass[i]);

          float gParticle_decay_vertex = sqrt(pow(MuonSystem->gParticle_decay_vertex_r[MuonSystem->nGenParticles], 2) + pow(MuonSystem->gParticle_decay_vertex_z[MuonSystem->nGenParticles], 2));

          MuonSystem->gParticle_ctau[MuonSystem->nGenParticles] = gParticle_decay_vertex / (particle.Beta() * particle.Gamma());
          MuonSystem->gParticle_beta[MuonSystem->nGenParticles] = particle.Beta();
          // if (abs(MuonSystem->gLLP_eta[MuonSystem->nGLLP]) < 2.4
          //   && abs(MuonSystem->gLLP_decay_vertex_z[MuonSystem->nGLLP])<1100 && abs(MuonSystem->gLLP_decay_vertex_z[MuonSystem->nGLLP])>400
          //   && MuonSystem->gLLP_decay_vertex_r[MuonSystem->nGLLP] < 695.5) MuonSystem->gLLP_csc[MuonSystem->nGLLP] = true;
          // if (abs(MuonSystem->gLLP_decay_vertex_z[MuonSystem->nGLLP])< 661.0
          //   && MuonSystem->gLLP_decay_vertex_r[MuonSystem->nGLLP] < 800
          //    && MuonSystem->gLLP_decay_vertex_r[MuonSystem->nGLLP] > 200.0) MuonSystem->gLLP_dt[MuonSystem->nGLLP] = true;
          MuonSystem->nGenParticles++;
        }
      }

      //for Twin higgs model
      MuonSystem->nGLLP = 0;
      for (int i = 0; i < nGenPart; i++) {
        if (abs(GenPart_pdgId[i]) == 9000006) {
          // cout<<eventNum<<","<<gParticleId[i]<<","<<MuonSystem->nGLLP<<endl;
          MuonSystem->gLLP_eta[MuonSystem->nGLLP] = GenPart_eta[i];
          MuonSystem->gLLP_phi[MuonSystem->nGLLP] = GenPart_phi[i];
          MuonSystem->gLLP_pt[MuonSystem->nGLLP] = GenPart_pt[i];
          float decay_x;
          float decay_y;
          float decay_z;
          bool foundDaughter = false;
          for (int j = 0; j < nGenPart; j++) {
            if (GenPart_genPartIdxMother[j] == i) {
              decay_x = GenPart_vx[j];
              decay_y = GenPart_vy[j];
              decay_z = GenPart_vz[j];
              foundDaughter = true;
              break;
            }
          }
          if (!foundDaughter)
            cout << "didn't find LLP daughter " << eventNum << endl;
          MuonSystem->gLLP_decay_vertex_r[MuonSystem->nGLLP] = sqrt(decay_x * decay_x + decay_y * decay_y);
          MuonSystem->gLLP_decay_vertex_x[MuonSystem->nGLLP] = decay_x;
          MuonSystem->gLLP_decay_vertex_y[MuonSystem->nGLLP] = decay_y;
          MuonSystem->gLLP_decay_vertex_z[MuonSystem->nGLLP] = decay_z;

          TLorentzVector LLP = makeTLorentzVectorPtEtaPhiM(GenPart_pt[i], GenPart_eta[i], GenPart_phi[i], GenPart_mass[i]);
          MuonSystem->gLLP_e[MuonSystem->nGLLP] = LLP.E();

          float gLLP_decay_vertex = sqrt(pow(MuonSystem->gLLP_decay_vertex_r[MuonSystem->nGLLP], 2) + pow(MuonSystem->gLLP_decay_vertex_z[MuonSystem->nGLLP], 2));

          MuonSystem->gLLP_ctau[MuonSystem->nGLLP] = gLLP_decay_vertex / (LLP.Beta() * LLP.Gamma());
          MuonSystem->gLLP_beta[MuonSystem->nGLLP] = LLP.Beta();
          if (abs(MuonSystem->gLLP_eta[MuonSystem->nGLLP]) < 2.4 && abs(MuonSystem->gLLP_decay_vertex_z[MuonSystem->nGLLP]) < 1100 && abs(MuonSystem->gLLP_decay_vertex_z[MuonSystem->nGLLP]) > 400 && MuonSystem->gLLP_decay_vertex_r[MuonSystem->nGLLP] < 695.5)
            MuonSystem->gLLP_csc[MuonSystem->nGLLP] = true;
          if (abs(MuonSystem->gLLP_decay_vertex_z[MuonSystem->nGLLP]) < 661.0 && MuonSystem->gLLP_decay_vertex_r[MuonSystem->nGLLP] < 800 && MuonSystem->gLLP_decay_vertex_r[MuonSystem->nGLLP] > 200.0)
            MuonSystem->gLLP_dt[MuonSystem->nGLLP] = true;
          MuonSystem->nGLLP++;
        }
      }

      MuonSystem->npu = Pileup_nTrueInt;
      MuonSystem->pileupWeight = helper->getPileupWeight(Pileup_nTrueInt);
      MuonSystem->pileupWeightUp = helper->getPileupWeightUp(Pileup_nTrueInt) / MuonSystem->pileupWeight;
      MuonSystem->pileupWeightDown = helper->getPileupWeightDown(Pileup_nTrueInt) / MuonSystem->pileupWeight;

      for (unsigned int i = 0; i < 9; i++) {
        MuonSystem->LHEScaleWeight[i] = LHEScaleWeight[i];
      }

    } //end of isData
    /* #endregion */

    // [7] Event-level observables and acceptance bookkeeping
    /* #region */
    //get NPU
    MuonSystem->npv = nPV;
    MuonSystem->rho = fixedGridRhoFastjetAll;
    MuonSystem->met = metType1Pt;
    MuonSystem->metPhi = metType1Phi;

    MuonSystem->Puppimet = PuppiMET_pt;
    MuonSystem->PuppimetPhi = PuppiMET_phi;

    if (signalScan && !isData)
      Total2D[make_pair(MuonSystem->mX, MuonSystem->ctau)]->Fill(1.0, genWeight * MuonSystem->pileupWeight);
    if (signalScan && !isData) {
      accep2D[make_pair(MuonSystem->mX, MuonSystem->ctau)]->Fill(1.0, genWeight * MuonSystem->pileupWeight);
    } else if (!isData) {
      if (MuonSystem->gLLP_csc[0] && MuonSystem->gLLP_csc[1])
        accep_csccsc->Fill(1.0, genWeight * MuonSystem->pileupWeight);
      if ((MuonSystem->gLLP_dt[0] && MuonSystem->gLLP_csc[1]) || (MuonSystem->gLLP_dt[1] && MuonSystem->gLLP_csc[0]))
        accep_cscdt->Fill(1.0, genWeight * MuonSystem->pileupWeight);
    }
    /* #endregion */

    // [8] Event filters, trigger bits, and jet veto maps
    /* #region */
    // noise filters
    MuonSystem->Flag_goodVertices = Flag_goodVertices;
    MuonSystem->Flag_globalSuperTightHalo2016Filter = Flag_globalSuperTightHalo2016Filter;
    MuonSystem->Flag_EcalDeadCellTriggerPrimitiveFilter = Flag_EcalDeadCellTriggerPrimitiveFilter;
    MuonSystem->Flag_BadPFMuonFilter = Flag_BadPFMuonFilter;
    MuonSystem->Flag_BadPFMuonDzFilter = Flag_BadPFMuonDzFilter;
    MuonSystem->Flag_hfNoisyHitsFilter = Flag_hfNoisyHitsFilter;
    MuonSystem->Flag_eeBadScFilter = Flag_eeBadScFilter;
    MuonSystem->Flag_all = Flag_eeBadScFilter && Flag_hfNoisyHitsFilter && Flag_BadPFMuonDzFilter && Flag_BadPFMuonFilter && Flag_EcalDeadCellTriggerPrimitiveFilter && Flag_globalSuperTightHalo2016Filter && Flag_goodVertices;
    if (analysisTag == "Summer24")
      MuonSystem->Flag_ecalBadCalibFilter = Flag_ecalBadCalibFilter;

    // Flag_ecalBadCalibFilter for nanoAOD: https://twiki.cern.ch/twiki/bin/view/CMS/MissingETOptionalFiltersRun2#ECal_BadCalibration_Filter_Flag
    if (analysisTag == "Summer24")
      MuonSystem->Flag_ecalBadCalibFilter = Flag_ecalBadCalibFilter;
    else {
      MuonSystem->Flag_ecalBadCalibFilter = true;
      if (isData && runNum >= 362433 && runNum <= 367144) {
        if (PuppiMET_pt > 100) {
          for (int i = 0; i < nJets; i++) {
            if (jetPt[i] < 50)
              continue;
            if (!(jetEta[i] <= -0.1 && jetEta[i] >= -0.5 && jetPhi[i] < -1.8 && jetPhi[i] > -2.1))
              continue;
            if (!(Jet_neEmEF[i] > 0.9 || Jet_chEmEF[i] > 0.9))
              continue;
            if (deltaPhi(PuppiMET_phi, jetPhi[i]) < 2.9)
              continue;
            Flag_ecalBadCalibFilter = false;
          }
        }
      }
    }

    // jet veto map, following selections here: https://cms-jerc.web.cern.ch/Recommendations/#jet-veto-maps
    MuonSystem->jetVeto = true;
    for (int i = 0; i < nJets; i++) {
      if (jetPt[i] <= 15)
        continue;
      if (Jet_neEmEF[i] + Jet_chEmEF[i] >= 0.9)
        continue;
      if (!jetPassIDTight[i])
        continue;
      //remove overlaps
      bool overlap = false;
      for (int j = 0; j < nMuons; j++) {
        if (!Muon_isPFcand[j])
          continue;
        if (RazorAnalyzerMerged::deltaR(jetEta[i], jetPhi[i], muonEta[j], muonPhi[j]) < 0.2)
          overlap = true;
      }
      if (overlap)
        continue;
      helper->getJetVetoMap(0, 1);
      if (helper->getJetVetoMap(jetEta[i], jetPhi[i]) > 0.0)
        MuonSystem->jetVeto = false;
      if (analysisTag == "Summer24" && helper->getJetVetoFpixMap(jetEta[i], jetPhi[i]) > 0.0)
        MuonSystem->jetVeto = false;
    }

    MuonSystem->HLT_CSCCSC = HLT_CscCluster_Loose;
    MuonSystem->HLT_CSCDT = HLT_L1CSCShower_DTCluster50;
    /* #endregion */

    // [9] Lepton object selection and lepton branch fill
    /* #region */
    //*************************************************************************
    //Start Object Selection
    //*************************************************************************
    std::vector<leptons> Leptons;
    //-------------------------------
    //Muons
    //-------------------------------
    for (int i = 0; i < nMuons; i++) {
      if (!muonIsLoose[i])
        continue;
      if (muonPt[i] < 25)
        continue;
      if (fabs(muonEta[i]) > 2.4)
        continue;

      //remove overlaps
      bool overlap = false;
      for (auto& lep : Leptons) {
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

      tmpMuon.passLooseIso = muon_pfRelIso04_all[i] < 0.25;
      tmpMuon.passTightIso = muon_pfRelIso04_all[i] < 0.15;
      tmpMuon.passVTightIso = muon_pfRelIso04_all[i] < 0.10;
      tmpMuon.passVVTightIso = muon_pfRelIso04_all[i] < 0.05;

      tmpMuon.passVetoId = false;
      Leptons.push_back(tmpMuon);
    }

    //-------------------------------
    //Electrons
    //-------------------------------
    for (int i = 0; i < nElectrons; i++) {
      if (!ele_passCutBasedIDVeto[i])
        continue;
      if (elePt[i] < 35)
        continue;
      if (fabs(eleEta[i]) > 2.5)
        continue;

      //remove overlaps
      bool overlap = false;
      for (auto& lep : Leptons) {
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

    for (auto& tmp : Leptons) {
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
    /* #endregion */

    // [10] Jet selection, JES propagation, and jet branch fill
    /* #region */
    //-----------------------------------------------
    //Select Jets
    //-----------------------------------------------

    std::vector<jets> Jets;
    float MetX_JESDown = 0;
    float MetY_JESDown = 0;

    float MetX_JESUp = 0;
    float MetY_JESUp = 0;

    for (int i = 0; i < nJets; i++) {
      if (fabs(jetEta[i]) >= 3.0)
        continue;
      if (jetPt[i] < 20)
        continue;
      if (!jetPassIDTight[i] && !jetPassIDTightLepVeto[i])
        continue;
      //------------------------------------------------------------
      //exclude selected muons and electrons from the jet collection
      //------------------------------------------------------------
      double deltaR = -1;
      for (auto& lep : Leptons) {
        double thisDR = RazorAnalyzerMerged::deltaR(jetEta[i], jetPhi[i], lep.lepton.Eta(), lep.lepton.Phi());
        if (deltaR < 0 || thisDR < deltaR)
          deltaR = thisDR;
      }
      if (deltaR > 0 && deltaR < 0.4)
        continue; //jet matches a selected lepton

      TLorentzVector thisJet = makeTLorentzVector(jetPt[i], jetEta[i], jetPhi[i], jetE[i]);

      jets tmpJet;
      tmpJet.jet = thisJet;
      tmpJet.passId = jetPassIDTightLepVeto[i];

      // calculate jet energy scale uncertainty
      double unc = helper->getJecUnc(jetPt[i], jetEta[i], runNum); //use run=999 as default
      tmpJet.jetPtJESUp = jetPt[i] * (1 + unc);
      tmpJet.jetPtJESDown = jetPt[i] * (1 - unc);
      tmpJet.jetEJESUp = jetE[i] * (1 + unc);
      tmpJet.jetEJESDown = jetE[i] * (1 - unc);
      // cout<<jetPt[i]<<","<<jetEta[i]<<","<<unc<<endl;
      tmpJet.JecUnc = unc;
      TLorentzVector thisJetJESUp = makeTLorentzVector(tmpJet.jetPtJESUp, jetEta[i], jetPhi[i], tmpJet.jetEJESUp);
      TLorentzVector thisJetJESDown = makeTLorentzVector(tmpJet.jetPtJESDown, jetEta[i], jetPhi[i], tmpJet.jetEJESDown);

      MetX_JESUp += -1 * (thisJetJESUp.Px() - thisJet.Px());
      MetY_JESUp += -1 * (thisJetJESUp.Py() - thisJet.Py());

      MetX_JESDown += -1 * (thisJetJESDown.Px() - thisJet.Px());
      MetY_JESDown += -1 * (thisJetJESDown.Py() - thisJet.Py());
      // done propogating JES uncertainty

      Jets.push_back(tmpJet);
    }

    sort(Jets.begin(), Jets.end(), my_largest_pt_jet);

    for (auto& tmp : Jets) {
      if (tmp.jet.Pt() < 30)
        continue;

      MuonSystem->jetE[MuonSystem->nJets] = tmp.jet.E();
      MuonSystem->jetPt[MuonSystem->nJets] = tmp.jet.Pt();
      MuonSystem->jetPtJESUp[MuonSystem->nJets] = tmp.jetPtJESUp;
      MuonSystem->jetPtJESDown[MuonSystem->nJets] = tmp.jetPtJESDown;
      MuonSystem->jetEta[MuonSystem->nJets] = tmp.jet.Eta();
      MuonSystem->jetPhi[MuonSystem->nJets] = tmp.jet.Phi();
      MuonSystem->jetTightPassId[MuonSystem->nJets] = tmp.passId;

      MuonSystem->nJets++;
    }

    TLorentzVector met = makeTLorentzVectorPtEtaPhiM(MuonSystem->Puppimet, 0, MuonSystem->PuppimetPhi, 0);

    //JES up
    float MetXJESUp = met.Px() + MetX_JESUp;
    float MetYJESUp = met.Py() + MetY_JESUp;
    MuonSystem->PuppimetJESUp = sqrt(pow(MetXJESUp, 2) + pow(MetYJESUp, 2));
    MuonSystem->PuppimetPhiJESUp = atan(MetYJESUp / MetXJESUp);
    if (MetXJESUp < 0.0)
      MuonSystem->PuppimetPhiJESUp = RazorAnalyzerMerged::deltaPhi(TMath::Pi() + MuonSystem->PuppimetPhiJESUp, 0.0);

    //JES down
    float MetXJESDown = met.Px() + MetX_JESDown;
    float MetYJESDown = met.Py() + MetY_JESDown;
    MuonSystem->PuppimetJESDown = sqrt(pow(MetXJESDown, 2) + pow(MetYJESDown, 2));
    MuonSystem->PuppimetPhiJESDown = atan(MetYJESDown / MetXJESDown);
    if (MetXJESDown < 0.0)
      MuonSystem->PuppimetPhiJESDown = RazorAnalyzerMerged::deltaPhi(TMath::Pi() + MuonSystem->PuppimetPhiJESDown, 0.0);
    /* #endregion */

    // [11] Rechit-level pass-through fill (CSC/DT/RPC + optional aliases)
    /* #region */
    MuonSystem->nDTRechits = nDtRechits;
    MuonSystem->nRpcRechits = nRpc;
    MuonSystem->nCscRechits = ncscRechits;
    if (MuonSystem->nCscRechits > N_MAX_CSCRECHITS)
      MuonSystem->nCscRechits = N_MAX_CSCRECHITS;
    if (MuonSystem->nDTRechits > N_MAX_DTRECHITS)
      MuonSystem->nDTRechits = N_MAX_DTRECHITS;
    if (MuonSystem->nRpcRechits > N_MAX_RPCRECHITS)
      MuonSystem->nRpcRechits = N_MAX_RPCRECHITS;

    for (int i = 0; i < MuonSystem->nCscRechits; i++) {
      MuonSystem->CscRechitsClusterId[i] = -1;
      MuonSystem->CscRechitsQuality[i] = cscRechitsQuality[i];
      MuonSystem->CscRechitsChamber[i] = cscRechitsChamber[i];
      MuonSystem->CscRechitsStation[i] = cscRechitsStation[i];
      MuonSystem->CscRechitsEta[i] = cscRechitsEta[i];
      MuonSystem->CscRechitsPhi[i] = cscRechitsPhi[i];
      MuonSystem->CscRechitsX[i] = cscRechitsX[i];
      MuonSystem->CscRechitsY[i] = cscRechitsY[i];
      MuonSystem->CscRechitsZ[i] = cscRechitsZ[i];
      MuonSystem->CscRechitsTpeak[i] = cscRechitsTpeak[i];
      MuonSystem->CscRechitsTwire[i] = cscRechitsTwire[i];

      if (hasCscRechitsIChamber)
        MuonSystem->CscRechitsIChamber[i] = cscRechitsIChamber_opt[i];
      if (hasCscRechitsNStrips)
        MuonSystem->CscRechitsNStrips[i] = cscRechitsNStrips_opt[i];
      if (hasCscRechitsWGroupsBX)
        MuonSystem->CscRechitsWGroupsBX[i] = cscRechitsWGroupsBX_opt[i];
      if (hasCscRechitsHitWire)
        MuonSystem->CscRechitsHitWire[i] = cscRechitsHitWire_opt[i];
      if (hasCscRechitsNWireGroups)
        MuonSystem->CscRechitsNWireGroups[i] = cscRechitsNWireGroups_opt[i];
      if (hasCscRechitsE)
        MuonSystem->CscRechitsE[i] = cscRechitsE_opt[i];
    }
    for (int i = 0; i < MuonSystem->nDTRechits; i++) {
      MuonSystem->DtRechitsClusterId[i] = -1;
      MuonSystem->DtRechitsLayer[i] = dtRechitLayer[i];
      MuonSystem->DtRechitsSuperLayer[i] = dtRechitSuperLayer[i];
      MuonSystem->DtRechitsStation[i] = dtRechitStation[i];
      MuonSystem->DtRechitsWheel[i] = dtRechitWheel[i];
      MuonSystem->DtRechitsEta[i] = dtRechitCorrectEta[i];
      MuonSystem->DtRechitsPhi[i] = dtRechitCorrectPhi[i];
      MuonSystem->DtRechitsX[i] = dtRechitCorrectX[i];
      MuonSystem->DtRechitsY[i] = dtRechitCorrectY[i];
      MuonSystem->DtRechitsZ[i] = dtRechitCorrectZ[i];
      if (hasDtRechitsSector)
        MuonSystem->DtRechitsSector[i] = dtRechitsSector_opt[i];
    }
    for (int i = 0; i < MuonSystem->nRpcRechits; i++) {
      MuonSystem->RpcRecHitsClusterId[i] = -1;
      MuonSystem->RpcRecHitsBx[i] = rpcBx[i];
      MuonSystem->RpcRecHitsRegion[i] = rpcRegion[i];
      MuonSystem->RpcRecHitsRing[i] = rpcRing[i];
      MuonSystem->RpcRecHitsLayer[i] = rpcLayer[i];
      MuonSystem->RpcRecHitsStation[i] = rpcStation[i];
      MuonSystem->RpcRecHitsSector[i] = rpcSector[i];
      MuonSystem->RpcRecHitsX[i] = rpcX[i];
      MuonSystem->RpcRecHitsY[i] = rpcY[i];
      MuonSystem->RpcRecHitsZ[i] = rpcZ[i];
      MuonSystem->RpcRecHitsPhi[i] = rpcPhi[i];
      MuonSystem->RpcRecHitsEta[i] = rpcEta[i];
      MuonSystem->RpcRecHitsTime[i] = rpcTime[i];
      MuonSystem->RpcRecHitsTimeError[i] = rpcTimeError[i];
    }
    /* #endregion */

    // [12] Ring occupancy summaries from raw rechits
    /* #region */
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

    for (int i = 0; i < nDtRechits; i++) {
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
    /* #endregion */

    // [13] CSC rechit clustering, feature fill, and matching/veto variables
    /* #region */
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

    for (int i = 0; i < ncscRechits; i++) {

      int layer = 0;
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
    //Do DBSCAN Clustering

    int min_point = 50; //minimum number of Rechitss to call it a cluster
    float epsilon = 0.4; //cluster radius parameter
    CACluster ds(min_point, epsilon, points);
    ds.run();
    ds.clusterProperties();
    ds.sort_clusters();
    const auto& cscPointsClustered = ds.points();
    int nCscPointsOut = MuonSystem->nCscRechits;
    if (nCscPointsOut > (int)cscPointsClustered.size())
      nCscPointsOut = cscPointsClustered.size();
    for (int i = 0; i < nCscPointsOut; i++) {
      MuonSystem->CscRechitsClusterId[i] = cscPointsClustered[i].clusterID;
    }

    MuonSystem->nCscRechitClusters = 0;
    MuonSystem->nCscRechitClusters_nocut = 0;
    for (auto& tmp : ds.clusters) {
      MuonSystem->nCscRechitClusters_nocut++;
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
      MuonSystem->cscRechitClusterNRechitME1112[MuonSystem->nCscRechitClusters] = tmp.nCscRechitsChamberMinus11 + tmp.nCscRechitsChamberMinus12 + tmp.nCscRechitsChamberPlus11 + tmp.nCscRechitsChamberPlus12;
      MuonSystem->cscRechitClusterNRechitStation12[MuonSystem->nCscRechitClusters] = MuonSystem->cscRechitClusterNRechitME1112[MuonSystem->nCscRechitClusters] + tmp.nCscRechitsChamberMinus13 + tmp.nCscRechitsChamberPlus13 + tmp.nCscRechitsChamberMinus21 + tmp.nCscRechitsChamberPlus21 + tmp.nCscRechitsChamberMinus22 + tmp.nCscRechitsChamberPlus22;

      float efficiencies[] = {
          helper->getHMTTriggerEff(13, tmp.nCscRechitsChamberPlus13),
          helper->getHMTTriggerEff(21, tmp.nCscRechitsChamberPlus21),
          helper->getHMTTriggerEff(22, tmp.nCscRechitsChamberPlus22),
          helper->getHMTTriggerEff(31, tmp.nCscRechitsChamberPlus31),
          helper->getHMTTriggerEff(32, tmp.nCscRechitsChamberPlus32),
          helper->getHMTTriggerEff(41, tmp.nCscRechitsChamberPlus41),
          helper->getHMTTriggerEff(42, tmp.nCscRechitsChamberPlus42),
          helper->getHMTTriggerEff(13, tmp.nCscRechitsChamberMinus13),
          helper->getHMTTriggerEff(21, tmp.nCscRechitsChamberMinus21),
          helper->getHMTTriggerEff(22, tmp.nCscRechitsChamberMinus22),
          helper->getHMTTriggerEff(31, tmp.nCscRechitsChamberMinus31),
          helper->getHMTTriggerEff(32, tmp.nCscRechitsChamberMinus32),
          helper->getHMTTriggerEff(41, tmp.nCscRechitsChamberMinus41),
          helper->getHMTTriggerEff(42, tmp.nCscRechitsChamberPlus42), //! Should this be Minus42 ???
      };

      MuonSystem->cscRechitClusterHMTEfficiency[MuonSystem->nCscRechitClusters] = 1.0;
      for (const float& eff : efficiencies) {
        MuonSystem->cscRechitClusterHMTEfficiency[MuonSystem->nCscRechitClusters] *= (1 - eff);
      }
      MuonSystem->cscRechitClusterHMTEfficiency[MuonSystem->nCscRechitClusters] = 1 - MuonSystem->cscRechitClusterHMTEfficiency[MuonSystem->nCscRechitClusters];

      MuonSystem->cscRechitClusterMaxChamber[MuonSystem->nCscRechitClusters] = tmp.maxChamber;
      MuonSystem->cscRechitClusterMaxChamberRatio[MuonSystem->nCscRechitClusters] = 1.0 * tmp.maxChamberRechits / tmp.nhits;
      MuonSystem->cscRechitClusterNChamber[MuonSystem->nCscRechitClusters] = tmp.nChamber;
      MuonSystem->cscRechitClusterMaxStation[MuonSystem->nCscRechitClusters] = tmp.maxStation;
      MuonSystem->cscRechitClusterMaxStationRatio[MuonSystem->nCscRechitClusters] = 1.0 * tmp.maxStationRechits / tmp.nhits;

      MuonSystem->cscRechitClusterNStation10[MuonSystem->nCscRechitClusters] = tmp.nStation10;
      MuonSystem->cscRechitClusterAvgStation10[MuonSystem->nCscRechitClusters] = tmp.avgStation10;

      //Jet veto/ muon veto
      MuonSystem->cscRechitClusterJetVetoPt[MuonSystem->nCscRechitClusters] = 0.0;
      MuonSystem->cscRechitClusterJetVetoE[MuonSystem->nCscRechitClusters] = 0.0;
      MuonSystem->cscRechitClusterMuonVetoPt[MuonSystem->nCscRechitClusters] = 0.0;
      MuonSystem->cscRechitClusterMuonVetoE[MuonSystem->nCscRechitClusters] = 0.0;

      // jet veto
      for (int i = 0; i < nJets; i++) {
        if (fabs(jetEta[i]) > 3.0)
          continue;
        if (RazorAnalyzerMerged::deltaR(jetEta[i], jetPhi[i], MuonSystem->cscRechitClusterEta[MuonSystem->nCscRechitClusters], MuonSystem->cscRechitClusterPhi[MuonSystem->nCscRechitClusters]) < 0.4 && jetPt[i] > MuonSystem->cscRechitClusterJetVetoPt[MuonSystem->nCscRechitClusters]) {
          MuonSystem->cscRechitClusterJetVetoPt[MuonSystem->nCscRechitClusters] = jetPt[i];
          MuonSystem->cscRechitClusterJetVetoE[MuonSystem->nCscRechitClusters] = jetE[i];
          MuonSystem->cscRechitClusterJetVetoTightId[MuonSystem->nCscRechitClusters] = jetPassIDTight[i];

          double unc = helper->getJecUnc(jetPt[i], jetEta[i], runNum); //use run=999 as default
          MuonSystem->cscRechitClusterJetVetoPtJESUp[MuonSystem->nCscRechitClusters] = jetPt[i] * (1 + unc);
          MuonSystem->cscRechitClusterJetVetoPtJESDown[MuonSystem->nCscRechitClusters] = jetPt[i] * (1 - unc);
        }
      }
      float min_deltaR = 15.;
      int index = 999;

      for (int i = 0; i < nMuons; i++) {
        if (fabs(muonEta[i]) > 3.0)
          continue;
        if (RazorAnalyzerMerged::deltaR(muonEta[i], muonPhi[i], MuonSystem->cscRechitClusterEta[MuonSystem->nCscRechitClusters], MuonSystem->cscRechitClusterPhi[MuonSystem->nCscRechitClusters]) < 0.4 && muonPt[i] > MuonSystem->cscRechitClusterMuonVetoPt[MuonSystem->nCscRechitClusters]) {
          MuonSystem->cscRechitClusterMuonVetoPt[MuonSystem->nCscRechitClusters] = muonPt[i];
          MuonSystem->cscRechitClusterMuonVetoE[MuonSystem->nCscRechitClusters] = muonE[i];
          MuonSystem->cscRechitClusterMuonVetoGlobal[MuonSystem->nCscRechitClusters] = muon_isGlobal[i];
          MuonSystem->cscRechitClusterMuonVetoLooseId[MuonSystem->nCscRechitClusters] = muonIsLoose[i];
        }
        if (RazorAnalyzerMerged::deltaR(muonEta[i], muonPhi[i], MuonSystem->cscRechitClusterEta[MuonSystem->nCscRechitClusters], MuonSystem->cscRechitClusterPhi[MuonSystem->nCscRechitClusters]) < 0.8 && muonPt[i] > MuonSystem->cscRechitClusterMuonVetoPt0p8Thresh[MuonSystem->nCscRechitClusters]) {
          MuonSystem->cscRechitClusterMuonVetoPt0p8Thresh[MuonSystem->nCscRechitClusters] = muonPt[i];
          MuonSystem->cscRechitClusterMuonVetoE0p8Thresh[MuonSystem->nCscRechitClusters] = muonE[i];
          MuonSystem->cscRechitClusterMuonVetoGlobal0p8Thresh[MuonSystem->nCscRechitClusters] = muon_isGlobal[i];
          MuonSystem->cscRechitClusterMuonVetoLooseId0p8Thresh[MuonSystem->nCscRechitClusters] = muonIsLoose[i];
        }
      }
      if (!isData) {
        // match to gen level LLP
        min_deltaR = 15.;
        index = 999;
        for (int j = 0; j < MuonSystem->nGLLP; j++) {
          double current_delta_r = RazorAnalyzerMerged::deltaR(MuonSystem->cscRechitClusterEta[MuonSystem->nCscRechitClusters], MuonSystem->cscRechitClusterPhi[MuonSystem->nCscRechitClusters], MuonSystem->gLLP_eta[j], MuonSystem->gLLP_phi[j]);
          if (current_delta_r < min_deltaR) {
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

      //match to MB1 DT segments
      for (int i = 0; i < nDtSeg; i++) {
        if (RazorAnalyzerMerged::deltaR(dtSegEta[i], dtSegPhi[i], MuonSystem->cscRechitClusterEta[MuonSystem->nCscRechitClusters], MuonSystem->cscRechitClusterPhi[MuonSystem->nCscRechitClusters]) < 0.4) {
          MuonSystem->cscRechitCluster_match_dtSeg_0p4[MuonSystem->nCscRechitClusters]++;
          if (dtSegStation[i] == 1)
            MuonSystem->cscRechitCluster_match_MB1Seg_0p4[MuonSystem->nCscRechitClusters]++;
        }
      }

      //match to RPC hits in RE1/2
      for (int i = 0; i < nRpc; i++) {
        float rpcR = sqrt(rpcX[i] * rpcX[i] + rpcY[i] * rpcY[i]);
        if (RazorAnalyzerMerged::deltaR(rpcEta[i], rpcPhi[i], MuonSystem->cscRechitClusterEta[MuonSystem->nCscRechitClusters], MuonSystem->cscRechitClusterPhi[MuonSystem->nCscRechitClusters]) < 0.4) {
          if (rpcR < 461.0 && rpcR > 275 && abs(rpcZ[i]) > 663 && abs(rpcZ[i]) < 730)
            MuonSystem->cscRechitCluster_match_RE12_0p4[MuonSystem->nCscRechitClusters]++;
          if (rpcR < 470 && rpcR > 380 && abs(rpcZ[i]) < 661)
            MuonSystem->cscRechitCluster_match_RB1_0p4[MuonSystem->nCscRechitClusters]++;
        }
      }

      MuonSystem->cscRechitClusterMet_dPhi[MuonSystem->nCscRechitClusters] = RazorAnalyzerMerged::deltaPhi(MuonSystem->cscRechitClusterPhi[MuonSystem->nCscRechitClusters], MuonSystem->metPhi);
      MuonSystem->cscRechitClusterPuppiMet_dPhi[MuonSystem->nCscRechitClusters] = RazorAnalyzerMerged::deltaPhi(MuonSystem->cscRechitClusterPhi[MuonSystem->nCscRechitClusters], MuonSystem->PuppimetPhi);
      if (MuonSystem->nTaus == 1) {
        MuonSystem->cscRechitClusterPromptTauDeltaEta[MuonSystem->nCscRechitClusters] = MuonSystem->tauEta[0] - MuonSystem->cscRechitClusterEta[MuonSystem->nCscRechitClusters];
        MuonSystem->cscRechitClusterPromptTauDeltaPhi[MuonSystem->nCscRechitClusters] = RazorAnalyzerMerged::deltaPhi(MuonSystem->tauPhi[0], MuonSystem->cscRechitClusterPhi[MuonSystem->nCscRechitClusters]);
        MuonSystem->cscRechitClusterPromptTauDeltaR[MuonSystem->nCscRechitClusters] = RazorAnalyzerMerged::deltaR(MuonSystem->tauEta[0], MuonSystem->tauPhi[0], MuonSystem->cscRechitClusterEta[MuonSystem->nCscRechitClusters], MuonSystem->cscRechitClusterPhi[MuonSystem->nCscRechitClusters]);
      }
      if (MuonSystem->nLeptons > 0) {
        if (abs(MuonSystem->lepPdgId[0]) == 11) {
          MuonSystem->cscRechitClusterPromptEleDeltaEta[MuonSystem->nCscRechitClusters] = MuonSystem->lepEta[0] - MuonSystem->cscRechitClusterEta[MuonSystem->nCscRechitClusters];
          MuonSystem->cscRechitClusterPromptEleDeltaPhi[MuonSystem->nCscRechitClusters] = RazorAnalyzerMerged::deltaPhi(MuonSystem->lepPhi[0], MuonSystem->cscRechitClusterPhi[MuonSystem->nCscRechitClusters]);
          MuonSystem->cscRechitClusterPromptEleDeltaR[MuonSystem->nCscRechitClusters] = RazorAnalyzerMerged::deltaR(MuonSystem->lepEta[0], MuonSystem->lepPhi[0], MuonSystem->cscRechitClusterEta[MuonSystem->nCscRechitClusters], MuonSystem->cscRechitClusterPhi[MuonSystem->nCscRechitClusters]);
        }
        if (abs(MuonSystem->lepPdgId[0]) == 13) {
          MuonSystem->cscRechitClusterPromptMuDeltaEta[MuonSystem->nCscRechitClusters] = MuonSystem->lepEta[0] - MuonSystem->cscRechitClusterEta[MuonSystem->nCscRechitClusters];
          MuonSystem->cscRechitClusterPromptMuDeltaPhi[MuonSystem->nCscRechitClusters] = RazorAnalyzerMerged::deltaPhi(MuonSystem->lepPhi[0], MuonSystem->cscRechitClusterPhi[MuonSystem->nCscRechitClusters]);
          MuonSystem->cscRechitClusterPromptMuDeltaR[MuonSystem->nCscRechitClusters] = RazorAnalyzerMerged::deltaR(MuonSystem->lepEta[0], MuonSystem->lepPhi[0], MuonSystem->cscRechitClusterEta[MuonSystem->nCscRechitClusters], MuonSystem->cscRechitClusterPhi[MuonSystem->nCscRechitClusters]);
        }
      }

      MuonSystem->nCscRechitClusters++;
    }
    /* #endregion */

    // [14] DT/RPC clustering, DT feature fill, and matching/veto variables
    /* #region */
    // DT cluster

    points.clear();

    for (int i = 0; i < nDtRechits; i++) {
      Rechits p;

      p.phi = dtRechitCorrectPhi[i];
      p.eta = dtRechitCorrectEta[i];
      p.x = dtRechitCorrectX[i];
      p.y = dtRechitCorrectY[i];
      p.z = dtRechitCorrectZ[i];
      p.t = -999.;
      p.twire = -999.;
      p.station = dtRechitStation[i];
      p.chamber = dtRechitWheel[i];
      p.superlayer = dtRechitSuperLayer[i];
      p.wheel = dtRechitWheel[i];
      p.clusterID = UNCLASSIFIED;
      points.push_back(p);
    }

    //Do DBSCAN Clustering
    int min_point_dt = 50; //minimum number of segments to call it a cluster
    float epsilon_dt = 0.2; //cluster radius parameter
    CACluster ds_dtRechit(min_point_dt, epsilon_dt, points);
    ds_dtRechit.run();
    ds_dtRechit.clusterProperties();
    ds_dtRechit.sort_clusters();
    const auto& dtPointsClustered = ds_dtRechit.points();
    int nDtPointsOut = MuonSystem->nDTRechits;
    if (nDtPointsOut > (int)dtPointsClustered.size())
      nDtPointsOut = dtPointsClustered.size();
    for (int i = 0; i < nDtPointsOut; i++) {
      MuonSystem->DtRechitsClusterId[i] = dtPointsClustered[i].clusterID;
    }

    // RPC cluster IDs
    points.clear();
    for (int i = 0; i < MuonSystem->nRpcRechits; i++) {
      Rechits p;
      p.phi = rpcPhi[i];
      p.eta = rpcEta[i];
      p.x = rpcX[i];
      p.y = rpcY[i];
      p.z = rpcZ[i];
      p.t = rpcTime[i];
      p.twire = rpcTime[i];
      p.station = rpcStation[i];
      p.chamber = rpcSector[i];
      p.layer = rpcLayer[i];
      p.superlayer = 0;
      p.wheel = rpcRing[i];
      p.clusterID = UNCLASSIFIED;
      points.push_back(p);
    }

    int min_point_rpc = 5;
    float epsilon_rpc = 0.2;
    CACluster ds_rpcRechit(min_point_rpc, epsilon_rpc, points);
    ds_rpcRechit.run();
    ds_rpcRechit.clusterProperties();
    ds_rpcRechit.sort_clusters();
    const auto& rpcPointsClustered = ds_rpcRechit.points();
    int nRpcPointsOut = MuonSystem->nRpcRechits;
    if (nRpcPointsOut > (int)rpcPointsClustered.size())
      nRpcPointsOut = rpcPointsClustered.size();
    for (int i = 0; i < nRpcPointsOut; i++) {
      MuonSystem->RpcRecHitsClusterId[i] = rpcPointsClustered[i].clusterID;
    }

    MuonSystem->nDtRechitClusters = 0;
    MuonSystem->nDtRechitClusters_nocut = 0;
    for (auto& tmp : ds_dtRechit.clusters) {
      //remove overlaps
      bool overlap = false;
      for (int i = 0; i < MuonSystem->nCscRechitClusters; i++) {
        if (RazorAnalyzerMerged::deltaR(MuonSystem->cscRechitClusterEta[i], MuonSystem->cscRechitClusterPhi[i], tmp.eta, tmp.phi) < 0.4)
          overlap = true;
      }
      if (overlap)
        continue;
      MuonSystem->nDtRechitClusters_nocut++;
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
      for (int i = 0; i < 12; ++i) {
        if (distribution(generator) < prob)
          MuonSystem->dtRechitClusterNoiseHitStation1[MuonSystem->nDtRechitClusters]++;
      }
      for (int i = 0; i < 12; ++i) {
        if (distribution(generator) < prob)
          MuonSystem->dtRechitClusterNoiseHitStation2[MuonSystem->nDtRechitClusters]++;
      }
      for (int i = 0; i < 12; ++i) {
        if (distribution(generator) < prob)
          MuonSystem->dtRechitClusterNoiseHitStation3[MuonSystem->nDtRechitClusters]++;
      }
      for (int i = 0; i < 8; ++i) {
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

      //Jet veto/ muon veto
      MuonSystem->dtRechitClusterJetVetoPt[MuonSystem->nDtRechitClusters] = 0.0;
      MuonSystem->dtRechitClusterJetVetoE[MuonSystem->nDtRechitClusters] = 0.0;
      MuonSystem->dtRechitClusterMuonVetoPt[MuonSystem->nDtRechitClusters] = 0.0;
      MuonSystem->dtRechitClusterMuonVetoE[MuonSystem->nDtRechitClusters] = 0.0;

      // jet veto
      for (int i = 0; i < nJets; i++) {
        if (fabs(jetEta[i]) > 3.0)
          continue;
        if (RazorAnalyzerMerged::deltaR(jetEta[i], jetPhi[i], MuonSystem->dtRechitClusterEta[MuonSystem->nDtRechitClusters], MuonSystem->dtRechitClusterPhi[MuonSystem->nDtRechitClusters]) < 0.4 && jetPt[i] > MuonSystem->dtRechitClusterJetVetoPt[MuonSystem->nDtRechitClusters]) {
          MuonSystem->dtRechitClusterJetVetoPt[MuonSystem->nDtRechitClusters] = jetPt[i];
          MuonSystem->dtRechitClusterJetVetoE[MuonSystem->nDtRechitClusters] = jetE[i];
          MuonSystem->dtRechitClusterJetVetoTightId[MuonSystem->nDtRechitClusters] = jetPassIDTight[i];

          double unc = helper->getJecUnc(jetPt[i], jetEta[i], runNum); //use run=999 as default
          MuonSystem->dtRechitClusterJetVetoPtJESUp[MuonSystem->nDtRechitClusters] = jetPt[i] * (1 + unc);
          MuonSystem->dtRechitClusterJetVetoPtJESDown[MuonSystem->nDtRechitClusters] = jetPt[i] * (1 - unc);
        }
      }

      for (int i = 0; i < nMuons; i++) {
        if (fabs(muonEta[i]) > 3.0)
          continue;
        if (RazorAnalyzerMerged::deltaR(muonEta[i], muonPhi[i], MuonSystem->dtRechitClusterEta[MuonSystem->nDtRechitClusters], MuonSystem->dtRechitClusterPhi[MuonSystem->nDtRechitClusters]) < 0.4 && muonPt[i] > MuonSystem->dtRechitClusterMuonVetoPt[MuonSystem->nDtRechitClusters]) {
          MuonSystem->dtRechitClusterMuonVetoPt[MuonSystem->nDtRechitClusters] = muonPt[i];
          MuonSystem->dtRechitClusterMuonVetoE[MuonSystem->nDtRechitClusters] = muonE[i];
          MuonSystem->dtRechitClusterMuonVetoGlobal[MuonSystem->nDtRechitClusters] = muon_isGlobal[i];
          MuonSystem->dtRechitClusterMuonVetoLooseId[MuonSystem->nDtRechitClusters] = muonIsLoose[i];
        }
      }

      for (int i = 0; i < nDtSeg; i++) {
        if (RazorAnalyzerMerged::deltaR(dtSegEta[i], dtSegPhi[i], MuonSystem->dtRechitClusterEta[MuonSystem->nDtRechitClusters], MuonSystem->dtRechitClusterPhi[MuonSystem->nDtRechitClusters]) < 0.4) {
          if (dtSegStation[i] == 1)
            MuonSystem->dtRechitClusterNSegStation1[MuonSystem->nDtRechitClusters] += 1;
          if (dtSegStation[i] == 2)
            MuonSystem->dtRechitClusterNSegStation2[MuonSystem->nDtRechitClusters] += 1;
          if (dtSegStation[i] == 3)
            MuonSystem->dtRechitClusterNSegStation3[MuonSystem->nDtRechitClusters] += 1;
          if (dtSegStation[i] == 4)
            MuonSystem->dtRechitClusterNSegStation4[MuonSystem->nDtRechitClusters] += 1;
        }
        if (abs(RazorAnalyzerMerged::deltaPhi(dtSegPhi[i], MuonSystem->dtRechitClusterPhi[MuonSystem->nDtRechitClusters])) > 2) {
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
      if (!isData) {
        for (int j = 0; j < MuonSystem->nGLLP; j++) {
          double current_delta_r = RazorAnalyzerMerged::deltaR(tmp.eta, tmp.phi, MuonSystem->gLLP_eta[j], MuonSystem->gLLP_phi[j]);
          if (current_delta_r < min_deltaR) {
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

      //match to MB1 DT segments
      for (int i = 0; i < nDtRechits; i++) {
        if (RazorAnalyzerMerged::deltaR(dtRechitCorrectEta[i], dtRechitCorrectPhi[i], MuonSystem->dtRechitClusterEta[MuonSystem->nDtRechitClusters], MuonSystem->dtRechitClusterPhi[MuonSystem->nDtRechitClusters]) < 0.5) {
          if (dtRechitStation[i] == 1)
            MuonSystem->dtRechitCluster_match_MB1hits_0p5[MuonSystem->nDtRechitClusters]++;
        }
        if (RazorAnalyzerMerged::deltaR(dtRechitCorrectEta[i], dtRechitCorrectPhi[i], MuonSystem->dtRechitClusterEta[MuonSystem->nDtRechitClusters], MuonSystem->dtRechitClusterPhi[MuonSystem->nDtRechitClusters]) < 0.4) {
          if (dtRechitStation[i] == 1)
            MuonSystem->dtRechitCluster_match_MB1hits_0p4[MuonSystem->nDtRechitClusters]++;
        }
        if (abs(dtRechitWheel[i] - MuonSystem->dtRechitClusterWheel[MuonSystem->nDtRechitClusters]) == 1 && dtRechitStation[i] == 1) {
          if (abs(RazorAnalyzerMerged::deltaPhi(dtRechitCorrectPhi[i], MuonSystem->dtRechitClusterPhi[MuonSystem->nDtRechitClusters])) < TMath::Pi() / 4.0) {
            if (dtRechitWheel[i] - MuonSystem->dtRechitClusterWheel[MuonSystem->nDtRechitClusters] == 1)
              MuonSystem->dtRechitCluster_match_MB1hits_cosmics_plus[MuonSystem->nDtRechitClusters]++;
            else
              MuonSystem->dtRechitCluster_match_MB1hits_cosmics_minus[MuonSystem->nDtRechitClusters]++;
          }
        }
      }

      std::vector<int> dtRechitCluster_match_rpcBx;

      //match to RPC hits with dPhi<0.5 and same wheel in DT
      for (int i = 0; i < nRpc; i++) {
        float rpcR = sqrt(rpcX[i] * rpcX[i] + rpcY[i] * rpcY[i]);
        if (rpcRegion[i] != 0)
          continue;
        if (abs(RazorAnalyzerMerged::deltaPhi(rpcPhi[i], MuonSystem->dtRechitClusterPhi[MuonSystem->nDtRechitClusters])) < 0.5) {
          if (rpcRing[i] == MuonSystem->dtRechitClusterWheel[MuonSystem->nDtRechitClusters]) {
            dtRechitCluster_match_rpcBx.push_back(rpcBx[i]);
            MuonSystem->dtRechitCluster_match_RPChits_dPhi0p5[MuonSystem->nDtRechitClusters]++;
            if (rpcR < 470 && rpcR > 380 && abs(rpcZ[i]) < 661)
              MuonSystem->dtRechitCluster_match_RB1_dPhi0p5[MuonSystem->nDtRechitClusters]++;
          }
        }
        if (RazorAnalyzerMerged::deltaR(rpcEta[i], rpcPhi[i], MuonSystem->dtRechitClusterEta[MuonSystem->nDtRechitClusters], MuonSystem->dtRechitClusterPhi[MuonSystem->nDtRechitClusters]) < 0.4) {
          if (rpcR < 470 && rpcR > 380 && abs(rpcZ[i]) < 661)
            MuonSystem->dtRechitCluster_match_RB1_0p4[MuonSystem->nDtRechitClusters]++;
        }
      }
      int max_occurence = 0;
      int max_bx = -999;
      for (unsigned int l = 0; l < dtRechitCluster_match_rpcBx.size(); l++) {
        int counter = 0;
        for (unsigned int j = 0; j < dtRechitCluster_match_rpcBx.size(); j++) {
          if (dtRechitCluster_match_rpcBx[j] == dtRechitCluster_match_rpcBx[l])
            counter++;
        }
        if (counter > max_occurence) {
          max_occurence = counter;
          max_bx = dtRechitCluster_match_rpcBx[l];
        }
      }
      MuonSystem->dtRechitCluster_match_RPCBx_dPhi0p5[MuonSystem->nDtRechitClusters] = max_bx;

      MuonSystem->dtRechitClusterMet_dPhi[MuonSystem->nDtRechitClusters] = RazorAnalyzerMerged::deltaPhi(MuonSystem->dtRechitClusterPhi[MuonSystem->nDtRechitClusters], MuonSystem->metPhi);
      MuonSystem->dtRechitClusterPuppiMet_dPhi[MuonSystem->nDtRechitClusters] = RazorAnalyzerMerged::deltaPhi(MuonSystem->dtRechitClusterPhi[MuonSystem->nDtRechitClusters], MuonSystem->PuppimetPhi);

      MuonSystem->nDtRechitClusters++;
    }
    /* #endregion */

    // [15] Final per-event tree fill
    /* #region */
    if (!isData && signalScan) {
      pair<int, int> smsPair = make_pair(MuonSystem->mX, MuonSystem->ctau);
      Trees2D[smsPair]->Fill();
    } else {
      MuonSystem->tree_->Fill();
    }
    /* #endregion */
  }
  /* #endregion */

  // [16] End-of-job writeout for trees and histograms
  /* #region */
  if (!isData && signalScan) {
    for (auto& filePtr : Files2D) {
      cout << "Writing output tree (" << filePtr.second->GetName() << ")" << endl;
      filePtr.second->cd();
      Trees2D[filePtr.first]->Write();
      NEvents2D[filePtr.first]->Write("NEvents");
      Total2D[filePtr.first]->Write("Total");
      accep2D[filePtr.first]->Write("acceptance");
      accep_met2D[filePtr.first]->Write("acceptance_met");
      filePtr.second->Close();
    }
  } else if (!isData) {
    cout << "Filled Total of " << NEvents->GetBinContent(1) << " Events\n";
    cout << "Writing output trees..." << endl;
    outFile->cd();
    MuonSystem->tree_->Write();
    NEvents->Write();
    accep->Write("acceptance");
    accep_csccsc->Write("acceptance_csccsc");
    accep_cscdt->Write("acceptance_cscdt");
    accep_met->Write("acceptance_met");
    outFile->Close();
  }

  else {
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
  /* #endregion */
}
