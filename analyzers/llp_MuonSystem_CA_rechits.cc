//#include "fastjet/ClusterSequence.hh"
//#include "fastjet/ClusterSequence.hh"
//#include "fastjet/Selector.hh"
#include "llp_MuonSystem_CA_rechits.h"
#include "RazorHelper.h"
#include "TreeMuonSystemCombinationRechits.h"

#include "CACluster.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include "TMath.h"
#include "Math/Vector3D.h"
#include <iostream>
#include <random>
//C++ includes
#include "assert.h"

//ROOT includes
#include "TH1F.h"

//using namespace fastjet;
using namespace std::chrono;
using namespace std;
using namespace ROOT::Math;

struct greater_than_pt
{
  inline bool operator() (const TLorentzVector& p1, const TLorentzVector& p2){return p1.Pt() > p2.Pt();}
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

//lepton highest pt comparator
struct largest_pt
{
  inline bool operator() (const leptons& p1, const leptons& p2){return p1.lepton.Pt() > p2.lepton.Pt();}
} my_largest_pt;

//jet highest pt comparator
struct largest_pt_jet
{
  inline bool operator() (const jets& p1, const jets& p2){return p1.jet.Pt() > p2.jet.Pt();}
} my_largest_pt_jet;


void llp_MuonSystem_CA_rechits::Analyze(bool isData, int options, string outputfilename, string analysisTag)
{
  //initialization: create one TTree for each analysis box
  cout << "Initializing..." << endl;
  cout << "IsData = " << isData << "\n";
  cout << "options = " << options << "\n";

  //---------------------------
  //options format: MH/MX/ctau/condor: 1000/300/0/1
  // mh can be 3-4 digits, mx is always 3 digits, ctau is one digit(number of zeros), last digit is condor option
  // mh can be 3-4 digits, mx is always 3 digits, ctau is 2 digit(number of zeros), last digit is condor option
  //
  //
  // int mx = int(options/1000)%1000;
  // int mh = options/1000000;
  // int ctau = pow(10, int(options/10)%10) * int(int(options/100)%10);
  //
  // cout<<"mh "<<mh<<", mx "<<mx<<", ctau "<<ctau<<endl;

  bool signalScan = int(options/10) == 1;
  int option = options%10;
  // if (options % 1){
  //   option = 1; // used when running condor
  // }
  // else{
  //   option = 0;// used when running locally
  // }

  if( isData )
  {
    std::cout << "[INFO]: running on data with option: " << option << std::endl;
  }
  else
  {
    std::cout << "[INFO]: running on MC with option: " << option << std::endl;
  }
  if( signalScan )
  {
    std::cout << "[INFO]: running with Signal scan" << std::endl;
  }
  else
  {
    std::cout << "[INFO]: running without Signal scan " << option << std::endl;
  }



  const float ELE_MASS = 0.000511;
  const float MU_MASS  = 0.105658;
  const float Z_MASS   = 91.2;

  if (analysisTag == ""){
    analysisTag = "Razor2016_80X";

  }



  //-----------------------------------------------
  //Set up Output File
  //-----------------------------------------------
  string outfilename = outputfilename;
  if (outfilename == "") outfilename = "MuonSystem_Tree.root";
  TFile *outFile;
  if (isData || !signalScan) outFile = new TFile(outfilename.c_str(), "RECREATE");


  TreeMuonSystemCombinationRechits *MuonSystem = new TreeMuonSystemCombinationRechits;
  MuonSystem->CreateTree();
  MuonSystem->tree_->SetAutoFlush(0);
  MuonSystem->InitTree();

  // for signals, need one output file for each signal point
  map<pair<int,int>, TFile*> Files2D;
  map<pair<int,int>, TTree*>Trees2D;
  map<pair<int,int>, TH1F*> NEvents2D;
  map<pair<int,int>, TH1F*> accep2D;
  map<pair<int,int>, TH1F*> accep_met2D;
  map<pair<int,int>, TH1F*> Total2D;



  //histogram containing total number of processed events (for normalization)
  TH1F *NEvents = new TH1F("NEvents", "NEvents", 1, 1, 2);
  TH1F *Total = new TH1F("Total", "Total", 1, 1, 2);

  TH1F *accep = new TH1F("accep", "acceptance", 1, 1, 2);
  TH1F *accep_met = new TH1F("accep_met", "acceptance_met", 1, 1, 2);

  TH1F *Nmet200 = new TH1F("Nmet200", "Nmet200", 1, 1, 2);
  TH1F *NmetFilter = new TH1F("NmetFilter", "NmetFilter", 1, 1, 2);
  TH1F *Nlep0 = new TH1F("Nlep0", "Nlep0", 1, 1, 2);
  TH1F *Njet1 = new TH1F("Njet1", "Njet1", 1, 1, 2);
  TH1F *NcosmicVeto = new TH1F("NcosmicVeto", "NcosmicVeto", 1, 1, 2);

  //JetDefinition jet_def( antikt_algorithm, .4 );
  //fastjet::JetDefinition jet_def(fastjet::cambridge_algorithm, 0.4);

  //vector<fastjet::PseudoJet> input_particles;


  char* cmsswPath;
  cmsswPath = getenv("CMSSW_BASE");
  string pathname;
  if(cmsswPath != NULL) pathname = string(cmsswPath) + "/src/llp_analyzer/data/JEC/";
  if(cmsswPath != NULL and option == 1) pathname = "JEC/"; //run on condor if option == 1




  //--------------------------------
  //Initialize helper
  //--------------------------------
  RazorHelper *helper = 0;
  helper = new RazorHelper(analysisTag, isData);



  //*************************************************************************
  //Look over Input File Events
  //*************************************************************************
  if (fChain == 0) return;
  cout << "Total Events: " << fChain->GetEntries() << "\n";
  Long64_t nbytes = 0, nb = 0;
  clock_t start, end;
  start = clock();
  for (Long64_t jentry=0; jentry<fChain->GetEntries();jentry++) {
    //begin event
    if(jentry % 1000 == 0)
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
    if (ientry < 0) break;
    //GetEntry(ientry);
    nb = fChain->GetEntry(jentry); nbytes += nb;

    //fill normalization histogram
    MuonSystem->InitVariables();
    // std::cout << "deb1 " << jentry << std::endl;




    if (!isData && signalScan)
    {

      string mh_substring = lheComments->substr(lheComments->find("MH-")+3);
      int mh = stoi(mh_substring.substr(0,mh_substring.find('_')));
      string mx_substring = lheComments->substr(lheComments->find("MS-")+3);
      int mx = stoi(mx_substring.substr(0,mx_substring.find('_')));
      string ctau_substring = lheComments->substr(lheComments->find("ctauS-")+6);
      int ctau = stoi(ctau_substring.substr(0,ctau_substring.find('_')));
      MuonSystem->mH = mh;
      MuonSystem->mX = mx;
      MuonSystem->ctau = ctau;

      // if (mh2 != mh || mx2!=mx || ctau2!=ctau) continue;
      // cout<<*lheComments<<endl;

      pair<int,int> signalPair = make_pair(mx, ctau);

      if (Files2D.count(signalPair) == 0){ //create file and tree
        //format file name
        string thisFileName = outfilename;
        thisFileName.erase(thisFileName.end()-5, thisFileName.end());
        thisFileName += "_" + to_string(mx) + "_" + to_string(ctau) + ".root";

        Files2D[signalPair] = new TFile(thisFileName.c_str(), "recreate");
        Trees2D[signalPair] =  MuonSystem->tree_->CloneTree(0);
        NEvents2D[signalPair] = new TH1F(Form("NEvents%d%d", mx, ctau), "NEvents", 1,0.5,1.5);
        Total2D[signalPair] = new TH1F(Form("Total%d%d", mx, ctau), "Total", 1,0.5,1.5);
        accep2D[signalPair] = new TH1F(Form("accep2D%d%d", mx, ctau), "acceptance", 1,0.5,1.5);
        accep_met2D[signalPair] = new TH1F(Form("accep_met2D%d%d", mx, ctau), "acceptance_met", 1,0.5,1.5);



        cout << "Created new output file " << thisFileName << endl;
      }
          //Fill NEvents hist
      NEvents2D[signalPair]->Fill(1.0, genWeight);


    }
    //event info
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
       for(int i = 0; i < 2;i++)
       {
	       MuonSystem->gLLP_eta[MuonSystem->nGLLP] = gLLP_eta[i];
         MuonSystem->gLLP_phi[MuonSystem->nGLLP] = gLLP_phi[i];
         MuonSystem->gLLP_e[MuonSystem->nGLLP] = gLLP_e[i];
         MuonSystem->gLLP_pt[MuonSystem->nGLLP] = gLLP_pt[i];

         MuonSystem->gLLP_decay_vertex_r[MuonSystem->nGLLP] = sqrt(gLLP_decay_vertex_x[i]*gLLP_decay_vertex_x[i]+gLLP_decay_vertex_y[i]*gLLP_decay_vertex_y[i]);
         MuonSystem->gLLP_decay_vertex_x[MuonSystem->nGLLP] = gLLP_decay_vertex_x[i];
         MuonSystem->gLLP_decay_vertex_y[MuonSystem->nGLLP] = gLLP_decay_vertex_y[i];
         MuonSystem->gLLP_decay_vertex_z[MuonSystem->nGLLP] = gLLP_decay_vertex_z[i];
         float beta = gLLP_beta[i];
         float gLLP_decay_vertex = sqrt(pow(MuonSystem->gLLP_decay_vertex_r[i], 2) + pow(MuonSystem->gLLP_decay_vertex_z[i],2));
         float gamma = 1.0/sqrt(1-beta*beta);
         MuonSystem->gLLP_ctau[MuonSystem->nGLLP] = gLLP_decay_vertex/(beta * gamma);
         MuonSystem->gLLP_beta[MuonSystem->nGLLP] = gLLP_beta[i];


           if (abs(MuonSystem->gLLP_eta[i]) < 2.4
             && abs(MuonSystem->gLLP_decay_vertex_z[i])<1100 && abs(MuonSystem->gLLP_decay_vertex_z[i])>400
             && MuonSystem->gLLP_decay_vertex_r[i] < 695.5) MuonSystem->gLLP_csc[MuonSystem->nGLLP] = true;
           if (abs(MuonSystem->gLLP_decay_vertex_z[i])< 661.0
             && MuonSystem->gLLP_decay_vertex_r[i] < 800
              && MuonSystem->gLLP_decay_vertex_r[i] > 200.0) MuonSystem->gLLP_dt[MuonSystem->nGLLP] = true;

            MuonSystem->nGLLP++;
       }


    }//end of isData

      //get NPU
      MuonSystem->npv = nPV;
      MuonSystem->rho = fixedGridRhoFastjetAll;
      MuonSystem->met = metType1Pt;
      MuonSystem->metPhi = metType1Phi;

      if(signalScan && !isData)Total2D[make_pair(MuonSystem->mX, MuonSystem->ctau)]->Fill(1.0, genWeight*MuonSystem->pileupWeight);
      if(signalScan && !isData)
      {
        accep2D[make_pair(MuonSystem->mX, MuonSystem->ctau)]->Fill(1.0, genWeight*MuonSystem->pileupWeight);
      }
      else if (!isData)
      {
        accep->Fill(1.0, genWeight*MuonSystem->pileupWeight);

      }


      //Triggers
      for(int i = 0; i < NTriggersMAX; i++){
        MuonSystem->HLTDecision[i] = HLTDecision[i];
      }

      //*************************************************************************
      //Start Object Selection
      //*************************************************************************

      std::vector<leptons> Leptons;
      //-------------------------------
      //Muons
      //-------------------------------
      for( int i = 0; i < nMuons; i++ )
      {

        if(!muonIsLoose[i]) continue;
        if(muonPt[i] < 25) continue;
        if(fabs(muonEta[i]) > 2.4) continue;

        //remove overlaps
        bool overlap = false;
        for(auto& lep : Leptons)
        {
          if (RazorAnalyzer::deltaR(muonEta[i],muonPhi[i],lep.lepton.Eta(),lep.lepton.Phi()) < 0.3) overlap = true;
        }
        if(overlap) continue;

        leptons tmpMuon;
        tmpMuon.lepton.SetPtEtaPhiM(muonPt[i],muonEta[i], muonPhi[i], MU_MASS);
        tmpMuon.pdgId = 13 * -1 * muonCharge[i];
        tmpMuon.dZ = muon_dZ[i];
        tmpMuon.passId = muonIsTight[i];
        float muonIso = (muon_chargedIso[i] + fmax(0.0,  muon_photonIso[i] + muon_neutralHadIso[i] - 0.5*muon_pileupIso[i])) / muonPt[i];

        tmpMuon.passLooseIso = muonIso<0.25;
        tmpMuon.passTightIso = muonIso<0.15;
        tmpMuon.passVTightIso = muonIso<0.10;
        tmpMuon.passVVTightIso = muonIso<0.05;

        tmpMuon.passVetoId = false;
        Leptons.push_back(tmpMuon);
      }

      //-------------------------------
      //Electrons
      //-------------------------------
      for( int i = 0; i < nElectrons; i++ )
      {
        if (!ele_passCutBasedIDVeto[i]) continue;
        if(elePt[i] < 35) continue;
        if(fabs(eleEta[i]) > 2.5) continue;

        //remove overlaps
        bool overlap = false;
        for(auto& lep : Leptons)
        {
          if (RazorAnalyzer::deltaR(eleEta[i],elePhi[i],lep.lepton.Eta(),lep.lepton.Phi()) < 0.3) overlap = true;
        }
        if(overlap) continue;
        leptons tmpElectron;
        tmpElectron.lepton.SetPtEtaPhiM(elePt[i],eleEta[i], elePhi[i], ELE_MASS);
        tmpElectron.pdgId = 11 * -1 * eleCharge[i];
        tmpElectron.dZ = ele_dZ[i];
        tmpElectron.passId = ele_passCutBasedIDTight[i];
        Leptons.push_back(tmpElectron);
      }

      sort(Leptons.begin(), Leptons.end(), my_largest_pt);


      for ( auto &tmp : Leptons )
      {
        MuonSystem->lepE[MuonSystem->nLeptons]      = tmp.lepton.E();
        MuonSystem->lepPt[MuonSystem->nLeptons]     = tmp.lepton.Pt();
        MuonSystem->lepEta[MuonSystem->nLeptons]    = tmp.lepton.Eta();
        MuonSystem->lepPhi[MuonSystem->nLeptons]    = tmp.lepton.Phi();
        MuonSystem->lepPdgId[MuonSystem->nLeptons]  = tmp.pdgId;
        MuonSystem->lepDZ[MuonSystem->nLeptons]     = tmp.dZ;
        MuonSystem->lepTightId[MuonSystem->nLeptons] = tmp.passId;
        MuonSystem->lepPassLooseIso[MuonSystem->nLeptons] = tmp.passLooseIso;
        MuonSystem->lepPassTightIso[MuonSystem->nLeptons] = tmp.passTightIso;
        MuonSystem->lepPassVTightIso[MuonSystem->nLeptons] = tmp.passVTightIso;
        MuonSystem->lepPassVVTightIso[MuonSystem->nLeptons] = tmp.passVVTightIso;
        MuonSystem->nLeptons++;
      }



    //-----------------------------------------------
    //Select Jets
    //-----------------------------------------------

    std::vector<jets> Jets;

    for(int i = 0; i < nJets; i++)
    {
      if (fabs(jetEta[i]) >= 3.0)continue;
      if( jetPt[i] < 20 ) continue;
      if (!jetPassIDLoose[i]) continue;
      //------------------------------------------------------------
      //exclude selected muons and electrons from the jet collection
      //------------------------------------------------------------
      double deltaR = -1;
      for(auto& lep : Leptons){
        double thisDR = RazorAnalyzer::deltaR(jetEta[i],jetPhi[i],lep.lepton.Eta(),lep.lepton.Phi());
        if(deltaR < 0 || thisDR < deltaR) deltaR = thisDR;
      }
      if(deltaR > 0 && deltaR < 0.4) continue; //jet matches a selected lepton

      TLorentzVector thisJet = makeTLorentzVector( jetPt[i], jetEta[i], jetPhi[i], jetE[i] );

      jets tmpJet;
      tmpJet.jet    = thisJet;
      tmpJet.passId = jetPassIDTight[i];

      Jets.push_back(tmpJet);

      }

      sort(Jets.begin(), Jets.end(), my_largest_pt_jet);



      for ( auto &tmp : Jets )
      {
        if(tmp.jet.Pt()<30)continue;

        MuonSystem->jetE[MuonSystem->nJets] = tmp.jet.E();
        MuonSystem->jetPt[MuonSystem->nJets] = tmp.jet.Pt();
        MuonSystem->jetEta[MuonSystem->nJets] = tmp.jet.Eta();
        MuonSystem->jetPhi[MuonSystem->nJets] = tmp.jet.Phi();
        MuonSystem->jetTightPassId[MuonSystem->nJets] = tmp.passId;


        MuonSystem->nJets++;
      }


      MuonSystem->nDTRechits  = 0;

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

        if (dtRechitStation[i] == 1 && dtRechitWheel[i] == -2) nDTRechitsChamberMinus12++;
        if (dtRechitStation[i] == 1 && dtRechitWheel[i] == -1) nDTRechitsChamberMinus11++;
        if (dtRechitStation[i] == 1 && dtRechitWheel[i] == 0) nDTRechitsChamber10++;
        if (dtRechitStation[i] == 1 && dtRechitWheel[i] == 1) nDTRechitsChamberPlus11++;
        if (dtRechitStation[i] == 1 && dtRechitWheel[i] == 2) nDTRechitsChamberPlus12++;
        if (dtRechitStation[i] == 2 && dtRechitWheel[i] == -2) nDTRechitsChamberMinus22++;
        if (dtRechitStation[i] == 2 && dtRechitWheel[i] == -1) nDTRechitsChamberMinus21++;
        if (dtRechitStation[i] == 2 && dtRechitWheel[i] == 0) nDTRechitsChamber20++;
        if (dtRechitStation[i] == 2 && dtRechitWheel[i] == 1) nDTRechitsChamberPlus21++;
        if (dtRechitStation[i] == 2 && dtRechitWheel[i] == 2) nDTRechitsChamberPlus22++;
        if (dtRechitStation[i] == 3 && dtRechitWheel[i] == -2) nDTRechitsChamberMinus32++;
        if (dtRechitStation[i] == 3 && dtRechitWheel[i] == -1) nDTRechitsChamberMinus31++;
        if (dtRechitStation[i] == 3 && dtRechitWheel[i] == 0) nDTRechitsChamber30++;
        if (dtRechitStation[i] == 3 && dtRechitWheel[i] == 1) nDTRechitsChamberPlus31++;
        if (dtRechitStation[i] == 3 && dtRechitWheel[i] == 2) nDTRechitsChamberPlus32++;
        if (dtRechitStation[i] == 4 && dtRechitWheel[i] == -2) nDTRechitsChamberMinus42++;
        if (dtRechitStation[i] == 4 && dtRechitWheel[i] == -1) nDTRechitsChamberMinus41++;
        if (dtRechitStation[i] == 4 && dtRechitWheel[i] == 0) nDTRechitsChamber40++;
        if (dtRechitStation[i] == 4 && dtRechitWheel[i] == 1) nDTRechitsChamberPlus41++;
        if (dtRechitStation[i] == 4 && dtRechitWheel[i] == 2) nDTRechitsChamberPlus42++;




      }


      if ( nDTRechitsChamberMinus12 > 50) MuonSystem->nDtRings++;
      if ( nDTRechitsChamberMinus11 > 50) MuonSystem->nDtRings++;
      if ( nDTRechitsChamber10 > 50) MuonSystem->nDtRings++;
      if ( nDTRechitsChamberPlus11 > 50) MuonSystem->nDtRings++;
      if ( nDTRechitsChamberPlus12 > 50) MuonSystem->nDtRings++;
      if ( nDTRechitsChamberMinus22 > 50) MuonSystem->nDtRings++;
      if ( nDTRechitsChamberMinus21 > 50) MuonSystem->nDtRings++;
      if ( nDTRechitsChamber20 > 50) MuonSystem->nDtRings++;
      if ( nDTRechitsChamberPlus21 > 50) MuonSystem->nDtRings++;
      if ( nDTRechitsChamberPlus22 > 50) MuonSystem->nDtRings++;
      if ( nDTRechitsChamberMinus32 > 50) MuonSystem->nDtRings++;
      if ( nDTRechitsChamberMinus31 > 50) MuonSystem->nDtRings++;
      if ( nDTRechitsChamber30 > 50) MuonSystem->nDtRings++;
      if ( nDTRechitsChamberPlus31 > 50) MuonSystem->nDtRings++;
      if ( nDTRechitsChamberPlus32 > 50) MuonSystem->nDtRings++;
      if ( nDTRechitsChamberMinus42 > 50) MuonSystem->nDtRings++;
      if ( nDTRechitsChamberMinus41 > 50) MuonSystem->nDtRings++;
      if ( nDTRechitsChamber40 > 50) MuonSystem->nDtRings++;
      if ( nDTRechitsChamberPlus41 > 50) MuonSystem->nDtRings++;
      if ( nDTRechitsChamberPlus42 > 50) MuonSystem->nDtRings++;



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

	//Some printouts --> these look all empty
	//if(i==0)
	//{
	//std::cout<< "cscRechitsPhi[0] " << cscRechitsPhi[i]  <<std::endl;
	//std::cout<< "cscRechitsE[0] " << cscRechitsE[i]  <<std::endl;
	//std::cout<< "cscRechitsQuality[0] " << cscRechitsQuality[i]  <<std::endl;
	//std::cout<< "cscRechitsChannels[0] " << cscRechitsChannels[i]  <<std::endl;
	//std::cout<< "cscRechitsNStrips[0] " << cscRechitsNStrips[i]  <<std::endl;
	//std::cout<< "cscRechitsHitWire[0] " << cscRechitsHitWire[i]  <<std::endl;
	//std::cout<< "cscRechitsWGroupsBX[0] " << cscRechitsWGroupsBX[i]  <<std::endl;
	//std::cout<< "cscRechitsNWireGroups[0] " << cscRechitsNWireGroups[i]  <<std::endl;
	//std::cout<< "cscRechitsDetId[0] " << cscRechitsDetId[i]  <<std::endl;
	//}

        //pick out the right bits for chamber
        int chamber = ((cscRechitsDetId[i] >> 3) & 077); //https://github.com/cms-sw/cmssw/blob/master/DataFormats/MuonDetId/interface/CSCDetId.h#L147

        int layer = (cscRechitsDetId[i] & 07);
        Rechits p;
        p.phi = cscRechitsPhi[i];
        p.eta = cscRechitsEta[i];
        //p.e = cscRechitsE[i];
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

        if (cscRechitsChamber[i] == 11)  nCscRechitsChamberPlus11++;
        if (cscRechitsChamber[i] == 12)  nCscRechitsChamberPlus12++;
        if (cscRechitsChamber[i] == 13)  nCscRechitsChamberPlus13++;
        if (cscRechitsChamber[i] == 21)  nCscRechitsChamberPlus21++;
        if (cscRechitsChamber[i] == 22)  nCscRechitsChamberPlus22++;
        if (cscRechitsChamber[i] == 31)  nCscRechitsChamberPlus31++;
        if (cscRechitsChamber[i] == 32)  nCscRechitsChamberPlus32++;
        if (cscRechitsChamber[i] == 41)  nCscRechitsChamberPlus41++;
        if (cscRechitsChamber[i] == 42)  nCscRechitsChamberPlus42++;
        if (cscRechitsChamber[i] == -11)  nCscRechitsChamberMinus11++;
        if (cscRechitsChamber[i] == -12)  nCscRechitsChamberMinus12++;
        if (cscRechitsChamber[i] == -13)  nCscRechitsChamberMinus13++;
        if (cscRechitsChamber[i] == -21)  nCscRechitsChamberMinus21++;
        if (cscRechitsChamber[i] == -22)  nCscRechitsChamberMinus22++;
        if (cscRechitsChamber[i] == -31)  nCscRechitsChamberMinus31++;
        if (cscRechitsChamber[i] == -32)  nCscRechitsChamberMinus32++;
        if (cscRechitsChamber[i] == -41)  nCscRechitsChamberMinus41++;
        if (cscRechitsChamber[i] == -42)  nCscRechitsChamberMinus42++;
      }
      MuonSystem->nCscRings = 0;
      if ( nCscRechitsChamberPlus11 > 50) MuonSystem->nCscRings++;
      if ( nCscRechitsChamberPlus12 > 50) MuonSystem->nCscRings++;
      if ( nCscRechitsChamberPlus13 > 50) MuonSystem->nCscRings++;
      if ( nCscRechitsChamberPlus21 > 50) MuonSystem->nCscRings++;
      if ( nCscRechitsChamberPlus22 > 50) MuonSystem->nCscRings++;
      if ( nCscRechitsChamberPlus31 > 50) MuonSystem->nCscRings++;
      if ( nCscRechitsChamberPlus32 > 50) MuonSystem->nCscRings++;
      if ( nCscRechitsChamberPlus41 > 50) MuonSystem->nCscRings++;
      if ( nCscRechitsChamberPlus42 > 50) MuonSystem->nCscRings++;
      if ( nCscRechitsChamberMinus11 > 50) MuonSystem->nCscRings++;
      if ( nCscRechitsChamberMinus12 > 50) MuonSystem->nCscRings++;
      if ( nCscRechitsChamberMinus13 > 50) MuonSystem->nCscRings++;
      if ( nCscRechitsChamberMinus21 > 50) MuonSystem->nCscRings++;
      if ( nCscRechitsChamberMinus22 > 50) MuonSystem->nCscRings++;
      if ( nCscRechitsChamberMinus31 > 50) MuonSystem->nCscRings++;
      if ( nCscRechitsChamberMinus32 > 50) MuonSystem->nCscRings++;
      if ( nCscRechitsChamberMinus41 > 50) MuonSystem->nCscRings++;
      if ( nCscRechitsChamberMinus42 > 50) MuonSystem->nCscRings++;
      //Do DBSCAN Clustering

      int min_point = 50;  //minimum number of Rechitss to call it a cluster
      float epsilon = 0.2; //cluster radius parameter
      CACluster ds(min_point, epsilon, points);
      vector<Rechits> clustered_csc_rechits = ds.run_rechits();
      ds.clusterProperties();
      //ds.merge_clusters();
      //ds.clusterProperties();
      ds.sort_clusters();

      MuonSystem->nCscRechitClusters = 0;
      for ( auto &tmp : ds.clusters  ) {
          MuonSystem->cscRechitClusterX[MuonSystem->nCscRechitClusters] =tmp.x;
          MuonSystem->cscRechitClusterY[MuonSystem->nCscRechitClusters] =tmp.y;
          MuonSystem->cscRechitClusterZ[MuonSystem->nCscRechitClusters] =tmp.z;
          MuonSystem->cscRechitClusterTimeWeighted[MuonSystem->nCscRechitClusters] = tmp.tWeighted;
          MuonSystem->cscRechitClusterTimeSpreadWeightedAll[MuonSystem->nCscRechitClusters] = tmp.TSpreadWeightedAll;
          
          MuonSystem->cscRechitClusterTime[MuonSystem->nCscRechitClusters] = tmp.tTotal;
          MuonSystem->cscRechitClusterTimeSpread[MuonSystem->nCscRechitClusters] = tmp.TSpread;

          MuonSystem->cscRechitClusterEta[MuonSystem->nCscRechitClusters] =tmp.eta;
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
          MuonSystem->cscRechitClusterMaxChamberRatio[MuonSystem->nCscRechitClusters] = 1.0*tmp.maxChamberRechits/tmp.nhits;
          MuonSystem->cscRechitClusterNChamber[MuonSystem->nCscRechitClusters] = tmp.nChamber;
          MuonSystem->cscRechitClusterMaxStation[MuonSystem->nCscRechitClusters] = tmp.maxStation;
          MuonSystem->cscRechitClusterMaxStationRatio[MuonSystem->nCscRechitClusters] = 1.0*tmp.maxStationRechits/tmp.nhits;

          MuonSystem->cscRechitClusterNStation10[MuonSystem->nCscRechitClusters] = tmp.nStation10;
          MuonSystem->cscRechitClusterAvgStation10[MuonSystem->nCscRechitClusters] = tmp.avgStation10;





          //Jet veto/ muon veto
          MuonSystem->cscRechitClusterJetVetoPt[MuonSystem->nCscRechitClusters] = 0.0;
          MuonSystem->cscRechitClusterJetVetoE[MuonSystem->nCscRechitClusters] = 0.0;
          MuonSystem->cscRechitClusterMuonVetoPt[MuonSystem->nCscRechitClusters] = 0.0;
          MuonSystem->cscRechitClusterMuonVetoE[MuonSystem->nCscRechitClusters] = 0.0;


          // jet veto
          for(int i = 0; i < nJets; i++)
          {
            if (fabs(jetEta[i]>3.0)) continue;
            if (RazorAnalyzer::deltaR(jetEta[i], jetPhi[i], MuonSystem->cscRechitClusterEta[MuonSystem->nCscRechitClusters],MuonSystem->cscRechitClusterPhi[MuonSystem->nCscRechitClusters]) < 0.4 && jetPt[i] > MuonSystem->cscRechitClusterJetVetoPt[MuonSystem->nCscRechitClusters] ) {
              MuonSystem->cscRechitClusterJetVetoPt[MuonSystem->nCscRechitClusters]  = jetPt[i];
              MuonSystem->cscRechitClusterJetVetoE[MuonSystem->nCscRechitClusters]  = jetE[i];
              MuonSystem->cscRechitClusterJetVetoTightId[MuonSystem->nCscRechitClusters]  = jetPassIDTight[i];
              MuonSystem->cscRechitClusterJetVetoLooseId[MuonSystem->nCscRechitClusters]  = jetPassIDLoose[i];

            }

          }
          float min_deltaR = 15.;
          int index = 999;


          for(int i = 0; i < nMuons; i++)
          {
            if (fabs(muonEta[i]>3.0)) continue;
            float muonIso = (muon_chargedIso[i] + fmax(0.0,  muon_photonIso[i] + muon_neutralHadIso[i] - 0.5*muon_pileupIso[i])) / muonPt[i];
            if (RazorAnalyzer::deltaR(muonEta[i], muonPhi[i], MuonSystem->cscRechitClusterEta[MuonSystem->nCscRechitClusters],MuonSystem->cscRechitClusterPhi[MuonSystem->nCscRechitClusters]) < 0.4 && muonPt[i] > MuonSystem->cscRechitClusterMuonVetoPt[MuonSystem->nCscRechitClusters] ) {
              MuonSystem->cscRechitClusterMuonVetoPt[MuonSystem->nCscRechitClusters]  = muonPt[i];
              MuonSystem->cscRechitClusterMuonVetoE[MuonSystem->nCscRechitClusters]  = muonE[i];
              MuonSystem->cscRechitClusterMuonVetoGlobal[MuonSystem->nCscRechitClusters]  = muon_isGlobal[i];
              MuonSystem->cscRechitClusterMuonVetoLooseId[MuonSystem->nCscRechitClusters]  = muonIsLoose[i];


            }
          }
          if(!isData)
          {
            // match to gen level LLP
            min_deltaR = 15.;
            index = 999;
            for(int j = 0; j < MuonSystem->nGLLP;j++)
            {

              double current_delta_r = RazorAnalyzer::deltaR(MuonSystem->cscRechitClusterEta[MuonSystem->nCscRechitClusters], MuonSystem->cscRechitClusterPhi[MuonSystem->nCscRechitClusters], MuonSystem->gLLP_eta[j], MuonSystem->gLLP_phi[j]);
              if (current_delta_r < min_deltaR)
              {
                min_deltaR = current_delta_r;
                index = j;
              }
            }
            if (min_deltaR < 0.4)MuonSystem->cscRechitCluster_match_gLLP[MuonSystem->nCscRechitClusters] = true;
            else MuonSystem->cscRechitCluster_match_gLLP[MuonSystem->nCscRechitClusters] = false;


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
            if (RazorAnalyzer::deltaR(dtSegEta[i], dtSegPhi[i], MuonSystem->cscRechitClusterEta[MuonSystem->nCscRechitClusters],MuonSystem->cscRechitClusterPhi[MuonSystem->nCscRechitClusters]) < 0.4 )
            {
              MuonSystem->cscRechitCluster_match_dtSeg_0p4[MuonSystem->nCscRechitClusters] ++;
              if (dtSegStation[i] == 1) MuonSystem->cscRechitCluster_match_MB1Seg_0p4[MuonSystem->nCscRechitClusters] ++;
            }
          }

          //match to RPC hits in RE1/2
          for (int i = 0; i < nRpc; i++) {
            float rpcR = sqrt(rpcX[i]*rpcX[i] + rpcY[i]*rpcY[i]);
            if (RazorAnalyzer::deltaR(rpcEta[i], rpcPhi[i], MuonSystem->cscRechitClusterEta[MuonSystem->nCscRechitClusters],MuonSystem->cscRechitClusterPhi[MuonSystem->nCscRechitClusters]) < 0.4 )
            {
              if (rpcR < 461.0 && rpcR > 275 && abs(rpcZ[i]) > 663 && abs(rpcZ[i]) < 730) MuonSystem->cscRechitCluster_match_RE12_0p4[MuonSystem->nCscRechitClusters] ++;
              if (rpcR < 470 && rpcR > 380 && abs(rpcZ[i]) < 661) MuonSystem->cscRechitCluster_match_RB1_0p4[MuonSystem->nCscRechitClusters] ++;
            }
          }

          MuonSystem->cscRechitClusterMet_dPhi[MuonSystem->nCscRechitClusters] =  RazorAnalyzer::deltaPhi(MuonSystem->cscRechitClusterPhi[MuonSystem->nCscRechitClusters],MuonSystem->metPhi);
          MuonSystem->cscRechitClusterIndex[MuonSystem->nCscRechitClusters] =  MuonSystem->nCscRechitClusters;
          MuonSystem->nCscRechitClusters++;
      }

      //L//std::cout << " n. of clusters in event after loop: " << MuonSystem->nCscRechitClusters << std::endl;
      //Here loop on points
      //MuonSystem->nCscClusteredRechits = clustered_csc_rechits.size();

      if (MuonSystem->nCscRechitClusters>0)
	{
	  for (unsigned int i = 0; i < clustered_csc_rechits.size(); i++)
	    {

	      MuonSystem->cscClusteredRechitsPhi[i] = clustered_csc_rechits.at(i).phi;
	      MuonSystem->cscClusteredRechitsEta[i] = clustered_csc_rechits.at(i).eta;
	      MuonSystem->cscClusteredRechitsX[i] = clustered_csc_rechits.at(i).x;
	      MuonSystem->cscClusteredRechitsY[i] = clustered_csc_rechits.at(i).y;
	      MuonSystem->cscClusteredRechitsZ[i] = clustered_csc_rechits.at(i).z;
	      //MuonSystem->cscClusteredRechitsE[i] = clustered_csc_rechits.at(i).e;
	      MuonSystem->cscClusteredRechitsTpeak[i] = clustered_csc_rechits.at(i).t;
	      MuonSystem->cscClusteredRechitsTwire[i] = clustered_csc_rechits.at(i).twire;
	      //MuonSystem->cscRechitsQuality[i] = cscRechitsQuality[i];
	      MuonSystem->cscClusteredRechitsChamber[i] = clustered_csc_rechits.at(i).chamber;
	      MuonSystem->cscClusteredRechitsStation[i] = clustered_csc_rechits.at(i).station;
	      //MuonSystem->cscRechitsChannels[i] = cscRechitsChannels[i];
	      //MuonSystem->cscRechitsNStrips[i] = cscRechitsNStrips[i];
	      //MuonSystem->cscRechitsHitWire[i] = cscRechitsHitWire[i];
	      //MuonSystem->cscRechitsWGroupsBX[i] = cscRechitsWGroupsBX[i];
	      //MuonSystem->cscRechitsNWireGroups[i] = cscRechitsNWireGroups[i];
	      //MuonSystem->cscRechitsDetId[i] = cscRechitsDetId[i];

	      //std::cout<<"point cluster ID "<< clustered_csc_rechits.at(i).clusterID << std::endl;
	      MuonSystem->cscClusteredRechitsIndex[i] = clustered_csc_rechits.at(i).clusterID;
	      MuonSystem->nCscClusteredRechits++;

	    }
	}


      // DT cluster

      points.clear();
      // cout<<"here"<<endl;

      for (int i = 0; i < nDtRechits; i++) {
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

      //Do DBSCAN Clustering
      int min_point_dt = 50;  //minimum number of segments to call it a cluster
      float epsilon_dt = 0.2; //cluster radius parameter
      CACluster ds_dtRechit(min_point_dt, epsilon_dt, points);
      vector<Rechits> clustered_dt_rechits = ds_dtRechit.run_rechits();
      //ds_dtRechit.run();
      ds_dtRechit.clusterProperties();
      //ds_dtRechit.merge_clusters();
      //ds_dtRechit.clusterProperties();

      ds_dtRechit.sort_clusters();


      MuonSystem->nDtRechitClusters = 0;

      for ( auto &tmp : ds_dtRechit.clusters  ) {

        //remove overlaps
        bool overlap = false;
        for(int i = 0; i < MuonSystem->nCscRechitClusters; i++)
        {
          if (RazorAnalyzer::deltaR(MuonSystem->cscRechitClusterEta[i],MuonSystem->cscRechitClusterPhi[i],tmp.eta, tmp.phi)<0.4) overlap = true;
        }
        if (overlap) continue;


          MuonSystem->dtRechitClusterX[MuonSystem->nDtRechitClusters] =tmp.x;
          MuonSystem->dtRechitClusterY[MuonSystem->nDtRechitClusters] =tmp.y;
          MuonSystem->dtRechitClusterZ[MuonSystem->nDtRechitClusters] =tmp.z;
          if (abs(tmp.z) < 126.8) MuonSystem->dtRechitClusterWheel[MuonSystem->nDtRechitClusters] = 0;
          else if (tmp.z > 126.8 && tmp.z < 395.4) MuonSystem->dtRechitClusterWheel[MuonSystem->nDtRechitClusters] = 1;
          else if (tmp.z < -126.8 && tmp.z > -395.4)MuonSystem->dtRechitClusterWheel[MuonSystem->nDtRechitClusters] = -1;
          else if (tmp.z<0) MuonSystem->dtRechitClusterWheel[MuonSystem->nDtRechitClusters] = -2;
          else MuonSystem->dtRechitClusterWheel[MuonSystem->nDtRechitClusters] = 2;
          MuonSystem->dtRechitClusterEta[MuonSystem->nDtRechitClusters] =tmp.eta;
          MuonSystem->dtRechitClusterPhi[MuonSystem->nDtRechitClusters] =tmp.phi;
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
          default_random_engine generator (seed);

	        // default_random_engine generator;
    	   uniform_real_distribution<double> distribution(0.0,1.0);
	       float prob = 0.03;
         for (int i=0; i<12; ++i) {
           if ( distribution(generator) < prob) MuonSystem->dtRechitClusterNoiseHitStation1[MuonSystem->nDtRechitClusters]++;
         }
         for (int i=0; i<12; ++i) {
           if ( distribution(generator) < prob) MuonSystem->dtRechitClusterNoiseHitStation2[MuonSystem->nDtRechitClusters]++;
      	 }
         for (int i=0; i<12; ++i) {
           if ( distribution(generator) < prob) MuonSystem->dtRechitClusterNoiseHitStation3[MuonSystem->nDtRechitClusters]++;
      	 }
         for (int i=0; i<8; ++i) {
           if ( distribution(generator) < prob) MuonSystem->dtRechitClusterNoiseHitStation4[MuonSystem->nDtRechitClusters]++;
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
        	MuonSystem->dtRechitClusterMaxChamberRatio[MuonSystem->nDtRechitClusters] = 1.0*tmp.maxChamberRechits/tmp.nhits;
        	MuonSystem->dtRechitClusterNChamber[MuonSystem->nDtRechitClusters] = tmp.nChamber;
        	MuonSystem->dtRechitClusterMaxStation[MuonSystem->nDtRechitClusters] = tmp.maxStation;
        	MuonSystem->dtRechitClusterMaxStationRatio[MuonSystem->nDtRechitClusters] = 1.0*tmp.maxStationRechits/tmp.nhits;
          MuonSystem->dtRechitClusterNStation10[MuonSystem->nDtRechitClusters] = tmp.nStation10;
          MuonSystem->dtRechitClusterAvgStation10[MuonSystem->nDtRechitClusters] = tmp.avgStation10;



          //Jet veto/ muon veto
          MuonSystem->dtRechitClusterJetVetoPt[MuonSystem->nDtRechitClusters] = 0.0;
          MuonSystem->dtRechitClusterJetVetoE[MuonSystem->nDtRechitClusters] = 0.0;
          MuonSystem->dtRechitClusterMuonVetoPt[MuonSystem->nDtRechitClusters] = 0.0;
          MuonSystem->dtRechitClusterMuonVetoE[MuonSystem->nDtRechitClusters] = 0.0;


          // jet veto
          for(int i = 0; i < nJets; i++)
          {
            if (fabs(jetEta[i]>3.0)) continue;
            if (RazorAnalyzer::deltaR(jetEta[i], jetPhi[i], MuonSystem->dtRechitClusterEta[MuonSystem->nDtRechitClusters],MuonSystem->dtRechitClusterPhi[MuonSystem->nDtRechitClusters]) < 0.4 && jetPt[i] > MuonSystem->dtRechitClusterJetVetoPt[MuonSystem->nDtRechitClusters] ) {
              MuonSystem->dtRechitClusterJetVetoPt[MuonSystem->nDtRechitClusters]  = jetPt[i];
              MuonSystem->dtRechitClusterJetVetoE[MuonSystem->nDtRechitClusters]  = jetE[i];
              MuonSystem->dtRechitClusterJetVetoLooseId[MuonSystem->nDtRechitClusters]  = jetPassIDLoose[i];
              MuonSystem->dtRechitClusterJetVetoTightId[MuonSystem->nDtRechitClusters]  = jetPassIDTight[i];

            }


          }


          for(int i = 0; i < nMuons; i++)
          {
            if (fabs(muonEta[i]>3.0)) continue;
            float muonIso = (muon_chargedIso[i] + fmax(0.0,  muon_photonIso[i] + muon_neutralHadIso[i] - 0.5*muon_pileupIso[i])) / muonPt[i];
            if (RazorAnalyzer::deltaR(muonEta[i], muonPhi[i], MuonSystem->dtRechitClusterEta[MuonSystem->nDtRechitClusters],MuonSystem->dtRechitClusterPhi[MuonSystem->nDtRechitClusters]) < 0.4 && muonPt[i] > MuonSystem->dtRechitClusterMuonVetoPt[MuonSystem->nDtRechitClusters] ) {
              MuonSystem->dtRechitClusterMuonVetoPt[MuonSystem->nDtRechitClusters]  = muonPt[i];
              MuonSystem->dtRechitClusterMuonVetoE[MuonSystem->nDtRechitClusters]  = muonE[i];
              MuonSystem->dtRechitClusterMuonVetoGlobal[MuonSystem->nDtRechitClusters]  = muon_isGlobal[i];
              MuonSystem->dtRechitClusterMuonVetoLooseId[MuonSystem->nDtRechitClusters]  = muonIsLoose[i];


            }
          }

          for (int i = 0; i < nDtSeg; i++) {
              if (RazorAnalyzer::deltaR(dtSegEta[i], dtSegPhi[i], MuonSystem->dtRechitClusterEta[MuonSystem->nDtRechitClusters],MuonSystem->dtRechitClusterPhi[MuonSystem->nDtRechitClusters]) < 0.4) {
                if (dtSegStation[i] == 1) MuonSystem->dtRechitClusterNSegStation1[MuonSystem->nDtRechitClusters]  +=1;
                if (dtSegStation[i] == 2) MuonSystem->dtRechitClusterNSegStation2[MuonSystem->nDtRechitClusters]  +=1;
                if (dtSegStation[i] == 3) MuonSystem->dtRechitClusterNSegStation3[MuonSystem->nDtRechitClusters]  +=1;
                if (dtSegStation[i] == 4) MuonSystem->dtRechitClusterNSegStation4[MuonSystem->nDtRechitClusters]  +=1;
              }
              if (abs(RazorAnalyzer::deltaPhi(dtSegPhi[i],MuonSystem->dtRechitClusterPhi[MuonSystem->nDtRechitClusters]))>2) {
                if (dtSegStation[i] == 1) MuonSystem->dtRechitClusterNOppositeSegStation1[MuonSystem->nDtRechitClusters]  +=1;
                if (dtSegStation[i] == 2) MuonSystem->dtRechitClusterNOppositeSegStation2[MuonSystem->nDtRechitClusters]  +=1;
                if (dtSegStation[i] == 3) MuonSystem->dtRechitClusterNOppositeSegStation3[MuonSystem->nDtRechitClusters]  +=1;
                if (dtSegStation[i] == 4) MuonSystem->dtRechitClusterNOppositeSegStation4[MuonSystem->nDtRechitClusters]  +=1;
              }
         }



          // match to gen-level LLP
          float min_deltaR = 15.;
          int index = 999;
          if (!isData)
          {
            for(int j = 0; j < MuonSystem->nGLLP;j++)
            {
              double current_delta_r = RazorAnalyzer::deltaR(tmp.eta, tmp.phi, MuonSystem->gLLP_eta[j], MuonSystem->gLLP_phi[j]);
              if (current_delta_r < min_deltaR)
              {
                min_deltaR = current_delta_r;
                index = j;
              }
            }
            if (min_deltaR < 0.4)MuonSystem->dtRechitCluster_match_gLLP[MuonSystem->nDtRechitClusters] = true;
            else MuonSystem->dtRechitCluster_match_gLLP[MuonSystem->nDtRechitClusters] = false;

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
          MuonSystem->nCscRechits = ncscRechits;

          for (int i = 0; i < nDtRechits; i++) {
            if (RazorAnalyzer::deltaR(dtRechitCorrectEta[i], dtRechitCorrectPhi[i], MuonSystem->dtRechitClusterEta[MuonSystem->nDtRechitClusters],MuonSystem->dtRechitClusterPhi[MuonSystem->nDtRechitClusters]) < 0.5 )
            {
              if (dtRechitStation[i] == 1) MuonSystem->dtRechitCluster_match_MB1hits_0p5[MuonSystem->nDtRechitClusters] ++;
            }
            if (RazorAnalyzer::deltaR(dtRechitCorrectEta[i], dtRechitCorrectPhi[i], MuonSystem->dtRechitClusterEta[MuonSystem->nDtRechitClusters],MuonSystem->dtRechitClusterPhi[MuonSystem->nDtRechitClusters]) < 0.4 )
            {
              if (dtRechitStation[i] == 1) MuonSystem->dtRechitCluster_match_MB1hits_0p4[MuonSystem->nDtRechitClusters] ++;
            }
            if(abs(dtRechitWheel[i]-MuonSystem->dtRechitClusterWheel[MuonSystem->nDtRechitClusters])==1 && dtRechitStation[i] == 1)
            {
              if (abs(RazorAnalyzer::deltaPhi(dtRechitCorrectPhi[i], MuonSystem->dtRechitClusterPhi[MuonSystem->nDtRechitClusters])) < TMath::Pi()/4.0 )
              {
                if (dtRechitWheel[i]-MuonSystem->dtRechitClusterWheel[MuonSystem->nDtRechitClusters] == 1) MuonSystem->dtRechitCluster_match_MB1hits_cosmics_plus[MuonSystem->nDtRechitClusters] ++;
                else MuonSystem->dtRechitCluster_match_MB1hits_cosmics_minus[MuonSystem->nDtRechitClusters] ++;
              }
            }

          }

         std::vector<int> dtRechitCluster_match_rpcBx;

         //match to RPC hits with dPhi<0.5 and same wheel in DT
         for (int i = 0; i < nRpc; i++) {
           float rpcR = sqrt(rpcX[i]*rpcX[i] + rpcY[i]*rpcY[i]);
           if (rpcRegion[i]!=0) continue;
           if (abs(RazorAnalyzer::deltaPhi(rpcPhi[i], MuonSystem->dtRechitClusterPhi[MuonSystem->nDtRechitClusters])) < 0.5 )
           {
             if (rpcRing[i] == MuonSystem->dtRechitClusterWheel[MuonSystem->nDtRechitClusters])
             {
               dtRechitCluster_match_rpcBx.push_back(rpcBx[i]);
               MuonSystem->dtRechitCluster_match_RPChits_dPhi0p5[MuonSystem->nDtRechitClusters]++;
               if (rpcR < 470 && rpcR > 380 && abs(rpcZ[i]) < 661)MuonSystem->dtRechitCluster_match_RB1_dPhi0p5[MuonSystem->nDtRechitClusters] ++;

             }
           }
           if(RazorAnalyzer::deltaR(rpcEta[i], rpcPhi[i], MuonSystem->dtRechitClusterEta[MuonSystem->nDtRechitClusters], MuonSystem->dtRechitClusterPhi[MuonSystem->nDtRechitClusters]) < 0.4 )
           {
             if (rpcR < 470 && rpcR > 380 && abs(rpcZ[i]) < 661)MuonSystem->dtRechitCluster_match_RB1_0p4[MuonSystem->nDtRechitClusters] ++;
           }
         }
         int max_occurence = 0;
         int max_bx = -999;
         for (unsigned int l = 0; l < dtRechitCluster_match_rpcBx.size(); l++)
         {
           int counter = 0;
           for(unsigned int j = 0; j < dtRechitCluster_match_rpcBx.size(); j ++)
           {
             if (dtRechitCluster_match_rpcBx[j] == dtRechitCluster_match_rpcBx[l]) counter++;
           }
           if (counter>max_occurence)
           {
             max_occurence = counter;
             max_bx = dtRechitCluster_match_rpcBx[l];
           }
         }
          MuonSystem->dtRechitCluster_match_RPCBx_dPhi0p5[MuonSystem->nDtRechitClusters] = max_bx;

          MuonSystem->dtRechitClusterMet_dPhi[MuonSystem->nDtRechitClusters] =  RazorAnalyzer::deltaPhi(MuonSystem->dtRechitClusterPhi[MuonSystem->nDtRechitClusters],MuonSystem->metPhi);

          MuonSystem->dtRechitClusterIndex[MuonSystem->nDtRechitClusters] =  MuonSystem->nDtRechitClusters;

          MuonSystem->nDtRechitClusters++;
        }

      //if (isData && MuonSystem->nDtRechitClusters + MuonSystem->nCscRechitClusters < 1) continue;
      //L
      if (MuonSystem->nDtRechitClusters>0)
	{
	  for (unsigned int i = 0; i < clustered_dt_rechits.size(); i++)
	    {

	      MuonSystem->dtClusteredRechitCorrectPhi[i] = clustered_dt_rechits.at(i).phi;
	      MuonSystem->dtClusteredRechitCorrectEta[i] = clustered_dt_rechits.at(i).eta;
	      MuonSystem->dtClusteredRechitCorrectX[i] = clustered_dt_rechits.at(i).x;
	      MuonSystem->dtClusteredRechitCorrectY[i] = clustered_dt_rechits.at(i).y;
	      MuonSystem->dtClusteredRechitCorrectZ[i] = clustered_dt_rechits.at(i).z;
	      MuonSystem->dtClusteredRechitTime[i] = clustered_dt_rechits.at(i).t;
	      MuonSystem->dtClusteredRechitStation[i] = clustered_dt_rechits.at(i).station;
	      MuonSystem->dtClusteredRechitWheel[i] = clustered_dt_rechits.at(i).wheel;
	      MuonSystem->dtClusteredRechitLayer[i] = clustered_dt_rechits.at(i).layer;
	      MuonSystem->dtClusteredRechitSuperLayer[i] = clustered_dt_rechits.at(i).superlayer;
	      MuonSystem->dtClusteredRechitIndex[i] = clustered_dt_rechits.at(i).clusterID;
	      MuonSystem->nDtClusteredRechits++;

	    }
	}
      //L

      if(!isData && signalScan)
      {
        pair<int,int> smsPair = make_pair(MuonSystem->mX, MuonSystem->ctau);
        Trees2D[smsPair]->Fill();
      }
      else
      {
        MuonSystem->tree_->Fill();
      }


    }
      if(!isData && signalScan)
      {
        for(auto &filePtr : Files2D)
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
