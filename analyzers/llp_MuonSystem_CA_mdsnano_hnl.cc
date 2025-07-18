//#include "fastjet/ClusterSequence.hh"
//#include "fastjet/ClusterSequence.hh"
//#include "fastjet/Selector.hh"
#include "llp_MuonSystem_CA_mdsnano_hnl.h"
#include "RazorHelper.h"
#include "TreeMuonSystemCombination.h"

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


void llp_MuonSystem_CA_mdsnano_hnl::Analyze(bool isData, int options, string outputfilename, string analysisTag)
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


  TreeMuonSystemCombination *MuonSystem = new TreeMuonSystemCombination;
  MuonSystem->CreateTree();
  MuonSystem->tree_->SetAutoFlush(0);
  MuonSystem->InitTree();





  //histogram containing total number of processed events (for normalization)
  TH1F *NEvents = new TH1F("NEvents", "NEvents", 1, 1, 2);
  TH1F *Total = new TH1F("Total", "Total", 1, 1, 2);

  TH1F *accep = new TH1F("accep", "acceptance", 1, 1, 2);
  TH1F *accep_csccsc = new TH1F("accep_csccsc", "accep_csccsc", 1, 1, 2);
  TH1F *accep_cscdt = new TH1F("accep_cscdt", "accep_cscdt", 1, 1, 2);
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
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry); nbytes += nb;

    //fill normalization histogram
    MuonSystem->InitVariables();
    
    auto nCscSeg = ncscSegments;
    auto nDtSeg = ndtSegments;
    auto nDtRechits = ndtRecHits;
    auto nRpc = nrpcRecHits;
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

    Float_t jetE[150] = {0}; // jetE
    for (int i = 0; i < nJet; ++i)
    {
      auto eta = Jet_eta[i];
      auto pt = Jet_pt[i];
      auto pz = pt * TMath::SinH(eta);
      auto mass = Jet_mass[i];
      jetE[i] = TMath::Sqrt(mass * mass + pt * pt + pz * pz);
    }
    auto jetEta = Jet_eta;         // jetEta
    bool jetPassIDTightLepVeto[150] = {0}; // jetPassIDLoose
    bool jetPassIDTight[150] = {0}; // jetPassIDTight
    for (int i = 0; i < nJet; ++i)
    {
      jetPassIDTight[i] = helper->jetTightLepVeto(analysisTag, false, Jet_neHEF[i], Jet_neEmEF[i], Jet_chEmEF[i], Jet_muEF[i], Jet_chHEF[i], Jet_chMultiplicity[i], Jet_neMultiplicity[i], Jet_eta[i], Jet_jetId[i]);
      jetPassIDTightLepVeto[i] = helper->jetTightLepVeto(analysisTag, true, Jet_neHEF[i], Jet_neEmEF[i], Jet_chEmEF[i], Jet_muEF[i], Jet_chHEF[i], Jet_chMultiplicity[i], Jet_neMultiplicity[i], Jet_eta[i], Jet_jetId[i]);
      // cout<<jetPassIDTight[i]<<","<<jetPassIDTightLepVeto[i]<<", "<< Jet_neHEF[i]<<", "<<Jet_neEmEF[i]<<", "<<Jet_chEmEF[i]<<", "<<Jet_muEF[i]<<", "<<Jet_chHEF[i]<<", "<<Jet_chMultiplicity[i]<<", "<<Jet_neMultiplicity[i]<<", "<< Jet_eta[i]<<", "<< Jet_jetId[i]<<endl;
    }
    auto jetPhi = Jet_phi;          // jetPhi
    auto jetPt = Jet_pt;            // jetPt
    auto lumiNum = luminosityBlock; // lumiNum
    auto metType1Phi = MET_phi;     // metType1Phi
    auto metType1Pt = MET_pt;       // metType1Pt
    if (analysisTag == "Summer24"){
      metType1Phi = PFMET_phi;     // metType1Phi
      metType1Pt = PFMET_pt;       // metType1Pt
    }
    
    auto muonCharge = Muon_charge;  // muonCharge
    Float_t muonE[50] = {0};        // muonE
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
    auto muon_pfRelIso04_all = Muon_pfRelIso04_all; // muon_chargedIso
    auto muon_dZ = Muon_dz;                     // muon_dZ
    auto muon_isGlobal = Muon_isGlobal;         // muon_isGlobal
    auto nElectrons = nElectron;                // nElectrons
    auto nJets = nJet;                          // nJets
    auto nMuons = nMuon;                        // nMuons
    int nPV = PV_npvs;                          // nPV
    auto runNum = run;                          // runNum

    auto *lheComments = (string *)"123";
    auto genWeight = Generator_weight; // genWeight


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

    bool llp_flag = false;
    if (!isData)
    {

       //for DS model
       MuonSystem->nGenParticles = 0;
       int nPromptTau=0;
      for(int i = 0;i < nGenPart; i++)
      {
        if (abs(GenPart_pdgId[i])>=999999)
        {
          // cout<<eventNum<<","<<gParticleId[i]<<","<<MuonSystem->nGLLP<<endl;
          MuonSystem->gParticleEta[MuonSystem->nGenParticles] = GenPart_eta[i];
          MuonSystem->gParticlePhi[MuonSystem->nGenParticles] = GenPart_phi[i];
          MuonSystem->gParticlePt[MuonSystem->nGenParticles] = GenPart_pt[i];
          MuonSystem->gParticleId[MuonSystem->nGenParticles] = GenPart_pdgId[i];
          float decay_x;
          float decay_y;
          float decay_z;
          bool foundDaughter = false;
          for (int j=0; j < nGenPart; j++)
          {
            if(GenPart_genPartIdxMother[j]==i)
            {
              decay_x = GenPart_vx[j];
              decay_y = GenPart_vx[j];
              decay_z = GenPart_vx[j];
              foundDaughter = true;
              break;
            }
          }
          if (!foundDaughter)cout<<"didn't find HV LLP daughter "<<eventNum<<endl;
          MuonSystem->gParticle_decay_vertex_r[MuonSystem->nGenParticles] = sqrt(decay_x*decay_x+decay_y*decay_y);
          MuonSystem->gParticle_decay_vertex_x[MuonSystem->nGenParticles] = decay_x;
          MuonSystem->gParticle_decay_vertex_y[MuonSystem->nGenParticles] = decay_y;
          MuonSystem->gParticle_decay_vertex_z[MuonSystem->nGenParticles] = decay_z;


          TLorentzVector particle = makeTLorentzVectorPtEtaPhiM( GenPart_pt[i], GenPart_eta[i], GenPart_phi[i], GenPart_mass[i] );

          float gParticle_decay_vertex = sqrt(pow(MuonSystem->gParticle_decay_vertex_r[MuonSystem->nGenParticles], 2) + pow(MuonSystem->gParticle_decay_vertex_z[MuonSystem->nGenParticles],2));

          
          MuonSystem->gParticle_ctau[MuonSystem->nGenParticles] = gParticle_decay_vertex/(particle.Beta()*particle.Gamma());
          MuonSystem->gParticle_beta[MuonSystem->nGenParticles] = particle.Beta();
            // if (abs(MuonSystem->gLLP_eta[MuonSystem->nGLLP]) < 2.4
            //   && abs(MuonSystem->gLLP_decay_vertex_z[MuonSystem->nGLLP])<1100 && abs(MuonSystem->gLLP_decay_vertex_z[MuonSystem->nGLLP])>400
            //   && MuonSystem->gLLP_decay_vertex_r[MuonSystem->nGLLP] < 695.5) MuonSystem->gLLP_csc[MuonSystem->nGLLP] = true;
            // if (abs(MuonSystem->gLLP_decay_vertex_z[MuonSystem->nGLLP])< 661.0
            //   && MuonSystem->gLLP_decay_vertex_r[MuonSystem->nGLLP] < 800
            //    && MuonSystem->gLLP_decay_vertex_r[MuonSystem->nGLLP] > 200.0) MuonSystem->gLLP_dt[MuonSystem->nGLLP] = true;
            MuonSystem->nGenParticles++;
          }
          
        if (abs(GenPart_pdgId[i])==15){
          int motherIdx = GenPart_genPartIdxMother[i];
          if (abs(GenPart_pdgId[motherIdx])!=24) continue;
          //cout<<"Tau mother id: "<<GenPart_pdgId[motherIdx]<<endl;
          for (int j=0; j < nGenPart; j++){
            if (GenPart_pdgId[j]==9900012){
              //cout<<"HNL Mother id: "<<GenPart_pdgId[GenPart_genPartIdxMother[j]]<<endl;
              if (GenPart_genPartIdxMother[j]==motherIdx){
              MuonSystem->gParticleEta[MuonSystem->nGenParticles] = GenPart_eta[i];
              MuonSystem->gParticlePhi[MuonSystem->nGenParticles] = GenPart_phi[i];
              MuonSystem->gParticlePt[MuonSystem->nGenParticles] = GenPart_pt[i];
              MuonSystem->gParticleId[MuonSystem->nGenParticles] = GenPart_pdgId[i];
              nPromptTau++;
              MuonSystem->nGenParticles++;
            }
        }
          }		
            
      }

}
      if (nPromptTau!=1){
        cout<<"Warning: found "<<nPromptTau<<" prompt taus in event "<<eventNum<<endl;
        cout<<"Skipping this event"<<endl;
        continue;
      }
      //for HNL
       MuonSystem->nGLLP = 0;
      for(int i = 0;i < nGenPart; i++)
      {
        if (abs(GenPart_pdgId[i])==9900012)
        {
          // cout<<eventNum<<","<<gParticleId[i]<<","<<MuonSystem->nGLLP<<endl;
          MuonSystem->gLLP_eta[MuonSystem->nGLLP] = GenPart_eta[i];
          MuonSystem->gLLP_phi[MuonSystem->nGLLP] = GenPart_phi[i];
          MuonSystem->gLLP_pt[MuonSystem->nGLLP] = GenPart_pt[i];
          // MuonSystem->gLLP_e[MuonSystem->nGenParticles] = GenPart_pdgId[i];
          float decay_x;
          float decay_y;
          float decay_z;
          bool foundDaughter = false;
          for (int j=0; j < nGenPart; j++)
          {
            if(GenPart_genPartIdxMother[j]==i)
            {
              decay_x = GenPart_vx[j];
              decay_y = GenPart_vy[j];
              decay_z = GenPart_vz[j];
              foundDaughter = true;
              break;
            }
          }
          if (!foundDaughter)cout<<"didn't find LLP daughter "<<eventNum<<endl;
          MuonSystem->gLLP_decay_vertex_r[MuonSystem->nGLLP] = sqrt(decay_x*decay_x+decay_y*decay_y);
          MuonSystem->gLLP_decay_vertex_x[MuonSystem->nGLLP] = decay_x;
          MuonSystem->gLLP_decay_vertex_y[MuonSystem->nGLLP] = decay_y;
          MuonSystem->gLLP_decay_vertex_z[MuonSystem->nGLLP] = decay_z;


          TLorentzVector LLP = makeTLorentzVectorPtEtaPhiM( GenPart_pt[i], GenPart_eta[i], GenPart_phi[i], GenPart_mass[i] );

          float gLLP_decay_vertex = sqrt(pow(MuonSystem->gLLP_decay_vertex_r[MuonSystem->nGLLP], 2) + pow(MuonSystem->gLLP_decay_vertex_z[MuonSystem->nGLLP],2));

          
          MuonSystem->gLLP_ctau[MuonSystem->nGLLP] = gLLP_decay_vertex/(LLP.Beta()*LLP.Gamma());
          MuonSystem->gLLP_beta[MuonSystem->nGLLP] = LLP.Beta();
            if (abs(MuonSystem->gLLP_eta[MuonSystem->nGLLP]) < 2.4
              && abs(MuonSystem->gLLP_decay_vertex_z[MuonSystem->nGLLP])<1100 && abs(MuonSystem->gLLP_decay_vertex_z[MuonSystem->nGLLP])>400
              && MuonSystem->gLLP_decay_vertex_r[MuonSystem->nGLLP] < 695.5) MuonSystem->gLLP_csc[MuonSystem->nGLLP] = true;
            if (abs(MuonSystem->gLLP_decay_vertex_z[MuonSystem->nGLLP])< 661.0
              && MuonSystem->gLLP_decay_vertex_r[MuonSystem->nGLLP] < 800
               && MuonSystem->gLLP_decay_vertex_r[MuonSystem->nGLLP] > 200.0) MuonSystem->gLLP_dt[MuonSystem->nGLLP] = true;
            MuonSystem->nGLLP++;
          }
	
      }



    //    for(int i = 0; i < 2;i++)
    //    {
	  //      MuonSystem->gLLP_eta[MuonSystem->nGLLP] = gLLP_eta[i];
    //      MuonSystem->gLLP_phi[MuonSystem->nGLLP] = gLLP_phi[i];
    //      MuonSystem->gLLP_e[MuonSystem->nGLLP] = gLLP_e[i];
    //      MuonSystem->gLLP_pt[MuonSystem->nGLLP] = gLLP_pt[i];

    //      MuonSystem->gLLP_decay_vertex_r[MuonSystem->nGLLP] = sqrt(gLLP_decay_vertex_x[i]*gLLP_decay_vertex_x[i]+gLLP_decay_vertex_y[i]*gLLP_decay_vertex_y[i]);
    //      MuonSystem->gLLP_decay_vertex_x[MuonSystem->nGLLP] = gLLP_decay_vertex_x[i];
    //      MuonSystem->gLLP_decay_vertex_y[MuonSystem->nGLLP] = gLLP_decay_vertex_y[i];
    //      MuonSystem->gLLP_decay_vertex_z[MuonSystem->nGLLP] = gLLP_decay_vertex_z[i];
    //      float beta = gLLP_beta[i];
    //      float gLLP_decay_vertex = sqrt(pow(MuonSystem->gLLP_decay_vertex_r[i], 2) + pow(MuonSystem->gLLP_decay_vertex_z[i],2));
    //      float gamma = 1.0/sqrt(1-beta*beta);
    //      MuonSystem->gLLP_ctau[MuonSystem->nGLLP] = gLLP_decay_vertex/(beta * gamma);
    //      MuonSystem->gLLP_beta[MuonSystem->nGLLP] = gLLP_beta[i];


    //        if (abs(MuonSystem->gLLP_eta[i]) < 2.4
    //          && abs(MuonSystem->gLLP_decay_vertex_z[i])<1100 && abs(MuonSystem->gLLP_decay_vertex_z[i])>400
    //          && MuonSystem->gLLP_decay_vertex_r[i] < 695.5) MuonSystem->gLLP_csc[MuonSystem->nGLLP] = true;
    //        if (abs(MuonSystem->gLLP_decay_vertex_z[i])< 661.0
    //          && MuonSystem->gLLP_decay_vertex_r[i] < 800
    //           && MuonSystem->gLLP_decay_vertex_r[i] > 200.0) MuonSystem->gLLP_dt[MuonSystem->nGLLP] = true;

	  //  if (abs(MuonSystem->gLLP_decay_vertex_z[i])< 1500 && 
		// 	   MuonSystem->gLLP_decay_vertex_r[i] < 1000 &&
		// 	   MuonSystem->gLLP_decay_vertex_r[i]>100 )llp_flag = true;
    //         MuonSystem->nGLLP++;
    //    }


      MuonSystem->npu = Pileup_nTrueInt;
      MuonSystem->pileupWeight = helper->getPileupWeight(Pileup_nTrueInt);
      MuonSystem->pileupWeightUp = helper->getPileupWeightUp(Pileup_nTrueInt) / MuonSystem->pileupWeight;
      MuonSystem->pileupWeightDown = helper->getPileupWeightDown(Pileup_nTrueInt) / MuonSystem->pileupWeight;

      for (unsigned int i = 0; i < 9; i++)
      {
       MuonSystem->LHEScaleWeight[i]= LHEScaleWeight[i];
      }
     



    }//end of isData

      //get NPU
      MuonSystem->npv = nPV;
      MuonSystem->rho = fixedGridRhoFastjetAll;
      MuonSystem->met = metType1Pt;
      MuonSystem->metPhi = metType1Phi;


      MuonSystem->Puppimet = PuppiMET_pt;
      MuonSystem->PuppimetPhi = PuppiMET_phi;

      if (!isData)
      {
        if (MuonSystem->gLLP_csc[0] && MuonSystem->gLLP_csc[1]) accep_csccsc->Fill(1.0, genWeight*MuonSystem->pileupWeight);
        if ((MuonSystem->gLLP_dt[0] && MuonSystem->gLLP_csc[1]) || (MuonSystem->gLLP_dt[1] && MuonSystem->gLLP_csc[0])) accep_cscdt->Fill(1.0, genWeight*MuonSystem->pileupWeight);

      }

      // noise filters

      MuonSystem->Flag_goodVertices = Flag_goodVertices;
      MuonSystem->Flag_globalSuperTightHalo2016Filter = Flag_globalSuperTightHalo2016Filter;
      MuonSystem->Flag_EcalDeadCellTriggerPrimitiveFilter = Flag_EcalDeadCellTriggerPrimitiveFilter;
      MuonSystem->Flag_BadPFMuonFilter = Flag_BadPFMuonFilter;
      MuonSystem->Flag_BadPFMuonDzFilter = Flag_BadPFMuonDzFilter;
      MuonSystem->Flag_hfNoisyHitsFilter = Flag_hfNoisyHitsFilter;
      MuonSystem->Flag_eeBadScFilter = Flag_eeBadScFilter;
      MuonSystem->Flag_all = Flag_eeBadScFilter && Flag_hfNoisyHitsFilter && Flag_BadPFMuonDzFilter && Flag_BadPFMuonFilter && Flag_EcalDeadCellTriggerPrimitiveFilter
                              && Flag_globalSuperTightHalo2016Filter && Flag_goodVertices;
      if (analysisTag == "Summer24") MuonSystem->Flag_ecalBadCalibFilter = Flag_ecalBadCalibFilter;

      // Flag_ecalBadCalibFilter for nanoAOD: https://twiki.cern.ch/twiki/bin/view/CMS/MissingETOptionalFiltersRun2#ECal_BadCalibration_Filter_Flag
      if (analysisTag == "Summer24") MuonSystem->Flag_ecalBadCalibFilter = Flag_ecalBadCalibFilter;
      else{
        MuonSystem->Flag_ecalBadCalibFilter = true;
        if (isData && runNum >= 362433 && runNum<=367144)
        {
          if (PuppiMET_pt > 100)
          {
            for(int i = 0; i < nJets; i++)
            {
              if (jetPt[i]<50) continue;
              if (!(jetEta[i] <= -0.1 && jetEta[i]>=-0.5 && jetPhi[i] <-1.8 && jetPhi[i]> -2.1)) continue;
              if (!(Jet_neEmEF[i] >0.9 || Jet_chEmEF[i]>0.9)) continue;
              if (deltaPhi(PuppiMET_phi, jetPhi[i])<2.9) continue;
              Flag_ecalBadCalibFilter = false;
            }
          }
        }
      }


      // jet veto map, following selections here: https://cms-jerc.web.cern.ch/Recommendations/#jet-veto-maps
      MuonSystem->jetVeto = true;
      for(int i = 0; i < nJets; i++)
      {
        if (jetPt[i] <= 15) continue;
        if (Jet_neEmEF[i] + Jet_chEmEF[i] >= 0.9) continue;
        if (!jetPassIDTight[i]) continue;
        //remove overlaps
        bool overlap = false;
        for(int j = 0; j < nMuons; j++)
        {
          if (!Muon_isPFcand[j])continue;
          if (RazorAnalyzerMerged::deltaR(jetEta[i],jetPhi[i],muonEta[j], muonPhi[j]) < 0.2) overlap = true;
        }
        if (overlap) continue;
        helper->getJetVetoMap(0,1);
        if (helper->getJetVetoMap(jetEta[i],jetPhi[i])>0.0) MuonSystem->jetVeto = false;
        if (analysisTag == "Summer24" && helper->getJetVetoFpixMap(jetEta[i],jetPhi[i])>0.0) MuonSystem->jetVeto = false;
      } 
      
      MuonSystem->HLT_CSCCSC = HLT_CscCluster_Loose;
      // MuonSystem->HLT_CSCDT = HLT_L1CSCShower_DTCluster50;
      MuonSystem->HLT_CscCluster100_PNetTauhPFJet10_Loose = HLT_CscCluster100_PNetTauhPFJet10_Loose;
      MuonSystem->HLT_CscCluster100_Ele5 = HLT_CscCluster100_Ele5;
      MuonSystem->HLT_CscCluster100_Mu5 = HLT_CscCluster100_Mu5;
      MuonSystem->HLT_CscCluster50_Photon30Unseeded = HLT_CscCluster50_Photon30Unseeded;
      MuonSystem->HLT_CscCluster50_Photon20Unseeded = HLT_CscCluster50_Photon20Unseeded;


      //if (!HLT_CscCluster100_PNetTauhPFJet10_Loose) continue; //UNCOMMENT WHEN WE HAVE CLUSTERS
      //cout<<"Pass trigger"<<endl;
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
          if (RazorAnalyzerMerged::deltaR(muonEta[i],muonPhi[i],lep.lepton.Eta(),lep.lepton.Phi()) < 0.3) overlap = true;
        }
        if(overlap) continue;

        leptons tmpMuon;
        tmpMuon.lepton.SetPtEtaPhiM(muonPt[i],muonEta[i], muonPhi[i], MU_MASS);
        tmpMuon.pdgId = 13 * -1 * muonCharge[i];
        tmpMuon.dZ = muon_dZ[i];
        tmpMuon.passId = muonIsTight[i];

        tmpMuon.passLooseIso = muon_pfRelIso04_all[i]<0.25;
        tmpMuon.passTightIso = muon_pfRelIso04_all[i]<0.15;
        tmpMuon.passVTightIso = muon_pfRelIso04_all[i]<0.10;
        tmpMuon.passVVTightIso = muon_pfRelIso04_all[i]<0.05;

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
          if (RazorAnalyzerMerged::deltaR(eleEta[i],elePhi[i],lep.lepton.Eta(),lep.lepton.Phi()) < 0.3) overlap = true;
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


      //********************************
      //Tau
      //********************************
      //cout<<"nTau: "<<nTau<<endl;
      for(int i = 0; i < nTau;  i++){
        // if(Tau_pt[i] < 20) continue;
        if(abs(Tau_eta[i]) > 2.3) continue;
        if(abs(Tau_dz[i]) > 0.2) continue;
        if(Tau_decayMode[i]==5 || Tau_decayMode[i]==6) continue;
        
        
        if ( Tau_idDeepTau2018v2p5VSjet[i] >= 1  && Tau_idDeepTau2018v2p5VSe[i] >= 2 && Tau_idDeepTau2018v2p5VSmu[i] == 4 ) MuonSystem->tauIsVVVLoose[MuonSystem->nTaus] = true; 
        if ( Tau_idDeepTau2018v2p5VSjet[i] >= 2  && Tau_idDeepTau2018v2p5VSe[i] >= 2 && Tau_idDeepTau2018v2p5VSmu[i] == 4 ) MuonSystem->tauIsVVLoose[MuonSystem->nTaus]  = true;
        if ( Tau_idDeepTau2018v2p5VSjet[i] >= 3  && Tau_idDeepTau2018v2p5VSe[i] >= 2 && Tau_idDeepTau2018v2p5VSmu[i] == 4 ) MuonSystem->tauIsVLoose[MuonSystem->nTaus]   = true;
        if ( Tau_idDeepTau2018v2p5VSjet[i] >= 4  && Tau_idDeepTau2018v2p5VSe[i] >= 2 && Tau_idDeepTau2018v2p5VSmu[i] == 4 ) MuonSystem->tauIsLoose[MuonSystem->nTaus]    = true;
        if ( Tau_idDeepTau2018v2p5VSjet[i] >= 5 && Tau_idDeepTau2018v2p5VSe[i] >= 2 && Tau_idDeepTau2018v2p5VSmu[i] == 4 ) MuonSystem->tauIsMedium[MuonSystem->nTaus]   = true;
        if ( Tau_idDeepTau2018v2p5VSjet[i] >= 6 && Tau_idDeepTau2018v2p5VSe[i] >= 2 && Tau_idDeepTau2018v2p5VSmu[i] == 4 ) MuonSystem->tauIsTight[MuonSystem->nTaus]    = true;
        if ( Tau_idDeepTau2018v2p5VSjet[i] >= 7 && Tau_idDeepTau2018v2p5VSe[i] >= 2 && Tau_idDeepTau2018v2p5VSmu[i] == 4 ) MuonSystem->tauIsVTight[MuonSystem->nTaus]   = true;
        if ( Tau_idDeepTau2018v2p5VSjet[i] >=8 && Tau_idDeepTau2018v2p5VSe[i] >= 2 && Tau_idDeepTau2018v2p5VSmu[i] == 4 ) MuonSystem->tauIsVVTight[MuonSystem->nTaus]  = true;

 
        MuonSystem->tauM[MuonSystem->nTaus]            = Tau_mass[i];
        MuonSystem->tauPt[MuonSystem->nTaus]           = Tau_pt[i];
        MuonSystem->tauEta[MuonSystem->nTaus]          = Tau_eta[i];
        MuonSystem->tauPhi[MuonSystem->nTaus]          = Tau_phi[i];
        MuonSystem->tauDz[MuonSystem->nTaus]           = Tau_dz[i];
        MuonSystem->tauGenPartFlav[MuonSystem->nTaus]  = Tau_genPartFlav[i];
        MuonSystem->tauDecayMode[MuonSystem->nTaus]    = Tau_decayMode[i];

        if (!isData){
          for (int j = 0; j < nGenPart; j++){
            if (abs(MuonSystem->gParticleId[j])!=15) continue;
            MuonSystem->deltaR_GenTauRecoTau[MuonSystem->nTaus] = (float) RazorAnalyzerMerged::deltaR(Tau_eta[i], Tau_phi[i], MuonSystem->gParticleEta[j], MuonSystem->gParticlePhi[j]);
            }
          }
        //load discriminator values
        MuonSystem->tauIdDeepTau2018v2p5VSjet[MuonSystem->nTaus] = Tau_idDeepTau2018v2p5VSjet[i];
        MuonSystem->tauIdDeepTau2018v2p5VSe[MuonSystem->nTaus] = Tau_idDeepTau2018v2p5VSe[i];
        MuonSystem->tauIdDeepTau2018v2p5VSmu[MuonSystem->nTaus] = Tau_idDeepTau2018v2p5VSmu[i];
        MuonSystem->nTaus++;
        // TLorentzVector thisTau; thisTau.SetPtEtaPhiM(Tau_pt[i], Tau_eta[i], Tau_phi[i], Tau_mass[i]);
      }
      if (MuonSystem->nTaus == 0)continue;
      //cout<<"Found Event with Tau: "<<eventNum<<", nTaus: "<<MuonSystem->nTaus<<endl;
    //-----------------------------------------------
    //Select Jets
    //-----------------------------------------------

    std::vector<jets> Jets;

    for(int i = 0; i < nJets; i++)
    {
      if (fabs(jetEta[i]) >= 3.0)continue;
      if( jetPt[i] < 20 ) continue;
      if (!jetPassIDTight[i] && !jetPassIDTightLepVeto[i]) continue;
      //------------------------------------------------------------
      //exclude selected muons and electrons from the jet collection
      //------------------------------------------------------------
      double deltaR = -1;
      for(auto& lep : Leptons){
        double thisDR = RazorAnalyzerMerged::deltaR(jetEta[i],jetPhi[i],lep.lepton.Eta(),lep.lepton.Phi());
        if(deltaR < 0 || thisDR < deltaR) deltaR = thisDR;
      }
      if(deltaR > 0 && deltaR < 0.4) continue; //jet matches a selected lepton

      TLorentzVector thisJet = makeTLorentzVector( jetPt[i], jetEta[i], jetPhi[i], jetE[i] );
      
      jets tmpJet;
      tmpJet.jet    = thisJet;
      tmpJet.passId = jetPassIDTightLepVeto[i];

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

      TLorentzVector met = makeTLorentzVectorPtEtaPhiM(MuonSystem->Puppimet, 0, MuonSystem->PuppimetPhi, 0);



      MuonSystem->nDtRechits  = 0;

      int nDtRechitsChamberMinus12 = 0;
      int nDtRechitsChamberMinus11 = 0;
      int nDtRechitsChamber10 = 0;
      int nDtRechitsChamberPlus11 = 0;
      int nDtRechitsChamberPlus12 = 0;
      int nDtRechitsChamberMinus22 = 0;
      int nDtRechitsChamberMinus21 = 0;
      int nDtRechitsChamber20 = 0;
      int nDtRechitsChamberPlus21 = 0;
      int nDtRechitsChamberPlus22 = 0;
      int nDtRechitsChamberMinus32 = 0;
      int nDtRechitsChamberMinus31 = 0;
      int nDtRechitsChamber30 = 0;
      int nDtRechitsChamberPlus31 = 0;
      int nDtRechitsChamberPlus32 = 0;
      int nDtRechitsChamberMinus42 = 0;
      int nDtRechitsChamberMinus41 = 0;
      int nDtRechitsChamber40 = 0;
      int nDtRechitsChamberPlus41 = 0;
      int nDtRechitsChamberPlus42 = 0;

      for (int i = 0; i < nDtRechits; i++) {

        if (dtRechitStation[i] == 1 && dtRechitWheel[i] == -2) nDtRechitsChamberMinus12++;
        if (dtRechitStation[i] == 1 && dtRechitWheel[i] == -1) nDtRechitsChamberMinus11++;
        if (dtRechitStation[i] == 1 && dtRechitWheel[i] == 0) nDtRechitsChamber10++;
        if (dtRechitStation[i] == 1 && dtRechitWheel[i] == 1) nDtRechitsChamberPlus11++;
        if (dtRechitStation[i] == 1 && dtRechitWheel[i] == 2) nDtRechitsChamberPlus12++;
        if (dtRechitStation[i] == 2 && dtRechitWheel[i] == -2) nDtRechitsChamberMinus22++;
        if (dtRechitStation[i] == 2 && dtRechitWheel[i] == -1) nDtRechitsChamberMinus21++;
        if (dtRechitStation[i] == 2 && dtRechitWheel[i] == 0) nDtRechitsChamber20++;
        if (dtRechitStation[i] == 2 && dtRechitWheel[i] == 1) nDtRechitsChamberPlus21++;
        if (dtRechitStation[i] == 2 && dtRechitWheel[i] == 2) nDtRechitsChamberPlus22++;
        if (dtRechitStation[i] == 3 && dtRechitWheel[i] == -2) nDtRechitsChamberMinus32++;
        if (dtRechitStation[i] == 3 && dtRechitWheel[i] == -1) nDtRechitsChamberMinus31++;
        if (dtRechitStation[i] == 3 && dtRechitWheel[i] == 0) nDtRechitsChamber30++;
        if (dtRechitStation[i] == 3 && dtRechitWheel[i] == 1) nDtRechitsChamberPlus31++;
        if (dtRechitStation[i] == 3 && dtRechitWheel[i] == 2) nDtRechitsChamberPlus32++;
        if (dtRechitStation[i] == 4 && dtRechitWheel[i] == -2) nDtRechitsChamberMinus42++;
        if (dtRechitStation[i] == 4 && dtRechitWheel[i] == -1) nDtRechitsChamberMinus41++;
        if (dtRechitStation[i] == 4 && dtRechitWheel[i] == 0) nDtRechitsChamber40++;
        if (dtRechitStation[i] == 4 && dtRechitWheel[i] == 1) nDtRechitsChamberPlus41++;
        if (dtRechitStation[i] == 4 && dtRechitWheel[i] == 2) nDtRechitsChamberPlus42++;




      }


      if ( nDtRechitsChamberMinus12 > 50) MuonSystem->nDtRings++;
      if ( nDtRechitsChamberMinus11 > 50) MuonSystem->nDtRings++;
      if ( nDtRechitsChamber10 > 50) MuonSystem->nDtRings++;
      if ( nDtRechitsChamberPlus11 > 50) MuonSystem->nDtRings++;
      if ( nDtRechitsChamberPlus12 > 50) MuonSystem->nDtRings++;
      if ( nDtRechitsChamberMinus22 > 50) MuonSystem->nDtRings++;
      if ( nDtRechitsChamberMinus21 > 50) MuonSystem->nDtRings++;
      if ( nDtRechitsChamber20 > 50) MuonSystem->nDtRings++;
      if ( nDtRechitsChamberPlus21 > 50) MuonSystem->nDtRings++;
      if ( nDtRechitsChamberPlus22 > 50) MuonSystem->nDtRings++;
      if ( nDtRechitsChamberMinus32 > 50) MuonSystem->nDtRings++;
      if ( nDtRechitsChamberMinus31 > 50) MuonSystem->nDtRings++;
      if ( nDtRechitsChamber30 > 50) MuonSystem->nDtRings++;
      if ( nDtRechitsChamberPlus31 > 50) MuonSystem->nDtRings++;
      if ( nDtRechitsChamberPlus32 > 50) MuonSystem->nDtRings++;
      if ( nDtRechitsChamberMinus42 > 50) MuonSystem->nDtRings++;
      if ( nDtRechitsChamberMinus41 > 50) MuonSystem->nDtRings++;
      if ( nDtRechitsChamber40 > 50) MuonSystem->nDtRings++;
      if ( nDtRechitsChamberPlus41 > 50) MuonSystem->nDtRings++;
      if ( nDtRechitsChamberPlus42 > 50) MuonSystem->nDtRings++;



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
        // cout<<cscRechitsChamber[i]<<endl;
        //pick out the right bits for chamber
        // int chamber = ((cscRechitsDetId[i] >> 3) & 077); //https://github.com/cms-sw/cmssw/blob/master/DataFormats/MuonDetId/interface/CSCDetId.h#L147

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
      float epsilon = 0.4; //cluster radius parameter
      CACluster ds(min_point, epsilon, points);
      ds.run();
      ds.clusterProperties();
      //ds.merge_clusters();
      //ds.clusterProperties();
      ds.sort_clusters();
      MuonSystem->nCscRechitClusters = 0;
      MuonSystem->nCscRechitClusters_nocut = 0;
      for ( auto &tmp : ds.clusters  ) {
	  MuonSystem->nCscRechitClusters_nocut++;
	  // if (tmp.nCscRechitsChamberPlus11 + tmp.nCscRechitsChamberPlus12 + tmp.nCscRechitsChamberMinus11 + tmp.nCscRechitsChamberMinus12 > 0)continue;
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


          // MuonSystem->cscRechitClusterHMTEffPlus11[MuonSystem->nCscRechitClusters] = helper->getHMTTriggerEff(11, tmp.nCscRechitsChamberPlus11);
          // MuonSystem->cscRechitClusterHMTEffPlus12[MuonSystem->nCscRechitClusters] = helper->getHMTTriggerEff(12, tmp.nCscRechitsChamberPlus12);
          // MuonSystem->cscRechitClusterHMTEffPlus13[MuonSystem->nCscRechitClusters] = helper->getHMTTriggerEff(13, tmp.nCscRechitsChamberPlus13);
          // MuonSystem->cscRechitClusterHMTEffPlus21[MuonSystem->nCscRechitClusters] = helper->getHMTTriggerEff(21, tmp.nCscRechitsChamberPlus21);
          // MuonSystem->cscRechitClusterHMTEffPlus22[MuonSystem->nCscRechitClusters] = helper->getHMTTriggerEff(22, tmp.nCscRechitsChamberPlus22);
          // MuonSystem->cscRechitClusterHMTEffPlus31[MuonSystem->nCscRechitClusters] = helper->getHMTTriggerEff(31, tmp.nCscRechitsChamberPlus31);
          // MuonSystem->cscRechitClusterHMTEffPlus32[MuonSystem->nCscRechitClusters] = helper->getHMTTriggerEff(32, tmp.nCscRechitsChamberPlus32);
          // MuonSystem->cscRechitClusterHMTEffPlus41[MuonSystem->nCscRechitClusters] = helper->getHMTTriggerEff(41, tmp.nCscRechitsChamberPlus41);
          // MuonSystem->cscRechitClusterHMTEffPlus42[MuonSystem->nCscRechitClusters] = helper->getHMTTriggerEff(42, tmp.nCscRechitsChamberPlus42);
          // MuonSystem->cscRechitClusterHMTEffMinus11[MuonSystem->nCscRechitClusters] = helper->getHMTTriggerEff(11, tmp.nCscRechitsChamberMinus11);
          // MuonSystem->cscRechitClusterHMTEffMinus12[MuonSystem->nCscRechitClusters] = helper->getHMTTriggerEff(12, tmp.nCscRechitsChamberMinus12);
          // MuonSystem->cscRechitClusterHMTEffMinus13[MuonSystem->nCscRechitClusters] = helper->getHMTTriggerEff(13, tmp.nCscRechitsChamberMinus13);
          // MuonSystem->cscRechitClusterHMTEffMinus21[MuonSystem->nCscRechitClusters] = helper->getHMTTriggerEff(21, tmp.nCscRechitsChamberMinus21);
          // MuonSystem->cscRechitClusterHMTEffMinus22[MuonSystem->nCscRechitClusters] = helper->getHMTTriggerEff(22, tmp.nCscRechitsChamberMinus22);
          // MuonSystem->cscRechitClusterHMTEffMinus31[MuonSystem->nCscRechitClusters] = helper->getHMTTriggerEff(31, tmp.nCscRechitsChamberMinus31);
          // MuonSystem->cscRechitClusterHMTEffMinus32[MuonSystem->nCscRechitClusters] = helper->getHMTTriggerEff(32, tmp.nCscRechitsChamberMinus32);
          // MuonSystem->cscRechitClusterHMTEffMinus41[MuonSystem->nCscRechitClusters] = helper->getHMTTriggerEff(41, tmp.nCscRechitsChamberMinus41);
          // MuonSystem->cscRechitClusterHMTEffMinus42[MuonSystem->nCscRechitClusters] = helper->getHMTTriggerEff(42, tmp.nCscRechitsChamberPlus42);


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
            helper->getHMTTriggerEff(42, tmp.nCscRechitsChamberPlus42),
          };

          MuonSystem->cscRechitClusterHMTEfficiency[MuonSystem->nCscRechitClusters] = 1.0;
          for(const float &eff : efficiencies){
              MuonSystem->cscRechitClusterHMTEfficiency[MuonSystem->nCscRechitClusters] *= (1-eff);
          }
          MuonSystem->cscRechitClusterHMTEfficiency[MuonSystem->nCscRechitClusters] = 1-MuonSystem->cscRechitClusterHMTEfficiency[MuonSystem->nCscRechitClusters];




          // MuonSystem->cscRechitClusterHMTEffMinus42[MuonSystem->nCscRechitClusters] = getHMTTriggerEff(11, tmp.nCscRechitsChamberPlus11);

          // cout<<tmp.nCscRechitsChamberPlus42<<", "<<helper->getHMTTriggerEff(42, tmp.nCscRechitsChamberPlus42)<<endl;;





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
            if (RazorAnalyzerMerged::deltaR(jetEta[i], jetPhi[i], MuonSystem->cscRechitClusterEta[MuonSystem->nCscRechitClusters],MuonSystem->cscRechitClusterPhi[MuonSystem->nCscRechitClusters]) < 0.4 && jetPt[i] > MuonSystem->cscRechitClusterJetVetoPt[MuonSystem->nCscRechitClusters] ) {
              MuonSystem->cscRechitClusterJetVetoPt[MuonSystem->nCscRechitClusters]  = jetPt[i];
              MuonSystem->cscRechitClusterJetVetoE[MuonSystem->nCscRechitClusters]  = jetE[i];
              MuonSystem->cscRechitClusterJetVetoTightId[MuonSystem->nCscRechitClusters]  = jetPassIDTight[i];

            }

          }
          float min_deltaR = 15.;
          int index = 999;


          for(int i = 0; i < nMuons; i++)
          {
            if (fabs(muonEta[i]>3.0)) continue;
            // float muonIso = (muon_chargedIso[i] + fmax(0.0,  muon_photonIso[i] + muon_neutralHadIso[i] - 0.5*muon_pileupIso[i])) / muonPt[i];
            if (RazorAnalyzerMerged::deltaR(muonEta[i], muonPhi[i], MuonSystem->cscRechitClusterEta[MuonSystem->nCscRechitClusters],MuonSystem->cscRechitClusterPhi[MuonSystem->nCscRechitClusters]) < 0.4 && muonPt[i] > MuonSystem->cscRechitClusterMuonVetoPt[MuonSystem->nCscRechitClusters] ) {
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

              double current_delta_r = RazorAnalyzerMerged::deltaR(MuonSystem->cscRechitClusterEta[MuonSystem->nCscRechitClusters], MuonSystem->cscRechitClusterPhi[MuonSystem->nCscRechitClusters], MuonSystem->gLLP_eta[j], MuonSystem->gLLP_phi[j]);
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
            if (RazorAnalyzerMerged::deltaR(dtSegEta[i], dtSegPhi[i], MuonSystem->cscRechitClusterEta[MuonSystem->nCscRechitClusters],MuonSystem->cscRechitClusterPhi[MuonSystem->nCscRechitClusters]) < 0.4 )
            {
              MuonSystem->cscRechitCluster_match_dtSeg_0p4[MuonSystem->nCscRechitClusters] ++;
              if (dtSegStation[i] == 1) MuonSystem->cscRechitCluster_match_MB1Seg_0p4[MuonSystem->nCscRechitClusters] ++;
            }
          }

          //match to RPC hits in RE1/2
          for (int i = 0; i < nRpc; i++) {
            float rpcR = sqrt(rpcX[i]*rpcX[i] + rpcY[i]*rpcY[i]);
            if (RazorAnalyzerMerged::deltaR(rpcEta[i], rpcPhi[i], MuonSystem->cscRechitClusterEta[MuonSystem->nCscRechitClusters],MuonSystem->cscRechitClusterPhi[MuonSystem->nCscRechitClusters]) < 0.4 )
            {
              if (rpcR < 461.0 && rpcR > 275 && abs(rpcZ[i]) > 663 && abs(rpcZ[i]) < 730) MuonSystem->cscRechitCluster_match_RE12_0p4[MuonSystem->nCscRechitClusters] ++;
              if (rpcR < 470 && rpcR > 380 && abs(rpcZ[i]) < 661) MuonSystem->cscRechitCluster_match_RB1_0p4[MuonSystem->nCscRechitClusters] ++;
            }
          }

          MuonSystem->cscRechitClusterMet_dPhi[MuonSystem->nCscRechitClusters] =  RazorAnalyzerMerged::deltaPhi(MuonSystem->cscRechitClusterPhi[MuonSystem->nCscRechitClusters],MuonSystem->metPhi);
          MuonSystem->cscRechitClusterPuppiMet_dPhi[MuonSystem->nCscRechitClusters] =  RazorAnalyzerMerged::deltaPhi(MuonSystem->cscRechitClusterPhi[MuonSystem->nCscRechitClusters],MuonSystem->PuppimetPhi);

          MuonSystem->nCscRechitClusters++;
      }



      //if (MuonSystem->nCscRechitClusters_nocut == 0) continue; //commented out for now, checking tau kinematics with prompt HNL


      MuonSystem->tree_->Fill();


    }

      if (!isData)
      {
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
