/* #region: includes and usings */
#include "main.h"
#include "MuonSystemCuttingPhase.h"
#include "MuonSystemFillingPhase.h"
#include "MuonSystemPhaseContext.h"
#include "MuonSystemPhaseTypes.h"
#include "MuonSystemSignalScanManager.h"
#include "MuonSystemSynthesisPhase.h"
#include "RazorHelper.h"
#include "TreeMuonSystemCombination.h"

#include "TLorentzVector.h"
#include <chrono>
#include <cmath>
#include <iostream>
#include <memory>

//using namespace fastjet;
using namespace std::chrono;
using namespace std;
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
  bool signalScan = int(options / 10) == 1;
  int option = options % 10;
  if (analysisTag.empty()) {
    analysisTag = "Summer24";
    std::cout << "[WARNING]: empty analysisTag received, defaulting to " << analysisTag << '\n';
  }
  /* #endregion */

  /* #region: run-mode logging*/
  std::cout << "[INFO]: running on " << (isData ? "data" : "MC")
            << " with option: " << option << '\n';
  std::cout << "[INFO]: running "
            << (signalScan ? "with Signal scan" : "without Signal scan ")
            << (signalScan ? "" : std::to_string(option))
            << '\n';
  /* #endregion */

  /* #region: setting up output file/tree */
  std::unique_ptr<TFile> outFile;
  if (isData || !signalScan) {
    string outfilename = outputfilename.empty() ? "MuonSystem_Tree.root" : outputfilename;
    outFile = std::make_unique<TFile>(outfilename.c_str(), "RECREATE");
    if (!outFile || outFile->IsZombie()) {
      cerr << "[ERROR]: failed to create output file: " << outfilename << endl;
      return;
    }
    outFile->cd();
  }

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
  MuonSystemSignalScanManager signalScanManager(
      !isData && signalScan,
      outputfilename,
      MuonSystem);
  // [1] Event loop
  /* #region: [1] event loop */
  cout << "[INFO]: Loop Starting" << endl;
  auto lastReport = steady_clock::now(); // mr. timer
  for (Long64_t jentry = 0; jentry < nEntries; jentry++) {

    // progress logging
    if (jentry % 1000 == 0) {
      const auto now = steady_clock::now();
      const double dt = duration<double>(now - lastReport).count();
      cout << "Processing entry " << jentry << " | dt=" << dt << " s\n";
      lastReport = now;
    }

    /* #region: [2] read event entry and reset output event variables */
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0)
      break;
    fChain->GetEntry(jentry);

    //fill normalization histogram
    MuonSystem->InitVariables();
    /* #endregion */

    /* #region: [3] build local aliases and per-event helper inputs */
    EventTransientState eventState;
    eventState.synth = buildEventSynthesis(
        helper.get(),
        analysisTag,
        nElectron,
        Electron_cutBased,
        nJet,
        Jet_eta,
        Jet_pt,
        Jet_mass,
        Jet_neHEF,
        Jet_neEmEF,
        Jet_chEmEF,
        Jet_muEF,
        Jet_chHEF,
        Jet_chMultiplicity,
        Jet_neMultiplicity,
        Jet_jetId,
        nMuon,
        Muon_eta,
        Muon_pt);

    /* #endregion */

    // [4] Signal-scan routing (split output by model point)
    /* #region: [4] signal-scan routing (split output by model point) */
    SignalEventState signalState;
    if (!isData && signalScan) {
      signalState = signalScanManager.beginEvent(*this, Generator_weight);
      if (!signalState.hasKey) {
        std::cerr << "[ERROR]: signalScan requested but event has no parsed signal-point key. Exiting." << std::endl;
        return;
      }
    }
    /* #endregion */

    // [5] Event metadata and per-event nominal weight
    /* #region: [5] event metadata and per-event nominal weight */
    //event info
    if (isData) { //? can this be better organized?
      NEvents->Fill(1);
      MuonSystem->weight = 1;
    } else {
      MuonSystem->weight = Generator_weight;
      NEvents->Fill(1, Generator_weight);
    }
    MuonSystem->run = run; //? is there any particular reason these need to be here?
    MuonSystem->luminosityBlock = luminosityBlock;
    MuonSystem->event = event;

    if (isData && run < 360019) continue; //? should this be moved somewhere higher to save compute?
    /* #endregion */

    // [6] MC-only truth matching inputs and MC weights
    /* #region: [6] MC-only truth matching inputs and MC weights */
    if (!isData) {
      //for DS model
      /* #region: [6a] DS model truth fill */
      MuonSystem->nGenParticles = 0;
      for (int i = 0; i < nGenPart; i++) {
        if (abs(GenPart_pdgId[i]) >= 999999) {
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
          if(foundDaughter){
            MuonSystem->gParticle_decay_vertex_r[MuonSystem->nGenParticles] = sqrt(decay_x * decay_x + decay_y * decay_y);
            MuonSystem->gParticle_decay_vertex_x[MuonSystem->nGenParticles] = decay_x;
            MuonSystem->gParticle_decay_vertex_y[MuonSystem->nGenParticles] = decay_y;
            MuonSystem->gParticle_decay_vertex_z[MuonSystem->nGenParticles] = decay_z;
          } else {
            cout << "didn't find HV LLP daughter " << event << endl;
          }

          TLorentzVector particle = makeTLorentzVectorPtEtaPhiM(GenPart_pt[i], GenPart_eta[i], GenPart_phi[i], GenPart_mass[i]);

          float gParticle_decay_vertex = sqrt(pow(MuonSystem->gParticle_decay_vertex_r[MuonSystem->nGenParticles], 2) + pow(MuonSystem->gParticle_decay_vertex_z[MuonSystem->nGenParticles], 2));

          MuonSystem->gParticle_ctau[MuonSystem->nGenParticles] = gParticle_decay_vertex / (particle.Beta() * particle.Gamma());
          MuonSystem->gParticle_beta[MuonSystem->nGenParticles] = particle.Beta();
          MuonSystem->nGenParticles++;
        }
      }
      /* #endregion */

      //for Twin higgs model
      /* #region: [6b] Twin Higgs truth fill */
      MuonSystem->nGLLP = 0;
      for (int i = 0; i < nGenPart; i++) {
        if (abs(GenPart_pdgId[i]) == 9000006) {
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
          if (foundDaughter){
            MuonSystem->gLLP_decay_vertex_r[MuonSystem->nGLLP] = sqrt(decay_x * decay_x + decay_y * decay_y);
            MuonSystem->gLLP_decay_vertex_x[MuonSystem->nGLLP] = decay_x;
            MuonSystem->gLLP_decay_vertex_y[MuonSystem->nGLLP] = decay_y;
            MuonSystem->gLLP_decay_vertex_z[MuonSystem->nGLLP] = decay_z;
          } else {
            cout << "didn't find LLP daughter " << event << endl;
          }

          TLorentzVector LLP = makeTLorentzVectorPtEtaPhiM(GenPart_pt[i], GenPart_eta[i], GenPart_phi[i], GenPart_mass[i]);
          MuonSystem->gLLP_e[MuonSystem->nGLLP] = LLP.E();

          float gLLP_decay_vertex = sqrt(pow(MuonSystem->gLLP_decay_vertex_r[MuonSystem->nGLLP], 2) + pow(MuonSystem->gLLP_decay_vertex_z[MuonSystem->nGLLP], 2));

          MuonSystem->gLLP_ctau[MuonSystem->nGLLP] = gLLP_decay_vertex / (LLP.Beta() * LLP.Gamma());
          MuonSystem->gLLP_beta[MuonSystem->nGLLP] = LLP.Beta();
          MuonSystem->gLLP_csc[MuonSystem->nGLLP] =
              abs(MuonSystem->gLLP_eta[MuonSystem->nGLLP]) < 2.4 &&
              abs(MuonSystem->gLLP_decay_vertex_z[MuonSystem->nGLLP]) < 1100 &&
              abs(MuonSystem->gLLP_decay_vertex_z[MuonSystem->nGLLP]) > 400 &&
              MuonSystem->gLLP_decay_vertex_r[MuonSystem->nGLLP] < 695.5;
          MuonSystem->gLLP_dt[MuonSystem->nGLLP] =
              abs(MuonSystem->gLLP_decay_vertex_z[MuonSystem->nGLLP]) < 661.0 &&
              MuonSystem->gLLP_decay_vertex_r[MuonSystem->nGLLP] < 800 &&
              MuonSystem->gLLP_decay_vertex_r[MuonSystem->nGLLP] > 200.0;
          MuonSystem->nGLLP++;
        }
      }
      /* #endregion */

      MuonSystem->Pileup_nTrueInt = Pileup_nTrueInt;
      MuonSystem->pileupWeight = helper->getPileupWeight(Pileup_nTrueInt);
      MuonSystem->pileupWeightUp = helper->getPileupWeightUp(Pileup_nTrueInt) / MuonSystem->pileupWeight;
      MuonSystem->pileupWeightDown = helper->getPileupWeightDown(Pileup_nTrueInt) / MuonSystem->pileupWeight;

      for (unsigned int i = 0; i < 9; i++) {
        MuonSystem->LHEScaleWeight[i] = LHEScaleWeight[i];
      }

    }
    /* #endregion */

    // [7] Event-level observables and acceptance bookkeeping
    /* #region: [7] event-level observables and acceptance bookkeeping */
    //get NPU
    MuonSystem->PV_npvs = PV_npvs;
    MuonSystem->Rho_fixedGridRhoFastjetAll = Rho_fixedGridRhoFastjetAll;
    MuonSystem->PFMET_pt = PFMET_pt;
    MuonSystem->PFMET_phi = PFMET_phi;

    MuonSystem->PuppiMET_pt = PuppiMET_pt;
    MuonSystem->PuppiMET_phi = PuppiMET_phi;

    if (signalScan && !isData && signalState.hasKey)
      signalScanManager.fillTotal(signalState, Generator_weight * MuonSystem->pileupWeight);
    if (signalScan && !isData && signalState.hasKey) {
      signalScanManager.fillAccep(signalState, Generator_weight * MuonSystem->pileupWeight);
    } else if (!isData) {
      if (MuonSystem->gLLP_csc[0] && MuonSystem->gLLP_csc[1])
        accep_csccsc->Fill(1.0, Generator_weight * MuonSystem->pileupWeight);
      if ((MuonSystem->gLLP_dt[0] && MuonSystem->gLLP_csc[1]) || (MuonSystem->gLLP_dt[1] && MuonSystem->gLLP_csc[0]))
        accep_cscdt->Fill(1.0, Generator_weight * MuonSystem->pileupWeight);
    }
    /* #endregion */

    // [8] Event filters, trigger bits, and jet veto maps
    /* #region: [8] event filters, trigger bits, and jet veto maps */
    const EventCutState cuts = buildEventCutState(
        *this,
        helper.get(),
        analysisTag,
        isData,
        run,
        PuppiMET_pt,
        PuppiMET_phi,
        nJet,
        Jet_pt,
        Jet_eta,
        Jet_phi,
        Jet_neEmEF,
        Jet_chEmEF,
        nMuon,
        Muon_isPFcand,
        Muon_eta,
        Muon_phi,
        eventState.synth,
        Flag_goodVertices,
        Flag_globalSuperTightHalo2016Filter,
        Flag_EcalDeadCellTriggerPrimitiveFilter,
        Flag_BadPFMuonFilter,
        Flag_BadPFMuonDzFilter,
        Flag_hfNoisyHitsFilter,
        Flag_eeBadScFilter,
        Flag_ecalBadCalibFilter,
        HLT_CscCluster_Loose,
        HLT_L1CSCShower_DTCluster50);
    writeEventCutState(MuonSystem, cuts);
    /* #endregion */

    // [9] Lepton object selection and lepton branch fill
    /* #region: [9] lepton object selection and lepton branch fill */
    //*************************************************************************
    //Start Object Selection
    //*************************************************************************
    eventState.leptons = buildSelectedLeptons(
        *this,
        eventState.synth,
        nMuon,
        Muon_looseId,
        Muon_pt,
        Muon_eta,
        Muon_phi,
        Muon_charge,
        Muon_dz,
        Muon_tightId,
        Muon_pfRelIso04_all,
        nElectron,
        Electron_pt,
        Electron_eta,
        Electron_phi,
        Electron_charge,
        Electron_dz);
    fillLeptonBranches(MuonSystem, eventState.leptons);
    /* #endregion */

    // [10] Jet selection, JES propagation, and jet branch fill
    /* #region: [10] jet selection, JES propagation, and jet branch fill */
    eventState.jetStage = buildJetStageResult(
        *this,
        helper.get(),
        run,
        nJet,
        Jet_eta,
        Jet_phi,
        Jet_pt,
        eventState.synth,
        eventState.leptons);

    fillJetBranches(MuonSystem, eventState.jetStage.selectedJets);
    fillPuppiMetJesFromShift(*this, MuonSystem, eventState.jetStage);
    /* #endregion */

    // [11] Rechit-level pass-through fill (CSC/DT/RPC + optional aliases)
    /* #region: [11] rechit-level pass-through fill (CSC/DT/RPC + optional aliases) */
    fillRawRechits(
        MuonSystem,
        ncscRechits,
        ndtRecHits,
        nrpcRecHits,
        cscRechits_Quality,
        cscRechits_Chamber,
        cscRechits_Station,
        cscRechits_Eta,
        cscRechits_Phi,
        cscRechits_X,
        cscRechits_Y,
        cscRechits_Z,
        cscRechits_Tpeak,
        cscRechits_Twire,
        cscRechits_IChamber,
        cscRechits_NStrips,
        cscRechits_WGroupsBX,
        cscRechits_HitWire,
        cscRechits_NWireGroups,
        cscRechits_E,
        dtRecHits_Layer,
        dtRecHits_SuperLayer,
        dtRecHits_Station,
        dtRecHits_Wheel,
        dtRecHits_Eta,
        dtRecHits_Phi,
        dtRecHits_X,
        dtRecHits_Y,
        dtRecHits_Z,
        dtRecHits_Sector,
        rpcRecHits_Bx,
        rpcRecHits_Region,
        rpcRecHits_Ring,
        rpcRecHits_Layer,
        rpcRecHits_Station,
        rpcRecHits_Sector,
        rpcRecHits_X,
        rpcRecHits_Y,
        rpcRecHits_Z,
        rpcRecHits_Phi,
        rpcRecHits_Eta,
        rpcRecHits_Time,
        rpcRecHits_TimeError);
    /* #endregion */

    // [12] Ring occupancy summaries from raw rechits
    /* #region: [12] ring occupancy summaries from raw rechits */
    MuonSystem->nDtRings = countDtRingsFromRecHits(
        ndtRecHits,
        dtRecHits_Station,
        dtRecHits_Wheel);
    /* #endregion */

    // [13] CSC rechit clustering, feature fill, and matching/veto variables
    /* #region: [13] CSC rechit clustering, feature fill, and matching/veto variables */
    processCscClusterStage(
        *this,
        helper.get(),
        MuonSystem,
        eventState.synth,
        isData,
        run);
    /* #endregion */

    // [14] DT/RPC clustering, DT feature fill, and matching/veto variables
    /* #region: [14] DT/RPC clustering, DT feature fill, and matching/veto variables */
    processDtRpcClusterStage(
        *this,
        helper.get(),
        MuonSystem,
        eventState.synth,
        isData,
        run);
    /* #endregion */

    // [15] Final per-event tree fill
    /* #region: [15] final per-event tree fill */
    if (!isData && signalScan) {
      signalScanManager.fillEventTree(signalState);
    } else {
      MuonSystem->tree_->Fill();
    }
    /* #endregion */
  }
  /* #endregion */

  // [16] End-of-job writeout for trees and histograms
  /* #region: [16] end-of-job writeout for trees and histograms */
  if (!isData && signalScan) {
    signalScanManager.finalize();
  } else if (!isData) {
    cout << "Filled Total of " << NEvents->GetBinContent(1) << " Events\n";
    cout << "Writing output trees..." << endl;
    outFile->cd();
    MuonSystem->tree_->Write();
    NEvents->Write();
    accep->Write("accep");
    accep_csccsc->Write("accep_csccsc");
    accep_cscdt->Write("accep_cscdt");
    accep_met->Write("accep_met");
    outFile->Close();
  } else {
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
          
