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

#include <chrono>
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
  // for (Long64_t jentry = 0; jentry < nEntries; jentry++) {
  for (Long64_t jentry = 0; jentry < 1000; jentry++) {

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
        *this,
        helper.get(),
        analysisTag);

    /* #endregion */

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

    /* #region: [5] event identity and nominal event weight */
    fillEventIdentityAndNominalWeight(
        MuonSystem,
        NEvents.get(),
        isData,
        Generator_weight,
        run,
        luminosityBlock,
        event);
    if (shouldSkipDataRun(isData, run))
      continue;
    /* #endregion */

    /* #region: [6] MC-only truth and per-event MC weights */
    if (!isData)
      fillMcTruthAndPileupWeights(*this, helper.get(), MuonSystem);
    /* #endregion */

    /* #region: [7] event observables and acceptance counters */
    fillEventObservables(
        MuonSystem,
        PV_npvs,
        Rho_fixedGridRhoFastjetAll,
        PFMET_pt,
        PFMET_phi,
        PuppiMET_pt,
        PuppiMET_phi);
    fillAcceptanceCounters(
        isData,
        signalScan,
        signalState,
        signalScanManager,
        MuonSystem,
        accep_csccsc.get(),
        accep_cscdt.get(),
        Generator_weight);
    /* #endregion */

    /* #region: [8] event filters, trigger bits, and jet veto maps */
    const EventCutState cuts = buildEventCutState(
        *this,
        helper.get(),
        analysisTag,
        isData,
        run,
        eventState.synth);
    writeEventCutState(MuonSystem, cuts);
    /* #endregion */

    /* #region: [9] lepton object selection and lepton branch fill */
    eventState.leptons = buildSelectedLeptons(
        *this,
        eventState.synth);
    fillLeptonBranches(MuonSystem, eventState.leptons);
    /* #endregion */

    /* #region: [10] jet selection, JES propagation, and jet branch fill */
    eventState.jetStage = buildJetStageResult(
        *this,
        helper.get(),
        run,
        eventState.synth,
        eventState.leptons);

    fillJetBranches(MuonSystem, eventState.jetStage.selectedJets);
    fillPuppiMetJesFromShift(*this, MuonSystem, eventState.jetStage);
    /* #endregion */

    /* #region: [11] rechit-level pass-through fill (CSC/DT/RPC + optional aliases) */
    fillRawRechits(
        MuonSystem,
        *this);
    /* #endregion */

    /* #region: [12] ring occupancy summaries from raw rechits */
    MuonSystem->nDtRings = countDtRingsFromRecHits(
        ndtRecHits,
        dtRecHits_Station,
        dtRecHits_Wheel);
    /* #endregion */

    /* #region: [13] CSC rechit clustering, feature fill, and matching/veto variables */
    processCscClusterStage(
        *this,
        helper.get(),
        MuonSystem,
        eventState.synth,
        isData,
        run);
    /* #endregion */

    /* #region: [14] DT/RPC clustering, DT feature fill, and matching/veto variables */
    processDtRpcClusterStage(
        *this,
        helper.get(),
        MuonSystem,
        eventState.synth,
        isData,
        run);
    /* #endregion */

    /* #region: [15] final per-event tree write */
    fillEventTreeForCurrentMode(
        isData,
        signalScan,
        MuonSystem,
        signalScanManager,
        signalState);
    /* #endregion */
  }
  /* #endregion */

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
          
