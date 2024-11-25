#include "RazorHelper.h"
#include "TFile.h"
#include "TH1F.h"
#include "TKey.h"
#include "TTree.h"
#include "TTreeFormula.h"
#include "TreeMuonSystemCombination.h"
#include "llp_event.h"
#include "nano_events.h"
#include <TRandom3.h>
#include <assert.h>
#include <fstream>
#include <iostream>
#include <iterator>
#include <regex>
#include <sstream>

#include "TString.h"
#include <vector>

#define N_MAX_RECHITS 20000
#define N_MAX_SEGMENT 1000

using namespace std;

// get list of files to open, add normalization branch to the tree in each file
int main(int argc, char *argv[]) {
    // parse input list to get names of ROOT files
    if (argc < 3) {
        cerr << "usage MergeNtuple [match_object.root] [output.root]" << endl;
        return -1;
    }
    TFile *matchFile = new TFile(argv[1]);
    TTree *inp_nano = (TTree *)matchFile->Get("inp1");
    TTree *inp_ntuple = (TTree *)matchFile->Get("inp2");
    TTree *matched_idx_tree = (TTree *)matchFile->Get("idx");
    // inp1
    //  |-- str: paths
    //  `-- i64: lengths
    // inp2
    //  |-- str: paths
    //  `-- i64: lengths
    // idx
    //  |-- i64: 1
    //  `-- i64: 2

    string output_fname(argv[2]);

    // load ntuple
    TChain *ntuple_chain = new TChain();

    char ntuple_path[1000];
    inp_ntuple->SetBranchAddress("paths", ntuple_path);

    for (int i = 0; i < inp_ntuple->GetEntries(); i++) {
        inp_ntuple->GetEntry(i);
        cout << "curFileName = " << ntuple_path << endl;
        // cout << "wtf" << endl;
        if (i == 0) {
            // checks root file structure and add first file
            cout << "[INFO]: loading file: " << ntuple_path << endl;
            TFile *f_0 = TFile::Open(ntuple_path, "READ");
            if (f_0->GetDirectory("ntuples")) {
                ntuple_chain->SetName("ntuples/llp");
                cout << "[INFO]: default configuration for tchain" << endl;
            } else {
                ntuple_chain->SetName("llp");
                cout << "[INFO]: alternative configuration for tchain" << endl;
            }
            ntuple_chain->Add(ntuple_path);
            delete f_0;
        } else {
            // Addind remaining files after file structure is decided
            ntuple_chain->Add(ntuple_path);
        }
    }

    if (ntuple_chain == NULL)
        return -1;

    bool is_mc = ntuple_chain->GetBranch("gLLP_eta") != NULL;
    cout << "Loaded ntuples: " << inp_ntuple->GetEntries() << " files and " << ntuple_chain->GetEntries() << " events. is_mc = " << is_mc << endl;

    ////////////////////////
    ////// Load nanoAOD files
    ////////////////////////
    TChain *nano_chain = new TChain();

    char nano_path[1000];
    inp_nano->SetBranchAddress("paths", nano_path);

    for (int i = 0; i < inp_nano->GetEntries(); i++) {
        inp_nano->GetEntry(i);
        cout << "curFileName = " << nano_path << endl;
        if (i == 0) {
            // checks root file structure and add first file
            cout << "[INFO]: loading file: " << nano_path << endl;
            TFile *f_0 = TFile::Open(nano_path, "READ");
            f_0->ls();
            nano_chain->SetName("Events");
            nano_chain->Add(nano_path);
            delete f_0;
        } else {
            // Addind remaining files after file structure is decided
            nano_chain->Add(nano_path);
        }
    }
    cout << "Loaded nanoAOD: " << inp_nano->GetEntries() << " files and " << nano_chain->GetEntries() << " events" << endl;
    if (nano_chain == NULL)
        return -1;

    //*****************************************************************************************
    // Make map of event number in tree 1 to event number in tree 2
    //*****************************************************************************************

    nano_events nano = nano_events(nano_chain);
    llp_event ntuple = llp_event(ntuple_chain);

    // find intersection of the two maps
    Long64_t idx_nano, idx_ntuple;
    matched_idx_tree->SetBranchAddress("1", &idx_nano);
    matched_idx_tree->SetBranchAddress("2", &idx_ntuple);

    std::vector<std::pair<int64_t, int64_t>> idx_pairs(matched_idx_tree->GetEntries());
    for (int i = 0; i < matched_idx_tree->GetEntries(); i++) {
        matched_idx_tree->GetEntry(i);
        // std::cout << idx_nano << " --- " << idx_ntuple << std::endl;
        idx_pairs[i] = make_pair(idx_nano, idx_ntuple);
    }

    // cout << "Sorting ... " << endl;

    // sort(idx_pairs.begin(), idx_pairs.end(), [](const pair<ulong, ulong> &left, const pair<ulong, ulong> &right) { return left.second > right.second; });

    ulong nMatchedEvents = idx_pairs.size();
    ulong NEventsTreeNTuple = ntuple_chain->GetEntries();

    cout << "Matched events    = " << nMatchedEvents << " / " << NEventsTreeNTuple << " \n";
    cout << "Un-Matched events = " << NEventsTreeNTuple - nMatchedEvents << " / " << NEventsTreeNTuple << " \n";
    if (nMatchedEvents == 0) {
        cout << "No matched events found. Exiting..." << endl;
        return 0;
    }
    // for (int i = 0; i < idx_pairs.size(); i++) {
    //     auto pair = idx_pairs[i];
    //     nano_chain->GetEntry(pair.first);
    //     ntuple_chain->GetEntry(pair.second);
    //     std::cout << pair.second << " --- " << ntuple.runNum << ":" << ntuple.lumiNum << ":" << ntuple.eventNum << std::endl;
    //     std::cout << pair.first << " ---" << nano.run << " : " << nano.luminosityBlock << " : " << nano.event << std::endl;
    // }
    //*****************************************************************************************
    // Produce Output Tree
    //*****************************************************************************************

    // clone tree with 0 entries, copy all the branch addresses only

    TFile *outputFile = new TFile(output_fname.c_str(), "RECREATE", "", 202);
    outputFile->cd();
    TTree *outputTree = nano_chain->CloneTree(0);

    ////////////////////////////////////
    // declare new branches to be added
    //////////////////////////////////
    Int_t nCscSeg, nCscRechits, nDtRechits, nDtSeg, nRpc;
    outputTree->Branch("nCscSeg", &nCscSeg, "nCscSeg/I");
    outputTree->Branch("nCscRechits", &nCscRechits, "nCscRechits/I");
    outputTree->Branch("nDtSeg", &nDtSeg, "nDtSeg/I");
    outputTree->Branch("nDtRechits", &nDtRechits, "nDtRechits/I");
    outputTree->Branch("nRpc", &nRpc, "nRpc/I");

    outputTree->SetBranchAddress("nCscSeg", &nCscSeg);
    outputTree->SetBranchAddress("nCscRechits", &nCscRechits);
    outputTree->SetBranchAddress("nDtSeg", &nDtSeg);
    outputTree->SetBranchAddress("nDtRechits", &nDtRechits);
    outputTree->SetBranchAddress("nRpc", &nRpc);
    ////////////////////////////////////////////////
    ///// Float for rechits related variables
    ////////////////////////////////////////////////

    int numFloatBranches = 12;

    const char *addBranchNamesFloat[numFloatBranches]{
        "cscRechitsPhi",
        "cscRechitsEta",
        "cscRechitsX",
        "cscRechitsY",
        "cscRechitsZ",
        "cscRechitsTpeak",
        "cscRechitsTwire",
        "dtRechitCorrectX",
        "dtRechitCorrectY",
        "dtRechitCorrectZ",
        "dtRechitCorrectEta",
        "dtRechitCorrectPhi",
    };
    Float_t *addBranchesRazorVarFloat[numFloatBranches]{
        ntuple.cscRechitsPhi,
        ntuple.cscRechitsEta,
        ntuple.cscRechitsX,
        ntuple.cscRechitsY,
        ntuple.cscRechitsZ,
        ntuple.cscRechitsTpeak,
        ntuple.cscRechitsTwire,
        ntuple.dtRechitCorrectX,
        ntuple.dtRechitCorrectY,
        ntuple.dtRechitCorrectZ,
        ntuple.dtRechitCorrectEta,
        ntuple.dtRechitCorrectPhi,
    };
    Float_t addBranchesInputVarFloat[numFloatBranches][N_MAX_RECHITS];

    for (int i = 0; i < numFloatBranches; i++) {
        cout << "Adding Branch: " << addBranchNamesFloat[i] << "\n";
        if (string(addBranchNamesFloat[i]).find("cscRechits") != string::npos)
            outputTree->Branch(addBranchNamesFloat[i], addBranchesInputVarFloat[i], TString::Format("%s[nCscRechits]/F", addBranchNamesFloat[i]));
        if (string(addBranchNamesFloat[i]).find("dtRechit") != string::npos)
            outputTree->Branch(addBranchNamesFloat[i], addBranchesInputVarFloat[i], TString::Format("%s[nDtRechits]/F", addBranchNamesFloat[i]));
        outputTree->SetBranchAddress(addBranchNamesFloat[i], addBranchesInputVarFloat[i]);
    }

    ////////////////////////////////////////////////
    ///// Float for segment related variables
    ////////////////////////////////////////////////

    int numFloatSegBranches = 11;

    const char *addBranchNamesFloatSeg[numFloatSegBranches]{
        "cscSegPhi", "cscSegEta",
        "dtSegPhi", "dtSegEta",
        "rpcPhi", "rpcEta", "rpcX", "rpcY", "rpcZ", "rpcT", "rpcTError"};
    Float_t *addBranchesRazorVarFloatSeg[numFloatSegBranches]{
        ntuple.cscSegPhi, ntuple.cscSegEta,
        ntuple.dtSegPhi, ntuple.dtSegEta,
        ntuple.rpcPhi, ntuple.rpcEta, ntuple.rpcX, ntuple.rpcY, ntuple.rpcZ, ntuple.rpcT, ntuple.rpcTError};
    Float_t addBranchesInputVarFloatSeg[numFloatSegBranches][N_MAX_SEGMENT];

    for (int i = 0; i < numFloatSegBranches; i++) {
        cout << "Adding Branch: " << addBranchNamesFloatSeg[i] << "\n";
        if (string(addBranchNamesFloatSeg[i]).find("cscSeg") != string::npos)
            outputTree->Branch(addBranchNamesFloatSeg[i], addBranchesInputVarFloatSeg[i], TString::Format("%s[nCscSeg]/F", addBranchNamesFloatSeg[i]));
        if (string(addBranchNamesFloatSeg[i]).find("dtSeg") != string::npos)
            outputTree->Branch(addBranchNamesFloatSeg[i], addBranchesInputVarFloatSeg[i], TString::Format("%s[nDtSeg]/F", addBranchNamesFloatSeg[i]));
        if (string(addBranchNamesFloatSeg[i]).find("rpc") != string::npos)
            outputTree->Branch(addBranchNamesFloatSeg[i], addBranchesInputVarFloatSeg[i], TString::Format("%s[nRpc]/F", addBranchNamesFloatSeg[i]));
        outputTree->SetBranchAddress(addBranchNamesFloatSeg[i], addBranchesInputVarFloatSeg[i]);
    }
    ////////////////////////////////////////////////
    /////integer for rechits related variables
    ////////////////////////////////////////////////

    int numIntBranches = 6;
    const char *addBranchNamesInt[numIntBranches]{"cscRechitsChamber", "cscRechitsStation", "cscRechitsDetId", "dtRechitStation", "dtRechitWheel", "dtRechitSuperLayer"};

    Int_t *addBranchesRazorVarInt[numIntBranches]{
        ntuple.cscRechitsChamber, ntuple.cscRechitsStation, ntuple.cscRechitsDetId,
        ntuple.dtRechitStation, ntuple.dtRechitWheel, ntuple.dtRechitSuperLayer};

    Int_t addBranchesInputVarInt[numIntBranches][N_MAX_RECHITS];

    for (int i = 0; i < numIntBranches; i++) {
        cout << "Adding Branch: " << addBranchNamesInt[i] << "\n";
        if (string(addBranchNamesInt[i]).find("cscRechits") != string::npos)
            outputTree->Branch(addBranchNamesInt[i], addBranchesInputVarInt[i], TString::Format("%s[nCscRechits]/I", addBranchNamesInt[i]));
        else if (string(addBranchNamesInt[i]).find("dtRechit") != string::npos)
            outputTree->Branch(addBranchNamesInt[i], addBranchesInputVarInt[i], TString::Format("%s[nDtRechits]/I", addBranchNamesInt[i]));
        outputTree->SetBranchAddress(addBranchNamesInt[i], addBranchesInputVarInt[i]);
    }

    // ////////////////////////////////////////////////
    // /////integer for segments related variables
    // ////////////////////////////////////////////////

    int numIntSegBranches = 11;
    const char *addBranchNamesIntSeg[numIntSegBranches]{
        "cscSegChamber", "cscSegStation", "cscSegNRecHits", "dtSegStation", "dtSegWheel",
        "rpcBx", "rpcRegion", "rpcRing", "rpcSector", "rpcStation", "rpcLayer"};

    Int_t *addBranchesRazorVarIntSeg[numIntSegBranches]{
        ntuple.cscSegChamber, ntuple.cscSegStation, ntuple.cscSegNRecHits, ntuple.dtSegStation, ntuple.dtSegWheel,
        ntuple.rpcBx, ntuple.rpcRegion, ntuple.rpcRing, ntuple.rpcSector, ntuple.rpcStation, ntuple.rpcLayer

    };
    Int_t addBranchesInputVarIntSeg[numIntSegBranches][N_MAX_SEGMENT];

    for (int i = 0; i < numIntSegBranches; i++) {
        cout << "Adding Branch: " << addBranchNamesIntSeg[i] << "\n";
        if (string(addBranchNamesIntSeg[i]).find("cscSeg") != string::npos)
            outputTree->Branch(addBranchNamesIntSeg[i], addBranchesInputVarIntSeg[i], TString::Format("%s[nCscSeg]/I", addBranchNamesIntSeg[i]));
        else if (string(addBranchNamesIntSeg[i]).find("dtSeg") != string::npos)
            outputTree->Branch(addBranchNamesIntSeg[i], addBranchesInputVarIntSeg[i], TString::Format("%s[nDtSeg]/I", addBranchNamesIntSeg[i]));
        else if (string(addBranchNamesIntSeg[i]).find("rpc") != string::npos)
            outputTree->Branch(addBranchNamesIntSeg[i], addBranchesInputVarIntSeg[i], TString::Format("%s[nRpc]/I", addBranchNamesIntSeg[i]));
        outputTree->SetBranchAddress(addBranchNamesIntSeg[i], addBranchesInputVarIntSeg[i]);
    }

    int nGLLP;
    float gLLP_eta[2], gLLP_phi[2], gLLP_csc[2], gLLP_dt[2], gLLP_beta[2], gLLP_e[2], gLLP_pt[2], gLLP_decay_vertex_r[2], gLLP_decay_vertex_x[2], gLLP_decay_vertex_y[2], gLLP_decay_vertex_z[2];
    if (is_mc) {
        outputTree->Branch("nGLLP", &nGLLP, "nGLLP/I");
        outputTree->Branch("gLLP_eta", gLLP_eta, "gLLP_eta[nGLLP]/F");
        outputTree->Branch("gLLP_phi", gLLP_phi, "gLLP_phi[nGLLP]/F");
        outputTree->Branch("gLLP_csc", gLLP_csc, "gLLP_csc[nGLLP]/F");
        outputTree->Branch("gLLP_dt", gLLP_dt, "gLLP_dt[nGLLP]/F");
        outputTree->Branch("gLLP_beta", gLLP_beta, "gLLP_beta[nGLLP]/F");
        outputTree->Branch("gLLP_e", gLLP_e, "gLLP_e[nGLLP]/F");
        outputTree->Branch("gLLP_pt", gLLP_pt, "gLLP_pt[nGLLP]/F");
        outputTree->Branch("gLLP_decay_vertex_r", gLLP_decay_vertex_r, "gLLP_decay_vertex_r[nGLLP]/F");
        outputTree->Branch("gLLP_decay_vertex_x", gLLP_decay_vertex_x, "gLLP_decay_vertex_x[nGLLP]/F");
        outputTree->Branch("gLLP_decay_vertex_y", gLLP_decay_vertex_y, "gLLP_decay_vertex_y[nGLLP]/F");
        outputTree->Branch("gLLP_decay_vertex_z", gLLP_decay_vertex_z, "gLLP_decay_vertex_z[nGLLP]/F");
        outputTree->SetBranchAddress("nGLLP", &nGLLP);
        outputTree->SetBranchAddress("gLLP_eta", gLLP_eta);
        outputTree->SetBranchAddress("gLLP_phi", gLLP_phi);
        outputTree->SetBranchAddress("gLLP_csc", gLLP_csc);
        outputTree->SetBranchAddress("gLLP_dt", gLLP_dt);
        outputTree->SetBranchAddress("gLLP_beta", gLLP_beta);
        outputTree->SetBranchAddress("gLLP_e", gLLP_e);
        outputTree->SetBranchAddress("gLLP_pt", gLLP_pt);
        outputTree->SetBranchAddress("gLLP_decay_vertex_r", gLLP_decay_vertex_r);
        outputTree->SetBranchAddress("gLLP_decay_vertex_x", gLLP_decay_vertex_x);
        outputTree->SetBranchAddress("gLLP_decay_vertex_y", gLLP_decay_vertex_y);
        outputTree->SetBranchAddress("gLLP_decay_vertex_z", gLLP_decay_vertex_z);
    }

    /////////////////////////////////////////////////////
    /////// FILL OUTPUT TREE
    /////////////////////////////////////////////////////

    int idx = 0;
    for (auto it : idx_pairs) {
        auto i_nano = it.first;
        auto i_ntuple = it.second;
        if (idx % 1000 == 0)
            cout << "Processing entry " << idx << "\n";

        nano_chain->GetEntry(i_nano);
        ntuple_chain->GetEntry(i_ntuple);
        if (nano.event != ntuple.eventNum || nano.luminosityBlock != ntuple.lumiNum || nano.run != ntuple.runNum) {
            std::cerr << "Mismatched event numbers: nano: " << nano.run << ":" << nano.luminosityBlock << ":" << nano.event << std::endl;
            std::cerr << "Mismatched event numbers: ntuple: " << ntuple.runNum << ":" << ntuple.lumiNum << ":" << ntuple.eventNum << std::endl;
            throw std::runtime_error("Mismatched event numbers");
        }

        // fill the new added branches

        nCscRechits = ntuple.ncscRechits;
        nCscSeg = ntuple.nCscSeg;
        nDtRechits = ntuple.nDtRechits;
        nDtSeg = ntuple.nDtSeg;
        nRpc = ntuple.nRpc;

        for (int i = 0; i < numFloatBranches; i++) {
            int temp_nhits = 0;
            if (string(addBranchNamesFloat[i]).find("cscRechits") != string::npos)
                temp_nhits = nCscRechits;
            if (string(addBranchNamesFloat[i]).find("dtRechit") != string::npos)
                temp_nhits = nDtRechits;
            for (int j = 0; j < min(temp_nhits, N_MAX_RECHITS); j++)
                addBranchesInputVarFloat[i][j] = addBranchesRazorVarFloat[i][j];
        }

        for (int i = 0; i < numFloatSegBranches; i++) {
            int temp_nhits = 0;
            if (string(addBranchNamesFloatSeg[i]).find("cscSeg") != string::npos)
                temp_nhits = nCscSeg;
            if (string(addBranchNamesFloatSeg[i]).find("dtSeg") != string::npos)
                temp_nhits = nDtSeg;
            if (string(addBranchNamesFloatSeg[i]).find("rpc") != string::npos)
                temp_nhits = nRpc;
            for (int j = 0; j < min(temp_nhits, N_MAX_SEGMENT); j++)
                addBranchesInputVarFloatSeg[i][j] = addBranchesRazorVarFloatSeg[i][j];
        }

        for (int i = 0; i < numIntBranches; i++) {
            int temp_nhits = 0;
            if (string(addBranchNamesInt[i]).find("cscRechits") != string::npos)
                temp_nhits = nCscRechits;
            if (string(addBranchNamesInt[i]).find("dtRechit") != string::npos)
                temp_nhits = nDtRechits;
            for (int j = 0; j < min(temp_nhits, N_MAX_RECHITS); j++)
                addBranchesInputVarInt[i][j] = addBranchesRazorVarInt[i][j];
        }

        for (int i = 0; i < numIntSegBranches; i++) {
            int temp_nhits = 0;
            if (string(addBranchNamesIntSeg[i]).find("cscSeg") != string::npos)
                temp_nhits = nCscSeg;
            if (string(addBranchNamesIntSeg[i]).find("dtSeg") != string::npos)
                temp_nhits = nDtSeg;
            if (string(addBranchNamesIntSeg[i]).find("rpc") != string::npos)
                temp_nhits = nRpc;
            for (int j = 0; j < min(temp_nhits, N_MAX_SEGMENT); j++)
                addBranchesInputVarIntSeg[i][j] = addBranchesRazorVarIntSeg[i][j];
        }

        if (is_mc) {
            nGLLP = 2;
            for (int i = 0; i < nGLLP; i++) {
                gLLP_eta[i] = ntuple.gLLP_eta[i];
                gLLP_phi[i] = ntuple.gLLP_phi[i];
                gLLP_csc[i] = ntuple.gLLP_csc[i];
                gLLP_dt[i] = ntuple.gLLP_dt[i];
                gLLP_beta[i] = ntuple.gLLP_beta[i];
                gLLP_e[i] = ntuple.gLLP_e[i];
                gLLP_pt[i] = ntuple.gLLP_pt[i];
                double x = ntuple.gLLP_decay_vertex_x[i];
                double y = ntuple.gLLP_decay_vertex_y[i];
                gLLP_decay_vertex_x[i] = x;
                gLLP_decay_vertex_y[i] = y;
                gLLP_decay_vertex_z[i] = ntuple.gLLP_decay_vertex_z[i];
                gLLP_decay_vertex_r[i] = sqrt(x * x + y * y);
            }
        }

        outputTree->Fill();
        idx++;
    }
    // save information
    cout << "Filled Total of " << outputTree->GetEntries() << " Matched Events\n";
    cout << "Writing output trees...\n";
    outputFile->cd();
    cout << "CD done\n";

    outputTree->Write();
    cout << "Write done\n";
    outputFile->Close();
    cout << "Close done\n";
    delete outputFile;
    return 0;
}