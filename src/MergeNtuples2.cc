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
#include <sstream>

#include "TString.h"
#include <vector>
using namespace std;

#define N_MAX_RECHITS 20000
#define N_MAX_SEGMENT 1000

std::string ParseCommandLine(int argc, char *argv[], std::string opt) {
    for (int i = 1; i < argc; i++) {
        std::string tmp(argv[i]);
        if (tmp.find(opt) != std::string::npos) {
            if (tmp.find("=") != std::string::npos)
                return tmp.substr(tmp.find_last_of("=") + 1);
            if (tmp.find("--") != std::string::npos)
                return "yes";
        }
    }

    return "";
};

// get list of files to open, add normalization branch to the tree in each file
int main(int argc, char *argv[]) {
    // parse input list to get names of ROOT files
    if (argc < 4) {
        cerr << "usage MergeNtuple [inputfile1] [inputfile2] [outputfile] [analysisTag]" << endl;
        return -1;
    }
    // string filenameNTuplers(argv[1]);
    string filenameNTuplerList(argv[1]);
    string filenameNanoAODList(argv[2]);
    string outputfilename(argv[3]);
    string analysisTag(argv[4]); // analysisTag

    if (analysisTag == "") {
        analysisTag = "Razor2016_80X";
    }
    // create output file
    TFile *outputFile = new TFile(outputfilename.c_str(), "RECREATE");

    ////////////////////////
    ////// Load ntuples
    ////////////////////////
    TChain *ntupleChain = new TChain();

    // TFile* f_0 = TFile::Open( filenameNTuplers.c_str() );
    // ntupleChain->SetName("llp");
    // ntupleChain->Add(filenameNTuplers.c_str());
    // delete f_0;
    string curNtupleName;
    cout << "filenameNTuplerList.c_str() = " << filenameNTuplerList.c_str() << endl;
    ifstream inputNtuple(filenameNTuplerList.c_str());
    int NtuplesLoaded = 0;
    if (!inputNtuple) {
        cerr << "Error: input ntuple file list not found!" << endl;
        return -1;
    }
    while (getline(inputNtuple, curNtupleName)) {
        cout << "curFileName = " << curNtupleName << endl;
        if (NtuplesLoaded == 0) {
            // checks root file structure and add first file
            std::cout << "[INFO]: loading file: " << curNtupleName.c_str() << std::endl;
            TFile *f_0 = TFile::Open(curNtupleName.c_str());
            if (f_0->GetDirectory("ntuples")) {
                ntupleChain->SetName("ntuples/llp");
                std::cout << "[INFO]: default configuration for tchain" << std::endl;
            } else {
                ntupleChain->SetName("llp");
                std::cout << "[INFO]: alternative configuration for tchain" << std::endl;
            }
            ntupleChain->Add(curNtupleName.c_str());
            delete f_0;

        } else {
            // Addind remaining files after file structure is decided
            ntupleChain->Add(curNtupleName.c_str());
        }
        NtuplesLoaded++;
    }
    if (ntupleChain == NULL)
        return -1;

    bool is_mc = ntupleChain->GetBranch("gLLP_eta") != NULL;
    is_mc = false;
    cout << "Loaded ntuples: " << NtuplesLoaded << " files and " << ntupleChain->GetEntries() << " events. is_mc = " << is_mc << endl;

    ////////////////////////
    ////// Load nanoAOD files
    ////////////////////////
    TChain *nanoChain = new TChain();
    string curFileName;
    ifstream inputFile(filenameNanoAODList.c_str());
    int NFilesLoaded = 0;
    if (!inputFile) {
        cerr << "Error: input nanoAOD file list not found!" << endl;
        return -1;
    }

    while (getline(inputFile, curFileName)) {
        cout << "curFileName = " << curFileName << endl;
        if (NFilesLoaded == 0) {
            // checks root file structure and add first file
            std::cout << "[INFO]: loading file: " << curFileName.c_str() << std::endl;
            TFile *f_0 = TFile::Open(curFileName.c_str());
            f_0->ls();
            nanoChain->SetName("Events");
            nanoChain->Add(curFileName.c_str());
            delete f_0;
        } else {
            // Addind remaining files after file structure is decided
            nanoChain->Add(curFileName.c_str());
        }
        NFilesLoaded++;
    }
    std::cout << "Loaded nanoAOD: " << NFilesLoaded << " files and " << nanoChain->GetEntries() << " events" << endl;
    if (nanoChain == NULL)
        return -1;

    uint NEventsTree1 = ntupleChain->GetEntries();
    uint NEventsTree2 = nanoChain->GetEntries();

    //*****************************************************************************************
    // Make map of event number in tree 1 to event number in tree 2
    //*****************************************************************************************
    int nMatchedEvents = 0;
    std::map<uint, uint> EventIndexToEventIndexMap;

    nano_events nano = nano_events(nanoChain);
    llp_event ntuple = llp_event(ntupleChain);

    // loop over tree2
    std::vector<std::pair<uint, ulong>> eventList2;
    std::vector<bool> matchedevent;

    cout << "building nanoAOD tree index map" << endl;
    for (UInt_t m = 0; m < NEventsTree2; m++) {
        if (m % 10000 == 0)
            cout << "Processing entry " << m << endl;
        nano.fChain->GetEntry(m);
        std::pair<uint, ulong> p(nano.run, nano.event);
        eventList2.push_back(p);
    }

    cout << "Looping over ntuple tree to find matched events" << endl;
    ;
    // loop over tree1
    for (uint n = 0; n < NEventsTree1; n++) {
        ntupleChain->GetEntry(n);
        if (n % 10000 == 0)
            cout << "Event " << n << "\n";

        bool matchFound = false;
        for (uint m = 0; m < eventList2.size(); m++) {
            if (eventList2[m].first != ntuple.runNum)
                continue;
            if (eventList2[m].second != ntuple.eventNum)
                continue;
            EventIndexToEventIndexMap[n] = m;
            matchFound = true;
            nMatchedEvents++;
        }

        if (matchFound)
            matchedevent.push_back(true);
        else
            matchedevent.push_back(false);
    }

    cout << "Matched events    = " << nMatchedEvents << " / " << NEventsTree1 << " \n";
    cout << "Un-Matched events = " << NEventsTree1 - nMatchedEvents << " / " << NEventsTree1 << " \n";
    //*****************************************************************************************
    // Produce Output Tree
    //*****************************************************************************************

    // clone tree with 0 entries, copy all the branch addresses only
    outputFile->cd();
    TTree *outputTree = nanoChain->CloneTree(0);

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
        if (string(addBranchNamesFloat[i]).find("cscRechits") != std::string::npos)
            outputTree->Branch(addBranchNamesFloat[i], addBranchesInputVarFloat[i], TString::Format("%s[nCscRechits]/F", addBranchNamesFloat[i]));
        if (string(addBranchNamesFloat[i]).find("dtRechit") != std::string::npos)
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
        if (string(addBranchNamesFloatSeg[i]).find("cscSeg") != std::string::npos)
            outputTree->Branch(addBranchNamesFloatSeg[i], addBranchesInputVarFloatSeg[i], TString::Format("%s[nCscSeg]/F", addBranchNamesFloatSeg[i]));
        if (string(addBranchNamesFloatSeg[i]).find("dtSeg") != std::string::npos)
            outputTree->Branch(addBranchNamesFloatSeg[i], addBranchesInputVarFloatSeg[i], TString::Format("%s[nDtSeg]/F", addBranchNamesFloatSeg[i]));
        if (string(addBranchNamesFloatSeg[i]).find("rpc") != std::string::npos)
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
        if (string(addBranchNamesInt[i]).find("cscRechits") != std::string::npos)
            outputTree->Branch(addBranchNamesInt[i], addBranchesInputVarInt[i], TString::Format("%s[nCscRechits]/I", addBranchNamesInt[i]));
        else if (string(addBranchNamesInt[i]).find("dtRechit") != std::string::npos)
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
        if (string(addBranchNamesIntSeg[i]).find("cscSeg") != std::string::npos)
            outputTree->Branch(addBranchNamesIntSeg[i], addBranchesInputVarIntSeg[i], TString::Format("%s[nCscSeg]/I", addBranchNamesIntSeg[i]));
        else if (string(addBranchNamesIntSeg[i]).find("dtSeg") != std::string::npos)
            outputTree->Branch(addBranchNamesIntSeg[i], addBranchesInputVarIntSeg[i], TString::Format("%s[nDtSeg]/I", addBranchNamesIntSeg[i]));
        else if (string(addBranchNamesIntSeg[i]).find("rpc") != std::string::npos)
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

    for (uint n = 0; n < NEventsTree1; n++) {
        if (n % 10000 == 0)
            cout << "Processing entry " << n << "\n";

        // Check if found a match
        if (!matchedevent[n])
            continue;
        // if (n == 6604)cout<< n << ", " << nano.event<< endl;
        // if (ntuple.ncscRechits >= 20000) cout<< n << ", " << ntuple.ncscRechits<< endl;
        // if (ntuple.nCscSeg >= 300) cout<< n << ", " << ntuple.nCscSeg<< endl;
        nanoChain->GetEntry(EventIndexToEventIndexMap[n]);
        ntupleChain->GetEntry(n);

        /////////////////////////////////////////////////////
        //////// cosmic shower veto ///////////
        /////////////////////////////////////////////////////

        int nCscRechitsChamberPlus11 = 0, nCscRechitsChamberPlus12 = 0, nCscRechitsChamberPlus13 = 0,
            nCscRechitsChamberPlus21 = 0, nCscRechitsChamberPlus22 = 0, nCscRechitsChamberPlus31 = 0, nCscRechitsChamberPlus32 = 0,
            nCscRechitsChamberPlus41 = 0, nCscRechitsChamberPlus42 = 0,
            nCscRechitsChamberMinus11 = 0, nCscRechitsChamberMinus12 = 0, nCscRechitsChamberMinus13 = 0,
            nCscRechitsChamberMinus21 = 0, nCscRechitsChamberMinus22 = 0, nCscRechitsChamberMinus31, nCscRechitsChamberMinus32 = 0,
            nCscRechitsChamberMinus41 = 0, nCscRechitsChamberMinus42 = 0, nCscRings = 0;
        for (int i = 0; i < ntuple.ncscRechits; i++) {
            if (ntuple.cscRechitsChamber[i] == 11)
                nCscRechitsChamberPlus11++;
            if (ntuple.cscRechitsChamber[i] == 12)
                nCscRechitsChamberPlus12++;
            if (ntuple.cscRechitsChamber[i] == 13)
                nCscRechitsChamberPlus13++;
            if (ntuple.cscRechitsChamber[i] == 21)
                nCscRechitsChamberPlus21++;
            if (ntuple.cscRechitsChamber[i] == 22)
                nCscRechitsChamberPlus22++;
            if (ntuple.cscRechitsChamber[i] == 31)
                nCscRechitsChamberPlus31++;
            if (ntuple.cscRechitsChamber[i] == 32)
                nCscRechitsChamberPlus32++;
            if (ntuple.cscRechitsChamber[i] == 41)
                nCscRechitsChamberPlus41++;
            if (ntuple.cscRechitsChamber[i] == 42)
                nCscRechitsChamberPlus42++;
            if (ntuple.cscRechitsChamber[i] == -11)
                nCscRechitsChamberMinus11++;
            if (ntuple.cscRechitsChamber[i] == -12)
                nCscRechitsChamberMinus12++;
            if (ntuple.cscRechitsChamber[i] == -13)
                nCscRechitsChamberMinus13++;
            if (ntuple.cscRechitsChamber[i] == -21)
                nCscRechitsChamberMinus21++;
            if (ntuple.cscRechitsChamber[i] == -22)
                nCscRechitsChamberMinus22++;
            if (ntuple.cscRechitsChamber[i] == -31)
                nCscRechitsChamberMinus31++;
            if (ntuple.cscRechitsChamber[i] == -32)
                nCscRechitsChamberMinus32++;
            if (ntuple.cscRechitsChamber[i] == -41)
                nCscRechitsChamberMinus41++;
            if (ntuple.cscRechitsChamber[i] == -42)
                nCscRechitsChamberMinus42++;
        }
        if (nCscRechitsChamberPlus11 > 50)
            nCscRings++;
        if (nCscRechitsChamberPlus12 > 50)
            nCscRings++;
        if (nCscRechitsChamberPlus13 > 50)
            nCscRings++;
        if (nCscRechitsChamberPlus21 > 50)
            nCscRings++;
        if (nCscRechitsChamberPlus22 > 50)
            nCscRings++;
        if (nCscRechitsChamberPlus31 > 50)
            nCscRings++;
        if (nCscRechitsChamberPlus32 > 50)
            nCscRings++;
        if (nCscRechitsChamberPlus41 > 50)
            nCscRings++;
        if (nCscRechitsChamberPlus42 > 50)
            nCscRings++;
        if (nCscRechitsChamberMinus11 > 50)
            nCscRings++;
        if (nCscRechitsChamberMinus12 > 50)
            nCscRings++;
        if (nCscRechitsChamberMinus13 > 50)
            nCscRings++;
        if (nCscRechitsChamberMinus21 > 50)
            nCscRings++;
        if (nCscRechitsChamberMinus22 > 50)
            nCscRings++;
        if (nCscRechitsChamberMinus31 > 50)
            nCscRings++;
        if (nCscRechitsChamberMinus32 > 50)
            nCscRings++;
        if (nCscRechitsChamberMinus41 > 50)
            nCscRings++;
        if (nCscRechitsChamberMinus42 > 50)
            nCscRings++;

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
        int nDtRings = 0;

        for (int i = 0; i < ntuple.nDtRechits; i++) {

            if (ntuple.dtRechitStation[i] == 1 && ntuple.dtRechitWheel[i] == -2)
                nDTRechitsChamberMinus12++;
            if (ntuple.dtRechitStation[i] == 1 && ntuple.dtRechitWheel[i] == -1)
                nDTRechitsChamberMinus11++;
            if (ntuple.dtRechitStation[i] == 1 && ntuple.dtRechitWheel[i] == 0)
                nDTRechitsChamber10++;
            if (ntuple.dtRechitStation[i] == 1 && ntuple.dtRechitWheel[i] == 1)
                nDTRechitsChamberPlus11++;
            if (ntuple.dtRechitStation[i] == 1 && ntuple.dtRechitWheel[i] == 2)
                nDTRechitsChamberPlus12++;
            if (ntuple.dtRechitStation[i] == 2 && ntuple.dtRechitWheel[i] == -2)
                nDTRechitsChamberMinus22++;
            if (ntuple.dtRechitStation[i] == 2 && ntuple.dtRechitWheel[i] == -1)
                nDTRechitsChamberMinus21++;
            if (ntuple.dtRechitStation[i] == 2 && ntuple.dtRechitWheel[i] == 0)
                nDTRechitsChamber20++;
            if (ntuple.dtRechitStation[i] == 2 && ntuple.dtRechitWheel[i] == 1)
                nDTRechitsChamberPlus21++;
            if (ntuple.dtRechitStation[i] == 2 && ntuple.dtRechitWheel[i] == 2)
                nDTRechitsChamberPlus22++;
            if (ntuple.dtRechitStation[i] == 3 && ntuple.dtRechitWheel[i] == -2)
                nDTRechitsChamberMinus32++;
            if (ntuple.dtRechitStation[i] == 3 && ntuple.dtRechitWheel[i] == -1)
                nDTRechitsChamberMinus31++;
            if (ntuple.dtRechitStation[i] == 3 && ntuple.dtRechitWheel[i] == 0)
                nDTRechitsChamber30++;
            if (ntuple.dtRechitStation[i] == 3 && ntuple.dtRechitWheel[i] == 1)
                nDTRechitsChamberPlus31++;
            if (ntuple.dtRechitStation[i] == 3 && ntuple.dtRechitWheel[i] == 2)
                nDTRechitsChamberPlus32++;
            if (ntuple.dtRechitStation[i] == 4 && ntuple.dtRechitWheel[i] == -2)
                nDTRechitsChamberMinus42++;
            if (ntuple.dtRechitStation[i] == 4 && ntuple.dtRechitWheel[i] == -1)
                nDTRechitsChamberMinus41++;
            if (ntuple.dtRechitStation[i] == 4 && ntuple.dtRechitWheel[i] == 0)
                nDTRechitsChamber40++;
            if (ntuple.dtRechitStation[i] == 4 && ntuple.dtRechitWheel[i] == 1)
                nDTRechitsChamberPlus41++;
            if (ntuple.dtRechitStation[i] == 4 && ntuple.dtRechitWheel[i] == 2)
                nDTRechitsChamberPlus42++;
        }

        if (nDTRechitsChamberMinus12 > 50)
            nDtRings++;
        if (nDTRechitsChamberMinus11 > 50)
            nDtRings++;
        if (nDTRechitsChamber10 > 50)
            nDtRings++;
        if (nDTRechitsChamberPlus11 > 50)
            nDtRings++;
        if (nDTRechitsChamberPlus12 > 50)
            nDtRings++;
        if (nDTRechitsChamberMinus22 > 50)
            nDtRings++;
        if (nDTRechitsChamberMinus21 > 50)
            nDtRings++;
        if (nDTRechitsChamber20 > 50)
            nDtRings++;
        if (nDTRechitsChamberPlus21 > 50)
            nDtRings++;
        if (nDTRechitsChamberPlus22 > 50)
            nDtRings++;
        if (nDTRechitsChamberMinus32 > 50)
            nDtRings++;
        if (nDTRechitsChamberMinus31 > 50)
            nDtRings++;
        if (nDTRechitsChamber30 > 50)
            nDtRings++;
        if (nDTRechitsChamberPlus31 > 50)
            nDtRings++;
        if (nDTRechitsChamberPlus32 > 50)
            nDtRings++;
        if (nDTRechitsChamberMinus42 > 50)
            nDtRings++;
        if (nDTRechitsChamberMinus41 > 50)
            nDtRings++;
        if (nDTRechitsChamber40 > 50)
            nDtRings++;
        if (nDTRechitsChamberPlus41 > 50)
            nDtRings++;
        if (nDTRechitsChamberPlus42 > 50)
            nDtRings++;

        if (nDtRings + nCscRings >= 10)
            continue;

        /////////////////////////////////////////////////////
        //////// end of cosmic shower veto ///////////
        /////////////////////////////////////////////////////

        // fill the new added branches
        nCscRechits = ntuple.ncscRechits;
        nCscSeg = ntuple.nCscSeg;
        nDtRechits = ntuple.nDtRechits;
        nDtSeg = ntuple.nDtSeg;
        nRpc = ntuple.nRpc;
        // if (nCscRechits > 20000) continue;

        for (int i = 0; i < numFloatBranches; i++) {
            int temp_nhits = 0;
            if (string(addBranchNamesFloat[i]).find("cscRechits") != std::string::npos)
                temp_nhits = nCscRechits;
            if (string(addBranchNamesFloat[i]).find("dtRechit") != std::string::npos)
                temp_nhits = nDtRechits;
            for (int j = 0; j < temp_nhits; j++)
                addBranchesInputVarFloat[i][j] = addBranchesRazorVarFloat[i][j];
        }

        for (int i = 0; i < numFloatSegBranches; i++) {
            int temp_nhits = 0;
            if (string(addBranchNamesFloatSeg[i]).find("cscSeg") != std::string::npos)
                temp_nhits = nCscSeg;
            if (string(addBranchNamesFloatSeg[i]).find("dtSeg") != std::string::npos)
                temp_nhits = nDtSeg;
            if (string(addBranchNamesFloatSeg[i]).find("rpc") != std::string::npos)
                temp_nhits = nRpc;
            for (int j = 0; j < temp_nhits; j++)
                addBranchesInputVarFloatSeg[i][j] = addBranchesRazorVarFloatSeg[i][j];
        }

        for (int i = 0; i < numIntBranches; i++) {
            int temp_nhits = 0;
            if (string(addBranchNamesInt[i]).find("cscRechits") != std::string::npos)
                temp_nhits = nCscRechits;
            if (string(addBranchNamesInt[i]).find("dtRechit") != std::string::npos)
                temp_nhits = nDtRechits;
            for (int j = 0; j < temp_nhits; j++)
                addBranchesInputVarInt[i][j] = addBranchesRazorVarInt[i][j];
        }

        for (int i = 0; i < numIntSegBranches; i++) {
            int temp_nhits = 0;
            if (string(addBranchNamesIntSeg[i]).find("cscSeg") != std::string::npos)
                temp_nhits = nCscSeg;
            if (string(addBranchNamesIntSeg[i]).find("dtSeg") != std::string::npos)
                temp_nhits = nDtSeg;
            if (string(addBranchNamesIntSeg[i]).find("rpc") != std::string::npos)
                temp_nhits = nRpc;
            for (int j = 0; j < temp_nhits; j++)
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
    }
    // save information
    cout << "Filled Total of " << outputTree->GetEntries() << " Matched Events\n";
    cout << "Writing output trees..." << endl;
    outputFile->cd();

    outputTree->Write();
    outputFile->Close();
    delete outputFile;
}
