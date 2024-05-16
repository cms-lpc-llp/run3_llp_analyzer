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
#include <regex>
#include "nano_events.h"
#include "llp_event.h"
#include "TreeMuonSystemCombination.h"
#include "RazorHelper.h"

#include <vector>
#include "TString.h"

#define N_MAX_RECHITS 20000
#define N_MAX_SEGMENT 1000

using namespace std;

string get_cache_path(string path_inp, string dir_out)
{
  regex e("^.*/(.*)/(.*)/(.*)/(.*\\.root)$");
  string result = regex_replace(path_inp, e, "$1_$2_$3_$4");
  string path_out = dir_out + "/" + result;
  return path_out;
}

map<pair<uint, ulong>, ulong> get_event_index_map(TChain *chain)
{
  map<pair<uint, ulong>, ulong> eventIdxMap;
  uint32_t run;
  ULong64_t event;
  chain->SetBranchAddress("run", &run);
  chain->SetBranchAddress("event", &event);
  for (UInt_t m = 0; m < chain->GetEntries(); m++)
  {
    chain->GetEntry(m);
    eventIdxMap[make_pair(run, event)] = m;
  }
  return eventIdxMap;
}

// get list of files to open, add normalization branch to the tree in each file
int main(int argc, char *argv[])
{
  // parse input list to get names of ROOT files
  if (argc < 4)
  {
    cerr << "usage MergeNtuple [inp_ntuple_list] [inp_nano_list] [output_file] [cache_dir]" << endl;
    return -1;
  }
  // string filenameNTuplers(argv[1]);
  string filenameNTuplerList(argv[1]);
  string filenameNanoAODList(argv[2]);
  string outputFileName(argv[3]);
  string cacheDir(argv[4]);

  ////////////////////////
  ////// Load ntuples
  ////////////////////////
  TChain *ntupleChain = new TChain();
  TChain *ntupleCacheChain = new TChain();

  string curNtupleName;
  string curNtupleCacheName;
  cout << "filenameNTuplerList.c_str() = " << filenameNTuplerList.c_str() << endl;
  ifstream inputNtuple(filenameNTuplerList.c_str());
  int NtuplesLoaded = 0;
  if (!inputNtuple)
  {
    cerr << "Error: input ntuple file list not found!" << endl;
    return -1;
  }
  while (getline(inputNtuple, curNtupleName))
  {
    curNtupleCacheName = get_cache_path(curNtupleName, cacheDir);
    cout << "curFileName = " << curNtupleName << endl;
    if (NtuplesLoaded == 0)
    {
      // checks root file structure and add first file
      cout << "[INFO]: loading file: " << curNtupleName.c_str() << endl;
      TFile *f_0 = TFile::Open(curNtupleName.c_str(), "READ");
      if (f_0->GetDirectory("ntuples"))
      {
        ntupleChain->SetName("ntuples/llp");
        cout << "[INFO]: default configuration for tchain" << endl;
      }
      else
      {
        ntupleChain->SetName("llp");
        cout << "[INFO]: alternative configuration for tchain" << endl;
      }
      ntupleChain->Add(curNtupleName.c_str());
      ntupleCacheChain->SetName("events");
      ntupleCacheChain->Add(curNtupleCacheName.c_str());
      delete f_0;
    }
    else
    {
      // Addind remaining files after file structure is decided
      ntupleChain->Add(curNtupleName.c_str());
      ntupleCacheChain->Add(curNtupleCacheName.c_str());
    }
    NtuplesLoaded++;
  }
  cout << "Loaded ntuples: " << NtuplesLoaded << " files and " << ntupleChain->GetEntries() << " events" << endl;
  if (ntupleChain == NULL)
    return -1;

  ////////////////////////
  ////// Load nanoAOD files
  ////////////////////////
  TChain *nanoChain = new TChain();
  TChain *nanoCacheChain = new TChain();

  string curFileName;
  string curCacheName;
  ifstream inputFile(filenameNanoAODList.c_str());
  int NFilesLoaded = 0;
  if (!inputFile)
  {
    cerr << "Error: input nanoAOD file list not found!" << endl;
    return -1;
  }

  while (getline(inputFile, curFileName))
  {
    curCacheName = get_cache_path(curFileName, cacheDir);
    cout << "curFileName = " << curFileName << endl;
    if (NFilesLoaded == 0)
    {
      // checks root file structure and add first file
      cout << "[INFO]: loading file: " << curFileName.c_str() << endl;
      TFile *f_0 = TFile::Open(curFileName.c_str(), "READ");
      f_0->ls();
      nanoChain->SetName("Events");
      nanoChain->Add(curFileName.c_str());
      delete f_0;
      nanoCacheChain->SetName("events");
      nanoCacheChain->Add(curCacheName.c_str());
    }
    else
    {
      // Addind remaining files after file structure is decided
      nanoChain->Add(curFileName.c_str());
      nanoCacheChain->Add(curCacheName.c_str());
    }
    NFilesLoaded++;
  }
  cout << "Loaded nanoAOD: " << NFilesLoaded << " files and " << nanoChain->GetEntries() << " events" << endl;
  if (nanoChain == NULL)
    return -1;

  //*****************************************************************************************
  // Make map of event number in tree 1 to event number in tree 2
  //*****************************************************************************************

  vector<pair<ulong, ulong>> indexPair; // nTuple - AOD idx pair
  nano_events nano = nano_events(nanoChain);
  llp_event ntuple = llp_event(ntupleChain);

  // loop over tree2
  cout << "building AOD tree index map" << endl;
  auto nTupleIdxMap = get_event_index_map(ntupleCacheChain);

  cout << "building NTuple tree index map" << endl;
  auto AODIdxMap = get_event_index_map(nanoCacheChain);

  // find intersection of the two maps
  for (auto it = nTupleIdxMap.begin(); it != nTupleIdxMap.end(); it++)
  {
    if (AODIdxMap.count(it->first) > 0)
    {
      indexPair.push_back(make_pair(it->second, AODIdxMap[it->first]));
    }
  }

  cout << "Sorting ... " << endl;

  sort(indexPair.begin(), indexPair.end(), [](const pair<ulong, ulong> &left, const pair<ulong, ulong> &right)
       { return left.second < right.second; });

  ulong nMatchedEvents = indexPair.size();
  ulong NEventsTreeNTuple = ntupleChain->GetEntries();

  cout << "Matched events    = " << nMatchedEvents << " / " << NEventsTreeNTuple << " \n";
  cout << "Un-Matched events = " << NEventsTreeNTuple - nMatchedEvents << " / " << NEventsTreeNTuple << " \n";
  if (nMatchedEvents == 0)
  {
    cout << "No matched events found. Exiting..." << endl;
    return 0;
  }
  //*****************************************************************************************
  // Produce Output Tree
  //*****************************************************************************************

  // clone tree with 0 entries, copy all the branch addresses only

  TFile *outputFile = new TFile(outputFileName.c_str(), "RECREATE");
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

  for (int i = 0; i < numFloatBranches; i++)
  {
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

  for (int i = 0; i < numFloatSegBranches; i++)
  {
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

  for (int i = 0; i < numIntBranches; i++)
  {
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

  for (int i = 0; i < numIntSegBranches; i++)
  {
    cout << "Adding Branch: " << addBranchNamesIntSeg[i] << "\n";
    if (string(addBranchNamesIntSeg[i]).find("cscSeg") != string::npos)
      outputTree->Branch(addBranchNamesIntSeg[i], addBranchesInputVarIntSeg[i], TString::Format("%s[nCscSeg]/I", addBranchNamesIntSeg[i]));
    else if (string(addBranchNamesIntSeg[i]).find("dtSeg") != string::npos)
      outputTree->Branch(addBranchNamesIntSeg[i], addBranchesInputVarIntSeg[i], TString::Format("%s[nDtSeg]/I", addBranchNamesIntSeg[i]));
    else if (string(addBranchNamesIntSeg[i]).find("rpc") != string::npos)
      outputTree->Branch(addBranchNamesIntSeg[i], addBranchesInputVarIntSeg[i], TString::Format("%s[nRpc]/I", addBranchNamesIntSeg[i]));
    outputTree->SetBranchAddress(addBranchNamesIntSeg[i], addBranchesInputVarIntSeg[i]);
  }

  /////////////////////////////////////////////////////
  /////// FILL OUTPUT TREE
  /////////////////////////////////////////////////////

  for (auto it : indexPair)
  {
    auto n = it.first;
    auto m = it.second;
    if (n % 10000 == 0)
      cout << "Processing entry " << n << "\n";

    nanoChain->GetEntry(m);
    ntupleChain->GetEntry(n);

    // fill the new added branches

    nCscRechits = ntuple.ncscRechits;
    nCscSeg = ntuple.nCscSeg;
    nDtRechits = ntuple.nDtRechits;
    nDtSeg = ntuple.nDtSeg;
    nRpc = ntuple.nRpc;

    for (int i = 0; i < numFloatBranches; i++)
    {
      int temp_nhits = 0;
      if (string(addBranchNamesFloat[i]).find("cscRechits") != string::npos)
        temp_nhits = nCscRechits;
      if (string(addBranchNamesFloat[i]).find("dtRechit") != string::npos)
        temp_nhits = nDtRechits;
      for (int j = 0; j < temp_nhits; j++)
        addBranchesInputVarFloat[i][j] = addBranchesRazorVarFloat[i][j];
    }

    for (int i = 0; i < numFloatSegBranches; i++)
    {
      int temp_nhits = 0;
      if (string(addBranchNamesFloatSeg[i]).find("cscSeg") != string::npos)
        temp_nhits = nCscSeg;
      if (string(addBranchNamesFloatSeg[i]).find("dtSeg") != string::npos)
        temp_nhits = nDtSeg;
      if (string(addBranchNamesFloatSeg[i]).find("rpc") != string::npos)
        temp_nhits = nRpc;
      for (int j = 0; j < temp_nhits; j++)
        addBranchesInputVarFloatSeg[i][j] = addBranchesRazorVarFloatSeg[i][j];
    }

    for (int i = 0; i < numIntBranches; i++)
    {
      int temp_nhits = 0;
      if (string(addBranchNamesInt[i]).find("cscRechits") != string::npos)
        temp_nhits = nCscRechits;
      if (string(addBranchNamesInt[i]).find("dtRechit") != string::npos)
        temp_nhits = nDtRechits;
      for (int j = 0; j < temp_nhits; j++)
        addBranchesInputVarInt[i][j] = addBranchesRazorVarInt[i][j];
    }

    for (int i = 0; i < numIntSegBranches; i++)
    {
      int temp_nhits = 0;
      if (string(addBranchNamesIntSeg[i]).find("cscSeg") != string::npos)
        temp_nhits = nCscSeg;
      if (string(addBranchNamesIntSeg[i]).find("dtSeg") != string::npos)
        temp_nhits = nDtSeg;
      if (string(addBranchNamesIntSeg[i]).find("rpc") != string::npos)
        temp_nhits = nRpc;
      for (int j = 0; j < temp_nhits; j++)
        addBranchesInputVarIntSeg[i][j] = addBranchesRazorVarIntSeg[i][j];
    }

    outputTree->Fill();
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