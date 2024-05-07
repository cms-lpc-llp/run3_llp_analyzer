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
#include "nano_events.h"
#include "llp_event.h"
#include "TreeMuonSystemCombination.h"
#include "RazorHelper.h"

#include <vector>
#include "TString.h"
using namespace std;


#define N_MAX_RECHITS 20000
#define N_MAX_SEGMENT 1000



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



// get list of files to open, add normalization branch to the tree in each file
int main(int argc, char *argv[])
{
  // parse input list to get names of ROOT files
  if (argc < 4)
  {
    cerr << "usage MergeNtuple [inputfile1] [inputfile2] [outputfile] [isData]" << endl;
    return -1;
  }
  // string filenameNTuplers(argv[1]);
  string filenameNTuplerList(argv[1]);
  string filenameNanoAODList(argv[2]);
  string outputfilename(argv[3]);
  string analysisTag(argv[4]); // analysisTag

  if (analysisTag == "")
  {
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
  if (!inputNtuple)
  {
    cerr << "Error: input ntuple file list not found!" << endl;
    return -1;
  }
  while (getline(inputNtuple, curNtupleName))
  {
    cout << "curFileName = " << curNtupleName << endl;
    if (NtuplesLoaded == 0)
    {
      // checks root file structure and add first file
      std::cout << "[INFO]: loading file: " << curNtupleName.c_str() << std::endl;
      TFile *f_0 = TFile::Open(curNtupleName.c_str());
      if( f_0->GetDirectory("ntuples") )
      {
        ntupleChain->SetName("ntuples/llp");
        std::cout << "[INFO]: default configuration for tchain" << std::endl;
      }
      else
      {
        ntupleChain->SetName("llp");
        std::cout << "[INFO]: alternative configuration for tchain"<< std::endl;
      }
      ntupleChain->Add(curNtupleName.c_str());
      delete f_0;

    }
    else 
    {
      // Addind remaining files after file structure is decided
      ntupleChain->Add(curNtupleName.c_str());
    }
    NtuplesLoaded++;
  }
  std::cout << "Loaded ntuples: " << NtuplesLoaded << " files and " << ntupleChain->GetEntries() << " events" <<endl;
  if ( ntupleChain == NULL ) return -1;


  ////////////////////////
  ////// Load nanoAOD files
  ////////////////////////
  TChain *nanoChain = new TChain();
  string curFileName;
  ifstream inputFile(filenameNanoAODList.c_str());
  int NFilesLoaded = 0;
  if (!inputFile)
  {
    cerr << "Error: input nanoAOD file list not found!" << endl;
    return -1;
  }

  while (getline(inputFile, curFileName))
  {
    cout << "curFileName = " << curFileName << endl;
    if (NFilesLoaded == 0)
    {
      // checks root file structure and add first file
      std::cout << "[INFO]: loading file: " << curFileName.c_str() << std::endl;
      TFile *f_0 = TFile::Open(curFileName.c_str());
      f_0->ls();
      nanoChain->SetName("Events");
      nanoChain->Add(curFileName.c_str());
      delete f_0;
    }
    else 
    {
      // Addind remaining files after file structure is decided
      nanoChain->Add(curFileName.c_str());
    }
    NFilesLoaded++;
  }
  std::cout << "Loaded nanoAOD: " << NFilesLoaded << " files and " << nanoChain->GetEntries() << " events" <<endl;
  if (nanoChain == NULL) return -1;
    


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

  cout << "building nanoAOD tree index map" <<endl;
  for (UInt_t m = 0; m < NEventsTree2; m++)
  {
    if (m % 10000 == 0) cout << "Processing entry " << m << endl;
    nano.fChain->GetEntry(m);
    std::pair<uint, ulong> p(nano.run, nano.event);
    eventList2.push_back(p);
  }

  cout << "Looping over ntuple tree to find matched events" <<endl;;
  // loop over tree1
  for (uint n = 0; n < NEventsTree1; n++)
  {
    ntupleChain->GetEntry(n);
    if (n % 10000 == 0) cout << "Event " << n << "\n";

    bool matchFound = false;
    for (uint m = 0; m < eventList2.size(); m++)
    {
      if (eventList2[m].first != ntuple.runNum)continue;
      if (eventList2[m].second != ntuple.eventNum) continue;
      EventIndexToEventIndexMap[n] = m;
      matchFound = true;
      nMatchedEvents++;
    }

    if (matchFound)matchedevent.push_back(true);
    else matchedevent.push_back(false);
  }

  cout << "Matched events    = " << nMatchedEvents << " / " << NEventsTree1 << " \n";
  cout << "Un-Matched events = " << NEventsTree1-nMatchedEvents << " / " << NEventsTree1 << " \n";
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

  outputTree->SetBranchAddress("nCscSeg",&nCscSeg);
  outputTree->SetBranchAddress("nCscRechits",&nCscRechits);
  outputTree->SetBranchAddress("nDtSeg",&nDtSeg);
  outputTree->SetBranchAddress("nDtRechits",&nDtRechits);
  outputTree->SetBranchAddress("nRpc",&nRpc);
  ////////////////////////////////////////////////
  ///// Float for rechits related variables
  ////////////////////////////////////////////////

  int numFloatBranches = 12;

  const char *addBranchNamesFloat[numFloatBranches]{
    "cscRechitsPhi","cscRechitsEta","cscRechitsX","cscRechitsY","cscRechitsZ","cscRechitsTpeak","cscRechitsTwire",
    "dtRechitCorrectX","dtRechitCorrectY","dtRechitCorrectZ","dtRechitCorrectEta","dtRechitCorrectPhi",
  };
  Float_t *addBranchesRazorVarFloat[numFloatBranches]{
    ntuple.cscRechitsPhi,ntuple.cscRechitsEta,ntuple.cscRechitsX,ntuple.cscRechitsY,ntuple.cscRechitsZ,ntuple.cscRechitsTpeak,ntuple.cscRechitsTwire,
    ntuple.dtRechitCorrectX,ntuple.dtRechitCorrectY,ntuple.dtRechitCorrectZ,ntuple.dtRechitCorrectEta,ntuple.dtRechitCorrectPhi,
  };
  Float_t addBranchesInputVarFloat[numFloatBranches][N_MAX_RECHITS];

  for (int i = 0; i < numFloatBranches; i++)
  {
    cout << "Adding Branch: " << addBranchNamesFloat[i] << "\n";
    if (string(addBranchNamesFloat[i]).find("cscRechits") != std::string::npos) outputTree->Branch(addBranchNamesFloat[i], addBranchesInputVarFloat[i], TString::Format("%s[nCscRechits]/F", addBranchNamesFloat[i]));
    if (string(addBranchNamesFloat[i]).find("dtRechit") != std::string::npos) outputTree->Branch(addBranchNamesFloat[i], addBranchesInputVarFloat[i], TString::Format("%s[nDtRechits]/F", addBranchNamesFloat[i]));
    outputTree->SetBranchAddress(addBranchNamesFloat[i],addBranchesInputVarFloat[i]);

  }


  ////////////////////////////////////////////////
  ///// Float for segment related variables
  ////////////////////////////////////////////////

  int numFloatSegBranches = 11;

  const char *addBranchNamesFloatSeg[numFloatSegBranches]{
  "cscSegPhi", "cscSegEta", 
  "dtSegPhi","dtSegEta", 
    "rpcPhi","rpcEta","rpcX","rpcY","rpcZ","rpcT","rpcTError"
  };
  Float_t *addBranchesRazorVarFloatSeg[numFloatSegBranches]{
  ntuple.cscSegPhi,ntuple.cscSegEta,
  ntuple.dtSegPhi,ntuple.dtSegEta,
    ntuple.rpcPhi,ntuple.rpcEta,ntuple.rpcX,ntuple.rpcY,ntuple.rpcZ, ntuple.rpcT, ntuple.rpcTError
  };
  Float_t addBranchesInputVarFloatSeg[numFloatSegBranches][N_MAX_SEGMENT];

  for (int i = 0; i < numFloatSegBranches; i++)
  {
    cout << "Adding Branch: " << addBranchNamesFloatSeg[i] << "\n";
    if (string(addBranchNamesFloatSeg[i]).find("cscSeg") != std::string::npos) outputTree->Branch(addBranchNamesFloatSeg[i], addBranchesInputVarFloatSeg[i], TString::Format("%s[nCscSeg]/F", addBranchNamesFloatSeg[i]));
    if (string(addBranchNamesFloatSeg[i]).find("dtSeg") != std::string::npos) outputTree->Branch(addBranchNamesFloatSeg[i], addBranchesInputVarFloatSeg[i], TString::Format("%s[nDtSeg]/F", addBranchNamesFloatSeg[i]));
    if (string(addBranchNamesFloatSeg[i]).find("rpc") != std::string::npos) outputTree->Branch(addBranchNamesFloatSeg[i], addBranchesInputVarFloatSeg[i], TString::Format("%s[nRpc]/F", addBranchNamesFloatSeg[i]));
    outputTree->SetBranchAddress(addBranchNamesFloatSeg[i],addBranchesInputVarFloatSeg[i]);

  }
  ////////////////////////////////////////////////
  /////integer for rechits related variables
  ////////////////////////////////////////////////

  int numIntBranches = 6;
  const char *addBranchNamesInt[numIntBranches]{ "cscRechitsChamber","cscRechitsStation","cscRechitsDetId","dtRechitStation","dtRechitWheel","dtRechitSuperLayer"};
    
  Int_t *addBranchesRazorVarInt[numIntBranches]{
    ntuple.cscRechitsChamber,ntuple.cscRechitsStation,ntuple.cscRechitsDetId,
    ntuple.dtRechitStation,ntuple.dtRechitWheel,ntuple.dtRechitSuperLayer};
  
  Int_t addBranchesInputVarInt[numIntBranches][N_MAX_RECHITS];

  for (int i = 0; i < numIntBranches; i++)
  {
    cout << "Adding Branch: " << addBranchNamesInt[i] << "\n";
    if (string(addBranchNamesInt[i]).find("cscRechits") != std::string::npos) outputTree->Branch(addBranchNamesInt[i], addBranchesInputVarInt[i], TString::Format("%s[nCscRechits]/I", addBranchNamesInt[i]));
    else if (string(addBranchNamesInt[i]).find("dtRechit") != std::string::npos) outputTree->Branch(addBranchNamesInt[i], addBranchesInputVarInt[i], TString::Format("%s[nDtRechits]/I", addBranchNamesInt[i]));
    outputTree->SetBranchAddress(addBranchNamesInt[i],addBranchesInputVarInt[i]);
  }

  // ////////////////////////////////////////////////
  // /////integer for segments related variables
  // ////////////////////////////////////////////////

  int numIntSegBranches = 11;
  const char *addBranchNamesIntSeg[numIntSegBranches]{
    "cscSegChamber","cscSegStation","cscSegNRecHits","dtSegStation","dtSegWheel",
    "rpcBx","rpcRegion","rpcRing","rpcSector","rpcStation","rpcLayer"
  };

  Int_t *addBranchesRazorVarIntSeg[numIntSegBranches]{
    ntuple.cscSegChamber,ntuple.cscSegStation,ntuple.cscSegNRecHits,ntuple.dtSegStation,ntuple.dtSegWheel,
    ntuple.rpcBx,ntuple.rpcRegion,ntuple.rpcRing,ntuple.rpcSector,ntuple.rpcStation, ntuple.rpcLayer
  
  };
  Int_t addBranchesInputVarIntSeg[numIntSegBranches][N_MAX_SEGMENT];

  for (int i = 0; i < numIntSegBranches; i++)
  {
    cout << "Adding Branch: " << addBranchNamesIntSeg[i] << "\n";
    if (string(addBranchNamesIntSeg[i]).find("cscSeg") != std::string::npos) outputTree->Branch(addBranchNamesIntSeg[i], addBranchesInputVarIntSeg[i], TString::Format("%s[nCscSeg]/I", addBranchNamesIntSeg[i]));
    else if (string(addBranchNamesIntSeg[i]).find("dtSeg") != std::string::npos) outputTree->Branch(addBranchNamesIntSeg[i], addBranchesInputVarIntSeg[i], TString::Format("%s[nDtSeg]/I", addBranchNamesIntSeg[i]));
    else if (string(addBranchNamesIntSeg[i]).find("rpc") != std::string::npos) outputTree->Branch(addBranchNamesIntSeg[i], addBranchesInputVarIntSeg[i], TString::Format("%s[nRpc]/I", addBranchNamesIntSeg[i]));
    outputTree->SetBranchAddress(addBranchNamesIntSeg[i],addBranchesInputVarIntSeg[i]);

  }


  /////////////////////////////////////////////////////
  /////// FILL OUTPUT TREE
  /////////////////////////////////////////////////////

  for (uint n = 0; n < NEventsTree1; n++)
  {
    if (n % 10000 == 0) cout << "Processing entry " << n << "\n";

    // Check if found a match
    if (!matchedevent[n])  continue; 
    nanoChain->GetEntry(EventIndexToEventIndexMap[n]); 
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
      if (string(addBranchNamesFloat[i]).find("cscRechits") != std::string::npos) temp_nhits = nCscRechits;
      if (string(addBranchNamesFloat[i]).find("dtRechit") != std::string::npos) temp_nhits = nDtRechits;
      for (int j = 0; j < temp_nhits; j++) addBranchesInputVarFloat[i][j] = addBranchesRazorVarFloat[i][j];

    }

    for (int i = 0; i < numFloatSegBranches; i++)
    {
      int temp_nhits = 0;
      if (string(addBranchNamesFloatSeg[i]).find("cscSeg") != std::string::npos) temp_nhits = nCscSeg;
      if (string(addBranchNamesFloatSeg[i]).find("dtSeg") != std::string::npos) temp_nhits = nDtSeg;
      if (string(addBranchNamesFloatSeg[i]).find("rpc") != std::string::npos) temp_nhits = nRpc;
      for (int j = 0; j < temp_nhits; j++) addBranchesInputVarFloatSeg[i][j] = addBranchesRazorVarFloatSeg[i][j];

    }

   for (int i = 0; i < numIntBranches; i++)
    {
      int temp_nhits = 0;
      if (string(addBranchNamesInt[i]).find("cscRechits") != std::string::npos) temp_nhits = nCscRechits;
      if (string(addBranchNamesInt[i]).find("dtRechit") != std::string::npos) temp_nhits = nDtRechits;
      for (int j = 0; j < temp_nhits; j++) addBranchesInputVarInt[i][j] = addBranchesRazorVarInt[i][j];

    }

    for (int i = 0; i < numIntSegBranches; i++)
    {
      int temp_nhits = 0;
      if (string(addBranchNamesIntSeg[i]).find("cscSeg") != std::string::npos) temp_nhits = nCscSeg;
      if (string(addBranchNamesIntSeg[i]).find("dtSeg") != std::string::npos) temp_nhits = nDtSeg;
      if (string(addBranchNamesIntSeg[i]).find("rpc") != std::string::npos) temp_nhits = nRpc;
      for (int j = 0; j < temp_nhits; j++) addBranchesInputVarIntSeg[i][j] = addBranchesRazorVarIntSeg[i][j];
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