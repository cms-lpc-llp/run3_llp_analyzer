#include <iterator>
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TKey.h"
#include <assert.h>
#include <iostream>
#include <vector>
#include <regex>
#include <sys/stat.h>

#include "nano_events.h"
#include "llp_event.h"

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

int main(int argc, char **argv)
{
  using namespace std;

  if (argc < 3)
  {
    cerr << "usage MergeNtuple [inp] [out_dir]" << endl;
    return -1;
  }

  string path_inp(argv[1]);
  string path_out(argv[2]);

  struct stat buffer;
  if (stat(path_out.c_str(), &buffer) != 0)
    mkdir(path_out.c_str(), 0750);

  std::regex e("^.*/(.*)/(.*)/(.*)/(.*\\.root)$");
  std::string result = std::regex_replace(path_inp, e, "$1_$2_$3_$4");
  path_out = path_out + "/" + result;

  cout << "cache path: " << path_out << endl;

  TFile *inp_file = TFile::Open(path_inp.c_str(), "READ");
  TFile *out_file = new TFile(path_out.c_str(), "RECREATE", "", 101);

  uint32_t run;
  ULong64_t event;
  const char *run_key = "runNum";
  const char *event_key = "eventNum";

  auto _tree = inp_file->Get("ntuples/llp");
  // ntuple
  if (_tree == nullptr)
    // ntuple alt
    _tree = inp_file->Get("llp");
  if (_tree == nullptr)
  {
    // AOD
    _tree = inp_file->Get("Events");
    event_key = "event";
    run_key = "run";
  }

  TTree *tree = (TTree *)_tree;

  tree->SetBranchAddress(run_key, &run);
  tree->SetBranchAddress(event_key, &event);

  out_file->cd();
  auto out_tree = new TTree("events", "events");
  out_tree->Branch("run", &run);
  out_tree->Branch("event", &event);
  out_tree->SetDirectory(out_file);

  auto len = tree->GetEntries();
  cout << "# events: " << len << endl;
  for (int i = 0; i < len; i++)
  {
    tree->GetEntry(i);
    if (i % 1000 == 0)
      cout << "Processing event " << i << endl;
    out_tree->Fill();
  }

  out_tree->Write();
  inp_file->Close();
  out_file->Close();
}