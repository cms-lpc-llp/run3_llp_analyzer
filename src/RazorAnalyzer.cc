#include "RazorAnalyzer.h"
#include "TLorentzVector.h"

#include <cassert>
#include <cmath>

using namespace std;

RazorAnalyzerMerged::RazorAnalyzerMerged(TTree* tree)
    : mdsnano_event(tree) {
  // turn off all branches
  //  fChain->SetBranchStatus("*", 1);
  fChain->SetBranchStatus("*", 0);
}

RazorAnalyzerMerged::~RazorAnalyzerMerged() {
}

void RazorAnalyzerMerged::EnableAll() {
  fChain->SetBranchStatus("*", 1);
}

void RazorAnalyzerMerged::Analyze(bool isData, int option, string outputFileName, string label) {
  std::cerr << "Virtual RazorAnalyzerMerged::Analyze" << std::endl;
  assert(0);
}

double RazorAnalyzerMerged::deltaPhi(double phi1, double phi2) {
  double dphi = phi1 - phi2;
  while (dphi > TMath::Pi())
    dphi -= TMath::TwoPi();
  while (dphi <= -TMath::Pi())
    dphi += TMath::TwoPi();
  return dphi;
}

double RazorAnalyzerMerged::deltaR(double eta1, double phi1, double eta2, double phi2) {
  double dphi = deltaPhi(phi1, phi2);
  double deta = eta1 - eta2;
  return std::sqrt(dphi * dphi + deta * deta);
}

TLorentzVector RazorAnalyzerMerged::makeTLorentzVector(double pt, double eta, double phi, double energy) {
  TLorentzVector vec;
  vec.SetPtEtaPhiE(pt, eta, phi, energy);
  return vec;
}

TLorentzVector RazorAnalyzerMerged::makeTLorentzVectorPtEtaPhiM(double pt, double eta, double phi, double mass) {
  TLorentzVector vec;
  vec.SetPtEtaPhiM(pt, eta, phi, mass);
  return vec;
}
