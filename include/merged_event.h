//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Thu May 23 17:59:07 2024 by ROOT version 6.30/06
// from TTree Events/Events
// found on file: 04bdd0bc-916b-4779-9bbd-b74e001e5423.root
//////////////////////////////////////////////////////////

#ifndef merged_event_h
#define merged_event_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

class merged_event
{
public:
   TTree *fChain;  //! pointer to the analyzed TTree or TChain
   Int_t fCurrent; //! current Tree number in a TChain

   // Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   // ===================== LLP Ntuple =====================
   Float_t gLLP_eta[2];
   Float_t gLLP_phi[2];
   Float_t gLLP_csc[2];
   Float_t gLLP_dt[2];
   Float_t gLLP_beta[2];
   Float_t gLLP_e[2];
   Float_t gLLP_pt[2];
   Float_t gLLP_decay_vertex_x[2];
   Float_t gLLP_decay_vertex_y[2];
   Float_t gLLP_decay_vertex_z[2];
   Float_t gLLP_decay_vertex_r[2];

   // ===================== MC =====================

   Int_t nsoftActivityVH;
   Int_t nGenJetAK8;
   Float_t GenJetAK8_eta[6 * 10];  //[nGenJetAK8]
   Float_t GenJetAK8_mass[6 * 10]; //[nGenJetAK8]
   Float_t GenJetAK8_phi[6 * 10];  //[nGenJetAK8]
   Float_t GenJetAK8_pt[6 * 10];   //[nGenJetAK8]
   Int_t nGenJet;
   Float_t GenJet_eta[19 * 10];  //[nGenJet]
   Float_t GenJet_mass[19 * 10]; //[nGenJet]
   Float_t GenJet_phi[19 * 10];  //[nGenJet]
   Float_t GenJet_pt[19 * 10];   //[nGenJet]
   Int_t nGenPart;
   Short_t GenPart_genPartIdxMother[125 * 10]; //[nGenPart]
   UShort_t GenPart_statusFlags[125 * 10];     //[nGenPart]
   Int_t GenPart_pdgId[125 * 10];              //[nGenPart]
   Int_t GenPart_status[125 * 10];             //[nGenPart]
   Float_t GenPart_eta[125 * 10];              //[nGenPart]
   Float_t GenPart_mass[125 * 10];             //[nGenPart]
   Float_t GenPart_phi[125 * 10];              //[nGenPart]
   Float_t GenPart_pt[125 * 10];               //[nGenPart]
   Int_t nGenProton;
   Bool_t GenProton_isPU[11 * 10]; //[nGenProton]
   Float_t GenProton_px[11 * 10];  //[nGenProton]
   Float_t GenProton_py[11 * 10];  //[nGenProton]
   Float_t GenProton_pz[11 * 10];  //[nGenProton]
   Float_t GenProton_vz[11 * 10];  //[nGenProton]
   Int_t nSubGenJetAK8;
   Float_t SubGenJetAK8_eta[12 * 10];  //[nSubGenJetAK8]
   Float_t SubGenJetAK8_mass[12 * 10]; //[nSubGenJetAK8]
   Float_t SubGenJetAK8_phi[12 * 10];  //[nSubGenJetAK8]
   Float_t SubGenJetAK8_pt[12 * 10];   //[nSubGenJetAK8]
   Int_t Generator_id1;
   Int_t Generator_id2;
   Float_t Generator_binvar;
   Float_t Generator_scalePDF;
   Float_t Generator_weight = -1;
   Float_t Generator_x1;
   Float_t Generator_x2;
   Float_t Generator_xpdf1;
   Float_t Generator_xpdf2;
   Float_t GenVtx_x;
   Float_t GenVtx_y;
   Float_t GenVtx_z;
   Int_t nGenVisTau;
   UChar_t GenVisTau_status[2 * 10];           //[nGenVisTau]
   Short_t GenVisTau_charge[2 * 10];           //[nGenVisTau]
   Short_t GenVisTau_genPartIdxMother[2 * 10]; //[nGenVisTau]
   Float_t GenVisTau_eta[2 * 10];              //[nGenVisTau]
   Float_t GenVisTau_mass[2 * 10];             //[nGenVisTau]
   Float_t GenVisTau_phi[2 * 10];              //[nGenVisTau]
   Float_t GenVisTau_pt[2 * 10];               //[nGenVisTau]
   Float_t genWeight;
   Float_t GenMET_phi;
   Float_t GenMET_pt;
   Int_t nGenDressedLepton;
   Bool_t GenDressedLepton_hasTauAnc[1 * 10]; //[nGenDressedLepton]
   Int_t GenDressedLepton_pdgId[1 * 10];      //[nGenDressedLepton]
   Float_t GenDressedLepton_eta[1 * 10];      //[nGenDressedLepton]
   Float_t GenDressedLepton_mass[1 * 10];     //[nGenDressedLepton]
   Float_t GenDressedLepton_phi[1 * 10];      //[nGenDressedLepton]
   Float_t GenDressedLepton_pt[1 * 10];       //[nGenDressedLepton]
   Int_t nGenIsolatedPhoton;
   Float_t GenIsolatedPhoton_eta[2 * 10];  //[nGenIsolatedPhoton]
   Float_t GenIsolatedPhoton_mass[2 * 10]; //[nGenIsolatedPhoton]
   Float_t GenIsolatedPhoton_phi[2 * 10];  //[nGenIsolatedPhoton]
   Float_t GenIsolatedPhoton_pt[2 * 10];   //[nGenIsolatedPhoton]
   Int_t genTtbarId;
   UChar_t boostedTau_genPartFlav[2 * 10];  //[nboostedTau]
   Short_t boostedTau_genPartIdx[2 * 10];   //[nboostedTau]
   UChar_t Electron_genPartFlav[5 * 10];    //[nElectron]
   Short_t Electron_genPartIdx[5 * 10];     //[nElectron]
   Short_t FatJet_genJetAK8Idx[4 * 10];     //[nFatJet]
   UChar_t GenJetAK8_hadronFlavour[6 * 10]; //[nGenJetAK8]
   Short_t GenJetAK8_partonFlavour[6 * 10]; //[nGenJetAK8]
   UChar_t GenJet_hadronFlavour[19 * 10];   //[nGenJet]
   Short_t GenJet_partonFlavour[19 * 10];   //[nGenJet]
   Float_t GenVtx_t0;
   Short_t Jet_genJetIdx[16 * 10];            //[nJet]
   UChar_t LowPtElectron_genPartFlav[6 * 10]; //[nLowPtElectron]
   Short_t LowPtElectron_genPartIdx[6 * 10];  //[nLowPtElectron]
   UChar_t Muon_genPartFlav[8 * 10];          //[nMuon]
   Short_t Muon_genPartIdx[8 * 10];           //[nMuon]
   UChar_t Photon_genPartFlav[6 * 10];        //[nPhoton]
   Short_t Photon_genPartIdx[6 * 10];         //[nPhoton]
   Float_t MET_fiducialGenPhi;
   Float_t MET_fiducialGenPt;
   UChar_t Tau_genPartFlav[5 * 10]; //[nTau]
   Short_t Tau_genPartIdx[5 * 10];  //[nTau]
   // ===================== Data =====================

   UInt_t run;
   UInt_t luminosityBlock;
   ULong64_t event;
   UInt_t bunchCrossing;
   Char_t BeamSpot_type;
   Float_t BeamSpot_sigmaZ;
   Float_t BeamSpot_sigmaZError;
   Float_t BeamSpot_z;
   Float_t BeamSpot_zError;
   Int_t nboostedTau;
   UChar_t boostedTau_idAntiEle2018[4];            //[nboostedTau]
   UChar_t boostedTau_idAntiMu[4];                 //[nboostedTau]
   UChar_t boostedTau_idMVAnewDM2017v2[4];         //[nboostedTau]
   UChar_t boostedTau_idMVAoldDM2017v2[4];         //[nboostedTau]
   Short_t boostedTau_jetIdx[4];                   //[nboostedTau]
   Short_t boostedTau_rawAntiEleCat2018[4];        //[nboostedTau]
   Int_t boostedTau_charge[4];                     //[nboostedTau]
   Int_t boostedTau_decayMode[4];                  //[nboostedTau]
   Float_t boostedTau_chargedIso[4];               //[nboostedTau]
   Float_t boostedTau_eta[4];                      //[nboostedTau]
   Float_t boostedTau_leadTkDeltaEta[4];           //[nboostedTau]
   Float_t boostedTau_leadTkDeltaPhi[4];           //[nboostedTau]
   Float_t boostedTau_leadTkPtOverTauPt[4];        //[nboostedTau]
   Float_t boostedTau_mass[4];                     //[nboostedTau]
   Float_t boostedTau_neutralIso[4];               //[nboostedTau]
   Float_t boostedTau_phi[4];                      //[nboostedTau]
   Float_t boostedTau_photonsOutsideSignalCone[4]; //[nboostedTau]
   Float_t boostedTau_pt[4];                       //[nboostedTau]
   Float_t boostedTau_puCorr[4];                   //[nboostedTau]
   Float_t boostedTau_rawAntiEle2018[4];           //[nboostedTau]
   Float_t boostedTau_rawIso[4];                   //[nboostedTau]
   Float_t boostedTau_rawIsodR03[4];               //[nboostedTau]
   Float_t boostedTau_rawMVAnewDM2017v2[4];        //[nboostedTau]
   Float_t boostedTau_rawMVAoldDM2017v2[4];        //[nboostedTau]
   Float_t CaloMET_phi;
   Float_t CaloMET_pt;
   Float_t CaloMET_sumEt;
   Float_t ChsMET_phi;
   Float_t ChsMET_pt;
   Float_t ChsMET_sumEt;
   Int_t nCorrT1METJet;
   Float_t CorrT1METJet_area[50];            //[nCorrT1METJet]
   Float_t CorrT1METJet_eta[50];             //[nCorrT1METJet]
   Float_t CorrT1METJet_muonSubtrFactor[50]; //[nCorrT1METJet]
   Float_t CorrT1METJet_phi[50];             //[nCorrT1METJet]
   Float_t CorrT1METJet_rawPt[50];           //[nCorrT1METJet]
   Float_t DeepMETResolutionTune_phi;
   Float_t DeepMETResolutionTune_pt;
   Float_t DeepMETResponseTune_phi;
   Float_t DeepMETResponseTune_pt;
   Int_t nElectron;
   Char_t Electron_seediEtaOriX[10];              //[nElectron]
   Bool_t Electron_convVeto[10];                  //[nElectron]
   UChar_t Electron_cutBased[10];                 //[nElectron]
   Bool_t Electron_cutBased_HEEP[10];             //[nElectron]
   Bool_t Electron_isPFcand[10];                  //[nElectron]
   UChar_t Electron_jetNDauCharged[10];           //[nElectron]
   UChar_t Electron_lostHits[10];                 //[nElectron]
   Bool_t Electron_mvaIso_WP80[10];               //[nElectron]
   Bool_t Electron_mvaIso_WP90[10];               //[nElectron]
   Bool_t Electron_mvaNoIso_WP80[10];             //[nElectron]
   Bool_t Electron_mvaNoIso_WP90[10];             //[nElectron]
   UChar_t Electron_seedGain[10];                 //[nElectron]
   UChar_t Electron_tightCharge[10];              //[nElectron]
   Short_t Electron_jetIdx[10];                   //[nElectron]
   Short_t Electron_photonIdx[10];                //[nElectron]
   Short_t Electron_svIdx[10];                    //[nElectron]
   Short_t Electron_fsrPhotonIdx[10];             //[nElectron]
   Int_t Electron_charge[10];                     //[nElectron]
   Int_t Electron_pdgId[10];                      //[nElectron]
   Int_t Electron_seediPhiOriY[10];               //[nElectron]
   Int_t Electron_vidNestedWPBitmap[10];          //[nElectron]
   Int_t Electron_vidNestedWPBitmapHEEP[10];      //[nElectron]
   Float_t Electron_deltaEtaSC[10];               //[nElectron]
   Float_t Electron_dr03EcalRecHitSumEt[10];      //[nElectron]
   Float_t Electron_dr03HcalDepth1TowerSumEt[10]; //[nElectron]
   Float_t Electron_dr03TkSumPt[10];              //[nElectron]
   Float_t Electron_dr03TkSumPtHEEP[10];          //[nElectron]
   Float_t Electron_dxy[10];                      //[nElectron]
   Float_t Electron_dxyErr[10];                   //[nElectron]
   Float_t Electron_dz[10];                       //[nElectron]
   Float_t Electron_dzErr[10];                    //[nElectron]
   Float_t Electron_eInvMinusPInv[10];            //[nElectron]
   Float_t Electron_energyErr[10];                //[nElectron]
   Float_t Electron_eta[10];                      //[nElectron]
   Float_t Electron_hoe[10];                      //[nElectron]
   Float_t Electron_ip3d[10];                     //[nElectron]
   Float_t Electron_jetPtRelv2[10];               //[nElectron]
   Float_t Electron_jetRelIso[10];                //[nElectron]
   Float_t Electron_mass[10];                     //[nElectron]
   Float_t Electron_miniPFRelIso_all[10];         //[nElectron]
   Float_t Electron_miniPFRelIso_chg[10];         //[nElectron]
   Float_t Electron_mvaHZZIso[10];                //[nElectron]
   Float_t Electron_mvaIso[10];                   //[nElectron]
   Float_t Electron_mvaNoIso[10];                 //[nElectron]
   Float_t Electron_pfRelIso03_all[10];           //[nElectron]
   Float_t Electron_pfRelIso03_chg[10];           //[nElectron]
   Float_t Electron_phi[10];                      //[nElectron]
   Float_t Electron_pt[10];                       //[nElectron]
   Float_t Electron_r9[10];                       //[nElectron]
   Float_t Electron_scEtOverPt[10];               //[nElectron]
   Float_t Electron_sieie[10];                    //[nElectron]
   Float_t Electron_sip3d[10];                    //[nElectron]
   Float_t Electron_mvaTTH[10];                   //[nElectron]
   Int_t nFatJet;
   UChar_t FatJet_jetId[15];                        //[nFatJet]
   UChar_t FatJet_nConstituents[15];                //[nFatJet]
   Short_t FatJet_subJetIdx1[15];                   //[nFatJet]
   Short_t FatJet_subJetIdx2[15];                   //[nFatJet]
   Short_t FatJet_electronIdx3SJ[15];               //[nFatJet]
   Short_t FatJet_muonIdx3SJ[15];                   //[nFatJet]
   Float_t FatJet_area[15];                         //[nFatJet]
   Float_t FatJet_btagDDBvLV2[15];                  //[nFatJet]
   Float_t FatJet_btagDDCvBV2[15];                  //[nFatJet]
   Float_t FatJet_btagDDCvLV2[15];                  //[nFatJet]
   Float_t FatJet_btagDeepB[15];                    //[nFatJet]
   Float_t FatJet_btagHbb[15];                      //[nFatJet]
   Float_t FatJet_eta[15];                          //[nFatJet]
   Float_t FatJet_mass[15];                         //[nFatJet]
   Float_t FatJet_msoftdrop[15];                    //[nFatJet]
   Float_t FatJet_n2b1[15];                         //[nFatJet]
   Float_t FatJet_n3b1[15];                         //[nFatJet]
   Float_t FatJet_particleNetWithMass_H4qvsQCD[15]; //[nFatJet]
   Float_t FatJet_particleNetWithMass_HbbvsQCD[15]; //[nFatJet]
   Float_t FatJet_particleNetWithMass_HccvsQCD[15]; //[nFatJet]
   Float_t FatJet_particleNetWithMass_QCD[15];      //[nFatJet]
   Float_t FatJet_particleNetWithMass_TvsQCD[15];   //[nFatJet]
   Float_t FatJet_particleNetWithMass_WvsQCD[15];   //[nFatJet]
   Float_t FatJet_particleNetWithMass_ZvsQCD[15];   //[nFatJet]
   Float_t FatJet_particleNet_QCD[15];              //[nFatJet]
   Float_t FatJet_particleNet_QCD0HF[15];           //[nFatJet]
   Float_t FatJet_particleNet_QCD1HF[15];           //[nFatJet]
   Float_t FatJet_particleNet_QCD2HF[15];           //[nFatJet]
   Float_t FatJet_particleNet_XbbVsQCD[15];         //[nFatJet]
   Float_t FatJet_particleNet_XccVsQCD[15];         //[nFatJet]
   Float_t FatJet_particleNet_XggVsQCD[15];         //[nFatJet]
   Float_t FatJet_particleNet_XqqVsQCD[15];         //[nFatJet]
   Float_t FatJet_particleNet_XteVsQCD[15];         //[nFatJet]
   Float_t FatJet_particleNet_XtmVsQCD[15];         //[nFatJet]
   Float_t FatJet_particleNet_XttVsQCD[15];         //[nFatJet]
   Float_t FatJet_particleNet_massCorr[15];         //[nFatJet]
   Float_t FatJet_phi[15];                          //[nFatJet]
   Float_t FatJet_pt[15];                           //[nFatJet]
   Float_t FatJet_rawFactor[15];                    //[nFatJet]
   Float_t FatJet_tau1[15];                         //[nFatJet]
   Float_t FatJet_tau2[15];                         //[nFatJet]
   Float_t FatJet_tau3[15];                         //[nFatJet]
   Float_t FatJet_tau4[15];                         //[nFatJet]
   Float_t FatJet_lsf3[15];                         //[nFatJet]
   Int_t nFsrPhoton;
   Short_t FsrPhoton_electronIdx[4]; //[nFsrPhoton]
   Short_t FsrPhoton_muonIdx[4];     //[nFsrPhoton]
   Float_t FsrPhoton_dROverEt2[4];   //[nFsrPhoton]
   Float_t FsrPhoton_eta[4];         //[nFsrPhoton]
   Float_t FsrPhoton_phi[4];         //[nFsrPhoton]
   Float_t FsrPhoton_pt[4];          //[nFsrPhoton]
   Float_t FsrPhoton_relIso03[4];    //[nFsrPhoton]
   Int_t nIsoTrack;
   Bool_t IsoTrack_isHighPurityTrack[20]; //[nIsoTrack]
   Bool_t IsoTrack_isPFcand[20];          //[nIsoTrack]
   Bool_t IsoTrack_isFromLostTrack[20];   //[nIsoTrack]
   Short_t IsoTrack_charge[20];           //[nIsoTrack]
   Short_t IsoTrack_fromPV[20];           //[nIsoTrack]
   Int_t IsoTrack_pdgId[20];              //[nIsoTrack]
   Float_t IsoTrack_dxy[20];              //[nIsoTrack]
   Float_t IsoTrack_dz[20];               //[nIsoTrack]
   Float_t IsoTrack_eta[20];              //[nIsoTrack]
   Float_t IsoTrack_pfRelIso03_all[20];   //[nIsoTrack]
   Float_t IsoTrack_pfRelIso03_chg[20];   //[nIsoTrack]
   Float_t IsoTrack_phi[20];              //[nIsoTrack]
   Float_t IsoTrack_pt[20];               //[nIsoTrack]
   Float_t IsoTrack_miniPFRelIso_all[20]; //[nIsoTrack]
   Float_t IsoTrack_miniPFRelIso_chg[20]; //[nIsoTrack]
   Int_t nJet;
   UChar_t Jet_jetId[150];                    //[nJet]
   UChar_t Jet_nConstituents[150];            //[nJet]
   UChar_t Jet_nElectrons[150];               //[nJet]
   UChar_t Jet_nMuons[150];                   //[nJet]
   UChar_t Jet_nSVs[150];                     //[nJet]
   Short_t Jet_electronIdx1[150];             //[nJet]
   Short_t Jet_electronIdx2[150];             //[nJet]
   Short_t Jet_muonIdx1[150];                 //[nJet]
   Short_t Jet_muonIdx2[150];                 //[nJet]
   Short_t Jet_svIdx1[150];                   //[nJet]
   Short_t Jet_svIdx2[150];                   //[nJet]
   Int_t Jet_hfadjacentEtaStripsSize[150];    //[nJet]
   Int_t Jet_hfcentralEtaStripSize[150];      //[nJet]
   Float_t Jet_PNetRegPtRawCorr[150];         //[nJet]
   Float_t Jet_PNetRegPtRawCorrNeutrino[150]; //[nJet]
   Float_t Jet_PNetRegPtRawRes[150];          //[nJet]
   Float_t Jet_area[150];                     //[nJet]
   Float_t Jet_btagDeepFlavB[150];            //[nJet]
   Float_t Jet_btagDeepFlavCvB[150];          //[nJet]
   Float_t Jet_btagDeepFlavCvL[150];          //[nJet]
   Float_t Jet_btagDeepFlavQG[150];           //[nJet]
   Float_t Jet_btagPNetB[150];                //[nJet]
   Float_t Jet_btagPNetCvB[150];              //[nJet]
   Float_t Jet_btagPNetCvL[150];              //[nJet]
   Float_t Jet_btagPNetQvG[150];              //[nJet]
   Float_t Jet_btagPNetTauVJet[150];          //[nJet]
   Float_t Jet_btagRobustParTAK4B[150];       //[nJet]
   Float_t Jet_btagRobustParTAK4CvB[150];     //[nJet]
   Float_t Jet_btagRobustParTAK4CvL[150];     //[nJet]
   Float_t Jet_btagRobustParTAK4QG[150];      //[nJet]
   Float_t Jet_chEmEF[150];                   //[nJet]
   Float_t Jet_chHEF[150];                    //[nJet]
   Float_t Jet_eta[150];                      //[nJet]
   Float_t Jet_hfsigmaEtaEta[150];            //[nJet]
   Float_t Jet_hfsigmaPhiPhi[150];            //[nJet]
   Float_t Jet_mass[150];                     //[nJet]
   Float_t Jet_muEF[150];                     //[nJet]
   Float_t Jet_muonSubtrFactor[150];          //[nJet]
   Float_t Jet_neEmEF[150];                   //[nJet]
   Float_t Jet_neHEF[150];                    //[nJet]
   Float_t Jet_phi[150];                      //[nJet]
   Float_t Jet_pt[150];                       //[nJet]
   Float_t Jet_rawFactor[150];                //[nJet]
   Int_t nLowPtElectron;
   Bool_t LowPtElectron_convVeto[11];          //[nLowPtElectron]
   UChar_t LowPtElectron_convWP[11];           //[nLowPtElectron]
   UChar_t LowPtElectron_lostHits[11];         //[nLowPtElectron]
   Short_t LowPtElectron_electronIdx[11];      //[nLowPtElectron]
   Short_t LowPtElectron_photonIdx[11];        //[nLowPtElectron]
   Int_t LowPtElectron_charge[11];             //[nLowPtElectron]
   Int_t LowPtElectron_pdgId[11];              //[nLowPtElectron]
   Float_t LowPtElectron_ID[11];               //[nLowPtElectron]
   Float_t LowPtElectron_convVtxRadius[11];    //[nLowPtElectron]
   Float_t LowPtElectron_deltaEtaSC[11];       //[nLowPtElectron]
   Float_t LowPtElectron_dxy[11];              //[nLowPtElectron]
   Float_t LowPtElectron_dxyErr[11];           //[nLowPtElectron]
   Float_t LowPtElectron_dz[11];               //[nLowPtElectron]
   Float_t LowPtElectron_dzErr[11];            //[nLowPtElectron]
   Float_t LowPtElectron_eInvMinusPInv[11];    //[nLowPtElectron]
   Float_t LowPtElectron_energyErr[11];        //[nLowPtElectron]
   Float_t LowPtElectron_eta[11];              //[nLowPtElectron]
   Float_t LowPtElectron_hoe[11];              //[nLowPtElectron]
   Float_t LowPtElectron_mass[11];             //[nLowPtElectron]
   Float_t LowPtElectron_miniPFRelIso_all[11]; //[nLowPtElectron]
   Float_t LowPtElectron_miniPFRelIso_chg[11]; //[nLowPtElectron]
   Float_t LowPtElectron_phi[11];              //[nLowPtElectron]
   Float_t LowPtElectron_pt[11];               //[nLowPtElectron]
   Float_t LowPtElectron_ptbiased[11];         //[nLowPtElectron]
   Float_t LowPtElectron_r9[11];               //[nLowPtElectron]
   Float_t LowPtElectron_scEtOverPt[11];       //[nLowPtElectron]
   Float_t LowPtElectron_sieie[11];            //[nLowPtElectron]
   Float_t LowPtElectron_unbiased[11];         //[nLowPtElectron]
   Float_t MET_MetUnclustEnUpDeltaX;
   Float_t MET_MetUnclustEnUpDeltaY;
   Float_t MET_covXX;
   Float_t MET_covXY;
   Float_t MET_covYY;
   Float_t MET_phi;
   Float_t MET_pt;
   Float_t MET_significance;
   Float_t MET_sumEt;
   Float_t MET_sumPtUnclustered;
   Int_t nProton_multiRP;
   UChar_t Proton_multiRP_arm[5];     //[nProton_multiRP]
   Float_t Proton_multiRP_t[5];       //[nProton_multiRP]
   Float_t Proton_multiRP_thetaX[5];  //[nProton_multiRP]
   Float_t Proton_multiRP_thetaY[5];  //[nProton_multiRP]
   Float_t Proton_multiRP_time[5];    //[nProton_multiRP]
   Float_t Proton_multiRP_timeUnc[5]; //[nProton_multiRP]
   Float_t Proton_multiRP_xi[5];      //[nProton_multiRP]
   Int_t nMuon;
   UChar_t Muon_highPtId[50];           //[nMuon]
   Bool_t Muon_highPurity[50];          //[nMuon]
   Bool_t Muon_inTimeMuon[50];          //[nMuon]
   Bool_t Muon_isGlobal[50];            //[nMuon]
   Bool_t Muon_isPFcand[50];            //[nMuon]
   Bool_t Muon_isStandalone[50];        //[nMuon]
   Bool_t Muon_isTracker[50];           //[nMuon]
   UChar_t Muon_jetNDauCharged[50];     //[nMuon]
   Bool_t Muon_looseId[50];             //[nMuon]
   Bool_t Muon_mediumId[50];            //[nMuon]
   Bool_t Muon_mediumPromptId[50];      //[nMuon]
   UChar_t Muon_miniIsoId[50];          //[nMuon]
   UChar_t Muon_multiIsoId[50];         //[nMuon]
   UChar_t Muon_mvaMuID_WP[50];         //[nMuon]
   UChar_t Muon_nStations[50];          //[nMuon]
   UChar_t Muon_nTrackerLayers[50];     //[nMuon]
   UChar_t Muon_pfIsoId[50];            //[nMuon]
   UChar_t Muon_puppiIsoId[50];         //[nMuon]
   Bool_t Muon_softId[50];              //[nMuon]
   Bool_t Muon_softMvaId[50];           //[nMuon]
   UChar_t Muon_tightCharge[50];        //[nMuon]
   Bool_t Muon_tightId[50];             //[nMuon]
   UChar_t Muon_tkIsoId[50];            //[nMuon]
   Bool_t Muon_triggerIdLoose[50];      //[nMuon]
   Short_t Muon_jetIdx[50];             //[nMuon]
   Short_t Muon_svIdx[50];              //[nMuon]
   Short_t Muon_fsrPhotonIdx[50];       //[nMuon]
   Int_t Muon_charge[50];               //[nMuon]
   Int_t Muon_pdgId[50];                //[nMuon]
   Float_t Muon_dxy[50];                //[nMuon]
   Float_t Muon_dxyErr[50];             //[nMuon]
   Float_t Muon_dxybs[50];              //[nMuon]
   Float_t Muon_dz[50];                 //[nMuon]
   Float_t Muon_dzErr[50];              //[nMuon]
   Float_t Muon_eta[50];                //[nMuon]
   Float_t Muon_ip3d[50];               //[nMuon]
   Float_t Muon_jetPtRelv2[50];         //[nMuon]
   Float_t Muon_jetRelIso[50];          //[nMuon]
   Float_t Muon_mass[50];               //[nMuon]
   Float_t Muon_miniPFRelIso_all[50];   //[nMuon]
   Float_t Muon_miniPFRelIso_chg[50];   //[nMuon]
   Float_t Muon_mvaMuID[50];            //[nMuon]
   Float_t Muon_pfRelIso03_all[50];     //[nMuon]
   Float_t Muon_pfRelIso03_chg[50];     //[nMuon]
   Float_t Muon_pfRelIso04_all[50];     //[nMuon]
   Float_t Muon_phi[50];                //[nMuon]
   Float_t Muon_pt[50];                 //[nMuon]
   Float_t Muon_ptErr[50];              //[nMuon]
   Float_t Muon_segmentComp[50];        //[nMuon]
   Float_t Muon_sip3d[50];              //[nMuon]
   Float_t Muon_softMva[50];            //[nMuon]
   Float_t Muon_tkRelIso[50];           //[nMuon]
   Float_t Muon_tunepRelPt[50];         //[nMuon]
   Float_t Muon_bsConstrainedChi2[50];  //[nMuon]
   Float_t Muon_bsConstrainedPt[50];    //[nMuon]
   Float_t Muon_bsConstrainedPtErr[50]; //[nMuon]
   Float_t Muon_mvaLowPt[50];           //[nMuon]
   Float_t Muon_mvaTTH[50];             //[nMuon]
   Float_t         PFMET_covXX;
   Float_t         PFMET_covXY;
   Float_t         PFMET_covYY;
   Float_t         PFMET_phi;
   Float_t         PFMET_phiUnclusteredDown;
   Float_t         PFMET_phiUnclusteredUp;
   Float_t         PFMET_pt;
   Float_t         PFMET_ptUnclusteredDown;
   Float_t         PFMET_ptUnclusteredUp;
   Float_t         PFMET_significance;
   Float_t         PFMET_sumEt;
   Float_t         PFMET_sumPtUnclustered;
   Int_t nPhoton;
   Char_t Photon_seediEtaOriX[15];              //[nPhoton]
   UChar_t Photon_cutBased[15];                 //[nPhoton]
   Bool_t Photon_electronVeto[15];              //[nPhoton]
   Bool_t Photon_hasConversionTracks[15];       //[nPhoton]
   Bool_t Photon_isScEtaEB[15];                 //[nPhoton]
   Bool_t Photon_isScEtaEE[15];                 //[nPhoton]
   Bool_t Photon_mvaID_WP80[15];                //[nPhoton]
   Bool_t Photon_mvaID_WP90[15];                //[nPhoton]
   Bool_t Photon_pixelSeed[15];                 //[nPhoton]
   UChar_t Photon_seedGain[15];                 //[nPhoton]
   Short_t Photon_electronIdx[15];              //[nPhoton]
   Short_t Photon_jetIdx[15];                   //[nPhoton]
   Int_t Photon_seediPhiOriY[15];               //[nPhoton]
   Int_t Photon_vidNestedWPBitmap[15];          //[nPhoton]
   Float_t Photon_ecalPFClusterIso[15];         //[nPhoton]
   Float_t Photon_energyErr[15];                //[nPhoton]
   Float_t Photon_energyRaw[15];                //[nPhoton]
   Float_t Photon_esEffSigmaRR[15];             //[nPhoton]
   Float_t Photon_esEnergyOverRawE[15];         //[nPhoton]
   Float_t Photon_eta[15];                      //[nPhoton]
   Float_t Photon_etaWidth[15];                 //[nPhoton]
   Float_t Photon_haloTaggerMVAVal[15];         //[nPhoton]
   Float_t Photon_hcalPFClusterIso[15];         //[nPhoton]
   Float_t Photon_hoe[15];                      //[nPhoton]
   Float_t Photon_hoe_PUcorr[15];               //[nPhoton]
   Float_t Photon_mvaID[15];                    //[nPhoton]
   Float_t Photon_pfChargedIso[15];             //[nPhoton]
   Float_t Photon_pfChargedIsoPFPV[15];         //[nPhoton]
   Float_t Photon_pfChargedIsoWorstVtx[15];     //[nPhoton]
   Float_t Photon_pfPhoIso03[15];               //[nPhoton]
   Float_t Photon_pfRelIso03_all_quadratic[15]; //[nPhoton]
   Float_t Photon_pfRelIso03_chg_quadratic[15]; //[nPhoton]
   Float_t Photon_phi[15];                      //[nPhoton]
   Float_t Photon_phiWidth[15];                 //[nPhoton]
   Float_t Photon_pt[15];                       //[nPhoton]
   Float_t Photon_r9[15];                       //[nPhoton]
   Float_t Photon_s4[15];                       //[nPhoton]
   Float_t Photon_sieie[15];                    //[nPhoton]
   Float_t Photon_sieip[15];                    //[nPhoton]
   Float_t Photon_sipip[15];                    //[nPhoton]
   Float_t Photon_trkSumPtHollowConeDR03[15];   //[nPhoton]
   Float_t Photon_trkSumPtSolidConeDR04[15];    //[nPhoton]
   Float_t Photon_x_calo[15];                   //[nPhoton]
   Float_t Photon_y_calo[15];                   //[nPhoton]
   Float_t Photon_z_calo[15];                   //[nPhoton]
   Int_t nPPSLocalTrack;
   Int_t PPSLocalTrack_multiRPProtonIdx[5];  //[nPPSLocalTrack]
   Int_t PPSLocalTrack_singleRPProtonIdx[5]; //[nPPSLocalTrack]
   Int_t PPSLocalTrack_decRPId[5];           //[nPPSLocalTrack]
   Int_t PPSLocalTrack_rpType[5];            //[nPPSLocalTrack]
   Float_t PPSLocalTrack_x[5];               //[nPPSLocalTrack]
   Float_t PPSLocalTrack_y[5];               //[nPPSLocalTrack]
   Float_t PPSLocalTrack_time[5];            //[nPPSLocalTrack]
   Float_t PPSLocalTrack_timeUnc[5];         //[nPPSLocalTrack]
   Int_t           Pileup_nPU;
   Int_t           Pileup_sumEOOT;
   Int_t           Pileup_sumLOOT;
   Float_t         Pileup_nTrueInt;
   Float_t         Pileup_pudensity;
   Float_t         Pileup_gpudensity;
   Float_t PuppiMET_phi;
   Float_t PuppiMET_phiJERDown;
   Float_t PuppiMET_phiJERUp;
   Float_t PuppiMET_phiJESDown;
   Float_t PuppiMET_phiJESUp;
   Float_t PuppiMET_phiUnclusteredDown;
   Float_t PuppiMET_phiUnclusteredUp;
   Float_t PuppiMET_pt;
   Float_t PuppiMET_ptJERDown;
   Float_t PuppiMET_ptJERUp;
   Float_t PuppiMET_ptJESDown;
   Float_t PuppiMET_ptJESUp;
   Float_t PuppiMET_ptUnclusteredDown;
   Float_t PuppiMET_ptUnclusteredUp;
   Float_t PuppiMET_sumEt;
   Float_t RawMET_phi;
   Float_t RawMET_pt;
   Float_t RawMET_sumEt;
   Float_t RawPuppiMET_phi;
   Float_t RawPuppiMET_pt;
   Float_t RawPuppiMET_sumEt;
   Float_t Rho_fixedGridRhoAll;
   Float_t Rho_fixedGridRhoFastjetAll;
   Float_t Rho_fixedGridRhoFastjetCentral;
   Float_t Rho_fixedGridRhoFastjetCentralCalo;
   Float_t Rho_fixedGridRhoFastjetCentralChargedPileUp;
   Float_t Rho_fixedGridRhoFastjetCentralNeutral;
   Int_t nSoftActivityJet;
   Float_t SoftActivityJet_eta[6]; //[nSoftActivityJet]
   Float_t SoftActivityJet_phi[6]; //[nSoftActivityJet]
   Float_t SoftActivityJet_pt[6];  //[nSoftActivityJet]
   Int_t SoftActivityJetNjets10;
   Int_t SoftActivityJetNjets2;
   Int_t SoftActivityJetNjets5;
   Float_t SoftActivityJetHT;
   Float_t SoftActivityJetHT10;
   Float_t SoftActivityJetHT2;
   Float_t SoftActivityJetHT5;
   Int_t nProton_singleRP;
   Short_t Proton_singleRP_decRPId[30]; //[nProton_singleRP]
   Float_t Proton_singleRP_thetaY[30];  //[nProton_singleRP]
   Float_t Proton_singleRP_xi[30];      //[nProton_singleRP]
   Int_t nSubJet;
   Float_t SubJet_btagDeepB[30]; //[nSubJet]
   Float_t SubJet_eta[30];       //[nSubJet]
   Float_t SubJet_mass[30];      //[nSubJet]
   Float_t SubJet_n2b1[30];      //[nSubJet]
   Float_t SubJet_n3b1[30];      //[nSubJet]
   Float_t SubJet_phi[30];       //[nSubJet]
   Float_t SubJet_pt[30];        //[nSubJet]
   Float_t SubJet_rawFactor[30]; //[nSubJet]
   Float_t SubJet_tau1[30];      //[nSubJet]
   Float_t SubJet_tau2[30];      //[nSubJet]
   Float_t SubJet_tau3[30];      //[nSubJet]
   Float_t SubJet_tau4[30];      //[nSubJet]
   Int_t nTau;
   UChar_t Tau_decayMode[10];                //[nTau]
   Bool_t Tau_idAntiEleDeadECal[10];         //[nTau]
   UChar_t Tau_idAntiMu[10];                 //[nTau]
   Bool_t Tau_idDecayModeNewDMs[10];         //[nTau]
   Bool_t Tau_idDecayModeOldDMs[10];         //[nTau]
   UChar_t Tau_idDeepTau2017v2p1VSe[10];     //[nTau]
   UChar_t Tau_idDeepTau2017v2p1VSjet[10];   //[nTau]
   UChar_t Tau_idDeepTau2017v2p1VSmu[10];    //[nTau]
   UChar_t Tau_idDeepTau2018v2p5VSe[10];     //[nTau]
   UChar_t Tau_idDeepTau2018v2p5VSjet[10];   //[nTau]
   UChar_t Tau_idDeepTau2018v2p5VSmu[10];    //[nTau]
   UChar_t Tau_nSVs[10];                     //[nTau]
   Short_t Tau_charge[10];                   //[nTau]
   Short_t Tau_decayModePNet[10];            //[nTau]
   Short_t Tau_eleIdx[10];                   //[nTau]
   Short_t Tau_jetIdx[10];                   //[nTau]
   Short_t Tau_muIdx[10];                    //[nTau]
   Short_t Tau_svIdx1[10];                   //[nTau]
   Short_t Tau_svIdx2[10];                   //[nTau]
   Float_t Tau_chargedIso[10];               //[nTau]
   Float_t Tau_dxy[10];                      //[nTau]
   Float_t Tau_dz[10];                       //[nTau]
   Float_t Tau_eta[10];                      //[nTau]
   Float_t Tau_leadTkDeltaEta[10];           //[nTau]
   Float_t Tau_leadTkDeltaPhi[10];           //[nTau]
   Float_t Tau_leadTkPtOverTauPt[10];        //[nTau]
   Float_t Tau_mass[10];                     //[nTau]
   Float_t Tau_neutralIso[10];               //[nTau]
   Float_t Tau_phi[10];                      //[nTau]
   Float_t Tau_photonsOutsideSignalCone[10]; //[nTau]
   Float_t Tau_probDM0PNet[10];              //[nTau]
   Float_t Tau_probDM10PNet[10];             //[nTau]
   Float_t Tau_probDM11PNet[10];             //[nTau]
   Float_t Tau_probDM1PNet[10];              //[nTau]
   Float_t Tau_probDM2PNet[10];              //[nTau]
   Float_t Tau_pt[10];                       //[nTau]
   Float_t Tau_ptCorrPNet[10];               //[nTau]
   Float_t Tau_puCorr[10];                   //[nTau]
   Float_t Tau_qConfPNet[10];                //[nTau]
   Float_t Tau_rawDeepTau2017v2p1VSe[10];    //[nTau]
   Float_t Tau_rawDeepTau2017v2p1VSjet[10];  //[nTau]
   Float_t Tau_rawDeepTau2017v2p1VSmu[10];   //[nTau]
   Float_t Tau_rawDeepTau2018v2p5VSe[10];    //[nTau]
   Float_t Tau_rawDeepTau2018v2p5VSjet[10];  //[nTau]
   Float_t Tau_rawDeepTau2018v2p5VSmu[10];   //[nTau]
   Float_t Tau_rawIso[10];                   //[nTau]
   Float_t Tau_rawIsodR03[10];               //[nTau]
   Float_t Tau_rawPNetVSe[10];               //[nTau]
   Float_t Tau_rawPNetVSjet[10];             //[nTau]
   Float_t Tau_rawPNetVSmu[7];              //[nTau]
   Float_t TkMET_phi;
   Float_t TkMET_pt;
   Float_t TkMET_sumEt;
   Int_t nTrigObj;
   Short_t TrigObj_l1charge[120]; //[nTrigObj]
   UShort_t TrigObj_id[120];      //[nTrigObj]
   Int_t TrigObj_l1iso[120];      //[nTrigObj]
   Int_t TrigObj_filterBits[120]; //[nTrigObj]
   Float_t TrigObj_pt[120];       //[nTrigObj]
   Float_t TrigObj_eta[120];      //[nTrigObj]
   Float_t TrigObj_phi[120];      //[nTrigObj]
   Float_t TrigObj_l1pt[120];     //[nTrigObj]
   Float_t TrigObj_l1pt_2[120];   //[nTrigObj]
   Float_t TrigObj_l2pt[120];     //[nTrigObj]
   Int_t nOtherPV;
   Float_t OtherPV_z[3];     //[nOtherPV]
   Float_t OtherPV_score[3]; //[nOtherPV]
   UChar_t PV_npvs;
   UChar_t PV_npvsGood;
   Float_t PV_ndof;
   Float_t PV_x;
   Float_t PV_y;
   Float_t PV_z;
   Float_t PV_chi2;
   Float_t PV_score;
   Int_t nSV;
   Short_t SV_charge[30];  //[nSV]
   Float_t SV_dlen[30];    //[nSV]
   Float_t SV_dlenSig[30]; //[nSV]
   Float_t SV_dxy[30];     //[nSV]
   Float_t SV_dxySig[30];  //[nSV]
   Float_t SV_pAngle[30];  //[nSV]
   UChar_t SV_ntracks[30]; //[nSV]
   Float_t SV_chi2[30];    //[nSV]
   Float_t SV_eta[30];     //[nSV]
   Float_t SV_mass[30];    //[nSV]
   Float_t SV_ndof[30];    //[nSV]
   Float_t SV_phi[30];     //[nSV]
   Float_t SV_pt[30];      //[nSV]
   Float_t SV_x[30];       //[nSV]
   Float_t SV_y[30];       //[nSV]
   Float_t SV_z[30];       //[nSV]
   Bool_t Flag_HBHENoiseFilter;
   Bool_t Flag_HBHENoiseIsoFilter;
   Bool_t Flag_CSCTightHaloFilter;
   Bool_t Flag_CSCTightHaloTrkMuUnvetoFilter;
   Bool_t Flag_CSCTightHalo2015Filter;
   Bool_t Flag_globalTightHalo2016Filter;
   Bool_t Flag_globalSuperTightHalo2016Filter;
   Bool_t Flag_HcalStripHaloFilter;
   Bool_t Flag_hcalLaserEventFilter;
   Bool_t Flag_EcalDeadCellTriggerPrimitiveFilter;
   Bool_t Flag_EcalDeadCellBoundaryEnergyFilter;
   Bool_t Flag_ecalBadCalibFilter;
   Bool_t Flag_goodVertices;
   Bool_t Flag_eeBadScFilter;
   Bool_t Flag_ecalLaserCorrFilter;
   Bool_t Flag_trkPOGFilters;
   Bool_t Flag_chargedHadronTrackResolutionFilter;
   Bool_t Flag_muonBadTrackFilter;
   Bool_t Flag_BadChargedCandidateFilter;
   Bool_t Flag_BadPFMuonFilter;
   Bool_t Flag_BadPFMuonDzFilter;
   Bool_t Flag_hfNoisyHitsFilter;
   Bool_t Flag_BadChargedCandidateSummer16Filter;
   Bool_t Flag_BadPFMuonSummer16Filter;
   Bool_t Flag_trkPOG_manystripclus53X;
   Bool_t Flag_trkPOG_toomanystripclus53X;
   Bool_t Flag_trkPOG_logErrorTooManyClusters;
   Bool_t Flag_METFilters;
   Bool_t L1_AlwaysTrue;
   Bool_t L1_BPTX_AND_Ref1_VME;
   Bool_t L1_BPTX_AND_Ref3_VME;
   Bool_t L1_BPTX_AND_Ref4_VME;
   Bool_t L1_BPTX_BeamGas_B1_VME;
   Bool_t L1_BPTX_BeamGas_B2_VME;
   Bool_t L1_BPTX_BeamGas_Ref1_VME;
   Bool_t L1_BPTX_BeamGas_Ref2_VME;
   Bool_t L1_BPTX_NotOR_VME;
   Bool_t L1_BPTX_OR_Ref3_VME;
   Bool_t L1_BPTX_OR_Ref4_VME;
   Bool_t L1_BPTX_RefAND_VME;
   Bool_t L1_BptxMinus;
   Bool_t L1_BptxOR;
   Bool_t L1_BptxPlus;
   Bool_t L1_BptxXOR;
   Bool_t L1_CDC_SingleMu_3_er1p2_TOP120_DPHI2p618_3p142;
   Bool_t L1_DoubleEG10_er1p2_dR_Max0p6;
   Bool_t L1_DoubleEG10p5_er1p2_dR_Max0p6;
   Bool_t L1_DoubleEG11_er1p2_dR_Max0p6;
   Bool_t L1_DoubleEG4_er1p2_dR_Max0p9;
   Bool_t L1_DoubleEG4p5_er1p2_dR_Max0p9;
   Bool_t L1_DoubleEG5_er1p2_dR_Max0p9;
   Bool_t L1_DoubleEG5p5_er1p2_dR_Max0p8;
   Bool_t L1_DoubleEG6_er1p2_dR_Max0p8;
   Bool_t L1_DoubleEG6p5_er1p2_dR_Max0p8;
   Bool_t L1_DoubleEG7_er1p2_dR_Max0p8;
   Bool_t L1_DoubleEG7p5_er1p2_dR_Max0p7;
   Bool_t L1_DoubleEG8_er1p2_dR_Max0p7;
   Bool_t L1_DoubleEG8er2p5_HTT260er;
   Bool_t L1_DoubleEG8er2p5_HTT280er;
   Bool_t L1_DoubleEG8er2p5_HTT300er;
   Bool_t L1_DoubleEG8er2p5_HTT320er;
   Bool_t L1_DoubleEG8er2p5_HTT340er;
   Bool_t L1_DoubleEG8p5_er1p2_dR_Max0p7;
   Bool_t L1_DoubleEG9_er1p2_dR_Max0p7;
   Bool_t L1_DoubleEG9p5_er1p2_dR_Max0p6;
   Bool_t L1_DoubleEG_15_10_er2p5;
   Bool_t L1_DoubleEG_20_10_er2p5;
   Bool_t L1_DoubleEG_22_10_er2p5;
   Bool_t L1_DoubleEG_25_12_er2p5;
   Bool_t L1_DoubleEG_25_14_er2p5;
   Bool_t L1_DoubleEG_27_14_er2p5;
   Bool_t L1_DoubleEG_LooseIso16_LooseIso12_er1p5;
   Bool_t L1_DoubleEG_LooseIso18_LooseIso12_er1p5;
   Bool_t L1_DoubleEG_LooseIso20_10_er2p5;
   Bool_t L1_DoubleEG_LooseIso20_LooseIso12_er1p5;
   Bool_t L1_DoubleEG_LooseIso22_10_er2p5;
   Bool_t L1_DoubleEG_LooseIso22_12_er2p5;
   Bool_t L1_DoubleEG_LooseIso22_LooseIso12_er1p5;
   Bool_t L1_DoubleEG_LooseIso25_12_er2p5;
   Bool_t L1_DoubleEG_LooseIso25_LooseIso12_er1p5;
   Bool_t L1_DoubleIsoTau26er2p1_Jet55_RmOvlp_dR0p5;
   Bool_t L1_DoubleIsoTau26er2p1_Jet70_RmOvlp_dR0p5;
   Bool_t L1_DoubleIsoTau28er2p1;
   Bool_t L1_DoubleIsoTau28er2p1_Mass_Max80;
   Bool_t L1_DoubleIsoTau28er2p1_Mass_Max90;
   Bool_t L1_DoubleIsoTau30er2p1;
   Bool_t L1_DoubleIsoTau30er2p1_Mass_Max80;
   Bool_t L1_DoubleIsoTau30er2p1_Mass_Max90;
   Bool_t L1_DoubleIsoTau32er2p1;
   Bool_t L1_DoubleIsoTau34er2p1;
   Bool_t L1_DoubleIsoTau35er2p1;
   Bool_t L1_DoubleIsoTau36er2p1;
   Bool_t L1_DoubleJet100er2p3_dEta_Max1p6;
   Bool_t L1_DoubleJet100er2p5;
   Bool_t L1_DoubleJet112er2p3_dEta_Max1p6;
   Bool_t L1_DoubleJet120er2p5;
   Bool_t L1_DoubleJet150er2p5;
   Bool_t L1_DoubleJet30er2p5_Mass_Min150_dEta_Max1p5;
   Bool_t L1_DoubleJet30er2p5_Mass_Min200_dEta_Max1p5;
   Bool_t L1_DoubleJet30er2p5_Mass_Min250_dEta_Max1p5;
   Bool_t L1_DoubleJet30er2p5_Mass_Min300_dEta_Max1p5;
   Bool_t L1_DoubleJet30er2p5_Mass_Min330_dEta_Max1p5;
   Bool_t L1_DoubleJet30er2p5_Mass_Min360_dEta_Max1p5;
   Bool_t L1_DoubleJet35_Mass_Min450_IsoTau45_RmOvlp;
   Bool_t L1_DoubleJet35_Mass_Min450_IsoTau45er2p1_RmOvlp_dR0p5;
   Bool_t L1_DoubleJet40er2p5;
   Bool_t L1_DoubleJet_100_30_DoubleJet30_Mass_Min620;
   Bool_t L1_DoubleJet_110_35_DoubleJet35_Mass_Min620;
   Bool_t L1_DoubleJet_115_40_DoubleJet40_Mass_Min620;
   Bool_t L1_DoubleJet_115_40_DoubleJet40_Mass_Min620_Jet60TT28;
   Bool_t L1_DoubleJet_120_45_DoubleJet45_Mass_Min620;
   Bool_t L1_DoubleJet_120_45_DoubleJet45_Mass_Min620_Jet60TT28;
   Bool_t L1_DoubleJet_80_30_Mass_Min420_DoubleMu0_SQ;
   Bool_t L1_DoubleJet_80_30_Mass_Min420_IsoTau40_RmOvlp;
   Bool_t L1_DoubleJet_80_30_Mass_Min420_Mu8;
   Bool_t L1_DoubleJet_90_30_DoubleJet30_Mass_Min620;
   Bool_t L1_DoubleLLPJet40;
   Bool_t L1_DoubleLooseIsoEG22er2p1;
   Bool_t L1_DoubleLooseIsoEG24er2p1;
   Bool_t L1_DoubleMu0;
   Bool_t L1_DoubleMu0_Mass_Min1;
   Bool_t L1_DoubleMu0_OQ;
   Bool_t L1_DoubleMu0_SQ;
   Bool_t L1_DoubleMu0_SQ_OS;
   Bool_t L1_DoubleMu0_Upt15_Upt7;
   Bool_t L1_DoubleMu0_Upt5_Upt5;
   Bool_t L1_DoubleMu0_Upt6_IP_Min1_Upt4;
   Bool_t L1_DoubleMu0_dR_Max1p6_Jet90er2p5_dR_Max0p8;
   Bool_t L1_DoubleMu0er1p4_OQ_OS_dEta_Max1p6;
   Bool_t L1_DoubleMu0er1p4_SQ_OS_dR_Max1p4;
   Bool_t L1_DoubleMu0er1p5_SQ;
   Bool_t L1_DoubleMu0er1p5_SQ_OS;
   Bool_t L1_DoubleMu0er1p5_SQ_OS_dR_Max1p4;
   Bool_t L1_DoubleMu0er1p5_SQ_dR_Max1p4;
   Bool_t L1_DoubleMu0er2p0_SQ_OS_dEta_Max1p5;
   Bool_t L1_DoubleMu0er2p0_SQ_OS_dEta_Max1p6;
   Bool_t L1_DoubleMu0er2p0_SQ_OS_dR_Max1p4;
   Bool_t L1_DoubleMu0er2p0_SQ_dEta_Max1p5;
   Bool_t L1_DoubleMu0er2p0_SQ_dEta_Max1p6;
   Bool_t L1_DoubleMu0er2p0_SQ_dR_Max1p4;
   Bool_t L1_DoubleMu18er2p1_SQ;
   Bool_t L1_DoubleMu3_OS_er2p3_Mass_Max14_DoubleEG7p5_er2p1_Mass_Max20;
   Bool_t L1_DoubleMu3_SQ_ETMHF50_HTT60er;
   Bool_t L1_DoubleMu3_SQ_ETMHF50_Jet60er2p5;
   Bool_t L1_DoubleMu3_SQ_ETMHF50_Jet60er2p5_OR_DoubleJet40er2p5;
   Bool_t L1_DoubleMu3_SQ_ETMHF60_Jet60er2p5;
   Bool_t L1_DoubleMu3_SQ_HTT220er;
   Bool_t L1_DoubleMu3_SQ_HTT240er;
   Bool_t L1_DoubleMu3_SQ_HTT260er;
   Bool_t L1_DoubleMu3_dR_Max1p6_Jet90er2p5_dR_Max0p8;
   Bool_t L1_DoubleMu3er2p0_SQ_OS_dR_Max1p4;
   Bool_t L1_DoubleMu4_SQ_EG9er2p5;
   Bool_t L1_DoubleMu4_SQ_OS;
   Bool_t L1_DoubleMu4_SQ_OS_dR_Max1p2;
   Bool_t L1_DoubleMu4p5_SQ_OS;
   Bool_t L1_DoubleMu4p5_SQ_OS_dR_Max1p2;
   Bool_t L1_DoubleMu4p5er2p0_SQ_OS;
   Bool_t L1_DoubleMu4p5er2p0_SQ_OS_Mass_7to18;
   Bool_t L1_DoubleMu4p5er2p0_SQ_OS_Mass_Min7;
   Bool_t L1_DoubleMu5_OS_er2p3_Mass_8to14_DoubleEG3er2p1_Mass_Max20;
   Bool_t L1_DoubleMu5_SQ_EG9er2p5;
   Bool_t L1_DoubleMu8_SQ;
   Bool_t L1_DoubleMu9_SQ;
   Bool_t L1_DoubleMu_12_5;
   Bool_t L1_DoubleMu_15_5_SQ;
   Bool_t L1_DoubleMu_15_7;
   Bool_t L1_DoubleMu_15_7_Mass_Min1;
   Bool_t L1_DoubleMu_15_7_SQ;
   Bool_t L1_DoubleTau70er2p1;
   Bool_t L1_ETM120;
   Bool_t L1_ETM150;
   Bool_t L1_ETMHF100;
   Bool_t L1_ETMHF100_HTT60er;
   Bool_t L1_ETMHF110;
   Bool_t L1_ETMHF110_HTT60er;
   Bool_t L1_ETMHF110_HTT60er_NotSecondBunchInTrain;
   Bool_t L1_ETMHF120;
   Bool_t L1_ETMHF120_HTT60er;
   Bool_t L1_ETMHF120_NotSecondBunchInTrain;
   Bool_t L1_ETMHF130;
   Bool_t L1_ETMHF130_HTT60er;
   Bool_t L1_ETMHF140;
   Bool_t L1_ETMHF150;
   Bool_t L1_ETMHF70;
   Bool_t L1_ETMHF70_HTT60er;
   Bool_t L1_ETMHF80;
   Bool_t L1_ETMHF80_HTT60er;
   Bool_t L1_ETMHF90;
   Bool_t L1_ETMHF90_HTT60er;
   Bool_t L1_ETMHF90_SingleJet60er2p5_dPhi_Min2p1;
   Bool_t L1_ETMHF90_SingleJet60er2p5_dPhi_Min2p6;
   Bool_t L1_ETMHF90_SingleJet80er2p5_dPhi_Min2p1;
   Bool_t L1_ETMHF90_SingleJet80er2p5_dPhi_Min2p6;
   Bool_t L1_ETT1200;
   Bool_t L1_ETT1600;
   Bool_t L1_ETT2000;
   Bool_t L1_FirstBunchAfterTrain;
   Bool_t L1_FirstBunchBeforeTrain;
   Bool_t L1_FirstBunchInTrain;
   Bool_t L1_FirstCollisionInOrbit;
   Bool_t L1_FirstCollisionInTrain;
   Bool_t L1_HCAL_LaserMon_Trig;
   Bool_t L1_HCAL_LaserMon_Veto;
   Bool_t L1_HTT120_SingleLLPJet40;
   Bool_t L1_HTT120er;
   Bool_t L1_HTT160_SingleLLPJet50;
   Bool_t L1_HTT160er;
   Bool_t L1_HTT200_SingleLLPJet60;
   Bool_t L1_HTT200er;
   Bool_t L1_HTT240_SingleLLPJet70;
   Bool_t L1_HTT255er;
   Bool_t L1_HTT280er;
   Bool_t L1_HTT280er_QuadJet_70_55_40_35_er2p5;
   Bool_t L1_HTT320er;
   Bool_t L1_HTT320er_QuadJet_70_55_40_40_er2p5;
   Bool_t L1_HTT320er_QuadJet_80_60_er2p1_45_40_er2p3;
   Bool_t L1_HTT320er_QuadJet_80_60_er2p1_50_45_er2p3;
   Bool_t L1_HTT360er;
   Bool_t L1_HTT400er;
   Bool_t L1_HTT450er;
   Bool_t L1_IsoEG32er2p5_Mt40;
   Bool_t L1_IsoTau52er2p1_QuadJet36er2p5;
   Bool_t L1_IsolatedBunch;
   Bool_t L1_LastBunchInTrain;
   Bool_t L1_LastCollisionInTrain;
   Bool_t L1_LooseIsoEG22er2p1_IsoTau26er2p1_dR_Min0p3;
   Bool_t L1_LooseIsoEG22er2p1_Tau70er2p1_dR_Min0p3;
   Bool_t L1_LooseIsoEG24er2p1_HTT100er;
   Bool_t L1_LooseIsoEG24er2p1_IsoTau27er2p1_dR_Min0p3;
   Bool_t L1_LooseIsoEG26er2p1_HTT100er;
   Bool_t L1_LooseIsoEG26er2p1_Jet34er2p5_dR_Min0p3;
   Bool_t L1_LooseIsoEG28er2p1_HTT100er;
   Bool_t L1_LooseIsoEG28er2p1_Jet34er2p5_dR_Min0p3;
   Bool_t L1_LooseIsoEG30er2p1_HTT100er;
   Bool_t L1_LooseIsoEG30er2p1_Jet34er2p5_dR_Min0p3;
   Bool_t L1_MinimumBiasHF0;
   Bool_t L1_MinimumBiasHF0_AND_BptxAND;
   Bool_t L1_Mu10er2p3_Jet32er2p3_dR_Max0p4_DoubleJet32er2p3_dEta_Max1p6;
   Bool_t L1_Mu12er2p3_Jet40er2p1_dR_Max0p4_DoubleJet40er2p1_dEta_Max1p6;
   Bool_t L1_Mu12er2p3_Jet40er2p3_dR_Max0p4_DoubleJet40er2p3_dEta_Max1p6;
   Bool_t L1_Mu18er2p1_Tau24er2p1;
   Bool_t L1_Mu18er2p1_Tau26er2p1;
   Bool_t L1_Mu18er2p1_Tau26er2p1_Jet55;
   Bool_t L1_Mu18er2p1_Tau26er2p1_Jet70;
   Bool_t L1_Mu20_EG10er2p5;
   Bool_t L1_Mu22er2p1_IsoTau28er2p1;
   Bool_t L1_Mu22er2p1_IsoTau30er2p1;
   Bool_t L1_Mu22er2p1_IsoTau32er2p1;
   Bool_t L1_Mu22er2p1_IsoTau34er2p1;
   Bool_t L1_Mu22er2p1_IsoTau36er2p1;
   Bool_t L1_Mu22er2p1_IsoTau40er2p1;
   Bool_t L1_Mu22er2p1_Tau70er2p1;
   Bool_t L1_Mu3_Jet120er2p5_dR_Max0p4;
   Bool_t L1_Mu3_Jet120er2p5_dR_Max0p8;
   Bool_t L1_Mu3_Jet16er2p5_dR_Max0p4;
   Bool_t L1_Mu3_Jet30er2p5;
   Bool_t L1_Mu3_Jet35er2p5_dR_Max0p4;
   Bool_t L1_Mu3_Jet60er2p5_dR_Max0p4;
   Bool_t L1_Mu3_Jet80er2p5_dR_Max0p4;
   Bool_t L1_Mu3er1p5_Jet100er2p5_ETMHF40;
   Bool_t L1_Mu3er1p5_Jet100er2p5_ETMHF50;
   Bool_t L1_Mu5_EG23er2p5;
   Bool_t L1_Mu5_LooseIsoEG20er2p5;
   Bool_t L1_Mu6_DoubleEG10er2p5;
   Bool_t L1_Mu6_DoubleEG12er2p5;
   Bool_t L1_Mu6_DoubleEG15er2p5;
   Bool_t L1_Mu6_DoubleEG17er2p5;
   Bool_t L1_Mu6_HTT240er;
   Bool_t L1_Mu6_HTT250er;
   Bool_t L1_Mu7_EG20er2p5;
   Bool_t L1_Mu7_EG23er2p5;
   Bool_t L1_Mu7_LooseIsoEG20er2p5;
   Bool_t L1_Mu7_LooseIsoEG23er2p5;
   Bool_t L1_NotBptxOR;
   Bool_t L1_QuadJet60er2p5;
   Bool_t L1_QuadJet_95_75_65_20_DoubleJet_75_65_er2p5_Jet20_FWD3p0;
   Bool_t L1_QuadMu0;
   Bool_t L1_QuadMu0_OQ;
   Bool_t L1_QuadMu0_SQ;
   Bool_t L1_SecondBunchInTrain;
   Bool_t L1_SecondLastBunchInTrain;
   Bool_t L1_SingleEG10er2p5;
   Bool_t L1_SingleEG15er2p5;
   Bool_t L1_SingleEG26er2p5;
   Bool_t L1_SingleEG28_FWD2p5;
   Bool_t L1_SingleEG28er1p5;
   Bool_t L1_SingleEG28er2p1;
   Bool_t L1_SingleEG28er2p5;
   Bool_t L1_SingleEG34er2p5;
   Bool_t L1_SingleEG36er2p5;
   Bool_t L1_SingleEG38er2p5;
   Bool_t L1_SingleEG40er2p5;
   Bool_t L1_SingleEG42er2p5;
   Bool_t L1_SingleEG45er2p5;
   Bool_t L1_SingleEG50;
   Bool_t L1_SingleEG60;
   Bool_t L1_SingleEG8er2p5;
   Bool_t L1_SingleIsoEG24er1p5;
   Bool_t L1_SingleIsoEG24er2p1;
   Bool_t L1_SingleIsoEG26er1p5;
   Bool_t L1_SingleIsoEG26er2p1;
   Bool_t L1_SingleIsoEG26er2p5;
   Bool_t L1_SingleIsoEG28_FWD2p5;
   Bool_t L1_SingleIsoEG28er1p5;
   Bool_t L1_SingleIsoEG28er2p1;
   Bool_t L1_SingleIsoEG28er2p5;
   Bool_t L1_SingleIsoEG30er2p1;
   Bool_t L1_SingleIsoEG30er2p5;
   Bool_t L1_SingleIsoEG32er2p1;
   Bool_t L1_SingleIsoEG32er2p5;
   Bool_t L1_SingleIsoEG34er2p5;
   Bool_t L1_SingleIsoTau32er2p1;
   Bool_t L1_SingleJet10erHE;
   Bool_t L1_SingleJet120;
   Bool_t L1_SingleJet120_FWD3p0;
   Bool_t L1_SingleJet120er2p5;
   Bool_t L1_SingleJet12erHE;
   Bool_t L1_SingleJet140er2p5;
   Bool_t L1_SingleJet140er2p5_ETMHF70;
   Bool_t L1_SingleJet140er2p5_ETMHF80;
   Bool_t L1_SingleJet140er2p5_ETMHF90;
   Bool_t L1_SingleJet160er2p5;
   Bool_t L1_SingleJet180;
   Bool_t L1_SingleJet180er2p5;
   Bool_t L1_SingleJet200;
   Bool_t L1_SingleJet20er2p5_NotBptxOR;
   Bool_t L1_SingleJet20er2p5_NotBptxOR_3BX;
   Bool_t L1_SingleJet35;
   Bool_t L1_SingleJet35_FWD3p0;
   Bool_t L1_SingleJet35er2p5;
   Bool_t L1_SingleJet43er2p5_NotBptxOR_3BX;
   Bool_t L1_SingleJet46er2p5_NotBptxOR_3BX;
   Bool_t L1_SingleJet60;
   Bool_t L1_SingleJet60_FWD3p0;
   Bool_t L1_SingleJet60er2p5;
   Bool_t L1_SingleJet8erHE;
   Bool_t L1_SingleJet90;
   Bool_t L1_SingleJet90_FWD3p0;
   Bool_t L1_SingleJet90er2p5;
   Bool_t L1_SingleLooseIsoEG26er1p5;
   Bool_t L1_SingleLooseIsoEG26er2p5;
   Bool_t L1_SingleLooseIsoEG28_FWD2p5;
   Bool_t L1_SingleLooseIsoEG28er1p5;
   Bool_t L1_SingleLooseIsoEG28er2p1;
   Bool_t L1_SingleLooseIsoEG28er2p5;
   Bool_t L1_SingleLooseIsoEG30er1p5;
   Bool_t L1_SingleLooseIsoEG30er2p5;
   Bool_t L1_SingleMu0_BMTF;
   Bool_t L1_SingleMu0_DQ;
   Bool_t L1_SingleMu0_EMTF;
   Bool_t L1_SingleMu0_OMTF;
   Bool_t L1_SingleMu10er1p5;
   Bool_t L1_SingleMu12_DQ_BMTF;
   Bool_t L1_SingleMu12_DQ_EMTF;
   Bool_t L1_SingleMu12_DQ_OMTF;
   Bool_t L1_SingleMu12er1p5;
   Bool_t L1_SingleMu14er1p5;
   Bool_t L1_SingleMu15_DQ;
   Bool_t L1_SingleMu16er1p5;
   Bool_t L1_SingleMu18;
   Bool_t L1_SingleMu18er1p5;
   Bool_t L1_SingleMu20;
   Bool_t L1_SingleMu22;
   Bool_t L1_SingleMu22_BMTF;
   Bool_t L1_SingleMu22_DQ;
   Bool_t L1_SingleMu22_EMTF;
   Bool_t L1_SingleMu22_OMTF;
   Bool_t L1_SingleMu22_OQ;
   Bool_t L1_SingleMu25;
   Bool_t L1_SingleMu3;
   Bool_t L1_SingleMu5;
   Bool_t L1_SingleMu6er1p5;
   Bool_t L1_SingleMu7;
   Bool_t L1_SingleMu7_DQ;
   Bool_t L1_SingleMu7er1p5;
   Bool_t L1_SingleMu8er1p5;
   Bool_t L1_SingleMu9er1p5;
   Bool_t L1_SingleMuCosmics;
   Bool_t L1_SingleMuCosmics_BMTF;
   Bool_t L1_SingleMuCosmics_EMTF;
   Bool_t L1_SingleMuCosmics_OMTF;
   Bool_t L1_SingleMuOpen;
   Bool_t L1_SingleMuOpen_NotBptxOR;
   Bool_t L1_SingleMuOpen_er1p1_NotBptxOR_3BX;
   Bool_t L1_SingleMuOpen_er1p4_NotBptxOR_3BX;
   Bool_t L1_SingleMuShower_Nominal;
   Bool_t L1_SingleMuShower_Tight;
   Bool_t L1_SingleTau120er2p1;
   Bool_t L1_SingleTau130er2p1;
   Bool_t L1_SingleTau70er2p1;
   Bool_t L1_TOTEM_1;
   Bool_t L1_TOTEM_2;
   Bool_t L1_TOTEM_3;
   Bool_t L1_TOTEM_4;
   Bool_t L1_TripleEG16er2p5;
   Bool_t L1_TripleEG_16_12_8_er2p5;
   Bool_t L1_TripleEG_16_15_8_er2p5;
   Bool_t L1_TripleEG_18_17_8_er2p5;
   Bool_t L1_TripleEG_18_18_12_er2p5;
   Bool_t L1_TripleJet_100_80_70_DoubleJet_80_70_er2p5;
   Bool_t L1_TripleJet_105_85_75_DoubleJet_85_75_er2p5;
   Bool_t L1_TripleJet_95_75_65_DoubleJet_75_65_er2p5;
   Bool_t L1_TripleMu0;
   Bool_t L1_TripleMu0_OQ;
   Bool_t L1_TripleMu0_SQ;
   Bool_t L1_TripleMu3;
   Bool_t L1_TripleMu3_SQ;
   Bool_t L1_TripleMu_2SQ_1p5SQ_0OQ;
   Bool_t L1_TripleMu_2SQ_1p5SQ_0OQ_Mass_Max12;
   Bool_t L1_TripleMu_3SQ_2p5SQ_0OQ_Mass_Max12;
   Bool_t L1_TripleMu_5SQ_3SQ_0OQ;
   Bool_t L1_TripleMu_5SQ_3SQ_0OQ_DoubleMu_5_3_SQ_OS_Mass_Max9;
   Bool_t L1_TripleMu_5SQ_3SQ_0_DoubleMu_5_3_SQ_OS_Mass_Max9;
   Bool_t L1_TripleMu_5_3_3;
   Bool_t L1_TripleMu_5_3_3_SQ;
   Bool_t L1_TripleMu_5_3p5_2p5;
   Bool_t L1_TripleMu_5_3p5_2p5_DoubleMu_5_2p5_OS_Mass_5to17;
   Bool_t L1_TripleMu_5_3p5_2p5_OQ_DoubleMu_5_2p5_OQ_OS_Mass_5to17;
   Bool_t L1_TripleMu_5_4_2p5_DoubleMu_5_2p5_OS_Mass_5to17;
   Bool_t L1_TripleMu_5_5_3;
   Bool_t L1_UnpairedBunchBptxMinus;
   Bool_t L1_UnpairedBunchBptxPlus;
   Bool_t L1_ZeroBias;
   Bool_t L1_ZeroBias_copy;
   Bool_t L1_UnprefireableEvent;
   Bool_t L1Reco_step;
   Bool_t Flag_HBHENoiseFilter_pRECO;
   Bool_t Flag_HBHENoiseIsoFilter_pRECO;
   Bool_t Flag_CSCTightHaloFilter_pRECO;
   Bool_t Flag_CSCTightHaloTrkMuUnvetoFilter_pRECO;
   Bool_t Flag_CSCTightHalo2015Filter_pRECO;
   Bool_t Flag_globalTightHalo2016Filter_pRECO;
   Bool_t Flag_globalSuperTightHalo2016Filter_pRECO;
   Bool_t Flag_HcalStripHaloFilter_pRECO;
   Bool_t Flag_hcalLaserEventFilter_pRECO;
   Bool_t Flag_EcalDeadCellTriggerPrimitiveFilter_pRECO;
   Bool_t Flag_EcalDeadCellBoundaryEnergyFilter_pRECO;
   Bool_t Flag_ecalBadCalibFilter_pRECO;
   Bool_t Flag_goodVertices_pRECO;
   Bool_t Flag_eeBadScFilter_pRECO;
   Bool_t Flag_ecalLaserCorrFilter_pRECO;
   Bool_t Flag_trkPOGFilters_pRECO;
   Bool_t Flag_chargedHadronTrackResolutionFilter_pRECO;
   Bool_t Flag_muonBadTrackFilter_pRECO;
   Bool_t Flag_BadChargedCandidateFilter_pRECO;
   Bool_t Flag_BadPFMuonFilter_pRECO;
   Bool_t Flag_BadPFMuonDzFilter_pRECO;
   Bool_t Flag_hfNoisyHitsFilter_pRECO;
   Bool_t Flag_BadChargedCandidateSummer16Filter_pRECO;
   Bool_t Flag_BadPFMuonSummer16Filter_pRECO;
   Bool_t Flag_trkPOG_manystripclus53X_pRECO;
   Bool_t Flag_trkPOG_toomanystripclus53X_pRECO;
   Bool_t Flag_trkPOG_logErrorTooManyClusters_pRECO;
   Bool_t Flag_METFilters_pRECO;
   Bool_t HLTriggerFirstPath;
   Bool_t HLT_AK8PFJet360_TrimMass30;
   Bool_t HLT_AK8PFJet380_TrimMass30;
   Bool_t HLT_AK8PFJet400_TrimMass30;
   Bool_t HLT_AK8PFJet420_TrimMass30;
   Bool_t HLT_AK8PFJet400_MassSD30;
   Bool_t HLT_AK8PFJet420_MassSD30;
   Bool_t HLT_AK8PFJet450_MassSD30;
   Bool_t HLT_AK8DiPFJet250_250_MassSD30;
   Bool_t HLT_AK8DiPFJet250_250_MassSD50;
   Bool_t HLT_AK8DiPFJet260_260_MassSD30;
   Bool_t HLT_AK8DiPFJet270_270_MassSD30;
   Bool_t HLT_AK8PFHT750_TrimMass50;
   Bool_t HLT_AK8PFHT800_TrimMass50;
   Bool_t HLT_AK8PFHT850_TrimMass50;
   Bool_t HLT_AK8PFHT900_TrimMass50;
   Bool_t HLT_CaloJet500_NoJetID;
   Bool_t HLT_CaloJet550_NoJetID;
   Bool_t HLT_DoubleMu5_Upsilon_DoubleEle3_CaloIdL_TrackIdL;
   Bool_t HLT_DoubleMu3_DoubleEle7p5_CaloIdL_TrackIdL_Upsilon;
   Bool_t HLT_Trimuon5_3p5_2_Upsilon_Muon;
   Bool_t HLT_TrimuonOpen_5_3p5_2_Upsilon_Muon;
   Bool_t HLT_DoubleEle25_CaloIdL_MW;
   Bool_t HLT_DoubleEle27_CaloIdL_MW;
   Bool_t HLT_DoubleEle33_CaloIdL_MW;
   Bool_t HLT_DoubleEle24_eta2p1_WPTight_Gsf;
   Bool_t HLT_DoubleEle8_CaloIdM_TrackIdM_Mass8_DZ_PFHT350;
   Bool_t HLT_DoubleEle8_CaloIdM_TrackIdM_Mass8_PFHT350;
   Bool_t HLT_Mu27_Ele37_CaloIdL_MW;
   Bool_t HLT_Mu37_Ele27_CaloIdL_MW;
   Bool_t HLT_Mu37_TkMu27;
   Bool_t HLT_DoubleMu4_3_Bs;
   Bool_t HLT_DoubleMu4_3_Jpsi;
   Bool_t HLT_DoubleMu4_3_LowMass;
   Bool_t HLT_DoubleMu4_LowMass_Displaced;
   Bool_t HLT_Mu0_L1DoubleMu;
   Bool_t HLT_Mu4_L1DoubleMu;
   Bool_t HLT_DoubleMu4_3_Photon4_BsToMMG;
   Bool_t HLT_DoubleMu4_3_Displaced_Photon4_BsToMMG;
   Bool_t HLT_DoubleMu3_Trk_Tau3mu;
   Bool_t HLT_DoubleMu3_TkMu_DsTau3Mu;
   Bool_t HLT_DoubleMu4_Mass3p8_DZ_PFHT350;
   Bool_t HLT_DoubleMu4_MuMuTrk_Displaced;
   Bool_t HLT_Mu3_PFJet40;
   Bool_t HLT_Mu7p5_L2Mu2_Jpsi;
   Bool_t HLT_Mu7p5_L2Mu2_Upsilon;
   Bool_t HLT_Mu3_L1SingleMu5orSingleMu7;
   Bool_t HLT_DoublePhoton33_CaloIdL;
   Bool_t HLT_DoublePhoton70;
   Bool_t HLT_DoublePhoton85;
   Bool_t HLT_Ele15_WPLoose_Gsf;
   Bool_t HLT_Ele20_WPLoose_Gsf;
   Bool_t HLT_DiEle27_WPTightCaloOnly_L1DoubleEG;
   Bool_t HLT_Ele27_WPTight_Gsf;
   Bool_t HLT_Ele28_WPTight_Gsf;
   Bool_t HLT_Ele30_WPTight_Gsf;
   Bool_t HLT_Ele32_WPTight_Gsf;
   Bool_t HLT_Ele35_WPTight_Gsf;
   Bool_t HLT_Ele35_WPTight_Gsf_L1EGMT;
   Bool_t HLT_Ele38_WPTight_Gsf;
   Bool_t HLT_Ele40_WPTight_Gsf;
   Bool_t HLT_Ele32_WPTight_Gsf_L1DoubleEG;
   Bool_t HLT_HT300_Beamspot;
   Bool_t HLT_ZeroBias_Beamspot;
   Bool_t HLT_IsoMu20_eta2p1_TightChargedIsoPFTauHPS27_eta2p1_TightID_CrossL1;
   Bool_t HLT_IsoMu27_MediumDeepTauPFTauHPS20_eta2p1_SingleL1;
   Bool_t HLT_IsoMu20;
   Bool_t HLT_IsoMu24;
   Bool_t HLT_IsoMu24_eta2p1;
   Bool_t HLT_IsoMu27;
   Bool_t HLT_UncorrectedJetE30_NoBPTX;
   Bool_t HLT_UncorrectedJetE30_NoBPTX3BX;
   Bool_t HLT_UncorrectedJetE60_NoBPTX3BX;
   Bool_t HLT_UncorrectedJetE70_NoBPTX3BX;
   Bool_t HLT_L1SingleMu18;
   Bool_t HLT_L1SingleMu25;
   Bool_t HLT_L1SingleMuCosmics;
   Bool_t HLT_L2Mu10_NoVertex_NoBPTX3BX;
   Bool_t HLT_L2Mu10_NoVertex_NoBPTX;
   Bool_t HLT_L2Mu45_NoVertex_3Sta_NoBPTX3BX;
   Bool_t HLT_L2Mu40_NoVertex_3Sta_NoBPTX3BX;
   Bool_t HLT_L2Mu23NoVtx_2Cha;
   Bool_t HLT_L2Mu23NoVtx_2Cha_CosmicSeed;
   Bool_t HLT_DoubleL2Mu30NoVtx_2Cha_CosmicSeed_Eta2p4;
   Bool_t HLT_DoubleL2Mu30NoVtx_2Cha_Eta2p4;
   Bool_t HLT_DoubleL2Mu50;
   Bool_t HLT_DoubleL2Mu23NoVtx_2Cha_CosmicSeed;
   Bool_t HLT_DoubleL2Mu25NoVtx_2Cha_CosmicSeed;
   Bool_t HLT_DoubleL2Mu25NoVtx_2Cha_CosmicSeed_Eta2p4;
   Bool_t HLT_DoubleL2Mu23NoVtx_2Cha;
   Bool_t HLT_DoubleL2Mu25NoVtx_2Cha;
   Bool_t HLT_DoubleL2Mu25NoVtx_2Cha_Eta2p4;
   Bool_t HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL;
   Bool_t HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL;
   Bool_t HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ;
   Bool_t HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ;
   Bool_t HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8;
   Bool_t HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_Mass8;
   Bool_t HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8;
   Bool_t HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_Mass3p8;
   Bool_t HLT_Mu25_TkMu0_Onia;
   Bool_t HLT_Mu30_TkMu0_Psi;
   Bool_t HLT_Mu30_TkMu0_Upsilon;
   Bool_t HLT_Mu20_TkMu0_Phi;
   Bool_t HLT_Mu25_TkMu0_Phi;
   Bool_t HLT_Mu15;
   Bool_t HLT_Mu20;
   Bool_t HLT_Mu27;
   Bool_t HLT_Mu50;
   Bool_t HLT_Mu55;
   Bool_t HLT_CascadeMu100;
   Bool_t HLT_HighPtTkMu100;
   Bool_t HLT_DiPFJetAve40;
   Bool_t HLT_DiPFJetAve60;
   Bool_t HLT_DiPFJetAve80;
   Bool_t HLT_DiPFJetAve140;
   Bool_t HLT_DiPFJetAve200;
   Bool_t HLT_DiPFJetAve260;
   Bool_t HLT_DiPFJetAve320;
   Bool_t HLT_DiPFJetAve400;
   Bool_t HLT_DiPFJetAve500;
   Bool_t HLT_DiPFJetAve60_HFJEC;
   Bool_t HLT_DiPFJetAve80_HFJEC;
   Bool_t HLT_DiPFJetAve100_HFJEC;
   Bool_t HLT_DiPFJetAve160_HFJEC;
   Bool_t HLT_DiPFJetAve220_HFJEC;
   Bool_t HLT_DiPFJetAve300_HFJEC;
   Bool_t HLT_AK8PFJet40;
   Bool_t HLT_AK8PFJet60;
   Bool_t HLT_AK8PFJet80;
   Bool_t HLT_AK8PFJet140;
   Bool_t HLT_AK8PFJet200;
   Bool_t HLT_AK8PFJet260;
   Bool_t HLT_AK8PFJet320;
   Bool_t HLT_AK8PFJet400;
   Bool_t HLT_AK8PFJet450;
   Bool_t HLT_AK8PFJet500;
   Bool_t HLT_AK8PFJet550;
   Bool_t HLT_PFJet40;
   Bool_t HLT_PFJet60;
   Bool_t HLT_PFJet80;
   Bool_t HLT_PFJet110;
   Bool_t HLT_PFJet140;
   Bool_t HLT_PFJet200;
   Bool_t HLT_PFJet260;
   Bool_t HLT_PFJet320;
   Bool_t HLT_PFJet400;
   Bool_t HLT_PFJet450;
   Bool_t HLT_PFJet500;
   Bool_t HLT_PFJet550;
   Bool_t HLT_PFJetFwd15;
   Bool_t HLT_PFJetFwd25;
   Bool_t HLT_PFJetFwd40;
   Bool_t HLT_PFJetFwd60;
   Bool_t HLT_PFJetFwd80;
   Bool_t HLT_PFJetFwd140;
   Bool_t HLT_PFJetFwd200;
   Bool_t HLT_PFJetFwd260;
   Bool_t HLT_PFJetFwd320;
   Bool_t HLT_PFJetFwd400;
   Bool_t HLT_PFJetFwd450;
   Bool_t HLT_PFJetFwd500;
   Bool_t HLT_AK8PFJetFwd15;
   Bool_t HLT_AK8PFJetFwd25;
   Bool_t HLT_AK8PFJetFwd40;
   Bool_t HLT_AK8PFJetFwd60;
   Bool_t HLT_AK8PFJetFwd80;
   Bool_t HLT_AK8PFJetFwd140;
   Bool_t HLT_AK8PFJetFwd200;
   Bool_t HLT_AK8PFJetFwd260;
   Bool_t HLT_AK8PFJetFwd320;
   Bool_t HLT_AK8PFJetFwd400;
   Bool_t HLT_AK8PFJetFwd450;
   Bool_t HLT_AK8PFJetFwd500;
   Bool_t HLT_PFHT180;
   Bool_t HLT_PFHT250;
   Bool_t HLT_PFHT370;
   Bool_t HLT_PFHT430;
   Bool_t HLT_PFHT510;
   Bool_t HLT_PFHT590;
   Bool_t HLT_PFHT680;
   Bool_t HLT_PFHT780;
   Bool_t HLT_PFHT890;
   Bool_t HLT_PFHT1050;
   Bool_t HLT_PFHT500_PFMET100_PFMHT100_IDTight;
   Bool_t HLT_PFHT500_PFMET110_PFMHT110_IDTight;
   Bool_t HLT_PFHT700_PFMET85_PFMHT85_IDTight;
   Bool_t HLT_PFHT800_PFMET75_PFMHT75_IDTight;
   Bool_t HLT_PFMET110_PFMHT110_IDTight;
   Bool_t HLT_PFMET120_PFMHT120_IDTight;
   Bool_t HLT_PFMET130_PFMHT130_IDTight;
   Bool_t HLT_PFMET140_PFMHT140_IDTight;
   Bool_t HLT_PFMET120_PFMHT120_IDTight_PFHT60;
   Bool_t HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60;
   Bool_t HLT_PFMETTypeOne120_PFMHT120_IDTight_PFHT60;
   Bool_t HLT_PFMETTypeOne110_PFMHT110_IDTight;
   Bool_t HLT_PFMETTypeOne120_PFMHT120_IDTight;
   Bool_t HLT_PFMETTypeOne130_PFMHT130_IDTight;
   Bool_t HLT_PFMETTypeOne140_PFMHT140_IDTight;
   Bool_t HLT_PFMETNoMu110_PFMHTNoMu110_IDTight;
   Bool_t HLT_PFMETNoMu120_PFMHTNoMu120_IDTight;
   Bool_t HLT_PFMETNoMu130_PFMHTNoMu130_IDTight;
   Bool_t HLT_PFMETNoMu140_PFMHTNoMu140_IDTight;
   Bool_t HLT_PFMETNoMu110_PFMHTNoMu110_IDTight_FilterHF;
   Bool_t HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_FilterHF;
   Bool_t HLT_PFMETNoMu130_PFMHTNoMu130_IDTight_FilterHF;
   Bool_t HLT_PFMETNoMu140_PFMHTNoMu140_IDTight_FilterHF;
   Bool_t HLT_L1ETMHadSeeds;
   Bool_t HLT_CaloMHT90;
   Bool_t HLT_CaloMET90_NotCleaned;
   Bool_t HLT_CaloMET350_NotCleaned;
   Bool_t HLT_PFMET200_NotCleaned;
   Bool_t HLT_PFMET250_NotCleaned;
   Bool_t HLT_PFMET300_NotCleaned;
   Bool_t HLT_PFMET200_BeamHaloCleaned;
   Bool_t HLT_PFMETTypeOne200_BeamHaloCleaned;
   Bool_t HLT_MET105_IsoTrk50;
   Bool_t HLT_MET120_IsoTrk50;
   Bool_t HLT_SingleJet30_Mu12_SinglePFJet40;
   Bool_t HLT_Mu12eta2p3;
   Bool_t HLT_Mu12eta2p3_PFJet40;
   Bool_t HLT_Mu12_DoublePFJets40_PFBTagDeepCSV_p71;
   Bool_t HLT_Mu12_DoublePFJets100_PFBTagDeepCSV_p71;
   Bool_t HLT_Mu12_DoublePFJets200_PFBTagDeepCSV_p71;
   Bool_t HLT_Mu12_DoublePFJets350_PFBTagDeepCSV_p71;
   Bool_t HLT_Mu12_DoublePFJets40MaxDeta1p6_DoublePFBTagDeepCSV_p71;
   Bool_t HLT_Mu12_DoublePFJets54MaxDeta1p6_DoublePFBTagDeepCSV_p71;
   Bool_t HLT_DoublePFJets40_PFBTagDeepCSV_p71;
   Bool_t HLT_DoublePFJets100_PFBTagDeepCSV_p71;
   Bool_t HLT_DoublePFJets200_PFBTagDeepCSV_p71;
   Bool_t HLT_DoublePFJets350_PFBTagDeepCSV_p71;
   Bool_t HLT_DoublePFJets116MaxDeta1p6_DoublePFBTagDeepCSV_p71;
   Bool_t HLT_DoublePFJets128MaxDeta1p6_DoublePFBTagDeepCSV_p71;
   Bool_t HLT_Mu12_DoublePFJets40_PFBTagDeepJet_p71;
   Bool_t HLT_Mu12_DoublePFJets100_PFBTagDeepJet_p71;
   Bool_t HLT_Mu12_DoublePFJets200_PFBTagDeepJet_p71;
   Bool_t HLT_Mu12_DoublePFJets350_PFBTagDeepJet_p71;
   Bool_t HLT_Mu12_DoublePFJets40MaxDeta1p6_DoublePFBTagDeepJet_p71;
   Bool_t HLT_Mu12_DoublePFJets54MaxDeta1p6_DoublePFBTagDeepJet_p71;
   Bool_t HLT_DoublePFJets40_PFBTagDeepJet_p71;
   Bool_t HLT_DoublePFJets100_PFBTagDeepJet_p71;
   Bool_t HLT_DoublePFJets200_PFBTagDeepJet_p71;
   Bool_t HLT_DoublePFJets350_PFBTagDeepJet_p71;
   Bool_t HLT_DoublePFJets116MaxDeta1p6_DoublePFBTagDeepJet_p71;
   Bool_t HLT_DoublePFJets128MaxDeta1p6_DoublePFBTagDeepJet_p71;
   Bool_t HLT_Photon300_NoHE;
   Bool_t HLT_Mu8_TrkIsoVVL;
   Bool_t HLT_Mu8_DiEle12_CaloIdL_TrackIdL_DZ;
   Bool_t HLT_Mu8_DiEle12_CaloIdL_TrackIdL;
   Bool_t HLT_Mu8_Ele8_CaloIdM_TrackIdM_Mass8_PFHT350_DZ;
   Bool_t HLT_Mu8_Ele8_CaloIdM_TrackIdM_Mass8_PFHT350;
   Bool_t HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ;
   Bool_t HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_PFDiJet30;
   Bool_t HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_CaloDiJet30;
   Bool_t HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_PFDiJet30_PFBtagDeepCSV_1p5;
   Bool_t HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_CaloDiJet30_CaloBtagDeepCSV_1p5;
   Bool_t HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL;
   Bool_t HLT_Mu17_TrkIsoVVL;
   Bool_t HLT_Mu19_TrkIsoVVL;
   Bool_t HLT_BTagMu_AK4DiJet20_Mu5;
   Bool_t HLT_BTagMu_AK4DiJet40_Mu5;
   Bool_t HLT_BTagMu_AK4DiJet70_Mu5;
   Bool_t HLT_BTagMu_AK4DiJet110_Mu5;
   Bool_t HLT_BTagMu_AK4DiJet170_Mu5;
   Bool_t HLT_BTagMu_AK4Jet300_Mu5;
   Bool_t HLT_BTagMu_AK8DiJet170_Mu5;
   Bool_t HLT_BTagMu_AK8Jet170_DoubleMu5;
   Bool_t HLT_BTagMu_AK8Jet300_Mu5;
   Bool_t HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ;
   Bool_t HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL;
   Bool_t HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ;
   Bool_t HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL;
   Bool_t HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL;
   Bool_t HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ;
   Bool_t HLT_Photon20;
   Bool_t HLT_Photon33;
   Bool_t HLT_Photon50;
   Bool_t HLT_Photon75;
   Bool_t HLT_Photon90;
   Bool_t HLT_Photon120;
   Bool_t HLT_Photon150;
   Bool_t HLT_Photon175;
   Bool_t HLT_Photon200;
   Bool_t HLT_Photon30EB_TightID_TightIso;
   Bool_t HLT_Photon110EB_TightID_TightIso;
   Bool_t HLT_Photon100EBHE10;
   Bool_t HLT_Photon50_R9Id90_HE10_IsoM;
   Bool_t HLT_Photon75_R9Id90_HE10_IsoM;
   Bool_t HLT_Photon90_R9Id90_HE10_IsoM;
   Bool_t HLT_Photon120_R9Id90_HE10_IsoM;
   Bool_t HLT_Photon165_R9Id90_HE10_IsoM;
   Bool_t HLT_Diphoton30_22_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass90;
   Bool_t HLT_Diphoton30_22_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass95;
   Bool_t HLT_Photon35_TwoProngs35;
   Bool_t HLT_IsoMu24_TwoProngs35;
   Bool_t HLT_Dimuon0_Jpsi_L1_NoOS;
   Bool_t HLT_Dimuon0_Jpsi_NoVertexing_NoOS;
   Bool_t HLT_Dimuon0_Jpsi;
   Bool_t HLT_Dimuon0_Jpsi_NoVertexing;
   Bool_t HLT_Dimuon0_Jpsi_L1_4R_0er1p5R;
   Bool_t HLT_Dimuon0_Jpsi_NoVertexing_L1_4R_0er1p5R;
   Bool_t HLT_Dimuon0_Jpsi3p5_Muon2;
   Bool_t HLT_Dimuon0_Upsilon_L1_4p5;
   Bool_t HLT_Dimuon0_Upsilon_L1_5;
   Bool_t HLT_Dimuon0_Upsilon_L1_4p5NoOS;
   Bool_t HLT_Dimuon0_Upsilon_L1_4p5er2p0;
   Bool_t HLT_Dimuon0_Upsilon_L1_4p5er2p0M;
   Bool_t HLT_Dimuon0_Upsilon_NoVertexing;
   Bool_t HLT_Dimuon0_Upsilon_L1_5M;
   Bool_t HLT_Dimuon0_LowMass_L1_0er1p5R;
   Bool_t HLT_Dimuon0_LowMass_L1_0er1p5;
   Bool_t HLT_Dimuon0_LowMass;
   Bool_t HLT_Dimuon0_LowMass_L1_4;
   Bool_t HLT_Dimuon0_LowMass_L1_4R;
   Bool_t HLT_Dimuon0_LowMass_L1_TM530;
   Bool_t HLT_Dimuon0_Upsilon_Muon_L1_TM0;
   Bool_t HLT_Dimuon0_Upsilon_Muon_NoL1Mass;
   Bool_t HLT_TripleMu_5_3_3_Mass3p8_DZ;
   Bool_t HLT_TripleMu_10_5_5_DZ;
   Bool_t HLT_TripleMu_12_10_5;
   Bool_t HLT_Tau3Mu_Mu7_Mu1_TkMu1_Tau15;
   Bool_t HLT_Tau3Mu_Mu7_Mu1_TkMu1_Tau15_Charge1;
   Bool_t HLT_Tau3Mu_Mu7_Mu1_TkMu1_IsoTau15;
   Bool_t HLT_Tau3Mu_Mu7_Mu1_TkMu1_IsoTau15_Charge1;
   Bool_t HLT_DoubleMu3_DZ_PFMET50_PFMHT60;
   Bool_t HLT_DoubleMu3_DZ_PFMET70_PFMHT70;
   Bool_t HLT_DoubleMu3_DZ_PFMET90_PFMHT90;
   Bool_t HLT_DoubleMu3_Trk_Tau3mu_NoL1Mass;
   Bool_t HLT_DoubleMu4_Jpsi_Displaced;
   Bool_t HLT_DoubleMu4_Jpsi_NoVertexing;
   Bool_t HLT_DoubleMu4_JpsiTrkTrk_Displaced;
   Bool_t HLT_DoubleMu4_JpsiTrk_Bc;
   Bool_t HLT_DoubleMu43NoFiltersNoVtx;
   Bool_t HLT_DoubleMu48NoFiltersNoVtx;
   Bool_t HLT_Mu43NoFiltersNoVtx_Photon43_CaloIdL;
   Bool_t HLT_Mu48NoFiltersNoVtx_Photon48_CaloIdL;
   Bool_t HLT_Mu38NoFiltersNoVtxDisplaced_Photon38_CaloIdL;
   Bool_t HLT_Mu43NoFiltersNoVtxDisplaced_Photon43_CaloIdL;
   Bool_t HLT_HT425;
   Bool_t HLT_HT430_DisplacedDijet40_DisplacedTrack;
   Bool_t HLT_HT500_DisplacedDijet40_DisplacedTrack;
   Bool_t HLT_HT430_DisplacedDijet60_DisplacedTrack;
   Bool_t HLT_HT400_DisplacedDijet40_DisplacedTrack;
   Bool_t HLT_HT650_DisplacedDijet60_Inclusive;
   Bool_t HLT_HT550_DisplacedDijet60_Inclusive;
   Bool_t HLT_DiJet110_35_Mjj650_PFMET110;
   Bool_t HLT_DiJet110_35_Mjj650_PFMET120;
   Bool_t HLT_DiJet110_35_Mjj650_PFMET130;
   Bool_t HLT_TripleJet110_35_35_Mjj650_PFMET110;
   Bool_t HLT_TripleJet110_35_35_Mjj650_PFMET120;
   Bool_t HLT_TripleJet110_35_35_Mjj650_PFMET130;
   Bool_t HLT_Ele30_eta2p1_WPTight_Gsf_CentralPFJet35_EleCleaned;
   Bool_t HLT_Ele28_eta2p1_WPTight_Gsf_HT150;
   Bool_t HLT_Ele28_HighEta_SC20_Mass55;
   Bool_t HLT_Ele15_IsoVVVL_PFHT450_CaloBTagDeepCSV_4p5;
   Bool_t HLT_Ele15_IsoVVVL_PFHT450_PFMET50;
   Bool_t HLT_Ele15_IsoVVVL_PFHT450;
   Bool_t HLT_Ele50_IsoVVVL_PFHT450;
   Bool_t HLT_Ele15_IsoVVVL_PFHT600;
   Bool_t HLT_Mu4_TrkIsoVVL_DiPFJet90_40_DEta3p5_MJJ750_HTT300_PFMETNoMu60;
   Bool_t HLT_Mu8_TrkIsoVVL_DiPFJet40_DEta3p5_MJJ750_HTT300_PFMETNoMu60;
   Bool_t HLT_Mu10_TrkIsoVVL_DiPFJet40_DEta3p5_MJJ750_HTT350_PFMETNoMu60;
   Bool_t HLT_Mu15_IsoVVVL_PFHT450_CaloBTagDeepCSV_4p5;
   Bool_t HLT_Mu15_IsoVVVL_PFHT450_PFMET50;
   Bool_t HLT_Mu15_IsoVVVL_PFHT450;
   Bool_t HLT_Mu50_IsoVVVL_PFHT450;
   Bool_t HLT_Mu15_IsoVVVL_PFHT600;
   Bool_t HLT_Mu3er1p5_PFJet100er2p5_PFMET70_PFMHT70_IDTight;
   Bool_t HLT_Mu3er1p5_PFJet100er2p5_PFMET80_PFMHT80_IDTight;
   Bool_t HLT_Mu3er1p5_PFJet100er2p5_PFMET90_PFMHT90_IDTight;
   Bool_t HLT_Mu3er1p5_PFJet100er2p5_PFMET100_PFMHT100_IDTight;
   Bool_t HLT_Mu3er1p5_PFJet100er2p5_PFMETNoMu70_PFMHTNoMu70_IDTight;
   Bool_t HLT_Mu3er1p5_PFJet100er2p5_PFMETNoMu80_PFMHTNoMu80_IDTight;
   Bool_t HLT_Mu3er1p5_PFJet100er2p5_PFMETNoMu90_PFMHTNoMu90_IDTight;
   Bool_t HLT_Mu3er1p5_PFJet100er2p5_PFMETNoMu100_PFMHTNoMu100_IDTight;
   Bool_t HLT_Dimuon10_PsiPrime_Barrel_Seagulls;
   Bool_t HLT_Dimuon20_Jpsi_Barrel_Seagulls;
   Bool_t HLT_Dimuon10_Upsilon_y1p4;
   Bool_t HLT_Dimuon12_Upsilon_y1p4;
   Bool_t HLT_Dimuon14_Phi_Barrel_Seagulls;
   Bool_t HLT_Dimuon25_Jpsi;
   Bool_t HLT_Dimuon14_PsiPrime;
   Bool_t HLT_Dimuon14_PsiPrime_noCorrL1;
   Bool_t HLT_Dimuon18_PsiPrime;
   Bool_t HLT_Dimuon18_PsiPrime_noCorrL1;
   Bool_t HLT_Dimuon24_Upsilon_noCorrL1;
   Bool_t HLT_Dimuon24_Phi_noCorrL1;
   Bool_t HLT_Dimuon25_Jpsi_noCorrL1;
   Bool_t HLT_DiMu4_Ele9_CaloIdL_TrackIdL_DZ_Mass3p8;
   Bool_t HLT_DiMu9_Ele9_CaloIdL_TrackIdL_DZ;
   Bool_t HLT_DiMu9_Ele9_CaloIdL_TrackIdL;
   Bool_t HLT_DoubleIsoMu20_eta2p1;
   Bool_t HLT_TrkMu12_DoubleTrkMu5NoFiltersNoVtx;
   Bool_t HLT_Mu8;
   Bool_t HLT_Mu17;
   Bool_t HLT_Mu19;
   Bool_t HLT_Mu17_Photon30_IsoCaloId;
   Bool_t HLT_Ele8_CaloIdL_TrackIdL_IsoVL_PFJet30;
   Bool_t HLT_Ele12_CaloIdL_TrackIdL_IsoVL_PFJet30;
   Bool_t HLT_Ele23_CaloIdL_TrackIdL_IsoVL_PFJet30;
   Bool_t HLT_Ele8_CaloIdM_TrackIdM_PFJet30;
   Bool_t HLT_Ele17_CaloIdM_TrackIdM_PFJet30;
   Bool_t HLT_Ele23_CaloIdM_TrackIdM_PFJet30;
   Bool_t HLT_Ele50_CaloIdVT_GsfTrkIdT_PFJet165;
   Bool_t HLT_Ele115_CaloIdVT_GsfTrkIdT;
   Bool_t HLT_Ele135_CaloIdVT_GsfTrkIdT;
   Bool_t HLT_PFHT330PT30_QuadPFJet_75_60_45_40_TriplePFBTagDeepCSV_4p5;
   Bool_t HLT_PFHT330PT30_QuadPFJet_75_60_45_40;
   Bool_t HLT_PFHT400_SixPFJet32_DoublePFBTagDeepCSV_2p94;
   Bool_t HLT_PFHT400_SixPFJet32;
   Bool_t HLT_PFHT450_SixPFJet36_PFBTagDeepCSV_1p59;
   Bool_t HLT_PFHT450_SixPFJet36;
   Bool_t HLT_PFHT400_FivePFJet_100_100_60_30_30;
   Bool_t HLT_PFHT400_FivePFJet_100_100_60_30_30_DoublePFBTagDeepCSV_4p5;
   Bool_t HLT_PFHT400_FivePFJet_120_120_60_30_30_DoublePFBTagDeepCSV_4p5;
   Bool_t HLT_PFHT350;
   Bool_t HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_PFHT350MinPFJet15;
   Bool_t HLT_ECALHT800;
   Bool_t HLT_DiSC30_18_EIso_AND_HE_Mass70;
   Bool_t HLT_Physics;
   Bool_t HLT_Random;
   Bool_t HLT_ZeroBias;
   Bool_t HLT_ZeroBias_Alignment;
   Bool_t HLT_Photon20_HoverELoose;
   Bool_t HLT_Photon30_HoverELoose;
   Bool_t HLT_EcalCalibration;
   Bool_t HLT_HcalCalibration;
   Bool_t HLT_L1UnpairedBunchBptxMinus;
   Bool_t HLT_L1UnpairedBunchBptxPlus;
   Bool_t HLT_L1NotBptxOR;
   Bool_t HLT_L1_CDC_SingleMu_3_er1p2_TOP120_DPHI2p618_3p142;
   Bool_t HLT_CDC_L2cosmic_10_er1p0;
   Bool_t HLT_CDC_L2cosmic_5p5_er1p0;
   Bool_t HLT_HcalNZS;
   Bool_t HLT_HcalPhiSym;
   Bool_t HLT_HcalIsolatedbunch;
   Bool_t HLT_IsoTrackHB;
   Bool_t HLT_IsoTrackHE;
   Bool_t HLT_ZeroBias_FirstCollisionAfterAbortGap;
   Bool_t HLT_ZeroBias_IsolatedBunches;
   Bool_t HLT_ZeroBias_FirstCollisionInTrain;
   Bool_t HLT_ZeroBias_LastCollisionInTrain;
   Bool_t HLT_ZeroBias_FirstBXAfterTrain;
   Bool_t HLT_Ele16_Ele12_Ele8_CaloIdL_TrackIdL;
   Bool_t HLT_VBF_DoubleMediumDeepTauPFTauHPS20_eta2p1;
   Bool_t HLT_Photon75_R9Id90_HE10_IsoM_EBOnly_PFJetsMJJ300DEta3;
   Bool_t HLT_PFMET100_PFMHT100_IDTight_PFHT60;
   Bool_t HLT_PFMETNoMu100_PFMHTNoMu100_IDTight_PFHT60;
   Bool_t HLT_PFMETTypeOne100_PFMHT100_IDTight_PFHT60;
   Bool_t HLT_Mu18_Mu9_SameSign;
   Bool_t HLT_DoubleMu2_Jpsi_DoubleTrk1_Phi1p05;
   Bool_t HLT_DoubleMu3_DCA_PFMET50_PFMHT60;
   Bool_t HLT_TripleMu_5_3_3_Mass3p8_DCA;
   Bool_t HLT_QuadPFJet103_88_75_15_DoublePFBTagDeepCSV_1p3_7p7_VBF1;
   Bool_t HLT_QuadPFJet111_90_80_15_DoublePFBTagDeepCSV_1p3_7p7_VBF1;
   Bool_t HLT_QuadPFJet103_88_75_15_PFBTagDeepCSV_1p3_VBF2;
   Bool_t HLT_QuadPFJet105_88_76_15_PFBTagDeepCSV_1p3_VBF2;
   Bool_t HLT_QuadPFJet111_90_80_15_PFBTagDeepCSV_1p3_VBF2;
   Bool_t HLT_QuadPFJet103_88_75_15;
   Bool_t HLT_QuadPFJet105_88_76_15;
   Bool_t HLT_QuadPFJet111_90_80_15;
   Bool_t HLT_AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_p02;
   Bool_t HLT_AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_np2;
   Bool_t HLT_AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_np4;
   Bool_t HLT_Diphoton30_18_R9IdL_AND_HE_AND_IsoCaloId_Mass55;
   Bool_t HLT_Diphoton30_18_R9IdL_AND_HE_AND_IsoCaloId;
   Bool_t HLT_Mu12_IP6;
   Bool_t HLT_QuadPFJet105_88_76_15_DoublePFBTagDeepCSV_1p3_7p7_VBF1;
   Bool_t HLT_DoubleMediumDeepTauPFTauHPS35_L2NN_eta2p1;
   Bool_t HLT_Ele24_eta2p1_WPTight_Gsf_LooseDeepTauPFTauHPS30_eta2p1_CrossL1;
   Bool_t HLT_IsoMu20_eta2p1_LooseDeepTauPFTauHPS27_eta2p1_CrossL1;
   Bool_t HLT_IsoMu24_eta2p1_LooseDeepTauPFTauHPS180_eta2p1;
   Bool_t HLT_IsoMu24_eta2p1_LooseDeepTauPFTauHPS30_eta2p1_CrossL1;
   Bool_t HLT_IsoMu24_eta2p1_MediumDeepTauPFTauHPS35_L2NN_eta2p1_CrossL1;
   Bool_t HLT_LooseDeepTauPFTauHPS180_L2NN_eta2p1;
   Bool_t HLT_PFHT330PT30_QuadPFJet_75_60_45_40_TriplePFBTagDeepJet_4p5;
   Bool_t HLT_PFHT400_FivePFJet_100_100_60_30_30_DoublePFBTagDeepJet_4p5;
   Bool_t HLT_PFHT400_FivePFJet_120_120_60_30_30_DoublePFBTagDeepJet_4p5;
   Bool_t HLT_PFHT400_SixPFJet32_DoublePFBTagDeepJet_2p94;
   Bool_t HLT_PFHT450_SixPFJet36_PFBTagDeepJet_1p59;
   Bool_t HLT_QuadPFJet103_88_75_15_DoublePFBTagDeepJet_1p3_7p7_VBF1;
   Bool_t HLT_QuadPFJet103_88_75_15_PFBTagDeepJet_1p3_VBF2;
   Bool_t HLT_QuadPFJet105_88_76_15_DoublePFBTagDeepJet_1p3_7p7_VBF1;
   Bool_t HLT_QuadPFJet105_88_76_15_PFBTagDeepJet_1p3_VBF2;
   Bool_t HLT_QuadPFJet111_90_80_15_DoublePFBTagDeepJet_1p3_7p7_VBF1;
   Bool_t HLT_QuadPFJet111_90_80_15_PFBTagDeepJet_1p3_VBF2;
   Bool_t HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_PFDiJet30_PFBtagDeepJet_1p5;
   Bool_t HLT_QuadPFJet70_50_40_30;
   Bool_t HLT_QuadPFJet70_50_40_30_PFBTagParticleNet_2BTagSum0p65;
   Bool_t HLT_QuadPFJet70_50_40_35_PFBTagParticleNet_2BTagSum0p65;
   Bool_t HLT_QuadPFJet70_50_45_35_PFBTagParticleNet_2BTagSum0p65;
   Bool_t HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_PFDiJet30_PFBTagParticleNet_2BTagSum0p65;
   Bool_t HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_QuadPFJet70_50_40_30;
   Bool_t HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_QuadPFJet70_50_40_30_PFBTagParticleNet_2BTagSum0p65;
   Bool_t HLT_AK8PFJet230_SoftDropMass40;
   Bool_t HLT_AK8PFJet230_SoftDropMass40_PFAK8ParticleNetBB0p35;
   Bool_t HLT_AK8PFJet250_SoftDropMass40_PFAK8ParticleNetBB0p35;
   Bool_t HLT_AK8PFJet275_SoftDropMass40_PFAK8ParticleNetBB0p35;
   Bool_t HLT_AK8PFJet400_SoftDropMass40;
   Bool_t HLT_AK8PFJet425_SoftDropMass40;
   Bool_t HLT_AK8PFJet450_SoftDropMass40;
   Bool_t HLT_AK8PFJet230_SoftDropMass40_PFAK8ParticleNetTauTau0p30;
   Bool_t HLT_AK8PFJet250_SoftDropMass40_PFAK8ParticleNetTauTau0p30;
   Bool_t HLT_AK8PFJet275_SoftDropMass40_PFAK8ParticleNetTauTau0p30;
   Bool_t HLT_IsoMu50_AK8PFJet230_SoftDropMass40;
   Bool_t HLT_IsoMu50_AK8PFJet230_SoftDropMass40_PFAK8ParticleNetBB0p35;
   Bool_t HLT_Ele50_CaloIdVT_GsfTrkIdT_AK8PFJet230_SoftDropMass40;
   Bool_t HLT_Ele50_CaloIdVT_GsfTrkIdT_AK8PFJet230_SoftDropMass40_PFAK8ParticleNetBB0p35;
   Bool_t HLT_DoubleMediumDeepTauPFTauHPS30_L2NN_eta2p1_PFJet60;
   Bool_t HLT_DoubleMediumDeepTauPFTauHPS30_L2NN_eta2p1_PFJet75;
   Bool_t HLT_IsoMu24_eta2p1_MediumDeepTauPFTauHPS30_L2NN_eta2p1_CrossL1;
   Bool_t HLT_IsoMu24_eta2p1_MediumDeepTauPFTauHPS30_L2NN_eta2p1_PFJet60_CrossL1;
   Bool_t HLT_IsoMu24_eta2p1_MediumDeepTauPFTauHPS30_L2NN_eta2p1_PFJet75_CrossL1;
   Bool_t HLT_DoublePFJets40_Mass500_MediumDeepTauPFTauHPS45_L2NN_MediumDeepTauPFTauHPS20_eta2p1;
   Bool_t HLT_IsoMu24_eta2p1_MediumDeepTauPFTauHPS20_eta2p1_SingleL1;
   Bool_t HLT_IsoMu24_eta2p1_MediumDeepTauPFTauHPS45_L2NN_eta2p1_CrossL1;
   Bool_t HLT_DoubleL2Mu10NoVtx_2Cha_VetoL3Mu0DxyMax1cm;
   Bool_t HLT_DoubleL2Mu12NoVtx_2Cha_VetoL3Mu0DxyMax1cm;
   Bool_t HLT_DoubleL2Mu14NoVtx_2Cha_VetoL3Mu0DxyMax1cm;
   Bool_t HLT_DoubleL3Mu16_10NoVtx_DxyMin0p01cm;
   Bool_t HLT_DoubleL3Mu18_10NoVtx_DxyMin0p01cm;
   Bool_t HLT_DoubleL3Mu20_10NoVtx_DxyMin0p01cm;
   Bool_t HLT_L2Mu10NoVtx_2Cha;
   Bool_t HLT_L2Mu10NoVtx_2Cha_VetoL3Mu0DxyMax1cm;
   Bool_t HLT_L3Mu10NoVtx;
   Bool_t HLT_L3Mu10NoVtx_DxyMin0p01cm;
   Bool_t HLT_DoubleL2Mu_L3Mu16NoVtx_VetoL3Mu0DxyMax0p1cm;
   Bool_t HLT_DoubleL2Mu_L3Mu18NoVtx_VetoL3Mu0DxyMax0p1cm;
   Bool_t HLT_DoubleL2Mu10NoVtx_2Cha_CosmicSeed_VetoL3Mu0DxyMax1cm;
   Bool_t HLT_DoubleL2Mu12NoVtx_2Cha_CosmicSeed_VetoL3Mu0DxyMax1cm;
   Bool_t HLT_L2Mu10NoVtx_2Cha_CosmicSeed;
   Bool_t HLT_L2Mu10NoVtx_2Cha_CosmicSeed_VetoL3Mu0DxyMax1cm;
   Bool_t HLT_DoubleL3dTksMu16_10NoVtx_DxyMin0p01cm;
   Bool_t HLT_L3dTksMu10_NoVtx_DxyMin0p01cm;
   Bool_t HLT_Mu20NoFiltersNoVtxDisplaced_Photon20_CaloCustomId;
   Bool_t HLT_DoubleMediumChargedIsoDisplacedPFTauHPS32_Trk1_eta2p1;
   Bool_t HLT_HT430_DelayedJet40_DoubleDelay0p5nsTrackless;
   Bool_t HLT_HT430_DelayedJet40_DoubleDelay1nsInclusive;
   Bool_t HLT_HT430_DelayedJet40_SingleDelay1nsTrackless;
   Bool_t HLT_HT430_DelayedJet40_SingleDelay2nsInclusive;
   Bool_t HLT_L1Mu6HT240;
   Bool_t HLT_Mu6HT240_DisplacedDijet30_Inclusive0PtrkShortSig5;
   Bool_t HLT_Mu6HT240_DisplacedDijet30_Inclusive1PtrkShortSig5_DisplacedLoose;
   Bool_t HLT_Mu6HT240_DisplacedDijet35_Inclusive0PtrkShortSig5;
   Bool_t HLT_Mu6HT240_DisplacedDijet35_Inclusive1PtrkShortSig5_DisplacedLoose;
   Bool_t HLT_Mu6HT240_DisplacedDijet40_Inclusive0PtrkShortSig5;
   Bool_t HLT_Mu6HT240_DisplacedDijet40_Inclusive1PtrkShortSig5_DisplacedLoose;
   Bool_t HLT_HT430_DisplacedDijet30_Inclusive1PtrkShortSig5;
   Bool_t HLT_HT430_DisplacedDijet35_Inclusive1PtrkShortSig5;
   Bool_t HLT_HT430_DisplacedDijet40_Inclusive1PtrkShortSig5;
   Bool_t HLT_CaloMET60_DTCluster50;
   Bool_t HLT_CaloMET60_DTClusterNoMB1S50;
   Bool_t HLT_L1MET_DTCluster50;
   Bool_t HLT_L1MET_DTClusterNoMB1S50;
   Bool_t HLT_CscCluster_Loose;
   Bool_t HLT_CscCluster_Medium;
   Bool_t HLT_CscCluster_Tight;
   Bool_t HLT_L1CSCShower_DTCluster50;
   Bool_t HLT_L1CSCShower_DTCluster75;
   Bool_t HLT_DoubleCscCluster75;
   Bool_t  HLT_DoubleCscCluster100;
   Bool_t HLT_PFMET105_IsoTrk50;
   Bool_t HLT_PFMET110_PFJet100;
   Bool_t HLT_HT170_L1SingleLLPJet_DisplacedDijet40_DisplacedTrack;
   Bool_t HLT_HT200_L1SingleLLPJet_DisplacedDijet40_DisplacedTrack;
   Bool_t HLT_HT200_L1SingleLLPJet_DisplacedDijet60_DisplacedTrack;
   Bool_t HLT_HT270_L1SingleLLPJet_DisplacedDijet40_DisplacedTrack;
   Bool_t HLT_HT320_L1SingleLLPJet_DisplacedDijet60_Inclusive;
   Bool_t HLT_HT420_L1SingleLLPJet_DisplacedDijet60_Inclusive;
   Bool_t HLT_HT200_L1SingleLLPJet_DelayedJet40_SingleDelay1nsTrackless;
   Bool_t HLT_HT200_L1SingleLLPJet_DelayedJet40_SingleDelay2nsInclusive;
   Bool_t HLT_HT200_L1SingleLLPJet_DelayedJet40_DoubleDelay0p5nsTrackless;
   Bool_t HLT_HT200_L1SingleLLPJet_DelayedJet40_DoubleDelay1nsInclusive;
   Bool_t HLT_HT200_L1SingleLLPJet_DisplacedDijet30_Inclusive1PtrkShortSig5;
   Bool_t HLT_HT200_L1SingleLLPJet_DisplacedDijet35_Inclusive1PtrkShortSig5;
   Bool_t HLT_HT200_L1SingleLLPJet_DisplacedDijet40_Inclusive1PtrkShortSig5;
   Bool_t HLT_DiPhoton10Time1p4ns;
   Bool_t HLT_DiPhoton10Time1p6ns;
   Bool_t HLT_DiPhoton10Time1p8ns;
   Bool_t HLT_DiPhoton10Time2ns;
   Bool_t HLT_DiPhoton10sminlt0p1;
   Bool_t HLT_DiPhoton10sminlt0p12;
   Bool_t HLT_DiPhoton10_CaloIdL;
   Bool_t HLT_DoubleEle4_eta1p22_mMax6;
   Bool_t HLT_DoubleEle4p5_eta1p22_mMax6;
   Bool_t HLT_DoubleEle5_eta1p22_mMax6;
   Bool_t HLT_DoubleEle5p5_eta1p22_mMax6;
   Bool_t HLT_DoubleEle6_eta1p22_mMax6;
   Bool_t HLT_DoubleEle6p5_eta1p22_mMax6;
   Bool_t HLT_DoubleEle7_eta1p22_mMax6;
   Bool_t HLT_DoubleEle7p5_eta1p22_mMax6;
   Bool_t HLT_DoubleEle8_eta1p22_mMax6;
   Bool_t HLT_DoubleEle8p5_eta1p22_mMax6;
   Bool_t HLT_DoubleEle9_eta1p22_mMax6;
   Bool_t HLT_DoubleEle9p5_eta1p22_mMax6;
   Bool_t HLT_DoubleEle10_eta1p22_mMax6;
   Bool_t HLT_Diphoton20_14_eta1p5_R9IdL_AND_HE_AND_IsoTCaloIdT;
   Bool_t HLT_Diphoton20_14_eta1p5_R9IdL_AND_HET_AND_IsoTCaloIdT;
   Bool_t HLT_Diphoton22_14_eta1p5_R9IdL_AND_HE_AND_IsoTCaloIdT;
   Bool_t HLT_Diphoton22_14_eta1p5_R9IdL_AND_HET_AND_IsoTCaloIdT;
   Bool_t HLT_ExpressMuons;
   Bool_t HLT_OnlineMonitorGroup;
   Bool_t HLT_PPSMaxTracksPerArm1;
   Bool_t HLT_PPSMaxTracksPerRP4;
   Bool_t HLTriggerFinalPath;
   Bool_t HLT_EphemeralPhysics;
   Bool_t HLT_EphemeralZeroBias;
   Bool_t L1_DoubleMu3_SQ_ETMHF30_HTT60er;
   Bool_t L1_DoubleMu3_SQ_ETMHF30_Jet60er2p5_OR_DoubleJet40er2p5;
   Bool_t L1_DoubleMu3_SQ_ETMHF40_HTT60er;
   Bool_t L1_DoubleMu3_SQ_ETMHF40_Jet60er2p5_OR_DoubleJet40er2p5;
   Bool_t L1_ETMHF80_SingleJet55er2p5_dPhi_Min2p1;
   Bool_t L1_ETMHF80_SingleJet55er2p5_dPhi_Min2p6;
   Bool_t L1_Mu3er1p5_Jet100er2p5_ETMHF30;
   Int_t nCscSeg;
   Int_t nCscRechits;
   Int_t nDtSeg;
   Int_t nDtRechits;
   Int_t nRpc;
   Float_t cscRechitsPhi[90000];      //[nCscRechits]
   Float_t cscRechitsEta[90000];      //[nCscRechits]
   Float_t cscRechitsX[90000];        //[nCscRechits]
   Float_t cscRechitsY[90000];        //[nCscRechits]
   Float_t cscRechitsZ[90000];        //[nCscRechits]
   Float_t cscRechitsTpeak[90000];    //[nCscRechits]
   Float_t cscRechitsTwire[90000];    //[nCscRechits]
   Float_t dtRechitCorrectX[20000];   //[nDtRechits]
   Float_t dtRechitCorrectY[20000];   //[nDtRechits]
   Float_t dtRechitCorrectZ[20000];   //[nDtRechits]
   Float_t dtRechitCorrectEta[20000]; //[nDtRechits]
   Float_t dtRechitCorrectPhi[20000]; //[nDtRechits]
   Float_t cscSegPhi[512];            //[nCscSeg]
   Float_t cscSegEta[512];            //[nCscSeg]
   Float_t dtSegPhi[1500];             //[nDtSeg]
   Float_t dtSegEta[1500];             //[nDtSeg]
   Float_t rpcPhi[2000];              //[nRpc]
   Float_t rpcEta[2000];              //[nRpc]
   Float_t rpcX[2000];                //[nRpc]
   Float_t rpcY[2000];                //[nRpc]
   Float_t rpcZ[2000];                //[nRpc]
   Float_t rpcT[2000];                //[nRpc]
   Float_t rpcTError[2000];           //[nRpc]
   Int_t cscRechitsChamber[90000];    //[nCscRechits]
   Int_t cscRechitsStation[90000];    //[nCscRechits]
   Int_t cscRechitsDetId[90000];      //[nCscRechits]
   Int_t dtRechitStation[20000];      //[nDtRechits]
   Int_t dtRechitWheel[20000];        //[nDtRechits]
   Int_t dtRechitSuperLayer[20000];   //[nDtRechits]
   Int_t cscSegChamber[512];          //[nCscSeg]
   Int_t cscSegStation[512];          //[nCscSeg]
   Int_t cscSegNRecHits[512];         //[nCscSeg]
   Int_t dtSegStation[1500];           //[nDtSeg]
   Int_t dtSegWheel[1500];             //[nDtSeg]
   Int_t rpcBx[2000];                 //[nRpc]
   Int_t rpcRegion[2000];             //[nRpc]
   Int_t rpcRing[2000];               //[nRpc]
   Int_t rpcSector[2000];             //[nRpc]
   Int_t rpcStation[2000];            //[nRpc]
   Int_t rpcLayer[2000];              //[nRpc]

   // List of branches
   TBranch *b_nsoftActivityVH;            //!
   TBranch *b_nGenJetAK8;                 //!
   TBranch *b_GenJetAK8_eta;              //!
   TBranch *b_GenJetAK8_mass;             //!
   TBranch *b_GenJetAK8_phi;              //!
   TBranch *b_GenJetAK8_pt;               //!
   TBranch *b_nGenJet;                    //!
   TBranch *b_GenJet_eta;                 //!
   TBranch *b_GenJet_mass;                //!
   TBranch *b_GenJet_phi;                 //!
   TBranch *b_GenJet_pt;                  //!
   TBranch *b_nGenPart;                   //!
   TBranch *b_GenPart_genPartIdxMother;   //!
   TBranch *b_GenPart_statusFlags;        //!
   TBranch *b_GenPart_pdgId;              //!
   TBranch *b_GenPart_status;             //!
   TBranch *b_GenPart_eta;                //!
   TBranch *b_GenPart_mass;               //!
   TBranch *b_GenPart_phi;                //!
   TBranch *b_GenPart_pt;                 //!
   TBranch *b_nGenProton;                 //!
   TBranch *b_GenProton_isPU;             //!
   TBranch *b_GenProton_px;               //!
   TBranch *b_GenProton_py;               //!
   TBranch *b_GenProton_pz;               //!
   TBranch *b_GenProton_vz;               //!
   TBranch *b_nSubGenJetAK8;              //!
   TBranch *b_SubGenJetAK8_eta;           //!
   TBranch *b_SubGenJetAK8_mass;          //!
   TBranch *b_SubGenJetAK8_phi;           //!
   TBranch *b_SubGenJetAK8_pt;            //!
   TBranch *b_Generator_id1;              //!
   TBranch *b_Generator_id2;              //!
   TBranch *b_Generator_binvar;           //!
   TBranch *b_Generator_scalePDF;         //!
   TBranch *b_Generator_weight;           //!
   TBranch *b_Generator_x1;               //!
   TBranch *b_Generator_x2;               //!
   TBranch *b_Generator_xpdf1;            //!
   TBranch *b_Generator_xpdf2;            //!
   TBranch *b_GenVtx_x;                   //!
   TBranch *b_GenVtx_y;                   //!
   TBranch *b_GenVtx_z;                   //!
   TBranch *b_nGenVisTau;                 //!
   TBranch *b_GenVisTau_status;           //!
   TBranch *b_GenVisTau_charge;           //!
   TBranch *b_GenVisTau_genPartIdxMother; //!
   TBranch *b_GenVisTau_eta;              //!
   TBranch *b_GenVisTau_mass;             //!
   TBranch *b_GenVisTau_phi;              //!
   TBranch *b_GenVisTau_pt;               //!
   TBranch *b_genWeight;                  //!
   TBranch *b_GenMET_phi;                 //!
   TBranch *b_GenMET_pt;                  //!
   TBranch *b_nGenDressedLepton;          //!
   TBranch *b_GenDressedLepton_hasTauAnc; //!
   TBranch *b_GenDressedLepton_pdgId;     //!
   TBranch *b_GenDressedLepton_eta;       //!
   TBranch *b_GenDressedLepton_mass;      //!
   TBranch *b_GenDressedLepton_phi;       //!
   TBranch *b_GenDressedLepton_pt;        //!
   TBranch *b_nGenIsolatedPhoton;         //!
   TBranch *b_GenIsolatedPhoton_eta;      //!
   TBranch *b_GenIsolatedPhoton_mass;     //!
   TBranch *b_GenIsolatedPhoton_phi;      //!
   TBranch *b_GenIsolatedPhoton_pt;       //!
   TBranch *b_genTtbarId;                 //!
   TBranch *b_boostedTau_genPartFlav;     //!
   TBranch *b_boostedTau_genPartIdx;      //!
   TBranch *b_Electron_genPartFlav;       //!
   TBranch *b_Electron_genPartIdx;        //!
   TBranch *b_FatJet_genJetAK8Idx;        //!
   TBranch *b_GenJetAK8_hadronFlavour;    //!
   TBranch *b_GenJetAK8_partonFlavour;    //!
   TBranch *b_GenJet_hadronFlavour;       //!
   TBranch *b_GenJet_partonFlavour;       //!
   TBranch *b_GenVtx_t0;                  //!
   TBranch *b_Jet_genJetIdx;              //!
   TBranch *b_LowPtElectron_genPartFlav;  //!
   TBranch *b_LowPtElectron_genPartIdx;   //!
   TBranch *b_Muon_genPartFlav;           //!
   TBranch *b_Muon_genPartIdx;            //!
   TBranch *b_Photon_genPartFlav;         //!
   TBranch *b_Photon_genPartIdx;          //!
   TBranch *b_MET_fiducialGenPhi;         //!
   TBranch *b_MET_fiducialGenPt;          //!
   TBranch *b_Tau_genPartFlav;            //!
   TBranch *b_Tau_genPartIdx;             //!

   TBranch *b_gLLP_eta;
   TBranch *b_gLLP_phi;
   TBranch *b_gLLP_csc;
   TBranch *b_gLLP_dt;
   TBranch *b_gLLP_beta;
   TBranch *b_gLLP_e;
   TBranch *b_gLLP_pt;
   TBranch *b_gLLP_decay_vertex_x;
   TBranch *b_gLLP_decay_vertex_y;
   TBranch *b_gLLP_decay_vertex_z;
   TBranch *b_gLLP_decay_vertex_r;

   TBranch *b_run;                                                                                                   //!
   TBranch *b_luminosityBlock;                                                                                       //!
   TBranch *b_event;                                                                                                 //!
   TBranch *b_bunchCrossing;                                                                                         //!
   TBranch *b_BeamSpot_type;                                                                                         //!
   TBranch *b_BeamSpot_sigmaZ;                                                                                       //!
   TBranch *b_BeamSpot_sigmaZError;                                                                                  //!
   TBranch *b_BeamSpot_z;                                                                                            //!
   TBranch *b_BeamSpot_zError;                                                                                       //!
   TBranch *b_nboostedTau;                                                                                           //!
   TBranch *b_boostedTau_idAntiEle2018;                                                                              //!
   TBranch *b_boostedTau_idAntiMu;                                                                                   //!
   TBranch *b_boostedTau_idMVAnewDM2017v2;                                                                           //!
   TBranch *b_boostedTau_idMVAoldDM2017v2;                                                                           //!
   TBranch *b_boostedTau_jetIdx;                                                                                     //!
   TBranch *b_boostedTau_rawAntiEleCat2018;                                                                          //!
   TBranch *b_boostedTau_charge;                                                                                     //!
   TBranch *b_boostedTau_decayMode;                                                                                  //!
   TBranch *b_boostedTau_chargedIso;                                                                                 //!
   TBranch *b_boostedTau_eta;                                                                                        //!
   TBranch *b_boostedTau_leadTkDeltaEta;                                                                             //!
   TBranch *b_boostedTau_leadTkDeltaPhi;                                                                             //!
   TBranch *b_boostedTau_leadTkPtOverTauPt;                                                                          //!
   TBranch *b_boostedTau_mass;                                                                                       //!
   TBranch *b_boostedTau_neutralIso;                                                                                 //!
   TBranch *b_boostedTau_phi;                                                                                        //!
   TBranch *b_boostedTau_photonsOutsideSignalCone;                                                                   //!
   TBranch *b_boostedTau_pt;                                                                                         //!
   TBranch *b_boostedTau_puCorr;                                                                                     //!
   TBranch *b_boostedTau_rawAntiEle2018;                                                                             //!
   TBranch *b_boostedTau_rawIso;                                                                                     //!
   TBranch *b_boostedTau_rawIsodR03;                                                                                 //!
   TBranch *b_boostedTau_rawMVAnewDM2017v2;                                                                          //!
   TBranch *b_boostedTau_rawMVAoldDM2017v2;                                                                          //!
   TBranch *b_CaloMET_phi;                                                                                           //!
   TBranch *b_CaloMET_pt;                                                                                            //!
   TBranch *b_CaloMET_sumEt;                                                                                         //!
   TBranch *b_ChsMET_phi;                                                                                            //!
   TBranch *b_ChsMET_pt;                                                                                             //!
   TBranch *b_ChsMET_sumEt;                                                                                          //!
   TBranch *b_nCorrT1METJet;                                                                                         //!
   TBranch *b_CorrT1METJet_area;                                                                                     //!
   TBranch *b_CorrT1METJet_eta;                                                                                      //!
   TBranch *b_CorrT1METJet_muonSubtrFactor;                                                                          //!
   TBranch *b_CorrT1METJet_phi;                                                                                      //!
   TBranch *b_CorrT1METJet_rawPt;                                                                                    //!
   TBranch *b_DeepMETResolutionTune_phi;                                                                             //!
   TBranch *b_DeepMETResolutionTune_pt;                                                                              //!
   TBranch *b_DeepMETResponseTune_phi;                                                                               //!
   TBranch *b_DeepMETResponseTune_pt;                                                                                //!
   TBranch *b_nElectron;                                                                                             //!
   TBranch *b_Electron_seediEtaOriX;                                                                                 //!
   TBranch *b_Electron_convVeto;                                                                                     //!
   TBranch *b_Electron_cutBased;                                                                                     //!
   TBranch *b_Electron_cutBased_HEEP;                                                                                //!
   TBranch *b_Electron_isPFcand;                                                                                     //!
   TBranch *b_Electron_jetNDauCharged;                                                                               //!
   TBranch *b_Electron_lostHits;                                                                                     //!
   TBranch *b_Electron_mvaIso_WP80;                                                                                  //!
   TBranch *b_Electron_mvaIso_WP90;                                                                                  //!
   TBranch *b_Electron_mvaNoIso_WP80;                                                                                //!
   TBranch *b_Electron_mvaNoIso_WP90;                                                                                //!
   TBranch *b_Electron_seedGain;                                                                                     //!
   TBranch *b_Electron_tightCharge;                                                                                  //!
   TBranch *b_Electron_jetIdx;                                                                                       //!
   TBranch *b_Electron_photonIdx;                                                                                    //!
   TBranch *b_Electron_svIdx;                                                                                        //!
   TBranch *b_Electron_fsrPhotonIdx;                                                                                 //!
   TBranch *b_Electron_charge;                                                                                       //!
   TBranch *b_Electron_pdgId;                                                                                        //!
   TBranch *b_Electron_seediPhiOriY;                                                                                 //!
   TBranch *b_Electron_vidNestedWPBitmap;                                                                            //!
   TBranch *b_Electron_vidNestedWPBitmapHEEP;                                                                        //!
   TBranch *b_Electron_deltaEtaSC;                                                                                   //!
   TBranch *b_Electron_dr03EcalRecHitSumEt;                                                                          //!
   TBranch *b_Electron_dr03HcalDepth1TowerSumEt;                                                                     //!
   TBranch *b_Electron_dr03TkSumPt;                                                                                  //!
   TBranch *b_Electron_dr03TkSumPtHEEP;                                                                              //!
   TBranch *b_Electron_dxy;                                                                                          //!
   TBranch *b_Electron_dxyErr;                                                                                       //!
   TBranch *b_Electron_dz;                                                                                           //!
   TBranch *b_Electron_dzErr;                                                                                        //!
   TBranch *b_Electron_eInvMinusPInv;                                                                                //!
   TBranch *b_Electron_energyErr;                                                                                    //!
   TBranch *b_Electron_eta;                                                                                          //!
   TBranch *b_Electron_hoe;                                                                                          //!
   TBranch *b_Electron_ip3d;                                                                                         //!
   TBranch *b_Electron_jetPtRelv2;                                                                                   //!
   TBranch *b_Electron_jetRelIso;                                                                                    //!
   TBranch *b_Electron_mass;                                                                                         //!
   TBranch *b_Electron_miniPFRelIso_all;                                                                             //!
   TBranch *b_Electron_miniPFRelIso_chg;                                                                             //!
   TBranch *b_Electron_mvaHZZIso;                                                                                    //!
   TBranch *b_Electron_mvaIso;                                                                                       //!
   TBranch *b_Electron_mvaNoIso;                                                                                     //!
   TBranch *b_Electron_pfRelIso03_all;                                                                               //!
   TBranch *b_Electron_pfRelIso03_chg;                                                                               //!
   TBranch *b_Electron_phi;                                                                                          //!
   TBranch *b_Electron_pt;                                                                                           //!
   TBranch *b_Electron_r9;                                                                                           //!
   TBranch *b_Electron_scEtOverPt;                                                                                   //!
   TBranch *b_Electron_sieie;                                                                                        //!
   TBranch *b_Electron_sip3d;                                                                                        //!
   TBranch *b_Electron_mvaTTH;                                                                                       //!
   TBranch *b_nFatJet;                                                                                               //!
   TBranch *b_FatJet_jetId;                                                                                          //!
   TBranch *b_FatJet_nConstituents;                                                                                  //!
   TBranch *b_FatJet_subJetIdx1;                                                                                     //!
   TBranch *b_FatJet_subJetIdx2;                                                                                     //!
   TBranch *b_FatJet_electronIdx3SJ;                                                                                 //!
   TBranch *b_FatJet_muonIdx3SJ;                                                                                     //!
   TBranch *b_FatJet_area;                                                                                           //!
   TBranch *b_FatJet_btagDDBvLV2;                                                                                    //!
   TBranch *b_FatJet_btagDDCvBV2;                                                                                    //!
   TBranch *b_FatJet_btagDDCvLV2;                                                                                    //!
   TBranch *b_FatJet_btagDeepB;                                                                                      //!
   TBranch *b_FatJet_btagHbb;                                                                                        //!
   TBranch *b_FatJet_eta;                                                                                            //!
   TBranch *b_FatJet_mass;                                                                                           //!
   TBranch *b_FatJet_msoftdrop;                                                                                      //!
   TBranch *b_FatJet_n2b1;                                                                                           //!
   TBranch *b_FatJet_n3b1;                                                                                           //!
   TBranch *b_FatJet_particleNetWithMass_H4qvsQCD;                                                                   //!
   TBranch *b_FatJet_particleNetWithMass_HbbvsQCD;                                                                   //!
   TBranch *b_FatJet_particleNetWithMass_HccvsQCD;                                                                   //!
   TBranch *b_FatJet_particleNetWithMass_QCD;                                                                        //!
   TBranch *b_FatJet_particleNetWithMass_TvsQCD;                                                                     //!
   TBranch *b_FatJet_particleNetWithMass_WvsQCD;                                                                     //!
   TBranch *b_FatJet_particleNetWithMass_ZvsQCD;                                                                     //!
   TBranch *b_FatJet_particleNet_QCD;                                                                                //!
   TBranch *b_FatJet_particleNet_QCD0HF;                                                                             //!
   TBranch *b_FatJet_particleNet_QCD1HF;                                                                             //!
   TBranch *b_FatJet_particleNet_QCD2HF;                                                                             //!
   TBranch *b_FatJet_particleNet_XbbVsQCD;                                                                           //!
   TBranch *b_FatJet_particleNet_XccVsQCD;                                                                           //!
   TBranch *b_FatJet_particleNet_XggVsQCD;                                                                           //!
   TBranch *b_FatJet_particleNet_XqqVsQCD;                                                                           //!
   TBranch *b_FatJet_particleNet_XteVsQCD;                                                                           //!
   TBranch *b_FatJet_particleNet_XtmVsQCD;                                                                           //!
   TBranch *b_FatJet_particleNet_XttVsQCD;                                                                           //!
   TBranch *b_FatJet_particleNet_massCorr;                                                                           //!
   TBranch *b_FatJet_phi;                                                                                            //!
   TBranch *b_FatJet_pt;                                                                                             //!
   TBranch *b_FatJet_rawFactor;                                                                                      //!
   TBranch *b_FatJet_tau1;                                                                                           //!
   TBranch *b_FatJet_tau2;                                                                                           //!
   TBranch *b_FatJet_tau3;                                                                                           //!
   TBranch *b_FatJet_tau4;                                                                                           //!
   TBranch *b_FatJet_lsf3;                                                                                           //!
   TBranch *b_nFsrPhoton;                                                                                            //!
   TBranch *b_FsrPhoton_electronIdx;                                                                                 //!
   TBranch *b_FsrPhoton_muonIdx;                                                                                     //!
   TBranch *b_FsrPhoton_dROverEt2;                                                                                   //!
   TBranch *b_FsrPhoton_eta;                                                                                         //!
   TBranch *b_FsrPhoton_phi;                                                                                         //!
   TBranch *b_FsrPhoton_pt;                                                                                          //!
   TBranch *b_FsrPhoton_relIso03;                                                                                    //!
   TBranch *b_nIsoTrack;                                                                                             //!
   TBranch *b_IsoTrack_isHighPurityTrack;                                                                            //!
   TBranch *b_IsoTrack_isPFcand;                                                                                     //!
   TBranch *b_IsoTrack_isFromLostTrack;                                                                              //!
   TBranch *b_IsoTrack_charge;                                                                                       //!
   TBranch *b_IsoTrack_fromPV;                                                                                       //!
   TBranch *b_IsoTrack_pdgId;                                                                                        //!
   TBranch *b_IsoTrack_dxy;                                                                                          //!
   TBranch *b_IsoTrack_dz;                                                                                           //!
   TBranch *b_IsoTrack_eta;                                                                                          //!
   TBranch *b_IsoTrack_pfRelIso03_all;                                                                               //!
   TBranch *b_IsoTrack_pfRelIso03_chg;                                                                               //!
   TBranch *b_IsoTrack_phi;                                                                                          //!
   TBranch *b_IsoTrack_pt;                                                                                           //!
   TBranch *b_IsoTrack_miniPFRelIso_all;                                                                             //!
   TBranch *b_IsoTrack_miniPFRelIso_chg;                                                                             //!
   TBranch *b_nJet;                                                                                                  //!
   TBranch *b_Jet_jetId;                                                                                             //!
   TBranch *b_Jet_nConstituents;                                                                                     //!
   TBranch *b_Jet_nElectrons;                                                                                        //!
   TBranch *b_Jet_nMuons;                                                                                            //!
   TBranch *b_Jet_nSVs;                                                                                              //!
   TBranch *b_Jet_electronIdx1;                                                                                      //!
   TBranch *b_Jet_electronIdx2;                                                                                      //!
   TBranch *b_Jet_muonIdx1;                                                                                          //!
   TBranch *b_Jet_muonIdx2;                                                                                          //!
   TBranch *b_Jet_svIdx1;                                                                                            //!
   TBranch *b_Jet_svIdx2;                                                                                            //!
   TBranch *b_Jet_hfadjacentEtaStripsSize;                                                                           //!
   TBranch *b_Jet_hfcentralEtaStripSize;                                                                             //!
   TBranch *b_Jet_PNetRegPtRawCorr;                                                                                  //!
   TBranch *b_Jet_PNetRegPtRawCorrNeutrino;                                                                          //!
   TBranch *b_Jet_PNetRegPtRawRes;                                                                                   //!
   TBranch *b_Jet_area;                                                                                              //!
   TBranch *b_Jet_btagDeepFlavB;                                                                                     //!
   TBranch *b_Jet_btagDeepFlavCvB;                                                                                   //!
   TBranch *b_Jet_btagDeepFlavCvL;                                                                                   //!
   TBranch *b_Jet_btagDeepFlavQG;                                                                                    //!
   TBranch *b_Jet_btagPNetB;                                                                                         //!
   TBranch *b_Jet_btagPNetCvB;                                                                                       //!
   TBranch *b_Jet_btagPNetCvL;                                                                                       //!
   TBranch *b_Jet_btagPNetQvG;                                                                                       //!
   TBranch *b_Jet_btagPNetTauVJet;                                                                                   //!
   TBranch *b_Jet_btagRobustParTAK4B;                                                                                //!
   TBranch *b_Jet_btagRobustParTAK4CvB;                                                                              //!
   TBranch *b_Jet_btagRobustParTAK4CvL;                                                                              //!
   TBranch *b_Jet_btagRobustParTAK4QG;                                                                               //!
   TBranch *b_Jet_chEmEF;                                                                                            //!
   TBranch *b_Jet_chHEF;                                                                                             //!
   TBranch *b_Jet_eta;                                                                                               //!
   TBranch *b_Jet_hfsigmaEtaEta;                                                                                     //!
   TBranch *b_Jet_hfsigmaPhiPhi;                                                                                     //!
   TBranch *b_Jet_mass;                                                                                              //!
   TBranch *b_Jet_muEF;                                                                                              //!
   TBranch *b_Jet_muonSubtrFactor;                                                                                   //!
   TBranch *b_Jet_neEmEF;                                                                                            //!
   TBranch *b_Jet_neHEF;                                                                                             //!
   TBranch *b_Jet_phi;                                                                                               //!
   TBranch *b_Jet_pt;                                                                                                //!
   TBranch *b_Jet_rawFactor;                                                                                         //!
   TBranch *b_nLowPtElectron;                                                                                        //!
   TBranch *b_LowPtElectron_convVeto;                                                                                //!
   TBranch *b_LowPtElectron_convWP;                                                                                  //!
   TBranch *b_LowPtElectron_lostHits;                                                                                //!
   TBranch *b_LowPtElectron_electronIdx;                                                                             //!
   TBranch *b_LowPtElectron_photonIdx;                                                                               //!
   TBranch *b_LowPtElectron_charge;                                                                                  //!
   TBranch *b_LowPtElectron_pdgId;                                                                                   //!
   TBranch *b_LowPtElectron_ID;                                                                                      //!
   TBranch *b_LowPtElectron_convVtxRadius;                                                                           //!
   TBranch *b_LowPtElectron_deltaEtaSC;                                                                              //!
   TBranch *b_LowPtElectron_dxy;                                                                                     //!
   TBranch *b_LowPtElectron_dxyErr;                                                                                  //!
   TBranch *b_LowPtElectron_dz;                                                                                      //!
   TBranch *b_LowPtElectron_dzErr;                                                                                   //!
   TBranch *b_LowPtElectron_eInvMinusPInv;                                                                           //!
   TBranch *b_LowPtElectron_energyErr;                                                                               //!
   TBranch *b_LowPtElectron_eta;                                                                                     //!
   TBranch *b_LowPtElectron_hoe;                                                                                     //!
   TBranch *b_LowPtElectron_mass;                                                                                    //!
   TBranch *b_LowPtElectron_miniPFRelIso_all;                                                                        //!
   TBranch *b_LowPtElectron_miniPFRelIso_chg;                                                                        //!
   TBranch *b_LowPtElectron_phi;                                                                                     //!
   TBranch *b_LowPtElectron_pt;                                                                                      //!
   TBranch *b_LowPtElectron_ptbiased;                                                                                //!
   TBranch *b_LowPtElectron_r9;                                                                                      //!
   TBranch *b_LowPtElectron_scEtOverPt;                                                                              //!
   TBranch *b_LowPtElectron_sieie;                                                                                   //!
   TBranch *b_LowPtElectron_unbiased;                                                                                //!
   TBranch *b_MET_MetUnclustEnUpDeltaX;                                                                              //!
   TBranch *b_MET_MetUnclustEnUpDeltaY;                                                                              //!
   TBranch *b_MET_covXX;                                                                                             //!
   TBranch *b_MET_covXY;                                                                                             //!
   TBranch *b_MET_covYY;                                                                                             //!
   TBranch *b_MET_phi;                                                                                               //!
   TBranch *b_MET_pt;                                                                                                //!
   TBranch *b_MET_significance;                                                                                      //!
   TBranch *b_MET_sumEt;                                                                                             //!
   TBranch *b_MET_sumPtUnclustered;                                                                                  //!
   TBranch *b_nProton_multiRP;                                                                                       //!
   TBranch *b_Proton_multiRP_arm;                                                                                    //!
   TBranch *b_Proton_multiRP_t;                                                                                      //!
   TBranch *b_Proton_multiRP_thetaX;                                                                                 //!
   TBranch *b_Proton_multiRP_thetaY;                                                                                 //!
   TBranch *b_Proton_multiRP_time;                                                                                   //!
   TBranch *b_Proton_multiRP_timeUnc;                                                                                //!
   TBranch *b_Proton_multiRP_xi;                                                                                     //!
   TBranch *b_nMuon;                                                                                                 //!
   TBranch *b_Muon_highPtId;                                                                                         //!
   TBranch *b_Muon_highPurity;                                                                                       //!
   TBranch *b_Muon_inTimeMuon;                                                                                       //!
   TBranch *b_Muon_isGlobal;                                                                                         //!
   TBranch *b_Muon_isPFcand;                                                                                         //!
   TBranch *b_Muon_isStandalone;                                                                                     //!
   TBranch *b_Muon_isTracker;                                                                                        //!
   TBranch *b_Muon_jetNDauCharged;                                                                                   //!
   TBranch *b_Muon_looseId;                                                                                          //!
   TBranch *b_Muon_mediumId;                                                                                         //!
   TBranch *b_Muon_mediumPromptId;                                                                                   //!
   TBranch *b_Muon_miniIsoId;                                                                                        //!
   TBranch *b_Muon_multiIsoId;                                                                                       //!
   TBranch *b_Muon_mvaMuID_WP;                                                                                       //!
   TBranch *b_Muon_nStations;                                                                                        //!
   TBranch *b_Muon_nTrackerLayers;                                                                                   //!
   TBranch *b_Muon_pfIsoId;                                                                                          //!
   TBranch *b_Muon_puppiIsoId;                                                                                       //!
   TBranch *b_Muon_softId;                                                                                           //!
   TBranch *b_Muon_softMvaId;                                                                                        //!
   TBranch *b_Muon_tightCharge;                                                                                      //!
   TBranch *b_Muon_tightId;                                                                                          //!
   TBranch *b_Muon_tkIsoId;                                                                                          //!
   TBranch *b_Muon_triggerIdLoose;                                                                                   //!
   TBranch *b_Muon_jetIdx;                                                                                           //!
   TBranch *b_Muon_svIdx;                                                                                            //!
   TBranch *b_Muon_fsrPhotonIdx;                                                                                     //!
   TBranch *b_Muon_charge;                                                                                           //!
   TBranch *b_Muon_pdgId;                                                                                            //!
   TBranch *b_Muon_dxy;                                                                                              //!
   TBranch *b_Muon_dxyErr;                                                                                           //!
   TBranch *b_Muon_dxybs;                                                                                            //!
   TBranch *b_Muon_dz;                                                                                               //!
   TBranch *b_Muon_dzErr;                                                                                            //!
   TBranch *b_Muon_eta;                                                                                              //!
   TBranch *b_Muon_ip3d;                                                                                             //!
   TBranch *b_Muon_jetPtRelv2;                                                                                       //!
   TBranch *b_Muon_jetRelIso;                                                                                        //!
   TBranch *b_Muon_mass;                                                                                             //!
   TBranch *b_Muon_miniPFRelIso_all;                                                                                 //!
   TBranch *b_Muon_miniPFRelIso_chg;                                                                                 //!
   TBranch *b_Muon_mvaMuID;                                                                                          //!
   TBranch *b_Muon_pfRelIso03_all;                                                                                   //!
   TBranch *b_Muon_pfRelIso03_chg;                                                                                   //!
   TBranch *b_Muon_pfRelIso04_all;                                                                                   //!
   TBranch *b_Muon_phi;                                                                                              //!
   TBranch *b_Muon_pt;                                                                                               //!
   TBranch *b_Muon_ptErr;                                                                                            //!
   TBranch *b_Muon_segmentComp;                                                                                      //!
   TBranch *b_Muon_sip3d;                                                                                            //!
   TBranch *b_Muon_softMva;                                                                                          //!
   TBranch *b_Muon_tkRelIso;                                                                                         //!
   TBranch *b_Muon_tunepRelPt;                                                                                       //!
   TBranch *b_Muon_bsConstrainedChi2;                                                                                //!
   TBranch *b_Muon_bsConstrainedPt;                                                                                  //!
   TBranch *b_Muon_bsConstrainedPtErr;                                                                               //!
   TBranch *b_Muon_mvaLowPt;                                                                                         //!
   TBranch *b_Muon_mvaTTH;   
   TBranch        *b_PFMET_covXX;   //!
   TBranch        *b_PFMET_covXY;   //!
   TBranch        *b_PFMET_covYY;   //!
   TBranch        *b_PFMET_phi;   //!
   TBranch        *b_PFMET_phiUnclusteredDown;   //!
   TBranch        *b_PFMET_phiUnclusteredUp;   //!
   TBranch        *b_PFMET_pt;   //!
   TBranch        *b_PFMET_ptUnclusteredDown;   //!
   TBranch        *b_PFMET_ptUnclusteredUp;   //!
   TBranch        *b_PFMET_significance;   //!
   TBranch        *b_PFMET_sumEt;   //!
   TBranch        *b_PFMET_sumPtUnclustered;   //!                                                                                        //!
   TBranch *b_nPhoton;                                                                                               //!
   TBranch *b_Photon_seediEtaOriX;                                                                                   //!
   TBranch *b_Photon_cutBased;                                                                                       //!
   TBranch *b_Photon_electronVeto;                                                                                   //!
   TBranch *b_Photon_hasConversionTracks;                                                                            //!
   TBranch *b_Photon_isScEtaEB;                                                                                      //!
   TBranch *b_Photon_isScEtaEE;                                                                                      //!
   TBranch *b_Photon_mvaID_WP80;                                                                                     //!
   TBranch *b_Photon_mvaID_WP90;                                                                                     //!
   TBranch *b_Photon_pixelSeed;                                                                                      //!
   TBranch *b_Photon_seedGain;                                                                                       //!
   TBranch *b_Photon_electronIdx;                                                                                    //!
   TBranch *b_Photon_jetIdx;                                                                                         //!
   TBranch *b_Photon_seediPhiOriY;                                                                                   //!
   TBranch *b_Photon_vidNestedWPBitmap;                                                                              //!
   TBranch *b_Photon_ecalPFClusterIso;                                                                               //!
   TBranch *b_Photon_energyErr;                                                                                      //!
   TBranch *b_Photon_energyRaw;                                                                                      //!
   TBranch *b_Photon_esEffSigmaRR;                                                                                   //!
   TBranch *b_Photon_esEnergyOverRawE;                                                                               //!
   TBranch *b_Photon_eta;                                                                                            //!
   TBranch *b_Photon_etaWidth;                                                                                       //!
   TBranch *b_Photon_haloTaggerMVAVal;                                                                               //!
   TBranch *b_Photon_hcalPFClusterIso;                                                                               //!
   TBranch *b_Photon_hoe;                                                                                            //!
   TBranch *b_Photon_hoe_PUcorr;                                                                                     //!
   TBranch *b_Photon_mvaID;                                                                                          //!
   TBranch *b_Photon_pfChargedIso;                                                                                   //!
   TBranch *b_Photon_pfChargedIsoPFPV;                                                                               //!
   TBranch *b_Photon_pfChargedIsoWorstVtx;                                                                           //!
   TBranch *b_Photon_pfPhoIso03;                                                                                     //!
   TBranch *b_Photon_pfRelIso03_all_quadratic;                                                                       //!
   TBranch *b_Photon_pfRelIso03_chg_quadratic;                                                                       //!
   TBranch *b_Photon_phi;                                                                                            //!
   TBranch *b_Photon_phiWidth;                                                                                       //!
   TBranch *b_Photon_pt;                                                                                             //!
   TBranch *b_Photon_r9;                                                                                             //!
   TBranch *b_Photon_s4;                                                                                             //!
   TBranch *b_Photon_sieie;                                                                                          //!
   TBranch *b_Photon_sieip;                                                                                          //!
   TBranch *b_Photon_sipip;                                                                                          //!
   TBranch *b_Photon_trkSumPtHollowConeDR03;                                                                         //!
   TBranch *b_Photon_trkSumPtSolidConeDR04;                                                                          //!
   TBranch *b_Photon_x_calo;                                                                                         //!
   TBranch *b_Photon_y_calo;                                                                                         //!
   TBranch *b_Photon_z_calo;                                                                                         //!
   TBranch *b_nPPSLocalTrack;                                                                                        //!
   TBranch *b_PPSLocalTrack_multiRPProtonIdx;                                                                        //!
   TBranch *b_PPSLocalTrack_singleRPProtonIdx;                                                                       //!
   TBranch *b_PPSLocalTrack_decRPId;                                                                                 //!
   TBranch *b_PPSLocalTrack_rpType;                                                                                  //!
   TBranch *b_PPSLocalTrack_x;                                                                                       //!
   TBranch *b_PPSLocalTrack_y;                                                                                       //!
   TBranch *b_PPSLocalTrack_time;                                                                                    //!
   TBranch *b_PPSLocalTrack_timeUnc;                                                                                 //!
   TBranch        *b_Pileup_nPU;   //!
   TBranch        *b_Pileup_sumEOOT;   //!
   TBranch        *b_Pileup_sumLOOT;   //!
   TBranch        *b_Pileup_nTrueInt;   //!
   TBranch        *b_Pileup_pudensity;   //!
   TBranch        *b_Pileup_gpudensity;   //!
   TBranch *b_PuppiMET_phi;                                                                                          //!
   TBranch *b_PuppiMET_phiJERDown;                                                                                   //!
   TBranch *b_PuppiMET_phiJERUp;                                                                                     //!
   TBranch *b_PuppiMET_phiJESDown;                                                                                   //!
   TBranch *b_PuppiMET_phiJESUp;                                                                                     //!
   TBranch *b_PuppiMET_phiUnclusteredDown;                                                                           //!
   TBranch *b_PuppiMET_phiUnclusteredUp;                                                                             //!
   TBranch *b_PuppiMET_pt;                                                                                           //!
   TBranch *b_PuppiMET_ptJERDown;                                                                                    //!
   TBranch *b_PuppiMET_ptJERUp;                                                                                      //!
   TBranch *b_PuppiMET_ptJESDown;                                                                                    //!
   TBranch *b_PuppiMET_ptJESUp;                                                                                      //!
   TBranch *b_PuppiMET_ptUnclusteredDown;                                                                            //!
   TBranch *b_PuppiMET_ptUnclusteredUp;                                                                              //!
   TBranch *b_PuppiMET_sumEt;                                                                                        //!
   TBranch *b_RawMET_phi;                                                                                            //!
   TBranch *b_RawMET_pt;                                                                                             //!
   TBranch *b_RawMET_sumEt;                                                                                          //!
   TBranch *b_RawPuppiMET_phi;                                                                                       //!
   TBranch *b_RawPuppiMET_pt;                                                                                        //!
   TBranch *b_RawPuppiMET_sumEt;                                                                                     //!
   TBranch *b_Rho_fixedGridRhoAll;                                                                                   //!
   TBranch *b_Rho_fixedGridRhoFastjetAll;                                                                            //!
   TBranch *b_Rho_fixedGridRhoFastjetCentral;                                                                        //!
   TBranch *b_Rho_fixedGridRhoFastjetCentralCalo;                                                                    //!
   TBranch *b_Rho_fixedGridRhoFastjetCentralChargedPileUp;                                                           //!
   TBranch *b_Rho_fixedGridRhoFastjetCentralNeutral;                                                                 //!
   TBranch *b_nSoftActivityJet;                                                                                      //!
   TBranch *b_SoftActivityJet_eta;                                                                                   //!
   TBranch *b_SoftActivityJet_phi;                                                                                   //!
   TBranch *b_SoftActivityJet_pt;                                                                                    //!
   TBranch *b_SoftActivityJetNjets10;                                                                                //!
   TBranch *b_SoftActivityJetNjets2;                                                                                 //!
   TBranch *b_SoftActivityJetNjets5;                                                                                 //!
   TBranch *b_SoftActivityJetHT;                                                                                     //!
   TBranch *b_SoftActivityJetHT10;                                                                                   //!
   TBranch *b_SoftActivityJetHT2;                                                                                    //!
   TBranch *b_SoftActivityJetHT5;                                                                                    //!
   TBranch *b_nProton_singleRP;                                                                                      //!
   TBranch *b_Proton_singleRP_decRPId;                                                                               //!
   TBranch *b_Proton_singleRP_thetaY;                                                                                //!
   TBranch *b_Proton_singleRP_xi;                                                                                    //!
   TBranch *b_nSubJet;                                                                                               //!
   TBranch *b_SubJet_btagDeepB;                                                                                      //!
   TBranch *b_SubJet_eta;                                                                                            //!
   TBranch *b_SubJet_mass;                                                                                           //!
   TBranch *b_SubJet_n2b1;                                                                                           //!
   TBranch *b_SubJet_n3b1;                                                                                           //!
   TBranch *b_SubJet_phi;                                                                                            //!
   TBranch *b_SubJet_pt;                                                                                             //!
   TBranch *b_SubJet_rawFactor;                                                                                      //!
   TBranch *b_SubJet_tau1;                                                                                           //!
   TBranch *b_SubJet_tau2;                                                                                           //!
   TBranch *b_SubJet_tau3;                                                                                           //!
   TBranch *b_SubJet_tau4;                                                                                           //!
   TBranch *b_nTau;                                                                                                  //!
   TBranch *b_Tau_decayMode;                                                                                         //!
   TBranch *b_Tau_idAntiEleDeadECal;                                                                                 //!
   TBranch *b_Tau_idAntiMu;                                                                                          //!
   TBranch *b_Tau_idDecayModeNewDMs;                                                                                 //!
   TBranch *b_Tau_idDecayModeOldDMs;                                                                                 //!
   TBranch *b_Tau_idDeepTau2017v2p1VSe;                                                                              //!
   TBranch *b_Tau_idDeepTau2017v2p1VSjet;                                                                            //!
   TBranch *b_Tau_idDeepTau2017v2p1VSmu;                                                                             //!
   TBranch *b_Tau_idDeepTau2018v2p5VSe;                                                                              //!
   TBranch *b_Tau_idDeepTau2018v2p5VSjet;                                                                            //!
   TBranch *b_Tau_idDeepTau2018v2p5VSmu;                                                                             //!
   TBranch *b_Tau_nSVs;                                                                                              //!
   TBranch *b_Tau_charge;                                                                                            //!
   TBranch *b_Tau_decayModePNet;                                                                                     //!
   TBranch *b_Tau_eleIdx;                                                                                            //!
   TBranch *b_Tau_jetIdx;                                                                                            //!
   TBranch *b_Tau_muIdx;                                                                                             //!
   TBranch *b_Tau_svIdx1;                                                                                            //!
   TBranch *b_Tau_svIdx2;                                                                                            //!
   TBranch *b_Tau_chargedIso;                                                                                        //!
   TBranch *b_Tau_dxy;                                                                                               //!
   TBranch *b_Tau_dz;                                                                                                //!
   TBranch *b_Tau_eta;                                                                                               //!
   TBranch *b_Tau_leadTkDeltaEta;                                                                                    //!
   TBranch *b_Tau_leadTkDeltaPhi;                                                                                    //!
   TBranch *b_Tau_leadTkPtOverTauPt;                                                                                 //!
   TBranch *b_Tau_mass;                                                                                              //!
   TBranch *b_Tau_neutralIso;                                                                                        //!
   TBranch *b_Tau_phi;                                                                                               //!
   TBranch *b_Tau_photonsOutsideSignalCone;                                                                          //!
   TBranch *b_Tau_probDM0PNet;                                                                                       //!
   TBranch *b_Tau_probDM10PNet;                                                                                      //!
   TBranch *b_Tau_probDM11PNet;                                                                                      //!
   TBranch *b_Tau_probDM1PNet;                                                                                       //!
   TBranch *b_Tau_probDM2PNet;                                                                                       //!
   TBranch *b_Tau_pt;                                                                                                //!
   TBranch *b_Tau_ptCorrPNet;                                                                                        //!
   TBranch *b_Tau_puCorr;                                                                                            //!
   TBranch *b_Tau_qConfPNet;                                                                                         //!
   TBranch *b_Tau_rawDeepTau2017v2p1VSe;                                                                             //!
   TBranch *b_Tau_rawDeepTau2017v2p1VSjet;                                                                           //!
   TBranch *b_Tau_rawDeepTau2017v2p1VSmu;                                                                            //!
   TBranch *b_Tau_rawDeepTau2018v2p5VSe;                                                                             //!
   TBranch *b_Tau_rawDeepTau2018v2p5VSjet;                                                                           //!
   TBranch *b_Tau_rawDeepTau2018v2p5VSmu;                                                                            //!
   TBranch *b_Tau_rawIso;                                                                                            //!
   TBranch *b_Tau_rawIsodR03;                                                                                        //!
   TBranch *b_Tau_rawPNetVSe;                                                                                        //!
   TBranch *b_Tau_rawPNetVSjet;                                                                                      //!
   TBranch *b_Tau_rawPNetVSmu;                                                                                       //!
   TBranch *b_TkMET_phi;                                                                                             //!
   TBranch *b_TkMET_pt;                                                                                              //!
   TBranch *b_TkMET_sumEt;                                                                                           //!
   TBranch *b_nTrigObj;                                                                                              //!
   TBranch *b_TrigObj_l1charge;                                                                                      //!
   TBranch *b_TrigObj_id;                                                                                            //!
   TBranch *b_TrigObj_l1iso;                                                                                         //!
   TBranch *b_TrigObj_filterBits;                                                                                    //!
   TBranch *b_TrigObj_pt;                                                                                            //!
   TBranch *b_TrigObj_eta;                                                                                           //!
   TBranch *b_TrigObj_phi;                                                                                           //!
   TBranch *b_TrigObj_l1pt;                                                                                          //!
   TBranch *b_TrigObj_l1pt_2;                                                                                        //!
   TBranch *b_TrigObj_l2pt;                                                                                          //!
   TBranch *b_nOtherPV;                                                                                              //!
   TBranch *b_OtherPV_z;                                                                                             //!
   TBranch *b_OtherPV_score;                                                                                         //!
   TBranch *b_PV_npvs;                                                                                               //!
   TBranch *b_PV_npvsGood;                                                                                           //!
   TBranch *b_PV_ndof;                                                                                               //!
   TBranch *b_PV_x;                                                                                                  //!
   TBranch *b_PV_y;                                                                                                  //!
   TBranch *b_PV_z;                                                                                                  //!
   TBranch *b_PV_chi2;                                                                                               //!
   TBranch *b_PV_score;                                                                                              //!
   TBranch *b_nSV;                                                                                                   //!
   TBranch *b_SV_charge;                                                                                             //!
   TBranch *b_SV_dlen;                                                                                               //!
   TBranch *b_SV_dlenSig;                                                                                            //!
   TBranch *b_SV_dxy;                                                                                                //!
   TBranch *b_SV_dxySig;                                                                                             //!
   TBranch *b_SV_pAngle;                                                                                             //!
   TBranch *b_SV_ntracks;                                                                                            //!
   TBranch *b_SV_chi2;                                                                                               //!
   TBranch *b_SV_eta;                                                                                                //!
   TBranch *b_SV_mass;                                                                                               //!
   TBranch *b_SV_ndof;                                                                                               //!
   TBranch *b_SV_phi;                                                                                                //!
   TBranch *b_SV_pt;                                                                                                 //!
   TBranch *b_SV_x;                                                                                                  //!
   TBranch *b_SV_y;                                                                                                  //!
   TBranch *b_SV_z;                                                                                                  //!
   TBranch *b_Flag_HBHENoiseFilter;                                                                                  //!
   TBranch *b_Flag_HBHENoiseIsoFilter;                                                                               //!
   TBranch *b_Flag_CSCTightHaloFilter;                                                                               //!
   TBranch *b_Flag_CSCTightHaloTrkMuUnvetoFilter;                                                                    //!
   TBranch *b_Flag_CSCTightHalo2015Filter;                                                                           //!
   TBranch *b_Flag_globalTightHalo2016Filter;                                                                        //!
   TBranch *b_Flag_globalSuperTightHalo2016Filter;                                                                   //!
   TBranch *b_Flag_HcalStripHaloFilter;                                                                              //!
   TBranch *b_Flag_hcalLaserEventFilter;                                                                             //!
   TBranch *b_Flag_EcalDeadCellTriggerPrimitiveFilter;                                                               //!
   TBranch *b_Flag_EcalDeadCellBoundaryEnergyFilter;                                                                 //!
   TBranch *b_Flag_ecalBadCalibFilter;                                                                               //!
   TBranch *b_Flag_goodVertices;                                                                                     //!
   TBranch *b_Flag_eeBadScFilter;                                                                                    //!
   TBranch *b_Flag_ecalLaserCorrFilter;                                                                              //!
   TBranch *b_Flag_trkPOGFilters;                                                                                    //!
   TBranch *b_Flag_chargedHadronTrackResolutionFilter;                                                               //!
   TBranch *b_Flag_muonBadTrackFilter;                                                                               //!
   TBranch *b_Flag_BadChargedCandidateFilter;                                                                        //!
   TBranch *b_Flag_BadPFMuonFilter;                                                                                  //!
   TBranch *b_Flag_BadPFMuonDzFilter;                                                                                //!
   TBranch *b_Flag_hfNoisyHitsFilter;                                                                                //!
   TBranch *b_Flag_BadChargedCandidateSummer16Filter;                                                                //!
   TBranch *b_Flag_BadPFMuonSummer16Filter;                                                                          //!
   TBranch *b_Flag_trkPOG_manystripclus53X;                                                                          //!
   TBranch *b_Flag_trkPOG_toomanystripclus53X;                                                                       //!
   TBranch *b_Flag_trkPOG_logErrorTooManyClusters;                                                                   //!
   TBranch *b_Flag_METFilters;                                                                                       //!
   TBranch *b_L1_AlwaysTrue;                                                                                         //!
   TBranch *b_L1_BPTX_AND_Ref1_VME;                                                                                  //!
   TBranch *b_L1_BPTX_AND_Ref3_VME;                                                                                  //!
   TBranch *b_L1_BPTX_AND_Ref4_VME;                                                                                  //!
   TBranch *b_L1_BPTX_BeamGas_B1_VME;                                                                                //!
   TBranch *b_L1_BPTX_BeamGas_B2_VME;                                                                                //!
   TBranch *b_L1_BPTX_BeamGas_Ref1_VME;                                                                              //!
   TBranch *b_L1_BPTX_BeamGas_Ref2_VME;                                                                              //!
   TBranch *b_L1_BPTX_NotOR_VME;                                                                                     //!
   TBranch *b_L1_BPTX_OR_Ref3_VME;                                                                                   //!
   TBranch *b_L1_BPTX_OR_Ref4_VME;                                                                                   //!
   TBranch *b_L1_BPTX_RefAND_VME;                                                                                    //!
   TBranch *b_L1_BptxMinus;                                                                                          //!
   TBranch *b_L1_BptxOR;                                                                                             //!
   TBranch *b_L1_BptxPlus;                                                                                           //!
   TBranch *b_L1_BptxXOR;                                                                                            //!
   TBranch *b_L1_CDC_SingleMu_3_er1p2_TOP120_DPHI2p618_3p142;                                                        //!
   TBranch *b_L1_DoubleEG10_er1p2_dR_Max0p6;                                                                         //!
   TBranch *b_L1_DoubleEG10p5_er1p2_dR_Max0p6;                                                                       //!
   TBranch *b_L1_DoubleEG11_er1p2_dR_Max0p6;                                                                         //!
   TBranch *b_L1_DoubleEG4_er1p2_dR_Max0p9;                                                                          //!
   TBranch *b_L1_DoubleEG4p5_er1p2_dR_Max0p9;                                                                        //!
   TBranch *b_L1_DoubleEG5_er1p2_dR_Max0p9;                                                                          //!
   TBranch *b_L1_DoubleEG5p5_er1p2_dR_Max0p8;                                                                        //!
   TBranch *b_L1_DoubleEG6_er1p2_dR_Max0p8;                                                                          //!
   TBranch *b_L1_DoubleEG6p5_er1p2_dR_Max0p8;                                                                        //!
   TBranch *b_L1_DoubleEG7_er1p2_dR_Max0p8;                                                                          //!
   TBranch *b_L1_DoubleEG7p5_er1p2_dR_Max0p7;                                                                        //!
   TBranch *b_L1_DoubleEG8_er1p2_dR_Max0p7;                                                                          //!
   TBranch *b_L1_DoubleEG8er2p5_HTT260er;                                                                            //!
   TBranch *b_L1_DoubleEG8er2p5_HTT280er;                                                                            //!
   TBranch *b_L1_DoubleEG8er2p5_HTT300er;                                                                            //!
   TBranch *b_L1_DoubleEG8er2p5_HTT320er;                                                                            //!
   TBranch *b_L1_DoubleEG8er2p5_HTT340er;                                                                            //!
   TBranch *b_L1_DoubleEG8p5_er1p2_dR_Max0p7;                                                                        //!
   TBranch *b_L1_DoubleEG9_er1p2_dR_Max0p7;                                                                          //!
   TBranch *b_L1_DoubleEG9p5_er1p2_dR_Max0p6;                                                                        //!
   TBranch *b_L1_DoubleEG_15_10_er2p5;                                                                               //!
   TBranch *b_L1_DoubleEG_20_10_er2p5;                                                                               //!
   TBranch *b_L1_DoubleEG_22_10_er2p5;                                                                               //!
   TBranch *b_L1_DoubleEG_25_12_er2p5;                                                                               //!
   TBranch *b_L1_DoubleEG_25_14_er2p5;                                                                               //!
   TBranch *b_L1_DoubleEG_27_14_er2p5;                                                                               //!
   TBranch *b_L1_DoubleEG_LooseIso16_LooseIso12_er1p5;                                                               //!
   TBranch *b_L1_DoubleEG_LooseIso18_LooseIso12_er1p5;                                                               //!
   TBranch *b_L1_DoubleEG_LooseIso20_10_er2p5;                                                                       //!
   TBranch *b_L1_DoubleEG_LooseIso20_LooseIso12_er1p5;                                                               //!
   TBranch *b_L1_DoubleEG_LooseIso22_10_er2p5;                                                                       //!
   TBranch *b_L1_DoubleEG_LooseIso22_12_er2p5;                                                                       //!
   TBranch *b_L1_DoubleEG_LooseIso22_LooseIso12_er1p5;                                                               //!
   TBranch *b_L1_DoubleEG_LooseIso25_12_er2p5;                                                                       //!
   TBranch *b_L1_DoubleEG_LooseIso25_LooseIso12_er1p5;                                                               //!
   TBranch *b_L1_DoubleIsoTau26er2p1_Jet55_RmOvlp_dR0p5;                                                             //!
   TBranch *b_L1_DoubleIsoTau26er2p1_Jet70_RmOvlp_dR0p5;                                                             //!
   TBranch *b_L1_DoubleIsoTau28er2p1;                                                                                //!
   TBranch *b_L1_DoubleIsoTau28er2p1_Mass_Max80;                                                                     //!
   TBranch *b_L1_DoubleIsoTau28er2p1_Mass_Max90;                                                                     //!
   TBranch *b_L1_DoubleIsoTau30er2p1;                                                                                //!
   TBranch *b_L1_DoubleIsoTau30er2p1_Mass_Max80;                                                                     //!
   TBranch *b_L1_DoubleIsoTau30er2p1_Mass_Max90;                                                                     //!
   TBranch *b_L1_DoubleIsoTau32er2p1;                                                                                //!
   TBranch *b_L1_DoubleIsoTau34er2p1;                                                                                //!
   TBranch *b_L1_DoubleIsoTau35er2p1;                                                                                //!
   TBranch *b_L1_DoubleIsoTau36er2p1;                                                                                //!
   TBranch *b_L1_DoubleJet100er2p3_dEta_Max1p6;                                                                      //!
   TBranch *b_L1_DoubleJet100er2p5;                                                                                  //!
   TBranch *b_L1_DoubleJet112er2p3_dEta_Max1p6;                                                                      //!
   TBranch *b_L1_DoubleJet120er2p5;                                                                                  //!
   TBranch *b_L1_DoubleJet150er2p5;                                                                                  //!
   TBranch *b_L1_DoubleJet30er2p5_Mass_Min150_dEta_Max1p5;                                                           //!
   TBranch *b_L1_DoubleJet30er2p5_Mass_Min200_dEta_Max1p5;                                                           //!
   TBranch *b_L1_DoubleJet30er2p5_Mass_Min250_dEta_Max1p5;                                                           //!
   TBranch *b_L1_DoubleJet30er2p5_Mass_Min300_dEta_Max1p5;                                                           //!
   TBranch *b_L1_DoubleJet30er2p5_Mass_Min330_dEta_Max1p5;                                                           //!
   TBranch *b_L1_DoubleJet30er2p5_Mass_Min360_dEta_Max1p5;                                                           //!
   TBranch *b_L1_DoubleJet35_Mass_Min450_IsoTau45_RmOvlp;                                                            //!
   TBranch *b_L1_DoubleJet35_Mass_Min450_IsoTau45er2p1_RmOvlp_dR0p5;                                                 //!
   TBranch *b_L1_DoubleJet40er2p5;                                                                                   //!
   TBranch *b_L1_DoubleJet_100_30_DoubleJet30_Mass_Min620;                                                           //!
   TBranch *b_L1_DoubleJet_110_35_DoubleJet35_Mass_Min620;                                                           //!
   TBranch *b_L1_DoubleJet_115_40_DoubleJet40_Mass_Min620;                                                           //!
   TBranch *b_L1_DoubleJet_115_40_DoubleJet40_Mass_Min620_Jet60TT28;                                                 //!
   TBranch *b_L1_DoubleJet_120_45_DoubleJet45_Mass_Min620;                                                           //!
   TBranch *b_L1_DoubleJet_120_45_DoubleJet45_Mass_Min620_Jet60TT28;                                                 //!
   TBranch *b_L1_DoubleJet_80_30_Mass_Min420_DoubleMu0_SQ;                                                           //!
   TBranch *b_L1_DoubleJet_80_30_Mass_Min420_IsoTau40_RmOvlp;                                                        //!
   TBranch *b_L1_DoubleJet_80_30_Mass_Min420_Mu8;                                                                    //!
   TBranch *b_L1_DoubleJet_90_30_DoubleJet30_Mass_Min620;                                                            //!
   TBranch *b_L1_DoubleLLPJet40;                                                                                     //!
   TBranch *b_L1_DoubleLooseIsoEG22er2p1;                                                                            //!
   TBranch *b_L1_DoubleLooseIsoEG24er2p1;                                                                            //!
   TBranch *b_L1_DoubleMu0;                                                                                          //!
   TBranch *b_L1_DoubleMu0_Mass_Min1;                                                                                //!
   TBranch *b_L1_DoubleMu0_OQ;                                                                                       //!
   TBranch *b_L1_DoubleMu0_SQ;                                                                                       //!
   TBranch *b_L1_DoubleMu0_SQ_OS;                                                                                    //!
   TBranch *b_L1_DoubleMu0_Upt15_Upt7;                                                                               //!
   TBranch *b_L1_DoubleMu0_Upt5_Upt5;                                                                                //!
   TBranch *b_L1_DoubleMu0_Upt6_IP_Min1_Upt4;                                                                        //!
   TBranch *b_L1_DoubleMu0_dR_Max1p6_Jet90er2p5_dR_Max0p8;                                                           //!
   TBranch *b_L1_DoubleMu0er1p4_OQ_OS_dEta_Max1p6;                                                                   //!
   TBranch *b_L1_DoubleMu0er1p4_SQ_OS_dR_Max1p4;                                                                     //!
   TBranch *b_L1_DoubleMu0er1p5_SQ;                                                                                  //!
   TBranch *b_L1_DoubleMu0er1p5_SQ_OS;                                                                               //!
   TBranch *b_L1_DoubleMu0er1p5_SQ_OS_dR_Max1p4;                                                                     //!
   TBranch *b_L1_DoubleMu0er1p5_SQ_dR_Max1p4;                                                                        //!
   TBranch *b_L1_DoubleMu0er2p0_SQ_OS_dEta_Max1p5;                                                                   //!
   TBranch *b_L1_DoubleMu0er2p0_SQ_OS_dEta_Max1p6;                                                                   //!
   TBranch *b_L1_DoubleMu0er2p0_SQ_OS_dR_Max1p4;                                                                     //!
   TBranch *b_L1_DoubleMu0er2p0_SQ_dEta_Max1p5;                                                                      //!
   TBranch *b_L1_DoubleMu0er2p0_SQ_dEta_Max1p6;                                                                      //!
   TBranch *b_L1_DoubleMu0er2p0_SQ_dR_Max1p4;                                                                        //!
   TBranch *b_L1_DoubleMu18er2p1_SQ;                                                                                 //!
   TBranch *b_L1_DoubleMu3_OS_er2p3_Mass_Max14_DoubleEG7p5_er2p1_Mass_Max20;                                         //!
   TBranch *b_L1_DoubleMu3_SQ_ETMHF50_HTT60er;                                                                       //!
   TBranch *b_L1_DoubleMu3_SQ_ETMHF50_Jet60er2p5;                                                                    //!
   TBranch *b_L1_DoubleMu3_SQ_ETMHF50_Jet60er2p5_OR_DoubleJet40er2p5;                                                //!
   TBranch *b_L1_DoubleMu3_SQ_ETMHF60_Jet60er2p5;                                                                    //!
   TBranch *b_L1_DoubleMu3_SQ_HTT220er;                                                                              //!
   TBranch *b_L1_DoubleMu3_SQ_HTT240er;                                                                              //!
   TBranch *b_L1_DoubleMu3_SQ_HTT260er;                                                                              //!
   TBranch *b_L1_DoubleMu3_dR_Max1p6_Jet90er2p5_dR_Max0p8;                                                           //!
   TBranch *b_L1_DoubleMu3er2p0_SQ_OS_dR_Max1p4;                                                                     //!
   TBranch *b_L1_DoubleMu4_SQ_EG9er2p5;                                                                              //!
   TBranch *b_L1_DoubleMu4_SQ_OS;                                                                                    //!
   TBranch *b_L1_DoubleMu4_SQ_OS_dR_Max1p2;                                                                          //!
   TBranch *b_L1_DoubleMu4p5_SQ_OS;                                                                                  //!
   TBranch *b_L1_DoubleMu4p5_SQ_OS_dR_Max1p2;                                                                        //!
   TBranch *b_L1_DoubleMu4p5er2p0_SQ_OS;                                                                             //!
   TBranch *b_L1_DoubleMu4p5er2p0_SQ_OS_Mass_7to18;                                                                  //!
   TBranch *b_L1_DoubleMu4p5er2p0_SQ_OS_Mass_Min7;                                                                   //!
   TBranch *b_L1_DoubleMu5_OS_er2p3_Mass_8to14_DoubleEG3er2p1_Mass_Max20;                                            //!
   TBranch *b_L1_DoubleMu5_SQ_EG9er2p5;                                                                              //!
   TBranch *b_L1_DoubleMu8_SQ;                                                                                       //!
   TBranch *b_L1_DoubleMu9_SQ;                                                                                       //!
   TBranch *b_L1_DoubleMu_12_5;                                                                                      //!
   TBranch *b_L1_DoubleMu_15_5_SQ;                                                                                   //!
   TBranch *b_L1_DoubleMu_15_7;                                                                                      //!
   TBranch *b_L1_DoubleMu_15_7_Mass_Min1;                                                                            //!
   TBranch *b_L1_DoubleMu_15_7_SQ;                                                                                   //!
   TBranch *b_L1_DoubleTau70er2p1;                                                                                   //!
   TBranch *b_L1_ETM120;                                                                                             //!
   TBranch *b_L1_ETM150;                                                                                             //!
   TBranch *b_L1_ETMHF100;                                                                                           //!
   TBranch *b_L1_ETMHF100_HTT60er;                                                                                   //!
   TBranch *b_L1_ETMHF110;                                                                                           //!
   TBranch *b_L1_ETMHF110_HTT60er;                                                                                   //!
   TBranch *b_L1_ETMHF110_HTT60er_NotSecondBunchInTrain;                                                             //!
   TBranch *b_L1_ETMHF120;                                                                                           //!
   TBranch *b_L1_ETMHF120_HTT60er;                                                                                   //!
   TBranch *b_L1_ETMHF120_NotSecondBunchInTrain;                                                                     //!
   TBranch *b_L1_ETMHF130;                                                                                           //!
   TBranch *b_L1_ETMHF130_HTT60er;                                                                                   //!
   TBranch *b_L1_ETMHF140;                                                                                           //!
   TBranch *b_L1_ETMHF150;                                                                                           //!
   TBranch *b_L1_ETMHF70;                                                                                            //!
   TBranch *b_L1_ETMHF70_HTT60er;                                                                                    //!
   TBranch *b_L1_ETMHF80;                                                                                            //!
   TBranch *b_L1_ETMHF80_HTT60er;                                                                                    //!
   TBranch *b_L1_ETMHF90;                                                                                            //!
   TBranch *b_L1_ETMHF90_HTT60er;                                                                                    //!
   TBranch *b_L1_ETMHF90_SingleJet60er2p5_dPhi_Min2p1;                                                               //!
   TBranch *b_L1_ETMHF90_SingleJet60er2p5_dPhi_Min2p6;                                                               //!
   TBranch *b_L1_ETMHF90_SingleJet80er2p5_dPhi_Min2p1;                                                               //!
   TBranch *b_L1_ETMHF90_SingleJet80er2p5_dPhi_Min2p6;                                                               //!
   TBranch *b_L1_ETT1200;                                                                                            //!
   TBranch *b_L1_ETT1600;                                                                                            //!
   TBranch *b_L1_ETT2000;                                                                                            //!
   TBranch *b_L1_FirstBunchAfterTrain;                                                                               //!
   TBranch *b_L1_FirstBunchBeforeTrain;                                                                              //!
   TBranch *b_L1_FirstBunchInTrain;                                                                                  //!
   TBranch *b_L1_FirstCollisionInOrbit;                                                                              //!
   TBranch *b_L1_FirstCollisionInTrain;                                                                              //!
   TBranch *b_L1_HCAL_LaserMon_Trig;                                                                                 //!
   TBranch *b_L1_HCAL_LaserMon_Veto;                                                                                 //!
   TBranch *b_L1_HTT120_SingleLLPJet40;                                                                              //!
   TBranch *b_L1_HTT120er;                                                                                           //!
   TBranch *b_L1_HTT160_SingleLLPJet50;                                                                              //!
   TBranch *b_L1_HTT160er;                                                                                           //!
   TBranch *b_L1_HTT200_SingleLLPJet60;                                                                              //!
   TBranch *b_L1_HTT200er;                                                                                           //!
   TBranch *b_L1_HTT240_SingleLLPJet70;                                                                              //!
   TBranch *b_L1_HTT255er;                                                                                           //!
   TBranch *b_L1_HTT280er;                                                                                           //!
   TBranch *b_L1_HTT280er_QuadJet_70_55_40_35_er2p5;                                                                 //!
   TBranch *b_L1_HTT320er;                                                                                           //!
   TBranch *b_L1_HTT320er_QuadJet_70_55_40_40_er2p5;                                                                 //!
   TBranch *b_L1_HTT320er_QuadJet_80_60_er2p1_45_40_er2p3;                                                           //!
   TBranch *b_L1_HTT320er_QuadJet_80_60_er2p1_50_45_er2p3;                                                           //!
   TBranch *b_L1_HTT360er;                                                                                           //!
   TBranch *b_L1_HTT400er;                                                                                           //!
   TBranch *b_L1_HTT450er;                                                                                           //!
   TBranch *b_L1_IsoEG32er2p5_Mt40;                                                                                  //!
   TBranch *b_L1_IsoTau52er2p1_QuadJet36er2p5;                                                                       //!
   TBranch *b_L1_IsolatedBunch;                                                                                      //!
   TBranch *b_L1_LastBunchInTrain;                                                                                   //!
   TBranch *b_L1_LastCollisionInTrain;                                                                               //!
   TBranch *b_L1_LooseIsoEG22er2p1_IsoTau26er2p1_dR_Min0p3;                                                          //!
   TBranch *b_L1_LooseIsoEG22er2p1_Tau70er2p1_dR_Min0p3;                                                             //!
   TBranch *b_L1_LooseIsoEG24er2p1_HTT100er;                                                                         //!
   TBranch *b_L1_LooseIsoEG24er2p1_IsoTau27er2p1_dR_Min0p3;                                                          //!
   TBranch *b_L1_LooseIsoEG26er2p1_HTT100er;                                                                         //!
   TBranch *b_L1_LooseIsoEG26er2p1_Jet34er2p5_dR_Min0p3;                                                             //!
   TBranch *b_L1_LooseIsoEG28er2p1_HTT100er;                                                                         //!
   TBranch *b_L1_LooseIsoEG28er2p1_Jet34er2p5_dR_Min0p3;                                                             //!
   TBranch *b_L1_LooseIsoEG30er2p1_HTT100er;                                                                         //!
   TBranch *b_L1_LooseIsoEG30er2p1_Jet34er2p5_dR_Min0p3;                                                             //!
   TBranch *b_L1_MinimumBiasHF0;                                                                                     //!
   TBranch *b_L1_MinimumBiasHF0_AND_BptxAND;                                                                         //!
   TBranch *b_L1_Mu10er2p3_Jet32er2p3_dR_Max0p4_DoubleJet32er2p3_dEta_Max1p6;                                        //!
   TBranch *b_L1_Mu12er2p3_Jet40er2p1_dR_Max0p4_DoubleJet40er2p1_dEta_Max1p6;                                        //!
   TBranch *b_L1_Mu12er2p3_Jet40er2p3_dR_Max0p4_DoubleJet40er2p3_dEta_Max1p6;                                        //!
   TBranch *b_L1_Mu18er2p1_Tau24er2p1;                                                                               //!
   TBranch *b_L1_Mu18er2p1_Tau26er2p1;                                                                               //!
   TBranch *b_L1_Mu18er2p1_Tau26er2p1_Jet55;                                                                         //!
   TBranch *b_L1_Mu18er2p1_Tau26er2p1_Jet70;                                                                         //!
   TBranch *b_L1_Mu20_EG10er2p5;                                                                                     //!
   TBranch *b_L1_Mu22er2p1_IsoTau28er2p1;                                                                            //!
   TBranch *b_L1_Mu22er2p1_IsoTau30er2p1;                                                                            //!
   TBranch *b_L1_Mu22er2p1_IsoTau32er2p1;                                                                            //!
   TBranch *b_L1_Mu22er2p1_IsoTau34er2p1;                                                                            //!
   TBranch *b_L1_Mu22er2p1_IsoTau36er2p1;                                                                            //!
   TBranch *b_L1_Mu22er2p1_IsoTau40er2p1;                                                                            //!
   TBranch *b_L1_Mu22er2p1_Tau70er2p1;                                                                               //!
   TBranch *b_L1_Mu3_Jet120er2p5_dR_Max0p4;                                                                          //!
   TBranch *b_L1_Mu3_Jet120er2p5_dR_Max0p8;                                                                          //!
   TBranch *b_L1_Mu3_Jet16er2p5_dR_Max0p4;                                                                           //!
   TBranch *b_L1_Mu3_Jet30er2p5;                                                                                     //!
   TBranch *b_L1_Mu3_Jet35er2p5_dR_Max0p4;                                                                           //!
   TBranch *b_L1_Mu3_Jet60er2p5_dR_Max0p4;                                                                           //!
   TBranch *b_L1_Mu3_Jet80er2p5_dR_Max0p4;                                                                           //!
   TBranch *b_L1_Mu3er1p5_Jet100er2p5_ETMHF40;                                                                       //!
   TBranch *b_L1_Mu3er1p5_Jet100er2p5_ETMHF50;                                                                       //!
   TBranch *b_L1_Mu5_EG23er2p5;                                                                                      //!
   TBranch *b_L1_Mu5_LooseIsoEG20er2p5;                                                                              //!
   TBranch *b_L1_Mu6_DoubleEG10er2p5;                                                                                //!
   TBranch *b_L1_Mu6_DoubleEG12er2p5;                                                                                //!
   TBranch *b_L1_Mu6_DoubleEG15er2p5;                                                                                //!
   TBranch *b_L1_Mu6_DoubleEG17er2p5;                                                                                //!
   TBranch *b_L1_Mu6_HTT240er;                                                                                       //!
   TBranch *b_L1_Mu6_HTT250er;                                                                                       //!
   TBranch *b_L1_Mu7_EG20er2p5;                                                                                      //!
   TBranch *b_L1_Mu7_EG23er2p5;                                                                                      //!
   TBranch *b_L1_Mu7_LooseIsoEG20er2p5;                                                                              //!
   TBranch *b_L1_Mu7_LooseIsoEG23er2p5;                                                                              //!
   TBranch *b_L1_NotBptxOR;                                                                                          //!
   TBranch *b_L1_QuadJet60er2p5;                                                                                     //!
   TBranch *b_L1_QuadJet_95_75_65_20_DoubleJet_75_65_er2p5_Jet20_FWD3p0;                                             //!
   TBranch *b_L1_QuadMu0;                                                                                            //!
   TBranch *b_L1_QuadMu0_OQ;                                                                                         //!
   TBranch *b_L1_QuadMu0_SQ;                                                                                         //!
   TBranch *b_L1_SecondBunchInTrain;                                                                                 //!
   TBranch *b_L1_SecondLastBunchInTrain;                                                                             //!
   TBranch *b_L1_SingleEG10er2p5;                                                                                    //!
   TBranch *b_L1_SingleEG15er2p5;                                                                                    //!
   TBranch *b_L1_SingleEG26er2p5;                                                                                    //!
   TBranch *b_L1_SingleEG28_FWD2p5;                                                                                  //!
   TBranch *b_L1_SingleEG28er1p5;                                                                                    //!
   TBranch *b_L1_SingleEG28er2p1;                                                                                    //!
   TBranch *b_L1_SingleEG28er2p5;                                                                                    //!
   TBranch *b_L1_SingleEG34er2p5;                                                                                    //!
   TBranch *b_L1_SingleEG36er2p5;                                                                                    //!
   TBranch *b_L1_SingleEG38er2p5;                                                                                    //!
   TBranch *b_L1_SingleEG40er2p5;                                                                                    //!
   TBranch *b_L1_SingleEG42er2p5;                                                                                    //!
   TBranch *b_L1_SingleEG45er2p5;                                                                                    //!
   TBranch *b_L1_SingleEG50;                                                                                         //!
   TBranch *b_L1_SingleEG60;                                                                                         //!
   TBranch *b_L1_SingleEG8er2p5;                                                                                     //!
   TBranch *b_L1_SingleIsoEG24er1p5;                                                                                 //!
   TBranch *b_L1_SingleIsoEG24er2p1;                                                                                 //!
   TBranch *b_L1_SingleIsoEG26er1p5;                                                                                 //!
   TBranch *b_L1_SingleIsoEG26er2p1;                                                                                 //!
   TBranch *b_L1_SingleIsoEG26er2p5;                                                                                 //!
   TBranch *b_L1_SingleIsoEG28_FWD2p5;                                                                               //!
   TBranch *b_L1_SingleIsoEG28er1p5;                                                                                 //!
   TBranch *b_L1_SingleIsoEG28er2p1;                                                                                 //!
   TBranch *b_L1_SingleIsoEG28er2p5;                                                                                 //!
   TBranch *b_L1_SingleIsoEG30er2p1;                                                                                 //!
   TBranch *b_L1_SingleIsoEG30er2p5;                                                                                 //!
   TBranch *b_L1_SingleIsoEG32er2p1;                                                                                 //!
   TBranch *b_L1_SingleIsoEG32er2p5;                                                                                 //!
   TBranch *b_L1_SingleIsoEG34er2p5;                                                                                 //!
   TBranch *b_L1_SingleIsoTau32er2p1;                                                                                //!
   TBranch *b_L1_SingleJet10erHE;                                                                                    //!
   TBranch *b_L1_SingleJet120;                                                                                       //!
   TBranch *b_L1_SingleJet120_FWD3p0;                                                                                //!
   TBranch *b_L1_SingleJet120er2p5;                                                                                  //!
   TBranch *b_L1_SingleJet12erHE;                                                                                    //!
   TBranch *b_L1_SingleJet140er2p5;                                                                                  //!
   TBranch *b_L1_SingleJet140er2p5_ETMHF70;                                                                          //!
   TBranch *b_L1_SingleJet140er2p5_ETMHF80;                                                                          //!
   TBranch *b_L1_SingleJet140er2p5_ETMHF90;                                                                          //!
   TBranch *b_L1_SingleJet160er2p5;                                                                                  //!
   TBranch *b_L1_SingleJet180;                                                                                       //!
   TBranch *b_L1_SingleJet180er2p5;                                                                                  //!
   TBranch *b_L1_SingleJet200;                                                                                       //!
   TBranch *b_L1_SingleJet20er2p5_NotBptxOR;                                                                         //!
   TBranch *b_L1_SingleJet20er2p5_NotBptxOR_3BX;                                                                     //!
   TBranch *b_L1_SingleJet35;                                                                                        //!
   TBranch *b_L1_SingleJet35_FWD3p0;                                                                                 //!
   TBranch *b_L1_SingleJet35er2p5;                                                                                   //!
   TBranch *b_L1_SingleJet43er2p5_NotBptxOR_3BX;                                                                     //!
   TBranch *b_L1_SingleJet46er2p5_NotBptxOR_3BX;                                                                     //!
   TBranch *b_L1_SingleJet60;                                                                                        //!
   TBranch *b_L1_SingleJet60_FWD3p0;                                                                                 //!
   TBranch *b_L1_SingleJet60er2p5;                                                                                   //!
   TBranch *b_L1_SingleJet8erHE;                                                                                     //!
   TBranch *b_L1_SingleJet90;                                                                                        //!
   TBranch *b_L1_SingleJet90_FWD3p0;                                                                                 //!
   TBranch *b_L1_SingleJet90er2p5;                                                                                   //!
   TBranch *b_L1_SingleLooseIsoEG26er1p5;                                                                            //!
   TBranch *b_L1_SingleLooseIsoEG26er2p5;                                                                            //!
   TBranch *b_L1_SingleLooseIsoEG28_FWD2p5;                                                                          //!
   TBranch *b_L1_SingleLooseIsoEG28er1p5;                                                                            //!
   TBranch *b_L1_SingleLooseIsoEG28er2p1;                                                                            //!
   TBranch *b_L1_SingleLooseIsoEG28er2p5;                                                                            //!
   TBranch *b_L1_SingleLooseIsoEG30er1p5;                                                                            //!
   TBranch *b_L1_SingleLooseIsoEG30er2p5;                                                                            //!
   TBranch *b_L1_SingleMu0_BMTF;                                                                                     //!
   TBranch *b_L1_SingleMu0_DQ;                                                                                       //!
   TBranch *b_L1_SingleMu0_EMTF;                                                                                     //!
   TBranch *b_L1_SingleMu0_OMTF;                                                                                     //!
   TBranch *b_L1_SingleMu10er1p5;                                                                                    //!
   TBranch *b_L1_SingleMu12_DQ_BMTF;                                                                                 //!
   TBranch *b_L1_SingleMu12_DQ_EMTF;                                                                                 //!
   TBranch *b_L1_SingleMu12_DQ_OMTF;                                                                                 //!
   TBranch *b_L1_SingleMu12er1p5;                                                                                    //!
   TBranch *b_L1_SingleMu14er1p5;                                                                                    //!
   TBranch *b_L1_SingleMu15_DQ;                                                                                      //!
   TBranch *b_L1_SingleMu16er1p5;                                                                                    //!
   TBranch *b_L1_SingleMu18;                                                                                         //!
   TBranch *b_L1_SingleMu18er1p5;                                                                                    //!
   TBranch *b_L1_SingleMu20;                                                                                         //!
   TBranch *b_L1_SingleMu22;                                                                                         //!
   TBranch *b_L1_SingleMu22_BMTF;                                                                                    //!
   TBranch *b_L1_SingleMu22_DQ;                                                                                      //!
   TBranch *b_L1_SingleMu22_EMTF;                                                                                    //!
   TBranch *b_L1_SingleMu22_OMTF;                                                                                    //!
   TBranch *b_L1_SingleMu22_OQ;                                                                                      //!
   TBranch *b_L1_SingleMu25;                                                                                         //!
   TBranch *b_L1_SingleMu3;                                                                                          //!
   TBranch *b_L1_SingleMu5;                                                                                          //!
   TBranch *b_L1_SingleMu6er1p5;                                                                                     //!
   TBranch *b_L1_SingleMu7;                                                                                          //!
   TBranch *b_L1_SingleMu7_DQ;                                                                                       //!
   TBranch *b_L1_SingleMu7er1p5;                                                                                     //!
   TBranch *b_L1_SingleMu8er1p5;                                                                                     //!
   TBranch *b_L1_SingleMu9er1p5;                                                                                     //!
   TBranch *b_L1_SingleMuCosmics;                                                                                    //!
   TBranch *b_L1_SingleMuCosmics_BMTF;                                                                               //!
   TBranch *b_L1_SingleMuCosmics_EMTF;                                                                               //!
   TBranch *b_L1_SingleMuCosmics_OMTF;                                                                               //!
   TBranch *b_L1_SingleMuOpen;                                                                                       //!
   TBranch *b_L1_SingleMuOpen_NotBptxOR;                                                                             //!
   TBranch *b_L1_SingleMuOpen_er1p1_NotBptxOR_3BX;                                                                   //!
   TBranch *b_L1_SingleMuOpen_er1p4_NotBptxOR_3BX;                                                                   //!
   TBranch *b_L1_SingleMuShower_Nominal;                                                                             //!
   TBranch *b_L1_SingleMuShower_Tight;                                                                               //!
   TBranch *b_L1_SingleTau120er2p1;                                                                                  //!
   TBranch *b_L1_SingleTau130er2p1;                                                                                  //!
   TBranch *b_L1_SingleTau70er2p1;                                                                                   //!
   TBranch *b_L1_TOTEM_1;                                                                                            //!
   TBranch *b_L1_TOTEM_2;                                                                                            //!
   TBranch *b_L1_TOTEM_3;                                                                                            //!
   TBranch *b_L1_TOTEM_4;                                                                                            //!
   TBranch *b_L1_TripleEG16er2p5;                                                                                    //!
   TBranch *b_L1_TripleEG_16_12_8_er2p5;                                                                             //!
   TBranch *b_L1_TripleEG_16_15_8_er2p5;                                                                             //!
   TBranch *b_L1_TripleEG_18_17_8_er2p5;                                                                             //!
   TBranch *b_L1_TripleEG_18_18_12_er2p5;                                                                            //!
   TBranch *b_L1_TripleJet_100_80_70_DoubleJet_80_70_er2p5;                                                          //!
   TBranch *b_L1_TripleJet_105_85_75_DoubleJet_85_75_er2p5;                                                          //!
   TBranch *b_L1_TripleJet_95_75_65_DoubleJet_75_65_er2p5;                                                           //!
   TBranch *b_L1_TripleMu0;                                                                                          //!
   TBranch *b_L1_TripleMu0_OQ;                                                                                       //!
   TBranch *b_L1_TripleMu0_SQ;                                                                                       //!
   TBranch *b_L1_TripleMu3;                                                                                          //!
   TBranch *b_L1_TripleMu3_SQ;                                                                                       //!
   TBranch *b_L1_TripleMu_2SQ_1p5SQ_0OQ;                                                                             //!
   TBranch *b_L1_TripleMu_2SQ_1p5SQ_0OQ_Mass_Max12;                                                                  //!
   TBranch *b_L1_TripleMu_3SQ_2p5SQ_0OQ_Mass_Max12;                                                                  //!
   TBranch *b_L1_TripleMu_5SQ_3SQ_0OQ;                                                                               //!
   TBranch *b_L1_TripleMu_5SQ_3SQ_0OQ_DoubleMu_5_3_SQ_OS_Mass_Max9;                                                  //!
   TBranch *b_L1_TripleMu_5SQ_3SQ_0_DoubleMu_5_3_SQ_OS_Mass_Max9;                                                    //!
   TBranch *b_L1_TripleMu_5_3_3;                                                                                     //!
   TBranch *b_L1_TripleMu_5_3_3_SQ;                                                                                  //!
   TBranch *b_L1_TripleMu_5_3p5_2p5;                                                                                 //!
   TBranch *b_L1_TripleMu_5_3p5_2p5_DoubleMu_5_2p5_OS_Mass_5to17;                                                    //!
   TBranch *b_L1_TripleMu_5_3p5_2p5_OQ_DoubleMu_5_2p5_OQ_OS_Mass_5to17;                                              //!
   TBranch *b_L1_TripleMu_5_4_2p5_DoubleMu_5_2p5_OS_Mass_5to17;                                                      //!
   TBranch *b_L1_TripleMu_5_5_3;                                                                                     //!
   TBranch *b_L1_UnpairedBunchBptxMinus;                                                                             //!
   TBranch *b_L1_UnpairedBunchBptxPlus;                                                                              //!
   TBranch *b_L1_ZeroBias;                                                                                           //!
   TBranch *b_L1_ZeroBias_copy;                                                                                      //!
   TBranch *b_L1_UnprefireableEvent;                                                                                 //!
   TBranch *b_L1Reco_step;                                                                                           //!
   TBranch *b_Flag_HBHENoiseFilter_pRECO;                                                                            //!
   TBranch *b_Flag_HBHENoiseIsoFilter_pRECO;                                                                         //!
   TBranch *b_Flag_CSCTightHaloFilter_pRECO;                                                                         //!
   TBranch *b_Flag_CSCTightHaloTrkMuUnvetoFilter_pRECO;                                                              //!
   TBranch *b_Flag_CSCTightHalo2015Filter_pRECO;                                                                     //!
   TBranch *b_Flag_globalTightHalo2016Filter_pRECO;                                                                  //!
   TBranch *b_Flag_globalSuperTightHalo2016Filter_pRECO;                                                             //!
   TBranch *b_Flag_HcalStripHaloFilter_pRECO;                                                                        //!
   TBranch *b_Flag_hcalLaserEventFilter_pRECO;                                                                       //!
   TBranch *b_Flag_EcalDeadCellTriggerPrimitiveFilter_pRECO;                                                         //!
   TBranch *b_Flag_EcalDeadCellBoundaryEnergyFilter_pRECO;                                                           //!
   TBranch *b_Flag_ecalBadCalibFilter_pRECO;                                                                         //!
   TBranch *b_Flag_goodVertices_pRECO;                                                                               //!
   TBranch *b_Flag_eeBadScFilter_pRECO;                                                                              //!
   TBranch *b_Flag_ecalLaserCorrFilter_pRECO;                                                                        //!
   TBranch *b_Flag_trkPOGFilters_pRECO;                                                                              //!
   TBranch *b_Flag_chargedHadronTrackResolutionFilter_pRECO;                                                         //!
   TBranch *b_Flag_muonBadTrackFilter_pRECO;                                                                         //!
   TBranch *b_Flag_BadChargedCandidateFilter_pRECO;                                                                  //!
   TBranch *b_Flag_BadPFMuonFilter_pRECO;                                                                            //!
   TBranch *b_Flag_BadPFMuonDzFilter_pRECO;                                                                          //!
   TBranch *b_Flag_hfNoisyHitsFilter_pRECO;                                                                          //!
   TBranch *b_Flag_BadChargedCandidateSummer16Filter_pRECO;                                                          //!
   TBranch *b_Flag_BadPFMuonSummer16Filter_pRECO;                                                                    //!
   TBranch *b_Flag_trkPOG_manystripclus53X_pRECO;                                                                    //!
   TBranch *b_Flag_trkPOG_toomanystripclus53X_pRECO;                                                                 //!
   TBranch *b_Flag_trkPOG_logErrorTooManyClusters_pRECO;                                                             //!
   TBranch *b_Flag_METFilters_pRECO;                                                                                 //!
   TBranch *b_HLTriggerFirstPath;                                                                                    //!
   TBranch *b_HLT_AK8PFJet360_TrimMass30;                                                                            //!
   TBranch *b_HLT_AK8PFJet380_TrimMass30;                                                                            //!
   TBranch *b_HLT_AK8PFJet400_TrimMass30;                                                                            //!
   TBranch *b_HLT_AK8PFJet420_TrimMass30;                                                                            //!
   TBranch *b_HLT_AK8PFJet400_MassSD30;                                                                              //!
   TBranch *b_HLT_AK8PFJet420_MassSD30;                                                                              //!
   TBranch *b_HLT_AK8PFJet450_MassSD30;                                                                              //!
   TBranch *b_HLT_AK8DiPFJet250_250_MassSD30;                                                                        //!
   TBranch *b_HLT_AK8DiPFJet250_250_MassSD50;                                                                        //!
   TBranch *b_HLT_AK8DiPFJet260_260_MassSD30;                                                                        //!
   TBranch *b_HLT_AK8DiPFJet270_270_MassSD30;                                                                        //!
   TBranch *b_HLT_AK8PFHT750_TrimMass50;                                                                             //!
   TBranch *b_HLT_AK8PFHT800_TrimMass50;                                                                             //!
   TBranch *b_HLT_AK8PFHT850_TrimMass50;                                                                             //!
   TBranch *b_HLT_AK8PFHT900_TrimMass50;                                                                             //!
   TBranch *b_HLT_CaloJet500_NoJetID;                                                                                //!
   TBranch *b_HLT_CaloJet550_NoJetID;                                                                                //!
   TBranch *b_HLT_DoubleMu5_Upsilon_DoubleEle3_CaloIdL_TrackIdL;                                                     //!
   TBranch *b_HLT_DoubleMu3_DoubleEle7p5_CaloIdL_TrackIdL_Upsilon;                                                   //!
   TBranch *b_HLT_Trimuon5_3p5_2_Upsilon_Muon;                                                                       //!
   TBranch *b_HLT_TrimuonOpen_5_3p5_2_Upsilon_Muon;                                                                  //!
   TBranch *b_HLT_DoubleEle25_CaloIdL_MW;                                                                            //!
   TBranch *b_HLT_DoubleEle27_CaloIdL_MW;                                                                            //!
   TBranch *b_HLT_DoubleEle33_CaloIdL_MW;                                                                            //!
   TBranch *b_HLT_DoubleEle24_eta2p1_WPTight_Gsf;                                                                    //!
   TBranch *b_HLT_DoubleEle8_CaloIdM_TrackIdM_Mass8_DZ_PFHT350;                                                      //!
   TBranch *b_HLT_DoubleEle8_CaloIdM_TrackIdM_Mass8_PFHT350;                                                         //!
   TBranch *b_HLT_Mu27_Ele37_CaloIdL_MW;                                                                             //!
   TBranch *b_HLT_Mu37_Ele27_CaloIdL_MW;                                                                             //!
   TBranch *b_HLT_Mu37_TkMu27;                                                                                       //!
   TBranch *b_HLT_DoubleMu4_3_Bs;                                                                                    //!
   TBranch *b_HLT_DoubleMu4_3_Jpsi;                                                                                  //!
   TBranch *b_HLT_DoubleMu4_3_LowMass;                                                                               //!
   TBranch *b_HLT_DoubleMu4_LowMass_Displaced;                                                                       //!
   TBranch *b_HLT_Mu0_L1DoubleMu;                                                                                    //!
   TBranch *b_HLT_Mu4_L1DoubleMu;                                                                                    //!
   TBranch *b_HLT_DoubleMu4_3_Photon4_BsToMMG;                                                                       //!
   TBranch *b_HLT_DoubleMu4_3_Displaced_Photon4_BsToMMG;                                                             //!
   TBranch *b_HLT_DoubleMu3_Trk_Tau3mu;                                                                              //!
   TBranch *b_HLT_DoubleMu3_TkMu_DsTau3Mu;                                                                           //!
   TBranch *b_HLT_DoubleMu4_Mass3p8_DZ_PFHT350;                                                                      //!
   TBranch *b_HLT_DoubleMu4_MuMuTrk_Displaced;                                                                       //!
   TBranch *b_HLT_Mu3_PFJet40;                                                                                       //!
   TBranch *b_HLT_Mu7p5_L2Mu2_Jpsi;                                                                                  //!
   TBranch *b_HLT_Mu7p5_L2Mu2_Upsilon;                                                                               //!
   TBranch *b_HLT_Mu3_L1SingleMu5orSingleMu7;                                                                        //!
   TBranch *b_HLT_DoublePhoton33_CaloIdL;                                                                            //!
   TBranch *b_HLT_DoublePhoton70;                                                                                    //!
   TBranch *b_HLT_DoublePhoton85;                                                                                    //!
   TBranch *b_HLT_Ele15_WPLoose_Gsf;                                                                                 //!
   TBranch *b_HLT_Ele20_WPLoose_Gsf;                                                                                 //!
   TBranch *b_HLT_DiEle27_WPTightCaloOnly_L1DoubleEG;                                                                //!
   TBranch *b_HLT_Ele27_WPTight_Gsf;                                                                                 //!
   TBranch *b_HLT_Ele28_WPTight_Gsf;                                                                                 //!
   TBranch *b_HLT_Ele30_WPTight_Gsf;                                                                                 //!
   TBranch *b_HLT_Ele32_WPTight_Gsf;                                                                                 //!
   TBranch *b_HLT_Ele35_WPTight_Gsf;                                                                                 //!
   TBranch *b_HLT_Ele35_WPTight_Gsf_L1EGMT;                                                                          //!
   TBranch *b_HLT_Ele38_WPTight_Gsf;                                                                                 //!
   TBranch *b_HLT_Ele40_WPTight_Gsf;                                                                                 //!
   TBranch *b_HLT_Ele32_WPTight_Gsf_L1DoubleEG;                                                                      //!
   TBranch *b_HLT_HT300_Beamspot;                                                                                    //!
   TBranch *b_HLT_ZeroBias_Beamspot;                                                                                 //!
   TBranch *b_HLT_IsoMu20_eta2p1_TightChargedIsoPFTauHPS27_eta2p1_TightID_CrossL1;                                   //!
   TBranch *b_HLT_IsoMu27_MediumDeepTauPFTauHPS20_eta2p1_SingleL1;                                                   //!
   TBranch *b_HLT_IsoMu20;                                                                                           //!
   TBranch *b_HLT_IsoMu24;                                                                                           //!
   TBranch *b_HLT_IsoMu24_eta2p1;                                                                                    //!
   TBranch *b_HLT_IsoMu27;                                                                                           //!
   TBranch *b_HLT_UncorrectedJetE30_NoBPTX;                                                                          //!
   TBranch *b_HLT_UncorrectedJetE30_NoBPTX3BX;                                                                       //!
   TBranch *b_HLT_UncorrectedJetE60_NoBPTX3BX;                                                                       //!
   TBranch *b_HLT_UncorrectedJetE70_NoBPTX3BX;                                                                       //!
   TBranch *b_HLT_L1SingleMu18;                                                                                      //!
   TBranch *b_HLT_L1SingleMu25;                                                                                      //!
   TBranch *b_HLT_L1SingleMuCosmics;                                                                                 //!
   TBranch *b_HLT_L2Mu10_NoVertex_NoBPTX3BX;                                                                         //!
   TBranch *b_HLT_L2Mu10_NoVertex_NoBPTX;                                                                            //!
   TBranch *b_HLT_L2Mu45_NoVertex_3Sta_NoBPTX3BX;                                                                    //!
   TBranch *b_HLT_L2Mu40_NoVertex_3Sta_NoBPTX3BX;                                                                    //!
   TBranch *b_HLT_L2Mu23NoVtx_2Cha;                                                                                  //!
   TBranch *b_HLT_L2Mu23NoVtx_2Cha_CosmicSeed;                                                                       //!
   TBranch *b_HLT_DoubleL2Mu30NoVtx_2Cha_CosmicSeed_Eta2p4;                                                          //!
   TBranch *b_HLT_DoubleL2Mu30NoVtx_2Cha_Eta2p4;                                                                     //!
   TBranch *b_HLT_DoubleL2Mu50;                                                                                      //!
   TBranch *b_HLT_DoubleL2Mu23NoVtx_2Cha_CosmicSeed;                                                                 //!
   TBranch *b_HLT_DoubleL2Mu25NoVtx_2Cha_CosmicSeed;                                                                 //!
   TBranch *b_HLT_DoubleL2Mu25NoVtx_2Cha_CosmicSeed_Eta2p4;                                                          //!
   TBranch *b_HLT_DoubleL2Mu23NoVtx_2Cha;                                                                            //!
   TBranch *b_HLT_DoubleL2Mu25NoVtx_2Cha;                                                                            //!
   TBranch *b_HLT_DoubleL2Mu25NoVtx_2Cha_Eta2p4;                                                                     //!
   TBranch *b_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL;                                                                      //!
   TBranch *b_HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL;                                                                      //!
   TBranch *b_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ;                                                                   //!
   TBranch *b_HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ;                                                                   //!
   TBranch *b_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8;                                                             //!
   TBranch *b_HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_Mass8;                                                             //!
   TBranch *b_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8;                                                           //!
   TBranch *b_HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_Mass3p8;                                                           //!
   TBranch *b_HLT_Mu25_TkMu0_Onia;                                                                                   //!
   TBranch *b_HLT_Mu30_TkMu0_Psi;                                                                                    //!
   TBranch *b_HLT_Mu30_TkMu0_Upsilon;                                                                                //!
   TBranch *b_HLT_Mu20_TkMu0_Phi;                                                                                    //!
   TBranch *b_HLT_Mu25_TkMu0_Phi;                                                                                    //!
   TBranch *b_HLT_Mu15;                                                                                              //!
   TBranch *b_HLT_Mu20;                                                                                              //!
   TBranch *b_HLT_Mu27;                                                                                              //!
   TBranch *b_HLT_Mu50;                                                                                              //!
   TBranch *b_HLT_Mu55;                                                                                              //!
   TBranch *b_HLT_CascadeMu100;                                                                                      //!
   TBranch *b_HLT_HighPtTkMu100;                                                                                     //!
   TBranch *b_HLT_DiPFJetAve40;                                                                                      //!
   TBranch *b_HLT_DiPFJetAve60;                                                                                      //!
   TBranch *b_HLT_DiPFJetAve80;                                                                                      //!
   TBranch *b_HLT_DiPFJetAve140;                                                                                     //!
   TBranch *b_HLT_DiPFJetAve200;                                                                                     //!
   TBranch *b_HLT_DiPFJetAve260;                                                                                     //!
   TBranch *b_HLT_DiPFJetAve320;                                                                                     //!
   TBranch *b_HLT_DiPFJetAve400;                                                                                     //!
   TBranch *b_HLT_DiPFJetAve500;                                                                                     //!
   TBranch *b_HLT_DiPFJetAve60_HFJEC;                                                                                //!
   TBranch *b_HLT_DiPFJetAve80_HFJEC;                                                                                //!
   TBranch *b_HLT_DiPFJetAve100_HFJEC;                                                                               //!
   TBranch *b_HLT_DiPFJetAve160_HFJEC;                                                                               //!
   TBranch *b_HLT_DiPFJetAve220_HFJEC;                                                                               //!
   TBranch *b_HLT_DiPFJetAve300_HFJEC;                                                                               //!
   TBranch *b_HLT_AK8PFJet40;                                                                                        //!
   TBranch *b_HLT_AK8PFJet60;                                                                                        //!
   TBranch *b_HLT_AK8PFJet80;                                                                                        //!
   TBranch *b_HLT_AK8PFJet140;                                                                                       //!
   TBranch *b_HLT_AK8PFJet200;                                                                                       //!
   TBranch *b_HLT_AK8PFJet260;                                                                                       //!
   TBranch *b_HLT_AK8PFJet320;                                                                                       //!
   TBranch *b_HLT_AK8PFJet400;                                                                                       //!
   TBranch *b_HLT_AK8PFJet450;                                                                                       //!
   TBranch *b_HLT_AK8PFJet500;                                                                                       //!
   TBranch *b_HLT_AK8PFJet550;                                                                                       //!
   TBranch *b_HLT_PFJet40;                                                                                           //!
   TBranch *b_HLT_PFJet60;                                                                                           //!
   TBranch *b_HLT_PFJet80;                                                                                           //!
   TBranch *b_HLT_PFJet110;                                                                                          //!
   TBranch *b_HLT_PFJet140;                                                                                          //!
   TBranch *b_HLT_PFJet200;                                                                                          //!
   TBranch *b_HLT_PFJet260;                                                                                          //!
   TBranch *b_HLT_PFJet320;                                                                                          //!
   TBranch *b_HLT_PFJet400;                                                                                          //!
   TBranch *b_HLT_PFJet450;                                                                                          //!
   TBranch *b_HLT_PFJet500;                                                                                          //!
   TBranch *b_HLT_PFJet550;                                                                                          //!
   TBranch *b_HLT_PFJetFwd15;                                                                                        //!
   TBranch *b_HLT_PFJetFwd25;                                                                                        //!
   TBranch *b_HLT_PFJetFwd40;                                                                                        //!
   TBranch *b_HLT_PFJetFwd60;                                                                                        //!
   TBranch *b_HLT_PFJetFwd80;                                                                                        //!
   TBranch *b_HLT_PFJetFwd140;                                                                                       //!
   TBranch *b_HLT_PFJetFwd200;                                                                                       //!
   TBranch *b_HLT_PFJetFwd260;                                                                                       //!
   TBranch *b_HLT_PFJetFwd320;                                                                                       //!
   TBranch *b_HLT_PFJetFwd400;                                                                                       //!
   TBranch *b_HLT_PFJetFwd450;                                                                                       //!
   TBranch *b_HLT_PFJetFwd500;                                                                                       //!
   TBranch *b_HLT_AK8PFJetFwd15;                                                                                     //!
   TBranch *b_HLT_AK8PFJetFwd25;                                                                                     //!
   TBranch *b_HLT_AK8PFJetFwd40;                                                                                     //!
   TBranch *b_HLT_AK8PFJetFwd60;                                                                                     //!
   TBranch *b_HLT_AK8PFJetFwd80;                                                                                     //!
   TBranch *b_HLT_AK8PFJetFwd140;                                                                                    //!
   TBranch *b_HLT_AK8PFJetFwd200;                                                                                    //!
   TBranch *b_HLT_AK8PFJetFwd260;                                                                                    //!
   TBranch *b_HLT_AK8PFJetFwd320;                                                                                    //!
   TBranch *b_HLT_AK8PFJetFwd400;                                                                                    //!
   TBranch *b_HLT_AK8PFJetFwd450;                                                                                    //!
   TBranch *b_HLT_AK8PFJetFwd500;                                                                                    //!
   TBranch *b_HLT_PFHT180;                                                                                           //!
   TBranch *b_HLT_PFHT250;                                                                                           //!
   TBranch *b_HLT_PFHT370;                                                                                           //!
   TBranch *b_HLT_PFHT430;                                                                                           //!
   TBranch *b_HLT_PFHT510;                                                                                           //!
   TBranch *b_HLT_PFHT590;                                                                                           //!
   TBranch *b_HLT_PFHT680;                                                                                           //!
   TBranch *b_HLT_PFHT780;                                                                                           //!
   TBranch *b_HLT_PFHT890;                                                                                           //!
   TBranch *b_HLT_PFHT1050;                                                                                          //!
   TBranch *b_HLT_PFHT500_PFMET100_PFMHT100_IDTight;                                                                 //!
   TBranch *b_HLT_PFHT500_PFMET110_PFMHT110_IDTight;                                                                 //!
   TBranch *b_HLT_PFHT700_PFMET85_PFMHT85_IDTight;                                                                   //!
   TBranch *b_HLT_PFHT800_PFMET75_PFMHT75_IDTight;                                                                   //!
   TBranch *b_HLT_PFMET110_PFMHT110_IDTight;                                                                         //!
   TBranch *b_HLT_PFMET120_PFMHT120_IDTight;                                                                         //!
   TBranch *b_HLT_PFMET130_PFMHT130_IDTight;                                                                         //!
   TBranch *b_HLT_PFMET140_PFMHT140_IDTight;                                                                         //!
   TBranch *b_HLT_PFMET120_PFMHT120_IDTight_PFHT60;                                                                  //!
   TBranch *b_HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60;                                                          //!
   TBranch *b_HLT_PFMETTypeOne120_PFMHT120_IDTight_PFHT60;                                                           //!
   TBranch *b_HLT_PFMETTypeOne110_PFMHT110_IDTight;                                                                  //!
   TBranch *b_HLT_PFMETTypeOne120_PFMHT120_IDTight;                                                                  //!
   TBranch *b_HLT_PFMETTypeOne130_PFMHT130_IDTight;                                                                  //!
   TBranch *b_HLT_PFMETTypeOne140_PFMHT140_IDTight;                                                                  //!
   TBranch *b_HLT_PFMETNoMu110_PFMHTNoMu110_IDTight;                                                                 //!
   TBranch *b_HLT_PFMETNoMu120_PFMHTNoMu120_IDTight;                                                                 //!
   TBranch *b_HLT_PFMETNoMu130_PFMHTNoMu130_IDTight;                                                                 //!
   TBranch *b_HLT_PFMETNoMu140_PFMHTNoMu140_IDTight;                                                                 //!
   TBranch *b_HLT_PFMETNoMu110_PFMHTNoMu110_IDTight_FilterHF;                                                        //!
   TBranch *b_HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_FilterHF;                                                        //!
   TBranch *b_HLT_PFMETNoMu130_PFMHTNoMu130_IDTight_FilterHF;                                                        //!
   TBranch *b_HLT_PFMETNoMu140_PFMHTNoMu140_IDTight_FilterHF;                                                        //!
   TBranch *b_HLT_L1ETMHadSeeds;                                                                                     //!
   TBranch *b_HLT_CaloMHT90;                                                                                         //!
   TBranch *b_HLT_CaloMET90_NotCleaned;                                                                              //!
   TBranch *b_HLT_CaloMET350_NotCleaned;                                                                             //!
   TBranch *b_HLT_PFMET200_NotCleaned;                                                                               //!
   TBranch *b_HLT_PFMET250_NotCleaned;                                                                               //!
   TBranch *b_HLT_PFMET300_NotCleaned;                                                                               //!
   TBranch *b_HLT_PFMET200_BeamHaloCleaned;                                                                          //!
   TBranch *b_HLT_PFMETTypeOne200_BeamHaloCleaned;                                                                   //!
   TBranch *b_HLT_MET105_IsoTrk50;                                                                                   //!
   TBranch *b_HLT_MET120_IsoTrk50;                                                                                   //!
   TBranch *b_HLT_SingleJet30_Mu12_SinglePFJet40;                                                                    //!
   TBranch *b_HLT_Mu12eta2p3;                                                                                        //!
   TBranch *b_HLT_Mu12eta2p3_PFJet40;                                                                                //!
   TBranch *b_HLT_Mu12_DoublePFJets40_PFBTagDeepCSV_p71;                                                             //!
   TBranch *b_HLT_Mu12_DoublePFJets100_PFBTagDeepCSV_p71;                                                            //!
   TBranch *b_HLT_Mu12_DoublePFJets200_PFBTagDeepCSV_p71;                                                            //!
   TBranch *b_HLT_Mu12_DoublePFJets350_PFBTagDeepCSV_p71;                                                            //!
   TBranch *b_HLT_Mu12_DoublePFJets40MaxDeta1p6_DoublePFBTagDeepCSV_p71;                                             //!
   TBranch *b_HLT_Mu12_DoublePFJets54MaxDeta1p6_DoublePFBTagDeepCSV_p71;                                             //!
   TBranch *b_HLT_DoublePFJets40_PFBTagDeepCSV_p71;                                                                  //!
   TBranch *b_HLT_DoublePFJets100_PFBTagDeepCSV_p71;                                                                 //!
   TBranch *b_HLT_DoublePFJets200_PFBTagDeepCSV_p71;                                                                 //!
   TBranch *b_HLT_DoublePFJets350_PFBTagDeepCSV_p71;                                                                 //!
   TBranch *b_HLT_DoublePFJets116MaxDeta1p6_DoublePFBTagDeepCSV_p71;                                                 //!
   TBranch *b_HLT_DoublePFJets128MaxDeta1p6_DoublePFBTagDeepCSV_p71;                                                 //!
   TBranch *b_HLT_Mu12_DoublePFJets40_PFBTagDeepJet_p71;                                                             //!
   TBranch *b_HLT_Mu12_DoublePFJets100_PFBTagDeepJet_p71;                                                            //!
   TBranch *b_HLT_Mu12_DoublePFJets200_PFBTagDeepJet_p71;                                                            //!
   TBranch *b_HLT_Mu12_DoublePFJets350_PFBTagDeepJet_p71;                                                            //!
   TBranch *b_HLT_Mu12_DoublePFJets40MaxDeta1p6_DoublePFBTagDeepJet_p71;                                             //!
   TBranch *b_HLT_Mu12_DoublePFJets54MaxDeta1p6_DoublePFBTagDeepJet_p71;                                             //!
   TBranch *b_HLT_DoublePFJets40_PFBTagDeepJet_p71;                                                                  //!
   TBranch *b_HLT_DoublePFJets100_PFBTagDeepJet_p71;                                                                 //!
   TBranch *b_HLT_DoublePFJets200_PFBTagDeepJet_p71;                                                                 //!
   TBranch *b_HLT_DoublePFJets350_PFBTagDeepJet_p71;                                                                 //!
   TBranch *b_HLT_DoublePFJets116MaxDeta1p6_DoublePFBTagDeepJet_p71;                                                 //!
   TBranch *b_HLT_DoublePFJets128MaxDeta1p6_DoublePFBTagDeepJet_p71;                                                 //!
   TBranch *b_HLT_Photon300_NoHE;                                                                                    //!
   TBranch *b_HLT_Mu8_TrkIsoVVL;                                                                                     //!
   TBranch *b_HLT_Mu8_DiEle12_CaloIdL_TrackIdL_DZ;                                                                   //!
   TBranch *b_HLT_Mu8_DiEle12_CaloIdL_TrackIdL;                                                                      //!
   TBranch *b_HLT_Mu8_Ele8_CaloIdM_TrackIdM_Mass8_PFHT350_DZ;                                                        //!
   TBranch *b_HLT_Mu8_Ele8_CaloIdM_TrackIdM_Mass8_PFHT350;                                                           //!
   TBranch *b_HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ;                                                     //!
   TBranch *b_HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_PFDiJet30;                                           //!
   TBranch *b_HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_CaloDiJet30;                                         //!
   TBranch *b_HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_PFDiJet30_PFBtagDeepCSV_1p5;                         //!
   TBranch *b_HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_CaloDiJet30_CaloBtagDeepCSV_1p5;                     //!
   TBranch *b_HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL;                                                        //!
   TBranch *b_HLT_Mu17_TrkIsoVVL;                                                                                    //!
   TBranch *b_HLT_Mu19_TrkIsoVVL;                                                                                    //!
   TBranch *b_HLT_BTagMu_AK4DiJet20_Mu5;                                                                             //!
   TBranch *b_HLT_BTagMu_AK4DiJet40_Mu5;                                                                             //!
   TBranch *b_HLT_BTagMu_AK4DiJet70_Mu5;                                                                             //!
   TBranch *b_HLT_BTagMu_AK4DiJet110_Mu5;                                                                            //!
   TBranch *b_HLT_BTagMu_AK4DiJet170_Mu5;                                                                            //!
   TBranch *b_HLT_BTagMu_AK4Jet300_Mu5;                                                                              //!
   TBranch *b_HLT_BTagMu_AK8DiJet170_Mu5;                                                                            //!
   TBranch *b_HLT_BTagMu_AK8Jet170_DoubleMu5;                                                                        //!
   TBranch *b_HLT_BTagMu_AK8Jet300_Mu5;                                                                              //!
   TBranch *b_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ;                                                             //!
   TBranch *b_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL;                                                                //!
   TBranch *b_HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ;                                                    //!
   TBranch *b_HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL;                                                       //!
   TBranch *b_HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL;                                                       //!
   TBranch *b_HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ;                                                    //!
   TBranch *b_HLT_Photon20;                                                                                          //!
   TBranch *b_HLT_Photon33;                                                                                          //!
   TBranch *b_HLT_Photon50;                                                                                          //!
   TBranch *b_HLT_Photon75;                                                                                          //!
   TBranch *b_HLT_Photon90;                                                                                          //!
   TBranch *b_HLT_Photon120;                                                                                         //!
   TBranch *b_HLT_Photon150;                                                                                         //!
   TBranch *b_HLT_Photon175;                                                                                         //!
   TBranch *b_HLT_Photon200;                                                                                         //!
   TBranch *b_HLT_Photon30EB_TightID_TightIso;                                                                       //!
   TBranch *b_HLT_Photon110EB_TightID_TightIso;                                                                      //!
   TBranch *b_HLT_Photon100EBHE10;                                                                                   //!
   TBranch *b_HLT_Photon50_R9Id90_HE10_IsoM;                                                                         //!
   TBranch *b_HLT_Photon75_R9Id90_HE10_IsoM;                                                                         //!
   TBranch *b_HLT_Photon90_R9Id90_HE10_IsoM;                                                                         //!
   TBranch *b_HLT_Photon120_R9Id90_HE10_IsoM;                                                                        //!
   TBranch *b_HLT_Photon165_R9Id90_HE10_IsoM;                                                                        //!
   TBranch *b_HLT_Diphoton30_22_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass90;                                                //!
   TBranch *b_HLT_Diphoton30_22_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass95;                                                //!
   TBranch *b_HLT_Photon35_TwoProngs35;                                                                              //!
   TBranch *b_HLT_IsoMu24_TwoProngs35;                                                                               //!
   TBranch *b_HLT_Dimuon0_Jpsi_L1_NoOS;                                                                              //!
   TBranch *b_HLT_Dimuon0_Jpsi_NoVertexing_NoOS;                                                                     //!
   TBranch *b_HLT_Dimuon0_Jpsi;                                                                                      //!
   TBranch *b_HLT_Dimuon0_Jpsi_NoVertexing;                                                                          //!
   TBranch *b_HLT_Dimuon0_Jpsi_L1_4R_0er1p5R;                                                                        //!
   TBranch *b_HLT_Dimuon0_Jpsi_NoVertexing_L1_4R_0er1p5R;                                                            //!
   TBranch *b_HLT_Dimuon0_Jpsi3p5_Muon2;                                                                             //!
   TBranch *b_HLT_Dimuon0_Upsilon_L1_4p5;                                                                            //!
   TBranch *b_HLT_Dimuon0_Upsilon_L1_5;                                                                              //!
   TBranch *b_HLT_Dimuon0_Upsilon_L1_4p5NoOS;                                                                        //!
   TBranch *b_HLT_Dimuon0_Upsilon_L1_4p5er2p0;                                                                       //!
   TBranch *b_HLT_Dimuon0_Upsilon_L1_4p5er2p0M;                                                                      //!
   TBranch *b_HLT_Dimuon0_Upsilon_NoVertexing;                                                                       //!
   TBranch *b_HLT_Dimuon0_Upsilon_L1_5M;                                                                             //!
   TBranch *b_HLT_Dimuon0_LowMass_L1_0er1p5R;                                                                        //!
   TBranch *b_HLT_Dimuon0_LowMass_L1_0er1p5;                                                                         //!
   TBranch *b_HLT_Dimuon0_LowMass;                                                                                   //!
   TBranch *b_HLT_Dimuon0_LowMass_L1_4;                                                                              //!
   TBranch *b_HLT_Dimuon0_LowMass_L1_4R;                                                                             //!
   TBranch *b_HLT_Dimuon0_LowMass_L1_TM530;                                                                          //!
   TBranch *b_HLT_Dimuon0_Upsilon_Muon_L1_TM0;                                                                       //!
   TBranch *b_HLT_Dimuon0_Upsilon_Muon_NoL1Mass;                                                                     //!
   TBranch *b_HLT_TripleMu_5_3_3_Mass3p8_DZ;                                                                         //!
   TBranch *b_HLT_TripleMu_10_5_5_DZ;                                                                                //!
   TBranch *b_HLT_TripleMu_12_10_5;                                                                                  //!
   TBranch *b_HLT_Tau3Mu_Mu7_Mu1_TkMu1_Tau15;                                                                        //!
   TBranch *b_HLT_Tau3Mu_Mu7_Mu1_TkMu1_Tau15_Charge1;                                                                //!
   TBranch *b_HLT_Tau3Mu_Mu7_Mu1_TkMu1_IsoTau15;                                                                     //!
   TBranch *b_HLT_Tau3Mu_Mu7_Mu1_TkMu1_IsoTau15_Charge1;                                                             //!
   TBranch *b_HLT_DoubleMu3_DZ_PFMET50_PFMHT60;                                                                      //!
   TBranch *b_HLT_DoubleMu3_DZ_PFMET70_PFMHT70;                                                                      //!
   TBranch *b_HLT_DoubleMu3_DZ_PFMET90_PFMHT90;                                                                      //!
   TBranch *b_HLT_DoubleMu3_Trk_Tau3mu_NoL1Mass;                                                                     //!
   TBranch *b_HLT_DoubleMu4_Jpsi_Displaced;                                                                          //!
   TBranch *b_HLT_DoubleMu4_Jpsi_NoVertexing;                                                                        //!
   TBranch *b_HLT_DoubleMu4_JpsiTrkTrk_Displaced;                                                                    //!
   TBranch *b_HLT_DoubleMu4_JpsiTrk_Bc;                                                                              //!
   TBranch *b_HLT_DoubleMu43NoFiltersNoVtx;                                                                          //!
   TBranch *b_HLT_DoubleMu48NoFiltersNoVtx;                                                                          //!
   TBranch *b_HLT_Mu43NoFiltersNoVtx_Photon43_CaloIdL;                                                               //!
   TBranch *b_HLT_Mu48NoFiltersNoVtx_Photon48_CaloIdL;                                                               //!
   TBranch *b_HLT_Mu38NoFiltersNoVtxDisplaced_Photon38_CaloIdL;                                                      //!
   TBranch *b_HLT_Mu43NoFiltersNoVtxDisplaced_Photon43_CaloIdL;                                                      //!
   TBranch *b_HLT_HT425;                                                                                             //!
   TBranch *b_HLT_HT430_DisplacedDijet40_DisplacedTrack;                                                             //!
   TBranch *b_HLT_HT500_DisplacedDijet40_DisplacedTrack;                                                             //!
   TBranch *b_HLT_HT430_DisplacedDijet60_DisplacedTrack;                                                             //!
   TBranch *b_HLT_HT400_DisplacedDijet40_DisplacedTrack;                                                             //!
   TBranch *b_HLT_HT650_DisplacedDijet60_Inclusive;                                                                  //!
   TBranch *b_HLT_HT550_DisplacedDijet60_Inclusive;                                                                  //!
   TBranch *b_HLT_DiJet110_35_Mjj650_PFMET110;                                                                       //!
   TBranch *b_HLT_DiJet110_35_Mjj650_PFMET120;                                                                       //!
   TBranch *b_HLT_DiJet110_35_Mjj650_PFMET130;                                                                       //!
   TBranch *b_HLT_TripleJet110_35_35_Mjj650_PFMET110;                                                                //!
   TBranch *b_HLT_TripleJet110_35_35_Mjj650_PFMET120;                                                                //!
   TBranch *b_HLT_TripleJet110_35_35_Mjj650_PFMET130;                                                                //!
   TBranch *b_HLT_Ele30_eta2p1_WPTight_Gsf_CentralPFJet35_EleCleaned;                                                //!
   TBranch *b_HLT_Ele28_eta2p1_WPTight_Gsf_HT150;                                                                    //!
   TBranch *b_HLT_Ele28_HighEta_SC20_Mass55;                                                                         //!
   TBranch *b_HLT_Ele15_IsoVVVL_PFHT450_CaloBTagDeepCSV_4p5;                                                         //!
   TBranch *b_HLT_Ele15_IsoVVVL_PFHT450_PFMET50;                                                                     //!
   TBranch *b_HLT_Ele15_IsoVVVL_PFHT450;                                                                             //!
   TBranch *b_HLT_Ele50_IsoVVVL_PFHT450;                                                                             //!
   TBranch *b_HLT_Ele15_IsoVVVL_PFHT600;                                                                             //!
   TBranch *b_HLT_Mu4_TrkIsoVVL_DiPFJet90_40_DEta3p5_MJJ750_HTT300_PFMETNoMu60;                                      //!
   TBranch *b_HLT_Mu8_TrkIsoVVL_DiPFJet40_DEta3p5_MJJ750_HTT300_PFMETNoMu60;                                         //!
   TBranch *b_HLT_Mu10_TrkIsoVVL_DiPFJet40_DEta3p5_MJJ750_HTT350_PFMETNoMu60;                                        //!
   TBranch *b_HLT_Mu15_IsoVVVL_PFHT450_CaloBTagDeepCSV_4p5;                                                          //!
   TBranch *b_HLT_Mu15_IsoVVVL_PFHT450_PFMET50;                                                                      //!
   TBranch *b_HLT_Mu15_IsoVVVL_PFHT450;                                                                              //!
   TBranch *b_HLT_Mu50_IsoVVVL_PFHT450;                                                                              //!
   TBranch *b_HLT_Mu15_IsoVVVL_PFHT600;                                                                              //!
   TBranch *b_HLT_Mu3er1p5_PFJet100er2p5_PFMET70_PFMHT70_IDTight;                                                    //!
   TBranch *b_HLT_Mu3er1p5_PFJet100er2p5_PFMET80_PFMHT80_IDTight;                                                    //!
   TBranch *b_HLT_Mu3er1p5_PFJet100er2p5_PFMET90_PFMHT90_IDTight;                                                    //!
   TBranch *b_HLT_Mu3er1p5_PFJet100er2p5_PFMET100_PFMHT100_IDTight;                                                  //!
   TBranch *b_HLT_Mu3er1p5_PFJet100er2p5_PFMETNoMu70_PFMHTNoMu70_IDTight;                                            //!
   TBranch *b_HLT_Mu3er1p5_PFJet100er2p5_PFMETNoMu80_PFMHTNoMu80_IDTight;                                            //!
   TBranch *b_HLT_Mu3er1p5_PFJet100er2p5_PFMETNoMu90_PFMHTNoMu90_IDTight;                                            //!
   TBranch *b_HLT_Mu3er1p5_PFJet100er2p5_PFMETNoMu100_PFMHTNoMu100_IDTight;                                          //!
   TBranch *b_HLT_Dimuon10_PsiPrime_Barrel_Seagulls;                                                                 //!
   TBranch *b_HLT_Dimuon20_Jpsi_Barrel_Seagulls;                                                                     //!
   TBranch *b_HLT_Dimuon10_Upsilon_y1p4;                                                                             //!
   TBranch *b_HLT_Dimuon12_Upsilon_y1p4;                                                                             //!
   TBranch *b_HLT_Dimuon14_Phi_Barrel_Seagulls;                                                                      //!
   TBranch *b_HLT_Dimuon25_Jpsi;                                                                                     //!
   TBranch *b_HLT_Dimuon14_PsiPrime;                                                                                 //!
   TBranch *b_HLT_Dimuon14_PsiPrime_noCorrL1;                                                                        //!
   TBranch *b_HLT_Dimuon18_PsiPrime;                                                                                 //!
   TBranch *b_HLT_Dimuon18_PsiPrime_noCorrL1;                                                                        //!
   TBranch *b_HLT_Dimuon24_Upsilon_noCorrL1;                                                                         //!
   TBranch *b_HLT_Dimuon24_Phi_noCorrL1;                                                                             //!
   TBranch *b_HLT_Dimuon25_Jpsi_noCorrL1;                                                                            //!
   TBranch *b_HLT_DiMu4_Ele9_CaloIdL_TrackIdL_DZ_Mass3p8;                                                            //!
   TBranch *b_HLT_DiMu9_Ele9_CaloIdL_TrackIdL_DZ;                                                                    //!
   TBranch *b_HLT_DiMu9_Ele9_CaloIdL_TrackIdL;                                                                       //!
   TBranch *b_HLT_DoubleIsoMu20_eta2p1;                                                                              //!
   TBranch *b_HLT_TrkMu12_DoubleTrkMu5NoFiltersNoVtx;                                                                //!
   TBranch *b_HLT_Mu8;                                                                                               //!
   TBranch *b_HLT_Mu17;                                                                                              //!
   TBranch *b_HLT_Mu19;                                                                                              //!
   TBranch *b_HLT_Mu17_Photon30_IsoCaloId;                                                                           //!
   TBranch *b_HLT_Ele8_CaloIdL_TrackIdL_IsoVL_PFJet30;                                                               //!
   TBranch *b_HLT_Ele12_CaloIdL_TrackIdL_IsoVL_PFJet30;                                                              //!
   TBranch *b_HLT_Ele23_CaloIdL_TrackIdL_IsoVL_PFJet30;                                                              //!
   TBranch *b_HLT_Ele8_CaloIdM_TrackIdM_PFJet30;                                                                     //!
   TBranch *b_HLT_Ele17_CaloIdM_TrackIdM_PFJet30;                                                                    //!
   TBranch *b_HLT_Ele23_CaloIdM_TrackIdM_PFJet30;                                                                    //!
   TBranch *b_HLT_Ele50_CaloIdVT_GsfTrkIdT_PFJet165;                                                                 //!
   TBranch *b_HLT_Ele115_CaloIdVT_GsfTrkIdT;                                                                         //!
   TBranch *b_HLT_Ele135_CaloIdVT_GsfTrkIdT;                                                                         //!
   TBranch *b_HLT_PFHT330PT30_QuadPFJet_75_60_45_40_TriplePFBTagDeepCSV_4p5;                                         //!
   TBranch *b_HLT_PFHT330PT30_QuadPFJet_75_60_45_40;                                                                 //!
   TBranch *b_HLT_PFHT400_SixPFJet32_DoublePFBTagDeepCSV_2p94;                                                       //!
   TBranch *b_HLT_PFHT400_SixPFJet32;                                                                                //!
   TBranch *b_HLT_PFHT450_SixPFJet36_PFBTagDeepCSV_1p59;                                                             //!
   TBranch *b_HLT_PFHT450_SixPFJet36;                                                                                //!
   TBranch *b_HLT_PFHT400_FivePFJet_100_100_60_30_30;                                                                //!
   TBranch *b_HLT_PFHT400_FivePFJet_100_100_60_30_30_DoublePFBTagDeepCSV_4p5;                                        //!
   TBranch *b_HLT_PFHT400_FivePFJet_120_120_60_30_30_DoublePFBTagDeepCSV_4p5;                                        //!
   TBranch *b_HLT_PFHT350;                                                                                           //!
   TBranch *b_HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_PFHT350MinPFJet15;                                       //!
   TBranch *b_HLT_ECALHT800;                                                                                         //!
   TBranch *b_HLT_DiSC30_18_EIso_AND_HE_Mass70;                                                                      //!
   TBranch *b_HLT_Physics;                                                                                           //!
   TBranch *b_HLT_Random;                                                                                            //!
   TBranch *b_HLT_ZeroBias;                                                                                          //!
   TBranch *b_HLT_ZeroBias_Alignment;                                                                                //!
   TBranch *b_HLT_Photon20_HoverELoose;                                                                              //!
   TBranch *b_HLT_Photon30_HoverELoose;                                                                              //!
   TBranch *b_HLT_EcalCalibration;                                                                                   //!
   TBranch *b_HLT_HcalCalibration;                                                                                   //!
   TBranch *b_HLT_L1UnpairedBunchBptxMinus;                                                                          //!
   TBranch *b_HLT_L1UnpairedBunchBptxPlus;                                                                           //!
   TBranch *b_HLT_L1NotBptxOR;                                                                                       //!
   TBranch *b_HLT_L1_CDC_SingleMu_3_er1p2_TOP120_DPHI2p618_3p142;                                                    //!
   TBranch *b_HLT_CDC_L2cosmic_10_er1p0;                                                                             //!
   TBranch *b_HLT_CDC_L2cosmic_5p5_er1p0;                                                                            //!
   TBranch *b_HLT_HcalNZS;                                                                                           //!
   TBranch *b_HLT_HcalPhiSym;                                                                                        //!
   TBranch *b_HLT_HcalIsolatedbunch;                                                                                 //!
   TBranch *b_HLT_IsoTrackHB;                                                                                        //!
   TBranch *b_HLT_IsoTrackHE;                                                                                        //!
   TBranch *b_HLT_ZeroBias_FirstCollisionAfterAbortGap;                                                              //!
   TBranch *b_HLT_ZeroBias_IsolatedBunches;                                                                          //!
   TBranch *b_HLT_ZeroBias_FirstCollisionInTrain;                                                                    //!
   TBranch *b_HLT_ZeroBias_LastCollisionInTrain;                                                                     //!
   TBranch *b_HLT_ZeroBias_FirstBXAfterTrain;                                                                        //!
   TBranch *b_HLT_Ele16_Ele12_Ele8_CaloIdL_TrackIdL;                                                                 //!
   TBranch *b_HLT_VBF_DoubleMediumDeepTauPFTauHPS20_eta2p1;                                                          //!
   TBranch *b_HLT_Photon75_R9Id90_HE10_IsoM_EBOnly_PFJetsMJJ300DEta3;                                                //!
   TBranch *b_HLT_PFMET100_PFMHT100_IDTight_PFHT60;                                                                  //!
   TBranch *b_HLT_PFMETNoMu100_PFMHTNoMu100_IDTight_PFHT60;                                                          //!
   TBranch *b_HLT_PFMETTypeOne100_PFMHT100_IDTight_PFHT60;                                                           //!
   TBranch *b_HLT_Mu18_Mu9_SameSign;                                                                                 //!
   TBranch *b_HLT_DoubleMu2_Jpsi_DoubleTrk1_Phi1p05;                                                                 //!
   TBranch *b_HLT_DoubleMu3_DCA_PFMET50_PFMHT60;                                                                     //!
   TBranch *b_HLT_TripleMu_5_3_3_Mass3p8_DCA;                                                                        //!
   TBranch *b_HLT_QuadPFJet103_88_75_15_DoublePFBTagDeepCSV_1p3_7p7_VBF1;                                            //!
   TBranch *b_HLT_QuadPFJet111_90_80_15_DoublePFBTagDeepCSV_1p3_7p7_VBF1;                                            //!
   TBranch *b_HLT_QuadPFJet103_88_75_15_PFBTagDeepCSV_1p3_VBF2;                                                      //!
   TBranch *b_HLT_QuadPFJet105_88_76_15_PFBTagDeepCSV_1p3_VBF2;                                                      //!
   TBranch *b_HLT_QuadPFJet111_90_80_15_PFBTagDeepCSV_1p3_VBF2;                                                      //!
   TBranch *b_HLT_QuadPFJet103_88_75_15;                                                                             //!
   TBranch *b_HLT_QuadPFJet105_88_76_15;                                                                             //!
   TBranch *b_HLT_QuadPFJet111_90_80_15;                                                                             //!
   TBranch *b_HLT_AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_p02;                                                    //!
   TBranch *b_HLT_AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_np2;                                                    //!
   TBranch *b_HLT_AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_np4;                                                    //!
   TBranch *b_HLT_Diphoton30_18_R9IdL_AND_HE_AND_IsoCaloId_Mass55;                                                   //!
   TBranch *b_HLT_Diphoton30_18_R9IdL_AND_HE_AND_IsoCaloId;                                                          //!
   TBranch *b_HLT_Mu12_IP6;                                                                                          //!
   TBranch *b_HLT_QuadPFJet105_88_76_15_DoublePFBTagDeepCSV_1p3_7p7_VBF1;                                            //!
   TBranch *b_HLT_DoubleMediumDeepTauPFTauHPS35_L2NN_eta2p1;                                                         //!
   TBranch *b_HLT_Ele24_eta2p1_WPTight_Gsf_LooseDeepTauPFTauHPS30_eta2p1_CrossL1;                                    //!
   TBranch *b_HLT_IsoMu20_eta2p1_LooseDeepTauPFTauHPS27_eta2p1_CrossL1;                                              //!
   TBranch *b_HLT_IsoMu24_eta2p1_LooseDeepTauPFTauHPS180_eta2p1;                                                     //!
   TBranch *b_HLT_IsoMu24_eta2p1_LooseDeepTauPFTauHPS30_eta2p1_CrossL1;                                              //!
   TBranch *b_HLT_IsoMu24_eta2p1_MediumDeepTauPFTauHPS35_L2NN_eta2p1_CrossL1;                                        //!
   TBranch *b_HLT_LooseDeepTauPFTauHPS180_L2NN_eta2p1;                                                               //!
   TBranch *b_HLT_PFHT330PT30_QuadPFJet_75_60_45_40_TriplePFBTagDeepJet_4p5;                                         //!
   TBranch *b_HLT_PFHT400_FivePFJet_100_100_60_30_30_DoublePFBTagDeepJet_4p5;                                        //!
   TBranch *b_HLT_PFHT400_FivePFJet_120_120_60_30_30_DoublePFBTagDeepJet_4p5;                                        //!
   TBranch *b_HLT_PFHT400_SixPFJet32_DoublePFBTagDeepJet_2p94;                                                       //!
   TBranch *b_HLT_PFHT450_SixPFJet36_PFBTagDeepJet_1p59;                                                             //!
   TBranch *b_HLT_QuadPFJet103_88_75_15_DoublePFBTagDeepJet_1p3_7p7_VBF1;                                            //!
   TBranch *b_HLT_QuadPFJet103_88_75_15_PFBTagDeepJet_1p3_VBF2;                                                      //!
   TBranch *b_HLT_QuadPFJet105_88_76_15_DoublePFBTagDeepJet_1p3_7p7_VBF1;                                            //!
   TBranch *b_HLT_QuadPFJet105_88_76_15_PFBTagDeepJet_1p3_VBF2;                                                      //!
   TBranch *b_HLT_QuadPFJet111_90_80_15_DoublePFBTagDeepJet_1p3_7p7_VBF1;                                            //!
   TBranch *b_HLT_QuadPFJet111_90_80_15_PFBTagDeepJet_1p3_VBF2;                                                      //!
   TBranch *b_HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_PFDiJet30_PFBtagDeepJet_1p5;                         //!
   TBranch *b_HLT_QuadPFJet70_50_40_30;                                                                              //!
   TBranch *b_HLT_QuadPFJet70_50_40_30_PFBTagParticleNet_2BTagSum0p65;                                               //!
   TBranch *b_HLT_QuadPFJet70_50_40_35_PFBTagParticleNet_2BTagSum0p65;                                               //!
   TBranch *b_HLT_QuadPFJet70_50_45_35_PFBTagParticleNet_2BTagSum0p65;                                               //!
   TBranch *b_HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_PFDiJet30_PFBTagParticleNet_2BTagSum0p65;            //!
   TBranch *b_HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_QuadPFJet70_50_40_30;                                //!
   TBranch *b_HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_QuadPFJet70_50_40_30_PFBTagParticleNet_2BTagSum0p65; //!
   TBranch *b_HLT_AK8PFJet230_SoftDropMass40;                                                                        //!
   TBranch *b_HLT_AK8PFJet230_SoftDropMass40_PFAK8ParticleNetBB0p35;                                                 //!
   TBranch *b_HLT_AK8PFJet250_SoftDropMass40_PFAK8ParticleNetBB0p35;                                                 //!
   TBranch *b_HLT_AK8PFJet275_SoftDropMass40_PFAK8ParticleNetBB0p35;                                                 //!
   TBranch *b_HLT_AK8PFJet400_SoftDropMass40;                                                                        //!
   TBranch *b_HLT_AK8PFJet425_SoftDropMass40;                                                                        //!
   TBranch *b_HLT_AK8PFJet450_SoftDropMass40;                                                                        //!
   TBranch *b_HLT_AK8PFJet230_SoftDropMass40_PFAK8ParticleNetTauTau0p30;                                             //!
   TBranch *b_HLT_AK8PFJet250_SoftDropMass40_PFAK8ParticleNetTauTau0p30;                                             //!
   TBranch *b_HLT_AK8PFJet275_SoftDropMass40_PFAK8ParticleNetTauTau0p30;                                             //!
   TBranch *b_HLT_IsoMu50_AK8PFJet230_SoftDropMass40;                                                                //!
   TBranch *b_HLT_IsoMu50_AK8PFJet230_SoftDropMass40_PFAK8ParticleNetBB0p35;                                         //!
   TBranch *b_HLT_Ele50_CaloIdVT_GsfTrkIdT_AK8PFJet230_SoftDropMass40;                                               //!
   TBranch *b_HLT_Ele50_CaloIdVT_GsfTrkIdT_AK8PFJet230_SoftDropMass40_PFAK8ParticleNetBB0p35;                        //!
   TBranch *b_HLT_DoubleMediumDeepTauPFTauHPS30_L2NN_eta2p1_PFJet60;                                                 //!
   TBranch *b_HLT_DoubleMediumDeepTauPFTauHPS30_L2NN_eta2p1_PFJet75;                                                 //!
   TBranch *b_HLT_IsoMu24_eta2p1_MediumDeepTauPFTauHPS30_L2NN_eta2p1_CrossL1;                                        //!
   TBranch *b_HLT_IsoMu24_eta2p1_MediumDeepTauPFTauHPS30_L2NN_eta2p1_PFJet60_CrossL1;                                //!
   TBranch *b_HLT_IsoMu24_eta2p1_MediumDeepTauPFTauHPS30_L2NN_eta2p1_PFJet75_CrossL1;                                //!
   TBranch *b_HLT_DoublePFJets40_Mass500_MediumDeepTauPFTauHPS45_L2NN_MediumDeepTauPFTauHPS20_eta2p1;                //!
   TBranch *b_HLT_IsoMu24_eta2p1_MediumDeepTauPFTauHPS20_eta2p1_SingleL1;                                            //!
   TBranch *b_HLT_IsoMu24_eta2p1_MediumDeepTauPFTauHPS45_L2NN_eta2p1_CrossL1;                                        //!
   TBranch *b_HLT_DoubleL2Mu10NoVtx_2Cha_VetoL3Mu0DxyMax1cm;                                                         //!
   TBranch *b_HLT_DoubleL2Mu12NoVtx_2Cha_VetoL3Mu0DxyMax1cm;                                                         //!
   TBranch *b_HLT_DoubleL2Mu14NoVtx_2Cha_VetoL3Mu0DxyMax1cm;                                                         //!
   TBranch *b_HLT_DoubleL3Mu16_10NoVtx_DxyMin0p01cm;                                                                 //!
   TBranch *b_HLT_DoubleL3Mu18_10NoVtx_DxyMin0p01cm;                                                                 //!
   TBranch *b_HLT_DoubleL3Mu20_10NoVtx_DxyMin0p01cm;                                                                 //!
   TBranch *b_HLT_L2Mu10NoVtx_2Cha;                                                                                  //!
   TBranch *b_HLT_L2Mu10NoVtx_2Cha_VetoL3Mu0DxyMax1cm;                                                               //!
   TBranch *b_HLT_L3Mu10NoVtx;                                                                                       //!
   TBranch *b_HLT_L3Mu10NoVtx_DxyMin0p01cm;                                                                          //!
   TBranch *b_HLT_DoubleL2Mu_L3Mu16NoVtx_VetoL3Mu0DxyMax0p1cm;                                                       //!
   TBranch *b_HLT_DoubleL2Mu_L3Mu18NoVtx_VetoL3Mu0DxyMax0p1cm;                                                       //!
   TBranch *b_HLT_DoubleL2Mu10NoVtx_2Cha_CosmicSeed_VetoL3Mu0DxyMax1cm;                                              //!
   TBranch *b_HLT_DoubleL2Mu12NoVtx_2Cha_CosmicSeed_VetoL3Mu0DxyMax1cm;                                              //!
   TBranch *b_HLT_L2Mu10NoVtx_2Cha_CosmicSeed;                                                                       //!
   TBranch *b_HLT_L2Mu10NoVtx_2Cha_CosmicSeed_VetoL3Mu0DxyMax1cm;                                                    //!
   TBranch *b_HLT_DoubleL3dTksMu16_10NoVtx_DxyMin0p01cm;                                                             //!
   TBranch *b_HLT_L3dTksMu10_NoVtx_DxyMin0p01cm;                                                                     //!
   TBranch *b_HLT_Mu20NoFiltersNoVtxDisplaced_Photon20_CaloCustomId;                                                 //!
   TBranch *b_HLT_DoubleMediumChargedIsoDisplacedPFTauHPS32_Trk1_eta2p1;                                             //!
   TBranch *b_HLT_HT430_DelayedJet40_DoubleDelay0p5nsTrackless;                                                      //!
   TBranch *b_HLT_HT430_DelayedJet40_DoubleDelay1nsInclusive;                                                        //!
   TBranch *b_HLT_HT430_DelayedJet40_SingleDelay1nsTrackless;                                                        //!
   TBranch *b_HLT_HT430_DelayedJet40_SingleDelay2nsInclusive;                                                        //!
   TBranch *b_HLT_L1Mu6HT240;                                                                                        //!
   TBranch *b_HLT_Mu6HT240_DisplacedDijet30_Inclusive0PtrkShortSig5;                                                 //!
   TBranch *b_HLT_Mu6HT240_DisplacedDijet30_Inclusive1PtrkShortSig5_DisplacedLoose;                                  //!
   TBranch *b_HLT_Mu6HT240_DisplacedDijet35_Inclusive0PtrkShortSig5;                                                 //!
   TBranch *b_HLT_Mu6HT240_DisplacedDijet35_Inclusive1PtrkShortSig5_DisplacedLoose;                                  //!
   TBranch *b_HLT_Mu6HT240_DisplacedDijet40_Inclusive0PtrkShortSig5;                                                 //!
   TBranch *b_HLT_Mu6HT240_DisplacedDijet40_Inclusive1PtrkShortSig5_DisplacedLoose;                                  //!
   TBranch *b_HLT_HT430_DisplacedDijet30_Inclusive1PtrkShortSig5;                                                    //!
   TBranch *b_HLT_HT430_DisplacedDijet35_Inclusive1PtrkShortSig5;                                                    //!
   TBranch *b_HLT_HT430_DisplacedDijet40_Inclusive1PtrkShortSig5;                                                    //!
   TBranch *b_HLT_CaloMET60_DTCluster50;                                                                             //!
   TBranch *b_HLT_CaloMET60_DTClusterNoMB1S50;                                                                       //!
   TBranch *b_HLT_L1MET_DTCluster50;                                                                                 //!
   TBranch *b_HLT_L1MET_DTClusterNoMB1S50;                                                                           //!
   TBranch *b_HLT_CscCluster_Loose;                                                                                  //!
   TBranch *b_HLT_CscCluster_Medium;                                                                                 //!
   TBranch *b_HLT_CscCluster_Tight;                                                                                  //!
   TBranch *b_HLT_L1CSCShower_DTCluster50;                                                                           //!
   TBranch *b_HLT_L1CSCShower_DTCluster75;                                                                           //!
   TBranch *b_HLT_DoubleCscCluster75;
   TBranch *b_HLT_DoubleCscCluster100;
   TBranch *b_HLT_PFMET105_IsoTrk50;                                                                                 //!
   TBranch *b_HLT_PFMET110_PFJet100;                                                                                 //!
   TBranch *b_HLT_HT170_L1SingleLLPJet_DisplacedDijet40_DisplacedTrack;                                              //!
   TBranch *b_HLT_HT200_L1SingleLLPJet_DisplacedDijet40_DisplacedTrack;                                              //!
   TBranch *b_HLT_HT200_L1SingleLLPJet_DisplacedDijet60_DisplacedTrack;                                              //!
   TBranch *b_HLT_HT270_L1SingleLLPJet_DisplacedDijet40_DisplacedTrack;                                              //!
   TBranch *b_HLT_HT320_L1SingleLLPJet_DisplacedDijet60_Inclusive;                                                   //!
   TBranch *b_HLT_HT420_L1SingleLLPJet_DisplacedDijet60_Inclusive;                                                   //!
   TBranch *b_HLT_HT200_L1SingleLLPJet_DelayedJet40_SingleDelay1nsTrackless;                                         //!
   TBranch *b_HLT_HT200_L1SingleLLPJet_DelayedJet40_SingleDelay2nsInclusive;                                         //!
   TBranch *b_HLT_HT200_L1SingleLLPJet_DelayedJet40_DoubleDelay0p5nsTrackless;                                       //!
   TBranch *b_HLT_HT200_L1SingleLLPJet_DelayedJet40_DoubleDelay1nsInclusive;                                         //!
   TBranch *b_HLT_HT200_L1SingleLLPJet_DisplacedDijet30_Inclusive1PtrkShortSig5;                                     //!
   TBranch *b_HLT_HT200_L1SingleLLPJet_DisplacedDijet35_Inclusive1PtrkShortSig5;                                     //!
   TBranch *b_HLT_HT200_L1SingleLLPJet_DisplacedDijet40_Inclusive1PtrkShortSig5;                                     //!
   TBranch *b_HLT_DiPhoton10Time1p4ns;                                                                               //!
   TBranch *b_HLT_DiPhoton10Time1p6ns;                                                                               //!
   TBranch *b_HLT_DiPhoton10Time1p8ns;                                                                               //!
   TBranch *b_HLT_DiPhoton10Time2ns;                                                                                 //!
   TBranch *b_HLT_DiPhoton10sminlt0p1;                                                                               //!
   TBranch *b_HLT_DiPhoton10sminlt0p12;                                                                              //!
   TBranch *b_HLT_DiPhoton10_CaloIdL;                                                                                //!
   TBranch *b_HLT_DoubleEle4_eta1p22_mMax6;                                                                          //!
   TBranch *b_HLT_DoubleEle4p5_eta1p22_mMax6;                                                                        //!
   TBranch *b_HLT_DoubleEle5_eta1p22_mMax6;                                                                          //!
   TBranch *b_HLT_DoubleEle5p5_eta1p22_mMax6;                                                                        //!
   TBranch *b_HLT_DoubleEle6_eta1p22_mMax6;                                                                          //!
   TBranch *b_HLT_DoubleEle6p5_eta1p22_mMax6;                                                                        //!
   TBranch *b_HLT_DoubleEle7_eta1p22_mMax6;                                                                          //!
   TBranch *b_HLT_DoubleEle7p5_eta1p22_mMax6;                                                                        //!
   TBranch *b_HLT_DoubleEle8_eta1p22_mMax6;                                                                          //!
   TBranch *b_HLT_DoubleEle8p5_eta1p22_mMax6;                                                                        //!
   TBranch *b_HLT_DoubleEle9_eta1p22_mMax6;                                                                          //!
   TBranch *b_HLT_DoubleEle9p5_eta1p22_mMax6;                                                                        //!
   TBranch *b_HLT_DoubleEle10_eta1p22_mMax6;                                                                         //!
   TBranch *b_HLT_Diphoton20_14_eta1p5_R9IdL_AND_HE_AND_IsoTCaloIdT;                                                 //!
   TBranch *b_HLT_Diphoton20_14_eta1p5_R9IdL_AND_HET_AND_IsoTCaloIdT;                                                //!
   TBranch *b_HLT_Diphoton22_14_eta1p5_R9IdL_AND_HE_AND_IsoTCaloIdT;                                                 //!
   TBranch *b_HLT_Diphoton22_14_eta1p5_R9IdL_AND_HET_AND_IsoTCaloIdT;                                                //!
   TBranch *b_HLT_ExpressMuons;                                                                                      //!
   TBranch *b_HLT_OnlineMonitorGroup;                                                                                //!
   TBranch *b_HLT_PPSMaxTracksPerArm1;                                                                               //!
   TBranch *b_HLT_PPSMaxTracksPerRP4;                                                                                //!
   TBranch *b_HLTriggerFinalPath;                                                                                    //!
   TBranch *b_HLT_EphemeralPhysics;                                                                                  //!
   TBranch *b_HLT_EphemeralZeroBias;                                                                                 //!
   TBranch *b_L1_DoubleMu3_SQ_ETMHF30_HTT60er;                                                                       //!
   TBranch *b_L1_DoubleMu3_SQ_ETMHF30_Jet60er2p5_OR_DoubleJet40er2p5;                                                //!
   TBranch *b_L1_DoubleMu3_SQ_ETMHF40_HTT60er;                                                                       //!
   TBranch *b_L1_DoubleMu3_SQ_ETMHF40_Jet60er2p5_OR_DoubleJet40er2p5;                                                //!
   TBranch *b_L1_ETMHF80_SingleJet55er2p5_dPhi_Min2p1;                                                               //!
   TBranch *b_L1_ETMHF80_SingleJet55er2p5_dPhi_Min2p6;                                                               //!
   TBranch *b_L1_Mu3er1p5_Jet100er2p5_ETMHF30;                                                                       //!
   TBranch *b_nCscSeg;                                                                                               //!
   TBranch *b_nCscRechits;                                                                                           //!
   TBranch *b_nDtSeg;                                                                                                //!
   TBranch *b_nDtRechits;                                                                                            //!
   TBranch *b_nRpc;                                                                                                  //!
   TBranch *b_cscRechitsPhi;                                                                                         //!
   TBranch *b_cscRechitsEta;                                                                                         //!
   TBranch *b_cscRechitsX;                                                                                           //!
   TBranch *b_cscRechitsY;                                                                                           //!
   TBranch *b_cscRechitsZ;                                                                                           //!
   TBranch *b_cscRechitsTpeak;                                                                                       //!
   TBranch *b_cscRechitsTwire;                                                                                       //!
   TBranch *b_dtRechitCorrectX;                                                                                      //!
   TBranch *b_dtRechitCorrectY;                                                                                      //!
   TBranch *b_dtRechitCorrectZ;                                                                                      //!
   TBranch *b_dtRechitCorrectEta;                                                                                    //!
   TBranch *b_dtRechitCorrectPhi;                                                                                    //!
   TBranch *b_cscSegPhi;                                                                                             //!
   TBranch *b_cscSegEta;                                                                                             //!
   TBranch *b_dtSegPhi;                                                                                              //!
   TBranch *b_dtSegEta;                                                                                              //!
   TBranch *b_rpcPhi;                                                                                                //!
   TBranch *b_rpcEta;                                                                                                //!
   TBranch *b_rpcX;                                                                                                  //!
   TBranch *b_rpcY;                                                                                                  //!
   TBranch *b_rpcZ;                                                                                                  //!
   TBranch *b_rpcT;                                                                                                  //!
   TBranch *b_rpcTError;                                                                                             //!
   TBranch *b_cscRechitsChamber;                                                                                     //!
   TBranch *b_cscRechitsStation;                                                                                     //!
   TBranch *b_cscRechitsDetId;                                                                                       //!
   TBranch *b_dtRechitStation;                                                                                       //!
   TBranch *b_dtRechitWheel;                                                                                         //!
   TBranch *b_dtRechitSuperLayer;                                                                                    //!
   TBranch *b_cscSegChamber;                                                                                         //!
   TBranch *b_cscSegStation;                                                                                         //!
   TBranch *b_cscSegNRecHits;                                                                                        //!
   TBranch *b_dtSegStation;                                                                                          //!
   TBranch *b_dtSegWheel;                                                                                            //!
   TBranch *b_rpcBx;                                                                                                 //!
   TBranch *b_rpcRegion;                                                                                             //!
   TBranch *b_rpcRing;                                                                                               //!
   TBranch *b_rpcSector;                                                                                             //!
   TBranch *b_rpcStation;                                                                                            //!
   TBranch *b_rpcLayer;                                                                                              //!

   merged_event(TTree *tree = 0);
   virtual ~merged_event();
   virtual Int_t Cut(Long64_t entry);
   virtual Int_t GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void Init(TTree *tree);
   virtual void Loop();
   virtual Bool_t Notify();
   virtual void Show(Long64_t entry = -1);
};

#endif

#ifdef merged_event_cxx
merged_event::merged_event(TTree *tree) : fChain(0)
{
   // if parameter tree is not specified (or zero), connect the file
   // used to generate this class and read the Tree.
   if (tree == 0)
   {
      TFile *f = (TFile *)gROOT->GetListOfFiles()->FindObject("04bdd0bc-916b-4779-9bbd-b74e001e5423.root");
      if (!f || !f->IsOpen())
      {
         f = new TFile("04bdd0bc-916b-4779-9bbd-b74e001e5423.root");
      }
      f->GetObject("Events", tree);
   }
   Init(tree);
}

merged_event::~merged_event()
{
   if (!fChain)
      return;
   delete fChain->GetCurrentFile();
}

Int_t merged_event::GetEntry(Long64_t entry)
{
   // Read contents of entry.
   if (!fChain)
      return 0;
   return fChain->GetEntry(entry);
}
Long64_t merged_event::LoadTree(Long64_t entry)
{
   // Set the environment to read one entry
   if (!fChain)
      return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0)
      return centry;
   if (fChain->GetTreeNumber() != fCurrent)
   {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void merged_event::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree)
      return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);
   fChain->SetBranchAddress("nsoftActivityVH", &nsoftActivityVH, &b_nsoftActivityVH);

   fChain->SetBranchAddress("nGenJetAK8", &nGenJetAK8, &b_nGenJetAK8);
   fChain->SetBranchAddress("GenJetAK8_eta", GenJetAK8_eta, &b_GenJetAK8_eta);
   fChain->SetBranchAddress("GenJetAK8_mass", GenJetAK8_mass, &b_GenJetAK8_mass);
   fChain->SetBranchAddress("GenJetAK8_phi", GenJetAK8_phi, &b_GenJetAK8_phi);
   fChain->SetBranchAddress("GenJetAK8_pt", GenJetAK8_pt, &b_GenJetAK8_pt);
   fChain->SetBranchAddress("nGenJet", &nGenJet, &b_nGenJet);
   fChain->SetBranchAddress("GenJet_eta", GenJet_eta, &b_GenJet_eta);
   fChain->SetBranchAddress("GenJet_mass", GenJet_mass, &b_GenJet_mass);
   fChain->SetBranchAddress("GenJet_phi", GenJet_phi, &b_GenJet_phi);
   fChain->SetBranchAddress("GenJet_pt", GenJet_pt, &b_GenJet_pt);
   fChain->SetBranchAddress("nGenPart", &nGenPart, &b_nGenPart);
   fChain->SetBranchAddress("GenPart_genPartIdxMother", GenPart_genPartIdxMother, &b_GenPart_genPartIdxMother);
   fChain->SetBranchAddress("GenPart_statusFlags", GenPart_statusFlags, &b_GenPart_statusFlags);
   fChain->SetBranchAddress("GenPart_pdgId", GenPart_pdgId, &b_GenPart_pdgId);
   fChain->SetBranchAddress("GenPart_status", GenPart_status, &b_GenPart_status);
   fChain->SetBranchAddress("GenPart_eta", GenPart_eta, &b_GenPart_eta);
   fChain->SetBranchAddress("GenPart_mass", GenPart_mass, &b_GenPart_mass);
   fChain->SetBranchAddress("GenPart_phi", GenPart_phi, &b_GenPart_phi);
   fChain->SetBranchAddress("GenPart_pt", GenPart_pt, &b_GenPart_pt);
   fChain->SetBranchAddress("nGenProton", &nGenProton, &b_nGenProton);
   fChain->SetBranchAddress("GenProton_isPU", GenProton_isPU, &b_GenProton_isPU);
   fChain->SetBranchAddress("GenProton_px", GenProton_px, &b_GenProton_px);
   fChain->SetBranchAddress("GenProton_py", GenProton_py, &b_GenProton_py);
   fChain->SetBranchAddress("GenProton_pz", GenProton_pz, &b_GenProton_pz);
   fChain->SetBranchAddress("GenProton_vz", GenProton_vz, &b_GenProton_vz);
   fChain->SetBranchAddress("nSubGenJetAK8", &nSubGenJetAK8, &b_nSubGenJetAK8);
   fChain->SetBranchAddress("SubGenJetAK8_eta", SubGenJetAK8_eta, &b_SubGenJetAK8_eta);
   fChain->SetBranchAddress("SubGenJetAK8_mass", SubGenJetAK8_mass, &b_SubGenJetAK8_mass);
   fChain->SetBranchAddress("SubGenJetAK8_phi", SubGenJetAK8_phi, &b_SubGenJetAK8_phi);
   fChain->SetBranchAddress("SubGenJetAK8_pt", SubGenJetAK8_pt, &b_SubGenJetAK8_pt);
   fChain->SetBranchAddress("Generator_id1", &Generator_id1, &b_Generator_id1);
   fChain->SetBranchAddress("Generator_id2", &Generator_id2, &b_Generator_id2);
   fChain->SetBranchAddress("Generator_binvar", &Generator_binvar, &b_Generator_binvar);
   fChain->SetBranchAddress("Generator_scalePDF", &Generator_scalePDF, &b_Generator_scalePDF);
   fChain->SetBranchAddress("Generator_weight", &Generator_weight, &b_Generator_weight);
   fChain->SetBranchAddress("Generator_x1", &Generator_x1, &b_Generator_x1);
   fChain->SetBranchAddress("Generator_x2", &Generator_x2, &b_Generator_x2);
   fChain->SetBranchAddress("Generator_xpdf1", &Generator_xpdf1, &b_Generator_xpdf1);
   fChain->SetBranchAddress("Generator_xpdf2", &Generator_xpdf2, &b_Generator_xpdf2);
   fChain->SetBranchAddress("GenVtx_x", &GenVtx_x, &b_GenVtx_x);
   fChain->SetBranchAddress("GenVtx_y", &GenVtx_y, &b_GenVtx_y);
   fChain->SetBranchAddress("GenVtx_z", &GenVtx_z, &b_GenVtx_z);
   fChain->SetBranchAddress("nGenVisTau", &nGenVisTau, &b_nGenVisTau);
   fChain->SetBranchAddress("GenVisTau_status", GenVisTau_status, &b_GenVisTau_status);
   fChain->SetBranchAddress("GenVisTau_charge", GenVisTau_charge, &b_GenVisTau_charge);
   fChain->SetBranchAddress("GenVisTau_genPartIdxMother", GenVisTau_genPartIdxMother, &b_GenVisTau_genPartIdxMother);
   fChain->SetBranchAddress("GenVisTau_eta", GenVisTau_eta, &b_GenVisTau_eta);
   fChain->SetBranchAddress("GenVisTau_mass", GenVisTau_mass, &b_GenVisTau_mass);
   fChain->SetBranchAddress("GenVisTau_phi", GenVisTau_phi, &b_GenVisTau_phi);
   fChain->SetBranchAddress("GenVisTau_pt", GenVisTau_pt, &b_GenVisTau_pt);
   fChain->SetBranchAddress("genWeight", &genWeight, &b_genWeight);
   fChain->SetBranchAddress("GenMET_phi", &GenMET_phi, &b_GenMET_phi);
   fChain->SetBranchAddress("GenMET_pt", &GenMET_pt, &b_GenMET_pt);
   fChain->SetBranchAddress("nGenDressedLepton", &nGenDressedLepton, &b_nGenDressedLepton);
   fChain->SetBranchAddress("GenDressedLepton_hasTauAnc", GenDressedLepton_hasTauAnc, &b_GenDressedLepton_hasTauAnc);
   fChain->SetBranchAddress("GenDressedLepton_pdgId", GenDressedLepton_pdgId, &b_GenDressedLepton_pdgId);
   fChain->SetBranchAddress("GenDressedLepton_eta", GenDressedLepton_eta, &b_GenDressedLepton_eta);
   fChain->SetBranchAddress("GenDressedLepton_mass", GenDressedLepton_mass, &b_GenDressedLepton_mass);
   fChain->SetBranchAddress("GenDressedLepton_phi", GenDressedLepton_phi, &b_GenDressedLepton_phi);
   fChain->SetBranchAddress("GenDressedLepton_pt", GenDressedLepton_pt, &b_GenDressedLepton_pt);
   fChain->SetBranchAddress("nGenIsolatedPhoton", &nGenIsolatedPhoton, &b_nGenIsolatedPhoton);
   fChain->SetBranchAddress("GenIsolatedPhoton_eta", GenIsolatedPhoton_eta, &b_GenIsolatedPhoton_eta);
   fChain->SetBranchAddress("GenIsolatedPhoton_mass", GenIsolatedPhoton_mass, &b_GenIsolatedPhoton_mass);
   fChain->SetBranchAddress("GenIsolatedPhoton_phi", GenIsolatedPhoton_phi, &b_GenIsolatedPhoton_phi);
   fChain->SetBranchAddress("GenIsolatedPhoton_pt", GenIsolatedPhoton_pt, &b_GenIsolatedPhoton_pt);
   fChain->SetBranchAddress("genTtbarId", &genTtbarId, &b_genTtbarId);
   fChain->SetBranchAddress("boostedTau_genPartFlav", boostedTau_genPartFlav, &b_boostedTau_genPartFlav);
   fChain->SetBranchAddress("boostedTau_genPartIdx", boostedTau_genPartIdx, &b_boostedTau_genPartIdx);
   fChain->SetBranchAddress("Electron_genPartFlav", Electron_genPartFlav, &b_Electron_genPartFlav);
   fChain->SetBranchAddress("Electron_genPartIdx", Electron_genPartIdx, &b_Electron_genPartIdx);
   fChain->SetBranchAddress("FatJet_genJetAK8Idx", FatJet_genJetAK8Idx, &b_FatJet_genJetAK8Idx);
   fChain->SetBranchAddress("GenJetAK8_hadronFlavour", GenJetAK8_hadronFlavour, &b_GenJetAK8_hadronFlavour);
   fChain->SetBranchAddress("GenJetAK8_partonFlavour", GenJetAK8_partonFlavour, &b_GenJetAK8_partonFlavour);
   fChain->SetBranchAddress("GenJet_hadronFlavour", GenJet_hadronFlavour, &b_GenJet_hadronFlavour);
   fChain->SetBranchAddress("GenJet_partonFlavour", GenJet_partonFlavour, &b_GenJet_partonFlavour);
   fChain->SetBranchAddress("GenVtx_t0", &GenVtx_t0, &b_GenVtx_t0);
   fChain->SetBranchAddress("Jet_genJetIdx", Jet_genJetIdx, &b_Jet_genJetIdx);
   fChain->SetBranchAddress("LowPtElectron_genPartFlav", LowPtElectron_genPartFlav, &b_LowPtElectron_genPartFlav);
   fChain->SetBranchAddress("LowPtElectron_genPartIdx", LowPtElectron_genPartIdx, &b_LowPtElectron_genPartIdx);
   fChain->SetBranchAddress("Muon_genPartFlav", Muon_genPartFlav, &b_Muon_genPartFlav);
   fChain->SetBranchAddress("Muon_genPartIdx", Muon_genPartIdx, &b_Muon_genPartIdx);
   fChain->SetBranchAddress("Photon_genPartFlav", Photon_genPartFlav, &b_Photon_genPartFlav);
   fChain->SetBranchAddress("Photon_genPartIdx", Photon_genPartIdx, &b_Photon_genPartIdx);
   fChain->SetBranchAddress("MET_fiducialGenPhi", &MET_fiducialGenPhi, &b_MET_fiducialGenPhi);
   fChain->SetBranchAddress("MET_fiducialGenPt", &MET_fiducialGenPt, &b_MET_fiducialGenPt);
   fChain->SetBranchAddress("Tau_genPartFlav", Tau_genPartFlav, &b_Tau_genPartFlav);
   fChain->SetBranchAddress("Tau_genPartIdx", Tau_genPartIdx, &b_Tau_genPartIdx);

   fChain->SetBranchAddress("gLLP_eta", gLLP_eta, &b_gLLP_eta);
   fChain->SetBranchAddress("gLLP_phi", gLLP_phi, &b_gLLP_phi);
   fChain->SetBranchAddress("gLLP_csc", gLLP_csc, &b_gLLP_csc);
   fChain->SetBranchAddress("gLLP_dt", gLLP_dt, &b_gLLP_dt);
   fChain->SetBranchAddress("gLLP_beta", gLLP_beta, &b_gLLP_beta);
   fChain->SetBranchAddress("gLLP_e", gLLP_e, &b_gLLP_e);
   fChain->SetBranchAddress("gLLP_pt", gLLP_pt, &b_gLLP_pt);
   fChain->SetBranchAddress("gLLP_decay_vertex_x", gLLP_decay_vertex_x, &b_gLLP_decay_vertex_x);
   fChain->SetBranchAddress("gLLP_decay_vertex_y", gLLP_decay_vertex_y, &b_gLLP_decay_vertex_y);
   fChain->SetBranchAddress("gLLP_decay_vertex_z", gLLP_decay_vertex_z, &b_gLLP_decay_vertex_z);
   fChain->SetBranchAddress("gLLP_decay_vertex_r", gLLP_decay_vertex_r, &b_gLLP_decay_vertex_r);

   fChain->SetBranchAddress("run", &run, &b_run);
   fChain->SetBranchAddress("luminosityBlock", &luminosityBlock, &b_luminosityBlock);
   fChain->SetBranchAddress("event", &event, &b_event);
   fChain->SetBranchAddress("bunchCrossing", &bunchCrossing, &b_bunchCrossing);
   fChain->SetBranchAddress("BeamSpot_type", &BeamSpot_type, &b_BeamSpot_type);
   fChain->SetBranchAddress("BeamSpot_sigmaZ", &BeamSpot_sigmaZ, &b_BeamSpot_sigmaZ);
   fChain->SetBranchAddress("BeamSpot_sigmaZError", &BeamSpot_sigmaZError, &b_BeamSpot_sigmaZError);
   fChain->SetBranchAddress("BeamSpot_z", &BeamSpot_z, &b_BeamSpot_z);
   fChain->SetBranchAddress("BeamSpot_zError", &BeamSpot_zError, &b_BeamSpot_zError);
   fChain->SetBranchAddress("nboostedTau", &nboostedTau, &b_nboostedTau);
   fChain->SetBranchAddress("boostedTau_idAntiEle2018", boostedTau_idAntiEle2018, &b_boostedTau_idAntiEle2018);
   fChain->SetBranchAddress("boostedTau_idAntiMu", boostedTau_idAntiMu, &b_boostedTau_idAntiMu);
   fChain->SetBranchAddress("boostedTau_idMVAnewDM2017v2", boostedTau_idMVAnewDM2017v2, &b_boostedTau_idMVAnewDM2017v2);
   fChain->SetBranchAddress("boostedTau_idMVAoldDM2017v2", boostedTau_idMVAoldDM2017v2, &b_boostedTau_idMVAoldDM2017v2);
   fChain->SetBranchAddress("boostedTau_jetIdx", boostedTau_jetIdx, &b_boostedTau_jetIdx);
   fChain->SetBranchAddress("boostedTau_rawAntiEleCat2018", boostedTau_rawAntiEleCat2018, &b_boostedTau_rawAntiEleCat2018);
   fChain->SetBranchAddress("boostedTau_charge", boostedTau_charge, &b_boostedTau_charge);
   fChain->SetBranchAddress("boostedTau_decayMode", boostedTau_decayMode, &b_boostedTau_decayMode);
   fChain->SetBranchAddress("boostedTau_chargedIso", boostedTau_chargedIso, &b_boostedTau_chargedIso);
   fChain->SetBranchAddress("boostedTau_eta", boostedTau_eta, &b_boostedTau_eta);
   fChain->SetBranchAddress("boostedTau_leadTkDeltaEta", boostedTau_leadTkDeltaEta, &b_boostedTau_leadTkDeltaEta);
   fChain->SetBranchAddress("boostedTau_leadTkDeltaPhi", boostedTau_leadTkDeltaPhi, &b_boostedTau_leadTkDeltaPhi);
   fChain->SetBranchAddress("boostedTau_leadTkPtOverTauPt", boostedTau_leadTkPtOverTauPt, &b_boostedTau_leadTkPtOverTauPt);
   fChain->SetBranchAddress("boostedTau_mass", boostedTau_mass, &b_boostedTau_mass);
   fChain->SetBranchAddress("boostedTau_neutralIso", boostedTau_neutralIso, &b_boostedTau_neutralIso);
   fChain->SetBranchAddress("boostedTau_phi", boostedTau_phi, &b_boostedTau_phi);
   fChain->SetBranchAddress("boostedTau_photonsOutsideSignalCone", boostedTau_photonsOutsideSignalCone, &b_boostedTau_photonsOutsideSignalCone);
   fChain->SetBranchAddress("boostedTau_pt", boostedTau_pt, &b_boostedTau_pt);
   fChain->SetBranchAddress("boostedTau_puCorr", boostedTau_puCorr, &b_boostedTau_puCorr);
   fChain->SetBranchAddress("boostedTau_rawAntiEle2018", boostedTau_rawAntiEle2018, &b_boostedTau_rawAntiEle2018);
   fChain->SetBranchAddress("boostedTau_rawIso", boostedTau_rawIso, &b_boostedTau_rawIso);
   fChain->SetBranchAddress("boostedTau_rawIsodR03", boostedTau_rawIsodR03, &b_boostedTau_rawIsodR03);
   fChain->SetBranchAddress("boostedTau_rawMVAnewDM2017v2", boostedTau_rawMVAnewDM2017v2, &b_boostedTau_rawMVAnewDM2017v2);
   fChain->SetBranchAddress("boostedTau_rawMVAoldDM2017v2", boostedTau_rawMVAoldDM2017v2, &b_boostedTau_rawMVAoldDM2017v2);
   fChain->SetBranchAddress("CaloMET_phi", &CaloMET_phi, &b_CaloMET_phi);
   fChain->SetBranchAddress("CaloMET_pt", &CaloMET_pt, &b_CaloMET_pt);
   fChain->SetBranchAddress("CaloMET_sumEt", &CaloMET_sumEt, &b_CaloMET_sumEt);
   fChain->SetBranchAddress("ChsMET_phi", &ChsMET_phi, &b_ChsMET_phi);
   fChain->SetBranchAddress("ChsMET_pt", &ChsMET_pt, &b_ChsMET_pt);
   fChain->SetBranchAddress("ChsMET_sumEt", &ChsMET_sumEt, &b_ChsMET_sumEt);
   fChain->SetBranchAddress("nCorrT1METJet", &nCorrT1METJet, &b_nCorrT1METJet);
   fChain->SetBranchAddress("CorrT1METJet_area", CorrT1METJet_area, &b_CorrT1METJet_area);
   fChain->SetBranchAddress("CorrT1METJet_eta", CorrT1METJet_eta, &b_CorrT1METJet_eta);
   fChain->SetBranchAddress("CorrT1METJet_muonSubtrFactor", CorrT1METJet_muonSubtrFactor, &b_CorrT1METJet_muonSubtrFactor);
   fChain->SetBranchAddress("CorrT1METJet_phi", CorrT1METJet_phi, &b_CorrT1METJet_phi);
   fChain->SetBranchAddress("CorrT1METJet_rawPt", CorrT1METJet_rawPt, &b_CorrT1METJet_rawPt);
   fChain->SetBranchAddress("DeepMETResolutionTune_phi", &DeepMETResolutionTune_phi, &b_DeepMETResolutionTune_phi);
   fChain->SetBranchAddress("DeepMETResolutionTune_pt", &DeepMETResolutionTune_pt, &b_DeepMETResolutionTune_pt);
   fChain->SetBranchAddress("DeepMETResponseTune_phi", &DeepMETResponseTune_phi, &b_DeepMETResponseTune_phi);
   fChain->SetBranchAddress("DeepMETResponseTune_pt", &DeepMETResponseTune_pt, &b_DeepMETResponseTune_pt);
   fChain->SetBranchAddress("nElectron", &nElectron, &b_nElectron);
   fChain->SetBranchAddress("Electron_seediEtaOriX", Electron_seediEtaOriX, &b_Electron_seediEtaOriX);
   fChain->SetBranchAddress("Electron_convVeto", Electron_convVeto, &b_Electron_convVeto);
   fChain->SetBranchAddress("Electron_cutBased", Electron_cutBased, &b_Electron_cutBased);
   fChain->SetBranchAddress("Electron_cutBased_HEEP", Electron_cutBased_HEEP, &b_Electron_cutBased_HEEP);
   fChain->SetBranchAddress("Electron_isPFcand", Electron_isPFcand, &b_Electron_isPFcand);
   fChain->SetBranchAddress("Electron_jetNDauCharged", Electron_jetNDauCharged, &b_Electron_jetNDauCharged);
   fChain->SetBranchAddress("Electron_lostHits", Electron_lostHits, &b_Electron_lostHits);
   fChain->SetBranchAddress("Electron_mvaIso_WP80", Electron_mvaIso_WP80, &b_Electron_mvaIso_WP80);
   fChain->SetBranchAddress("Electron_mvaIso_WP90", Electron_mvaIso_WP90, &b_Electron_mvaIso_WP90);
   fChain->SetBranchAddress("Electron_mvaNoIso_WP80", Electron_mvaNoIso_WP80, &b_Electron_mvaNoIso_WP80);
   fChain->SetBranchAddress("Electron_mvaNoIso_WP90", Electron_mvaNoIso_WP90, &b_Electron_mvaNoIso_WP90);
   fChain->SetBranchAddress("Electron_seedGain", Electron_seedGain, &b_Electron_seedGain);
   fChain->SetBranchAddress("Electron_tightCharge", Electron_tightCharge, &b_Electron_tightCharge);
   fChain->SetBranchAddress("Electron_jetIdx", Electron_jetIdx, &b_Electron_jetIdx);
   fChain->SetBranchAddress("Electron_photonIdx", Electron_photonIdx, &b_Electron_photonIdx);
   fChain->SetBranchAddress("Electron_svIdx", Electron_svIdx, &b_Electron_svIdx);
   fChain->SetBranchAddress("Electron_fsrPhotonIdx", Electron_fsrPhotonIdx, &b_Electron_fsrPhotonIdx);
   fChain->SetBranchAddress("Electron_charge", Electron_charge, &b_Electron_charge);
   fChain->SetBranchAddress("Electron_pdgId", Electron_pdgId, &b_Electron_pdgId);
   fChain->SetBranchAddress("Electron_seediPhiOriY", Electron_seediPhiOriY, &b_Electron_seediPhiOriY);
   fChain->SetBranchAddress("Electron_vidNestedWPBitmap", Electron_vidNestedWPBitmap, &b_Electron_vidNestedWPBitmap);
   fChain->SetBranchAddress("Electron_vidNestedWPBitmapHEEP", Electron_vidNestedWPBitmapHEEP, &b_Electron_vidNestedWPBitmapHEEP);
   fChain->SetBranchAddress("Electron_deltaEtaSC", Electron_deltaEtaSC, &b_Electron_deltaEtaSC);
   fChain->SetBranchAddress("Electron_dr03EcalRecHitSumEt", Electron_dr03EcalRecHitSumEt, &b_Electron_dr03EcalRecHitSumEt);
   fChain->SetBranchAddress("Electron_dr03HcalDepth1TowerSumEt", Electron_dr03HcalDepth1TowerSumEt, &b_Electron_dr03HcalDepth1TowerSumEt);
   fChain->SetBranchAddress("Electron_dr03TkSumPt", Electron_dr03TkSumPt, &b_Electron_dr03TkSumPt);
   fChain->SetBranchAddress("Electron_dr03TkSumPtHEEP", Electron_dr03TkSumPtHEEP, &b_Electron_dr03TkSumPtHEEP);
   fChain->SetBranchAddress("Electron_dxy", Electron_dxy, &b_Electron_dxy);
   fChain->SetBranchAddress("Electron_dxyErr", Electron_dxyErr, &b_Electron_dxyErr);
   fChain->SetBranchAddress("Electron_dz", Electron_dz, &b_Electron_dz);
   fChain->SetBranchAddress("Electron_dzErr", Electron_dzErr, &b_Electron_dzErr);
   fChain->SetBranchAddress("Electron_eInvMinusPInv", Electron_eInvMinusPInv, &b_Electron_eInvMinusPInv);
   fChain->SetBranchAddress("Electron_energyErr", Electron_energyErr, &b_Electron_energyErr);
   fChain->SetBranchAddress("Electron_eta", Electron_eta, &b_Electron_eta);
   fChain->SetBranchAddress("Electron_hoe", Electron_hoe, &b_Electron_hoe);
   fChain->SetBranchAddress("Electron_ip3d", Electron_ip3d, &b_Electron_ip3d);
   fChain->SetBranchAddress("Electron_jetPtRelv2", Electron_jetPtRelv2, &b_Electron_jetPtRelv2);
   fChain->SetBranchAddress("Electron_jetRelIso", Electron_jetRelIso, &b_Electron_jetRelIso);
   fChain->SetBranchAddress("Electron_mass", Electron_mass, &b_Electron_mass);
   fChain->SetBranchAddress("Electron_miniPFRelIso_all", Electron_miniPFRelIso_all, &b_Electron_miniPFRelIso_all);
   fChain->SetBranchAddress("Electron_miniPFRelIso_chg", Electron_miniPFRelIso_chg, &b_Electron_miniPFRelIso_chg);
   fChain->SetBranchAddress("Electron_mvaHZZIso", Electron_mvaHZZIso, &b_Electron_mvaHZZIso);
   fChain->SetBranchAddress("Electron_mvaIso", Electron_mvaIso, &b_Electron_mvaIso);
   fChain->SetBranchAddress("Electron_mvaNoIso", Electron_mvaNoIso, &b_Electron_mvaNoIso);
   fChain->SetBranchAddress("Electron_pfRelIso03_all", Electron_pfRelIso03_all, &b_Electron_pfRelIso03_all);
   fChain->SetBranchAddress("Electron_pfRelIso03_chg", Electron_pfRelIso03_chg, &b_Electron_pfRelIso03_chg);
   fChain->SetBranchAddress("Electron_phi", Electron_phi, &b_Electron_phi);
   fChain->SetBranchAddress("Electron_pt", Electron_pt, &b_Electron_pt);
   fChain->SetBranchAddress("Electron_r9", Electron_r9, &b_Electron_r9);
   fChain->SetBranchAddress("Electron_scEtOverPt", Electron_scEtOverPt, &b_Electron_scEtOverPt);
   fChain->SetBranchAddress("Electron_sieie", Electron_sieie, &b_Electron_sieie);
   fChain->SetBranchAddress("Electron_sip3d", Electron_sip3d, &b_Electron_sip3d);
   fChain->SetBranchAddress("Electron_mvaTTH", Electron_mvaTTH, &b_Electron_mvaTTH);
   fChain->SetBranchAddress("nFatJet", &nFatJet, &b_nFatJet);
   fChain->SetBranchAddress("FatJet_jetId", FatJet_jetId, &b_FatJet_jetId);
   fChain->SetBranchAddress("FatJet_nConstituents", FatJet_nConstituents, &b_FatJet_nConstituents);
   fChain->SetBranchAddress("FatJet_subJetIdx1", FatJet_subJetIdx1, &b_FatJet_subJetIdx1);
   fChain->SetBranchAddress("FatJet_subJetIdx2", FatJet_subJetIdx2, &b_FatJet_subJetIdx2);
   fChain->SetBranchAddress("FatJet_electronIdx3SJ", FatJet_electronIdx3SJ, &b_FatJet_electronIdx3SJ);
   fChain->SetBranchAddress("FatJet_muonIdx3SJ", FatJet_muonIdx3SJ, &b_FatJet_muonIdx3SJ);
   fChain->SetBranchAddress("FatJet_area", FatJet_area, &b_FatJet_area);
   fChain->SetBranchAddress("FatJet_btagDDBvLV2", FatJet_btagDDBvLV2, &b_FatJet_btagDDBvLV2);
   fChain->SetBranchAddress("FatJet_btagDDCvBV2", FatJet_btagDDCvBV2, &b_FatJet_btagDDCvBV2);
   fChain->SetBranchAddress("FatJet_btagDDCvLV2", FatJet_btagDDCvLV2, &b_FatJet_btagDDCvLV2);
   fChain->SetBranchAddress("FatJet_btagDeepB", FatJet_btagDeepB, &b_FatJet_btagDeepB);
   fChain->SetBranchAddress("FatJet_btagHbb", FatJet_btagHbb, &b_FatJet_btagHbb);
   fChain->SetBranchAddress("FatJet_eta", FatJet_eta, &b_FatJet_eta);
   fChain->SetBranchAddress("FatJet_mass", FatJet_mass, &b_FatJet_mass);
   fChain->SetBranchAddress("FatJet_msoftdrop", FatJet_msoftdrop, &b_FatJet_msoftdrop);
   fChain->SetBranchAddress("FatJet_n2b1", FatJet_n2b1, &b_FatJet_n2b1);
   fChain->SetBranchAddress("FatJet_n3b1", FatJet_n3b1, &b_FatJet_n3b1);
   fChain->SetBranchAddress("FatJet_particleNetWithMass_H4qvsQCD", FatJet_particleNetWithMass_H4qvsQCD, &b_FatJet_particleNetWithMass_H4qvsQCD);
   fChain->SetBranchAddress("FatJet_particleNetWithMass_HbbvsQCD", FatJet_particleNetWithMass_HbbvsQCD, &b_FatJet_particleNetWithMass_HbbvsQCD);
   fChain->SetBranchAddress("FatJet_particleNetWithMass_HccvsQCD", FatJet_particleNetWithMass_HccvsQCD, &b_FatJet_particleNetWithMass_HccvsQCD);
   fChain->SetBranchAddress("FatJet_particleNetWithMass_QCD", FatJet_particleNetWithMass_QCD, &b_FatJet_particleNetWithMass_QCD);
   fChain->SetBranchAddress("FatJet_particleNetWithMass_TvsQCD", FatJet_particleNetWithMass_TvsQCD, &b_FatJet_particleNetWithMass_TvsQCD);
   fChain->SetBranchAddress("FatJet_particleNetWithMass_WvsQCD", FatJet_particleNetWithMass_WvsQCD, &b_FatJet_particleNetWithMass_WvsQCD);
   fChain->SetBranchAddress("FatJet_particleNetWithMass_ZvsQCD", FatJet_particleNetWithMass_ZvsQCD, &b_FatJet_particleNetWithMass_ZvsQCD);
   fChain->SetBranchAddress("FatJet_particleNet_QCD", FatJet_particleNet_QCD, &b_FatJet_particleNet_QCD);
   fChain->SetBranchAddress("FatJet_particleNet_QCD0HF", FatJet_particleNet_QCD0HF, &b_FatJet_particleNet_QCD0HF);
   fChain->SetBranchAddress("FatJet_particleNet_QCD1HF", FatJet_particleNet_QCD1HF, &b_FatJet_particleNet_QCD1HF);
   fChain->SetBranchAddress("FatJet_particleNet_QCD2HF", FatJet_particleNet_QCD2HF, &b_FatJet_particleNet_QCD2HF);
   fChain->SetBranchAddress("FatJet_particleNet_XbbVsQCD", FatJet_particleNet_XbbVsQCD, &b_FatJet_particleNet_XbbVsQCD);
   fChain->SetBranchAddress("FatJet_particleNet_XccVsQCD", FatJet_particleNet_XccVsQCD, &b_FatJet_particleNet_XccVsQCD);
   fChain->SetBranchAddress("FatJet_particleNet_XggVsQCD", FatJet_particleNet_XggVsQCD, &b_FatJet_particleNet_XggVsQCD);
   fChain->SetBranchAddress("FatJet_particleNet_XqqVsQCD", FatJet_particleNet_XqqVsQCD, &b_FatJet_particleNet_XqqVsQCD);
   fChain->SetBranchAddress("FatJet_particleNet_XteVsQCD", FatJet_particleNet_XteVsQCD, &b_FatJet_particleNet_XteVsQCD);
   fChain->SetBranchAddress("FatJet_particleNet_XtmVsQCD", FatJet_particleNet_XtmVsQCD, &b_FatJet_particleNet_XtmVsQCD);
   fChain->SetBranchAddress("FatJet_particleNet_XttVsQCD", FatJet_particleNet_XttVsQCD, &b_FatJet_particleNet_XttVsQCD);
   fChain->SetBranchAddress("FatJet_particleNet_massCorr", FatJet_particleNet_massCorr, &b_FatJet_particleNet_massCorr);
   fChain->SetBranchAddress("FatJet_phi", FatJet_phi, &b_FatJet_phi);
   fChain->SetBranchAddress("FatJet_pt", FatJet_pt, &b_FatJet_pt);
   fChain->SetBranchAddress("FatJet_rawFactor", FatJet_rawFactor, &b_FatJet_rawFactor);
   fChain->SetBranchAddress("FatJet_tau1", FatJet_tau1, &b_FatJet_tau1);
   fChain->SetBranchAddress("FatJet_tau2", FatJet_tau2, &b_FatJet_tau2);
   fChain->SetBranchAddress("FatJet_tau3", FatJet_tau3, &b_FatJet_tau3);
   fChain->SetBranchAddress("FatJet_tau4", FatJet_tau4, &b_FatJet_tau4);
   fChain->SetBranchAddress("FatJet_lsf3", FatJet_lsf3, &b_FatJet_lsf3);
   fChain->SetBranchAddress("nFsrPhoton", &nFsrPhoton, &b_nFsrPhoton);
   fChain->SetBranchAddress("FsrPhoton_electronIdx", FsrPhoton_electronIdx, &b_FsrPhoton_electronIdx);
   fChain->SetBranchAddress("FsrPhoton_muonIdx", FsrPhoton_muonIdx, &b_FsrPhoton_muonIdx);
   fChain->SetBranchAddress("FsrPhoton_dROverEt2", FsrPhoton_dROverEt2, &b_FsrPhoton_dROverEt2);
   fChain->SetBranchAddress("FsrPhoton_eta", FsrPhoton_eta, &b_FsrPhoton_eta);
   fChain->SetBranchAddress("FsrPhoton_phi", FsrPhoton_phi, &b_FsrPhoton_phi);
   fChain->SetBranchAddress("FsrPhoton_pt", FsrPhoton_pt, &b_FsrPhoton_pt);
   fChain->SetBranchAddress("FsrPhoton_relIso03", FsrPhoton_relIso03, &b_FsrPhoton_relIso03);
   fChain->SetBranchAddress("nIsoTrack", &nIsoTrack, &b_nIsoTrack);
   fChain->SetBranchAddress("IsoTrack_isHighPurityTrack", IsoTrack_isHighPurityTrack, &b_IsoTrack_isHighPurityTrack);
   fChain->SetBranchAddress("IsoTrack_isPFcand", IsoTrack_isPFcand, &b_IsoTrack_isPFcand);
   fChain->SetBranchAddress("IsoTrack_isFromLostTrack", IsoTrack_isFromLostTrack, &b_IsoTrack_isFromLostTrack);
   fChain->SetBranchAddress("IsoTrack_charge", IsoTrack_charge, &b_IsoTrack_charge);
   fChain->SetBranchAddress("IsoTrack_fromPV", IsoTrack_fromPV, &b_IsoTrack_fromPV);
   fChain->SetBranchAddress("IsoTrack_pdgId", IsoTrack_pdgId, &b_IsoTrack_pdgId);
   fChain->SetBranchAddress("IsoTrack_dxy", IsoTrack_dxy, &b_IsoTrack_dxy);
   fChain->SetBranchAddress("IsoTrack_dz", IsoTrack_dz, &b_IsoTrack_dz);
   fChain->SetBranchAddress("IsoTrack_eta", IsoTrack_eta, &b_IsoTrack_eta);
   fChain->SetBranchAddress("IsoTrack_pfRelIso03_all", IsoTrack_pfRelIso03_all, &b_IsoTrack_pfRelIso03_all);
   fChain->SetBranchAddress("IsoTrack_pfRelIso03_chg", IsoTrack_pfRelIso03_chg, &b_IsoTrack_pfRelIso03_chg);
   fChain->SetBranchAddress("IsoTrack_phi", IsoTrack_phi, &b_IsoTrack_phi);
   fChain->SetBranchAddress("IsoTrack_pt", IsoTrack_pt, &b_IsoTrack_pt);
   fChain->SetBranchAddress("IsoTrack_miniPFRelIso_all", IsoTrack_miniPFRelIso_all, &b_IsoTrack_miniPFRelIso_all);
   fChain->SetBranchAddress("IsoTrack_miniPFRelIso_chg", IsoTrack_miniPFRelIso_chg, &b_IsoTrack_miniPFRelIso_chg);
   fChain->SetBranchAddress("nJet", &nJet, &b_nJet);
   fChain->SetBranchAddress("Jet_jetId", Jet_jetId, &b_Jet_jetId);
   fChain->SetBranchAddress("Jet_nConstituents", Jet_nConstituents, &b_Jet_nConstituents);
   fChain->SetBranchAddress("Jet_nElectrons", Jet_nElectrons, &b_Jet_nElectrons);
   fChain->SetBranchAddress("Jet_nMuons", Jet_nMuons, &b_Jet_nMuons);
   fChain->SetBranchAddress("Jet_nSVs", Jet_nSVs, &b_Jet_nSVs);
   fChain->SetBranchAddress("Jet_electronIdx1", Jet_electronIdx1, &b_Jet_electronIdx1);
   fChain->SetBranchAddress("Jet_electronIdx2", Jet_electronIdx2, &b_Jet_electronIdx2);
   fChain->SetBranchAddress("Jet_muonIdx1", Jet_muonIdx1, &b_Jet_muonIdx1);
   fChain->SetBranchAddress("Jet_muonIdx2", Jet_muonIdx2, &b_Jet_muonIdx2);
   fChain->SetBranchAddress("Jet_svIdx1", Jet_svIdx1, &b_Jet_svIdx1);
   fChain->SetBranchAddress("Jet_svIdx2", Jet_svIdx2, &b_Jet_svIdx2);
   fChain->SetBranchAddress("Jet_hfadjacentEtaStripsSize", Jet_hfadjacentEtaStripsSize, &b_Jet_hfadjacentEtaStripsSize);
   fChain->SetBranchAddress("Jet_hfcentralEtaStripSize", Jet_hfcentralEtaStripSize, &b_Jet_hfcentralEtaStripSize);
   fChain->SetBranchAddress("Jet_PNetRegPtRawCorr", Jet_PNetRegPtRawCorr, &b_Jet_PNetRegPtRawCorr);
   fChain->SetBranchAddress("Jet_PNetRegPtRawCorrNeutrino", Jet_PNetRegPtRawCorrNeutrino, &b_Jet_PNetRegPtRawCorrNeutrino);
   fChain->SetBranchAddress("Jet_PNetRegPtRawRes", Jet_PNetRegPtRawRes, &b_Jet_PNetRegPtRawRes);
   fChain->SetBranchAddress("Jet_area", Jet_area, &b_Jet_area);
   fChain->SetBranchAddress("Jet_btagDeepFlavB", Jet_btagDeepFlavB, &b_Jet_btagDeepFlavB);
   fChain->SetBranchAddress("Jet_btagDeepFlavCvB", Jet_btagDeepFlavCvB, &b_Jet_btagDeepFlavCvB);
   fChain->SetBranchAddress("Jet_btagDeepFlavCvL", Jet_btagDeepFlavCvL, &b_Jet_btagDeepFlavCvL);
   fChain->SetBranchAddress("Jet_btagDeepFlavQG", Jet_btagDeepFlavQG, &b_Jet_btagDeepFlavQG);
   fChain->SetBranchAddress("Jet_btagPNetB", Jet_btagPNetB, &b_Jet_btagPNetB);
   fChain->SetBranchAddress("Jet_btagPNetCvB", Jet_btagPNetCvB, &b_Jet_btagPNetCvB);
   fChain->SetBranchAddress("Jet_btagPNetCvL", Jet_btagPNetCvL, &b_Jet_btagPNetCvL);
   fChain->SetBranchAddress("Jet_btagPNetQvG", Jet_btagPNetQvG, &b_Jet_btagPNetQvG);
   fChain->SetBranchAddress("Jet_btagPNetTauVJet", Jet_btagPNetTauVJet, &b_Jet_btagPNetTauVJet);
   fChain->SetBranchAddress("Jet_btagRobustParTAK4B", Jet_btagRobustParTAK4B, &b_Jet_btagRobustParTAK4B);
   fChain->SetBranchAddress("Jet_btagRobustParTAK4CvB", Jet_btagRobustParTAK4CvB, &b_Jet_btagRobustParTAK4CvB);
   fChain->SetBranchAddress("Jet_btagRobustParTAK4CvL", Jet_btagRobustParTAK4CvL, &b_Jet_btagRobustParTAK4CvL);
   fChain->SetBranchAddress("Jet_btagRobustParTAK4QG", Jet_btagRobustParTAK4QG, &b_Jet_btagRobustParTAK4QG);
   fChain->SetBranchAddress("Jet_chEmEF", Jet_chEmEF, &b_Jet_chEmEF);
   fChain->SetBranchAddress("Jet_chHEF", Jet_chHEF, &b_Jet_chHEF);
   fChain->SetBranchAddress("Jet_eta", Jet_eta, &b_Jet_eta);
   fChain->SetBranchAddress("Jet_hfsigmaEtaEta", Jet_hfsigmaEtaEta, &b_Jet_hfsigmaEtaEta);
   fChain->SetBranchAddress("Jet_hfsigmaPhiPhi", Jet_hfsigmaPhiPhi, &b_Jet_hfsigmaPhiPhi);
   fChain->SetBranchAddress("Jet_mass", Jet_mass, &b_Jet_mass);
   fChain->SetBranchAddress("Jet_muEF", Jet_muEF, &b_Jet_muEF);
   fChain->SetBranchAddress("Jet_muonSubtrFactor", Jet_muonSubtrFactor, &b_Jet_muonSubtrFactor);
   fChain->SetBranchAddress("Jet_neEmEF", Jet_neEmEF, &b_Jet_neEmEF);
   fChain->SetBranchAddress("Jet_neHEF", Jet_neHEF, &b_Jet_neHEF);
   fChain->SetBranchAddress("Jet_phi", Jet_phi, &b_Jet_phi);
   fChain->SetBranchAddress("Jet_pt", Jet_pt, &b_Jet_pt);
   fChain->SetBranchAddress("Jet_rawFactor", Jet_rawFactor, &b_Jet_rawFactor);
   fChain->SetBranchAddress("nLowPtElectron", &nLowPtElectron, &b_nLowPtElectron);
   fChain->SetBranchAddress("LowPtElectron_convVeto", LowPtElectron_convVeto, &b_LowPtElectron_convVeto);
   fChain->SetBranchAddress("LowPtElectron_convWP", LowPtElectron_convWP, &b_LowPtElectron_convWP);
   fChain->SetBranchAddress("LowPtElectron_lostHits", LowPtElectron_lostHits, &b_LowPtElectron_lostHits);
   fChain->SetBranchAddress("LowPtElectron_electronIdx", LowPtElectron_electronIdx, &b_LowPtElectron_electronIdx);
   fChain->SetBranchAddress("LowPtElectron_photonIdx", LowPtElectron_photonIdx, &b_LowPtElectron_photonIdx);
   fChain->SetBranchAddress("LowPtElectron_charge", LowPtElectron_charge, &b_LowPtElectron_charge);
   fChain->SetBranchAddress("LowPtElectron_pdgId", LowPtElectron_pdgId, &b_LowPtElectron_pdgId);
   fChain->SetBranchAddress("LowPtElectron_ID", LowPtElectron_ID, &b_LowPtElectron_ID);
   fChain->SetBranchAddress("LowPtElectron_convVtxRadius", LowPtElectron_convVtxRadius, &b_LowPtElectron_convVtxRadius);
   fChain->SetBranchAddress("LowPtElectron_deltaEtaSC", LowPtElectron_deltaEtaSC, &b_LowPtElectron_deltaEtaSC);
   fChain->SetBranchAddress("LowPtElectron_dxy", LowPtElectron_dxy, &b_LowPtElectron_dxy);
   fChain->SetBranchAddress("LowPtElectron_dxyErr", LowPtElectron_dxyErr, &b_LowPtElectron_dxyErr);
   fChain->SetBranchAddress("LowPtElectron_dz", LowPtElectron_dz, &b_LowPtElectron_dz);
   fChain->SetBranchAddress("LowPtElectron_dzErr", LowPtElectron_dzErr, &b_LowPtElectron_dzErr);
   fChain->SetBranchAddress("LowPtElectron_eInvMinusPInv", LowPtElectron_eInvMinusPInv, &b_LowPtElectron_eInvMinusPInv);
   fChain->SetBranchAddress("LowPtElectron_energyErr", LowPtElectron_energyErr, &b_LowPtElectron_energyErr);
   fChain->SetBranchAddress("LowPtElectron_eta", LowPtElectron_eta, &b_LowPtElectron_eta);
   fChain->SetBranchAddress("LowPtElectron_hoe", LowPtElectron_hoe, &b_LowPtElectron_hoe);
   fChain->SetBranchAddress("LowPtElectron_mass", LowPtElectron_mass, &b_LowPtElectron_mass);
   fChain->SetBranchAddress("LowPtElectron_miniPFRelIso_all", LowPtElectron_miniPFRelIso_all, &b_LowPtElectron_miniPFRelIso_all);
   fChain->SetBranchAddress("LowPtElectron_miniPFRelIso_chg", LowPtElectron_miniPFRelIso_chg, &b_LowPtElectron_miniPFRelIso_chg);
   fChain->SetBranchAddress("LowPtElectron_phi", LowPtElectron_phi, &b_LowPtElectron_phi);
   fChain->SetBranchAddress("LowPtElectron_pt", LowPtElectron_pt, &b_LowPtElectron_pt);
   fChain->SetBranchAddress("LowPtElectron_ptbiased", LowPtElectron_ptbiased, &b_LowPtElectron_ptbiased);
   fChain->SetBranchAddress("LowPtElectron_r9", LowPtElectron_r9, &b_LowPtElectron_r9);
   fChain->SetBranchAddress("LowPtElectron_scEtOverPt", LowPtElectron_scEtOverPt, &b_LowPtElectron_scEtOverPt);
   fChain->SetBranchAddress("LowPtElectron_sieie", LowPtElectron_sieie, &b_LowPtElectron_sieie);
   fChain->SetBranchAddress("LowPtElectron_unbiased", LowPtElectron_unbiased, &b_LowPtElectron_unbiased);
   fChain->SetBranchAddress("MET_MetUnclustEnUpDeltaX", &MET_MetUnclustEnUpDeltaX, &b_MET_MetUnclustEnUpDeltaX);
   fChain->SetBranchAddress("MET_MetUnclustEnUpDeltaY", &MET_MetUnclustEnUpDeltaY, &b_MET_MetUnclustEnUpDeltaY);
   fChain->SetBranchAddress("MET_covXX", &MET_covXX, &b_MET_covXX);
   fChain->SetBranchAddress("MET_covXY", &MET_covXY, &b_MET_covXY);
   fChain->SetBranchAddress("MET_covYY", &MET_covYY, &b_MET_covYY);
   fChain->SetBranchAddress("MET_phi", &MET_phi, &b_MET_phi);
   fChain->SetBranchAddress("MET_pt", &MET_pt, &b_MET_pt);
   fChain->SetBranchAddress("MET_significance", &MET_significance, &b_MET_significance);
   fChain->SetBranchAddress("MET_sumEt", &MET_sumEt, &b_MET_sumEt);
   fChain->SetBranchAddress("MET_sumPtUnclustered", &MET_sumPtUnclustered, &b_MET_sumPtUnclustered);
   fChain->SetBranchAddress("nProton_multiRP", &nProton_multiRP, &b_nProton_multiRP);
   fChain->SetBranchAddress("Proton_multiRP_arm", &Proton_multiRP_arm, &b_Proton_multiRP_arm);
   fChain->SetBranchAddress("Proton_multiRP_t", &Proton_multiRP_t, &b_Proton_multiRP_t);
   fChain->SetBranchAddress("Proton_multiRP_thetaX", &Proton_multiRP_thetaX, &b_Proton_multiRP_thetaX);
   fChain->SetBranchAddress("Proton_multiRP_thetaY", &Proton_multiRP_thetaY, &b_Proton_multiRP_thetaY);
   fChain->SetBranchAddress("Proton_multiRP_time", &Proton_multiRP_time, &b_Proton_multiRP_time);
   fChain->SetBranchAddress("Proton_multiRP_timeUnc", &Proton_multiRP_timeUnc, &b_Proton_multiRP_timeUnc);
   fChain->SetBranchAddress("Proton_multiRP_xi", &Proton_multiRP_xi, &b_Proton_multiRP_xi);
   fChain->SetBranchAddress("nMuon", &nMuon, &b_nMuon);
   fChain->SetBranchAddress("Muon_highPtId", Muon_highPtId, &b_Muon_highPtId);
   fChain->SetBranchAddress("Muon_highPurity", Muon_highPurity, &b_Muon_highPurity);
   fChain->SetBranchAddress("Muon_inTimeMuon", Muon_inTimeMuon, &b_Muon_inTimeMuon);
   fChain->SetBranchAddress("Muon_isGlobal", Muon_isGlobal, &b_Muon_isGlobal);
   fChain->SetBranchAddress("Muon_isPFcand", Muon_isPFcand, &b_Muon_isPFcand);
   fChain->SetBranchAddress("Muon_isStandalone", Muon_isStandalone, &b_Muon_isStandalone);
   fChain->SetBranchAddress("Muon_isTracker", Muon_isTracker, &b_Muon_isTracker);
   fChain->SetBranchAddress("Muon_jetNDauCharged", Muon_jetNDauCharged, &b_Muon_jetNDauCharged);
   fChain->SetBranchAddress("Muon_looseId", Muon_looseId, &b_Muon_looseId);
   fChain->SetBranchAddress("Muon_mediumId", Muon_mediumId, &b_Muon_mediumId);
   fChain->SetBranchAddress("Muon_mediumPromptId", Muon_mediumPromptId, &b_Muon_mediumPromptId);
   fChain->SetBranchAddress("Muon_miniIsoId", Muon_miniIsoId, &b_Muon_miniIsoId);
   fChain->SetBranchAddress("Muon_multiIsoId", Muon_multiIsoId, &b_Muon_multiIsoId);
   fChain->SetBranchAddress("Muon_mvaMuID_WP", Muon_mvaMuID_WP, &b_Muon_mvaMuID_WP);
   fChain->SetBranchAddress("Muon_nStations", Muon_nStations, &b_Muon_nStations);
   fChain->SetBranchAddress("Muon_nTrackerLayers", Muon_nTrackerLayers, &b_Muon_nTrackerLayers);
   fChain->SetBranchAddress("Muon_pfIsoId", Muon_pfIsoId, &b_Muon_pfIsoId);
   fChain->SetBranchAddress("Muon_puppiIsoId", Muon_puppiIsoId, &b_Muon_puppiIsoId);
   fChain->SetBranchAddress("Muon_softId", Muon_softId, &b_Muon_softId);
   fChain->SetBranchAddress("Muon_softMvaId", Muon_softMvaId, &b_Muon_softMvaId);
   fChain->SetBranchAddress("Muon_tightCharge", Muon_tightCharge, &b_Muon_tightCharge);
   fChain->SetBranchAddress("Muon_tightId", Muon_tightId, &b_Muon_tightId);
   fChain->SetBranchAddress("Muon_tkIsoId", Muon_tkIsoId, &b_Muon_tkIsoId);
   fChain->SetBranchAddress("Muon_triggerIdLoose", Muon_triggerIdLoose, &b_Muon_triggerIdLoose);
   fChain->SetBranchAddress("Muon_jetIdx", Muon_jetIdx, &b_Muon_jetIdx);
   fChain->SetBranchAddress("Muon_svIdx", Muon_svIdx, &b_Muon_svIdx);
   fChain->SetBranchAddress("Muon_fsrPhotonIdx", Muon_fsrPhotonIdx, &b_Muon_fsrPhotonIdx);
   fChain->SetBranchAddress("Muon_charge", Muon_charge, &b_Muon_charge);
   fChain->SetBranchAddress("Muon_pdgId", Muon_pdgId, &b_Muon_pdgId);
   fChain->SetBranchAddress("Muon_dxy", Muon_dxy, &b_Muon_dxy);
   fChain->SetBranchAddress("Muon_dxyErr", Muon_dxyErr, &b_Muon_dxyErr);
   fChain->SetBranchAddress("Muon_dxybs", Muon_dxybs, &b_Muon_dxybs);
   fChain->SetBranchAddress("Muon_dz", Muon_dz, &b_Muon_dz);
   fChain->SetBranchAddress("Muon_dzErr", Muon_dzErr, &b_Muon_dzErr);
   fChain->SetBranchAddress("Muon_eta", Muon_eta, &b_Muon_eta);
   fChain->SetBranchAddress("Muon_ip3d", Muon_ip3d, &b_Muon_ip3d);
   fChain->SetBranchAddress("Muon_jetPtRelv2", Muon_jetPtRelv2, &b_Muon_jetPtRelv2);
   fChain->SetBranchAddress("Muon_jetRelIso", Muon_jetRelIso, &b_Muon_jetRelIso);
   fChain->SetBranchAddress("Muon_mass", Muon_mass, &b_Muon_mass);
   fChain->SetBranchAddress("Muon_miniPFRelIso_all", Muon_miniPFRelIso_all, &b_Muon_miniPFRelIso_all);
   fChain->SetBranchAddress("Muon_miniPFRelIso_chg", Muon_miniPFRelIso_chg, &b_Muon_miniPFRelIso_chg);
   fChain->SetBranchAddress("Muon_mvaMuID", Muon_mvaMuID, &b_Muon_mvaMuID);
   fChain->SetBranchAddress("Muon_pfRelIso03_all", Muon_pfRelIso03_all, &b_Muon_pfRelIso03_all);
   fChain->SetBranchAddress("Muon_pfRelIso03_chg", Muon_pfRelIso03_chg, &b_Muon_pfRelIso03_chg);
   fChain->SetBranchAddress("Muon_pfRelIso04_all", Muon_pfRelIso04_all, &b_Muon_pfRelIso04_all);
   fChain->SetBranchAddress("Muon_phi", Muon_phi, &b_Muon_phi);
   fChain->SetBranchAddress("Muon_pt", Muon_pt, &b_Muon_pt);
   fChain->SetBranchAddress("Muon_ptErr", Muon_ptErr, &b_Muon_ptErr);
   fChain->SetBranchAddress("Muon_segmentComp", Muon_segmentComp, &b_Muon_segmentComp);
   fChain->SetBranchAddress("Muon_sip3d", Muon_sip3d, &b_Muon_sip3d);
   fChain->SetBranchAddress("Muon_softMva", Muon_softMva, &b_Muon_softMva);
   fChain->SetBranchAddress("Muon_tkRelIso", Muon_tkRelIso, &b_Muon_tkRelIso);
   fChain->SetBranchAddress("Muon_tunepRelPt", Muon_tunepRelPt, &b_Muon_tunepRelPt);
   fChain->SetBranchAddress("Muon_bsConstrainedChi2", Muon_bsConstrainedChi2, &b_Muon_bsConstrainedChi2);
   fChain->SetBranchAddress("Muon_bsConstrainedPt", Muon_bsConstrainedPt, &b_Muon_bsConstrainedPt);
   fChain->SetBranchAddress("Muon_bsConstrainedPtErr", Muon_bsConstrainedPtErr, &b_Muon_bsConstrainedPtErr);
   fChain->SetBranchAddress("Muon_mvaLowPt", Muon_mvaLowPt, &b_Muon_mvaLowPt);
   fChain->SetBranchAddress("Muon_mvaTTH", Muon_mvaTTH, &b_Muon_mvaTTH);
   fChain->SetBranchAddress("PFMET_covXX", &PFMET_covXX, &b_PFMET_covXX);
   fChain->SetBranchAddress("PFMET_covXY", &PFMET_covXY, &b_PFMET_covXY);
   fChain->SetBranchAddress("PFMET_covYY", &PFMET_covYY, &b_PFMET_covYY);
   fChain->SetBranchAddress("PFMET_phi", &PFMET_phi, &b_PFMET_phi);
   fChain->SetBranchAddress("PFMET_phiUnclusteredDown", &PFMET_phiUnclusteredDown, &b_PFMET_phiUnclusteredDown);
   fChain->SetBranchAddress("PFMET_phiUnclusteredUp", &PFMET_phiUnclusteredUp, &b_PFMET_phiUnclusteredUp);
   // fChain->SetBranchAddress("PFMET_pt", &PFMET_pt, &b_PFMET_pt);
   fChain->SetBranchAddress("PFMET_ptUnclusteredDown", &PFMET_ptUnclusteredDown, &b_PFMET_ptUnclusteredDown);
   fChain->SetBranchAddress("PFMET_ptUnclusteredUp", &PFMET_ptUnclusteredUp, &b_PFMET_ptUnclusteredUp);
   fChain->SetBranchAddress("PFMET_significance", &PFMET_significance, &b_PFMET_significance);
   fChain->SetBranchAddress("PFMET_sumEt", &PFMET_sumEt, &b_PFMET_sumEt);
   fChain->SetBranchAddress("PFMET_sumPtUnclustered", &PFMET_sumPtUnclustered, &b_PFMET_sumPtUnclustered);

   fChain->SetBranchAddress("nPhoton", &nPhoton, &b_nPhoton);
   fChain->SetBranchAddress("Photon_seediEtaOriX", Photon_seediEtaOriX, &b_Photon_seediEtaOriX);
   fChain->SetBranchAddress("Photon_cutBased", Photon_cutBased, &b_Photon_cutBased);
   fChain->SetBranchAddress("Photon_electronVeto", Photon_electronVeto, &b_Photon_electronVeto);
   fChain->SetBranchAddress("Photon_hasConversionTracks", Photon_hasConversionTracks, &b_Photon_hasConversionTracks);
   fChain->SetBranchAddress("Photon_isScEtaEB", Photon_isScEtaEB, &b_Photon_isScEtaEB);
   fChain->SetBranchAddress("Photon_isScEtaEE", Photon_isScEtaEE, &b_Photon_isScEtaEE);
   fChain->SetBranchAddress("Photon_mvaID_WP80", Photon_mvaID_WP80, &b_Photon_mvaID_WP80);
   fChain->SetBranchAddress("Photon_mvaID_WP90", Photon_mvaID_WP90, &b_Photon_mvaID_WP90);
   fChain->SetBranchAddress("Photon_pixelSeed", Photon_pixelSeed, &b_Photon_pixelSeed);
   fChain->SetBranchAddress("Photon_seedGain", Photon_seedGain, &b_Photon_seedGain);
   fChain->SetBranchAddress("Photon_electronIdx", Photon_electronIdx, &b_Photon_electronIdx);
   fChain->SetBranchAddress("Photon_jetIdx", Photon_jetIdx, &b_Photon_jetIdx);
   fChain->SetBranchAddress("Photon_seediPhiOriY", Photon_seediPhiOriY, &b_Photon_seediPhiOriY);
   fChain->SetBranchAddress("Photon_vidNestedWPBitmap", Photon_vidNestedWPBitmap, &b_Photon_vidNestedWPBitmap);
   fChain->SetBranchAddress("Photon_ecalPFClusterIso", Photon_ecalPFClusterIso, &b_Photon_ecalPFClusterIso);
   fChain->SetBranchAddress("Photon_energyErr", Photon_energyErr, &b_Photon_energyErr);
   fChain->SetBranchAddress("Photon_energyRaw", Photon_energyRaw, &b_Photon_energyRaw);
   fChain->SetBranchAddress("Photon_esEffSigmaRR", Photon_esEffSigmaRR, &b_Photon_esEffSigmaRR);
   fChain->SetBranchAddress("Photon_esEnergyOverRawE", Photon_esEnergyOverRawE, &b_Photon_esEnergyOverRawE);
   fChain->SetBranchAddress("Photon_eta", Photon_eta, &b_Photon_eta);
   fChain->SetBranchAddress("Photon_etaWidth", Photon_etaWidth, &b_Photon_etaWidth);
   fChain->SetBranchAddress("Photon_haloTaggerMVAVal", Photon_haloTaggerMVAVal, &b_Photon_haloTaggerMVAVal);
   fChain->SetBranchAddress("Photon_hcalPFClusterIso", Photon_hcalPFClusterIso, &b_Photon_hcalPFClusterIso);
   fChain->SetBranchAddress("Photon_hoe", Photon_hoe, &b_Photon_hoe);
   fChain->SetBranchAddress("Photon_hoe_PUcorr", Photon_hoe_PUcorr, &b_Photon_hoe_PUcorr);
   fChain->SetBranchAddress("Photon_mvaID", Photon_mvaID, &b_Photon_mvaID);
   fChain->SetBranchAddress("Photon_pfChargedIso", Photon_pfChargedIso, &b_Photon_pfChargedIso);
   fChain->SetBranchAddress("Photon_pfChargedIsoPFPV", Photon_pfChargedIsoPFPV, &b_Photon_pfChargedIsoPFPV);
   fChain->SetBranchAddress("Photon_pfChargedIsoWorstVtx", Photon_pfChargedIsoWorstVtx, &b_Photon_pfChargedIsoWorstVtx);
   fChain->SetBranchAddress("Photon_pfPhoIso03", Photon_pfPhoIso03, &b_Photon_pfPhoIso03);
   fChain->SetBranchAddress("Photon_pfRelIso03_all_quadratic", Photon_pfRelIso03_all_quadratic, &b_Photon_pfRelIso03_all_quadratic);
   fChain->SetBranchAddress("Photon_pfRelIso03_chg_quadratic", Photon_pfRelIso03_chg_quadratic, &b_Photon_pfRelIso03_chg_quadratic);
   fChain->SetBranchAddress("Photon_phi", Photon_phi, &b_Photon_phi);
   fChain->SetBranchAddress("Photon_phiWidth", Photon_phiWidth, &b_Photon_phiWidth);
   fChain->SetBranchAddress("Photon_pt", Photon_pt, &b_Photon_pt);
   fChain->SetBranchAddress("Photon_r9", Photon_r9, &b_Photon_r9);
   fChain->SetBranchAddress("Photon_s4", Photon_s4, &b_Photon_s4);
   fChain->SetBranchAddress("Photon_sieie", Photon_sieie, &b_Photon_sieie);
   fChain->SetBranchAddress("Photon_sieip", Photon_sieip, &b_Photon_sieip);
   fChain->SetBranchAddress("Photon_sipip", Photon_sipip, &b_Photon_sipip);
   fChain->SetBranchAddress("Photon_trkSumPtHollowConeDR03", Photon_trkSumPtHollowConeDR03, &b_Photon_trkSumPtHollowConeDR03);
   fChain->SetBranchAddress("Photon_trkSumPtSolidConeDR04", Photon_trkSumPtSolidConeDR04, &b_Photon_trkSumPtSolidConeDR04);
   fChain->SetBranchAddress("Photon_x_calo", Photon_x_calo, &b_Photon_x_calo);
   fChain->SetBranchAddress("Photon_y_calo", Photon_y_calo, &b_Photon_y_calo);
   fChain->SetBranchAddress("Photon_z_calo", Photon_z_calo, &b_Photon_z_calo);
   fChain->SetBranchAddress("nPPSLocalTrack", &nPPSLocalTrack, &b_nPPSLocalTrack);
   fChain->SetBranchAddress("PPSLocalTrack_multiRPProtonIdx", &PPSLocalTrack_multiRPProtonIdx, &b_PPSLocalTrack_multiRPProtonIdx);
   fChain->SetBranchAddress("PPSLocalTrack_singleRPProtonIdx", &PPSLocalTrack_singleRPProtonIdx, &b_PPSLocalTrack_singleRPProtonIdx);
   fChain->SetBranchAddress("PPSLocalTrack_decRPId", &PPSLocalTrack_decRPId, &b_PPSLocalTrack_decRPId);
   fChain->SetBranchAddress("PPSLocalTrack_rpType", &PPSLocalTrack_rpType, &b_PPSLocalTrack_rpType);
   fChain->SetBranchAddress("PPSLocalTrack_x", &PPSLocalTrack_x, &b_PPSLocalTrack_x);
   fChain->SetBranchAddress("PPSLocalTrack_y", &PPSLocalTrack_y, &b_PPSLocalTrack_y);
   fChain->SetBranchAddress("PPSLocalTrack_time", &PPSLocalTrack_time, &b_PPSLocalTrack_time);
   fChain->SetBranchAddress("PPSLocalTrack_timeUnc", &PPSLocalTrack_timeUnc, &b_PPSLocalTrack_timeUnc);
    fChain->SetBranchAddress("Pileup_nPU", &Pileup_nPU, &b_Pileup_nPU);
   fChain->SetBranchAddress("Pileup_sumEOOT", &Pileup_sumEOOT, &b_Pileup_sumEOOT);
   fChain->SetBranchAddress("Pileup_sumLOOT", &Pileup_sumLOOT, &b_Pileup_sumLOOT);
   fChain->SetBranchAddress("Pileup_nTrueInt", &Pileup_nTrueInt, &b_Pileup_nTrueInt);
   fChain->SetBranchAddress("Pileup_pudensity", &Pileup_pudensity, &b_Pileup_pudensity);
   fChain->SetBranchAddress("Pileup_gpudensity", &Pileup_gpudensity, &b_Pileup_gpudensity);
   fChain->SetBranchAddress("PuppiMET_phi", &PuppiMET_phi, &b_PuppiMET_phi);
   fChain->SetBranchAddress("PuppiMET_phiJERDown", &PuppiMET_phiJERDown, &b_PuppiMET_phiJERDown);
   fChain->SetBranchAddress("PuppiMET_phiJERUp", &PuppiMET_phiJERUp, &b_PuppiMET_phiJERUp);
   fChain->SetBranchAddress("PuppiMET_phiJESDown", &PuppiMET_phiJESDown, &b_PuppiMET_phiJESDown);
   fChain->SetBranchAddress("PuppiMET_phiJESUp", &PuppiMET_phiJESUp, &b_PuppiMET_phiJESUp);
   fChain->SetBranchAddress("PuppiMET_phiUnclusteredDown", &PuppiMET_phiUnclusteredDown, &b_PuppiMET_phiUnclusteredDown);
   fChain->SetBranchAddress("PuppiMET_phiUnclusteredUp", &PuppiMET_phiUnclusteredUp, &b_PuppiMET_phiUnclusteredUp);
   fChain->SetBranchAddress("PuppiMET_pt", &PuppiMET_pt, &b_PuppiMET_pt);
   fChain->SetBranchAddress("PuppiMET_ptJERDown", &PuppiMET_ptJERDown, &b_PuppiMET_ptJERDown);
   fChain->SetBranchAddress("PuppiMET_ptJERUp", &PuppiMET_ptJERUp, &b_PuppiMET_ptJERUp);
   fChain->SetBranchAddress("PuppiMET_ptJESDown", &PuppiMET_ptJESDown, &b_PuppiMET_ptJESDown);
   fChain->SetBranchAddress("PuppiMET_ptJESUp", &PuppiMET_ptJESUp, &b_PuppiMET_ptJESUp);
   fChain->SetBranchAddress("PuppiMET_ptUnclusteredDown", &PuppiMET_ptUnclusteredDown, &b_PuppiMET_ptUnclusteredDown);
   fChain->SetBranchAddress("PuppiMET_ptUnclusteredUp", &PuppiMET_ptUnclusteredUp, &b_PuppiMET_ptUnclusteredUp);
   fChain->SetBranchAddress("PuppiMET_sumEt", &PuppiMET_sumEt, &b_PuppiMET_sumEt);
   fChain->SetBranchAddress("RawMET_phi", &RawMET_phi, &b_RawMET_phi);
   fChain->SetBranchAddress("RawMET_pt", &RawMET_pt, &b_RawMET_pt);
   fChain->SetBranchAddress("RawMET_sumEt", &RawMET_sumEt, &b_RawMET_sumEt);
   fChain->SetBranchAddress("RawPuppiMET_phi", &RawPuppiMET_phi, &b_RawPuppiMET_phi);
   fChain->SetBranchAddress("RawPuppiMET_pt", &RawPuppiMET_pt, &b_RawPuppiMET_pt);
   fChain->SetBranchAddress("RawPuppiMET_sumEt", &RawPuppiMET_sumEt, &b_RawPuppiMET_sumEt);
   fChain->SetBranchAddress("Rho_fixedGridRhoAll", &Rho_fixedGridRhoAll, &b_Rho_fixedGridRhoAll);
   fChain->SetBranchAddress("Rho_fixedGridRhoFastjetAll", &Rho_fixedGridRhoFastjetAll, &b_Rho_fixedGridRhoFastjetAll);
   fChain->SetBranchAddress("Rho_fixedGridRhoFastjetCentral", &Rho_fixedGridRhoFastjetCentral, &b_Rho_fixedGridRhoFastjetCentral);
   fChain->SetBranchAddress("Rho_fixedGridRhoFastjetCentralCalo", &Rho_fixedGridRhoFastjetCentralCalo, &b_Rho_fixedGridRhoFastjetCentralCalo);
   fChain->SetBranchAddress("Rho_fixedGridRhoFastjetCentralChargedPileUp", &Rho_fixedGridRhoFastjetCentralChargedPileUp, &b_Rho_fixedGridRhoFastjetCentralChargedPileUp);
   fChain->SetBranchAddress("Rho_fixedGridRhoFastjetCentralNeutral", &Rho_fixedGridRhoFastjetCentralNeutral, &b_Rho_fixedGridRhoFastjetCentralNeutral);
   fChain->SetBranchAddress("nSoftActivityJet", &nSoftActivityJet, &b_nSoftActivityJet);
   fChain->SetBranchAddress("SoftActivityJet_eta", SoftActivityJet_eta, &b_SoftActivityJet_eta);
   fChain->SetBranchAddress("SoftActivityJet_phi", SoftActivityJet_phi, &b_SoftActivityJet_phi);
   fChain->SetBranchAddress("SoftActivityJet_pt", SoftActivityJet_pt, &b_SoftActivityJet_pt);
   fChain->SetBranchAddress("SoftActivityJetNjets10", &SoftActivityJetNjets10, &b_SoftActivityJetNjets10);
   fChain->SetBranchAddress("SoftActivityJetNjets2", &SoftActivityJetNjets2, &b_SoftActivityJetNjets2);
   fChain->SetBranchAddress("SoftActivityJetNjets5", &SoftActivityJetNjets5, &b_SoftActivityJetNjets5);
   fChain->SetBranchAddress("SoftActivityJetHT", &SoftActivityJetHT, &b_SoftActivityJetHT);
   fChain->SetBranchAddress("SoftActivityJetHT10", &SoftActivityJetHT10, &b_SoftActivityJetHT10);
   fChain->SetBranchAddress("SoftActivityJetHT2", &SoftActivityJetHT2, &b_SoftActivityJetHT2);
   fChain->SetBranchAddress("SoftActivityJetHT5", &SoftActivityJetHT5, &b_SoftActivityJetHT5);
   fChain->SetBranchAddress("nProton_singleRP", &nProton_singleRP, &b_nProton_singleRP);
   fChain->SetBranchAddress("Proton_singleRP_decRPId", &Proton_singleRP_decRPId, &b_Proton_singleRP_decRPId);
   fChain->SetBranchAddress("Proton_singleRP_thetaY", &Proton_singleRP_thetaY, &b_Proton_singleRP_thetaY);
   fChain->SetBranchAddress("Proton_singleRP_xi", &Proton_singleRP_xi, &b_Proton_singleRP_xi);
   fChain->SetBranchAddress("nSubJet", &nSubJet, &b_nSubJet);
   fChain->SetBranchAddress("SubJet_btagDeepB", SubJet_btagDeepB, &b_SubJet_btagDeepB);
   fChain->SetBranchAddress("SubJet_eta", SubJet_eta, &b_SubJet_eta);
   fChain->SetBranchAddress("SubJet_mass", SubJet_mass, &b_SubJet_mass);
   fChain->SetBranchAddress("SubJet_n2b1", SubJet_n2b1, &b_SubJet_n2b1);
   fChain->SetBranchAddress("SubJet_n3b1", SubJet_n3b1, &b_SubJet_n3b1);
   fChain->SetBranchAddress("SubJet_phi", SubJet_phi, &b_SubJet_phi);
   fChain->SetBranchAddress("SubJet_pt", SubJet_pt, &b_SubJet_pt);
   fChain->SetBranchAddress("SubJet_rawFactor", SubJet_rawFactor, &b_SubJet_rawFactor);
   fChain->SetBranchAddress("SubJet_tau1", SubJet_tau1, &b_SubJet_tau1);
   fChain->SetBranchAddress("SubJet_tau2", SubJet_tau2, &b_SubJet_tau2);
   fChain->SetBranchAddress("SubJet_tau3", SubJet_tau3, &b_SubJet_tau3);
   fChain->SetBranchAddress("SubJet_tau4", SubJet_tau4, &b_SubJet_tau4);
   fChain->SetBranchAddress("nTau", &nTau, &b_nTau);
   fChain->SetBranchAddress("Tau_decayMode", Tau_decayMode, &b_Tau_decayMode);
   fChain->SetBranchAddress("Tau_idAntiEleDeadECal", Tau_idAntiEleDeadECal, &b_Tau_idAntiEleDeadECal);
   fChain->SetBranchAddress("Tau_idAntiMu", Tau_idAntiMu, &b_Tau_idAntiMu);
   fChain->SetBranchAddress("Tau_idDecayModeNewDMs", Tau_idDecayModeNewDMs, &b_Tau_idDecayModeNewDMs);
   fChain->SetBranchAddress("Tau_idDecayModeOldDMs", Tau_idDecayModeOldDMs, &b_Tau_idDecayModeOldDMs);
   fChain->SetBranchAddress("Tau_idDeepTau2017v2p1VSe", Tau_idDeepTau2017v2p1VSe, &b_Tau_idDeepTau2017v2p1VSe);
   fChain->SetBranchAddress("Tau_idDeepTau2017v2p1VSjet", Tau_idDeepTau2017v2p1VSjet, &b_Tau_idDeepTau2017v2p1VSjet);
   fChain->SetBranchAddress("Tau_idDeepTau2017v2p1VSmu", Tau_idDeepTau2017v2p1VSmu, &b_Tau_idDeepTau2017v2p1VSmu);
   fChain->SetBranchAddress("Tau_idDeepTau2018v2p5VSe", Tau_idDeepTau2018v2p5VSe, &b_Tau_idDeepTau2018v2p5VSe);
   fChain->SetBranchAddress("Tau_idDeepTau2018v2p5VSjet", Tau_idDeepTau2018v2p5VSjet, &b_Tau_idDeepTau2018v2p5VSjet);
   fChain->SetBranchAddress("Tau_idDeepTau2018v2p5VSmu", Tau_idDeepTau2018v2p5VSmu, &b_Tau_idDeepTau2018v2p5VSmu);
   fChain->SetBranchAddress("Tau_nSVs", Tau_nSVs, &b_Tau_nSVs);
   fChain->SetBranchAddress("Tau_charge", Tau_charge, &b_Tau_charge);
   fChain->SetBranchAddress("Tau_decayModePNet", Tau_decayModePNet, &b_Tau_decayModePNet);
   fChain->SetBranchAddress("Tau_eleIdx", Tau_eleIdx, &b_Tau_eleIdx);
   fChain->SetBranchAddress("Tau_jetIdx", Tau_jetIdx, &b_Tau_jetIdx);
   fChain->SetBranchAddress("Tau_muIdx", Tau_muIdx, &b_Tau_muIdx);
   fChain->SetBranchAddress("Tau_svIdx1", Tau_svIdx1, &b_Tau_svIdx1);
   fChain->SetBranchAddress("Tau_svIdx2", Tau_svIdx2, &b_Tau_svIdx2);
   fChain->SetBranchAddress("Tau_chargedIso", Tau_chargedIso, &b_Tau_chargedIso);
   fChain->SetBranchAddress("Tau_dxy", Tau_dxy, &b_Tau_dxy);
   fChain->SetBranchAddress("Tau_dz", Tau_dz, &b_Tau_dz);
   fChain->SetBranchAddress("Tau_eta", Tau_eta, &b_Tau_eta);
   fChain->SetBranchAddress("Tau_leadTkDeltaEta", Tau_leadTkDeltaEta, &b_Tau_leadTkDeltaEta);
   fChain->SetBranchAddress("Tau_leadTkDeltaPhi", Tau_leadTkDeltaPhi, &b_Tau_leadTkDeltaPhi);
   fChain->SetBranchAddress("Tau_leadTkPtOverTauPt", Tau_leadTkPtOverTauPt, &b_Tau_leadTkPtOverTauPt);
   fChain->SetBranchAddress("Tau_mass", Tau_mass, &b_Tau_mass);
   fChain->SetBranchAddress("Tau_neutralIso", Tau_neutralIso, &b_Tau_neutralIso);
   fChain->SetBranchAddress("Tau_phi", Tau_phi, &b_Tau_phi);
   fChain->SetBranchAddress("Tau_photonsOutsideSignalCone", Tau_photonsOutsideSignalCone, &b_Tau_photonsOutsideSignalCone);
   fChain->SetBranchAddress("Tau_probDM0PNet", Tau_probDM0PNet, &b_Tau_probDM0PNet);
   fChain->SetBranchAddress("Tau_probDM10PNet", Tau_probDM10PNet, &b_Tau_probDM10PNet);
   fChain->SetBranchAddress("Tau_probDM11PNet", Tau_probDM11PNet, &b_Tau_probDM11PNet);
   fChain->SetBranchAddress("Tau_probDM1PNet", Tau_probDM1PNet, &b_Tau_probDM1PNet);
   fChain->SetBranchAddress("Tau_probDM2PNet", Tau_probDM2PNet, &b_Tau_probDM2PNet);
   fChain->SetBranchAddress("Tau_pt", Tau_pt, &b_Tau_pt);
   fChain->SetBranchAddress("Tau_ptCorrPNet", Tau_ptCorrPNet, &b_Tau_ptCorrPNet);
   fChain->SetBranchAddress("Tau_puCorr", Tau_puCorr, &b_Tau_puCorr);
   fChain->SetBranchAddress("Tau_qConfPNet", Tau_qConfPNet, &b_Tau_qConfPNet);
   fChain->SetBranchAddress("Tau_rawDeepTau2017v2p1VSe", Tau_rawDeepTau2017v2p1VSe, &b_Tau_rawDeepTau2017v2p1VSe);
   fChain->SetBranchAddress("Tau_rawDeepTau2017v2p1VSjet", Tau_rawDeepTau2017v2p1VSjet, &b_Tau_rawDeepTau2017v2p1VSjet);
   fChain->SetBranchAddress("Tau_rawDeepTau2017v2p1VSmu", Tau_rawDeepTau2017v2p1VSmu, &b_Tau_rawDeepTau2017v2p1VSmu);
   fChain->SetBranchAddress("Tau_rawDeepTau2018v2p5VSe", Tau_rawDeepTau2018v2p5VSe, &b_Tau_rawDeepTau2018v2p5VSe);
   fChain->SetBranchAddress("Tau_rawDeepTau2018v2p5VSjet", Tau_rawDeepTau2018v2p5VSjet, &b_Tau_rawDeepTau2018v2p5VSjet);
   fChain->SetBranchAddress("Tau_rawDeepTau2018v2p5VSmu", Tau_rawDeepTau2018v2p5VSmu, &b_Tau_rawDeepTau2018v2p5VSmu);
   fChain->SetBranchAddress("Tau_rawIso", Tau_rawIso, &b_Tau_rawIso);
   fChain->SetBranchAddress("Tau_rawIsodR03", Tau_rawIsodR03, &b_Tau_rawIsodR03);
   fChain->SetBranchAddress("Tau_rawPNetVSe", Tau_rawPNetVSe, &b_Tau_rawPNetVSe);
   fChain->SetBranchAddress("Tau_rawPNetVSjet", Tau_rawPNetVSjet, &b_Tau_rawPNetVSjet);
   fChain->SetBranchAddress("Tau_rawPNetVSmu", Tau_rawPNetVSmu, &b_Tau_rawPNetVSmu);
   fChain->SetBranchAddress("TkMET_phi", &TkMET_phi, &b_TkMET_phi);
   fChain->SetBranchAddress("TkMET_pt", &TkMET_pt, &b_TkMET_pt);
   fChain->SetBranchAddress("TkMET_sumEt", &TkMET_sumEt, &b_TkMET_sumEt);
   fChain->SetBranchAddress("nTrigObj", &nTrigObj, &b_nTrigObj);
   fChain->SetBranchAddress("TrigObj_l1charge", TrigObj_l1charge, &b_TrigObj_l1charge);
   fChain->SetBranchAddress("TrigObj_id", TrigObj_id, &b_TrigObj_id);
   fChain->SetBranchAddress("TrigObj_l1iso", TrigObj_l1iso, &b_TrigObj_l1iso);
   fChain->SetBranchAddress("TrigObj_filterBits", TrigObj_filterBits, &b_TrigObj_filterBits);
   fChain->SetBranchAddress("TrigObj_pt", TrigObj_pt, &b_TrigObj_pt);
   fChain->SetBranchAddress("TrigObj_eta", TrigObj_eta, &b_TrigObj_eta);
   fChain->SetBranchAddress("TrigObj_phi", TrigObj_phi, &b_TrigObj_phi);
   fChain->SetBranchAddress("TrigObj_l1pt", TrigObj_l1pt, &b_TrigObj_l1pt);
   fChain->SetBranchAddress("TrigObj_l1pt_2", TrigObj_l1pt_2, &b_TrigObj_l1pt_2);
   fChain->SetBranchAddress("TrigObj_l2pt", TrigObj_l2pt, &b_TrigObj_l2pt);
   fChain->SetBranchAddress("nOtherPV", &nOtherPV, &b_nOtherPV);
   fChain->SetBranchAddress("OtherPV_z", OtherPV_z, &b_OtherPV_z);
   fChain->SetBranchAddress("OtherPV_score", OtherPV_score, &b_OtherPV_score);
   fChain->SetBranchAddress("PV_npvs", &PV_npvs, &b_PV_npvs);
   fChain->SetBranchAddress("PV_npvsGood", &PV_npvsGood, &b_PV_npvsGood);
   fChain->SetBranchAddress("PV_ndof", &PV_ndof, &b_PV_ndof);
   fChain->SetBranchAddress("PV_x", &PV_x, &b_PV_x);
   fChain->SetBranchAddress("PV_y", &PV_y, &b_PV_y);
   fChain->SetBranchAddress("PV_z", &PV_z, &b_PV_z);
   fChain->SetBranchAddress("PV_chi2", &PV_chi2, &b_PV_chi2);
   fChain->SetBranchAddress("PV_score", &PV_score, &b_PV_score);
   fChain->SetBranchAddress("nSV", &nSV, &b_nSV);
   fChain->SetBranchAddress("SV_charge", SV_charge, &b_SV_charge);
   fChain->SetBranchAddress("SV_dlen", SV_dlen, &b_SV_dlen);
   fChain->SetBranchAddress("SV_dlenSig", SV_dlenSig, &b_SV_dlenSig);
   fChain->SetBranchAddress("SV_dxy", SV_dxy, &b_SV_dxy);
   fChain->SetBranchAddress("SV_dxySig", SV_dxySig, &b_SV_dxySig);
   fChain->SetBranchAddress("SV_pAngle", SV_pAngle, &b_SV_pAngle);
   fChain->SetBranchAddress("SV_ntracks", SV_ntracks, &b_SV_ntracks);
   fChain->SetBranchAddress("SV_chi2", SV_chi2, &b_SV_chi2);
   fChain->SetBranchAddress("SV_eta", SV_eta, &b_SV_eta);
   fChain->SetBranchAddress("SV_mass", SV_mass, &b_SV_mass);
   fChain->SetBranchAddress("SV_ndof", SV_ndof, &b_SV_ndof);
   fChain->SetBranchAddress("SV_phi", SV_phi, &b_SV_phi);
   fChain->SetBranchAddress("SV_pt", SV_pt, &b_SV_pt);
   fChain->SetBranchAddress("SV_x", SV_x, &b_SV_x);
   fChain->SetBranchAddress("SV_y", SV_y, &b_SV_y);
   fChain->SetBranchAddress("SV_z", SV_z, &b_SV_z);
   fChain->SetBranchAddress("Flag_HBHENoiseFilter", &Flag_HBHENoiseFilter, &b_Flag_HBHENoiseFilter);
   fChain->SetBranchAddress("Flag_HBHENoiseIsoFilter", &Flag_HBHENoiseIsoFilter, &b_Flag_HBHENoiseIsoFilter);
   fChain->SetBranchAddress("Flag_CSCTightHaloFilter", &Flag_CSCTightHaloFilter, &b_Flag_CSCTightHaloFilter);
   fChain->SetBranchAddress("Flag_CSCTightHaloTrkMuUnvetoFilter", &Flag_CSCTightHaloTrkMuUnvetoFilter, &b_Flag_CSCTightHaloTrkMuUnvetoFilter);
   fChain->SetBranchAddress("Flag_CSCTightHalo2015Filter", &Flag_CSCTightHalo2015Filter, &b_Flag_CSCTightHalo2015Filter);
   fChain->SetBranchAddress("Flag_globalTightHalo2016Filter", &Flag_globalTightHalo2016Filter, &b_Flag_globalTightHalo2016Filter);
   fChain->SetBranchAddress("Flag_globalSuperTightHalo2016Filter", &Flag_globalSuperTightHalo2016Filter, &b_Flag_globalSuperTightHalo2016Filter);
   fChain->SetBranchAddress("Flag_HcalStripHaloFilter", &Flag_HcalStripHaloFilter, &b_Flag_HcalStripHaloFilter);
   fChain->SetBranchAddress("Flag_hcalLaserEventFilter", &Flag_hcalLaserEventFilter, &b_Flag_hcalLaserEventFilter);
   fChain->SetBranchAddress("Flag_EcalDeadCellTriggerPrimitiveFilter", &Flag_EcalDeadCellTriggerPrimitiveFilter, &b_Flag_EcalDeadCellTriggerPrimitiveFilter);
   fChain->SetBranchAddress("Flag_EcalDeadCellBoundaryEnergyFilter", &Flag_EcalDeadCellBoundaryEnergyFilter, &b_Flag_EcalDeadCellBoundaryEnergyFilter);
   fChain->SetBranchAddress("Flag_ecalBadCalibFilter", &Flag_ecalBadCalibFilter, &b_Flag_ecalBadCalibFilter);
   fChain->SetBranchAddress("Flag_goodVertices", &Flag_goodVertices, &b_Flag_goodVertices);
   fChain->SetBranchAddress("Flag_eeBadScFilter", &Flag_eeBadScFilter, &b_Flag_eeBadScFilter);
   fChain->SetBranchAddress("Flag_ecalLaserCorrFilter", &Flag_ecalLaserCorrFilter, &b_Flag_ecalLaserCorrFilter);
   fChain->SetBranchAddress("Flag_trkPOGFilters", &Flag_trkPOGFilters, &b_Flag_trkPOGFilters);
   fChain->SetBranchAddress("Flag_chargedHadronTrackResolutionFilter", &Flag_chargedHadronTrackResolutionFilter, &b_Flag_chargedHadronTrackResolutionFilter);
   fChain->SetBranchAddress("Flag_muonBadTrackFilter", &Flag_muonBadTrackFilter, &b_Flag_muonBadTrackFilter);
   fChain->SetBranchAddress("Flag_BadChargedCandidateFilter", &Flag_BadChargedCandidateFilter, &b_Flag_BadChargedCandidateFilter);
   fChain->SetBranchAddress("Flag_BadPFMuonFilter", &Flag_BadPFMuonFilter, &b_Flag_BadPFMuonFilter);
   fChain->SetBranchAddress("Flag_BadPFMuonDzFilter", &Flag_BadPFMuonDzFilter, &b_Flag_BadPFMuonDzFilter);
   fChain->SetBranchAddress("Flag_hfNoisyHitsFilter", &Flag_hfNoisyHitsFilter, &b_Flag_hfNoisyHitsFilter);
   fChain->SetBranchAddress("Flag_BadChargedCandidateSummer16Filter", &Flag_BadChargedCandidateSummer16Filter, &b_Flag_BadChargedCandidateSummer16Filter);
   fChain->SetBranchAddress("Flag_BadPFMuonSummer16Filter", &Flag_BadPFMuonSummer16Filter, &b_Flag_BadPFMuonSummer16Filter);
   fChain->SetBranchAddress("Flag_trkPOG_manystripclus53X", &Flag_trkPOG_manystripclus53X, &b_Flag_trkPOG_manystripclus53X);
   fChain->SetBranchAddress("Flag_trkPOG_toomanystripclus53X", &Flag_trkPOG_toomanystripclus53X, &b_Flag_trkPOG_toomanystripclus53X);
   fChain->SetBranchAddress("Flag_trkPOG_logErrorTooManyClusters", &Flag_trkPOG_logErrorTooManyClusters, &b_Flag_trkPOG_logErrorTooManyClusters);
   fChain->SetBranchAddress("Flag_METFilters", &Flag_METFilters, &b_Flag_METFilters);
   fChain->SetBranchAddress("L1_AlwaysTrue", &L1_AlwaysTrue, &b_L1_AlwaysTrue);
   fChain->SetBranchAddress("L1_BPTX_AND_Ref1_VME", &L1_BPTX_AND_Ref1_VME, &b_L1_BPTX_AND_Ref1_VME);
   fChain->SetBranchAddress("L1_BPTX_AND_Ref3_VME", &L1_BPTX_AND_Ref3_VME, &b_L1_BPTX_AND_Ref3_VME);
   fChain->SetBranchAddress("L1_BPTX_AND_Ref4_VME", &L1_BPTX_AND_Ref4_VME, &b_L1_BPTX_AND_Ref4_VME);
   fChain->SetBranchAddress("L1_BPTX_BeamGas_B1_VME", &L1_BPTX_BeamGas_B1_VME, &b_L1_BPTX_BeamGas_B1_VME);
   fChain->SetBranchAddress("L1_BPTX_BeamGas_B2_VME", &L1_BPTX_BeamGas_B2_VME, &b_L1_BPTX_BeamGas_B2_VME);
   fChain->SetBranchAddress("L1_BPTX_BeamGas_Ref1_VME", &L1_BPTX_BeamGas_Ref1_VME, &b_L1_BPTX_BeamGas_Ref1_VME);
   fChain->SetBranchAddress("L1_BPTX_BeamGas_Ref2_VME", &L1_BPTX_BeamGas_Ref2_VME, &b_L1_BPTX_BeamGas_Ref2_VME);
   fChain->SetBranchAddress("L1_BPTX_NotOR_VME", &L1_BPTX_NotOR_VME, &b_L1_BPTX_NotOR_VME);
   fChain->SetBranchAddress("L1_BPTX_OR_Ref3_VME", &L1_BPTX_OR_Ref3_VME, &b_L1_BPTX_OR_Ref3_VME);
   fChain->SetBranchAddress("L1_BPTX_OR_Ref4_VME", &L1_BPTX_OR_Ref4_VME, &b_L1_BPTX_OR_Ref4_VME);
   fChain->SetBranchAddress("L1_BPTX_RefAND_VME", &L1_BPTX_RefAND_VME, &b_L1_BPTX_RefAND_VME);
   fChain->SetBranchAddress("L1_BptxMinus", &L1_BptxMinus, &b_L1_BptxMinus);
   fChain->SetBranchAddress("L1_BptxOR", &L1_BptxOR, &b_L1_BptxOR);
   fChain->SetBranchAddress("L1_BptxPlus", &L1_BptxPlus, &b_L1_BptxPlus);
   fChain->SetBranchAddress("L1_BptxXOR", &L1_BptxXOR, &b_L1_BptxXOR);
   fChain->SetBranchAddress("L1_CDC_SingleMu_3_er1p2_TOP120_DPHI2p618_3p142", &L1_CDC_SingleMu_3_er1p2_TOP120_DPHI2p618_3p142, &b_L1_CDC_SingleMu_3_er1p2_TOP120_DPHI2p618_3p142);
   fChain->SetBranchAddress("L1_DoubleEG10_er1p2_dR_Max0p6", &L1_DoubleEG10_er1p2_dR_Max0p6, &b_L1_DoubleEG10_er1p2_dR_Max0p6);
   fChain->SetBranchAddress("L1_DoubleEG10p5_er1p2_dR_Max0p6", &L1_DoubleEG10p5_er1p2_dR_Max0p6, &b_L1_DoubleEG10p5_er1p2_dR_Max0p6);
   fChain->SetBranchAddress("L1_DoubleEG11_er1p2_dR_Max0p6", &L1_DoubleEG11_er1p2_dR_Max0p6, &b_L1_DoubleEG11_er1p2_dR_Max0p6);
   fChain->SetBranchAddress("L1_DoubleEG4_er1p2_dR_Max0p9", &L1_DoubleEG4_er1p2_dR_Max0p9, &b_L1_DoubleEG4_er1p2_dR_Max0p9);
   fChain->SetBranchAddress("L1_DoubleEG4p5_er1p2_dR_Max0p9", &L1_DoubleEG4p5_er1p2_dR_Max0p9, &b_L1_DoubleEG4p5_er1p2_dR_Max0p9);
   fChain->SetBranchAddress("L1_DoubleEG5_er1p2_dR_Max0p9", &L1_DoubleEG5_er1p2_dR_Max0p9, &b_L1_DoubleEG5_er1p2_dR_Max0p9);
   fChain->SetBranchAddress("L1_DoubleEG5p5_er1p2_dR_Max0p8", &L1_DoubleEG5p5_er1p2_dR_Max0p8, &b_L1_DoubleEG5p5_er1p2_dR_Max0p8);
   fChain->SetBranchAddress("L1_DoubleEG6_er1p2_dR_Max0p8", &L1_DoubleEG6_er1p2_dR_Max0p8, &b_L1_DoubleEG6_er1p2_dR_Max0p8);
   fChain->SetBranchAddress("L1_DoubleEG6p5_er1p2_dR_Max0p8", &L1_DoubleEG6p5_er1p2_dR_Max0p8, &b_L1_DoubleEG6p5_er1p2_dR_Max0p8);
   fChain->SetBranchAddress("L1_DoubleEG7_er1p2_dR_Max0p8", &L1_DoubleEG7_er1p2_dR_Max0p8, &b_L1_DoubleEG7_er1p2_dR_Max0p8);
   fChain->SetBranchAddress("L1_DoubleEG7p5_er1p2_dR_Max0p7", &L1_DoubleEG7p5_er1p2_dR_Max0p7, &b_L1_DoubleEG7p5_er1p2_dR_Max0p7);
   fChain->SetBranchAddress("L1_DoubleEG8_er1p2_dR_Max0p7", &L1_DoubleEG8_er1p2_dR_Max0p7, &b_L1_DoubleEG8_er1p2_dR_Max0p7);
   fChain->SetBranchAddress("L1_DoubleEG8er2p5_HTT260er", &L1_DoubleEG8er2p5_HTT260er, &b_L1_DoubleEG8er2p5_HTT260er);
   fChain->SetBranchAddress("L1_DoubleEG8er2p5_HTT280er", &L1_DoubleEG8er2p5_HTT280er, &b_L1_DoubleEG8er2p5_HTT280er);
   fChain->SetBranchAddress("L1_DoubleEG8er2p5_HTT300er", &L1_DoubleEG8er2p5_HTT300er, &b_L1_DoubleEG8er2p5_HTT300er);
   fChain->SetBranchAddress("L1_DoubleEG8er2p5_HTT320er", &L1_DoubleEG8er2p5_HTT320er, &b_L1_DoubleEG8er2p5_HTT320er);
   fChain->SetBranchAddress("L1_DoubleEG8er2p5_HTT340er", &L1_DoubleEG8er2p5_HTT340er, &b_L1_DoubleEG8er2p5_HTT340er);
   fChain->SetBranchAddress("L1_DoubleEG8p5_er1p2_dR_Max0p7", &L1_DoubleEG8p5_er1p2_dR_Max0p7, &b_L1_DoubleEG8p5_er1p2_dR_Max0p7);
   fChain->SetBranchAddress("L1_DoubleEG9_er1p2_dR_Max0p7", &L1_DoubleEG9_er1p2_dR_Max0p7, &b_L1_DoubleEG9_er1p2_dR_Max0p7);
   fChain->SetBranchAddress("L1_DoubleEG9p5_er1p2_dR_Max0p6", &L1_DoubleEG9p5_er1p2_dR_Max0p6, &b_L1_DoubleEG9p5_er1p2_dR_Max0p6);
   fChain->SetBranchAddress("L1_DoubleEG_15_10_er2p5", &L1_DoubleEG_15_10_er2p5, &b_L1_DoubleEG_15_10_er2p5);
   fChain->SetBranchAddress("L1_DoubleEG_20_10_er2p5", &L1_DoubleEG_20_10_er2p5, &b_L1_DoubleEG_20_10_er2p5);
   fChain->SetBranchAddress("L1_DoubleEG_22_10_er2p5", &L1_DoubleEG_22_10_er2p5, &b_L1_DoubleEG_22_10_er2p5);
   fChain->SetBranchAddress("L1_DoubleEG_25_12_er2p5", &L1_DoubleEG_25_12_er2p5, &b_L1_DoubleEG_25_12_er2p5);
   fChain->SetBranchAddress("L1_DoubleEG_25_14_er2p5", &L1_DoubleEG_25_14_er2p5, &b_L1_DoubleEG_25_14_er2p5);
   fChain->SetBranchAddress("L1_DoubleEG_27_14_er2p5", &L1_DoubleEG_27_14_er2p5, &b_L1_DoubleEG_27_14_er2p5);
   fChain->SetBranchAddress("L1_DoubleEG_LooseIso16_LooseIso12_er1p5", &L1_DoubleEG_LooseIso16_LooseIso12_er1p5, &b_L1_DoubleEG_LooseIso16_LooseIso12_er1p5);
   fChain->SetBranchAddress("L1_DoubleEG_LooseIso18_LooseIso12_er1p5", &L1_DoubleEG_LooseIso18_LooseIso12_er1p5, &b_L1_DoubleEG_LooseIso18_LooseIso12_er1p5);
   fChain->SetBranchAddress("L1_DoubleEG_LooseIso20_10_er2p5", &L1_DoubleEG_LooseIso20_10_er2p5, &b_L1_DoubleEG_LooseIso20_10_er2p5);
   fChain->SetBranchAddress("L1_DoubleEG_LooseIso20_LooseIso12_er1p5", &L1_DoubleEG_LooseIso20_LooseIso12_er1p5, &b_L1_DoubleEG_LooseIso20_LooseIso12_er1p5);
   fChain->SetBranchAddress("L1_DoubleEG_LooseIso22_10_er2p5", &L1_DoubleEG_LooseIso22_10_er2p5, &b_L1_DoubleEG_LooseIso22_10_er2p5);
   fChain->SetBranchAddress("L1_DoubleEG_LooseIso22_12_er2p5", &L1_DoubleEG_LooseIso22_12_er2p5, &b_L1_DoubleEG_LooseIso22_12_er2p5);
   fChain->SetBranchAddress("L1_DoubleEG_LooseIso22_LooseIso12_er1p5", &L1_DoubleEG_LooseIso22_LooseIso12_er1p5, &b_L1_DoubleEG_LooseIso22_LooseIso12_er1p5);
   fChain->SetBranchAddress("L1_DoubleEG_LooseIso25_12_er2p5", &L1_DoubleEG_LooseIso25_12_er2p5, &b_L1_DoubleEG_LooseIso25_12_er2p5);
   fChain->SetBranchAddress("L1_DoubleEG_LooseIso25_LooseIso12_er1p5", &L1_DoubleEG_LooseIso25_LooseIso12_er1p5, &b_L1_DoubleEG_LooseIso25_LooseIso12_er1p5);
   fChain->SetBranchAddress("L1_DoubleIsoTau26er2p1_Jet55_RmOvlp_dR0p5", &L1_DoubleIsoTau26er2p1_Jet55_RmOvlp_dR0p5, &b_L1_DoubleIsoTau26er2p1_Jet55_RmOvlp_dR0p5);
   fChain->SetBranchAddress("L1_DoubleIsoTau26er2p1_Jet70_RmOvlp_dR0p5", &L1_DoubleIsoTau26er2p1_Jet70_RmOvlp_dR0p5, &b_L1_DoubleIsoTau26er2p1_Jet70_RmOvlp_dR0p5);
   fChain->SetBranchAddress("L1_DoubleIsoTau28er2p1", &L1_DoubleIsoTau28er2p1, &b_L1_DoubleIsoTau28er2p1);
   fChain->SetBranchAddress("L1_DoubleIsoTau28er2p1_Mass_Max80", &L1_DoubleIsoTau28er2p1_Mass_Max80, &b_L1_DoubleIsoTau28er2p1_Mass_Max80);
   fChain->SetBranchAddress("L1_DoubleIsoTau28er2p1_Mass_Max90", &L1_DoubleIsoTau28er2p1_Mass_Max90, &b_L1_DoubleIsoTau28er2p1_Mass_Max90);
   fChain->SetBranchAddress("L1_DoubleIsoTau30er2p1", &L1_DoubleIsoTau30er2p1, &b_L1_DoubleIsoTau30er2p1);
   fChain->SetBranchAddress("L1_DoubleIsoTau30er2p1_Mass_Max80", &L1_DoubleIsoTau30er2p1_Mass_Max80, &b_L1_DoubleIsoTau30er2p1_Mass_Max80);
   fChain->SetBranchAddress("L1_DoubleIsoTau30er2p1_Mass_Max90", &L1_DoubleIsoTau30er2p1_Mass_Max90, &b_L1_DoubleIsoTau30er2p1_Mass_Max90);
   fChain->SetBranchAddress("L1_DoubleIsoTau32er2p1", &L1_DoubleIsoTau32er2p1, &b_L1_DoubleIsoTau32er2p1);
   fChain->SetBranchAddress("L1_DoubleIsoTau34er2p1", &L1_DoubleIsoTau34er2p1, &b_L1_DoubleIsoTau34er2p1);
   fChain->SetBranchAddress("L1_DoubleIsoTau35er2p1", &L1_DoubleIsoTau35er2p1, &b_L1_DoubleIsoTau35er2p1);
   fChain->SetBranchAddress("L1_DoubleIsoTau36er2p1", &L1_DoubleIsoTau36er2p1, &b_L1_DoubleIsoTau36er2p1);
   fChain->SetBranchAddress("L1_DoubleJet100er2p3_dEta_Max1p6", &L1_DoubleJet100er2p3_dEta_Max1p6, &b_L1_DoubleJet100er2p3_dEta_Max1p6);
   fChain->SetBranchAddress("L1_DoubleJet100er2p5", &L1_DoubleJet100er2p5, &b_L1_DoubleJet100er2p5);
   fChain->SetBranchAddress("L1_DoubleJet112er2p3_dEta_Max1p6", &L1_DoubleJet112er2p3_dEta_Max1p6, &b_L1_DoubleJet112er2p3_dEta_Max1p6);
   fChain->SetBranchAddress("L1_DoubleJet120er2p5", &L1_DoubleJet120er2p5, &b_L1_DoubleJet120er2p5);
   fChain->SetBranchAddress("L1_DoubleJet150er2p5", &L1_DoubleJet150er2p5, &b_L1_DoubleJet150er2p5);
   fChain->SetBranchAddress("L1_DoubleJet30er2p5_Mass_Min150_dEta_Max1p5", &L1_DoubleJet30er2p5_Mass_Min150_dEta_Max1p5, &b_L1_DoubleJet30er2p5_Mass_Min150_dEta_Max1p5);
   fChain->SetBranchAddress("L1_DoubleJet30er2p5_Mass_Min200_dEta_Max1p5", &L1_DoubleJet30er2p5_Mass_Min200_dEta_Max1p5, &b_L1_DoubleJet30er2p5_Mass_Min200_dEta_Max1p5);
   fChain->SetBranchAddress("L1_DoubleJet30er2p5_Mass_Min250_dEta_Max1p5", &L1_DoubleJet30er2p5_Mass_Min250_dEta_Max1p5, &b_L1_DoubleJet30er2p5_Mass_Min250_dEta_Max1p5);
   fChain->SetBranchAddress("L1_DoubleJet30er2p5_Mass_Min300_dEta_Max1p5", &L1_DoubleJet30er2p5_Mass_Min300_dEta_Max1p5, &b_L1_DoubleJet30er2p5_Mass_Min300_dEta_Max1p5);
   fChain->SetBranchAddress("L1_DoubleJet30er2p5_Mass_Min330_dEta_Max1p5", &L1_DoubleJet30er2p5_Mass_Min330_dEta_Max1p5, &b_L1_DoubleJet30er2p5_Mass_Min330_dEta_Max1p5);
   fChain->SetBranchAddress("L1_DoubleJet30er2p5_Mass_Min360_dEta_Max1p5", &L1_DoubleJet30er2p5_Mass_Min360_dEta_Max1p5, &b_L1_DoubleJet30er2p5_Mass_Min360_dEta_Max1p5);
   fChain->SetBranchAddress("L1_DoubleJet35_Mass_Min450_IsoTau45_RmOvlp", &L1_DoubleJet35_Mass_Min450_IsoTau45_RmOvlp, &b_L1_DoubleJet35_Mass_Min450_IsoTau45_RmOvlp);
   fChain->SetBranchAddress("L1_DoubleJet35_Mass_Min450_IsoTau45er2p1_RmOvlp_dR0p5", &L1_DoubleJet35_Mass_Min450_IsoTau45er2p1_RmOvlp_dR0p5, &b_L1_DoubleJet35_Mass_Min450_IsoTau45er2p1_RmOvlp_dR0p5);
   fChain->SetBranchAddress("L1_DoubleJet40er2p5", &L1_DoubleJet40er2p5, &b_L1_DoubleJet40er2p5);
   fChain->SetBranchAddress("L1_DoubleJet_100_30_DoubleJet30_Mass_Min620", &L1_DoubleJet_100_30_DoubleJet30_Mass_Min620, &b_L1_DoubleJet_100_30_DoubleJet30_Mass_Min620);
   fChain->SetBranchAddress("L1_DoubleJet_110_35_DoubleJet35_Mass_Min620", &L1_DoubleJet_110_35_DoubleJet35_Mass_Min620, &b_L1_DoubleJet_110_35_DoubleJet35_Mass_Min620);
   fChain->SetBranchAddress("L1_DoubleJet_115_40_DoubleJet40_Mass_Min620", &L1_DoubleJet_115_40_DoubleJet40_Mass_Min620, &b_L1_DoubleJet_115_40_DoubleJet40_Mass_Min620);
   fChain->SetBranchAddress("L1_DoubleJet_115_40_DoubleJet40_Mass_Min620_Jet60TT28", &L1_DoubleJet_115_40_DoubleJet40_Mass_Min620_Jet60TT28, &b_L1_DoubleJet_115_40_DoubleJet40_Mass_Min620_Jet60TT28);
   fChain->SetBranchAddress("L1_DoubleJet_120_45_DoubleJet45_Mass_Min620", &L1_DoubleJet_120_45_DoubleJet45_Mass_Min620, &b_L1_DoubleJet_120_45_DoubleJet45_Mass_Min620);
   fChain->SetBranchAddress("L1_DoubleJet_120_45_DoubleJet45_Mass_Min620_Jet60TT28", &L1_DoubleJet_120_45_DoubleJet45_Mass_Min620_Jet60TT28, &b_L1_DoubleJet_120_45_DoubleJet45_Mass_Min620_Jet60TT28);
   fChain->SetBranchAddress("L1_DoubleJet_80_30_Mass_Min420_DoubleMu0_SQ", &L1_DoubleJet_80_30_Mass_Min420_DoubleMu0_SQ, &b_L1_DoubleJet_80_30_Mass_Min420_DoubleMu0_SQ);
   fChain->SetBranchAddress("L1_DoubleJet_80_30_Mass_Min420_IsoTau40_RmOvlp", &L1_DoubleJet_80_30_Mass_Min420_IsoTau40_RmOvlp, &b_L1_DoubleJet_80_30_Mass_Min420_IsoTau40_RmOvlp);
   fChain->SetBranchAddress("L1_DoubleJet_80_30_Mass_Min420_Mu8", &L1_DoubleJet_80_30_Mass_Min420_Mu8, &b_L1_DoubleJet_80_30_Mass_Min420_Mu8);
   fChain->SetBranchAddress("L1_DoubleJet_90_30_DoubleJet30_Mass_Min620", &L1_DoubleJet_90_30_DoubleJet30_Mass_Min620, &b_L1_DoubleJet_90_30_DoubleJet30_Mass_Min620);
   fChain->SetBranchAddress("L1_DoubleLLPJet40", &L1_DoubleLLPJet40, &b_L1_DoubleLLPJet40);
   fChain->SetBranchAddress("L1_DoubleLooseIsoEG22er2p1", &L1_DoubleLooseIsoEG22er2p1, &b_L1_DoubleLooseIsoEG22er2p1);
   fChain->SetBranchAddress("L1_DoubleLooseIsoEG24er2p1", &L1_DoubleLooseIsoEG24er2p1, &b_L1_DoubleLooseIsoEG24er2p1);
   fChain->SetBranchAddress("L1_DoubleMu0", &L1_DoubleMu0, &b_L1_DoubleMu0);
   fChain->SetBranchAddress("L1_DoubleMu0_Mass_Min1", &L1_DoubleMu0_Mass_Min1, &b_L1_DoubleMu0_Mass_Min1);
   fChain->SetBranchAddress("L1_DoubleMu0_OQ", &L1_DoubleMu0_OQ, &b_L1_DoubleMu0_OQ);
   fChain->SetBranchAddress("L1_DoubleMu0_SQ", &L1_DoubleMu0_SQ, &b_L1_DoubleMu0_SQ);
   fChain->SetBranchAddress("L1_DoubleMu0_SQ_OS", &L1_DoubleMu0_SQ_OS, &b_L1_DoubleMu0_SQ_OS);
   fChain->SetBranchAddress("L1_DoubleMu0_Upt15_Upt7", &L1_DoubleMu0_Upt15_Upt7, &b_L1_DoubleMu0_Upt15_Upt7);
   fChain->SetBranchAddress("L1_DoubleMu0_Upt5_Upt5", &L1_DoubleMu0_Upt5_Upt5, &b_L1_DoubleMu0_Upt5_Upt5);
   fChain->SetBranchAddress("L1_DoubleMu0_Upt6_IP_Min1_Upt4", &L1_DoubleMu0_Upt6_IP_Min1_Upt4, &b_L1_DoubleMu0_Upt6_IP_Min1_Upt4);
   fChain->SetBranchAddress("L1_DoubleMu0_dR_Max1p6_Jet90er2p5_dR_Max0p8", &L1_DoubleMu0_dR_Max1p6_Jet90er2p5_dR_Max0p8, &b_L1_DoubleMu0_dR_Max1p6_Jet90er2p5_dR_Max0p8);
   fChain->SetBranchAddress("L1_DoubleMu0er1p4_OQ_OS_dEta_Max1p6", &L1_DoubleMu0er1p4_OQ_OS_dEta_Max1p6, &b_L1_DoubleMu0er1p4_OQ_OS_dEta_Max1p6);
   fChain->SetBranchAddress("L1_DoubleMu0er1p4_SQ_OS_dR_Max1p4", &L1_DoubleMu0er1p4_SQ_OS_dR_Max1p4, &b_L1_DoubleMu0er1p4_SQ_OS_dR_Max1p4);
   fChain->SetBranchAddress("L1_DoubleMu0er1p5_SQ", &L1_DoubleMu0er1p5_SQ, &b_L1_DoubleMu0er1p5_SQ);
   fChain->SetBranchAddress("L1_DoubleMu0er1p5_SQ_OS", &L1_DoubleMu0er1p5_SQ_OS, &b_L1_DoubleMu0er1p5_SQ_OS);
   fChain->SetBranchAddress("L1_DoubleMu0er1p5_SQ_OS_dR_Max1p4", &L1_DoubleMu0er1p5_SQ_OS_dR_Max1p4, &b_L1_DoubleMu0er1p5_SQ_OS_dR_Max1p4);
   fChain->SetBranchAddress("L1_DoubleMu0er1p5_SQ_dR_Max1p4", &L1_DoubleMu0er1p5_SQ_dR_Max1p4, &b_L1_DoubleMu0er1p5_SQ_dR_Max1p4);
   fChain->SetBranchAddress("L1_DoubleMu0er2p0_SQ_OS_dEta_Max1p5", &L1_DoubleMu0er2p0_SQ_OS_dEta_Max1p5, &b_L1_DoubleMu0er2p0_SQ_OS_dEta_Max1p5);
   fChain->SetBranchAddress("L1_DoubleMu0er2p0_SQ_OS_dEta_Max1p6", &L1_DoubleMu0er2p0_SQ_OS_dEta_Max1p6, &b_L1_DoubleMu0er2p0_SQ_OS_dEta_Max1p6);
   fChain->SetBranchAddress("L1_DoubleMu0er2p0_SQ_OS_dR_Max1p4", &L1_DoubleMu0er2p0_SQ_OS_dR_Max1p4, &b_L1_DoubleMu0er2p0_SQ_OS_dR_Max1p4);
   fChain->SetBranchAddress("L1_DoubleMu0er2p0_SQ_dEta_Max1p5", &L1_DoubleMu0er2p0_SQ_dEta_Max1p5, &b_L1_DoubleMu0er2p0_SQ_dEta_Max1p5);
   fChain->SetBranchAddress("L1_DoubleMu0er2p0_SQ_dEta_Max1p6", &L1_DoubleMu0er2p0_SQ_dEta_Max1p6, &b_L1_DoubleMu0er2p0_SQ_dEta_Max1p6);
   fChain->SetBranchAddress("L1_DoubleMu0er2p0_SQ_dR_Max1p4", &L1_DoubleMu0er2p0_SQ_dR_Max1p4, &b_L1_DoubleMu0er2p0_SQ_dR_Max1p4);
   fChain->SetBranchAddress("L1_DoubleMu18er2p1_SQ", &L1_DoubleMu18er2p1_SQ, &b_L1_DoubleMu18er2p1_SQ);
   fChain->SetBranchAddress("L1_DoubleMu3_OS_er2p3_Mass_Max14_DoubleEG7p5_er2p1_Mass_Max20", &L1_DoubleMu3_OS_er2p3_Mass_Max14_DoubleEG7p5_er2p1_Mass_Max20, &b_L1_DoubleMu3_OS_er2p3_Mass_Max14_DoubleEG7p5_er2p1_Mass_Max20);
   fChain->SetBranchAddress("L1_DoubleMu3_SQ_ETMHF50_HTT60er", &L1_DoubleMu3_SQ_ETMHF50_HTT60er, &b_L1_DoubleMu3_SQ_ETMHF50_HTT60er);
   fChain->SetBranchAddress("L1_DoubleMu3_SQ_ETMHF50_Jet60er2p5", &L1_DoubleMu3_SQ_ETMHF50_Jet60er2p5, &b_L1_DoubleMu3_SQ_ETMHF50_Jet60er2p5);
   fChain->SetBranchAddress("L1_DoubleMu3_SQ_ETMHF50_Jet60er2p5_OR_DoubleJet40er2p5", &L1_DoubleMu3_SQ_ETMHF50_Jet60er2p5_OR_DoubleJet40er2p5, &b_L1_DoubleMu3_SQ_ETMHF50_Jet60er2p5_OR_DoubleJet40er2p5);
   fChain->SetBranchAddress("L1_DoubleMu3_SQ_ETMHF60_Jet60er2p5", &L1_DoubleMu3_SQ_ETMHF60_Jet60er2p5, &b_L1_DoubleMu3_SQ_ETMHF60_Jet60er2p5);
   fChain->SetBranchAddress("L1_DoubleMu3_SQ_HTT220er", &L1_DoubleMu3_SQ_HTT220er, &b_L1_DoubleMu3_SQ_HTT220er);
   fChain->SetBranchAddress("L1_DoubleMu3_SQ_HTT240er", &L1_DoubleMu3_SQ_HTT240er, &b_L1_DoubleMu3_SQ_HTT240er);
   fChain->SetBranchAddress("L1_DoubleMu3_SQ_HTT260er", &L1_DoubleMu3_SQ_HTT260er, &b_L1_DoubleMu3_SQ_HTT260er);
   fChain->SetBranchAddress("L1_DoubleMu3_dR_Max1p6_Jet90er2p5_dR_Max0p8", &L1_DoubleMu3_dR_Max1p6_Jet90er2p5_dR_Max0p8, &b_L1_DoubleMu3_dR_Max1p6_Jet90er2p5_dR_Max0p8);
   fChain->SetBranchAddress("L1_DoubleMu3er2p0_SQ_OS_dR_Max1p4", &L1_DoubleMu3er2p0_SQ_OS_dR_Max1p4, &b_L1_DoubleMu3er2p0_SQ_OS_dR_Max1p4);
   fChain->SetBranchAddress("L1_DoubleMu4_SQ_EG9er2p5", &L1_DoubleMu4_SQ_EG9er2p5, &b_L1_DoubleMu4_SQ_EG9er2p5);
   fChain->SetBranchAddress("L1_DoubleMu4_SQ_OS", &L1_DoubleMu4_SQ_OS, &b_L1_DoubleMu4_SQ_OS);
   fChain->SetBranchAddress("L1_DoubleMu4_SQ_OS_dR_Max1p2", &L1_DoubleMu4_SQ_OS_dR_Max1p2, &b_L1_DoubleMu4_SQ_OS_dR_Max1p2);
   fChain->SetBranchAddress("L1_DoubleMu4p5_SQ_OS", &L1_DoubleMu4p5_SQ_OS, &b_L1_DoubleMu4p5_SQ_OS);
   fChain->SetBranchAddress("L1_DoubleMu4p5_SQ_OS_dR_Max1p2", &L1_DoubleMu4p5_SQ_OS_dR_Max1p2, &b_L1_DoubleMu4p5_SQ_OS_dR_Max1p2);
   fChain->SetBranchAddress("L1_DoubleMu4p5er2p0_SQ_OS", &L1_DoubleMu4p5er2p0_SQ_OS, &b_L1_DoubleMu4p5er2p0_SQ_OS);
   fChain->SetBranchAddress("L1_DoubleMu4p5er2p0_SQ_OS_Mass_7to18", &L1_DoubleMu4p5er2p0_SQ_OS_Mass_7to18, &b_L1_DoubleMu4p5er2p0_SQ_OS_Mass_7to18);
   fChain->SetBranchAddress("L1_DoubleMu4p5er2p0_SQ_OS_Mass_Min7", &L1_DoubleMu4p5er2p0_SQ_OS_Mass_Min7, &b_L1_DoubleMu4p5er2p0_SQ_OS_Mass_Min7);
   fChain->SetBranchAddress("L1_DoubleMu5_OS_er2p3_Mass_8to14_DoubleEG3er2p1_Mass_Max20", &L1_DoubleMu5_OS_er2p3_Mass_8to14_DoubleEG3er2p1_Mass_Max20, &b_L1_DoubleMu5_OS_er2p3_Mass_8to14_DoubleEG3er2p1_Mass_Max20);
   fChain->SetBranchAddress("L1_DoubleMu5_SQ_EG9er2p5", &L1_DoubleMu5_SQ_EG9er2p5, &b_L1_DoubleMu5_SQ_EG9er2p5);
   fChain->SetBranchAddress("L1_DoubleMu8_SQ", &L1_DoubleMu8_SQ, &b_L1_DoubleMu8_SQ);
   fChain->SetBranchAddress("L1_DoubleMu9_SQ", &L1_DoubleMu9_SQ, &b_L1_DoubleMu9_SQ);
   fChain->SetBranchAddress("L1_DoubleMu_12_5", &L1_DoubleMu_12_5, &b_L1_DoubleMu_12_5);
   fChain->SetBranchAddress("L1_DoubleMu_15_5_SQ", &L1_DoubleMu_15_5_SQ, &b_L1_DoubleMu_15_5_SQ);
   fChain->SetBranchAddress("L1_DoubleMu_15_7", &L1_DoubleMu_15_7, &b_L1_DoubleMu_15_7);
   fChain->SetBranchAddress("L1_DoubleMu_15_7_Mass_Min1", &L1_DoubleMu_15_7_Mass_Min1, &b_L1_DoubleMu_15_7_Mass_Min1);
   fChain->SetBranchAddress("L1_DoubleMu_15_7_SQ", &L1_DoubleMu_15_7_SQ, &b_L1_DoubleMu_15_7_SQ);
   fChain->SetBranchAddress("L1_DoubleTau70er2p1", &L1_DoubleTau70er2p1, &b_L1_DoubleTau70er2p1);
   fChain->SetBranchAddress("L1_ETM120", &L1_ETM120, &b_L1_ETM120);
   fChain->SetBranchAddress("L1_ETM150", &L1_ETM150, &b_L1_ETM150);
   fChain->SetBranchAddress("L1_ETMHF100", &L1_ETMHF100, &b_L1_ETMHF100);
   fChain->SetBranchAddress("L1_ETMHF100_HTT60er", &L1_ETMHF100_HTT60er, &b_L1_ETMHF100_HTT60er);
   fChain->SetBranchAddress("L1_ETMHF110", &L1_ETMHF110, &b_L1_ETMHF110);
   fChain->SetBranchAddress("L1_ETMHF110_HTT60er", &L1_ETMHF110_HTT60er, &b_L1_ETMHF110_HTT60er);
   fChain->SetBranchAddress("L1_ETMHF110_HTT60er_NotSecondBunchInTrain", &L1_ETMHF110_HTT60er_NotSecondBunchInTrain, &b_L1_ETMHF110_HTT60er_NotSecondBunchInTrain);
   fChain->SetBranchAddress("L1_ETMHF120", &L1_ETMHF120, &b_L1_ETMHF120);
   fChain->SetBranchAddress("L1_ETMHF120_HTT60er", &L1_ETMHF120_HTT60er, &b_L1_ETMHF120_HTT60er);
   fChain->SetBranchAddress("L1_ETMHF120_NotSecondBunchInTrain", &L1_ETMHF120_NotSecondBunchInTrain, &b_L1_ETMHF120_NotSecondBunchInTrain);
   fChain->SetBranchAddress("L1_ETMHF130", &L1_ETMHF130, &b_L1_ETMHF130);
   fChain->SetBranchAddress("L1_ETMHF130_HTT60er", &L1_ETMHF130_HTT60er, &b_L1_ETMHF130_HTT60er);
   fChain->SetBranchAddress("L1_ETMHF140", &L1_ETMHF140, &b_L1_ETMHF140);
   fChain->SetBranchAddress("L1_ETMHF150", &L1_ETMHF150, &b_L1_ETMHF150);
   fChain->SetBranchAddress("L1_ETMHF70", &L1_ETMHF70, &b_L1_ETMHF70);
   fChain->SetBranchAddress("L1_ETMHF70_HTT60er", &L1_ETMHF70_HTT60er, &b_L1_ETMHF70_HTT60er);
   fChain->SetBranchAddress("L1_ETMHF80", &L1_ETMHF80, &b_L1_ETMHF80);
   fChain->SetBranchAddress("L1_ETMHF80_HTT60er", &L1_ETMHF80_HTT60er, &b_L1_ETMHF80_HTT60er);
   fChain->SetBranchAddress("L1_ETMHF90", &L1_ETMHF90, &b_L1_ETMHF90);
   fChain->SetBranchAddress("L1_ETMHF90_HTT60er", &L1_ETMHF90_HTT60er, &b_L1_ETMHF90_HTT60er);
   fChain->SetBranchAddress("L1_ETMHF90_SingleJet60er2p5_dPhi_Min2p1", &L1_ETMHF90_SingleJet60er2p5_dPhi_Min2p1, &b_L1_ETMHF90_SingleJet60er2p5_dPhi_Min2p1);
   fChain->SetBranchAddress("L1_ETMHF90_SingleJet60er2p5_dPhi_Min2p6", &L1_ETMHF90_SingleJet60er2p5_dPhi_Min2p6, &b_L1_ETMHF90_SingleJet60er2p5_dPhi_Min2p6);
   fChain->SetBranchAddress("L1_ETMHF90_SingleJet80er2p5_dPhi_Min2p1", &L1_ETMHF90_SingleJet80er2p5_dPhi_Min2p1, &b_L1_ETMHF90_SingleJet80er2p5_dPhi_Min2p1);
   fChain->SetBranchAddress("L1_ETMHF90_SingleJet80er2p5_dPhi_Min2p6", &L1_ETMHF90_SingleJet80er2p5_dPhi_Min2p6, &b_L1_ETMHF90_SingleJet80er2p5_dPhi_Min2p6);
   fChain->SetBranchAddress("L1_ETT1200", &L1_ETT1200, &b_L1_ETT1200);
   fChain->SetBranchAddress("L1_ETT1600", &L1_ETT1600, &b_L1_ETT1600);
   fChain->SetBranchAddress("L1_ETT2000", &L1_ETT2000, &b_L1_ETT2000);
   fChain->SetBranchAddress("L1_FirstBunchAfterTrain", &L1_FirstBunchAfterTrain, &b_L1_FirstBunchAfterTrain);
   fChain->SetBranchAddress("L1_FirstBunchBeforeTrain", &L1_FirstBunchBeforeTrain, &b_L1_FirstBunchBeforeTrain);
   fChain->SetBranchAddress("L1_FirstBunchInTrain", &L1_FirstBunchInTrain, &b_L1_FirstBunchInTrain);
   fChain->SetBranchAddress("L1_FirstCollisionInOrbit", &L1_FirstCollisionInOrbit, &b_L1_FirstCollisionInOrbit);
   fChain->SetBranchAddress("L1_FirstCollisionInTrain", &L1_FirstCollisionInTrain, &b_L1_FirstCollisionInTrain);
   fChain->SetBranchAddress("L1_HCAL_LaserMon_Trig", &L1_HCAL_LaserMon_Trig, &b_L1_HCAL_LaserMon_Trig);
   fChain->SetBranchAddress("L1_HCAL_LaserMon_Veto", &L1_HCAL_LaserMon_Veto, &b_L1_HCAL_LaserMon_Veto);
   fChain->SetBranchAddress("L1_HTT120_SingleLLPJet40", &L1_HTT120_SingleLLPJet40, &b_L1_HTT120_SingleLLPJet40);
   fChain->SetBranchAddress("L1_HTT120er", &L1_HTT120er, &b_L1_HTT120er);
   fChain->SetBranchAddress("L1_HTT160_SingleLLPJet50", &L1_HTT160_SingleLLPJet50, &b_L1_HTT160_SingleLLPJet50);
   fChain->SetBranchAddress("L1_HTT160er", &L1_HTT160er, &b_L1_HTT160er);
   fChain->SetBranchAddress("L1_HTT200_SingleLLPJet60", &L1_HTT200_SingleLLPJet60, &b_L1_HTT200_SingleLLPJet60);
   fChain->SetBranchAddress("L1_HTT200er", &L1_HTT200er, &b_L1_HTT200er);
   fChain->SetBranchAddress("L1_HTT240_SingleLLPJet70", &L1_HTT240_SingleLLPJet70, &b_L1_HTT240_SingleLLPJet70);
   fChain->SetBranchAddress("L1_HTT255er", &L1_HTT255er, &b_L1_HTT255er);
   fChain->SetBranchAddress("L1_HTT280er", &L1_HTT280er, &b_L1_HTT280er);
   fChain->SetBranchAddress("L1_HTT280er_QuadJet_70_55_40_35_er2p5", &L1_HTT280er_QuadJet_70_55_40_35_er2p5, &b_L1_HTT280er_QuadJet_70_55_40_35_er2p5);
   fChain->SetBranchAddress("L1_HTT320er", &L1_HTT320er, &b_L1_HTT320er);
   fChain->SetBranchAddress("L1_HTT320er_QuadJet_70_55_40_40_er2p5", &L1_HTT320er_QuadJet_70_55_40_40_er2p5, &b_L1_HTT320er_QuadJet_70_55_40_40_er2p5);
   fChain->SetBranchAddress("L1_HTT320er_QuadJet_80_60_er2p1_45_40_er2p3", &L1_HTT320er_QuadJet_80_60_er2p1_45_40_er2p3, &b_L1_HTT320er_QuadJet_80_60_er2p1_45_40_er2p3);
   fChain->SetBranchAddress("L1_HTT320er_QuadJet_80_60_er2p1_50_45_er2p3", &L1_HTT320er_QuadJet_80_60_er2p1_50_45_er2p3, &b_L1_HTT320er_QuadJet_80_60_er2p1_50_45_er2p3);
   fChain->SetBranchAddress("L1_HTT360er", &L1_HTT360er, &b_L1_HTT360er);
   fChain->SetBranchAddress("L1_HTT400er", &L1_HTT400er, &b_L1_HTT400er);
   fChain->SetBranchAddress("L1_HTT450er", &L1_HTT450er, &b_L1_HTT450er);
   fChain->SetBranchAddress("L1_IsoEG32er2p5_Mt40", &L1_IsoEG32er2p5_Mt40, &b_L1_IsoEG32er2p5_Mt40);
   fChain->SetBranchAddress("L1_IsoTau52er2p1_QuadJet36er2p5", &L1_IsoTau52er2p1_QuadJet36er2p5, &b_L1_IsoTau52er2p1_QuadJet36er2p5);
   fChain->SetBranchAddress("L1_IsolatedBunch", &L1_IsolatedBunch, &b_L1_IsolatedBunch);
   fChain->SetBranchAddress("L1_LastBunchInTrain", &L1_LastBunchInTrain, &b_L1_LastBunchInTrain);
   fChain->SetBranchAddress("L1_LastCollisionInTrain", &L1_LastCollisionInTrain, &b_L1_LastCollisionInTrain);
   fChain->SetBranchAddress("L1_LooseIsoEG22er2p1_IsoTau26er2p1_dR_Min0p3", &L1_LooseIsoEG22er2p1_IsoTau26er2p1_dR_Min0p3, &b_L1_LooseIsoEG22er2p1_IsoTau26er2p1_dR_Min0p3);
   fChain->SetBranchAddress("L1_LooseIsoEG22er2p1_Tau70er2p1_dR_Min0p3", &L1_LooseIsoEG22er2p1_Tau70er2p1_dR_Min0p3, &b_L1_LooseIsoEG22er2p1_Tau70er2p1_dR_Min0p3);
   fChain->SetBranchAddress("L1_LooseIsoEG24er2p1_HTT100er", &L1_LooseIsoEG24er2p1_HTT100er, &b_L1_LooseIsoEG24er2p1_HTT100er);
   fChain->SetBranchAddress("L1_LooseIsoEG24er2p1_IsoTau27er2p1_dR_Min0p3", &L1_LooseIsoEG24er2p1_IsoTau27er2p1_dR_Min0p3, &b_L1_LooseIsoEG24er2p1_IsoTau27er2p1_dR_Min0p3);
   fChain->SetBranchAddress("L1_LooseIsoEG26er2p1_HTT100er", &L1_LooseIsoEG26er2p1_HTT100er, &b_L1_LooseIsoEG26er2p1_HTT100er);
   fChain->SetBranchAddress("L1_LooseIsoEG26er2p1_Jet34er2p5_dR_Min0p3", &L1_LooseIsoEG26er2p1_Jet34er2p5_dR_Min0p3, &b_L1_LooseIsoEG26er2p1_Jet34er2p5_dR_Min0p3);
   fChain->SetBranchAddress("L1_LooseIsoEG28er2p1_HTT100er", &L1_LooseIsoEG28er2p1_HTT100er, &b_L1_LooseIsoEG28er2p1_HTT100er);
   fChain->SetBranchAddress("L1_LooseIsoEG28er2p1_Jet34er2p5_dR_Min0p3", &L1_LooseIsoEG28er2p1_Jet34er2p5_dR_Min0p3, &b_L1_LooseIsoEG28er2p1_Jet34er2p5_dR_Min0p3);
   fChain->SetBranchAddress("L1_LooseIsoEG30er2p1_HTT100er", &L1_LooseIsoEG30er2p1_HTT100er, &b_L1_LooseIsoEG30er2p1_HTT100er);
   fChain->SetBranchAddress("L1_LooseIsoEG30er2p1_Jet34er2p5_dR_Min0p3", &L1_LooseIsoEG30er2p1_Jet34er2p5_dR_Min0p3, &b_L1_LooseIsoEG30er2p1_Jet34er2p5_dR_Min0p3);
   fChain->SetBranchAddress("L1_MinimumBiasHF0", &L1_MinimumBiasHF0, &b_L1_MinimumBiasHF0);
   fChain->SetBranchAddress("L1_MinimumBiasHF0_AND_BptxAND", &L1_MinimumBiasHF0_AND_BptxAND, &b_L1_MinimumBiasHF0_AND_BptxAND);
   fChain->SetBranchAddress("L1_Mu10er2p3_Jet32er2p3_dR_Max0p4_DoubleJet32er2p3_dEta_Max1p6", &L1_Mu10er2p3_Jet32er2p3_dR_Max0p4_DoubleJet32er2p3_dEta_Max1p6, &b_L1_Mu10er2p3_Jet32er2p3_dR_Max0p4_DoubleJet32er2p3_dEta_Max1p6);
   fChain->SetBranchAddress("L1_Mu12er2p3_Jet40er2p1_dR_Max0p4_DoubleJet40er2p1_dEta_Max1p6", &L1_Mu12er2p3_Jet40er2p1_dR_Max0p4_DoubleJet40er2p1_dEta_Max1p6, &b_L1_Mu12er2p3_Jet40er2p1_dR_Max0p4_DoubleJet40er2p1_dEta_Max1p6);
   fChain->SetBranchAddress("L1_Mu12er2p3_Jet40er2p3_dR_Max0p4_DoubleJet40er2p3_dEta_Max1p6", &L1_Mu12er2p3_Jet40er2p3_dR_Max0p4_DoubleJet40er2p3_dEta_Max1p6, &b_L1_Mu12er2p3_Jet40er2p3_dR_Max0p4_DoubleJet40er2p3_dEta_Max1p6);
   fChain->SetBranchAddress("L1_Mu18er2p1_Tau24er2p1", &L1_Mu18er2p1_Tau24er2p1, &b_L1_Mu18er2p1_Tau24er2p1);
   fChain->SetBranchAddress("L1_Mu18er2p1_Tau26er2p1", &L1_Mu18er2p1_Tau26er2p1, &b_L1_Mu18er2p1_Tau26er2p1);
   fChain->SetBranchAddress("L1_Mu18er2p1_Tau26er2p1_Jet55", &L1_Mu18er2p1_Tau26er2p1_Jet55, &b_L1_Mu18er2p1_Tau26er2p1_Jet55);
   fChain->SetBranchAddress("L1_Mu18er2p1_Tau26er2p1_Jet70", &L1_Mu18er2p1_Tau26er2p1_Jet70, &b_L1_Mu18er2p1_Tau26er2p1_Jet70);
   fChain->SetBranchAddress("L1_Mu20_EG10er2p5", &L1_Mu20_EG10er2p5, &b_L1_Mu20_EG10er2p5);
   fChain->SetBranchAddress("L1_Mu22er2p1_IsoTau28er2p1", &L1_Mu22er2p1_IsoTau28er2p1, &b_L1_Mu22er2p1_IsoTau28er2p1);
   fChain->SetBranchAddress("L1_Mu22er2p1_IsoTau30er2p1", &L1_Mu22er2p1_IsoTau30er2p1, &b_L1_Mu22er2p1_IsoTau30er2p1);
   fChain->SetBranchAddress("L1_Mu22er2p1_IsoTau32er2p1", &L1_Mu22er2p1_IsoTau32er2p1, &b_L1_Mu22er2p1_IsoTau32er2p1);
   fChain->SetBranchAddress("L1_Mu22er2p1_IsoTau34er2p1", &L1_Mu22er2p1_IsoTau34er2p1, &b_L1_Mu22er2p1_IsoTau34er2p1);
   fChain->SetBranchAddress("L1_Mu22er2p1_IsoTau36er2p1", &L1_Mu22er2p1_IsoTau36er2p1, &b_L1_Mu22er2p1_IsoTau36er2p1);
   fChain->SetBranchAddress("L1_Mu22er2p1_IsoTau40er2p1", &L1_Mu22er2p1_IsoTau40er2p1, &b_L1_Mu22er2p1_IsoTau40er2p1);
   fChain->SetBranchAddress("L1_Mu22er2p1_Tau70er2p1", &L1_Mu22er2p1_Tau70er2p1, &b_L1_Mu22er2p1_Tau70er2p1);
   fChain->SetBranchAddress("L1_Mu3_Jet120er2p5_dR_Max0p4", &L1_Mu3_Jet120er2p5_dR_Max0p4, &b_L1_Mu3_Jet120er2p5_dR_Max0p4);
   fChain->SetBranchAddress("L1_Mu3_Jet120er2p5_dR_Max0p8", &L1_Mu3_Jet120er2p5_dR_Max0p8, &b_L1_Mu3_Jet120er2p5_dR_Max0p8);
   fChain->SetBranchAddress("L1_Mu3_Jet16er2p5_dR_Max0p4", &L1_Mu3_Jet16er2p5_dR_Max0p4, &b_L1_Mu3_Jet16er2p5_dR_Max0p4);
   fChain->SetBranchAddress("L1_Mu3_Jet30er2p5", &L1_Mu3_Jet30er2p5, &b_L1_Mu3_Jet30er2p5);
   fChain->SetBranchAddress("L1_Mu3_Jet35er2p5_dR_Max0p4", &L1_Mu3_Jet35er2p5_dR_Max0p4, &b_L1_Mu3_Jet35er2p5_dR_Max0p4);
   fChain->SetBranchAddress("L1_Mu3_Jet60er2p5_dR_Max0p4", &L1_Mu3_Jet60er2p5_dR_Max0p4, &b_L1_Mu3_Jet60er2p5_dR_Max0p4);
   fChain->SetBranchAddress("L1_Mu3_Jet80er2p5_dR_Max0p4", &L1_Mu3_Jet80er2p5_dR_Max0p4, &b_L1_Mu3_Jet80er2p5_dR_Max0p4);
   fChain->SetBranchAddress("L1_Mu3er1p5_Jet100er2p5_ETMHF40", &L1_Mu3er1p5_Jet100er2p5_ETMHF40, &b_L1_Mu3er1p5_Jet100er2p5_ETMHF40);
   fChain->SetBranchAddress("L1_Mu3er1p5_Jet100er2p5_ETMHF50", &L1_Mu3er1p5_Jet100er2p5_ETMHF50, &b_L1_Mu3er1p5_Jet100er2p5_ETMHF50);
   fChain->SetBranchAddress("L1_Mu5_EG23er2p5", &L1_Mu5_EG23er2p5, &b_L1_Mu5_EG23er2p5);
   fChain->SetBranchAddress("L1_Mu5_LooseIsoEG20er2p5", &L1_Mu5_LooseIsoEG20er2p5, &b_L1_Mu5_LooseIsoEG20er2p5);
   fChain->SetBranchAddress("L1_Mu6_DoubleEG10er2p5", &L1_Mu6_DoubleEG10er2p5, &b_L1_Mu6_DoubleEG10er2p5);
   fChain->SetBranchAddress("L1_Mu6_DoubleEG12er2p5", &L1_Mu6_DoubleEG12er2p5, &b_L1_Mu6_DoubleEG12er2p5);
   fChain->SetBranchAddress("L1_Mu6_DoubleEG15er2p5", &L1_Mu6_DoubleEG15er2p5, &b_L1_Mu6_DoubleEG15er2p5);
   fChain->SetBranchAddress("L1_Mu6_DoubleEG17er2p5", &L1_Mu6_DoubleEG17er2p5, &b_L1_Mu6_DoubleEG17er2p5);
   fChain->SetBranchAddress("L1_Mu6_HTT240er", &L1_Mu6_HTT240er, &b_L1_Mu6_HTT240er);
   fChain->SetBranchAddress("L1_Mu6_HTT250er", &L1_Mu6_HTT250er, &b_L1_Mu6_HTT250er);
   fChain->SetBranchAddress("L1_Mu7_EG20er2p5", &L1_Mu7_EG20er2p5, &b_L1_Mu7_EG20er2p5);
   fChain->SetBranchAddress("L1_Mu7_EG23er2p5", &L1_Mu7_EG23er2p5, &b_L1_Mu7_EG23er2p5);
   fChain->SetBranchAddress("L1_Mu7_LooseIsoEG20er2p5", &L1_Mu7_LooseIsoEG20er2p5, &b_L1_Mu7_LooseIsoEG20er2p5);
   fChain->SetBranchAddress("L1_Mu7_LooseIsoEG23er2p5", &L1_Mu7_LooseIsoEG23er2p5, &b_L1_Mu7_LooseIsoEG23er2p5);
   fChain->SetBranchAddress("L1_NotBptxOR", &L1_NotBptxOR, &b_L1_NotBptxOR);
   fChain->SetBranchAddress("L1_QuadJet60er2p5", &L1_QuadJet60er2p5, &b_L1_QuadJet60er2p5);
   fChain->SetBranchAddress("L1_QuadJet_95_75_65_20_DoubleJet_75_65_er2p5_Jet20_FWD3p0", &L1_QuadJet_95_75_65_20_DoubleJet_75_65_er2p5_Jet20_FWD3p0, &b_L1_QuadJet_95_75_65_20_DoubleJet_75_65_er2p5_Jet20_FWD3p0);
   fChain->SetBranchAddress("L1_QuadMu0", &L1_QuadMu0, &b_L1_QuadMu0);
   fChain->SetBranchAddress("L1_QuadMu0_OQ", &L1_QuadMu0_OQ, &b_L1_QuadMu0_OQ);
   fChain->SetBranchAddress("L1_QuadMu0_SQ", &L1_QuadMu0_SQ, &b_L1_QuadMu0_SQ);
   fChain->SetBranchAddress("L1_SecondBunchInTrain", &L1_SecondBunchInTrain, &b_L1_SecondBunchInTrain);
   fChain->SetBranchAddress("L1_SecondLastBunchInTrain", &L1_SecondLastBunchInTrain, &b_L1_SecondLastBunchInTrain);
   fChain->SetBranchAddress("L1_SingleEG10er2p5", &L1_SingleEG10er2p5, &b_L1_SingleEG10er2p5);
   fChain->SetBranchAddress("L1_SingleEG15er2p5", &L1_SingleEG15er2p5, &b_L1_SingleEG15er2p5);
   fChain->SetBranchAddress("L1_SingleEG26er2p5", &L1_SingleEG26er2p5, &b_L1_SingleEG26er2p5);
   fChain->SetBranchAddress("L1_SingleEG28_FWD2p5", &L1_SingleEG28_FWD2p5, &b_L1_SingleEG28_FWD2p5);
   fChain->SetBranchAddress("L1_SingleEG28er1p5", &L1_SingleEG28er1p5, &b_L1_SingleEG28er1p5);
   fChain->SetBranchAddress("L1_SingleEG28er2p1", &L1_SingleEG28er2p1, &b_L1_SingleEG28er2p1);
   fChain->SetBranchAddress("L1_SingleEG28er2p5", &L1_SingleEG28er2p5, &b_L1_SingleEG28er2p5);
   fChain->SetBranchAddress("L1_SingleEG34er2p5", &L1_SingleEG34er2p5, &b_L1_SingleEG34er2p5);
   fChain->SetBranchAddress("L1_SingleEG36er2p5", &L1_SingleEG36er2p5, &b_L1_SingleEG36er2p5);
   fChain->SetBranchAddress("L1_SingleEG38er2p5", &L1_SingleEG38er2p5, &b_L1_SingleEG38er2p5);
   fChain->SetBranchAddress("L1_SingleEG40er2p5", &L1_SingleEG40er2p5, &b_L1_SingleEG40er2p5);
   fChain->SetBranchAddress("L1_SingleEG42er2p5", &L1_SingleEG42er2p5, &b_L1_SingleEG42er2p5);
   fChain->SetBranchAddress("L1_SingleEG45er2p5", &L1_SingleEG45er2p5, &b_L1_SingleEG45er2p5);
   fChain->SetBranchAddress("L1_SingleEG50", &L1_SingleEG50, &b_L1_SingleEG50);
   fChain->SetBranchAddress("L1_SingleEG60", &L1_SingleEG60, &b_L1_SingleEG60);
   fChain->SetBranchAddress("L1_SingleEG8er2p5", &L1_SingleEG8er2p5, &b_L1_SingleEG8er2p5);
   fChain->SetBranchAddress("L1_SingleIsoEG24er1p5", &L1_SingleIsoEG24er1p5, &b_L1_SingleIsoEG24er1p5);
   fChain->SetBranchAddress("L1_SingleIsoEG24er2p1", &L1_SingleIsoEG24er2p1, &b_L1_SingleIsoEG24er2p1);
   fChain->SetBranchAddress("L1_SingleIsoEG26er1p5", &L1_SingleIsoEG26er1p5, &b_L1_SingleIsoEG26er1p5);
   fChain->SetBranchAddress("L1_SingleIsoEG26er2p1", &L1_SingleIsoEG26er2p1, &b_L1_SingleIsoEG26er2p1);
   fChain->SetBranchAddress("L1_SingleIsoEG26er2p5", &L1_SingleIsoEG26er2p5, &b_L1_SingleIsoEG26er2p5);
   fChain->SetBranchAddress("L1_SingleIsoEG28_FWD2p5", &L1_SingleIsoEG28_FWD2p5, &b_L1_SingleIsoEG28_FWD2p5);
   fChain->SetBranchAddress("L1_SingleIsoEG28er1p5", &L1_SingleIsoEG28er1p5, &b_L1_SingleIsoEG28er1p5);
   fChain->SetBranchAddress("L1_SingleIsoEG28er2p1", &L1_SingleIsoEG28er2p1, &b_L1_SingleIsoEG28er2p1);
   fChain->SetBranchAddress("L1_SingleIsoEG28er2p5", &L1_SingleIsoEG28er2p5, &b_L1_SingleIsoEG28er2p5);
   fChain->SetBranchAddress("L1_SingleIsoEG30er2p1", &L1_SingleIsoEG30er2p1, &b_L1_SingleIsoEG30er2p1);
   fChain->SetBranchAddress("L1_SingleIsoEG30er2p5", &L1_SingleIsoEG30er2p5, &b_L1_SingleIsoEG30er2p5);
   fChain->SetBranchAddress("L1_SingleIsoEG32er2p1", &L1_SingleIsoEG32er2p1, &b_L1_SingleIsoEG32er2p1);
   fChain->SetBranchAddress("L1_SingleIsoEG32er2p5", &L1_SingleIsoEG32er2p5, &b_L1_SingleIsoEG32er2p5);
   fChain->SetBranchAddress("L1_SingleIsoEG34er2p5", &L1_SingleIsoEG34er2p5, &b_L1_SingleIsoEG34er2p5);
   fChain->SetBranchAddress("L1_SingleIsoTau32er2p1", &L1_SingleIsoTau32er2p1, &b_L1_SingleIsoTau32er2p1);
   fChain->SetBranchAddress("L1_SingleJet10erHE", &L1_SingleJet10erHE, &b_L1_SingleJet10erHE);
   fChain->SetBranchAddress("L1_SingleJet120", &L1_SingleJet120, &b_L1_SingleJet120);
   fChain->SetBranchAddress("L1_SingleJet120_FWD3p0", &L1_SingleJet120_FWD3p0, &b_L1_SingleJet120_FWD3p0);
   fChain->SetBranchAddress("L1_SingleJet120er2p5", &L1_SingleJet120er2p5, &b_L1_SingleJet120er2p5);
   fChain->SetBranchAddress("L1_SingleJet12erHE", &L1_SingleJet12erHE, &b_L1_SingleJet12erHE);
   fChain->SetBranchAddress("L1_SingleJet140er2p5", &L1_SingleJet140er2p5, &b_L1_SingleJet140er2p5);
   fChain->SetBranchAddress("L1_SingleJet140er2p5_ETMHF70", &L1_SingleJet140er2p5_ETMHF70, &b_L1_SingleJet140er2p5_ETMHF70);
   fChain->SetBranchAddress("L1_SingleJet140er2p5_ETMHF80", &L1_SingleJet140er2p5_ETMHF80, &b_L1_SingleJet140er2p5_ETMHF80);
   fChain->SetBranchAddress("L1_SingleJet140er2p5_ETMHF90", &L1_SingleJet140er2p5_ETMHF90, &b_L1_SingleJet140er2p5_ETMHF90);
   fChain->SetBranchAddress("L1_SingleJet160er2p5", &L1_SingleJet160er2p5, &b_L1_SingleJet160er2p5);
   fChain->SetBranchAddress("L1_SingleJet180", &L1_SingleJet180, &b_L1_SingleJet180);
   fChain->SetBranchAddress("L1_SingleJet180er2p5", &L1_SingleJet180er2p5, &b_L1_SingleJet180er2p5);
   fChain->SetBranchAddress("L1_SingleJet200", &L1_SingleJet200, &b_L1_SingleJet200);
   fChain->SetBranchAddress("L1_SingleJet20er2p5_NotBptxOR", &L1_SingleJet20er2p5_NotBptxOR, &b_L1_SingleJet20er2p5_NotBptxOR);
   fChain->SetBranchAddress("L1_SingleJet20er2p5_NotBptxOR_3BX", &L1_SingleJet20er2p5_NotBptxOR_3BX, &b_L1_SingleJet20er2p5_NotBptxOR_3BX);
   fChain->SetBranchAddress("L1_SingleJet35", &L1_SingleJet35, &b_L1_SingleJet35);
   fChain->SetBranchAddress("L1_SingleJet35_FWD3p0", &L1_SingleJet35_FWD3p0, &b_L1_SingleJet35_FWD3p0);
   fChain->SetBranchAddress("L1_SingleJet35er2p5", &L1_SingleJet35er2p5, &b_L1_SingleJet35er2p5);
   fChain->SetBranchAddress("L1_SingleJet43er2p5_NotBptxOR_3BX", &L1_SingleJet43er2p5_NotBptxOR_3BX, &b_L1_SingleJet43er2p5_NotBptxOR_3BX);
   fChain->SetBranchAddress("L1_SingleJet46er2p5_NotBptxOR_3BX", &L1_SingleJet46er2p5_NotBptxOR_3BX, &b_L1_SingleJet46er2p5_NotBptxOR_3BX);
   fChain->SetBranchAddress("L1_SingleJet60", &L1_SingleJet60, &b_L1_SingleJet60);
   fChain->SetBranchAddress("L1_SingleJet60_FWD3p0", &L1_SingleJet60_FWD3p0, &b_L1_SingleJet60_FWD3p0);
   fChain->SetBranchAddress("L1_SingleJet60er2p5", &L1_SingleJet60er2p5, &b_L1_SingleJet60er2p5);
   fChain->SetBranchAddress("L1_SingleJet8erHE", &L1_SingleJet8erHE, &b_L1_SingleJet8erHE);
   fChain->SetBranchAddress("L1_SingleJet90", &L1_SingleJet90, &b_L1_SingleJet90);
   fChain->SetBranchAddress("L1_SingleJet90_FWD3p0", &L1_SingleJet90_FWD3p0, &b_L1_SingleJet90_FWD3p0);
   fChain->SetBranchAddress("L1_SingleJet90er2p5", &L1_SingleJet90er2p5, &b_L1_SingleJet90er2p5);
   fChain->SetBranchAddress("L1_SingleLooseIsoEG26er1p5", &L1_SingleLooseIsoEG26er1p5, &b_L1_SingleLooseIsoEG26er1p5);
   fChain->SetBranchAddress("L1_SingleLooseIsoEG26er2p5", &L1_SingleLooseIsoEG26er2p5, &b_L1_SingleLooseIsoEG26er2p5);
   fChain->SetBranchAddress("L1_SingleLooseIsoEG28_FWD2p5", &L1_SingleLooseIsoEG28_FWD2p5, &b_L1_SingleLooseIsoEG28_FWD2p5);
   fChain->SetBranchAddress("L1_SingleLooseIsoEG28er1p5", &L1_SingleLooseIsoEG28er1p5, &b_L1_SingleLooseIsoEG28er1p5);
   fChain->SetBranchAddress("L1_SingleLooseIsoEG28er2p1", &L1_SingleLooseIsoEG28er2p1, &b_L1_SingleLooseIsoEG28er2p1);
   fChain->SetBranchAddress("L1_SingleLooseIsoEG28er2p5", &L1_SingleLooseIsoEG28er2p5, &b_L1_SingleLooseIsoEG28er2p5);
   fChain->SetBranchAddress("L1_SingleLooseIsoEG30er1p5", &L1_SingleLooseIsoEG30er1p5, &b_L1_SingleLooseIsoEG30er1p5);
   fChain->SetBranchAddress("L1_SingleLooseIsoEG30er2p5", &L1_SingleLooseIsoEG30er2p5, &b_L1_SingleLooseIsoEG30er2p5);
   fChain->SetBranchAddress("L1_SingleMu0_BMTF", &L1_SingleMu0_BMTF, &b_L1_SingleMu0_BMTF);
   fChain->SetBranchAddress("L1_SingleMu0_DQ", &L1_SingleMu0_DQ, &b_L1_SingleMu0_DQ);
   fChain->SetBranchAddress("L1_SingleMu0_EMTF", &L1_SingleMu0_EMTF, &b_L1_SingleMu0_EMTF);
   fChain->SetBranchAddress("L1_SingleMu0_OMTF", &L1_SingleMu0_OMTF, &b_L1_SingleMu0_OMTF);
   fChain->SetBranchAddress("L1_SingleMu10er1p5", &L1_SingleMu10er1p5, &b_L1_SingleMu10er1p5);
   fChain->SetBranchAddress("L1_SingleMu12_DQ_BMTF", &L1_SingleMu12_DQ_BMTF, &b_L1_SingleMu12_DQ_BMTF);
   fChain->SetBranchAddress("L1_SingleMu12_DQ_EMTF", &L1_SingleMu12_DQ_EMTF, &b_L1_SingleMu12_DQ_EMTF);
   fChain->SetBranchAddress("L1_SingleMu12_DQ_OMTF", &L1_SingleMu12_DQ_OMTF, &b_L1_SingleMu12_DQ_OMTF);
   fChain->SetBranchAddress("L1_SingleMu12er1p5", &L1_SingleMu12er1p5, &b_L1_SingleMu12er1p5);
   fChain->SetBranchAddress("L1_SingleMu14er1p5", &L1_SingleMu14er1p5, &b_L1_SingleMu14er1p5);
   fChain->SetBranchAddress("L1_SingleMu15_DQ", &L1_SingleMu15_DQ, &b_L1_SingleMu15_DQ);
   fChain->SetBranchAddress("L1_SingleMu16er1p5", &L1_SingleMu16er1p5, &b_L1_SingleMu16er1p5);
   fChain->SetBranchAddress("L1_SingleMu18", &L1_SingleMu18, &b_L1_SingleMu18);
   fChain->SetBranchAddress("L1_SingleMu18er1p5", &L1_SingleMu18er1p5, &b_L1_SingleMu18er1p5);
   fChain->SetBranchAddress("L1_SingleMu20", &L1_SingleMu20, &b_L1_SingleMu20);
   fChain->SetBranchAddress("L1_SingleMu22", &L1_SingleMu22, &b_L1_SingleMu22);
   fChain->SetBranchAddress("L1_SingleMu22_BMTF", &L1_SingleMu22_BMTF, &b_L1_SingleMu22_BMTF);
   fChain->SetBranchAddress("L1_SingleMu22_DQ", &L1_SingleMu22_DQ, &b_L1_SingleMu22_DQ);
   fChain->SetBranchAddress("L1_SingleMu22_EMTF", &L1_SingleMu22_EMTF, &b_L1_SingleMu22_EMTF);
   fChain->SetBranchAddress("L1_SingleMu22_OMTF", &L1_SingleMu22_OMTF, &b_L1_SingleMu22_OMTF);
   fChain->SetBranchAddress("L1_SingleMu22_OQ", &L1_SingleMu22_OQ, &b_L1_SingleMu22_OQ);
   fChain->SetBranchAddress("L1_SingleMu25", &L1_SingleMu25, &b_L1_SingleMu25);
   fChain->SetBranchAddress("L1_SingleMu3", &L1_SingleMu3, &b_L1_SingleMu3);
   fChain->SetBranchAddress("L1_SingleMu5", &L1_SingleMu5, &b_L1_SingleMu5);
   fChain->SetBranchAddress("L1_SingleMu6er1p5", &L1_SingleMu6er1p5, &b_L1_SingleMu6er1p5);
   fChain->SetBranchAddress("L1_SingleMu7", &L1_SingleMu7, &b_L1_SingleMu7);
   fChain->SetBranchAddress("L1_SingleMu7_DQ", &L1_SingleMu7_DQ, &b_L1_SingleMu7_DQ);
   fChain->SetBranchAddress("L1_SingleMu7er1p5", &L1_SingleMu7er1p5, &b_L1_SingleMu7er1p5);
   fChain->SetBranchAddress("L1_SingleMu8er1p5", &L1_SingleMu8er1p5, &b_L1_SingleMu8er1p5);
   fChain->SetBranchAddress("L1_SingleMu9er1p5", &L1_SingleMu9er1p5, &b_L1_SingleMu9er1p5);
   fChain->SetBranchAddress("L1_SingleMuCosmics", &L1_SingleMuCosmics, &b_L1_SingleMuCosmics);
   fChain->SetBranchAddress("L1_SingleMuCosmics_BMTF", &L1_SingleMuCosmics_BMTF, &b_L1_SingleMuCosmics_BMTF);
   fChain->SetBranchAddress("L1_SingleMuCosmics_EMTF", &L1_SingleMuCosmics_EMTF, &b_L1_SingleMuCosmics_EMTF);
   fChain->SetBranchAddress("L1_SingleMuCosmics_OMTF", &L1_SingleMuCosmics_OMTF, &b_L1_SingleMuCosmics_OMTF);
   fChain->SetBranchAddress("L1_SingleMuOpen", &L1_SingleMuOpen, &b_L1_SingleMuOpen);
   fChain->SetBranchAddress("L1_SingleMuOpen_NotBptxOR", &L1_SingleMuOpen_NotBptxOR, &b_L1_SingleMuOpen_NotBptxOR);
   fChain->SetBranchAddress("L1_SingleMuOpen_er1p1_NotBptxOR_3BX", &L1_SingleMuOpen_er1p1_NotBptxOR_3BX, &b_L1_SingleMuOpen_er1p1_NotBptxOR_3BX);
   fChain->SetBranchAddress("L1_SingleMuOpen_er1p4_NotBptxOR_3BX", &L1_SingleMuOpen_er1p4_NotBptxOR_3BX, &b_L1_SingleMuOpen_er1p4_NotBptxOR_3BX);
   fChain->SetBranchAddress("L1_SingleMuShower_Nominal", &L1_SingleMuShower_Nominal, &b_L1_SingleMuShower_Nominal);
   fChain->SetBranchAddress("L1_SingleMuShower_Tight", &L1_SingleMuShower_Tight, &b_L1_SingleMuShower_Tight);
   fChain->SetBranchAddress("L1_SingleTau120er2p1", &L1_SingleTau120er2p1, &b_L1_SingleTau120er2p1);
   fChain->SetBranchAddress("L1_SingleTau130er2p1", &L1_SingleTau130er2p1, &b_L1_SingleTau130er2p1);
   fChain->SetBranchAddress("L1_SingleTau70er2p1", &L1_SingleTau70er2p1, &b_L1_SingleTau70er2p1);
   fChain->SetBranchAddress("L1_TOTEM_1", &L1_TOTEM_1, &b_L1_TOTEM_1);
   fChain->SetBranchAddress("L1_TOTEM_2", &L1_TOTEM_2, &b_L1_TOTEM_2);
   fChain->SetBranchAddress("L1_TOTEM_3", &L1_TOTEM_3, &b_L1_TOTEM_3);
   fChain->SetBranchAddress("L1_TOTEM_4", &L1_TOTEM_4, &b_L1_TOTEM_4);
   fChain->SetBranchAddress("L1_TripleEG16er2p5", &L1_TripleEG16er2p5, &b_L1_TripleEG16er2p5);
   fChain->SetBranchAddress("L1_TripleEG_16_12_8_er2p5", &L1_TripleEG_16_12_8_er2p5, &b_L1_TripleEG_16_12_8_er2p5);
   fChain->SetBranchAddress("L1_TripleEG_16_15_8_er2p5", &L1_TripleEG_16_15_8_er2p5, &b_L1_TripleEG_16_15_8_er2p5);
   fChain->SetBranchAddress("L1_TripleEG_18_17_8_er2p5", &L1_TripleEG_18_17_8_er2p5, &b_L1_TripleEG_18_17_8_er2p5);
   fChain->SetBranchAddress("L1_TripleEG_18_18_12_er2p5", &L1_TripleEG_18_18_12_er2p5, &b_L1_TripleEG_18_18_12_er2p5);
   fChain->SetBranchAddress("L1_TripleJet_100_80_70_DoubleJet_80_70_er2p5", &L1_TripleJet_100_80_70_DoubleJet_80_70_er2p5, &b_L1_TripleJet_100_80_70_DoubleJet_80_70_er2p5);
   fChain->SetBranchAddress("L1_TripleJet_105_85_75_DoubleJet_85_75_er2p5", &L1_TripleJet_105_85_75_DoubleJet_85_75_er2p5, &b_L1_TripleJet_105_85_75_DoubleJet_85_75_er2p5);
   fChain->SetBranchAddress("L1_TripleJet_95_75_65_DoubleJet_75_65_er2p5", &L1_TripleJet_95_75_65_DoubleJet_75_65_er2p5, &b_L1_TripleJet_95_75_65_DoubleJet_75_65_er2p5);
   fChain->SetBranchAddress("L1_TripleMu0", &L1_TripleMu0, &b_L1_TripleMu0);
   fChain->SetBranchAddress("L1_TripleMu0_OQ", &L1_TripleMu0_OQ, &b_L1_TripleMu0_OQ);
   fChain->SetBranchAddress("L1_TripleMu0_SQ", &L1_TripleMu0_SQ, &b_L1_TripleMu0_SQ);
   fChain->SetBranchAddress("L1_TripleMu3", &L1_TripleMu3, &b_L1_TripleMu3);
   fChain->SetBranchAddress("L1_TripleMu3_SQ", &L1_TripleMu3_SQ, &b_L1_TripleMu3_SQ);
   fChain->SetBranchAddress("L1_TripleMu_2SQ_1p5SQ_0OQ", &L1_TripleMu_2SQ_1p5SQ_0OQ, &b_L1_TripleMu_2SQ_1p5SQ_0OQ);
   fChain->SetBranchAddress("L1_TripleMu_2SQ_1p5SQ_0OQ_Mass_Max12", &L1_TripleMu_2SQ_1p5SQ_0OQ_Mass_Max12, &b_L1_TripleMu_2SQ_1p5SQ_0OQ_Mass_Max12);
   fChain->SetBranchAddress("L1_TripleMu_3SQ_2p5SQ_0OQ_Mass_Max12", &L1_TripleMu_3SQ_2p5SQ_0OQ_Mass_Max12, &b_L1_TripleMu_3SQ_2p5SQ_0OQ_Mass_Max12);
   fChain->SetBranchAddress("L1_TripleMu_5SQ_3SQ_0OQ", &L1_TripleMu_5SQ_3SQ_0OQ, &b_L1_TripleMu_5SQ_3SQ_0OQ);
   fChain->SetBranchAddress("L1_TripleMu_5SQ_3SQ_0OQ_DoubleMu_5_3_SQ_OS_Mass_Max9", &L1_TripleMu_5SQ_3SQ_0OQ_DoubleMu_5_3_SQ_OS_Mass_Max9, &b_L1_TripleMu_5SQ_3SQ_0OQ_DoubleMu_5_3_SQ_OS_Mass_Max9);
   fChain->SetBranchAddress("L1_TripleMu_5SQ_3SQ_0_DoubleMu_5_3_SQ_OS_Mass_Max9", &L1_TripleMu_5SQ_3SQ_0_DoubleMu_5_3_SQ_OS_Mass_Max9, &b_L1_TripleMu_5SQ_3SQ_0_DoubleMu_5_3_SQ_OS_Mass_Max9);
   fChain->SetBranchAddress("L1_TripleMu_5_3_3", &L1_TripleMu_5_3_3, &b_L1_TripleMu_5_3_3);
   fChain->SetBranchAddress("L1_TripleMu_5_3_3_SQ", &L1_TripleMu_5_3_3_SQ, &b_L1_TripleMu_5_3_3_SQ);
   fChain->SetBranchAddress("L1_TripleMu_5_3p5_2p5", &L1_TripleMu_5_3p5_2p5, &b_L1_TripleMu_5_3p5_2p5);
   fChain->SetBranchAddress("L1_TripleMu_5_3p5_2p5_DoubleMu_5_2p5_OS_Mass_5to17", &L1_TripleMu_5_3p5_2p5_DoubleMu_5_2p5_OS_Mass_5to17, &b_L1_TripleMu_5_3p5_2p5_DoubleMu_5_2p5_OS_Mass_5to17);
   fChain->SetBranchAddress("L1_TripleMu_5_3p5_2p5_OQ_DoubleMu_5_2p5_OQ_OS_Mass_5to17", &L1_TripleMu_5_3p5_2p5_OQ_DoubleMu_5_2p5_OQ_OS_Mass_5to17, &b_L1_TripleMu_5_3p5_2p5_OQ_DoubleMu_5_2p5_OQ_OS_Mass_5to17);
   fChain->SetBranchAddress("L1_TripleMu_5_4_2p5_DoubleMu_5_2p5_OS_Mass_5to17", &L1_TripleMu_5_4_2p5_DoubleMu_5_2p5_OS_Mass_5to17, &b_L1_TripleMu_5_4_2p5_DoubleMu_5_2p5_OS_Mass_5to17);
   fChain->SetBranchAddress("L1_TripleMu_5_5_3", &L1_TripleMu_5_5_3, &b_L1_TripleMu_5_5_3);
   fChain->SetBranchAddress("L1_UnpairedBunchBptxMinus", &L1_UnpairedBunchBptxMinus, &b_L1_UnpairedBunchBptxMinus);
   fChain->SetBranchAddress("L1_UnpairedBunchBptxPlus", &L1_UnpairedBunchBptxPlus, &b_L1_UnpairedBunchBptxPlus);
   fChain->SetBranchAddress("L1_ZeroBias", &L1_ZeroBias, &b_L1_ZeroBias);
   fChain->SetBranchAddress("L1_ZeroBias_copy", &L1_ZeroBias_copy, &b_L1_ZeroBias_copy);
   fChain->SetBranchAddress("L1_UnprefireableEvent", &L1_UnprefireableEvent, &b_L1_UnprefireableEvent);
   fChain->SetBranchAddress("L1Reco_step", &L1Reco_step, &b_L1Reco_step);
   fChain->SetBranchAddress("Flag_HBHENoiseFilter_pRECO", &Flag_HBHENoiseFilter_pRECO, &b_Flag_HBHENoiseFilter_pRECO);
   fChain->SetBranchAddress("Flag_HBHENoiseIsoFilter_pRECO", &Flag_HBHENoiseIsoFilter_pRECO, &b_Flag_HBHENoiseIsoFilter_pRECO);
   fChain->SetBranchAddress("Flag_CSCTightHaloFilter_pRECO", &Flag_CSCTightHaloFilter_pRECO, &b_Flag_CSCTightHaloFilter_pRECO);
   fChain->SetBranchAddress("Flag_CSCTightHaloTrkMuUnvetoFilter_pRECO", &Flag_CSCTightHaloTrkMuUnvetoFilter_pRECO, &b_Flag_CSCTightHaloTrkMuUnvetoFilter_pRECO);
   fChain->SetBranchAddress("Flag_CSCTightHalo2015Filter_pRECO", &Flag_CSCTightHalo2015Filter_pRECO, &b_Flag_CSCTightHalo2015Filter_pRECO);
   fChain->SetBranchAddress("Flag_globalTightHalo2016Filter_pRECO", &Flag_globalTightHalo2016Filter_pRECO, &b_Flag_globalTightHalo2016Filter_pRECO);
   fChain->SetBranchAddress("Flag_globalSuperTightHalo2016Filter_pRECO", &Flag_globalSuperTightHalo2016Filter_pRECO, &b_Flag_globalSuperTightHalo2016Filter_pRECO);
   fChain->SetBranchAddress("Flag_HcalStripHaloFilter_pRECO", &Flag_HcalStripHaloFilter_pRECO, &b_Flag_HcalStripHaloFilter_pRECO);
   fChain->SetBranchAddress("Flag_hcalLaserEventFilter_pRECO", &Flag_hcalLaserEventFilter_pRECO, &b_Flag_hcalLaserEventFilter_pRECO);
   fChain->SetBranchAddress("Flag_EcalDeadCellTriggerPrimitiveFilter_pRECO", &Flag_EcalDeadCellTriggerPrimitiveFilter_pRECO, &b_Flag_EcalDeadCellTriggerPrimitiveFilter_pRECO);
   fChain->SetBranchAddress("Flag_EcalDeadCellBoundaryEnergyFilter_pRECO", &Flag_EcalDeadCellBoundaryEnergyFilter_pRECO, &b_Flag_EcalDeadCellBoundaryEnergyFilter_pRECO);
   fChain->SetBranchAddress("Flag_ecalBadCalibFilter_pRECO", &Flag_ecalBadCalibFilter_pRECO, &b_Flag_ecalBadCalibFilter_pRECO);
   fChain->SetBranchAddress("Flag_goodVertices_pRECO", &Flag_goodVertices_pRECO, &b_Flag_goodVertices_pRECO);
   fChain->SetBranchAddress("Flag_eeBadScFilter_pRECO", &Flag_eeBadScFilter_pRECO, &b_Flag_eeBadScFilter_pRECO);
   fChain->SetBranchAddress("Flag_ecalLaserCorrFilter_pRECO", &Flag_ecalLaserCorrFilter_pRECO, &b_Flag_ecalLaserCorrFilter_pRECO);
   fChain->SetBranchAddress("Flag_trkPOGFilters_pRECO", &Flag_trkPOGFilters_pRECO, &b_Flag_trkPOGFilters_pRECO);
   fChain->SetBranchAddress("Flag_chargedHadronTrackResolutionFilter_pRECO", &Flag_chargedHadronTrackResolutionFilter_pRECO, &b_Flag_chargedHadronTrackResolutionFilter_pRECO);
   fChain->SetBranchAddress("Flag_muonBadTrackFilter_pRECO", &Flag_muonBadTrackFilter_pRECO, &b_Flag_muonBadTrackFilter_pRECO);
   fChain->SetBranchAddress("Flag_BadChargedCandidateFilter_pRECO", &Flag_BadChargedCandidateFilter_pRECO, &b_Flag_BadChargedCandidateFilter_pRECO);
   fChain->SetBranchAddress("Flag_BadPFMuonFilter_pRECO", &Flag_BadPFMuonFilter_pRECO, &b_Flag_BadPFMuonFilter_pRECO);
   fChain->SetBranchAddress("Flag_BadPFMuonDzFilter_pRECO", &Flag_BadPFMuonDzFilter_pRECO, &b_Flag_BadPFMuonDzFilter_pRECO);
   fChain->SetBranchAddress("Flag_hfNoisyHitsFilter_pRECO", &Flag_hfNoisyHitsFilter_pRECO, &b_Flag_hfNoisyHitsFilter_pRECO);
   fChain->SetBranchAddress("Flag_BadChargedCandidateSummer16Filter_pRECO", &Flag_BadChargedCandidateSummer16Filter_pRECO, &b_Flag_BadChargedCandidateSummer16Filter_pRECO);
   fChain->SetBranchAddress("Flag_BadPFMuonSummer16Filter_pRECO", &Flag_BadPFMuonSummer16Filter_pRECO, &b_Flag_BadPFMuonSummer16Filter_pRECO);
   fChain->SetBranchAddress("Flag_trkPOG_manystripclus53X_pRECO", &Flag_trkPOG_manystripclus53X_pRECO, &b_Flag_trkPOG_manystripclus53X_pRECO);
   fChain->SetBranchAddress("Flag_trkPOG_toomanystripclus53X_pRECO", &Flag_trkPOG_toomanystripclus53X_pRECO, &b_Flag_trkPOG_toomanystripclus53X_pRECO);
   fChain->SetBranchAddress("Flag_trkPOG_logErrorTooManyClusters_pRECO", &Flag_trkPOG_logErrorTooManyClusters_pRECO, &b_Flag_trkPOG_logErrorTooManyClusters_pRECO);
   fChain->SetBranchAddress("Flag_METFilters_pRECO", &Flag_METFilters_pRECO, &b_Flag_METFilters_pRECO);
   fChain->SetBranchAddress("HLTriggerFirstPath", &HLTriggerFirstPath, &b_HLTriggerFirstPath);
   fChain->SetBranchAddress("HLT_AK8PFJet360_TrimMass30", &HLT_AK8PFJet360_TrimMass30, &b_HLT_AK8PFJet360_TrimMass30);
   fChain->SetBranchAddress("HLT_AK8PFJet380_TrimMass30", &HLT_AK8PFJet380_TrimMass30, &b_HLT_AK8PFJet380_TrimMass30);
   fChain->SetBranchAddress("HLT_AK8PFJet400_TrimMass30", &HLT_AK8PFJet400_TrimMass30, &b_HLT_AK8PFJet400_TrimMass30);
   fChain->SetBranchAddress("HLT_AK8PFJet420_TrimMass30", &HLT_AK8PFJet420_TrimMass30, &b_HLT_AK8PFJet420_TrimMass30);
   fChain->SetBranchAddress("HLT_AK8PFJet400_MassSD30", &HLT_AK8PFJet400_MassSD30, &b_HLT_AK8PFJet400_MassSD30);
   fChain->SetBranchAddress("HLT_AK8PFJet420_MassSD30", &HLT_AK8PFJet420_MassSD30, &b_HLT_AK8PFJet420_MassSD30);
   fChain->SetBranchAddress("HLT_AK8PFJet450_MassSD30", &HLT_AK8PFJet450_MassSD30, &b_HLT_AK8PFJet450_MassSD30);
   fChain->SetBranchAddress("HLT_AK8DiPFJet250_250_MassSD30", &HLT_AK8DiPFJet250_250_MassSD30, &b_HLT_AK8DiPFJet250_250_MassSD30);
   fChain->SetBranchAddress("HLT_AK8DiPFJet250_250_MassSD50", &HLT_AK8DiPFJet250_250_MassSD50, &b_HLT_AK8DiPFJet250_250_MassSD50);
   fChain->SetBranchAddress("HLT_AK8DiPFJet260_260_MassSD30", &HLT_AK8DiPFJet260_260_MassSD30, &b_HLT_AK8DiPFJet260_260_MassSD30);
   fChain->SetBranchAddress("HLT_AK8DiPFJet270_270_MassSD30", &HLT_AK8DiPFJet270_270_MassSD30, &b_HLT_AK8DiPFJet270_270_MassSD30);
   fChain->SetBranchAddress("HLT_AK8PFHT750_TrimMass50", &HLT_AK8PFHT750_TrimMass50, &b_HLT_AK8PFHT750_TrimMass50);
   fChain->SetBranchAddress("HLT_AK8PFHT800_TrimMass50", &HLT_AK8PFHT800_TrimMass50, &b_HLT_AK8PFHT800_TrimMass50);
   fChain->SetBranchAddress("HLT_AK8PFHT850_TrimMass50", &HLT_AK8PFHT850_TrimMass50, &b_HLT_AK8PFHT850_TrimMass50);
   fChain->SetBranchAddress("HLT_AK8PFHT900_TrimMass50", &HLT_AK8PFHT900_TrimMass50, &b_HLT_AK8PFHT900_TrimMass50);
   fChain->SetBranchAddress("HLT_CaloJet500_NoJetID", &HLT_CaloJet500_NoJetID, &b_HLT_CaloJet500_NoJetID);
   fChain->SetBranchAddress("HLT_CaloJet550_NoJetID", &HLT_CaloJet550_NoJetID, &b_HLT_CaloJet550_NoJetID);
   fChain->SetBranchAddress("HLT_DoubleMu5_Upsilon_DoubleEle3_CaloIdL_TrackIdL", &HLT_DoubleMu5_Upsilon_DoubleEle3_CaloIdL_TrackIdL, &b_HLT_DoubleMu5_Upsilon_DoubleEle3_CaloIdL_TrackIdL);
   fChain->SetBranchAddress("HLT_DoubleMu3_DoubleEle7p5_CaloIdL_TrackIdL_Upsilon", &HLT_DoubleMu3_DoubleEle7p5_CaloIdL_TrackIdL_Upsilon, &b_HLT_DoubleMu3_DoubleEle7p5_CaloIdL_TrackIdL_Upsilon);
   fChain->SetBranchAddress("HLT_Trimuon5_3p5_2_Upsilon_Muon", &HLT_Trimuon5_3p5_2_Upsilon_Muon, &b_HLT_Trimuon5_3p5_2_Upsilon_Muon);
   fChain->SetBranchAddress("HLT_TrimuonOpen_5_3p5_2_Upsilon_Muon", &HLT_TrimuonOpen_5_3p5_2_Upsilon_Muon, &b_HLT_TrimuonOpen_5_3p5_2_Upsilon_Muon);
   fChain->SetBranchAddress("HLT_DoubleEle25_CaloIdL_MW", &HLT_DoubleEle25_CaloIdL_MW, &b_HLT_DoubleEle25_CaloIdL_MW);
   fChain->SetBranchAddress("HLT_DoubleEle27_CaloIdL_MW", &HLT_DoubleEle27_CaloIdL_MW, &b_HLT_DoubleEle27_CaloIdL_MW);
   fChain->SetBranchAddress("HLT_DoubleEle33_CaloIdL_MW", &HLT_DoubleEle33_CaloIdL_MW, &b_HLT_DoubleEle33_CaloIdL_MW);
   fChain->SetBranchAddress("HLT_DoubleEle24_eta2p1_WPTight_Gsf", &HLT_DoubleEle24_eta2p1_WPTight_Gsf, &b_HLT_DoubleEle24_eta2p1_WPTight_Gsf);
   fChain->SetBranchAddress("HLT_DoubleEle8_CaloIdM_TrackIdM_Mass8_DZ_PFHT350", &HLT_DoubleEle8_CaloIdM_TrackIdM_Mass8_DZ_PFHT350, &b_HLT_DoubleEle8_CaloIdM_TrackIdM_Mass8_DZ_PFHT350);
   fChain->SetBranchAddress("HLT_DoubleEle8_CaloIdM_TrackIdM_Mass8_PFHT350", &HLT_DoubleEle8_CaloIdM_TrackIdM_Mass8_PFHT350, &b_HLT_DoubleEle8_CaloIdM_TrackIdM_Mass8_PFHT350);
   fChain->SetBranchAddress("HLT_Mu27_Ele37_CaloIdL_MW", &HLT_Mu27_Ele37_CaloIdL_MW, &b_HLT_Mu27_Ele37_CaloIdL_MW);
   fChain->SetBranchAddress("HLT_Mu37_Ele27_CaloIdL_MW", &HLT_Mu37_Ele27_CaloIdL_MW, &b_HLT_Mu37_Ele27_CaloIdL_MW);
   fChain->SetBranchAddress("HLT_Mu37_TkMu27", &HLT_Mu37_TkMu27, &b_HLT_Mu37_TkMu27);
   fChain->SetBranchAddress("HLT_DoubleMu4_3_Bs", &HLT_DoubleMu4_3_Bs, &b_HLT_DoubleMu4_3_Bs);
   fChain->SetBranchAddress("HLT_DoubleMu4_3_Jpsi", &HLT_DoubleMu4_3_Jpsi, &b_HLT_DoubleMu4_3_Jpsi);
   fChain->SetBranchAddress("HLT_DoubleMu4_3_LowMass", &HLT_DoubleMu4_3_LowMass, &b_HLT_DoubleMu4_3_LowMass);
   fChain->SetBranchAddress("HLT_DoubleMu4_LowMass_Displaced", &HLT_DoubleMu4_LowMass_Displaced, &b_HLT_DoubleMu4_LowMass_Displaced);
   fChain->SetBranchAddress("HLT_Mu0_L1DoubleMu", &HLT_Mu0_L1DoubleMu, &b_HLT_Mu0_L1DoubleMu);
   fChain->SetBranchAddress("HLT_Mu4_L1DoubleMu", &HLT_Mu4_L1DoubleMu, &b_HLT_Mu4_L1DoubleMu);
   fChain->SetBranchAddress("HLT_DoubleMu4_3_Photon4_BsToMMG", &HLT_DoubleMu4_3_Photon4_BsToMMG, &b_HLT_DoubleMu4_3_Photon4_BsToMMG);
   fChain->SetBranchAddress("HLT_DoubleMu4_3_Displaced_Photon4_BsToMMG", &HLT_DoubleMu4_3_Displaced_Photon4_BsToMMG, &b_HLT_DoubleMu4_3_Displaced_Photon4_BsToMMG);
   fChain->SetBranchAddress("HLT_DoubleMu3_Trk_Tau3mu", &HLT_DoubleMu3_Trk_Tau3mu, &b_HLT_DoubleMu3_Trk_Tau3mu);
   fChain->SetBranchAddress("HLT_DoubleMu3_TkMu_DsTau3Mu", &HLT_DoubleMu3_TkMu_DsTau3Mu, &b_HLT_DoubleMu3_TkMu_DsTau3Mu);
   fChain->SetBranchAddress("HLT_DoubleMu4_Mass3p8_DZ_PFHT350", &HLT_DoubleMu4_Mass3p8_DZ_PFHT350, &b_HLT_DoubleMu4_Mass3p8_DZ_PFHT350);
   fChain->SetBranchAddress("HLT_DoubleMu4_MuMuTrk_Displaced", &HLT_DoubleMu4_MuMuTrk_Displaced, &b_HLT_DoubleMu4_MuMuTrk_Displaced);
   fChain->SetBranchAddress("HLT_Mu3_PFJet40", &HLT_Mu3_PFJet40, &b_HLT_Mu3_PFJet40);
   fChain->SetBranchAddress("HLT_Mu7p5_L2Mu2_Jpsi", &HLT_Mu7p5_L2Mu2_Jpsi, &b_HLT_Mu7p5_L2Mu2_Jpsi);
   fChain->SetBranchAddress("HLT_Mu7p5_L2Mu2_Upsilon", &HLT_Mu7p5_L2Mu2_Upsilon, &b_HLT_Mu7p5_L2Mu2_Upsilon);
   fChain->SetBranchAddress("HLT_Mu3_L1SingleMu5orSingleMu7", &HLT_Mu3_L1SingleMu5orSingleMu7, &b_HLT_Mu3_L1SingleMu5orSingleMu7);
   fChain->SetBranchAddress("HLT_DoublePhoton33_CaloIdL", &HLT_DoublePhoton33_CaloIdL, &b_HLT_DoublePhoton33_CaloIdL);
   fChain->SetBranchAddress("HLT_DoublePhoton70", &HLT_DoublePhoton70, &b_HLT_DoublePhoton70);
   fChain->SetBranchAddress("HLT_DoublePhoton85", &HLT_DoublePhoton85, &b_HLT_DoublePhoton85);
   fChain->SetBranchAddress("HLT_Ele15_WPLoose_Gsf", &HLT_Ele15_WPLoose_Gsf, &b_HLT_Ele15_WPLoose_Gsf);
   fChain->SetBranchAddress("HLT_Ele20_WPLoose_Gsf", &HLT_Ele20_WPLoose_Gsf, &b_HLT_Ele20_WPLoose_Gsf);
   fChain->SetBranchAddress("HLT_DiEle27_WPTightCaloOnly_L1DoubleEG", &HLT_DiEle27_WPTightCaloOnly_L1DoubleEG, &b_HLT_DiEle27_WPTightCaloOnly_L1DoubleEG);
   fChain->SetBranchAddress("HLT_Ele27_WPTight_Gsf", &HLT_Ele27_WPTight_Gsf, &b_HLT_Ele27_WPTight_Gsf);
   fChain->SetBranchAddress("HLT_Ele28_WPTight_Gsf", &HLT_Ele28_WPTight_Gsf, &b_HLT_Ele28_WPTight_Gsf);
   fChain->SetBranchAddress("HLT_Ele30_WPTight_Gsf", &HLT_Ele30_WPTight_Gsf, &b_HLT_Ele30_WPTight_Gsf);
   fChain->SetBranchAddress("HLT_Ele32_WPTight_Gsf", &HLT_Ele32_WPTight_Gsf, &b_HLT_Ele32_WPTight_Gsf);
   fChain->SetBranchAddress("HLT_Ele35_WPTight_Gsf", &HLT_Ele35_WPTight_Gsf, &b_HLT_Ele35_WPTight_Gsf);
   fChain->SetBranchAddress("HLT_Ele35_WPTight_Gsf_L1EGMT", &HLT_Ele35_WPTight_Gsf_L1EGMT, &b_HLT_Ele35_WPTight_Gsf_L1EGMT);
   fChain->SetBranchAddress("HLT_Ele38_WPTight_Gsf", &HLT_Ele38_WPTight_Gsf, &b_HLT_Ele38_WPTight_Gsf);
   fChain->SetBranchAddress("HLT_Ele40_WPTight_Gsf", &HLT_Ele40_WPTight_Gsf, &b_HLT_Ele40_WPTight_Gsf);
   fChain->SetBranchAddress("HLT_Ele32_WPTight_Gsf_L1DoubleEG", &HLT_Ele32_WPTight_Gsf_L1DoubleEG, &b_HLT_Ele32_WPTight_Gsf_L1DoubleEG);
   fChain->SetBranchAddress("HLT_HT300_Beamspot", &HLT_HT300_Beamspot, &b_HLT_HT300_Beamspot);
   fChain->SetBranchAddress("HLT_ZeroBias_Beamspot", &HLT_ZeroBias_Beamspot, &b_HLT_ZeroBias_Beamspot);
   fChain->SetBranchAddress("HLT_IsoMu20_eta2p1_TightChargedIsoPFTauHPS27_eta2p1_TightID_CrossL1", &HLT_IsoMu20_eta2p1_TightChargedIsoPFTauHPS27_eta2p1_TightID_CrossL1, &b_HLT_IsoMu20_eta2p1_TightChargedIsoPFTauHPS27_eta2p1_TightID_CrossL1);
   fChain->SetBranchAddress("HLT_IsoMu27_MediumDeepTauPFTauHPS20_eta2p1_SingleL1", &HLT_IsoMu27_MediumDeepTauPFTauHPS20_eta2p1_SingleL1, &b_HLT_IsoMu27_MediumDeepTauPFTauHPS20_eta2p1_SingleL1);
   fChain->SetBranchAddress("HLT_IsoMu20", &HLT_IsoMu20, &b_HLT_IsoMu20);
   fChain->SetBranchAddress("HLT_IsoMu24", &HLT_IsoMu24, &b_HLT_IsoMu24);
   fChain->SetBranchAddress("HLT_IsoMu24_eta2p1", &HLT_IsoMu24_eta2p1, &b_HLT_IsoMu24_eta2p1);
   fChain->SetBranchAddress("HLT_IsoMu27", &HLT_IsoMu27, &b_HLT_IsoMu27);
   fChain->SetBranchAddress("HLT_UncorrectedJetE30_NoBPTX", &HLT_UncorrectedJetE30_NoBPTX, &b_HLT_UncorrectedJetE30_NoBPTX);
   fChain->SetBranchAddress("HLT_UncorrectedJetE30_NoBPTX3BX", &HLT_UncorrectedJetE30_NoBPTX3BX, &b_HLT_UncorrectedJetE30_NoBPTX3BX);
   fChain->SetBranchAddress("HLT_UncorrectedJetE60_NoBPTX3BX", &HLT_UncorrectedJetE60_NoBPTX3BX, &b_HLT_UncorrectedJetE60_NoBPTX3BX);
   fChain->SetBranchAddress("HLT_UncorrectedJetE70_NoBPTX3BX", &HLT_UncorrectedJetE70_NoBPTX3BX, &b_HLT_UncorrectedJetE70_NoBPTX3BX);
   fChain->SetBranchAddress("HLT_L1SingleMu18", &HLT_L1SingleMu18, &b_HLT_L1SingleMu18);
   fChain->SetBranchAddress("HLT_L1SingleMu25", &HLT_L1SingleMu25, &b_HLT_L1SingleMu25);
   fChain->SetBranchAddress("HLT_L1SingleMuCosmics", &HLT_L1SingleMuCosmics, &b_HLT_L1SingleMuCosmics);
   fChain->SetBranchAddress("HLT_L2Mu10_NoVertex_NoBPTX3BX", &HLT_L2Mu10_NoVertex_NoBPTX3BX, &b_HLT_L2Mu10_NoVertex_NoBPTX3BX);
   fChain->SetBranchAddress("HLT_L2Mu10_NoVertex_NoBPTX", &HLT_L2Mu10_NoVertex_NoBPTX, &b_HLT_L2Mu10_NoVertex_NoBPTX);
   fChain->SetBranchAddress("HLT_L2Mu45_NoVertex_3Sta_NoBPTX3BX", &HLT_L2Mu45_NoVertex_3Sta_NoBPTX3BX, &b_HLT_L2Mu45_NoVertex_3Sta_NoBPTX3BX);
   fChain->SetBranchAddress("HLT_L2Mu40_NoVertex_3Sta_NoBPTX3BX", &HLT_L2Mu40_NoVertex_3Sta_NoBPTX3BX, &b_HLT_L2Mu40_NoVertex_3Sta_NoBPTX3BX);
   fChain->SetBranchAddress("HLT_L2Mu23NoVtx_2Cha", &HLT_L2Mu23NoVtx_2Cha, &b_HLT_L2Mu23NoVtx_2Cha);
   fChain->SetBranchAddress("HLT_L2Mu23NoVtx_2Cha_CosmicSeed", &HLT_L2Mu23NoVtx_2Cha_CosmicSeed, &b_HLT_L2Mu23NoVtx_2Cha_CosmicSeed);
   fChain->SetBranchAddress("HLT_DoubleL2Mu30NoVtx_2Cha_CosmicSeed_Eta2p4", &HLT_DoubleL2Mu30NoVtx_2Cha_CosmicSeed_Eta2p4, &b_HLT_DoubleL2Mu30NoVtx_2Cha_CosmicSeed_Eta2p4);
   fChain->SetBranchAddress("HLT_DoubleL2Mu30NoVtx_2Cha_Eta2p4", &HLT_DoubleL2Mu30NoVtx_2Cha_Eta2p4, &b_HLT_DoubleL2Mu30NoVtx_2Cha_Eta2p4);
   fChain->SetBranchAddress("HLT_DoubleL2Mu50", &HLT_DoubleL2Mu50, &b_HLT_DoubleL2Mu50);
   fChain->SetBranchAddress("HLT_DoubleL2Mu23NoVtx_2Cha_CosmicSeed", &HLT_DoubleL2Mu23NoVtx_2Cha_CosmicSeed, &b_HLT_DoubleL2Mu23NoVtx_2Cha_CosmicSeed);
   fChain->SetBranchAddress("HLT_DoubleL2Mu25NoVtx_2Cha_CosmicSeed", &HLT_DoubleL2Mu25NoVtx_2Cha_CosmicSeed, &b_HLT_DoubleL2Mu25NoVtx_2Cha_CosmicSeed);
   fChain->SetBranchAddress("HLT_DoubleL2Mu25NoVtx_2Cha_CosmicSeed_Eta2p4", &HLT_DoubleL2Mu25NoVtx_2Cha_CosmicSeed_Eta2p4, &b_HLT_DoubleL2Mu25NoVtx_2Cha_CosmicSeed_Eta2p4);
   fChain->SetBranchAddress("HLT_DoubleL2Mu23NoVtx_2Cha", &HLT_DoubleL2Mu23NoVtx_2Cha, &b_HLT_DoubleL2Mu23NoVtx_2Cha);
   fChain->SetBranchAddress("HLT_DoubleL2Mu25NoVtx_2Cha", &HLT_DoubleL2Mu25NoVtx_2Cha, &b_HLT_DoubleL2Mu25NoVtx_2Cha);
   fChain->SetBranchAddress("HLT_DoubleL2Mu25NoVtx_2Cha_Eta2p4", &HLT_DoubleL2Mu25NoVtx_2Cha_Eta2p4, &b_HLT_DoubleL2Mu25NoVtx_2Cha_Eta2p4);
   fChain->SetBranchAddress("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL", &HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL, &b_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL);
   fChain->SetBranchAddress("HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL", &HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL, &b_HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL);
   fChain->SetBranchAddress("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ", &HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ, &b_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ);
   fChain->SetBranchAddress("HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ", &HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ, &b_HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ);
   fChain->SetBranchAddress("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8", &HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8, &b_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8);
   fChain->SetBranchAddress("HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_Mass8", &HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_Mass8, &b_HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_Mass8);
   fChain->SetBranchAddress("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8", &HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8, &b_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8);
   fChain->SetBranchAddress("HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_Mass3p8", &HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_Mass3p8, &b_HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_Mass3p8);
   fChain->SetBranchAddress("HLT_Mu25_TkMu0_Onia", &HLT_Mu25_TkMu0_Onia, &b_HLT_Mu25_TkMu0_Onia);
   fChain->SetBranchAddress("HLT_Mu30_TkMu0_Psi", &HLT_Mu30_TkMu0_Psi, &b_HLT_Mu30_TkMu0_Psi);
   fChain->SetBranchAddress("HLT_Mu30_TkMu0_Upsilon", &HLT_Mu30_TkMu0_Upsilon, &b_HLT_Mu30_TkMu0_Upsilon);
   fChain->SetBranchAddress("HLT_Mu20_TkMu0_Phi", &HLT_Mu20_TkMu0_Phi, &b_HLT_Mu20_TkMu0_Phi);
   fChain->SetBranchAddress("HLT_Mu25_TkMu0_Phi", &HLT_Mu25_TkMu0_Phi, &b_HLT_Mu25_TkMu0_Phi);
   fChain->SetBranchAddress("HLT_Mu15", &HLT_Mu15, &b_HLT_Mu15);
   fChain->SetBranchAddress("HLT_Mu20", &HLT_Mu20, &b_HLT_Mu20);
   fChain->SetBranchAddress("HLT_Mu27", &HLT_Mu27, &b_HLT_Mu27);
   fChain->SetBranchAddress("HLT_Mu50", &HLT_Mu50, &b_HLT_Mu50);
   fChain->SetBranchAddress("HLT_Mu55", &HLT_Mu55, &b_HLT_Mu55);
   fChain->SetBranchAddress("HLT_CascadeMu100", &HLT_CascadeMu100, &b_HLT_CascadeMu100);
   fChain->SetBranchAddress("HLT_HighPtTkMu100", &HLT_HighPtTkMu100, &b_HLT_HighPtTkMu100);
   fChain->SetBranchAddress("HLT_DiPFJetAve40", &HLT_DiPFJetAve40, &b_HLT_DiPFJetAve40);
   fChain->SetBranchAddress("HLT_DiPFJetAve60", &HLT_DiPFJetAve60, &b_HLT_DiPFJetAve60);
   fChain->SetBranchAddress("HLT_DiPFJetAve80", &HLT_DiPFJetAve80, &b_HLT_DiPFJetAve80);
   fChain->SetBranchAddress("HLT_DiPFJetAve140", &HLT_DiPFJetAve140, &b_HLT_DiPFJetAve140);
   fChain->SetBranchAddress("HLT_DiPFJetAve200", &HLT_DiPFJetAve200, &b_HLT_DiPFJetAve200);
   fChain->SetBranchAddress("HLT_DiPFJetAve260", &HLT_DiPFJetAve260, &b_HLT_DiPFJetAve260);
   fChain->SetBranchAddress("HLT_DiPFJetAve320", &HLT_DiPFJetAve320, &b_HLT_DiPFJetAve320);
   fChain->SetBranchAddress("HLT_DiPFJetAve400", &HLT_DiPFJetAve400, &b_HLT_DiPFJetAve400);
   fChain->SetBranchAddress("HLT_DiPFJetAve500", &HLT_DiPFJetAve500, &b_HLT_DiPFJetAve500);
   fChain->SetBranchAddress("HLT_DiPFJetAve60_HFJEC", &HLT_DiPFJetAve60_HFJEC, &b_HLT_DiPFJetAve60_HFJEC);
   fChain->SetBranchAddress("HLT_DiPFJetAve80_HFJEC", &HLT_DiPFJetAve80_HFJEC, &b_HLT_DiPFJetAve80_HFJEC);
   fChain->SetBranchAddress("HLT_DiPFJetAve100_HFJEC", &HLT_DiPFJetAve100_HFJEC, &b_HLT_DiPFJetAve100_HFJEC);
   fChain->SetBranchAddress("HLT_DiPFJetAve160_HFJEC", &HLT_DiPFJetAve160_HFJEC, &b_HLT_DiPFJetAve160_HFJEC);
   fChain->SetBranchAddress("HLT_DiPFJetAve220_HFJEC", &HLT_DiPFJetAve220_HFJEC, &b_HLT_DiPFJetAve220_HFJEC);
   fChain->SetBranchAddress("HLT_DiPFJetAve300_HFJEC", &HLT_DiPFJetAve300_HFJEC, &b_HLT_DiPFJetAve300_HFJEC);
   fChain->SetBranchAddress("HLT_AK8PFJet40", &HLT_AK8PFJet40, &b_HLT_AK8PFJet40);
   fChain->SetBranchAddress("HLT_AK8PFJet60", &HLT_AK8PFJet60, &b_HLT_AK8PFJet60);
   fChain->SetBranchAddress("HLT_AK8PFJet80", &HLT_AK8PFJet80, &b_HLT_AK8PFJet80);
   fChain->SetBranchAddress("HLT_AK8PFJet140", &HLT_AK8PFJet140, &b_HLT_AK8PFJet140);
   fChain->SetBranchAddress("HLT_AK8PFJet200", &HLT_AK8PFJet200, &b_HLT_AK8PFJet200);
   fChain->SetBranchAddress("HLT_AK8PFJet260", &HLT_AK8PFJet260, &b_HLT_AK8PFJet260);
   fChain->SetBranchAddress("HLT_AK8PFJet320", &HLT_AK8PFJet320, &b_HLT_AK8PFJet320);
   fChain->SetBranchAddress("HLT_AK8PFJet400", &HLT_AK8PFJet400, &b_HLT_AK8PFJet400);
   fChain->SetBranchAddress("HLT_AK8PFJet450", &HLT_AK8PFJet450, &b_HLT_AK8PFJet450);
   fChain->SetBranchAddress("HLT_AK8PFJet500", &HLT_AK8PFJet500, &b_HLT_AK8PFJet500);
   fChain->SetBranchAddress("HLT_AK8PFJet550", &HLT_AK8PFJet550, &b_HLT_AK8PFJet550);
   fChain->SetBranchAddress("HLT_PFJet40", &HLT_PFJet40, &b_HLT_PFJet40);
   fChain->SetBranchAddress("HLT_PFJet60", &HLT_PFJet60, &b_HLT_PFJet60);
   fChain->SetBranchAddress("HLT_PFJet80", &HLT_PFJet80, &b_HLT_PFJet80);
   fChain->SetBranchAddress("HLT_PFJet110", &HLT_PFJet110, &b_HLT_PFJet110);
   fChain->SetBranchAddress("HLT_PFJet140", &HLT_PFJet140, &b_HLT_PFJet140);
   fChain->SetBranchAddress("HLT_PFJet200", &HLT_PFJet200, &b_HLT_PFJet200);
   fChain->SetBranchAddress("HLT_PFJet260", &HLT_PFJet260, &b_HLT_PFJet260);
   fChain->SetBranchAddress("HLT_PFJet320", &HLT_PFJet320, &b_HLT_PFJet320);
   fChain->SetBranchAddress("HLT_PFJet400", &HLT_PFJet400, &b_HLT_PFJet400);
   fChain->SetBranchAddress("HLT_PFJet450", &HLT_PFJet450, &b_HLT_PFJet450);
   fChain->SetBranchAddress("HLT_PFJet500", &HLT_PFJet500, &b_HLT_PFJet500);
   fChain->SetBranchAddress("HLT_PFJet550", &HLT_PFJet550, &b_HLT_PFJet550);
   fChain->SetBranchAddress("HLT_PFJetFwd15", &HLT_PFJetFwd15, &b_HLT_PFJetFwd15);
   fChain->SetBranchAddress("HLT_PFJetFwd25", &HLT_PFJetFwd25, &b_HLT_PFJetFwd25);
   fChain->SetBranchAddress("HLT_PFJetFwd40", &HLT_PFJetFwd40, &b_HLT_PFJetFwd40);
   fChain->SetBranchAddress("HLT_PFJetFwd60", &HLT_PFJetFwd60, &b_HLT_PFJetFwd60);
   fChain->SetBranchAddress("HLT_PFJetFwd80", &HLT_PFJetFwd80, &b_HLT_PFJetFwd80);
   fChain->SetBranchAddress("HLT_PFJetFwd140", &HLT_PFJetFwd140, &b_HLT_PFJetFwd140);
   fChain->SetBranchAddress("HLT_PFJetFwd200", &HLT_PFJetFwd200, &b_HLT_PFJetFwd200);
   fChain->SetBranchAddress("HLT_PFJetFwd260", &HLT_PFJetFwd260, &b_HLT_PFJetFwd260);
   fChain->SetBranchAddress("HLT_PFJetFwd320", &HLT_PFJetFwd320, &b_HLT_PFJetFwd320);
   fChain->SetBranchAddress("HLT_PFJetFwd400", &HLT_PFJetFwd400, &b_HLT_PFJetFwd400);
   fChain->SetBranchAddress("HLT_PFJetFwd450", &HLT_PFJetFwd450, &b_HLT_PFJetFwd450);
   fChain->SetBranchAddress("HLT_PFJetFwd500", &HLT_PFJetFwd500, &b_HLT_PFJetFwd500);
   fChain->SetBranchAddress("HLT_AK8PFJetFwd15", &HLT_AK8PFJetFwd15, &b_HLT_AK8PFJetFwd15);
   fChain->SetBranchAddress("HLT_AK8PFJetFwd25", &HLT_AK8PFJetFwd25, &b_HLT_AK8PFJetFwd25);
   fChain->SetBranchAddress("HLT_AK8PFJetFwd40", &HLT_AK8PFJetFwd40, &b_HLT_AK8PFJetFwd40);
   fChain->SetBranchAddress("HLT_AK8PFJetFwd60", &HLT_AK8PFJetFwd60, &b_HLT_AK8PFJetFwd60);
   fChain->SetBranchAddress("HLT_AK8PFJetFwd80", &HLT_AK8PFJetFwd80, &b_HLT_AK8PFJetFwd80);
   fChain->SetBranchAddress("HLT_AK8PFJetFwd140", &HLT_AK8PFJetFwd140, &b_HLT_AK8PFJetFwd140);
   fChain->SetBranchAddress("HLT_AK8PFJetFwd200", &HLT_AK8PFJetFwd200, &b_HLT_AK8PFJetFwd200);
   fChain->SetBranchAddress("HLT_AK8PFJetFwd260", &HLT_AK8PFJetFwd260, &b_HLT_AK8PFJetFwd260);
   fChain->SetBranchAddress("HLT_AK8PFJetFwd320", &HLT_AK8PFJetFwd320, &b_HLT_AK8PFJetFwd320);
   fChain->SetBranchAddress("HLT_AK8PFJetFwd400", &HLT_AK8PFJetFwd400, &b_HLT_AK8PFJetFwd400);
   fChain->SetBranchAddress("HLT_AK8PFJetFwd450", &HLT_AK8PFJetFwd450, &b_HLT_AK8PFJetFwd450);
   fChain->SetBranchAddress("HLT_AK8PFJetFwd500", &HLT_AK8PFJetFwd500, &b_HLT_AK8PFJetFwd500);
   fChain->SetBranchAddress("HLT_PFHT180", &HLT_PFHT180, &b_HLT_PFHT180);
   fChain->SetBranchAddress("HLT_PFHT250", &HLT_PFHT250, &b_HLT_PFHT250);
   fChain->SetBranchAddress("HLT_PFHT370", &HLT_PFHT370, &b_HLT_PFHT370);
   fChain->SetBranchAddress("HLT_PFHT430", &HLT_PFHT430, &b_HLT_PFHT430);
   fChain->SetBranchAddress("HLT_PFHT510", &HLT_PFHT510, &b_HLT_PFHT510);
   fChain->SetBranchAddress("HLT_PFHT590", &HLT_PFHT590, &b_HLT_PFHT590);
   fChain->SetBranchAddress("HLT_PFHT680", &HLT_PFHT680, &b_HLT_PFHT680);
   fChain->SetBranchAddress("HLT_PFHT780", &HLT_PFHT780, &b_HLT_PFHT780);
   fChain->SetBranchAddress("HLT_PFHT890", &HLT_PFHT890, &b_HLT_PFHT890);
   fChain->SetBranchAddress("HLT_PFHT1050", &HLT_PFHT1050, &b_HLT_PFHT1050);
   fChain->SetBranchAddress("HLT_PFHT500_PFMET100_PFMHT100_IDTight", &HLT_PFHT500_PFMET100_PFMHT100_IDTight, &b_HLT_PFHT500_PFMET100_PFMHT100_IDTight);
   fChain->SetBranchAddress("HLT_PFHT500_PFMET110_PFMHT110_IDTight", &HLT_PFHT500_PFMET110_PFMHT110_IDTight, &b_HLT_PFHT500_PFMET110_PFMHT110_IDTight);
   fChain->SetBranchAddress("HLT_PFHT700_PFMET85_PFMHT85_IDTight", &HLT_PFHT700_PFMET85_PFMHT85_IDTight, &b_HLT_PFHT700_PFMET85_PFMHT85_IDTight);
   fChain->SetBranchAddress("HLT_PFHT800_PFMET75_PFMHT75_IDTight", &HLT_PFHT800_PFMET75_PFMHT75_IDTight, &b_HLT_PFHT800_PFMET75_PFMHT75_IDTight);
   fChain->SetBranchAddress("HLT_PFMET110_PFMHT110_IDTight", &HLT_PFMET110_PFMHT110_IDTight, &b_HLT_PFMET110_PFMHT110_IDTight);
   fChain->SetBranchAddress("HLT_PFMET120_PFMHT120_IDTight", &HLT_PFMET120_PFMHT120_IDTight, &b_HLT_PFMET120_PFMHT120_IDTight);
   fChain->SetBranchAddress("HLT_PFMET130_PFMHT130_IDTight", &HLT_PFMET130_PFMHT130_IDTight, &b_HLT_PFMET130_PFMHT130_IDTight);
   fChain->SetBranchAddress("HLT_PFMET140_PFMHT140_IDTight", &HLT_PFMET140_PFMHT140_IDTight, &b_HLT_PFMET140_PFMHT140_IDTight);
   fChain->SetBranchAddress("HLT_PFMET120_PFMHT120_IDTight_PFHT60", &HLT_PFMET120_PFMHT120_IDTight_PFHT60, &b_HLT_PFMET120_PFMHT120_IDTight_PFHT60);
   fChain->SetBranchAddress("HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60", &HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60, &b_HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60);
   fChain->SetBranchAddress("HLT_PFMETTypeOne120_PFMHT120_IDTight_PFHT60", &HLT_PFMETTypeOne120_PFMHT120_IDTight_PFHT60, &b_HLT_PFMETTypeOne120_PFMHT120_IDTight_PFHT60);
   fChain->SetBranchAddress("HLT_PFMETTypeOne110_PFMHT110_IDTight", &HLT_PFMETTypeOne110_PFMHT110_IDTight, &b_HLT_PFMETTypeOne110_PFMHT110_IDTight);
   fChain->SetBranchAddress("HLT_PFMETTypeOne120_PFMHT120_IDTight", &HLT_PFMETTypeOne120_PFMHT120_IDTight, &b_HLT_PFMETTypeOne120_PFMHT120_IDTight);
   fChain->SetBranchAddress("HLT_PFMETTypeOne130_PFMHT130_IDTight", &HLT_PFMETTypeOne130_PFMHT130_IDTight, &b_HLT_PFMETTypeOne130_PFMHT130_IDTight);
   fChain->SetBranchAddress("HLT_PFMETTypeOne140_PFMHT140_IDTight", &HLT_PFMETTypeOne140_PFMHT140_IDTight, &b_HLT_PFMETTypeOne140_PFMHT140_IDTight);
   fChain->SetBranchAddress("HLT_PFMETNoMu110_PFMHTNoMu110_IDTight", &HLT_PFMETNoMu110_PFMHTNoMu110_IDTight, &b_HLT_PFMETNoMu110_PFMHTNoMu110_IDTight);
   fChain->SetBranchAddress("HLT_PFMETNoMu120_PFMHTNoMu120_IDTight", &HLT_PFMETNoMu120_PFMHTNoMu120_IDTight, &b_HLT_PFMETNoMu120_PFMHTNoMu120_IDTight);
   fChain->SetBranchAddress("HLT_PFMETNoMu130_PFMHTNoMu130_IDTight", &HLT_PFMETNoMu130_PFMHTNoMu130_IDTight, &b_HLT_PFMETNoMu130_PFMHTNoMu130_IDTight);
   fChain->SetBranchAddress("HLT_PFMETNoMu140_PFMHTNoMu140_IDTight", &HLT_PFMETNoMu140_PFMHTNoMu140_IDTight, &b_HLT_PFMETNoMu140_PFMHTNoMu140_IDTight);
   fChain->SetBranchAddress("HLT_PFMETNoMu110_PFMHTNoMu110_IDTight_FilterHF", &HLT_PFMETNoMu110_PFMHTNoMu110_IDTight_FilterHF, &b_HLT_PFMETNoMu110_PFMHTNoMu110_IDTight_FilterHF);
   fChain->SetBranchAddress("HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_FilterHF", &HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_FilterHF, &b_HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_FilterHF);
   fChain->SetBranchAddress("HLT_PFMETNoMu130_PFMHTNoMu130_IDTight_FilterHF", &HLT_PFMETNoMu130_PFMHTNoMu130_IDTight_FilterHF, &b_HLT_PFMETNoMu130_PFMHTNoMu130_IDTight_FilterHF);
   fChain->SetBranchAddress("HLT_PFMETNoMu140_PFMHTNoMu140_IDTight_FilterHF", &HLT_PFMETNoMu140_PFMHTNoMu140_IDTight_FilterHF, &b_HLT_PFMETNoMu140_PFMHTNoMu140_IDTight_FilterHF);
   fChain->SetBranchAddress("HLT_L1ETMHadSeeds", &HLT_L1ETMHadSeeds, &b_HLT_L1ETMHadSeeds);
   fChain->SetBranchAddress("HLT_CaloMHT90", &HLT_CaloMHT90, &b_HLT_CaloMHT90);
   fChain->SetBranchAddress("HLT_CaloMET90_NotCleaned", &HLT_CaloMET90_NotCleaned, &b_HLT_CaloMET90_NotCleaned);
   fChain->SetBranchAddress("HLT_CaloMET350_NotCleaned", &HLT_CaloMET350_NotCleaned, &b_HLT_CaloMET350_NotCleaned);
   fChain->SetBranchAddress("HLT_PFMET200_NotCleaned", &HLT_PFMET200_NotCleaned, &b_HLT_PFMET200_NotCleaned);
   fChain->SetBranchAddress("HLT_PFMET250_NotCleaned", &HLT_PFMET250_NotCleaned, &b_HLT_PFMET250_NotCleaned);
   fChain->SetBranchAddress("HLT_PFMET300_NotCleaned", &HLT_PFMET300_NotCleaned, &b_HLT_PFMET300_NotCleaned);
   fChain->SetBranchAddress("HLT_PFMET200_BeamHaloCleaned", &HLT_PFMET200_BeamHaloCleaned, &b_HLT_PFMET200_BeamHaloCleaned);
   fChain->SetBranchAddress("HLT_PFMETTypeOne200_BeamHaloCleaned", &HLT_PFMETTypeOne200_BeamHaloCleaned, &b_HLT_PFMETTypeOne200_BeamHaloCleaned);
   fChain->SetBranchAddress("HLT_MET105_IsoTrk50", &HLT_MET105_IsoTrk50, &b_HLT_MET105_IsoTrk50);
   fChain->SetBranchAddress("HLT_MET120_IsoTrk50", &HLT_MET120_IsoTrk50, &b_HLT_MET120_IsoTrk50);
   fChain->SetBranchAddress("HLT_SingleJet30_Mu12_SinglePFJet40", &HLT_SingleJet30_Mu12_SinglePFJet40, &b_HLT_SingleJet30_Mu12_SinglePFJet40);
   fChain->SetBranchAddress("HLT_Mu12eta2p3", &HLT_Mu12eta2p3, &b_HLT_Mu12eta2p3);
   fChain->SetBranchAddress("HLT_Mu12eta2p3_PFJet40", &HLT_Mu12eta2p3_PFJet40, &b_HLT_Mu12eta2p3_PFJet40);
   fChain->SetBranchAddress("HLT_Mu12_DoublePFJets40_PFBTagDeepCSV_p71", &HLT_Mu12_DoublePFJets40_PFBTagDeepCSV_p71, &b_HLT_Mu12_DoublePFJets40_PFBTagDeepCSV_p71);
   fChain->SetBranchAddress("HLT_Mu12_DoublePFJets100_PFBTagDeepCSV_p71", &HLT_Mu12_DoublePFJets100_PFBTagDeepCSV_p71, &b_HLT_Mu12_DoublePFJets100_PFBTagDeepCSV_p71);
   fChain->SetBranchAddress("HLT_Mu12_DoublePFJets200_PFBTagDeepCSV_p71", &HLT_Mu12_DoublePFJets200_PFBTagDeepCSV_p71, &b_HLT_Mu12_DoublePFJets200_PFBTagDeepCSV_p71);
   fChain->SetBranchAddress("HLT_Mu12_DoublePFJets350_PFBTagDeepCSV_p71", &HLT_Mu12_DoublePFJets350_PFBTagDeepCSV_p71, &b_HLT_Mu12_DoublePFJets350_PFBTagDeepCSV_p71);
   fChain->SetBranchAddress("HLT_Mu12_DoublePFJets40MaxDeta1p6_DoublePFBTagDeepCSV_p71", &HLT_Mu12_DoublePFJets40MaxDeta1p6_DoublePFBTagDeepCSV_p71, &b_HLT_Mu12_DoublePFJets40MaxDeta1p6_DoublePFBTagDeepCSV_p71);
   fChain->SetBranchAddress("HLT_Mu12_DoublePFJets54MaxDeta1p6_DoublePFBTagDeepCSV_p71", &HLT_Mu12_DoublePFJets54MaxDeta1p6_DoublePFBTagDeepCSV_p71, &b_HLT_Mu12_DoublePFJets54MaxDeta1p6_DoublePFBTagDeepCSV_p71);
   fChain->SetBranchAddress("HLT_DoublePFJets40_PFBTagDeepCSV_p71", &HLT_DoublePFJets40_PFBTagDeepCSV_p71, &b_HLT_DoublePFJets40_PFBTagDeepCSV_p71);
   fChain->SetBranchAddress("HLT_DoublePFJets100_PFBTagDeepCSV_p71", &HLT_DoublePFJets100_PFBTagDeepCSV_p71, &b_HLT_DoublePFJets100_PFBTagDeepCSV_p71);
   fChain->SetBranchAddress("HLT_DoublePFJets200_PFBTagDeepCSV_p71", &HLT_DoublePFJets200_PFBTagDeepCSV_p71, &b_HLT_DoublePFJets200_PFBTagDeepCSV_p71);
   fChain->SetBranchAddress("HLT_DoublePFJets350_PFBTagDeepCSV_p71", &HLT_DoublePFJets350_PFBTagDeepCSV_p71, &b_HLT_DoublePFJets350_PFBTagDeepCSV_p71);
   fChain->SetBranchAddress("HLT_DoublePFJets116MaxDeta1p6_DoublePFBTagDeepCSV_p71", &HLT_DoublePFJets116MaxDeta1p6_DoublePFBTagDeepCSV_p71, &b_HLT_DoublePFJets116MaxDeta1p6_DoublePFBTagDeepCSV_p71);
   fChain->SetBranchAddress("HLT_DoublePFJets128MaxDeta1p6_DoublePFBTagDeepCSV_p71", &HLT_DoublePFJets128MaxDeta1p6_DoublePFBTagDeepCSV_p71, &b_HLT_DoublePFJets128MaxDeta1p6_DoublePFBTagDeepCSV_p71);
   fChain->SetBranchAddress("HLT_Mu12_DoublePFJets40_PFBTagDeepJet_p71", &HLT_Mu12_DoublePFJets40_PFBTagDeepJet_p71, &b_HLT_Mu12_DoublePFJets40_PFBTagDeepJet_p71);
   fChain->SetBranchAddress("HLT_Mu12_DoublePFJets100_PFBTagDeepJet_p71", &HLT_Mu12_DoublePFJets100_PFBTagDeepJet_p71, &b_HLT_Mu12_DoublePFJets100_PFBTagDeepJet_p71);
   fChain->SetBranchAddress("HLT_Mu12_DoublePFJets200_PFBTagDeepJet_p71", &HLT_Mu12_DoublePFJets200_PFBTagDeepJet_p71, &b_HLT_Mu12_DoublePFJets200_PFBTagDeepJet_p71);
   fChain->SetBranchAddress("HLT_Mu12_DoublePFJets350_PFBTagDeepJet_p71", &HLT_Mu12_DoublePFJets350_PFBTagDeepJet_p71, &b_HLT_Mu12_DoublePFJets350_PFBTagDeepJet_p71);
   fChain->SetBranchAddress("HLT_Mu12_DoublePFJets40MaxDeta1p6_DoublePFBTagDeepJet_p71", &HLT_Mu12_DoublePFJets40MaxDeta1p6_DoublePFBTagDeepJet_p71, &b_HLT_Mu12_DoublePFJets40MaxDeta1p6_DoublePFBTagDeepJet_p71);
   fChain->SetBranchAddress("HLT_Mu12_DoublePFJets54MaxDeta1p6_DoublePFBTagDeepJet_p71", &HLT_Mu12_DoublePFJets54MaxDeta1p6_DoublePFBTagDeepJet_p71, &b_HLT_Mu12_DoublePFJets54MaxDeta1p6_DoublePFBTagDeepJet_p71);
   fChain->SetBranchAddress("HLT_DoublePFJets40_PFBTagDeepJet_p71", &HLT_DoublePFJets40_PFBTagDeepJet_p71, &b_HLT_DoublePFJets40_PFBTagDeepJet_p71);
   fChain->SetBranchAddress("HLT_DoublePFJets100_PFBTagDeepJet_p71", &HLT_DoublePFJets100_PFBTagDeepJet_p71, &b_HLT_DoublePFJets100_PFBTagDeepJet_p71);
   fChain->SetBranchAddress("HLT_DoublePFJets200_PFBTagDeepJet_p71", &HLT_DoublePFJets200_PFBTagDeepJet_p71, &b_HLT_DoublePFJets200_PFBTagDeepJet_p71);
   fChain->SetBranchAddress("HLT_DoublePFJets350_PFBTagDeepJet_p71", &HLT_DoublePFJets350_PFBTagDeepJet_p71, &b_HLT_DoublePFJets350_PFBTagDeepJet_p71);
   fChain->SetBranchAddress("HLT_DoublePFJets116MaxDeta1p6_DoublePFBTagDeepJet_p71", &HLT_DoublePFJets116MaxDeta1p6_DoublePFBTagDeepJet_p71, &b_HLT_DoublePFJets116MaxDeta1p6_DoublePFBTagDeepJet_p71);
   fChain->SetBranchAddress("HLT_DoublePFJets128MaxDeta1p6_DoublePFBTagDeepJet_p71", &HLT_DoublePFJets128MaxDeta1p6_DoublePFBTagDeepJet_p71, &b_HLT_DoublePFJets128MaxDeta1p6_DoublePFBTagDeepJet_p71);
   fChain->SetBranchAddress("HLT_Photon300_NoHE", &HLT_Photon300_NoHE, &b_HLT_Photon300_NoHE);
   fChain->SetBranchAddress("HLT_Mu8_TrkIsoVVL", &HLT_Mu8_TrkIsoVVL, &b_HLT_Mu8_TrkIsoVVL);
   fChain->SetBranchAddress("HLT_Mu8_DiEle12_CaloIdL_TrackIdL_DZ", &HLT_Mu8_DiEle12_CaloIdL_TrackIdL_DZ, &b_HLT_Mu8_DiEle12_CaloIdL_TrackIdL_DZ);
   fChain->SetBranchAddress("HLT_Mu8_DiEle12_CaloIdL_TrackIdL", &HLT_Mu8_DiEle12_CaloIdL_TrackIdL, &b_HLT_Mu8_DiEle12_CaloIdL_TrackIdL);
   fChain->SetBranchAddress("HLT_Mu8_Ele8_CaloIdM_TrackIdM_Mass8_PFHT350_DZ", &HLT_Mu8_Ele8_CaloIdM_TrackIdM_Mass8_PFHT350_DZ, &b_HLT_Mu8_Ele8_CaloIdM_TrackIdM_Mass8_PFHT350_DZ);
   fChain->SetBranchAddress("HLT_Mu8_Ele8_CaloIdM_TrackIdM_Mass8_PFHT350", &HLT_Mu8_Ele8_CaloIdM_TrackIdM_Mass8_PFHT350, &b_HLT_Mu8_Ele8_CaloIdM_TrackIdM_Mass8_PFHT350);
   fChain->SetBranchAddress("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ", &HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ, &b_HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ);
   fChain->SetBranchAddress("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_PFDiJet30", &HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_PFDiJet30, &b_HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_PFDiJet30);
   fChain->SetBranchAddress("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_CaloDiJet30", &HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_CaloDiJet30, &b_HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_CaloDiJet30);
   fChain->SetBranchAddress("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_PFDiJet30_PFBtagDeepCSV_1p5", &HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_PFDiJet30_PFBtagDeepCSV_1p5, &b_HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_PFDiJet30_PFBtagDeepCSV_1p5);
   fChain->SetBranchAddress("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_CaloDiJet30_CaloBtagDeepCSV_1p5", &HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_CaloDiJet30_CaloBtagDeepCSV_1p5, &b_HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_CaloDiJet30_CaloBtagDeepCSV_1p5);
   fChain->SetBranchAddress("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL", &HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL, &b_HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL);
   fChain->SetBranchAddress("HLT_Mu17_TrkIsoVVL", &HLT_Mu17_TrkIsoVVL, &b_HLT_Mu17_TrkIsoVVL);
   fChain->SetBranchAddress("HLT_Mu19_TrkIsoVVL", &HLT_Mu19_TrkIsoVVL, &b_HLT_Mu19_TrkIsoVVL);
   fChain->SetBranchAddress("HLT_BTagMu_AK4DiJet20_Mu5", &HLT_BTagMu_AK4DiJet20_Mu5, &b_HLT_BTagMu_AK4DiJet20_Mu5);
   fChain->SetBranchAddress("HLT_BTagMu_AK4DiJet40_Mu5", &HLT_BTagMu_AK4DiJet40_Mu5, &b_HLT_BTagMu_AK4DiJet40_Mu5);
   fChain->SetBranchAddress("HLT_BTagMu_AK4DiJet70_Mu5", &HLT_BTagMu_AK4DiJet70_Mu5, &b_HLT_BTagMu_AK4DiJet70_Mu5);
   fChain->SetBranchAddress("HLT_BTagMu_AK4DiJet110_Mu5", &HLT_BTagMu_AK4DiJet110_Mu5, &b_HLT_BTagMu_AK4DiJet110_Mu5);
   fChain->SetBranchAddress("HLT_BTagMu_AK4DiJet170_Mu5", &HLT_BTagMu_AK4DiJet170_Mu5, &b_HLT_BTagMu_AK4DiJet170_Mu5);
   fChain->SetBranchAddress("HLT_BTagMu_AK4Jet300_Mu5", &HLT_BTagMu_AK4Jet300_Mu5, &b_HLT_BTagMu_AK4Jet300_Mu5);
   fChain->SetBranchAddress("HLT_BTagMu_AK8DiJet170_Mu5", &HLT_BTagMu_AK8DiJet170_Mu5, &b_HLT_BTagMu_AK8DiJet170_Mu5);
   fChain->SetBranchAddress("HLT_BTagMu_AK8Jet170_DoubleMu5", &HLT_BTagMu_AK8Jet170_DoubleMu5, &b_HLT_BTagMu_AK8Jet170_DoubleMu5);
   fChain->SetBranchAddress("HLT_BTagMu_AK8Jet300_Mu5", &HLT_BTagMu_AK8Jet300_Mu5, &b_HLT_BTagMu_AK8Jet300_Mu5);
   fChain->SetBranchAddress("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ", &HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ, &b_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ);
   fChain->SetBranchAddress("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL", &HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL, &b_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL);
   fChain->SetBranchAddress("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ", &HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ, &b_HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ);
   fChain->SetBranchAddress("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL", &HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL, &b_HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL);
   fChain->SetBranchAddress("HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL", &HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL, &b_HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL);
   fChain->SetBranchAddress("HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ", &HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ, &b_HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ);
   fChain->SetBranchAddress("HLT_Photon20", &HLT_Photon20, &b_HLT_Photon20);
   fChain->SetBranchAddress("HLT_Photon33", &HLT_Photon33, &b_HLT_Photon33);
   fChain->SetBranchAddress("HLT_Photon50", &HLT_Photon50, &b_HLT_Photon50);
   fChain->SetBranchAddress("HLT_Photon75", &HLT_Photon75, &b_HLT_Photon75);
   fChain->SetBranchAddress("HLT_Photon90", &HLT_Photon90, &b_HLT_Photon90);
   fChain->SetBranchAddress("HLT_Photon120", &HLT_Photon120, &b_HLT_Photon120);
   fChain->SetBranchAddress("HLT_Photon150", &HLT_Photon150, &b_HLT_Photon150);
   fChain->SetBranchAddress("HLT_Photon175", &HLT_Photon175, &b_HLT_Photon175);
   fChain->SetBranchAddress("HLT_Photon200", &HLT_Photon200, &b_HLT_Photon200);
   fChain->SetBranchAddress("HLT_Photon30EB_TightID_TightIso", &HLT_Photon30EB_TightID_TightIso, &b_HLT_Photon30EB_TightID_TightIso);
   fChain->SetBranchAddress("HLT_Photon110EB_TightID_TightIso", &HLT_Photon110EB_TightID_TightIso, &b_HLT_Photon110EB_TightID_TightIso);
   fChain->SetBranchAddress("HLT_Photon100EBHE10", &HLT_Photon100EBHE10, &b_HLT_Photon100EBHE10);
   fChain->SetBranchAddress("HLT_Photon50_R9Id90_HE10_IsoM", &HLT_Photon50_R9Id90_HE10_IsoM, &b_HLT_Photon50_R9Id90_HE10_IsoM);
   fChain->SetBranchAddress("HLT_Photon75_R9Id90_HE10_IsoM", &HLT_Photon75_R9Id90_HE10_IsoM, &b_HLT_Photon75_R9Id90_HE10_IsoM);
   fChain->SetBranchAddress("HLT_Photon90_R9Id90_HE10_IsoM", &HLT_Photon90_R9Id90_HE10_IsoM, &b_HLT_Photon90_R9Id90_HE10_IsoM);
   fChain->SetBranchAddress("HLT_Photon120_R9Id90_HE10_IsoM", &HLT_Photon120_R9Id90_HE10_IsoM, &b_HLT_Photon120_R9Id90_HE10_IsoM);
   fChain->SetBranchAddress("HLT_Photon165_R9Id90_HE10_IsoM", &HLT_Photon165_R9Id90_HE10_IsoM, &b_HLT_Photon165_R9Id90_HE10_IsoM);
   fChain->SetBranchAddress("HLT_Diphoton30_22_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass90", &HLT_Diphoton30_22_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass90, &b_HLT_Diphoton30_22_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass90);
   fChain->SetBranchAddress("HLT_Diphoton30_22_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass95", &HLT_Diphoton30_22_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass95, &b_HLT_Diphoton30_22_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass95);
   fChain->SetBranchAddress("HLT_Photon35_TwoProngs35", &HLT_Photon35_TwoProngs35, &b_HLT_Photon35_TwoProngs35);
   fChain->SetBranchAddress("HLT_IsoMu24_TwoProngs35", &HLT_IsoMu24_TwoProngs35, &b_HLT_IsoMu24_TwoProngs35);
   fChain->SetBranchAddress("HLT_Dimuon0_Jpsi_L1_NoOS", &HLT_Dimuon0_Jpsi_L1_NoOS, &b_HLT_Dimuon0_Jpsi_L1_NoOS);
   fChain->SetBranchAddress("HLT_Dimuon0_Jpsi_NoVertexing_NoOS", &HLT_Dimuon0_Jpsi_NoVertexing_NoOS, &b_HLT_Dimuon0_Jpsi_NoVertexing_NoOS);
   fChain->SetBranchAddress("HLT_Dimuon0_Jpsi", &HLT_Dimuon0_Jpsi, &b_HLT_Dimuon0_Jpsi);
   fChain->SetBranchAddress("HLT_Dimuon0_Jpsi_NoVertexing", &HLT_Dimuon0_Jpsi_NoVertexing, &b_HLT_Dimuon0_Jpsi_NoVertexing);
   fChain->SetBranchAddress("HLT_Dimuon0_Jpsi_L1_4R_0er1p5R", &HLT_Dimuon0_Jpsi_L1_4R_0er1p5R, &b_HLT_Dimuon0_Jpsi_L1_4R_0er1p5R);
   fChain->SetBranchAddress("HLT_Dimuon0_Jpsi_NoVertexing_L1_4R_0er1p5R", &HLT_Dimuon0_Jpsi_NoVertexing_L1_4R_0er1p5R, &b_HLT_Dimuon0_Jpsi_NoVertexing_L1_4R_0er1p5R);
   fChain->SetBranchAddress("HLT_Dimuon0_Jpsi3p5_Muon2", &HLT_Dimuon0_Jpsi3p5_Muon2, &b_HLT_Dimuon0_Jpsi3p5_Muon2);
   fChain->SetBranchAddress("HLT_Dimuon0_Upsilon_L1_4p5", &HLT_Dimuon0_Upsilon_L1_4p5, &b_HLT_Dimuon0_Upsilon_L1_4p5);
   fChain->SetBranchAddress("HLT_Dimuon0_Upsilon_L1_5", &HLT_Dimuon0_Upsilon_L1_5, &b_HLT_Dimuon0_Upsilon_L1_5);
   fChain->SetBranchAddress("HLT_Dimuon0_Upsilon_L1_4p5NoOS", &HLT_Dimuon0_Upsilon_L1_4p5NoOS, &b_HLT_Dimuon0_Upsilon_L1_4p5NoOS);
   fChain->SetBranchAddress("HLT_Dimuon0_Upsilon_L1_4p5er2p0", &HLT_Dimuon0_Upsilon_L1_4p5er2p0, &b_HLT_Dimuon0_Upsilon_L1_4p5er2p0);
   fChain->SetBranchAddress("HLT_Dimuon0_Upsilon_L1_4p5er2p0M", &HLT_Dimuon0_Upsilon_L1_4p5er2p0M, &b_HLT_Dimuon0_Upsilon_L1_4p5er2p0M);
   fChain->SetBranchAddress("HLT_Dimuon0_Upsilon_NoVertexing", &HLT_Dimuon0_Upsilon_NoVertexing, &b_HLT_Dimuon0_Upsilon_NoVertexing);
   fChain->SetBranchAddress("HLT_Dimuon0_Upsilon_L1_5M", &HLT_Dimuon0_Upsilon_L1_5M, &b_HLT_Dimuon0_Upsilon_L1_5M);
   fChain->SetBranchAddress("HLT_Dimuon0_LowMass_L1_0er1p5R", &HLT_Dimuon0_LowMass_L1_0er1p5R, &b_HLT_Dimuon0_LowMass_L1_0er1p5R);
   fChain->SetBranchAddress("HLT_Dimuon0_LowMass_L1_0er1p5", &HLT_Dimuon0_LowMass_L1_0er1p5, &b_HLT_Dimuon0_LowMass_L1_0er1p5);
   fChain->SetBranchAddress("HLT_Dimuon0_LowMass", &HLT_Dimuon0_LowMass, &b_HLT_Dimuon0_LowMass);
   fChain->SetBranchAddress("HLT_Dimuon0_LowMass_L1_4", &HLT_Dimuon0_LowMass_L1_4, &b_HLT_Dimuon0_LowMass_L1_4);
   fChain->SetBranchAddress("HLT_Dimuon0_LowMass_L1_4R", &HLT_Dimuon0_LowMass_L1_4R, &b_HLT_Dimuon0_LowMass_L1_4R);
   fChain->SetBranchAddress("HLT_Dimuon0_LowMass_L1_TM530", &HLT_Dimuon0_LowMass_L1_TM530, &b_HLT_Dimuon0_LowMass_L1_TM530);
   fChain->SetBranchAddress("HLT_Dimuon0_Upsilon_Muon_L1_TM0", &HLT_Dimuon0_Upsilon_Muon_L1_TM0, &b_HLT_Dimuon0_Upsilon_Muon_L1_TM0);
   fChain->SetBranchAddress("HLT_Dimuon0_Upsilon_Muon_NoL1Mass", &HLT_Dimuon0_Upsilon_Muon_NoL1Mass, &b_HLT_Dimuon0_Upsilon_Muon_NoL1Mass);
   fChain->SetBranchAddress("HLT_TripleMu_5_3_3_Mass3p8_DZ", &HLT_TripleMu_5_3_3_Mass3p8_DZ, &b_HLT_TripleMu_5_3_3_Mass3p8_DZ);
   fChain->SetBranchAddress("HLT_TripleMu_10_5_5_DZ", &HLT_TripleMu_10_5_5_DZ, &b_HLT_TripleMu_10_5_5_DZ);
   fChain->SetBranchAddress("HLT_TripleMu_12_10_5", &HLT_TripleMu_12_10_5, &b_HLT_TripleMu_12_10_5);
   fChain->SetBranchAddress("HLT_Tau3Mu_Mu7_Mu1_TkMu1_Tau15", &HLT_Tau3Mu_Mu7_Mu1_TkMu1_Tau15, &b_HLT_Tau3Mu_Mu7_Mu1_TkMu1_Tau15);
   fChain->SetBranchAddress("HLT_Tau3Mu_Mu7_Mu1_TkMu1_Tau15_Charge1", &HLT_Tau3Mu_Mu7_Mu1_TkMu1_Tau15_Charge1, &b_HLT_Tau3Mu_Mu7_Mu1_TkMu1_Tau15_Charge1);
   fChain->SetBranchAddress("HLT_Tau3Mu_Mu7_Mu1_TkMu1_IsoTau15", &HLT_Tau3Mu_Mu7_Mu1_TkMu1_IsoTau15, &b_HLT_Tau3Mu_Mu7_Mu1_TkMu1_IsoTau15);
   fChain->SetBranchAddress("HLT_Tau3Mu_Mu7_Mu1_TkMu1_IsoTau15_Charge1", &HLT_Tau3Mu_Mu7_Mu1_TkMu1_IsoTau15_Charge1, &b_HLT_Tau3Mu_Mu7_Mu1_TkMu1_IsoTau15_Charge1);
   fChain->SetBranchAddress("HLT_DoubleMu3_DZ_PFMET50_PFMHT60", &HLT_DoubleMu3_DZ_PFMET50_PFMHT60, &b_HLT_DoubleMu3_DZ_PFMET50_PFMHT60);
   fChain->SetBranchAddress("HLT_DoubleMu3_DZ_PFMET70_PFMHT70", &HLT_DoubleMu3_DZ_PFMET70_PFMHT70, &b_HLT_DoubleMu3_DZ_PFMET70_PFMHT70);
   fChain->SetBranchAddress("HLT_DoubleMu3_DZ_PFMET90_PFMHT90", &HLT_DoubleMu3_DZ_PFMET90_PFMHT90, &b_HLT_DoubleMu3_DZ_PFMET90_PFMHT90);
   fChain->SetBranchAddress("HLT_DoubleMu3_Trk_Tau3mu_NoL1Mass", &HLT_DoubleMu3_Trk_Tau3mu_NoL1Mass, &b_HLT_DoubleMu3_Trk_Tau3mu_NoL1Mass);
   fChain->SetBranchAddress("HLT_DoubleMu4_Jpsi_Displaced", &HLT_DoubleMu4_Jpsi_Displaced, &b_HLT_DoubleMu4_Jpsi_Displaced);
   fChain->SetBranchAddress("HLT_DoubleMu4_Jpsi_NoVertexing", &HLT_DoubleMu4_Jpsi_NoVertexing, &b_HLT_DoubleMu4_Jpsi_NoVertexing);
   fChain->SetBranchAddress("HLT_DoubleMu4_JpsiTrkTrk_Displaced", &HLT_DoubleMu4_JpsiTrkTrk_Displaced, &b_HLT_DoubleMu4_JpsiTrkTrk_Displaced);
   fChain->SetBranchAddress("HLT_DoubleMu4_JpsiTrk_Bc", &HLT_DoubleMu4_JpsiTrk_Bc, &b_HLT_DoubleMu4_JpsiTrk_Bc);
   fChain->SetBranchAddress("HLT_DoubleMu43NoFiltersNoVtx", &HLT_DoubleMu43NoFiltersNoVtx, &b_HLT_DoubleMu43NoFiltersNoVtx);
   fChain->SetBranchAddress("HLT_DoubleMu48NoFiltersNoVtx", &HLT_DoubleMu48NoFiltersNoVtx, &b_HLT_DoubleMu48NoFiltersNoVtx);
   fChain->SetBranchAddress("HLT_Mu43NoFiltersNoVtx_Photon43_CaloIdL", &HLT_Mu43NoFiltersNoVtx_Photon43_CaloIdL, &b_HLT_Mu43NoFiltersNoVtx_Photon43_CaloIdL);
   fChain->SetBranchAddress("HLT_Mu48NoFiltersNoVtx_Photon48_CaloIdL", &HLT_Mu48NoFiltersNoVtx_Photon48_CaloIdL, &b_HLT_Mu48NoFiltersNoVtx_Photon48_CaloIdL);
   fChain->SetBranchAddress("HLT_Mu38NoFiltersNoVtxDisplaced_Photon38_CaloIdL", &HLT_Mu38NoFiltersNoVtxDisplaced_Photon38_CaloIdL, &b_HLT_Mu38NoFiltersNoVtxDisplaced_Photon38_CaloIdL);
   fChain->SetBranchAddress("HLT_Mu43NoFiltersNoVtxDisplaced_Photon43_CaloIdL", &HLT_Mu43NoFiltersNoVtxDisplaced_Photon43_CaloIdL, &b_HLT_Mu43NoFiltersNoVtxDisplaced_Photon43_CaloIdL);
   fChain->SetBranchAddress("HLT_HT425", &HLT_HT425, &b_HLT_HT425);
   fChain->SetBranchAddress("HLT_HT430_DisplacedDijet40_DisplacedTrack", &HLT_HT430_DisplacedDijet40_DisplacedTrack, &b_HLT_HT430_DisplacedDijet40_DisplacedTrack);
   fChain->SetBranchAddress("HLT_HT500_DisplacedDijet40_DisplacedTrack", &HLT_HT500_DisplacedDijet40_DisplacedTrack, &b_HLT_HT500_DisplacedDijet40_DisplacedTrack);
   fChain->SetBranchAddress("HLT_HT430_DisplacedDijet60_DisplacedTrack", &HLT_HT430_DisplacedDijet60_DisplacedTrack, &b_HLT_HT430_DisplacedDijet60_DisplacedTrack);
   fChain->SetBranchAddress("HLT_HT400_DisplacedDijet40_DisplacedTrack", &HLT_HT400_DisplacedDijet40_DisplacedTrack, &b_HLT_HT400_DisplacedDijet40_DisplacedTrack);
   fChain->SetBranchAddress("HLT_HT650_DisplacedDijet60_Inclusive", &HLT_HT650_DisplacedDijet60_Inclusive, &b_HLT_HT650_DisplacedDijet60_Inclusive);
   fChain->SetBranchAddress("HLT_HT550_DisplacedDijet60_Inclusive", &HLT_HT550_DisplacedDijet60_Inclusive, &b_HLT_HT550_DisplacedDijet60_Inclusive);
   fChain->SetBranchAddress("HLT_DiJet110_35_Mjj650_PFMET110", &HLT_DiJet110_35_Mjj650_PFMET110, &b_HLT_DiJet110_35_Mjj650_PFMET110);
   fChain->SetBranchAddress("HLT_DiJet110_35_Mjj650_PFMET120", &HLT_DiJet110_35_Mjj650_PFMET120, &b_HLT_DiJet110_35_Mjj650_PFMET120);
   fChain->SetBranchAddress("HLT_DiJet110_35_Mjj650_PFMET130", &HLT_DiJet110_35_Mjj650_PFMET130, &b_HLT_DiJet110_35_Mjj650_PFMET130);
   fChain->SetBranchAddress("HLT_TripleJet110_35_35_Mjj650_PFMET110", &HLT_TripleJet110_35_35_Mjj650_PFMET110, &b_HLT_TripleJet110_35_35_Mjj650_PFMET110);
   fChain->SetBranchAddress("HLT_TripleJet110_35_35_Mjj650_PFMET120", &HLT_TripleJet110_35_35_Mjj650_PFMET120, &b_HLT_TripleJet110_35_35_Mjj650_PFMET120);
   fChain->SetBranchAddress("HLT_TripleJet110_35_35_Mjj650_PFMET130", &HLT_TripleJet110_35_35_Mjj650_PFMET130, &b_HLT_TripleJet110_35_35_Mjj650_PFMET130);
   fChain->SetBranchAddress("HLT_Ele30_eta2p1_WPTight_Gsf_CentralPFJet35_EleCleaned", &HLT_Ele30_eta2p1_WPTight_Gsf_CentralPFJet35_EleCleaned, &b_HLT_Ele30_eta2p1_WPTight_Gsf_CentralPFJet35_EleCleaned);
   fChain->SetBranchAddress("HLT_Ele28_eta2p1_WPTight_Gsf_HT150", &HLT_Ele28_eta2p1_WPTight_Gsf_HT150, &b_HLT_Ele28_eta2p1_WPTight_Gsf_HT150);
   fChain->SetBranchAddress("HLT_Ele28_HighEta_SC20_Mass55", &HLT_Ele28_HighEta_SC20_Mass55, &b_HLT_Ele28_HighEta_SC20_Mass55);
   fChain->SetBranchAddress("HLT_Ele15_IsoVVVL_PFHT450_CaloBTagDeepCSV_4p5", &HLT_Ele15_IsoVVVL_PFHT450_CaloBTagDeepCSV_4p5, &b_HLT_Ele15_IsoVVVL_PFHT450_CaloBTagDeepCSV_4p5);
   fChain->SetBranchAddress("HLT_Ele15_IsoVVVL_PFHT450_PFMET50", &HLT_Ele15_IsoVVVL_PFHT450_PFMET50, &b_HLT_Ele15_IsoVVVL_PFHT450_PFMET50);
   fChain->SetBranchAddress("HLT_Ele15_IsoVVVL_PFHT450", &HLT_Ele15_IsoVVVL_PFHT450, &b_HLT_Ele15_IsoVVVL_PFHT450);
   fChain->SetBranchAddress("HLT_Ele50_IsoVVVL_PFHT450", &HLT_Ele50_IsoVVVL_PFHT450, &b_HLT_Ele50_IsoVVVL_PFHT450);
   fChain->SetBranchAddress("HLT_Ele15_IsoVVVL_PFHT600", &HLT_Ele15_IsoVVVL_PFHT600, &b_HLT_Ele15_IsoVVVL_PFHT600);
   fChain->SetBranchAddress("HLT_Mu4_TrkIsoVVL_DiPFJet90_40_DEta3p5_MJJ750_HTT300_PFMETNoMu60", &HLT_Mu4_TrkIsoVVL_DiPFJet90_40_DEta3p5_MJJ750_HTT300_PFMETNoMu60, &b_HLT_Mu4_TrkIsoVVL_DiPFJet90_40_DEta3p5_MJJ750_HTT300_PFMETNoMu60);
   fChain->SetBranchAddress("HLT_Mu8_TrkIsoVVL_DiPFJet40_DEta3p5_MJJ750_HTT300_PFMETNoMu60", &HLT_Mu8_TrkIsoVVL_DiPFJet40_DEta3p5_MJJ750_HTT300_PFMETNoMu60, &b_HLT_Mu8_TrkIsoVVL_DiPFJet40_DEta3p5_MJJ750_HTT300_PFMETNoMu60);
   fChain->SetBranchAddress("HLT_Mu10_TrkIsoVVL_DiPFJet40_DEta3p5_MJJ750_HTT350_PFMETNoMu60", &HLT_Mu10_TrkIsoVVL_DiPFJet40_DEta3p5_MJJ750_HTT350_PFMETNoMu60, &b_HLT_Mu10_TrkIsoVVL_DiPFJet40_DEta3p5_MJJ750_HTT350_PFMETNoMu60);
   fChain->SetBranchAddress("HLT_Mu15_IsoVVVL_PFHT450_CaloBTagDeepCSV_4p5", &HLT_Mu15_IsoVVVL_PFHT450_CaloBTagDeepCSV_4p5, &b_HLT_Mu15_IsoVVVL_PFHT450_CaloBTagDeepCSV_4p5);
   fChain->SetBranchAddress("HLT_Mu15_IsoVVVL_PFHT450_PFMET50", &HLT_Mu15_IsoVVVL_PFHT450_PFMET50, &b_HLT_Mu15_IsoVVVL_PFHT450_PFMET50);
   fChain->SetBranchAddress("HLT_Mu15_IsoVVVL_PFHT450", &HLT_Mu15_IsoVVVL_PFHT450, &b_HLT_Mu15_IsoVVVL_PFHT450);
   fChain->SetBranchAddress("HLT_Mu50_IsoVVVL_PFHT450", &HLT_Mu50_IsoVVVL_PFHT450, &b_HLT_Mu50_IsoVVVL_PFHT450);
   fChain->SetBranchAddress("HLT_Mu15_IsoVVVL_PFHT600", &HLT_Mu15_IsoVVVL_PFHT600, &b_HLT_Mu15_IsoVVVL_PFHT600);
   fChain->SetBranchAddress("HLT_Mu3er1p5_PFJet100er2p5_PFMET70_PFMHT70_IDTight", &HLT_Mu3er1p5_PFJet100er2p5_PFMET70_PFMHT70_IDTight, &b_HLT_Mu3er1p5_PFJet100er2p5_PFMET70_PFMHT70_IDTight);
   fChain->SetBranchAddress("HLT_Mu3er1p5_PFJet100er2p5_PFMET80_PFMHT80_IDTight", &HLT_Mu3er1p5_PFJet100er2p5_PFMET80_PFMHT80_IDTight, &b_HLT_Mu3er1p5_PFJet100er2p5_PFMET80_PFMHT80_IDTight);
   fChain->SetBranchAddress("HLT_Mu3er1p5_PFJet100er2p5_PFMET90_PFMHT90_IDTight", &HLT_Mu3er1p5_PFJet100er2p5_PFMET90_PFMHT90_IDTight, &b_HLT_Mu3er1p5_PFJet100er2p5_PFMET90_PFMHT90_IDTight);
   fChain->SetBranchAddress("HLT_Mu3er1p5_PFJet100er2p5_PFMET100_PFMHT100_IDTight", &HLT_Mu3er1p5_PFJet100er2p5_PFMET100_PFMHT100_IDTight, &b_HLT_Mu3er1p5_PFJet100er2p5_PFMET100_PFMHT100_IDTight);
   fChain->SetBranchAddress("HLT_Mu3er1p5_PFJet100er2p5_PFMETNoMu70_PFMHTNoMu70_IDTight", &HLT_Mu3er1p5_PFJet100er2p5_PFMETNoMu70_PFMHTNoMu70_IDTight, &b_HLT_Mu3er1p5_PFJet100er2p5_PFMETNoMu70_PFMHTNoMu70_IDTight);
   fChain->SetBranchAddress("HLT_Mu3er1p5_PFJet100er2p5_PFMETNoMu80_PFMHTNoMu80_IDTight", &HLT_Mu3er1p5_PFJet100er2p5_PFMETNoMu80_PFMHTNoMu80_IDTight, &b_HLT_Mu3er1p5_PFJet100er2p5_PFMETNoMu80_PFMHTNoMu80_IDTight);
   fChain->SetBranchAddress("HLT_Mu3er1p5_PFJet100er2p5_PFMETNoMu90_PFMHTNoMu90_IDTight", &HLT_Mu3er1p5_PFJet100er2p5_PFMETNoMu90_PFMHTNoMu90_IDTight, &b_HLT_Mu3er1p5_PFJet100er2p5_PFMETNoMu90_PFMHTNoMu90_IDTight);
   fChain->SetBranchAddress("HLT_Mu3er1p5_PFJet100er2p5_PFMETNoMu100_PFMHTNoMu100_IDTight", &HLT_Mu3er1p5_PFJet100er2p5_PFMETNoMu100_PFMHTNoMu100_IDTight, &b_HLT_Mu3er1p5_PFJet100er2p5_PFMETNoMu100_PFMHTNoMu100_IDTight);
   fChain->SetBranchAddress("HLT_Dimuon10_PsiPrime_Barrel_Seagulls", &HLT_Dimuon10_PsiPrime_Barrel_Seagulls, &b_HLT_Dimuon10_PsiPrime_Barrel_Seagulls);
   fChain->SetBranchAddress("HLT_Dimuon20_Jpsi_Barrel_Seagulls", &HLT_Dimuon20_Jpsi_Barrel_Seagulls, &b_HLT_Dimuon20_Jpsi_Barrel_Seagulls);
   fChain->SetBranchAddress("HLT_Dimuon10_Upsilon_y1p4", &HLT_Dimuon10_Upsilon_y1p4, &b_HLT_Dimuon10_Upsilon_y1p4);
   fChain->SetBranchAddress("HLT_Dimuon12_Upsilon_y1p4", &HLT_Dimuon12_Upsilon_y1p4, &b_HLT_Dimuon12_Upsilon_y1p4);
   fChain->SetBranchAddress("HLT_Dimuon14_Phi_Barrel_Seagulls", &HLT_Dimuon14_Phi_Barrel_Seagulls, &b_HLT_Dimuon14_Phi_Barrel_Seagulls);
   fChain->SetBranchAddress("HLT_Dimuon25_Jpsi", &HLT_Dimuon25_Jpsi, &b_HLT_Dimuon25_Jpsi);
   fChain->SetBranchAddress("HLT_Dimuon14_PsiPrime", &HLT_Dimuon14_PsiPrime, &b_HLT_Dimuon14_PsiPrime);
   fChain->SetBranchAddress("HLT_Dimuon14_PsiPrime_noCorrL1", &HLT_Dimuon14_PsiPrime_noCorrL1, &b_HLT_Dimuon14_PsiPrime_noCorrL1);
   fChain->SetBranchAddress("HLT_Dimuon18_PsiPrime", &HLT_Dimuon18_PsiPrime, &b_HLT_Dimuon18_PsiPrime);
   fChain->SetBranchAddress("HLT_Dimuon18_PsiPrime_noCorrL1", &HLT_Dimuon18_PsiPrime_noCorrL1, &b_HLT_Dimuon18_PsiPrime_noCorrL1);
   fChain->SetBranchAddress("HLT_Dimuon24_Upsilon_noCorrL1", &HLT_Dimuon24_Upsilon_noCorrL1, &b_HLT_Dimuon24_Upsilon_noCorrL1);
   fChain->SetBranchAddress("HLT_Dimuon24_Phi_noCorrL1", &HLT_Dimuon24_Phi_noCorrL1, &b_HLT_Dimuon24_Phi_noCorrL1);
   fChain->SetBranchAddress("HLT_Dimuon25_Jpsi_noCorrL1", &HLT_Dimuon25_Jpsi_noCorrL1, &b_HLT_Dimuon25_Jpsi_noCorrL1);
   fChain->SetBranchAddress("HLT_DiMu4_Ele9_CaloIdL_TrackIdL_DZ_Mass3p8", &HLT_DiMu4_Ele9_CaloIdL_TrackIdL_DZ_Mass3p8, &b_HLT_DiMu4_Ele9_CaloIdL_TrackIdL_DZ_Mass3p8);
   fChain->SetBranchAddress("HLT_DiMu9_Ele9_CaloIdL_TrackIdL_DZ", &HLT_DiMu9_Ele9_CaloIdL_TrackIdL_DZ, &b_HLT_DiMu9_Ele9_CaloIdL_TrackIdL_DZ);
   fChain->SetBranchAddress("HLT_DiMu9_Ele9_CaloIdL_TrackIdL", &HLT_DiMu9_Ele9_CaloIdL_TrackIdL, &b_HLT_DiMu9_Ele9_CaloIdL_TrackIdL);
   fChain->SetBranchAddress("HLT_DoubleIsoMu20_eta2p1", &HLT_DoubleIsoMu20_eta2p1, &b_HLT_DoubleIsoMu20_eta2p1);
   fChain->SetBranchAddress("HLT_TrkMu12_DoubleTrkMu5NoFiltersNoVtx", &HLT_TrkMu12_DoubleTrkMu5NoFiltersNoVtx, &b_HLT_TrkMu12_DoubleTrkMu5NoFiltersNoVtx);
   fChain->SetBranchAddress("HLT_Mu8", &HLT_Mu8, &b_HLT_Mu8);
   fChain->SetBranchAddress("HLT_Mu17", &HLT_Mu17, &b_HLT_Mu17);
   fChain->SetBranchAddress("HLT_Mu19", &HLT_Mu19, &b_HLT_Mu19);
   fChain->SetBranchAddress("HLT_Mu17_Photon30_IsoCaloId", &HLT_Mu17_Photon30_IsoCaloId, &b_HLT_Mu17_Photon30_IsoCaloId);
   fChain->SetBranchAddress("HLT_Ele8_CaloIdL_TrackIdL_IsoVL_PFJet30", &HLT_Ele8_CaloIdL_TrackIdL_IsoVL_PFJet30, &b_HLT_Ele8_CaloIdL_TrackIdL_IsoVL_PFJet30);
   fChain->SetBranchAddress("HLT_Ele12_CaloIdL_TrackIdL_IsoVL_PFJet30", &HLT_Ele12_CaloIdL_TrackIdL_IsoVL_PFJet30, &b_HLT_Ele12_CaloIdL_TrackIdL_IsoVL_PFJet30);
   fChain->SetBranchAddress("HLT_Ele23_CaloIdL_TrackIdL_IsoVL_PFJet30", &HLT_Ele23_CaloIdL_TrackIdL_IsoVL_PFJet30, &b_HLT_Ele23_CaloIdL_TrackIdL_IsoVL_PFJet30);
   fChain->SetBranchAddress("HLT_Ele8_CaloIdM_TrackIdM_PFJet30", &HLT_Ele8_CaloIdM_TrackIdM_PFJet30, &b_HLT_Ele8_CaloIdM_TrackIdM_PFJet30);
   fChain->SetBranchAddress("HLT_Ele17_CaloIdM_TrackIdM_PFJet30", &HLT_Ele17_CaloIdM_TrackIdM_PFJet30, &b_HLT_Ele17_CaloIdM_TrackIdM_PFJet30);
   fChain->SetBranchAddress("HLT_Ele23_CaloIdM_TrackIdM_PFJet30", &HLT_Ele23_CaloIdM_TrackIdM_PFJet30, &b_HLT_Ele23_CaloIdM_TrackIdM_PFJet30);
   fChain->SetBranchAddress("HLT_Ele50_CaloIdVT_GsfTrkIdT_PFJet165", &HLT_Ele50_CaloIdVT_GsfTrkIdT_PFJet165, &b_HLT_Ele50_CaloIdVT_GsfTrkIdT_PFJet165);
   fChain->SetBranchAddress("HLT_Ele115_CaloIdVT_GsfTrkIdT", &HLT_Ele115_CaloIdVT_GsfTrkIdT, &b_HLT_Ele115_CaloIdVT_GsfTrkIdT);
   fChain->SetBranchAddress("HLT_Ele135_CaloIdVT_GsfTrkIdT", &HLT_Ele135_CaloIdVT_GsfTrkIdT, &b_HLT_Ele135_CaloIdVT_GsfTrkIdT);
   fChain->SetBranchAddress("HLT_PFHT330PT30_QuadPFJet_75_60_45_40_TriplePFBTagDeepCSV_4p5", &HLT_PFHT330PT30_QuadPFJet_75_60_45_40_TriplePFBTagDeepCSV_4p5, &b_HLT_PFHT330PT30_QuadPFJet_75_60_45_40_TriplePFBTagDeepCSV_4p5);
   fChain->SetBranchAddress("HLT_PFHT330PT30_QuadPFJet_75_60_45_40", &HLT_PFHT330PT30_QuadPFJet_75_60_45_40, &b_HLT_PFHT330PT30_QuadPFJet_75_60_45_40);
   fChain->SetBranchAddress("HLT_PFHT400_SixPFJet32_DoublePFBTagDeepCSV_2p94", &HLT_PFHT400_SixPFJet32_DoublePFBTagDeepCSV_2p94, &b_HLT_PFHT400_SixPFJet32_DoublePFBTagDeepCSV_2p94);
   fChain->SetBranchAddress("HLT_PFHT400_SixPFJet32", &HLT_PFHT400_SixPFJet32, &b_HLT_PFHT400_SixPFJet32);
   fChain->SetBranchAddress("HLT_PFHT450_SixPFJet36_PFBTagDeepCSV_1p59", &HLT_PFHT450_SixPFJet36_PFBTagDeepCSV_1p59, &b_HLT_PFHT450_SixPFJet36_PFBTagDeepCSV_1p59);
   fChain->SetBranchAddress("HLT_PFHT450_SixPFJet36", &HLT_PFHT450_SixPFJet36, &b_HLT_PFHT450_SixPFJet36);
   fChain->SetBranchAddress("HLT_PFHT400_FivePFJet_100_100_60_30_30", &HLT_PFHT400_FivePFJet_100_100_60_30_30, &b_HLT_PFHT400_FivePFJet_100_100_60_30_30);
   fChain->SetBranchAddress("HLT_PFHT400_FivePFJet_100_100_60_30_30_DoublePFBTagDeepCSV_4p5", &HLT_PFHT400_FivePFJet_100_100_60_30_30_DoublePFBTagDeepCSV_4p5, &b_HLT_PFHT400_FivePFJet_100_100_60_30_30_DoublePFBTagDeepCSV_4p5);
   fChain->SetBranchAddress("HLT_PFHT400_FivePFJet_120_120_60_30_30_DoublePFBTagDeepCSV_4p5", &HLT_PFHT400_FivePFJet_120_120_60_30_30_DoublePFBTagDeepCSV_4p5, &b_HLT_PFHT400_FivePFJet_120_120_60_30_30_DoublePFBTagDeepCSV_4p5);
   fChain->SetBranchAddress("HLT_PFHT350", &HLT_PFHT350, &b_HLT_PFHT350);
   fChain->SetBranchAddress("HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_PFHT350MinPFJet15", &HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_PFHT350MinPFJet15, &b_HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_PFHT350MinPFJet15);
   fChain->SetBranchAddress("HLT_ECALHT800", &HLT_ECALHT800, &b_HLT_ECALHT800);
   fChain->SetBranchAddress("HLT_DiSC30_18_EIso_AND_HE_Mass70", &HLT_DiSC30_18_EIso_AND_HE_Mass70, &b_HLT_DiSC30_18_EIso_AND_HE_Mass70);
   fChain->SetBranchAddress("HLT_Physics", &HLT_Physics, &b_HLT_Physics);
   fChain->SetBranchAddress("HLT_Random", &HLT_Random, &b_HLT_Random);
   fChain->SetBranchAddress("HLT_ZeroBias", &HLT_ZeroBias, &b_HLT_ZeroBias);
   fChain->SetBranchAddress("HLT_ZeroBias_Alignment", &HLT_ZeroBias_Alignment, &b_HLT_ZeroBias_Alignment);
   fChain->SetBranchAddress("HLT_Photon20_HoverELoose", &HLT_Photon20_HoverELoose, &b_HLT_Photon20_HoverELoose);
   fChain->SetBranchAddress("HLT_Photon30_HoverELoose", &HLT_Photon30_HoverELoose, &b_HLT_Photon30_HoverELoose);
   fChain->SetBranchAddress("HLT_EcalCalibration", &HLT_EcalCalibration, &b_HLT_EcalCalibration);
   fChain->SetBranchAddress("HLT_HcalCalibration", &HLT_HcalCalibration, &b_HLT_HcalCalibration);
   fChain->SetBranchAddress("HLT_L1UnpairedBunchBptxMinus", &HLT_L1UnpairedBunchBptxMinus, &b_HLT_L1UnpairedBunchBptxMinus);
   fChain->SetBranchAddress("HLT_L1UnpairedBunchBptxPlus", &HLT_L1UnpairedBunchBptxPlus, &b_HLT_L1UnpairedBunchBptxPlus);
   fChain->SetBranchAddress("HLT_L1NotBptxOR", &HLT_L1NotBptxOR, &b_HLT_L1NotBptxOR);
   fChain->SetBranchAddress("HLT_L1_CDC_SingleMu_3_er1p2_TOP120_DPHI2p618_3p142", &HLT_L1_CDC_SingleMu_3_er1p2_TOP120_DPHI2p618_3p142, &b_HLT_L1_CDC_SingleMu_3_er1p2_TOP120_DPHI2p618_3p142);
   fChain->SetBranchAddress("HLT_CDC_L2cosmic_10_er1p0", &HLT_CDC_L2cosmic_10_er1p0, &b_HLT_CDC_L2cosmic_10_er1p0);
   fChain->SetBranchAddress("HLT_CDC_L2cosmic_5p5_er1p0", &HLT_CDC_L2cosmic_5p5_er1p0, &b_HLT_CDC_L2cosmic_5p5_er1p0);
   fChain->SetBranchAddress("HLT_HcalNZS", &HLT_HcalNZS, &b_HLT_HcalNZS);
   fChain->SetBranchAddress("HLT_HcalPhiSym", &HLT_HcalPhiSym, &b_HLT_HcalPhiSym);
   fChain->SetBranchAddress("HLT_HcalIsolatedbunch", &HLT_HcalIsolatedbunch, &b_HLT_HcalIsolatedbunch);
   fChain->SetBranchAddress("HLT_IsoTrackHB", &HLT_IsoTrackHB, &b_HLT_IsoTrackHB);
   fChain->SetBranchAddress("HLT_IsoTrackHE", &HLT_IsoTrackHE, &b_HLT_IsoTrackHE);
   fChain->SetBranchAddress("HLT_ZeroBias_FirstCollisionAfterAbortGap", &HLT_ZeroBias_FirstCollisionAfterAbortGap, &b_HLT_ZeroBias_FirstCollisionAfterAbortGap);
   fChain->SetBranchAddress("HLT_ZeroBias_IsolatedBunches", &HLT_ZeroBias_IsolatedBunches, &b_HLT_ZeroBias_IsolatedBunches);
   fChain->SetBranchAddress("HLT_ZeroBias_FirstCollisionInTrain", &HLT_ZeroBias_FirstCollisionInTrain, &b_HLT_ZeroBias_FirstCollisionInTrain);
   fChain->SetBranchAddress("HLT_ZeroBias_LastCollisionInTrain", &HLT_ZeroBias_LastCollisionInTrain, &b_HLT_ZeroBias_LastCollisionInTrain);
   fChain->SetBranchAddress("HLT_ZeroBias_FirstBXAfterTrain", &HLT_ZeroBias_FirstBXAfterTrain, &b_HLT_ZeroBias_FirstBXAfterTrain);
   fChain->SetBranchAddress("HLT_Ele16_Ele12_Ele8_CaloIdL_TrackIdL", &HLT_Ele16_Ele12_Ele8_CaloIdL_TrackIdL, &b_HLT_Ele16_Ele12_Ele8_CaloIdL_TrackIdL);
   fChain->SetBranchAddress("HLT_VBF_DoubleMediumDeepTauPFTauHPS20_eta2p1", &HLT_VBF_DoubleMediumDeepTauPFTauHPS20_eta2p1, &b_HLT_VBF_DoubleMediumDeepTauPFTauHPS20_eta2p1);
   fChain->SetBranchAddress("HLT_Photon75_R9Id90_HE10_IsoM_EBOnly_PFJetsMJJ300DEta3", &HLT_Photon75_R9Id90_HE10_IsoM_EBOnly_PFJetsMJJ300DEta3, &b_HLT_Photon75_R9Id90_HE10_IsoM_EBOnly_PFJetsMJJ300DEta3);
   fChain->SetBranchAddress("HLT_PFMET100_PFMHT100_IDTight_PFHT60", &HLT_PFMET100_PFMHT100_IDTight_PFHT60, &b_HLT_PFMET100_PFMHT100_IDTight_PFHT60);
   fChain->SetBranchAddress("HLT_PFMETNoMu100_PFMHTNoMu100_IDTight_PFHT60", &HLT_PFMETNoMu100_PFMHTNoMu100_IDTight_PFHT60, &b_HLT_PFMETNoMu100_PFMHTNoMu100_IDTight_PFHT60);
   fChain->SetBranchAddress("HLT_PFMETTypeOne100_PFMHT100_IDTight_PFHT60", &HLT_PFMETTypeOne100_PFMHT100_IDTight_PFHT60, &b_HLT_PFMETTypeOne100_PFMHT100_IDTight_PFHT60);
   fChain->SetBranchAddress("HLT_Mu18_Mu9_SameSign", &HLT_Mu18_Mu9_SameSign, &b_HLT_Mu18_Mu9_SameSign);
   fChain->SetBranchAddress("HLT_DoubleMu2_Jpsi_DoubleTrk1_Phi1p05", &HLT_DoubleMu2_Jpsi_DoubleTrk1_Phi1p05, &b_HLT_DoubleMu2_Jpsi_DoubleTrk1_Phi1p05);
   fChain->SetBranchAddress("HLT_DoubleMu3_DCA_PFMET50_PFMHT60", &HLT_DoubleMu3_DCA_PFMET50_PFMHT60, &b_HLT_DoubleMu3_DCA_PFMET50_PFMHT60);
   fChain->SetBranchAddress("HLT_TripleMu_5_3_3_Mass3p8_DCA", &HLT_TripleMu_5_3_3_Mass3p8_DCA, &b_HLT_TripleMu_5_3_3_Mass3p8_DCA);
   fChain->SetBranchAddress("HLT_QuadPFJet103_88_75_15_DoublePFBTagDeepCSV_1p3_7p7_VBF1", &HLT_QuadPFJet103_88_75_15_DoublePFBTagDeepCSV_1p3_7p7_VBF1, &b_HLT_QuadPFJet103_88_75_15_DoublePFBTagDeepCSV_1p3_7p7_VBF1);
   fChain->SetBranchAddress("HLT_QuadPFJet111_90_80_15_DoublePFBTagDeepCSV_1p3_7p7_VBF1", &HLT_QuadPFJet111_90_80_15_DoublePFBTagDeepCSV_1p3_7p7_VBF1, &b_HLT_QuadPFJet111_90_80_15_DoublePFBTagDeepCSV_1p3_7p7_VBF1);
   fChain->SetBranchAddress("HLT_QuadPFJet103_88_75_15_PFBTagDeepCSV_1p3_VBF2", &HLT_QuadPFJet103_88_75_15_PFBTagDeepCSV_1p3_VBF2, &b_HLT_QuadPFJet103_88_75_15_PFBTagDeepCSV_1p3_VBF2);
   fChain->SetBranchAddress("HLT_QuadPFJet105_88_76_15_PFBTagDeepCSV_1p3_VBF2", &HLT_QuadPFJet105_88_76_15_PFBTagDeepCSV_1p3_VBF2, &b_HLT_QuadPFJet105_88_76_15_PFBTagDeepCSV_1p3_VBF2);
   fChain->SetBranchAddress("HLT_QuadPFJet111_90_80_15_PFBTagDeepCSV_1p3_VBF2", &HLT_QuadPFJet111_90_80_15_PFBTagDeepCSV_1p3_VBF2, &b_HLT_QuadPFJet111_90_80_15_PFBTagDeepCSV_1p3_VBF2);
   fChain->SetBranchAddress("HLT_QuadPFJet103_88_75_15", &HLT_QuadPFJet103_88_75_15, &b_HLT_QuadPFJet103_88_75_15);
   fChain->SetBranchAddress("HLT_QuadPFJet105_88_76_15", &HLT_QuadPFJet105_88_76_15, &b_HLT_QuadPFJet105_88_76_15);
   fChain->SetBranchAddress("HLT_QuadPFJet111_90_80_15", &HLT_QuadPFJet111_90_80_15, &b_HLT_QuadPFJet111_90_80_15);
   fChain->SetBranchAddress("HLT_AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_p02", &HLT_AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_p02, &b_HLT_AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_p02);
   fChain->SetBranchAddress("HLT_AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_np2", &HLT_AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_np2, &b_HLT_AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_np2);
   fChain->SetBranchAddress("HLT_AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_np4", &HLT_AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_np4, &b_HLT_AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_np4);
   fChain->SetBranchAddress("HLT_Diphoton30_18_R9IdL_AND_HE_AND_IsoCaloId_Mass55", &HLT_Diphoton30_18_R9IdL_AND_HE_AND_IsoCaloId_Mass55, &b_HLT_Diphoton30_18_R9IdL_AND_HE_AND_IsoCaloId_Mass55);
   fChain->SetBranchAddress("HLT_Diphoton30_18_R9IdL_AND_HE_AND_IsoCaloId", &HLT_Diphoton30_18_R9IdL_AND_HE_AND_IsoCaloId, &b_HLT_Diphoton30_18_R9IdL_AND_HE_AND_IsoCaloId);
   fChain->SetBranchAddress("HLT_Mu12_IP6", &HLT_Mu12_IP6, &b_HLT_Mu12_IP6);
   fChain->SetBranchAddress("HLT_QuadPFJet105_88_76_15_DoublePFBTagDeepCSV_1p3_7p7_VBF1", &HLT_QuadPFJet105_88_76_15_DoublePFBTagDeepCSV_1p3_7p7_VBF1, &b_HLT_QuadPFJet105_88_76_15_DoublePFBTagDeepCSV_1p3_7p7_VBF1);
   fChain->SetBranchAddress("HLT_DoubleMediumDeepTauPFTauHPS35_L2NN_eta2p1", &HLT_DoubleMediumDeepTauPFTauHPS35_L2NN_eta2p1, &b_HLT_DoubleMediumDeepTauPFTauHPS35_L2NN_eta2p1);
   fChain->SetBranchAddress("HLT_Ele24_eta2p1_WPTight_Gsf_LooseDeepTauPFTauHPS30_eta2p1_CrossL1", &HLT_Ele24_eta2p1_WPTight_Gsf_LooseDeepTauPFTauHPS30_eta2p1_CrossL1, &b_HLT_Ele24_eta2p1_WPTight_Gsf_LooseDeepTauPFTauHPS30_eta2p1_CrossL1);
   fChain->SetBranchAddress("HLT_IsoMu20_eta2p1_LooseDeepTauPFTauHPS27_eta2p1_CrossL1", &HLT_IsoMu20_eta2p1_LooseDeepTauPFTauHPS27_eta2p1_CrossL1, &b_HLT_IsoMu20_eta2p1_LooseDeepTauPFTauHPS27_eta2p1_CrossL1);
   fChain->SetBranchAddress("HLT_IsoMu24_eta2p1_LooseDeepTauPFTauHPS180_eta2p1", &HLT_IsoMu24_eta2p1_LooseDeepTauPFTauHPS180_eta2p1, &b_HLT_IsoMu24_eta2p1_LooseDeepTauPFTauHPS180_eta2p1);
   fChain->SetBranchAddress("HLT_IsoMu24_eta2p1_LooseDeepTauPFTauHPS30_eta2p1_CrossL1", &HLT_IsoMu24_eta2p1_LooseDeepTauPFTauHPS30_eta2p1_CrossL1, &b_HLT_IsoMu24_eta2p1_LooseDeepTauPFTauHPS30_eta2p1_CrossL1);
   fChain->SetBranchAddress("HLT_IsoMu24_eta2p1_MediumDeepTauPFTauHPS35_L2NN_eta2p1_CrossL1", &HLT_IsoMu24_eta2p1_MediumDeepTauPFTauHPS35_L2NN_eta2p1_CrossL1, &b_HLT_IsoMu24_eta2p1_MediumDeepTauPFTauHPS35_L2NN_eta2p1_CrossL1);
   fChain->SetBranchAddress("HLT_LooseDeepTauPFTauHPS180_L2NN_eta2p1", &HLT_LooseDeepTauPFTauHPS180_L2NN_eta2p1, &b_HLT_LooseDeepTauPFTauHPS180_L2NN_eta2p1);
   fChain->SetBranchAddress("HLT_PFHT330PT30_QuadPFJet_75_60_45_40_TriplePFBTagDeepJet_4p5", &HLT_PFHT330PT30_QuadPFJet_75_60_45_40_TriplePFBTagDeepJet_4p5, &b_HLT_PFHT330PT30_QuadPFJet_75_60_45_40_TriplePFBTagDeepJet_4p5);
   fChain->SetBranchAddress("HLT_PFHT400_FivePFJet_100_100_60_30_30_DoublePFBTagDeepJet_4p5", &HLT_PFHT400_FivePFJet_100_100_60_30_30_DoublePFBTagDeepJet_4p5, &b_HLT_PFHT400_FivePFJet_100_100_60_30_30_DoublePFBTagDeepJet_4p5);
   fChain->SetBranchAddress("HLT_PFHT400_FivePFJet_120_120_60_30_30_DoublePFBTagDeepJet_4p5", &HLT_PFHT400_FivePFJet_120_120_60_30_30_DoublePFBTagDeepJet_4p5, &b_HLT_PFHT400_FivePFJet_120_120_60_30_30_DoublePFBTagDeepJet_4p5);
   fChain->SetBranchAddress("HLT_PFHT400_SixPFJet32_DoublePFBTagDeepJet_2p94", &HLT_PFHT400_SixPFJet32_DoublePFBTagDeepJet_2p94, &b_HLT_PFHT400_SixPFJet32_DoublePFBTagDeepJet_2p94);
   fChain->SetBranchAddress("HLT_PFHT450_SixPFJet36_PFBTagDeepJet_1p59", &HLT_PFHT450_SixPFJet36_PFBTagDeepJet_1p59, &b_HLT_PFHT450_SixPFJet36_PFBTagDeepJet_1p59);
   fChain->SetBranchAddress("HLT_QuadPFJet103_88_75_15_DoublePFBTagDeepJet_1p3_7p7_VBF1", &HLT_QuadPFJet103_88_75_15_DoublePFBTagDeepJet_1p3_7p7_VBF1, &b_HLT_QuadPFJet103_88_75_15_DoublePFBTagDeepJet_1p3_7p7_VBF1);
   fChain->SetBranchAddress("HLT_QuadPFJet103_88_75_15_PFBTagDeepJet_1p3_VBF2", &HLT_QuadPFJet103_88_75_15_PFBTagDeepJet_1p3_VBF2, &b_HLT_QuadPFJet103_88_75_15_PFBTagDeepJet_1p3_VBF2);
   fChain->SetBranchAddress("HLT_QuadPFJet105_88_76_15_DoublePFBTagDeepJet_1p3_7p7_VBF1", &HLT_QuadPFJet105_88_76_15_DoublePFBTagDeepJet_1p3_7p7_VBF1, &b_HLT_QuadPFJet105_88_76_15_DoublePFBTagDeepJet_1p3_7p7_VBF1);
   fChain->SetBranchAddress("HLT_QuadPFJet105_88_76_15_PFBTagDeepJet_1p3_VBF2", &HLT_QuadPFJet105_88_76_15_PFBTagDeepJet_1p3_VBF2, &b_HLT_QuadPFJet105_88_76_15_PFBTagDeepJet_1p3_VBF2);
   fChain->SetBranchAddress("HLT_QuadPFJet111_90_80_15_DoublePFBTagDeepJet_1p3_7p7_VBF1", &HLT_QuadPFJet111_90_80_15_DoublePFBTagDeepJet_1p3_7p7_VBF1, &b_HLT_QuadPFJet111_90_80_15_DoublePFBTagDeepJet_1p3_7p7_VBF1);
   fChain->SetBranchAddress("HLT_QuadPFJet111_90_80_15_PFBTagDeepJet_1p3_VBF2", &HLT_QuadPFJet111_90_80_15_PFBTagDeepJet_1p3_VBF2, &b_HLT_QuadPFJet111_90_80_15_PFBTagDeepJet_1p3_VBF2);
   fChain->SetBranchAddress("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_PFDiJet30_PFBtagDeepJet_1p5", &HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_PFDiJet30_PFBtagDeepJet_1p5, &b_HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_PFDiJet30_PFBtagDeepJet_1p5);
   fChain->SetBranchAddress("HLT_QuadPFJet70_50_40_30", &HLT_QuadPFJet70_50_40_30, &b_HLT_QuadPFJet70_50_40_30);
   fChain->SetBranchAddress("HLT_QuadPFJet70_50_40_30_PFBTagParticleNet_2BTagSum0p65", &HLT_QuadPFJet70_50_40_30_PFBTagParticleNet_2BTagSum0p65, &b_HLT_QuadPFJet70_50_40_30_PFBTagParticleNet_2BTagSum0p65);
   fChain->SetBranchAddress("HLT_QuadPFJet70_50_40_35_PFBTagParticleNet_2BTagSum0p65", &HLT_QuadPFJet70_50_40_35_PFBTagParticleNet_2BTagSum0p65, &b_HLT_QuadPFJet70_50_40_35_PFBTagParticleNet_2BTagSum0p65);
   fChain->SetBranchAddress("HLT_QuadPFJet70_50_45_35_PFBTagParticleNet_2BTagSum0p65", &HLT_QuadPFJet70_50_45_35_PFBTagParticleNet_2BTagSum0p65, &b_HLT_QuadPFJet70_50_45_35_PFBTagParticleNet_2BTagSum0p65);
   fChain->SetBranchAddress("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_PFDiJet30_PFBTagParticleNet_2BTagSum0p65", &HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_PFDiJet30_PFBTagParticleNet_2BTagSum0p65, &b_HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_PFDiJet30_PFBTagParticleNet_2BTagSum0p65);
   fChain->SetBranchAddress("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_QuadPFJet70_50_40_30", &HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_QuadPFJet70_50_40_30, &b_HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_QuadPFJet70_50_40_30);
   fChain->SetBranchAddress("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_QuadPFJet70_50_40_30_PFBTagParticleNet_2BTagSum0p65", &HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_QuadPFJet70_50_40_30_PFBTagParticleNet_2BTagSum0p65, &b_HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_QuadPFJet70_50_40_30_PFBTagParticleNet_2BTagSum0p65);
   fChain->SetBranchAddress("HLT_AK8PFJet230_SoftDropMass40", &HLT_AK8PFJet230_SoftDropMass40, &b_HLT_AK8PFJet230_SoftDropMass40);
   fChain->SetBranchAddress("HLT_AK8PFJet230_SoftDropMass40_PFAK8ParticleNetBB0p35", &HLT_AK8PFJet230_SoftDropMass40_PFAK8ParticleNetBB0p35, &b_HLT_AK8PFJet230_SoftDropMass40_PFAK8ParticleNetBB0p35);
   fChain->SetBranchAddress("HLT_AK8PFJet250_SoftDropMass40_PFAK8ParticleNetBB0p35", &HLT_AK8PFJet250_SoftDropMass40_PFAK8ParticleNetBB0p35, &b_HLT_AK8PFJet250_SoftDropMass40_PFAK8ParticleNetBB0p35);
   fChain->SetBranchAddress("HLT_AK8PFJet275_SoftDropMass40_PFAK8ParticleNetBB0p35", &HLT_AK8PFJet275_SoftDropMass40_PFAK8ParticleNetBB0p35, &b_HLT_AK8PFJet275_SoftDropMass40_PFAK8ParticleNetBB0p35);
   fChain->SetBranchAddress("HLT_AK8PFJet400_SoftDropMass40", &HLT_AK8PFJet400_SoftDropMass40, &b_HLT_AK8PFJet400_SoftDropMass40);
   fChain->SetBranchAddress("HLT_AK8PFJet425_SoftDropMass40", &HLT_AK8PFJet425_SoftDropMass40, &b_HLT_AK8PFJet425_SoftDropMass40);
   fChain->SetBranchAddress("HLT_AK8PFJet450_SoftDropMass40", &HLT_AK8PFJet450_SoftDropMass40, &b_HLT_AK8PFJet450_SoftDropMass40);
   fChain->SetBranchAddress("HLT_AK8PFJet230_SoftDropMass40_PFAK8ParticleNetTauTau0p30", &HLT_AK8PFJet230_SoftDropMass40_PFAK8ParticleNetTauTau0p30, &b_HLT_AK8PFJet230_SoftDropMass40_PFAK8ParticleNetTauTau0p30);
   fChain->SetBranchAddress("HLT_AK8PFJet250_SoftDropMass40_PFAK8ParticleNetTauTau0p30", &HLT_AK8PFJet250_SoftDropMass40_PFAK8ParticleNetTauTau0p30, &b_HLT_AK8PFJet250_SoftDropMass40_PFAK8ParticleNetTauTau0p30);
   fChain->SetBranchAddress("HLT_AK8PFJet275_SoftDropMass40_PFAK8ParticleNetTauTau0p30", &HLT_AK8PFJet275_SoftDropMass40_PFAK8ParticleNetTauTau0p30, &b_HLT_AK8PFJet275_SoftDropMass40_PFAK8ParticleNetTauTau0p30);
   fChain->SetBranchAddress("HLT_IsoMu50_AK8PFJet230_SoftDropMass40", &HLT_IsoMu50_AK8PFJet230_SoftDropMass40, &b_HLT_IsoMu50_AK8PFJet230_SoftDropMass40);
   fChain->SetBranchAddress("HLT_IsoMu50_AK8PFJet230_SoftDropMass40_PFAK8ParticleNetBB0p35", &HLT_IsoMu50_AK8PFJet230_SoftDropMass40_PFAK8ParticleNetBB0p35, &b_HLT_IsoMu50_AK8PFJet230_SoftDropMass40_PFAK8ParticleNetBB0p35);
   fChain->SetBranchAddress("HLT_Ele50_CaloIdVT_GsfTrkIdT_AK8PFJet230_SoftDropMass40", &HLT_Ele50_CaloIdVT_GsfTrkIdT_AK8PFJet230_SoftDropMass40, &b_HLT_Ele50_CaloIdVT_GsfTrkIdT_AK8PFJet230_SoftDropMass40);
   fChain->SetBranchAddress("HLT_Ele50_CaloIdVT_GsfTrkIdT_AK8PFJet230_SoftDropMass40_PFAK8ParticleNetBB0p35", &HLT_Ele50_CaloIdVT_GsfTrkIdT_AK8PFJet230_SoftDropMass40_PFAK8ParticleNetBB0p35, &b_HLT_Ele50_CaloIdVT_GsfTrkIdT_AK8PFJet230_SoftDropMass40_PFAK8ParticleNetBB0p35);
   fChain->SetBranchAddress("HLT_DoubleMediumDeepTauPFTauHPS30_L2NN_eta2p1_PFJet60", &HLT_DoubleMediumDeepTauPFTauHPS30_L2NN_eta2p1_PFJet60, &b_HLT_DoubleMediumDeepTauPFTauHPS30_L2NN_eta2p1_PFJet60);
   fChain->SetBranchAddress("HLT_DoubleMediumDeepTauPFTauHPS30_L2NN_eta2p1_PFJet75", &HLT_DoubleMediumDeepTauPFTauHPS30_L2NN_eta2p1_PFJet75, &b_HLT_DoubleMediumDeepTauPFTauHPS30_L2NN_eta2p1_PFJet75);
   fChain->SetBranchAddress("HLT_IsoMu24_eta2p1_MediumDeepTauPFTauHPS30_L2NN_eta2p1_CrossL1", &HLT_IsoMu24_eta2p1_MediumDeepTauPFTauHPS30_L2NN_eta2p1_CrossL1, &b_HLT_IsoMu24_eta2p1_MediumDeepTauPFTauHPS30_L2NN_eta2p1_CrossL1);
   fChain->SetBranchAddress("HLT_IsoMu24_eta2p1_MediumDeepTauPFTauHPS30_L2NN_eta2p1_PFJet60_CrossL1", &HLT_IsoMu24_eta2p1_MediumDeepTauPFTauHPS30_L2NN_eta2p1_PFJet60_CrossL1, &b_HLT_IsoMu24_eta2p1_MediumDeepTauPFTauHPS30_L2NN_eta2p1_PFJet60_CrossL1);
   fChain->SetBranchAddress("HLT_IsoMu24_eta2p1_MediumDeepTauPFTauHPS30_L2NN_eta2p1_PFJet75_CrossL1", &HLT_IsoMu24_eta2p1_MediumDeepTauPFTauHPS30_L2NN_eta2p1_PFJet75_CrossL1, &b_HLT_IsoMu24_eta2p1_MediumDeepTauPFTauHPS30_L2NN_eta2p1_PFJet75_CrossL1);
   fChain->SetBranchAddress("HLT_DoublePFJets40_Mass500_MediumDeepTauPFTauHPS45_L2NN_MediumDeepTauPFTauHPS20_eta2p1", &HLT_DoublePFJets40_Mass500_MediumDeepTauPFTauHPS45_L2NN_MediumDeepTauPFTauHPS20_eta2p1, &b_HLT_DoublePFJets40_Mass500_MediumDeepTauPFTauHPS45_L2NN_MediumDeepTauPFTauHPS20_eta2p1);
   fChain->SetBranchAddress("HLT_IsoMu24_eta2p1_MediumDeepTauPFTauHPS20_eta2p1_SingleL1", &HLT_IsoMu24_eta2p1_MediumDeepTauPFTauHPS20_eta2p1_SingleL1, &b_HLT_IsoMu24_eta2p1_MediumDeepTauPFTauHPS20_eta2p1_SingleL1);
   fChain->SetBranchAddress("HLT_IsoMu24_eta2p1_MediumDeepTauPFTauHPS45_L2NN_eta2p1_CrossL1", &HLT_IsoMu24_eta2p1_MediumDeepTauPFTauHPS45_L2NN_eta2p1_CrossL1, &b_HLT_IsoMu24_eta2p1_MediumDeepTauPFTauHPS45_L2NN_eta2p1_CrossL1);
   fChain->SetBranchAddress("HLT_DoubleL2Mu10NoVtx_2Cha_VetoL3Mu0DxyMax1cm", &HLT_DoubleL2Mu10NoVtx_2Cha_VetoL3Mu0DxyMax1cm, &b_HLT_DoubleL2Mu10NoVtx_2Cha_VetoL3Mu0DxyMax1cm);
   fChain->SetBranchAddress("HLT_DoubleL2Mu12NoVtx_2Cha_VetoL3Mu0DxyMax1cm", &HLT_DoubleL2Mu12NoVtx_2Cha_VetoL3Mu0DxyMax1cm, &b_HLT_DoubleL2Mu12NoVtx_2Cha_VetoL3Mu0DxyMax1cm);
   fChain->SetBranchAddress("HLT_DoubleL2Mu14NoVtx_2Cha_VetoL3Mu0DxyMax1cm", &HLT_DoubleL2Mu14NoVtx_2Cha_VetoL3Mu0DxyMax1cm, &b_HLT_DoubleL2Mu14NoVtx_2Cha_VetoL3Mu0DxyMax1cm);
   fChain->SetBranchAddress("HLT_DoubleL3Mu16_10NoVtx_DxyMin0p01cm", &HLT_DoubleL3Mu16_10NoVtx_DxyMin0p01cm, &b_HLT_DoubleL3Mu16_10NoVtx_DxyMin0p01cm);
   fChain->SetBranchAddress("HLT_DoubleL3Mu18_10NoVtx_DxyMin0p01cm", &HLT_DoubleL3Mu18_10NoVtx_DxyMin0p01cm, &b_HLT_DoubleL3Mu18_10NoVtx_DxyMin0p01cm);
   fChain->SetBranchAddress("HLT_DoubleL3Mu20_10NoVtx_DxyMin0p01cm", &HLT_DoubleL3Mu20_10NoVtx_DxyMin0p01cm, &b_HLT_DoubleL3Mu20_10NoVtx_DxyMin0p01cm);
   fChain->SetBranchAddress("HLT_L2Mu10NoVtx_2Cha", &HLT_L2Mu10NoVtx_2Cha, &b_HLT_L2Mu10NoVtx_2Cha);
   fChain->SetBranchAddress("HLT_L2Mu10NoVtx_2Cha_VetoL3Mu0DxyMax1cm", &HLT_L2Mu10NoVtx_2Cha_VetoL3Mu0DxyMax1cm, &b_HLT_L2Mu10NoVtx_2Cha_VetoL3Mu0DxyMax1cm);
   fChain->SetBranchAddress("HLT_L3Mu10NoVtx", &HLT_L3Mu10NoVtx, &b_HLT_L3Mu10NoVtx);
   fChain->SetBranchAddress("HLT_L3Mu10NoVtx_DxyMin0p01cm", &HLT_L3Mu10NoVtx_DxyMin0p01cm, &b_HLT_L3Mu10NoVtx_DxyMin0p01cm);
   fChain->SetBranchAddress("HLT_DoubleL2Mu_L3Mu16NoVtx_VetoL3Mu0DxyMax0p1cm", &HLT_DoubleL2Mu_L3Mu16NoVtx_VetoL3Mu0DxyMax0p1cm, &b_HLT_DoubleL2Mu_L3Mu16NoVtx_VetoL3Mu0DxyMax0p1cm);
   fChain->SetBranchAddress("HLT_DoubleL2Mu_L3Mu18NoVtx_VetoL3Mu0DxyMax0p1cm", &HLT_DoubleL2Mu_L3Mu18NoVtx_VetoL3Mu0DxyMax0p1cm, &b_HLT_DoubleL2Mu_L3Mu18NoVtx_VetoL3Mu0DxyMax0p1cm);
   fChain->SetBranchAddress("HLT_DoubleL2Mu10NoVtx_2Cha_CosmicSeed_VetoL3Mu0DxyMax1cm", &HLT_DoubleL2Mu10NoVtx_2Cha_CosmicSeed_VetoL3Mu0DxyMax1cm, &b_HLT_DoubleL2Mu10NoVtx_2Cha_CosmicSeed_VetoL3Mu0DxyMax1cm);
   fChain->SetBranchAddress("HLT_DoubleL2Mu12NoVtx_2Cha_CosmicSeed_VetoL3Mu0DxyMax1cm", &HLT_DoubleL2Mu12NoVtx_2Cha_CosmicSeed_VetoL3Mu0DxyMax1cm, &b_HLT_DoubleL2Mu12NoVtx_2Cha_CosmicSeed_VetoL3Mu0DxyMax1cm);
   fChain->SetBranchAddress("HLT_L2Mu10NoVtx_2Cha_CosmicSeed", &HLT_L2Mu10NoVtx_2Cha_CosmicSeed, &b_HLT_L2Mu10NoVtx_2Cha_CosmicSeed);
   fChain->SetBranchAddress("HLT_L2Mu10NoVtx_2Cha_CosmicSeed_VetoL3Mu0DxyMax1cm", &HLT_L2Mu10NoVtx_2Cha_CosmicSeed_VetoL3Mu0DxyMax1cm, &b_HLT_L2Mu10NoVtx_2Cha_CosmicSeed_VetoL3Mu0DxyMax1cm);
   fChain->SetBranchAddress("HLT_DoubleL3dTksMu16_10NoVtx_DxyMin0p01cm", &HLT_DoubleL3dTksMu16_10NoVtx_DxyMin0p01cm, &b_HLT_DoubleL3dTksMu16_10NoVtx_DxyMin0p01cm);
   fChain->SetBranchAddress("HLT_L3dTksMu10_NoVtx_DxyMin0p01cm", &HLT_L3dTksMu10_NoVtx_DxyMin0p01cm, &b_HLT_L3dTksMu10_NoVtx_DxyMin0p01cm);
   fChain->SetBranchAddress("HLT_Mu20NoFiltersNoVtxDisplaced_Photon20_CaloCustomId", &HLT_Mu20NoFiltersNoVtxDisplaced_Photon20_CaloCustomId, &b_HLT_Mu20NoFiltersNoVtxDisplaced_Photon20_CaloCustomId);
   fChain->SetBranchAddress("HLT_DoubleMediumChargedIsoDisplacedPFTauHPS32_Trk1_eta2p1", &HLT_DoubleMediumChargedIsoDisplacedPFTauHPS32_Trk1_eta2p1, &b_HLT_DoubleMediumChargedIsoDisplacedPFTauHPS32_Trk1_eta2p1);
   fChain->SetBranchAddress("HLT_HT430_DelayedJet40_DoubleDelay0p5nsTrackless", &HLT_HT430_DelayedJet40_DoubleDelay0p5nsTrackless, &b_HLT_HT430_DelayedJet40_DoubleDelay0p5nsTrackless);
   fChain->SetBranchAddress("HLT_HT430_DelayedJet40_DoubleDelay1nsInclusive", &HLT_HT430_DelayedJet40_DoubleDelay1nsInclusive, &b_HLT_HT430_DelayedJet40_DoubleDelay1nsInclusive);
   fChain->SetBranchAddress("HLT_HT430_DelayedJet40_SingleDelay1nsTrackless", &HLT_HT430_DelayedJet40_SingleDelay1nsTrackless, &b_HLT_HT430_DelayedJet40_SingleDelay1nsTrackless);
   fChain->SetBranchAddress("HLT_HT430_DelayedJet40_SingleDelay2nsInclusive", &HLT_HT430_DelayedJet40_SingleDelay2nsInclusive, &b_HLT_HT430_DelayedJet40_SingleDelay2nsInclusive);
   fChain->SetBranchAddress("HLT_L1Mu6HT240", &HLT_L1Mu6HT240, &b_HLT_L1Mu6HT240);
   fChain->SetBranchAddress("HLT_Mu6HT240_DisplacedDijet30_Inclusive0PtrkShortSig5", &HLT_Mu6HT240_DisplacedDijet30_Inclusive0PtrkShortSig5, &b_HLT_Mu6HT240_DisplacedDijet30_Inclusive0PtrkShortSig5);
   fChain->SetBranchAddress("HLT_Mu6HT240_DisplacedDijet30_Inclusive1PtrkShortSig5_DisplacedLoose", &HLT_Mu6HT240_DisplacedDijet30_Inclusive1PtrkShortSig5_DisplacedLoose, &b_HLT_Mu6HT240_DisplacedDijet30_Inclusive1PtrkShortSig5_DisplacedLoose);
   fChain->SetBranchAddress("HLT_Mu6HT240_DisplacedDijet35_Inclusive0PtrkShortSig5", &HLT_Mu6HT240_DisplacedDijet35_Inclusive0PtrkShortSig5, &b_HLT_Mu6HT240_DisplacedDijet35_Inclusive0PtrkShortSig5);
   fChain->SetBranchAddress("HLT_Mu6HT240_DisplacedDijet35_Inclusive1PtrkShortSig5_DisplacedLoose", &HLT_Mu6HT240_DisplacedDijet35_Inclusive1PtrkShortSig5_DisplacedLoose, &b_HLT_Mu6HT240_DisplacedDijet35_Inclusive1PtrkShortSig5_DisplacedLoose);
   fChain->SetBranchAddress("HLT_Mu6HT240_DisplacedDijet40_Inclusive0PtrkShortSig5", &HLT_Mu6HT240_DisplacedDijet40_Inclusive0PtrkShortSig5, &b_HLT_Mu6HT240_DisplacedDijet40_Inclusive0PtrkShortSig5);
   fChain->SetBranchAddress("HLT_Mu6HT240_DisplacedDijet40_Inclusive1PtrkShortSig5_DisplacedLoose", &HLT_Mu6HT240_DisplacedDijet40_Inclusive1PtrkShortSig5_DisplacedLoose, &b_HLT_Mu6HT240_DisplacedDijet40_Inclusive1PtrkShortSig5_DisplacedLoose);
   fChain->SetBranchAddress("HLT_HT430_DisplacedDijet30_Inclusive1PtrkShortSig5", &HLT_HT430_DisplacedDijet30_Inclusive1PtrkShortSig5, &b_HLT_HT430_DisplacedDijet30_Inclusive1PtrkShortSig5);
   fChain->SetBranchAddress("HLT_HT430_DisplacedDijet35_Inclusive1PtrkShortSig5", &HLT_HT430_DisplacedDijet35_Inclusive1PtrkShortSig5, &b_HLT_HT430_DisplacedDijet35_Inclusive1PtrkShortSig5);
   fChain->SetBranchAddress("HLT_HT430_DisplacedDijet40_Inclusive1PtrkShortSig5", &HLT_HT430_DisplacedDijet40_Inclusive1PtrkShortSig5, &b_HLT_HT430_DisplacedDijet40_Inclusive1PtrkShortSig5);
   fChain->SetBranchAddress("HLT_CaloMET60_DTCluster50", &HLT_CaloMET60_DTCluster50, &b_HLT_CaloMET60_DTCluster50);
   fChain->SetBranchAddress("HLT_CaloMET60_DTClusterNoMB1S50", &HLT_CaloMET60_DTClusterNoMB1S50, &b_HLT_CaloMET60_DTClusterNoMB1S50);
   fChain->SetBranchAddress("HLT_L1MET_DTCluster50", &HLT_L1MET_DTCluster50, &b_HLT_L1MET_DTCluster50);
   fChain->SetBranchAddress("HLT_L1MET_DTClusterNoMB1S50", &HLT_L1MET_DTClusterNoMB1S50, &b_HLT_L1MET_DTClusterNoMB1S50);
   fChain->SetBranchAddress("HLT_CscCluster_Loose", &HLT_CscCluster_Loose, &b_HLT_CscCluster_Loose);
   fChain->SetBranchAddress("HLT_CscCluster_Medium", &HLT_CscCluster_Medium, &b_HLT_CscCluster_Medium);
   fChain->SetBranchAddress("HLT_CscCluster_Tight", &HLT_CscCluster_Tight, &b_HLT_CscCluster_Tight);
   fChain->SetBranchAddress("HLT_L1CSCShower_DTCluster50", &HLT_L1CSCShower_DTCluster50, &b_HLT_L1CSCShower_DTCluster50);
   fChain->SetBranchAddress("HLT_L1CSCShower_DTCluster75", &HLT_L1CSCShower_DTCluster75, &b_HLT_L1CSCShower_DTCluster75);
   fChain->SetBranchAddress("HLT_DoubleCscCluster75", &HLT_DoubleCscCluster75, &b_HLT_DoubleCscCluster75);
   fChain->SetBranchAddress("HLT_DoubleCscCluster100", &HLT_DoubleCscCluster100, &b_HLT_DoubleCscCluster100);
   fChain->SetBranchAddress("HLT_PFMET105_IsoTrk50", &HLT_PFMET105_IsoTrk50, &b_HLT_PFMET105_IsoTrk50);
   fChain->SetBranchAddress("HLT_PFMET110_PFJet100", &HLT_PFMET110_PFJet100, &b_HLT_PFMET110_PFJet100);
   fChain->SetBranchAddress("HLT_HT170_L1SingleLLPJet_DisplacedDijet40_DisplacedTrack", &HLT_HT170_L1SingleLLPJet_DisplacedDijet40_DisplacedTrack, &b_HLT_HT170_L1SingleLLPJet_DisplacedDijet40_DisplacedTrack);
   fChain->SetBranchAddress("HLT_HT200_L1SingleLLPJet_DisplacedDijet40_DisplacedTrack", &HLT_HT200_L1SingleLLPJet_DisplacedDijet40_DisplacedTrack, &b_HLT_HT200_L1SingleLLPJet_DisplacedDijet40_DisplacedTrack);
   fChain->SetBranchAddress("HLT_HT200_L1SingleLLPJet_DisplacedDijet60_DisplacedTrack", &HLT_HT200_L1SingleLLPJet_DisplacedDijet60_DisplacedTrack, &b_HLT_HT200_L1SingleLLPJet_DisplacedDijet60_DisplacedTrack);
   fChain->SetBranchAddress("HLT_HT270_L1SingleLLPJet_DisplacedDijet40_DisplacedTrack", &HLT_HT270_L1SingleLLPJet_DisplacedDijet40_DisplacedTrack, &b_HLT_HT270_L1SingleLLPJet_DisplacedDijet40_DisplacedTrack);
   fChain->SetBranchAddress("HLT_HT320_L1SingleLLPJet_DisplacedDijet60_Inclusive", &HLT_HT320_L1SingleLLPJet_DisplacedDijet60_Inclusive, &b_HLT_HT320_L1SingleLLPJet_DisplacedDijet60_Inclusive);
   fChain->SetBranchAddress("HLT_HT420_L1SingleLLPJet_DisplacedDijet60_Inclusive", &HLT_HT420_L1SingleLLPJet_DisplacedDijet60_Inclusive, &b_HLT_HT420_L1SingleLLPJet_DisplacedDijet60_Inclusive);
   fChain->SetBranchAddress("HLT_HT200_L1SingleLLPJet_DelayedJet40_SingleDelay1nsTrackless", &HLT_HT200_L1SingleLLPJet_DelayedJet40_SingleDelay1nsTrackless, &b_HLT_HT200_L1SingleLLPJet_DelayedJet40_SingleDelay1nsTrackless);
   fChain->SetBranchAddress("HLT_HT200_L1SingleLLPJet_DelayedJet40_SingleDelay2nsInclusive", &HLT_HT200_L1SingleLLPJet_DelayedJet40_SingleDelay2nsInclusive, &b_HLT_HT200_L1SingleLLPJet_DelayedJet40_SingleDelay2nsInclusive);
   fChain->SetBranchAddress("HLT_HT200_L1SingleLLPJet_DelayedJet40_DoubleDelay0p5nsTrackless", &HLT_HT200_L1SingleLLPJet_DelayedJet40_DoubleDelay0p5nsTrackless, &b_HLT_HT200_L1SingleLLPJet_DelayedJet40_DoubleDelay0p5nsTrackless);
   fChain->SetBranchAddress("HLT_HT200_L1SingleLLPJet_DelayedJet40_DoubleDelay1nsInclusive", &HLT_HT200_L1SingleLLPJet_DelayedJet40_DoubleDelay1nsInclusive, &b_HLT_HT200_L1SingleLLPJet_DelayedJet40_DoubleDelay1nsInclusive);
   fChain->SetBranchAddress("HLT_HT200_L1SingleLLPJet_DisplacedDijet30_Inclusive1PtrkShortSig5", &HLT_HT200_L1SingleLLPJet_DisplacedDijet30_Inclusive1PtrkShortSig5, &b_HLT_HT200_L1SingleLLPJet_DisplacedDijet30_Inclusive1PtrkShortSig5);
   fChain->SetBranchAddress("HLT_HT200_L1SingleLLPJet_DisplacedDijet35_Inclusive1PtrkShortSig5", &HLT_HT200_L1SingleLLPJet_DisplacedDijet35_Inclusive1PtrkShortSig5, &b_HLT_HT200_L1SingleLLPJet_DisplacedDijet35_Inclusive1PtrkShortSig5);
   fChain->SetBranchAddress("HLT_HT200_L1SingleLLPJet_DisplacedDijet40_Inclusive1PtrkShortSig5", &HLT_HT200_L1SingleLLPJet_DisplacedDijet40_Inclusive1PtrkShortSig5, &b_HLT_HT200_L1SingleLLPJet_DisplacedDijet40_Inclusive1PtrkShortSig5);
   fChain->SetBranchAddress("HLT_DiPhoton10Time1p4ns", &HLT_DiPhoton10Time1p4ns, &b_HLT_DiPhoton10Time1p4ns);
   fChain->SetBranchAddress("HLT_DiPhoton10Time1p6ns", &HLT_DiPhoton10Time1p6ns, &b_HLT_DiPhoton10Time1p6ns);
   fChain->SetBranchAddress("HLT_DiPhoton10Time1p8ns", &HLT_DiPhoton10Time1p8ns, &b_HLT_DiPhoton10Time1p8ns);
   fChain->SetBranchAddress("HLT_DiPhoton10Time2ns", &HLT_DiPhoton10Time2ns, &b_HLT_DiPhoton10Time2ns);
   fChain->SetBranchAddress("HLT_DiPhoton10sminlt0p1", &HLT_DiPhoton10sminlt0p1, &b_HLT_DiPhoton10sminlt0p1);
   fChain->SetBranchAddress("HLT_DiPhoton10sminlt0p12", &HLT_DiPhoton10sminlt0p12, &b_HLT_DiPhoton10sminlt0p12);
   fChain->SetBranchAddress("HLT_DiPhoton10_CaloIdL", &HLT_DiPhoton10_CaloIdL, &b_HLT_DiPhoton10_CaloIdL);
   fChain->SetBranchAddress("HLT_DoubleEle4_eta1p22_mMax6", &HLT_DoubleEle4_eta1p22_mMax6, &b_HLT_DoubleEle4_eta1p22_mMax6);
   fChain->SetBranchAddress("HLT_DoubleEle4p5_eta1p22_mMax6", &HLT_DoubleEle4p5_eta1p22_mMax6, &b_HLT_DoubleEle4p5_eta1p22_mMax6);
   fChain->SetBranchAddress("HLT_DoubleEle5_eta1p22_mMax6", &HLT_DoubleEle5_eta1p22_mMax6, &b_HLT_DoubleEle5_eta1p22_mMax6);
   fChain->SetBranchAddress("HLT_DoubleEle5p5_eta1p22_mMax6", &HLT_DoubleEle5p5_eta1p22_mMax6, &b_HLT_DoubleEle5p5_eta1p22_mMax6);
   fChain->SetBranchAddress("HLT_DoubleEle6_eta1p22_mMax6", &HLT_DoubleEle6_eta1p22_mMax6, &b_HLT_DoubleEle6_eta1p22_mMax6);
   fChain->SetBranchAddress("HLT_DoubleEle6p5_eta1p22_mMax6", &HLT_DoubleEle6p5_eta1p22_mMax6, &b_HLT_DoubleEle6p5_eta1p22_mMax6);
   fChain->SetBranchAddress("HLT_DoubleEle7_eta1p22_mMax6", &HLT_DoubleEle7_eta1p22_mMax6, &b_HLT_DoubleEle7_eta1p22_mMax6);
   fChain->SetBranchAddress("HLT_DoubleEle7p5_eta1p22_mMax6", &HLT_DoubleEle7p5_eta1p22_mMax6, &b_HLT_DoubleEle7p5_eta1p22_mMax6);
   fChain->SetBranchAddress("HLT_DoubleEle8_eta1p22_mMax6", &HLT_DoubleEle8_eta1p22_mMax6, &b_HLT_DoubleEle8_eta1p22_mMax6);
   fChain->SetBranchAddress("HLT_DoubleEle8p5_eta1p22_mMax6", &HLT_DoubleEle8p5_eta1p22_mMax6, &b_HLT_DoubleEle8p5_eta1p22_mMax6);
   fChain->SetBranchAddress("HLT_DoubleEle9_eta1p22_mMax6", &HLT_DoubleEle9_eta1p22_mMax6, &b_HLT_DoubleEle9_eta1p22_mMax6);
   fChain->SetBranchAddress("HLT_DoubleEle9p5_eta1p22_mMax6", &HLT_DoubleEle9p5_eta1p22_mMax6, &b_HLT_DoubleEle9p5_eta1p22_mMax6);
   fChain->SetBranchAddress("HLT_DoubleEle10_eta1p22_mMax6", &HLT_DoubleEle10_eta1p22_mMax6, &b_HLT_DoubleEle10_eta1p22_mMax6);
   fChain->SetBranchAddress("HLT_Diphoton20_14_eta1p5_R9IdL_AND_HE_AND_IsoTCaloIdT", &HLT_Diphoton20_14_eta1p5_R9IdL_AND_HE_AND_IsoTCaloIdT, &b_HLT_Diphoton20_14_eta1p5_R9IdL_AND_HE_AND_IsoTCaloIdT);
   fChain->SetBranchAddress("HLT_Diphoton20_14_eta1p5_R9IdL_AND_HET_AND_IsoTCaloIdT", &HLT_Diphoton20_14_eta1p5_R9IdL_AND_HET_AND_IsoTCaloIdT, &b_HLT_Diphoton20_14_eta1p5_R9IdL_AND_HET_AND_IsoTCaloIdT);
   fChain->SetBranchAddress("HLT_Diphoton22_14_eta1p5_R9IdL_AND_HE_AND_IsoTCaloIdT", &HLT_Diphoton22_14_eta1p5_R9IdL_AND_HE_AND_IsoTCaloIdT, &b_HLT_Diphoton22_14_eta1p5_R9IdL_AND_HE_AND_IsoTCaloIdT);
   fChain->SetBranchAddress("HLT_Diphoton22_14_eta1p5_R9IdL_AND_HET_AND_IsoTCaloIdT", &HLT_Diphoton22_14_eta1p5_R9IdL_AND_HET_AND_IsoTCaloIdT, &b_HLT_Diphoton22_14_eta1p5_R9IdL_AND_HET_AND_IsoTCaloIdT);
   fChain->SetBranchAddress("HLT_ExpressMuons", &HLT_ExpressMuons, &b_HLT_ExpressMuons);
   fChain->SetBranchAddress("HLT_OnlineMonitorGroup", &HLT_OnlineMonitorGroup, &b_HLT_OnlineMonitorGroup);
   fChain->SetBranchAddress("HLT_PPSMaxTracksPerArm1", &HLT_PPSMaxTracksPerArm1, &b_HLT_PPSMaxTracksPerArm1);
   fChain->SetBranchAddress("HLT_PPSMaxTracksPerRP4", &HLT_PPSMaxTracksPerRP4, &b_HLT_PPSMaxTracksPerRP4);
   fChain->SetBranchAddress("HLTriggerFinalPath", &HLTriggerFinalPath, &b_HLTriggerFinalPath);
   fChain->SetBranchAddress("HLT_EphemeralPhysics", &HLT_EphemeralPhysics, &b_HLT_EphemeralPhysics);
   fChain->SetBranchAddress("HLT_EphemeralZeroBias", &HLT_EphemeralZeroBias, &b_HLT_EphemeralZeroBias);
   fChain->SetBranchAddress("L1_DoubleMu3_SQ_ETMHF30_HTT60er", &L1_DoubleMu3_SQ_ETMHF30_HTT60er, &b_L1_DoubleMu3_SQ_ETMHF30_HTT60er);
   fChain->SetBranchAddress("L1_DoubleMu3_SQ_ETMHF30_Jet60er2p5_OR_DoubleJet40er2p5", &L1_DoubleMu3_SQ_ETMHF30_Jet60er2p5_OR_DoubleJet40er2p5, &b_L1_DoubleMu3_SQ_ETMHF30_Jet60er2p5_OR_DoubleJet40er2p5);
   fChain->SetBranchAddress("L1_DoubleMu3_SQ_ETMHF40_HTT60er", &L1_DoubleMu3_SQ_ETMHF40_HTT60er, &b_L1_DoubleMu3_SQ_ETMHF40_HTT60er);
   fChain->SetBranchAddress("L1_DoubleMu3_SQ_ETMHF40_Jet60er2p5_OR_DoubleJet40er2p5", &L1_DoubleMu3_SQ_ETMHF40_Jet60er2p5_OR_DoubleJet40er2p5, &b_L1_DoubleMu3_SQ_ETMHF40_Jet60er2p5_OR_DoubleJet40er2p5);
   fChain->SetBranchAddress("L1_ETMHF80_SingleJet55er2p5_dPhi_Min2p1", &L1_ETMHF80_SingleJet55er2p5_dPhi_Min2p1, &b_L1_ETMHF80_SingleJet55er2p5_dPhi_Min2p1);
   fChain->SetBranchAddress("L1_ETMHF80_SingleJet55er2p5_dPhi_Min2p6", &L1_ETMHF80_SingleJet55er2p5_dPhi_Min2p6, &b_L1_ETMHF80_SingleJet55er2p5_dPhi_Min2p6);
   fChain->SetBranchAddress("L1_Mu3er1p5_Jet100er2p5_ETMHF30", &L1_Mu3er1p5_Jet100er2p5_ETMHF30, &b_L1_Mu3er1p5_Jet100er2p5_ETMHF30);
   fChain->SetBranchAddress("nCscSeg", &nCscSeg, &b_nCscSeg);
   fChain->SetBranchAddress("nCscRechits", &nCscRechits, &b_nCscRechits);
   fChain->SetBranchAddress("nDtSeg", &nDtSeg, &b_nDtSeg);
   fChain->SetBranchAddress("nDtRechits", &nDtRechits, &b_nDtRechits);
   fChain->SetBranchAddress("nRpc", &nRpc, &b_nRpc);
   fChain->SetBranchAddress("cscRechitsPhi", cscRechitsPhi, &b_cscRechitsPhi);
   fChain->SetBranchAddress("cscRechitsEta", cscRechitsEta, &b_cscRechitsEta);
   fChain->SetBranchAddress("cscRechitsX", cscRechitsX, &b_cscRechitsX);
   fChain->SetBranchAddress("cscRechitsY", cscRechitsY, &b_cscRechitsY);
   fChain->SetBranchAddress("cscRechitsZ", cscRechitsZ, &b_cscRechitsZ);
   fChain->SetBranchAddress("cscRechitsTpeak", cscRechitsTpeak, &b_cscRechitsTpeak);
   fChain->SetBranchAddress("cscRechitsTwire", cscRechitsTwire, &b_cscRechitsTwire);
   fChain->SetBranchAddress("dtRechitCorrectX", dtRechitCorrectX, &b_dtRechitCorrectX);
   fChain->SetBranchAddress("dtRechitCorrectY", dtRechitCorrectY, &b_dtRechitCorrectY);
   fChain->SetBranchAddress("dtRechitCorrectZ", dtRechitCorrectZ, &b_dtRechitCorrectZ);
   fChain->SetBranchAddress("dtRechitCorrectEta", dtRechitCorrectEta, &b_dtRechitCorrectEta);
   fChain->SetBranchAddress("dtRechitCorrectPhi", dtRechitCorrectPhi, &b_dtRechitCorrectPhi);
   fChain->SetBranchAddress("cscSegPhi", cscSegPhi, &b_cscSegPhi);
   fChain->SetBranchAddress("cscSegEta", cscSegEta, &b_cscSegEta);
   fChain->SetBranchAddress("dtSegPhi", dtSegPhi, &b_dtSegPhi);
   fChain->SetBranchAddress("dtSegEta", dtSegEta, &b_dtSegEta);
   fChain->SetBranchAddress("rpcPhi", rpcPhi, &b_rpcPhi);
   fChain->SetBranchAddress("rpcEta", rpcEta, &b_rpcEta);
   fChain->SetBranchAddress("rpcX", rpcX, &b_rpcX);
   fChain->SetBranchAddress("rpcY", rpcY, &b_rpcY);
   fChain->SetBranchAddress("rpcZ", rpcZ, &b_rpcZ);
   fChain->SetBranchAddress("rpcT", rpcT, &b_rpcT);
   fChain->SetBranchAddress("rpcTError", rpcTError, &b_rpcTError);
   fChain->SetBranchAddress("cscRechitsChamber", cscRechitsChamber, &b_cscRechitsChamber);
   fChain->SetBranchAddress("cscRechitsStation", cscRechitsStation, &b_cscRechitsStation);
   fChain->SetBranchAddress("cscRechitsDetId", cscRechitsDetId, &b_cscRechitsDetId);
   fChain->SetBranchAddress("dtRechitStation", dtRechitStation, &b_dtRechitStation);
   fChain->SetBranchAddress("dtRechitWheel", dtRechitWheel, &b_dtRechitWheel);
   fChain->SetBranchAddress("dtRechitSuperLayer", dtRechitSuperLayer, &b_dtRechitSuperLayer);
   fChain->SetBranchAddress("cscSegChamber", cscSegChamber, &b_cscSegChamber);
   fChain->SetBranchAddress("cscSegStation", cscSegStation, &b_cscSegStation);
   fChain->SetBranchAddress("cscSegNRecHits", cscSegNRecHits, &b_cscSegNRecHits);
   fChain->SetBranchAddress("dtSegStation", dtSegStation, &b_dtSegStation);
   fChain->SetBranchAddress("dtSegWheel", dtSegWheel, &b_dtSegWheel);
   fChain->SetBranchAddress("rpcBx", rpcBx, &b_rpcBx);
   fChain->SetBranchAddress("rpcRegion", rpcRegion, &b_rpcRegion);
   fChain->SetBranchAddress("rpcRing", rpcRing, &b_rpcRing);
   fChain->SetBranchAddress("rpcSector", rpcSector, &b_rpcSector);
   fChain->SetBranchAddress("rpcStation", rpcStation, &b_rpcStation);
   fChain->SetBranchAddress("rpcLayer", rpcLayer, &b_rpcLayer);
   Notify();
}

Bool_t merged_event::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void merged_event::Show(Long64_t entry)
{
   // Print contents of entry.
   // If entry is not specified, print current entry
   if (!fChain)
      return;
   fChain->Show(entry);
}
Int_t merged_event::Cut(Long64_t entry)
{
   // This function may be called from Loop.
   // returns  1 if entry is accepted.
   // returns -1 otherwise.
   return 1;
}
#endif // #ifdef merged_event_cxx
