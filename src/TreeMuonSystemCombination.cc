#include "RazorHelper.h"
#include "TreeMuonSystemCombination.h"
#include "assert.h"
#include "TTree.h"
#include "DBSCAN.h"

// Constructor
TreeMuonSystemCombination::TreeMuonSystemCombination()
{
  InitVariables();
};
TreeMuonSystemCombination::~TreeMuonSystemCombination()
{
  if (f_) f_->Close();
};
void TreeMuonSystemCombination::InitVariables()
{
  runNum=0; lumiSec=0; evtNum=0; MC_condition = 0;
  npv=0; npu=0;
  pileupWeight = 0;
  pileupWeightDown = 0;
  pileupWeightUp = 0;
  weight=-1.0;rho=-1;
  met=-1; metPhi=-1;
  Flag_HBHENoiseFilter = false; Flag_HBHEIsoNoiseFilter = false; Flag_BadPFMuonFilter = false; Flag_CSCTightHaloFilter = false; Flag_goodVertices = false;
  Flag_ecalBadCalibFilter = false; Flag_all = false; Flag_globalSuperTightHalo2016Filter = false; Flag_BadChargedCandidateFilter = false; Flag_eeBadScFilter = false;

  Flag2_HBHENoiseFilter = false; Flag2_HBHEIsoNoiseFilter = false; Flag2_BadPFMuonFilter = false; Flag2_globalSuperTightHalo2016Filter = false;
  Flag2_globalTightHalo2016Filter = false; Flag2_BadChargedCandidateFilter = false; Flag2_EcalDeadCellTriggerPrimitiveFilter = false; Flag2_ecalBadCalibFilter = false;
  Flag2_eeBadScFilter = false;
  Flag2_all = false;
  mH = 0; mX = 0; ctau = 0;




  gHiggsPt = 0.; gHiggsPhi = 0.; gHiggsEta = 0.; gHiggsE = 0.;
  //CSC

  nCscRechitClusters = 0;
  nDtRechitClusters = 0;
  nDtRings = 0;
  nCscRings = 0;



  for( int i = 0; i < N_MAX_CSC; i++ )
  {

        cscRechitCluster_match_gLLP[i] = false;
        cscRechitCluster_match_gLLP_minDeltaR[i] = 999;
        cscRechitCluster_match_gLLP_index[i] = 999;
        cscRechitCluster_match_gLLP_eta[i] = 999.;
        cscRechitCluster_match_gLLP_phi[i] = 999.;
        cscRechitCluster_match_gLLP_decay_r[i] = 999.;
        cscRechitCluster_match_gLLP_decay_z[i] = 999.;
        cscRechitCluster_match_gLLP_csc[i] = false;
        cscRechitCluster_match_gLLP_dt[i] = false;
        cscRechitCluster_match_gLLP_e[i] = 999.;


        cscRechitClusterSize[i] = -999;
        cscRechitClusterX[i] = -999.;
        cscRechitClusterY[i] = -999.;
        cscRechitClusterZ[i] = -999.;
        cscRechitClusterTimeWeighted[i] = -999.;
        cscRechitCluster_match_dtSeg_0p4[i] = 0;
        cscRechitCluster_match_MB1Seg_0p4[i] = 0;
        cscRechitCluster_match_RE12_0p4[i] = 0;
        cscRechitCluster_match_RB1_0p4[i] = 0;
        cscRechitClusterTimeSpreadWeightedAll[i] = -999.;
        cscRechitClusterTime[i] = -999.;
        cscRechitClusterTimeSpread[i] = -999.;
        cscRechitClusternXY[i] = -999;
        cscRechitClusternZ[i] = -999;
        cscRechitClusterXSpread[i] = -999.;
        cscRechitClusterYSpread[i] = -999.;
        cscRechitClusterZSpread[i] = -999.;
        cscRechitClusterXYSpread[i] = -999.;
        cscRechitClusterRSpread[i] = -999.;
        cscRechitClusterEtaPhiSpread[i] = -999.;
        cscRechitClusterEtaSpread[i] = -999.;
        cscRechitClusterPhiSpread[i] = -999.;
        cscRechitClusterDeltaRSpread[i] = -999.;
        cscRechitClusterMajorAxis[i] = -999.;
        cscRechitClusterMinorAxis[i] = -999.;
        cscRechitClusterSkewX[i] = -999.;
        cscRechitClusterSkewY[i] = -999.;
        cscRechitClusterSkewZ[i] = -999.;
        cscRechitClusterKurtX[i] = -999.;
        cscRechitClusterKurtY[i] = -999.;
        cscRechitClusterKurtZ[i] = -999.;
        cscRechitClusterEta[i] = -999.;
        cscRechitClusterPhi[i] = -999.;
        cscRechitClusterJetVetoPt[i] = 0.0;
        cscRechitClusterJetVetoLooseId[i] = false;
        cscRechitClusterJetVetoTightId[i] = false;
        cscRechitClusterJetVetoE[i] = 0.0;
        cscRechitClusterMuonVetoPt[i] = 0.0;
        cscRechitClusterMuonVetoE[i] = 0.0;
        cscRechitClusterMuonVetoLooseId[i] = false;
        cscRechitClusterMuonVetoGlobal[i] = false;

        cscRechitClusterNChamber[i] = -999;
        cscRechitClusterMaxChamberRatio[i] = -999.;
        cscRechitClusterMaxChamber[i] = -999;
        cscRechitClusterNStation10[i] = -999;
        cscRechitClusterAvgStation10[i] = -999.;
        cscRechitClusterMaxStationRatio[i] = -999.;
        cscRechitClusterMaxStation[i] = -999;
        cscRechitClusterNRechitChamberPlus11[i] = -999;
        cscRechitClusterNRechitChamberPlus12[i] = -999;
        cscRechitClusterNRechitChamberPlus13[i] = -999;
        cscRechitClusterNRechitChamberPlus21[i] = -999;
        cscRechitClusterNRechitChamberPlus22[i] = -999;
        cscRechitClusterNRechitChamberPlus31[i] = -999;
        cscRechitClusterNRechitChamberPlus32[i] = -999;
        cscRechitClusterNRechitChamberPlus41[i] = -999;
        cscRechitClusterNRechitChamberPlus42[i] = -999;
        cscRechitClusterNRechitChamberMinus11[i] = -999;
        cscRechitClusterNRechitChamberMinus12[i] = -999;
        cscRechitClusterNRechitChamberMinus13[i] = -999;
        cscRechitClusterNRechitChamberMinus21[i] = -999;
        cscRechitClusterNRechitChamberMinus22[i] = -999;
        cscRechitClusterNRechitChamberMinus31[i] = -999;
        cscRechitClusterNRechitChamberMinus32[i] = -999;
        cscRechitClusterNRechitChamberMinus41[i] = -999;
        cscRechitClusterNRechitChamberMinus42[i] = -999;
        cscRechitClusterMet_dPhi[i] = 999.;

        dtRechitClusterNSegStation1[i] = 0;
        dtRechitClusterNSegStation2[i] = 0;
        dtRechitClusterNSegStation3[i] = 0;
        dtRechitClusterNSegStation4[i] = 0;

        dtRechitClusterNOppositeSegStation1[i] = 0;
        dtRechitClusterNOppositeSegStation2[i] = 0;
        dtRechitClusterNOppositeSegStation3[i] = 0;
        dtRechitClusterNOppositeSegStation4[i] = 0;

        dtRechitCluster_match_MB1hits_0p4[i] = 0;
        dtRechitCluster_match_MB1hits_0p5[i] = 0;
        dtRechitCluster_match_MB1hits_cosmics_plus[i] = 0;
        dtRechitCluster_match_MB1hits_cosmics_minus[i] = 0;
        dtRechitCluster_match_MB1Seg_0p4[i] = 0;
        dtRechitCluster_match_MB1Seg_0p5[i] = 0;
        dtRechitCluster_match_RPChits_dPhi0p5[i] = 0;
        dtRechitCluster_match_RPCBx_dPhi0p5[i] = 0;
        dtRechitCluster_match_RB1_0p4[i] = 0;
        dtRechitCluster_match_RB1_dPhi0p5[i] = 0;


          dtRechitCluster_match_gLLP[i] = false;
          dtRechitCluster_match_gLLP_minDeltaR[i] = 999;
          dtRechitCluster_match_gLLP_index[i] = 999;
          dtRechitCluster_match_gLLP_eta[i] = 999.;
          dtRechitCluster_match_gLLP_phi[i] = 999.;
          dtRechitCluster_match_gLLP_decay_r[i] = 999.;
          dtRechitCluster_match_gLLP_decay_z[i] = 999.;
          dtRechitCluster_match_gLLP_csc[i] = false;
          dtRechitCluster_match_gLLP_dt[i] = false;



          dtRechitCluster_match_gLLP_e[i] = 999.;


          dtRechitClusterSize[i] = -999;
          dtRechitClusterOverlap[i] = false;
          dtRechitClusterNoiseHit[i] = 0;
          dtRechitClusterNoiseHitStation1[i] = 0;
          dtRechitClusterNoiseHitStation2[i] = 0;
          dtRechitClusterNoiseHitStation3[i] = 0;
          dtRechitClusterNoiseHitStation4[i] = 0;
          dtRechitClusterX[i] = -999.;
          dtRechitClusterY[i] = -999.;
          dtRechitClusterZ[i] = -999.;

          dtRechitClusterWheel[i] = -999;



          dtRechitClusterEta[i] = -999.;
          dtRechitClusterPhi[i] = -999.;
          dtRechitClusterJetVetoPt[i] = 0.0;
          dtRechitClusterJetVetoLooseId[i] = false;
          dtRechitClusterJetVetoTightId[i] = false;


          dtRechitClusterJetVetoE[i] = 0.0;
          dtRechitClusterMuonVetoPt[i] = 0.0;
          dtRechitClusterMuonVetoE[i] = 0.0;


          dtRechitClusterMuonVetoTightId[i] = false;
          dtRechitClusterMuonVetoLooseId[i] = false;
          dtRechitClusterMuonVetoGlobal[i] = false;



          dtRechitClusterNChamber[i] = -999;
          dtRechitClusterMaxChamberRatio[i] = -999.;
          dtRechitClusterMaxChamber[i] = -999;
          dtRechitClusterNStation10[i] = -999;
          dtRechitClusterAvgStation10[i] = -999.;
          dtRechitClusterMaxStationRatio[i] = -999.;
          dtRechitClusterMaxStation[i] = -999;

          dtRechitClusterNHitStation1[i] = -999;
          dtRechitClusterNHitStation2[i] = -999;
          dtRechitClusterNHitStation3[i] = -999;
          dtRechitClusterNHitStation4[i] = -999;
          dtRechitClusterNHitWheel0[i] = -999;
          dtRechitClusterNHitWheel1[i] = -999;
          dtRechitClusterNHitWheel2[i] = -999;
          dtRechitClusterMet_dPhi[i] = 999.;
	  dtRechitClusternXY[i] = -999;
	  dtRechitClusternZ[i] = -999;
	  dtRechitClusterXSpread[i] = -999.;
	  dtRechitClusterYSpread[i] = -999.;
	  dtRechitClusterZSpread[i] = -999.;
	  dtRechitClusterXYSpread[i] = -999.;
	  dtRechitClusterRSpread[i] = -999.;
	  dtRechitClusterEtaPhiSpread[i] = -999.;
	  dtRechitClusterEtaSpread[i] = -999.;
	  dtRechitClusterPhiSpread[i] = -999.;
	  dtRechitClusterDeltaRSpread[i] = -999.;
	  dtRechitClusterMajorAxis[i] = -999.;
	  dtRechitClusterMinorAxis[i] = -999.;
	  dtRechitClusterSkewX[i] = -999.;
	  dtRechitClusterSkewY[i] = -999.;
	  dtRechitClusterSkewZ[i] = -999.;
	  dtRechitClusterKurtX[i] = -999.;
	  dtRechitClusterKurtY[i] = -999.;
	  dtRechitClusterKurtZ[i] = -999.;


  }

   nGLLP = 0;
  for(int i = 0;i<N_MAX_LLP;i++)
  {


    gLLP_eta[i] = 0.0;
    gLLP_phi[i] = 0.0;
    gLLP_beta[i] = 0.0;
    gLLP_e[i] = 0.0;
    gLLP_pt[i] = 0.0;
    gLLP_csc[i] = 0.0;
    gLLP_dt[i] = 0.0;
    gLLP_ctau[i] = 0.0;
    gLLP_decay_vertex_r[i] = 0.0;
    gLLP_decay_vertex_x[i] = 0.0;
    gLLP_decay_vertex_y[i] = 0.0;
    gLLP_decay_vertex_z[i] = 0.0;


  }


  //leptons
  nLeptons = 0;
  for( int i = 0; i < N_MAX_LEPTONS; i++ )
  {
    lepE[i]      = -999.;
    lepPt[i]     = -999.;
    lepEta[i]    = -999.;
    lepPhi[i]    = -999.;
    lepPdgId[i]  = -999;
    lepDZ[i]     = -999.;
    lepTightId[i] = false;
    lepPassLooseIso[i] = false;
    lepPassTightIso[i] = false;
    lepPassVTightIso[i] = false;
    lepPassVTightIso[i] = false;


  }

  //jets
  nJets = 0;
  for( int i = 0; i < N_MAX_JETS; i++ )
  {
    jetE[i]      = -999.;
    jetPt[i]     = -999.;
    jetEta[i]    = -999.;
    jetPhi[i]    = -999.;
    jetTightPassId[i] = false;
  }

  for(int i = 0; i <NTriggersMAX; i++){
    HLTDecision[i] = false;
  }

};

void TreeMuonSystemCombination::InitTree()
{
  assert(tree_);
  InitVariables();

  tree_->SetBranchAddress("runNum",      &runNum);
  tree_->SetBranchAddress("MC_condition",      &MC_condition);
  tree_->SetBranchAddress("lumiSec",     &lumiSec);
  tree_->SetBranchAddress("evtNum",      &evtNum);
  tree_->SetBranchAddress("mX",      &mX);
  tree_->SetBranchAddress("mH",      &mH);
  tree_->SetBranchAddress("ctau",      &ctau);

  tree_->SetBranchAddress("npv",         &npv);
  tree_->SetBranchAddress("npu",         &npu);
  tree_->SetBranchAddress("weight",      &weight);


  tree_->SetBranchAddress("pileupWeight",      &pileupWeight);
  tree_->SetBranchAddress("pileupWeightUp",      &pileupWeightUp);
  tree_->SetBranchAddress("pileupWeightDown",      &pileupWeightDown);
  tree_->SetBranchAddress("Flag_HBHENoiseFilter",      &Flag_HBHENoiseFilter);
  tree_->SetBranchAddress("Flag_HBHEIsoNoiseFilter",      &Flag_HBHEIsoNoiseFilter);
  tree_->SetBranchAddress("Flag_BadPFMuonFilter",      &Flag_BadPFMuonFilter);
  tree_->SetBranchAddress("Flag_CSCTightHaloFilter",      &Flag_CSCTightHaloFilter);
  tree_->SetBranchAddress("Flag_BadChargedCandidateFilter",      &Flag_BadChargedCandidateFilter);
  tree_->SetBranchAddress("Flag_eeBadScFilter",      &Flag_eeBadScFilter);
  tree_->SetBranchAddress("Flag_globalSuperTightHalo2016Filter",      &Flag_globalSuperTightHalo2016Filter);
  tree_->SetBranchAddress("Flag_goodVertices",      &Flag_goodVertices);
  tree_->SetBranchAddress("Flag_ecalBadCalibFilter",      &Flag_ecalBadCalibFilter);
  // tree_->SetBranchAddress("Flag_all",      &Flag_all);

  tree_->SetBranchAddress("Flag2_HBHENoiseFilter",      &Flag2_HBHENoiseFilter);
  tree_->SetBranchAddress("Flag2_HBHEIsoNoiseFilter",      &Flag2_HBHEIsoNoiseFilter);
  tree_->SetBranchAddress("Flag2_BadPFMuonFilter",      &Flag2_BadPFMuonFilter);
  tree_->SetBranchAddress("Flag2_globalSuperTightHalo2016Filter",      &Flag2_globalSuperTightHalo2016Filter);
  tree_->SetBranchAddress("Flag2_globalTightHalo2016Filter",      &Flag2_globalTightHalo2016Filter);
  tree_->SetBranchAddress("Flag2_BadChargedCandidateFilter",      &Flag2_BadChargedCandidateFilter);
  tree_->SetBranchAddress("Flag2_EcalDeadCellTriggerPrimitiveFilter",      &Flag2_EcalDeadCellTriggerPrimitiveFilter);
  tree_->SetBranchAddress("Flag2_ecalBadCalibFilter",      &Flag2_ecalBadCalibFilter);
  tree_->SetBranchAddress("Flag2_eeBadScFilter",      &Flag2_eeBadScFilter);
  tree_->SetBranchAddress("Flag2_all",      &Flag2_all);


  tree_->SetBranchAddress("rho",         &rho);
  tree_->SetBranchAddress("met",         &met);
  tree_->SetBranchAddress("metPhi",      &metPhi);




  tree_->SetBranchAddress("gHiggsPt",      &gHiggsPt);
  tree_->SetBranchAddress("gHiggsPhi",      &gHiggsPhi);
  tree_->SetBranchAddress("gHiggsE",      &gHiggsE);
  tree_->SetBranchAddress("gHiggsEta",      &gHiggsEta);



  tree_->SetBranchAddress("nCscRings",             &nCscRings);
  tree_->SetBranchAddress("nDtRings",             &nDtRings);



  tree_->SetBranchAddress("nDtRechitClusters",             &nDtRechitClusters);


  tree_->SetBranchAddress("dtRechitClusterNSegStation1",             &dtRechitClusterNSegStation1);
  tree_->SetBranchAddress("dtRechitClusterNSegStation2",             &dtRechitClusterNSegStation2);
  tree_->SetBranchAddress("dtRechitClusterNSegStation3",             &dtRechitClusterNSegStation3);
  tree_->SetBranchAddress("dtRechitClusterNSegStation4",             &dtRechitClusterNSegStation4);

  tree_->SetBranchAddress("dtRechitClusterNOppositeSegStation1",             &dtRechitClusterNOppositeSegStation1);
  tree_->SetBranchAddress("dtRechitClusterNOppositeSegStation2",             &dtRechitClusterNOppositeSegStation2);
  tree_->SetBranchAddress("dtRechitClusterNOppositeSegStation3",             &dtRechitClusterNOppositeSegStation3);
  tree_->SetBranchAddress("dtRechitClusterNOppositeSegStation4",             &dtRechitClusterNOppositeSegStation4);

  tree_->SetBranchAddress("dtRechitCluster_match_gLLP",             &dtRechitCluster_match_gLLP);
  tree_->SetBranchAddress("dtRechitCluster_match_gLLP_index",             &dtRechitCluster_match_gLLP_index);
  tree_->SetBranchAddress("dtRechitCluster_match_gLLP_minDeltaR",             &dtRechitCluster_match_gLLP_minDeltaR);
  tree_->SetBranchAddress("dtRechitCluster_match_gLLP_eta",             &dtRechitCluster_match_gLLP_eta);
  tree_->SetBranchAddress("dtRechitCluster_match_gLLP_phi",             &dtRechitCluster_match_gLLP_phi);
  tree_->SetBranchAddress("dtRechitCluster_match_gLLP_decay_r",             &dtRechitCluster_match_gLLP_decay_r);
  tree_->SetBranchAddress("dtRechitCluster_match_gLLP_decay_z",             &dtRechitCluster_match_gLLP_decay_z);
  tree_->SetBranchAddress("dtRechitCluster_match_gLLP_csc",             &dtRechitCluster_match_gLLP_csc);
  tree_->SetBranchAddress("dtRechitCluster_match_gLLP_dt",             &dtRechitCluster_match_gLLP_dt);

  tree_->SetBranchAddress("dtRechitCluster_match_gLLP_e",             &dtRechitCluster_match_gLLP_e);

  tree_->SetBranchAddress("dtRechitClusterMaxStation",             &dtRechitClusterMaxStation);
  tree_->SetBranchAddress("dtRechitClusterMaxStationRatio",             &dtRechitClusterMaxStationRatio);
  tree_->SetBranchAddress("dtRechitClusterNStation10",             &dtRechitClusterNStation10);
  tree_->SetBranchAddress("dtRechitClusterAvgStation10",             &dtRechitClusterAvgStation10);
  tree_->SetBranchAddress("dtRechitClusterMaxChamber",             &dtRechitClusterMaxChamber);
  tree_->SetBranchAddress("dtRechitClusterMaxChamberRatio",             &dtRechitClusterMaxChamberRatio);
  tree_->SetBranchAddress("dtRechitClusterNChamber",             &dtRechitClusterNChamber);
  tree_->SetBranchAddress("dtRechitClusterX",             dtRechitClusterX);
  tree_->SetBranchAddress("dtRechitClusterY",             dtRechitClusterY);
  tree_->SetBranchAddress("dtRechitClusterZ",             dtRechitClusterZ);
  tree_->SetBranchAddress("dtRechitClusterWheel",             dtRechitClusterWheel);


  tree_->SetBranchAddress("dtRechitClusterEta",             dtRechitClusterEta);
  tree_->SetBranchAddress("dtRechitClusterPhi",             dtRechitClusterPhi);

  tree_->SetBranchAddress("dtRechitClusterJetVetoPt",             dtRechitClusterJetVetoPt);
  tree_->SetBranchAddress("dtRechitClusterJetVetoLooseId",             dtRechitClusterJetVetoLooseId);
  tree_->SetBranchAddress("dtRechitClusterJetVetoTightId",             dtRechitClusterJetVetoTightId);
  tree_->SetBranchAddress("dtRechitClusterJetVetoE",             dtRechitClusterJetVetoE);
  tree_->SetBranchAddress("dtRechitClusterMuonVetoPt",             dtRechitClusterMuonVetoPt);
  tree_->SetBranchAddress("dtRechitClusterMuonVetoE",             dtRechitClusterMuonVetoE);

    tree_->SetBranchAddress("dtRechitClusterMuonVetoTightId",             dtRechitClusterMuonVetoTightId);
    tree_->SetBranchAddress("dtRechitClusterMuonVetoLooseId",             dtRechitClusterMuonVetoLooseId);
    tree_->SetBranchAddress("dtRechitClusterMuonVetoGlobal",             dtRechitClusterMuonVetoGlobal);

  tree_->SetBranchAddress("dtRechitClusterSize",             dtRechitClusterSize);
  tree_->SetBranchAddress("dtRechitClusterOverlap",             dtRechitClusterOverlap);
  tree_->SetBranchAddress("dtRechitClusterNoiseHit",             dtRechitClusterNoiseHit);
  tree_->SetBranchAddress("dtRechitClusterNoiseHitStation1",             dtRechitClusterNoiseHitStation1);
  tree_->SetBranchAddress("dtRechitClusterNoiseHitStation2",             dtRechitClusterNoiseHitStation2);
  tree_->SetBranchAddress("dtRechitClusterNoiseHitStation3",             dtRechitClusterNoiseHitStation3);
  tree_->SetBranchAddress("dtRechitClusterNoiseHitStation4",             dtRechitClusterNoiseHitStation4);

  tree_->SetBranchAddress("dtRechitClusterNHitStation1",             dtRechitClusterNHitStation1);
  tree_->SetBranchAddress("dtRechitClusterNHitStation2",             dtRechitClusterNHitStation2);
  tree_->SetBranchAddress("dtRechitClusterNHitStation3",             dtRechitClusterNHitStation3);
  tree_->SetBranchAddress("dtRechitClusterNHitStation4",             dtRechitClusterNHitStation4);

  tree_->SetBranchAddress("dtRechitClusterNHitWheel0",             dtRechitClusterNHitWheel0);
  tree_->SetBranchAddress("dtRechitClusterNHitWheel1",             dtRechitClusterNHitWheel1);
  tree_->SetBranchAddress("dtRechitClusterNHitWheel2",             dtRechitClusterNHitWheel2);

  tree_->SetBranchAddress("dtRechitCluster_match_MB1hits_0p4",             dtRechitCluster_match_MB1hits_0p4);
  tree_->SetBranchAddress("dtRechitCluster_match_MB1hits_0p5",             dtRechitCluster_match_MB1hits_0p5);
  tree_->SetBranchAddress("dtRechitCluster_match_MB1hits_cosmics_plus",             dtRechitCluster_match_MB1hits_cosmics_plus);
  tree_->SetBranchAddress("dtRechitCluster_match_MB1hits_cosmics_minus",             dtRechitCluster_match_MB1hits_cosmics_minus);
  tree_->SetBranchAddress("dtRechitCluster_match_MB1Seg_0p4",             dtRechitCluster_match_MB1Seg_0p4);
  tree_->SetBranchAddress("dtRechitCluster_match_MB1Seg_0p5",             dtRechitCluster_match_MB1Seg_0p5);
  tree_->SetBranchAddress("dtRechitCluster_match_RPChits_dPhi0p5",             dtRechitCluster_match_RPChits_dPhi0p5);
  tree_->SetBranchAddress("dtRechitCluster_match_RPCBx_dPhi0p5",             dtRechitCluster_match_RPCBx_dPhi0p5);
  tree_->SetBranchAddress("dtRechitCluster_match_RB1_0p4",             dtRechitCluster_match_RB1_0p4);
  tree_->SetBranchAddress("dtRechitCluster_match_RB1_dPhi0p5",             dtRechitCluster_match_RB1_dPhi0p5);



    tree_->SetBranchAddress("dtRechitClusterMet_dPhi",             dtRechitClusterMet_dPhi);
    tree_->SetBranchAddress("dtRechitClusternXY",             dtRechitClusternXY);
    tree_->SetBranchAddress("dtRechitClusternZ",             dtRechitClusternZ);
    tree_->SetBranchAddress("dtRechitClusterXSpread",             dtRechitClusterXSpread);
    tree_->SetBranchAddress("dtRechitClusterYSpread",             dtRechitClusterYSpread);
    tree_->SetBranchAddress("dtRechitClusterZSpread",             dtRechitClusterZSpread);
    tree_->SetBranchAddress("dtRechitClusterXYSpread",             dtRechitClusterXYSpread);
    tree_->SetBranchAddress("dtRechitClusterRSpread",             dtRechitClusterRSpread);
    tree_->SetBranchAddress("dtRechitClusterEtaPhiSpread",             dtRechitClusterEtaPhiSpread);
    tree_->SetBranchAddress("dtRechitClusterEtaSpread",             dtRechitClusterEtaSpread);
    tree_->SetBranchAddress("dtRechitClusterPhiSpread",             dtRechitClusterPhiSpread);
    tree_->SetBranchAddress("dtRechitClusterDeltaRSpread",             dtRechitClusterDeltaRSpread);
    tree_->SetBranchAddress("dtRechitClusterMajorAxis",             dtRechitClusterMajorAxis);
    tree_->SetBranchAddress("dtRechitClusterMinorAxis",             dtRechitClusterMinorAxis);
    tree_->SetBranchAddress("dtRechitClusterSkewX",             dtRechitClusterSkewX);
    tree_->SetBranchAddress("dtRechitClusterSkewY",             dtRechitClusterSkewY);
    tree_->SetBranchAddress("dtRechitClusterSkewZ",             dtRechitClusterSkewZ);
    tree_->SetBranchAddress("dtRechitClusterKurtX",             dtRechitClusterKurtX);
    tree_->SetBranchAddress("dtRechitClusterKurtY",             dtRechitClusterKurtY);
    tree_->SetBranchAddress("dtRechitClusterKurtZ",             dtRechitClusterKurtZ);






  tree_->SetBranchAddress("nCscRechitClusters",             &nCscRechitClusters);
  tree_->SetBranchAddress("cscRechitCluster_match_gLLP",             &cscRechitCluster_match_gLLP);
  tree_->SetBranchAddress("cscRechitCluster_match_gLLP_index",             &cscRechitCluster_match_gLLP_index);
  tree_->SetBranchAddress("cscRechitCluster_match_gLLP_minDeltaR",             &cscRechitCluster_match_gLLP_minDeltaR);
  tree_->SetBranchAddress("cscRechitCluster_match_gLLP_eta",             &cscRechitCluster_match_gLLP_eta);
  tree_->SetBranchAddress("cscRechitCluster_match_gLLP_phi",             &cscRechitCluster_match_gLLP_phi);
  tree_->SetBranchAddress("cscRechitCluster_match_gLLP_decay_r",             &cscRechitCluster_match_gLLP_decay_r);
  tree_->SetBranchAddress("cscRechitCluster_match_gLLP_decay_z",             &cscRechitCluster_match_gLLP_decay_z);
  tree_->SetBranchAddress("cscRechitCluster_match_gLLP_csc",             &cscRechitCluster_match_gLLP_csc);
  tree_->SetBranchAddress("cscRechitCluster_match_gLLP_dt",             &cscRechitCluster_match_gLLP_dt);

  tree_->SetBranchAddress("cscRechitCluster_match_gLLP_e",             &cscRechitCluster_match_gLLP_e);

  tree_->SetBranchAddress("cscRechitClusterMaxStation",             &cscRechitClusterMaxStation);
  tree_->SetBranchAddress("cscRechitClusterMaxStationRatio",             &cscRechitClusterMaxStationRatio);
  tree_->SetBranchAddress("cscRechitClusterNStation10",             &cscRechitClusterNStation10);
  tree_->SetBranchAddress("cscRechitClusterAvgStation10",             &cscRechitClusterAvgStation10);
  tree_->SetBranchAddress("cscRechitClusterMaxChamber",             &cscRechitClusterMaxChamber);
  tree_->SetBranchAddress("cscRechitClusterMaxChamberRatio",             &cscRechitClusterMaxChamberRatio);
  tree_->SetBranchAddress("cscRechitClusterNChamber",             &cscRechitClusterNChamber);
  tree_->SetBranchAddress("cscRechitClusterX",             cscRechitClusterX);
  tree_->SetBranchAddress("cscRechitClusterY",             cscRechitClusterY);
  tree_->SetBranchAddress("cscRechitClusterZ",             cscRechitClusterZ);
  tree_->SetBranchAddress("cscRechitClusterTimeWeighted",             cscRechitClusterTimeWeighted);


  tree_->SetBranchAddress("cscRechitClusterTimeSpreadWeightedAll",             cscRechitClusterTimeSpreadWeightedAll);


  tree_->SetBranchAddress("cscRechitClusterTime",             cscRechitClusterTime);
  tree_->SetBranchAddress("cscRechitClusterTimeSpread",             cscRechitClusterTimeSpread);

  tree_->SetBranchAddress("cscRechitClusternXY",             cscRechitClusternXY);
  tree_->SetBranchAddress("cscRechitClusternZ",             cscRechitClusternZ);
  tree_->SetBranchAddress("cscRechitClusterXSpread",             cscRechitClusterXSpread);
  tree_->SetBranchAddress("cscRechitClusterYSpread",             cscRechitClusterYSpread);
  tree_->SetBranchAddress("cscRechitClusterZSpread",             cscRechitClusterZSpread);
  tree_->SetBranchAddress("cscRechitClusterXYSpread",             cscRechitClusterXYSpread);
  tree_->SetBranchAddress("cscRechitClusterRSpread",             cscRechitClusterRSpread);
  tree_->SetBranchAddress("cscRechitClusterEtaPhiSpread",             cscRechitClusterEtaPhiSpread);
  tree_->SetBranchAddress("cscRechitClusterEtaSpread",             cscRechitClusterEtaSpread);
  tree_->SetBranchAddress("cscRechitClusterPhiSpread",             cscRechitClusterPhiSpread);
  tree_->SetBranchAddress("cscRechitClusterDeltaRSpread",             cscRechitClusterDeltaRSpread);
  tree_->SetBranchAddress("cscRechitClusterMajorAxis",             cscRechitClusterMajorAxis);
  tree_->SetBranchAddress("cscRechitClusterMinorAxis",             cscRechitClusterMinorAxis);
  tree_->SetBranchAddress("cscRechitClusterSkewX",             cscRechitClusterSkewX);
  tree_->SetBranchAddress("cscRechitClusterSkewY",             cscRechitClusterSkewY);
  tree_->SetBranchAddress("cscRechitClusterSkewZ",             cscRechitClusterSkewZ);
  tree_->SetBranchAddress("cscRechitClusterKurtX",             cscRechitClusterKurtX);
  tree_->SetBranchAddress("cscRechitClusterKurtY",             cscRechitClusterKurtY);
  tree_->SetBranchAddress("cscRechitClusterKurtZ",             cscRechitClusterKurtZ);


  tree_->SetBranchAddress("cscRechitClusterEta",             cscRechitClusterEta);
  tree_->SetBranchAddress("cscRechitClusterPhi",             cscRechitClusterPhi);


  tree_->SetBranchAddress("cscRechitClusterJetVetoPt",             cscRechitClusterJetVetoPt);


  tree_->SetBranchAddress("cscRechitClusterJetVetoLooseId",             cscRechitClusterJetVetoLooseId);
  tree_->SetBranchAddress("cscRechitClusterJetVetoTightId",             cscRechitClusterJetVetoTightId);

  tree_->SetBranchAddress("cscRechitClusterJetVetoE",             cscRechitClusterJetVetoE);
  tree_->SetBranchAddress("cscRechitClusterMuonVetoPt",             cscRechitClusterMuonVetoPt);
  tree_->SetBranchAddress("cscRechitClusterMuonVetoE",             cscRechitClusterMuonVetoE);


    tree_->SetBranchAddress("cscRechitClusterMuonVetoLooseId",             cscRechitClusterMuonVetoLooseId);
    tree_->SetBranchAddress("cscRechitClusterMuonVetoGlobal",             cscRechitClusterMuonVetoGlobal);



  tree_->SetBranchAddress("cscRechitClusterSize",             cscRechitClusterSize);
  tree_->SetBranchAddress("cscRechitCluster_match_dtSeg_0p4",             cscRechitCluster_match_dtSeg_0p4);
  tree_->SetBranchAddress("cscRechitCluster_match_MB1Seg_0p4",             cscRechitCluster_match_MB1Seg_0p4);
  tree_->SetBranchAddress("cscRechitCluster_match_RE12_0p4",             cscRechitCluster_match_RE12_0p4);
  tree_->SetBranchAddress("cscRechitCluster_match_RB1_0p4",             cscRechitCluster_match_RB1_0p4);



  tree_->SetBranchAddress("cscRechitClusterNRechitChamberPlus11",             cscRechitClusterNRechitChamberPlus11);
  tree_->SetBranchAddress("cscRechitClusterNRechitChamberPlus12",             cscRechitClusterNRechitChamberPlus12);
  tree_->SetBranchAddress("cscRechitClusterNRechitChamberPlus13",             cscRechitClusterNRechitChamberPlus13);
  tree_->SetBranchAddress("cscRechitClusterNRechitChamberPlus21",             cscRechitClusterNRechitChamberPlus21);
  tree_->SetBranchAddress("cscRechitClusterNRechitChamberPlus22",             cscRechitClusterNRechitChamberPlus22);
  tree_->SetBranchAddress("cscRechitClusterNRechitChamberPlus31",             cscRechitClusterNRechitChamberPlus31);
  tree_->SetBranchAddress("cscRechitClusterNRechitChamberPlus32",             cscRechitClusterNRechitChamberPlus32);
  tree_->SetBranchAddress("cscRechitClusterNRechitChamberPlus41",             cscRechitClusterNRechitChamberPlus41);
  tree_->SetBranchAddress("cscRechitClusterNRechitChamberPlus42",             cscRechitClusterNRechitChamberPlus42);

  tree_->SetBranchAddress("cscRechitClusterNRechitChamberMinus11",             cscRechitClusterNRechitChamberMinus11);
  tree_->SetBranchAddress("cscRechitClusterNRechitChamberMinus12",             cscRechitClusterNRechitChamberMinus12);
  tree_->SetBranchAddress("cscRechitClusterNRechitChamberMinus13",             cscRechitClusterNRechitChamberMinus13);
  tree_->SetBranchAddress("cscRechitClusterNRechitChamberMinus21",             cscRechitClusterNRechitChamberMinus21);
  tree_->SetBranchAddress("cscRechitClusterNRechitChamberMinus22",             cscRechitClusterNRechitChamberMinus22);
  tree_->SetBranchAddress("cscRechitClusterNRechitChamberMinus31",             cscRechitClusterNRechitChamberMinus31);
  tree_->SetBranchAddress("cscRechitClusterNRechitChamberMinus32",             cscRechitClusterNRechitChamberMinus32);
  tree_->SetBranchAddress("cscRechitClusterNRechitChamberMinus41",             cscRechitClusterNRechitChamberMinus41);
  tree_->SetBranchAddress("cscRechitClusterNRechitChamberMinus42",             cscRechitClusterNRechitChamberMinus42);

  tree_->SetBranchAddress("cscRechitClusterMet_dPhi",             cscRechitClusterMet_dPhi);


  tree_->SetBranchAddress("nGLLP",    &nGLLP);

  tree_->SetBranchAddress("gLLP_eta",    gLLP_eta);
  tree_->SetBranchAddress("gLLP_phi",    gLLP_phi);
  tree_->SetBranchAddress("gLLP_csc",    gLLP_csc);
  tree_->SetBranchAddress("gLLP_dt",    gLLP_dt);
  tree_->SetBranchAddress("gLLP_ctau",    gLLP_ctau);
  tree_->SetBranchAddress("gLLP_beta",    gLLP_beta);


  tree_->SetBranchAddress("gLLP_e",    gLLP_e);
  tree_->SetBranchAddress("gLLP_pt",    gLLP_pt);

  tree_->SetBranchAddress("gLLP_decay_vertex_r",    gLLP_decay_vertex_r);
  tree_->SetBranchAddress("gLLP_decay_vertex_x",    gLLP_decay_vertex_x);
  tree_->SetBranchAddress("gLLP_decay_vertex_y",    gLLP_decay_vertex_y);
  tree_->SetBranchAddress("gLLP_decay_vertex_z",    gLLP_decay_vertex_z);



  //Leptons


  tree_->SetBranchAddress("nLeptons",    &nLeptons);
  tree_->SetBranchAddress("lepE",        lepE);
  tree_->SetBranchAddress("lepPt",       lepPt);
  tree_->SetBranchAddress("lepEta",      lepEta);
  tree_->SetBranchAddress("lepPhi",      lepPhi);
  tree_->SetBranchAddress("lepPdgId",  lepPdgId);
  tree_->SetBranchAddress("lepDZ",     lepDZ);
  tree_->SetBranchAddress("lepTightId", lepTightId);
  tree_->SetBranchAddress("lepPassLooseIso", lepPassLooseIso);
  tree_->SetBranchAddress("lepPassTightIso", lepPassTightIso);
  tree_->SetBranchAddress("lepPassVTightIso", lepPassVTightIso);
  tree_->SetBranchAddress("lepPassVTightIso", lepPassVTightIso);


  //Z-candidate

  //jets
  tree_->SetBranchAddress("nJets",     &nJets);
  tree_->SetBranchAddress("jetE",      jetE);
  tree_->SetBranchAddress("jetPt",     jetPt);
  tree_->SetBranchAddress("jetEta",    jetEta);
  tree_->SetBranchAddress("jetPhi",    jetPhi);
  tree_->SetBranchAddress("jetTightPassId", jetTightPassId);
  // triggers
  tree_->SetBranchAddress("HLTDecision",   HLTDecision);


};

void TreeMuonSystemCombination::LoadTree(const char* file)
{
  f_ = TFile::Open(file);
  assert(f_);
  tree_ = dynamic_cast<TTree*>(f_->Get("MuonSystem"));
  InitTree();
  assert(tree_);
};

void TreeMuonSystemCombination::CreateTree()
{
  tree_ = new TTree("MuonSystem","MuonSystem");
  f_ = 0;

  tree_->Branch("runNum",      &runNum,     "runNum/i");      // event run number
  tree_->Branch("MC_condition",      &MC_condition,     "MC_condition/i");      // event run number
  tree_->Branch("lumiSec",     &lumiSec,    "lumiSec/i");     // event lumi section
  tree_->Branch("evtNum",      &evtNum,     "evtNum/i");      // event number
  tree_->Branch("mH",      &mH,     "mH/I");      // event number
  tree_->Branch("mX",      &mX,     "mX/I");      // event number
  tree_->Branch("ctau",      &ctau,     "ctau/I");      // event number

  tree_->Branch("npv",         &npv,        "npv/i");         // number of primary vertices
  tree_->Branch("npu",         &npu,        "npu/i");         // number of in-time PU events (MC)
  tree_->Branch("weight",      &weight,     "weight/F");

  tree_->Branch("pileupWeight",      &pileupWeight,     "pileupWeight/F");
  tree_->Branch("pileupWeightUp",      &pileupWeightUp,     "pileupWeightUp/F");
  tree_->Branch("pileupWeightDown",      &pileupWeightDown,     "pileupWeightDown/F");
  tree_->Branch("Flag_HBHENoiseFilter",      &Flag_HBHENoiseFilter,     "Flag_HBHENoiseFilter/O");
  tree_->Branch("Flag_BadPFMuonFilter",      &Flag_BadPFMuonFilter,     "Flag_BadPFMuonFilter/O");
  tree_->Branch("Flag_HBHEIsoNoiseFilter",      &Flag_HBHEIsoNoiseFilter,     "Flag_HBHEIsoNoiseFilter/O");
  tree_->Branch("Flag_CSCTightHaloFilter",      &Flag_CSCTightHaloFilter,     "Flag_CSCTightHaloFilter/O");
  tree_->Branch("Flag_globalSuperTightHalo2016Filter",      &Flag_globalSuperTightHalo2016Filter,     "Flag_globalSuperTightHalo2016Filter/O");
  tree_->Branch("Flag_goodVertices",      &Flag_goodVertices,     "Flag_goodVertices/O");
  tree_->Branch("Flag_ecalBadCalibFilter",      &Flag_ecalBadCalibFilter,     "Flag_ecalBadCalibFilter/O");
  tree_->Branch("Flag_BadChargedCandidateFilter",      &Flag_BadChargedCandidateFilter,     "Flag_BadChargedCandidateFilter/O");
  tree_->Branch("Flag_eeBadScFilter",      &Flag_eeBadScFilter,     "Flag_eeBadScFilter/O");
  tree_->Branch("Flag_all",      &Flag_all,     "Flag_all/O");

  tree_->Branch("Flag2_HBHENoiseFilter",      &Flag2_HBHENoiseFilter,     "Flag2_HBHENoiseFilter/O");
  tree_->Branch("Flag2_HBHEIsoNoiseFilter",      &Flag2_HBHEIsoNoiseFilter,     "Flag2_HBHEIsoNoiseFilter/O");
  tree_->Branch("Flag2_BadPFMuonFilter",      &Flag2_BadPFMuonFilter,     "Flag2_BadPFMuonFilter/O");
  tree_->Branch("Flag2_globalSuperTightHalo2016Filter",      &Flag2_globalSuperTightHalo2016Filter,     "Flag2_globalSuperTightHalo2016Filter/O");
  tree_->Branch("Flag2_globalTightHalo2016Filter",      &Flag2_globalTightHalo2016Filter,     "Flag2_globalTightHalo2016Filter/O");
  tree_->Branch("Flag2_BadChargedCandidateFilter",      &Flag2_BadChargedCandidateFilter,     "Flag2_BadChargedCandidateFilter/O");
  tree_->Branch("Flag2_EcalDeadCellTriggerPrimitiveFilter",      &Flag2_EcalDeadCellTriggerPrimitiveFilter,     "Flag2_EcalDeadCellTriggerPrimitiveFilter/O");
  tree_->Branch("Flag2_ecalBadCalibFilter",      &Flag2_ecalBadCalibFilter,     "Flag2_ecalBadCalibFilter/O");
  tree_->Branch("Flag2_eeBadScFilter",      &Flag2_eeBadScFilter,     "Flag2_eeBadScFilter/O");
  tree_->Branch("Flag2_all",      &Flag2_all,     "Flag2_all/O");



  tree_->Branch("rho",         &rho,        "rho/F");
  tree_->Branch("met",         &met,        "met/F");         // MET
  tree_->Branch("metPhi",      &metPhi,     "metPhi/F");      // phi(MET)



  tree_->Branch("gHiggsPt",      &gHiggsPt,     "gHiggsPt/F");      // phi(MET)
  tree_->Branch("gHiggsE",      &gHiggsE,     "gHiggsE/F");      // phi(MET)
  tree_->Branch("gHiggsEta",      &gHiggsEta,     "gHiggsEta/F");      // phi(MET)
  tree_->Branch("gHiggsPhi",      &gHiggsPhi,     "gHiggsPhi/F");      // phi(MET)




  tree_->Branch("nCscRings",             &nCscRings, "nCscRings/I");
  tree_->Branch("nDtRings",             &nDtRings, "nDtRings/I");

    tree_->Branch("nCscRechitClusters",             &nCscRechitClusters, "nCscRechitClusters/I");

    tree_->Branch("cscRechitCluster_match_gLLP",             cscRechitCluster_match_gLLP,             "cscRechitCluster_match_gLLP[nCscRechitClusters]/O");
    tree_->Branch("cscRechitCluster_match_gLLP_minDeltaR",             cscRechitCluster_match_gLLP_minDeltaR,             "cscRechitCluster_match_gLLP_minDeltaR[nCscRechitClusters]/F");
    tree_->Branch("cscRechitCluster_match_gLLP_index",             cscRechitCluster_match_gLLP_index,             "cscRechitCluster_match_gLLP_index[nCscRechitClusters]/I");
    tree_->Branch("cscRechitCluster_match_gLLP_eta",             cscRechitCluster_match_gLLP_eta, "cscRechitCluster_match_gLLP_eta[nCscRechitClusters]/F");
    tree_->Branch("cscRechitCluster_match_gLLP_phi",             cscRechitCluster_match_gLLP_phi, "cscRechitCluster_match_gLLP_phi[nCscRechitClusters]/F");
    tree_->Branch("cscRechitCluster_match_gLLP_decay_r",             cscRechitCluster_match_gLLP_decay_r, "cscRechitCluster_match_gLLP_decay_r[nCscRechitClusters]/F");
    tree_->Branch("cscRechitCluster_match_gLLP_decay_z",             cscRechitCluster_match_gLLP_decay_z, "cscRechitCluster_match_gLLP_decay_z[nCscRechitClusters]/F");
    tree_->Branch("cscRechitCluster_match_gLLP_csc",             cscRechitCluster_match_gLLP_csc, "cscRechitCluster_match_gLLP_csc[nCscRechitClusters]/O");
    tree_->Branch("cscRechitCluster_match_gLLP_dt",             cscRechitCluster_match_gLLP_dt, "cscRechitCluster_match_gLLP_dt[nCscRechitClusters]/O");
    tree_->Branch("cscRechitCluster_match_gLLP_e",             cscRechitCluster_match_gLLP_e, "cscRechitCluster_match_gLLP_e[nCscRechitClusters]/F");


    tree_->Branch("cscRechitClusterX",             cscRechitClusterX,             "cscRechitClusterX[nCscRechitClusters]/F");
    tree_->Branch("cscRechitClusterY",             cscRechitClusterY,             "cscRechitClusterY[nCscRechitClusters]/F");
    tree_->Branch("cscRechitClusterZ",             cscRechitClusterZ,             "cscRechitClusterZ[nCscRechitClusters]/F");
    tree_->Branch("cscRechitClusterTimeWeighted",             cscRechitClusterTimeWeighted,             "cscRechitClusterTimeWeighted[nCscRechitClusters]/F");

    tree_->Branch("cscRechitClusterTimeSpreadWeightedAll",             cscRechitClusterTimeSpreadWeightedAll,             "cscRechitClusterTimeSpreadWeightedAll[nCscRechitClusters]/F");
    tree_->Branch("cscRechitClusterTime",             cscRechitClusterTime,             "cscRechitClusterTime[nCscRechitClusters]/F");
    tree_->Branch("cscRechitClusterTimeSpread",             cscRechitClusterTimeSpread,             "cscRechitClusterTimeSpread[nCscRechitClusters]/F");
    tree_->Branch("cscRechitClusternXY",             cscRechitClusternXY,             "cscRechitClusternXY[nCscRechitClusters]/I");
    tree_->Branch("cscRechitClusternZ",             cscRechitClusternZ,             "cscRechitClusternZ[nCscRechitClusters]/I");
    tree_->Branch("cscRechitClusterXSpread",             cscRechitClusterXSpread,             "cscRechitClusterXSpread[nCscRechitClusters]/F");
    tree_->Branch("cscRechitClusterYSpread",             cscRechitClusterYSpread,             "cscRechitClusterYSpread[nCscRechitClusters]/F");
    tree_->Branch("cscRechitClusterZSpread",             cscRechitClusterZSpread,             "cscRechitClusterZSpread[nCscRechitClusters]/F");
    tree_->Branch("cscRechitClusterXYSpread",             cscRechitClusterXYSpread,             "cscRechitClusterXYSpread[nCscRechitClusters]/F");
    tree_->Branch("cscRechitClusterRSpread",             cscRechitClusterRSpread,             "cscRechitClusterRSpread[nCscRechitClusters]/F");
    tree_->Branch("cscRechitClusterEtaPhiSpread",             cscRechitClusterEtaPhiSpread,             "cscRechitClusterEtaPhiSpread[nCscRechitClusters]/F");
    tree_->Branch("cscRechitClusterEtaSpread",             cscRechitClusterEtaSpread,             "cscRechitClusterEtaSpread[nCscRechitClusters]/F");
    tree_->Branch("cscRechitClusterPhiSpread",             cscRechitClusterPhiSpread,             "cscRechitClusterPhiSpread[nCscRechitClusters]/F");
    tree_->Branch("cscRechitClusterDeltaRSpread",             cscRechitClusterDeltaRSpread,             "cscRechitClusterDeltaRSpread[nCscRechitClusters]/F");
    tree_->Branch("cscRechitClusterMajorAxis",             cscRechitClusterMajorAxis,             "cscRechitClusterMajorAxis[nCscRechitClusters]/F");
    tree_->Branch("cscRechitClusterMinorAxis",             cscRechitClusterMinorAxis,             "cscRechitClusterMinorAxis[nCscRechitClusters]/F");
    tree_->Branch("cscRechitClusterSkewX",             cscRechitClusterSkewX,             "cscRechitClusterSkewX[nCscRechitClusters]/F");
    tree_->Branch("cscRechitClusterSkewY",             cscRechitClusterSkewY,             "cscRechitClusterSkewY[nCscRechitClusters]/F");
    tree_->Branch("cscRechitClusterSkewZ",             cscRechitClusterSkewZ,             "cscRechitClusterSkewZ[nCscRechitClusters]/F");
    tree_->Branch("cscRechitClusterKurtX",             cscRechitClusterKurtX,             "cscRechitClusterKurtX[nCscRechitClusters]/F");
    tree_->Branch("cscRechitClusterKurtY",             cscRechitClusterKurtY,             "cscRechitClusterKurtY[nCscRechitClusters]/F");
    tree_->Branch("cscRechitClusterKurtZ",             cscRechitClusterKurtZ,             "cscRechitClusterKurtZ[nCscRechitClusters]/F");


    tree_->Branch("cscRechitClusterPhi",             cscRechitClusterPhi,             "cscRechitClusterPhi[nCscRechitClusters]/F");
    tree_->Branch("cscRechitClusterEta",             cscRechitClusterEta,             "cscRechitClusterEta[nCscRechitClusters]/F");
    tree_->Branch("cscRechitClusterJetVetoPt",             cscRechitClusterJetVetoPt,             "cscRechitClusterJetVetoPt[nCscRechitClusters]/F");
    tree_->Branch("cscRechitClusterJetVetoLooseId",             cscRechitClusterJetVetoLooseId,             "cscRechitClusterJetVetoLooseId[nCscRechitClusters]/O");
    tree_->Branch("cscRechitClusterJetVetoTightId",             cscRechitClusterJetVetoTightId,             "cscRechitClusterJetVetoTightId[nCscRechitClusters]/O");
    tree_->Branch("cscRechitClusterJetVetoE",             cscRechitClusterJetVetoE,             "cscRechitClusterJetVetoE[nCscRechitClusters]/F");
    tree_->Branch("cscRechitClusterMuonVetoPt",             cscRechitClusterMuonVetoPt,             "cscRechitClusterMuonVetoPt[nCscRechitClusters]/F");
    tree_->Branch("cscRechitClusterMuonVetoE",             cscRechitClusterMuonVetoE,             "cscRechitClusterMuonVetoE[nCscRechitClusters]/F");

    tree_->Branch("cscRechitClusterMuonVetoLooseId",             cscRechitClusterMuonVetoLooseId,             "cscRechitClusterMuonVetoLooseId[nCscRechitClusters]/O");
    tree_->Branch("cscRechitClusterMuonVetoGlobal",             cscRechitClusterMuonVetoGlobal,             "cscRechitClusterMuonVetoGlobal[nCscRechitClusters]/O");



      tree_->Branch("cscRechitCluster_match_dtSeg_0p4",             cscRechitCluster_match_dtSeg_0p4,             "cscRechitCluster_match_dtSeg_0p4[nCscRechitClusters]/I");
      tree_->Branch("cscRechitCluster_match_MB1Seg_0p4",             cscRechitCluster_match_MB1Seg_0p4,             "cscRechitCluster_match_MB1Seg_0p4[nCscRechitClusters]/I");
      tree_->Branch("cscRechitCluster_match_RE12_0p4",             cscRechitCluster_match_RE12_0p4,             "cscRechitCluster_match_RE12_0p4[nCscRechitClusters]/I");
      tree_->Branch("cscRechitCluster_match_RB1_0p4",             cscRechitCluster_match_RB1_0p4,             "cscRechitCluster_match_RB1_0p4[nCscRechitClusters]/I");

    tree_->Branch("cscRechitClusterSize",             cscRechitClusterSize,             "cscRechitClusterSize[nCscRechitClusters]/I");
    tree_->Branch("cscRechitClusterNStation10",             cscRechitClusterNStation10,             "cscRechitClusterNStation10[nCscRechitClusters]/I");
    tree_->Branch("cscRechitClusterAvgStation10",             cscRechitClusterAvgStation10,             "cscRechitClusterAvgStation10[nCscRechitClusters]/F");


    tree_->Branch("cscRechitClusterMaxStation",             cscRechitClusterMaxStation,             "cscRechitClusterMaxStation[nCscRechitClusters]/I");
    tree_->Branch("cscRechitClusterMaxStationRatio",             cscRechitClusterMaxStationRatio,             "cscRechitClusterMaxStationRatio[nCscRechitClusters]/F");
    tree_->Branch("cscRechitClusterNChamber",             cscRechitClusterNChamber,             "cscRechitClusterNChamber[nCscRechitClusters]/I");
    tree_->Branch("cscRechitClusterMaxChamber",             cscRechitClusterMaxChamber,             "cscRechitClusterMaxChamber[nCscRechitClusters]/I");
    tree_->Branch("cscRechitClusterMaxChamberRatio",             cscRechitClusterMaxChamberRatio,             "cscRechitClusterMaxChamberRatio[nCscRechitClusters]/F");
    tree_->Branch("cscRechitClusterNRechitChamberPlus11",             cscRechitClusterNRechitChamberPlus11,             "cscRechitClusterNRechitChamberPlus11[nCscRechitClusters]/I");
    tree_->Branch("cscRechitClusterNRechitChamberPlus12",             cscRechitClusterNRechitChamberPlus12,             "cscRechitClusterNRechitChamberPlus12[nCscRechitClusters]/I");
    tree_->Branch("cscRechitClusterNRechitChamberPlus13",             cscRechitClusterNRechitChamberPlus13,             "cscRechitClusterNRechitChamberPlus13[nCscRechitClusters]/I");
    tree_->Branch("cscRechitClusterNRechitChamberPlus21",             cscRechitClusterNRechitChamberPlus21,             "cscRechitClusterNRechitChamberPlus21[nCscRechitClusters]/I");
    tree_->Branch("cscRechitClusterNRechitChamberPlus22",             cscRechitClusterNRechitChamberPlus22,             "cscRechitClusterNRechitChamberPlus22[nCscRechitClusters]/I");
    tree_->Branch("cscRechitClusterNRechitChamberPlus31",             cscRechitClusterNRechitChamberPlus31,             "cscRechitClusterNRechitChamberPlus31[nCscRechitClusters]/I");
    tree_->Branch("cscRechitClusterNRechitChamberPlus32",             cscRechitClusterNRechitChamberPlus32,             "cscRechitClusterNRechitChamberPlus32[nCscRechitClusters]/I");
    tree_->Branch("cscRechitClusterNRechitChamberPlus41",             cscRechitClusterNRechitChamberPlus41,             "cscRechitClusterNRechitChamberPlus41[nCscRechitClusters]/I");
    tree_->Branch("cscRechitClusterNRechitChamberPlus42",             cscRechitClusterNRechitChamberPlus42,             "cscRechitClusterNRechitChamberPlus42[nCscRechitClusters]/I");
    tree_->Branch("cscRechitClusterNRechitChamberMinus11",             cscRechitClusterNRechitChamberMinus11,             "cscRechitClusterNRechitChamberMinus11[nCscRechitClusters]/I");
    tree_->Branch("cscRechitClusterNRechitChamberMinus12",             cscRechitClusterNRechitChamberMinus12,             "cscRechitClusterNRechitChamberMinus12[nCscRechitClusters]/I");
    tree_->Branch("cscRechitClusterNRechitChamberMinus13",             cscRechitClusterNRechitChamberMinus13,             "cscRechitClusterNRechitChamberMinus13[nCscRechitClusters]/I");
    tree_->Branch("cscRechitClusterNRechitChamberMinus21",             cscRechitClusterNRechitChamberMinus21,             "cscRechitClusterNRechitChamberMinus21[nCscRechitClusters]/I");
    tree_->Branch("cscRechitClusterNRechitChamberMinus22",             cscRechitClusterNRechitChamberMinus22,             "cscRechitClusterNRechitChamberMinus22[nCscRechitClusters]/I");
    tree_->Branch("cscRechitClusterNRechitChamberMinus31",             cscRechitClusterNRechitChamberMinus31,             "cscRechitClusterNRechitChamberMinus31[nCscRechitClusters]/I");
    tree_->Branch("cscRechitClusterNRechitChamberMinus32",             cscRechitClusterNRechitChamberMinus32,             "cscRechitClusterNRechitChamberMinus32[nCscRechitClusters]/I");
    tree_->Branch("cscRechitClusterNRechitChamberMinus41",             cscRechitClusterNRechitChamberMinus41,             "cscRechitClusterNRechitChamberMinus41[nCscRechitClusters]/I");
    tree_->Branch("cscRechitClusterNRechitChamberMinus42",             cscRechitClusterNRechitChamberMinus42,             "cscRechitClusterNRechitChamberMinus42[nCscRechitClusters]/I");
    tree_->Branch("cscRechitClusterMet_dPhi",             cscRechitClusterMet_dPhi,             "cscRechitClusterMet_dPhi[nCscRechitClusters]/F");




    tree_->Branch("nDtRechitClusters",             &nDtRechitClusters, "nDtRechitClusters/I");

    tree_->Branch("dtRechitClusterNSegStation1",             dtRechitClusterNSegStation1,             "dtRechitClusterNSegStation1[nDtRechitClusters]/I");
    tree_->Branch("dtRechitClusterNSegStation2",             dtRechitClusterNSegStation2,             "dtRechitClusterNSegStation2[nDtRechitClusters]/I");
    tree_->Branch("dtRechitClusterNSegStation3",             dtRechitClusterNSegStation3,             "dtRechitClusterNSegStation3[nDtRechitClusters]/I");
    tree_->Branch("dtRechitClusterNSegStation4",             dtRechitClusterNSegStation4,             "dtRechitClusterNSegStation4[nDtRechitClusters]/I");

    tree_->Branch("dtRechitClusterNOppositeSegStation1",             dtRechitClusterNOppositeSegStation1,             "dtRechitClusterNOppositeSegStation1[nDtRechitClusters]/I");
    tree_->Branch("dtRechitClusterNOppositeSegStation2",             dtRechitClusterNOppositeSegStation2,             "dtRechitClusterNOppositeSegStation2[nDtRechitClusters]/I");
    tree_->Branch("dtRechitClusterNOppositeSegStation3",             dtRechitClusterNOppositeSegStation3,             "dtRechitClusterNOppositeSegStation3[nDtRechitClusters]/I");
    tree_->Branch("dtRechitClusterNOppositeSegStation4",             dtRechitClusterNOppositeSegStation4,             "dtRechitClusterNOppositeSegStation4[nDtRechitClusters]/I");

        tree_->Branch("dtRechitCluster_match_gLLP",             dtRechitCluster_match_gLLP,             "dtRechitCluster_match_gLLP[nDtRechitClusters]/O");
        tree_->Branch("dtRechitCluster_match_gLLP_minDeltaR",             dtRechitCluster_match_gLLP_minDeltaR,             "dtRechitCluster_match_gLLP_minDeltaR[nDtRechitClusters]/F");
        tree_->Branch("dtRechitCluster_match_gLLP_index",             dtRechitCluster_match_gLLP_index,             "dtRechitCluster_match_gLLP_index[nDtRechitClusters]/I");
        tree_->Branch("dtRechitCluster_match_gLLP_eta",             dtRechitCluster_match_gLLP_eta, "dtRechitCluster_match_gLLP_eta[nDtRechitClusters]/F");
        tree_->Branch("dtRechitCluster_match_gLLP_phi",             dtRechitCluster_match_gLLP_phi, "dtRechitCluster_match_gLLP_phi[nDtRechitClusters]/F");
        tree_->Branch("dtRechitCluster_match_gLLP_decay_r",             dtRechitCluster_match_gLLP_decay_r, "dtRechitCluster_match_gLLP_decay_r[nDtRechitClusters]/F");
        tree_->Branch("dtRechitCluster_match_gLLP_decay_z",             dtRechitCluster_match_gLLP_decay_z, "dtRechitCluster_match_gLLP_decay_z[nDtRechitClusters]/F");
        tree_->Branch("dtRechitCluster_match_gLLP_csc",             dtRechitCluster_match_gLLP_csc, "dtRechitCluster_match_gLLP_csc[nDtRechitClusters]/O");
        tree_->Branch("dtRechitCluster_match_gLLP_dt",             dtRechitCluster_match_gLLP_dt, "dtRechitCluster_match_gLLP_dt[nDtRechitClusters]/O");
        tree_->Branch("dtRechitCluster_match_gLLP_e",             dtRechitCluster_match_gLLP_e, "dtRechitCluster_match_gLLP_e[nDtRechitClusters]/F");


        tree_->Branch("dtRechitClusterX",             dtRechitClusterX,             "dtRechitClusterX[nDtRechitClusters]/F");
        tree_->Branch("dtRechitClusterY",             dtRechitClusterY,             "dtRechitClusterY[nDtRechitClusters]/F");
        tree_->Branch("dtRechitClusterZ",             dtRechitClusterZ,             "dtRechitClusterZ[nDtRechitClusters]/F");

        tree_->Branch("dtRechitClusterWheel",             dtRechitClusterWheel,             "dtRechitClusterWheel[nDtRechitClusters]/I");


        tree_->Branch("dtRechitClusterPhi",             dtRechitClusterPhi,             "dtRechitClusterPhi[nDtRechitClusters]/F");
        tree_->Branch("dtRechitClusterEta",             dtRechitClusterEta,             "dtRechitClusterEta[nDtRechitClusters]/F");


        tree_->Branch("dtRechitClusterJetVetoPt",             dtRechitClusterJetVetoPt,             "dtRechitClusterJetVetoPt[nDtRechitClusters]/F");
        tree_->Branch("dtRechitClusterJetVetoLooseId",             dtRechitClusterJetVetoLooseId,             "dtRechitClusterJetVetoLooseId[nDtRechitClusters]/O");
        tree_->Branch("dtRechitClusterJetVetoTightId",             dtRechitClusterJetVetoTightId,             "dtRechitClusterJetVetoTightId[nDtRechitClusters]/O");

        tree_->Branch("dtRechitClusterJetVetoE",             dtRechitClusterJetVetoE,             "dtRechitClusterJetVetoE[nDtRechitClusters]/F");
        tree_->Branch("dtRechitClusterMuonVetoPt",             dtRechitClusterMuonVetoPt,             "dtRechitClusterMuonVetoPt[nDtRechitClusters]/F");
        tree_->Branch("dtRechitClusterMuonVetoE",             dtRechitClusterMuonVetoE,             "dtRechitClusterMuonVetoE[nDtRechitClusters]/F");


        tree_->Branch("dtRechitClusterMuonVetoTightId",             dtRechitClusterMuonVetoTightId,             "dtRechitClusterMuonVetoTightId[nDtRechitClusters]/O");
        tree_->Branch("dtRechitClusterMuonVetoLooseId",             dtRechitClusterMuonVetoLooseId,             "dtRechitClusterMuonVetoLooseId[nDtRechitClusters]/O");
        tree_->Branch("dtRechitClusterMuonVetoGlobal",             dtRechitClusterMuonVetoGlobal,             "dtRechitClusterMuonVetoGlobal[nDtRechitClusters]/O");





        tree_->Branch("dtRechitClusterOverlap",             dtRechitClusterOverlap,             "dtRechitClusterOverlap[nDtRechitClusters]/O");

        tree_->Branch("dtRechitClusterSize",             dtRechitClusterSize,             "dtRechitClusterSize[nDtRechitClusters]/I");
        tree_->Branch("dtRechitClusterNoiseHit",             dtRechitClusterNoiseHit,             "dtRechitClusterNoiseHit[nDtRechitClusters]/I");
        tree_->Branch("dtRechitClusterNoiseHitStation1",             dtRechitClusterNoiseHitStation1,             "dtRechitClusterNoiseHitStation1[nDtRechitClusters]/I");
        tree_->Branch("dtRechitClusterNoiseHitStation2",             dtRechitClusterNoiseHitStation2,             "dtRechitClusterNoiseHitStation2[nDtRechitClusters]/I");
        tree_->Branch("dtRechitClusterNoiseHitStation3",             dtRechitClusterNoiseHitStation3,             "dtRechitClusterNoiseHitStation3[nDtRechitClusters]/I");
        tree_->Branch("dtRechitClusterNoiseHitStation4",             dtRechitClusterNoiseHitStation4,             "dtRechitClusterNoiseHitStation4[nDtRechitClusters]/I");
        tree_->Branch("dtRechitClusterNStation10",             dtRechitClusterNStation10,             "dtRechitClusterNStation10[nDtRechitClusters]/I");
        tree_->Branch("dtRechitClusterAvgStation10",             dtRechitClusterAvgStation10,             "dtRechitClusterAvgStation10[nDtRechitClusters]/F");


        tree_->Branch("dtRechitClusterMaxStation",             dtRechitClusterMaxStation,             "dtRechitClusterMaxStation[nDtRechitClusters]/I");
        tree_->Branch("dtRechitClusterMaxStationRatio",             dtRechitClusterMaxStationRatio,             "dtRechitClusterMaxStationRatio[nDtRechitClusters]/F");
        tree_->Branch("dtRechitClusterNChamber",             dtRechitClusterNChamber,             "dtRechitClusterNChamber[nDtRechitClusters]/I");
        tree_->Branch("dtRechitClusterMaxChamber",             dtRechitClusterMaxChamber,             "dtRechitClusterMaxChamber[nDtRechitClusters]/I");
        tree_->Branch("dtRechitClusterMaxChamberRatio",             dtRechitClusterMaxChamberRatio,             "dtRechitClusterMaxChamberRatio[nDtRechitClusters]/F");
        tree_->Branch("dtRechitClusterNHitStation1",             dtRechitClusterNHitStation1,             "dtRechitClusterNHitStation1[nDtRechitClusters]/I");
        tree_->Branch("dtRechitClusterNHitStation2",             dtRechitClusterNHitStation2,             "dtRechitClusterNHitStation2[nDtRechitClusters]/I");
        tree_->Branch("dtRechitClusterNHitStation3",             dtRechitClusterNHitStation3,             "dtRechitClusterNHitStation3[nDtRechitClusters]/I");
        tree_->Branch("dtRechitClusterNHitStation4",             dtRechitClusterNHitStation4,             "dtRechitClusterNHitStation4[nDtRechitClusters]/I");
        tree_->Branch("dtRechitClusterNHitWheel0",             dtRechitClusterNHitWheel0,             "dtRechitClusterNHitWheel0[nDtRechitClusters]/I");
        tree_->Branch("dtRechitClusterNHitWheel1",             dtRechitClusterNHitWheel0,             "dtRechitClusterNHitWheel1[nDtRechitClusters]/I");
        tree_->Branch("dtRechitClusterNHitWheel2",             dtRechitClusterNHitWheel0,             "dtRechitClusterNHitWheel2[nDtRechitClusters]/I");
        tree_->Branch("dtRechitClusterMet_dPhi",             dtRechitClusterMet_dPhi,             "dtRechitClusterMet_dPhi[nDtRechitClusters]/F");



        tree_->Branch("dtRechitCluster_match_RPChits_dPhi0p5",             dtRechitCluster_match_RPChits_dPhi0p5,             "dtRechitCluster_match_RPChits_dPhi0p5[nDtRechitClusters]/I");
        tree_->Branch("dtRechitCluster_match_RPCBx_dPhi0p5",             dtRechitCluster_match_RPCBx_dPhi0p5,             "dtRechitCluster_match_RPCBx_dPhi0p5[nDtRechitClusters]/I");
        tree_->Branch("dtRechitCluster_match_RB1_0p4",             dtRechitCluster_match_RB1_0p4,             "dtRechitCluster_match_RB1_0p4[nDtRechitClusters]/I");
        tree_->Branch("dtRechitCluster_match_RB1_dPhi0p5",             dtRechitCluster_match_RB1_dPhi0p5,             "dtRechitCluster_match_RB1_dPhi0p5[nDtRechitClusters]/I");
        tree_->Branch("dtRechitCluster_match_MB1Seg_0p4",             dtRechitCluster_match_MB1Seg_0p4,             "dtRechitCluster_match_MB1Seg_0p4[nDtRechitClusters]/I");
        tree_->Branch("dtRechitCluster_match_MB1Seg_0p5",             dtRechitCluster_match_MB1Seg_0p5,             "dtRechitCluster_match_MB1Seg_0p5[nDtRechitClusters]/I");
        tree_->Branch("dtRechitCluster_match_MB1hits_0p4",             dtRechitCluster_match_MB1hits_0p4,             "dtRechitCluster_match_MB1hits_0p4[nDtRechitClusters]/I");
        tree_->Branch("dtRechitCluster_match_MB1hits_0p5",             dtRechitCluster_match_MB1hits_0p5,             "dtRechitCluster_match_MB1hits_0p5[nDtRechitClusters]/I");
        tree_->Branch("dtRechitCluster_match_MB1hits_cosmics_plus",             dtRechitCluster_match_MB1hits_cosmics_plus,             "dtRechitCluster_match_MB1hits_cosmics_plus[nDtRechitClusters]/I");
        tree_->Branch("dtRechitCluster_match_MB1hits_cosmics_minus",             dtRechitCluster_match_MB1hits_cosmics_minus,             "dtRechitCluster_match_MB1hits_cosmics_minus[nDtRechitClusters]/I");

    tree_->Branch("dtRechitClusternXY",             dtRechitClusternXY,             "dtRechitClusternXY[nDtRechitClusters]/I");
    tree_->Branch("dtRechitClusternZ",             dtRechitClusternZ,             "dtRechitClusternZ[nDtRechitClusters]/I");
    tree_->Branch("dtRechitClusterXSpread",             dtRechitClusterXSpread,             "dtRechitClusterXSpread[nDtRechitClusters]/F");
    tree_->Branch("dtRechitClusterYSpread",             dtRechitClusterYSpread,             "dtRechitClusterYSpread[nDtRechitClusters]/F");
    tree_->Branch("dtRechitClusterZSpread",             dtRechitClusterZSpread,             "dtRechitClusterZSpread[nDtRechitClusters]/F");
    tree_->Branch("dtRechitClusterXYSpread",             dtRechitClusterXYSpread,             "dtRechitClusterXYSpread[nDtRechitClusters]/F");
    tree_->Branch("dtRechitClusterRSpread",             dtRechitClusterRSpread,             "dtRechitClusterRSpread[nDtRechitClusters]/F");
    tree_->Branch("dtRechitClusterEtaPhiSpread",             dtRechitClusterEtaPhiSpread,             "dtRechitClusterEtaPhiSpread[nDtRechitClusters]/F");
    tree_->Branch("dtRechitClusterEtaSpread",             dtRechitClusterEtaSpread,             "dtRechitClusterEtaSpread[nDtRechitClusters]/F");
    tree_->Branch("dtRechitClusterPhiSpread",             dtRechitClusterPhiSpread,             "dtRechitClusterPhiSpread[nDtRechitClusters]/F");
    tree_->Branch("dtRechitClusterDeltaRSpread",             dtRechitClusterDeltaRSpread,             "dtRechitClusterDeltaRSpread[nDtRechitClusters]/F");
    tree_->Branch("dtRechitClusterMajorAxis",             dtRechitClusterMajorAxis,             "dtRechitClusterMajorAxis[nDtRechitClusters]/F");
    tree_->Branch("dtRechitClusterMinorAxis",             dtRechitClusterMinorAxis,             "dtRechitClusterMinorAxis[nDtRechitClusters]/F");
    tree_->Branch("dtRechitClusterSkewX",             dtRechitClusterSkewX,             "dtRechitClusterSkewX[nDtRechitClusters]/F");
    tree_->Branch("dtRechitClusterSkewY",             dtRechitClusterSkewY,             "dtRechitClusterSkewY[nDtRechitClusters]/F");
    tree_->Branch("dtRechitClusterSkewZ",             dtRechitClusterSkewZ,             "dtRechitClusterSkewZ[nDtRechitClusters]/F");
    tree_->Branch("dtRechitClusterKurtX",             dtRechitClusterKurtX,             "dtRechitClusterKurtX[nDtRechitClusters]/F");
    tree_->Branch("dtRechitClusterKurtY",             dtRechitClusterKurtY,             "dtRechitClusterKurtY[nDtRechitClusters]/F");
    tree_->Branch("dtRechitClusterKurtZ",             dtRechitClusterKurtZ,             "dtRechitClusterKurtZ[nDtRechitClusters]/F");


  //gLLP branches
  tree_->Branch("nGLLP",          &nGLLP,          "nGLLP/I");


  tree_->Branch("gLLP_eta",          gLLP_eta,          "gLLP_eta[nGLLP]/F");
  tree_->Branch("gLLP_phi",          gLLP_phi,          "gLLP_phi[nGLLP]/F");
  tree_->Branch("gLLP_csc",          gLLP_csc,          "gLLP_csc[nGLLP]/F");
  tree_->Branch("gLLP_dt",          gLLP_dt,          "gLLP_dt[nGLLP]/F");
  tree_->Branch("gLLP_beta",          gLLP_beta,          "gLLP_beta[nGLLP]/F");



  tree_->Branch("gLLP_e",          gLLP_e,          "gLLP_e[nGLLP]/F");
  tree_->Branch("gLLP_pt",          gLLP_pt,          "gLLP_pt[nGLLP]/F");


  tree_->Branch("gLLP_ctau",          gLLP_ctau,          "gLLP_ctau[nGLLP]/F");

  tree_->Branch("gLLP_decay_vertex_r",          gLLP_decay_vertex_r,          "gLLP_decay_vertex_r[nGLLP]/F");
  tree_->Branch("gLLP_decay_vertex_x",          gLLP_decay_vertex_x,          "gLLP_decay_vertex_x[nGLLP]/F");
  tree_->Branch("gLLP_decay_vertex_y",          gLLP_decay_vertex_y,          "gLLP_decay_vertex_y[nGLLP]/F");
  tree_->Branch("gLLP_decay_vertex_z",          gLLP_decay_vertex_z,          "gLLP_decay_vertex_z[nGLLP]/F");


  //leptons
  tree_->Branch("nLeptons",  &nLeptons, "nLeptons/I");
  tree_->Branch("lepE",      lepE,      "lepE[nLeptons]/F");
  tree_->Branch("lepPt",     lepPt,     "lepPt[nLeptons]/F");
  tree_->Branch("lepEta",    lepEta,    "lepEta[nLeptons]/F");
  tree_->Branch("lepPhi",    lepPhi,    "lepPhi[nLeptons]/F");
  tree_->Branch("lepPdgId",  lepPdgId,  "lepPdgId[nLeptons]/I");
  tree_->Branch("lepDZ",     lepDZ,     "lepDZ[nLeptons]/F");
  tree_->Branch("lepTightId", lepTightId, "lepTightId[nLeptons]/O");

  tree_->Branch("lepPassLooseIso", lepPassLooseIso, "lepPassLooseIso[nLeptons]/O");
  tree_->Branch("lepPassTightIso", lepPassTightIso, "lepPassTightIso[nLeptons]/O");
  tree_->Branch("lepPassVTightIso", lepPassVTightIso, "lepPassVTightIso[nLeptons]/O");
  tree_->Branch("lepPassVVTightIso", lepPassVVTightIso, "lepPassVVTightIso[nLeptons]/O");


  //jets
  tree_->Branch("nJets",     &nJets,    "nJets/I");
  tree_->Branch("jetE",      jetE,      "jetE[nJets]/F");
  tree_->Branch("jetPt",     jetPt,     "jetPt[nJets]/F");
  tree_->Branch("jetEta",    jetEta,    "jetEta[nJets]/F");
  tree_->Branch("jetPhi",    jetPhi,    "jetPhi[nJets]/F");
  tree_->Branch("jetTightPassId", jetTightPassId, "jetTightPassId[nJets]/O");
  tree_->Branch("HLTDecision", HLTDecision, "HLTDecision[1201]/O"); //hardcoded

};
