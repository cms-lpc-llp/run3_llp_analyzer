#include "CACluster.h"
#include "TMath.h"
#include <iostream>
#include "TVector3.h"
#include "TGraph.h"
#include "TF1.h"
struct largest_nhit_cluster_
{
  inline bool operator() (const cluster& c1, const cluster& c2){return c1.nhits > c2.nhits;}
} largest_nhit_cluster;

struct hits
{
  float time;
  float error;
  bool strip;
};
const double theWireError_ = 8.6;
const double theStripError_ = 7.0;
const double thePruneCut_ = 9.0;
const int nRechitMin_ = 50;
const int nStationThres_ = 10;
int CACluster::run()
{
  fastjet::JetDefinition jet_def(fastjet::cambridge_algorithm, 0.4);
  std::vector<fastjet::PseudoJet> fjInput;
  vector<Rechits>::iterator iter;

  int recIt = 0;
  for(iter = m_points.begin(); iter != m_points.end(); ++iter)
  {      
    
    float x = (*iter).x;
    float y = (*iter).y;
    float z = (*iter).z;
    float mag = sqrt(x*x + y*y + z*z);
    fjInput.push_back(fastjet::PseudoJet(x, y, z, mag));
    fjInput.back().set_user_index(recIt);
    recIt++;
  }
  fastjet::ClusterSequence clus_seq(fjInput, jet_def);


  //keep all the clusters
  double ptmin = 0.0;
  std::vector<fastjet::PseudoJet> fjJets = clus_seq.inclusive_jets(ptmin);

  int cluster_index = 0;
  // auto clusters = std::make_unique<RecHitClusterCollection>();
  for (auto const& fjJet : fjJets) {
    // skip if the cluster has too few rechits
    if (int(fjJet.constituents().size()) < nRechitMin_) continue;
    cluster tmpCluster;
    // get the constituents from fastjet
    vector<Rechits> rechits;
    for (auto const& constituent : fjJet.constituents()) {
      auto index = constituent.user_index();
      if (index >= 0 && static_cast<unsigned int>(index) < m_points.size()) {
        rechits.push_back(m_points[index]);
      }
    }
   
 

    vector<float> wireTimes;
    vector<float> stripTimes;
    std::vector<hits> cscHits;
    float avg_x(0.0), avg_y(0.0), avg_z(0.0), avg_tWire(0.0), avg_tWirePruned(0.0), avg_t(0.0), avg_tTotal(0.0),tTotalSpreadPruned(0.0);
    float avg_x_sl2(0.0), avg_y_sl2(0.0), avg_z_sl2(0.0);
    float avg_eta(0.0), avg_phi(0.0);
    int size(0), size_z(0), size_xy(0);
    
        // cluster position is the average position of the constituent rechits

    tmpCluster.nhits = rechits.size();

    tmpCluster.tTotal = 0.0;
    for (auto const& rechit : rechits) {
      tmpCluster.tTotal += rechit.twire;
      tmpCluster.tTotal += rechit.t;
    }
    tmpCluster.tTotal /= (2.0 * rechits.size());


    //new timing calculation, error weighted
    // https://github.com/cms-sw/cmssw/blob/master/RecoMuon/MuonIdentification/src/CSCTimingExtractor.cc
    bool modified = false;
    double totalWeightTimeVtx = 0;
    double timeVtx = 0;
    for (auto const& rechit : rechits) 
    {
      if ( rechit.superlayer == 2) //for DT rechits that only have coordinates in Z
      {
        avg_x_sl2 +=  rechit.x;
        avg_y_sl2 +=  rechit.y;
        avg_z_sl2 +=  rechit.z;
        size_z++;
      }
      else if ( rechit.superlayer == 1 ||  rechit.superlayer == 3)
      {
        avg_x +=  rechit.x;
        avg_y +=  rechit.y;
        avg_z +=  rechit.z;
        avg_t +=  rechit.t;
        size_xy ++;
      }
      else //csc or for DT "wrong" rechit coordinates
      {
        avg_x +=  rechit.x;
        avg_y +=  rechit.y;
        avg_z +=  rechit.z;
        avg_t +=  rechit.t;
        avg_tWire +=  rechit.twire;
        wireTimes.push_back( rechit.twire);
        stripTimes.push_back( rechit.t);
        hits thisHit;
        thisHit.time =  rechit.twire;
        thisHit.error = 1./(theWireError_*theWireError_);
        thisHit.strip = false;
        cscHits.push_back(thisHit);
        thisHit.time =  rechit.t;
        thisHit.error = 1./(theStripError_*theStripError_);
        thisHit.strip = true;
        cscHits.push_back(thisHit);

      }
      size ++;

    }

    if (size_xy > 0 && size_z > 0) //for DT correct position, calculate average Z using sl2 and average XY using sl1/3
    {
      avg_x = avg_x/size_xy;
      avg_y = avg_y/size_xy;
      avg_z = avg_z_sl2/size_z;
    }
    else if (size_xy == 0 && size_z == 0) //csc or DT wrong position
    {
      avg_x = avg_x/size;
      avg_y = avg_y/size;
      avg_z = avg_z/size;
      // cout<<avg_x<<","<<avg_y<<","<<avg_z<<endl;
    }
    else if (size_xy > 0 && size_z == 0)
    {
      avg_x = avg_x/size_xy;
      avg_y = avg_y/size_xy;
      avg_z = avg_z/size_xy;

    }
    else
    {
      avg_x = avg_x_sl2/size_z;
      avg_y = avg_y_sl2/size_z;
      avg_z = avg_z_sl2/size_z;

    }

    tmpCluster.x = avg_x;
    tmpCluster.y = avg_y;
    tmpCluster.z = avg_z;

    // calculate cluster eta and phi
    avg_phi = atan(tmpCluster.y/tmpCluster.x);
    if  (tmpCluster.x < 0.0) avg_phi = TMath::Pi() + avg_phi;

    avg_phi = deltaPhi(avg_phi,0.0);
    avg_eta = atan(sqrt(pow(tmpCluster.x,2)+pow(tmpCluster.y,2))/abs(tmpCluster.z));
    avg_eta = -1.0*TMath::Sign(1.0, tmpCluster.z)*log(tan(avg_eta/2));

    tmpCluster.eta = avg_eta;
    tmpCluster.phi = avg_phi;

    do {
      modified = false;
      totalWeightTimeVtx = 0;
      timeVtx = 0;
      for (std::vector<hits>::iterator it = cscHits.begin(); it != cscHits.end(); ++it) {
        timeVtx += it->time * it->error;
        totalWeightTimeVtx += it->error;
      }
      timeVtx /= totalWeightTimeVtx;

      // cut away outliers
      double diff_tvtx;
      double chimax = 0.0;
      int tmmax;
      for (unsigned int i = 0; i < cscHits.size(); i++) {
        diff_tvtx = (cscHits[i].time - timeVtx) * (cscHits[i].time - timeVtx) * cscHits[i].error;

        if (diff_tvtx > chimax) {
          tmmax =  i;
          chimax = diff_tvtx;
        }
      }
      // cut away the outliers
      if (chimax > thePruneCut_) {
        cscHits.erase(cscHits.begin()+tmmax);
        modified = true;
      }
    } while (modified);

    tmpCluster.tWeighted = timeVtx;

    // time spread calculation 
    tmpCluster.TSpread = 0.0;
    tmpCluster.TSpreadWeightedAll = 0.0;
    for (auto const& rechit : rechits) 
    {
      tmpCluster.TSpreadWeightedAll += (rechit.t - tmpCluster.tWeighted) * (rechit.t - tmpCluster.tWeighted);
      tmpCluster.TSpread += (rechit.t - tmpCluster.tTotal) * (rechit.t - tmpCluster.tTotal);
    }
    tmpCluster.TSpread = sqrt(tmpCluster.TSpread/tmpCluster.nhits);
    tmpCluster.TSpreadWeightedAll = sqrt(tmpCluster.TSpreadWeightedAll/tmpCluster.nhits);
    tmpCluster.nCscRechitsChamberPlus11 = 0;
    tmpCluster.nCscRechitsChamberPlus12 = 0;
    tmpCluster.nCscRechitsChamberPlus13 = 0;
    tmpCluster.nCscRechitsChamberPlus21 = 0;
    tmpCluster.nCscRechitsChamberPlus22 = 0;
    tmpCluster.nCscRechitsChamberPlus31 = 0;
    tmpCluster.nCscRechitsChamberPlus32 = 0;
    tmpCluster.nCscRechitsChamberPlus41 = 0;
    tmpCluster.nCscRechitsChamberPlus42 = 0;
    tmpCluster.nCscRechitsChamberMinus11 = 0;
    tmpCluster.nCscRechitsChamberMinus12 = 0;
    tmpCluster.nCscRechitsChamberMinus13 = 0;
    tmpCluster.nCscRechitsChamberMinus21 = 0;
    tmpCluster.nCscRechitsChamberMinus22 = 0;
    tmpCluster.nCscRechitsChamberMinus31 = 0;
    tmpCluster.nCscRechitsChamberMinus32 = 0;
    tmpCluster.nCscRechitsChamberMinus41 = 0;
    tmpCluster.nCscRechitsChamberMinus42 = 0;

    tmpCluster.nDtRechitsStation1 = 0;
    tmpCluster.nDtRechitsStation2 = 0;
    tmpCluster.nDtRechitsStation3 = 0;
    tmpCluster.nDtRechitsStation4 = 0;

  // number of rechits in each chamber/station 
    for (auto const& rechit : rechits) {
      if (rechit.chamber == 11) tmpCluster.nCscRechitsChamberPlus11++;
      if (rechit.chamber == 12) tmpCluster.nCscRechitsChamberPlus12++;
      if (rechit.chamber == 13) tmpCluster.nCscRechitsChamberPlus13++;
      if (rechit.chamber == 21) tmpCluster.nCscRechitsChamberPlus21++;
      if (rechit.chamber == 22) tmpCluster.nCscRechitsChamberPlus22++;
      if (rechit.chamber == 31) tmpCluster.nCscRechitsChamberPlus31++;
      if (rechit.chamber == 32) tmpCluster.nCscRechitsChamberPlus32++;
      if (rechit.chamber == 41) tmpCluster.nCscRechitsChamberPlus41++;
      if (rechit.chamber == 42) tmpCluster.nCscRechitsChamberPlus42++;
      if (rechit.chamber == -11) tmpCluster.nCscRechitsChamberMinus11++;
      if (rechit.chamber == -12) tmpCluster.nCscRechitsChamberMinus12++;
      if (rechit.chamber == -13) tmpCluster.nCscRechitsChamberMinus13++;
      if (rechit.chamber == -21) tmpCluster.nCscRechitsChamberMinus21++;
      if (rechit.chamber == -22) tmpCluster.nCscRechitsChamberMinus22++;
      if (rechit.chamber == -31) tmpCluster.nCscRechitsChamberMinus31++;
      if (rechit.chamber == -32) tmpCluster.nCscRechitsChamberMinus32++;
      if (rechit.chamber == -41) tmpCluster.nCscRechitsChamberMinus41++;
      if (rechit.chamber == -42) tmpCluster.nCscRechitsChamberMinus42++;

      if (abs(rechit.station) == 1) tmpCluster.nDtRechitsStation1++;
      if (abs(rechit.station) == 2) tmpCluster.nDtRechitsStation2++;
      if (abs(rechit.station) == 3) tmpCluster.nDtRechitsStation3++;
      if (abs(rechit.station) == 4) tmpCluster.nDtRechitsStation4++;

    }

    // Nstation avg station
    tmpCluster.nStation10 = 0;
    tmpCluster.avgStation10 = 0;
    int counter = 0;
    int max_count = -999;
    std::map<int, int> station_count_map;
    for (auto const& rechit : rechits) {
      station_count_map[rechit.station]++;
    }
    //station statistics
    std::map<int, int>::iterator it;
    for (auto const& [station, count] : station_count_map) {
      if (count >= nStationThres_) {
        tmpCluster.nStation10++;
        tmpCluster.avgStation10 += station * count;
        counter += count;
      }
      if (count > max_count)
      {
        tmpCluster.maxStation  = station;
        tmpCluster.maxStationRechits = count;
        max_count = count;
      }
    }
    if (counter != 0) tmpCluster.avgStation10 = tmpCluster.avgStation10 / counter;

    clusters.push_back(tmpCluster);

  }

  sort(clusters.begin(), clusters.end(), largest_nhit_cluster);

};




double CACluster::deltaPhi(double phi1, double phi2)
{
  double dphi = phi1-phi2;
  while (dphi > TMath::Pi())
  {
    dphi -= TMath::TwoPi();
  }
  while (dphi <= -TMath::Pi())
  {
    dphi += TMath::TwoPi();
  }
  return dphi;
};
