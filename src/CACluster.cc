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
  fastjet::JetDefinition jet_def(fastjet::cambridge_algorithm, 0.5);
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

  nClusters = 0;
  // auto clusters = std::make_unique<RecHitClusterCollection>();
  for (auto const& fjJet : fjJets) {
    // skip if the cluster has too few rechits
    if (int(fjJet.constituents().size()) < nRechitMin_) continue;
    // get the constituents from fastjet
    vector<Rechits> rechits;
    for (auto const& constituent : fjJet.constituents()) {
      auto index = constituent.user_index();
      if (index >= 0 && static_cast<unsigned int>(index) < m_points.size()) {
        m_points[index].clusterID = nClusters;
        rechits.push_back(m_points[index]);
      }
    }
   nClusters++;
  }
}
void CACluster::clusterProperties()
{
  

  for (int i = 0; i < nClusters; i++){
    cluster tmpCluster;
    vector<Rechits> rechits;
    for(int j = 0; j < m_points.size(); j++)
    {
      if (m_points[j].clusterID != i) continue;
      rechits.push_back(m_points[j]);
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

    // time spread and momentum calculation 
    float m11(0.0), m12(0.0), m22(0.0), p3_x(0.0), p4_x(0.0), p3_y(0.0), p4_y(0.0), p3_z(0.0), p4_z(0.0);
    float XSpread(0.0), YSpread(0.0), ZSpread(0.0), TSpread(0.0),  TSpreadAll(0.0), XYSpread(0.0), RSpread(0.0), DeltaRSpread(0.0);
    int nXY = 0;
    int nZ = 0;

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


    for (auto const& rechit : rechits) 
    {

      m11 += (rechit.eta-tmpCluster.eta)*(rechit.eta-tmpCluster.eta);
      m12 += (rechit.eta-tmpCluster.eta)* deltaPhi(rechit.phi,tmpCluster.phi);
      m22 += deltaPhi(rechit.phi,tmpCluster.phi)*deltaPhi(rechit.phi,tmpCluster.phi);
      DeltaRSpread +=  pow(deltaR(tmpCluster.eta, tmpCluster.phi, rechit.eta, rechit.phi),2);
      
      // DT superlayer 2 -> measures Z only
      // DT superlayers 1 & 3 -> measure XY only
      // CSCs -> take only one Z value per each station
      
      if(rechit.superlayer > 0){ //hit in DT
	if(rechit.superlayer == 1 || rechit.superlayer == 3){ //XY information
	  XYSpread += (rechit.x - tmpCluster.x)*(rechit.y - tmpCluster.y);
	  XSpread += (rechit.x - tmpCluster.x) * (rechit.x - tmpCluster.x);
	  YSpread += (rechit.y - tmpCluster.y) * (rechit.y - tmpCluster.y);
	  float radius = sqrt(pow(rechit.x, 2) + pow(rechit.y, 2));
	  RSpread += pow(radius-sqrt(tmpCluster.x*tmpCluster.x+tmpCluster.y*tmpCluster.y),2);
	  p3_x += (rechit.x - tmpCluster.x) * (rechit.x - tmpCluster.x) * (rechit.x - tmpCluster.x);
	  p4_x += (rechit.x - tmpCluster.x) * (rechit.x - tmpCluster.x) * (rechit.x - tmpCluster.x) * (rechit.x - tmpCluster.x);
	  p3_y += (rechit.y - tmpCluster.y) * (rechit.y - tmpCluster.y) * (rechit.y - tmpCluster.y);
	  p4_y += (rechit.y - tmpCluster.y) * (rechit.y - tmpCluster.y) * (rechit.y - tmpCluster.y) * (rechit.y - tmpCluster.y);
	  nXY++;
	}
	else if(rechit.superlayer == 2){ //Z information
	  ZSpread += (rechit.z - tmpCluster.z) * (rechit.z - tmpCluster.z);
	  p3_z += (rechit.z - tmpCluster.z) * (rechit.z - tmpCluster.z) * (rechit.z - tmpCluster.z);
	  p4_z += (rechit.z - tmpCluster.z) * (rechit.z - tmpCluster.z) * (rechit.z - tmpCluster.z) * (rechit.z - tmpCluster.z);
	  nZ++;
	}
      }else{ //hit in CSC
	XYSpread += (rechit.x - tmpCluster.x)*(rechit.y - tmpCluster.y);
	XSpread += (rechit.x - tmpCluster.x) * (rechit.x - tmpCluster.x);
	YSpread += (rechit.y - tmpCluster.y) * (rechit.y - tmpCluster.y);
	float radius = sqrt(pow(rechit.x, 2) + pow(rechit.y, 2));
	RSpread += pow(radius-sqrt(tmpCluster.x*tmpCluster.x+tmpCluster.y*tmpCluster.y),2);
	p3_x += (rechit.x - tmpCluster.x) * (rechit.x - tmpCluster.x) * (rechit.x - tmpCluster.x);
	p4_x += (rechit.x - tmpCluster.x) * (rechit.x - tmpCluster.x) * (rechit.x - tmpCluster.x) * (rechit.x - tmpCluster.x);
	p3_y += (rechit.y - tmpCluster.y) * (rechit.y - tmpCluster.y) * (rechit.y - tmpCluster.y);
	p4_y += (rechit.y - tmpCluster.y) * (rechit.y - tmpCluster.y) * (rechit.y - tmpCluster.y) * (rechit.y - tmpCluster.y);
	nXY++;
	ZSpread += (rechit.z - tmpCluster.z) * (rechit.z - tmpCluster.z);
	p3_z += (rechit.z - tmpCluster.z) * (rechit.z - tmpCluster.z) * (rechit.z - tmpCluster.z);
	p4_z += (rechit.z - tmpCluster.z) * (rechit.z - tmpCluster.z) * (rechit.z - tmpCluster.z) * (rechit.z - tmpCluster.z);
	nZ++;
      }

      TSpread += (rechit.t - tmpCluster.tTotal) * (rechit.t - tmpCluster.tTotal);
      TSpreadAll += (rechit.t - tmpCluster.tWeighted) * (rechit.t - tmpCluster.tWeighted);

      // number of rechits in each chamber/station 
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

    float a = (m11+m22)/2;
    float b = 0.5*sqrt((m11+m22)*(m11+m22)-4*(m11*m22-m12*m12));
    tmpCluster.XSpread = sqrt(XSpread/(float)nXY);
    tmpCluster.YSpread = sqrt(YSpread/(float)nXY);
    tmpCluster.ZSpread = sqrt(ZSpread/(float)nZ);
    tmpCluster.RSpread = sqrt(RSpread/(float)nXY);
    tmpCluster.XYSpread = sqrt(abs(XYSpread)/nXY);
    tmpCluster.SkewX = (p3_x/(float)nXY) / pow(XSpread/(float)nXY, 3/2);
    tmpCluster.SkewY = (p3_y/(float)nXY) / pow(YSpread/(float)nXY, 3/2);
    tmpCluster.SkewZ = (p3_z/(float)nZ) / pow(ZSpread/(float)nZ, 3/2);
    tmpCluster.KurtX = (p4_x/(float)nXY) / (pow(XSpread/(float)nXY, 2)) - 3.0;
    tmpCluster.KurtY = (p4_y/(float)nXY) / (pow(YSpread/(float)nXY, 2)) - 3.0;
    tmpCluster.KurtZ = (p4_z/(float)nZ) / (pow(ZSpread/(float)nZ, 2)) - 3.0;
    tmpCluster.nXY = nXY;
    tmpCluster.nZ = nZ;
    tmpCluster.DeltaRSpread = sqrt(DeltaRSpread/(float)tmpCluster.nhits);
    tmpCluster.TSpread = sqrt(TSpread/(float)tmpCluster.nhits);
    tmpCluster.TSpreadWeightedAll = sqrt(TSpreadAll/(float)tmpCluster.nhits);
    tmpCluster.EtaSpread = sqrt(abs(m11)/(float)tmpCluster.nhits);
    tmpCluster.EtaPhiSpread = sqrt(abs(m12)/(float)tmpCluster.nhits);
    tmpCluster.PhiSpread = sqrt(abs(m22)/(float)tmpCluster.nhits);
    tmpCluster.MajorAxis = sqrt((a+b)/(float)tmpCluster.nhits);
    tmpCluster.MinorAxis = sqrt((a-b)/(float)tmpCluster.nhits);


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


};
void CACluster::sort_clusters() //only run sort after merg_clusters, or it will mess up the clusterID of the m_points
{
  sort(clusters.begin(), clusters.end(), largest_nhit_cluster);
}

void CACluster::merge_clusters()
{
  // clear all the cluster variables
  //change cluster ID of points

  // get the list of eta and phi of the clusters
  vector<float> clusterEta;
  vector<float> clusterPhi;
  for(unsigned int j = 0; j < nClusters; j++){
    clusterEta.push_back(clusters.at(j).eta);
    clusterPhi.push_back(clusters.at(j).phi);

  }
  bool modified = true;
  while(modified){
    modified = false;
    float mindR = 15;
    int cluster1 = 999;
    int cluster2 = 999;

    for(unsigned int i = 0; i < nClusters; i++){ //find the min_deltaR between any two clusters
      for(unsigned int j = i+1; j < nClusters; j++){
        float current_dR = deltaR(clusters[i].eta, clusters[i].phi, clusters[j].eta, clusters[j].phi);
        if(current_dR<mindR)
        {
          mindR = current_dR;
          cluster1 = i;
          cluster2 = j;
        }
      }
    }
    if (mindR < CA_MERGE_CLUSTER_DR){ //if min deltaR < the deltaR merging threshold, then merge the two clusters
      vector<Rechits>::iterator iter;
      float avg_x(0.0), avg_y(0.0), avg_z(0.0);
      float avg_x_sl2(0.0), avg_y_sl2(0.0), avg_z_sl2(0.0);
      float avg_eta(0.0), avg_phi(0.0);
      int size(0), size_z(0), size_xy(0);

      
      for(iter = m_points.begin(); iter != m_points.end(); ++iter)
      {
        if ( iter->clusterID == cluster2 ){ //change the clusterID from cluster 2 to cluster1
          iter->clusterID = cluster1;
        }
        if ( iter->clusterID > cluster2 )iter->clusterID = iter->clusterID-1; // change the cluster ID for all clusters beyond to ID-1

        if (iter->clusterID == cluster1) //recalculate the eta/phi position of the new cluster
        {

            if ( iter->superlayer == 2) //for DT rechits that only have coordinates in Z
            {
              avg_x_sl2 +=  iter->x;
              avg_y_sl2 +=  iter->y;
              avg_z_sl2 +=  iter->z;
              size_z++;
            }
            else if ( iter->superlayer == 1 ||  iter->superlayer == 3)
            {
              avg_x +=  iter->x;
              avg_y +=  iter->y;
              avg_z +=  iter->z;
              size_xy ++;
            }
            else //csc or for DT "wrong" rechit coordinates
            {
              avg_x +=  iter->x;
              avg_y +=  iter->y;
              avg_z +=  iter->z;
            
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

          // calculate cluster eta and phi
          avg_phi = atan(avg_y/avg_x);
          if  (avg_x < 0.0) avg_phi = TMath::Pi() + avg_phi;

          avg_phi = deltaPhi(avg_phi,0.0);
          avg_eta = atan(sqrt(pow(avg_x,2)+pow(avg_y,2))/abs(avg_z));
          avg_eta = -1.0*TMath::Sign(1.0, avg_z)*log(tan(avg_eta/2));

      }
      clusterEta.erase(clusterEta.begin() + cluster2);
      clusterPhi.erase(clusterPhi.begin() + cluster2);
      clusterEta[cluster1] = avg_eta;
      clusterPhi[cluster1] = avg_phi;

      nClusters--;
      modified = true;
    }
  }
  clusters.clear();
}


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

double CACluster::deltaR(double eta1, double phi1, double eta2, double phi2) {
  double dphi = deltaPhi(phi1,phi2);
  double deta = eta1 - eta2;
  return sqrt( dphi*dphi + deta*deta);
}
