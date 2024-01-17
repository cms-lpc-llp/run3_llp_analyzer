#ifndef CACluster_H
#define CACluster_H

#include <vector>
#include <cmath>
#include "fastjet/ClusterSequence.hh"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/Selector.hh"
using namespace std;
const double CA_MERGE_CLUSTER_DR = 0.6;



struct cluster
{
  float x, y, z, tTotal, tWeighted, eta, phi;//t is t_strip, tWire, tTotal is sum
  int nhits;
  int maxChamber, maxChamberRechits, nChamber;
  // vector<int> cscChambers;
  int maxStation, maxStationRechits, nStation10;
  float avgStation10;
  int nCscRechitsChamberPlus11, nCscRechitsChamberPlus12, nCscRechitsChamberPlus13, nCscRechitsChamberPlus21, nCscRechitsChamberPlus22, nCscRechitsChamberPlus31, nCscRechitsChamberPlus32, nCscRechitsChamberPlus41, nCscRechitsChamberPlus42;
  int nCscRechitsChamberMinus11, nCscRechitsChamberMinus12, nCscRechitsChamberMinus13, nCscRechitsChamberMinus21, nCscRechitsChamberMinus22, nCscRechitsChamberMinus31, nCscRechitsChamberMinus32, nCscRechitsChamberMinus41, nCscRechitsChamberMinus42;

  int nDtRechitsStation1,nDtRechitsStation2,nDtRechitsStation3,nDtRechitsStation4,nDtRechitsWheel0,nDtRechitsWheel1,nDtRechitsWheel2;

  float Me11Ratio, Me12Ratio;
  float XSpread, YSpread, ZSpread, RSpread, DeltaRSpread, XYSpread, TSpread, TSpreadWeightedAll, EtaSpread, EtaPhiSpread, PhiSpread, MajorAxis, MinorAxis, SkewX, SkewY, SkewZ, KurtX, KurtY, KurtZ;
  int nXY, nZ;

};

typedef struct Rechits_
{
    float x, y, z, t, twire;  // X, Y, Z position
    float eta,phi;
    float dirX, dirY, dirZ;
  int station, chamber, layer, superlayer, wheel; //superlayer and wheel exist only for DT
    int clusterID;  // clustered ID
}Rechits;

class CACluster {
public:
    CACluster(unsigned int minPts, float eps, vector<Rechits> rechits){
        m_minPoints = minPts;
        m_epsilon = eps;
        m_points = rechits;
        m_pointSize = rechits.size();
    }
    ~CACluster(){}
   
    vector<cluster> clusters;
    int nClusters;
    int run();
    vector<Rechits> run_rechits();
    void clusterProperties();
    void sort_clusters();
    void merge_clusters();

    double deltaPhi(double phi1, double phi2);
    double deltaR(double eta1, double phi1, double eta2, double phi2);

private:
    vector<Rechits> m_points;
    unsigned int m_pointSize;
    unsigned int m_minPoints;
    float m_epsilon;
};

#endif // CACluster_H
