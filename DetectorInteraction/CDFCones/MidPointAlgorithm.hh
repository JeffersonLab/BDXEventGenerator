#ifndef _MIDPOINT_ALGORITHM_HH_
#define _MIDPOINT_ALGORITHM_HH_

#include "PhysicsTower.hh"
#include "Cluster.hh"
#include <vector>

class MidPointAlgorithm
{
 private:
  double _seedThreshold;
  double _coneRadius;
  double _coneAreaFraction;
  int    _maxPairSize;
  int    _maxIterations;
  double _overlapThreshold;

 public:
  MidPointAlgorithm():
    _seedThreshold(1),
    _coneRadius(0.7),
    _coneAreaFraction(0.25),
    _maxPairSize(2),
    _maxIterations(100),
    _overlapThreshold(0.75)
  {}
  MidPointAlgorithm(double st, double cr, double caf, int mps, int mi, double ot):
    _seedThreshold(st),
    _coneRadius(cr),
    _coneAreaFraction(caf),
    _maxPairSize(mps),
    _maxIterations(mi),
    _overlapThreshold(ot)
  {}
  void findStableConesFromSeeds(std::vector<PhysicsTower>& particles, std::vector<Cluster>& stableCones);
  void findStableConesFromMidPoints(std::vector<PhysicsTower>& particles, std::vector<Cluster>& stableCones);
  void iterateCone(double startRapidity, double startPhi, double startPt, std::vector<PhysicsTower>& particles,
		   std::vector<Cluster>& stableCones, bool reduceConeSize);
  void addClustersToPairs(std::vector<int>& testPair, std::vector< std::vector<int> >& pairs,
			  std::vector< std::vector<bool> >& distanceOK, int maxClustersInPair);
  void splitAndMerge(std::vector<Cluster>& stableCones, std::vector<Cluster>& jets);
  void run(std::vector<PhysicsTower>& particles, std::vector<Cluster>& jets);
};

#endif
