#ifndef _CLUSTER_HH_
#define _CLUSTER_HH_

#include "PhysicsTower.hh"
#include "LorentzVector.hh"
#include "Centroid.hh"
#include <vector>

class Cluster
{
 public:
  std::vector<PhysicsTower> towerList;
  LorentzVector fourVector;
  Centroid centroid;

  Cluster()
  {
    clear();
  }
  void clear()
  {
    towerList.clear();
    fourVector = LorentzVector();
    centroid = Centroid();
  }
  void addTower(PhysicsTower p)
  {
    towerList.push_back(p);
    fourVector.add(p.fourVector);
    centroid.add(Centroid(p.Et(),p.eta(),p.phi()));
  }
  void removeTower(PhysicsTower p)
  {
    for(std::vector<PhysicsTower>::iterator towerIter = towerList.begin(); towerIter != towerList.end(); towerIter++)
      if(towerIter->isEqual(p)){
	fourVector.subtract(towerIter->fourVector);
	centroid.subtract(Centroid(towerIter->Et(),towerIter->eta(),towerIter->phi()));
	towerList.erase(towerIter);
	break;
      }
  }
  int size(){return towerList.size();}
};

#endif
