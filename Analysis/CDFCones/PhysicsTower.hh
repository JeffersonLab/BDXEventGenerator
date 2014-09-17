#ifndef _PHYSICS_TOWER_HH_
#define _PHYSICS_TOWER_HH_

#include "LorentzVector.hh"
#include "CalTower.hh"

class PhysicsTower
{
 public:

  LorentzVector fourVector;
  CalTower calTower;

  PhysicsTower(): fourVector(LorentzVector()), calTower(CalTower()) {}
  PhysicsTower(LorentzVector v, CalTower c): fourVector(v), calTower(c) {}
  PhysicsTower(const PhysicsTower& p): fourVector(p.fourVector), calTower(p.calTower) {}
  PhysicsTower(CalTower c):
    fourVector(LorentzVector(c.Et*cos(c.phi),c.Et*sin(c.phi),c.Et*sinh(c.eta),c.Et*cosh(c.eta))), calTower(c) {}
  PhysicsTower(LorentzVector v): fourVector(v), calTower(CalTower(v.Et(),v.eta(),v.phi())) {}
  double Et()   const {return calTower.Et;}
  double eta()  const {return calTower.eta;}
  double phi()  const {return calTower.phi;}
  int    iEta() const {return calTower.iEta;}
  int    iPhi() const {return calTower.iPhi;}
  bool isEqual(PhysicsTower p)
  {
    return fourVector.isEqual(p.fourVector) && calTower.isEqual(p.calTower);
  }
};

#endif
