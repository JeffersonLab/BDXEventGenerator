#ifndef MadGraphConeJetFinder_h
#define MadGraphConeJetFinder_h

#include "ExRootAnalysis/ExRootModule.h"

#include "CDFCones/PhysicsTower.hh"
#include "CDFCones/Cluster.hh"

#include <vector>

class TClonesArray;
class TObjArray;
class TIterator;

class MidPointAlgorithm;

class MadGraphConeJetFinder: public ExRootModule
{
public:

  MadGraphConeJetFinder();
  ~MadGraphConeJetFinder();

  void Init();
  void Event();
  void Finish();
  
private:

  Double_t fMinParticlePT, fMinJetPT;

  std::vector<PhysicsTower> fTowersList;
  std::vector<Cluster> fJetsList;

  MidPointAlgorithm *fJetAlgo; //!

  TIterator *fItParticle; //!

  TClonesArray *fBranchParticle; //!
  
  TObjArray *fOutputArray; //!

  ClassDef(MadGraphConeJetFinder, 1)
};

#endif
