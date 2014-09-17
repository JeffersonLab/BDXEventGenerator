#ifndef MadGraphKtJetFinder_h
#define MadGraphKtJetFinder_h

#include "ExRootAnalysis/ExRootModule.h"

#include "KtJet/KtLorentzVector.h"

#include <vector>

class TClonesArray;
class TObjArray;
class TIterator;

class MadGraphKtJetFinder: public ExRootModule
{
public:

  MadGraphKtJetFinder();
  ~MadGraphKtJetFinder();

  void Init();
  void Event();
  void Finish();

private:

  Double_t fMaxParticleEta, fMinParticlePT, fMinJetPT, fParameterR;

  Int_t fCollisionType, fDistanceScheme, fRecombinationScheme;

  std::vector<KtJet::KtLorentzVector> fTowersList; //!
  std::vector<KtJet::KtLorentzVector> fJetsList; //!

  TIterator *fItParticle; //!

  TClonesArray *fBranchParticle; //!
  
  TObjArray *fOutputArray; //!

  ClassDef(MadGraphKtJetFinder, 1)
};

#endif
