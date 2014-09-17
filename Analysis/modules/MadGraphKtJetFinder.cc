
#include "modules/MadGraphKtJetFinder.h"


#include "ExRootAnalysis/ExRootClasses.h"
#include "ExRootAnalysis/ExRootFactory.h"

#include "KtJet/KtEvent.h"
#include "KtJet/KtLorentzVector.h"

#include "TClonesArray.h"

#include "TH1.h"
#include "TH2.h"
#include "TString.h"
#include "TCanvas.h"
#include "TLorentzVector.h"

#include <iostream>
#include <vector>

using namespace std;
using namespace KtJet;

//------------------------------------------------------------------------------

MadGraphKtJetFinder::MadGraphKtJetFinder() :
  fItParticle(0)
{

}

//------------------------------------------------------------------------------

MadGraphKtJetFinder::~MadGraphKtJetFinder()
{
  if(fItParticle) delete fItParticle;
}

//------------------------------------------------------------------------------

void MadGraphKtJetFinder::Init()
{
  
  // define KtJet algorithm
  
  fMaxParticleEta = GetDouble("MaxParticleEta", 5.0);
  fMinParticlePT = GetDouble("ParticleThreshold", 0.5);
  fMinJetPT = GetDouble("JetThreshold", 5.0);

  fCollisionType = GetInt("CollisionType", 4); // PP
  fDistanceScheme = GetInt("DistanceScheme", 1); // Angular
  fRecombinationScheme = GetInt("RecombinationScheme", 1); // E
  fParameterR = GetDouble("ParameterR", 1.0);

  // import ROOT tree branch

  fBranchParticle = UseBranch("GenParticle");
  fItParticle = fBranchParticle->MakeIterator();
  
  // create output arrays

  fOutputArray = ExportArray("jets");

}

//------------------------------------------------------------------------------

void MadGraphKtJetFinder::Finish()
{

}

//------------------------------------------------------------------------------

void MadGraphKtJetFinder::Event()
{
  TRootGenJet *jet;
  TRootGenParticle *particle;
  
  ExRootFactory *factory = GetFactory();

  fTowersList.clear();

  // loop over all particles in event and select stable ones
  fItParticle->Reset();
  while((particle = (TRootGenParticle*) fItParticle->Next()))
  {
 	  if((particle->Status == 1) &&
	     (TMath::Abs(particle->PID) != 12) &&
	     (TMath::Abs(particle->PID) != 14) &&
	     (TMath::Abs(particle->PID) != 16) &&
	     (TMath::Abs(particle->Eta) < fMaxParticleEta) &&
       (particle->PT > fMinParticlePT))
	  {
      fTowersList.push_back(KtLorentzVector(particle->Px, particle->Py,
                                            particle->Pz, particle->E));
    }
  }

  // construct jets from a list of stable particles
  KtEvent event(fTowersList, fCollisionType, fDistanceScheme,
                fRecombinationScheme, fParameterR);

  fJetsList.clear();
  fJetsList = event.getJetsPt();

  Double_t signPz;

  // loop over all jets and export them
  vector<KtLorentzVector>::iterator itJet;
  for(itJet = fJetsList.begin(); itJet != fJetsList.end(); ++itJet)
  {
  	jet = factory->New<TRootGenJet>();

  	jet->E = itJet->e();
  	jet->Px = itJet->px();
  	jet->Py = itJet->py();
  	jet->Pz = itJet->pz();

  	jet->PT = itJet->perp();
    signPz = (jet->Pz >= 0.0) ? 1.0 : -1.0;
    jet->Eta = jet->PT == 0.0 ? signPz*999.9 : itJet->eta();
    jet->Phi = itJet->phi();

    jet->Rapidity = jet->PT == 0.0 ? signPz*999.9 : itJet->rapidity();

  	jet->Mass = itJet->m();

    fOutputArray->Add(jet);
  }
  
  // sort jets by PT
  TRootGenJet::fgCompare = TComparePT<TRootGenJet>::Instance();
  fOutputArray->Sort();
}

