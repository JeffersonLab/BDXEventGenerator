
#include "modules/MadGraphConeJetFinder.h"


#include "ExRootAnalysis/ExRootClasses.h"
#include "ExRootAnalysis/ExRootFactory.h"

#include "CDFCones/JetCluAlgorithm.hh"
#include "CDFCones/MidPointAlgorithm.hh"

#include "TClonesArray.h"

#include "TH1.h"
#include "TH2.h"
#include "TString.h"
#include "TCanvas.h"
#include "TLorentzVector.h"

#include <iostream>
#include <vector>

using namespace std;

//------------------------------------------------------------------------------

MadGraphConeJetFinder::MadGraphConeJetFinder() :
  fJetAlgo(0), fItParticle(0)
{

}

//------------------------------------------------------------------------------

MadGraphConeJetFinder::~MadGraphConeJetFinder()
{
  if(fJetAlgo) delete fJetAlgo;
  if(fItParticle) delete fItParticle;
}

//------------------------------------------------------------------------------

void MadGraphConeJetFinder::Init()
{
  
  // define MidPoint algorithm
  
  fMinParticlePT = GetDouble("ParticleThreshold", 0.5);
  fMinJetPT = GetDouble("JetThreshold", 5.0);

  double seedThreshold    = GetDouble("SeedThreshold", 1.0);
  double coneRadius       = GetDouble("ConeRadius", 0.5);
  double coneAreaFraction = GetDouble("ConeAreaFraction", 1.0);
  int    maxPairSize      = GetInt("MaxPairSize", 2);
  int    maxIterations    = GetInt("MaxIterations", 100);
  double overlapThreshold = GetDouble("OverlapThreshold", 0.75);

  fJetAlgo = new MidPointAlgorithm(seedThreshold, coneRadius, coneAreaFraction,
                                   maxPairSize, maxIterations, overlapThreshold);

  // import ROOT tree branch

  fBranchParticle = UseBranch("GenParticle");
  fItParticle = fBranchParticle->MakeIterator();
  
  // create output arrays

  fOutputArray = ExportArray("jets");

}

//------------------------------------------------------------------------------

void MadGraphConeJetFinder::Finish()
{

}

//------------------------------------------------------------------------------

void MadGraphConeJetFinder::Event()
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
       particle->PT > fMinParticlePT)
	  {
      fTowersList.push_back(PhysicsTower(LorentzVector(particle->Px, particle->Py,
                                                       particle->Pz, particle->E)));
    }
  }

  // construct jets from a list of stable particles
  fJetsList.clear();
  fJetAlgo->run(fTowersList, fJetsList);

  Double_t signPz;
  LorentzVector jetMomentum;
  TLorentzVector momentum;

  // loop over all jets and export them
  vector<Cluster>::iterator itJet;
  for(itJet = fJetsList.begin(); itJet != fJetsList.end(); ++itJet)
  {
  	jet = factory->New<TRootGenJet>();

  	jetMomentum = itJet->fourVector;

  	momentum.SetPxPyPzE(jetMomentum.px, jetMomentum.py, jetMomentum.pz, jetMomentum.E);

  	jet->E = jetMomentum.E;
  	jet->Px = jetMomentum.px;
  	jet->Py = jetMomentum.py;
  	jet->Pz = jetMomentum.pz;

  	jet->PT = momentum.Perp();
    signPz = (jet->Pz >= 0.0) ? 1.0 : -1.0;
    jet->Eta = jet->PT == 0.0 ? signPz*999.9 : momentum.Eta();
    jet->Phi = momentum.Phi();

    jet->Rapidity = jet->PT == 0.0 ? signPz*999.9 : momentum.Rapidity();

  	jet->Mass = momentum.M();

    fOutputArray->Add(jet);
  }
  
  // sort jets by PT
  TRootGenJet::fgCompare = TComparePT<TRootGenJet>::Instance();
  fOutputArray->Sort();
}

