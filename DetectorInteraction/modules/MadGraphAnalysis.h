#ifndef MadGraphAnalysis_h
#define MadGraphAnalysis_h

#include "ExRootAnalysis/ExRootModule.h"

#include "TString.h"

#include <map>

class TH1;
class TObjArray;
class TClonesArray;

class MadGraphAnalysis: public ExRootModule
{
public:
    
  MadGraphAnalysis();
  ~MadGraphAnalysis();

  void Init();
  void Event();
  void Finish();
  
private:

  struct ParticleHistograms
  {
    TH1 *fParticlePt;
    TH1 *fParticleRapidity;
  };

  struct PairHistograms
  {
    TH1 *fPairDeltaR;
    TH1 *fPairMass;
  };

  void BookParticleHistograms(ParticleHistograms *histograms,
                              const char *name, const char *title);
  void BookPairHistograms(PairHistograms *histograms,
                          const char *name, const char *title);

  ParticleHistograms *GetParticleHistograms(const char *module, Int_t number);
  PairHistograms *GetPairHistograms(const char *module1, Int_t number1,
                                    const char *module2, Int_t number2);

  TString fOutputFileName; //!

  const TObjArray *fInputArray; //!

  TClonesArray *fBranchEvent; //!

  std::map<TString, ParticleHistograms *> fParticleHistogramsMap; //!
  std::map<TString, PairHistograms *> fPairHistogramsMap; //!

  ClassDef(MadGraphAnalysis, 1)
};

#endif
