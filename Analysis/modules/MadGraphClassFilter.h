#ifndef MadGraphClassFilter_h
#define MadGraphClassFilter_h

#include "ExRootAnalysis/ExRootModule.h"

class TClonesArray;
class TObjArray;
class TIterator;

class ExRootFilter;
class MadGraphParticleClassifier;

class MadGraphClassFilter: public ExRootModule
{
public:
    
  MadGraphClassFilter();
  ~MadGraphClassFilter();

  void Init();
  void Event();
  void Finish();
  
private:

  ExRootFilter *fFilter; //!
  MadGraphParticleClassifier *fClassifier; //!

  TIterator *fItParticle; //!

  TClonesArray *fBranchParticle; //!

  TObjArray *fOutputArray; //!

  ClassDef(MadGraphClassFilter, 1)
};

#endif
