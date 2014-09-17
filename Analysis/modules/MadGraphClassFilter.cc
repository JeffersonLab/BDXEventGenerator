
#include "modules/MadGraphClassFilter.h"


#include "ExRootAnalysis/ExRootResult.h"
#include "ExRootAnalysis/ExRootClasses.h"

#include "ExRootAnalysis/ExRootFilter.h"
#include "ExRootAnalysis/ExRootClassifier.h"

#include "TClonesArray.h"

#include <map>
#include <set>
#include <deque>

#include <iostream>

using namespace std;

//------------------------------------------------------------------------------

class MadGraphParticleClassifier : public ExRootClassifier
{
public:
  MadGraphParticleClassifier() :
  fMaxCategories(0) { }
  void InsertParticleStatus(Int_t status);
  void InsertClassPID(const TString &className, Int_t pid);
  Int_t GetCategory(TObject *object);
  Int_t GetMaxCategories() const { return fMaxCategories; }
  TString GetCategoryClassName(Int_t category) const { return fClassNameArray[category]; }
private:

  map< Int_t, Int_t > fPIDMap;
  map< TString, Int_t > fClassNameMap;
  deque< TString > fClassNameArray;
  set< Int_t > fParticleStatusSet;

  Int_t fMaxCategories;
};

//------------------------------------------------------------------------------

void MadGraphParticleClassifier::InsertParticleStatus(Int_t status)
{
  fParticleStatusSet.insert(status);
}

//------------------------------------------------------------------------------

void MadGraphParticleClassifier::InsertClassPID(const TString &className, Int_t pid)
{
  Int_t category;
  map< TString, Int_t >::const_iterator itClassNameMap;

  itClassNameMap = fClassNameMap.find(className);
  
  if(itClassNameMap == fClassNameMap.end())
  {
    category = fMaxCategories;
    fClassNameMap[className] = category;
    fClassNameArray.push_back(className);
    ++fMaxCategories;
  }
  else
  {
    category = itClassNameMap->second;
  }
  fPIDMap[pid] = category;
}

//------------------------------------------------------------------------------

Int_t MadGraphParticleClassifier::GetCategory(TObject *object)
{
  TRootLHEFParticle *particle = static_cast<TRootLHEFParticle*>(object);

  map< Int_t, Int_t >::const_iterator itPIDMap;

  TString className;

  Int_t pidAbs, pid = particle->PID;
  Int_t result = -1;

  if(fParticleStatusSet.find(particle->Status) == fParticleStatusSet.end()) return -1;

  itPIDMap = fPIDMap.find(pid);

  if(itPIDMap != fPIDMap.end())
  {
    result = itPIDMap->second;
  }
  else
  {
    pidAbs = TMath::Abs(pid);
    className = Form("%d", pidAbs);
    result = fMaxCategories;
    InsertClassPID(className, pidAbs);
    InsertClassPID(className, -pidAbs);
  }

  return result;
}

//------------------------------------------------------------------------------

MadGraphClassFilter::MadGraphClassFilter()
{
}

//------------------------------------------------------------------------------

MadGraphClassFilter::~MadGraphClassFilter()
{
}

//------------------------------------------------------------------------------

void MadGraphClassFilter::Init()
{
  TString className;
  ExRootConfParam param, classParticles;

  Int_t i, j, status, pid, sizeParam, sizeParticles;

  // import ROOT tree branch

  fBranchParticle = UseBranch("Particle");

  fItParticle = fBranchParticle->MakeIterator();

  // create classifier and filter

  fClassifier = new MadGraphParticleClassifier();
  fFilter = new ExRootFilter(fBranchParticle);

  // read particle status from configuration file and setup classifier

  param = GetParam("ParticleStatus");
  sizeParam = param.GetSize();

  for(i = 0; i < sizeParam; ++i)
  {
    status = param[i].GetInt();
    fClassifier->InsertParticleStatus(status);
  }

  // read particle classes from configuration file and setup classifier

  param = GetParam("ClassParticles");
  sizeParam = param.GetSize();

  for(i = 0; i < sizeParam/2; ++i)
  {
    className = param[i*2].GetString();
    classParticles = param[i*2 + 1];
    sizeParticles = classParticles.GetSize();

    for(j = 0; j < sizeParticles; ++j)
    {
      pid = classParticles[j].GetInt();
      fClassifier->InsertClassPID(className, pid);
    }
  }

  // create output arrays

  fOutputArray = ExportArray("particles");
}

//------------------------------------------------------------------------------

void MadGraphClassFilter::Finish()
{
  if(fFilter) delete fFilter;
  if(fClassifier) delete fClassifier;
  if(fItParticle) delete fItParticle;
}

//------------------------------------------------------------------------------

void MadGraphClassFilter::Event()
{
  TObjArray *subarray;
  Int_t category;

  fFilter->Reset();

  // make filter classify particles and fill all subarrays
  // at this point classifier creates additional/missing classes
  fFilter->GetSubArray(fClassifier, 0);

  // loop over all classes and export class names and classified particles
  for(category = 0; category < fClassifier->GetMaxCategories(); ++category)
  {
    subarray = fFilter->GetSubArray(fClassifier, category);
    if(subarray)
    {
      subarray->SetName(fClassifier->GetCategoryClassName(category));
      fOutputArray->Add(subarray);
  
      // sort particles by PT
      TRootLHEFParticle::fgCompare = TComparePT<TRootLHEFParticle>::Instance();
      subarray->Sort();
    }
  }
}

