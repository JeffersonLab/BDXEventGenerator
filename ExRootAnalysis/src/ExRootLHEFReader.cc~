
/** \class ExRootLHEFReader
 *
 *  Reads LHEF file
 *
 *
 *  $Date$
 *  $Revision$
 *
 *
 *  \author P. Demin - UCL, Louvain-la-Neuve
 *
 */

#include "ExRootAnalysis/ExRootLHEFReader.h"

#include <stdexcept>
#include <iostream>
#include <sstream>

#include <stdio.h>

#include "TLorentzVector.h"

#include "ExRootAnalysis/ExRootClasses.h"
#include "ExRootAnalysis/ExRootFactory.h"
#include "ExRootAnalysis/ExRootStream.h"

#include "ExRootAnalysis/ExRootTreeBranch.h"

using namespace std;

static const int kBufferSize  = 1024;

//---------------------------------------------------------------------------

ExRootLHEFReader::ExRootLHEFReader() :
  fInputFile(0), fBuffer(0),
  fEventReady(kFALSE), fEventCounter(-1), fParticleCounter(-1)
{
  fBuffer = new char[kBufferSize];
}

//---------------------------------------------------------------------------

ExRootLHEFReader::~ExRootLHEFReader()
{
  if(fBuffer) delete[] fBuffer;
}

//---------------------------------------------------------------------------

void ExRootLHEFReader::SetInputFile(FILE *inputFile)
{
  fInputFile = inputFile;
}

//---------------------------------------------------------------------------

void ExRootLHEFReader::Clear()
{
  fEventReady = kFALSE;
  fEventCounter = -1;
  fParticleCounter = -1;
  fRwgtList.clear();
}

//---------------------------------------------------------------------------

bool ExRootLHEFReader::EventReady()
{
  return fEventReady;
}

//---------------------------------------------------------------------------

bool ExRootLHEFReader::ReadBlock(ExRootTreeBranch *branch)
{
  int rc;
  char *pch;
  double weight;

  if(!fgets(fBuffer, kBufferSize, fInputFile)) return kFALSE;

  if(strstr(fBuffer, "<event>"))
  {
    Clear();
    fEventCounter = 1;
  }
  else if(fEventCounter > 0)
  {
    ExRootStream bufferStream(fBuffer);

    rc = bufferStream.ReadInt(fNparticles)
      && bufferStream.ReadInt(fProcessID)
      && bufferStream.ReadDbl(fWeight)
      && bufferStream.ReadDbl(fScalePDF)
      && bufferStream.ReadDbl(fAlphaQED)
      && bufferStream.ReadDbl(fAlphaQCD);

    if(!rc)
    {
      cerr << "** ERROR: " << "invalid event format" << endl;
      return kFALSE;
    }
    fParticleCounter = fNparticles;
    --fEventCounter;
  }
  else if(fParticleCounter > 0)
  {
    ExRootStream bufferStream(fBuffer);

    rc = bufferStream.ReadInt(fPID)
      && bufferStream.ReadInt(fStatus)
      && bufferStream.ReadInt(fM1)
      && bufferStream.ReadInt(fM2)
      && bufferStream.ReadInt(fC1)
      && bufferStream.ReadInt(fC2)
      && bufferStream.ReadDbl(fPx)
      && bufferStream.ReadDbl(fPy)
      && bufferStream.ReadDbl(fPz)
      && bufferStream.ReadDbl(fE)
      && bufferStream.ReadDbl(fMass)
      && bufferStream.ReadDbl(fLifeTime)
      && bufferStream.ReadDbl(fSpin);

    if(!rc)
    {
      cerr << "** ERROR: " << "invalid particle format" << endl;
      return kFALSE;
    }

    AnalyzeParticle(branch);

    --fParticleCounter;
  }
  else if(strstr(fBuffer, "<wgt"))
  {
    pch = strstr(fBuffer, ">");
    if(!pch)
    {
      cerr << "** ERROR: " << "invalid weight format" << endl;
      return kFALSE;
    }

    ExRootStream bufferStream(pch + 1);
    rc = bufferStream.ReadDbl(weight);

    if(!rc)
    {
      cerr << "** ERROR: " << "invalid weight format" << endl;
      return kFALSE;
    }

    fRwgtList.push_back(weight);
  }
  else if(strstr(fBuffer, "</event>"))
  {
    fEventReady = kTRUE;
  }

  return kTRUE;
}

//---------------------------------------------------------------------------

void ExRootLHEFReader::AnalyzeEvent(ExRootTreeBranch *branch, long long eventNumber)
{
  TRootLHEFEvent *element;

  element = static_cast<TRootLHEFEvent *>(branch->NewEntry());

  element->Number = eventNumber;

  element->Nparticles = fNparticles;
  element->ProcessID = fProcessID;
  element->Weight = fWeight;
  element->ScalePDF = fScalePDF;
  element->CouplingQED = fAlphaQED;
  element->CouplingQCD = fAlphaQCD;
}

//---------------------------------------------------------------------------

void ExRootLHEFReader::AnalyzeRwgt(ExRootTreeBranch *branch)
{
  TRootWeight *element;
  vector<double>::const_iterator itRwgtList;

  for(itRwgtList = fRwgtList.begin(); itRwgtList != fRwgtList.end(); ++itRwgtList)
  {
    element = static_cast<TRootWeight *>(branch->NewEntry());

    element->Weight = *itRwgtList;
  }
}

//---------------------------------------------------------------------------

void ExRootLHEFReader::AnalyzeParticle(ExRootTreeBranch *branch)
{
  TRootLHEFParticle *element;

  Double_t signPz, cosTheta;
  TLorentzVector momentum;

  element = static_cast<TRootLHEFParticle*>(branch->NewEntry());

  element->PID = fPID;
  element->Status = fStatus;

  element->Mother1 = fM1 - 1;
  element->Mother2 = fM2 - 1;

  element->ColorLine1 = fC1;
  element->ColorLine2 = fC2;

  element->Px = fPx;
  element->Py = fPy;
  element->Pz = fPz;
  element->E = fE;
  element->M = fMass;

  momentum.SetPxPyPzE(fPx, fPy, fPz, fE);

  cosTheta = TMath::Abs(momentum.CosTheta());
  signPz = (momentum.Pz() >= 0.0) ? 1.0 : -1.0;

  element->PT = momentum.Perp();
  element->Phi = momentum.Phi();
  element->Eta = (cosTheta == 1.0 ? signPz*999.9 : momentum.Eta());
  element->Rapidity = (cosTheta == 1.0 ? signPz*999.9 : momentum.Rapidity());

  element->LifeTime = fLifeTime;
  element->Spin = fSpin;
}

//---------------------------------------------------------------------------
