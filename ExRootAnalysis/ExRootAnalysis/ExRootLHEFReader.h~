#ifndef ExRootLHEFReader_h
#define ExRootLHEFReader_h

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

#include <stdio.h>

#include <vector>

class ExRootTreeBranch;
class ExRootFactory;

class ExRootLHEFReader
{
public:

  ExRootLHEFReader();
  ~ExRootLHEFReader();

  void SetInputFile(FILE *inputFile);

  void Clear();
  bool EventReady();

  bool ReadBlock(ExRootTreeBranch *branch);

  void AnalyzeEvent(ExRootTreeBranch *branch, long long eventNumber);

  void AnalyzeRwgt(ExRootTreeBranch *branch);

private:

  void AnalyzeParticle(ExRootTreeBranch *branch);

  FILE *fInputFile;

  char *fBuffer;

  bool fEventReady;

  int fEventCounter;

  int fParticleCounter, fNparticles, fProcessID;
  double fWeight, fScalePDF, fAlphaQCD, fAlphaQED;

  int fPID, fStatus, fM1, fM2, fC1, fC2;
  double fPx, fPy, fPz, fE, fMass, fLifeTime, fSpin;
  
  std::vector<double> fRwgtList;
};

#endif // ExRootLHEFReader_h


