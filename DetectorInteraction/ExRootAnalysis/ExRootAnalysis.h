#ifndef ExRootAnalysis_h
#define ExRootAnalysis_h

/** \class ExRootAnalysis
 *
 *  Analysis steering class.
 *  Implements events loop and modules management.
 *
 *  $Date: 2006/09/21 13:08:01 $
 *  $Revision: 1.1 $
 *
 *
 *  \author P. Demin - UCL, Louvain-la-Neuve
 *
 */

#include "ExRootAnalysis/ExRootTask.h"

class TFile;
class TFolder;
class TObjArray;

class ExRootConfReader;
class ExRootTreeReader;
class ExRootTreeWriter;

class ExRootFactory;

class ExRootAnalysis: public ExRootTask 
{
public:

  ExRootAnalysis();
  ~ExRootAnalysis();

  void SetTclFileName(const char *name) { fTclFileName = name; }
  void SetPDGFileName(const char *name) { fPDGFileName = name; }

  Long64_t GetEntries() const;
  Bool_t ReadEvent(Long64_t entry);

  void Loop();

  virtual void ProcessEvent();

  virtual void Init();
  virtual void Event();
  virtual void Finish();

  virtual void Clear();

private:

  TFile *fTreeFile, *fInfoFile;

  TString fTclFileName, fPDGFileName;

  TObjArray *fChains;

  ExRootConfReader *fConfReader;
  ExRootTreeReader *fTreeReader;
  ExRootTreeWriter *fTreeWriter;

  ExRootFactory *fFactory;

  ClassDef(ExRootAnalysis, 1)
};

#endif /* ExRootAnalysis_h */

