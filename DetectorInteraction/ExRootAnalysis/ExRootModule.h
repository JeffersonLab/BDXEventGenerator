#ifndef ExRootModule_h
#define ExRootModule_h

/** \class ExRootModule
 *
 *  Base class for analysis modules
 *
 *  $Date: 2006/09/21 13:08:01 $
 *  $Revision: 1.1 $
 *
 *
 *  \author P. Demin - UCL, Louvain-la-Neuve
 *
 */

#include "ExRootAnalysis/ExRootTask.h"

class TClass;
class TObject;
class TFolder;
class TClonesArray;

class ExRootResult;
class ExRootTreeBranch;
class ExRootTreeReader;
class ExRootTreeWriter;

class ExRootFactory;

class ExRootModule: public ExRootTask 
{
public:

  ExRootModule();
  ~ExRootModule();

  virtual void Init() = 0;
  virtual void Event() = 0;
  virtual void Finish() = 0;

  const TObjArray *ImportArray(const char *name);
  TObjArray *ExportArray(const char *name);

  TClonesArray *UseBranch(const char *name);

  ExRootTreeBranch *NewBranch(const char *name, TClass *cl);

  ExRootResult *GetPlots();
  ExRootFactory *GetFactory();

private:

  TFolder *NewFolder(const char *name);
  TObject *GetObject(const char *name, TClass *cl);

  TFolder *fPlotFolder, *fExportFolder;

  ExRootTreeReader *fTreeReader;
  ExRootTreeWriter *fTreeWriter;

  ExRootFactory *fFactory;

  ExRootResult *fPlots;

  ClassDef(ExRootModule, 1)
};

#endif /* ExRootModule_h */

