
/** \class ExRootModule
 *
 *  Base class for analysis modules
 *
 *  $Date: 2006/09/21 13:10:52 $
 *  $Revision: 1.1 $
 *
 *
 *  \author P. Demin - UCL, Louvain-la-Neuve
 *
 */

#include "ExRootAnalysis/ExRootTreeReader.h"
#include "ExRootAnalysis/ExRootTreeBranch.h"
#include "ExRootAnalysis/ExRootTreeWriter.h"
#include "ExRootAnalysis/ExRootResult.h"

#include "ExRootAnalysis/ExRootModule.h"
#include "ExRootAnalysis/ExRootFactory.h"

#include "TROOT.h"
#include "TClass.h"
#include "TFolder.h"
#include "TObjArray.h"

#include <iostream>

using namespace std;

ExRootModule::ExRootModule() :
  fPlotFolder(0), fExportFolder(0),
  fTreeReader(0), fTreeWriter(0),
  fFactory(0), fPlots(0)
{
}

//------------------------------------------------------------------------------

ExRootModule::~ExRootModule()
{
}

//------------------------------------------------------------------------------

const TObjArray *ExRootModule::ImportArray(const char *name)
{
  TObjArray *object;

  object = static_cast<TObjArray *>(GetObject(name, TObjArray::Class()));
  if(!object)
  {
    cout << "** ERROR: cannot access input list '" << name << "'" << endl;
    return 0;
  }

  return object;
}

//------------------------------------------------------------------------------

TObjArray *ExRootModule::ExportArray(const char *name)
{
  TObjArray *array;
  if(!fExportFolder)
  {
    fExportFolder = NewFolder("Export");
  }

  array = GetFactory()->NewPermanentArray();

  array->SetName(name);
  fExportFolder->Add(array);

  return array;
}

//------------------------------------------------------------------------------

ExRootTreeBranch *ExRootModule::NewBranch(const char *name, TClass *cl)
{
  if(!fTreeWriter)
  {
    fTreeWriter = static_cast<ExRootTreeWriter *>(GetObject("TreeWriter", ExRootTreeWriter::Class()));
    if(!fTreeWriter)
    {
      cout << "** ERROR: cannot access tree writer" << endl;
      return 0;
    }
  }
  return fTreeWriter->NewBranch(name, cl);
}

//------------------------------------------------------------------------------

TClonesArray *ExRootModule::UseBranch(const char *name)
{
  if(!fTreeReader)
  {
    fTreeReader = static_cast<ExRootTreeReader *>(GetObject("TreeReader", ExRootTreeReader::Class()));
    if(!fTreeReader)
    {
      cout << "** ERROR: cannot access tree reader" << endl;
      return 0;
    }
  }
  return fTreeReader->UseBranch(name);
}

//------------------------------------------------------------------------------

TFolder *ExRootModule::NewFolder(const char *name)
{
  TFolder *folder;
  folder = static_cast<TFolder *>(GetObject(name, TFolder::Class()));
  if(!folder) folder = GetFolder()->AddFolder(name, "");
  if(!folder)
  {
    cout << "** ERROR: cannot create folder '" << name << "'" << endl;
    return 0;
  }
  folder = folder->AddFolder(GetName(), GetTitle());
  if(!folder)
  {
    cout << "** ERROR: cannot create folder '";
    cout << name << "/" << GetName() << "'" << endl;
    return 0;
  }
  return folder;
}

//------------------------------------------------------------------------------

TObject *ExRootModule::GetObject(const char *name, TClass *cl)
{
  TObject *object = GetFolder()->FindObjectAny(name);
  if(object && object->IsA() != cl)
  {
    cout << "** ERROR: object '" << name;
    cout << "' is not of class '" << cl->GetName() << "'" << endl;
    return 0;
  }
  return object;
}

//------------------------------------------------------------------------------

ExRootResult *ExRootModule::GetPlots()
{
  if(!fPlots)
  {
    fPlots = new ExRootResult();
    fPlots->SetFolder(GetFolder());
  }
  return fPlots;
}

//------------------------------------------------------------------------------

ExRootFactory *ExRootModule::GetFactory()
{
  if(!fFactory)
  {
    fFactory = static_cast<ExRootFactory *>(GetObject("ObjectFactory", ExRootFactory::Class()));
    if(!fFactory)
    {
      cout << "** ERROR: cannot access factory" << endl;
      return 0;
    }
  }
  return fFactory;
}


