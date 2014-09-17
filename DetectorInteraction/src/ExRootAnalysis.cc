
/** \class ExRootAnalysis
 *
 *  Analysis steering class.
 *  Implements events loop and modules management.
 *
 *  $Date: 2006/12/13 18:49:30 $
 *  $Revision: 1.3 $
 *
 *
 *  \author P. Demin - UCL, Louvain-la-Neuve
 *
 */

#include "ExRootAnalysis/ExRootAnalysis.h"
#include "ExRootAnalysis/ExRootFactory.h"

#include "ExRootAnalysis/ExRootConfReader.h"
#include "ExRootAnalysis/ExRootTreeReader.h"
#include "ExRootAnalysis/ExRootTreeWriter.h"

#include "ExRootAnalysis/ExRootUtilities.h"
#include "ExRootAnalysis/ExRootProgressBar.h"

#include "TROOT.h"
#include "TClass.h"
#include "TSystem.h"
#include "TFolder.h"
#include "TObjArray.h"

#include <iostream>

#include <string.h>
#include <stdio.h>

using namespace std;

ExRootAnalysis::ExRootAnalysis() :
  fTreeFile(0), fInfoFile(0)
{
  TFolder *folder = new TFolder("", "");
  SetFolder(folder);

  fChains = new TObjArray;
  fChains->SetOwner();

  fConfReader = new ExRootConfReader;

  fTreeReader = new ExRootTreeReader();

  fTreeWriter = new ExRootTreeWriter();

  fFactory = new ExRootFactory();

}

//------------------------------------------------------------------------------

ExRootAnalysis::~ExRootAnalysis()
{
  delete fFactory;
  delete fTreeWriter;
  delete fTreeReader;
  delete fConfReader;
  delete fChains;
  delete GetFolder();
}

//------------------------------------------------------------------------------

Long64_t ExRootAnalysis::GetEntries() const
{
  return fTreeReader ? fTreeReader->GetEntries() : 0;
}

//------------------------------------------------------------------------------

Bool_t ExRootAnalysis::ReadEvent(Long64_t entry)
{
  return fTreeReader ? fTreeReader->ReadEntry(entry) : kFALSE;
}

//------------------------------------------------------------------------------

void ExRootAnalysis::ProcessEvent()
{
  Clear();
  ExRootTask::ProcessEvent();
  if(fTreeWriter) fTreeWriter->Fill();
}

//------------------------------------------------------------------------------

void ExRootAnalysis::Loop()
{
  cout << "** Calculating number of events to process. Please wait..." << endl;
  Long64_t allEntries = GetEntries();
  cout << "** Chain contains " << allEntries << " events" << endl;
  Long64_t entry;

  if(allEntries > 0)
  {
    ExRootProgressBar progressBar(allEntries);
    // Loop over all events
    for(entry = 0; entry < allEntries; ++entry)
    {
      if(!ReadEvent(entry))
      {
        cout << "** ERROR: cannot read event " << entry << endl;
        break;
      }
  
      ProcessEvent();

      progressBar.Update(entry);
    }
    progressBar.Finish();
  }
}

//------------------------------------------------------------------------------

void ExRootAnalysis::Init()
{
  fConfReader->ReadFile(fTclFileName);

  TString name = fConfReader->GetString("AppName", "ExRootAnalysis");

  TFolder *folder = GetFolder();
  folder->SetName(name);
  gROOT->GetListOfBrowsables()->Add(folder);

  SetName(name);
  folder->Add(this);

  fConfReader->SetName("ConfReader");
  folder->Add(fConfReader);

  ExRootConfParam param = fConfReader->GetParam("::InputCollection");
  Long_t i, size;
  TChain *chain = 0, *firstChain = 0;
  size = param.GetSize();
  if(size > 0)
  {
    for(i = 0; i < size; ++i)
    {
      chain = new TChain("", "");
      fChains->Add(chain);
      name = param[i][0].GetString();
      chain->SetName(name);
      FillChain(chain, param[i][1].GetString());
      if(i == 0)
      {
        firstChain = chain;
      }
      else
      {
        firstChain->AddFriend(chain, name + i);
      }
    }
    fTreeReader->SetTree(firstChain);
  }
  fTreeReader->SetName("TreeReader");
  folder->Add(fTreeReader);

  name = fConfReader->GetString("OutputFile", "Analysis");
  name.ReplaceAll(".root", "");
  fTreeFile = TFile::Open(name + "Tree.root", "RECREATE");
  if(!fTreeFile)
  {
    cout << "** ERROR: cannot create output tree file" << endl;
    return;
  }

  fInfoFile = TFile::Open(name + "Info.root", "RECREATE");
  if(!fInfoFile)
  {
    cout << "** ERROR: cannot create output info file" << endl;
    return;
  }
  
  name = fConfReader->GetString("TreeName", "Analysis");

  fTreeWriter->SetTreeFile(fTreeFile);
  fTreeWriter->SetTreeName(name);
  fTreeWriter->SetName("TreeWriter");
  folder->Add(fTreeWriter);

  fFactory->SetName("ObjectFactory");
  folder->Add(fFactory);

  ExRootTask *task;
  const ExRootConfReader::ExRootTaskMap *modules = fConfReader->GetModules();
  ExRootConfReader::ExRootTaskMap::const_iterator itModules;

  param = fConfReader->GetParam("::ExecutionPath");
  size = param.GetSize();

  for(i = 0; i < size; ++i)
  {
    name = param[i].GetString();
    itModules = modules->find(name);
    if(itModules != modules->end())
    {
      cout << itModules->second << " \t " <<  itModules->first << endl;
      task = NewTask(itModules->second, itModules->first);
      if(task)
      {
        task->SetFolder(GetFolder());
        task->SetConfReader(fConfReader);
      }      
    }
    else
    {
      cout << "** ERROR: module '" << name;
      cout << "' is specified in ExecutionPath but not configured.";
      return;
    }

  }

}

//------------------------------------------------------------------------------

void ExRootAnalysis::Event()
{
}

//------------------------------------------------------------------------------

void ExRootAnalysis::Finish()
{
  if(fTreeWriter) fTreeWriter->Write();
}

//------------------------------------------------------------------------------

void ExRootAnalysis::Clear()
{
  if(fTreeWriter) fTreeWriter->Clear();
  if(fFactory) fFactory->Clear();
}

//------------------------------------------------------------------------------

