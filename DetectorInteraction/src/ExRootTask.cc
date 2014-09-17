
/** \class ExRootTask
 *
 *  Class handling output ROOT tree
 *
 *  $Date: 2006/09/21 13:10:52 $
 *  $Revision: 1.1 $
 *
 *
 *  \author P. Demin - UCL, Louvain-la-Neuve
 *
 */

#include "ExRootAnalysis/ExRootTask.h"
#include "ExRootAnalysis/ExRootConfReader.h"

#include "TROOT.h"
#include "TClass.h"
#include "TString.h"

#include <iostream>

const char *const kINIT = "0";
const char *const kPROCESS = "1";
const char *const kFINISH = "2";

using namespace std;

ExRootTask::ExRootTask() :
  TTask("", ""), fFolder(0), fConfReader(0)
{
}

//------------------------------------------------------------------------------

ExRootTask::~ExRootTask()
{
}

//------------------------------------------------------------------------------

void ExRootTask::Exec(Option_t *option)
{
  if(option == kINIT)
  {
    Init();
  }
  else if(option == kPROCESS)
  {
    Event();
  }
  else if(option == kFINISH)
  {
    Finish();
  }
}

//------------------------------------------------------------------------------

void ExRootTask::InitTask()
{
  ExecuteTask(kINIT);
}

//------------------------------------------------------------------------------

void ExRootTask::ProcessEvent()
{
  ExecuteTask(kPROCESS);
}

//------------------------------------------------------------------------------

void ExRootTask::FinishTask()
{
  ExecuteTask(kFINISH);
}

//------------------------------------------------------------------------------

void ExRootTask::Add(TTask *task)
{
  if(!task) return;

  if(!task->IsA()->InheritsFrom(ExRootTask::Class()))
  {
    cout << "** ERROR: task '" << task->IsA()->GetName();
    cout << "' does not inherit from ExRootTask" << endl;
    return;
  }

  TTask::Add(task);
}

//------------------------------------------------------------------------------

ExRootTask *ExRootTask::NewTask(TClass *cl, const char *taskName)
{
  if(!cl) return 0;
 
  if(!cl->InheritsFrom(ExRootTask::Class()))
  {
    cout << "** ERROR: task '" << cl->GetName();
    cout << "' does not inherit from ExRootTask" << endl;
    return 0;
  }

  ExRootTask *task = static_cast<ExRootTask *>(cl->New());
  task->SetName(taskName);
  task->SetFolder(fFolder);

  Add(task);

  return task;
}

//------------------------------------------------------------------------------

ExRootTask *ExRootTask::NewTask(const char *className, const char *taskName)
{
  TClass *cl = gROOT->GetClass(className);
  if(!cl)
  {
    cout << "** ERROR: cannot find class '" << className << "'" << endl;
    return 0;
  }

  return NewTask(cl, taskName);
}

//------------------------------------------------------------------------------

ExRootConfParam ExRootTask::GetParam(const char *name)
{
  if(fConfReader)
  {
    return fConfReader->GetParam(TString(GetName()) + "::" + name);
  }
  else
  {
    return ExRootConfParam(TString(GetName()) + "::" + name, 0, 0);
  }
}

//------------------------------------------------------------------------------

int ExRootTask::GetInt(const char *name, int defaultValue, int index)
{
  if(fConfReader)
  {
    return fConfReader->GetInt(TString(GetName()) + "::" + name, defaultValue, index);
  }
  else
  {
    return defaultValue;
  }
}

//------------------------------------------------------------------------------

long ExRootTask::GetLong(const char *name, long defaultValue, int index)
{
  if(fConfReader)
  {
    return fConfReader->GetLong(TString(GetName()) + "::" + name, defaultValue, index);
  }
  else
  {
    return defaultValue;
  }
}

//------------------------------------------------------------------------------

double ExRootTask::GetDouble(const char *name, double defaultValue, int index)
{
  if(fConfReader)
  {
    return fConfReader->GetDouble(TString(GetName()) + "::" + name, defaultValue, index);
  }
  else
  {
    return defaultValue;
  }
}

//------------------------------------------------------------------------------

bool ExRootTask::GetBool(const char *name, bool defaultValue, int index)
{
  if(fConfReader)
  {
    return fConfReader->GetBool(TString(GetName()) + "::" + name, defaultValue, index);
  }
  else
  {
    return defaultValue;
  }
}

//------------------------------------------------------------------------------

const char *ExRootTask::GetString(const char *name, const char *defaultValue, int index)
{
  if(fConfReader)
  {
    return fConfReader->GetString(TString(GetName()) + "::" + name, defaultValue, index);
  }
  else
  {
    return defaultValue;
  }
}


