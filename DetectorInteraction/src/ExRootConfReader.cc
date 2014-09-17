
/** \class ExRootConfReader
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

#include "ExRootAnalysis/ExRootConfReader.h"

#include "tcl/tcl.h"

#include <iostream>
#include <fstream>
#include <string>

using namespace std;

static Tcl_ObjCmdProc ModuleObjCmdProc;

//------------------------------------------------------------------------------

ExRootConfReader::ExRootConfReader() :
  fTclInterp(0)
{
  fTclInterp = Tcl_CreateInterp();

  Tcl_CreateObjCommand(fTclInterp, "module", ModuleObjCmdProc, this, 0);
}

//------------------------------------------------------------------------------

ExRootConfReader::~ExRootConfReader()
{
  Tcl_DeleteInterp(fTclInterp);
}

//------------------------------------------------------------------------------

void ExRootConfReader::ReadFile(const char *fileName)
{
/*
  ifstream infile(fileName);
  string cmdBuffer = string(istreambuf_iterator<char>(infile), istreambuf_iterator<char>());

  Tcl_Obj *cmdObjPtr = Tcl_NewObj();
  cmdObjPtr->bytes = const_cast<char *>(cmdBuffer.c_str());
  cmdObjPtr->length = cmdBuffer.size();
*/

  ifstream infile(fileName, ios::in | ios::ate);
  int file_length = infile.tellg();
  infile.seekg(0, ios::beg);
  infile.clear();
  char *cmdBuffer = new char[file_length];
  infile.read(cmdBuffer, file_length);

  Tcl_Obj *cmdObjPtr = Tcl_NewObj();
  cmdObjPtr->bytes = cmdBuffer;
  cmdObjPtr->length = file_length;

  Tcl_IncrRefCount(cmdObjPtr);

  if(Tcl_EvalObj(fTclInterp, cmdObjPtr) != TCL_OK)
  {
    cerr << endl;
    cerr << "** ERROR: " << endl;
    cerr << Tcl_GetStringResult(fTclInterp) << endl;
  }

  cmdObjPtr->bytes = 0;
  cmdObjPtr->length = 0;

  Tcl_DecrRefCount(cmdObjPtr);

  delete[] cmdBuffer;
}

//------------------------------------------------------------------------------

ExRootConfParam ExRootConfReader::GetParam(const char *name)
{
  Tcl_Obj *object;
  Tcl_Obj *variableName = Tcl_NewStringObj(const_cast<char *>(name),-1);
  object = Tcl_ObjGetVar2(fTclInterp, variableName, 0, TCL_GLOBAL_ONLY);
  return ExRootConfParam(name, object, fTclInterp);
}

//------------------------------------------------------------------------------

int ExRootConfReader::GetInt(const char *name, int defaultValue, int index)
{
  ExRootConfParam object = GetParam(name);
  if(index >= 0)
  {
    object = object[index];
  }

  return object.GetInt(defaultValue);
}

//------------------------------------------------------------------------------

long ExRootConfReader::GetLong(const char *name, long defaultValue, int index)
{
  ExRootConfParam object = GetParam(name);
  if(index >= 0)
  {
    object = object[index];
  }

  return object.GetLong(defaultValue);
}

//------------------------------------------------------------------------------

double ExRootConfReader::GetDouble(const char *name, double defaultValue, int index)
{
  ExRootConfParam object = GetParam(name);
  if(index >= 0)
  {
    object = object[index];
  }

  return object.GetDouble(defaultValue);
}

//------------------------------------------------------------------------------

bool ExRootConfReader::GetBool(const char *name, bool defaultValue, int index)
{
  ExRootConfParam object = GetParam(name);
  if(index >= 0)
  {
    object = object[index];
  }

  return object.GetBool(defaultValue);
}

//------------------------------------------------------------------------------

const char *ExRootConfReader::GetString(const char *name, const char *defaultValue, int index)
{
  ExRootConfParam object = GetParam(name);
  if(index >= 0)
  {
    object = object[index];
  }

  return object.GetString(defaultValue);
}

//------------------------------------------------------------------------------

void ExRootConfReader::AddModule(const char *className, const char *moduleName)
{
  ExRootTaskMap::iterator itMoudles = fModules.find(moduleName);

  if(itMoudles != fModules.end())
  {
    cout << "** WARNING: module '" << moduleName << "' is already configured.";
    cout << " Only first entry will be used." << endl;
  }
  else
  {
    fModules.insert(make_pair(moduleName, className));
    cout << "** INFO: adding module ";
    cout << className << " \t ";
    cout << moduleName << endl;
  }
}

//------------------------------------------------------------------------------

int ModuleObjCmdProc(ClientData clientData, Tcl_Interp *interp, int objc, Tcl_Obj *const objv[])
{
  if(objc < 3)
  {
    Tcl_SetResult(interp, "wrong # args: should be \"module className moduleName arg ?arg...?\"", 0);
    return TCL_ERROR;
  }

  ExRootConfReader *test = (ExRootConfReader*) clientData;

  // add module to a list of modules to be created

  test->AddModule(Tcl_GetStringFromObj(objv[1], 0), Tcl_GetStringFromObj(objv[2], 0));

  if(objc > 3)
  {
    Tcl_Obj *object = Tcl_NewListObj(0, 0);
    Tcl_ListObjAppendElement(interp, object, Tcl_NewStringObj("namespace", -1));
    Tcl_ListObjAppendElement(interp, object, Tcl_NewStringObj("eval", -1));
    Tcl_ListObjAppendList(interp, object, Tcl_NewListObj(objc-2, objv+2));

    return Tcl_GlobalEvalObj(interp, object);
  }

  return TCL_OK;
}

//------------------------------------------------------------------------------

ExRootConfParam::ExRootConfParam(const char *name, Tcl_Obj *object, Tcl_Interp *interp) :
  fName(name), fObject(object), fTclInterp(interp)
{
}

//------------------------------------------------------------------------------

int ExRootConfParam::GetInt(int defaultValue)
{
  int result = defaultValue;
  if(fObject && TCL_OK != Tcl_GetIntFromObj(fTclInterp, fObject, &result))
  {
    cerr << endl;
    cerr << "** ERROR: parameter '"<< fName << "' is not an integer." << endl;
    cerr << "** ERROR: " << fName << " = " << Tcl_GetStringFromObj(fObject, 0) << endl;
  }
  return result;
}

//------------------------------------------------------------------------------

long ExRootConfParam::GetLong(long defaultValue)
{
  long result = defaultValue;
  if(fObject && TCL_OK != Tcl_GetLongFromObj(fTclInterp, fObject, &result))
  {
    cerr << endl;
    cerr << "** ERROR: parameter '"<< fName << "' is not a long integer." << endl;
    cerr << "** ERROR: " << fName << " = " << Tcl_GetStringFromObj(fObject, 0) << endl;
  }
  return result;
}

//------------------------------------------------------------------------------

double ExRootConfParam::GetDouble(double defaultValue)
{
  double result = defaultValue;
  if(fObject && TCL_OK != Tcl_GetDoubleFromObj(fTclInterp, fObject, &result))
  {
    cerr << endl;
    cerr << "** ERROR: parameter '"<< fName << "' is not a number." << endl;
    cerr << "** ERROR: " << fName << " = " << Tcl_GetStringFromObj(fObject, 0) << endl;
  }
  return result;
}

//------------------------------------------------------------------------------

bool ExRootConfParam::GetBool(bool defaultValue)
{
  int result = defaultValue;
  if(fObject && TCL_OK != Tcl_GetBooleanFromObj(fTclInterp, fObject, &result))
  {
    cerr << endl;
    cerr << "** ERROR: parameter '"<< fName << "' is not a boolean." << endl;
    cerr << "** ERROR: " << fName << " = " << Tcl_GetStringFromObj(fObject, 0) << endl;
  }
  return result;
}

//------------------------------------------------------------------------------

const char *ExRootConfParam::GetString(const char *defaultValue)
{
  const char *result = defaultValue;
  if(fObject) result = Tcl_GetStringFromObj(fObject, 0);
  return result;  
}

//------------------------------------------------------------------------------

int ExRootConfParam::GetSize()
{
  int length = 0;
  if(fObject && TCL_OK != Tcl_ListObjLength(fTclInterp, fObject, &length))
  {
    cerr << endl;
    cerr << "** ERROR: parameter '"<< fName << "' is not a list." << endl;
    cerr << "** ERROR: " << fName << " = " << Tcl_GetStringFromObj(fObject, 0) << endl;
  }
  return length;  
}

//------------------------------------------------------------------------------

ExRootConfParam ExRootConfParam::operator[](int index)
{
  Tcl_Obj *object = 0;
  if(fObject && TCL_OK != Tcl_ListObjIndex(fTclInterp, fObject, index, &object))
  {
    cerr << endl;
    cerr << "** ERROR: parameter '"<< fName << "' is not a list." << endl;
    cerr << "** ERROR: " << fName << " = " << Tcl_GetStringFromObj(fObject, 0) << endl;
  }
  return ExRootConfParam(fName, object, fTclInterp);
}


