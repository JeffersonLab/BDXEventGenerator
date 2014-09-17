
#include "ExRootAnalysis/ExRootAnalysis.h"

#include "TROOT.h"
#include "TApplication.h"

#include <iostream>

using namespace std;

int main(int argc, char *argv[])
{
  char *appName = "ExRootMain";

  if(argc != 2)
  {
    cout << " Usage: " << appName << " input_file" << endl;
    cout << " input_file - configuration file in Tcl format." << endl;
    return 1;
  }

  gROOT->SetBatch();

  int appargc = 1;
  char *appargv[] = {appName};
  TApplication app(appName, &appargc, appargv);
  
  TString inputFile(argv[1]);
  
  ExRootAnalysis test;
  test.SetTclFileName(inputFile);
  test.InitTask();
  test.Loop();
  test.FinishTask();
}
