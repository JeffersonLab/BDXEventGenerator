#include <iostream>
#include <utility>
#include <vector>

#include "TROOT.h"
#include "TApplication.h"

#include "TChain.h"
#include "TString.h"

#include "TClonesArray.h"

#include "TH2.h"
#include "THStack.h"
#include "TLegend.h"
#include "TPaveText.h"
#include "TLorentzVector.h"

#include "ExRootAnalysis/ExRootClasses.h"

#include "ExRootAnalysis/ExRootTreeReader.h"
#include "ExRootAnalysis/ExRootTreeWriter.h"
#include "ExRootAnalysis/ExRootTreeBranch.h"
#include "ExRootAnalysis/ExRootResult.h"
#include "ExRootAnalysis/ExRootUtilities.h"

using namespace std;

//------------------------------------------------------------------------------

// Here you can put your analysis macro

#include "Example.C"

//------------------------------------------------------------------------------

int main(int argc, char *argv[])
{
  char appName[] = "JetsSim";

  if(argc != 2)
  {
    cout << " Usage: " << appName << " input_file_list" << " output_file" << endl;
    cout << " input_file_list - list of input files in ROOT format ('Analysis' tree)," << endl;
    return 1;
  }

  gROOT->SetBatch();

  int appargc = 1;
  char *appargv[] = {appName};
  TApplication app(appName, &appargc, appargv);

  TString inputFileList(argv[1]);

//------------------------------------------------------------------------------

// Here you call your macro's main function 

  Example(inputFileList);

//------------------------------------------------------------------------------

}

