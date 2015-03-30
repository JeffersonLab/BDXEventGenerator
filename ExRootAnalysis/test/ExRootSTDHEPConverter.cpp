#include <stdexcept>
#include <iostream>
#include <sstream>

#include <signal.h>

#include "TROOT.h"
#include "TApplication.h"

#include "TFile.h"
#include "TLorentzVector.h"

#include "ExRootAnalysis/ExRootClasses.h"
#include "ExRootAnalysis/ExRootSTDHEPReader.h"

#include "ExRootAnalysis/ExRootTreeWriter.h"
#include "ExRootAnalysis/ExRootTreeBranch.h"
#include "ExRootAnalysis/ExRootProgressBar.h"

using namespace std;

//---------------------------------------------------------------------------

static bool interrupted = false;

void SignalHandler(int sig)
{
  interrupted = true;
}

//---------------------------------------------------------------------------

int main(int argc, char *argv[])
{
  char appName[] = "ExRootSTDHEPConverter";
  stringstream message;
  FILE *inputFile = 0;
  TFile *outputFile = 0;
  ExRootTreeWriter *treeWriter = 0;
  ExRootTreeBranch *branchGenEvent = 0, *branchGenParticle = 0;
  ExRootSTDHEPReader *reader = 0;
  Long64_t length, eventCounter;

  if(argc != 3)
  {
    cout << " Usage: " << appName << " input_file" << " output_file" << endl;
    cout << " input_file - input file in STDHEP format," << endl;
    cout << " output_file - output file in ROOT format." << endl;
    return 1;
  }

  signal(SIGINT, SignalHandler);

  gROOT->SetBatch();

  int appargc = 1;
  char *appargv[] = {appName};
  TApplication app(appName, &appargc, appargv);

  try
  {
    outputFile = TFile::Open(argv[2], "CREATE");

    if(outputFile == NULL)
    {
      message << "can't create output file " << argv[2];
      throw runtime_error(message.str());
    }

    treeWriter = new ExRootTreeWriter(outputFile, "STDHEP");

    // information about generated event
    branchGenEvent = treeWriter->NewBranch("Event", TRootLHEFEvent::Class());
    // generated particles from HEPEVT
    branchGenParticle = treeWriter->NewBranch("GenParticle", TRootGenParticle::Class());

    reader = new ExRootSTDHEPReader;

    cout << "** Reading " << argv[1] << endl;
    inputFile = fopen(argv[1], "r");

    if(inputFile == NULL)
    {
      message << "can't open " << argv[1];
      throw runtime_error(message.str());
    }

    fseek(inputFile, 0L, SEEK_END);
    length = ftello(inputFile);
    fseek(inputFile, 0L, SEEK_SET);

    if(length > 0)
    {
      reader->SetInputFile(inputFile);

      ExRootProgressBar progressBar(length);

      // Loop over all objects
      eventCounter = 0;
      treeWriter->Clear();
      reader->Clear();
      while(reader->ReadBlock(branchGenParticle) && !interrupted)
      {
        if(reader->EventReady())
        {
          ++eventCounter;

          reader->AnalyzeEvent(branchGenEvent, eventCounter);

          treeWriter->Fill();

          treeWriter->Clear();

          reader->Clear();
        }
        progressBar.Update(ftello(inputFile), eventCounter);
      }

      fseek(inputFile, 0L, SEEK_END);
      progressBar.Update(ftello(inputFile), eventCounter, kTRUE);
      progressBar.Finish();
    }

    fclose(inputFile);

    treeWriter->Write();

    cout << "** Exiting..." << endl;

    delete reader;
    delete treeWriter;
    delete outputFile;

    return 0;
  }
  catch(runtime_error &e)
  {
    if(treeWriter) delete treeWriter;
    if(outputFile) delete outputFile;
    cerr << "** ERROR: " << e.what() << endl;
    return 1;
  }
}
