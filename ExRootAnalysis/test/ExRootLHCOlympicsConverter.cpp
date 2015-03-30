#include <stdexcept>
#include <iostream>
#include <sstream>

#include <stdlib.h>
#include <signal.h>
#include <stdio.h>

#include "TROOT.h"
#include "TApplication.h"

#include "TFile.h"
#include "TLorentzVector.h"

#include "ExRootAnalysis/ExRootStream.h"
#include "ExRootAnalysis/ExRootClasses.h"

#include "ExRootAnalysis/ExRootTreeWriter.h"
#include "ExRootAnalysis/ExRootTreeBranch.h"
#include "ExRootAnalysis/ExRootProgressBar.h"

using namespace std;

static const int kBufferSize  = 1024;

/*
LHC Olympics format discription from http://www.jthaler.net/olympicswiki/doku.php?id=lhc_olympics:data_file_format

    * The first column of each row is just a counter that labels the object.
    * The event begins with a row labelled "0"; this row contains the event number and the triggering information. The last row of the event is always the missing transverse momentum (MET).
    * The second column of each row gives the type of object being listed [0, 1, 2, 3, 4, 6 = photon, electron, muon, hadronically-decaying tau, jet, missing transverse energy].
    * The next three columns give the pseudorapidity, the azimuthal angle, and the transverse momentum of the object.
    * The sixth column gives the invariant mass of the object.
    * The seventh column gives the number of tracks associated with the object; in the case of a lepton, this number is multiplied by the charge of the lepton.
    * The eighth column is 1 or 2 for a jet that has been "tagged" as containing a b-quark (actually a heavy flavor tag that sometimes indicates c-quarks), otherwise it is 0. For muons, the integer part of this number is the identity of the jet (see column 1) that is closest ot this muon in Delta R.
    * The ninth column is the ratio of the hadronic versus electromagnetic energy deposited in the calorimeter cells associated with the object. For muons to the left of the decimal point is the summed pT in a R=0.4 cone (excluding the muon). To the right of the decimal point is etrat, which is a percentage between .00 and .99. It is the ratio of the transverse energy in a 3x3 grid surrounding the muon to the pT of the muon.
*/

//------------------------------------------------------------------------------

class LHCOConverter
{
public:
  LHCOConverter(TFile *outputFile);
  ~LHCOConverter();

  void Write();

  Bool_t ReadLine(FILE *inputFile);

private:

  void AddMissingEvents();

  void AnalyseEvent(ExRootTreeBranch *branch);

  void AnalysePhoton(ExRootTreeBranch *branch);
  void AnalyseElectron(ExRootTreeBranch *branch);
  void AnalyseMuon(ExRootTreeBranch *branch);
  void AnalyseTau(ExRootTreeBranch *branch);
  void AnalyseJet(ExRootTreeBranch *branch);
  void AnalyseMissingET(ExRootTreeBranch *branch);

  enum {kIntParamSize = 2, kDblParamSize = 7};
  Int_t fIntParam[kIntParamSize];
  Double_t fDblParam[kDblParamSize];

  Bool_t fIsReadyToFill;

  Int_t fTriggerWord, fEventNumber;

  char *fBuffer;

  ExRootTreeWriter *fTreeWriter;

  ExRootTreeBranch *fBranchEvent;
  ExRootTreeBranch *fBranchPhoton;
  ExRootTreeBranch *fBranchElectron;
  ExRootTreeBranch *fBranchMuon;
  ExRootTreeBranch *fBranchTau;
  ExRootTreeBranch *fBranchJet;
  ExRootTreeBranch *fBranchMissingET;
};

//------------------------------------------------------------------------------

LHCOConverter::LHCOConverter(TFile *outputFile) :
  fIsReadyToFill(kFALSE),
  fTriggerWord(0), fEventNumber(1),
  fBuffer(0), fTreeWriter(0)
{
  fBuffer = new char[kBufferSize];
  fTreeWriter = new ExRootTreeWriter(outputFile, "LHCO");

  // information about reconstructed event
  fBranchEvent = fTreeWriter->NewBranch("Event", TRootEvent::Class());
  // reconstructed photons
  fBranchPhoton = fTreeWriter->NewBranch("Photon", TRootPhoton::Class());
  // reconstructed electrons
  fBranchElectron = fTreeWriter->NewBranch("Electron", TRootElectron::Class());
  // reconstructed muons
  fBranchMuon = fTreeWriter->NewBranch("Muon", TRootMuon::Class());
  // reconstructed hadronically-decaying tau leptons
  fBranchTau = fTreeWriter->NewBranch("Tau", TRootTau::Class());
  // reconstructed jets
  fBranchJet = fTreeWriter->NewBranch("Jet", TRootJet::Class());
  // missing transverse energy
  fBranchMissingET = fTreeWriter->NewBranch("MissingET", TRootMissingET::Class());
}

//------------------------------------------------------------------------------

LHCOConverter::~LHCOConverter()
{
  if(fTreeWriter) delete fTreeWriter;
  if(fBuffer) delete[] fBuffer;
}

//------------------------------------------------------------------------------

Bool_t LHCOConverter::ReadLine(FILE *inputFile)
{
  int rc;

  if(!fgets(fBuffer, kBufferSize, inputFile)) return kFALSE;

  ExRootStream bufferStream(fBuffer);

  rc = bufferStream.ReadInt(fIntParam[0]);

  if(!rc)
  {
    return kTRUE;
  }

  if(fIntParam[0] == 0)
  {
    rc = bufferStream.ReadInt(fEventNumber)
      && bufferStream.ReadInt(fTriggerWord);

    if(!rc)
    {
      cerr << "** ERROR: " << "invalid event format" << endl;
      return kFALSE;
    }

    if(fIsReadyToFill && fTreeWriter)
    {
      fTreeWriter->Fill();
      fTreeWriter->Clear();
    }

    AnalyseEvent(fBranchEvent);
    fIsReadyToFill = kTRUE;
  }
  else
  {
    rc = bufferStream.ReadInt(fIntParam[1])
      && bufferStream.ReadDbl(fDblParam[0])
      && bufferStream.ReadDbl(fDblParam[1])
      && bufferStream.ReadDbl(fDblParam[2])
      && bufferStream.ReadDbl(fDblParam[3])
      && bufferStream.ReadDbl(fDblParam[4])
      && bufferStream.ReadDbl(fDblParam[5])
      && bufferStream.ReadDbl(fDblParam[6]);

    if(!rc)
    {
      cerr << "** ERROR: " << "invalid object format" << endl;
      return kFALSE;
    }

    switch(fIntParam[1])
    {
      case 0: AnalysePhoton(fBranchPhoton); break;
      case 1: AnalyseElectron(fBranchElectron); break;
      case 2: AnalyseMuon(fBranchMuon); break;
      case 3: AnalyseTau(fBranchTau); break;
      case 4: AnalyseJet(fBranchJet); break;
      case 6: AnalyseMissingET(fBranchMissingET); break;
    }
  }

  return kTRUE;
}

//---------------------------------------------------------------------------

void LHCOConverter::Write()
{
  if(fIsReadyToFill && fTreeWriter) fTreeWriter->Fill();
  if(fTreeWriter) fTreeWriter->Write();
  fIsReadyToFill = kFALSE;
}

//---------------------------------------------------------------------------

void LHCOConverter::AnalyseEvent(ExRootTreeBranch *branch)
{
  TRootEvent *element;

  element = static_cast<TRootEvent*>(branch->NewEntry());

  element->Number = fEventNumber;
  element->Trigger = fTriggerWord;
}

//---------------------------------------------------------------------------

void LHCOConverter::AnalysePhoton(ExRootTreeBranch *branch)
{
  TRootPhoton *element;

  element = static_cast<TRootPhoton*>(branch->NewEntry());

  element->Eta = fDblParam[0];
  element->Phi = fDblParam[1];
  element->PT = fDblParam[2];
  element->EhadOverEem = fDblParam[6];
}

//---------------------------------------------------------------------------

void LHCOConverter::AnalyseElectron(ExRootTreeBranch *branch)
{
  TRootElectron *element;

  element = static_cast<TRootElectron*>(branch->NewEntry());

  element->Eta = fDblParam[0];
  element->Phi = fDblParam[1];
  element->PT = fDblParam[2];

  element->Charge = fDblParam[4] < 0.0 ? -1 : 1;

  element->Ntrk = TMath::Abs(fDblParam[4]);

  element->EhadOverEem = fDblParam[6];
}

//---------------------------------------------------------------------------

void LHCOConverter::AnalyseMuon(ExRootTreeBranch *branch)
{
  TRootMuon *element;

  element = static_cast<TRootMuon*>(branch->NewEntry());

  element->Eta = fDblParam[0];
  element->Phi = fDblParam[1];
  element->PT = fDblParam[2];

  element->Charge = fDblParam[4] < 0.0 ? -1 : 1;

  element->Ntrk = TMath::Abs(fDblParam[4]);

  element->JetIndex = Int_t(fDblParam[5]);

  element->PTiso = Int_t(fDblParam[6]);
  element->ETiso = fDblParam[6] - element->PTiso;
}

//---------------------------------------------------------------------------

void LHCOConverter::AnalyseTau(ExRootTreeBranch *branch)
{
  TRootTau *element;

  element = static_cast<TRootTau*>(branch->NewEntry());

  element->Eta = fDblParam[0];
  element->Phi = fDblParam[1];
  element->PT = fDblParam[2];

  element->Charge = fDblParam[4] < 0 ? -1 : 1;

  element->Ntrk = TMath::Abs(fDblParam[4]);

  element->EhadOverEem = fDblParam[6];
}

//---------------------------------------------------------------------------

void LHCOConverter::AnalyseJet(ExRootTreeBranch *branch)
{
  TRootJet *element;

  element = static_cast<TRootJet*>(branch->NewEntry());

  element->Eta = fDblParam[0];
  element->Phi = fDblParam[1];
  element->PT = fDblParam[2];

  element->Mass = fDblParam[3];

  element->Ntrk = TMath::Abs(Int_t(fDblParam[4]));

  element->BTag = Int_t(fDblParam[5]);

  element->EhadOverEem = fDblParam[6];

  element->Index = fIntParam[0];
}

//---------------------------------------------------------------------------

void LHCOConverter::AnalyseMissingET(ExRootTreeBranch *branch)
{
  TRootMissingET *element;

  element = static_cast<TRootMissingET*>(branch->NewEntry());

  element->Phi = fDblParam[1];
  element->MET = fDblParam[2];
}

//---------------------------------------------------------------------------

static bool interrupted = false;

void SignalHandler(int sig)
{
  interrupted = true;
}

//---------------------------------------------------------------------------

int main(int argc, char *argv[])
{
  char appName[] = "ExRootLHCOConverter";
  stringstream message;
  FILE *inputFile = 0;
  TFile *outputFile = 0;
  LHCOConverter *converter = 0;
  Long64_t length, eventCounter;

  if(argc != 3)
  {
    cout << " Usage: " << appName << " input_file" << " output_file" << endl;
    cout << " input_file - input file in LHEF format," << endl;
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
      message << "can't open " << argv[2];
      throw runtime_error(message.str());
    }

    converter = new LHCOConverter(outputFile);

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
      eventCounter = 0;
      ExRootProgressBar progressBar(length);

      // Loop over all objects
      while(converter->ReadLine(inputFile) && !interrupted)
      {
        ++eventCounter;

        progressBar.Update(ftello(inputFile), eventCounter);
      }
      converter->Write();

      fseek(inputFile, 0L, SEEK_END);
      progressBar.Update(ftello(inputFile), eventCounter, kTRUE);
      progressBar.Finish();
    }

    fclose(inputFile);

    cout << "** Exiting..." << endl;

    delete converter;
    delete outputFile;

    return 0;
  }
  catch(runtime_error &e)
  {
    if(converter) delete converter;
    if(outputFile) delete outputFile;
    cerr << "** ERROR: " << e.what() << endl;
    return 1;
  }
}


