
#include <iostream>
#include <sstream>
#include <map>

#include "TROOT.h"
#include "TApplication.h"

#include "TFile.h"
#include "TChain.h"
#include "TString.h"

#include "TH2.h"
#include "THStack.h"
#include "TLegend.h"
#include "TPaveText.h"
#include "TLorentzVector.h"

#include "LHEF.h"

#include "ExRootAnalysis/ExRootClasses.h"

#include "ExRootAnalysis/ExRootTreeWriter.h"
#include "ExRootAnalysis/ExRootTreeBranch.h"

#include "ExRootAnalysis/ExRootUtilities.h"
#include "ExRootAnalysis/ExRootProgressBar.h"

using namespace std;

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

struct LHCOlympicsObject
{
  enum {maxIntParam = 2, maxDblParam = 7};

  Int_t intParam[maxIntParam];
  Double_t dblParam[maxDblParam];
};

//------------------------------------------------------------------------------

class LHCOlympicsConverter
{
public:
  LHCOlympicsConverter(const char *outputFileName);
  ~LHCOlympicsConverter();

  void ProcessObject();
  void Write();

  Long64_t GetNumberOfObjects(ifstream &inputFileStream);
  Bool_t ReadObject(ifstream &inputFileStream);

private:

  void AddMissingEvents();

  void AnalyseEvent(ExRootTreeBranch *branch,
                    Long64_t eventNumber, Int_t triggerWord);

  void AnalysePhoton(ExRootTreeBranch *branch);
  void AnalyseElectron(ExRootTreeBranch *branch);
  void AnalyseMuon(ExRootTreeBranch *branch);
  void AnalyseTau(ExRootTreeBranch *branch);
  void AnalyseJet(ExRootTreeBranch *branch);
  void AnalyseMissingET(ExRootTreeBranch *branch);

  istringstream fBufferStream;
  string fBuffer;

  LHCOlympicsObject fCurrentObject;

  Bool_t fIsFirstEvent, fIsNewEvent, fIsReadyToFill;
  
  Long64_t fPreviousObjectNumber, fTriggerWord, fEventNumber, fRecordNumber;

  TFile *fOutputFile;
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

LHCOlympicsConverter::LHCOlympicsConverter(const char *outputFileName) :
  fIsFirstEvent(kTRUE), fIsNewEvent(kFALSE), fIsReadyToFill(kFALSE),
  fPreviousObjectNumber(0), fTriggerWord(0), fEventNumber(1), fRecordNumber(1),
  fOutputFile(0), fTreeWriter(0)
{
  fOutputFile = TFile::Open(outputFileName, "RECREATE");
  fTreeWriter = new ExRootTreeWriter(fOutputFile, "LHCO");

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

LHCOlympicsConverter::~LHCOlympicsConverter()
{
  if(fTreeWriter) delete fTreeWriter;
  if(fOutputFile) delete fOutputFile;
}

//------------------------------------------------------------------------------

Long64_t LHCOlympicsConverter::GetNumberOfObjects(ifstream &inputFileStream)
{
  Long64_t counter = 0;
  Bool_t canReadNumber, canReadFile = kTRUE;
  Int_t number;
  int position = inputFileStream.tellg();
  inputFileStream.seekg(0, std::ios::beg);

  inputFileStream.clear();

  while(canReadFile)
  {
    do
    {
      getline(inputFileStream, fBuffer);
  
      if(!inputFileStream.good())
      {
        canReadFile = kFALSE;
        break;
      }

      fBufferStream.clear();
      fBufferStream.str(fBuffer);
      
      canReadNumber = (fBufferStream >> number);
    }
    while(!canReadNumber);

    ++counter;
  }

  inputFileStream.clear();

  inputFileStream.seekg(position, std::ios::beg);

  return (counter - 1);
}

//------------------------------------------------------------------------------

Bool_t LHCOlympicsConverter::ReadObject(ifstream &inputFileStream)
{
  Int_t i;
  Bool_t canReadNumber;

  do
  {
    getline(inputFileStream, fBuffer);

    if(!inputFileStream.good()) return kFALSE;

    fBufferStream.clear();
    fBufferStream.str(fBuffer);

    canReadNumber = kTRUE;

    for(i = 0; canReadNumber && i < LHCOlympicsObject::maxIntParam; ++i)
    {
      canReadNumber = (fBufferStream >> fCurrentObject.intParam[i]);
    }
    
    if(canReadNumber && fCurrentObject.intParam[0] == 0)
    {
      fEventNumber = fCurrentObject.intParam[1];
      canReadNumber = (fBufferStream >> fTriggerWord);
    }
    else
    {
      for(i = 0; canReadNumber && i < LHCOlympicsObject::maxDblParam; ++i)
      {
        canReadNumber = (fBufferStream >> fCurrentObject.dblParam[i]);
      }
    }
  }
  while(!canReadNumber);

  return kTRUE;
}

//---------------------------------------------------------------------------

void LHCOlympicsConverter::Write()
{
  if(fIsReadyToFill && fTreeWriter) fTreeWriter->Fill();
  if(fTreeWriter) fTreeWriter->Write();
  fIsReadyToFill = kFALSE;
}

//---------------------------------------------------------------------------
// add empty events for missing event numbers

void LHCOlympicsConverter::AddMissingEvents()
{
  while(fRecordNumber < fEventNumber)
  {
    fTreeWriter->Clear();
    AnalyseEvent(fBranchEvent, fRecordNumber, 0);
    fTreeWriter->Fill();

    ++fRecordNumber;
  }

  fTreeWriter->Clear();

  fIsReadyToFill = kFALSE;
}

//---------------------------------------------------------------------------

void LHCOlympicsConverter::ProcessObject()
{
  fIsNewEvent = (fCurrentObject.intParam[0] <= fPreviousObjectNumber);

  fPreviousObjectNumber = fCurrentObject.intParam[0];

  if(fIsNewEvent && fIsFirstEvent && fTreeWriter)
  {
    fIsFirstEvent = kFALSE;

    AddMissingEvents();
  }

  if(fIsNewEvent && fIsReadyToFill && fTreeWriter)
  {
    fIsReadyToFill = kFALSE;

    fTreeWriter->Fill();
    fTreeWriter->Clear();

    ++fRecordNumber;

    AddMissingEvents();
  }

  if(fCurrentObject.intParam[0] == 0)
  {
    AnalyseEvent(fBranchEvent, fEventNumber, fTriggerWord);
  }
  else
  {
    switch(fCurrentObject.intParam[1])
    {
      case 0: AnalysePhoton(fBranchPhoton); break;
      case 1: AnalyseElectron(fBranchElectron); break;
      case 2: AnalyseMuon(fBranchMuon); break;
      case 3: AnalyseTau(fBranchTau); break;
      case 4: AnalyseJet(fBranchJet); break;
      case 6: AnalyseMissingET(fBranchMissingET); break;
    }
  }

  fIsReadyToFill = kTRUE;
}

//---------------------------------------------------------------------------

void LHCOlympicsConverter::AnalyseEvent(ExRootTreeBranch *branch,
                                        Long64_t eventNumber, Int_t triggerWord)
{
  TRootEvent *element;

  element = static_cast<TRootEvent*>(branch->NewEntry());

  element->Number = eventNumber;
  element->Trigger = triggerWord;
}

//---------------------------------------------------------------------------

void LHCOlympicsConverter::AnalysePhoton(ExRootTreeBranch *branch)
{
  TRootPhoton *element;

  element = static_cast<TRootPhoton*>(branch->NewEntry());

  element->Eta = fCurrentObject.dblParam[0];
  element->Phi = fCurrentObject.dblParam[1];
  element->PT = fCurrentObject.dblParam[2];

  element->EhadOverEem = fCurrentObject.dblParam[6];
}

//---------------------------------------------------------------------------

void LHCOlympicsConverter::AnalyseElectron(ExRootTreeBranch *branch)
{
  TRootElectron *element;

  element = static_cast<TRootElectron*>(branch->NewEntry());

  element->Eta = fCurrentObject.dblParam[0];
  element->Phi = fCurrentObject.dblParam[1];
  element->PT = fCurrentObject.dblParam[2];

  element->Ntrk = TMath::Abs(fCurrentObject.dblParam[4]);

  element->Charge = fCurrentObject.dblParam[4] < 0.0 ? -1.0 : 1.0;

  element->EhadOverEem = fCurrentObject.dblParam[6];
}

//---------------------------------------------------------------------------

void LHCOlympicsConverter::AnalyseMuon(ExRootTreeBranch *branch)
{
  TRootMuon *element;

  element = static_cast<TRootMuon*>(branch->NewEntry());

  element->Eta = fCurrentObject.dblParam[0];
  element->Phi = fCurrentObject.dblParam[1];
  element->PT = fCurrentObject.dblParam[2];

  element->Ntrk = TMath::Abs(fCurrentObject.dblParam[4]);

  element->Charge = fCurrentObject.dblParam[4] < 0.0 ? -1.0 : 1.0;

  element->JetIndex = Int_t(fCurrentObject.dblParam[5]);

  element->PTiso = Int_t(fCurrentObject.dblParam[6]);
  element->ETiso = fCurrentObject.dblParam[6] - element->PTiso;
}

//---------------------------------------------------------------------------

void LHCOlympicsConverter::AnalyseTau(ExRootTreeBranch *branch)
{
  TRootTau *element;

  element = static_cast<TRootTau*>(branch->NewEntry());

  element->Eta = fCurrentObject.dblParam[0];
  element->Phi = fCurrentObject.dblParam[1];
  element->PT = fCurrentObject.dblParam[2];

  element->Ntrk = TMath::Abs(fCurrentObject.dblParam[4]);

  element->Charge = fCurrentObject.dblParam[4] < 0 ? -1.0 : 1.0;

  element->EhadOverEem = fCurrentObject.dblParam[6];
}

//---------------------------------------------------------------------------

void LHCOlympicsConverter::AnalyseJet(ExRootTreeBranch *branch)
{
  TRootJet *element;

  element = static_cast<TRootJet*>(branch->NewEntry());

  element->Eta = fCurrentObject.dblParam[0];
  element->Phi = fCurrentObject.dblParam[1];
  element->PT = fCurrentObject.dblParam[2];

  element->Mass = fCurrentObject.dblParam[3];

  element->Ntrk = TMath::Abs(fCurrentObject.dblParam[4]);

  element->BTag = fCurrentObject.dblParam[5];

  element->EhadOverEem = fCurrentObject.dblParam[6];

  element->Index = fCurrentObject.intParam[0];

}

//---------------------------------------------------------------------------

void LHCOlympicsConverter::AnalyseMissingET(ExRootTreeBranch *branch)
{
  TRootMissingET *element;

  element = static_cast<TRootMissingET*>(branch->NewEntry());

  element->Phi = fCurrentObject.dblParam[1];
  element->MET = fCurrentObject.dblParam[2];
}

//---------------------------------------------------------------------------

int main(int argc, char *argv[])
{
  char appName[] = "ExRootLHEFConverter";

  if(argc != 3)
  {
    cout << " Usage: " << appName << " input_file" << " output_file" << endl;
    cout << " input_file - input file in LHEF format," << endl;
    cout << " output_file - output file in ROOT format." << endl;
    return 1;
  }

  gROOT->SetBatch();

  int appargc = 1;
  char *appargv[] = {appName};
  TApplication app(appName, &appargc, appargv);

  // Open a stream connected to an event file:
  ifstream inputFileStream(argv[1]);

  if(!inputFileStream.is_open())
  {
    cerr << "** ERROR: Can't open '" << argv[1] << "' for input" << endl;
    return 1;
  }

  // Create LHC Olympics converter:
  LHCOlympicsConverter *converter = new LHCOlympicsConverter(argv[2]);

  cout << "** Calculating number of objects to process. Please wait..." << endl;
  Long64_t allEntries = converter->GetNumberOfObjects(inputFileStream);
  cout << "** Input file contains " << allEntries << " objects" << endl;

  if(allEntries > 0)
  {
    ExRootProgressBar progressBar(allEntries);

    // Loop over all objects
    Long64_t entry = 0;
    while(converter->ReadObject(inputFileStream))
    {
      converter->ProcessObject();

      progressBar.Update(entry);
      
      ++entry;
    }
    progressBar.Finish();

    converter->Write();
  }

  cout << "** Exiting..." << endl;

  delete converter;
}


