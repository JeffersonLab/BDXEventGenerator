
#include <iostream>
#include <utility>
#include <deque>

#include "TROOT.h"
#include "TApplication.h"

#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TLeaf.h"
#include "TString.h"
#include "TLorentzVector.h"

#include "ExRootAnalysis/ExRootClasses.h"

#include "ExRootAnalysis/ExRootTreeReader.h"
#include "ExRootAnalysis/ExRootTreeWriter.h"
#include "ExRootAnalysis/ExRootTreeBranch.h"
#include "ExRootAnalysis/ExRootUtilities.h"

using namespace std;

//------------------------------------------------------------------------------

struct HEPEvent
{
  Int_t           Nevhep;
  Int_t           Nhep;
  Int_t           *Idhep;    //[Nhep]
  Int_t           *Jsmhep;   //[Nhep]
  Int_t           *Jsdhep;   //[Nhep]
  Float_t         *Phep;    //[Nhep][5]
  Float_t         *Vhep;    //[Nhep][4]
  Int_t           Irun;
  Int_t           Ievt;
  Float_t         Weight;
  Float_t         Xsecn;
  Int_t           Ifilter;
  Int_t           Nparam;
  Float_t         *Param;    //[Nparam]
};

//------------------------------------------------------------------------------

class HEPTreeReader
{
public:

  HEPTreeReader(TTree *tree, HEPEvent *event);
  ~HEPTreeReader();

  Long64_t GetEntries() const { return fChain ? static_cast<Long64_t>(fChain->GetEntries()) : 0; }

  Bool_t ReadEntry(Long64_t element);

private:

  void Notify();

  TTree *fChain;  // pointer to the analyzed TTree or TChain
  Int_t fCurrentTree; // current Tree number in a TChain

  HEPEvent *fEvent;

  deque< pair<TString, TBranch*> > fBranches;

};

//------------------------------------------------------------------------------

HEPTreeReader::HEPTreeReader(TTree *tree, HEPEvent *event) : fChain(tree), fCurrentTree(-1), fEvent(event)
{
  if(!fChain) return;

  TBranch *branch;
  TLeaf *leaf;
  TString name;
  Int_t i;

  name = "Nhep";
  branch = fChain->GetBranch(name);
  branch->SetAddress(&(event->Nhep));
  fBranches.push_back(make_pair(name, branch));

  TString intNames[3] = {"Idhep", "Jsmhep", "Jsdhep"};
  Int_t **intData[3] = {&event->Idhep, &event->Jsmhep, &event->Jsdhep};

  for(i = 0; i < 3; ++i)
  {
    name = intNames[i];
    branch = fChain->GetBranch(name);
    leaf = branch->GetLeaf(name);
    *intData[i] = new Int_t[leaf->GetNdata()];
    branch->SetAddress(*intData[i]);
    fBranches.push_back(make_pair(name, branch));
  }

  TString floatNames[2] = {"Phep", "Vhep"};
  Float_t **floatData[2] = {&event->Phep, &event->Vhep};

  for(i = 0; i < 2; ++i)
  {
    name = floatNames[i];
    branch = fChain->GetBranch(name);
    leaf = branch->GetLeaf(name);
    *floatData[i] = new Float_t[leaf->GetNdata()];
    branch->SetAddress(*floatData[i]);
    fBranches.push_back(make_pair(name, branch));
  }
}

//------------------------------------------------------------------------------

HEPTreeReader::~HEPTreeReader()
{
  if(!fChain) return;

  delete[] fEvent->Idhep;
  delete[] fEvent->Jsmhep;
  delete[] fEvent->Jsdhep;
  delete[] fEvent->Phep;
  delete[] fEvent->Vhep;
}

//------------------------------------------------------------------------------

Bool_t HEPTreeReader::ReadEntry(Long64_t element)
{
  if(!fChain) return kFALSE;

  Int_t treeEntry = fChain->LoadTree(element);
  if(treeEntry < 0) return kFALSE;

  if(fChain->IsA() == TChain::Class())
  {
    TChain *chain = static_cast<TChain*>(fChain);
    if(chain->GetTreeNumber() != fCurrentTree)
    {
      fCurrentTree = chain->GetTreeNumber();
      Notify();
    }
  }

  deque< pair<TString, TBranch*> >::iterator it_deque;
  TBranch *branch;

  for(it_deque = fBranches.begin(); it_deque != fBranches.end(); ++it_deque)
  {
    branch = it_deque->second;
    if(branch)
    {
      branch->GetEntry(treeEntry);
    }
  }

  return kTRUE;
}

//------------------------------------------------------------------------------

void HEPTreeReader::Notify()
{
  // Called when loading a new file.
  // Get branch pointers.
  if(!fChain) return;

  deque< pair<TString, TBranch*> >::iterator it_deque;

  for(it_deque = fBranches.begin(); it_deque != fBranches.end(); ++it_deque)
  {
    it_deque->second = fChain->GetBranch(it_deque->first);
    if(!it_deque->second)
    {
      cout << "** WARNING: cannot get branch '" << it_deque->first << "'" << endl;
    }
  }
}


//------------------------------------------------------------------------------

int main(int argc, char *argv[])
{
  char appName[] = "ExRootHEPEVTConverter";

  if(argc != 3)
  {
    cout << " Usage: " << appName << " input_file" << " output_file" << endl;
    cout << " input_file - list of HEPEVT files in ROOT format ('h101' tree)," << endl;
    cout << " output_file - output file in ROOT format." << endl;
    return 1;
  }

  gROOT->SetBatch();

  int appargc = 1;
  char *appargv[] = {appName};
  TApplication app(appName, &appargc, appargv);  

  TString inputFileList(argv[1]);
  TString outputFileName(argv[2]);
  string buffer;

  TChain *chain = new TChain("h101");
  if(FillChain(chain, inputFileList))
  {

    HEPEvent event;
    HEPTreeReader *treeReader = new HEPTreeReader(chain, &event);
  
    TFile *outputFile = TFile::Open(outputFileName, "RECREATE");
    ExRootTreeWriter *treeWriter = new ExRootTreeWriter(outputFile, "Analysis");
  
    ExRootTreeBranch *branchGen = treeWriter->NewBranch("Gen", TRootGenParticle::Class());
  
    Long64_t entry, allEntries = treeReader->GetEntries();
  
    cout << "** Chain contains " << allEntries << " events" << endl;
  
    Int_t address, particle;
    Double_t signPz, cosTheta;
  
    TLorentzVector momentum;
  
    TRootGenParticle *element;
  
    // Loop over all events
    for(entry = 0; entry < allEntries; ++entry)
    {
      // Load selected branches with data from specified event
      treeReader->ReadEntry(entry);
      treeWriter->Clear();
  
      if((entry % 100) == 0 && entry > 0 )
      {
        cout << "** Processing element # " << entry << endl;
      }
  
      for(particle = 0; particle < event.Nhep; ++particle)
      {
        element = (TRootGenParticle*) branchGen->NewEntry();
  
        element->PID = event.Idhep[particle];
        element->Status = event.Jsmhep[particle]/16000000
                        + event.Jsdhep[particle]/16000000*100;
        element->M1 = (event.Jsmhep[particle]%16000000)/4000 - 1;
        element->M2 = event.Jsmhep[particle]%4000 - 1;
        element->D1 = (event.Jsdhep[particle]%16000000)/4000 - 1;
        element->D2 = event.Jsdhep[particle]%4000 - 1;
  
        address = particle*5;
        element->E = event.Phep[address + 3];
        element->Px = event.Phep[address + 0];
        element->Py = event.Phep[address + 1];
        element->Pz = event.Phep[address + 2];
  
        momentum.SetPxPyPzE(element->Px, element->Py, element->Pz, element->E);

        cosTheta = TMath::Abs(momentum.CosTheta());
        signPz = (momentum.Pz() >= 0.0) ? 1.0 : -1.0;

        element->PT = momentum.Perp();
        element->Phi = momentum.Phi();
        element->Eta = (cosTheta == 1.0 ? signPz*999.9 : momentum.Eta());
        element->Rapidity = (cosTheta == 1.0 ? signPz*999.9 : momentum.Rapidity());

        address = particle*4;
        element->T = event.Vhep[address + 3];
        element->X = event.Vhep[address + 0];
        element->Y = event.Vhep[address + 1];
        element->Z = event.Vhep[address + 2];
      }
      treeWriter->Fill();
    }
  
    treeWriter->Write();
  
    cout << "** Exiting..." << endl;
  
    delete treeWriter;
    delete outputFile;
    delete treeReader;
  }
  delete chain;
}



