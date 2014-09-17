
#include "modules/MadGraphAnalysis.h"


#include "ExRootAnalysis/ExRootResult.h"
#include "ExRootAnalysis/ExRootClasses.h"

#include "TClonesArray.h"

#include "TH1.h"
#include "TH2.h"
#include "TString.h"
#include "TCanvas.h"
#include "TLorentzVector.h"

#include <iostream>

using namespace std;

//------------------------------------------------------------------------------

MadGraphAnalysis::MadGraphAnalysis()
{
}

//------------------------------------------------------------------------------

MadGraphAnalysis::~MadGraphAnalysis()
{
}

//------------------------------------------------------------------------------

void MadGraphAnalysis::Init()
{
  fOutputFileName = GetString("OutputFile", "madgraph_plots.root");

  // import array with output from filter/classifier module

  fInputArray = ImportArray(GetString("InputArray", "classification/particles"));

  // import ROOT tree branch

  fBranchEvent = UseBranch("Event");
}

//------------------------------------------------------------------------------

void MadGraphAnalysis::Finish()
{
  GetPlots()->Write(fOutputFileName);

  GetPlots()->GetCanvas()->SetLogy(1);
  GetPlots()->Print();
  GetPlots()->GetCanvas()->SetLogy(0);
}

//------------------------------------------------------------------------------

void MadGraphAnalysis::Event()
{
  TRootLHEFParticle *particle1 = 0, *particle2 = 0;
  ParticleHistograms *histogramsParticle = 0;
  PairHistograms *histogramsPair = 0;
  const TObjArray *array1 = 0, *array2 = 0;
  TLorentzVector vector1, vector2;
  TString name;
  Int_t entry1, entry2, maxEntry1, maxEntry2;
  Int_t inputEntry1, inputEntry2, maxInputEntry;

  // read event weight

  Double_t weight = 0.0;
  Double_t pt1, pt2, dr, rapidity, signPz;

  TRootLHEFEvent *eventInfo = 0;

  if(fBranchEvent->GetEntriesFast() == 1)
  {
    eventInfo = static_cast<TRootLHEFEvent*>(fBranchEvent->At(0));

    weight = eventInfo->Weight;
  }

  // fill histograms
  
  maxInputEntry = fInputArray->GetEntriesFast();
  for(inputEntry1 = 0; inputEntry1 < maxInputEntry; ++inputEntry1)
  {
    array1 = static_cast<TObjArray*>(fInputArray->At(inputEntry1));

    maxEntry1 = array1->GetEntriesFast();
    for(entry1 = 0; entry1 < maxEntry1; ++entry1)
    {
      particle1 = static_cast<TRootLHEFParticle*>(array1->At(entry1));
      vector1.SetPxPyPzE(particle1->Px, particle1->Py, particle1->Pz, particle1->E);

      pt1 = vector1.Pt();
      signPz = (vector1.Pz() >= 0.0) ? 1.0 : -1.0;
      rapidity = (pt1 == 0.0 ? signPz*999.9 : vector1.Rapidity());

      // fill histograms for single particles

      histogramsParticle = GetParticleHistograms(array1->GetName(), entry1 + 1);

      histogramsParticle->fParticlePt->Fill(pt1, weight);

      histogramsParticle->fParticleRapidity->Fill(rapidity, weight);
      
      if(particle1->Status == 2) continue;

      for(entry2 = entry1 + 1; entry2 < maxEntry1; ++entry2)
      {
        particle2 = static_cast<TRootLHEFParticle*>(array1->At(entry2));
 
        if(particle2->Status == 2) continue;

        vector2.SetPxPyPzE(particle2->Px, particle2->Py,
                           particle2->Pz, particle2->E);

        pt2 = vector2.Pt();
        dr = (pt1 == 0.0 || pt2 == 0.0 ? 999.9 : vector1.DeltaR(vector2));

        // fill histograms for pairs of particles

        histogramsPair = GetPairHistograms(array1->GetName(), entry1 + 1,
                                           array1->GetName(), entry2 + 1);

        histogramsPair->fPairDeltaR->Fill(dr, weight);
        histogramsPair->fPairMass->Fill((vector1 + vector2).M(), weight);
      }

      for(inputEntry2 = inputEntry1 + 1; inputEntry2 < maxInputEntry; ++inputEntry2)
      {
        array2 = static_cast<TObjArray*>(fInputArray->At(inputEntry2));

        maxEntry2 = array2->GetEntriesFast();
        for(entry2 = 0; entry2 < maxEntry2; ++entry2)
        {
          particle2 = static_cast<TRootLHEFParticle*>(array2->At(entry2));

          if(particle2->Status == 2) continue;

          vector2.SetPxPyPzE(particle2->Px, particle2->Py, particle2->Pz, particle2->E);

          pt2 = vector2.Pt();
          dr = (pt1 == 0.0 || pt2 == 0.0 ? 999.9 : vector1.DeltaR(vector2));

          // fill histograms for pairs of particles

          histogramsPair = GetPairHistograms(array1->GetName(), entry1 + 1,
                                             array2->GetName(), entry2 + 1);

          histogramsPair->fPairDeltaR->Fill(dr, weight);
          histogramsPair->fPairMass->Fill((vector1 + vector2).M(), weight);
        }
      }
    }
  }
}

//------------------------------------------------------------------------------

void MadGraphAnalysis::BookParticleHistograms(MadGraphAnalysis::ParticleHistograms *histograms,
                                                 const char *name, const char *title)
{
  ExRootResult *result = GetPlots();
  histograms->fParticlePt = result->AddHist1D(Form("pt_%s", name),
                                              Form("P_{T}(%s)", title),
                                              Form("P_{T}(%s), GeV/c", title),
                                              "pb/bin",
                                              60, 0.0, 300.0, 0, 1);
  histograms->fParticlePt->SetStats(kTRUE);

  histograms->fParticleRapidity = result->AddHist1D(Form("y_%s", name),
                                                    Form("y(%s)", title),
                                                    Form("y(%s)", title),
                                                    "pb/bin",
                                                    100, -5.0, 5.0, 0, 0);
  histograms->fParticleRapidity->SetStats(kTRUE);
  
}

//------------------------------------------------------------------------------

void MadGraphAnalysis::BookPairHistograms(MadGraphAnalysis::PairHistograms *histograms,
                                             const char *name, const char *title)
{
  ExRootResult *result = GetPlots();
  histograms->fPairDeltaR = result->AddHist1D(Form("dr_%s", name),
                                              Form("#DeltaR(%s)", title),
                                              Form("#DeltaR(%s)", title),
                                              "pb/bin",
                                              70, 0.0, 7.0, 0, 1);
  histograms->fPairDeltaR->SetStats(kTRUE);

  histograms->fPairMass = result->AddHist1D(Form("mass_%s", name),
                                            Form("M_{inv}(%s)", title),
                                            Form("M_{inv}(%s), GeV/c^{2}", title),
                                            "pb/bin",
                                            120, 0.0, 600.0, 0, 1);
  histograms->fPairMass->SetStats(kTRUE);

}

//------------------------------------------------------------------------------

MadGraphAnalysis::ParticleHistograms *MadGraphAnalysis::GetParticleHistograms(const char *module, Int_t number)
{
  map<TString, ParticleHistograms *>::iterator itParticleHistogramsMap;
  ParticleHistograms *histograms = 0;
  TString name = Form("%s_%d", module, number);
  name.ReplaceAll("{", "");
  name.ReplaceAll("}", "");
  name.ReplaceAll("^", "");
  name.ReplaceAll("#bar", "anti_");
  name.ReplaceAll("#", "");
  TString title = Form("%s_{%d}", module, number);
  itParticleHistogramsMap = fParticleHistogramsMap.find(name);
  if(itParticleHistogramsMap == fParticleHistogramsMap.end())
  {
    histograms = new ParticleHistograms;

    BookParticleHistograms(histograms, name, title);

    fParticleHistogramsMap[name] = histograms;
  }
  else
  {
    histograms = itParticleHistogramsMap->second;
  }
  return histograms;
}

//------------------------------------------------------------------------------

MadGraphAnalysis::PairHistograms *MadGraphAnalysis::GetPairHistograms(const char *module1, Int_t number1,
                                                                            const char *module2, Int_t number2)
{
  map<TString, PairHistograms *>::iterator itPairHistogramsMap;
  PairHistograms *histograms = 0;
  TString name = Form("%s_%d_%s_%d", module1, number1, module2, number2);
  name.ReplaceAll("{", "");
  name.ReplaceAll("}", "");
  name.ReplaceAll("^", "");
  name.ReplaceAll("#bar", "anti_");
  name.ReplaceAll("#", "");
  TString title = Form("%s_{%d}, %s_{%d}", module1, number1, module2, number2);
  itPairHistogramsMap = fPairHistogramsMap.find(name);
  if(itPairHistogramsMap == fPairHistogramsMap.end())
  {
    histograms = new PairHistograms;

    BookPairHistograms(histograms, name, title);

    fPairHistogramsMap[name] = histograms;
  }
  else
  {
    histograms = itPairHistogramsMap->second;
  }
  return histograms;
}

