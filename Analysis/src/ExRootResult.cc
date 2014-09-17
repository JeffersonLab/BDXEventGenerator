
/** \class ExRootResult
 *
 *  Class simplifying work with histograms
 *
 *  $Date: 2007/07/23 14:36:04 $
 *  $Revision: 1.5 $
 *
 *
 *  \author P. Demin - UCL, Louvain-la-Neuve
 *
 */

#include "ExRootAnalysis/ExRootResult.h"

#include "ExRootAnalysis/ExRootUtilities.h"

#include "TROOT.h"
#include "TFile.h"
#include "TClass.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TPaveText.h"
#include "TPaveStats.h"
#include "TList.h"
#include "TH2.h"
#include "THStack.h"
#include "TProfile.h"
#include "TObjArray.h"
#include "TFolder.h"

#include <algorithm>
#include <iostream>

using namespace std;

const Font_t kExRootFont = 42;
const Float_t kExRootFontSize = 0.04;
const Color_t kExRootBackgroundColor = 10;

//------------------------------------------------------------------------------

static void DeleteTObjectPtr(TObject *x)
{
  delete x;
}

//------------------------------------------------------------------------------

ExRootResult::ExRootResult() : fCanvas(0), fFolder(0)
{

}

//------------------------------------------------------------------------------

ExRootResult::~ExRootResult()
{
  for_each(fPool.begin(), fPool.end(), DeleteTObjectPtr);

  if(fCanvas) delete fCanvas;
}

//------------------------------------------------------------------------------

void ExRootResult::Reset()
{

}

//------------------------------------------------------------------------------

void ExRootResult::Write(const char *fileName)
{
  TObject *object;
  TDirectory *currentDirectory = gDirectory; 
  TFile *file = new TFile(fileName, "RECREATE");
  file->cd();
  map<TObject*, TObjArray*>::iterator it_plots;
  for(it_plots = fPlots.begin(); it_plots != fPlots.end(); ++it_plots)
  {
    object = it_plots->first;
    object->Write();    
  }
  currentDirectory->cd();
  delete file;
}

//------------------------------------------------------------------------------

void ExRootResult::CreateCanvas()
{
  TDirectory *currentDirectory = gDirectory;

  // Graphics style parameters to avoid grey background on figures
  gStyle->SetCanvasColor(kExRootBackgroundColor);
  gStyle->SetStatColor(kExRootBackgroundColor);
  //  gStyle->SetTitleColor(kExRootBackgroundColor);
  gStyle->SetPadColor(kExRootBackgroundColor);

  gStyle->SetPadTopMargin(0.10);
  gStyle->SetPadRightMargin(0.05);
  gStyle->SetPadBottomMargin(0.15);
  gStyle->SetPadLeftMargin(0.15);

  gStyle->SetStatFont(kExRootFont);
  gStyle->SetStatFontSize(kExRootFontSize);

  gStyle->SetTitleFont(kExRootFont, "");
  gStyle->SetTitleFont(kExRootFont, "X");
  gStyle->SetTitleFont(kExRootFont, "Y");
  gStyle->SetTitleFont(kExRootFont, "Z");
  gStyle->SetTitleSize(kExRootFontSize, "");
  gStyle->SetTitleSize(kExRootFontSize, "X");
  gStyle->SetTitleSize(kExRootFontSize, "Y");
  gStyle->SetTitleSize(kExRootFontSize, "Z");

  gStyle->SetLabelFont(kExRootFont, "X");
  gStyle->SetLabelFont(kExRootFont, "Y");
  gStyle->SetLabelFont(kExRootFont, "Z");
  gStyle->SetLabelSize(kExRootFontSize, "X");
  gStyle->SetLabelSize(kExRootFontSize, "Y");
  gStyle->SetLabelSize(kExRootFontSize, "Z");

  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);

  gStyle->SetTextFont(kExRootFont);
  gStyle->SetTextSize(kExRootFontSize);

  gStyle->SetOptStat(111110);
  // gStyle->SetOptFit(101);
  
  fCanvas = static_cast<TCanvas*>(gROOT->FindObject("c1"));
  if(fCanvas)
  {
    fCanvas->Clear();
    fCanvas->UseCurrentStyle();
    fCanvas->SetWindowSize(800, 650);
  }
  else
  {
    fCanvas = new TCanvas("c1", "c1", 800, 650);
  }
  fCanvas->SetLogy(0);
  fCanvas->SetHighLightColor(kExRootBackgroundColor);

  currentDirectory->cd();
}

//------------------------------------------------------------------------------

TCanvas *ExRootResult::GetCanvas()
{
  if(!fCanvas) CreateCanvas();
  return fCanvas;
}

//------------------------------------------------------------------------------

void ExRootResult::Attach(TObject *plot, TObject *object)
{
  if(!plot) return;

  map<TObject*, TObjArray*>::iterator it_plots = fPlots.find(plot);
  if(it_plots != fPlots.end())
  {
    TObjArray *attachment = it_plots->second;
    if(!attachment)
    {
      attachment = new TObjArray();
      it_plots->second = attachment;
    }
    attachment->Add(object);
  }
}


//------------------------------------------------------------------------------

void ExRootResult::PrintPlot(TObject *plot, const char *sufix, const char *format)
{
  if(!plot) return;

  TCanvas *canvas = GetCanvas();
  TH1 *histogram = 0;

  if(plot->IsA()->InheritsFrom(TH1::Class()))
  {
    histogram = static_cast<TH1*>(plot);
  }

  map<TObject*, PlotSettings>::iterator it_settings = fSettings.find(plot);
  if(it_settings != fSettings.end())
  {
    canvas->SetLogx(it_settings->second.logx);
    if(histogram == 0 || histogram->Integral() > 0.0)
    {
      canvas->SetLogy(it_settings->second.logy);
    }
    else
    {
      canvas->SetLogy(0);
    }
  }
  
  map<TObject*, TObjArray*>::iterator it_plots = fPlots.find(plot);
  if(it_plots != fPlots.end())
  {
    TObjArray *attachment = it_plots->second;
    if(attachment)
    {
      TIter iterator(attachment);
      TObject *object;
      while((object = iterator()))
      {
        object->Draw();
      }
    }
  }

  TString name = plot->GetName();
  canvas->Print(name + sufix + "." + format);
}

//------------------------------------------------------------------------------

void ExRootResult::Print(const char *format)
{
  TObjArray *attachment;
  TObject *object;
  TH1 *histogram;
  TPaveStats *stats;
  TString name;

  TCanvas *canvas = GetCanvas();

  map<TObject*, TObjArray*>::iterator it_plots;
  map<TObject*, PlotSettings>::iterator it_settings;

  for(it_plots = fPlots.begin(); it_plots != fPlots.end(); ++it_plots)
  {
    object = it_plots->first;
    attachment = it_plots->second;
    name = object->GetName();
    histogram = 0;

    if(object->IsA()->InheritsFrom(TH1::Class()))
    {
      histogram = static_cast<TH1*>(object);
    }

    it_settings = fSettings.find(object);
    if(it_settings != fSettings.end())
    {
      canvas->SetLogx(it_settings->second.logx);
      if(histogram == 0 || histogram->Integral() > 0.0)
      {
        canvas->SetLogy(it_settings->second.logy);
      }
      else
      {
        canvas->SetLogy(0);
      }
    }

    if(object->IsA()->InheritsFrom(THStack::Class()))
    {
      object->Draw("nostack");
    }
    else
    {
      object->Draw();
    }

    canvas->Update();

    if(histogram)
    {
      stats = static_cast<TPaveStats*>(histogram->GetListOfFunctions()->FindObject("stats"));
      if(stats)
      {
        stats->SetX1NDC(0.67);
        stats->SetX2NDC(0.99);
        stats->SetY1NDC(0.77);
        stats->SetY2NDC(0.99);
        stats->SetTextFont(kExRootFont);
        stats->SetTextSize(kExRootFontSize);
        canvas->Draw();
      }
    }
    if(attachment)
    {
      TIter iterator(attachment);
      while((object = iterator()))
      {
        object->Draw();
      }
    }
    canvas->Print(name + "." + format);
  }
}

//------------------------------------------------------------------------------

TH1 *ExRootResult::AddHist1D(const char *name, const char *title,
                             const char *xlabel, const char *ylabel,
                             Int_t nxbins, Axis_t xmin, Axis_t xmax,
                             Int_t logx, Int_t logy)
{
  TH1F *hist = new TH1F(name, title, nxbins, xmin, xmax);
  PlotSettings settings;
  settings.logx = logx;
  settings.logy = logy;
  
  fPool.insert(hist);
  hist->GetXaxis()->SetTitle(xlabel);
  hist->GetYaxis()->SetTitle(ylabel);
  fPlots[hist] = 0;
  fSettings[hist] = settings;
  HistStyle(hist, kFALSE);
  return hist;
}

//------------------------------------------------------------------------------

TH1 *ExRootResult::AddHist1D(const char *name, const char *title,
                             const char *xlabel, const char *ylabel,
                             Int_t nxbins, const Float_t *bins,
                             Int_t logx, Int_t logy)
{
  TH1F *hist = new TH1F(name, title, nxbins, bins);
  PlotSettings settings;
  settings.logx = logx;
  settings.logy = logy;

  fPool.insert(hist);
  hist->GetXaxis()->SetTitle(xlabel);
  hist->GetYaxis()->SetTitle(ylabel);
  fPlots[hist] = 0;
  fSettings[hist] = settings;
  HistStyle(hist, kFALSE);
  if(fFolder) fFolder->Add(hist);
  return hist;
}

//------------------------------------------------------------------------------

TProfile *ExRootResult::AddProfile(const char *name, const char *title,
                                   const char *xlabel, const char *ylabel,
                                   Int_t nxbins, Axis_t xmin, Axis_t xmax,
                                   Int_t logx, Int_t logy)
{
  TProfile *profile = new TProfile(name, title, nxbins, xmin, xmax);
  PlotSettings settings;
  settings.logx = logx;
  settings.logy = logy;

  fPool.insert(profile);
  profile->GetXaxis()->SetTitle(xlabel);
  profile->GetYaxis()->SetTitle(ylabel);
  fPlots[profile] = 0;
  fSettings[profile] = settings;
  HistStyle(profile, kFALSE);
  if(fFolder) fFolder->Add(profile);
  return profile;
}

//------------------------------------------------------------------------------

TH2 *ExRootResult::AddHist2D(const char *name, const char *title,
                             const char *xlabel, const char *ylabel,
                             Int_t nxbins, Axis_t xmin, Axis_t xmax,
                             Int_t nybins, Axis_t ymin, Axis_t ymax,
                             Int_t logx, Int_t logy)
{
  TH2F *hist = new TH2F(name, title, nxbins, xmin, xmax, nybins, ymin, ymax);
  PlotSettings settings;
  settings.logx = logx;
  settings.logy = logy;

  fPool.insert(hist);
  hist->GetXaxis()->SetTitle(xlabel);
  hist->GetYaxis()->SetTitle(ylabel);
  fPlots[hist] = 0;
  fSettings[hist] = settings;
  HistStyle(hist, kFALSE);
  if(fFolder) fFolder->Add(hist);
  return hist;
}

//------------------------------------------------------------------------------

THStack *ExRootResult::AddHistStack(const char *name, const char *title)
{
  THStack *stack = new THStack(name, title);
//  segmentaion violation when deleting stack in ~ExRootResult()
//  fPool.insert(stack);
  fPlots[stack] = 0;
  if(fFolder) fFolder->Add(stack);
  return stack;
}

//------------------------------------------------------------------------------

TPaveText *ExRootResult::AddComment(Double_t x1, Double_t y1, Double_t x2, Double_t y2)
{
  TPaveText *comment = new TPaveText(x1, y1, x2, y2, "brNDC");
  fPool.insert(comment);
  comment->SetTextSize(kExRootFontSize);
  comment->SetTextFont(kExRootFont);
  comment->SetTextAlign(22);
  comment->SetFillColor(kExRootBackgroundColor);
  comment->SetBorderSize(2);
  return comment;
}

//------------------------------------------------------------------------------

TLegend *ExRootResult::AddLegend(Double_t x1, Double_t y1, Double_t x2, Double_t y2)
{
  TLegend *legend = new TLegend(x1, y1, x2, y2);
  fPool.insert(legend);
  legend->SetTextSize(kExRootFontSize);
  legend->SetTextFont(kExRootFont);
  legend->SetFillColor(kExRootBackgroundColor);
  legend->SetBorderSize(2);
  return legend;
}

//------------------------------------------------------------------------------


