
/** \class ExRootProgressBar
 *
 *  Class showing progress bar
 *
 *  $Date: 2006/12/13 18:49:30 $
 *  $Revision: 1.2 $
 *
 *
 *  \author P. Demin - UCL, Louvain-la-Neuve
 *
 */

#include "ExRootAnalysis/ExRootProgressBar.h"

#include "TSystem.h"

#include <iostream>

#include <string.h>
#include <stdio.h>

using namespace std;

ExRootProgressBar::ExRootProgressBar(Long64_t entries, Int_t width) :
  fEntries(entries), fWidth(width), fTime(0), fHashes(0), fBar(0)
{
  fBar = new char[width + 1];
  memset(fBar, '-', width);
  fBar[width] = 0;

}

//------------------------------------------------------------------------------

ExRootProgressBar::~ExRootProgressBar()
{
  if(fBar) delete[] fBar;
}

//------------------------------------------------------------------------------

void ExRootProgressBar::Update(Long64_t entry)
{
  ULong64_t time = gSystem->Now();

  if(time < fTime + 1000 && entry < fEntries - 1) return;

  fTime = time;

  Int_t hashes = Int_t((entry + 1.0)/fEntries*fWidth);

  if(hashes > fHashes)
  {
    memset(fBar + fHashes, '#', hashes - fHashes);
    fHashes = hashes;
  }

/*
  cerr << "[" << fBar << "] (";
  cerr.setf(ios::fixed);
  cerr.precision(2);
  cerr << (entry + 1.0)/fEntries*100.0 << "%) : ";
  cerr << entry + 1 << "/" << fEntries;
  cerr << " events processed\r" << flush;
*/

  fprintf(stderr, "[%s] (%.2f%%) : %lli/%lli entries processed\r", fBar,
          (entry + 1.0)/fEntries*100.0, entry + 1, fEntries);
  fflush(stderr);
}

//------------------------------------------------------------------------------

void ExRootProgressBar::Finish()
{
  fprintf(stderr, "\n");
  fflush(stderr);
}

//------------------------------------------------------------------------------

