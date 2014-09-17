#ifndef ExRootProgressBar_h
#define ExRootProgressBar_h

#include "Rtypes.h"

class ExRootProgressBar
{
public:

  ExRootProgressBar(Long64_t entries, Int_t width = 25);
  ~ExRootProgressBar();

  void Update(Long64_t entry);
  void Finish();

private:

  Long64_t fEntries;
  Int_t fWidth;

  ULong64_t fTime;
  Int_t fHashes;

  char *fBar;

};

#endif /* ExRootProgressBar */

