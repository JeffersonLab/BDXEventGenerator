#ifndef ExRootUtilities_h
#define ExRootUtilities_h

/** \class ExRootUtilities
 *
 *  Functions simplifying ROOT tree analysis
 *
 *  $Date: 2006/09/21 13:06:18 $
 *  $Revision: 1.2 $
 *
 *
 *  \author P. Demin - UCL, Louvain-la-Neuve
 *
 */

#include "Rtypes.h"

class TH1;
class TChain;

void HistStyle(TH1 *hist, Bool_t stats = kTRUE);

Bool_t FillChain(TChain *chain, const char *inputFileList);

#endif // ExRootUtilities_h
