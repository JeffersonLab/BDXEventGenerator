/** \class ExRootAnalysisLinkDef
 *
 *  Lists classes to be included in cint dicitonary
 *
 *  $Date: 2006/09/22 14:15:42 $
 *  $Revision: 1.1 $
 *
 *  
 *  \author P. Demin - UCL, Louvain-la-Neuve
 *
 */

#include "ExRootAnalysis/ExRootClasses.h"
#include "ExRootAnalysis/ExRootTreeWriter.h"
#include "ExRootAnalysis/ExRootTreeBranch.h"

#ifdef __CINT__

#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

#pragma link C++ class TSortableObject+;
#pragma link C++ class TRootLHEFEvent+;
#pragma link C++ class TRootLHEFParticle+;
#pragma link C++ class TRootGenParticle+;
#pragma link C++ class TRootTrack+;
#pragma link C++ class TRootCalTower+;
#pragma link C++ class TRootMissingET+;
#pragma link C++ class TRootCalCluster+;
#pragma link C++ class TRootPhoton+;
#pragma link C++ class TRootElectron+;
#pragma link C++ class TRootMuon+;
#pragma link C++ class TRootTau+;
#pragma link C++ class TRootJet+;
#pragma link C++ class TRootHeavy+;

#pragma link C++ class ExRootTreeBranch+;
#pragma link C++ class ExRootTreeWriter+;

#endif

