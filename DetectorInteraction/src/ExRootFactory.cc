
/** \class ExRootFactory
 *
 *  Class handling creation of ExRootCandidate,
 *  ExRootCandList and all other objects.
 *
 *  $Date: 2006/09/21 13:10:52 $
 *  $Revision: 1.1 $
 *
 *
 *  \author P. Demin - UCL, Louvain-la-Neuve
 *
 */

#include "ExRootAnalysis/ExRootTreeWriter.h"
#include "ExRootAnalysis/ExRootTreeBranch.h"

#include "ExRootAnalysis/ExRootFactory.h"

#include "TClass.h"
#include "TObjArray.h"

using namespace std;

//------------------------------------------------------------------------------

ExRootFactory::ExRootFactory() :
  fTreeWriter(0), fPermanentObjArrays(0)
{
  fTreeWriter = new ExRootTreeWriter();
  fPermanentObjArrays = fTreeWriter->NewFactory("PermanentObjArrays", TObjArray::Class());
}

//------------------------------------------------------------------------------

ExRootFactory::~ExRootFactory()
{
  if(fTreeWriter) delete fTreeWriter;
}

//------------------------------------------------------------------------------

void ExRootFactory::Clear()
{
  map<const TClass *, ExRootTreeBranch *>::iterator it_map;
  for(it_map = fMakers.begin(); it_map != fMakers.begin(); ++it_map)
  {
    it_map->second->Clear();
  }

  set<TObjArray *>::iterator it_set;
  for(it_set = fPool.begin(); it_set != fPool.end(); ++it_set)
  {
    (*it_set)->Clear();
  }
}

//------------------------------------------------------------------------------

TObjArray *ExRootFactory::NewArray()
{
  return New<TObjArray>();
}

//------------------------------------------------------------------------------

TObjArray *ExRootFactory::NewPermanentArray()
{
  TObjArray *array = static_cast<TObjArray *>(fPermanentObjArrays->NewEntry());
  fPool.insert(array);
  return array;
}

//------------------------------------------------------------------------------

TObject *ExRootFactory::New(TClass *cl)
{
  ExRootTreeBranch *maker = 0;
  map<const TClass *, ExRootTreeBranch *>::iterator it = fMakers.find(cl);
  
  if(it != fMakers.end())
  {
    maker = it->second;
  }
  else
  {
    maker = fTreeWriter->NewFactory(cl->GetName(), cl);
    fMakers.insert(make_pair(cl, maker));
  }

  return maker->NewEntry();
}


