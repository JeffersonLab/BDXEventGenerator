
/** \class ExRootFilter
 *
 *  Class simplifying classification and subarrays handling
 *
 *  $Date: 2010/05/05 15:32:18 $
 *  $Revision: 1.5 $
 *
 *
 *  \author P. Demin - UCL, Louvain-la-Neuve
 *
 */

#include "ExRootAnalysis/ExRootFilter.h"
#include "ExRootAnalysis/ExRootClassifier.h"

#include "TSeqCollection.h"
#include "TObjArray.h"

using namespace std;

typedef map<Int_t, TObjArray*> TCategoryMap;
typedef map<ExRootClassifier*, pair<Bool_t, TCategoryMap> > TClassifierMap;

ExRootFilter::ExRootFilter(const TSeqCollection *collection) :
  fCollection(collection)
{
  fIter = fCollection->MakeIterator();
}

//------------------------------------------------------------------------------

ExRootFilter::~ExRootFilter()
{
  TClassifierMap::iterator it_map;
  TCategoryMap::iterator it_submap;
  for(it_map = fMap.begin(); it_map != fMap.end(); ++it_map)
  {
    for(it_submap = it_map->second.second.begin();
        it_submap != it_map->second.second.end(); ++it_submap)
    {
      delete (it_submap->second);
    }
  }

  delete fIter;
}

//------------------------------------------------------------------------------

void ExRootFilter::Reset(ExRootClassifier *classifier)
{
  TClassifierMap::iterator it_map;
  TCategoryMap::iterator it_submap;
  if(classifier)
  {
    it_map = fMap.find(classifier);
    if(it_map != fMap.end())
    {
      it_map->second.first = kTRUE;
      for(it_submap = it_map->second.second.begin();
          it_submap != it_map->second.second.end(); ++it_submap)
      {
        it_submap->second->Clear();
      }
    }
  }
  else
  {
    for(it_map = fMap.begin(); it_map != fMap.end(); ++it_map)
    {
      it_map->second.first = kTRUE;
      for(it_submap = it_map->second.second.begin();
          it_submap != it_map->second.second.end(); ++it_submap)
      {
        it_submap->second->Clear();
      }
    }
  }  
}

//------------------------------------------------------------------------------

TObjArray *ExRootFilter::GetSubArray(ExRootClassifier *classifier, Int_t category)
{
  Int_t result;
  TObject *element;
  TObjArray *array;
  TCategoryMap::iterator it_submap;
  pair<TCategoryMap::iterator, bool> pair_submap;
  pair<TClassifierMap::iterator, bool> pair_map;

  TClassifierMap::iterator it_map = fMap.find(classifier);
  if(it_map == fMap.end())
  {
    pair_map = fMap.insert(make_pair(classifier, make_pair(kTRUE, TCategoryMap())));
    if(!pair_map.second) throw FilterExeption();

    it_map = pair_map.first;
  }

  if(it_map->second.first)
  {
    it_map->second.first = kFALSE;
    fIter->Reset();
    while((element = fIter->Next()) != 0)
    {
      result = classifier->GetCategory(element);
      if(result < 0) continue;
      it_submap = it_map->second.second.find(result);
      if(it_submap == it_map->second.second.end())
      {
        array = new TObjArray(fCollection->GetSize());
        pair_submap = it_map->second.second.insert(make_pair(result, array));
        if(!pair_submap.second) throw FilterExeption();

        it_submap = pair_submap.first;
      }
      it_submap->second->Add(element);
    }
  }

  it_submap = it_map->second.second.find(category);
  return (it_submap != it_map->second.second.end()) ? it_submap->second : 0;  
}

//------------------------------------------------------------------------------

