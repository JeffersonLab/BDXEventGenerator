#ifndef ExRootFactory_h
#define ExRootFactory_h

/** \class ExRootFactory
 *
 *  Class handling creation of ExRootCandidate,
 *  ExRootCandList and all other objects.
 *
 *  $Date: 2006/09/21 13:08:01 $
 *  $Revision: 1.1 $
 *
 *
 *  \author P. Demin - UCL, Louvain-la-Neuve
 *
 */

#include "TNamed.h"

#include <map>
#include <set>

class TObjArray;

class ExRootTreeBranch;
class ExRootTreeWriter;

class ExRootFactory: public TNamed
{
public:
  
  ExRootFactory();
  ~ExRootFactory();

  void Clear();

  TObjArray *NewArray();

  TObjArray *NewPermanentArray();

  TObject *New(TClass *cl);

  template<typename T>
  T *New() { return static_cast<T *>(New(T::Class())); }


private:

  ExRootTreeWriter *fTreeWriter; //!

  ExRootTreeBranch *fPermanentObjArrays; //!

  std::map<const TClass *, ExRootTreeBranch *> fMakers; //!
  std::set<TObjArray *> fPool; //!

  ClassDef(ExRootFactory, 1)
};

#endif /* ExRootFactory */

