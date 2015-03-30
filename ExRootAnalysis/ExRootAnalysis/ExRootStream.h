#ifndef ExRootStream_h
#define ExRootStream_h

/** \class ExRootStream
 *
 *  Provides an interface to manipulate c strings as if they were input streams
 *
 *
 *  $Date: 2012-11-15 13:57:55 $
 *  $Revision: 1.1 $
 *
 *
 *  \author P. Demin - UCL, Louvain-la-Neuve
 *
 */

class ExRootStream
{
public:

  ExRootStream(char *buffer);

  bool ReadDbl(double &value);
  bool ReadInt(int &value);

private:

  char *fBuffer;
  
  static bool fFirstLongMin;
  static bool fFirstLongMax;
  static bool fFirstHugePos;
  static bool fFirstHugeNeg;
  static bool fFirstZero;
};

#endif // ExRootStream_h


