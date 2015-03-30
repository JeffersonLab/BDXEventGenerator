
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

#include "ExRootAnalysis/ExRootStream.h"

#include <stdlib.h>
#include <limits.h>
#include <string.h>
#include <stdio.h>
#include <errno.h>
#include <math.h>

#include <iostream>

using namespace std;

//------------------------------------------------------------------------------

bool ExRootStream::fFirstLongMin = true;
bool ExRootStream::fFirstLongMax = true;
bool ExRootStream::fFirstHugePos = true;
bool ExRootStream::fFirstHugeNeg = true;
bool ExRootStream::fFirstZero = true;

//------------------------------------------------------------------------------

ExRootStream::ExRootStream(char *buffer) :
  fBuffer(buffer)
{
}

//------------------------------------------------------------------------------

bool ExRootStream::ReadDbl(double &value)
{
  char *start = fBuffer;
  errno = 0;
  value = strtod(start, &fBuffer);
  if(errno == ERANGE)
  {
    if(fFirstHugePos && value == HUGE_VAL)
    {
      fFirstHugePos = false; 
      cout << "** WARNING: too large positive value, return " << value << endl;
    }
    else if(fFirstHugeNeg && value == -HUGE_VAL)
    {
      fFirstHugeNeg = false; 
      cout << "** WARNING: too large negative value, return " << value << endl;
    }
    else if(fFirstZero)
    {
      fFirstZero = false; 
      value = 0.0;
      cout << "** WARNING: too small value, return " << value << endl;
    }
  }
  return start != fBuffer;
}

//------------------------------------------------------------------------------

bool ExRootStream::ReadInt(int &value)
{
  char *start = fBuffer;
  errno = 0;
  value = strtol(start, &fBuffer, 10);
  if(errno == ERANGE)
  {
    if(fFirstLongMin && value == LONG_MIN)
    {
      fFirstLongMin = false; 
      cout << "** WARNING: too large positive value, return " << value << endl;
    }
    else if(fFirstLongMax && value == LONG_MAX)
    {
      fFirstLongMax = false; 
      cout << "** WARNING: too large negative value, return " << value << endl;
    }
  }
  return start != fBuffer;
}

//------------------------------------------------------------------------------
