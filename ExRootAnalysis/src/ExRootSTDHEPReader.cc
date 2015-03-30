
/** \class ExRootSTDHEPReader
 *
 *  Reads STDHEP file
 *
 *
 *  $Date$
 *  $Revision$
 *
 *
 *  \author P. Demin - UCL, Louvain-la-Neuve
 *
 */

#include "ExRootAnalysis/ExRootSTDHEPReader.h"

#include <stdexcept>
#include <iostream>
#include <sstream>

#include <stdio.h>
#include <errno.h>
#include <rpc/types.h>
#include <rpc/xdr.h>

#include "TLorentzVector.h"

#include "ExRootAnalysis/ExRootClasses.h"
#include "ExRootAnalysis/ExRootFactory.h"

#include "ExRootAnalysis/ExRootTreeBranch.h"

using namespace std;

static const int kBufferSize  = 1000000;

//---------------------------------------------------------------------------

ExRootSTDHEPReader::ExRootSTDHEPReader() :
  fInputFile(0), fInputXDR(0), fBuffer(0), fBlockType(-1)
{
  fInputXDR = new XDR;
  fBuffer = new char[kBufferSize*96 + 24];
}

//---------------------------------------------------------------------------

ExRootSTDHEPReader::~ExRootSTDHEPReader()
{
  if(fBuffer) delete fBuffer;
  if(fInputXDR) delete fInputXDR;
}

//---------------------------------------------------------------------------

void ExRootSTDHEPReader::SetInputFile(FILE *inputFile)
{
  fInputFile = inputFile;
  xdrstdio_create(fInputXDR, inputFile, XDR_DECODE);
  ReadFileHeader();
}

//---------------------------------------------------------------------------

void ExRootSTDHEPReader::Clear()
{
  fBlockType = -1;
}

//---------------------------------------------------------------------------

bool ExRootSTDHEPReader::EventReady()
{
  return (fBlockType == MCFIO_STDHEP) || (fBlockType == MCFIO_STDHEP4);
}

//---------------------------------------------------------------------------

bool ExRootSTDHEPReader::ReadBlock(ExRootTreeBranch *branch)
{
  if(feof(fInputFile)) return kFALSE;

  xdr_int(fInputXDR, &fBlockType);

  SkipBytes(4);

  if(fBlockType == EVENTTABLE)
  {
    ReadEventTable();
  }
  else if(fBlockType == EVENTHEADER)
  {
    ReadEventHeader();
  }
  else if(fBlockType == MCFIO_STDHEPBEG ||
          fBlockType == MCFIO_STDHEPEND)
  {
    ReadSTDCM1();
  }
  else if(fBlockType == MCFIO_STDHEP)
  {
    ReadSTDHEP();
    AnalyzeParticles(branch);
  }
  else if(fBlockType == MCFIO_STDHEP4)
  {
    ReadSTDHEP();
    AnalyzeParticles(branch);
    ReadSTDHEP4();
  }
  else
  {
    throw runtime_error("Unsupported block type.");
  }

  return kTRUE;
}

//---------------------------------------------------------------------------

void ExRootSTDHEPReader::SkipBytes(u_int size)
{
  int rc;
  u_int rndup;

  rndup = size % 4;
  if(rndup > 0)
  {
    rndup = 4 - rndup;
  }

  rc = fseek(fInputFile, size + rndup, SEEK_CUR);

  if(rc != 0 && errno == ESPIPE)
  {
    xdr_opaque(fInputXDR, fBuffer, size);
  }
}

//---------------------------------------------------------------------------

void ExRootSTDHEPReader::SkipArray(u_int elsize)
{
  u_int size;
  xdr_u_int(fInputXDR, &size);
  SkipBytes(size*elsize);
}

//---------------------------------------------------------------------------

void ExRootSTDHEPReader::ReadFileHeader()
{
  u_int i;
  enum STDHEPVersion {UNKNOWN, V1, V2, V21} version;

  xdr_int(fInputXDR, &fBlockType);
  if (fBlockType != FILEHEADER)
  {
    throw runtime_error("Header block not found. File is probably corrupted.");
  }

  SkipBytes(4);

  // version
  xdr_string(fInputXDR, &fBuffer, 100);
  if(fBuffer[0] == '\0' || fBuffer[1] == '\0') version = UNKNOWN;
  else if(fBuffer[0] == '1') version = V1;
  else if(strncmp(fBuffer, "2.01", 4) == 0) version = V21;
  else if(fBuffer[0] == '2') version = V2;
  else version = UNKNOWN;

  if(version == UNKNOWN)
  {
    throw runtime_error("Unknown file format version.");
  }

  SkipArray(1);
  SkipArray(1);
  SkipArray(1);

  if(version == V21)
  {
    SkipArray(1);
  }

  // Expected number of events
  SkipBytes(4);

  // Number of events
  xdr_u_int(fInputXDR, &fEntries);

  SkipBytes(8);

  // Number of blocks
  u_int nBlocks = 0;
  xdr_u_int(fInputXDR, &nBlocks);

  // Number of NTuples
  u_int nNTuples = 0;
  if(version != V1)
  {
    xdr_u_int(fInputXDR, &nNTuples);
  }

  if(nNTuples != 0)
  {
    throw runtime_error("Files containing n-tuples are not supported.");
  }

  // Processing blocks extraction
  if(nBlocks != 0)
  {
    SkipArray(4);

    for(i = 0; i < nBlocks; i++)
    {
      SkipArray(1);
    }
  }
}

//---------------------------------------------------------------------------

void ExRootSTDHEPReader::ReadEventTable()
{
  // version
  xdr_string(fInputXDR, &fBuffer, 100);
  if(strncmp(fBuffer, "1.00", 4) == 0)
  {
    SkipBytes(8);

    SkipArray(4);
    SkipArray(4);
    SkipArray(4);
    SkipArray(4);
    SkipArray(4);
  }
  else if(strncmp(fBuffer, "2.00", 4) == 0)
  {
    SkipBytes(12);

    SkipArray(4);
    SkipArray(4);
    SkipArray(4);
    SkipArray(4);
    SkipArray(8);
  }
}

//---------------------------------------------------------------------------

void ExRootSTDHEPReader::ReadEventHeader()
{
  bool skipNTuples = false;
  u_int skipSize = 4;

  // version
  xdr_string(fInputXDR, &fBuffer, 100);
  if(strncmp(fBuffer, "2.00", 4) == 0)
  {
    skipNTuples = true;
  }
  else if(strncmp(fBuffer, "3.00", 4) == 0)
  {
    skipNTuples = true;
    skipSize = 8;
  }

  SkipBytes(20);

  u_int dimBlocks = 0;
  xdr_u_int(fInputXDR, &dimBlocks);

  u_int dimNTuples = 0;
  if(skipNTuples)
  {
    SkipBytes(4);
    xdr_u_int(fInputXDR, &dimNTuples);
  }

  // Processing blocks extraction
  if(dimBlocks > 0)
  {
    SkipArray(4);
    SkipArray(skipSize);
  }

  // Processing blocks extraction
  if(skipNTuples && dimNTuples > 0)
  {
    SkipArray(4);
    SkipArray(skipSize);
  }
}

//---------------------------------------------------------------------------

void ExRootSTDHEPReader::ReadSTDCM1()
{
  // version
  xdr_string(fInputXDR, &fBuffer, 100);

  // skip 5*4 + 2*8 = 36 bytes
  SkipBytes(36);

  if((strncmp(fBuffer, "1.", 2) == 0) || (strncmp(fBuffer, "2.", 2) == 0) ||
     (strncmp(fBuffer, "3.", 2) == 0) || (strncmp(fBuffer, "4.", 2) == 0) ||
     (strncmp(fBuffer, "5.00", 4) == 0))
  {
    return;
  }

  SkipArray(1);
  SkipArray(1);

  if(strncmp(fBuffer, "5.01", 4) == 0)
  {
    return;
  }

  SkipBytes(4);
}

//---------------------------------------------------------------------------

void ExRootSTDHEPReader::ReadSTDHEP()
{
  u_int idhepSize, isthepSize, jmohepSize, jdahepSize, phepSize, vhepSize;

  // version
  xdr_string(fInputXDR, &fBuffer, 100);

  // Extracting the event number
  xdr_int(fInputXDR, &fEventNumber);

  // Extracting the number of particles
  xdr_int(fInputXDR, &fEventSize);

  if(fEventSize >= kBufferSize)
  {
    throw runtime_error("too many particles in event");
  }

  // 4*n + 4*n + 8*n + 8*n + 40*n + 32*n +
  // 4 + 4 + 4 + 4 + 4 + 4 = 96*n + 24

  xdr_opaque(fInputXDR, fBuffer, 96*fEventSize + 24);

  idhepSize = ntohl(*(u_int*)(fBuffer));
  isthepSize = ntohl(*(u_int*)(fBuffer + 4*1 + 4*1*fEventSize));
  jmohepSize = ntohl(*(u_int*)(fBuffer + 4*2 + 4*2*fEventSize));
  jdahepSize = ntohl(*(u_int*)(fBuffer + 4*3 + 4*4*fEventSize));
  phepSize = ntohl(*(u_int*)(fBuffer + 4*4 + 4*6*fEventSize));
  vhepSize = ntohl(*(u_int*)(fBuffer + 4*5 + 4*16*fEventSize));

  if(fEventSize < 0 ||
     fEventSize != (int)idhepSize      || fEventSize != (int)isthepSize     ||
     (2*fEventSize) != (int)jmohepSize || (2*fEventSize) != (int)jdahepSize ||
     (5*fEventSize) != (int)phepSize   || (4*fEventSize) != (int)vhepSize)
  {
    throw runtime_error("Inconsistent size of arrays. File is probably corrupted.");
  }

  fWeight = 1.0;
  fAlphaQED = 0.0;
  fAlphaQCD = 0.0;
  fScaleSize = 0;
  memset(fScale, 0, 10*sizeof(double));
}

//---------------------------------------------------------------------------

void ExRootSTDHEPReader::ReadSTDHEP4()
{
  u_int number;

  // Extracting the event weight
  xdr_double(fInputXDR, &fWeight);

  // Extracting alpha QED
  xdr_double(fInputXDR, &fAlphaQED);

  // Extracting alpha QCD
  xdr_double(fInputXDR, &fAlphaQCD);

  // Extracting the event scale
  xdr_u_int(fInputXDR, &fScaleSize);
  for(number = 0; number < fScaleSize; ++number)
  {
    xdr_double(fInputXDR, &fScale[number]);
  }

  SkipArray(8);
  SkipArray(4);

  SkipBytes(4);
}

//---------------------------------------------------------------------------

void ExRootSTDHEPReader::AnalyzeEvent(ExRootTreeBranch *branch, long long eventNumber)
{
  TRootLHEFEvent *element;

  element = static_cast<TRootLHEFEvent *>(branch->NewEntry());

  element->Number = fEventNumber;

  element->ProcessID = 0;

  element->Weight = fWeight;
  element->ScalePDF = fScale[0];
  element->CouplingQED = fAlphaQED;
  element->CouplingQCD = fAlphaQCD;
}

//---------------------------------------------------------------------------

void ExRootSTDHEPReader::AnalyzeParticles(ExRootTreeBranch *branch)
{
  TRootGenParticle *element;

  Double_t signPz, cosTheta;
  TLorentzVector momentum;

  int number;
  int pid, status, m1, m2, d1, d2;
  double px, py, pz, e, mass;
  double x, y, z, t;

  XDR bufferXDR[6];
  xdrmem_create(&bufferXDR[0], fBuffer + 4*1, 4*fEventSize, XDR_DECODE);
  xdrmem_create(&bufferXDR[1], fBuffer + 4*2 + 4*1*fEventSize, 4*fEventSize, XDR_DECODE);
  xdrmem_create(&bufferXDR[2], fBuffer + 4*3 + 4*2*fEventSize, 8*fEventSize, XDR_DECODE);
  xdrmem_create(&bufferXDR[3], fBuffer + 4*4 + 4*4*fEventSize, 8*fEventSize, XDR_DECODE);
  xdrmem_create(&bufferXDR[4], fBuffer + 4*5 + 4*6*fEventSize, 40*fEventSize, XDR_DECODE);
  xdrmem_create(&bufferXDR[5], fBuffer + 4*6 + 4*16*fEventSize, 32*fEventSize, XDR_DECODE);

  for(number = 0; number < fEventSize; ++number)
  {
    xdr_int(&bufferXDR[0], &status);
    xdr_int(&bufferXDR[1], &pid);
    xdr_int(&bufferXDR[2], &m1);
    xdr_int(&bufferXDR[2], &m2);
    xdr_int(&bufferXDR[3], &d1);
    xdr_int(&bufferXDR[3], &d2);

    xdr_double(&bufferXDR[4], &px);
    xdr_double(&bufferXDR[4], &py);
    xdr_double(&bufferXDR[4], &pz);
    xdr_double(&bufferXDR[4], &e);
    xdr_double(&bufferXDR[4], &mass);

    xdr_double(&bufferXDR[5], &x);
    xdr_double(&bufferXDR[5], &y);
    xdr_double(&bufferXDR[5], &z);
    xdr_double(&bufferXDR[5], &t);

    element = static_cast<TRootGenParticle*>(branch->NewEntry());

    element->PID = pid;
    element->Status = status;

    element->M1 = m1 - 1;
    element->M2 = m2 - 1;

    element->D1 = d1 - 1;
    element->D2 = d2 - 1;

    element->E = e;
    element->Px = px;
    element->Py = py;
    element->Pz = pz;

    momentum.SetPxPyPzE(px, py, pz, e);

    cosTheta = TMath::Abs(momentum.CosTheta());
    signPz = (momentum.Pz() >= 0.0) ? 1.0 : -1.0;

    element->PT = momentum.Perp();
    element->Phi = momentum.Phi();
    element->Eta = (cosTheta == 1.0 ? signPz*999.9 : momentum.Eta());
    element->Rapidity = (cosTheta == 1.0 ? signPz*999.9 : momentum.Rapidity());

    element->T = t;
    element->X = x;
    element->Y = y;
    element->Z = z;
  }
}

//---------------------------------------------------------------------------
