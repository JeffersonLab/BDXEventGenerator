/*******************************************************************************
*									       *
* stdhep_internal_utils.c -- C version of stdhep internal utility routines     *
*									       *
* Copyright (c) 1995 Universities Research Association, Inc.		       *
* All rights reserved.							       *
* 									       *
* This material resulted from work developed under a Government Contract and   *
* is subject to the following license:  The Government retains a paid-up,      *
* nonexclusive, irrevocable worldwide license to reproduce, prepare derivative *
* works, perform publicly and display publicly by or for the Government,       *
* including the right to distribute to other Government contractors.  Neither  *
* the United States nor the United States Department of Energy, nor any of     *
* their employees, makes any warranty, express or implied, or assumes any      *
* legal liability or responsibility for the accuracy, completeness, or         *
* usefulness of any information, apparatus, product, or process disclosed, or  *
* represents that its use would not infringe privately owned rights.           *
*                                        				       *
*									       *
* Written by Lynn Garren    					       	       *
*									       *
*									       *
*******************************************************************************/
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#ifndef FALSE
#define FALSE 0
#endif
#ifndef TRUE
#define TRUE 1
#endif
/* 
*   StdHep definitions and include files
*/
#include "stdhep.h"
#include "hepev4.h"
#include "stdtmp.h"
#include "hepeup.h"

extern struct stdtmp stdtmp_;
extern struct tmpev4 tmpev4_;

extern struct hepevt myhepevt;
extern struct hepev2 hepev2_;
extern struct hepev3 hepev3_;
extern struct hepev4 hepev4_;
extern struct hepev5 hepev5_;
extern struct hepeup hepeup_;

/* Purpose: copy an event to/from the standard common block */
int StdHepTempCopy(int idir, int istr)
{
    int nh, i, k;
    if (idir == 1) {        /* copy from hepevt to stdtmp */
        stdtmp_.nevhept = myhepevt.nevhep;
        stdtmp_.nhept = myhepevt.nhep;
	tmpev4_.eventweightt = hepev4_.eventweightlh;
	tmpev4_.alphaqedt = hepev4_.alphaqedlh;
	tmpev4_.alphaqcdt = hepev4_.alphaqcdlh;
        for (i = 0; i < 10; i++) {
	    tmpev4_.scalet[i] = hepev4_.scalelh[i];
	}
	tmpev4_.idrupt = hepev4_.idruplh;
        for (i = 0; i < myhepevt.nhep; i++) {
            stdtmp_.isthept[i] = myhepevt.isthep[i];
            stdtmp_.idhept[i] = myhepevt.idhep[i];
            for (k = 0; k < 2; k++) {
                stdtmp_.jmohept[i][k] = myhepevt.jmohep[i][k];
                stdtmp_.jdahept[i][k] = myhepevt.jdahep[i][k];
		tmpev4_.icolorflowt[i][k] = hepev4_.icolorflowlh[i][k];
                }
            for (k = 0; k < 5; k++)
                stdtmp_.phept[i][k] = myhepevt.phep[i][k];
            for (k = 0; k < 4; k++)
                stdtmp_.vhept[i][k] = myhepevt.vhep[i][k];
            for (k = 0; k < 3; k++)
                tmpev4_.spint[i][k] = hepev4_.spinlh[i][k];
            }
    } else if (idir == 2) {    /* copy from stdtmp to hepevt */
        if (myhepevt.nhep + stdtmp_.nhept > NMXHEP) {
            fprintf(stderr,
              "     StdHepTempCopy: event would overflow HEPEVT array size\n");
            fprintf(stderr,"     StdHepTempCopy: event %d has been lost\n",
                   stdtmp_.nevhept);
            return 5;
            }
        myhepevt.nevhep = stdtmp_.nevhept;
        nh = myhepevt.nhep;
	hepev4_.eventweightlh = tmpev4_.eventweightt;
	hepev4_.alphaqedlh = tmpev4_.alphaqedt;
	hepev4_.alphaqcdlh = tmpev4_.alphaqcdt;
        for (i = 0; i < 10; i++) {
	    hepev4_.scalelh[i] = tmpev4_.scalet[i];
	}
	hepev4_.idruplh = tmpev4_.idrupt;
        for (i = 0; i < stdtmp_.nhept; i++) {
            myhepevt.isthep[nh+i] = stdtmp_.isthept[i];
            myhepevt.idhep[nh+i] = stdtmp_.idhept[i];
            for (k = 0; k < 2; k++) {
                myhepevt.jmohep[nh+i][k] = stdtmp_.jmohept[i][k];
                myhepevt.jdahep[nh+i][k] = stdtmp_.jdahept[i][k];
		hepev4_.icolorflowlh[nh+i][k] = tmpev4_.icolorflowt[i][k];
                }
            for (k = 0; k < 5; k++)
                myhepevt.phep[nh+i][k] = stdtmp_.phept[i][k];
            for (k = 0; k < 4; k++)
                myhepevt.vhep[nh+i][k] = stdtmp_.vhept[i][k];
            for (k = 0; k < 3; k++)
                hepev4_.spinlh[nh+i][k] = tmpev4_.spint[i][k];
            }
        hepev2_.nmulti += 1;
	if (hepev2_.nmulti <= NMXMLT ) {
	    hepev3_.nevmulti[hepev2_.nmulti] = stdtmp_.nevhept;
	    hepev3_.itrkmulti[hepev2_.nmulti] = stdtmp_.nhept + 1;
	    hepev3_.mltstr[hepev2_.nmulti] = istr;
	    hepev5_.eventweightmulti[i] = tmpev4_.eventweightt;
	    hepev5_.alphaqedmulti[i] = tmpev4_.alphaqedt;
	    hepev5_.alphaqcdmulti[i] = tmpev4_.alphaqcdt;
	    for( k = 0; k < 10; ++k) {
		hepev5_.scalemulti[i][k] = tmpev4_.scalet[k];
	    }
	    hepev5_.idrupmulti[i] = tmpev4_.idrupt;
	} else {
	    fprintf(stderr," StdHepTempCopy: %d multiple interactions in this event\n",
	         hepev2_.nmulti );  
	    fprintf(stderr," StdHepTempCopy: only %d multiple interactions are allowed\n",
	         NMXMLT );  
	}
        for (i = 0; i < stdtmp_.nhept; i++) {
            hepev2_.jmulti[nh+i] = hepev2_.nmulti;
            for (k = 0; k < 2; k++) {
                if (myhepevt.jmohep[nh+i][k] != 0) {
		   myhepevt.jmohep[nh+i][k] += myhepevt.nhep;
		   }
                if (myhepevt.jdahep[nh+i][k] != 0) {
		   myhepevt.jdahep[nh+i][k] += myhepevt.nhep;
		   }
                if (hepev4_.icolorflowlh[nh+i][k] != 0) {
		   hepev4_.icolorflowlh[nh+i][k] += myhepevt.nhep;
		   }
                }
        }
        myhepevt.nhep += stdtmp_.nhept;
    } else {
        fprintf(stderr," StdHepTempCopy: improper calling flag\n");
    }
    return 0;
}

void StdHepZero(void)
{
    int i, k;
    myhepevt.nhep = 0;
    hepev2_.nmulti = 0;
    for (i = 0; i < NMXHEP; i++) {
        myhepevt.isthep[i] = 0;
        myhepevt.idhep[i] = 0;
        hepev2_.jmulti[i] = 0;
        for (k = 0; k < 2; k++) {
            myhepevt.jmohep[i][k] = 0;
            myhepevt.jdahep[i][k] = 0;
	    hepev4_.icolorflowlh[i][k] = 0;
            }
        for (k = 0; k < 5; k++)
            myhepevt.phep[i][k] = 0.;
        for (k = 0; k < 4; k++)
            myhepevt.vhep[i][k] = 0.;
        for (k = 0; k < 3; k++)
            hepev4_.spinlh[i][k] = 0.;
        }
    for (i = 0; i < NMXMLT; i++) {
        hepev3_.nevmulti[i] = 0;
        hepev3_.itrkmulti[i] = 0;
        hepev3_.mltstr[i] = 0;
	hepev5_.eventweightmulti[i] = 0.;
	hepev5_.alphaqedmulti[i] = 0.;
	hepev5_.alphaqcdmulti[i] = 0.;
	for( k = 0; k < 10; ++k) {
	    hepev5_.scalemulti[i][k] = 0.;
	}
	hepev5_.idrupmulti[i] = 0;
    }
    hepev4_.eventweightlh = 0.;
    hepev4_.alphaqedlh = 0.;
    hepev4_.alphaqcdlh = 0.;
    for (i = 0; i < 10; i++) {
        hepev4_.scalelh[i] = 0.;
    }
    hepev4_.idruplh = 0;
}

void StdHepZeroHEPEUP(void)
{
   int i, k;
   hepeup_.nup = 0;
   hepeup_.idprup = 0;
   hepeup_.xwgtup = 0;
   hepeup_.scalup = 0;
   hepeup_.aqedup = 0;
   hepeup_.aqcdup = 0;
   for (i = 0; i < MAXNUP; ++i) {
      hepeup_.idup[i] = 0;
      hepeup_.istup[i] = 0;
      for (k = 0; k < 2; ++k) {
         hepeup_.mothup[i][k] = 0;
         hepeup_.icolup[i][k] = 0;
      }
      for (k = 0; k < 5; ++k) {
         hepeup_.pup[i][k] = 0;
      }
      hepeup_.vtimup[i] = 0;
      hepeup_.spinup[i] = 0;
   }
}
