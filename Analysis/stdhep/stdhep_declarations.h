#ifndef STDHEP_DECLARATIONS_H
#define STDHEP_DECLARATIONS_H

/* declare the struct instances */

#include "stdhep.h"
#include "hepev4.h"
#include "hepeup.h"
#include "heprup.h"
#include "stdtmp.h"
#include "stdhd.h"
#include "stdcnt.h"
#include "stdcm1.h"
#include "stdver.h"

#if defined(c_plusplus) || defined(__cplusplus)
extern "C" {
#endif

struct hepevt myhepevt;
struct hepev2 hepev2_;
struct hepev3 hepev3_;
struct hepev4 hepev4_;
struct hepev5 hepev5_;
struct hepeup hepeup_;
struct heprup heprup_;
struct stdcnt stdcnt_;
struct stdhd1 stdhd1_;
struct stdhd2 stdhd2_;
struct stdtmp stdtmp_;
struct tmpev4 tmpev4_;
struct stdcm1 stdcm1_;
struct stdcm2 stdcm2_;
struct stdver stdver_;

#if defined(c_plusplus) || defined(__cplusplus)
}
#endif

#endif /* STDHEP_DECLARATIONS_H */
