#ifndef STDHEP_H
#define STDHEP_H
/*
** Basic COMMON block from STDHEP: the HEPEVT COMMON block 
** See product StDhep
*/
/*  note that to avoid alignment problems, structures and common blocks
    should be in the order: double precision, real, integer.
*/

#if defined(c_plusplus) || defined(__cplusplus)
extern "C" {
#endif

#define NMXHEP 4000
#define NMXMLT 16
struct hepevt {
int nevhep;		/* The event number */
int nhep;		/* The number of entries in this event */
int isthep[NMXHEP]; 	/* The Particle id */
int idhep[NMXHEP];      /* The particle id */
int jmohep[NMXHEP][2];    /* The position of the mother particle */
int jdahep[NMXHEP][2];    /* Position of the first daughter... */
double phep[NMXHEP][5];    /* 4-Momentum, mass */
double vhep[NMXHEP][4];    /* Vertex information */
};
struct hepev2 {
int nmulti;		/* number of interactions in the list */
int jmulti[NMXHEP];     /* multiple interaction number */
};
struct hepev3 {
int nevmulti[NMXMLT];     /* event number of original interaction */
int itrkmulti[NMXMLT];     /* first particle in the original interaction */
int mltstr[NMXMLT];     /* stream this event is from */
};

#if defined(c_plusplus) || defined(__cplusplus)
}
#endif

#endif /* STDHEP_H */

