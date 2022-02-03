
#ifndef DATATYPES_H
#define DATATYPES_H


#include "double2.h"
#include "double3.h"
#include "double4.h"
#include "double5.h"
#include "float2.h"
#include "float3.h"
#include "float4.h"
#include "float5.h"
#include "int2.h"
#include "int3.h"
#include "int5.h"

static const int MAXMAT = 10;
 static const int NS = 10000;
static const int NB = 5000;
static const int NXG = 250;
static const int NS2M = 2 * NS;

typedef struct {
    double E, X, Y, Z, U, V, W, WGFT, SP1, SP2, SP3, PAGE;
   // bool lage = false;
    bool LAGE;
    int KPAR, IBODY, MAT, ILB[5], IPOL;
} TRACK_MOD;

typedef struct {
    double EABS[3][MAXMAT], C1[MAXMAT], C2[MAXMAT], WCC[MAXMAT], WCR[MAXMAT], DEN[MAXMAT], RDEN[MAXMAT], E0STEP, DESOFT, SSOFT;
   // int NMS = 1000;
//	int NEGP = 200;
   int NMS;
   int NEGP;
   int NMAT;
} PENELOPE_MOD;

//PENELOPE_MOD penelope_mod;


 typedef struct {
    char BALIAS[NB];
  //  bool lverb = false;
    bool LVERB;
    double dstot;
    int MATER[NB], KDET[NB], KSLAST, NBODY;
} PENGEOM_MOD;

//PENGEOM_MOD pengeom_mod;



//TRACK_MOD track_mod;


 typedef struct {
    double AXX[NS], AXY[NS], AXZ[NS], AYY[NS], AYZ[NS],
        AZZ[NS], AX[NS], AY[NS], AZ[NS], A0[NS];
    int NSURF, KPLANE[NS];
} QSURF;


typedef struct{
    int NBODYS, KMOTH[NB], KDGHT[NB][NXG], KSURF[NB][NXG], KFLAG[NB][NXG], KSP[NS], NWARN;
} QTREE;


typedef struct {
    int KBODY[NB][NXG];
    int KBOMO[NB];

} QBODY;


typedef struct {
    int s[NS2M];
    int is[NS2M];
    int mat0;
    double dsres;
    double dsp;
    int ks;
    int kb1;
    int ks1;
    int kf;
    int kflo;
    int ierr;
    int nerr;
    int sw;

   // int nsc = 0; // numero de cruzamentos da superficie à frente da particula
  //  int nsct = 0;
  //  int nst = 0;
    
	int nsc; // numero de cruzamentos da superficie à frente da particula
    int nsct;
    int nst;
    int matl;
    int ibodyl;
} STEP_MOD;

/*typedef struct {
    int kbody[NB][NXG];
    int KBOMO[NB];

} QBODY;*/




#endif // DATATYPES_H
