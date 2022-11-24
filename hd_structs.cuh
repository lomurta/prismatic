static const int MAXMAT = 10;
static const int NS = 10000;
static const int NB = 5000;
static const int NXG = 250;
static const int NS2M = 2 * NS;
static const int NEGP = 200;
static const int NMS = 1000;
static const int NRP = 8000;
static const int NTP = 12000;
static const int NRX = 60000;
static const int NM = 512;
static const int NR = 128;
static const double ZERO = 1.0e-75;
static const double ZEROT = 0.1 * ZERO;
static const double A0B = 5.2917721092e-9; //Bohr radius (cm)
static const double HREV = 27.21138505e0; //Hartree energy (eV)
static const double AVOG = 6.02214129e23; //Avogadro's number
static const double SL = 137.035999074e0; //Speed of light (1/alpha)
static const double PI = 3.1415926535897932e0;
static const double FOURPI = 4.0e0 * PI;
static const int NO = 512;
static const int NOCO = 512;
static const int NDIM = 12000;
static const double HBAR = 6.58211928e-16; //Constante de Planck eV*s
static const int NBE = 57;
static const int NBW = 32;
static const int NE = 96;
static const int NA = 606;
static const int NP = 128;
static const int NQ = 250;
static const int NEX = 1024;
static const int NP2 = 150;
static const int NSEM = 1000;
static const int NPSFM = 100;
static const int NIDM = 25;
static const int NBV = 5000;
static const int NBEM = 1500;
static const int NBTHM = 1800;
static const int NBPHM = 180;
static const int NBEM2 = 1000;
static const int NDXM = 201;
static const int NDYM = 201;
static const int NDZM = 201;

int panar2 = 0;

static const int pilhaPart = 3072; //64*64
static const int pilhaSec = 3072; //64*64

#define MAX_THREADS_PER_BLOCK 1024
#define MIN_BLOCKS_PER_MP 16

char LINHA[200];
char APOIO[200];
char SPCDIO[NIDM][20];
char SPCFLO[NIDM][20];
char SPCAGE[NIDM][20];
char SPCDEO[NIDM][20];
double S[NS2M];
int IS[NS2M];
int imprimiu = 0;
int wIPOLI = 0;
int h_N = 0; //mesma funcionalidade do CNTRL_.N
int h_vetN[pilhaPart];

static const int tamISEEDs = 10001;

int IS1[tamISEEDs];
int IS2[tamISEEDs];

int const blockSize = 128;
int tamMemShared = 8*blockSize*7;


__constant__ double d_PI = 3.1415926535897932e0;
__constant__ double d_REV = 5.10998928e5;
__constant__ double d_TREV = 2.0e0 * 5.10998928e5;
__constant__ double d_RREV = 1.0e0 / 5.10998928e5;
__constant__ double d_TWOPI = 2.0e0 * 3.1415926535897932e0;
__constant__ double d_RTREV = 1.0e0 / (2.0e0 * 5.10998928e5);

__constant__ double d_FUZZL = 1.0e-12;

__constant__ int d_NS2M = 2 * NS;

__constant__ int ROLL = 4;


clock_t start, end;

typedef struct
{
	double X[blockSize],Y[blockSize],Z[blockSize],U[blockSize],V[blockSize], W[blockSize];
} hd_TRACK_MOD_SHARED;


typedef struct
{
	int IEXIT[pilhaPart], KEn[pilhaPart], IBODYL[pilhaPart], NCROSS[pilhaPart], IDET[pilhaPart], MATL[pilhaPart], ICOL[pilhaPart], LEFT[pilhaPart];
	double DSEF[pilhaPart], DS[pilhaPart], DSMAX[pilhaPart], DEP[pilhaPart], XL[pilhaPart], YL[pilhaPart], ZL[pilhaPart], DECSD[pilhaPart], DSEFR[pilhaPart], XD[pilhaPart], YD[pilhaPart], ZD[pilhaPart], DE[pilhaPart], WS[pilhaPart], US[pilhaPart], VS[pilhaPart], SDTS[pilhaPart], DF[pilhaPart];
	bool LINTF[pilhaPart], CROSS[pilhaPart];
} hd_wSHOWERS;



typedef	struct {
	double* AXX, * AXY, * AXZ, * AYY, * AYZ,
		* AZZ, * AX, * AY, * AZ, * A0;
	int* NSURF, * KPLANE;
} QSURF;

typedef	struct {
	double AXX[NS],  AXY[NS],  AXZ[NS],  AYY[NS], AYZ[NS],
	       AZZ[NS],  AX[NS], AY[NS], AZ[NS], A0[NS];
	int NSURF, KPLANE[NS];
} hd_QSURF;

typedef struct {
	int* NBODYS, * KMOTH, (*KDGHT)[NB], (*KSURF)[NB], (*KFLAG)[NB], * KSP, * NWARN;
} QTREE;

typedef struct {
	int NBODYS, KMOTH[NB], KDGHT[NXG][NB], KSURF[NXG][NB], KFLAG[NXG][NB], KSP[NS][pilhaPart], NWARN[pilhaPart];
} hd_QTREE;

typedef struct {
	double* E, * X, * Y, * Z, * U, * V, * W, * WGHT, * SP1, * SP2, * SP3, * PAGE;
	int* KPAR, * IBODY, * MAT, * ILB, * IPOL;
	bool* LAGE;
} TRACK_MOD;

typedef struct {
	double X,  Y,  Z,  U,  V,  W;
	int IBODY,  MAT;
} preTRACK_MOD;

typedef struct {
	double E[pilhaPart], 
	X[pilhaPart], Y[pilhaPart], Z[pilhaPart], 
	U[pilhaPart], V[pilhaPart], W[pilhaPart], WGHT[pilhaPart], 							SP1[pilhaPart], SP2[pilhaPart], SP3[pilhaPart], PAGE[pilhaPart];
	int KPAR[pilhaPart], IBODY[pilhaPart], MAT[pilhaPart], ILB[5][pilhaPart], 			IPOL[pilhaPart], INDEX[pilhaPart], N[pilhaPart], STEP[pilhaPart], IEXIT[pilhaPart];
																					bool LAGE[pilhaPart];
} hd_TRACK_MOD;

typedef struct {
	double E[pilhaPart], 
		   X[pilhaPart], Y[pilhaPart], Z[pilhaPart], 
		   U[pilhaPart], V[pilhaPart], W[pilhaPart], 
		   WGHT[pilhaPart], 							
	int KPAR[pilhaPart], IBODY[pilhaPart],
		MAT[pilhaPart], ILB[5][pilhaPart], 																							
} TRACK_MOD;



static const int sec = 100;

typedef struct {
	double E[pilhaPart*sec], X[pilhaPart*sec], Y[pilhaPart*sec], Z[pilhaPart*sec], U[pilhaPart*sec], V[pilhaPart*sec], W[pilhaPart*sec], WGHT[pilhaPart*sec], SP1[pilhaPart*sec], SP2[pilhaPart*sec], SP3[pilhaPart*sec], PAGE[pilhaPart*sec];
	int KPAR[pilhaPart*sec], IBODY[pilhaPart*sec], MAT[pilhaPart*sec], ILB[5][pilhaPart*sec], IPOL[pilhaPart*sec], INDEX[pilhaPart*sec], N[pilhaPart*sec], STEP[pilhaPart*sec], IEXIT[pilhaPart*sec];
	bool LAGE[pilhaPart*sec];
} hd_TRACK_MOD_SEC;


typedef struct {
	double(*EABS)[3], * C1, * C2, * WCC, * WCR, * DEN, * RDEN, * E0STEP, * DESOFT, * SSOFT;
	int* NMS, * NEGP, * NMAT;
}PENELOPE_MOD;

typedef struct {
	double EABS[MAXMAT][3], C1[MAXMAT], C2[MAXMAT], WCC[MAXMAT], WCR[MAXMAT], DEN[MAXMAT], RDEN[MAXMAT], E0STEP[pilhaPart], DESOFT[pilhaPart], SSOFT[pilhaPart];
	int NMS, NEGP, NMAT;
}  hd_PENELOPE_MOD;


typedef struct {
	char(*BALIAS)[5];
	double* DSTOT;
	int* MATER, * KDET, * KSLAST, * NBODY;
	bool* LVERB;
} PENGEOM_MOD;

typedef struct {
	int NBODY, MATER[NB], KDET[NB], KSLAST[pilhaPart];
	double DSTOT[pilhaPart];
	char BALIAS[NB][5];
	bool LVERB;
} hd_PENGEOM_MOD;

typedef struct {
	int(*KBODY)[NB], * KBOMO;
} QBODY;

typedef struct {
	int KBODY[NXG][NB], KBOMO[NB];
} hd_QBODY;

typedef struct {
	double* ECUTR;
} CECUTR; 

typedef struct {
	double ECUTR[MAXMAT];
} hd_CECUTR; 

typedef struct {
	int* ISGAW;
} CSGAWR;

typedef struct {
	int* IERSEC;
} CERSEC;

typedef struct {
	int IERSEC;
} hd_CERSEC;

typedef struct {
	double* EMIN, * EL, * EU, * ET, * DLEMP, * DLEMP1, * DLFC, * XEL, * XE, * XEK;
	int* KE;
} CEGRID; //Rede de energia e constantes de interpolacao para a energia atual.

typedef struct {
	double EMIN, EL, EU, ET[NEGP], DLEMP[NEGP], DLEMP1, DLFC, XEL[pilhaPart], XE[pilhaPart], XEK[pilhaPart];
	int KE[pilhaPart];
} hd_CEGRID;


typedef struct {
	double(*XESI)[NRP];
	int* IESIF, * IESIL, * NSESI, * NCURE;
} CESI0; //Ioniza��o da camada interna por impacto de el�trons e p�sitrons.


typedef struct {
	double(*XPSI)[NRP];
	int* IPSIF, * IPSIL, * NSPSI, * NCURP;
} CPSI0; //Ioniza��o da camada interna por impacto de el�trons e p�sitrons.


typedef struct {
	double* ATW, * EPX, * RSCR, * ETA, (*EB)[99], (*ALW)[99], (*CP0)[99];
	int(*IFI)[99], (*IKS)[99], * NSHT;
	char(*LASYMB)[3];
} CADATA; //Dados elemento

typedef struct {
	double ATW[99], EPX[99], RSCR[99], ETA[99], EB[30][99], ALW[30][99], CP0[30][99];
	int IFI[30][99], IKS[30][99], NSHT[99];
	char LASYMB[99][3];
} hd_CADATA;

typedef struct {
	double* EPH, (*XPH)[NTP];
	int* IPHF, * IPHL, * NPHS, * NCUR;
} CGPH00;

typedef struct {
	double EPH[NTP], XPH[17][NTP];
	int IPHF[99], IPHL[99], NPHS[99], NCUR;
} hd_CGPH00;


typedef struct {
	double* P, * ET, * F;
	int* IS0, * IS1, * IS2, (*IFIRST)[99], (*ILAST)[99], * NCUR, * KS, * MODER;
} CRELAX;

typedef struct {
	double P[NRX], ET[NRX], F[NRX];
	int IS0[NRX], IS1[NRX], IS2[NRX], IFIRST[16][99], ILAST[16][99], NCUR, KS[pilhaPart], MODER;
} hd_CRELAX;


typedef struct {

	double* XA, * AA, * BA;
	double* FA;
	int* IA, * NPM1A;

} CRITAA;

typedef struct {
	double* X, * A, * B;
	double* F;
	int* KA, * NPM1;
} CRNDG3;


typedef struct {
	double X[NR], A[NR], B[NR], F[NR];
	int KA[NR], NPM1;
} hd_CRNDG3;


typedef struct {
	double* XT, * PAC, * DPAC, * A, * B;
	int* IL, * IU, * NPM1;

	double* QTI, * PACI, * DPACI, * AI, * BI;
	int* ITLI, * ITUI, * NPM1I;

	double* XTI;

} CRITA;

typedef struct {

	double* CNORM;

} CRITAN;


typedef struct {
	double(*STF)[MAXMAT], * ZT, * AT, * RHO, * VMOL;
	int(*IZ)[MAXMAT], * NELEM;
} COMPOS; //Dados Composicao

typedef struct {
	double STF[30][MAXMAT], ZT[MAXMAT], AT[MAXMAT], RHO[MAXMAT], VMOL[MAXMAT];
	int IZ[30][MAXMAT], NELEM[MAXMAT];
} hd_COMPOS;

typedef struct {

	double(*RANGE)[MAXMAT][3], (*RANGEL)[MAXMAT][3];

} CRANGE;

typedef struct {
	double* EXPOT, * OP2, (*F)[MAXMAT], (*UI)[MAXMAT], (*WRI)[MAXMAT];
	int(*KZ)[MAXMAT], (*KS)[MAXMAT], * NOSC;
}CEIN; //Colisoes Inelasticas

typedef struct {
	double EXPOT[MAXMAT], OP2[MAXMAT], F[NO][MAXMAT], UI[NO][MAXMAT], WRI[NO][MAXMAT];
	int KZ[NO][MAXMAT], KS[NO][MAXMAT], NOSC[MAXMAT];
}hd_CEIN;

typedef struct {
	double* T1EI, * T2EI, * T1PI, * T2PI;

} CEINTF;

typedef struct {
	double* SEH0, * SEH1, * SEH2, * SES0, * SES1, * SES2, * SET0, * SET1, * SET2;
} CEIN00; //Se��es transversais parciais de conchas / osciladores individuais.

typedef struct {

	double* SPH0, * SPH1, * SPH2, * SPS0, * SPS1, * SPS2, * SPT0, * SPT1, * SPT2;

} CPIN00;

typedef struct {
	double(*EINAC)[NEGP][MAXMAT];
	int(*IEIN)[MAXMAT], * NEIN;
} CEINAC; //Inelestica de eletrons. e tabelas de ionizacao de camada interna.

typedef struct {
	double EINAC[NO][NEGP][MAXMAT];
	int IEIN[NO][MAXMAT], NEIN[MAXMAT];
} hd_CEINAC;

typedef struct {
	double(*ESIAC)[NEGP][MAXMAT];
	int(*IESI)[MAXMAT], * NESI;
} CESIAC; //Inelastica de el�trons. e tabelas de ionizazaoo de camada interna.

typedef struct {
	double ESIAC[NO][NEGP][MAXMAT];
	int IESI[NO][MAXMAT], NESI[MAXMAT];
} hd_CESIAC;

typedef struct {

	double(*XSEIN)[NEGP], (*XSESI)[NEGP];
	int* ISIE;


} CESIN; //Inel�stica de el�trons. e tabelas de ioniza��o de camada interna.


typedef struct {
	double(*PINAC)[NEGP][MAXMAT];
	int(*IPIN)[MAXMAT], * NPIN;
} CPINAC; //Positron inel�stica coll. e tabelas de ioniza��o de camada interna.

typedef struct {
	double PINAC[NO][NEGP][MAXMAT];
	int IPIN[NO][MAXMAT], NPIN[MAXMAT];
} hd_CPINAC;

typedef struct {
	double(*PSIAC)[NEGP][MAXMAT];
	int(*IPSI)[MAXMAT], * NPSI;
} CPSIAC; //Positron inel�stica coll. e tabelas de ioniza��o de camada interna.

typedef struct {
	double PSIAC[NO][NEGP][MAXMAT];
	int IPSI[NO][MAXMAT], NPSI[MAXMAT];
} hd_CPSIAC; 

typedef struct {

	double(*XSPIN)[NEGP], (*XSPSI)[NEGP];
	int* ISIP;


} CPSIN; //Positron inel�stica coll. e tabelas de ioniza��o de camada interna.

//novos
typedef struct {
	double(*FCO)[MAXMAT], (*UICO)[MAXMAT], (*FJ0)[MAXMAT], (*PTRSH)[MAXMAT];
	int(*KZCO)[MAXMAT], (*KSCO)[MAXMAT], * NOSCCO;
} CGCO; //Espalhamento Compton

typedef struct {
	double FCO[NOCO][MAXMAT], UICO[NOCO][MAXMAT], FJ0[NOCO][MAXMAT], PTRSH[NOCO][MAXMAT];
	int KZCO[NOCO][MAXMAT], KSCO[NOCO][MAXMAT], NOSCCO[MAXMAT];
} hd_CGCO;

typedef struct {
	double(*SEHEL)[MAXMAT], (*SEHIN)[MAXMAT], (*SEISI)[MAXMAT], (*SEHBR)[MAXMAT], (*SEAUX)[MAXMAT],
		(*SETOT)[MAXMAT], (*CSTPE)[MAXMAT], (*RSTPE)[MAXMAT], (*DEL)[MAXMAT], (*W1E)[MAXMAT],
		(*W2E)[MAXMAT], (*DW1EL)[MAXMAT], (*DW2EL)[MAXMAT], (*RNDCE)[MAXMAT], (*AE)[MAXMAT], (*BE)[MAXMAT],
		(*T1E)[MAXMAT], (*T2E)[MAXMAT];

} CEIMFP; //Tabela e simula��o do Eletron

typedef struct {
	double SEHEL[NEGP][MAXMAT], SEHIN[NEGP][MAXMAT], SEISI[NEGP][MAXMAT], SEHBR[NEGP][MAXMAT], SEAUX[NEGP][MAXMAT],
		SETOT[NEGP][MAXMAT], CSTPE[NEGP][MAXMAT], RSTPE[NEGP][MAXMAT], DEL[NEGP][MAXMAT], W1E[NEGP][MAXMAT],
		W2E[NEGP][MAXMAT], DW1EL[NEGP][MAXMAT], DW2EL[NEGP][MAXMAT], RNDCE[NEGP][MAXMAT], AE[NEGP][MAXMAT], BE[NEGP][MAXMAT],
		T1E[NEGP][MAXMAT], T2E[NEGP][MAXMAT];
} hd_CEIMFP;


typedef struct {
	double(*TSTPE)[MAXMAT], (*TSTRE)[MAXMAT], (*TRL1E)[MAXMAT], (*TRL2E)[MAXMAT];

} CLAS1E;

typedef struct {
	double(*SPHEL)[MAXMAT], (*SPHIN)[MAXMAT], (*SPISI)[MAXMAT], (*SPHBR)[MAXMAT], (*SPAN)[MAXMAT],
		(*SPAUX)[MAXMAT], (*SPTOT)[MAXMAT], (*CSTPP)[MAXMAT], (*RSTPP)[MAXMAT], (*W1P)[MAXMAT],
		(*W2P)[MAXMAT], (*DW1PL)[MAXMAT], (*DW2PL)[MAXMAT], (*RNDCP)[MAXMAT], (*AP)[MAXMAT], (*BP)[MAXMAT],
		(*T1P)[MAXMAT], (*T2P)[MAXMAT];
} CPIMFP; // Tabelas de simula��o dos positrons

typedef struct {
	double SPHEL[NEGP][MAXMAT], SPHIN[NEGP][MAXMAT], SPISI[NEGP][MAXMAT], SPHBR[NEGP][MAXMAT], SPAN[NEGP][MAXMAT],
		SPAUX[NEGP][MAXMAT], SPTOT[NEGP][MAXMAT], CSTPP[NEGP][MAXMAT], RSTPP[NEGP][MAXMAT], W1P[NEGP][MAXMAT],
		W2P[NEGP][MAXMAT], DW1PL[NEGP][MAXMAT], DW2PL[NEGP][MAXMAT], RNDCP[NEGP][MAXMAT], AP[NEGP][MAXMAT], BP[NEGP][MAXMAT],
		T1P[NEGP][MAXMAT], T2P[NEGP][MAXMAT];
} hd_CPIMFP;


typedef struct {
	double(*TSTPP)[MAXMAT], (*TSTRP)[MAXMAT], (*TRL1P)[MAXMAT], (*TRL2P)[MAXMAT];

} CLAS1P;

//Espalhamento elastico de eletrons e positrons

typedef struct {

	double* EJT, * XE0, * XE1, * XE2, * XP0, * XP1, * XP2, * T1E0, * T2E0, * T1P0,
		* T2P0, * EJTL, * FJL, * A, * B, * C, * D;

} CEEL00;

//Rendimentos radiativos de el�trons e p�sitrons.

typedef struct {
	double(*EBRY)[MAXMAT], (*PBRY)[MAXMAT];
} CBRYLD;

//tabelas de simula��o dos fotons

typedef struct {
	double(*SGRA)[MAXMAT], (*SGCO)[MAXMAT], (*SGPH)[MAXMAT], (*SGPP)[MAXMAT], (*SGAUX)[MAXMAT];
} CGIMFP;

typedef struct {
	double SGRA[NEGP][MAXMAT], SGCO[NEGP][MAXMAT], SGPH[NEGP][MAXMAT], SGPP[NEGP][MAXMAT], SGAUX[NEGP][MAXMAT];
} hd_CGIMFP;

typedef struct {
	double* ER, * XSR;
	int* NPHD;

} CGPH01;

typedef struct {
	double(*TRIP)[MAXMAT];
} CGPP01;

typedef struct {
	double TRIP[NEGP][MAXMAT];
} hd_CGPP01;

typedef struct {
	double* WB, (*PBCUT)[MAXMAT], (*WBCUT)[MAXMAT], (*PDFB)[NEGP][MAXMAT], (*DPDFB)[NEGP][MAXMAT], (*PACB)[NEGP][MAXMAT], * ZBR2;
} CEBR;

typedef struct {
	double WB[NBW], PBCUT[NEGP][MAXMAT], WBCUT[NEGP][MAXMAT], PDFB[NBW][NEGP][MAXMAT], DPDFB[NBW][NEGP][MAXMAT], PACB[NBW][NEGP][MAXMAT], ZBR2[MAXMAT];
} hd_CEBR;


typedef struct {
	double* EBT, (*XS)[NBE], * TXS, * X, * Y;

} CEBR01;


typedef struct {
	double(*P0)[NEGP][MAXMAT];

}CEBR02;

typedef struct {
	double* BET, * BK, (*BP1)[21][6][MAXMAT], (*BP2)[21][6][MAXMAT], * ZBEQ;
}CBRANG;

typedef struct {
	double BET[6], BK[21], BP1[4][21][6][MAXMAT], BP2[4][21][6][MAXMAT], ZBEQ[MAXMAT];
} hd_CBRANG;

typedef struct {

	double* EI, * EE, * CPS, * AMOL;
	int* MOM;

} CEIN01;

typedef struct {
	int* IERGA, * NCALL;

} CSUMGA;

typedef struct {

	double* EI, * CPS, * BHA1, * BHA2, * BHA3, * BHA4;
	int* MOM;

} CPIN01;


//Elastic scattering simulation tables.
typedef struct {

	double* ETS, * ETL, * TH, * THR, * XMU, * XMUL, * ECS, * ETCS1, * ETCS2, (*EDCS)[NE],
		* PCS, * PTCS1, * PTCS2, (*PDCS)[NE], * DCSI, * DCSIL, * CSI, * TCS1I, * TCS2I;

} CDCSEP;

typedef struct {
	double(*XSE)[NEGP][NP], (*PSE)[NEGP][NP], (*ASE)[NEGP][NP], (*BSE)[NEGP][NP];
	int(*ITLE)[NEGP][NP], (*ITUE)[NEGP][NP];
}CEELDB;

typedef struct{
    double XSE[MAXMAT][NEGP][NP], PSE[MAXMAT][NEGP][NP], ASE[MAXMAT][NEGP][NP], BSE[MAXMAT][NEGP][NP];
	int ITLE[MAXMAT][NEGP][NP], ITUE[MAXMAT][NEGP][NP];
}hd_CEELDB;


typedef struct {
	double(*XSP)[NEGP][NP], (*PSP)[NEGP][NP], (*ASP)[NEGP][NP], (*BSP)[NEGP][NP];
	int(*ITLP)[NEGP][NP], (*ITUP)[NEGP][NP];
}CPELDB;

typedef struct {
	double XSP[MAXMAT][NEGP][NP], PSP[MAXMAT][NEGP][NP], ASP[MAXMAT][NEGP][NP], BSP[MAXMAT][NEGP][NP];
	int ITLP[MAXMAT][NEGP][NP], ITUP[MAXMAT][NEGP][NP];
} hd_CPELDB;

typedef struct {
	double* EELMAX, * PELMAX, (*RNDCED)[MAXMAT], (*RNDCPD)[MAXMAT];
}CELSEP;

typedef struct {
	double EELMAX[MAXMAT], PELMAX[MAXMAT], RNDCED[NEGP][MAXMAT], RNDCPD[NEGP][MAXMAT];
} hd_CELSEP;

typedef struct {
	double* FACTE, * Q2MAX;
	int* MM, * MOM;
} CGRA00;

typedef struct {
	double(*FF)[MAXMAT], * ERA, (*XSRA)[MAXMAT];
	int* IED, * IEU, * NE;
}CGRA01;

typedef struct {
	double FF[NQ][MAXMAT], ERA[NEX], XSRA[NEX][MAXMAT];
	int IED[NEGP], IEU[NEGP], NE;
} hd_CGRA01;


typedef struct {
	double* QQ, (*AR)[MAXMAT], (*BR)[MAXMAT], (*CR)[MAXMAT], (*DR)[MAXMAT], * FF0, * QQM;
}CGRA02;

typedef struct {
	double(*QRA)[NP2], (*PRA)[NP2], (*DPRA)[NP2], (*ARA)[NP2],
		(*BRA)[NP2], (*PMAX)[NEGP];
	int(*ITLRA)[NP2], (*ITURA)[NP2];
}CGRA03;


typedef struct {
	double QRA[MAXMAT][NP2], PRA[MAXMAT][NP2], DPRA[MAXMAT][NP2], ARA[MAXMAT][NP2],
		BRA[MAXMAT][NP2], PMAX[MAXMAT][NEGP];
	int ITLRA[MAXMAT][NP2], ITURA[MAXMAT][NP2];
} hd_CGRA03;

typedef struct {
	double* ZEQPP, (*F0)[MAXMAT], * BCB;
}CGPP00;

typedef struct {
	double ZEQPP[MAXMAT], F0[2][MAXMAT], BCB[MAXMAT];
}hd_CGPP00;


//novos commons

typedef struct {
	int* ISEED1, * ISEED2;
}RSEED;

typedef struct {
	int ISEED1[tamISEEDs], ISEED2[tamISEEDs];
} hd_RSEED;




typedef struct {
	char* TITLE, * TITLE2;

} CTITLE;

typedef struct {
	char* DATE23;

} CDATE;

typedef struct {
	double* DSMAX, (*EABSB)[3];
}CSPGEO;

typedef struct {
	double DSMAX[NB], EABSB[NB][3];
} hd_CSPGEO;

//Forçando interação, janelas de peso.
typedef struct {
	double(*WLOW)[NBV], (*WHIG)[NBV];
	bool(*LFORCE)[NBV];

}CFORCI;

//Divisão de raios-X
typedef struct {
	int* IXRSPL, * ILBA;
	bool* LXRSPL;
}CXRSPL;

typedef struct {
	int IXRSPL[NBV], ILBA[pilhaPart][5];
	bool LXRSPL[NBV];
} hd_CXRSPL;

//Definição de origem.
// Particulas primarias

typedef struct {

	double* CTHL, * DCTH, * PHIL, * DPHI;
	int* KPARP, * JOBEND;
	bool* LSCONE, * LGPOL, * LPSF;
} CSOUR0;

typedef struct {
	double* E0, * EPMAX, * SP10, * SP20, * SP30;

} CSOUR1;

//Espectro de Energia
typedef struct {
	double* ESRC, * PSRC, * FSRC;
	int* IASRC;
	bool* LSPEC;

}CSOUR2;

// Fonte estendida

typedef struct {
	double* SX0, * SY0, * SZ0, * SSX, * SSY, * SSZ;
	int* IXSBOD;
	bool* LEXSRC, * LEXBD;

} CSOUR3;

//Arquivo de espaço de fase de entrada
typedef struct {
	double* WGMIN, * RWGMIN, * WGMAX, * RLREAD;
	int* IPSFI, * NPSF, * NPSN, * NSPLIT, * KODEPS;
} CSOUR4;

typedef struct {
	char(*PSFI)[20];

}CSOUR5;

// contadores discretos

typedef struct {
	double* PRIM, * PRIM2, * DPRIM; //Numero de particulas IEXIT;
	double(*SEC)[3], (*SEC2)[3], (*DSEC)[3]; // Geradores de particulas secundarias.
	double* AVW, * AVW2, * DAVW; //Cosseno final do diretor polar.
	double* AVA, * AVA2, * DAVA; // Angulo final polar
	double* AVE, * AVE2, * DAVE; // Energia final
}CNT0;

typedef struct {
	double PRIM[3], PRIM2[3], DPRIM[3][pilhaPart]; //Numero de particulas IEXIT;
	double SEC[3][3], SEC2[3][3], DSEC[3][3][pilhaPart]; // Geradores de particulas secundarias.
	double AVW[2], AVW2[2], DAVW[2][pilhaPart]; //Cosseno final do diretor polar.
	double AVA[2], AVA2[2], DAVA[2][pilhaPart]; // Angulo final polar
	double AVE[2], AVE2[2], DAVE[2][pilhaPart]; // Energia final
} hd_CNT0;


// Energias depositadas em vários corpos.

typedef struct {
	double* TDEBO, * TDEBO2, *DEBO;

}CNT1;

typedef struct {
	double TDEBO[NB], TDEBO2[NB], DEBO[NB][pilhaPart];
} hd_CNT1;

//Distribuições contínuas.
//Espectro de energia da fonte.

typedef struct {
	double* SHIST; //definicao
	int* NSEB;
}CNT2;

typedef struct {
	double(*SEDS)[3], (*SEDS2)[3], * DSDE, * RDSDE;
	int* NSDE;
}CNT3;

typedef struct {
	double SEDS[NSEM][3][pilhaPart],  SEDS2[NSEM][3][pilhaPart], DSDE, RDSDE;
	int NSDE;
} hd_CNT3;


//Detectores (até diferentes detectores NIDM).
typedef struct {
	double* RLAST, * RWRITE;
	int* IDCUT, (*KKDI)[NIDM], * IPSF, * NID, * NPSFO, * IPSFO;
}CNT4;

typedef struct {
	double RLAST, RWRITE;
	int IDCUT[NIDM], KKDI[3][NIDM], IPSF[NIDM], NID, NPSFO, IPSFO;
} hd_CNT4;

typedef struct {
	double* DEDE;
	int* KBDE, * NED;
}CNT5;

typedef struct {
	double DEDE[NIDM][pilhaPart];
	int KBDE[NB], NED;
} hd_CNT5;

typedef struct {
	bool* LDOSEM;
}CNT6;

typedef struct {
	bool LDOSEM;
} hd_CNT6;

// Detalhes do trabalho

//Arquivo dump
typedef struct {
	bool* LDUMP;
	char* PFILED;
}CDUMP;

//controlador de tempo e contador de simulacoes.

typedef struct {
	double* TSIM, * TSEC, * TSECA, * TSECAD, * CPUT0, * DUMPP, * DSHN, * SHN;
	int* N;
}CNTRL;

typedef struct {
	double TSIM, TSEC, TSECA, TSECAD, CPUT0, DUMPP, DSHN, SHN;
	int N[pilhaPart];
} hd_CNTRL;

typedef struct {
	double* CPCT, * CPST, * SPCT, * SPST, * SPHI, * CPHI, * STHE, * CTHE, * CAPER;
}CGCONE;

typedef struct {
	double CPCT, CPST, SPCT, SPST, SPHI, CPHI, STHE, CTHE, CAPER;
}hd_CGCONE;

typedef struct {
	double* EL, * EU, * THL, * THU, * BSE, * RBSE, * BSTH, * RBSTH, * BSPH, * RBSPH;
	double(*PDE)[2][3], (*PDE2)[2][3], (*PDEP)[2][3];
	double(*PDA)[NBTHM][3], (*PDA2)[NBTHM][3], (*PDAP)[NBTHM][3];
	int(*LPDE)[2][3];
	int(*LPDA)[NBTHM][3];
	int* NE, * NTH, * NPH, *LLE;
	bool * LLTH;
}CENANG;


typedef struct {
	double EL, EU, THL, THU, BSE, RBSE, BSTH, RBSTH, BSPH, RBSPH;
	double PDE[NBEM][2][3], PDE2[NBEM][2][3], PDEP[NBEM][2][3];
	double PDA[NBPHM][NBTHM][3], PDA2[NBPHM][NBTHM][3], PDAP[NBPHM][NBTHM][3];
	int LPDE[NBEM][2][3];
	int LPDA[NBPHM][NBTHM][3];
	int NE, NTH, NPH, LLE;
	bool LLTH;
} hd_CENANG;


typedef struct {
	double* EL, * EU, * BSE, * RBSE;
	double(*ET)[NIDM], * EDEP, * EDEP2, * EDEPP;
	double(*DIT)[NIDM], (*DIT2)[NIDM], (*DITP)[NIDM];
	double(*DIP)[NBEM2][NIDM], (*DIP2)[NBEM2][NIDM], (*DIPP)[NBEM2][NIDM];
	double(*FLT)[NIDM], (*FLT2)[NIDM], (*FLTP)[NIDM];
	double(*FLP)[NBEM2][NIDM], (*FLP2)[NBEM2][NIDM], (*FLPP)[NBEM2][NIDM];
	double* AGEL, * AGEU, * BAGE, * RBAGE, (*AGE)[NIDM];
	double(*AGE2)[NIDM], (*AGEP)[NIDM];
	int* LEDEP;
	int(*LDIT)[NIDM];
	int(*LDIP)[NBEM2][NIDM];
	int(*LFLT)[NIDM];
	int(*LFLP)[NBEM2][NIDM];
	int(*LAGEA)[NIDM];
	int* IDCUT, * NE, *LLE;
	bool * LLAGE;
	int* NAGE, * NID;
}CIMDET;

typedef struct {
	double EL[NIDM], EU[NIDM], BSE[NIDM], RBSE[NIDM];
	double ET[NBEM2+1][NIDM], EDEP[NIDM], EDEP2[NIDM], EDEPP[NIDM];
	double DIT[NBEM2][NIDM], DIT2[NBEM2][NIDM], DITP[NBEM2][NIDM];
	double DIP[3][NBEM2][NIDM], DIP2[3][NBEM2][NIDM], DIPP[3][NBEM2][NIDM];
	double FLT[NBEM2][NIDM], FLT2[NBEM2][NIDM], FLTP[NBEM2][NIDM];
	double FLP[3][NBEM2][NIDM], FLP2[3][NBEM2][NIDM], FLPP[3][NBEM2][NIDM];
	double AGEL[NIDM], AGEU[NIDM], BAGE[NIDM], RBAGE[NIDM], AGE[NBEM2][NIDM];
	double AGE2[NBEM2][NIDM], AGEP[NBEM2][NIDM];
	int LEDEP[NIDM];
	int LDIT[NBEM2][NIDM];
	int LDIP[3][NBEM2][NIDM];
	int LFLT[NBEM2][NIDM];
	int LFLP[3][NBEM2][NIDM];
	int LAGEA[NBEM2][NIDM];
	int IDCUT[NIDM], NE[NIDM], LLE[NIDM];
	bool LLAGE[NIDM];
	int NAGE[NIDM], NID;
} hd_CIMDET;



typedef struct {
	double* EL, * EU, * BSE, * RBSE, * EDEP, * EDEP2, (*DET)[NIDM];
	int* NE, * NID;
	int* LLE;
}CENDET;


typedef struct {
	double  EL[NIDM], EU[NIDM], BSE[NIDM], RBSE[NIDM], EDEP[NIDM], EDEP2[NIDM], DET[NBEM2][NIDM];
	int NE[NIDM], NID[NIDM];
	int LLE[NIDM];
} hd_CENDET;

typedef struct {
	double(*DOSE)[NDYM][NDXM], (*DOSE2)[NDYM][NDXM], (*DOSEP)[NDYM][NDXM];
	int(*LDOSE)[NDYM][NDXM], * KDOSE;
} CDOSE1;

typedef struct {
	double DOSE[NDZM][NDYM][NDXM], DOSE2[NDZM][NDYM][NDXM], DOSEP[NDZM][NDYM][NDXM];
	int LDOSE[NDZM][NDYM][NDXM], KDOSE;
} hd_CDOSE1;

typedef struct {
	double* DDOSE, * DDOSE2, * DDOSEP;
	int* LDDOSE;
}CDOSE2;

typedef struct {
	double DDOSE[NDZM], DDOSE2[NDZM], DDOSEP[NDZM];
	int LDDOSE[NDZM];
} hd_CDOSE2;

typedef struct {
	double* DXL, * DXU, * BDOSE, * RBDOSE;
	int* NDB;
}CDOSE3;

typedef struct {
	double DXL[3], DXU[3], BDOSE[3], RBDOSE[3];
	int NDB[3];
} hd_CDOSE3;

typedef struct {
	double(*VMASS)[NDYM][NDXM];
}CDOSE4;

typedef struct {
	double VMASS[NDZM][NDYM][NDXM];
} hd_CDOSE4;


typedef struct {
	double* ELAST1, * ELAST2;
	int* MHINGE, * KSOFTE, * KSOFTI, * KDELTA;
}CJUMP1;

typedef struct {
	double ELAST1[pilhaPart], ELAST2[pilhaPart];
	int MHINGE[pilhaPart], KSOFTE[pilhaPart], KSOFTI[pilhaPart], KDELTA[pilhaPart];
}hd_CJUMP1;

typedef struct {
	double* ES, * XS, * YS, * ZS, * US, * VS, * WS, * WGHTS, * SP1S, * SP2S, * SP3S, * PAGES;
	int* KS, * IBODYS, * MS, (*ILBS)[5], * IPOLS, * NSEC;
}SECST;

typedef struct {
	double ES[NMS],  XS[NMS],  YS[NMS],  ZS[NMS],  US[NMS],  VS[NMS],  WS[NMS],  WGHTS[NMS],  SP1S[NMS],  SP2S[NMS],  SP3S[NMS],  PAGES[NMS];
	int KS[NMS],  IBODYS[NMS],  MS[NMS], ILBS[NMS][5],  IPOLS[NMS],  NSEC;
}hd_SECST;

typedef struct {
	double* P, * ST, * DST, * DSR, * W1, * W2, * T1, * T2;
}CJUMP0;

typedef struct {
	double P[8][pilhaPart], ST[pilhaPart], DST[pilhaPart], DSR[pilhaPart], W1[pilhaPart], W2[pilhaPart], T1[pilhaPart], T2[pilhaPart];
}hd_CJUMP0;

typedef struct {
	int* ILBA;
}CHIST;

typedef struct {
	int ILBA[pilhaPart][5];
} hd_CHIST;

typedef struct{
	int nPRITRACK, nSECTRACK_E, nSECTRACK_G, nSECTRACK_P, nFINISH, TIPO;
} hd_nTRACKS;


typedef struct{
	int step1, step2, step3, step4, step5, step6, step7_0, step7_1, step7, step8, step9, step10, 
		step11, step12, step13, step14, step15, step16, step17, step18, step19, step20;

} hd_steps;

PENELOPE_MOD PENELOPE_mod_;
PENGEOM_MOD PENGEOM_mod_;
TRACK_MOD TRACK_mod_;
QSURF QSURF_;
QTREE QTREE_;
QBODY QBODY_;
CECUTR CECUTR_;
CSGAWR CSGAWR_;
CERSEC CERSEC_;
CEGRID CEGRID_;
CESI0 CESI0_;
CPSI0 CPSI0_;
CADATA CADATA_;
CGPH00 CGPH00_;
CRELAX CRELAX_;
CRITAA CRITAA_;
CRNDG3 CRNDG3_;
CRITA CRITA_;
CRITAN CRITAN_;
COMPOS COMPOS_;
CRANGE CRANGE_;
CEIN CEIN_;
CEINTF CEINTF_;
CEIN00 CEIN00_;
CPIN00 CPIN00_;
CEINAC CEINAC_;
CESIAC CESIAC_;
CESIN CESIN_;
CPINAC CPINAC_;
CPSIAC CPSIAC_;
CPSIN CPSIN_;
CGCO CGCO_;
CEIMFP CEIMFP_;
CLAS1E CLAS1E_;
CPIMFP CPIMFP_;
CLAS1P CLAS1P_;
CEEL00 CEEL00_;
CBRYLD CBRYLD_;
CGIMFP CGIMFP_;
CGPH01 CGPH01_;
CGPP01 CGPP01_;
CEBR CEBR_;
CEBR01 CEBR01_;
CEBR02 CEBR02_;
CBRANG CBRANG_;
CEIN01 CEIN01_;
CSUMGA CSUMGA_;
CPIN01 CPIN01_;
CDCSEP CDCSEP_;
CEELDB CEELDB_;
CPELDB CPELDB_;
CELSEP CELSEP_;
CGRA00 CGRA00_;
CGRA01 CGRA01_;
CGRA02 CGRA02_;
CGRA03 CGRA03_;
CGPP00 CGPP00_;
RSEED RSEED_;
CTITLE CTITLE_;
CDATE CDATE_;
CSPGEO CSPGEO_;
CFORCI CFORCI_;
CXRSPL CXRSPL_;
CSOUR0 CSOUR0_;
CSOUR1 CSOUR1_;
CSOUR2 CSOUR2_;
CSOUR3 CSOUR3_;
CSOUR4 CSOUR4_;
CSOUR5 CSOUR5_;
CNT0 CNT0_;
CNT1 CNT1_;
CNT2 CNT2_;
CNT3 CNT3_;
CNT4 CNT4_;
CNT5 CNT5_;
CNT6 CNT6_;
CDUMP CDUMP_;
CNTRL CNTRL_;
CGCONE CGCONE_;
CENANG CENANG_;
CIMDET CIMDET_;
CENDET CENDET_;
CDOSE1 CDOSE1_;
CDOSE2 CDOSE2_;
CDOSE3 CDOSE3_;
CDOSE4 CDOSE4_;
CJUMP1 CJUMP1_;
SECST SECST_;
CJUMP0 CJUMP0_;
CHIST CHIST_;
hd_nTRACKS nTRACKS_;

/*hd_TRACK_MOD *PRITRACK; 
hd_TRACK_MOD_SEC *SECTRACK_G; 
hd_TRACK_MOD_SEC *SECTRACK_E;
hd_TRACK_MOD_SEC *SECTRACK_P; 

hd_TRACK_MOD_SEC *vTrack_Simular;*/

hd_TRACK_MOD PRITRACK; 
hd_TRACK_MOD_SEC SECTRACK_G; 
hd_TRACK_MOD_SEC SECTRACK_E;
hd_TRACK_MOD_SEC SECTRACK_P; 

hd_TRACK_MOD vTrack_Simular;

__device__ hd_wSHOWERS dg_wSHOWERS_;

extern __shared__ int shared[];

int tam_sISEEDS = 2*blockSize*sizeof(int);

//__shared__ hd_TRACK_MOD_SHARED TRACK_MOD_SHARED;

/*__shared__ int ds_ISEED1[blockSize];
__shared__ int ds_ISEED2[blockSize];*/

/*__shared__ int knock_e_eela;
__shared__ int knock_e_eina;
__shared__ int knock_e_ebra;
__shared__ int knock_e_esia;
__shared__ int knock_e_eaux;

__device__ int dg_knock_e_eela[pilhaPart];
__device__ int dg_knock_e_eina[pilhaPart];
__device__ int dg_knock_e_ebra[pilhaPart];
__device__ int dg_knock_e_esia[pilhaPart];
__device__ int dg_knock_e_eaux[pilhaPart];

__device__ int dg_knock_g_graa[pilhaPart];
__device__ int dg_knock_g_gcoa[pilhaPart];
__device__ int dg_knock_g_gpha[pilhaPart];
__device__ int dg_knock_g_gppa[pilhaPart];
__device__ int dg_knock_g_gaux[pilhaPart];

__device__ int dg_knock_p_eela[pilhaPart];
__device__ int dg_knock_p_pina[pilhaPart];
__device__ int dg_knock_p_ebra[pilhaPart];
__device__ int dg_knock_p_psia[pilhaPart];
__device__ int dg_knock_p_pana[pilhaPart];
__device__ int dg_knock_p_paux[pilhaPart];*/


__device__ curandState_t dg_rand[pilhaPart];
//__device__ curandStatePhilox4_32_10_t dg_rand[pilhaPart];
//__device__ curandStateMRG32k3a_t dg_rand[pilhaPart];



//DECLARACOES PARA COPIA NA GPU
__device__ hd_CEELDB dg_CEELDB_;
hd_CEELDB *d_CEELDB;

__device__ hd_SECST dg_SECST_;
hd_SECST *d_SECST;

__device__ hd_PENGEOM_MOD dg_PENGEOM_mod_;
hd_PENGEOM_MOD * d_PENGEOM_mod;

__device__ hd_CGCONE dg_CGCONE_;
hd_CGCONE * d_CGCONE;

__device__ hd_QSURF dg_QSURF_;
hd_QSURF* d_QSURF;

__device__ hd_QTREE dg_QTREE_;
hd_QTREE* d_QTREE;

//__device__ hd_RSEED dg_RSEED_;
//hd_RSEED* d_RSEED;

__device__ hd_RSEED dg_RSEED2_;
hd_RSEED* d_RSEED2_;
hd_RSEED h_RSEED2_;

__device__ hd_steps dg_steps;
hd_steps h_steps;

__device__ hd_CIMDET dg_CIMDET_;
hd_CIMDET* d_CIMDET;

__device__ hd_CERSEC dg_CERSEC_;
hd_CERSEC* d_CERSEC;

__device__ hd_CEGRID dg_CEGRID_;
hd_CEGRID* d_CEGRID;

//__device__ hd_CJUMP1 dg_CJUMP1_;
__device__ hd_CJUMP1 dg_CJUMP1_;
hd_CJUMP1* d_CJUMP1;

__device__ hd_CEIMFP dg_CEIMFP_;
hd_CEIMFP* d_CEIMFP;

__device__ hd_CPIMFP dg_CPIMFP_;
hd_CPIMFP* d_CPIMFP;

//__device__ hd_CJUMP0 dg_CJUMP0_;
__device__ hd_CJUMP0 dg_CJUMP0_;
hd_CJUMP0* d_CJUMP0;

__device__ hd_CGIMFP dg_CGIMFP_;
hd_CGIMFP* d_CGIMFP;

__device__ hd_CHIST dg_CHIST_;
hd_CHIST* d_CHIST;

__device__ hd_COMPOS dg_COMPOS_;
hd_COMPOS* d_COMPOS;

__device__ hd_CADATA dg_CADATA_;
hd_CADATA* d_CADATA;

__device__ hd_CEIN dg_CEIN_;
hd_CEIN* d_CEIN;

__device__ hd_CGCO dg_CGCO_;
hd_CGCO* d_CGCO;

__device__ hd_CEBR dg_CEBR_;
hd_CEBR* d_CEBR;

__device__ hd_CELSEP dg_CELSEP_;
hd_CELSEP* d_CELSEP;

__device__ hd_CECUTR dg_CECUTR_;
hd_CECUTR* d_CECUTR;

__device__ hd_CESIAC dg_CESIAC_;
hd_CESIAC* d_CESIAC;

__device__ hd_CEINAC dg_CEINAC_;
hd_CEINAC* d_CEINAC;

__device__ hd_CBRANG dg_CBRANG_;
hd_CBRANG* d_CBRANG;

__device__ hd_CRELAX dg_CRELAX_;
hd_CRELAX* d_CRELAX;

__device__ hd_CGRA01 dg_CGRA01_;
hd_CGRA01* d_CGRA01;

__device__ hd_CGRA03 dg_CGRA03_;
hd_CGRA03* d_CGRA03;

__device__ hd_CGPH00 dg_CGPH00_;
hd_CGPH00* d_CGPH00;

__device__ hd_CGPP00 dg_CGPP00_;
hd_CGPP00* d_CGPP00;

__device__ hd_CGPP01 dg_CGPP01_;
hd_CGPP01* d_CGPP01;

__device__ hd_CPELDB dg_CPELDB_;
hd_CPELDB* d_CPELDB;

__device__ hd_CPSIAC dg_CPSIAC_;
hd_CPSIAC* d_CPSIAC;

__device__ hd_CPINAC dg_CPINAC_;
hd_CPINAC* d_CPINAC;

__device__ hd_CENANG dg_CENANG_;
hd_CENANG* d_CENANG;

__device__ hd_CENDET dg_CENDET_;
hd_CENDET* d_CENDET;

__device__ hd_PENELOPE_MOD dg_PENELOPE_mod_;
hd_PENELOPE_MOD* d_PENELOPE_mod;

/*__device__ hd_TRACK_MOD dg_TRACK_mod_;
hd_TRACK_MOD* d_TRACK_mod;*/

__device__ hd_TRACK_MOD dg_TRACK_mod_;
//__device__ hd_TRACK_MOD *dg_TRACK_mod_;
hd_TRACK_MOD* d_TRACK_mod;

__device__ hd_QBODY dg_QBODY_;
hd_QBODY* d_QBODY;

__device__ hd_CDOSE1 dg_CDOSE1_;
hd_CDOSE1* d_CDOSE1;

__device__ hd_CDOSE2 dg_CDOSE2_;
hd_CDOSE2* d_CDOSE2;

__device__ hd_CDOSE3 dg_CDOSE3_;
hd_CDOSE3* d_CDOSE3;

__device__ hd_CDOSE4 dg_CDOSE4_;
hd_CDOSE4* d_CDOSE4;

__device__ hd_CRNDG3 dg_CRNDG3_;
hd_CRNDG3* d_CRNDG3;

__device__ hd_CNT0 dg_CNT0_;
hd_CNT0* d_CNT0;
hd_CNT0 h_CNT0_;

__device__ hd_CNT1 dg_CNT1_;
hd_CNT1* d_CNT1;
hd_CNT1 h_CNT1_;

__device__ hd_CNT3 dg_CNT3_;
hd_CNT3* d_CNT3;
hd_CNT3 h_CNT3_;

__device__ hd_CNT4 dg_CNT4_;
hd_CNT4* d_CNT4;

__device__ hd_CNT5 dg_CNT5_;
hd_CNT5* d_CNT5;

__device__ hd_CNT6 dg_CNT6_;
hd_CNT6* d_CNT6;

__device__ hd_CNTRL dg_CNTRL_;
hd_CNTRL* d_CNTRL;

__device__ hd_CXRSPL dg_CXRSPL_;
hd_CXRSPL* d_CXRSPL;

__device__ hd_CSPGEO dg_CSPGEO_;
hd_CSPGEO* d_CSPGEO;

__device__ hd_TRACK_MOD dg_PRITRACK_;
hd_TRACK_MOD* d_PRITRACK;

__device__ hd_TRACK_MOD_SEC dg_SECTRACK_G_;
hd_TRACK_MOD_SEC *d_SECTRACK_G;

__device__ hd_TRACK_MOD_SEC dg_SECTRACK_E_;
hd_TRACK_MOD_SEC *d_SECTRACK_E;

__device__ hd_TRACK_MOD_SEC dg_SECTRACK_P_;
hd_TRACK_MOD_SEC *d_SECTRACK_P;

__device__ hd_nTRACKS dg_nTRACKS_;
hd_nTRACKS* d_nTRACKS;

__device__ int size;
int * d_size;

//__device__ double d_S[NS2M];
//__device__ int d_IS[NS2M];
__device__ int d_wIPOLI = 0;



//__device__ double d_S[pilhaPart*NS2M/100];
//__device__ int d_IS[pilhaPart*NS2M/100];

__device__ double d_S[(NS2M/100)+48][pilhaPart]; //matrix (2048,4096)
__device__ int d_IS[(NS2M/100)+48][pilhaPart];






