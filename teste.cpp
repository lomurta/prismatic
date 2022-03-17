
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <cmath>
#include <math.h>
#include <algorithm>
#include <stdint.h>
//#include <chrono>
#include <time.h>
#include <iostream>
#include <fstream>
#include <string>

using namespace std;


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
static const int NM=512;
static const int NR=128;
static const double EPS=1e-10;
static const double ZERO=1.0e-75;
static const double ZEROT=0.1*ZERO;
static const int NIP=51;
static const double A0B = 5.2917721092e-9; //Bohr radius (cm)
static const double HREV= 27.21138505e0; //Hartree energy (eV)
static const double AVOG= 6.02214129e23; //Avogadro's number
static const double SL=137.035999074e0; //Speed of light (1/alpha)
static const double PI=3.1415926535897932e0; 
static const double FOURPI = 4.0e0*PI; 
static const int NO = 512; 
static const int NOCO=512;
static const int NDIM=12000;
static const double HBAR = 6.58211928e-16; //Constante de Planck eV*s
static const int NTRAN = 2500;

static const int NBE = 57;
static const int NBW = 32;

static const int NE=96;
static const int NA=606;
static const int NP=128;
static const int NQ=250;
static const int NEX=1024;
static const int NP2=150;
static const int NSEM=1000;
static const int NPSFM=100;
static const int NIDM=25;
static const int NBV=5000;

static const int NBEM=1500;
static const int NBTHM=1800;
static const int NBPHM=180;
static const int NBEM2=1000;


static const int NDXM=201;
static const int NDYM=201;
static const int NDZM=201;

char LINHA[200];
char APOIO[200];

char SPCDIO[NIDM][20];
char SPCFLO[NIDM][20];
char SPCAGE[NIDM][20];
char SPCDEO[NIDM][20];

double S[NS2M];
int IS[NS2M];

int imprimiu=0;
int wIPOLI=0;

clock_t start, end;




 typedef	struct {
    double *AXX, *AXY, *AXZ, *AYY, *AYZ,
        *AZZ, *AX, *AY, *AZ, *A0;
    int *NSURF, *KPLANE;
} QSURF;

typedef struct{
    int *NBODYS, *KMOTH, (*KDGHT)[NB], (*KSURF)[NB], (*KFLAG)[NB], *KSP, *NWARN;
}  QTREE;
		
typedef struct {
    double *E, *X, *Y, *Z, *U, *V, *W, *WGHT, *SP1, *SP2, *SP3, *PAGE;
    int *KPAR, *IBODY, *MAT, *ILB, *IPOL;
    bool *LAGE;
}TRACK_MOD;	
	
typedef struct {
   double (*EABS)[3], *C1, *C2, *WCC, *WCR, *DEN, *RDEN, *E0STEP, *DESOFT, *SSOFT;
   int *NMS, *NEGP, *NMAT;
}PENELOPE_MOD;

 typedef struct {
    char (*BALIAS)[5];
    double *DSTOT;
    int *MATER, *KDET, *KSLAST, *NBODY;
    bool *LVERB;
    
} PENGEOM_MOD;


 typedef struct {
 	int (*KBODY)[NB], *KBOMO;    
} QBODY;

 typedef struct {
 	double *ECUTR;    
} CECUTR; //parametros de simula��o

 typedef struct {
 	int *ISGAW;    
} CSGAWR;

 typedef struct {
 	int *IERSEC;   
} CERSEC;

typedef struct {
	double *EMIN, *EL, *EU, *ET, *DLEMP, *DLEMP1, *DLFC, *XEL, *XE, *XEK;
	int *KE;	
} CEGRID; //Rede de energia e constantes de interpolacao para a energia atual.


typedef struct{
	double (*XESI)[NRP];
	int *IESIF, *IESIL, *NSESI, *NCURE;
} CESI0; //Ioniza��o da camada interna por impacto de el�trons e p�sitrons.


typedef struct{
	double (*XPSI)[NRP];
	int *IPSIF, *IPSIL, *NSPSI, *NCURP;
} CPSI0; //Ioniza��o da camada interna por impacto de el�trons e p�sitrons.


typedef struct {
	double *ATW, *EPX, *RSCR, *ETA, (*EB)[99], (*ALW)[99], (*CP0)[99];
	int (*IFI)[99], (*IKS)[99], *NSHT;
	char (*LASYMB)[2];  
	
} CADATA; //Dados elemento

typedef struct {
	double *EPH, (*XPH)[NTP];
	int *IPHF, *IPHL, *NPHS, *NCUR;
	
} CGPH00;


typedef struct {
	double *P, *ET, *F;
	int *IS0, *IS1, *IS2, (*IFIRST)[99], (*ILAST)[99], *NCUR, *KS, *MODER;
	
} CRELAX;


typedef struct {
	
	double *XA, *AA, *BA;
	double *FA;
	int *IA, *NPM1A;
	
} CRITAA;

typedef struct {
	
	double *X, *A, *B;
	double *F;
	int *KA, *NPM1;
	
} CRNDG3;


typedef struct {
	double *XT, *PAC, *DPAC, *A, *B;
	int *IL, *IU, *NPM1;
	
	double *QTI, *PACI, *DPACI, *AI, *BI;
	int *ITLI, *ITUI, *NPM1I;
	
	double *XTI;
	
} CRITA;

typedef struct{
	
	double *CNORM;
	
} CRITAN;


typedef struct{
	double (*STF)[MAXMAT], *ZT, *AT, *RHO, *VMOL;
	int (*IZ)[MAXMAT], *NELEM;
	
} COMPOS; //Dados Composi��o


typedef struct{
	
	double (*RANGE)[MAXMAT][3], (*RANGEL)[MAXMAT][3];

} CRANGE;

typedef struct {
	double *EXPOT, *OP2, (*F)[MAXMAT], (*UI)[MAXMAT], (*WRI)[MAXMAT];
	int (*KZ)[MAXMAT], (*KS)[MAXMAT], *NOSC;
	
}CEIN; //Colisoes Inelasticas

typedef struct{
	double *T1EI, *T2EI , *T1PI, *T2PI;
	
} CEINTF;

typedef struct{
	double *SEH0, *SEH1, *SEH2 , *SES0, *SES1, *SES2, *SET0, *SET1, *SET2;
} CEIN00; //Se��es transversais parciais de conchas / osciladores individuais.

typedef struct{
	
	double *SPH0, *SPH1, *SPH2, *SPS0, *SPS1, *SPS2, *SPT0, *SPT1, *SPT2;
	
} CPIN00;

typedef struct{
	
	double (*EINAC)[NEGP][MAXMAT];
	int (*IEIN)[MAXMAT], *NEIN;
	
	
} CEINAC; //Inel�stica de el�trons. e tabelas de ioniza��o de camada interna.


typedef struct{
	
	double (*ESIAC)[NEGP][MAXMAT];
	int (*IESI)[MAXMAT], *NESI;
	
	
} CESIAC; //Inel�stica de el�trons. e tabelas de ioniza��o de camada interna.


typedef struct{
	
	double (*XSEIN)[NEGP], (*XSESI)[NEGP];
	int *ISIE;
	
	
} CESIN; //Inel�stica de el�trons. e tabelas de ioniza��o de camada interna.


typedef struct {
	
	double (*PINAC)[NEGP][MAXMAT];
	int (*IPIN)[MAXMAT], *NPIN;
	
} CPINAC; //Positron inel�stica coll. e tabelas de ioniza��o de camada interna.


typedef struct{
	
	double (*PSIAC)[NEGP][MAXMAT];
	int (*IPSI)[MAXMAT], *NPSI;
	
	
} CPSIAC; //Positron inel�stica coll. e tabelas de ioniza��o de camada interna.


typedef struct{
	
	double (*XSPIN)[NEGP], (*XSPSI)[NEGP];
	int *ISIP;
	
	
} CPSIN; //Positron inel�stica coll. e tabelas de ioniza��o de camada interna.

//novos
typedef struct{
	
	double (*FCO)[MAXMAT], (*UICO)[MAXMAT], (*FJ0)[MAXMAT], (*PTRSH)[MAXMAT];
	int (*KZCO)[MAXMAT], (*KSCO)[MAXMAT], *NOSCCO;
	
} CGCO; //Espalhamento Compton

typedef struct{
	double (*SEHEL)[MAXMAT], (*SEHIN)[MAXMAT], (*SEISI)[MAXMAT], (*SEHBR)[MAXMAT], (*SEAUX)[MAXMAT], 
           (*SETOT)[MAXMAT], (*CSTPE)[MAXMAT], (*RSTPE)[MAXMAT], (*DEL)[MAXMAT], (*W1E)[MAXMAT],
           (*W2E)[MAXMAT], (*DW1EL)[MAXMAT], (*DW2EL)[MAXMAT], (*RNDCE)[MAXMAT], (*AE)[MAXMAT], (*BE)[MAXMAT],
    	    (*T1E)[MAXMAT], (*T2E)[MAXMAT];
    	    
} CEIMFP; //Tabela e simula��o do Eletron


typedef struct {
	double (*TSTPE)[MAXMAT], (*TSTRE)[MAXMAT], (*TRL1E)[MAXMAT], (*TRL2E)[MAXMAT];
	
} CLAS1E;

typedef struct{
	double (*SPHEL)[MAXMAT], (*SPHIN)[MAXMAT], (*SPISI)[MAXMAT], (*SPHBR)[MAXMAT], (*SPAN)[MAXMAT], 
           (*SPAUX)[MAXMAT], (*SPTOT)[MAXMAT], (*CSTPP)[MAXMAT], (*RSTPP)[MAXMAT], (*W1P)[MAXMAT],
           (*W2P)[MAXMAT], (*DW1PL)[MAXMAT], (*DW2PL)[MAXMAT], (*RNDCP)[MAXMAT], (*AP)[MAXMAT], (*BP)[MAXMAT],
           (*T1P)[MAXMAT], (*T2P)[MAXMAT];		
} CPIMFP; // Tabelas de simula��o dos positrons


typedef struct{
	double (*TSTPP)[MAXMAT], (*TSTRP)[MAXMAT], (*TRL1P)[MAXMAT], (*TRL2P)[MAXMAT];

} CLAS1P;

//Espalhamento elastico de eletrons e positrons

typedef struct{
	
	double *EJT, *XE0, *XE1, *XE2, *XP0, *XP1, *XP2, *T1E0, *T2E0, *T1P0, 
           *T2P0, *EJTL, *FJL, *A, *B, *C, *D; 
	
} CEEL00;

//Rendimentos radiativos de el�trons e p�sitrons.

typedef struct{ 
	double (*EBRY)[MAXMAT], (*PBRY)[MAXMAT];
} CBRYLD;

//tabelas de simula��o dos fotons

typedef struct{
	double (*SGRA)[MAXMAT], (*SGCO)[MAXMAT], (*SGPH)[MAXMAT], (*SGPP)[MAXMAT], (*SGAUX)[MAXMAT];
	
} CGIMFP;

typedef struct{
	   	double *ER, *XSR;
		int *NPHD;
	
} CGPH01;

typedef struct{
	double (*TRIP)[MAXMAT];
	
} CGPP01;

typedef struct {
	double *WB, (*PBCUT)[MAXMAT], (*WBCUT)[MAXMAT], (*PDFB)[NEGP][MAXMAT], (*DPDFB)[NEGP][MAXMAT], (*PACB)[NEGP][MAXMAT], *ZBR2;
	
} CEBR;


typedef struct {
	double *EBT, (*XS)[NBE], *TXS, *X, *Y;
	
} CEBR01;


typedef struct{
	double (*P0)[NEGP][MAXMAT];
	
}CEBR02;

typedef struct{
	
//	double *BET, *BK, (*BP1)[4][21][6], (*BP2)[4][21][6], *ZBEQ;
    double *BET, *BK, (*BP1)[21][6][MAXMAT], (*BP2)[21][6][MAXMAT], *ZBEQ;
	
}CBRANG;

typedef struct{
	
	double *EI, *EE, *CPS, *AMOL;
	int *MOM;
	
} CEIN01;

typedef struct{
	int *IERGA, *NCALL;
	
} CSUMGA;

typedef struct{
	
	double *EI, *CPS, *BHA1, *BHA2, *BHA3, *BHA4;
	int *MOM;
	
} CPIN01;


//Elastic scattering simulation tables.
typedef struct {
	
	double *ETS, *ETL, *TH, *THR, *XMU , *XMUL, *ECS, *ETCS1 , *ETCS2, (*EDCS)[NE],
     	    *PCS, *PTCS1, *PTCS2, (*PDCS)[NE], *DCSI, *DCSIL, *CSI, *TCS1I, *TCS2I;
	
} CDCSEP;

typedef struct {
	
	double (*XSE)[NEGP][NP], (*PSE)[NEGP][NP],(*ASE)[NEGP][NP], (*BSE)[NEGP][NP];
    int (*ITLE)[NEGP][NP], (*ITUE)[NEGP][NP];

}CEELDB;


typedef struct{
	
	double (*XSP)[NEGP][NP], (*PSP)[NEGP][NP],(*ASP)[NEGP][NP], (*BSP)[NEGP][NP];
    int(*ITLP)[NEGP][NP], (*ITUP)[NEGP][NP];
	
}CPELDB;

typedef struct {
	double *EELMAX, *PELMAX, (*RNDCED)[MAXMAT], (*RNDCPD)[MAXMAT];
	
}CELSEP;

typedef struct {
	double *FACTE, *Q2MAX;
	int *MM, *MOM;
} CGRA00;

typedef struct{
	double (*FF)[MAXMAT], *ERA, (*XSRA)[MAXMAT];
    int    *IED, *IEU, *NE;
}CGRA01;

typedef struct{
	double *QQ, (*AR)[MAXMAT], (*BR)[MAXMAT], (*CR)[MAXMAT], (*DR)[MAXMAT], *FF0, *QQM;
}CGRA02;

typedef struct{
	double  (*QRA)[NP2], (*PRA)[NP2], (*DPRA)[NP2], (*ARA)[NP2],
			(*BRA)[NP2], (*PMAX)[NEGP];
	int		(*ITLRA)[NP2], (*ITURA)[NP2];
}CGRA03;

typedef struct {
	double *ZEQPP, (*F0)[MAXMAT], *BCB;
	
}CGPP00;


//novos commons

typedef struct {
	int *ISEED1, *ISEED2;
}RSEED;

typedef struct{
	char *TITLE, *TITLE2;

} CTITLE;

typedef struct{
	char *DATE23;

} CDATE;

typedef struct{
	double *DSMAX, (*EABSB)[3];

}CSPGEO;

//Forçando interação, janelas de peso.
typedef struct{
	double (*WLOW)[NBV], (*WHIG)[NBV];
	bool (*LFORCE)[NBV];

}CFORCI;

//Divisão de raios-X
typedef struct {
	int *IXRSPL, *ILBA;
	bool *LXRSPL;
}CXRSPL;

//Definição de origem.
// Particulas primarias

typedef struct {

	double *CTHL, *DCTH, *PHIL, *DPHI;
	int *KPARP, *JOBEND;
	bool *LSCONE, *LGPOL, *LPSF;
} CSOUR0;

typedef struct{
	double *E0, *EPMAX, *SP10, *SP20, *SP30;

} CSOUR1;

//Espectro de Energia
typedef struct{
	double *ESRC, *PSRC, *FSRC;
	int *IASRC;
	bool *LSPEC;

}CSOUR2;

// Fonte estendida

typedef struct {
	double *SX0, *SY0, *SZ0, *SSX, *SSY, *SSZ;
	int *IXSBOD;
	bool *LEXSRC, *LEXBD;

} CSOUR3;

//Arquivo de espaço de fase de entrada
typedef struct{
	double *WGMIN, *RWGMIN, *WGMAX, *RLREAD;
	int *IPSFI, *NPSF, *NPSN, *NSPLIT, *KODEPS;
} CSOUR4;

typedef struct{
	char (*PSFI)[20];

}CSOUR5;

// contadores discretos

typedef struct{
    double *PRIM ,*PRIM2 ,*DPRIM; //Numero de particulas IEXIT;
	double (*SEC)[3], (*SEC2)[3], (*DSEC)[3]; // Geradores de particulas secundarias.
	double *AVW, *AVW2, *DAVW; //Cosseno final do diretor polar.
	double *AVA, *AVA2, *DAVA; // Angulo final polar
	double *AVE, *AVE2, *DAVE; // Energia final
}CNT0;


// Energias depositadas em vários corpos.

typedef struct{
	double *TDEBO, *TDEBO2, *DEBO;

}CNT1;

//Distribuições contínuas.
//Espectro de energia da fonte.

typedef struct{
	double *SHIST; //definicao
	int *NSEB;
}CNT2;

typedef struct{
	double (*SEDS)[3], (*SEDS2)[3], *DSDE, *RDSDE;
	int *NSDE; 

}CNT3;

//Detectores (até diferentes detectores NIDM).
typedef struct{
	double *RLAST, *RWRITE;
	int *IDCUT, (*KKDI)[NIDM], *IPSF, *NID, *NPSFO, *IPSFO;

}CNT4;

typedef struct{
	double *DEDE;
	int *KBDE, *NED;
}CNT5;

typedef struct{
	bool *LDOSEM;

}CNT6;

// Detalhes do trabalho

//Arquivo dump
typedef struct{
	bool *LDUMP;
	char *PFILED;
}CDUMP;

//controlador de tempo e contador de simulacoes.

typedef struct{
	double *TSIM, *TSEC, *TSECA, *TSECAD, *CPUT0, *DUMPP, *DSHN, *SHN;
	int *N;

}CNTRL;

typedef struct{
	double *CPCT, *CPST, *SPCT, *SPST, *SPHI, *CPHI, *STHE, *CTHE, *CAPER;
}CGCONE;

typedef struct{
	double *EL, *EU, *THL, *THU, *BSE, *RBSE, *BSTH, *RBSTH, *BSPH, *RBSPH;
	double (*PDE)[2][3], (*PDE2)[2][3], (*PDEP)[2][3];
	double (*PDA)[NBTHM][3], (*PDA2)[NBTHM][3], (*PDAP)[NBTHM][3];
	bool (*LPDE)[2][3], (*LPDA)[NBTHM][3];
	int *NE, *NTH, *NPH;
	bool *LLE, *LLTH;
}CENANG;


typedef struct{
	double *EL, *EU, *BSE, *RBSE;
	double (*ET)[NIDM], *EDEP, *EDEP2, *EDEPP;
    double  (*DIT)[NIDM], (*DIT2)[NIDM], (*DITP)[NIDM];
    double  (*DIP)[NBEM2][NIDM], (*DIP2)[NBEM2][NIDM], (*DIPP)[NBEM2][NIDM];
    double  (*FLT)[NIDM], (*FLT2)[NIDM], (*FLTP)[NIDM];
    double  (*FLP)[NBEM2][NIDM], (*FLP2)[NBEM2][NIDM], (*FLPP)[NBEM2][NIDM];
    double  *AGEL, *AGEU, *BAGE, *RBAGE, (*AGE)[NIDM];
    double  (*AGE2)[NIDM],(*AGEP)[NIDM];
	bool *LEDEP, (*LDIT)[NIDM], (*LDIP)[NBEM2][NIDM], (*LFLT)[NIDM], (*LFLP)[NBEM2][NIDM], (*LAGEA)[NIDM];
	int *IDCUT, *NE;
	bool *LLE, *LLAGE;
    int *NAGE, *NID;

}CIMDET;

typedef struct{
	double *EL, *EU, *BSE, *RBSE, *EDEP, *EDEP2, (*DET)[NIDM];
    int  *NE, *NID;
	bool *LLE;
}CENDET;

typedef struct{
	double (*DOSE)[NDYM][NDXM], (*DOSE2)[NDYM][NDXM], (*DOSEP)[NDYM][NDXM];
	int (*LDOSE)[NDYM][NDXM], *KDOSE;
} CDOSE1;

typedef struct{
	double *DDOSE, *DDOSE2, *DDOSEP;
	int *LDDOSE;

}CDOSE2;

typedef struct{
	double *DXL, *DXU, *BDOSE, *RBDOSE;
	int *NDB;
}CDOSE3;

typedef struct{
	double (*VMASS)[NDYM][NDXM];

}CDOSE4;


typedef struct{
	double *ELAST1, *ELAST2;
	int *MHINGE, *KSOFTE, *KSOFTI, *KDELTA;
}CJUMP1;

typedef struct{
	double *ES, *XS, *YS, *ZS, *US , *VS, *WS, *WGHTS, *SP1S, *SP2S, *SP3S, *PAGES;
    int *KS, *IBODYS, *MS, (*ILBS)[5], *IPOLS, *NSEC;
}SECST;

typedef struct{
	double *P, *ST, *DST, *DSR, *W1, *W2, *T1, *T2;
}CJUMP0;

typedef struct{
	int *ILBA;
}CHIST;






void extrairString(char destino[], char origem[], int inicio, int qtde);
bool comentario(char BLINE[]);
void imprimirKSURF(FILE* IW, int &KB);
void imprimirKFLAG(FILE* IW, int &KB);
void imprimirKBODY(FILE* IW, int &KB);
void imprimirKDGHT(FILE* IW, int &KB);



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


extern "C" {
	
void transfqsurf_(double *AXX, double *AXY, double *AXZ, double *AYY, double *AYZ, double *AZZ, double *AX, double *AY, double *AZ, double *A0, int *NSURF, int *KPLANE);

void transfqtree_(int *NBODYS, int *KMOTH, int (*KDGHT)[NB], int (*KSURF)[NB], int (*KFLAG)[NB], int *KSP, int *NWARN);

void transftrack_mod_(double *E, double *X, double *Y, double *Z, double *U, double *V, double *W, double *WGHT, double *SP1, double *SP2, double *SP3, double *PAGE, int *KPAR, int *IBODY, int *MAT, int *ILB, int *IPOL, bool *LAGE);
  	
void transfpenelope_mod_( double (*EABS)[3], double *C1, double *C2, double *WCC, double *WCR, double *DEN, double *RDEN, double *E0STEP, double *DESOFT, double *SSOFT, int *NMS, int *NEGP, int *NMAT);
 	
void transfpengeom_mod_( char (*BALIAS)[5], double *DSTOT, int *MATER, int *KDET, int *KSLAST, int *NBODY, bool *LVERB);

void transfqbody_(int (*KBODY)[NB], int *KBOMO);

void transfcecutr_(double *ECUTR);

void transfcsgawr_(int *ISGAW);

void transfcersec_(int *IERSEC);

void transfcegrid_(double *EMIN, double *EL, double *EU, double *ET, double *DLEMP, double *DLEMP1, double *DLFC, double *XEL, double *XE, double *XEK, int *KE);

void transfcesi0_(double (*XESI)[NRP], int *IESIF, int *IESIL, int *NSESI, int *NCURE);

void transfcpsi0_(double (*XPSI)[NRP], int *IPSIF, int *IPSIL, int *NSPSI, int *NCURP);

void transfcadata_(double *ATW, double *EPX, double *RSCR, double *ETA, double (*EB)[99], double (*ALW)[99], double (*CP0)[99], int (*IFI)[99], int (*IKS)[99], int *NSHT,  char (*LASYMB)[2]);


void transfcgph00_(int *IPHF, int *IPHL, int *NPHS, int *NCUR, double *EPH, double (*XPH)[NTP]);


void transfcrelax_(double *P, double *ET, double *F, int *IS0, int *IS1, int *IS2, int (*IFIRST)[99], int (*ILAST)[99], int *NCUR, int *KS, int *MODER );

void transfcritaa_(	double *XA, double *AA, double *BA, double *FA, int *IA, int *NPM1A);

void transfcrndg3_(double *X, double *A, double *B, double *F, int *KA, int *NPM1);

void transfcrita_(double *XT, double *PAC, double *DPAC, double *A, double *B, int *IL, int *IU, int *NPM1);
// double *QTI, double *PACI, double *DPACI, double *AI, double *BI, int *ITLI, int *ITUI, int *NPM1I,
// double *XTI);

void transfcgra00_(double *FACTE, double *Q2MAX, int *MM, int *MOM);

void transfcgra01_(double (*FF)[MAXMAT], double *ERA, double (*XSRA)[MAXMAT], int *IED, int *IEU, int *NE);

void transfcgra02_(double *QQ, double (*AR)[MAXMAT], double (*BR)[MAXMAT], double (*CR)[MAXMAT], double (*DR)[MAXMAT], double *FF0, double *QQM);

void transfcgra03_(double  (*QRA)[NP2], double (*PRA)[NP2], double (*DPRA)[NP2], double (*ARA)[NP2],
			double (*BRA)[NP2],double (*PMAX)[NEGP], int (*ITLRA)[NP2],int (*ITURA)[NP2]);


void transfcritan_(double *CNORM);

void transfcompos_(double (*STF)[MAXMAT], double *ZT, double *AT, double *RHO, double *VMOL, int (*IZ)[MAXMAT], int *NELEM);

void transfcrange_(double (*RANGE)[MAXMAT][3], double (*RANGEL)[MAXMAT][3]);

void transfcein_(double *EXPOT, double *OP2, double (*F)[MAXMAT],double (*UI)[MAXMAT], double (*WRI)[MAXMAT], int (*KZ)[MAXMAT], int (*KS)[MAXMAT], int *NOSC);

void transfceintf_(double *T1EI, double *T2EI , double *T1PI, double *T2PI);

void transfcein00_(double *SEH0, double *SEH1, double *SEH2 , double *SES0, double *SES1, double *SES2, double *SET0, double *SET1, double *SET2);

void transfcpin00_(double *SPH0, double *SPH1, double *SPH2, double *SPS0, double *SPS1, double *SPS2, double *SPT0, double *SPT1, double *SPT2);

void transfceinac_(double (*EINAC)[NEGP][MAXMAT], int (*IEIN)[MAXMAT], int *NEIN);

void transfcesiac_(double (*ESIAC)[NEGP][MAXMAT], int (*IESI)[MAXMAT], int *NESI);

void transfcesin_(double (*XSEIN)[NEGP], double (*XSESI)[NEGP], int *ISIE);

void transfcpinac_(double (*PINAC)[NEGP][MAXMAT], int (*IPIN)[MAXMAT], int *NPIN);

void transfcpsiac_(double (*PSIAC)[NEGP][MAXMAT], int (*IPSI)[MAXMAT], int *NPSI);

void transfcpsin_(double (*XSPIN)[NEGP], double (*XSPSI)[NEGP], int *ISIP);

void transfcgco_(double (*FCO)[MAXMAT], double (*UICO)[MAXMAT], double (*FJ0)[MAXMAT], double (*PTRSH)[MAXMAT], int (*KZCO)[MAXMAT], int (*KSCO)[MAXMAT], int *NOSCCO);

void transfceimfp_(double (*SEHEL)[MAXMAT], double(*SEHIN)[MAXMAT], double(*SEISI)[MAXMAT], double(*SEHBR)[MAXMAT],double (*SEAUX)[MAXMAT], 
          double (*SETOT)[MAXMAT], double (*CSTPE)[MAXMAT], double (*RSTPE)[MAXMAT], double (*DEL)[MAXMAT], double (*W1E)[MAXMAT],
          double (*W2E)[MAXMAT], double (*DW1EL)[MAXMAT],  double(*DW2EL)[MAXMAT], double(*RNDCE)[MAXMAT],double (*AE)[MAXMAT], double (*BE)[MAXMAT],
    	   double (*T1E)[MAXMAT],double (*T2E)[MAXMAT]);

void transfclas1e_(double (*TSTPE)[MAXMAT], double(*TSTRE)[MAXMAT], double(*TRL1E)[MAXMAT], double(*TRL2E)[MAXMAT]);

void transfcpimfp_(double (*SPHEL)[MAXMAT], double(*SPHIN)[MAXMAT], double(*SPISI)[MAXMAT], double(*SPHBR)[MAXMAT],double (*SPAN)[MAXMAT], 
           double (*SPAUX)[MAXMAT], double (*SPTOT)[MAXMAT], double (*CSTPP)[MAXMAT], double (*RSTPP)[MAXMAT], double (*W1P)[MAXMAT],
           double (*W2P)[MAXMAT], double (*DW1PL)[MAXMAT], double (*DW2PL)[MAXMAT], double (*RNDCP)[MAXMAT], double (*AP)[MAXMAT], double (*BP)[MAXMAT],
           double (*T1P)[MAXMAT], double (*T2P)[MAXMAT]);

void transfclas1p_(double (*TSTPP)[MAXMAT], double (*TSTRP)[MAXMAT], double (*TRL1P)[MAXMAT], double (*TRL2P)[MAXMAT]);

void transfceel00_(double *EJT, double  *XE0,double  *XE1, double *XE2, double *XP0, double *XP1, double *XP2, double *T1E0, double *T2E0, double *T1P0, 
            double *T2P0, double *EJTL, double *FJL, double *A, double *B, double *C, double *D);

void transfcbryld_(double (*EBRY)[MAXMAT], double (*PBRY)[MAXMAT]);

void transfcgimfp_(double (*SGRA)[MAXMAT], double (*SGCO)[MAXMAT], double (*SGPH)[MAXMAT], double (*SGPP)[MAXMAT], double (*SGAUX)[MAXMAT]);

void transfcgph01_(double *ER, double *XSR, int *NPHD);

void transfcgpp01_(double (*TRIP)[MAXMAT]);

void transfcebr_(double *WB, double (*PBCUT)[MAXMAT], double (*WBCUT)[MAXMAT], double (*PDFB)[NEGP][MAXMAT], double (*DPDFB)[NEGP][MAXMAT], double (*PACB)[NEGP][MAXMAT], double *ZBR2);

void transfcebr01_(double *EBT, double (*XS)[NBE], double *TXS, double *X, double *Y);

void transfcebr02_(double (*P0)[NEGP][MAXMAT]);

void transfcbrang_(double *BET, double *BK, double(*BP1)[21][6][MAXMAT], double (*BP2)[21][6][MAXMAT], double *ZBEQ);

void transfcein01_(double *EI, double *EE, double *CPS, double *AMOL, int *MOM);
	
void transfcsumga_(int *IERGA, int *NCALL);

void transfcpin01_(double *EI, double *CPS, double *BHA1, double *BHA2, double *BHA3, double *BHA4, int *MOM);

void transfcdcsep_(double *ETS, double *ETL, double *TH, double *THR, double *XMU , double *XMUL, double *ECS, double *ETCS1 , double *ETCS2, double (*EDCS)[NE],
     	   double *PCS, double *PTCS1,  double *PTCS2, double (*PDCS)[NE], double *DCSI, double *DCSIL, double *CSI, double *TCS1I, double *TCS2I);

void transfceeldb_(	double (*XSE)[NEGP][NP], double (*PSE)[NEGP][NP], double(*ASE)[NEGP][NP], double (*BSE)[NEGP][NP],
    int (*ITLE)[NEGP][NP], int (*ITUE)[NEGP][NP]);

void transfcpeldb_(double (*XSP)[NEGP][NP],double  (*PSP)[NEGP][NP],double (*ASP)[NEGP][NP],double  (*BSP)[NEGP][NP],
    int(*ITLP)[NEGP][NP], int (*ITUP)[NEGP][NP]);

void transfcelsep_(double *EELMAX, double *PELMAX, double (*RNDCED)[MAXMAT], double (*RNDCPD)[MAXMAT]);

void transfcgpp00_(double *ZEQPP, double (*F0)[MAXMAT], double *BCB);

void transfrseed_(int *ISEED1, int *ISEED2);

void transfctitle_(char *TITLE, char *TITLE2);

void transfcdate_(char *DATE23);

void transfcspgeo_(double *DSMAX, double (*EABSB)[3]);

void transfcforci_(double (*WLOW)[NBV], double(*WHIG)[NBV], bool (*LFORCE)[NBV]);

void transfcxrspl_(int *IXRSPL, int *ILBA, bool *LXRSPL);

void transfcsour0_(double *CTHL, double *DCTH, double *PHIL, double *DPHI, int *KPARP, int *JOBEND, bool *LSCONE, bool *LGPOL, bool *LPSF);

void transfcsour1_(double *E0, double *EPMAX, double *SP10, double *SP20, double *SP30);

void transfcsour2_(double *ESRC, double *PSRC, int *IASRC, double *FSRC, bool *LSPEC);

void transfcsour3_(double *SX0, double *SY0, double *SZ0, double *SSX, double *SSY, double *SSZ, int *IXSBOD, bool *LEXSRC, bool *LEXBD);

void transfcsour4_(double *WGMIN, double *RWGMIN, double *WGMAX, double *RLREAD, int *IPSFI, int *NPSF, int *NPSN, int *NSPLIT, int *KODEPS);

void transfcsour5_(char (*PSFI)[20]);

void transfcnt0_(double *PRIM , double *PRIM2, double *DPRIM, double (*SEC)[3], double (*SEC2)[3], double (*DSEC)[3], 
			    double *AVW, double *AVW2, double *DAVW, double *AVA, double *AVA2, double *DAVA, double *AVE, double *AVE2, double *DAVE);

void transfcnt1_(double *TDEBO, double *TDEBO2, double *DEBO);

void transfcnt2_(double *SHIST, int *NSEB);

void transfcnt3_(double (*SEDS)[3], double (*SEDS2)[3], double *DSDE, double *RDSDE, int *NSDE);

void transfcnt4_(double *RLAST, double *RWRITE, int *IDCUT, int (*KKDI)[NIDM], int *IPSF, int *NID, int *NPSFO, int *IPSFO);

void transfcnt5_(double *DEDE, int *KBDE, int *NED);

void transfcnt6_(bool *LDOSEM);

void transfcdump_(bool *LDUMP, char *PFILED);

void transfcntrl_(double *TSIM, double *TSEC, double *TSECA, double *TSECAD, double *CPUT0, double *DUMPP, double *DSHN, double *SHN, int *N);

void transfcgcone_(double *CPCT, double *CPST, double *SPCT, double *SPST, double *SPHI, double *CPHI, double *STHE, double *CTHE, double *CAPER);

void transfcenang_(	double *EL, double *EU, double *THL, double *THU, double *BSE, double *RBSE, double *BSTH, double *RBSTH, double *BSPH, double *RBSPH,
	double (*PDE)[2][3], double (*PDE2)[2][3], double (*PDEP)[2][3],
	double (*PDA)[NBTHM][3], double (*PDA2)[NBTHM][3], double (*PDAP)[NBTHM][3],
	bool (*LPDE)[2][3], bool (*LPDA)[NBTHM][3],
	int *NE, int *NTH, int *NPH,
	bool *LLE, bool *LLTH);


void transfcimdet_(double *EL, double *EU, double *BSE, double *RBSE,
	double (*ET)[NIDM], double *EDEP, double *EDEP2, double *EDEPP,
    double  (*DIT)[NIDM], double (*DIT2)[NIDM], double (*DITP)[NIDM],
    double  (*DIP)[NBEM2][NIDM], double (*DIP2)[NBEM2][NIDM], double (*DIPP)[NBEM2][NIDM],
    double  (*FLT)[NIDM], double (*FLT2)[NIDM], double (*FLTP)[NIDM],
    double  (*FLP)[NBEM2][NIDM], double (*FLP2)[NBEM2][NIDM], double (*FLPP)[NBEM2][NIDM],
    double  *AGEL, double *AGEU, double *BAGE, double *RBAGE, double (*AGE)[NIDM],
    double  (*AGE2)[NIDM], double (*AGEP)[NIDM],
	bool *LEDEP, bool (*LDIT)[NIDM], bool (*LDIP)[NBEM2][NIDM], bool (*LFLT)[NIDM], bool (*LFLP)[NBEM2][NIDM], bool (*LAGEA)[NIDM],
	int *IDCUT, int *NE,
	bool *LLE, bool *LLAGE,
    int *NAGE, int *NID);

void transfcendet_(double *EL, double *EU, double *BSE, double *RBSE, double *EDEP, double *EDEP2, double (*DET)[NIDM],
    int  *NE, int *NID, bool *LLE);

void transfcdose1_(double (*DOSE)[NDYM][NDXM], double (*DOSE2)[NDYM][NDXM], double (*DOSEP)[NDYM][NDXM], int (*LDOSE)[NDYM][NDXM], int *KDOSE);

void transfcdose2_(double *DDOSE, double *DDOSE2, double *DDOSEP, int *LDDOSE);

void transfcdose3_(double *DXL, double *DXU, double *BDOSE, double *RBDOSE, int *NDB);

void transfcdose4_(double (*VMASS)[NDYM][NDXM]);

void transfcjump1_(double *ELAST1, double *ELAST2, int *MHINGE, int *KSOFTE, int *KSOFTI, int *KDELTA);

void transfsecst_(double *ES, double *XS, double *YS, double *ZS, double *US , double *VS, double *WS, double *WGHTS, double *SP1S, double *SP2S, double *SP3S, double *PAGES,
    int *KS, int *IBODYS, int *MS, int (*ILBS)[5], int *IPOLS, int *NSEC);

void transfcjump0_(double *P, double *ST, double *DST, double *DSR, double *W1, double *W2, double *T1, double *T2);

void transfchist_(int *ILBA);



void einat_(double &E, double &WCCM, double &XH0, double &XH1, double &XH2, double &XS0, double &XS1, double &XS2, double &XT1, double &XT2, double &DELTA, int *M); //REVISTO


void ebrar_(double EBT, int M, int IRD, int IWR, int INFO);


void locate2_();
    
void fsurf2_(int &KS, double &A, double &B, double &C);

void step2_(double *DS, double *DSEF, int *NCROSS);

void stepsi2_(int &KB, int &NSC);

//void stepsi2_(int &KB, double *S, int *IS, int &NSC);

void stepsi_(int &KB, double *S, int *IS, int &NSC);

void steplb2_(int &KB, int &IERR);

void geomin2_(double *PARINP, int *NPINP, int *NMATG, int *NBOD, FILE *IRD, FILE *IWR);

void rotshf2_(double &OMEGA, double &THETA, double &PHI, double &DX, double &DY, double &DZ, double &AXX, double &AXY, double &AXZ, double &AYY, double &AYZ, double &AZZ, double &AX, double &AY, double &AZ, double &A0);

void peinit2_(double *EMAX, int &NMATER, FILE *IWR, int *INFO, char (*PMFILE)[100]);

void egrid2_(double &EMINu, double *EMAXu);

void esia02_();

void psia02_();

void gpha02_();

void relax02_();

void rndg302_();

double rndg3f2_(double X);

void rita02_(int &PDF, double &XLOW, double &XHIGH, int &N, int &NU, double &ERRM, FILE *IWR);

void ritai02_(int &PDF, double &XLOW, double &XHIGH, int &N, int &NU, double &ERRM, FILE *IWR);

void pematr2_(int *M, FILE *IRD, FILE *IWR, int *INFO);


void irnd02_(double *W, double *F, int *K, int *N);


void pematr_(int *M, uintptr_t IRD, uintptr_t IWR, int *INFO);

void egrid_(double &EMINu, double *EMAXu, FILE *IW);

void esia0_();

void psia0_();

void gpha0_();

void relax0_();

void rndg30_();

double rndg3f_(double &X);

void rita0_(double &XLOW, double &XHIGH, int &N, int &NU, double &ERRM, int *IWR);

void ritai0_(double &XLOW, double &XHIGH, int &N, int &NU, double &ERRM, int *IWR);

void irnd0_(double *W, double *F, int *K, int *N);

void relaxr2_(FILE *IRD, FILE *IWR, int *INFO);

void esiar2_(int *M, FILE *IRD, FILE *IWR,int *INFO);

void findi2_(double *X, double *XC, int *N, int *J);

void findi_(double *X, double &XC, int &N, int &J);


void psiar2_(int *M, FILE *IRD, FILE *IWR, int *INFO);

void ebrar2_(double *WCRM, int *M, FILE *IRD, FILE *IWR, int *INFO);

void spline2_(double *X, double *Y, double *A, double *B, double *C, double *D, double S1, double SN, int const N);

void spline_(double *X, double *Y, double *A, double *B, double *C, double *D, double &S1, double &SN, int &N);


void rlpac2_(double *X, double *PDF, double *PAC, int *NP);


double rlmom2_(double *X, double *FCT, double XC, int NP, int MOM);

void braar2_(int *M, FILE *IRD, FILE *IWR, int *INFO);


void einat2_(double &E, double &WCCM, double &XH0, double &XH1, double &XH2, double &XS0, double &XS1, double &XS2, double &XT1, double &XT2, double &DELTA, int *M);


void einat12_(double &E, double &UK, double &WK, double &DELTA, double &WCCM, double &H0, double &H1, double &H2, double &S0, double &S1, double &S2, double &R0, double &R1, double &R2);

double sumga2_(int &funcao, double &XL, double &XU, double &TOL);

double einads2_(double RMU);

void pinat2_(double &E, double &WCCM, double &XH0, double &XH1, double &XH2, double &XS0, double &XS1, double &XS2, double &XT1, double &XT2, double &DELTA, int *M);

void pinat_(double &E, double &WCCM, double &XH0, double &XH1, double &XH2, double &XS0, double &XS1, double &XS2, double &XT1, double &XT2, double &DELTA, int *M);


void pinat12_(double &E, double &UK, double &WK, double &DELTA, double &WCCM, double &H0, double &H1, double &H2, double &S0, double &S1, double &S2, double &R0, double &R1, double &R2);

void pinat1_(double &E, double &UK, double &WK, double &DELTA, double &WCCM, double &H0, double &H1, double &H2, double &S0, double &S1, double &S2, double &R0, double &R1, double &R2);


double pinads2_(double RMU);

void ebrat2_(double &E, double &WCRM, double &XH0, double &XH1, double &XH2, double &XS1, double &XS2, int *M);

void sinteg2_(double *X, double *A, double *B, double *C, double *D,  double &XL, double &XU, double &SUM, int &N);

double rmomx2_(double *X, double *PDF, double XD, double XU, int NP, int MOM);

void pbrat2_(double &E, double &WCRM, double &XH0, double &XH1, double &XH2, double &XS1, double &XS2, int *M);

void panat2_(double &E, double &TXS);

void eelar2_(int *M, FILE *IRD, FILE *IWR, int *INFO);

void eela02_(double &XSO, double &XS1, double &XS2, double &XS0H, double &A, double &B, double &RNDC, double &XS1S, double &XS2S);

void eeldr2_(int *M, FILE *IRD, FILE *IWR, int *INFO);

void dcsel02_(double E, int IELEC);

double dcsel2_(double RMU);

void ritam2_(double XD, double XU, double &XM0, double &XM1, double &XM2);

void graar2_(int *M, FILE *IRD, FILE *IWR, int *INFO);

double graaf22_(double Q2);

void slag62_(double H, double *Y, double *S, int N);

void gphar2_(int *M, FILE *IRD, FILE *IWR, int *INFO);

void merge22_(double *X1, double *Y1, double *X2, double *Y2, double *XM, double *YM, int &N1, int &N2, int &N);

void sort22_(double *X, double *Y, int N);

void graati2_(double E, double &ECS, int *M);

void gppa02_(int *M);

void pmrdr2_();

void gcone02_(double THETA, double PHI, double ALPHA);

void enang02_(double &EMIN, double &EMAX, int &NBE, int &NBTH, int &NBPH, FILE *IWR);

void imdet02_(double &EMIN, double &EMAX, int &NBE, double &AGEMIN, double &AGEMAX, int &NBAGE, int &ICUT, char *FNSPC, char *FNFLU, char *FNAGE, int &ID, FILE *IWR);

void endet02_(double &EMIN, double &EMAX, int &NB, char *FNSPC, int &ID, FILE *IWR);

void dose02_(double &XL, double &XU,double &YL,double &YU,double &ZL,double &ZU, int &NBX, int &NBY, int &NBZ, int &IDOSE, FILE *IWR);

void rand02_(int N);

void enangr2_(ifstream &IRD, FILE *IWDUMP);

void imdetr2_(ifstream &IRD, FILE *IWDUMP);

void endetr2_(ifstream &IRD, FILE *IWDUMP);

void doser2_(ifstream &IRD, FILE *IWDUMP);

void tabelas_();

void shower2_();

void cleans2_();

double rand2_(double DUMMY);

void gcone2_(double &UF, double &VF, double &WF);

void panar2_(double &ECUT);

void direct2_(double &CDT, double &DF, double &U, double &V, double &W);

void stores2_(double &EI, double &XI, double &YI, double &ZI, double &UI, double &VI, double &WI, double &WGHTI, int &KPARI, int *ILBI, int &IPOLI);

void sdose2_(double &DEP, double &XD, double &YD, double &ZD, int &MATC, int &N);

void simdet2_(int &N, int &ID);

void jump2_(double &DSMAX, double &DS);

void start2_();

double rndg32_();

void eimfp2_(int &IEND);

void pimfp2_(int &IEND);

void gimfp2_();

void fimdet2_(int &N, int &ID, double &DSEF);

void fimdes2_(int &N, int &ID, double &EI, double &DECSD, double &DSEF );

void eela2_(double &A, double &B, double &RNDC, double &RMU);

void eeld2_(double &RNDC, double &RMU);

void eina2_(double &E, double &DELTA, double &DE, double &EP, double &CDT, double &ES, double &CDTS, int &M, int &IOSC);

void knock2_(double &DE, int &ICOL);

void ebra2_(double &E, double &W, int &M);

void ebraa2_(double &E, double &DE, double &CDT, int &M);

void esia2_(double &E, double &DELTA, double &DE, double &EP, double &CDT, double &ES, double &CDTS, int &M, int &IZZ, int &ISH);

void relax2_(int &IZ, int &IS);

void eaux2_();

void graa2_(double &E, double &CDT, int &IEFF, int &M);

void dirpol2_(double &CDT, double &DF, double &CONS, double &SP1, double &SP2, double &SP3, double &U, double &V, double &W);

void gcoa2_(double &E, double &DE, double &EP, double &CDT, double &ES, double &CDTS, int &M, int &IZZ, int &ISH);

void  gpha2_(double &ES, int &IZZ, int &ISH);

void sauter2_(double &ES, double &CDTS);

void schiff2_(double &B, double &G1, double &G2);

void gppa2_(double &EE, double &CDTE, double &EP, double &CDTP, int &IZZ, int &ISH);

void gaux2_();

void peld2_(double &RNDC, double &RMU);

void pina2_(double &E, double &DELTA, double &DE, double &EP, double &CDT, double &ES, double &CDTS, int &M, int &IOSC);

void psia2_(double &E, double &DELTA, double &DE, double &EP, double &CDT, double &ES, double &CDTS, int &M, int &IZZ, int &ISH);

void pana2_(double &E, double &E1, double &CDT1, double &E2, double &CDT2, int &M);

void paux2_();

void tenang2_(int &IEXIT, int &N);

void sendet2_(double &ED, int &ID);

void secpar2_(int &LEFT);

void pmwrt2_(int ICLOSE);

void enangd2_(FILE *IWR);

void imdetd2_(FILE *IWR);

void endetd2_(FILE *IWR);

void dosed2_(FILE *IWR);

void enangw2_(double &SHN);

void dosew2_(double &SHN, double &TSIM, FILE *IWR);

void timer2_(double &SEC);


double cputim2_();

void inicializarStructs();

void memoryFree();

}


void transfqsurf_(double *AXX, double *AXY, double *AXZ, double *AYY, double *AYZ, double *AZZ, double *AX, double *AY, double *AZ, double *A0, int *NSURF, int *KPLANE){
	QSURF_.AXX = AXX;
	QSURF_.AXY = AXY;
	QSURF_.AXZ = AXZ;
	QSURF_.AYY = AYY;
	QSURF_.AYZ = AYZ;
	QSURF_.AZZ = AZZ;
	QSURF_.AX = AX;
	QSURF_.AY = AY;
	QSURF_.AZ = AZ;
	QSURF_.A0 = A0;
	QSURF_.NSURF = NSURF;
	QSURF_.KPLANE = KPLANE;			
}
       
void transfqtree_(int *NBODYS, int *KMOTH, int (*KDGHT)[NB], int (*KSURF)[NB], int (*KFLAG)[NB], int *KSP, int *NWARN){
	QTREE_.NBODYS = NBODYS;
	QTREE_.KMOTH = KMOTH;
	QTREE_.KDGHT = KDGHT;
	QTREE_.KSURF = KSURF;
	QTREE_.KFLAG = KFLAG;
	QTREE_.KSP = KSP;
	QTREE_.NWARN = NWARN;
} 

void transftrack_mod_(double *E, double *X, double *Y, double *Z, double *U, double *V, double *W, double *WGHT, double *SP1, double *SP2, double *SP3, double *PAGE, int *KPAR, int *IBODY, int *MAT, int *ILB, int *IPOL, bool *LAGE){

/*

	M�DULO TRACK_mod
 
   **** Vari�veis de part�cula TRACK (a serem inicializadas antes de chamar
         sub-rotina START).
      SALVE ! Salva todos os itens do m�dulo.
 
   ---- Energia, posi��o, dire��o e peso.
      DUPLA PRECIS�O :: E, X, Y, Z, U, V, W, WGHT
   ---- Tipo de part�cula, corpo atual e material.
      INTEGER * 4 :: KPAR, IBODY, MAT
   ---- Sinalizadores de hist�rico de part�culas.
      INTEGER * 4, DIMENSION (5) :: ILB
 
   **** Polariza��o de f�tons.
   ---- F�tons polarizados se IPOL = 1, caso contr�rio f�tons n�o polarizados.
      INTEGER * 4 :: IPOL = 0
   ---- Par�metros de Stokes.
      DUPLA PRECIS�O :: SP1, SP2, SP3
 
   **** A idade das part�culas (tempo decorrido desde o in�cio da simula��o)
         � registrado quando LAGE = .TRUE.
      LOGICAL :: LAGE = .FALSE.
*/

	TRACK_mod_.E = E;
	TRACK_mod_.X = X;
	TRACK_mod_.Y = Y;
	TRACK_mod_.Z = Z;
	TRACK_mod_.U = U;
	TRACK_mod_.V = V;
	TRACK_mod_.W = W;
	TRACK_mod_.WGHT = WGHT;
	TRACK_mod_.SP1 = SP1;
	TRACK_mod_.SP2 = SP2;
	TRACK_mod_.SP3 = SP3;
	TRACK_mod_.PAGE = PAGE;
	TRACK_mod_.KPAR = KPAR;
	TRACK_mod_.IBODY = IBODY;
	TRACK_mod_.MAT = MAT;
	TRACK_mod_.ILB = ILB;
	TRACK_mod_.IPOL = IPOL;
	TRACK_mod_.LAGE = LAGE;
}
 	
 void transfpenelope_mod_( double (*EABS)[3], double *C1, double *C2, double *WCC, double *WCR, double *DEN, double *RDEN, double *E0STEP, double *DESOFT, double *SSOFT, int *NMS, int *NEGP, int *NMAT){
 		


/*
	M�DULO PENELOPE_mod
      SAVE ?  ! Salva todos os itens do m�dulo.
 
   **** N�mero m�ximo de materiais na geometria.
 
      INTEGER * 4, PARAMETER :: MAXMAT = 10
 
   **** Par�metros de simula��o (devem ser definidos antes de chamar o
         sub-rotina de inicializa��o PEINIT).
 
   ---- Energias de absor��o, EABS (KPAR, MAT).
      DUPLA PRECIS�O, DIMENS�O (3, MAXMAT) :: EABS = 50,0D0
   ---- Par�metros de transporte de el�tron / p�sitron.
      DUPLA PRECIS�O, DIMENS�O (MAXMAT) :: C1 = 0,01D0, C2 = 0,01D0
      DUPLA PRECIS�O, DIMENS�O (MAXMAT) :: WCC = 1.0D2, WCR = 1.0D2
 
   **** Tamanho da pilha secund�ria (n�mero m�ximo de part�culas que
         pode ser armazenado).
      INTEGER * 4, PARAMETER :: NMS = 1000
 
   **** Interpola��o de energia, n�mero de pontos da grade.
      INTEGER * 4, PARAMETER :: NEGP = 200
 
   **** Informa��es globais sobre o sistema de materiais (definido por
         sub-rotina PEINIT).
   ---- N�mero de materiais presentes.
      INTEGER * 4 :: NMAT
   ---- Densidades de materiais e seus rec�procos.
      DUPLA PRECIS�O, DIMENS�O (MAXMAT) :: DEN = 1.0D0, RDEN = 1.0D0
 
   **** Par�metros de desacelera��o de dobradi�a aleat�ria (sa�da de sub-rotinas
         JUMP e KNOCK).
           E0STEP ... energia no in�cio do segmento,
           DESOFT ... perda de energia ao longo da etapa,
           SSOFT .... pot�ncia de parada efetiva, = DESOFT / step_length.
      DUPLA PRECIS�O :: E0STEP, DESOFT, SSOFT
 
      END MODULE PENELOPE_mod

	  
	  **** Estado atual na simula��o de classe II.
        MHINGE = 0 (1) antes (depois) da dobradi�a.
     COMMON / CJUMP1 / ELAST1, ELAST2, MHINGE, KSOFTE, KSOFTI, KDELTA



*/




    PENELOPE_mod_.EABS = EABS;
	PENELOPE_mod_.C1 = C1;
	PENELOPE_mod_.C2 = C2;
	PENELOPE_mod_.WCC = WCC;
	PENELOPE_mod_.WCR = WCR;
	PENELOPE_mod_.DEN = DEN;
	PENELOPE_mod_.RDEN = RDEN;
	PENELOPE_mod_.E0STEP = E0STEP;
	PENELOPE_mod_.DESOFT = DESOFT;
	PENELOPE_mod_.SSOFT = SSOFT;
	PENELOPE_mod_.NMS = NMS;
	PENELOPE_mod_.NEGP = NEGP;
	PENELOPE_mod_.NMAT = NMAT;	
	
}

void transfpengeom_mod_( char (*BALIAS)[5], double *DSTOT, int *MATER, int *KDET, int *KSLAST, int *NBODY, bool *LVERB){
      	
	/*
	
	O m�dulo TRACK_mod cont�m vari�veis de trilha de part�culas (posi��o,
   dire��o, corpo atual e material), que s�o compartilhados com o
   programa principal de dire��o. Quando PENGEOM � usado com PENELOPE, este
   m�dulo � comentado, porque � definido na lista de fontes
   'penelope.f'. Para usar o PENGEOM com outros c�digos Monte Carlo, descomente
   as quatro linhas come�ando com 'c>' abaixo, incluem a instru��o 'USE
   TRACK_mod 'em todos os subprogramas que definem ou modificam a posi��o ou
   dire��o da part�cula, e compilar o arquivo 'pengeom.f' antes
   quaisquer outros arquivos de origem que usam as sub-rotinas atuais.
 
   >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> >>>>>>>>>>>>>>>>>>>>>>>> >>>>>>>>>>>>>>>>>>>>>
    M�DULO TRACK_mod! Vari�veis ??de part�cula TRACK.
    SALVE ?   ! Salva todos os itens do m�dulo.
    DUPLA PRECIS�O :: X, Y, Z, U, V, W! Posi��o e dire��o.
    INTEGER * 4 :: IBODY, MAT! Corpo e material atuais.
   FIM DO M�DULO TRACK_mod
   <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< <<<<<<<<<<<<<<<<<<<<<<


   >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> >>>>>>>>>>>>>>>>>>>>>>>> >>>>>>>>>>>>>>>>>>>>>
      M�DULO PENGEOM_mod
 
   **** Par�metros de defini��o de geometria e grandezas de E / S.
 
      SAVE Salva todos os itens do m�dulo.
   ---- Tamanhos de matriz de geometria.
         N�mero m�ximo de superf�cies, corpos e elementos limitantes.
      INTEGER * 4, PARAMETER :: NS = 10000, NB = 5000, NXG = 250
         N�mero de corpos no sistema de materiais (fornecido pelo PENGEOM).
      INTEGER * 4 :: NBODY
 
   ---- Aliases de corpo (r�tulos de usu�rio).
      PERSONAGEM * 4 :: BALIAS (NB) = ''
 
   ---- Materiais do corpo. MATER (KB) � o material no corpo KB.
      INTEGER * 4 :: MATER (NB) = 0
 
   ---- Defini��o do detector.
         KDET (KB) = ID se o corpo KB fizer parte do ID do detector.
      INTEGER * 4 :: KDET (NB) = 0
 
   **** Mensagens de aviso para disparos acidentais ou erros de arredondamento
         s�o emitidos quando LVERB = .TRUE.
 
      LOGICAL :: LVERB = .FALSE.
 
   **** Recursos da �ltima etapa (sa�da da sub-rotina STEP).
 
   ---- Comprimento do caminho percorrido, incluindo segmentos em volumes vazios.
      DUPLA PRECIS�O :: DSTOT
   ---- R�tulo da �ltima superf�cie cruzada pela part�cula antes
         inserir um corpo material (definido apenas quando NCROSS / = 0).
      INTEGER * 4 :: KSLAST
 
      END MODULE PENGEOM_mod
	
	
	*/
	
	
	
	PENGEOM_mod_.BALIAS = BALIAS;
	PENGEOM_mod_.DSTOT = DSTOT;
	PENGEOM_mod_.MATER = MATER;
	PENGEOM_mod_.KDET = KDET;
	PENGEOM_mod_.KSLAST = KSLAST;
	PENGEOM_mod_.NBODY = NBODY;
	PENGEOM_mod_.LVERB = LVERB;

}

void transfqbody_(int (*KBODY)[NB], int *KBOMO){
	
	QBODY_.KBODY = KBODY;
	QBODY_.KBOMO = KBOMO;
	
}

     
void transfcecutr_(double *ECUTR){
	
	CECUTR_.ECUTR = ECUTR;
	
	
	
}

void transfcsgawr_(int *ISGAW){
	
	CSGAWR_.ISGAW = ISGAW;
	
}

void transfcersec_(int *IERSEC){
	
	CERSEC_.IERSEC = IERSEC;
	
}     
     
  
void transfcegrid_(double *EMIN, double *EL, double *EU, double *ET, double *DLEMP, double *DLEMP1, double *DLFC, double *XEL, double *XE, double *XEK, int *KE){
	
	CEGRID_.EMIN = EMIN;
	CEGRID_.EL = EL; 
	CEGRID_.EU = EU;
	CEGRID_.ET = ET;
    CEGRID_.DLEMP = DLEMP;
	CEGRID_.DLEMP1 = DLEMP1;
	CEGRID_.DLFC = DLFC;
	CEGRID_.XEL = XEL;
	CEGRID_.XE = XE;
	CEGRID_.XEK = XEK;
	CEGRID_.KE = KE;
	
}  
    

void transfcesi0_(double (*XESI)[NRP], int *IESIF, int *IESIL, int *NSESI, int *NCURE){
	
	CESI0_.XESI = XESI;
	CESI0_.IESIF = IESIF;
	CESI0_.IESIL = IESIL;
	CESI0_.NSESI = NSESI;
	CESI0_.NCURE = NCURE;

}     
     
void transfcpsi0_(double (*XPSI)[NRP], int *IPSIF, int *IPSIL, int *NSPSI, int *NCURP){
	
	CPSI0_.XPSI = XPSI;
	CPSI0_.IPSIF = IPSIF;
	CPSI0_.IPSIL = IPSIL;
	CPSI0_.NSPSI = NSPSI;
	CPSI0_.NCURP = NCURP;

} 


void transfcadata_(double *ATW, double *EPX, double *RSCR, double *ETA, double (*EB)[99], double (*ALW)[99], double (*CP0)[99], int (*IFI)[99], int (*IKS)[99], int *NSHT,  char (*LASYMB)[2]){
	
	CADATA_.IFI = IFI;
	CADATA_.IKS = IKS;
	CADATA_.NSHT = NSHT;
	CADATA_.ATW = ATW;
	CADATA_.EPX = EPX;
	CADATA_.RSCR = RSCR;
	CADATA_.ETA = ETA;
	CADATA_.EB = EB;
	CADATA_.ALW = ALW;
	CADATA_.CP0 = CP0;	
	CADATA_.LASYMB = LASYMB;
	
}
     
void transfcgph00_(int *IPHF, int *IPHL, int *NPHS, int *NCUR, double *EPH, double (*XPH)[NTP]){
	
	CGPH00_.IPHF = IPHF;
	CGPH00_.IPHL = IPHL;
	CGPH00_.NPHS = NPHS;
	CGPH00_.NCUR = NCUR;
	CGPH00_.EPH = EPH;
	CGPH00_.XPH = XPH;
		
}
     
     
void transfcrelax_(double *P, double *ET, double *F, int *IS0, int *IS1, int *IS2, int (*IFIRST)[99], int (*ILAST)[99], int *NCUR, int *KS, int *MODER){
	
	CRELAX_.IS0 = IS0;
	CRELAX_.IS1 = IS1;
	CRELAX_.IS2 = IS2;
	CRELAX_.IFIRST = IFIRST;
	CRELAX_.ILAST = ILAST;
	CRELAX_.NCUR = NCUR;
	CRELAX_.KS = KS; 
	CRELAX_.MODER = MODER;
	CRELAX_.P = P;
	CRELAX_.ET = ET;
	CRELAX_.F = F;

}
     
 
 void transfcritaa_(	double *XA, double *AA, double *BA, double *FA, int *IA, int *NPM1A){
	 
	CRITAA_.XA = XA;
	CRITAA_.AA = AA;
	CRITAA_.BA = BA; 
	CRITAA_.FA = FA;
	CRITAA_.IA = IA;
	CRITAA_.NPM1A = NPM1A;
 
 }    
     

 void transfcrndg3_(double *X, double *A, double *B, double *F, int *KA, int *NPM1){
	 
	 CRNDG3_.X = X;
	 CRNDG3_.A = A;
	 CRNDG3_.B = B;
	 CRNDG3_.F = F;
	 CRNDG3_.KA = KA;
	 CRNDG3_.NPM1 = NPM1;
	 
	 
 }
  
void transfcrita_(double *XT, double *PAC, double *DPAC, double *A, double *B, int *IL, int *IU, int *NPM1){
 //double *QTI, double *PACI, double *DPACI, double *AI, double *BI, int *ITLI, int *ITUI, int *NPM1I,
 //double *XTI){
 	
	CRITA_.XT = XT;
	CRITA_.PAC = PAC;
	CRITA_.DPAC = DPAC;
	CRITA_.A = A;
	CRITA_.B = B;
	CRITA_.IL = IL;
	CRITA_.IU = IU;
	CRITA_.NPM1 = NPM1;
	
	CRITA_.QTI = CRITA_.XT;
	CRITA_.PACI = CRITA_.PAC;
	CRITA_.DPACI = CRITA_.DPAC;
	CRITA_.AI = CRITA_.A;
	CRITA_.BI = CRITA_.B; 
	CRITA_.ITLI = CRITA_.IL;
	CRITA_.ITUI = CRITA_.IU;
	CRITA_.NPM1I = CRITA_.NPM1;
	
	CRITA_.XTI = CRITA_.XT;
	
}

void transfcritan_(double *CNORM){
	
	CRITAN_.CNORM = CNORM;
		
} 
  
  
void transfcompos_(double (*STF)[MAXMAT], double *ZT, double *AT, double *RHO, double *VMOL, int (*IZ)[MAXMAT], int *NELEM){
	
	COMPOS_.STF = STF;
	COMPOS_.ZT = ZT;
	COMPOS_.AT = AT;
	COMPOS_.RHO = RHO;
	COMPOS_.VMOL = VMOL;
	COMPOS_.IZ = IZ;
	COMPOS_.NELEM = NELEM;
	
}

void transfcrange_(double (*RANGE)[MAXMAT][3], double (*RANGEL)[MAXMAT][3]){
	
	CRANGE_.RANGE = RANGE;
	CRANGE_.RANGEL = RANGEL;
	
}

void transfcein_(double *EXPOT, double *OP2, double (*F)[MAXMAT],double (*UI)[MAXMAT], double (*WRI)[MAXMAT], int (*KZ)[MAXMAT], int (*KS)[MAXMAT], int *NOSC){

	CEIN_.EXPOT = EXPOT;
	CEIN_.OP2 = OP2;
	CEIN_.F = F;
	CEIN_.UI = UI;
	CEIN_.WRI = WRI;
	CEIN_.KZ = KZ;
	CEIN_.KS = KS;
	CEIN_.NOSC = NOSC;
	
		
}

void transfceintf_(double *T1EI, double *T2EI , double *T1PI, double *T2PI){
	
	CEINTF_.T1EI = T1EI;
	CEINTF_.T2EI = T2EI;
	CEINTF_.T1PI = T1PI;
	CEINTF_.T2PI = T2PI;
	
}

void transfcein00_(double *SEH0, double *SEH1, double *SEH2 , double *SES0, double *SES1, double *SES2, double *SET0, double *SET1, double *SET2){
	
	CEIN00_.SEH0 = SEH0;
	CEIN00_.SEH1 = SEH1;
	CEIN00_.SEH2 = SEH2;
	CEIN00_.SES0 = SES0;
	CEIN00_.SES1 = SES1;
	CEIN00_.SES2 = SES2;
	CEIN00_.SET0 = SET0;
	CEIN00_.SET1 = SET1;
	CEIN00_.SET2 = SET2;
	
}

void transfcpin00_(double *SPH0, double *SPH1, double *SPH2, double *SPS0, double *SPS1, double *SPS2, double *SPT0, double *SPT1, double *SPT2){
	
	CPIN00_.SPH0 = SPH0;
	CPIN00_.SPH1 = SPH1; 
	CPIN00_.SPH2 = SPH2;
	CPIN00_.SPS0 = SPS0;
	CPIN00_.SPS1 = SPS1;
	CPIN00_.SPS2 = SPS2;
	CPIN00_.SPT0 = SPT0;
	CPIN00_.SPT1 = SPT1;
	CPIN00_.SPT2 = SPT2;
	
}

void transfceinac_(double (*EINAC)[NEGP][MAXMAT], int (*IEIN)[MAXMAT], int *NEIN){
	
	CEINAC_.EINAC = EINAC;
	CEINAC_.IEIN = IEIN;
	CEINAC_.NEIN = NEIN;
	
}

void transfcesiac_(double (*ESIAC)[NEGP][MAXMAT], int (*IESI)[MAXMAT], int *NESI){

	CESIAC_.ESIAC = ESIAC;
	CESIAC_.IESI = IESI;
	CESIAC_.NESI =NESI;
	
	
}

void transfcesin_(double (*XSEIN)[NEGP], double (*XSESI)[NEGP], int *ISIE){
	
	CESIN_.XSEIN = XSEIN;
	CESIN_.XSESI = XSESI;
	CESIN_.ISIE = ISIE;
	
}

void transfcpinac_(double (*PINAC)[NEGP][MAXMAT], int (*IPIN)[MAXMAT], int *NPIN){
	
	CPINAC_.PINAC = PINAC;
	CPINAC_.IPIN = IPIN;
	CPINAC_.NPIN = NPIN;
	
}

void transfcpsiac_(double (*PSIAC)[NEGP][MAXMAT], int (*IPSI)[MAXMAT], int *NPSI){
	
	CPSIAC_.PSIAC = PSIAC;
	CPSIAC_.IPSI = IPSI; 
	CPSIAC_.NPSI = NPSI;
	
}

void transfcpsin_(double (*XSPIN)[NEGP], double (*XSPSI)[NEGP], int *ISIP){
	CPSIN_.XSPIN = XSPIN;
	CPSIN_.XSPSI = XSPSI;
	CPSIN_.ISIP = ISIP;
	
}  
  
  
void transfcgco_(double (*FCO)[MAXMAT], double (*UICO)[MAXMAT], double (*FJ0)[MAXMAT], double (*PTRSH)[MAXMAT], int (*KZCO)[MAXMAT], int (*KSCO)[MAXMAT], int *NOSCCO){

	CGCO_.FCO   =  FCO;
	CGCO_.UICO  =  UICO;
	CGCO_.FJ0   =  FJ0;
	CGCO_.PTRSH =  PTRSH;
	CGCO_.KZCO  =  KZCO;
	CGCO_.KSCO  =  KSCO;
	CGCO_.NOSCCO = NOSCCO;
	
	
}
	

void transfceimfp_(double (*SEHEL)[MAXMAT], double(*SEHIN)[MAXMAT], double(*SEISI)[MAXMAT], double(*SEHBR)[MAXMAT],double (*SEAUX)[MAXMAT], 
          double (*SETOT)[MAXMAT], double (*CSTPE)[MAXMAT], double (*RSTPE)[MAXMAT], double (*DEL)[MAXMAT], double (*W1E)[MAXMAT],
          double (*W2E)[MAXMAT], double (*DW1EL)[MAXMAT],  double(*DW2EL)[MAXMAT], double(*RNDCE)[MAXMAT],double (*AE)[MAXMAT], double (*BE)[MAXMAT],
    	   double (*T1E)[MAXMAT],double (*T2E)[MAXMAT]){
			   
	CEIMFP_.SEHEL = SEHEL;
	CEIMFP_.SEHIN = SEHIN;
	CEIMFP_.SEISI = SEISI;
	CEIMFP_.SEHBR = SEHBR;
	CEIMFP_.SEAUX = SEAUX;
    CEIMFP_.SETOT = SETOT;
	CEIMFP_.CSTPE = CSTPE;
	CEIMFP_.RSTPE = RSTPE;
	CEIMFP_.DEL =   DEL;
	CEIMFP_.W1E =   W1E;
    CEIMFP_.W2E =   W2E;
	CEIMFP_.DW1EL = DW1EL;
	CEIMFP_.DW2EL = DW2EL;
	CEIMFP_.RNDCE = RNDCE;
	CEIMFP_.AE =    AE;
	CEIMFP_.BE =    BE;
    CEIMFP_.T1E =   T1E;
	CEIMFP_.T2E =   T2E;	   	   
}

void transfclas1e_(double (*TSTPE)[MAXMAT], double(*TSTRE)[MAXMAT], double(*TRL1E)[MAXMAT], double(*TRL2E)[MAXMAT]){
	
	CLAS1E_.TSTPE = TSTPE;
	CLAS1E_.TSTRE = TSTRE;
	CLAS1E_.TRL1E = TRL1E;
	CLAS1E_.TRL2E = TRL2E;
	
}

void transfcpimfp_(double (*SPHEL)[MAXMAT], double(*SPHIN)[MAXMAT], double(*SPISI)[MAXMAT], double(*SPHBR)[MAXMAT],double (*SPAN)[MAXMAT], 
           double (*SPAUX)[MAXMAT], double (*SPTOT)[MAXMAT], double (*CSTPP)[MAXMAT], double (*RSTPP)[MAXMAT], double (*W1P)[MAXMAT],
           double (*W2P)[MAXMAT], double (*DW1PL)[MAXMAT], double (*DW2PL)[MAXMAT], double (*RNDCP)[MAXMAT], double (*AP)[MAXMAT], double (*BP)[MAXMAT],
           double (*T1P)[MAXMAT], double (*T2P)[MAXMAT]){
			
	CPIMFP_.SPHEL = SPHEL;
	CPIMFP_.SPHIN = SPHIN;
	CPIMFP_.SPISI = SPISI;
	CPIMFP_.SPHBR = SPHBR;
	CPIMFP_.SPAN =  SPAN;
	CPIMFP_.SPAUX = SPAUX;
	CPIMFP_.SPTOT = SPTOT;
	CPIMFP_.CSTPP = CSTPP;
	CPIMFP_.RSTPP = RSTPP;
	CPIMFP_.W1P =   W1P;
	CPIMFP_.W2P =   W2P;
	CPIMFP_.DW1PL = DW1PL;
	CPIMFP_.DW2PL = DW2PL;
	CPIMFP_.RNDCP = RNDCP;
	CPIMFP_.AP =    AP;
	CPIMFP_.BP =    BP;
	CPIMFP_.T1P =   T1P;
	CPIMFP_.T2P =   T2P;			   
			   
}

void transfclas1p_(double (*TSTPP)[MAXMAT], double (*TSTRP)[MAXMAT], double (*TRL1P)[MAXMAT], double (*TRL2P)[MAXMAT]){

	CLAS1P_.TSTPP = TSTPP; 
	CLAS1P_.TSTRP = TSTRP; 
	CLAS1P_.TRL1P = TRL1P; 
	CLAS1P_.TRL2P = TRL2P; 
}

void transfceel00_(double *EJT, double  *XE0,double  *XE1, double *XE2, double *XP0, double *XP1, double *XP2, double *T1E0, double *T2E0, double *T1P0, 
            double *T2P0, double *EJTL, double *FJL, double *A, double *B, double *C, double *D){

	CEEL00_.EJT  = EJT; 
	CEEL00_.XE0  = XE0; 
	CEEL00_.XE1  = XE1; 
	CEEL00_.XE2  = XE2; 
	CEEL00_.XP0  = XP0; 
	CEEL00_.XP1  = XP1; 
	CEEL00_.XP2  = XP2; 
	CEEL00_.T1E0 = T1E0;
	CEEL00_.T2E0 = T2E0;
	CEEL00_.T1P0 = T1P0;
	CEEL00_.T2P0 = T2P0;
	CEEL00_.EJTL = EJTL;
	CEEL00_.FJL  = FJL; 
	CEEL00_.A    = A;   
	CEEL00_.B    = B;   
	CEEL00_.C    = C;   
	CEEL00_.D =    D;
					
}

void transfcbryld_(double (*EBRY)[MAXMAT], double (*PBRY)[MAXMAT]){

	CBRYLD_.EBRY = EBRY;

	CBRYLD_.PBRY = PBRY;
}

void transfcgimfp_(double (*SGRA)[MAXMAT], double (*SGCO)[MAXMAT], double (*SGPH)[MAXMAT], double (*SGPP)[MAXMAT], double (*SGAUX)[MAXMAT]){
	
	CGIMFP_.SGRA =  SGRA;
	CGIMFP_.SGCO =  SGCO;
	CGIMFP_.SGPH =  SGPH;
	CGIMFP_.SGPP =  SGPP;
	CGIMFP_.SGAUX = SGAUX;
	
	
}

void transfcgph01_(double *ER, double *XSR, int *NPHD){
	CGPH01_.ER = ER;
	CGPH01_.XSR = XSR;
	CGPH01_.NPHD = NPHD;
}

void transfcgpp01_(double (*TRIP)[MAXMAT]){
	
	CGPP01_.TRIP = TRIP;
	
}
  
  
  
void transfcebr_(double *WB, double (*PBCUT)[MAXMAT], double (*WBCUT)[MAXMAT], double (*PDFB)[NEGP][MAXMAT], double (*DPDFB)[NEGP][MAXMAT], double (*PACB)[NEGP][MAXMAT], double *ZBR2){

	CEBR_.WB = WB;
	CEBR_.PBCUT = PBCUT;
	CEBR_.WBCUT = WBCUT;
	CEBR_.PDFB = PDFB;
	CEBR_.DPDFB = DPDFB;
	CEBR_.PACB = PACB;
	CEBR_.ZBR2 = ZBR2;
}
	
void transfcebr01_(double *EBT, double (*XS)[NBE], double *TXS, double *X, double *Y){
	
	CEBR01_.EBT = EBT;
	CEBR01_.XS = XS;
	CEBR01_.TXS = TXS;
	CEBR01_.X = X;
	CEBR01_.Y = Y;
}



void transfcebr02_(double (*P0)[NEGP][MAXMAT]){
	
	CEBR02_.P0 = P0;
	
}
  
  
 void transfcbrang_(double *BET, double *BK, double(*BP1)[21][6][MAXMAT], double (*BP2)[21][6][MAXMAT], double *ZBEQ){
	CBRANG_.BET = BET;
	CBRANG_.BK = BK;
	CBRANG_.BP1 = BP1;
	CBRANG_.BP2 = BP2;
	CBRANG_.ZBEQ = ZBEQ;
 } 
  
  
 void transfcein01_(double *EI, double *EE, double *CPS, double *AMOL, int *MOM){
	
	CEIN01_.EI = EI;
	CEIN01_.EE = EE;
	CEIN01_.CPS = CPS; 
	CEIN01_.AMOL = AMOL; 
	CEIN01_.MOM = MOM;
	
}

void transfcsumga_(int *IERGA, int *NCALL){
	
	CSUMGA_.IERGA = IERGA;
	CSUMGA_.NCALL = NCALL;
	
}

void transfcpin01_(double *EI, double *CPS, double *BHA1, double *BHA2, double *BHA3, double *BHA4, int *MOM){
	CPIN01_.EI = EI;
	CPIN01_.CPS = CPS;
	CPIN01_.BHA1 = BHA1;
	CPIN01_.BHA2 = BHA2;
	CPIN01_.BHA3 = BHA3;
	CPIN01_.BHA4 = BHA4;
	CPIN01_.MOM = MOM;
}
  
  
void transfcdcsep_(double *ETS, double *ETL, double *TH, double *THR, double *XMU , double *XMUL, double *ECS, double *ETCS1 , double *ETCS2, double (*EDCS)[NE],
     	   double *PCS, double *PTCS1,  double *PTCS2, double (*PDCS)[NE], double *DCSI, double *DCSIL, double *CSI, double *TCS1I, double *TCS2I){
				
    CDCSEP_.ETS = ETS;
	CDCSEP_.ETL = ETL;
	CDCSEP_.TH = TH;
	CDCSEP_.THR = THR;
	CDCSEP_.XMU = XMU; 
	CDCSEP_.XMUL = XMUL;
	CDCSEP_.ECS = ECS; 
	CDCSEP_.ETCS1 = ETCS1;
	CDCSEP_.ETCS2 = ETCS2;
	CDCSEP_.EDCS = EDCS;
	CDCSEP_.PCS = PCS;
	CDCSEP_.PTCS1 = PTCS1;
	CDCSEP_.PTCS2 = PTCS2;
	CDCSEP_.PDCS = PDCS;
	CDCSEP_.DCSI =  DCSI;
	CDCSEP_.DCSIL = DCSIL;
	CDCSEP_.CSI = CSI;
	CDCSEP_.TCS1I = TCS1I;
	CDCSEP_.TCS2I =  TCS2I;
												
}

void transfceeldb_(	double (*XSE)[NEGP][NP], double (*PSE)[NEGP][NP], double(*ASE)[NEGP][NP], double (*BSE)[NEGP][NP],
    int (*ITLE)[NEGP][NP], int (*ITUE)[NEGP][NP]){
		
    CEELDB_.XSE = XSE;
	CEELDB_.PSE = PSE;
	CEELDB_.ASE = ASE;
	CEELDB_.BSE = BSE;
	CEELDB_.ITLE = ITLE;
	CEELDB_.ITUE = ITUE;
				
}

void transfcpeldb_(double (*XSP)[NEGP][NP],double  (*PSP)[NEGP][NP],double (*ASP)[NEGP][NP],double  (*BSP)[NEGP][NP],
    int(*ITLP)[NEGP][NP], int (*ITUP)[NEGP][NP]){
		
    CPELDB_.XSP = XSP;
	CPELDB_.PSP = PSP;
	CPELDB_.ASP = ASP;
	CPELDB_.BSP = BSP;
	CPELDB_.ITLP = ITLP;
	CPELDB_.ITUP = ITUP;
		
				
}

void transfcelsep_(double *EELMAX, double *PELMAX, double (*RNDCED)[MAXMAT], double (*RNDCPD)[MAXMAT]){
	
    CELSEP_.EELMAX = EELMAX;
	CELSEP_.PELMAX = PELMAX;
	CELSEP_.RNDCED = RNDCED;
	CELSEP_.RNDCPD = RNDCPD;
}
  
  
  
void transfcgra00_(double *FACTE, double *Q2MAX, int *MM, int *MOM){
	CGRA00_.FACTE = FACTE;
	CGRA00_.Q2MAX = Q2MAX;
	CGRA00_.MM = MM;
	CGRA00_.MOM = MOM;	
}

void transfcgra01_(double (*FF)[MAXMAT], double *ERA, double (*XSRA)[MAXMAT], int *IED, int *IEU, int *NE){
    CGRA01_.FF = FF;
	CGRA01_.ERA = ERA;
	CGRA01_.XSRA = XSRA;
	CGRA01_.IED = IED; 
	CGRA01_.IEU = IEU;
	CGRA01_.NE = NE;	
}

void transfcgra02_(double *QQ, double (*AR)[MAXMAT], double (*BR)[MAXMAT], double (*CR)[MAXMAT], double (*DR)[MAXMAT], double *FF0, double *QQM){
    CGRA02_.QQ = QQ;
	CGRA02_.AR = AR;
	CGRA02_.BR = BR;
	CGRA02_.CR = CR;
	CGRA02_.DR = DR;
	CGRA02_.FF0 = FF0;
	CGRA02_.QQM = QQM;	
}

void transfcgra03_(double  (*QRA)[NP2], double (*PRA)[NP2], double (*DPRA)[NP2], double (*ARA)[NP2],
			double (*BRA)[NP2],double (*PMAX)[NEGP], int		(*ITLRA)[NP2],int (*ITURA)[NP2]){
    CGRA03_.QRA = QRA;
	CGRA03_.PRA = PRA;
	CGRA03_.DPRA = DPRA;
	CGRA03_.ARA = ARA;
	CGRA03_.BRA = BRA;
	CGRA03_.PMAX = PMAX;
	CGRA03_.ITLRA = ITLRA;
	CGRA03_.ITURA = ITURA;	
}

void transfcgpp00_(double *ZEQPP, double (*F0)[MAXMAT], double *BCB){
	CGPP00_.ZEQPP = ZEQPP;
	CGPP00_.F0 = F0;
	CGPP00_.BCB = BCB;
	
}
  
  
  //implementar todas as transfs abaixo
  void transfrseed_(int *ISEED1, int *ISEED2){
	  RSEED_.ISEED1 = ISEED1;
	  RSEED_.ISEED2 = ISEED2;
  }

void transfctitle_(char *TITLE, char *TITLE2){
	CTITLE_.TITLE = TITLE;
	CTITLE_.TITLE2 = TITLE2;
}

void transfcdate_(char *DATE23){
	CDATE_.DATE23 = DATE23;
}

void transfcspgeo_(double *DSMAX, double (*EABSB)[3]){
	CSPGEO_.DSMAX = DSMAX; 
	CSPGEO_.EABSB = EABSB;
}

void transfcforci_(double (*WLOW)[NBV], double(*WHIG)[NBV], bool (*LFORCE)[NBV]){
	CFORCI_.WLOW = WLOW;
	CFORCI_.WHIG = WHIG;
	CFORCI_.LFORCE = LFORCE;
}

void transfcxrspl_(int *IXRSPL, int *ILBA, bool *LXRSPL){
	
	CXRSPL_.IXRSPL = IXRSPL;
	CXRSPL_.ILBA = ILBA;
	CXRSPL_.LXRSPL = LXRSPL;

}

void transfcsour0_(double *CTHL, double *DCTH, double *PHIL, double *DPHI, int *KPARP, int *JOBEND, bool *LSCONE, bool *LGPOL, bool *LPSF){
   CSOUR0_.CTHL =   CTHL;
   CSOUR0_.DCTH =   DCTH;
   CSOUR0_.PHIL =   PHIL;
   CSOUR0_.DPHI =   DPHI;
   CSOUR0_.KPARP =  KPARP;
   CSOUR0_.JOBEND = JOBEND;
   CSOUR0_.LSCONE = LSCONE;
   CSOUR0_.LGPOL =  LGPOL;
   CSOUR0_.LPSF =   LPSF;
 	
}

void transfcsour1_(double *E0, double *EPMAX, double *SP10, double *SP20, double *SP30){
    CSOUR1_.E0    = E0;    
	CSOUR1_.EPMAX = EPMAX;
	CSOUR1_.SP10  = SP10; 
	CSOUR1_.SP20  = SP20;
	CSOUR1_.SP30  = SP30; 

}

void transfcsour2_(double *ESRC, double *PSRC, int *IASRC, double *FSRC, bool *LSPEC){
	CSOUR2_.ESRC = ESRC; 
	CSOUR2_.PSRC = PSRC;
	CSOUR2_.IASRC = IASRC;
	CSOUR2_.FSRC = FSRC;
	CSOUR2_.LSPEC = LSPEC;
}

void transfcsour3_(double *SX0, double *SY0, double *SZ0, double *SSX, double *SSY, double *SSZ, int *IXSBOD, bool *LEXSRC, bool *LEXBD){
    CSOUR3_.SX0 = SX0; 
	CSOUR3_.SY0 = SY0; 
	CSOUR3_.SZ0 = SZ0; 
	CSOUR3_.SSX = SSX; 
	CSOUR3_.SSY = SSY; 
	CSOUR3_.SSZ = SSZ; 
	CSOUR3_.IXSBOD = IXSBOD;
	CSOUR3_.LEXSRC = LEXSRC;
	CSOUR3_.LEXBD = LEXBD; 
}

void transfcsour4_(double *WGMIN, double *RWGMIN, double *WGMAX, double *RLREAD, int *IPSFI, int *NPSF, int *NPSN, int *NSPLIT, int *KODEPS){
	CSOUR4_.WGMIN =   WGMIN;
	CSOUR4_.RWGMIN =	RWGMIN;
	CSOUR4_.WGMAX =	 WGMAX;
	CSOUR4_.RLREAD =	RLREAD;
	CSOUR4_.IPSFI =	IPSFI;
	CSOUR4_.NPSF =	NPSF;
	CSOUR4_.NPSN =	NPSN;
	CSOUR4_.NSPLIT =	NSPLIT;
	CSOUR4_.KODEPS =	KODEPS;

}

void transfcsour5_(char (*PSFI)[20]){
	CSOUR5_.PSFI = PSFI;

}

void transfcnt0_(double *PRIM , double *PRIM2, double *DPRIM, double (*SEC)[3], double (*SEC2)[3], double (*DSEC)[3], 
			    double *AVW, double *AVW2, double *DAVW, double *AVA, double *AVA2, double *DAVA, double *AVE, double *AVE2, double *DAVE){
	
   CNT0_.PRIM =    PRIM;
	CNT0_.PRIM2	= PRIM2;
	CNT0_.DPRIM	= DPRIM;
	CNT0_.SEC	= SEC;
	CNT0_.SEC2	= SEC2;
	CNT0_.DSEC	= DSEC;
	CNT0_.AVW	= AVW;
	CNT0_.AVW2	= AVW2;
	CNT0_.DAVW	= DAVW;
	CNT0_.AVA	= AVA;
	CNT0_.AVA2	= AVA2;
	CNT0_.DAVA	= DAVA;
	CNT0_.AVE	= AVE;
	CNT0_.AVE2	= AVE2;
	CNT0_.DAVE	= DAVE;

}

/*void transfcnt0_(double *PRIM , double *PRIM2, double *DPRIM, double (*SEC)[3], double (*SEC2)[3], double (*DSEC)[3], 
			    double *AVW, double *AVW2, double *DAVW, double *AVA, double *AVA2, double *DAVA, double *AVE, double *AVE2, double *DAVE){

    CNT0_.PRIM[0] =    PRIM;
	CNT0_.PRIM2[0]	= PRIM2;
	CNT0_.DPRIM[0]	= DPRIM;
	CNT0_.SEC[0][0]	= SEC;
	CNT0_.SEC2[0][0]	= SEC2;
	CNT0_.DSEC[0][0]	= DSEC;
	CNT0_.AVW[0]	= AVW;
	CNT0_.AVW2[0]	= AVW2;
	CNT0_.DAVW[0]	= DAVW;
	CNT0_.AVA[0]	= AVA;
	CNT0_.AVA2[0]	= AVA2;
	CNT0_.DAVA[0]	= DAVA;
	CNT0_.AVE[0]	= AVE;
	CNT0_.AVE2[0]	= AVE2;
	CNT0_.DAVE[0]	= DAVE;

}*/



void transfcnt1_(double *TDEBO, double *TDEBO2, double *DEBO){
	CNT1_.TDEBO = TDEBO;
	CNT1_.TDEBO2 = TDEBO2;
	CNT1_.DEBO = DEBO;
}

void transfcnt2_(double *SHIST, int *NSEB){
	CNT2_.SHIST = SHIST; 
	CNT2_.NSEB = NSEB;
}

void transfcnt3_(double (*SEDS)[3], double (*SEDS2)[3], double *DSDE, double *RDSDE, int *NSDE){
	CNT3_.SEDS = SEDS;
	CNT3_.SEDS2 = SEDS2;
	CNT3_.DSDE = DSDE;
	CNT3_.RDSDE = RDSDE;
	CNT3_.NSDE =  NSDE;
}

void transfcnt4_(double *RLAST, double *RWRITE, int *IDCUT, int (*KKDI)[NIDM], int *IPSF, int *NID, int *NPSFO, int *IPSFO){
	CNT4_.RLAST	= RLAST;
	CNT4_.RWRITE =	RWRITE;
	CNT4_.IDCUT	= IDCUT;
	CNT4_.KKDI	= KKDI;
	CNT4_.IPSF =	IPSF;
	CNT4_.NID =	NID;
	CNT4_.NPSFO	 = NPSFO;
	CNT4_.IPSFO	= IPSFO;
}

void transfcnt5_(double *DEDE, int *KBDE, int *NED){
	CNT5_.DEDE = DEDE;
	CNT5_.KBDE = KBDE;
	CNT5_.NED = NED;
}

void transfcnt6_(bool *LDOSEM){
	CNT6_.LDOSEM = LDOSEM;
}

void transfcdump_(bool *LDUMP, char *PFILED){
	CDUMP_.LDUMP = LDUMP;
	CDUMP_.PFILED = PFILED;
}

void transfcntrl_(double *TSIM, double *TSEC, double *TSECA, double *TSECAD, double *CPUT0, double *DUMPP, double *DSHN, double *SHN, int *N){
	CNTRL_.TSIM	= TSIM;
	CNTRL_.TSEC =	TSEC;
	CNTRL_.TSECA =	TSECA;
	CNTRL_.TSECAD =	TSECAD;
	CNTRL_.CPUT0 =	CPUT0;
	CNTRL_.DUMPP =	DUMPP;
	CNTRL_.DSHN =	DSHN;
	CNTRL_.SHN =	SHN;
	CNTRL_.N =	N;
}

//Parâmetros para direções de amostragem dentro de um cone.
void transfcgcone_(double *CPCT, double *CPST, double *SPCT, double *SPST, double *SPHI, double *CPHI, double *STHE, double *CTHE, double *CAPER){
    CGCONE_.CPCT = CPCT;
	CGCONE_.CPST = CPST;
	CGCONE_.SPCT = SPCT;
	CGCONE_.SPST = SPST;
	CGCONE_.SPHI = SPHI;
	CGCONE_.CPHI = CPHI;
	CGCONE_.STHE = STHE;
	CGCONE_.CTHE = CTHE;
	CGCONE_.CAPER = CAPER;
}


void transfcenang_(	double *EL, double *EU, double *THL, double *THU, double *BSE, double *RBSE, double *BSTH, double *RBSTH, double *BSPH, double *RBSPH,
	double (*PDE)[2][3], double (*PDE2)[2][3], double (*PDEP)[2][3],
	double (*PDA)[NBTHM][3], double (*PDA2)[NBTHM][3], double (*PDAP)[NBTHM][3],
	bool (*LPDE)[2][3], bool (*LPDA)[NBTHM][3],
	int *NE, int *NTH, int *NPH,
	bool *LLE, bool *LLTH){

			CENANG_.EL =	EL;
			CENANG_.EU =	EU;
			CENANG_.THL =	THL;
			CENANG_.THU	= THU;
			CENANG_.BSE	= BSE;
			CENANG_.RBSE	= RBSE;
			CENANG_.BSTH	= BSTH;
			CENANG_.RBSTH	= RBSTH;
			CENANG_.BSPH	= BSPH;
			CENANG_.RBSPH	= RBSPH;
            CENANG_.PDE    = PDE;
			CENANG_.PDE2	= PDE2;
			CENANG_.PDEP	= PDEP;
            CENANG_.PDA    = PDA;
			CENANG_.PDA2	= PDA2;
			CENANG_.PDAP	= PDAP;
            CENANG_.LPDE    = LPDE;
			CENANG_.LPDA	= LPDA;
			CENANG_.NE	= NE;
			CENANG_.NTH	= NTH;
			CENANG_.NPH	= NPH;
			CENANG_.LLE	= LLE;
			CENANG_.LLTH	= LLTH;

	}


void transfcimdet_(double *EL, double *EU, double *BSE, double *RBSE,
	double (*ET)[NIDM], double *EDEP, double *EDEP2, double *EDEPP,
    double  (*DIT)[NIDM], double (*DIT2)[NIDM], double (*DITP)[NIDM],
    double  (*DIP)[NBEM2][NIDM], double (*DIP2)[NBEM2][NIDM], double (*DIPP)[NBEM2][NIDM],
    double  (*FLT)[NIDM], double (*FLT2)[NIDM], double (*FLTP)[NIDM],
    double  (*FLP)[NBEM2][NIDM], double (*FLP2)[NBEM2][NIDM], double (*FLPP)[NBEM2][NIDM],
    double  *AGEL, double *AGEU, double *BAGE, double *RBAGE, double (*AGE)[NIDM],
    double  (*AGE2)[NIDM], double (*AGEP)[NIDM],
	bool *LEDEP, bool (*LDIT)[NIDM], bool (*LDIP)[NBEM2][NIDM], bool (*LFLT)[NIDM], bool (*LFLP)[NBEM2][NIDM], bool (*LAGEA)[NIDM],
	int *IDCUT, int *NE,
	bool *LLE, bool *LLAGE,
    int *NAGE, int *NID){

	CIMDET_.EL =	EL; 	
	CIMDET_.EU =	EU;
	CIMDET_.BSE =	BSE;
	CIMDET_.RBSE =	RBSE; 
	CIMDET_.ET =	ET; 
	CIMDET_.EDEP =	EDEP; 
	CIMDET_.EDEP2 =	EDEP2; 
	CIMDET_.EDEPP =	EDEPP; 
	CIMDET_.DIT =	DIT; 
	CIMDET_.DIT2 =	DIT2; 
	CIMDET_.DITP =	DITP; 
	CIMDET_.DIP =	DIP; 
	CIMDET_.DIP2 =	DIP2; 
	CIMDET_.DIPP =	DIPP; 
	CIMDET_.FLT =	FLT; 
	CIMDET_.FLT2 =	FLT2; 
	CIMDET_.FLTP =	FLTP;
	CIMDET_.FLP =	FLP; 
	CIMDET_.FLP2 =	FLP2; 
	CIMDET_.FLPP =	FLPP; 
	CIMDET_.AGEL =	AGEL; 
	CIMDET_.AGEU =	AGEU; 
	CIMDET_.BAGE =	BAGE; 
	CIMDET_.RBAGE =	RBAGE; 
	CIMDET_.AGE =	AGE; 
	CIMDET_.AGE2 =	AGE2; 
	CIMDET_.AGEP =	AGEP; 
	CIMDET_.LEDEP =	LEDEP; 
	CIMDET_.LDIT =	LDIT; 
	CIMDET_.LDIP =	LDIP; 
	CIMDET_.LFLT =	LFLT; 
	CIMDET_.LFLP =	LFLP; 
	CIMDET_.LAGEA =	LAGEA; 
	CIMDET_.IDCUT =	IDCUT; 
	CIMDET_.NE =	NE; 
	CIMDET_.LLE =	LLE; 
	CIMDET_.LLAGE =	LLAGE;
	CIMDET_.NAGE =	NAGE; 
	CIMDET_.NID =	NID;

	}

void transfcendet_(double *EL, double *EU, double *BSE, double *RBSE, double *EDEP, double *EDEP2, double (*DET)[NIDM],
    int  *NE, int *NID, bool *LLE){
	CENDET_.EL = EL;
	CENDET_.EU = EU;
	CENDET_.BSE = BSE;
	CENDET_.RBSE = RBSE;
	CENDET_.EDEP = EDEP;
	CENDET_.EDEP2 = EDEP2;
	CENDET_.DET = DET;
	CENDET_.NE = NE;
	CENDET_.NID = NID;
	CENDET_.LLE = LLE;
}  


void transfcdose1_(double (*DOSE)[NDYM][NDXM], double (*DOSE2)[NDYM][NDXM], double (*DOSEP)[NDYM][NDXM], int (*LDOSE)[NDYM][NDXM], int *KDOSE){
	CDOSE1_.DOSE = DOSE;
	CDOSE1_.DOSE2 = DOSE2;
	CDOSE1_.DOSEP = DOSEP;
	CDOSE1_.LDOSE = LDOSE;
	CDOSE1_.KDOSE = KDOSE;

}

void transfcdose2_(double *DDOSE, double *DDOSE2, double *DDOSEP, int *LDDOSE){
	CDOSE2_.DDOSE = DDOSE;
	CDOSE2_.DDOSE2 = DDOSE2;
	CDOSE2_.DDOSEP = DDOSEP;
	CDOSE2_.LDDOSE = LDDOSE;

}

void transfcdose3_(double *DXL, double *DXU, double *BDOSE, double *RBDOSE, int *NDB){
	CDOSE3_.DXL = DXL;
	CDOSE3_.DXU = DXU;
	CDOSE3_.BDOSE = BDOSE;
	CDOSE3_.RBDOSE = RBDOSE;
	CDOSE3_.NDB = NDB;

}

void transfcdose4_(double (*VMASS)[NDYM][NDXM]){
	CDOSE4_.VMASS = VMASS;
}


void transfcjump1_(double *ELAST1, double *ELAST2, int *MHINGE, int *KSOFTE, int *KSOFTI, int *KDELTA){

	CJUMP1_.ELAST1 = ELAST1;
	CJUMP1_.ELAST2 = ELAST2;
	CJUMP1_.MHINGE = MHINGE;
	CJUMP1_.KSOFTE = KSOFTE;
	CJUMP1_.KSOFTI = KSOFTI;
	CJUMP1_.KDELTA = KDELTA;
}

void transfsecst_(double *ES, double *XS, double *YS, double *ZS, double *US , double *VS, double *WS, double *WGHTS, double *SP1S, double *SP2S, double *SP3S, double *PAGES,
    int *KS, int *IBODYS, int *MS, int (*ILBS)[5], int *IPOLS, int *NSEC){
		SECST_.ES = ES;
		SECST_.XS = XS;
		SECST_.YS = YS;
		SECST_.ZS = ZS;
		SECST_.US = US;
		SECST_.VS = VS;
		SECST_.WS = WS;
		SECST_.WGHTS = WGHTS;
		SECST_.SP1S = SP1S;
		SECST_.SP2S = SP2S;
		SECST_.SP3S = SP3S;
		SECST_.PAGES = PAGES;
        SECST_.KS = KS;
		SECST_.IBODYS = IBODYS;
		SECST_.MS = MS;
		SECST_.ILBS = ILBS;
		SECST_.IPOLS = IPOLS;
		SECST_.NSEC = NSEC;
	}


void transfcjump0_(double *P, double *ST, double *DST, double *DSR, double *W1, double *W2, double *T1, double *T2){
	CJUMP0_.P = P;
	CJUMP0_.ST = ST;
	CJUMP0_.DST = DST;
	CJUMP0_.DSR = DSR;
	CJUMP0_.W1 = W1;
	CJUMP0_.W2 = W2;
	CJUMP0_.T1 = T1;
	CJUMP0_.T2 = T2;
}

void transfchist_(int *ILBA){
	CHIST_.ILBA = ILBA;
}


      
void fsurf2_(int& KS, double& A, double& B, double& C){
	
	
	/*
	Calcula os par�metros da fun��o mestre da superf�cie KS e o raio (X, Y, Z) + S * (U, V, W).
	*/
 
     
    if (QSURF_.KPLANE[KS-1] == 0) {
        A = *TRACK_mod_.U * (QSURF_.AXX[KS-1] * *TRACK_mod_.U + QSURF_.AXY[KS-1] * *TRACK_mod_.V + QSURF_.AXZ[KS-1] * *TRACK_mod_.W) +
            *TRACK_mod_.V * (QSURF_.AYY[KS-1] * *TRACK_mod_.V + QSURF_.AYZ[KS-1] * *TRACK_mod_.W) + *TRACK_mod_.W * QSURF_.AZZ[KS-1] * *TRACK_mod_.W;
		double XXX = QSURF_.AXX[KS-1] * *TRACK_mod_.X + QSURF_.AXY[KS-1] * *TRACK_mod_.Y + QSURF_.AXZ[KS-1] * *TRACK_mod_.Z + QSURF_.AX[KS-1];
        double YYY = QSURF_.AYY[KS-1] * *TRACK_mod_.Y + QSURF_.AYZ[KS-1] * *TRACK_mod_.Z + QSURF_.AY[KS-1];
        double ZZZ = QSURF_.AZZ[KS-1] * *TRACK_mod_.Z + QSURF_.AZ[KS-1];

        B = *TRACK_mod_.U * (QSURF_.AXX[KS-1] * *TRACK_mod_.X + XXX) + *TRACK_mod_.V * (QSURF_.AXY[KS-1] * *TRACK_mod_.X + QSURF_.AYY[KS-1] * *TRACK_mod_.Y + YYY) +
            *TRACK_mod_.W * (QSURF_.AXZ[KS-1] * *TRACK_mod_.X + QSURF_.AYZ[KS-1] * *TRACK_mod_.Y + QSURF_.AZZ[KS-1] * *TRACK_mod_.Z + ZZZ);

        C = *TRACK_mod_.X * XXX + *TRACK_mod_.Y * YYY + *TRACK_mod_.Z * ZZZ + QSURF_.A0[KS-1];
    }
    else{
        A = 0.0e0;
    	B = *TRACK_mod_.U * QSURF_.AX[KS-1] + *TRACK_mod_.V * QSURF_.AY[KS-1] + *TRACK_mod_.W * QSURF_.AZ[KS-1];
    	C = *TRACK_mod_.X * QSURF_.AX[KS-1] + *TRACK_mod_.Y * QSURF_.AY[KS-1] + *TRACK_mod_.Z * QSURF_.AZ[KS-1] + QSURF_.A0[KS-1];
	} 
     
     
}

void locate2_(){

/*	if (imprimiu == 0){
		printf("\n\nLOCATE2\n\n");
		imprimiu++;
	}*/
	
	
	/*
	
	Esta sub-rotina determina o corpo que cont�m o ponto com
   coordenadas (X, Y, Z). Os efeitos dos erros de arredondamento num�ricos s�o
   evitado considerando superf�cies difusas, que aumentam ou encolhem ligeiramente
   quando a part�cula os cruza.
 
   Valores de entrada (m�dulo TRACK_mod):
      X, Y, Z ... coordenadas da part�cula,
      U, V, W ... dire��o do movimento.
 
   Valores de sa�da (m�dulo TRACK_mod):
      IBODY ..... corpo onde a part�cula se move,
      MAT ...... material em IBODY,
                     MAT = 0, indica uma regi�o vazia
	
	
	*/
	

    double A = 0.0;
    double B = 0.0;
    double C = 0.0;
    double FUZZ = 0.0;
    const double FUZZL = 1.0e-12;;
    int KS;
    int KB;
    int KF;
    double ABSA = 0.0;
											
    for (int I = 1; I <= *QSURF_.NSURF; I++) { 
	    QTREE_.KSP[I-1] = 0;
    }
	int KB0 = *QTREE_.NBODYS;
	
d100:
	for (int KSS = 1; KSS <= QTREE_.KSURF[NXG-1][KB0-1]; KSS++) {
		KS = QTREE_.KSURF[KSS-1][KB0-1];
        if ((QTREE_.KSP[KS-1] != 0) || (QTREE_.KFLAG[KSS-1][KB0-1] > 4)) {
        	goto d101;
		}
		
		fsurf2_(KS, A, B, C);
		
	   	ABSA = fabs(A);

    	if (ABSA > 1.0e-36)
	        FUZZ = FUZZL * (B * B - 4.0e0 * A * C) / ABSA;
        else
   	   	   FUZZ = FUZZL * fabs(B);

    	if (C < (-FUZZ)) {
			QTREE_.KSP[KS-1] = 1;
		} else if (C > FUZZ){
				QTREE_.KSP[KS-1] = 2;
		} else {
			if (B < 0.0e0)	
				QTREE_.KSP[KS-1] = 1; //particula movendo-se para dentro
			else
				QTREE_.KSP[KS-1] = 2; // particula movendo-se para fora
		}
d101: ;
	}
	for (int KBB = 1; KBB <= QTREE_.KDGHT[NXG-1][KB0 - 1]; KBB++) {
		KB = QTREE_.KDGHT[KBB - 1][KB0 - 1];
		for (int KSS = 1; KSS <= QTREE_.KSURF[NXG-1][KB-1]; KSS++) {
			KS = QTREE_.KSURF[KSS-1][KB-1];
			KF = QTREE_.KFLAG[KSS-1][KB-1];
                
			if ((KF < 3) && (QTREE_.KSP[KS-1] != KF)){
				goto d102;
			}  	   
		}	
		if (KB == KB0) {
			*TRACK_mod_.IBODY = KB; // a particula está dentro do corpo ou modulo KB
			*TRACK_mod_.MAT = PENGEOM_mod_.MATER[KB-1];// a particula está dentro do MATERial KB    
			return; 
		}else if (QTREE_.KDGHT[NXG-1][KB-1] > 1) {
			KB0 = KB; //o ponto está dentro de um submodulo
			goto d100;
		}else{
			*TRACK_mod_.IBODY = KB; //a particula esta dentro de um corpo ou modulo irmão
			*TRACK_mod_.MAT = PENGEOM_mod_.MATER[KB-1];
			return;
		}				
d102: ;
	}
	*TRACK_mod_.IBODY = *QTREE_.NBODYS + 1;
	*TRACK_mod_.MAT = 0;
 
}




void step2_(double *DS, double *DSEF, int *NCROSS){

	/*if (imprimiu == 0){
		printf("\n\nStep2\n\n");
		imprimiu++;
	}*/
	
	/*
	
	Esta sub-rotina lida com a parte geom�trica da simula��o de pista
   ��o A part�cula come�a no ponto (X, Y, Z) e percorre um comprimento
   DS na dire��o (U, V, W) dentro do material para onde se move. Quando
   a trilha deixa o material inicial, a part�cula � parada apenas
   depois de entrar no pr�ximo corpo material (regi�es vazias com MAT = 0 s�o
   cruzado automaticamente). Al�m disso, quando a part�cula chega de
   uma regi�o vazia, ele � interrompido logo ap�s entrar no primeiro material
   corpo.
 
   Valores de entrada (m�dulo TRACK_mod):
      X, Y, Z ... coordenadas do ponto inicial,
      U, V, W ... cossenos de dire��o do deslocamento,
      IBODY ..... corpo onde o ponto inicial est� localizado,
      MAT ....... material no corpo IBODY.
   NB: Quando uma trilha de part�culas � iniciada, as vari�veis ??IBODY e MAT
   deve ser definido chamando a sub-rotina LOCATE.
 
   Argumento de entrada:
      DS ........ comprimento do caminho a percorrer.
 
   Argumentos de sa�da:
      DSEF ....... percorreu o comprimento do caminho antes de deixar o inicial
                  material ou completar o salto (menos do que DS se o
                  a trilha atravessa uma interface),
      NCROSS .... = 0 se toda a etapa estiver contida no inicial
                    material,
                  .gt.0 se a part�cula cruzou uma interface, ou seja,
                    se entrou em um novo material.
 
   Valores de sa�da (m�dulo TRACK_mod):
      X, Y, Z ... coordenadas da posi��o final,
      IBODY ..... corpo onde o ponto final est� localizado,
      MAT ....... material em IBODY. O valor MAT = 0 indica que o
                  part�cula escapou do sistema.
 
   Valores de sa�da (m�dulo PENGEOM_mod):
      DSTOT ..... comprimento do caminho percorrido, incluindo segmentos de caminho no vazio
                  volumes.
      KSLAST .... quando NCROSS.ne.0, o valor de sa�da de KSLAST � o
         r�tulo da �ltima superf�cie cruzada pela part�cula antes
         entrar em um corpo material. KSLAST � usado para renderiza��o em 3D.
	
	
	*/

	static const int NS = 10000;
	static const int NS2M = 2*NS;
	
	*DSEF = 0.0e0;
	*PENGEOM_mod_.DSTOT = 0.0e0;
	*NCROSS = 0;
	*PENGEOM_mod_.KSLAST = 0;
	double DSRES;
	int KB1;
	//double S[NS2M];
	//int IS[NS2M];

	//double *S = (double *)malloc(NS2M*sizeof(double));
	//int *IS = (int *)malloc(NS2M*sizeof(int));

	
	int NST;
	double DSP;
	int KS1;
	int KF;
	int IERR;
	int IBODYL;
	int NERR;
	int KS;
	int KFLO;
	int SW;
	int MATL;
	double A, B, C;

	int NSC = 0; //Número de cruzamentos da superfície à frente da partícula.
	int NSCT;
	int MAT0;
	
//	printf("1\n");
	
//	printf("LVERB %d\n", *PENGEOM_mod_.LVERB);
//	printf(*PENGEOM_mod_.LVERB ? "true\n" : "false\n")
	for (int I = 1; I <= *QSURF_.NSURF; I++){
	 	QTREE_.KSP[I - 1]  = 0; //Ponteiros laterais das superfícies avaliadas.
	}
	
	//	printf("2\n");
		
		
	MAT0 = *TRACK_mod_.MAT; //Material Inicial
	
	if (*TRACK_mod_.MAT == 0){
		DSRES = 1.0e35; //No vácuo, as partículas voam livremente.
	} else {
		DSRES = *DS; // comprimento do camimho residual
	}
	
//		printf("3\n");
	//A partícula entra de fora do recinto.
	
	if (*TRACK_mod_.IBODY > *QTREE_.NBODYS){
		KB1 = *QTREE_.NBODYS;
		//stepsi2_(KB1, S, IS, NSC);
		stepsi2_(KB1, NSC);
	//	printf("4\n");
		if (NSC == 0)
			goto L300;
		NSCT = NSC;
		NST = QTREE_.KSURF[NXG -1][KB1 - 1];
			
		for (int KI = NSCT; KI >= 1; KI--){
	//		printf("5\n");
		//	printf("KI: %d\n", KI);
	//		printf("*PENGEOM_mod_.KSLAST : %d\n" ,*PENGEOM_mod_.KSLAST); 
			// a particula atravessa uma superficie
			*PENGEOM_mod_.KSLAST = IS[KI - 1];
	//		printf("6\n");
			if (QTREE_.KSP[*PENGEOM_mod_.KSLAST - 1] == 1)
				QTREE_.KSP[*PENGEOM_mod_.KSLAST - 1] = 2;
			else 
				QTREE_.KSP[*PENGEOM_mod_.KSLAST - 1] = 1;
			
			DSP = S[KI - 1];
			*DSEF = *DSEF + DSP;
			*PENGEOM_mod_.DSTOT = *PENGEOM_mod_.DSTOT + DSP;
			*TRACK_mod_.X = *TRACK_mod_.X + DSP * *TRACK_mod_.U;
			*TRACK_mod_.Y = *TRACK_mod_.Y + DSP * *TRACK_mod_.V;
			*TRACK_mod_.Z = *TRACK_mod_.Z + DSP * *TRACK_mod_.W;
			NSC = NSC -1;
			
			if (NSC > 0){
				for (int I = 1; I <= NSC; I++){
					S[I-1] = S[I-1] - DSP;
				}
			}
			
			for (int KSS = 1; KSS <= NST; KSS++){
				KS1 = QTREE_.KSURF[KSS-1][KB1-1];
				KF = QTREE_.KFLAG[KSS-1][KB1 -1];
				if ((KF < 3) && (QTREE_.KSP[KS1 -1] != KF))
					goto L101;
			}
				// A partícula entra no invólucro.
			L100:;
			steplb2_(KB1, IERR);
			// A partícula entra em um submódulo.
			if (IERR == -1){
				KB1 = *TRACK_mod_.IBODY;
				//stepsi2_(KB1, S, IS, NSC);
				stepsi2_(KB1, NSC);
				goto L100;
			} else{
				//A particula entrou em um corpo material
				if (*TRACK_mod_.MAT != 0){
					*NCROSS = 1;
	//				free(S);
	//				free(IS);
					return;
				} else{
					KB1 = *TRACK_mod_.IBODY;
					//stepsi2_(KB1, S, IS, NSC);
					stepsi2_(KB1, NSC);
				
					goto L200;
				}
						
			}
				
				//Neste ponto, o programa saiu do ciclo DO.
			L101:;
		}
		//	printf("6\n");
			goto L300;		
	}
	
	//Cruzamentos de superfície.
	
	IBODYL = *TRACK_mod_.IBODY;
//	printf("teste\n");
//	printf("LVERB: %d\n", *PENGEOM_mod_.LVERB);
	
//	printf(*PENGEOM_mod_.LVERB ? "true\n" : "false\n");
	if (*PENGEOM_mod_.LVERB)
		NERR = 0;
L102:;
	KB1 = *TRACK_mod_.IBODY;
	//stepsi2_(KB1, S, IS, NSC);
	stepsi2_(KB1, NSC);
	steplb2_(KB1, IERR);
	
	//Evidência de erros de arredondamento.
//printf("7\n");
	if (IERR != 0){
		if (NSC > 0){
			//Quando uma superfície está muito próxima, movemos a partícula além dela.
			if (S[NSC - 1] < 1e-10){
				*PENGEOM_mod_.KSLAST = IS[NSC -1];
				if (QTREE_.KSP[*PENGEOM_mod_.KSLAST-1] == 1){
					QTREE_.KSP[*PENGEOM_mod_.KSLAST-1] = 2;
				} else{
					QTREE_.KSP[*PENGEOM_mod_.KSLAST-1] = 1;
				}
				
				DSP = S[NSC-1];
				*TRACK_mod_.X = *TRACK_mod_.X + DSP * *TRACK_mod_.U;
				*TRACK_mod_.Y = *TRACK_mod_.Y + DSP * *TRACK_mod_.V;
				*TRACK_mod_.Z = *TRACK_mod_.Z + DSP * *TRACK_mod_.W;
				if (*TRACK_mod_.MAT == MAT0){
					*DSEF  = *DSEF + DSP;
					DSRES = DSRES - DSP;
				}
				
				*PENGEOM_mod_.DSTOT = *PENGEOM_mod_.DSTOT + DSP;
				NSC = NSC - 1;
				
				if (*TRACK_mod_.IBODY <= *QTREE_.NBODYS)
					goto L102;		
			}
		}
	//	printf("8\n");
       // printf("LVERB %d\n", *PENGEOM_mod_.LVERB);
		if (*PENGEOM_mod_.LVERB){
			
			NERR = NERR +1;
			if ((*QTREE_.NWARN < 100) && (*TRACK_mod_.MAT != 0)){
				
				printf("WARNING, STEP: Accidental undershot or r");
				
				for (int KSS = 1; KSS <= QTREE_.KSURF[NXG-1][KB1-1]; KSS++){
					KS = QTREE_.KSURF[KSS -1][KB1 -1];
					KFLO = QTREE_.KFLAG[KSS -1][KB1 -1];
					if (KFLO < 3){
						for (int KI = NSC; KI >= 1; KI--){
							if (KS == IS[KI - 1]){
								SW = S[KI - 1];
								goto L103;
							}
						}
						
						SW = 0.0e0;
L103:;
						fsurf2_(KS, A, B, C);
						if (KFLO == QTREE_.KSP[KS -1]){
							printf("KS, KFLO, KSP");
						} else{
							printf("KS, KFLO, KSP");
							*PENGEOM_mod_.KSLAST = KS;
						}	
					}
				}
				*QTREE_.NWARN = *QTREE_.NWARN + 1;	
			}
		
		}
		if (*TRACK_mod_.IBODY <= *QTREE_.NBODYS)
			goto L102;
	}
//	printf("9\n");
	if ((PENGEOM_mod_.KDET[*TRACK_mod_.IBODY - 1] != PENGEOM_mod_.KDET[IBODYL -1]) || (*TRACK_mod_.MAT != MAT0)){
		*NCROSS = 1;
		*DSEF = 0.0e0;
//		free(S);
//		free(IS);
		return;
	}
		   	
	//A particula permanece no mesmo material
	
	if ((*TRACK_mod_.MAT != 0) && (DSRES < S[NSC -1])){
		if (*TRACK_mod_.MAT == MAT0)
			*DSEF = *DSEF + DSRES;
		*PENGEOM_mod_.DSTOT = *PENGEOM_mod_.DSTOT + DSRES;
		*TRACK_mod_.X = *TRACK_mod_.X + DSRES * *TRACK_mod_.U;
		*TRACK_mod_.Y = *TRACK_mod_.Y + DSRES * *TRACK_mod_.V;
		*TRACK_mod_.Z = *TRACK_mod_.Z + DSRES * *TRACK_mod_.W;
//		free(S);
//		free(IS);
		return;	
	}
	
	
	// Nova posição
	
	L200:;
//	printf("10\n");
	if (NSC == 0 ){
		if (*TRACK_mod_.MAT == MAT0)
			*DSEF = *DSEF + DSRES;
		*PENGEOM_mod_.DSTOT = *PENGEOM_mod_.DSTOT + DSRES;
		*TRACK_mod_.X = *TRACK_mod_.X + DSRES * *TRACK_mod_.U;
		*TRACK_mod_.Y = *TRACK_mod_.Y + DSRES * *TRACK_mod_.V;
		*TRACK_mod_.Z = *TRACK_mod_.Z + DSRES * *TRACK_mod_.W;
		return;
	}
	NSCT = NSC;
	MATL = *TRACK_mod_.MAT;
	IBODYL = *TRACK_mod_.IBODY;
	for (int KI = NSCT; KI >= 1; KI--){
		//A etapa termina dentro do corpo
		if (DSRES < S[KI -1]){
			if (*TRACK_mod_.MAT == MAT0)
				*DSEF = *DSEF + DSRES;
			*PENGEOM_mod_.DSTOT = *PENGEOM_mod_.DSTOT + DSRES;
			*TRACK_mod_.X = *TRACK_mod_.X + DSRES * *TRACK_mod_.U;
			*TRACK_mod_.Y = *TRACK_mod_.Y + DSRES * *TRACK_mod_.V;
			*TRACK_mod_.Z = *TRACK_mod_.Z + DSRES * *TRACK_mod_.W;
//			free(S);
//			free(IS);
			return;
		}
		
		// A particula atravessa uma superfice
		*PENGEOM_mod_.KSLAST = IS[KI -1];
		if (QTREE_.KSP[*PENGEOM_mod_.KSLAST-1] == 1)
			QTREE_.KSP[*PENGEOM_mod_.KSLAST-1] = 2;
		else
			QTREE_.KSP[*PENGEOM_mod_.KSLAST-1] = 1;
		
		DSP = S[KI -1];
		*TRACK_mod_.X = *TRACK_mod_.X + DSP * *TRACK_mod_.U;
		*TRACK_mod_.Y = *TRACK_mod_.Y + DSP * *TRACK_mod_.V;
		*TRACK_mod_.Z = *TRACK_mod_.Z + DSP * *TRACK_mod_.W;
		if (*TRACK_mod_.MAT == MAT0){
			*DSEF = *DSEF + DSP;
			DSRES = DSRES - DSP;
		}
		*PENGEOM_mod_.DSTOT = *PENGEOM_mod_.DSTOT + DSP;
		NSC = NSC -1;
		if (NSC > 0){
			for (int I = 1; I <= NSC; I++){
				S[I-1] = S[I-1] - DSP;
			}
		}
		steplb2_(KB1, IERR);
		
L201:;
		KB1 = *TRACK_mod_.IBODY;
		if (IERR == -1){
			// A particula entrou em um submodulo
			//stepsi2_(KB1, S, IS, NSC);
			stepsi2_(KB1, NSC);
			steplb2_(KB1, IERR);
			goto L201;
		} 
		else if (IERR == 1){
			//A partícula deixa o corpo ou módulo.
			if (*TRACK_mod_.IBODY <= *QTREE_.NBODYS){
				//stepsi2_(KB1, S, IS, NSC);
				stepsi2_(KB1, NSC);
				steplb2_(KB1, IERR);
				goto L201;
			} 
			else{
				//A partícula sai do recinto.
				if (*TRACK_mod_.MAT != MATL)
					*NCROSS = *NCROSS + 1;
				goto L300;
			}
		}
		
		//A partícula continua voando quando entra em uma região vazia 
		if (*TRACK_mod_.MAT == 0){
			if (MATL == MAT0)
				*NCROSS = *NCROSS +1;
			MATL=0;
			DSRES = 1.0e35;
			goto L202;
			//A partícula continua voando quando entra em um novo corpo do 
			//mesmo material que não faz parte de um detector diferente ...
		} 
		else if (*TRACK_mod_.MAT == MATL){
			if (PENGEOM_mod_.KDET[*TRACK_mod_.IBODY - 1] == PENGEOM_mod_.KDET[IBODYL -1]){
				goto L202;
			} 
			else{
				*NCROSS = *NCROSS +1;
				return;
			}
					//.. e para quando penetra um novo corpo material ou umDetector
		} else {
			*NCROSS = *NCROSS +1;
//			free(S);
//			free(IS);
			return;
		}
L202:;
		//stepsi2_(KB1, S, IS, NSC);
		stepsi2_(KB1, NSC);
		goto L200;
		//Neste ponto, o programa saiu do ciclo DO.	
L203:;	
	}
	//A particula sai do recinto.
	L300:;
	
	DSP = 1.0e36;
	*TRACK_mod_.IBODY = *QTREE_.NBODYS +1;
	*TRACK_mod_.MAT = 0;
	if (*TRACK_mod_.MAT == MAT0)
		*DSEF = *DSEF + DSP;
	*PENGEOM_mod_.DSTOT = *PENGEOM_mod_.DSTOT + DSP;
	*TRACK_mod_.X = *TRACK_mod_.X + DSP * *TRACK_mod_.U;
	*TRACK_mod_.Y = *TRACK_mod_.Y + DSP * *TRACK_mod_.V;
	*TRACK_mod_.Z = *TRACK_mod_.Z + DSP * *TRACK_mod_.W;
//	free(S);
//	free(IS);
}


//void stepsi2_(int &KB, double *S, int *IS, int &NSC){
//void stepsi2_(int &KB, double S[NS2M], int IS[NS2M], int &NSC){
void stepsi2_(int &KB,  int &NSC){
/*Calcula as interseções da trajetória com o limite
	Superfícies do corpo KB. Os cruzamentos são adicionados à lista e
	classificados em ordem decrescente. 
	Esta sub-rotina funciona apenas quando chamada de dentro da sub-rotina STEP.*/
	
//stepsi_(KB, S, IS, NSC);
	//return;

	int KFL;
	int KS;
	double ABSA;
	double ABSB;
	double A, B, C;
	double FUZZL = 1.0e-12;
	double T1;
	double DISCR;
	double FUZZ;
	int IAMBIG;
	double R2A;
	double DELTA;
	double SH;
	double T2;
	int KMAX;
	double SMAX;
	int KKMAX;
	
	
	
	//Determine cruzamentos de superfície.
	
	
	for ( int KSS = 1; KSS <= QTREE_.KSURF[NXG - 1][KB -1]; KSS++){
		//printf("\nKSURF: %d\n",  QTREE_.KSURF[NXG - 1][KB -1]);
		
		/*As interseções com uma determinada superfície são calculadas apenas uma vez.
		O ponteiro lateral de uma superfície deve ser alterado cada vez que o a superfície está cruzada.*/
		
		KFL = QTREE_.KFLAG[KSS-1][KB-1];
		if (KFL > 4)
			goto L100;
		KS = QTREE_.KSURF[KSS-1][KB-1];
		if (QTREE_.KSP[KS-1] != 0) 
			goto L100;
		fsurf2_(KS, A, B, C);
		ABSA = fabs(A);
		ABSB = fabs(B);
		
		// Plano, unica raiz
		if(ABSA < 1.0e-36){
			if (ABSB > 0.0e0){
				if (C < -FUZZL){
					QTREE_.KSP[KS-1] = 1;
				} 
				else if (C > FUZZL){
						QTREE_.KSP[KS -1] = 2;
				} else{
					if (B < 0.0e0){
						QTREE_.KSP[KS -1] = 1;
					} else{
						QTREE_.KSP[KS -1] = 2;
					}
					goto L100;	
				}
				T1 = -C/B;
				if (T1 > 0.0e0){
					NSC = NSC +1;
					IS[NSC -1] = KS;
					S [NSC -1] = T1;
				}
			}
			else {
				if (C < 0.0e0){
					QTREE_.KSP[KS - 1] = 1;
				}
				else{
					QTREE_.KSP[KS - 1] = 2;
				}
			}
				
			// Superficie não plana, duas raizes			
		} 
		else{
			DISCR = B*B-4.0e0*A*C;
			FUZZ = FUZZL*DISCR/ABSA;
			if (C < -FUZZ){
				IAMBIG = 0;
				QTREE_.KSP[KS -1] = 1;
			} else if (C > FUZZ){
				IAMBIG = 0;
				QTREE_.KSP[KS -1] = 2;
			} else{
				IAMBIG = 1;
				if (B < 0.0e0){
					QTREE_.KSP[KS -1] = 1;
				}else {
					QTREE_.KSP[KS -1] = 2;
				}		
			}
			
			if (DISCR < 1.0e-36)
				goto L100;
			
			if (IAMBIG == 0){
				R2A = 0.5e0/A;
				DELTA = sqrt(DISCR)*fabs(R2A);
				SH = -B*R2A;
				T1 = SH - DELTA;
				if (T1 > 0.0e0) {
					NSC = NSC +1;
					IS[NSC -1] = KS;
					S[NSC - 1] = T1;
				}
				T2 = SH + DELTA;
				if (T2 > 0.0e0){
					NSC = NSC +1;
					IS[NSC -1] = KS;
					S[NSC - 1] = T2;
				}
			} else{
				if (B*A < 0.0e0){
					R2A = 0.5e0/A;
				   	DELTA = sqrt(DISCR)*fabs(R2A);
				   	SH = -B*R2A;
				   	T2 = SH + DELTA;
				   	NSC = NSC +1;
				   	IS[NSC -1] = KS;
				   	S[NSC - 1] = fmax(T2, 0.0e0);	
				}	
			}	
		}		
L100:;		  		
	}

	
	
	//Classifique as distâncias da superfície em ordem decrescente.
   // printf("\nNSC %d\n", NSC);
	if (NSC > 1){
		for (int KI = 1; KI <= (NSC - 1); KI++){
			SMAX = S[KI -1];
			KMAX = KI;
			for ( int KJ = (KI+1); KJ <= NSC; KJ++){
				if (S[KJ -1] > SMAX){
					SMAX = S[KJ -1];
			   	   	KMAX = KJ;
				}	
			}
			if (KMAX != KI){
				SMAX = S[KI -1];
				S[KI -1] = S[KMAX -1];
				S[KMAX -1] = SMAX;
   	   			KKMAX = IS[KI -1];
   	   			IS[KI-1] = IS[KMAX -1];
				IS[KMAX -1] = KKMAX; 
			}	
		}
	}	
}


void steplb2_( int &KB, int &IERR){
	
	/*Ajuda a encontrar o corpo ou módulo que contém os ponteiros laterais fornecidos para as superfícies analisadas.
	 A subrotina STEPLB funciona apenas quando invocada de dentro da sub-rotina STEP. 
	 Ele se move através da árvore de módulos a única etapa*/
	 
	 
	 int NLBOD;
	 int KBS;
	 int KS;
	 int KF;
	 int KBD;
	 
	//Analisa o corpo ou o modulo atual. 
	
	if (QBODY_.KBOMO[KB -1] == 0){
		//Corpo
		NLBOD = QBODY_.KBODY[NXG - 1][KB - 1];
		if (NLBOD > 0) {
			for (int KBB = 1; KBB <= NLBOD; KBB++){
				KBS = QBODY_.KBODY[KBB - 1][KB - 1];
				for (int KSS = 1; KSS <= QTREE_.KSURF[NXG -1 ][KBS -1]; KSS++){
					KS = QTREE_.KSURF[KSS -1 ][KBS -1];
					KF = QTREE_.KFLAG[KSS -1][KBS -1];
					if ((KF < 3) && (QTREE_.KSP[KS -1] != KF))
						goto L100;		
				}
				*TRACK_mod_.IBODY = KBS;
				if (QTREE_.KDGHT[NXG - 1][*TRACK_mod_.IBODY - 1] > 1){
					IERR = -1; //A particula está dentro de um módulo irmã
				} else{
					IERR = 0; //A particula está dentro de um corpo irma	
					*TRACK_mod_.MAT = PENGEOM_mod_.MATER[*TRACK_mod_.IBODY - 1];
				}
				return;
L100:;
			}
		}
		
		for ( int KSS = 1; KSS <= QTREE_.KSURF[NXG -1][KB -1]; KSS++){
			KS = QTREE_.KSURF[KSS -1 ][KB -1];
			KF = QTREE_.KFLAG[KSS -1][KB -1];
			if ((KF < 3) && (QTREE_.KSP[KS -1] != KF))
				goto L300;
		}
		*TRACK_mod_.IBODY = KB;
		IERR = 0;
		*TRACK_mod_.MAT = PENGEOM_mod_.MATER[*TRACK_mod_.IBODY - 1];
		return;	
	} 
	else{	
		// Modulo
		for (int KBB = 1; KBB <= QTREE_.KDGHT[NXG -1][KB -1]; KBB++){
			KBD = QTREE_.KDGHT[KBB -1][KB -1];
			for (int KSS = 1; KSS <= QTREE_.KSURF[NXG -1 ][KBD -1]; KSS++){
				KS = QTREE_.KSURF[KSS -1 ][KBD -1];
				KF = QTREE_.KFLAG[KSS -1][KBD -1];
				if ((KF < 3) && (QTREE_.KSP[KS-1] != KF))
					goto L200;		
			}
			*TRACK_mod_.IBODY = KBD;
			if (KBD == KB){
				IERR = 0; //A partícula permanece dentro do módulo atual.
				*TRACK_mod_.MAT = PENGEOM_mod_.MATER[*TRACK_mod_.IBODY - 1];
			} else{
				if (QTREE_.KDGHT[NXG -1][KBD -1] > 1){
					IERR = -1; // A particula esta dentro de um submodulo
				} else{
					IERR = 0; //A particula esta dentro de um corpo simples;
					*TRACK_mod_.MAT = PENGEOM_mod_.MATER[*TRACK_mod_.IBODY - 1];
				}
			}
			return;
			L200:;		
		}	
	}
	
	//A partícula está fora do corpo ou módulo atual.
	L300:
	IERR=1;
	*TRACK_mod_.IBODY = QTREE_.KMOTH[KB - 1];
	if (*TRACK_mod_.IBODY == 0){
		*TRACK_mod_.IBODY = *QTREE_.NBODYS +1;
		*TRACK_mod_.MAT = 0;
	}
}
	










	



void geomin2_(double *PARINP, int *NPINP, int *NMATG, int *NBOD, FILE *IRD, FILE *IWR){
//	printf("\ngeomin2\n");
	/*
	Lê o arquivo de definição de geometria e configura os arrays usados para rastrear partículas através do sistema.
	
	Argumentos de entrada:
	PARINP .... array contendo parâmetros opcionais, que podem substitui os inseridos no arquivo de entrada.
    NPINP ..... índice do maior componente de PARINP que carrega um valor de parâmetro efetivo. 
				Quando NPINP é negativo ou zero, os valores numéricos do arquivo de geometria são mantidos inalterados.
    IRD ....... unidade de entrada do arquivo de definição de geometria (aberto em o programa principal).
    IWR ....... unidade de saída (aberta no programa principal). Na saída este arquivo contém uma duplicata do arquivo de definição
				com os valores dos parâmetros efetivos e com elementos dos vários tipos rotulados em estritamente crescente pedido. 
				Esta parte do arquivo de saída descreve o geometria real usada na simulação. Depois do fim Linha do bloco de definição de geometria, sub-rotina
				GEOMIN escreve um relatório detalhado com a estrutura de a árvore de módulos e com informações sobre redundantes Superfícies(duplicadas).
    
    Argumentos de saída C:
    NMATG ..... número de diferentes materiais em corpos inteiros (excluindo ing regiões vazias).
    NBOD ...... Número de corpos definidos.
    
	**** Alias, material, mãe, filhas, superfícies e ponteiros laterais de corpos. 
		O último segundo componente das matrizes duplas é o número de componentes usados. 
		Por exemplo, KSURF (KB, NXG) é o número de superfícies que limitam o corpo ou módulo KB. 
		O número de valores armazenado na matriz KFLAG (KB, -), entretanto, é fornecido por KSURF (KB, NXG).
		DO KB = 1, NB
	**** A matriz KDET (KB) é usada para rotular corpos que são partes de detectores de impacto.
		Os detectores devem ser definidos no programa principal, após a chamada à sub-rotina GEOMIN
	*/
	
	bool LOPCLO;
	
	char CHR[6];
	char CCR[2];
	char CNL[2];
	char C5[6];
	char C5C[6];
	char C4[5];
	char C1[2];
//	char CA[36] = {'0', '1', '2', '3', '4', '5', '6', '7', '8', '9', 'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 
  //                 'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P', 'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'X', 'Y', 'Z'};
    char CA[] = "0123456789ABCDEFBHIJKLMNOPQRSTUVWXYZ";
	
	char GFILE[13];
	char BLINE[100];
	char DEF[73];
	
	char LKEYW[15];
	char LANGLE[9];
	char APOIO[100];
	
	static const char LNUL[] ="00000000";
	static const char LSUR[]="SURFACE ";
	static const char LIND[]="INDICES=";
	static const char LBOD[]="BODY    ";
	static const char LMAT[]="MATERIAL";
	 static const char LMOD[]="MODULE  ";
     static const char LOPEN[]=")       ";
	 static const char LEND[]="END     ";
	 static const char LONE[]="11111111";
     static const char LXSC[]="X-SCALE=";
	 static const char LYSC[]="Y-SCALE=";
	 static const char LZSC[]="Z-SCALE=";
     static const char LTHE[]="  THETA=";
	 static const char LPHI[]="    PHI=";
	 static const char LOME[]="  OMEGA=";
     static const char LXSH[]="X-SHIFT=";
	 static const char LYSH[]="Y-SHIFT=";
	 static const char LZSH[]="Z-SHIFT=";
     static const char LAXX[]="    AXX=";
	 static const char LAXY[]="    AXY=";
	 static const char LAXZ[]="    AXZ=";
     static const char LAYY[]="    AYY=";
	 static const char LAYZ[]="    AYZ=";
	 static const char LAZZ[]="    AZZ=";
     static const char LAX[]="     AX=";
	 static const char LAY[]="     AY="; 
	 static const char LAZ[]="     AZ=";
     static const char LA0[]="     A0=";
	 static const char LDEG[]=") DEG   ";
	 static const char LRAD[]=") RAD   ";
     static const char LCLON[]="CLONE   ";
	 static const char LSUA[]="SURFACE*";
	 static const char LINC[]="INCLUDE ";
     static const char LINA[]="INCLUDE*";
	 static const char LFIL[]="   FILE=";
	 static const char ZEROS[]="0000000000000000000000000000000000000000000000000000000000000000\n";
	 
	 int LARRAY[9];
	 int KQ[5];
	 int KM[NS];
	 int IDESC[NB];
	 int IDONE[NB];
	 int ISCL[NS];
	 int IBCL[NB];
	 int IBOR[NB];
	 char DEFB[NB][73];
	 char DEFS[NS][73];
	 char ALIAS[NS][6];
	 char ALIAB[NB][6];
	 int KSFF[NS];
	 int NSFF;
	 int ICLONE;
	 int KEEPL;
	 int NINCL;
	 int IRI;
	 int IMODE;
	 int ICHPAR;
	 int IMAT;
	 int INDS;
	 int KS;
	 int KST;
	 int KST0;
	 int KN1;
	 int KN2;
	 int KN3;
	 int KSURF1;
	 int KSURF2;
	 int KDT;
	 int KMIN;
	 int KBMIN;
	 int KSAVE;
	 int KBD;
	 int IDGHT;
	 int K2M;
	 int KNS;
	 int K2;
	 int KBMOTH;
	 int KORIG;
	 int ND;
	 int NDC;
	 int KDG;
	 int IN;
	 int ID;
	 int KSD;
	 int KBO;
	 int MLESS;
	 int KBENC;
	 int NT;
	 int IWRITE;
	 int NSB;
	 int KSLP;
	 int KSL;
	 int KFL;
	 int KFLP;
	 int NBU;
	 int NSU;
	 int NSE;
	 int KB;
	 int KBB;
	 int KBI;
	 int KB1;
	 int KF;
	 
	 
	 double XSCALE;
	 double YSCALE;
	 double ZSCALE;
	 double OMEGA;
	 double THETA;
	 double PHI;
	 double XSHIFT;
	 double YSHIFT;
	 double ZSHIFT;
	 double VALUE;
	 
	 double QXX;
	 double QXY;
	 double QXZ;
	 double QYY;
	 double QYZ;
	 double QZZ;
	 double QX;
	 double QY;
	 double QZ;
	 double Q0;
	 
	 
	double TSTL ;
	double TSTU;
	double TOL;
	double F;
	double FM;
	double FP;
	double FFP;
	double RFFP;
	double TST;
	double ISPF;
	double TF;
	
	
	 
	CCR[1] = char(13);
	CNL[1] = char(10); 
	
	NSFF = 0;

	
	
 	FILE* IR = fopen("quadril.geo", "r");
	if (IR == NULL){
		printf("Não foi possível abrir o arquivo de geometria");
		exit(0);
	}
	
	FILE* IW = fopen("geometry2.geo", "w");
	if (IW == NULL){
		printf("Não foi possível abrir o arquivo de geometria-res");
		exit(0);
	}

	//FILE* IR = IRD;
//	FILE* IW = IWR;


	
	for (int KS = 1; KS <= NS; KS++){
		strcpy(ALIAS[KS -1], "    0");
		KSFF[KS - 1] = 0;
		QSURF_.AXX[KS - 1] = 0.0e0;
		QSURF_.AXY[KS - 1] = 0.0e0;
		QSURF_.AXZ[KS - 1] = 0.0e0;
		QSURF_.AYY[KS - 1] = 0.0e0;
		QSURF_.AYZ[KS - 1] = 0.0e0;
		QSURF_.AZZ[KS - 1] = 0.0e0;
		QSURF_.AX[KS - 1] = 0.0e0;
		QSURF_.AY[KS - 1] = 0.0e0;
		QSURF_.AZ[KS - 1] = 0.0e0;
		QSURF_.A0[KS - 1] = 0.0e0;
		QSURF_.KPLANE[KS - 1] = 0;		
	}
	
	for (int KB = 1; KB <= NB; KB++){
		PENGEOM_mod_.KDET[KB -1] = 0; //KDET (KB) .ne.0 se o corpo KB fizer parte de um detector.
		strcpy(ALIAB[KB-1], "    0"); //Aliases de corpo (rótulos de usuário).
		PENGEOM_mod_.MATER[KB -1] = 0; //Material no corpo KB.
		QBODY_.KBOMO[KB - 1] = 0; //0 para corpos, 1 para módulos.
		QTREE_.KMOTH[KB - 1] = 0; //Mãe do corpo KB (deve ser única).

		for (int K2 = 1; K2 <= NXG; K2++){
			QBODY_.KBODY[K2 - 1][KB - 1] = 0; //Limitando os corpos do corpo KB (mantenha-o pequeno)
			QTREE_.KDGHT[K2 -1][KB - 1] = 0; //Filhas do módulo KB (deve ser pequeno).
			QTREE_.KSURF[K2 - 1][KB -1] = 0; //Limitando as superfícies do corpo KB.
			QTREE_.KFLAG[K2 - 1][KB - 1] = 5;//Ponteiros laterais de superfície de corpos materiais.
			
			/*
			KFLAG (KB, KS) = 1, se KS é uma superfície limitante de KB e KB está dentro de KS (ou seja, ponteiro lateral = -1).
			KFLAG (KB, KS) = 2, se KS é uma superfície limitante de KB e KB está fora KS (ou seja, ponteiro lateral = + 1).
			KFLAG (KB, KS) = 3, se KB for um corpo e KS não limitar diretamente KB, mas KS aparece na definição de um corpo que limita KB.
			KFLAG (KB, KS) = 4, se KB for um módulo e KS limitar uma de suas filhas (corpos e submódulos), mas não aparece no Definiçãode KB.
			KFLAG (KB, KS) = 5, caso contrário.
			*/
		}	
	}
	
	//	printf("linha 1121\n");
	
	*NMATG = 0;
	*QSURF_.NSURF = 0;
	*QTREE_.NBODYS = 0;
	ICLONE = 0;
	KEEPL = 0;
	*QTREE_.NWARN = 0;
	
	
	
	//Lendo o arquivo de entrada de geometria.
	
	
	

	
/*	IR = IRD;
	IW = IWR;
	NINCL = 0;
	C1 = CA [NINCL+1];
	if (IW == IR){
		fprintf(IW, "SUBROUTINE GEOMIN. Input arguments.");
		fprintf(IW, "IRD = %d IWR = %d));
		fprintf(IW, "*** The input and output units must be different.");
		exit(0)(0);
	}
	
	 //Procure uma unidade de E / S que não esteja aberta, para ser conectada possíveis arquivos INCLUÍDOS.	
	 
	 IRI = max(IR, IW) + 1;*/
	 
	NINCL = 0;
	C1[0] = CA[NINCL+1-1];
	C1[1] = '\0';
	
	 

L9:;

	/*FILE* file = fopen("penmain-res.dat", "r");
	if (file == NULL){
		IRI = IRI + 1;
		goto L9;
		
	}*/
	
L1:;
	
	
	fgets(BLINE, sizeof(BLINE), IR);
//	printf(" linha 1165 \n ");
//	printf(" BLINE %s \n", BLINE);
	

//	printf("LKEYW: %s", BLINE.substr);

//	strcpy(LKEYW, extrairString(BLINE, 0, 8)); 
	//extrairString(LKEYW, BLINE, 0, 8);
	
	extrairString(LKEYW, BLINE, 0, 8);

	
	
//	printf("LKEYW: %s\n", LKEYW);
//	strcpy(LKEYW, "\0");
	
	
//	strcpy
	
//	printf(" linha 1169: LKEYW: %s \n ", strncpyy(LKEYW, BLINE+0, 8));

	if (strcmp(LKEYW, LNUL)){
		fputs(BLINE, IW);
		goto L1;
	} else{
		fputs(ZEROS, IW);
	}
	
		
	
L2:;
//printf(" linha 1178 ");

	fgets(BLINE, sizeof(BLINE), IR);
	
//	printf(" linha 1188 ");
	if (comentario(BLINE)){
	//	printf("linha 1197\n ");
		fputs(BLINE, IW);
		goto L2;   
	}
	//	printf("linha 1201\n ");
	extrairString(LKEYW, BLINE, 0, 8); 
	extrairString(C4, BLINE, 9, 4);
//	extrairString(C4, BLINE, 9, 4);
	
	 //printf(" linha 1198 LKEYW: %s , C4: %s \n", LKEYW, C4);
	
	if ((!strcmp(LKEYW, LSUR)) || (!strcmp(LKEYW, LSUA))){
	
		goto L100; //superficies
		
	}else if (!strcmp(LKEYW, LBOD))
		goto L200; //corpos
	else if (!strcmp(LKEYW, LMOD))
		goto L300; //modulos
	else if (!strcmp(LKEYW, LCLON))
		goto L400;
	else if ((!strcmp(LKEYW, LINC)) || (!strcmp(LKEYW, LINA)))
		goto L500;
	else if (!strcmp(LKEYW, LEND)){
	
		goto L600;
		
	}
	else {
		fputs(BLINE, IW);
		fputs("*** What do you mean?", IW);
		exit(0);
	}
//	char tmpC4[5];
	
	//Superficies
	
	
L100:;
	
	//if ((strcmp(strncpyy(APOIO, BLINE+8, 1), "(" )) || (strcmp(strncpyy(APOIO, BLINE+13, 1), ")"))){
	if ((BLINE[8] != '(' ) || (BLINE[13] != ')')){
		fputs(BLINE, IW);
		fputs("*** Incorrect label format.", IW);
		exit(0);
	}
 	extrairString(C4, BLINE, 9, 4);
 
 	if (KEEPL == 1){
		strcat(strcpy(C5, C4), "0");
		C5[5] = '\0';
 
	 }else {
 		strcat(strcpy(C5, C4), C1);
	 }

 	if (*QSURF_.NSURF > 0){
 		for (int KS0 = 1; KS0 <= *QSURF_.NSURF; KS0++){
 			extrairString(APOIO, ALIAS[KS0-1], 0, 5);
 			if (!strcmp(C5, APOIO)){
				fputs(BLINE, IW);
				fputs("*** Same label for two surfaces.", IW);
				exit(0); 
			}
		}
	}
	*QSURF_.NSURF = *QSURF_.NSURF + 1;
	KS = *QSURF_.NSURF;
	fprintf(IW, "%s(%4d) linha 340\n", LKEYW, KS);
	//fputs("Linha 342\n", IW);
	if (KS > NS){
		fputs("*** The parameter NS must be increased.", IW);
		exit(0);
	}
	//escrita no DEF fputs("Linha 345\n", IW);
	strcpy(ALIAS[KS-1], C5);
	
	//printf("C5: %s, ALIAS: %s", C5,ALIAS[KS-1] );
	
//	printf("antes KSFF %d    NSFF: %d\n", 	KSFF[KS-1], NSFF );
	
	if (!strcmp(LKEYW, LSUA)){
		KSFF[KS-1] = 1;
		NSFF = NSFF + 1;
	}
	
//	printf("depois KSFF %d    NSFF: %d\n", 	KSFF[KS-1], NSFF );
//	printf("1266\n");
	
	//Indices

L190:;

	fgets(BLINE, sizeof(BLINE), IR);
	if (comentario(BLINE)){

		   fputs(BLINE, IW);
		   goto L190;   
	}
	extrairString(LKEYW, BLINE, 0, 8); 
	CHR[0] = BLINE[8]; CHR[1] = BLINE[11]; CHR[2] = BLINE[14]; CHR[3] = BLINE[17]; CHR[4] = BLINE[20]; CHR[5] = BLINE[23];
//	KQ[0] = atoi(BLINE+8);	KQ[1] = atoi(BLINE+12);KQ[2] = atoi(BLINE+15);	KQ[3] = atoi(BLINE+18);	KQ[4] = atoi(BLINE+22);

//	memcpy(APOIO,BLINE+9,2);// KQ[1] = atoi(strncpyy(APOIO,BLINE+12,2));  KQ[2] =  atoi(strncpyy(APOIO,BLINE+15,2)); 
//	KQ[0] = atoi(APOIO);
	
//	memcpy(APOIO,BLINE+12,2);// KQ[1] = atoi(strncpyy(APOIO,BLINE+12,2));  KQ[2] =  atoi(strncpyy(APOIO,BLINE+15,2)); 
//	KQ[1] = atoi(APOIO);
//	memcpy(APOIO,BLINE+15,2);// KQ[1] = atoi(strncpyy(APOIO,BLINE+12,2));  KQ[2] =  atoi(strncpyy(APOIO,BLINE+15,2)); 
//	KQ[2] = atoi(APOIO);
//	memcpy(APOIO,BLINE+18,2);// KQ[1] = atoi(strncpyy(APOIO,BLINE+12,2));  KQ[2] =  atoi(strncpyy(APOIO,BLINE+15,2)); 
//	KQ[3] = atoi(APOIO);
//	memcpy(APOIO,BLINE+21,2);// KQ[1] = atoi(strncpyy(APOIO,BLINE+12,2));  KQ[2] =  atoi(strncpyy(APOIO,BLINE+15,2)); 
//	KQ[4] = atoi(APOIO);

	
	
	
	
	
//	KQ[3] = atoi(strncpyy(APOIO,BLINE+18,2)); KQ[4] = atoi(strncpyy(APOIO,BLINE+21,2));

/*	strncpyy(APOIO,BLINE+9,2);
	APOIO[2] = '\0';
	KQ[0] = atoi(APOIO);
	strncpyy(APOIO,BLINE+12,2);
	APOIO[2] = '\0';
	KQ[1] = atoi(APOIO);
	strncpyy(APOIO,BLINE+15,2);
	APOIO[2] = '\0';
	KQ[2] = atoi(APOIO);
	strncpyy(APOIO,BLINE+18,2);
	APOIO[2] = '\0';
	KQ[3] = atoi(APOIO);
	strncpyy(APOIO,BLINE+21,2);
	APOIO[2] = '\0';
	KQ[4] = atoi(APOIO);*/

		extrairString(APOIO,BLINE, 9,2); 
			KQ[0] = atoi(APOIO); 
	extrairString(APOIO,BLINE, 12,2);
	KQ[1] = atoi(APOIO);
	extrairString(APOIO,BLINE, 15,2);
	KQ[2] = atoi(APOIO); 
	extrairString(APOIO,BLINE, 18,2);
		KQ[3] = atoi(APOIO); 
	extrairString(APOIO,BLINE, 21,2);
	KQ[4] = atoi(APOIO); 




//	KQ[0] = atoi(strcpy(APOIO, strncpyy(APOIO,BLINE+9,2))); KQ[1] = atoi(strcpy(APOIO, strncpyy(APOIO,BLINE+12,2))); KQ[2] = atoi(strcpy(APOIO, strncpyy(APOIO,BLINE+15,2)));
//	KQ[3] = atoi(strcpy(APOIO, strncpyy(APOIO,BLINE+18,2))); KQ[4] = atoi(strcpy(APOIO, strncpyy(APOIO,BLINE+21,2)));
//	strcpy(&CHR[0], extrairString(BLINE, 9,9)); strcpy(&CHR[1], extrairString(BLINE, 12,12)); strcpy(&CHR[2], extrairString(BLINE, 15,15)); 
//	str	
	
//	KQ[0] = atoi(strcpy(APOIO, strncpyy(APOIO,BLINE+9,2))); KQ[1] = atoi(strncpyy(APOIO,BLINE+12,2));  KQ[2] =  atoi(strncpyy(APOIO,BLINE+15,2)); 
//	KQ[3] = atoi(strncpyy(APOIO,BLINE+18,2)); KQ[4] = atoi(strncpyy(APOIO,BLINE+21,2));
//	strcpy(&CHR[0], extrairString(BLINE, 9,9)); strcpy(&CHR[1], extrairString(BLINE, 12,12)); strcpy(&CHR[2], extrairString(BLINE, 15,15)); 
//	strcpy(&CHR[3], extrairString(BLINE, 18,18)); strcpy(&CHR[4], extrairString(BLINE, 21,21)); strcpy(&CHR[5], extrairString(BLINE, 24,24));
//	KQ[0] =  atoi(extrairString(BLINE, 10,11)); KQ[1] = atoi(extrairString(BLINE, 13,14));  KQ[2] =  atoi(extrairString(BLINE, 16,17)); 
//	KQ[3] = atoi(extrairString(BLINE, 19,20)); KQ[4] = atoi(extrairString(BLINE, 22,23));
	
	//Verificando os parenteses e virgulas
/*	if ((strcmp(&CHR[0], "(")) || (strcmp(&CHR[1], ",")) || (strcmp(&CHR[2], ",")) || (strcmp(&CHR[3], ",")) ||
		(strcmp(&CHR[4], ",")) || (strcmp(&CHR[5], ")"))){
		   fputs("*** Incorrect format of surface indices.", IW);
		   exit(0);   
	}*/
	
	
	
//	printf("KQ[0]: %d KQ[1]: %d KQ[2]: %d KQ[3]: %d KQ[4]: %d --- CHR: %s\n", KQ[0], KQ[1], KQ[2],KQ[3],KQ[4],CHR);
	
	fprintf(IW, "%s%c%2d%c%2d%c%2d%c%2d%c%2d%c\n", LKEYW,CHR[0],KQ[0],CHR[1],KQ[1],CHR[2],KQ[2],CHR[3],KQ[3],CHR[4],KQ[4],CHR[5]);
	
	if ((CHR[0] != '(') || (CHR[1] != ',') || (CHR[2] != ',') || (CHR[3] != ',') ||
		(CHR[4] != ',') || (CHR[5] != ')') ){
		   fputs("*** Incorrect format of surface indices.", IW);
		   exit(0);   
	}
	
	
	// Verificando os valores dos indices.
//	IMODE = max({abs(KQ[0]), abs(KQ[1]), abs(KQ[2]), abs(KQ[3]), abs(KQ[4])});
	IMODE = max(abs(KQ[0]), max(max(max(abs(KQ[1]), abs(KQ[2])), abs(KQ[3])), abs(KQ[4])));
//	IMODE = max(abs(KQ[0]), abs(KQ[1]));
//	IMODE = max(IMODE, abs(KQ[2]));
//	IMODE = max(IMODE, abs(KQ[3]));
//	IMODE = max(IMODE, abs(KQ[4]));
	if ((IMODE > 1) || (strcmp(LKEYW, LIND))) {
		fputs("*** Incorrect surface indices.", IW);
		exit(0);
	}
	if (IMODE == 0)
		goto L103;
	
	
	//forma reduzida
	 XSCALE = 1.0e0;
	 YSCALE = 1.0e0;
	 ZSCALE = 1.0e0;
	 OMEGA = 0.0e0;
	 THETA = 0.0e0;
	 PHI = 0.0e0;
	 XSHIFT = 0.0e0;
	 YSHIFT = 0.0e0;
	 ZSHIFT = 0.0e0;
	 
	 //Parametros da quadratica
	 
L101:;

	fgets(BLINE, sizeof(BLINE), IR);
	if (comentario(BLINE)){
		   fputs(BLINE, IW);
		   goto L101;   
	}
	
	extrairString(LKEYW, BLINE, 0, 8); 
	extrairString(APOIO, BLINE, 9,22);
	VALUE = atof(APOIO);
	extrairString(APOIO, BLINE, 32, 4);
	ICHPAR = atoi(APOIO);
    extrairString(LANGLE, BLINE, 36, 8);
	
//	printf("linha 1343: VALUE:%f  ICHPAR:%4d   LANGLE:%s\n", VALUE, ICHPAR, LANGLE);
	
//	printf("linha 1342\n");
	
	/*VALUE = atof(extrairString(BLINE, 9,30));
	ICHPAR = atoi(extrairString(BLINE, 32,35));
	strcpy(LANGLE, extrairString(BLINE, 36,43));*/
	
//	printf("linha 1348: LKEYW: %s\n", LKEYW);
	
	if (!strcmp(LKEYW, LNUL))
		goto L102;
	if (ICHPAR > 0){
		if (ICHPAR <= *NPINP){
			VALUE = PARINP[ICHPAR-1];
			ICHPAR = -ICHPAR; //Desativa a opção de alteração de parâmetro.
		} else{
			fprintf(IW, "%s(%f,%4d,%s)\n", LKEYW, VALUE, ICHPAR, LANGLE);
		//	fputs("linha 398", IW);
			fputs("*** NPINP is too small (check PARINP).", IW );
			exit(0);
		}
	}
	
	//Parâmetros de transformação.
	
	if (!strcmp(LKEYW, LXSC)){
		fprintf(IW, "%s(%.14E,%4d%s\n", LKEYW, VALUE, ICHPAR, LOPEN);
		//fputs("linha 406\n", IW);
		if (VALUE < 1e-15){
			fputs("*** Scale factor less than 1.0D-15.", IW);
			exit(0);
		}
		XSCALE = VALUE;
	} else if (!strcmp(LKEYW, LYSC)){
		//fputs("linha 414\n", IW);
		fprintf(IW, "%s(%.14E,%4d%s\n", LKEYW, VALUE, ICHPAR, LOPEN);
		if (VALUE < 1e-15){
			fputs("*** Scale factor less than 1.0D-15.", IW);
			exit(0);
		}
		YSCALE = VALUE;
	} else if (!strcmp(LKEYW, LZSC)){
		//fputs("linha 422\n", IW);
		fprintf(IW, "%s(%.14E,%4d%s\n", LKEYW, VALUE, ICHPAR, LOPEN);
		if (VALUE < 1e-15){
			fputs("*** Scale factor less than 1.0D-15.", IW);
			exit(0);
		}
		ZSCALE = VALUE;
	} else if (!strcmp(LKEYW, LOME)){
		if (!strcmp(LANGLE, LRAD)){
			//fputs("linha 431\n", IW);
			fprintf(IW, "%s(%.14E,%4d%s\n", LKEYW, VALUE, ICHPAR, LRAD);
			OMEGA = VALUE;
		} else{
			//fputs("linha 435\n", IW);
			fprintf(IW, "%s(%.14E,%4d%s\n", LKEYW, VALUE, ICHPAR, LDEG);
			OMEGA = VALUE * PI / 180.0e0;
		}
	} else if (!strcmp(LKEYW, LTHE)){
		if (!strcmp(LANGLE, LRAD)){
			//fputs("linha 441\n", IW);
			fprintf(IW, "%s(%.14E,%4d%s\n", LKEYW, VALUE, ICHPAR, LRAD);
			THETA = VALUE;
		} else{
			//fputs("linha 445\n", IW);
			fprintf(IW, "%s(%.14E,%4d%s\n", LKEYW, VALUE, ICHPAR, LDEG);
			THETA = VALUE * PI / 180.0e0;
		}
	} else if (!strcmp(LKEYW, LPHI)){
		if (!strcmp(LANGLE, LRAD)){
			//fputs("linha 451\n", IW);
			fprintf(IW, "%s(%.14E,%4d%s\n", LKEYW, VALUE, ICHPAR, LRAD);
			PHI = VALUE;
		} else{
			//fputs("linha 455\n", IW);
			fprintf(IW, "%s(%.14E,%4d%s\n", LKEYW, VALUE, ICHPAR, LDEG);
			PHI = VALUE * PI / 180.0;
		}
	} else if (!strcmp(LKEYW, LXSH)){
			//fputs("linha 460\n", IW);
			fprintf(IW, "%s(%.14E,%4d%s\n", LKEYW, VALUE, ICHPAR, LOPEN);
			XSHIFT = VALUE;
	} else if (!strcmp(LKEYW, LYSH)){
			//fputs("linha 464\n", IW);
			fprintf(IW, "%s(%.14E,%4d%s\n", LKEYW, VALUE, ICHPAR, LOPEN);
			YSHIFT = VALUE;
	} else if (!strcmp(LKEYW, LZSH)){
			//fputs("linha 468\n", IW)%s(%22.15f, %4d%s);
			fprintf(IW, "%s(%.14E,%4d%s\n", LKEYW, VALUE, ICHPAR, LOPEN);
			ZSHIFT = VALUE;
	} else {
		fputs(BLINE, IW);
		fputs("*** What do you mean?\n", IW);
		exit(0);	
	}
	goto L101;	
	
	//Quádrica expandida.
	
L102:; 
	QXX = KQ[1-1] / (XSCALE*XSCALE);
	QXY = 0.0e0;
	QXZ = 0.0e0;
	QYY = KQ[2-1] / (YSCALE*YSCALE);
	QYZ = 0.0e0;
	QZZ = KQ[3-1] / (ZSCALE*ZSCALE);
	QX = 0.0e0;
	QY = 0.0e0;
	QZ = KQ[4-1] / ZSCALE;
	Q0 = KQ[5-1];
	
	// Quádrica girada e deslocada.
	
///	printf("\n\nC++ ANTES KS: %d\nQXX: %f\nQXY: %f\nQXZ: %f\nQYY: %f\nQYZ: %f\nQZZ: %f\nQX: %f\nQY: %f\nQZ: %f\nQ0: %f\n", KS, QXX, QXY, QXZ, QYY, QYZ, QZZ, QX, QY, QZ, Q0 );
	
	
	rotshf2_(OMEGA, THETA, PHI, XSHIFT, YSHIFT, ZSHIFT, QXX, QXY, QXZ, QYY, QYZ, QZZ, QX, QY, QZ, Q0);
	
/*	TSTL = fmin(fmin(QXX,QXY),fmin(fmin(fmin(QXZ,QYY),QYZ),QZZ));
	TSTU = fmax(fmax(QXX,QXY),fmax(fmax(fmax(QXZ,QYY),QYZ),QZZ));
	
	if (fmax(fabs(TSTL),fabs(TSTU)) < 1.0e-30) 
		QSURF_.KPLANE[KS -1] = 1;	*/
	
	if (fmax(fabs(fmin(QXX, fmin(fmin(fmin(QXY, QXZ), QYY), QZZ))), 
			 fabs(fmax(QXX, fmax(fmax(fmax(QXY, QXZ), QYY), QZZ)))) < 1e-30){
		 QSURF_.KPLANE[KS -1] = 1;	
	}
			
   //  printf("\n\nC++ DEPOIS KS: %d KPLANE: %d\nQXX: %f20.8\nQXY: %f20.8\nQXZ: %f20.8\nQYY: %f20.8\nQYZ: %f20.8\nQZZ: %f20.8\nQX: %f20.8\nQY: %f20.8\nQZ: %f20.8\nQ0: %f20.8\n", KS, QSURF_.KPLANE[KS -1], QXX, QXY, QXZ, QYY, QYZ, QZZ, QX, QY, QZ, Q0 );
	
	
	QSURF_.AXX[KS - 1]=QXX;
	QSURF_.AXY[KS - 1]=QXY;
	QSURF_.AXZ[KS - 1]=QXZ;
	QSURF_.AYY[KS - 1]=QYY;
	QSURF_.AYZ[KS - 1]=QYZ;
	QSURF_.AZZ[KS - 1]=QZZ;
	QSURF_.AX[KS - 1]=QX;
	QSURF_.AY[KS - 1]=QY;
	QSURF_.AZ[KS - 1]=QZ;
	QSURF_.A0[KS - 1]=Q0;
	
	fputs(ZEROS, IW);
	strcpy(DEFS[KS - 1], DEF);
	
	goto L2;
	
	
	//Forma Implicita
	
L103:;

	QXX = 0.0e0;
	QXY = 0.0e0;
	QXZ = 0.0e0;
	QYY = 0.0e0;
	QYZ = 0.0e0;
	QZZ = 0.0e0;
	QX = 0.0e0;
	QY = 0.0e0;
	QZ = 0.0e0;
	Q0 =0.0e0;
	
L193:;

	fgets(BLINE, sizeof(BLINE), IR);
	if (comentario(BLINE)){
		   fputs(BLINE, IW);
		   goto L193;   
	}
	
	extrairString(LKEYW, BLINE, 0, 8); 
	extrairString(APOIO, BLINE, 9,22);
	VALUE = atof(APOIO);
	extrairString(APOIO, BLINE, 32, 4);
	ICHPAR = atoi(APOIO);
	extrairString(LANGLE, BLINE, 36, 8);
	
	if (!strcmp(LKEYW, LNUL)) goto L107;
	if (!strcmp(LKEYW, LONE)) goto L104;
	if (ICHPAR > 0){
		if (ICHPAR <= *NPINP){
			VALUE = PARINP[ICHPAR-1];
			ICHPAR = -ICHPAR; //Desativa a opção de alteração de parâmetro.
		} else{
			fprintf(IW, "%s(%f,%d,%s)\n", LKEYW, VALUE, ICHPAR, LANGLE);
			//fputs("linha 539", IW);
			fputs("*** NPINP is too small (check PARINP).", IW);
			exit(0);
		}
	}
	
	if (!strcmp(LKEYW, LAXX)){
		QXX = VALUE;
		//fputs("linha 547\n", IW);
		fprintf(IW, "%s(%.14E,%4d%s \n", LAXX, QXX, ICHPAR, LOPEN);
	} else if (!strcmp(LKEYW, LAXY)){
		QXY = VALUE;
		//fputs("linha 551\n", IW);
		fprintf(IW, "%s(%.14E,%4d%s\n", LAXY, QXY, ICHPAR, LOPEN);
	} else if (!strcmp(LKEYW, LAXZ)){
		QXZ = VALUE;
		fprintf(IW, "%s(%.14E,%4d%s\n", LAXZ, QXZ, ICHPAR, LOPEN);
		//fputs("linha 555\n", IW);
	} else if (!strcmp(LKEYW, LAYY)){
		QYY = VALUE;
		//fputs("linha 559\n", IW);
		fprintf(IW, "%s(%.14E,%4d%s\n", LAYY, QYY, ICHPAR, LOPEN);
	} else if (!strcmp(LKEYW, LAYZ)){
		QYZ = VALUE;
		fprintf(IW, "%s(%.14E,%4d%s\n", LAYZ, QYZ, ICHPAR, LOPEN);
		//fputs("linha 563\n", IW);
	} else if (!strcmp(LKEYW, LAZZ)){
		QZZ = VALUE;
		fprintf(IW, "%s(%.14E,%4d%s\n", LAZZ, QZZ, ICHPAR, LOPEN);
		//fputs("linha 567\n", IW);
	} else if (!strcmp(LKEYW, LAX)){
		QX = VALUE;
		fprintf(IW, "%s(%.14E,%4d%s\n", LAX, QX, ICHPAR, LOPEN);
		//fputs("linha 571\n", IW);
	} else if (!strcmp(LKEYW, LAY)){
		QY = VALUE;
		fprintf(IW, "%s(%.14E,%4d%s\n", LAY, QY, ICHPAR, LOPEN);
		//fputs("linha 575\n", IW);
	} else if (!strcmp(LKEYW, LAZ)){
		QZ = VALUE;
		fprintf(IW, "%s(%.14E,%4d%s\n", LAZ, QZ, ICHPAR, LOPEN);
		//fputs("linha 579\n", IW);
	} else if (!strcmp(LKEYW, LA0)){
		Q0 = VALUE;
		fprintf(IW, "%s(%.14E,%4d%s\n", LA0, Q0, ICHPAR, LOPEN);
		//fputs("linha 583\n", IW);
	} else {
		fputs(BLINE, IW);
		fputs("*** What do you mean?\n", IW);
		exit(0);
	}
	goto L193;
	
	
	//Parâmetros de transformação.
L104:;

	fputs(ZEROS, IW);
	OMEGA = 0.0e0;
	THETA = 0.0e0;
	PHI = 0.0e0;
	XSHIFT = 0.0e0;
	YSHIFT = 0.0e0;
	ZSHIFT = 0.0e0;
	
L105:;

	fgets(BLINE, sizeof(BLINE), IR);
if (comentario(BLINE)){

		   fputs(BLINE, IW);
		   goto L105;   
	}
	
	extrairString(LKEYW, BLINE, 0, 8); 
	extrairString(APOIO, BLINE, 9,22);
	VALUE = atof(APOIO);
	extrairString(APOIO, BLINE, 32, 4);
	ICHPAR = atoi(APOIO);
	extrairString(LANGLE, BLINE, 36, 8);
	
	if (!strcmp(LKEYW, LNUL)) 
		goto L106;
	if (ICHPAR > 0)
		if (ICHPAR <= *NPINP){
			VALUE = PARINP[ICHPAR-1];
			ICHPAR = -ICHPAR; //Desativa a opção de alteração de parâmetro.
	} else{
		fprintf(IW, "%s(%f,%d,%s)\n", LKEYW, VALUE, ICHPAR, LANGLE);
		fputs("*** NPINP is too small (check PARINP).\n", IW);
		exit(0);
	}
	
	if (!strcmp(LKEYW, LOME)){
		if (!strcmp(LANGLE, LRAD)){
			//fputs("linha 626\n", IW);
			fprintf(IW, "%s(%.14E,%4d%s\n", LKEYW, VALUE, ICHPAR, LRAD);
			OMEGA = VALUE;
		} else{
			//fputs("linha 630\n", IW);
			fprintf(IW, "%s(%.14E,%4d%s\n", LKEYW, VALUE, ICHPAR, LDEG);
			OMEGA = VALUE * PI / 180.0;
		}
	} else if (!strcmp(LKEYW, LTHE)){
		if (!strcmp(LANGLE, LRAD)){
			//fputs("linha 636\n", IW);
			fprintf(IW, "%s(%.14E,%4d%s\n", LKEYW, VALUE, ICHPAR, LRAD);
			THETA = VALUE;
		} else{
			//fputs("linha 640\n", IW);
			fprintf(IW, "%s(%.14E,%4d%s\n", LKEYW, VALUE, ICHPAR, LDEG);
			THETA = VALUE * PI / 180.0;
		}
	} else if (!strcmp(LKEYW, LPHI)){
		if (!strcmp(LANGLE, LRAD)){
		//	fputs("linha 646\n", IW);
			fprintf(IW, "%s(%.14E,%4d%s\n", LKEYW, VALUE, ICHPAR, LRAD);
			PHI = VALUE;
		} else{
			//fputs("linha 650\n", IW);
			fprintf(IW, "%s(%.14E,%4d%s\n", LKEYW, VALUE, ICHPAR, LDEG);
			PHI = VALUE * PI / 180.0;
		}
	} else if (!strcmp(LKEYW, LXSH)){
			   //	fputs("linha 655\n", IW);
			fprintf(IW, "%s(%.14E,%4d%s\n", LKEYW, VALUE, ICHPAR, LOPEN);
			XSHIFT = VALUE;
	} else if (!strcmp(LKEYW, LYSH)){
			//fputs("linha 659\n", IW);
			fprintf(IW, "%s(%.14E,%4d%s\n", LKEYW, VALUE, ICHPAR, LOPEN);
			YSHIFT = VALUE;
	} else if (!strcmp(LKEYW, LZSH)){
			//fputs("linha 663\n", IW);
			fprintf(IW, "%s(%.14E,%4d%s\n", LKEYW, VALUE, ICHPAR, LOPEN);
			ZSHIFT = VALUE;
	} else {
		fputs(BLINE, IW);
		fputs("*** What do you mean?\n", IW);
		exit(0);	
	}
	goto L105;	
	
	//Rotação e translação da superfície.
L106:;

//	printf("\n\nC++ ANTES KS: %d\nQXX: %f\nQXY: %f\nQXZ: %f\nQYY: %f\nQYZ: %f\nQZZ: %f\nQX: %f\nQY: %f\nQZ: %f\nQ0: %f\n", KS, QXX, QXY, QXZ, QYY, QYZ, QZZ, QX, QY, QZ, Q0 );
	
	
	rotshf2_(OMEGA, THETA, PHI, XSHIFT, YSHIFT, ZSHIFT, QXX, QXY, QXZ, QYY, QYZ, QZZ, QX, QY, QZ, Q0);

	
	 //   printf("\n\nC++ DEPOIS KS: %d\nQXX: %f20.8\nQXY: %f20.8\nQXZ: %f20.8\nQYY: %f20.8\nQYZ: %f20.8\nQZZ: %f20.8\nQX: %f20.8\nQY: %f20.8\nQZ: %f20.8\nQ0: %f20.8\n", KS, QXX, QXY, QXZ, QYY, QYZ, QZZ, QX, QY, QZ, Q0 );

L107:;

	TSTL = fmin(fmin(QXX,QXY),fmin(fmin(fmin(QXZ,QYY),QYZ),QZZ));
	TSTU = fmax(fmax(QXX,QXY),fmax(fmax(fmax(QXZ,QYY),QYZ),QZZ));
	
	if (fmax(fabs(TSTL),fabs(TSTU)) < 1.0e-30) 
		QSURF_.KPLANE[KS -1] = 1;
	  
  
   //  printf("\n\nC++ 107 KS: %d  KPLANE: %d\nQXX: %f20.8\nQXY: %f20.8\nQXZ: %f20.8\nQYY: %f20.8\nQYZ: %f20.8\nQZZ: %f20.8\nQX: %f20.8\nQY: %f20.8\nQZ: %f20.8\nQ0: %f20.8\n", KS, QSURF_.KPLANE[KS -1], QXX, QXY, QXZ, QYY, QYZ, QZZ, QX, QY, QZ, Q0 );
	
	QSURF_.AXX[KS - 1]=QXX;
	QSURF_.AXY[KS - 1]=QXY;
	QSURF_.AXZ[KS - 1]=QXZ;
	QSURF_.AYY[KS - 1]=QYY;
	QSURF_.AYZ[KS - 1]=QYZ;
	QSURF_.AZZ[KS - 1]=QZZ;
	QSURF_.AX[KS - 1]=QX;
	QSURF_.AY[KS - 1]=QY;
	QSURF_.AZ[KS - 1]=QZ;
	QSURF_.A0[KS - 1]=Q0;
	
	fputs(ZEROS, IW);
	strcpy(DEFS[KS - 1], DEF);
	
	goto L2;
	
	
	//Corpos
	
	
L200:;

	if ((BLINE[8] != '(' ) || (BLINE[13] != ')')){
		fputs(BLINE, IW);
		fputs("*** Incorrect label format.", IW);
		exit(0);
	}
 	extrairString(C4, BLINE, 9, 4);
 	if (KEEPL == 1){
 		strcat(strcpy(C5, C4), "0");
 		C5[5] = '\0';
	 }else 
 		strcat(strcpy(C5, C4), C1);
 	if (*QTREE_.NBODYS > 0){
 		for (int KB0 = 1; KB0 <= *QTREE_.NBODYS; KB0++){
 			extrairString(APOIO, ALIAB[KB0-1], 0, 5);
 			if (!strcmp(C5, APOIO)){
				fputs(BLINE, IW);
				fputs("*** Same label for two bodys (or modules).", IW);
				exit(0); 
			}
		}
	}
	*QTREE_.NBODYS = *QTREE_.NBODYS + 1;
	fprintf(IW, "%s(%4d)\n", LKEYW, *QTREE_.NBODYS);
	//fputs("linha 722\n", IW);
	if (*QTREE_.NBODYS > NB){
		fputs("*** The parameter NB must be increased.", IW);
		exit(0);
	}
	//escrito no def fputs("linha 727\n", IW);
	strcpy(ALIAB[*QTREE_.NBODYS -1], C5);
	
//	printf("linha 1705 ALIAB: %s --- C5: %s",ALIAB[*QTREE_.NBODYS -1], C5 );
	
	
L295:;

	fgets(BLINE, sizeof(BLINE), IR);
	if (comentario(BLINE)){
		   fputs(BLINE, IW);
		   goto L295;   
	}
	
////	char tmpIMAT[5];
//	strncpyy(tmpIMAT, BLINE+9, 4);
	extrairString(APOIO, BLINE, 9, 4);
	IMAT = atoi(APOIO);
	extrairString(LKEYW, BLINE, 0, 8); 
//	IMAT = atoi(strncpyy(APOIO, BLINE+9, 4));
//	printf("Linha 1704 LKEYW:%s IMAT: %d\n", LKEYW, IMAT);
	
	if (strcmp(LKEYW, LMAT)){
		fputs(BLINE, IW);
		fputs("*** Incorrect material definition line.", IW);
		exit(0);
	}
	fprintf(IW, "%s(%4d)\n", LKEYW, IMAT);
	//fputs("linha 751\n", IW);
	if (IMAT < 0) 
		IMAT = 0;
//	printf("linha 1710\n");
	PENGEOM_mod_.MATER[*QTREE_.NBODYS-1] = IMAT;
//	printf("linha 1710\n");
//	printf("Linha 1717 IMAT: %d    NMATG: %d\n", IMAT, *NMATG);

	*NMATG = max(*NMATG, IMAT);
	
	
	
	//	printf("linha 1709\n");
	
L201:;
	
	fgets(BLINE, sizeof(BLINE), IR);
	if (comentario(BLINE)){
		   fputs(BLINE, IW);
		   goto L201;   
	}
	
	extrairString(LKEYW, BLINE, 0, 8); 
	if (!strcmp(LKEYW, LNUL)) 
		goto L208;
	
	extrairString(C4, BLINE, 9, 4);
	if (KEEPL == 1){
 		strcat(strcpy(C5, C4), "0");
 		C5[5] = '\0';
	}else 
 		strcat(strcpy(C5, C4), C1);
 	
 	//	printf("linha 1730\n");
 	
 	if ((!strcmp(LKEYW, LSUR)) || (!strcmp(LKEYW, LSUA))){
 		extrairString(APOIO, BLINE, 30, 2);
		INDS = atoi(APOIO);
	//	printf("INDS CORPO: %d\n", INDS );
 		
 		//	printf("linha 1735 INDS: %d\n", INDS);
 		//Superficie
 		for (int KS0=1; KS0<= *QSURF_.NSURF; KS0++){
 		//	printf("linha 1738: C5: %s      ALIAS: %s \n",C5, strncpyy( APOIO, ALIAS[KS0-1], 5) );
 			 extrairString(APOIO, ALIAS[KS0-1], 0, 5 );
			 if (!strcmp(C5, APOIO)){
				 KS = KS0;
			//	 printf("linha 1740");
				 goto L202;
			 } 
		 }
		fputs(BLINE, IW);
		fputs("*** Undefined surface label.\n", IW);
		exit(0);
L202:;
		fprintf(IW, "%s(%4d), SIDE POINTER=(%2d)\n", LKEYW, KS, INDS);
	//	fputs("linha 780\n", IW);
		KST = QTREE_.KSURF[NXG-1][*QTREE_.NBODYS - 1];
		if (KST > 0){
			KST0 = KST;
			for (int K = 1; K<= KST0; K++){
				if (KS == QTREE_.KSURF[K-1][*QTREE_.NBODYS -1]){
					if (QTREE_.KFLAG[K-1][*QTREE_.NBODYS -1] < 3){
						fputs("*** The last limiting surface has been defined twice.\n", IW);
						exit(0);
					} else{
						KST = K; //KS limita um corpo limitador.
						goto L203;
					}
				}
			}
		}
		KST = KST+1;
		if (KST >= NXG){
			fputs("*** The number of limiting surfaces is too large.", IW);
			exit(0);
		}
		QTREE_.KSURF[NXG -1][*QTREE_.NBODYS -1] = KST;
		QTREE_.KSURF[KST -1][*QTREE_.NBODYS -1] = KS;
L203:;
		if (INDS == -1){
			QTREE_.KFLAG[KST-1][*QTREE_.NBODYS-1] = 1;
		}else if (INDS == 1){
			QTREE_.KFLAG[KST-1][*QTREE_.NBODYS-1] = 2;	
		} else{
			fputs(BLINE, IW);
			fputs("*** Check side pointer value.", IW);
			exit(0);	
		}	
	 } else if(!strcmp(LKEYW, LBOD)){
	 	//corpo
	 	if (*QTREE_.NBODYS > 1){
			 for (int KB0 = 1; KB0 <= *QTREE_.NBODYS - 1; KB0++){
			 	extrairString(APOIO, ALIAB[KB0 - 1], 0, 5 );
			 	if (!strcmp(C5, APOIO)){
			 		KB = KB0;
			 		goto L204;											 
				 }
			 }
			fputs(BLINE, IW);
			fputs("*** Undefined body label.", IW);
			exit(0);	
		 }
L204:;
		fprintf(IW, "%s(%d)\n", LKEYW, KB);
		//fputs("LInha 831\n", IW);
		if (QBODY_.KBOMO[KB - 1] != 0){
			fputs("*** This body is a module.", IW);
			exit(0);
		}
		
		KN1 = QTREE_.KSURF[NXG-1][KB-1];
	
		for (int KS1 = 1; KS1 <= KN1; KS1++){
			KSURF1 = QTREE_.KSURF[KS1-1][KB-1];
			KN2 = QTREE_.KSURF[NXG-1][*QTREE_.NBODYS-1];
			for (int KS2 = 1; KS2 <= KN2; KS2++){
				KSURF2 = QTREE_.KSURF[KS2-1][*QTREE_.NBODYS-1];
				if (KSURF2 == KSURF1) 
					goto L205;	
			}
			KST = KN2 + 1;
			if (KST >= NXG){
				fputs(BLINE, IW);
				fputs("*** The number of limiting surfaces is too large.", IW);
				exit(0);		
			}
			QTREE_.KSURF[NXG-1][*QTREE_.NBODYS-1] = KST;
			QTREE_.KSURF[KST-1][*QTREE_.NBODYS-1] = KSURF1;
			QTREE_.KFLAG[KST-1][*QTREE_.NBODYS-1] = 3;
L205:;	
		}
		QBODY_.KBODY[NXG-1][*QTREE_.NBODYS-1] = QBODY_.KBODY[NXG-1][*QTREE_.NBODYS-1] + 1;
		QBODY_.KBODY[QBODY_.KBODY[NXG-1][*QTREE_.NBODYS-1]-1][*QTREE_.NBODYS-1] = KB;
	 } else if(!strcmp(LKEYW, LMOD)){
	 	//Modulo.
	 	if (*QTREE_.NBODYS > 1){
			 for (int KB0 = 1; KB0 <= *QTREE_.NBODYS - 1; KB0++){
			 	extrairString(APOIO, ALIAB[KB0-1], 0, 5);
				 if (!strcmp(C5, APOIO)){
					 KB = KB0;
					 goto L206;
				 }
			 }
			 fputs(BLINE, IW);
		     fputs("*** Undefined body label.\n", IW);
		 	 exit(0);		 
		}
L206:;
		fprintf(IW, "%s(%d) linha 879\n", LKEYW, KB);
	//	fputs("LInha 872\n", IW);
		if (QBODY_.KBOMO[KB -1] != 1){
			fputs("*** This module is a body.\n", IW);
			exit(0);
		}
		
		KN1 = QTREE_.KSURF[NXG-1][KB-1];
		for (int KS1 = 1; KS1 <= KN1; KS1++){
			KSURF1 = QTREE_.KSURF[KS1-1][KB-1];
			if (QTREE_.KFLAG[KS1-1][KB-1] > 2)
				goto L207;
			KN2 = QTREE_.KSURF[NXG-1][*QTREE_.NBODYS-1];
			for (int KS2 = 1; KS2 <= KN2; KS2++){
				KSURF2 = QTREE_.KSURF[KS2-1][*QTREE_.NBODYS-1];
				if (KSURF2 == KSURF1)
					goto L207;
			}
			KST = KN2+1;
			if (KST >= NXG){
				fputs(BLINE, IW);
    			fputs("*** The number of limiting surfaces is too large.\n'", IW);
				exit(0);		
			}
			QTREE_.KSURF[NXG-1][*QTREE_.NBODYS-1] = KST;
			QTREE_.KSURF[KST-1][*QTREE_.NBODYS-1] = KSURF1;
			QTREE_.KFLAG[KST-1][*QTREE_.NBODYS-1] = 3;
L207:;	
		}
		QBODY_.KBODY[NXG-1][*QTREE_.NBODYS-1] = QBODY_.KBODY[NXG-1][*QTREE_.NBODYS-1] + 1;
		QBODY_.KBODY[QBODY_.KBODY[NXG-1][*QTREE_.NBODYS-1]-1][*QTREE_.NBODYS-1] = KB;	 
	 } else{
	 	fputs(BLINE, IW);
		fputs("*** What do you mean?\n", IW);	
		exit(0);
	 }
	 goto L201;
L208:;

 	fputs(ZEROS, IW);
//	printf("linha 1884: DEF: %s\n\n ", DEF);
	strcpy(DEFB[*QTREE_.NBODYS -1], DEF);
	goto L2;
	
	//Modulos

L300:;

	if ((BLINE[8] != '(' ) || (BLINE[13] != ')')){
		fputs(BLINE, IW);
		fputs("*** Incorrect label format.\n", IW);
		exit(0);
	}
 	extrairString(C4, BLINE, 9, 4);
 	if (KEEPL == 1){
 		strcat(strcpy(C5, C4), "0");
 		C5[5] = '\0';
 		
	 }else 
 		strcat(strcpy(C5, C4), C1);
 	if (*QTREE_.NBODYS > 0){
 		for (int KB0 = 1; KB0 <= *QTREE_.NBODYS; KB0++){
 			extrairString(APOIO, ALIAB[KB0-1], 0, 5);
 			if (!strcmp(C5, APOIO)){
				fputs(BLINE, IW);
				fputs("*** Same label for two bodys (or modules).\n", IW);
				exit(0); 
			}
		}
	}
	*QTREE_.NBODYS = *QTREE_.NBODYS + 1;
	fprintf(IW, "%s(%4d)\n", LKEYW, *QTREE_.NBODYS);
	
	//fputs("linha 939\n", IW);
	if (*QTREE_.NBODYS > NB){
		fputs("*** The parameter NB must be increased.\n", IW);
		exit(0);
	}
	//escrito no DEF fputs("linha 944\n", IW);
	strcpy(ALIAB[*QTREE_.NBODYS-1], C5);
	
	
L391:;

	fgets(BLINE, sizeof(BLINE), IR);
	

 	 if (comentario(BLINE)){

		   fputs(BLINE, IW);
		   goto L391;   
	}
	
	extrairString(LKEYW, BLINE, 0, 8); 
//	strncpyy(tmpIMAT, BLINE+9, 4);
	extrairString(APOIO, BLINE, 9, 4);
	IMAT = atoi(APOIO);
//	IMAT = atoi(strncpyy(APOIO, BLINE+9, 4));
	
	
	fprintf(IW, "%s(%4d)\n", LKEYW, IMAT);
	if (strcmp(LKEYW, LMAT)){
		fputs("*** Incorrect material definition line.", IW);
		exit(0);
	}
	
	
	if (IMAT < 0) 
		IMAT = 0;
	PENGEOM_mod_.MATER[*QTREE_.NBODYS -1] = IMAT;
	*NMATG = max(*NMATG, IMAT);
	
	QTREE_.KDGHT[NXG-1][*QTREE_.NBODYS-1] = 1;
	QTREE_.KDGHT[1-1][*QTREE_.NBODYS-1] = *QTREE_.NBODYS;
	
	
	//Limitando superfícies e componentes.

L301:;
	fgets(BLINE, sizeof(BLINE), IR);
	if (comentario(BLINE)){
		fputs(BLINE, IW);
		goto L301;   
	}
	
	extrairString(LKEYW, BLINE, 0, 8); 
	if ((!strcmp(LKEYW, LNUL)) || (!strcmp(LKEYW, LONE))){
		KDT = QTREE_.KDGHT[NXG-1][*QTREE_.NBODYS-1];
		
		//Classifique as filhas em ordem crescente.
		if (KDT > 1){
			for (int KI = 1; KI <= KDT-1; KI++){
				KBMIN = QTREE_.KDGHT[KI-1][*QTREE_.NBODYS-1];
				KMIN = KI;
				for (int KJ = KI+1; KJ <= KDT; KJ++){
					if (QTREE_.KDGHT[KJ-1][*QTREE_.NBODYS-1] < KBMIN){
						KBMIN = QTREE_.KDGHT[KJ-1][*QTREE_.NBODYS-1];
						KMIN = KJ;
					}
				}
				if (KMIN != KI){
					KSAVE = QTREE_.KDGHT[KI-1][*QTREE_.NBODYS-1];
					QTREE_.KDGHT[KI-1][*QTREE_.NBODYS-1] = QTREE_.KDGHT[KMIN-1][*QTREE_.NBODYS-1];
					QTREE_.KDGHT[KMIN-1][*QTREE_.NBODYS-1] = KSAVE;
				}
			}
		}
		
		if (!strcmp(LKEYW, LONE)) 
			goto L308;
		if (!strcmp(LKEYW, LNUL))
			goto L312;
	} 
	
	//Limitando Superficies
	if ((!strcmp(LKEYW, LSUR)) || (!strcmp(LKEYW, LSUA))){
		extrairString(C4, BLINE, 9, 4);
		extrairString(APOIO, BLINE, 30, 2);
		INDS = atoi(APOIO);
	//	printf("INDS: %d\n", INDS );
		
		if (KEEPL == 1){
 		   	strcat(strcpy(C5, C4), "0");
 		   	C5[5] = '\0'; 
		}
 		else 
 		   	strcat(strcpy(C5, C4), C1);
 	
 		
 		//Superficie
 		for (int KS0=1; KS0<= *QSURF_.NSURF; KS0++){
 			extrairString(APOIO, ALIAS[KS0-1], 0, 5);
			 if (!strcmp(C5, APOIO)){
				 KS = KS0;
				 goto L302;
			 } 
		 }
		 fprintf(IW, "%s(%4d), SIDE POINTER=(%2d)\n", LKEYW, KS, INDS);
		//fputs("linha 1020\n", IW);
		fputs("*** Undefined surface label.\n", IW);
		exit(0);
L302:;
		//fputs("linha 1025\n", IW);
		 fprintf(IW, "%s(%4d), SIDE POINTER=(%2d)\n", LKEYW, KS, INDS);
		
		KST = QTREE_.KSURF[NXG-1][*QTREE_.NBODYS - 1];
		if (KST > 0){
			KST0 = KST;
			for (int K = 1; K<= KST0; K++){
				if (KS == QTREE_.KSURF[K-1][*QTREE_.NBODYS -1]){
					if (QTREE_.KFLAG[K-1][*QTREE_.NBODYS -1] < 3){
						fputs("*** The last limiting surface has been defined twice.\n", IW);
						exit(0);
					} else{
						KST = K; //KS limita um corpo limitador.
						goto L303;
					}
				}
			}
		}
		KST = KST+1;
		if (KST >= NXG){
			fputs("*** The number of limiting surfaces is too large.\n", IW);
			exit(0);
		}
		QTREE_.KSURF[NXG -1][*QTREE_.NBODYS -1] = KST;
		QTREE_.KSURF[KST -1][*QTREE_.NBODYS -1] = KS;
L303:;
		if (INDS == -1){
			QTREE_.KFLAG[KST-1][*QTREE_.NBODYS-1] = 1;
		}else if (INDS == 1){
			QTREE_.KFLAG[KST-1][*QTREE_.NBODYS-1] = 2;	
		} else{
			fputs("*** Check side pointer value.\n", IW);
			exit(0);	
		}
		
		//Corpos
	}else if (!strcmp(LKEYW, LBOD)){
		extrairString(C4, BLINE, 9, 4);
		if (KEEPL == 1){
 		   	strcat(strcpy(C5, C4), "0");
 		   	C5[5] = '\0';
		}else 
 		   	strcat(strcpy(C5, C4), C1);
			 
		if (*QTREE_.NBODYS > 1){ //Verifique se KB está na lista.
			 for (int KB0 = 1; KB0 <= *QTREE_.NBODYS - 1; KB0++){
			 	extrairString(APOIO, ALIAB[KB0 - 1], 0, 5 );
			 	if (!strcmp(C5, APOIO)){
			 		KB = KB0;
			 		goto L304;											 
				 }
			 }
			//fputs("linha 1080\n", IW);
			fprintf(IW, "%s(%s)n", LKEYW, C4);
			fputs("*** Undefined body label.\n", IW);
			exit(0);	
		 }
L304:;
		//fputs("LInha 1085\n", IW);
		fprintf(IW, "%s(%4d)\n", LKEYW, KB);
		if (QBODY_.KBOMO[KB - 1] != 0){
			fputs("*** This body is a module.\n", IW);
			exit(0);
		}
		
		if ((QTREE_.KMOTH[KB - 1] > 0) && (QTREE_.KMOTH[KB - 1] != *QTREE_.NBODYS)){
			fputs("You are trying to assign two mothers to the last body.\n", IW);
			exit(0);
		}
		
		QTREE_.KMOTH[KB - 1] = *QTREE_.NBODYS;
		KDT = QTREE_.KDGHT[NXG-1][*QTREE_.NBODYS-1] + 1;
	   	QTREE_.KDGHT[NXG-1][*QTREE_.NBODYS-1] = KDT;
	   	QTREE_.KDGHT[KDT-1][*QTREE_.NBODYS-1] = KB;
	   	
	   	for (int K2 = 1; K2 <= QBODY_.KBODY[NXG-1][KB-1]; K2++){
			if (QTREE_.KMOTH[QBODY_.KBODY[K2-1][KB-1]-1] == 0)
				QTREE_.KMOTH[QBODY_.KBODY[K2-1][KB-1]-1] = *QTREE_.NBODYS;	   
		}
		
		//Atribua atributos genealógicos às irmãs do corpo KB.
		
		K2M = QBODY_.KBODY[NXG-1][KB-1];
		for (int K2 = 1; K2 <= K2M; K2++){
			KBD = QBODY_.KBODY[K2-1][KB-1];
			if (QTREE_.KMOTH[KBD - 1] == 0)
				QTREE_.KMOTH[KBD - 1] = *QTREE_.NBODYS;
			IDGHT = 0;
			for (int K = 1; K <= KDT; K++){
				if (QTREE_.KDGHT[K-1][*QTREE_.NBODYS-1] == KBD)
					IDGHT = K;	
			}
			if (IDGHT == 0){
				QTREE_.KDGHT[NXG-1][*QTREE_.NBODYS-1] = QTREE_.KDGHT[NXG-1][*QTREE_.NBODYS-1] + 1;
				QTREE_.KDGHT[QTREE_.KDGHT[NXG-1][*QTREE_.NBODYS-1]-1][*QTREE_.NBODYS-1] = KB;	
			}
		}
		
		// superficies dos corpos irmãos
		KN1 = QTREE_.KSURF[NXG-1][KB-1];
		for (int KS1 = 1; KS1 <= KN1; KS1++){
			if (QTREE_.KFLAG[KS1-1][KB-1] > 3)
				goto L317;
			KSURF1 = QTREE_.KSURF[KS1-1][KB-1];
			KN2 = QTREE_.KSURF[NXG-1][*QTREE_.NBODYS-1];
			for (int KS2 = 1; KS2 <= KN2; KS2++){
				KSURF2 = QTREE_.KSURF[KS2-1][*QTREE_.NBODYS-1];
				if (KSURF2 == KSURF1)
					goto L317;
			}
			KN2 = KN2 + 1; 
			if (KN2 >= NXG){
				fputs("*** The number of limiting surfaces is too large.\n", IW);
				exit(0);
			}
			QTREE_.KSURF[NXG-1][*QTREE_.NBODYS-1] = KN2;
			QTREE_.KSURF[KN2-1][*QTREE_.NBODYS-1] = KSURF1;
			QTREE_.KFLAG[KN2-1][*QTREE_.NBODYS-1] = 4;
L317:;
		}
	
		//Modulo
		   
	} else if (!strcmp(LKEYW, LMOD)){ //Verifique se KB está na lista.
		extrairString(C4, BLINE, 9, 4);
		if (KEEPL == 1){
 		   	strcat(strcpy(C5, C4), "0");
 		   	C5[5] = '\0';
 		   	
		}else 
 			strcat(strcpy(C5, C4), C1);
 		if (*QTREE_.NBODYS > 0){
			 for (int KB0 = 1; KB0 <= *QTREE_.NBODYS - 1; KB0++){
			 	extrairString(APOIO, ALIAB[KB0-1], 0, 5);
			 //	printf("C5: %s   ALIAB: %s\n", C5, strncpyy(APOIO, ALIAB[KB0-1], 5));
				 if (!strcmp(C5, APOIO)){
				 //	printf("\n\n\n\nKB = KBO\n\n\n\n");
					 KB = KB0;
					 goto L305;
				 }
			 }
			 //fputs("linha 1155", IW);
			 fprintf(IW, "%s(%s) linha 1062\n", LKEYW, C4);
		     fputs("*** Undefined body label.\n", IW);
		 	 exit(0);		 
		}
L305:;
		//fputs("LInha 1160\n", IW);
		fprintf(IW, "%s(%4d)\n", LKEYW, KB);
		if (QBODY_.KBOMO[KB-1] != 1){
			fputs("*** This module is a body.\n", IW);
			exit(0);
		}
		if ((QTREE_.KMOTH[KB - 1] > 0) && (QTREE_.KMOTH[KB - 1] != *QTREE_.NBODYS)){
			fputs("You are trying to assign two mothers to the last body.\n", IW);
			exit(0);
		}
		
		QTREE_.KMOTH[KB - 1] = *QTREE_.NBODYS;
		KDT = QTREE_.KDGHT[NXG-1][*QTREE_.NBODYS-1] + 1;
	   	QTREE_.KDGHT[NXG-1][*QTREE_.NBODYS-1] = KDT;
	   	QTREE_.KDGHT[KDT-1][*QTREE_.NBODYS-1] = KB;
	   	
	   	//Superficies Limitantes
		KN1 = QTREE_.KSURF[NXG-1][KB-1];
		for (int KS1 = 1; KS1 <= KN1; KS1++){
			if (QTREE_.KFLAG[KS1-1][KB-1] > 2)
				goto L307;
			KSURF1 = QTREE_.KSURF[KS1-1][KB-1];
			KN2 = QTREE_.KSURF[NXG-1][*QTREE_.NBODYS-1];
			for (int KS2 = 1; KS2 <= KN2; KS2++){
				KSURF2 = QTREE_.KSURF[KS2-1][*QTREE_.NBODYS-1];
				if (KSURF2 == KSURF1)
					goto L307;
			}
			KN2 = KN2 + 1;
			if (KN2 >= NXG){
				fputs("*** The number of limiting surfaces is too large.\n", IW);
				exit(0);
			}
			QTREE_.KSURF[NXG-1][*QTREE_.NBODYS-1] = KN2;
			QTREE_.KSURF[KN2-1][*QTREE_.NBODYS-1] = KSURF1;
			QTREE_.KFLAG[KN2-1][*QTREE_.NBODYS-1] = 4;
		//	printf("\n\n\n MODLUE KFLAG 4 \n\n\n");			
L307:;  		
		}
	} else{
		fputs(BLINE, IW);
		fputs("*** What do you mean?", IW);
		exit(0);	
	}
	
	goto L301;
	
	//Parametros de Transformação
	
L308:;

	fputs(ZEROS, IW);
	OMEGA=0.0e0;
	THETA=0.0e0;
	PHI=0.0e0;
	XSHIFT=0.0e0;
	YSHIFT=0.0e0;
	ZSHIFT=0.0e0;
	
L309:;

	fgets(BLINE, sizeof(BLINE), IR);
	if (comentario(BLINE)){

		   fputs(BLINE, IW);
		   goto L309;   
	}
	
	extrairString(LKEYW, BLINE, 0, 8); 
	extrairString(APOIO, BLINE, 9,22);
	VALUE = atof(APOIO);
	extrairString(APOIO, BLINE, 32, 4);
	ICHPAR = atoi(APOIO);
	extrairString(LANGLE, BLINE, 36, 8);
	
	if (!strcmp(LKEYW, LNUL))
		goto L310;
	if (ICHPAR > 0){
		if (ICHPAR <= *NPINP){
			VALUE = PARINP[ICHPAR-1];
			ICHPAR = -ICHPAR; //Desativa a opção de alteração de parâmetro.
		} else{
			//fputs("linha 1229", IW);
			fprintf(IW, "%s(%f,%d%s\n", LKEYW, VALUE, ICHPAR, LANGLE);
			fputs("*** NPINP is too small (check PARINP).\n", IW );
			exit(0);
		}
	}
	if (!strcmp(LKEYW, LOME)){
		if (!strcmp(LANGLE, LRAD)){
			//fputs("linha 1237\n", IW);
			fprintf(IW, "%s(%f,%d%s\n", LKEYW, VALUE, ICHPAR, LRAD);
			OMEGA = VALUE;
		} else{
			//fputs("linha 1241\n", IW);
			fprintf(IW, "%s(%f,%d%s\n", LKEYW, VALUE, ICHPAR, LDEG);
			OMEGA = VALUE * PI / 180.0;
		}
	} else if (!strcmp(LKEYW, LTHE)){
		if (!strcmp(LANGLE, LRAD)){
			//fputs("linha 1247\n", IW);
			fprintf(IW, "%s(%f,%d%s\n", LKEYW, VALUE, ICHPAR, LRAD);
			THETA = VALUE;
		} else{
			//fputs("linha 1251\n", IW);
			fprintf(IW, "%s(%f,%d%s\n", LKEYW, VALUE, ICHPAR, LDEG);
			THETA = VALUE * PI / 180.0;
		}
	} else if (!strcmp(LKEYW, LPHI)){
		if (!strcmp(LANGLE, LRAD)){
			//fputs("linha 1257\n", IW);
			fprintf(IW, "%s(%f,%d%s\n", LKEYW, VALUE, ICHPAR, LRAD);
			PHI = VALUE;
		} else{
			//fputs("linha 1261\n", IW);
			fprintf(IW, "%s(%f,%d%s\n", LKEYW, VALUE, ICHPAR, LDEG);
			PHI = VALUE * PI / 180.0;
		}
	} else if (!strcmp(LKEYW, LXSH)){
			//fputs("linha 1266\n", IW);
			fprintf(IW, "%s(%f,%d%s\n", LKEYW, VALUE, ICHPAR, LOPEN);
			XSHIFT = VALUE;
	} else if (!strcmp(LKEYW, LYSH)){
			//fputs("linha 1270\n", IW);
			fprintf(IW, "%s(%f,%d%s\n", LKEYW, VALUE, ICHPAR, LOPEN);
			YSHIFT = VALUE;
	} else if (!strcmp(LKEYW, LZSH)){
			//fputs("linha 1274\n", IW);
			fprintf(IW, "%s(%f,%d%s\n", LKEYW, VALUE, ICHPAR, LOPEN);
			ZSHIFT = VALUE;
	} else {
		fputs(BLINE, IW);
		fputs("*** What do you mean?\n", IW);
		exit(0);	
	}
	goto L309;	
	
	
	//Rotação e translação do conteúdo do módulo (suas superfícies).

L310:;
	for (int KS = 1; KS <= *QSURF_.NSURF; KS++){
		KM[KS-1] = 0;	
	}
	
	for (int KB = 1; KB <= *QTREE_.NBODYS; KB++){
		KBMOTH = KB; //Precisamos transformar todos os descendentes.
		
L311:;
		if (KBMOTH == *QTREE_.NBODYS){
			KNS = QTREE_.KSURF[NXG-1][KB-1];
			for (int KSS = 1;  KSS <= KNS; KSS++){
				KS = QTREE_.KSURF[KSS-1][KB-1];
				
				if ((QTREE_.KFLAG[KSS-1][KB-1] < 5) && (KM[KS-1] == 0) && KSFF[KS-1] == 0){
				//	printf("linha 2487 KS:%d\n ", KS);
					QXX = QSURF_.AXX[KS - 1];
					QXY = QSURF_.AXY[KS - 1];
					QXZ = QSURF_.AXZ[KS - 1];
					QYY = QSURF_.AYY[KS - 1];
					QYZ = QSURF_.AYZ[KS - 1];
					QZZ = QSURF_.AZZ[KS - 1];
					QX = QSURF_.AX[KS - 1];
					QY = QSURF_.AY[KS - 1];
					QZ = QSURF_.AZ[KS - 1];
					Q0 = QSURF_.A0[KS - 1];
					
					
   	   	   	//   printf("\n\nC++ ANTES KS: %d  KPLANE: %d\nQXX: %f20.8\nQXY: %f20.8\nQXZ: %f20.8\nQYY: %f20.8\nQYZ: %f20.8\nQZZ: %f20.8\nQX: %f20.8\nQY: %f20.8\nQZ: %f20.8\nQ0: %f20.8\n", KS, QSURF_.KPLANE[KS -1], QXX, QXY, QXZ, QYY, QYZ, QZZ, QX, QY, QZ, Q0 );
	
					
					rotshf2_(OMEGA, THETA, PHI, XSHIFT, YSHIFT, ZSHIFT, QXX, QXY, QXZ, QYY, QYZ, QZZ, QX, QY, QZ, Q0);
			
					
    // printf("\n\nC++ DEPOIS KS: %d  KPLANE: %d\nQXX: %f20.8\nQXY: %f20.8\nQXZ: %f20.8\nQYY: %f20.8\nQYZ: %f20.8\nQZZ: %f20.8\nQX: %f20.8\nQY: %f20.8\nQZ: %f20.8\nQ0: %f20.8\n", KS, QSURF_.KPLANE[KS -1], QXX, QXY, QXZ, QYY, QYZ, QZZ, QX, QY, QZ, Q0 );
			
					QSURF_.AXX[KS - 1] = QXX;
					QSURF_.AXY[KS - 1] = QXY;
					QSURF_.AXZ[KS - 1] = QXZ;
					QSURF_.AYY[KS - 1] = QYY;
					QSURF_.AYZ[KS - 1] = QYZ;
					QSURF_.AZZ[KS - 1] = QZZ;
					QSURF_.AX[KS - 1] = QX;
					QSURF_.AY[KS - 1] = QY;
					QSURF_.AZ[KS - 1] = QZ;
					QSURF_.A0[KS - 1] = Q0;
					KM[KS-1] = 1;

				}
			}
		} else{
			KBMOTH = QTREE_.KMOTH[KBMOTH-1];
			if (KBMOTH > 0)
				goto L311;	
		}	
	}
	
L312:;
	fputs(ZEROS, IW);
	
	QBODY_.KBOMO[*QTREE_.NBODYS -1] = 1;
	strcpy(DEFB[*QTREE_.NBODYS -1], DEF);
	goto L2;
	
	//Clonando um Modulo.
	
L400:;
	if ((BLINE[8] != '(' ) || (BLINE[13] != ')')){
		fputs(BLINE, IW);
		fputs("*** Incorrect label format.", IW);
		exit(0);
	}
	fputs("linha 1345", IW);
	fputs("linha 1346", IW);
	fputs("C ************  Cloned module:", IW);
	fputs("linha 1348", IW);
	ICLONE = ICLONE+1;
 	extrairString(C4, BLINE, 9, 4);
 	if (KEEPL == 1)
 		strcat(strcpy(C5, C4), "0");
 	else 
 		strcat(strcpy(C5, C4), C1);
 	
 	strcpy(C5C, C5);
    if (comentario(BLINE)){
		fputs(BLINE, IW);
	}else {
		fputs("linha 1364", IW);
	}
	
	if (*QTREE_.NBODYS > 0){
		for (int KB0=1; KB0 <= *QTREE_.NBODYS; KB0++){
			extrairString(APOIO, ALIAB[KB0-1], 0, 5);
			if (!strcmp(C5, APOIO)){
				fputs(BLINE, IW);
			   	fputs("*** Same label for two bodies or modules..", IW );
				exit(0);
			}
		}
	}
	
L401:;

	fgets(BLINE, sizeof(BLINE), IR);
	if (comentario(BLINE)){
		   fputs(BLINE, IW);
		   goto L401;   
	}
	extrairString(LKEYW, BLINE, 0, 8);
	//strcpy(LKEYW, extrairString(BLINE, 0, 7)); 
	extrairString(C4, BLINE, 9, 4);
 	if (KEEPL == 1)
 		strcat(strcpy(C5, C4), "0");
 	else 
 		strcat(strcpy(C5, C4), C1);
 	
 	fputs("linha 1393", IW);
 	
 	if(!strcmp(LKEYW, LMOD)){
 		fputs(BLINE, IW);
 		fputs("*** The cloned object must be a module.", IW);
 		exit(0);
	 }
	 
	 if (*QTREE_.NBODYS == 0){
	 	fputs(BLINE, IW);
 		fputs("*** This module is not defined.", IW);
 		exit(0);
	 }
	 
	for (int KB0 = 1; KB0 <= *QTREE_.NBODYS; KB0++){
		extrairString(APOIO, ALIAB[KB0-1], 0, 5);
		if (!strcmp(C5, APOIO)){
		 	 KORIG = KB0;
			if (QBODY_.KBOMO[KORIG-1] <= 1){
				fputs("*** The cloned object must be a module.", IW);
				fputs("*** The selected object is a body.", IW);
				exit(0);
			}
			goto L402;
		}	
	}
	fputs("*** The label does not correspond to a module.", IW);
	exit(0);
	 
L402:;

	fgets(BLINE, sizeof(BLINE), IR);
	if (comentario(BLINE)){
		   fputs(BLINE, IW);
		   goto L402;   
	}
	extrairString(LKEYW, BLINE, 0, 8); 
	if ((!strcmp(LKEYW, LONE)) || (!strcmp(LKEYW, LNUL))) {
		goto L403;
	} else{
		fputs(BLINE, IW);
		fputs("*** What do you mean?", IW);
		exit(0);	
	}
	
	//Parâmetros de transformação.

L403:;
	OMEGA = 0.0;
	THETA = 0.0;
	PHI = 0.0;
	XSHIFT = 0.0;
	YSHIFT = 0.0;
	ZSHIFT = 0.0;
	
	if (!strcmp(LKEYW, LNUL))
		goto L405;

L404:;

	fgets(BLINE, sizeof(BLINE), IR);
	if (comentario(BLINE)){

		   fputs(BLINE, IW);
		   goto L404;   
	}
	
	extrairString(LKEYW, BLINE, 0, 8); 
	extrairString(APOIO, BLINE, 9,22);
	VALUE = atof(APOIO);
	extrairString(APOIO, BLINE, 32, 4);
	ICHPAR = atoi(APOIO);
	extrairString(LANGLE, BLINE, 36, 8);
	
	if (!strcmp(LKEYW, LNUL))
		goto L405;
	if (ICHPAR > 0){
		if (ICHPAR <= *NPINP){
			VALUE = PARINP[ICHPAR-1];
			ICHPAR = -ICHPAR; //Desativa a opção de alteração de parâmetro.
		} else{
			fputs("linha 1462", IW);
			fputs("*** NPINP is too small (check PARINP).", IW );
			exit(0);
		}
	}
	
	if (!strcmp(LKEYW, LOME)){
		if (!strcmp(LANGLE, LRAD)){
			fputs("linha 1470", IW);
			OMEGA = VALUE;
		} else{
			fputs("linha 1474", IW);
			OMEGA = VALUE * PI / 180.0;
		}
	} else if (!strcmp(LKEYW, LTHE)){
		if (!strcmp(LANGLE, LRAD)){
			fputs("linha 1480", IW);
			THETA = VALUE;
		} else{
			fputs("linha 1484", IW);
			THETA = VALUE * PI / 180.0;
		}
	} else if (!strcmp(LKEYW, LPHI)){
		if (!strcmp(LANGLE, LRAD)){
			fputs("linha 1490", IW);
			PHI = VALUE;
		} else{
			fputs("linha 1494", IW);
			PHI = VALUE * PI / 180.0;
		}
	} else if (!strcmp(LKEYW, LXSH)){
			fputs("linha 1500", IW);
			XSHIFT = VALUE;
	} else if (!strcmp(LKEYW, LYSH)){
			fputs("linha 1504", IW);
			YSHIFT = VALUE;
	} else if (!strcmp(LKEYW, LZSH)){
			fputs("linha 1508", IW);
			ZSHIFT = VALUE;
	} else {
		fputs(BLINE, IW);
		fputs("*** What do you mean?", IW);
		exit(0);	
	}
	goto L404;
	
L405:;
//	printf("linha 2491\n");

	//Determine todos os descendentes do módulo KORIG.
	ND = 1;
	IDESC[1-1] = KORIG; //Descendentes do módulo clonado.
	IDONE[1-1] = 0; //Os descendentes ainda não foram identificados.

L406:;

	NDC = ND;
	KDG = 0;
	for (int I = 1; I <= NDC; I++){
		if (IDONE[I-1] == 0){
			KB = IDESC[I-1];
			if (QBODY_.KBODY[NXG-1][KB-1] > 0){
				for (int J = 1; J <= QBODY_.KBODY[NXG-1][KB-1]; J++ ){
					ND = ND+1;
					IDESC[ND-1] = QBODY_.KBODY[J-1][KB-1]; //Novo descendente
					IDONE[ND-1] = 0;
					KDG = KDG+1;
				}
			} else if (QTREE_.KDGHT[NXG-1][KB-1] > 0){
				for (int J = 1; J <= QBODY_.KBODY[NXG-1][KB-1]; J++ ){
					if (QTREE_.KDGHT[NXG-1][KB-1] != KB){
						ND = ND+1;
						IDESC[ND-1] = QTREE_.KDGHT[J-1][KB-1]; //Novo descendente
						IDONE[ND-1] = 0;
						KDG = KDG+1;
					}
				}
			}
			IDONE[I-1] = 1; //Os descendentes de KB = IDESC (I) foram listados.
		}
	}
	if (KDG > 0)
		goto L406;
	
	IN =0;
	for (int I=1; I<=NB; I++){
		IBCL[I-1] = 0; //Rótulo de um corpo ou módulo clonado.
		IBOR[I-1] = 0; //Label of the original body or module.
	}
	KSD = *QSURF_.NSURF;
	
	for (int I=1; I <= NS; I++){
		ISCL[I-1] = 0; //Label of a cloned surface.
	}
	
	KBD = *QTREE_.NBODYS;
	for (int KBB=1; KBB <= *QTREE_.NBODYS; KBB++){
		for (int ID = ND; ID >=1; ID--){
			KB = IDESC[ID -1];
			if (KBB == KB){
				KBD = KBD+1;
				IBCL[KB-1] = KBD;
				IBOR[KBD -1] = KB;
				PENGEOM_mod_.MATER[KBD-1] = PENGEOM_mod_.MATER[KB-1];
				//Clone as superfícies do módulo original e seus descendentes.
				QTREE_.KSURF[NXG-1][KBD-1] = QTREE_.KSURF[NXG-1][KB-1];
				for (int KSS = QTREE_.KSURF[NXG-1][KB-1]; KSS >=1; KSS--){
					KS = QTREE_.KSURF[KSS-1][KB-1];
					if ((QTREE_.KFLAG[KSS-1][KB-1] < 3) && (ISCL[KS-1] == 0)) {
						if (KSFF[KS-1] == 0){
							KSD = KSD +1;
							if (KSD > NS){
								fputs("*** The parameter NS must be increased.", IW);
								exit(0);
							}
							strcpy(DEFS[KSD-1], DEFS[KS-1]);
							ISCL[KS-1] = KSD;
							strcpy(BLINE, DEFS[KSD-1]);
							QXX=QSURF_.AXX[KS - 1];
							QXY=QSURF_.AXY[KS - 1];
							QXZ=QSURF_.AXZ[KS - 1];
							QYY=QSURF_.AYY[KS - 1];
							QYZ=QSURF_.AYZ[KS - 1];
							QZZ=QSURF_.AZZ[KS - 1];
							QX=QSURF_.AX[KS - 1];
							QY=QSURF_.AY[KS - 1];
							QZ=QSURF_.AZ[KS - 1];
							Q0=QSURF_.A0[KS - 1];
							
							rotshf2_(OMEGA, THETA, PHI, XSHIFT, YSHIFT, ZSHIFT, QXX, QXY, QXZ, QYY, QYZ, QZZ, QX, QY, QZ, Q0);
							
							QSURF_.AXX[KSD - 1]=QXX;
							QSURF_.AXY[KSD - 1]=QXY;
							QSURF_.AXZ[KSD - 1]=QXZ;
							QSURF_.AYY[KSD - 1]=QYY;
							QSURF_.AYZ[KSD - 1]=QYZ;
							QSURF_.AZZ[KSD - 1]=QZZ;
							QSURF_.AX[KSD - 1]=QX;
							QSURF_.AY[KSD - 1]=QY;
							QSURF_.AZ[KSD - 1]=QZ;
							QSURF_.A0[KSD - 1]=Q0;
							KSFF[KSD-1] = 0;
							fputs("linha 1605", IW);
							fputs("linha 1606", IW);
							
							if (fabs(QSURF_.AXX[KSD-1]) > 1e-20){
								fputs("linha 1609", IW);
							}
							if (fabs(QSURF_.AXY[KSD-1]) > 1e-20){
								fputs("linha 1612", IW);
							}
							if (fabs(QSURF_.AXZ[KSD-1]) > 1e-20){
								fputs("linha 1615", IW);
							}
							if (fabs(QSURF_.AYY[KSD-1]) > 1e-20){
								fputs("linha 1618", IW);
							}
							if (fabs(QSURF_.AYZ[KSD-1]) > 1e-20){
								fputs("linha 1621", IW);
							}
							if (fabs(QSURF_.AZZ[KSD-1]) > 1e-20){
								fputs("linha 1624", IW);
							}
							if (fabs(QSURF_.AX[KSD-1]) > 1e-20){
								fputs("linha 1627", IW);
							}
							if (fabs(QSURF_.AY[KSD-1]) > 1e-20){
								fputs("linha 1630", IW);
							}
							if (fabs(QSURF_.AZ[KSD-1]) > 1e-20){
								fputs("linha 1633", IW);
							}
							if (fabs(QSURF_.A0[KSD-1]) > 1e-20){
								fputs("linha 1636", IW);
							}
							fputs("linha 1638", IW);	
						} else{
							ISCL[KS-1]= KS;
						}
					}
				}
				goto L407;
			}
		}
L407:;
	}
	
	//Clone o módulo original e seus descendentes.
	
	for (int KB = *QTREE_.NBODYS+1; KB <= KBD; KB++){
		if (KB > NB){
			fputs("*** The parameter NB must be increased.", IW);
			exit(0);
		}
	
	
		KBO = IBOR[KB-1];
		if (QTREE_.KMOTH[KBO-1] > 0){
			QTREE_.KMOTH[KB-1] = IBCL[QTREE_.KMOTH[KBO-1]-1];
		}else {
			QTREE_.KMOTH[KB-1] = 0;	
		}
		
		QBODY_.KBOMO[KB-1] = QBODY_.KBOMO[KBO-1];
		strcpy(DEFB[KB-1], DEFB[KBO-1]);
		strcpy(BLINE, DEFB[KBO]);
		
		if (QBODY_.KBOMO[KB-1] == 0){
			strcpy(LKEYW, LBOD);
			fputs("linha 1665", IW);
		} else if (QBODY_.KBOMO[KB-1] == 1){ 
			strcpy(LKEYW, LMOD);
			fputs("linha 1668", IW);
		}else {
			fputs("linha 1670", IW);
			fputs("*** Something wrong...", IW);
			exit(0);
		}
		fputs("linha 1674", IW);
		
		for (int KS = 1; KS <= QTREE_.KSURF[NXG-1][KB-1]; KS++){
			QTREE_.KSURF[KS-1][KB-1] = ISCL[QTREE_.KSURF[KS-1][KBO-1]-1];
			QTREE_.KFLAG[KS-1][KB-1] = QTREE_.KFLAG[KS-1][KBO-1];
			if (QTREE_.KFLAG[KS-1][KB-1] < 3){
				if (QTREE_.KFLAG[KS-1][KB-1] == 1){
					INDS = -1;
				} else{
					INDS = +1;
				}
				fputs("LInha 1684", IW);
			}
		}
		
		if (QBODY_.KBOMO[KB-1] == 0){
			
			QBODY_.KBODY[NXG-1][KB-1] = QBODY_.KBODY[NXG-1][KBO-1];
			for (int I = 1; I <= QBODY_.KBODY[NXG-1][KB-1]; I++){
				QBODY_.KBODY[I-1][KB-1] = IBCL[QBODY_.KBODY[I-1][KBO-1]-1];
				KBB = QBODY_.KBODY[I-1][KB-1];
				if (QBODY_.KBOMO[KBB-1] == 0){
					strcpy(LKEYW, LBOD);
					fputs("linha 1695", IW);
				} else if (KBB != KB){
					strcpy(LKEYW, LMOD);
					fputs("linha 1698", IW);	
				}
				if (KBB >= KB){
					fputs("*** The limiting body or module is not yet defined", IW);
					exit(0);
				}
			} 
			
		} else{
			QTREE_.KDGHT[NXG-1][KB-1] = QTREE_.KDGHT[NXG-1][KBO-1];
			for (int I = 1; I <= QTREE_.KDGHT[NXG-1][KB-1]; I++){
				QTREE_.KDGHT[I-1][KB-1] = IBCL[QTREE_.KDGHT[I-1][KBO-1]-1];
				KBB = QTREE_.KDGHT[I-1][KB-1];
				if (KBB != KB){
					if (QBODY_.KBOMO[KBB-1] == 0){
						strcpy(LKEYW, LBOD);
						fputs("linha 1714", IW);
					} else if (KBB != KB){
						strcpy(LKEYW, LMOD);
						fputs("linha 1717", IW);	
					}
					if (KBB >= KB){
						fputs("*** The limiting body or module is not yet defined", IW);
						exit(0);
					}
				}
			}  
		}
		fputs("linha 1727", IW);
		if (KB == KBD){
			fputs("Linha 1729", IW);
			fputs("C ************  End of cloned module.9", IW);
			fputs("Linha 1731", IW);
			fputs("Linha 1732", IW);
			
		}
	}
	strcpy(ALIAB[KBD], C5C);
	QTREE_.KMOTH[KBD-1]=0;
	*QTREE_.NBODYS=KBD;
	*QSURF_.NSURF=KSD;
	
	goto L2;
	
	//Arquivo de geometria incluído.

L500:;
	
	if (!strcmp(LKEYW, LINA)){
		KEEPL=1; //Mantenha os rótulos do usuário dos elementos de arquivo incluídos.
	} else{
		KEEPL=0;
	}
	fputs("linha 1749", IW);
	fputs("linha 1750", IW);
	fputs("C ************  Included file:  ", IW);
	if (KEEPL == 1)
		fputs("The included elements keep their user labels", IW);
	fputs("linha 1754", IW);
	fputs("linha 1755", IW);
L501:;
	fgets(BLINE, sizeof(BLINE), IR);
	if (comentario(BLINE)){

		   fputs(BLINE, IW);
		   goto L501;   
	}
	
	if ((BLINE[8] != '(' ) || (BLINE[13] != ')')){
		fputs(BLINE, IW);
		fputs("*** Incorrect label format.", IW);
		exit(0);
	}
	extrairString(LKEYW, BLINE, 0, 8);
	extrairString(GFILE, BLINE, 9, 20);
	
	if (strcmp(LKEYW, LFIL)){
		fputs(BLINE, IW);
		fputs("*** What do you mean?", IW);
		exit(0);
	} else {
		fputs("Linha 1776", IW);
	}
	
L502:;

	fgets(BLINE, sizeof(BLINE), IR);
	if (comentario(BLINE)){

		   fputs(BLINE, IW);
		   goto L502;   
	}
	
	extrairString(LKEYW, BLINE, 0, 8);
	if (!strcmp(LKEYW, LNUL)){
		
		//codigo reduzido para abrir um arquivo incluido no arquivo de geometria linhas 1788 até 1817
		freopen(GFILE, "r", IR);		
	}
	
	goto L1;
	
	//Linha final no arquivo de entrada.
	
L600:;
//	fclose(IR);
//	fputs("C\n", IW);
//	fputs("C ************  End of included file.\n", IW);
//	fputs("linha 1825\n", IW);
//	fputs("linha 1826\n", IW);
	
//	C1[0] = '0';
	//goto L2;
	
	if (*QTREE_.NBODYS == 1) {
		QTREE_.KDGHT[NXG-1][1-1] = 1;
		QTREE_.KDGHT[1-1][1-1] = 1;	
	}
	
	for (int KB=1; KB <= *QTREE_.NBODYS; KB++){
		extrairString(C5, ALIAB[KB-1], 0, 5);
	//	strncpyy(APOIO, ALIAB[KB-1], 5);
	//	APOIO[5] = '\0';
		//strcpy(C5, APOIO);
		if (C5[4] == '0'){
		//	strcpy(PENGEOM_mod_.BALIAS[KB-1], 
			extrairString(PENGEOM_mod_.BALIAS[KB-1], C5, 0, 4);
			//strncpyy(PENGEOM_mod_.BALIAS[KB-1], C5, 4); 
		//	PENGEOM_mod_.BALIAS[KB-1][4] = '\0';
		//	printf("BALIAS %s   -   C5: %s\n", PENGEOM_mod_.BALIAS[KB-1], C5);
			

		}
			
	}
	
	//Verifique se há corpos ou módulos sem mãe.
	
	MLESS = 0;
	KBENC = 0;
	for (int KBB=1; KBB<= *QTREE_.NBODYS; KBB++){
		if (QTREE_.KMOTH[KBB-1] == 0){
			MLESS = MLESS +1;
			KBENC = KBB;
		}
	}
	if (MLESS == 1){
		if (QBODY_.KBOMO[KBENC-1] == 1){
		//	printf("\n\nPULOU l602\n\n");
			goto L602; //existe um modulo raiz.
		}
	}
	
	//Defina o gabinete.
	if (*QTREE_.NBODYS > (NB -1)){
		fputs("** The parameter NB must be increased.", IW);
		exit(0);
	}
	if (*QSURF_.NSURF > (NS-1)){
		fputs("** The parameter NS must be increased.", IW);
		exit(0);
	}
	
	//A próxima linha serve apenas para evitar um aviso emitido por certos compiladores.
	NT = 0;
	KS = *QSURF_.NSURF+1;
	
	//O recinto é uma esfera centrada na origem e com um Raio de unidades de comprimento 1.0E7.
	
	QSURF_.AXX[KS-1] = 1.0e0;
	QSURF_.AYY[KS-1] = 1.0e0;
	QSURF_.AZZ[KS-1] = 1.0e0;
	QSURF_.A0[KS-1] = -1.0e14;
	KB = *QTREE_.NBODYS + 1;
	QBODY_.KBOMO[KB-1] = 1;
	QTREE_.KSURF[NXG-1][KB-1] = 1;
	QTREE_.KSURF[1-1][KB-1] = KS;
	QTREE_.KFLAG[1-1][KB-1] = 1;
	QTREE_.KDGHT[NXG-1][KB-1] = 1;
	QTREE_.KDGHT[1-1][KB-1] = KB;
	
//	printf("\n\n  600 \n\n");
	
	for (int KBB=1; KBB <= *QTREE_.NBODYS; KBB++){
		if (QTREE_.KMOTH[KBB-1] == 0){
			NT = QTREE_.KDGHT[NXG-1][KB-1] +1;
			QTREE_.KDGHT[NXG-1][KB-1] = NT;
			QTREE_.KDGHT[NT-1][KB-1] = KBB;
			KN1 = QTREE_.KSURF[NXG-1][KBB-1];
			for (int KS1 = 1; KS1 <= KN1; KS1++){
				if (QTREE_.KFLAG[KS1-1][KBB-1] > 3){
				//	printf("\n\n\nPULOOOUUU\n\n\n");
					goto L601;
				}
					
				KSURF1= QTREE_.KSURF[KS1-1][KBB-1];
				KN2 = QTREE_.KSURF[NXG-1][KB-1];
				for (int KS2 = 1; KS2 <= KN2; KS2++){
					if (QTREE_.KSURF[KS2-1][KB-1] == KSURF1){
						//	printf("\n\n\nPULOOOUUU\n\n\n");
						goto L601;
					}
				}
				KN2 = KN2+1;
				if (KN2 >= NXG){
					fputs("*** The parameter NXG is too small.", IW);
					exit(0);
				}
				QTREE_.KSURF[NXG-1][KB-1] = KN2;
				QTREE_.KSURF[KN2-1][KB-1] = KSURF1;
				QTREE_.KFLAG[KN2-1][KB-1] = 4;
			//	printf("\n\n\nEXECUTOUUU\n\n\n");
L601:;			
			}
			QTREE_.KMOTH[KBB-1] = KB;
			KN3 = KN1+1;
			QTREE_.KSURF[NXG-1][KBB-1] = KN3;
			QTREE_.KSURF[KN3-1][KBB-1] = KS;
			QTREE_.KFLAG[KN3-1][KBB-1] = 1;
		//	printf("\n\n\nPULOU\n\n\n");
		}
	}
	
	*QSURF_.NSURF = KS;
	*QTREE_.NBODYS = KB;
	
	//Classifique as filhas em ordem crescente.
	
	
	if (NT > 1){
			for (int KI = 1; KI <= NT-1; KI++){
				KBMIN = QTREE_.KDGHT[KI-1][*QTREE_.NBODYS-1];
				KMIN = KI;
				for (int KJ = KI+1; KJ <= NT; KJ++){
					if (QTREE_.KDGHT[KJ-1][*QTREE_.NBODYS-1] < KBMIN){
						KBMIN = QTREE_.KDGHT[KJ-1][*QTREE_.NBODYS-1];
						KMIN = KJ;
					}
				}
				if (KMIN != KI){
					KSAVE = QTREE_.KDGHT[KI-1][*QTREE_.NBODYS-1];
					QTREE_.KDGHT[KI-1][*QTREE_.NBODYS-1] = QTREE_.KDGHT[KMIN-1][*QTREE_.NBODYS-1];
					QTREE_.KDGHT[KMIN-1][*QTREE_.NBODYS-1] = KSAVE;
				}
			}
		}
		
L602:;

	fprintf(IW,"%s 0000000000000000000000000000000000000000000000000000000\n\n\n", LKEYW);
	//fputs("linha 1933\n", IW);	
	*NBOD = *QTREE_.NBODYS;
	*PENGEOM_mod_.NBODY = *QTREE_.NBODYS;	
	
	//Superfícies duplicadas (dentro da tolerância de arredondamento) são removidas.
	
	fputs("*****************************************\n", IW);
	fputs("****     PENGEOM (version 2014)      ****\n", IW);
	fputs("****  Constructive Quadric Geometry  ****\n", IW);
	fputs("*****************************************\n", IW);
	
	if (NSFF > 0){		
		fputs("WARNING: The system contains fixed (starred) surfaces, which are'',/9X,''not affected by translations and rotations. Hence, any translation or rotation that does not leave these surfaces, invariant will distort the system.\n", IW);
	}
	
	if (*QSURF_.NSURF < 2)
		goto L704;
	IWRITE = 0;
	TOL = 1.0e-14;
	
	for (int KS = 1; KS <= *QSURF_.NSURF; KS++){
		KM[KS-1] = 0;
	}
	
//	printf(" linha 3178  NSURF: %d",*QSURF_.NSURF );
	
	
	for (int KS =1 ; KS <= *QSURF_.NSURF-1; KS++){
		if (KM[KS-1] != 0)
			goto L703;
		F = fmax(fmax(QSURF_.AXX[KS-1],QSURF_.AXY[KS-1]), fmax(fmax(fmax(fmax(fmax(fmax(fmax(QSURF_.AXZ[KS-1],QSURF_.AYY[KS-1]),QSURF_.AYZ[KS-1]),QSURF_.AZZ[KS-1]),QSURF_.AX[KS-1]),QSURF_.AY[KS-1]),QSURF_.AZ[KS-1]),QSURF_.A0[KS-1]));
		FM = fmin(fmin(QSURF_.AXX[KS-1],QSURF_.AXY[KS-1]), fmin(fmin(fmin(fmin(fmin(fmin(fmin(QSURF_.AXZ[KS-1],QSURF_.AYY[KS-1]),QSURF_.AYZ[KS-1]),QSURF_.AZZ[KS-1]),QSURF_.AX[KS-1]),QSURF_.AY[KS-1]),QSURF_.AZ[KS-1]),QSURF_.A0[KS-1]));
		if (fabs(FM) > fabs(F))
		   	F = FM;
		for (int KST = KS+1; KST <= *QSURF_.NSURF; KST++){
			if (KM[KST-1] != 0)
				goto L702;
			FP = fmax(fmax(QSURF_.AXX[KST-1],QSURF_.AXY[KST-1]), fmax(fmax(fmax(fmax(fmax(fmax(fmax(QSURF_.AXZ[KST-1],QSURF_.AYY[KST-1]),QSURF_.AYZ[KST-1]),QSURF_.AZZ[KST-1]),QSURF_.AX[KST-1]),QSURF_.AY[KST-1]),QSURF_.AZ[KST-1]),QSURF_.A0[KST-1]));
			FM = fmin(fmin(QSURF_.AXX[KST-1],QSURF_.AXY[KST-1]), fmin(fmin(fmin(fmin(fmin(fmin(fmin(QSURF_.AXZ[KST-1],QSURF_.AYY[KST-1]),QSURF_.AYZ[KST-1]),QSURF_.AZZ[KST-1]),QSURF_.AX[KST-1]),QSURF_.AY[KST-1]),QSURF_.AZ[KST-1]),QSURF_.A0[KST-1]));
			if (fabs(FM) > fabs(FP))
				FP = FM;
			FFP = F/FP;
			RFFP = 1.0e0/FFP;
			
			TST = 0.0e0;
			
			if (fabs(QSURF_.AXX[KS-1]) > 1.0e-16){
				TST = fmax(TST, fabs((QSURF_.AXX[KS-1] - QSURF_.AXX[KST-1] * FFP) / QSURF_.AXX[KS-1]));
			} else if (fabs(QSURF_.AXX[KST-1]) > 1.0e-16){
				TST = fmax(TST, fabs((QSURF_.AXX[KS-1] * RFFP - QSURF_.AXX[KST-1] ) / QSURF_.AXX[KST-1]));
			}
			if (fabs(QSURF_.AXY[KS-1]) > 1.0e-16){
				TST = fmax(TST, fabs((QSURF_.AXY[KS-1] - QSURF_.AXY[KST-1] * FFP) / QSURF_.AXY[KS-1]));
			} else if (fabs(QSURF_.AXY[KST-1]) > 1.0e-16){
				TST = fmax(TST, fabs((QSURF_.AXY[KS-1] * RFFP - QSURF_.AXY[KST-1]) / QSURF_.AXY[KST-1]));
			}
			if (fabs(QSURF_.AXZ[KS-1]) > 1.0e-16){
				TST = fmax(TST, fabs((QSURF_.AXZ[KS-1] - QSURF_.AXZ[KST-1] * FFP) / QSURF_.AXZ[KS-1]));
			} else if (fabs(QSURF_.AXZ[KST-1]) > 1.0e-16){
				TST = fmax(TST, fabs((QSURF_.AXZ[KS-1] * RFFP - QSURF_.AXZ[KST-1]) / QSURF_.AXZ[KST-1]));
			}
			if (fabs(QSURF_.AYY[KS-1]) > 1.0e-16){
				TST = fmax(TST, fabs((QSURF_.AYY[KS-1] - QSURF_.AYY[KST-1] * FFP) / QSURF_.AYY[KS-1]));
			} else if (fabs(QSURF_.AYY[KST-1]) > 1.0e-16){
				TST = fmax(TST, fabs((QSURF_.AYY[KS-1] * RFFP - QSURF_.AYY[KST-1]) / QSURF_.AYY[KST-1]));
			}
			if (fabs(QSURF_.AYZ[KS-1]) > 1.0e-16){
				TST = fmax(TST, fabs((QSURF_.AYZ[KS-1] - QSURF_.AYZ[KST-1] * FFP) / QSURF_.AYZ[KS-1]));
			} else if (fabs(QSURF_.AYZ[KST-1]) > 1.0e-16){
				TST = fmax(TST, fabs((QSURF_.AYZ[KS-1] * RFFP - QSURF_.AYZ[KST-1]) / QSURF_.AYZ[KST-1]));
			}
			if (fabs(QSURF_.AZZ[KS-1]) > 1.0e-16){
				TST = fmax(TST, fabs((QSURF_.AZZ[KS-1] - QSURF_.AZZ[KST-1] * FFP) / QSURF_.AZZ[KS-1]));
			} else if (fabs(QSURF_.AZZ[KST-1]) > 1.0e-16){
				TST = fmax(TST, fabs((QSURF_.AZZ[KS-1] * RFFP - QSURF_.AZZ[KST-1]) / QSURF_.AZZ[KST-1]));
			}
			if (fabs(QSURF_.AX[KS-1]) > 1.0e-16){
				TST = fmax(TST, fabs((QSURF_.AX[KS-1] - QSURF_.AX[KST-1] * FFP) / QSURF_.AX[KS-1]));
			} else if (fabs(QSURF_.AX[KST-1]) > 1.0e-16){
				TST = fmax(TST, fabs((QSURF_.AX[KS-1] * RFFP - QSURF_.AX[KST-1]) / QSURF_.AX[KST-1]));
			}
			if (fabs(QSURF_.AY[KS-1]) > 1.0e-16){
				TST = fmax(TST, fabs((QSURF_.AY[KS-1] - QSURF_.AY[KST-1] * FFP) / QSURF_.AY[KS-1]));
			} else if (fabs(QSURF_.AY[KST-1]) > 1.0e-16){
				TST = fmax(TST, fabs((QSURF_.AY[KS-1] * RFFP - QSURF_.AY[KST-1]) / QSURF_.AY[KST-1]));
			}
			if (fabs(QSURF_.AZ[KS-1]) > 1.0e-16){
				TST = fmax(TST, fabs((QSURF_.AZ[KS-1] - QSURF_.AZ[KST-1] * FFP) / QSURF_.AZ[KS-1]));
			} else if (fabs(QSURF_.AZ[KST-1]) > 1.0e-16){
				TST = fmax(TST, fabs((QSURF_.AZ[KS-1] * RFFP - QSURF_.AZ[KST-1]) / QSURF_.AZ[KST-1]));
			}
			if (fabs(QSURF_.A0[KS-1]) > 1.0e-16){
				TST = fmax(TST, fabs((QSURF_.A0[KS-1] - QSURF_.A0[KST-1] * FFP) / QSURF_.A0[KS-1]));
			} else if (fabs(QSURF_.A0[KST-1]) > 1.0e-16){
				TST = fmax(TST, fabs((QSURF_.A0[KS-1] * RFFP - QSURF_.A0[KST-1]) / QSURF_.A0[KST-1]));
			}
			
			if (TST < TOL){
				if (IWRITE == 0){
					fputs("************  Removal of duplicated user-defined surfaces.\n", IW);
					IWRITE = 1;
				}
				fprintf(IW, "(SURFACE (%d) is replaced by SURFACE (%d) linha 2040\n", KST, KS);
				
				//fputs("linha 2033 e 2034\n", IW);
				/*
				Verifique se as duas funções de superfície têm o mesmo global
				Sinal C (ISPF = 0) ou não (ISPF = 1).	
				*/
				ISPF = 0;
				TF = TOL*F;
				
				if (((fabs(QSURF_.AXX[KS-1])) > TF) && (QSURF_.AXX[KS-1]*QSURF_.AXX[KST-1] < 0.0)) 
					ISPF = 1;
				if (((fabs(QSURF_.AXY[KS-1])) > TF) && (QSURF_.AXY[KS-1]*QSURF_.AXY[KST-1] < 0.0))
					ISPF = 1;
				if (((fabs(QSURF_.AXZ[KS-1])) > TF) && (QSURF_.AXZ[KS-1]*QSURF_.AXZ[KST-1] < 0.0)) 
					ISPF = 1;
				if (((fabs(QSURF_.AYY[KS-1])) > TF) && (QSURF_.AYY[KS-1]*QSURF_.AYY[KST-1] < 0.0))
					ISPF = 1;
				if (((fabs(QSURF_.AYZ[KS-1])) > TF) && (QSURF_.AYZ[KS-1]*QSURF_.AYZ[KST-1] < 0.0)) 
					ISPF = 1;
				if (((fabs(QSURF_.AZZ[KS-1])) > TF) && (QSURF_.AZZ[KS-1]*QSURF_.AZZ[KST-1] < 0.0))
					ISPF = 1;
				if (((fabs(QSURF_.AX[KS-1])) > TF) && (QSURF_.AX[KS-1]*QSURF_.AX[KST-1] < 0.0))
					ISPF = 1;
				if (((fabs(QSURF_.AY[KS-1])) > TF) && (QSURF_.AY[KS-1]*QSURF_.AY[KST-1] < 0.0))
					ISPF = 1;
				if (((fabs(QSURF_.AZ[KS-1])) > TF) && (QSURF_.AZ[KS-1]*QSURF_.AZ[KST-1] < 0.0))
					ISPF = 1;
				if (((fabs(QSURF_.A0[KS-1])) > TF) && (QSURF_.A0[KS-1]*QSURF_.A0[KST-1] < 0.0))
					ISPF = 1;
				
				if (TST > 1.0e-15){
					fputs("linha 2051\n", IW);
					fputs("linha 2052\n", IW);
					fputs("linha 2054\n", IW);
					fputs("linha 2056\n", IW);
					fputs("linha 2058\n", IW);
					fputs("linha 2060\n", IW);
					fputs("linha 2062\n", IW);
					fputs("linha 2064\n", IW);
					fputs("linha 2066\n", IW);
					fputs("linha 2068\n", IW);
					fputs("linha 2070\n", IW);
					fputs("linha 2072\n", IW);
				}

				QSURF_.AXX[KST - 1] = QSURF_.AXX[KS - 1];
				QSURF_.AXY[KST - 1] = QSURF_.AXY[KS - 1];
				QSURF_.AXZ[KST - 1] = QSURF_.AXZ[KS - 1];
				QSURF_.AYY[KST - 1] = QSURF_.AYY[KS - 1];
				QSURF_.AYZ[KST - 1] = QSURF_.AYZ[KS - 1];
				QSURF_.AZZ[KST - 1] = QSURF_.AZZ[KS - 1];
				QSURF_.AX[KST - 1] = QSURF_.AX[KS - 1];
				QSURF_.AY[KST - 1] = QSURF_.AY[KS - 1];
				QSURF_.AZ[KST - 1] = QSURF_.AZ[KS - 1];
				QSURF_.A0[KST - 1] = QSURF_.A0[KS - 1];
				
				QSURF_.KPLANE[KST-1] = QSURF_.KPLANE[KS-1];
				KM[KST-1] = KS;
				
				KBI = 1;
				
L701:;		

				for (int KB = KBI; KB <= *QTREE_.NBODYS; KB++){
					KSL = 0;
					KSLP = 0;
					NSB = QTREE_.KSURF[NXG-1][KB-1];
					for ( int K = 1; K <= NSB; K++){
						if (QTREE_.KSURF[K-1][KB-1] == KS)
							KSL = K;
						if (QTREE_.KSURF[K-1][KB-1] == KST)
							KSLP = K;
					}
					if (KSLP > 0){
						QTREE_.KSURF[KSLP-1][KB-1] = KS;
						//Se as equações implícitas das superfícies KS e KST diferem por um Sinal global, o ponteiro lateral do KST deve ser invertido.
						if (ISPF == 1){
							if (QTREE_.KFLAG[KSLP-1][KB-1] == 1){
								QTREE_.KFLAG[KSLP-1][KB-1] = 2;
							}else if (QTREE_.KFLAG[KSLP-1][KB-1] == 2){
								QTREE_.KFLAG[KSLP-1][KB-1] = 1;
							}
						}
						if (KSL > 0){
							KFL = QTREE_.KFLAG[KSL-1][KB-1];
							KFLP = QTREE_.KFLAG[KSLP-1][KB-1];
							if ((min(KFL, KFLP)) < 3){
								QTREE_.KFLAG[KSL-1][KB-1] = min(KFL, KFLP);
								
								if ((KFL != KFLP) && (max(KFL, KFLP) < 3)){
									if (QBODY_.KBOMO[KB-1] == 0){
										fprintf(IW, "ERROR: BODY(%d) is limited by two equivalent surfaces. Probably\n  this body cannot be resolved because it is small and located,far from the origin.\n", KB);
										//fputs(" linha 2118, 2119, 2120\n", IW);
									} else{
										fprintf(IW, "ERROR: MODULE(%d) is limited by two equivalent surfaces. Probably\n  this modulo cannot be resolved because it is small and located,far from the origin.\n", KB);
								
										//fputs(" linha 2122, 2123, 2124, 2125\n", IW);	
									}
									exit(0);
									
									
								}
							} else if (QBODY_.KBOMO[KB-1] == 0){
								QTREE_.KFLAG[KSL-1][KB-1] = 3; 
							} else if (QBODY_.KBOMO[KB-1] == 1){
								QTREE_.KFLAG[KSL-1][KB-1] = 4;
							}else {
								QTREE_.KFLAG[KSL-1][KB-1] = 5;
							}
							
							for (int K = KSLP; K <= NSB-1; K++){
								QTREE_.KSURF[K-1][KB-1] = QTREE_.KSURF[K+1-1][KB-1];
								QTREE_.KFLAG[K-1][KB-1] = QTREE_.KFLAG[K+1-1][KB-1];
							}
							QTREE_.KSURF[NSB-1][KB-1] = 0;
							QTREE_.KFLAG[NSB-1][KB-1] = 5;

							QTREE_.KSURF[NXG-1][KB-1] = NSB-1;
							printf("\n\nKSURF NXG-1 NB-1 = %d\n\n", NSB-1);
							if (KB < *QTREE_.NBODYS) {
								KBI = KB;
								goto L701;
							}	
						}
						
					}
				}
			}
L702:;
		}
L703:;
	}
L704:;

	fputs("\n\n************  Genealogical tree.\n\n ", IW);
	
	for (int KB = 1; KB <= *QTREE_.NBODYS; KB++){
		if (QBODY_.KBOMO[KB-1] == 0){
			//fputs("linha 2159 e 2160\n", IW);
			fprintf(IW, "\n*** BODY   = %5d,  KMOTH = %5d,  MAT = %3d\n", KB, QTREE_.KMOTH[KB-1], PENGEOM_mod_.MATER[KB-1]);
			if (QBODY_.KBODY[NXG-1][KB-1] > 0){
				//fputs("linha 2162 e 2163\n", IW);
				imprimirKBODY(IW, KB);
			}
		} else if (QBODY_.KBOMO[KB-1] == 1){
			//fputs("linha 2165 e 2166\n", IW);
			fprintf(IW, "\n*** MODULE   = %5d,  KMOTH =%5d,  MAT =%3d\n", KB, QTREE_.KMOTH[KB-1], PENGEOM_mod_.MATER[KB-1]);
			imprimirKDGHT(IW, KB);
			//fputs("linha 2167 e 2168\n", IW);
		} else{
			fprintf(IW, "\n*** ERROR: the label %d does not correspond to a body.", KB);
		//	fputs("linha 2170\n", IW);
			exit(0);
		}
		//fprintf(IW, "KSURF = %d", QTREE_.KSURF);
		imprimirKSURF(IW, KB);
		imprimirKFLAG(IW, KB);

	
	
		//fputs("linha 2174\n", IW);
	//	fputs("linha 2176\n", IW);	
	}
		
	
	if (MLESS == 1){
		if (QBODY_.KBOMO[KBENC-1] == 1){
			fprintf(IW, "\nThe module %5d is the enclosure.\n\n", KBENC);
		//	fputs("The module'',I5,'' is the enclosure.\n", IW);
		}
	}

	
	//Teste de consistência de superfície (F. Tola).
	for (int KB = 1; KB <= *QTREE_.NBODYS-1; KB++){
		KB1 = QTREE_.KMOTH[KB-1];
		for (int I = 1; I <= QTREE_.KSURF[NXG-1][KB-1]; I++){
			KS = QTREE_.KSURF[I-1][KB-1];
			KF = QTREE_.KFLAG[I-1][KB-1];
			for (int J = 1; J <= QTREE_.KSURF[NXG-1][KB1-1]; J++){
				if (QTREE_.KSURF[J-1][KB1-1] == KS){
					//Surface tem filha e mãe em lados opostos.
					if (((KF == 1) && (QTREE_.KSURF[J-1][KB1-1] == 2)) || ((KF == 2) &&  (QTREE_.KSURF[J-1][KB1-1] == 1))) {
						if (QBODY_.KBOMO[KB-1] == 0){
							fputs("linha 2200 ERROR: the SURFACE  which limits BODY as inconsistent side pointers\n", IW);
						} else{
							fputs("linha 2204 ERROR: the SURFACE  which limits MODULE as inconsistent side pointers\n", IW);
								
						}
						exit(0);
					}
					goto L801;
				}
			}
L801:;
		}
	}
		
	
	//Teste de facilidade.
	
	NBU = 0;
	NSU =0;
	for (int KB = 1; KB <= *QTREE_.NBODYS; KB++ ){
		if (QBODY_.KBOMO[KB-1] == 0){
			NBU = max(NBU, QBODY_.KBODY[NXG-1][KB-1]);
		}else if (QBODY_.KBOMO[KB-1] == 1){
			NBU = max(NBU, QTREE_.KDGHT[NXG-1][KB-1]);
		}
		NSE = 0;
		for (int K = 1; K<= QTREE_.KSURF[NXG-1][KB-1]; K++){
			if (QTREE_.KFLAG[K-1][KB-1] < 5){
				NSE = NSE + 1;
			}
		}
		NSU = max(NSU, NSE);		
	}

	
	
	fputs("\n************ Adequacy of the geometry definition\n", IW);
	fputs("\nThe largest number of bodies in a module or\n", IW);
	fprintf(IW, "		bodies limiting a single body is ........%4d\n", NBU);
	fputs("The largest number of limiting surfaces for\n", IW);
	fprintf(IW, "		a single body or module is ..............%4d\n", NSU);
	
	if ((*QTREE_.NBODYS < 15) && (*QSURF_.NSURF < 15)){
		fputs("\nThe simulation of this geometry will be relatively fast\n", IW);
		fputs("     no further optimization seems to be required\n", IW);
	} else if ((NBU < 10) && (NSU < 10)){		
	   	fputs("\nThe simulation of this geometry will be relatively fast,\n", IW);
		fputs("no further optimization seems to be re\n", IW);
	} else if ((NBU < 15) && (NSU < 20)){		
	   	fputs("\nThe simulation of this geometry is expected be slow\n", IW);
		fputs(" try to split complex bodies into several modules\n", IW);
	} else if ((NBU < 25) && (NSU < 30)){		
	   	fputs("\nThe simulation of this geometry will be very slow, you should\n", IW);
		fputs(" try to optimize the structure of the tree of modules\n", IW);
	} else {
		fputs("\nSimulating this geometry will be extremely slow\n", IW);	
	}
	fputs("\n************  The end.\n", IW);
	
	return;
	
L900:;
	fputs(BLINE, IW);
	fputs("\n*** Wrong input format.\n", IW);
	exit(0);

}
	
void rotshf2_(double &OMEGA, double &THETA, double &PHI, double &DX, double &DY, double &DZ, double &AXX, double &AXY, double &AXZ, double &AYY, double &AYZ, double &AZZ, double &AX, double &AY, double &AZ, double &A0){
	
	/*
	
	Esta sub-rotina gira e desloca uma superfície quádrica.
	Parâmetros de entrada C:
	OMEGA, THETA, PHI ... ângulos de rotação de Euler,
	DX, DY, DZ .......... componentes do vetor de deslocamento,
	AXX, ..., A0 ........ coeficientes da quádrica inicial.
	
	Parâmetros de saída C:
	AXX, ..., A0 ........ coeficientes da quádrica transformada.

	*/
	
	double R[3][3];
	double A2[3][3];
	double B2[3][3];
	double A1[3];
	double B1[3];
	double D1[3];
	double B0;
	double A2D;

	
/*	B2[0][0] = AXX;
	B2[1][0] = 0.5e0*AXY;
	B2[2][0] = 0.5e0*AXZ;
	B2[0][1] = B2[1][0];
	B2[1][1] = AYY;
	B2[2][1] = 0.5*AYZ;
	B2[0][2] = B2[2][0];
	B2[1][2] = B2[2][1];
	B2[2][2] = AZZ;
	B1[0] = AX;
	B1[1] = AY;
	B1[2] = AZ;*/
	
	B2[0][0] = AXX;
	B2[0][1] = 0.5e0*AXY;
	B2[0][2] = 0.5e0*AXZ;
	B2[1][0] = B2[0][1];
	B2[1][1] = AYY;
	B2[1][2] = 0.5e0*AYZ;
	B2[2][0] = B2[0][2];
	B2[2][1] = B2[1][2];
	B2[2][2] = AZZ;
	B1[0] = AX;
	B1[1] = AY;
	B1[2] = AZ;
	B0 = A0;
	D1[0]= DX;
	D1[1] = DY;
	D1[2] = DZ;
	
	// Rotação da Matrix

	double STHETA=sin(THETA);
	double CTHETA=cos(THETA);
	double SPHI=sin(PHI);
	double CPHI=cos(PHI);
	double SOMEGA=sin(OMEGA);
	double COMEGA=cos(OMEGA);
	
/*	R[0][0]=CPHI*CTHETA*COMEGA-SPHI*SOMEGA;
	R[1][0]=-CPHI*CTHETA*SOMEGA-SPHI*COMEGA;
	R[2][0]=CPHI*STHETA;
	R[0][1]=SPHI*CTHETA*COMEGA+CPHI*SOMEGA;
	R[1][1]=-SPHI*CTHETA*SOMEGA+CPHI*COMEGA;
	R[2][1]=SPHI*STHETA;
	R[0][2]=-STHETA*COMEGA;
	R[1][2]=STHETA*SOMEGA;
	R[2][2]=CTHETA;	*/
	
	R[0][0]=CPHI*CTHETA*COMEGA-SPHI*SOMEGA;
	R[0][1]=-CPHI*CTHETA*SOMEGA-SPHI*COMEGA;
	R[0][2]=CPHI*STHETA;
	R[1][0]=SPHI*CTHETA*COMEGA+CPHI*SOMEGA;
	R[1][1]=-SPHI*CTHETA*SOMEGA+CPHI*COMEGA;
	R[1][2]=SPHI*STHETA;
	R[2][0]=-STHETA*COMEGA;
	R[2][1]=STHETA*SOMEGA;
	R[2][2]=CTHETA;	
	
	
	
	
	// Rotação quadratica
	
/*	for (int I = 1; I <=3; I++){
		A1[I-1] = 0.0;
		for (int J=1; J<=3; J++){
			A1[I-1] = A1[I-1] + R[J-1][I-1]*B1[J-1];
			A2[J-1][I-1] = 0.0e0;
			for (int M = 1; M <=3; M++){
				for (int K = 1; K<=3; K++){
					A2[J-1][I-1] = A2[J-1][I-1] + R[K-1][I-1]*B2[M-1][K-1]*R[M-1][J-1];
				}
			}
		}
	}*/
	
	for (int I = 1; I <=3; I++){
		A1[I-1] = 0.0e0;
		for (int J=1; J<=3; J++){
			A1[I-1] = A1[I-1] + R[I-1][J-1]*B1[J-1];
			A2[I-1][J-1] = 0.0e0;
			for (int M = 1; M <=3; M++){
				for (int K = 1; K<=3; K++){
					A2[I-1][J-1] = A2[I-1][J-1] + R[I-1][K-1]*B2[K-1][M-1]*R[J-1][M-1];
				}
			}
		}
	}
	
	
	
	//Quadratica com rotação deslocada
	
/*	for (int I = 1; I <= 3; I++){
		A2D = 0.0e0;
		for (int J = 1; J<= 3; J++){
			A2D = A2D + A2[J-1][I-1]*D1[J-1];
		}
		B1[I-1] = A1[I-1]-2.0*A2D;
		B0 = B0+D1[I-1]*(A2D-A1[I-1]);
	}*/
	
	
	for (int I = 1; I <= 3; I++){
		A2D = 0.0e0;
		for (int J = 1; J<= 3; J++){
			A2D = A2D + A2[I-1][J-1]*D1[J-1];
		}
		B1[I-1] = A1[I-1]-2.0e0*A2D;
		B0 = B0+D1[I-1]*(A2D-A1[I-1]);
	}
	
	
/*	AXX=A2[0][0];
	AXY=A2[1][0]+A2[0][1];
	AXZ=A2[2][0]+A2[0][2];
	AYY=A2[1][1];
	AYZ=A2[2][1]+A2[1][2];
	AZZ=A2[2][2];
	AX=B1[0];
	AY=B1[1];
	AZ=B1[2];*/
	
	AXX=A2[0][0];
	AXY=A2[0][1]+A2[1][0];
	AXZ=A2[0][2]+A2[2][0];
	AYY=A2[1][1];
	AYZ=A2[1][2]+A2[2][1];
	AZZ=A2[2][2];
	AX=B1[0];
	AY=B1[1];
	AZ=B1[2];
	
	
	
	A0=B0;
	if (fabs(AXX) < 1.0e-16) AXX=0.0e0;
	if(fabs(AXY) < 1.0e-16) AXY=0.0e0;
	if(fabs(AXZ) < 1.0e-16) AXZ=0.0e0;
	if(fabs(AYY) < 1.0e-16) AYY=0.0e0;
	if(fabs(AYZ) < 1.0e-16) AYZ=0.0e0;
	if(fabs(AZZ) < 1.0e-16) AZZ=0.0e0;
	if(fabs(AX) < 1.0e-16) AX=0.0e0;
	if(fabs(AY) < 1.0e-16) AY=0.0e0;
	if(fabs(AZ) < 1.0e-16) AZ=0.0e0;
	if(fabs(A0) < 1.0e-16) A0=0.0e0;		
	
}


bool comentario(char BLINE[]){
	
	char CCR = 13;
	char CNL = 10;
	char C[] = "C";
	char c[] = "c";
	char A[] = "#";
//	char APOIO[9]; 
	
	bool comentario = false;
	
	if (((BLINE[0] == C[0]) || (BLINE[0] == c[0])) && ((BLINE[1] == CCR) || (BLINE[1] == CNL) || (BLINE[1] == A[0]))){
		comentario = true;
	}
	
	
/*	if (((!strcmp(strncpyy(APOIO, BLINE+0, 1), "C")) || (!strcmp(strncpyy(APOIO, BLINE+0, 1), "c"))) && 
	   (!strcmp(strncpyy(APOIO, BLINE+1, 1), CCR)) || (!strcmp(strncpyy(APOIO, BLINE+1, 1), CNL)) || (!strcmp(strncpyy(APOIO, BLINE+1, 1), " ")) ||
       (!strcmp(strncpyy(APOIO, BLINE+0, 1), "#"))){
		   
		   comentario = true;
		   
	   }*/
	
	   
	return comentario;
	
	
}


void extrairString(char destino[], char origem[], int inicio, int qtde){
	
	strncpy(destino, origem+inicio, qtde);
	destino[qtde] = '\0';
	
}

void imprimirKSURF(FILE* IW, int &KB){
	fprintf(IW, "KSURF =");
	//printf("\n\nQTREE_.KSURF %d\n\n", QTREE_.KSURF[NXG-1][KB-1]);
	for (int KS = 1; KS <= QTREE_.KSURF[NXG-1][KB-1]; KS++){ 
		//if (QTREE_.KSURF[i-1][KB-1] != 0){
			fprintf(IW, "%5d", QTREE_.KSURF[KS-1][KB-1]);
			if ((KS == 15) || (KS == 30))
				fprintf(IW, "\n       ");
		//}
	}
	fprintf(IW, "\n");
	
}

void imprimirKFLAG(FILE* IW, int &KB){
	fprintf(IW, "KFLAG =");
	for (int i = 1; i <= QTREE_.KSURF[NXG-1][KB-1]; i++){
		if (QTREE_.KFLAG[i-1][KB-1] != 5){
		   	fprintf(IW, "%5d", QTREE_.KFLAG[i-1][KB-1]);
			if (i == 15)
			   	fprintf(IW, "\n       ");
		}
		
	}
	fprintf(IW, "\n");
	
}

void imprimirKBODY(FILE* IW, int &KB){
	fprintf(IW, "KBODY =");
	for (int i = 1; i <= QBODY_.KBODY[NXG-1][KB-1]; i++){
		if (QBODY_.KBODY[i-1][KB-1] != 0){
		   	fprintf(IW, "%5d", QBODY_.KBODY[i-1][KB-1]);
			if (i == 15)
				fprintf(IW, "\n       ");
		}
	}
	fprintf(IW, "\n");
	
}

void imprimirKDGHT(FILE* IW, int &KB){
	fprintf(IW, "KDGHT =");
	for (int i = 1; i <= QTREE_.KDGHT[NXG-1][KB-1]; i++){
		if (QTREE_.KDGHT[i-1][KB-1] != 0){
			fprintf(IW, "%5d", QTREE_.KDGHT[i-1][KB-1]);
			if (i == 15)
			   	fprintf(IW, "\n       ");
		}
	}
	fprintf(IW, "\n");
	
}

void peinit2_(double *EMAX, int &NMATER, FILE *IWR, int *INFO, char (*PMFILE)[100]){
	
	
/*Modulo Penelope.f

  A conven��o usada para nomear as sub-rotinas de simula��o de intera��o �
   a seguir:
   - A primeira letra indica a part�cula (E para el�trons, P para
     p�sitrons, G para f�tons).
   - A segunda e terceira letras denotam o mecanismo de intera��o
     (EL para el�stico, IN para inel�stico, SI para ioniza��o de camada interna,
     BR para bremsstrahlung, AN para aniquila��o, RA para Rayleigh, CO para
      Compton, PH para fotoel�trico e PP para produ��o de pares).
   - A quarta letra (min�scula) indica o modelo te�rico
     usado para descrever as intera��es. Isso serve para distinguir o
     modelo padr�o (denotado pela letra 'a') de modelos alternativos.
   - As rotinas de amostragem aleat�ria t�m nomes de quatro letras. Auxiliar
     rotinas, que realizam c�lculos espec�ficos, t�m nomes mais longos,
     com a quinta e subsequentes letras e / ou n�meros indicando
     o tipo de c�lculo (T para se��o x total, D para diferen�as
     se��o x inicial) ou a��o (W para gravar dados em um arquivo, R para ler
     dados de um arquivo, I para inicializa��o do algoritmo de simula��o).
 
   As presentes sub-rotinas podem imprimir mensagens de aviso e erro em
   a unidade de E / S 26. Esta � a unidade de sa�da padr�o no exemplo principal
   programas PENCYL e PENMAIN.
 
   A subrotina PEMATW conecta os arquivos �s unidades de E / S 3 (entrada) e 7
   (sa�da). No entanto, isso n�o entra em conflito com o programa principal,
   porque PEMATW n�o � chamado durante a simula��o. Esta sub-rotina �
   utilizado apenas pelo programa MATERIAL, para gerar arquivos de dados de materiais.
 
   A subrotina PEINIT conecta os arquivos de defini��o de material � entrada
   unidade 3; esta unidade � fechada antes de retornar ao programa principal


   SUBROTINA PEINIT
   Entrada de dados de materiais e inicializa��o de rotinas de simula��o.
 
   Cada material � definido por meio de um arquivo de entrada, que � criado por
   o programa MATERIAL utilizando informa��es contidas no banco de dados.
   Este arquivo pode ser modificado pelo usu�rio se houver intera��o mais precisa
   os dados est�o dispon�veis.
 
   Argumentos de entrada:
     EMAX ... energia m�xima da part�cula (energia cin�tica para el�trons e
              positrons) usados na simula��o. Nota: Positrons com
              a energia E pode produzir f�tons com energia E + 1.022E6.
     NMATER .... n�mero de materiais na geometria.
     PMFILE .... array de strings de MAXMAT de 20 caracteres. O primeiro NMATER
              elementos s�o os nomes dos arquivos de dados do material.
              O arquivo PMFILE (M) cont�m dados de intera��o de radia��o
              para o material M (ou seja, a ordem � importante!).
     IWR .... unidade de sa�da.
     INFO ... determina a quantidade de informa��es que est�o gravadas em
              o arquivo de sa�da,
                INFO = 1 (ou menos), m�nimo (apenas dados de composi��o).
                INFO = 2, m�dio (mesma informa��o que no material
                  arquivo de dados de defini��o, �til para verificar se a estrutura
                  tura do �ltimo est� correta).
                INFO = 3 ou maior, informa��es completas, incluindo tabelas de
                  propriedades de intera��o usadas na simula��o.
 
   Para os c�lculos preliminares, PEINIT precisa saber a absor��o
   energias EABS (KPAR, M) e os par�metros de simula��o C1 (M), C2 (M),
   WCC (M) e WCR (M). Esta informa��o � introduzida atrav�s do m�dulo
   PENELOPE_mod, que deve ser carregado antes de chamar PEINIT.
 
      USE PENELOPE_mod
      USE TRACK_mod
 
      PRECIS�O DUPLA IMPL�CITA (A-H, O-Z), INTEGER * 4 (I-N)
      CHARACTER * 3 LIT
      CHARACTER * 20 PMFILE
   **** Par�metros de simula��o.
      COMMON / CECUTR / ECUTR (MAXMAT)
      DIMENSION PMFILE (MAXMAT), EABS0 (3, MAXMAT)
 
      COMMON / CSGAWR / ISGAW! Controla mensagens de aviso de SUMGA.
      COMMON / CERSEC / IERSEC
*/


	char LIT[4];
	double EABS0[MAXMAT][3];
	char APOIO[100];
	
	
	//FILE* IWR2 = fopen("material2.dat", "w");
	//IWR = IWR2;
   // if (IWR == NULL){
 	//	printf("N�o foi possivel abrir o arquivo material2.dat");
 	//	exit(0);
 //	}


	fprintf(IWR, "\n **********************************\n");
	fprintf(IWR, " **   PENELOPE  (version 2014)   **\n");
	fprintf(IWR, " **********************************\n");
	
	*CSGAWR_.ISGAW = 0;
	*CERSEC_.IERSEC = 0;
	
	
	*PENELOPE_mod_.NMAT = NMATER;
//	printf("\n\nNMAT: %d\n\n", *PENELOPE_mod_.NMAT);
	for (int M = 1; M <= *PENELOPE_mod_.NMAT; M++){
		if (PENELOPE_mod_.EABS[M-1][1-1] < 49.999e0){
			fprintf(IWR, "EABS(1, %2d) = %.4E eV \n ERROR: electron absorption energy cannot be less than 50 eV\n", M, PENELOPE_mod_.EABS[M-1][1-1]);
			printf("Electron absorption energy less than 50 eV.");
			exit(0);
		}
		EABS0[M-1][1-1] = PENELOPE_mod_.EABS[M-1][1-1];
		
		if (PENELOPE_mod_.EABS[M-1][2-1] < 49.999e0){
			fprintf(IWR, "EABS(2, %2d) = %.4E eV \n ERROR: photon absorption energy cannot be less than 50 eV\n", M, PENELOPE_mod_.EABS[M-1][2-1]);
			printf("Photon absorption energy less than 50 eV.");
			exit(0);
		}
		EABS0[M-1][2-1] = PENELOPE_mod_.EABS[M-1][2-1];
		
		if (PENELOPE_mod_.EABS[M-1][3-1] < 49.999e0){
			fprintf(IWR, "EABS(3, %2d) = %.4E eV \n ERROR: positron absorption energy cannot be less than 50 eV.\n", M, PENELOPE_mod_.EABS[M-1][3-1]);
			printf("Positron absorption energy less than 50 eV.");
			exit(0);
		}
		EABS0[M-1][3-1] = PENELOPE_mod_.EABS[M-1][3-1];	
	}
	
	//Por padr�o, a polariza��o do f�ton n�o � considerada.
	
	*TRACK_mod_.IPOL = 0;
	*TRACK_mod_.SP1 = 0.0;
	*TRACK_mod_.SP2 = 0.0;
	*TRACK_mod_.SP3 = 0.0;
	
	
	//Limite inferior da rede el�trica.
	
	double EMIN = 1.0e35;
	
	for (int M = 1; M <= *PENELOPE_mod_.NMAT; M++){
		if (PENELOPE_mod_.EABS[M-1][1-1] >= *EMAX){
			PENELOPE_mod_.EABS[M-1][1-1] = *EMAX * 0.9999e0;
			fprintf(IWR, "WARNING: EABS(%1d,%2d) has been modified.\n", 1, M);
		}
		
		if (PENELOPE_mod_.EABS[M-1][2-1] >= *EMAX){
			PENELOPE_mod_.EABS[M-1][2-1] = *EMAX * 0.9999e0;
			fprintf(IWR, "WARNING: EABS(%1d,%2d) has been modified.\n", 2, M);
		}
		
		if (PENELOPE_mod_.EABS[M-1][3-1] >= *EMAX){
			PENELOPE_mod_.EABS[M-1][3-1] = *EMAX * 0.9999e0;
			fprintf(IWR, "WARNING: EABS(%1d,%2d) has been modified.\n", 3, M);
		}	
		EMIN = fmin(EMIN, fmin(fmin(PENELOPE_mod_.EABS[M-1][1-1], PENELOPE_mod_.EABS[M-1][2-1]), PENELOPE_mod_.EABS[M-1][3-1]));
	}
			
	
	if (EMIN < 50.0e0)
		EMIN = 50.0e0;

		
	
	fprintf(IWR, "EMIN = %.4E eV,  EMAX = %.4E eV\n", EMIN, *EMAX);
	if (*EMAX < EMIN + 10.0){
		printf("The energy interval is too narrow.");
		exit(0);
	}

	
	
	if (*PENELOPE_mod_.NMAT > MAXMAT){
		fprintf(IWR, "*** PENELOPE cannot handle %2d different materials. \nEdit the source file and change the parameter MAXMAT = %2d to MAXMAT = %2d\n\n", *PENELOPE_mod_.NMAT, MAXMAT, *PENELOPE_mod_.NMAT); 
		printf("PEINIT. Too many materials.");
		exit(0);
	}
	
	
	if (*INFO > 2) 
		fprintf(IWR, "NOTE: 1 mtu = 1 g/cm**2\n");
	
	; //Define a grade de energia da simula��o.
/*	ESIa02_(); //Inicializa rotinas de ioniza��o por impacto de el�trons.
	PSIa02_(); //Inicializa as rotinas de ioniza��o por impacto de p�sitrons.
	GPHa02_(); //Inicializa rotinas fotoel�tricas.
	RELAX02_(); //Inicializa rotinas de relaxamento at�mico.
	RNDG302_(); //Inicializa a rotina de amostragem gaussiana.*/
	

	egrid2_(EMIN, EMAX);
	

	esia02_(); //Inicializa rotinas de ioniza��o por impacto de el�trons.
	psia02_(); //Inicializa as rotinas de ioniza��o por impacto de p�sitrons.
	gpha02_(); //Inicializa rotinas fotoel�tricas.
	//relax0_();
	relax02_(); //Inicializa rotinas de relaxamento at�mico.
	//rndg30_();
	rndg302_(); //Inicializa a rotina de amostragem gaussiana.


	
	for (int M = 1; M <= *PENELOPE_mod_.NMAT; M++){
		
		if (M == 1) strcpy(LIT, "st");
		if (M == 2) strcpy(LIT, "nd");//LIT = "nd";
		if (M == 3) strcpy(LIT, "rd");//LIT = "rd";
		if (M == 4) strcpy(LIT, "th");//LIT = "th";
		
		fprintf(IWR, "\n *********************\n");
		fprintf(IWR, " **	%2d%s material  **\n", M, LIT);
		fprintf(IWR, " *********************\n");
        extrairString(APOIO, PMFILE[M-1], 0 , 20);
		//printf("\nPMFILE %s\n", PMFILE[M-1]);
		
		fprintf(IWR, "\nMaterial data file: %s\n\n", APOIO);
		
		FILE* IRD = fopen(APOIO, "r");
	
		if (IRD == NULL){
			fprintf(IWR, "Nao foi possivel abrir o arquivo %s\n", APOIO);
  	   		printf("Nao foi possivel abrir o arquivo %s\n", APOIO);
		   	exit(0);
		}
		
	//Limites e limiares de energia.
	
		fprintf(IWR, "*** Simulation parameters:\n");
		fprintf(IWR, "	Electron absorption energy = %.4E eV\n", PENELOPE_mod_.EABS[M-1][1-1]);
		fprintf(IWR, "	Photon absorption energy = %.4E eV\n", PENELOPE_mod_.EABS[M-1][2-1]);
		fprintf(IWR, "	Positron absorption energy = %.4E eV\n",PENELOPE_mod_.EABS[M-1][3-1]);
		
		PENELOPE_mod_.C1[M-1] = fmin(0.2e0, fabs(PENELOPE_mod_.C1[M-1]));
		PENELOPE_mod_.C2[M-1] = fmin(0.2e0, fabs(PENELOPE_mod_.C2[M-1]));
		
		PENELOPE_mod_.WCC[M-1] = fmin(fabs(PENELOPE_mod_.WCC[M-1]), *EMAX);
		if (PENELOPE_mod_.WCC[M-1] > PENELOPE_mod_.EABS[M-1][1-1])
			PENELOPE_mod_.WCC[M-1] = PENELOPE_mod_.EABS[M-1][1-1];
		
		if (PENELOPE_mod_.WCR[M-1] > PENELOPE_mod_.EABS[M-1][2-1]) 
			PENELOPE_mod_.WCR[M-1] = PENELOPE_mod_.EABS[M-1][2-1];
		
		if (PENELOPE_mod_.WCR[M-1] < 0.0e0)
			fprintf(IWR, "*** Warning: soft radiative losses are switched off\n");
		fprintf(IWR, "    C1 = %.4E,       C2 = %.4E \n    WCC = %.4E eV,       WCR = %.4E  eV\n", PENELOPE_mod_.C1[M-1], PENELOPE_mod_.C2[M-1],PENELOPE_mod_.WCC[M-1], fmax(PENELOPE_mod_.WCR[M-1],10.0e0  ));
		
		
		CECUTR_.ECUTR[M-1] = fmin(PENELOPE_mod_.EABS[M-1][1-1], PENELOPE_mod_.EABS[M-1][2-1]);
	
//		uintptr_t IRDD = (uintptr_t)IRD;
//		uintptr_t IWRR = (uintptr_t)IWR;
		//pematr_(&M, IRDD, IWRR, INFO);
		

		pematr2_(&M, IRD, IWR, INFO);
		fclose(IRD);
		
	}
	
	//Restaura os valores do usu�rio de EABS (.), Para evitar inconsist�ncias no programa principal.
	for (int M = 1; M <= *PENELOPE_mod_.NMAT; M++){
		PENELOPE_mod_.EABS[M-1][1-1] = EABS0[M-1][1-1];
		PENELOPE_mod_.EABS[M-1][2-1] = EABS0[M-1][2-1];
		PENELOPE_mod_.EABS[M-1][3-1] = EABS0[M-1][3-1];
	}
	
	
}


void egrid2_(double &EMINu, double *EMAXu){ //OK
	
	/*
	Esta sub-rotina define a rede de energia onde as fun��es de transporte s�o
   tabulado. A grade � espa�ada logaritmicamente e assumimos que ela
   � denso o suficiente para permitir a interpola��o log-log linear precisa de
   as fun��es tabuladas.
 
      USE PENELOPE_mod
 
      PRECIS�O DUPLA IMPL�CITA (A-H, O-Z), INTEGER * 4 (I-N)
   **** Rede de energia e constantes de interpola��o para a energia atual.
      COMMON / CEGRID / EMIN, EL, EU, ET (NEGP), DLEMP (NEGP), DLEMP1, DLFC,
     1 XEL, XE, XEK, KE
 
   **** Consist�ncia dos pontos finais do intervalo.
	
	*/
	
	if (EMINu < 50.0e0)
		EMINu = 50.0e0;
	if (EMINu > *EMAXu - 1.0e0){
	//	fprintf(IW, "EMIN = %E.4 eV,  EMAX = %E.4 eV\n");
		printf("EGRID. The energy interval is too narrow.\n");
		exit(0);
	}
	
	//Pontos da grade de energia.
	
	*CEGRID_.EMIN = EMINu;
	*CEGRID_.EL = 0.99999e0 * EMINu;
	*CEGRID_.EU = 1.00001e0 * *EMAXu;
	*CEGRID_.DLFC = log(*CEGRID_.EU / *CEGRID_.EL) / (NEGP - 1);
	*CEGRID_.DLEMP1 = log(*CEGRID_.EL);
	CEGRID_.DLEMP[1-1] = *CEGRID_.DLEMP1;
	CEGRID_.ET[1-1] = *CEGRID_.EL;

	
	for (int I = 2; I <= NEGP; I++){
		CEGRID_.DLEMP[I-1] = CEGRID_.DLEMP[I-1-1] + *CEGRID_.DLFC;
		CEGRID_.ET[I-1] = exp(CEGRID_.DLEMP[I-1]);
	} 
	
	*CEGRID_.DLFC = 1.0e0 / *CEGRID_.DLFC;
	//NOTA: Para determinar o intervalo KE onde a energia E est� localizada, n�s fa�a o seguinte,
	
/*	CEGRID_.XEL = log (TRACK_mod_.E);
	CEGRID_.XE = 1.0e0 + (CEGRID_.XEL - CEGRID_.DLEMP1) * CEGRID_.DFLC;
	CEGRID_.KE = CEGRID_.XE;
	CEGRID_.XEK = CEGRID_.XE - CEGRID_.KE;  //parte 'fracion�ria de XE (usada para interpola��o).  */ 
}


void esia02_(){ //OK
	
	/*
	Esta sub-rotina define todas as vari�veis em comum / CESI0 / para zero.
 
   Deve ser invocado antes de ler a primeira defini��o de material
   Arquivo.
   
   **** Tabelas de se��o transversal de ioniza��o.
	*/
	const static int NRP = 8000;
	for (int I = 1; I<= 99; I++){
		CESI0_.IESIF[I-1] = 0;
		CESI0_.IESIL[I-1] = 0;
		CESI0_.NSESI[I-1] = 0;
		
	}
	
	for (int I = 1; I <= NRP; I++){
		
		for (int J = 1; J <= 16; J++ ){
			CESI0_.XESI[J-1][I-1] = -80.6e0;
		}
	}
	
	*CESI0_.NCURE = 0;
}



void psia02_(){ //OK
	
	/*
	Esta sub-rotina define todas as vari�veis em comum / CPSI0 / para zero.
 
   Deve ser invocado antes de ler a primeira defini��o de material
   Arquivo.
 
   **** Tabelas de se��o transversal de ioniza��o.
	*/
	
	const static int NRP=8000;
	for (int I = 1; I<= 99; I++){
		CPSI0_.IPSIF[I-1] = 0;
		CPSI0_.IPSIL[I-1] = 0;
		CPSI0_.NSPSI[I-1] = 0;
		
	}
	
	for (int I = 1; I <= NRP; I++){
		
		for (int J = 1; J <= 16; J++ ){
			CPSI0_.XPSI[J-1][I-1] = -80.6e0;
		}
	}
	
	*CPSI0_.NCURP = 0;
		
}


void gpha02_(){ //OK
	
	/*
	Esta sub-rotina define todas as vari�veis em comum / CGPH00 / para zero.
   Deve ser invocado antes de ler a primeira defini��o de material Arquivo.
 
      PRECIS�O DUPLA IMPL�CITA (A-H, O-Z), INTEGER * 4 (I-N)
      CHARACTER * 2 LASYMB
   **** Dados do elemento.
      COMMON / CADATA / ATW (99), EPX (99), RSCR (99), ETA (99), EB (99,30),
     1 ALW (99,30), CP0 (99,30), IFI (99,30), IKS (99,30), NSHT (99), LASYMB (99)
   **** Se��es transversais fotoel�tricas.
      PAR�METRO (NTP = 12000)
      COMMON / CGPH00 / EPH (NTP), XPH (NTP, 17), IPHF (99), IPHL (99), NPHS (99), NCUR
	*/
	const static int NTP= 12000;
	
	for (int I = 1; I<= 99; I++){
		CGPH00_.NPHS[I-1] = 0;
		CGPH00_.IPHF[I-1] = 0;
		CGPH00_.IPHL[I-1] = 0;
		
	}
	
	for (int I = 1; I <= NTP; I++){
		CGPH00_.EPH[I-1] = 0.0;
		for (int J = 1; J <= 17; J++ ){
			CGPH00_.XPH[J-1][I-1] = 1.0e-35;
		}
	}
	
	*CGPH00_.NCUR = 0;
		
}


void relax02_(){//OK
	
	/*
	Esta sub-rotina define todas as vari�veis em comum / CRELAX / para zero.
 
      PRECIS�O DUPLA IMPL�CITA (A-H, O-Z), INTEGER * 4 (I-N)
      CHARACTER * 2 LASYMB
   **** Dados do elemento.
      COMMON / CADATA / ATW (99), EPX (99), RSCR (99), ETA (99), EB (99,30),
     1 ALW (99,30), CP0 (99,30), IFI (99,30), IKS (99,30), NSHT (99), LASYMB (99)
   **** Dados de relaxamento at�mico.
      PAR�METRO (NRX = 60.000)
      COMMON / CRELAX / P (NRX), ET (NRX), F (NRX), IS0 (NRX), IS1 (NRX), IS2 (NRX),
     1 IFIRST (99,16), ILAST (99,16), NCUR, KS, MODER
	
	*/
	const static int NRX=60000;
	
	for (int I = 1; I <= 99; I++){
		CADATA_.NSHT[I-1] = 0;
		CADATA_.ETA[I-1] = 0.0e0;;
		
		for (int J = 1; J <= 16; J++){
			CRELAX_.IFIRST[J-1][I-1] = 0;
			CRELAX_.ILAST[J-1][I-1] = 0;
		}
		for (int J = 1; J <= 30; J++){
			CADATA_.IFI[J-1][I-1] = 0;
			CADATA_.EB[J-1][I-1] = 0.0e0;
			CADATA_.IKS[J-1][I-1] = 0;
		}
	}
	
	for (int I = 1; I <= NRX; I++){

		CRELAX_.IS0[I-1] = 0;
		CRELAX_.IS1[I-1] = 0;
		CRELAX_.IS2[I-1] = 0;
		CRELAX_.ET[I-1] = 0;
		CRELAX_.P[I-1] = 0;
		CRELAX_.F[I-1] = 0;	
	}
	
	*CRELAX_.NCUR = 0;
	
}


void rndg302_(){ //OK
	
	//Inicializa��o da fun��o de amostragem RNDG3.

	static const int NR = 128;
	
	int N = NR;
	int NU = N / 4;
	double XMIN = -3.0e0;
	double XMAX = 3.0e0;
	double ERRM;
	FILE *F;
	int PDF = 1;

	
	rita02_(PDF, XMIN, XMAX, N, NU, ERRM, F);

	if (N != NR){
		printf ("RNDG30: initialisation error (1).");
		exit(0);
	}
	
	if (ERRM > 1.0e-6){
		printf ("RNDG30: initialisation error (2).");
		exit(0);	
	}
	
	for (int I = 1; I <= NR; I++){
		CRNDG3_.X[I-1] = CRITAA_.XA[I-1];
		CRNDG3_.A[I-1] = CRITAA_.AA[I-1];
		CRNDG3_.B[I-1] = CRITAA_.BA[I-1];
		CRNDG3_.F[I-1] = CRITAA_.FA[I-1];
		CRNDG3_.KA[I-1] = CRITAA_.IA[I-1];
//		printf("\n %f \n %f \n %f \n %f \n %d\n", CRNDG3_.X[I-1], CRNDG3_.A[I-1], CRNDG3_.B[I-1], CRNDG3_.F[I-1], CRNDG3_.KA[I-1]);
	}
	
	*CRNDG3_.NPM1 = *CRITAA_.NPM1A;
//	printf("\nNPM1: %d\n", *CRNDG3_.NPM1);
	
	
	
}

double rndg3f2_(double X){ //OK
	
	//Distribui��o gaussiana truncada, restrita ao intervalo (-3,0, 3,0), com m�dia zero e vari�ncia unit�ria.
	double resultado;
	if (fabs(X) > 3.00001e0){
		resultado = 0.0e0;
	} else{
		resultado = exp(-X*X*0.5e0/pow(1.01538698e0,2));
	}
	return resultado;
	
}

void rita02_(int &PDF, double &XLOW, double &XHIGH, int &N, int &NU, double &ERRM, FILE *IWR){ //OK
	
	/*
	
	Inicializa��o do algoritmo RITA para amostragem aleat�ria de um
   vari�vel aleat�ria cont�nua X de uma fun��o de distribui��o de probabilidade
   PDF (X) definido no intervalo (XLOW, XHIGH). A fun��o externa
   PDF (X) --n�o necessariamente normalizado-- deve ser fornecido pelo usu�rio.
   N � o n�mero de pontos na grade de amostragem. Esses pontos s�o
   determinado por meio de uma estrat�gia adaptativa que minimiza local
   erros de interpola��o. Os primeiros pontos da grade NU s�o uniformemente espa�ados
   em (XLOW, XHIGH); quando NU � negativo, a grade inicial consiste em
   -NU pontos espa�ados logaritmicamente (neste caso, XLOW deve ser
   n�o negativo).
 
   ERRM � uma medida do erro de interpola��o (o maior valor de
   o erro absoluto da interpola��o racional integrada sobre cada
   intervalo de grade).
 
   **** Coeficientes de interpola��o e tabelas PDF s�o impressos em
         arquivos separados (UNIT = IWR) se IWR for maior que zero.
 
   Outros subprogramas necess�rios: fun��o EXTERNA PDF,
                             sub-rotinas RITAI0 e IRND0.
	
	*/
//	printf("\n\nPDF: %d\nXLOW: %f, XHIGH: %f\n",PDF, XLOW, XHIGH );
	const static int NM = 512;
	
	
	ritai02_(PDF, XLOW, XHIGH, N, NU, ERRM, IWR);


	
	*CRITAA_.NPM1A = *CRITA_.NPM1;
	
	for (int I = 1; I <= *CRITA_.NPM1; I++){
		CRITAA_.XA[I-1] = CRITA_.XT[I-1];
		CRITAA_.AA[I-1] = CRITA_.A[I-1];
		CRITAA_.BA[I-1] = CRITA_.B[I-1];	
	}
	
	CRITAA_.XA[*CRITA_.NPM1+1-1] = CRITA_.XT[*CRITA_.NPM1+1-1];
	double SAVE = *CRITAN_.CNORM;
	irnd02_(CRITA_.DPAC, CRITAA_.FA, CRITAA_.IA, CRITAA_.NPM1A);
	*CRITAN_.CNORM = SAVE;
	CRITAA_.FA[*CRITA_.NPM1+1-1] = 1.0e0;
	CRITAA_.IA[*CRITA_.NPM1+1-1] = *CRITAA_.NPM1A;
	
/*	for (int I = 1; I <= *CRITA_.NPM1; I++){
		printf("\nCRITAA\n");
		printf("XA: %f, AA: %f, BA: %f,\n FA: %f, IA: %d, NPM1A: %d\n", CRITAA_.XA[I-1], CRITAA_.AA[I-1], CRITAA_.BA[I-1], CRITAA_.FA[I-1], CRITAA_.IA[I-1], *CRITAA_.NPM1A);
	}*/
	
}


void pematr2_(int *M, FILE *IRD, FILE *IWR, int *INFO){	
	
	/*
	
	Esta sub-rotina l� o arquivo de defini��o do material M (unidade IRD) 
	e inicializa as rotinas de simula��o para este material. 
	A informa��o � escrita na unidade IWR, a quantidade de informa��o escrita 
	� determinada pelo valor de INFO.
	
	*/

    static const int NO = 512;
	static const int NRP = 8000;
	static const int NOCO=512;
	static const int NDIM=12000;
	static const int NEGP = 200;
	
	//arrays auxiliares
/*	double STFI[NO];
	double STFO[NO];
	double INOUT[NO];
	double EIT[NEGP];
	double EITL[NEGP];
	double FL[NEGP];
	double F1[NEGP];
	double F2[NEGP];
    double F3[NEGP];
	double F4[NEGP];
	double RADY[NEGP];
	double RADN[NEGP];*/
	


    
	
	
	double *STFI = (double *) malloc(NO * sizeof(double));
	double *STFO = (double *) malloc(NO * sizeof(double));
	double *INOUT = (double *) malloc(NO * sizeof(double));
	double *EIT = (double *) malloc(NEGP * sizeof(double));
	double *EITL = (double *) malloc(NEGP * sizeof(double));
	double *FL = (double *) malloc(NEGP * sizeof(double));
	double *F1 = (double *) malloc(NEGP * sizeof(double));
	double *F2 = (double *) malloc(NEGP * sizeof(double));
	double *F3 = (double *) malloc(NEGP * sizeof(double));
	double *F4 = (double *) malloc(NEGP * sizeof(double));
	double *RADY = (double *) malloc(NEGP * sizeof(double));
	double *RADN = (double *) malloc(NEGP * sizeof(double));
	
	double WCCM;
	double DIFMAX;
	char NAME[100];
	int IZZ;
	int ISH;
	double OMEGA;
	double EXPT;
	double WCRM;
	double IBREMS;
	double XH0, XH1, XH2, XS0, XS1, XS2, XT1, XT2, DELTA;
	double H0, H1, H2, S0, S1, S2, R0, R1, R2;
	double XS1SI, XS2SI, XT1SI, XT2SI, STPSI, STRSI, XS1IN, XS2IN, XT1IN, XT2IN, STPIN, STRIN;
	double STPI, STPC;
	double EC, TXAN;
	int J, INDC, ICOR, NP; 
	double EBIN, EE, FACT, PCSI, STOT, DFERMI, RNFO, SSIT, FNORM, SINT, STPR, XL, XU, DR;
	double FSOFTE, FSOFTP, W1EW, W2EW, W1PW, W2PW, FMTU;
	double FPEL, FPSI, FPIN, FPBR, FPAN, FPTOT, VMOLL, TRIPL, PAIR, ECS, PCO, PRA, PPP, PT, PPH;
	
	/*
	Reescalonamento de MFPs calculados.
  Quando ISCALE � definido igual a 1, o programa redimensiona o valor calculado
  se��es transversais totais e MFPs para reproduzir as se��es transversais e
  poderes de parada lidos do arquivo de dados do material de entrada. Para evitar isso
  redimensionando, defina ISCALE = 0.
	*/
	
	int ISCALE = 1;
	
	// Caracteristicas dos Materiais
	
	if (*M > MAXMAT){
		printf("PEMATR. Too many materials.");
		exit(0);
	}
	
	char LNAME[58] = " PENELOPE (v. 2014)  Material data file ...............";
	
	fgets(LINHA, sizeof(LINHA), IRD);
	
	extrairString(NAME, LINHA, 0, 55);
	if (strcmp(NAME, LNAME)){
		fprintf(IWR, "I/O error. Corrupt material data file.\n");
		fprintf(IWR, "     The first line is: %s\n", NAME);
		fprintf(IWR, "     ... and should be: %s\n", LNAME);
		printf("I/O error. Corrupt material data file.\n");
		printf("     The first line is: %s\n\n", NAME);
		printf("     ... and should be: %s\n", LNAME);
		printf("PEMATR. Corrupt material data file\n");
		exit(0);
	}
	
	fprintf(IWR, "\n%s\n", NAME);
	
	
	fgets(LINHA, sizeof(LINHA), IRD);
	extrairString(LNAME, LINHA, 10, 62);
	
	fprintf(IWR, " Material: %s\n", LNAME);
	
	fgets(LINHA, sizeof(LINHA), IRD);
	extrairString(APOIO, LINHA, 15, 15);
	COMPOS_.RHO[*M-1] = atof(APOIO);
	fprintf(IWR, " Mass density: %.8E g/cm^3\n", COMPOS_.RHO[*M-1]);
	
	fgets(LINHA, sizeof(LINHA), IRD);
	extrairString(APOIO, LINHA, 37, 3);
	COMPOS_.NELEM[*M-1] = atoi(APOIO);
	fprintf(IWR, " Number of elements in the molecule = %2d\n", COMPOS_.NELEM[*M-1]);
	
	if (COMPOS_.NELEM[*M-1] > 30){
		printf("PEMATR. Too many elements.\n");
		exit(0);
	}
	
	
	COMPOS_.ZT[*M-1] = 0.0e0;
	COMPOS_.AT[*M-1] = 0.0e0;
	
	for (int I = 1; I <= COMPOS_.NELEM[*M-1]; I++){
		fgets(LINHA, sizeof(LINHA), IRD);
		extrairString(APOIO, LINHA, 18, 3);
		COMPOS_.IZ[I-1][*M-1] = atoi(APOIO);
	
	
		extrairString(APOIO, LINHA, 40, 15);
		COMPOS_.STF[I-1][*M-1] = atof(APOIO);
		char teste[3];
		strcpy(teste,CADATA_.LASYMB[COMPOS_.IZ[I-1][*M-1]-1]); 
		teste[2] = '\0';
		
	//	printf("\n\nELEMENTO %s\n\n", CADATA_.LASYMB[COMPOS_.IZ[I-1][*M-1]-1]);
		
		fprintf(IWR, " Element: %s (Z=%2d), atoms/molecule = %.8E\n", teste, COMPOS_.IZ[I-1][*M-1], COMPOS_.STF[I-1][*M-1]);
		COMPOS_.ZT[*M-1] = COMPOS_.ZT[*M-1] + COMPOS_.STF[I-1][*M-1] * COMPOS_.IZ[I-1][*M-1];
		IZZ = COMPOS_.IZ[I-1][*M-1];
		COMPOS_.AT[*M-1] = COMPOS_.AT[*M-1] + CADATA_.ATW[IZZ-1] * COMPOS_.STF[I-1][*M-1];
	}
	
	COMPOS_.VMOL[*M-1] = AVOG*COMPOS_.RHO[*M-1] / COMPOS_.AT[*M-1];
	if (*INFO >= 2){
		fprintf(IWR, "\n Molecular density = %.8E\n", COMPOS_.VMOL[*M-1]);	
	}
	CEIN_.OP2[*M-1]= FOURPI * COMPOS_.ZT[*M-1] * COMPOS_.VMOL[*M-1] * pow(A0B,3) * pow(HREV,2);
	OMEGA=sqrt(CEIN_.OP2[*M-1]);
	
	PENELOPE_mod_.DEN[*M-1] = COMPOS_.RHO[*M-1];
	PENELOPE_mod_.RDEN[*M-1] = 1.0e0 / PENELOPE_mod_.DEN[*M-1];
	
	if (*INFO >= 2){
		fprintf(IWR, "\n *** Electron/positron inelastic scattering.\n");
		fprintf(IWR," Plasma energy = %.8E eV\n", OMEGA );
	}
	
	fgets(LINHA, sizeof(LINHA), IRD);
	extrairString(APOIO, LINHA, 25, 15);
	CEIN_.EXPOT[*M-1] = atof(APOIO);
	fprintf(IWR, " Mean excitation energy = %.8E eV\n", CEIN_.EXPOT[*M-1]);
	
	//Colisoes Inelasticas E/P
	
	fgets(LINHA, sizeof(LINHA), IRD);
	extrairString(APOIO, LINHA, 24, 3);
	CEIN_.NOSC[*M-1] = atoi(APOIO);
		
	if ((*INFO >= 2 ) || (CEIN_.NOSC[*M-1] > NO)) 
		fprintf(IWR," Number of oscillators = %4d\n",CEIN_.NOSC[*M-1] );
	if (CEIN_.NOSC[*M-1] > NO) {
		printf("PEMATR. Too many oscillators.\n");
		exit(0);
	}
	if (*INFO >= 2 ){
		fprintf(IWR, "\n           Fi            Ui (eV)         Wi (eV)       KZ  KS\n");
	    fprintf(IWR, "---------------------------------------------------------------\n");
	}
	
	EXPT = 0.0;
	
	for ( int I = 1; I <= CEIN_.NOSC[*M-1]; I++){
		fgets(LINHA, sizeof(LINHA), IRD);
		extrairString(APOIO, LINHA, 4, 16);
		CEIN_.F[I-1][*M-1] = atof(APOIO);
		extrairString(APOIO, LINHA, 20, 16);
		CEIN_.UI[I-1][*M-1] = atof(APOIO);
		extrairString(APOIO, LINHA,36, 16);
		CEIN_.WRI[I-1][*M-1] = atof(APOIO);
		extrairString(APOIO, LINHA, 52, 4);
		CEIN_.KZ[I-1][*M-1] = atoi(APOIO);
		extrairString(APOIO, LINHA, 56, 4);
		CEIN_.KS[I-1][*M-1] = atoi(APOIO);
		
		if (CEIN_.UI[I-1][*M-1] < 1.0e-3){
			CEIN_.UI[I-1][*M-1] = 0.0e0;	
		}
		
		if (*INFO >= 2) {
			fprintf(IWR, "%4d %.8E %.8E %.8E %4d %4d\n", I, CEIN_.F[I-1][*M-1], CEIN_.UI[I-1][*M-1], CEIN_.WRI[I-1][*M-1], CEIN_.KZ[I-1][*M-1], CEIN_.KS[I-1][*M-1] );
		}
		
		EXPT = EXPT + CEIN_.F[I-1][*M-1] * log(CEIN_.WRI[I-1][*M-1]);

	}
	
	EXPT = exp(EXPT / COMPOS_.ZT[*M-1]);
	
	
	if (fabs(EXPT - CEIN_.EXPOT[*M-1]) > 1.0e-6*CEIN_.EXPOT[*M-1]){
		printf("EXPT      = %f\n", EXPT);
		printf("EXPOT (M) = %f\n", CEIN_.EXPOT[*M-1]);
		printf("Inconsistent oscillator data.");
		exit(0);	
	}
		
	//Espalhamento Compton
	
	fgets(LINHA, sizeof(LINHA), IRD);
	extrairString(APOIO, LINHA, 19, 3);
	CGCO_.NOSCCO[*M-1] = atoi(APOIO);
	
	if ((*INFO >= 2 ) || (CGCO_.NOSCCO[*M-1] > NOCO)) {
		fprintf(IWR, "\n *** Compton scattering (Impulse Approximation)\n");
		fprintf(IWR," Number of shells = %4d\n",CGCO_.NOSCCO[*M-1] );
	
		if (CGCO_.NOSCCO[*M-1] > NOCO) {
			printf("PEMATR. Too many shells.\n");
			exit(0);
		}
	}
	
	if (*INFO >= 2 ){
		fprintf(IWR, "\n           Fi            Ui (eV)         Ji (0)       KZ  KS\n");
    	fprintf(IWR, "--------------------------------------------------------------\n");
	}
	
	for ( int I = 1; I <= CGCO_.NOSCCO[*M-1]; I++){
		fgets(LINHA, sizeof(LINHA), IRD);
		extrairString(APOIO, LINHA, 4, 16);
		CGCO_.FCO[I-1][*M-1] = atof(APOIO);
		extrairString(APOIO, LINHA, 20, 16);
		CGCO_.UICO[I-1][*M-1] = atof(APOIO);
		extrairString(APOIO, LINHA,36, 16);
		CGCO_.FJ0[I-1][*M-1] = atof(APOIO);
		extrairString(APOIO, LINHA, 52, 4);
		CGCO_.KZCO[I-1][*M-1] = atoi(APOIO);
		extrairString(APOIO, LINHA, 56, 4);
		CGCO_.KSCO[I-1][*M-1] = atoi(APOIO);
		
		if (*INFO >= 2) {
			fprintf(IWR, "%4d %.8E %.8E %.8E %4d %4d\n", I, CGCO_.FCO[I-1][*M-1], CGCO_.UICO[I-1][*M-1], CGCO_.FJ0[I-1][*M-1], CGCO_.KZCO[I-1][*M-1], CGCO_.KSCO[I-1][*M-1] );
		}
		
		CGCO_.FJ0[I-1][*M-1] = CGCO_.FJ0[I-1][*M-1] * SL;

	}
	
	//Dados do relaxamento atomico
	
	for (int I = 1; I <= COMPOS_.NELEM[*M-1]; I++){
         relaxr2_(IRD, IWR, INFO);
	}	
	
	
	
	//Ioniza��o da camada interna por impacto de el�trons e p�sitrons.
		
	esiar2_(M, IRD, IWR, INFO);
	
	psiar2_(M, IRD, IWR, INFO); 
	
  

	//propriedades de intera��o do eletron e positron
	
	//Se��o x de bremsstrahlung em escala.
	
	if (PENELOPE_mod_.WCR[*M-1] >= 0.0e0){
		PENELOPE_mod_.WCR[*M-1] = fmax(PENELOPE_mod_.WCR[*M-1], 10.0e0);
		WCRM = PENELOPE_mod_.WCR[*M-1];
		IBREMS = 1;
	} else{
		WCRM = 10.0e0;
		PENELOPE_mod_.WCR[*M-1] = 10.0e0;
		IBREMS = 0;
	}

	ebrar2_(&WCRM, M, IRD, IWR, INFO);
	
	
	
	//Distribui��o angular de Bremsstrahlung.
	
	braar2_(M, IRD, IWR, INFO);
	
	
	
	//Parando poderes.
	
		
	
	fgets(LINHA, sizeof(LINHA), IRD);
	extrairString(APOIO, LINHA, 58, 4);
	int NDATA = atoi(APOIO);
	
	
	if (*INFO >= 2){
		fprintf(IWR, "\n*** Stopping powers for electrons and positrons,  NDATA =%4d\n", NDATA);
	}
	if (NDATA > NEGP){
		printf("PEMATR. Too many data points (1)\n");
		exit(0);
	}
	
	if (*INFO >= 2){
		fprintf(IWR, "\n  Energy      Scol,e-      Srad,e-      Scol,e+      Srad,e+      \n    (eV)      (MeV/mtu)    (MeV/mtu)    (MeV/mtu)    (MeV/mtu)\n\n");
	    fprintf(IWR, " --------------------------------------------------------------\n");
	}
	
	for (int I = 1; I <= NDATA; I++){
		fgets(LINHA, sizeof(LINHA), IRD);
		extrairString(APOIO, LINHA, 0, 10);
		EIT[I-1] = atof(APOIO);
		extrairString(APOIO, LINHA, 10, 12);
		F1[I-1] = atof(APOIO);
		extrairString(APOIO, LINHA, 22, 12);
		F2[I-1] = atof(APOIO);
		extrairString(APOIO, LINHA, 34, 12);
		F3[I-1] = atof(APOIO);
		extrairString(APOIO, LINHA, 46, 12);
		F4[I-1] = atof(APOIO);	
		if (*INFO >= 2){
			fprintf(IWR, " %.3E %.5E %.5E %.5E %.5E\n", EIT[I-1], F1[I-1], F2[I-1], F3[I-1], F4[I-1] );
		}
		EITL[I-1] = log(EIT[I-1]);
		
	}
	
//	  exit(0);
	
	//Se��o x inel�stica para el�trons.
	
	WCCM = PENELOPE_mod_.WCC[*M-1];
	DIFMAX = 0.0e0;

	for (int I = 1; I <= NDATA; I++){
		if ((EIT[I-1] >= *CEGRID_.EL) && (EIT[I-1] <= *CEGRID_.EU)){
	        STPI = F1[I-1] * COMPOS_.RHO[*M-1] * 1.0e6; //! Collision stopping power.
	       	
		    einat2_(EIT[I-1], WCCM, XH0, XH1, XH2, XS0, XS1, XS2, XT1, XT2, DELTA, M);
		  //  einat_(EIT[I-1], WCCM, XH0, XH1, XH2, XS0, XS1, XS2, XT1, XT2, DELTA, M);
            
		    STPC = (XS1+XH1)* COMPOS_.VMOL[*M-1];
		    DIFMAX = fmax(DIFMAX, 100.0e0 * fabs(STPC - STPI) / STPI);	
		}		
	}
	
	if ((DIFMAX > 1.0e-3) && (ISCALE == 1)){
		ICOR = 1;
		for (int I = 1; I <= NDATA; I++){
			FL[I-1] = log(F1[I-1] * COMPOS_.RHO[*M-1] * 1.0e6);	
		}
		spline2_(EITL, FL, CEEL00_.A, CEEL00_.B, CEEL00_.C, CEEL00_.D, 0.0, 0.0, NDATA);
		
	} else{
		ICOR = 0;
	}
	
	//Classifique as cascas internas e externas.
	
	int NI = 0;
	int NS = 0;
	
	for (int IO = 1; IO <= NO; IO++){
		STFI[IO-1] = 0.0e0;
        STFO[IO-1] = 0.0e0;
	}
	
	for (int KO = 1; KO <= CEIN_.NOSC[*M-1]; KO++){
    	IZZ = CEIN_.KZ[KO-1][*M-1];
        ISH = CEIN_.KS[KO-1][*M-1];
        if ((IZZ < 3) || (ISH > 16) || (CEIN_.UI[KO-1][*M-1] < 50.0e0)){
        	NI=NI+1;
            CEINAC_.IEIN[NI-1][*M-1] = KO;
            CESIN_.ISIE[KO-1] = -NI;
            INOUT[NI-1] = 0;
            for (int IEL = 1; IEL <= COMPOS_.NELEM[*M-1]; IEL++){
            	if (IZZ == COMPOS_.IZ[IEL-1][*M-1])
            		STFO[NI-1] = COMPOS_.STF[IEL-1][*M-1];
			}   	
		}else if(ISH <= CESI0_.NSESI[IZZ-1]){
			EBIN = CADATA_.EB[ISH-1][IZZ-1];
			if (EBIN > CECUTR_.ECUTR[*M-1]){
				NS = NS + 1;
            	CESIAC_.IESI[NS-1][*M-1] = KO;
            	CESIN_.ISIE[KO-1] = NS;
            	for (int IEL = 1; IEL <= COMPOS_.NELEM[*M-1]; IEL++){
    		   		if (IZZ == COMPOS_.IZ[IEL-1][*M-1])
            		   	STFI[NS-1] = COMPOS_.STF[IEL-1][*M-1];
				}   		
			} else{
				NI=NI+1;
            	CEINAC_.IEIN[NI-1][*M-1] = KO;
            	CESIN_.ISIE[KO-1] = -NI;
            	INOUT[NI-1] = 1;
            	for (int IEL = 1; IEL <= COMPOS_.NELEM[*M-1]; IEL++){
            		if (IZZ == COMPOS_.IZ[IEL-1][*M-1])
            		   	STFO[NI-1] = COMPOS_.STF[IEL-1][*M-1];
				}   	
					
			}
				
		} else{
			NI=NI+1;
			CEINAC_.IEIN[NI-1][*M-1] = KO;
   			CESIN_.ISIE[KO-1] = -NI;
   			INOUT[NI-1] = 0;
            for (int IEL = 1; IEL <= COMPOS_.NELEM[*M-1]; IEL++){
            	if (IZZ == COMPOS_.IZ[IEL-1][*M-1])
            	   	STFO[NI-1] = COMPOS_.STF[IEL-1][*M-1];
			}   		
		}	
	}
	
	CEINAC_.NEIN[*M-1] = NI;
	CESIAC_.NESI[*M-1] = NS;
	
	//Arrays de Simula��o.
	
	for (int I = 1; I <= NEGP; I++){
		EE = CEGRID_.ET[I-1];
		
		einat2_(EE, WCCM, XH0, XH1, XH2, XS0, XS1, XS2, XT1, XT2, DELTA, M);
		//einat_(EE, WCCM, XH0, XH1, XH2, XS0, XS1, XS2, XT1, XT2, DELTA, M);
		
		STPC = (XS1+XH1)* COMPOS_.VMOL[*M-1];
		if (ICOR == 1){
			EC = CEGRID_.DLEMP[I-1];
		
			findi2_(EITL, &EC, &NDATA, &J);
			
			STPI = exp(CEEL00_.A[J-1] + EC*( CEEL00_.B[J-1] + EC * (CEEL00_.C[J-1]+ EC * CEEL00_.D[J-1])));
			FACT = STPI / STPC;
		} else{
			FACT = 1.0e0;
		}
		
		CEIMFP_.DEL[I-1][*M-1] = DELTA;
        CEIMFP_.CSTPE[I-1][*M-1] = STPC * FACT;
        XS1SI = 0.0e0;
        XS2SI = 0.0e0;
        XT1SI = 0.0e0;
        XT2SI = 0.0e0;
        STPSI = 0.0e0;
        STRSI = 0.0e0;
        XS1IN = 0.0e0;
        XS2IN = 0.0e0;
        XT1IN = 0.0e0;
        XT2IN = 0.0e0;
        STPIN = 0.0e0;
        STRIN = 0.0e0;
        
        for( int IO = 1; IO <= CEINAC_.NEIN[*M-1]; IO++){
        	CESIN_.XSEIN[IO-1][I-1] = 0.0e0;
		}
		
		for( int IO = 1; IO <= CESIAC_.NESI[*M-1]; IO++){
        	CESIN_.XSESI[IO-1][I-1] = 0.0e0;
		}
		
		for (int KO = 1; KO <= CEIN_.NOSC[*M-1]; KO++){
			int IO = CESIN_.ISIE[KO-1];
			//Cascas internas
			
			if (IO > 0){
				IZZ = CEIN_.KZ[KO-1][*M-1];
        		ISH = CEIN_.KS[KO-1][*M-1];
        		EBIN = CADATA_.EB[ISH-1][IZZ-1];
        		
        		if (EBIN < EE){
        			INDC = CESI0_.IESIF[IZZ-1] - 1;
  	  	        	PCSI = exp(CESI0_.XESI[ISH-1][INDC+I-1]) * STFI[IO-1];
  	  	        	if (PCSI > 1.0e-35){
  	  	        		double AUX = 0.0;
  	  	        		einat12_(EE, CEIN_.UI[KO-1][*M-1], CEIN_.WRI[KO-1][*M-1], AUX, WCCM, H0, H1, H2, S0, S1, S2, R0, R1, R2);
						
						
						STOT = CEIN_.F[KO-1][*M-1] * (S0 + H0);	
						
						if (STOT > 1.0e-35){
							DFERMI = (CEIN00_.SES0[KO-1] + CEIN00_.SEH0[KO-1]) / STOT;
                  	    	PCSI = PCSI * DFERMI;
                    		CESIN_.XSESI[IO-1][I-1] = PCSI;
                    		RNFO = PCSI / (CEIN00_.SES0[KO-1] + CEIN00_.SEH0[KO-1]);
                    		STPSI = STPSI +(CEIN00_.SES1[KO-1] + CEIN00_.SEH1[KO-1]) * RNFO;
                    		STRSI = STRSI +(CEIN00_.SES2[KO-1] + CEIN00_.SEH2[KO-1]) * RNFO;	
						} else{
							CESIN_.XSESI[IO-1][I-1] = PCSI;
                    		STPSI = STPSI + PCSI * EBIN;
                    		STRSI = STRSI + PCSI * pow(EBIN,2);
							
						}
						
					}	
				}	
			}else{
				//outras cascas
				IO = -IO;
				
				if (INOUT[IO-1] == 1){
					IZZ = CEIN_.KZ[KO-1][*M-1];
        			ISH = CEIN_.KS[KO-1][*M-1];
        			EBIN = CADATA_.EB[ISH-1][IZZ-1];
        			INDC = CESI0_.IESIF[IZZ-1] - 1;
  	  	        	PCSI = exp(CESI0_.XESI[ISH-1][INDC+I-1]) * STFO[IO-1];
  	  	        	double AUX = 0.0;
  	  	        	
  	       			einat12_(EE, CEIN_.UI[KO-1][*M-1], CEIN_.WRI[KO-1][*M-1], AUX, WCCM, H0, H1, H2, S0, S1, S2, R0, R1, R2);
                    STOT = CEIN_.F[KO-1][*M-1] * (S0 + H0);
                    if (STOT > 1.0e-35){
						DFERMI = (CEIN00_.SES0[KO-1] + CEIN00_.SEH0[KO-1]) / STOT;
                  	   	PCSI = PCSI * DFERMI;
                    	RNFO = PCSI / (CEIN00_.SES0[KO-1] + CEIN00_.SEH0[KO-1]);
						XS1SI = XS1SI + CEIN00_.SES1[KO-1] * RNFO;
                		XS2SI = XS2SI + CEIN00_.SES2[KO-1] * RNFO;
                		XT1SI = XT1SI + CEIN00_.SET1[KO-1] * RNFO;
                		XT2SI = XT2SI + CEIN00_.SET2[KO-1] * RNFO;	
						CESIN_.XSEIN[IO-1][I-1] = CEIN00_.SEH0[KO-1] * RNFO;
                    	STPSI = STPSI +(CEIN00_.SES1[KO-1] + CEIN00_.SEH1[KO-1]) * RNFO;
                    	STRSI = STRSI +(CEIN00_.SES2[KO-1] + CEIN00_.SEH2[KO-1]) * RNFO;	
					} else{
						XS1SI = XS1SI + CEIN00_.SES1[KO-1];
                		XS2SI = XS2SI + CEIN00_.SES2[KO-1];
                		XT1SI = XT1SI + CEIN00_.SET1[KO-1];
                		XT2SI = XT2SI + CEIN00_.SET2[KO-1];	
						CESIN_.XSEIN[IO-1][I-1] = CEIN00_.SEH0[KO-1];
                    	STPSI = STPSI +(CEIN00_.SES1[KO-1] + CEIN00_.SEH1[KO-1]);
                    	STRSI = STRSI +(CEIN00_.SES2[KO-1] + CEIN00_.SEH2[KO-1]);	
					}
					
				} else{
					XS1IN = XS1IN + CEIN00_.SES1[KO-1];
                	XS2IN = XS2IN + CEIN00_.SES2[KO-1];
                	XT1IN = XT1IN + CEIN00_.SET1[KO-1];
                	XT2IN = XT2IN + CEIN00_.SET2[KO-1];
                	CESIN_.XSEIN[IO-1][I-1] = CEIN00_.SEH0[KO-1];
                	STPIN = STPIN + (CEIN00_.SES1[KO-1] + CEIN00_.SEH1[KO-1]);
                	STRIN = STRIN + (CEIN00_.SES2[KO-1] + CEIN00_.SEH2[KO-1]);
				}
				
			}	
		}
		
		SSIT = 0.0;
		if (CESIAC_.NESI[*M-1] > 0){
			for (int IO =1; IO <= CESIAC_.NESI[*M-1]; IO++){
				CESIAC_.ESIAC[IO-1][I-1][*M-1] = SSIT;
				SSIT = SSIT + CESIN_.XSESI[IO-1][I-1];	
			}
		    CEIMFP_.SEISI[I-1][*M-1] = log(fmax(SSIT * COMPOS_.VMOL[*M-1], 1.0e-35));	
		} else{
			CEIMFP_.SEISI[I-1][*M-1] = log(1.0e-35);	
		}

		if (SSIT > 1.0e-35){
			for (int IO = 1; IO <= CESIAC_.NESI[*M-1]; IO++){
				CESIAC_.ESIAC[IO-1][I-1][*M-1] = CESIAC_.ESIAC[IO-1][I-1][*M-1] / SSIT;	
			}
		} else{
			for (int IO = 1; IO <= CESIAC_.NESI[*M-1]; IO++){
				CESIAC_.ESIAC[IO-1][I-1][*M-1] = 1.0e0;	
			}
		}
		
		STPSI = STPSI * COMPOS_.VMOL[*M-1];
        STPIN = STPIN * COMPOS_.VMOL[*M-1];
        FNORM = (CEIMFP_.CSTPE[I-1][*M-1]  - STPSI) / STPIN;
        
        SINT = 0.0e0;
        
        
        if (CEINAC_.NEIN[*M-1] == 0){
			printf("PHMATR. No outer shells?\n");
			exit(0);
		}
		
		for (int IO = 1; IO <= CEINAC_.NEIN[*M-1]; IO++){
			if (INOUT[IO-1] == 0)
				CESIN_.XSEIN[IO-1][I-1] = FNORM * CESIN_.XSEIN[IO-1][I-1];
			CEINAC_.EINAC[IO-1][I-1][*M-1] = SINT;
			SINT  = SINT + CESIN_.XSEIN[IO-1][I-1];	
		}
		
		if (SINT > 1.0e-35){
			for (int IO = 1; IO <= CEINAC_.NEIN[*M-1]; IO++){
				CEINAC_.EINAC[IO-1][I-1][*M-1] = CEINAC_.EINAC[IO-1][I-1][*M-1] / SINT;
			}
		} else{
			for (int IO = 1; IO <= CEINAC_.NEIN[*M-1]; IO++){
				CEINAC_.EINAC[IO-1][I-1][*M-1] = 1.0e0;
			}
		}
		CEIMFP_.SEHIN[I-1][*M-1] = log(fmax(SINT * COMPOS_.VMOL[*M-1], 1.0e-35));
        CLAS1E_.TSTRE[I-1][*M-1] = (STRSI + FNORM * STRIN) * COMPOS_.VMOL[*M-1];
        CEIMFP_.W1E[I-1][*M-1] = (XS1SI + FNORM * XS1IN) * COMPOS_.VMOL[*M-1];
        CEIMFP_.W2E[I-1][*M-1] = (XS2SI + FNORM * XS2IN) * COMPOS_.VMOL[*M-1];
        CEINTF_.T1EI[I-1] = (XT1SI + FNORM * XT1IN) * COMPOS_.VMOL[*M-1];
        CEINTF_.T2EI[I-1] = (XT2SI + FNORM * XT2IN) * COMPOS_.VMOL[*M-1];	
	}
	
	
	//Se��o x inel�stica para p�sitrons.
	
		
	DIFMAX = 0.0e0;

	for (int I = 1; I <= NDATA; I++){
		if ((EIT[I-1] >= *CEGRID_.EL) && (EIT[I-1] <= *CEGRID_.EU)){
	        STPI = F3[I-1] * COMPOS_.RHO[*M-1] * 1.0e6; //! Collision stopping power.
	       	
		    pinat2_(EIT[I-1], WCCM, XH0, XH1, XH2, XS0, XS1, XS2, XT1, XT2, DELTA, M);

		    STPC = (XS1+XH1)* COMPOS_.VMOL[*M-1];
		    DIFMAX = fmax(DIFMAX, 100.0e0 * fabs(STPC - STPI) / STPI);	
		}		
	}
	

	
	if ((DIFMAX > 1.0e-3) && (ISCALE == 1)){
		ICOR = 1;
		for (int I = 1; I <= NDATA; I++){
			FL[I-1] = log(F3[I-1] * COMPOS_.RHO[*M-1] * 1.0e6);	
		}
		spline2_(EITL, FL, CEEL00_.A, CEEL00_.B, CEEL00_.C, CEEL00_.D, 0.0, 0.0, NDATA);
	} else{
		ICOR = 0;
	}
	
	//Classifique as cascas internas e externas.
	
	NI = 0;
	NS = 0;
	
	for (int IO = 1; IO <= NO; IO++){
		STFI[IO-1] = 0.0e0;
        STFO[IO-1] = 0.0e0;
	}
	
	for (int KO = 1; KO <= CEIN_.NOSC[*M-1]; KO++){
    	IZZ = CEIN_.KZ[KO-1][*M-1];
        ISH = CEIN_.KS[KO-1][*M-1];
        if ((IZZ < 3) || (ISH > 16) || (CEIN_.UI[KO-1][*M-1] < 50.0e0)){
        	NI=NI+1;
            CPINAC_.IPIN[NI-1][*M-1] = KO;
            CPSIN_.ISIP[KO-1] = -NI;
	//		printf("\n\nCPSIN1 %d\n\n", CPSIN_.ISIP[KO-1]);
            INOUT[NI-1] = 0;
            for (int IEL = 1; IEL <= COMPOS_.NELEM[*M-1]; IEL++){
            	if (IZZ == COMPOS_.IZ[IEL-1][*M-1])
            		STFO[NI-1] = COMPOS_.STF[IEL-1][*M-1];
			}   	
		}else if(ISH <= CPSI0_.NSPSI[IZZ-1]){
			EBIN = CADATA_.EB[ISH-1][IZZ-1];
			if (EBIN > CECUTR_.ECUTR[*M-1]){
				NS = NS + 1;
            	CPSIAC_.IPSI[NS-1][*M-1] = KO;
            	CPSIN_.ISIP[KO-1] = NS;
		//		printf("\n\nCPSIN2 %d\n\n", CPSIN_.ISIP[KO-1]);
            	for (int IEL = 1; IEL <= COMPOS_.NELEM[*M-1]; IEL++){
    		   		if (IZZ == COMPOS_.IZ[IEL-1][*M-1])
            		   	STFI[NS-1] = COMPOS_.STF[IEL-1][*M-1];
				}   		
			} else{
				NI=NI+1;
            	CPINAC_.IPIN[NI-1][*M-1] = KO;
            	CPSIN_.ISIP[KO-1] = -NI;
			//	printf("\n\nCPSIN3: %d\n\n", CPSIN_.ISIP[KO-1]);
            	INOUT[NI-1] = 1;
            	for (int IEL = 1; IEL <= COMPOS_.NELEM[*M-1]; IEL++){
            		if (IZZ == COMPOS_.IZ[IEL-1][*M-1])
            		   	STFO[NI-1] = COMPOS_.STF[IEL-1][*M-1];
				}   	
					
			}
				
		} else{
			NI=NI+1;
			CPINAC_.IPIN[NI-1][*M-1] = KO;
   			CPSIN_.ISIP[KO-1] = -NI;
		//	   printf("\n\nCPSIN4 %d\n\n", CPSIN_.ISIP[KO-1]);
   			INOUT[NI-1] = 0;
            for (int IEL = 1; IEL <= COMPOS_.NELEM[*M-1]; IEL++){
            	if (IZZ == COMPOS_.IZ[IEL-1][*M-1])
            	   	STFO[NI-1] = COMPOS_.STF[IEL-1][*M-1];
			}   		
		}	
	}
	
	CPINAC_.NPIN[*M-1] = NI;
	CPSIAC_.NPSI[*M-1] = NS;
	
	//Arrays de Simula��o.
	
	

	
	for (int I = 1; I <= NEGP; I++){
		EE = CEGRID_.ET[I-1];
		pinat2_(EE, WCCM, XH0, XH1, XH2, XS0, XS1, XS2, XT1, XT2, DELTA, M);
		STPC = (XS1+XH1)* COMPOS_.VMOL[*M-1];
		if (ICOR == 1){
			EC = CEGRID_.DLEMP[I-1];
			findi2_(EITL, &EC, &NDATA, &J);
			STPI = exp(CEEL00_.A[J-1] + EC*( CEEL00_.B[J-1] + EC * (CEEL00_.C[J-1]+ EC * CEEL00_.D[J-1])));
			FACT = STPI / STPC;
		} else{
			FACT = 1.0e0;
		}
		
	
        
		CPIMFP_.CSTPP[I-1][*M-1] = STPC * FACT;
        XS1SI = 0.0e0;
        XS2SI = 0.0e0;
        XT1SI = 0.0e0;
        XT2SI = 0.0e0;
        STPSI = 0.0e0;
        STRSI = 0.0e0;
        XS1IN = 0.0e0;
        XS2IN = 0.0e0;
        XT1IN = 0.0e0;
        XT2IN = 0.0e0;
        STPIN = 0.0e0;
        STRIN = 0.0e0;
        
        for( int IO = 1; IO <= CPINAC_.NPIN[*M-1]; IO++){
        	CPSIN_.XSPIN[IO-1][I-1] = 0.0e0;
		}
		
		for( int IO = 1; IO <= CPSIAC_.NPSI[*M-1]; IO++){
        	CPSIN_.XSPSI[IO-1][I-1] = 0.0e0;
		}
		
		for (int KO = 1; KO <= CEIN_.NOSC[*M-1]; KO++){
			int IO = CPSIN_.ISIP[KO-1];
		//	printf("\n\nIO: %d\n\n", CPSIN_.ISIP[KO-1]);
			//Cascas internas
			if (IO > 0){
				IZZ = CEIN_.KZ[KO-1][*M-1];
        		ISH = CEIN_.KS[KO-1][*M-1];
        		EBIN = CADATA_.EB[ISH-1][IZZ-1];
        		if (EBIN < EE){
        			INDC = CPSI0_.IPSIF[IZZ-1] - 1;
  	  	        	PCSI = exp(CPSI0_.XPSI[ISH-1][INDC+I-1]) * STFI[IO-1];
  	  	        	if (PCSI > 1.0e-35){
  	  	        		double AUX = 0.0;
  	  	        		pinat12_(EE, CEIN_.UI[KO-1][*M-1], CEIN_.WRI[KO-1][*M-1], AUX, WCCM, H0, H1, H2, S0, S1, S2, R0, R1, R2);
						STOT = CEIN_.F[KO-1][*M-1] * (S0 + H0);	
						
						if (STOT > 1.0e-35){
							DFERMI = (CPIN00_.SPS0[KO-1] + CPIN00_.SPH0[KO-1]) / STOT;
                  	    	PCSI = PCSI * DFERMI;
                    		CPSIN_.XSPSI[IO-1][I-1] = PCSI;
                    		RNFO = PCSI / (CPIN00_.SPS0[KO-1] + CPIN00_.SPH0[KO-1]);
                    		STPSI = STPSI +(CPIN00_.SPS1[KO-1] + CPIN00_.SPH1[KO-1]) * RNFO;
                    		STRSI = STRSI +(CPIN00_.SPS2[KO-1] + CPIN00_.SPH2[KO-1]) * RNFO;	
						} else{
							CPSIN_.XSPSI[IO-1][I-1] = PCSI;
                    		STPSI = STPSI + PCSI * EBIN;
                    		STRSI = STRSI + PCSI * pow(EBIN,2);	
						}	
					}	
				}	
			}else{
				//outras cascas
				IO = -IO;
				if (INOUT[IO-1] == 1){
					IZZ = CEIN_.KZ[KO-1][*M-1];
        			ISH = CEIN_.KS[KO-1][*M-1];
        			EBIN = CADATA_.EB[ISH-1][IZZ-1];
        			INDC = CPSI0_.IPSIF[IZZ-1] - 1;
  	  	        	PCSI = exp(CPSI0_.XPSI[ISH-1][INDC+I-1]) * STFO[IO-1];
  	  	        	double AUX = 0.0e0;
  	       			pinat12_(EE, CEIN_.UI[KO-1][*M-1], CEIN_.WRI[KO-1][*M-1], AUX, WCCM, H0, H1, H2, S0, S1, S2, R0, R1, R2);
                    STOT = CEIN_.F[KO-1][*M-1] * (S0 + H0);
                    if (STOT > 1.0e-35){
						DFERMI = (CPIN00_.SPS0[KO-1] + CPIN00_.SPH0[KO-1]) / STOT;
                  	   	PCSI = PCSI * DFERMI;
                    	RNFO = PCSI / (CPIN00_.SPS0[KO-1] + CPIN00_.SPH0[KO-1]);
						XS1SI = XS1SI + CPIN00_.SPS1[KO-1] * RNFO;
                		XS2SI = XS2SI + CPIN00_.SPS2[KO-1] * RNFO;
                		XT1SI = XT1SI + CPIN00_.SPT1[KO-1] * RNFO;
                		XT2SI = XT2SI + CPIN00_.SPT2[KO-1] * RNFO;	
						CPSIN_.XSPIN[IO-1][I-1] = CPIN00_.SPH0[KO-1] * RNFO;
                    	STPSI = STPSI +(CPIN00_.SPS1[KO-1] + CPIN00_.SPH1[KO-1]) * RNFO;
                    	STRSI = STRSI +(CPIN00_.SPS2[KO-1] + CPIN00_.SPH2[KO-1]) * RNFO;	
					} else{
						XS1SI = XS1SI + CPIN00_.SPS1[KO-1];
                		XS2SI = XS2SI + CPIN00_.SPS2[KO-1];
                		XT1SI = XT1SI + CPIN00_.SPT1[KO-1];
                		XT2SI = XT2SI + CPIN00_.SPT2[KO-1];	
						CPSIN_.XSPIN[IO-1][I-1] = CPIN00_.SPH0[KO-1];
                    	STPSI = STPSI +(CPIN00_.SPS1[KO-1] + CPIN00_.SPH1[KO-1]);
                    	STRSI = STRSI +(CPIN00_.SPS2[KO-1] + CPIN00_.SPH2[KO-1]);	
					}
					
				} else{
					XS1IN = XS1IN + CPIN00_.SPS1[KO-1];
                	XS2IN = XS2IN + CPIN00_.SPS2[KO-1];
                	XT1IN = XT1IN + CPIN00_.SPT1[KO-1];
                	XT2IN = XT2IN + CPIN00_.SPT2[KO-1];
                	CPSIN_.XSPIN[IO-1][I-1] = CPIN00_.SPH0[KO-1];
                	STPIN = STPIN + (CPIN00_.SPS1[KO-1] + CPIN00_.SPH1[KO-1]);
                	STRIN = STRIN + (CPIN00_.SPS2[KO-1] + CPIN00_.SPH2[KO-1]);	
				}
				
			}	
		}
		
		SSIT = 0.0;
		if (CPSIAC_.NPSI[*M-1] > 0){
			for (int IO =1; IO <= CPSIAC_.NPSI[*M-1]; IO++){
				CPSIAC_.PSIAC[IO-1][I-1][*M-1] = SSIT;
				SSIT = SSIT + CPSIN_.XSPSI[IO-1][I-1];	
			}
		    CPIMFP_.SPISI[I-1][*M-1] = log(fmax(SSIT * COMPOS_.VMOL[*M-1], 1.0e-35));	
		} else{
			CPIMFP_.SPISI[I-1][*M-1] = log(1.0e-35);	
		}
		
		if (SSIT > 1.0e-35){
			for (int IO = 1; IO <= CPSIAC_.NPSI[*M-1]; IO++){
				CPSIAC_.PSIAC[IO-1][I-1][*M-1] = CPSIAC_.PSIAC[IO-1][I-1][*M-1] / SSIT;	
			}
		} else{
			for (int IO = 1; IO <= CPSIAC_.NPSI[*M-1]; IO++){
				CPSIAC_.PSIAC[IO-1][I-1][*M-1] = 1.0e0;	
			}
		}
		
		STPSI = STPSI * COMPOS_.VMOL[*M-1];
        STPIN = STPIN * COMPOS_.VMOL[*M-1];
        FNORM = (CPIMFP_.CSTPP[I-1][*M-1]  - STPSI) / STPIN;
        
        SINT = 0.0e0;
        if (CPINAC_.NPIN[*M-1] == 0){
			printf("PHMATR. No outer shells?\n");
			exit(0);
		}
		
		for (int IO = 1; IO <= CPINAC_.NPIN[*M-1]; IO++){
			if (INOUT[IO-1] == 0)
				CPSIN_.XSPIN[IO-1][I-1] = FNORM * CPSIN_.XSPIN[IO-1][I-1];
			CPINAC_.PINAC[IO-1][I-1][*M-1] = SINT;
			SINT  = SINT + CPSIN_.XSPIN[IO-1][I-1];	
		}
		
		if (SINT > 1.0e-35){
			for (int IO = 1; IO <= CPINAC_.NPIN[*M-1]; IO++){
				CPINAC_.PINAC[IO-1][I-1][*M-1] = CPINAC_.PINAC[IO-1][I-1][*M-1] / SINT;
			}
		} else{
			for (int IO = 1; IO <= CPINAC_.NPIN[*M-1]; IO++){
				CPINAC_.PINAC[IO-1][I-1][*M-1] = 1.0e0;
			}
		}
		CPIMFP_.SPHIN[I-1][*M-1] = log(fmax(SINT * COMPOS_.VMOL[*M-1], 1.0e-35));
        CLAS1P_.TSTRP[I-1][*M-1] = (STRSI + FNORM * STRIN) * COMPOS_.VMOL[*M-1];
        CPIMFP_.W1P[I-1][*M-1] = (XS1SI + FNORM * XS1IN) * COMPOS_.VMOL[*M-1];
        CPIMFP_.W2P[I-1][*M-1] = (XS2SI + FNORM * XS2IN) * COMPOS_.VMOL[*M-1];
        CEINTF_.T1PI[I-1] = (XT1SI + FNORM * XT1IN) * COMPOS_.VMOL[*M-1];
        CEINTF_.T2PI[I-1] = (XT2SI + FNORM * XT2IN) * COMPOS_.VMOL[*M-1];	
	}
	
	//Se��o x de Bremsstrahlung para el�trons.
	

	DIFMAX = 0.0e0;

	for (int I = 1; I <= NDATA; I++){
		if ((EIT[I-1] >= *CEGRID_.EL) && (EIT[I-1] < *CEGRID_.EU)){
	        STPI = F2[I-1] * COMPOS_.RHO[*M-1] * 1.0e6; //! Collision stopping power.
	       	
		    ebrat2_(EIT[I-1], WCRM, XH0, XH1, XH2, XS1, XS2, M);
   
		    STPR = (XS1+XH1)* COMPOS_.VMOL[*M-1];
		    DIFMAX = fmax(DIFMAX, fabs(STPR - STPI) / (0.01e0 * STPI));	
		}		
	}


	
	if ((DIFMAX > 1.0e-3) && (ISCALE == 1)){
		ICOR = 1;
		for (int I = 1; I <= NDATA; I++){
			FL[I-1] = log(F2[I-1] * COMPOS_.RHO[*M-1] * 1.0e6);	
		}
		spline2_(EITL, FL, CEEL00_.A, CEEL00_.B, CEEL00_.C, CEEL00_.D, 0.0e0, 0.0e0, NDATA);
	} else{
		ICOR = 0;
	}
	
	
	for (int I = 1; I <= NEGP; I++){
		ebrat2_(CEGRID_.ET[I-1], WCRM, XH0, XH1, XH2, XS1, XS2, M);
//		printf("\nI: %d, ET: %.15f WCRM: %.15f XH0: %.15e XH1: %.15e XH2 %.15e XS1 %.15e XS2 %.15e, M: %d\n", I, CEGRID_.ET[I-1], WCRM, XH0, XH1, XH2, XS1, XS2, *M);
		STPR = (XS1+XH1)* COMPOS_.VMOL[*M-1];
		if (ICOR == 1){
			EC = CEGRID_.DLEMP[I-1];
			findi2_(EITL, &EC, &NDATA, &J);
			STPI = exp(CEEL00_.A[J-1] + EC*( CEEL00_.B[J-1] + EC * (CEEL00_.C[J-1]+ EC * CEEL00_.D[J-1])));
			FACT = STPI / STPR;
		} else{
			FACT = 1.0e0;
		}
		
		CEIMFP_.RSTPE[I-1][*M-1]=STPR*FACT;
        CLAS1E_.TSTPE[I-1][*M-1]=log(CEIMFP_.CSTPE[I-1][*M-1]+CEIMFP_.RSTPE[I-1][*M-1]);
        CLAS1E_.TSTRE[I-1][*M-1]=log(CLAS1E_.TSTRE[I-1][*M-1]+(XS2+XH2)*COMPOS_.VMOL[*M-1]*FACT);
     //   printf("\n\nXH0: %.15e, VOL: %.15e, FACT: %.15f\n\n", XH0, COMPOS_.VMOL[*M-1], FACT );
        CEIMFP_.SEHBR[I-1][*M-1]=log(fmax(XH0*COMPOS_.VMOL[*M-1]*FACT,1.0e-35));
        
        if (IBREMS == 1){
        	CEIMFP_.W1E[I-1][*M-1]=CEIMFP_.W1E[I-1][*M-1]+XS1*COMPOS_.VMOL[*M-1]*FACT;
            CEIMFP_.W2E[I-1][*M-1]=CEIMFP_.W2E[I-1][*M-1]+XS2*COMPOS_.VMOL[*M-1]*FACT;	
		}
	}
	

	
	
	//Alcance dos el�trons em fun��o da energia.
	
	for (int I = 1; I <= NEGP; I++){
		F1[I-1] = 1.0e0/(CEIMFP_.CSTPE[I-1][*M-1]+CEIMFP_.RSTPE[I-1][*M-1]);	
	}
	
	spline2_(CEGRID_.ET,F1, CEEL00_.A, CEEL00_.B, CEEL00_.C, CEEL00_.D ,0.0e0,0.0e0,NEGP);
	CRANGE_.RANGE[1-1][*M-1][1-1] = 1.0e-8;
	CRANGE_.RANGEL[1-1][*M-1][1-1] = log(CRANGE_.RANGE[1-1][*M-1][1-1]);
	for (int I = 2; I <= NEGP; I++){
		XL= CEGRID_.ET[I-1-1];
        XU= CEGRID_.ET[I-1];
        int w = NEGP;
        sinteg2_(CEGRID_.ET, CEEL00_.A, CEEL00_.B, CEEL00_.C, CEEL00_.D, XL, XU, DR, w);
        CRANGE_.RANGE[I-1][*M-1][1-1]=CRANGE_.RANGE[I-1-1][*M-1][1-1]+DR;
        CRANGE_.RANGEL[I-1][*M-1][1-1]= log(CRANGE_.RANGE[I-1][*M-1][1-1]);
	}		
	
	//Rendimento radiativo de el�trons em fun��o da energia.
	
	for (int I = 1; I <= NEGP; I++){
		F1[I-1] = CEIMFP_.RSTPE[I-1][*M-1] / (CEIMFP_.CSTPE[I-1][*M-1] + CEIMFP_.RSTPE[I-1][*M-1]);
	}
	
	RADY[1-1] = 1.0e-35;
	CBRYLD_.EBRY[1-1][*M-1] = log(RADY[1-1] / CEGRID_.ET[1-1]);
	
	for (int I = 2; I <= NEGP; I++){
		XL= CEGRID_.ET[I-1-1];
        XU= CEGRID_.ET[I-1];
        RADY[I-1] = RADY[I-1-1] + rmomx2_(CEGRID_.ET,F1,XL,XU,NEGP,0);
        CBRYLD_.EBRY[I-1][*M-1]= log(RADY[I-1] / CEGRID_.ET[I-1]);
	}
	
		
	//Bremss de el�trons. rendimento do n�mero de f�tons em fun��o da energia.
	
	for (int I = 1; I <= NEGP; I++){
		F1[I-1] = exp(CEIMFP_.SEHBR[I-1][*M-1]) / (CEIMFP_.CSTPE[I-1][*M-1] + CEIMFP_.RSTPE[I-1][*M-1]);	
	}
	RADN[1-1] = 0.0e0;
	for (int I = 2; I <= NEGP; I++){
		XL= CEGRID_.ET[I-1-1];
        XU= CEGRID_.ET[I-1];
        RADN[I-1] = RADN[I-1-1] + rmomx2_(CEGRID_.ET,F1,XL,XU,NEGP,0);
	}
	
	//Imprimir tabelas de pot�ncia de interrup��o de el�trons.
	
	if (*INFO >= 3){
		fprintf(IWR,"\n PENELOPE >>>  Stopping powers for electrons\n");
		fprintf(IWR, "\n   Energy         Scol         Srad          range      Rad. Yield     PhotonYield      delta\n    (eV)        (eV/mtu)      (eV/mtu)        (mtu)                     (W>WCR)\n");
        fprintf(IWR, "---------------------------------------------------------------------------------------\n");
		
		for (int I = 1; I <= NEGP; I++){
			fprintf(IWR, " %.5E  %.5E  %.5E  %.5E  %.5E  %.5E  %.5E\n", CEGRID_.ET[I-1], CEIMFP_.CSTPE[I-1][*M-1]/COMPOS_.RHO[*M-1],
          	  	    CEIMFP_.RSTPE[I-1][*M-1]/COMPOS_.RHO[*M-1], CRANGE_.RANGE[I-1][*M-1][1-1]*COMPOS_.RHO[*M-1], RADY[I-1]/CEGRID_.ET[I-1],
          	  	  	RADN[I-1], CEIMFP_.DEL[I-1][*M-1]);
		}
	}
	
	
	
	//Se��o x de Bremss para p�sitrons.
	
	DIFMAX = 0.0;

	for (int I = 1; I <= NDATA; I++){
		if ((EIT[I-1] >= *CEGRID_.EL) && (EIT[I-1] < *CEGRID_.EU)){
	        STPI = F4[I-1] * COMPOS_.RHO[*M-1] * 1.0e6; //! Collision stopping power.
	       	
		    pbrat2_(EIT[I-1], WCRM, XH0, XH1, XH2, XS1, XS2, M);
   
		    STPR = (XS1+XH1)* COMPOS_.VMOL[*M-1];
		    DIFMAX = fmax(DIFMAX, fabs(STPR - STPI) / (0.01e0 * STPI));	
		}		
	}
	
	if ((DIFMAX > 1.0e-3) && (ISCALE == 1)){
		ICOR = 1;
		for (int I = 1; I <= NDATA; I++){
			FL[I-1] = log(F4[I-1] * COMPOS_.RHO[*M-1] * 1.0e6);	
		}
		spline2_(EITL, FL, CEEL00_.A, CEEL00_.B, CEEL00_.C, CEEL00_.D, 0.0, 0.0, NDATA);
	} else{
		ICOR = 0;
	}
	
	
	for (int I = 1; I <= NEGP; I++){
		pbrat2_(CEGRID_.ET[I-1], WCRM, XH0, XH1, XH2, XS1, XS2, M);
		STPR = (XS1+XH1)* COMPOS_.VMOL[*M-1];
		if (ICOR == 1){
			EC = CEGRID_.DLEMP[I-1];
			findi2_(EITL, &EC, &NDATA, &J);
			STPI = exp(CEEL00_.A[J-1] + EC*(CEEL00_.B[J-1] + EC * (CEEL00_.C[J-1]+ EC * CEEL00_.D[J-1])));
			FACT = STPI / STPR;
		} else{
			FACT = 1.0e0;
		}
		
		CPIMFP_.RSTPP[I-1][*M-1]=STPR*FACT;
        CLAS1P_.TSTPP[I-1][*M-1]=log(CPIMFP_.CSTPP[I-1][*M-1]+CPIMFP_.RSTPP[I-1][*M-1]);
        CLAS1P_.TSTRP[I-1][*M-1]=log(CLAS1P_.TSTRP[I-1][*M-1]+(XS2+XH2)*COMPOS_.VMOL[*M-1]*FACT);
        CPIMFP_.SPHBR[I-1][*M-1]=log(fmax(XH0*COMPOS_.VMOL[*M-1]*FACT,1.0e-35));
        
        if (IBREMS == 1){
        	CPIMFP_.W1P[I-1][*M-1]=CPIMFP_.W1P[I-1][*M-1]+XS1*COMPOS_.VMOL[*M-1]*FACT;
            CPIMFP_.W2P[I-1][*M-1]=CPIMFP_.W2P[I-1][*M-1]+XS2*COMPOS_.VMOL[*M-1]*FACT;	
		}
	}
	

	//Aniquila��o de Positron
	
	for (int I = 1; I <= NEGP; I++){
		panat2_(CEGRID_.ET[I-1], TXAN);
		CPIMFP_.SPAN[I-1][*M-1] = log(TXAN * COMPOS_.ZT[*M-1]  * COMPOS_.VMOL[*M-1] );
	}
	
	//Faixa de p�sitrons em fun��o da energia.
	
	for (int I = 1; I <= NEGP; I++){
		F3[I-1] = 1.0e0/(CPIMFP_.CSTPP[I-1][*M-1]+CPIMFP_.RSTPP[I-1][*M-1]);
			
	}
	
	spline2_(CEGRID_.ET,F3, CEEL00_.A, CEEL00_.B, CEEL00_.C, CEEL00_.D ,0.0e0,0.0e0,NEGP);
	CRANGE_.RANGE[1-1][*M-1][3-1] = 1.0e-8;
	CRANGE_.RANGEL[1-1][*M-1][3-1] = log(CRANGE_.RANGE[1-1][*M-1][3-1]);
	for (int I = 2; I <= NEGP; I++){
		XL= CEGRID_.ET[I-1-1];
        XU= CEGRID_.ET[I-1];
        int w = NEGP;
        sinteg2_(CEGRID_.ET, CEEL00_.A, CEEL00_.B, CEEL00_.C, CEEL00_.D, XL, XU, DR, w);
        CRANGE_.RANGE[I-1][*M-1][3-1]=CRANGE_.RANGE[I-1-1][*M-1][3-1]+DR;
        CRANGE_.RANGEL[I-1][*M-1][3-1]= log(CRANGE_.RANGE[I-1][*M-1][3-1]);
	}
		
	
	//Rendimento radiativo de positron em fun��o da energia.
	
	for (int I = 1; I <= NEGP; I++){
		F3[I-1] = CPIMFP_.RSTPP[I-1][*M-1] / (CPIMFP_.CSTPP[I-1][*M-1] + CPIMFP_.RSTPP[I-1][*M-1]);
	}
	
	RADY[1-1] = 1.0e-35;
	CBRYLD_.PBRY[1-1][*M-1] = log(RADY[1-1] / CEGRID_.ET[1-1]);
	
	
	for (int I = 2; I <= NEGP; I++){
		XL= CEGRID_.ET[I-1-1];
        XU= CEGRID_.ET[I-1];
        RADY[I-1] = RADY[I-1-1] + rmomx2_(CEGRID_.ET,F3,XL,XU,NEGP,0);
        CBRYLD_.PBRY[I-1][*M-1]= log(RADY[I-1] / CEGRID_.ET[I-1]);
	}
		
		
	//Bremss de positrons. rendimento do n�mero de f�tons em fun��o da energia.
	
	for (int I = 1; I <= NEGP; I++){
		F3[I-1] = exp(CPIMFP_.SPHBR[I-1][*M-1]) / (CPIMFP_.CSTPP[I-1][*M-1] + CPIMFP_.RSTPP[I-1][*M-1]);	
	}
	RADN[1-1] = 1.0e-35;
	for (int I = 2; I <= NEGP; I++){
		XL= CEGRID_.ET[I-1-1];
        XU= CEGRID_.ET[I-1];
        RADN[I-1] = RADN[I-1-1] + rmomx2_(CEGRID_.ET,F3,XL,XU,NEGP,0);
	}
	
	//Imprimir tabelas de pot�ncia de interrup��o de positrons.
	
	if (*INFO >= 3){
		fprintf(IWR,"\n PENELOPE >>>  Stopping powers for positrons\n");
		fprintf(IWR, "\n   Energy        Scol         Srad         range     Rad. Yield   PhotonYield  annih. mfp\n    (eV)       (eV/mtu)     (eV/mtu)       (mtu)                    (W>WCR)      (mtu)\n");
		fprintf(IWR, " ---------------------------------------------------------------------------------------\n");
		for (int I = 1; I <= NEGP; I++){
			fprintf(IWR, " %.5E  %.5E  %.5E  %.5E  %.5E  %.5E  %.5E\n", CEGRID_.ET[I-1],CPIMFP_.CSTPP[I-1][*M-1]/COMPOS_.RHO[*M-1],
          	  	    CPIMFP_.RSTPP[I-1][*M-1]/COMPOS_.RHO[*M-1], CRANGE_.RANGE[I-1][*M-1][3-1]*COMPOS_.RHO[*M-1],RADY[I-1]/CEGRID_.ET[I-1],
          	  	  	RADN[I-1],COMPOS_.RHO[*M-1]/exp(CPIMFP_.SPAN[I-1][*M-1]));
		}
	}

	
	if (*INFO >= 3){
		fprintf(IWR, "\nPENELOPE >>>  Soft stopping power and energy straggling\n");
        fprintf(IWR, "  \n    Energy       SStp,e-      SStr,e-      STP,e-       SStp,e+      SStr,e+     Stp,e+\n     (eV)        (eV/mtu)     (eV**2/mtu)   soft/tot    (eV/mtu)     (eV**2/mtu)   soft/tot\n");
		fprintf(IWR, " ---------------------------------------------------------------------------------------\n");

	}	
	
	
	FMTU=1.0e0/COMPOS_.RHO[*M-1];
	for (int I = 1; I <= NEGP; I++){
		//Os eventos de perda de energia suave s�o desligados quando W1 � muito pequeno.
		FSOFTE=CEIMFP_.W1E[I-1][*M-1]/(CEIMFP_.CSTPE[I-1][*M-1]+CEIMFP_.RSTPE[I-1][*M-1]);
        FSOFTP=CPIMFP_.W1P[I-1][*M-1]/(CPIMFP_.CSTPP[I-1][*M-1]+CPIMFP_.RSTPP[I-1][*M-1]);
		if (CEIMFP_.W1E[I-1][*M-1] > 1.0e-5*(CEIMFP_.CSTPE[I-1][*M-1]+CEIMFP_.RSTPE[I-1][*M-1])){
			W1EW=CEIMFP_.W1E[I-1][*M-1];
            W2EW=CEIMFP_.W2E[I-1][*M-1];
            CEIMFP_.W1E[I-1][*M-1]=log(fmax(CEIMFP_.W1E[I-1][*M-1],1.0e-35));
            CEIMFP_.W2E[I-1][*M-1]=log(fmax(CEIMFP_.W2E[I-1][*M-1],1.0e-35));	
		}else{
			W1EW=0.0e0;
            W2EW=0.0e0;
            CEIMFP_.W1E[I-1][*M-1]=-100.0e0;
            CEIMFP_.W2E[I-1][*M-1]=-100.0e0;
		}
		if (CPIMFP_.W1P[I-1][*M-1] > 1.0e-5*(CPIMFP_.CSTPP[I-1][*M-1]+CPIMFP_.RSTPP[I-1][*M-1])){
			W1PW=CPIMFP_.W1P[I-1][*M-1];
            W2PW=CPIMFP_.W2P[I-1][*M-1];
            CPIMFP_.W1P[I-1][*M-1]=log(fmax(CPIMFP_.W1P[I-1][*M-1],1.0e-35));
            CPIMFP_.W2P[I-1][*M-1]=log(fmax(CPIMFP_.W2P[I-1][*M-1],1.0e-35));	
		}else{
			W1PW=0.0e0;
            W2PW=0.0e0;
            CPIMFP_.W1P[I-1][*M-1]=-100.0e0;
            CPIMFP_.W2P[I-1][*M-1]=-100.0e0;
		}
		if (*INFO >= 3){
			fprintf(IWR, " %.5E  %.5E  %.5E  %.2E  %.5E  %.5E  %.2E\n", CEGRID_.ET[I-1], W1EW*FMTU, W2EW*FMTU, FSOFTE, W1PW*FMTU, W2PW*FMTU, FSOFTP);
		}	
	}

	//Espalhamento el�stico de el�trons e p�sitrons.
	
	eelar2_(M, IRD, IWR, INFO);

	for (int I = 1; I <= NEGP; I++){
		CLAS1E_.TRL1E[I-1][*M-1]=-log(CEEL00_.XE1[I-1]*COMPOS_.VMOL[*M-1]);
        CLAS1E_.TRL2E[I-1][*M-1]=-log(CEEL00_.XE2[I-1]*COMPOS_.VMOL[*M-1]);
        CLAS1P_.TRL1P[I-1][*M-1]=-log(CEEL00_.XP1[I-1]*COMPOS_.VMOL[*M-1]);
        CLAS1P_.TRL2P[I-1][*M-1]=-log(CEEL00_.XP2[I-1]*COMPOS_.VMOL[*M-1]);	
	}
	
	eeldr2_(M, IRD, IWR, INFO); // ! Usa o banco de dados ELSEPA.
	
	if (*INFO >= 3){
		fprintf(IWR, "\n PENELOPE >>>  Soft angular transport coefficients\n");
		fprintf(IWR, "\n   Energy       SITMFP1,e-    SITMFP2,e-    SITMFP1,e+    SITMFP2,e+\n    (eV)         (1/mtu)        (1/mtu)        (1/mtu)        (1/mtu)\n");
        fprintf(IWR, "-----------------------------------------------------------------------\n");
		for (int I = 1; I <= NEGP; I++){
			fprintf(IWR," %.5E  %.5E  %.5E  %.5E  %.5E\n", 
	        CEGRID_.ET[I-1],exp(CEIMFP_.T1E[I-1][*M-1])*FMTU,exp(CEIMFP_.T2E[I-1][*M-1])*FMTU, exp(CPIMFP_.T1P[I-1][*M-1])*FMTU,exp(CPIMFP_.T2P[I-1][*M-1])*FMTU);
		}
	}
	
	
	
	//Imprimir tabelas de caminho livre de m�dia de el�trons.
	
	if (*INFO >= 3){
		fprintf(IWR, "\n PENELOPE >>>  Electron mean free paths (hard events)\n");
		fprintf(IWR, " *** NOTE: The MFP for inner-shell ionisation(isi) is listed only for\n");
		fprintf(IWR, "           completeness. The MFP or inelastic collisions (in) has been\n");
		fprintf(IWR, "           calculated by considering all inelastic events, including isi.\n");
		fprintf(IWR, "\n   Energy        MFPel        MFPin        MFPbr       MFPtot       MFPisi\n    (eV)         (mtu)        (mtu)        (mtu)        (mtu)        (mtu)\n");
	    fprintf(IWR, " -----------------------------------------------------------------------\n");
	}
	
	for (int I = 1; I <= NEGP; I++){
		FPEL=exp(CEIMFP_.SEHEL[I-1][*M-1]);
        FPSI=exp(CEIMFP_.SEISI[I-1][*M-1]);
        FPIN=exp(CEIMFP_.SEHIN[I-1][*M-1])+FPSI;
        FPBR=exp(CEIMFP_.SEHBR[I-1][*M-1]);
        FPTOT=FPEL+FPIN+FPBR;
        if (*INFO >= 3){
			fprintf(IWR, " %.5E  %.5E  %.5E  %.5E  %.5E  %.5E\n", 
			CEGRID_.ET[I-1],COMPOS_.RHO[*M-1]/FPEL,COMPOS_.RHO[*M-1]/FPIN,COMPOS_.RHO[*M-1]/FPBR,COMPOS_.RHO[*M-1]/FPTOT,COMPOS_.RHO[*M-1]/FPSI);
		}
    	CEIMFP_.SEAUX[I-1][*M-1]=log(1.0e-35);
        CEIMFP_.SETOT[I-1][*M-1]=log(FPTOT);
	//	printf("\n\nSETOT %f\n\n", CEIMFP_.SETOT[I-1][*M-1])	;
	}
	
	for (int I = 2; I <= NEGP-1; I++){
		if ((exp(CEIMFP_.SETOT[I-1][*M-1]) > 1.005e0*exp(CEIMFP_.SETOT[I-1-1][*M-1])) &&
            (exp(CEIMFP_.SETOT[I-1][*M-1]) > 1.005e0*exp(CEIMFP_.SETOT[I+1-1][*M-1])) &&
            (CEGRID_.ET[I-1] > PENELOPE_mod_.EABS[*M-1][1-1]) && (CEGRID_.ET[I-1] < 1.0e6)){
            	fprintf(IWR, "WARNING: The electron hard IMFP has a maximum at E = %.5E eV --- I: %d\n", CEGRID_.ET[I-1], I);	
			}
	}
	
	//Imprimir tabelas de caminho livre de m�dia de positrons.
	
	if (*INFO >= 3){
		fprintf(IWR, "\n PENELOPE >>>  Positron mean free paths (hard events)\n");
		fprintf(IWR, " *** NOTE: The MFP for inner-shell ionisation(isi) is listed only for\n");
		fprintf(IWR, "           completeness. The MFP or inelastic collisions (in) has been\n");
		fprintf(IWR, "           calculated by considering all inelastic events, including isi.\n");
		fprintf(IWR, "\n   Energy         MFPel         MFPin         MFPbr         MFPan        MFPtot        MFPisi\n    (eV)          (mtu)         (mtu)         (mtu)         (mtu)          (mtu)        (mtu)\n");
    	fprintf(IWR, " ---------------------------------------------------------------------------------------------\n");
	}
	
	for (int I = 1; I <= NEGP; I++){
		FPEL=exp(CPIMFP_.SPHEL[I-1][*M-1]);
        FPSI=exp(CPIMFP_.SPISI[I-1][*M-1]);
        FPIN=exp(CPIMFP_.SPHIN[I-1][*M-1])+FPSI;
        FPBR=exp(CPIMFP_.SPHBR[I-1][*M-1]);
        FPAN=exp(CPIMFP_.SPAN[I-1][*M-1]);
        FPTOT=FPEL+FPIN+FPBR+FPAN;
        if (*INFO >= 3){
			fprintf(IWR, " %.5E  %.5E  %.5E  %.5E  %.5E  %.5E  %.5E \n", 
			CEGRID_.ET[I-1],COMPOS_.RHO[*M-1]/FPEL,COMPOS_.RHO[*M-1]/FPIN,COMPOS_.RHO[*M-1]/FPBR,COMPOS_.RHO[*M-1]/FPAN,COMPOS_.RHO[*M-1]/FPTOT,COMPOS_.RHO[*M-1]/FPSI);
		}
    	CPIMFP_.SPAUX[I-1][*M-1]=log(1.0e-35);
        CPIMFP_.SPTOT[I-1][*M-1]=log(FPTOT);	
	}
	
	for (int I = 2; I <= NEGP-1; I++){
		if ((exp(CPIMFP_.SPTOT[I-1][*M-1]) > 1.005e0*exp(CPIMFP_.SPTOT[I-1-1][*M-1])) &&
            (exp(CPIMFP_.SPTOT[I-1][*M-1]) > 1.005e0*exp(CPIMFP_.SPTOT[I+1-1][*M-1])) &&
            (CEGRID_.ET[I-1] > PENELOPE_mod_.EABS[*M-1][3-1]) && (CEGRID_.ET[I-1] < 1.0e6)){
            	fprintf(IWR, "WARNING: The positron hard IMFP has a maximum  at E = %.5E eV\n ", CEGRID_.ET[I-1]);	
			}
	}

	
	//Corre��o para a depend�ncia energ�tica do poder de parada e o par�metro de dispers�o de energia de intera��es suaves.
	
	
	for (int I = 1; I <= NEGP -1; I++){
		CEIMFP_.DW1EL[I-1][*M-1]=(CEIMFP_.W1E[I+1-1][*M-1]-CEIMFP_.W1E[I-1][*M-1])/(CEGRID_.ET[I+1-1]-CEGRID_.ET[I-1]);
        CEIMFP_.DW2EL[I-1][*M-1]=0.5e0*(CEIMFP_.W2E[I+1-1][*M-1]-CEIMFP_.W2E[I-1][*M-1])/(CEGRID_.ET[I+1-1]-CEGRID_.ET[I-1])+CEIMFP_.DW1EL[I-1][*M-1];
		CEIMFP_.DW1EL[I-1][*M-1]=0.5e0*CEIMFP_.DW1EL[I-1][*M-1];
		
		CPIMFP_.DW1PL[I-1][*M-1]=(CPIMFP_.W1P[I+1-1][*M-1]-CPIMFP_.W1P[I-1][*M-1])/(CEGRID_.ET[I+1-1]-CEGRID_.ET[I-1]);
        CPIMFP_.DW2PL[I-1][*M-1]=0.5e0*(CPIMFP_.W2P[I+1-1][*M-1]-CPIMFP_.W2P[I-1][*M-1])/(CEGRID_.ET[I+1-1]-CEGRID_.ET[I-1])+CPIMFP_.DW1PL[I-1][*M-1];
		CPIMFP_.DW1PL[I-1][*M-1]=0.5e0*CPIMFP_.DW1PL[I-1][*M-1];	
	}
 	CEIMFP_.DW1EL[NEGP-1][*M-1]=0.0e0;
    CEIMFP_.DW2EL[NEGP-1][*M-1]=0.0e0;
    CPIMFP_.DW1PL[NEGP-1][*M-1]=0.0e0;
    CPIMFP_.DW2PL[NEGP-1][*M-1]=0.0e0;
    
    
    //Propriedades de Interacoes dos Fotons
       
    //Espalhamento Rayleigh

    graar2_(M,IRD,IWR,INFO);
  
    
    //espalhamento Compton e producao de pares
    
    fgets(LINHA, sizeof(LINHA), IRD);
	extrairString(APOIO, LINHA, 57, 4);
	NDATA = atoi(APOIO);
//	printf("\n\nNDATA: %d\n\n", NDATA);
//	printf("\n\nLINHA %s\n\n", LINHA);
    
    if (*INFO >= 2){
		fprintf(IWR, "\n*** Compton and pair-production cross sections,  NDATA = %4d\n", NDATA);
	}
	if (NDATA > NEGP){
		printf("PEMATR. Too many data points (2).\n");
		exit(0);
	}
	if (*INFO >=2){
		fprintf(IWR, "\n  Energy      CS-Comp      CS-pair      CS-triplet\n   (eV)       (cm**2)      (cm**2)      (cm**2)\n");
		fprintf(IWR, " ---------------------------------------------------------------\n");
	}
	for (int I = 1; I <= NDATA; I++){
		fgets(LINHA, sizeof(LINHA), IRD);
		extrairString(APOIO, LINHA, 0, 10);
		EIT[I-1] = atof(APOIO);
		extrairString(APOIO, LINHA, 10, 12);
		F2[I-1] = atof(APOIO);
		extrairString(APOIO, LINHA, 22, 12);
		F3[I-1] = atof(APOIO);
		extrairString(APOIO, LINHA, 34, 12);
		F4[I-1] = atof(APOIO);
		if (*INFO >= 2){
			fprintf(IWR, "%.3E %.5E %.5E %.5E\n",EIT[I-1],F2[I-1], F3[I-1], F4[I-1] );
		}
		EITL[I-1]=log(EIT[I-1]);		
	}
	
	//espalhamento compton
	VMOLL=log(COMPOS_.VMOL[*M-1]);
	for (int I = 1; I <= NDATA; I++){
		FL[I-1]=log(fmax(F2[I-1],1.0e-35));
	}
//	printf("\n\nNP: %d \n\n", NDATA);
	spline2_(EITL,FL,CEEL00_.A,CEEL00_.B,CEEL00_.C,CEEL00_.D,0.0e0,0.0e0,NDATA);
	for (int I = 1; I <= NEGP; I++){
		EC=CEGRID_.DLEMP[I-1];
        findi2_(EITL,&EC,&NDATA,&J);
        CGIMFP_.SGCO[I-1][*M-1]=CEEL00_.A[J-1]+EC*(CEEL00_.B[J-1]+EC*(CEEL00_.C[J-1]+EC*CEEL00_.D[J-1]))+VMOLL;	
	}
	
	//Par el�tron-p�sitron e produ��o de trig�meos.
	
	//Pares
	
	NP = 0;
	for (int I = 1; I <= NDATA; I++){
		if (EIT[I-1] < 1.023e6)
			goto L1;
		NP=NP+1;
        F1[NP-1]=EITL[I-1];
        FACT= pow(EIT[I-1]/(EIT[I-1]-1.022e6),3);
        FL[NP-1]=log(FACT*fmax(F3[I-1],1.0e-35));
L1:;
	}
//	printf("\n\nNP: %d \n\n", NP);
	spline2_(F1,FL,CEEL00_.A,CEEL00_.B,CEEL00_.C,CEEL00_.D,0.0e0,0.0e0,NP);
	for (int I = 1; I <= NEGP; I++){
		if (CEGRID_.ET[I-1] < 1.023e6)
			CGIMFP_.SGPP[I-1][*M-1]=-80.6e0;
		else{
			FACT=pow(CEGRID_.ET[I-1]/(CEGRID_.ET[I-1]-1.022e6),3);
            EC=CEGRID_.DLEMP[I-1];
            findi2_(F1,&EC,&NP,&J);
            CGIMFP_.SGPP[I-1][*M-1]=CEEL00_.A[J-1]+EC*(CEEL00_.B[J-1]+EC*(CEEL00_.C[J-1]+EC*CEEL00_.D[J-1]))+VMOLL-log(FACT);	
		}		
	}
	
	//Tripletos
	
	NP = 0;
	for (int I = 1; I <= NDATA; I++){
		if (EIT[I-1] < 2.045e6)
			goto L2;
		NP=NP+1;
        FACT=pow(EIT[I-1]/(EIT[I-1]-2.044e6),3);
        F1[NP-1]=EITL[I-1];
        FL[NP-1]=log(FACT*fmax(F4[I-1],1.0e-35));
L2:;
	}

	spline2_(F1,FL,CEEL00_.A,CEEL00_.B,CEEL00_.C,CEEL00_.D,0.0e0,0.0e0,NP);
	for (int I = 1; I <= NEGP; I++){
		if (CEGRID_.ET[I-1] < 2.045e6)
			CGPP01_.TRIP[I-1][*M-1]=-80.6e0;
		else{
			FACT=pow(CEGRID_.ET[I-1]/(CEGRID_.ET[I-1]-2.044e6),3);
            EC=CEGRID_.DLEMP[I-1];
            findi2_(F1,&EC,&NP,&J);
            TRIPL=exp(CEEL00_.A[J-1]+EC*(CEEL00_.B[J-1]+EC*(CEEL00_.C[J-1]+EC*CEEL00_.D[J-1]))+VMOLL-log(FACT));
            PAIR=exp(CGIMFP_.SGPP[I-1][*M-1]);
            CGIMFP_.SGPP[I-1][*M-1]=log(PAIR+TRIPL);
        	CGPP01_.TRIP[I-1][*M-1]=TRIPL/(PAIR+TRIPL);	
		}		
	}
	
	if (CGCO_.NOSCCO[*M-1] > 1){
		CGCO_.PTRSH[1-1][*M-1]=CGCO_.FCO[1-1][*M-1];
		for (int I = 2; I <= CGCO_.NOSCCO[*M-1]; I++){
			CGCO_.PTRSH[I-1][*M-1]=CGCO_.PTRSH[I-1-1][*M-1]+CGCO_.FCO[I-1][*M-1];
		}
		
		for (int I = 1; I <= CGCO_.NOSCCO[*M-1]; I++){
			CGCO_.PTRSH[I-1][*M-1]=CGCO_.PTRSH[I-1][*M-1]/CGCO_.PTRSH[CGCO_.NOSCCO[*M-1]-1][*M-1];
		}
	} else{
		CGCO_.PTRSH[1-1][*M-1]=1.0e0;	
	}
	
	

	//Absorcao Fotoeletrica
	
	gphar2_(M,IRD,IWR,INFO);

	
	//'Alcance' do f�ton (= caminho livre m�dio, ligeiramente subestimado).
	
	for (int KE = 1; KE <= NEGP; KE++){
		graati2_(CEGRID_.ET[KE-1],ECS,M);
        PRA=ECS*COMPOS_.VMOL[*M-1];
        PCO=exp(CGIMFP_.SGCO[KE-1][*M-1]);
        if (CEGRID_.ET[KE-1] < 1.023e6)
        	PPP=1.0e-35;
        else
        	PPP=exp(CGIMFP_.SGPP[KE-1][*M-1]);
		PPH=CGIMFP_.SGPH[KE-1][*M-1];
        PT=(PRA+PCO+PPP+PPH);
        CRANGE_.RANGE[KE-1][*M-1][2-1]=1.0e0/PT;
        CRANGE_.RANGEL[KE-1][*M-1][2-1]=log(CRANGE_.RANGE[KE-1][*M-1][2-1]);
	}
	

	
	if (*INFO >= 3){
		fprintf(IWR, "\n PENELOPE >>>  Photon mass attenuation coefficients'\n");
		fprintf(IWR, "\n   Energy       Rayleigh       Compton     Photoelect      Pair          Total\n    (eV)         (1/mtu)       (1/mtu)       (1/mtu)       (1/mtu)       (1/mtu)\n");
	    fprintf(IWR, " ---------------------------------------------------------------------------------\n");
	}
	
	for (int I = 1; I <= *CGPH01_.NPHD; I++){
		if ((CGPH01_.ER[I-1] >= *CEGRID_.EL) && (CGPH01_.ER[I-1] <= *CEGRID_.EU)){
			*CEGRID_.XEL=log(CGPH01_.ER[I-1]);
            *CEGRID_.XE=1.0e0+(*CEGRID_.XEL-*CEGRID_.DLEMP1)* *CEGRID_.DLFC;
            *CEGRID_.KE=*CEGRID_.XE;
            *CEGRID_.XEK=*CEGRID_.XE-*CEGRID_.KE;
            graati2_(CGPH01_.ER[I-1],ECS,M);
            PRA=ECS*COMPOS_.VMOL[*M-1]/COMPOS_.RHO[*M-1];
            PCO=exp(CGIMFP_.SGCO[*CEGRID_.KE-1][*M-1]+(CGIMFP_.SGCO[*CEGRID_.KE+1-1][*M-1]-CGIMFP_.SGCO[*CEGRID_.KE-1][*M-1])* *CEGRID_.XEK)/COMPOS_.RHO[*M-1];
			if (CGPH01_.ER[I-1] < 1.022e6)
				PPP=1.0e-35;
			else
				PPP=exp(CGIMFP_.SGPP[*CEGRID_.KE-1][*M-1]+(CGIMFP_.SGPP[*CEGRID_.KE+1-1][*M-1]-CGIMFP_.SGPP[*CEGRID_.KE-1][*M-1])* *CEGRID_.XEK)/COMPOS_.RHO[*M-1];
			PPH=CGPH01_.XSR[I-1]*COMPOS_.VMOL[*M-1]/COMPOS_.RHO[*M-1];
            PT=PRA+PCO+PPP+PPH;	
            if (*INFO >= 3){
				fprintf(IWR, " %.5E  %.5E  %.5E  %.5E  %.5E  %.5E\n", CGPH01_.ER[I-1],PRA,PCO,PPH,PPP,PT);
			}
				
		}
	}
	
	//Produ��o de pares. Rotina de inicializa��o.

	gppa02_(M);
	
	for (int IE = 1; IE <= NEGP; IE++){
		CGIMFP_.SGAUX[IE-1][*M-1]=log(1.0e-35);
	}
	
	strcpy(LNAME, "\n PENELOPE (v. 2014)  End of material data file ........");
	fgets(NAME, sizeof(NAME), IRD);
	if (!strcmp(LNAME, NAME)){
		printf("I/O error. Corrupt material data file.\n");
		printf("      The last line is: %s '\n", NAME);
		printf("     ... and should be: %s \n", LNAME);
		fprintf(IWR, "I/O error. Corrupt material data file.\n");
		fprintf(IWR, "      The last line is: %s '\n", NAME);
		fprintf(IWR, "     ... and should be: %s \n", LNAME);
		exit(0);	
	}
	
	fprintf(IWR, "%s\n\n", NAME);
	
	free(STFI);
	free(STFO);
	free(INOUT);
	free(EIT);
	free(EITL);
	free(FL);
	free(F1);
	free(F2);
	free(F3);
	free(F4);
	free(RADY);
	free(RADN);
	
//	printf("FIM PMATR\n");
//	 exit(0);
}

void irnd02_(double *W, double *F, int *K, int *N){ //OK
	

//	printf("\nirnd02\n");
	/*
	Inicializa��o do algoritmo de aliasing de Walker para amostragem aleat�ria
   de distribui��es de probabilidade discretas.
 
   Argumentos de entrada:
     N ........ n�mero de diferentes valores da vari�vel aleat�ria.
     W (1: N) ... probabilidades de pontos correspondentes (n�o necessariamente
                normalizado para a unidade).
   Argumentos de sa�da:
     F (1: N) ... valores de corte.
     K (1: N) ... valores de alias.
	*/
	
	
	//Renormalizacao
	
	double HLOW;
	double HIGH;
	int ILOW;
	int IHIGH;
	double FACT;
	
	*CRITAN_.CNORM = 0.0e0;
	for (int I = 1; I <= *N; I++){
		if (W[I-1] < 0.0e0){
			printf("IRND0. Negative point probability.");
			exit(0);
		}
		*CRITAN_.CNORM = *CRITAN_.CNORM + W[I-1];
	}
	
	*CRITAN_.CNORM = 1.0e0 / *CRITAN_.CNORM;
	FACT = double(*N) * (*CRITAN_.CNORM);
//	printf("\n\nN: %d, CNORM: %f, valor FACT: %f", *N, *CRITAN_.CNORM, FACT);
	for (int I = 1; I <= *N; I++){
		K[I-1] = I;
		F[I-1] = W[I-1] * FACT;	
	}
	
	if (*N == 1){
		return;
	}
	
	for (int I = 1; I <= (*N - 1); I++){
		HLOW = 1.0e0;
		HIGH = 1.0e0;
		ILOW = 0;
		IHIGH = 0;
		for (int J = 1; J <= *N; J++){
			if (K[J-1] == J){
				if (F[J-1] < HLOW){
					HLOW = F[J-1];
					ILOW = J;
				//	printf("\n\nILOW C++: %d, HLOW: %17.17f\n\n", ILOW, HLOW);
				} else if (F[J-1] > HIGH){
					HIGH = F[J-1];
					IHIGH = J;
				//	printf("\n\nIHIGH C++: %d, HIGH: %17.17f\n\n", IHIGH, HIGH);
				}
			}
		}
		if ((ILOW == 0) || (IHIGH == 0))
			return;
		K[ILOW-1] = IHIGH;
		F[IHIGH-1] = HIGH+HLOW-1.0e0;	
	}	
	
//	printf("\nIRND02\n");
}


void ritai02_(int &PDF, double &XLOW, double &XHIGH, int &N, int &NU, double &ERRM, FILE *IWR){ //OK
	
	/*
	Inicializa��o do algoritmo RITA para amostragem aleat�ria de um
   vari�vel aleat�ria cont�nua X de uma fun��o de distribui��o de probabilidade
   PDF (X) definido no intervalo (XLOW, XHIGH). A fun��o externa
   PDF (X) --n�o necessariamente normalizado-- deve ser fornecido pelo usu�rio.
   N � o n�mero de pontos na grade de amostragem. Esses pontos s�o
   determinado por meio de uma estrat�gia adaptativa que minimiza local
   erros de interpola��o. Os primeiros pontos da grade NU s�o uniformemente espa�ados
   em (XLOW, XHIGH); quando NU � negativo, a grade inicial consiste em
   -NU pontos espa�ados logaritmicamente (neste caso, XLOW deve ser
   n�o negativo).
 
   ERRM � uma medida do erro de interpola��o (o maior valor de
   o erro absoluto da interpola��o racional integrada sobre cada
   intervalo de grade).
 
   **** Coeficientes de interpola��o e tabelas PDF s�o impressos em
         arquivos separados (UNIT = IWR) se IWR for maior que zero.
 
   Outros subprogramas necess�rios: Fun��o EXTERNA PDF.
 
      PRECIS�O DUPLA IMPL�CITA (A-H, O-Z), INTEIRO * 4 (I-N)
      PAR�METRO (EPS = 1.0D-10, ZERO = 1.0D-75, ZEROT = 0.1D0 * ZERO)
      PAR�METRO (NM = 512)
 
   A informa��o usada pela fun��o de amostragem RITAI � exportada
   atrav�s do seguinte bloco comum,
      COMMON / CRITA / XT (NM), PAC (NM), DPAC (NM), A (NM), B (NM), IL (NM), IU (NM),
     1 NPM1
   Onde
     XT (I) ..... pontos da grade, em ordem crescente.
     PAC (I) .... valor da pdf cumulativa em XT (I).
     DPAC (I) ... probabilidade do intervalo I-�simo.
     A (I), B (I) ... par�metros de distribui��o cumulativa inversa racional.
     IL (I) .... maior J para o qual PAC (J) <(I-1) / (NP-1).
     IU (I) .... menor J para o qual PAC (J)> I / (NP-1).
     NPM1 ..... n�mero de pontos de grade menos um, NP-1 (8.LE.NP.LE.NM).
 
      COMUM / CRITAN / CNORM
     CNORM .... constante de normaliza��o do PDF externo, sa�da.
 
      ERR DE DIMENS�O (NM), (NM)
      PAR�METRO (NIP = 51)
      DIMENS�O XS (NIP), YS (NIP), SUMI (NIP)
      PDF EXTERNO
	*/
	static const int NM = 512;
	double EPS = 1.0e-10;
	double ZERO = 1.0e-75;
	double ZEROT = 0.1e0*ZERO;
	static const int NIP = 51;
	
	
	double XS[NIP];
	double YS[NIP];
	double SUMI[NIP];
	double ERR[NM];
	double C[NM];
	
/*	double *XS = (double *) malloc(NIP * sizeof(double));
	double *YS = (double *) malloc(NIP * sizeof(double));
	double *SUMI = (double *) malloc(NIP * sizeof(double));
	double *ERR = (double *) malloc(NM * sizeof(double));
	double *C = (double *) malloc(NM * sizeof(double));*/
	
	
	
	
	int NUNIF;
	double DX;
	double DXI;
	int NP;
	int I1;
	double FX;
	double PDFMAX;
	double CONS;
	double XIH;
	double YSH;
	double FACT;
	double PLI;
	double PUI;
	int ICASE;
	double RR;
	double PAP;
	int LMAX;
	double BIN;
	double PTST;
	double PDFE;
	double XTAU, P1, TAU, CON1, CON2, ETA, P2;
	
	
	
	
	if (N < 9){
		printf(" Error in RITAI0: N must be larger than 8 N=%d11\n", N);
		exit(0);	
	}
	
	if (N > NM) {
		printf("Error in RITAI0: N must be less than NM=512. N=%11d\n", N);
		exit(0);
	}
	
	if (XLOW > XHIGH-EPS) {
		printf("Error in RITAI0: XLOW must be larger than XHIGH. XLOW = %.E6, XHIGH = %.E6\n", XLOW, XHIGH);
		exit(0);
	}
		
	//Come�amos com uma grade de pontos NUNIF uniformemente espa�ados no Intervalo (XLOW, XHIGH).	
	
	
	if (NU >= 0){
		NUNIF=min(max(8,NU),N/2);
        NP=NUNIF;
        
        DX=(XHIGH-XLOW)/(double)(NP-1);
        CRITA_.XT[1-1]=XLOW;
        for (int I = 1; I <= NP-1; I++){
        	CRITA_.XT[I+1-1]=XLOW+I*DX;
		}
		CRITA_.XT[NP-1]=XHIGH;
	} else{
		// Se NU.LT.0, os pontos NUNIF s�o espa�ados logaritmicamente. XLOW deve ser maior ou igual a zero
		NUNIF=min(max(8,-NU),N/2);
        NP=NUNIF;
        if (XLOW < 0.0e0){
        	printf("Error in RITAI0: XLOW and NU are negative. XLOW= %.E7, NU=%11d", XLOW, NU);
			exit(0);
		}
		CRITA_.XT[1-1]=XLOW;
		if (XLOW < 1.0e-16){
		 	 CRITA_.XT[2-1]=XLOW+1.0e-6*(XHIGH-XLOW);
			 I1=2;
		} else{
			I1=1;
		}
		FX=exp(log(XHIGH/CRITA_.XT[I1-1])/(double)(NP-I1));
		for (int I = I1; I <= NP-1; I++){
			CRITA_.XT[I+1-1]=CRITA_.XT[I-1]*FX;
		}
		CRITA_.XT[NP-1]=XHIGH;
	}
	
//	printf("\n\nNP222: %4d\n\n", NP);
	
	for (int I = 1; I <= NP-1; I++){
		DX=CRITA_.XT[I+1-1]-CRITA_.XT[I-1];
        DXI=DX/(double)(NIP-1);
        PDFMAX=0.0e0;
        for (int K = 1; K <= NIP; K++){
  		    XS[K-1]=CRITA_.XT[I-1]+(double)(K-1)*DXI;
  		    switch (PDF){
				case 1: //rndg3f2
					YS[K-1] = fmax(rndg3f2_(XS[K-1]), ZEROT);
					break;
				case 2: //dcsel2
					YS[K-1] = fmax(dcsel2_(XS[K-1]), ZEROT);
					break;	
				case 3: //dcsel2
					YS[K-1] = fmax(graaf22_(XS[K-1]), ZEROT);
					break;	
			}
		//	printf("\n\nYS %d: %f\n\n", K, YS[K-1]);
            PDFMAX=fmax(PDFMAX,YS[K-1]);	
		}
		
		//Integra��o de Simpson.
		CONS=DXI * 3.3333333333333333e-1 * 0.5e0;
		SUMI[1-1]=0.0e0;
		
		for (int K = 2; K <= NIP; K++){
		 	XIH=XS[K-1]-0.5e0*DXI;
		 	switch (PDF){
				case 1: //rndg3f2
					YSH = fmax(rndg3f2_(XIH), ZEROT);
					break;
				case 2: //dcsel2
					YSH = fmax(dcsel2_(XIH), ZEROT);
					break;
				case 3: //dcsel2
					YSH = fmax(graaf22_(XIH), ZEROT);
					break;
			 }
		//	 printf("\n\nYSH %d: %f\n\n", K, YSH);
            PDFMAX=fmax(PDFMAX,YSH);
            SUMI[K-1]=SUMI[K-1-1]+CONS*(YS[K-1-1]+4.0e0*YSH+YS[K-1]);	
		}
		
		CRITA_.DPAC[I-1]=SUMI[NIP-1];
        FACT=1.0e0/CRITA_.DPAC[I-1];
        
        for (int K = 1; K <= NIP; K++){
        	SUMI[K-1]=FACT*SUMI[K-1];
		}
		
		//Quando o PDF desaparece em um dos pontos finais do intervalo, � O valor � modificado.
		
		if (YS[1-1] < ZERO)
			YS[1-1]=1.0e-5*PDFMAX;
		
		if (YS[NIP-1] < ZERO) 
			YS[NIP-1]=1.0e-5*PDFMAX;
		
		PLI=YS[1-1]*FACT;
        PUI=YS[NIP-1]*FACT;
        CRITA_.B[I-1]=1.0e0-1.0e0/(PLI*PUI*DX*DX);
        CRITA_.A[I-1]=(1.0e0/(PLI*DX))-1.0e0-CRITA_.B[I-1];
        C[I-1]=1.0e0+CRITA_.A[I-1]+CRITA_.B[I-1];
        
        if (C[I-1] < ZERO){
        	CRITA_.A[I-1]=0.0e0;
            CRITA_.B[I-1]=0.0e0;
            C[I-1]=1.0e0;	
		}
		
		/*
		ERR (I) � o integral da diferen�a absoluta entre o
		Interpola��o racional e o verdadeiro PDF, estendido ao longo do intervalo
		(XT (I), XT (I + 1)). Calculado usando a regra trapezoidal.
		*/
		
		ICASE=1;

		

L100:;
		ERR[I-1]=0.0e0;
		for (int K = 1; K <= NIP; K++){
			RR=SUMI[K-1];
            PAP=CRITA_.DPAC[I-1]*pow(1.0e0+(CRITA_.A[I-1]+CRITA_.B[I-1]*RR)*RR,2)/
		        ((1.0e0-CRITA_.B[I-1]*RR*RR)*C[I-1]*(CRITA_.XT[I+1-1]-CRITA_.XT[I-1]));
			if ((K == 1) || (K == NIP))
				ERR[I-1]=ERR[I-1]+0.5e0*fabs(PAP-YS[K-1]);
			else
				ERR[I-1]=ERR[I-1]+fabs(PAP-YS[K-1]);	
		}
		
		ERR[I-1]=ERR[I-1]*DXI;
		
		// Se ERR (I) for muito grande, o PDF � aproximado por um uniforme Distribui��o.
		
		if ((ERR[I-1] > 0.10e0*CRITA_.DPAC[I-1]) && (ICASE == 1)){
			CRITA_.B[I-1]=0.0e0;
			CRITA_.A[I-1]=0.0e0;
            C[I-1]=1.0e0;
            ICASE=2;
            goto L100;
		}
	}
	
	CRITA_.XT[NP-1]=XHIGH;
    CRITA_.A[NP-1]=0.0e0;
    CRITA_.B[NP-1]=0.0e0;
    C[NP-1]=0.0e0;
    ERR[NP-1]=0.0e0;
    CRITA_.DPAC[NP-1]=0.0e0;
    
    //Novos pontos de grade s�o adicionados reduzindo � metade o subintervalo com o maior erro absoluto.
    
L200:;

	ERRM=0.0e0;

	LMAX=1;
	
	for (int I = 1; I <= NP-1; I++){
		//ERRM � o maior dos erros de intervalo ERR (I).
		if (ERR[I-1] > ERRM){
			ERRM=ERR[I-1];
            LMAX=I;	
		}
	}
	
	NP=NP+1;
	for (int I = NP; I >= LMAX+1; I--){
		CRITA_.XT[I-1]=CRITA_.XT[I-1-1];
        CRITA_.A[I-1]=CRITA_.A[I-1-1];
        CRITA_.B[I-1]=CRITA_.B[I-1-1];
        C[I-1]=C[I-1-1];
        ERR[I-1]=ERR[I-1-1];
        CRITA_.DPAC[I-1]=CRITA_.DPAC[I-1-1];
	}
	
	CRITA_.XT[LMAX+1-1]=0.5e0*(CRITA_.XT[LMAX-1]+CRITA_.XT[LMAX+2-1]);
	for (int I = LMAX; I <= LMAX+1; I++){
		DX=CRITA_.XT[I+1-1]-CRITA_.XT[I-1];
        DXI=(CRITA_.XT[I+1-1]-CRITA_.XT[I-1])/(double)(NIP-1);
        PDFMAX=0.0e0;
        
        for (int K = 1; K <= NIP; K++){
        	XS[K-1]=CRITA_.XT[I-1]+(double)(K-1)*DXI;
        	switch (PDF){
				case 1: //rndg3f2
					YS[K-1] = fmax(rndg3f2_(XS[K-1]), ZEROT);
					break;
				case 2: //dcsel2
					YS[K-1] = fmax(dcsel2_(XS[K-1]), ZEROT);
					break;		
				case 3: //dcsel2
					YS[K-1] = fmax(graaf22_(XS[K-1]), ZEROT);
					break;
			}
		//	printf("\n\nYS %d: %f\n\n", K, YS[K-1]);
        	PDFMAX=fmax(PDFMAX,YS[K-1]);
		}
		
        //Integra��o de Simpson.
		CONS=DXI*3.3333333333333333e-1*0.5e0;
		SUMI[1-1]=0.0e0;
		
		for (int K = 2; K <= NIP; K++){
		 	XIH=XS[K-1]-0.5e0*DXI;
		 	switch (PDF){
				case 1: //rndg3f2
					YSH = fmax(rndg3f2_(XIH), ZEROT);
					break;
				case 2: //dcsel2
					YSH = fmax(dcsel2_(XIH), ZEROT);
					break;
				case 3: //
					YSH = fmax(graaf22_(XIH), ZEROT);
					break;
			 }
			// printf("\n\nYSH %d: %f\n\n", K, YSH);
            SUMI[K-1]=SUMI[K-1-1]+CONS*(YS[K-1-1]+4.0e0*YSH+YS[K-1]);	
		}
		
		CRITA_.DPAC[I-1]=SUMI[NIP-1];
        FACT=1.0e0/CRITA_.DPAC[I-1];
        
        for (int K = 1; K <= NIP; K++){
        	SUMI[K-1]=FACT*SUMI[K-1];
		}
		
		//Quando o PDF desaparece em um dos pontos finais do intervalo, � O valor � modificado.
		
		if (YS[1-1] < ZERO)
			YS[1-1]=1.0e-5*PDFMAX;
		
		if (YS[NIP-1] < ZERO) 
			YS[NIP-1]=1.0e-5*PDFMAX;
		
		PLI=YS[1-1]*FACT;
        PUI=YS[NIP-1]*FACT;
        CRITA_.B[I-1]=1.0e0-1.0e0/(PLI*PUI*DX*DX);
        CRITA_.A[I-1]=(1.0e0/(PLI*DX))-1.0e0-CRITA_.B[I-1];
        C[I-1]=1.0e0+CRITA_.A[I-1]+CRITA_.B[I-1];
        
        if (C[I-1] < ZERO){
        	CRITA_.A[I-1]=0.0e0;
            CRITA_.B[I-1]=0.0e0;
            C[I-1]=1.0e0;	
		}
		
		/*
		ERR (I) � o integral da diferen�a absoluta entre o
		Interpola��o racional e o verdadeiro PDF, estendido ao longo do intervalo
		(XT (I), XT (I + 1)). Calculado usando a regra trapezoidal.
		*/
		
		ICASE=1;
		
L300:;
		
		ERR[I-1]=0.0e0;
		for (int K = 1; K <= NIP; K++){
			RR=SUMI[K-1];
            PAP=CRITA_.DPAC[I-1]*pow(1.0e0+(CRITA_.A[I-1]+CRITA_.B[I-1]*RR)*RR,2)/
		        ((1.0e0-CRITA_.B[I-1]*RR*RR)*C[I-1]*(CRITA_.XT[I+1-1]-CRITA_.XT[I-1]));
			if ((K == 1) || (K == NIP))
				ERR[I-1]=ERR[I-1]+0.5e0*fabs(PAP-YS[K-1]);
			else
				ERR[I-1]=ERR[I-1]+fabs(PAP-YS[K-1]);	
		}
		
		ERR[I-1]=ERR[I-1]*DXI;
		
		// Se ERR (I) for muito grande, o PDF � aproximado por um uniforme Distribui��o.
		
		if ((ERR[I-1] > 0.10e0*CRITA_.DPAC[I-1]) && (ICASE == 1)){
			CRITA_.B[I-1]=0.0e0;
			CRITA_.A[I-1]=0.0e0;
            C[I-1]=1.0e0;
            ICASE=2;
            goto L300;
		}	
	}
	
	if (NP < N)
		goto L200;
	*CRITA_.NPM1=NP-1;
	
	//renormalizacao
	
	*CRITAN_.CNORM=0.0e0;
	for (int I = 1; I <= *CRITA_.NPM1; I++){
		*CRITAN_.CNORM = *CRITAN_.CNORM + CRITA_.DPAC[I-1];	
	}
	
	*CRITAN_.CNORM=1.0e0 / *CRITAN_.CNORM;
	ERRM=0.0e0;
	
	for (int I = 1; I <= *CRITA_.NPM1; I++){
		CRITA_.DPAC[I-1]=CRITA_.DPAC[I-1] * *CRITAN_.CNORM;
        ERR[I-1]=ERR[I-1] * *CRITAN_.CNORM;
        ERRM=fmax(ERRM,ERR[I-1]);
	}	
	
	CRITA_.PAC[1-1]=0.0e0;
	
	
	for (int I = 1; I <= *CRITA_.NPM1; I++){
		CRITA_.PAC[I+1-1]=CRITA_.PAC[I-1]+CRITA_.DPAC[I-1];
	}
	CRITA_.PAC[NP-1]=1.0e0;
	
	//Limites pr�-calculados para a pesquisa bin�ria inicial em Sub-rotina RITAI.
	
	BIN= 1.0e0 / (double)*CRITA_.NPM1;
	CRITA_.IL[1-1]=1;
	
	for (int I = 2; I <= *CRITA_.NPM1; I++){
		PTST= double(I-1) * BIN;
		for (int J = CRITA_.IL[I-1-1]; J <= NP; J++){
			if (CRITA_.PAC[J-1] > PTST){
				CRITA_.IL[I-1]=J-1;
                CRITA_.IU[I-1-1]=J;
                goto L400;
			}
		}
L400:;
	}
	CRITA_.IU[*CRITA_.NPM1-1]=NP;
    CRITA_.IL[NP-1]=NP-1;
    CRITA_.IU[NP-1]=NP;

	

    
    
    //Imprimir tabelas de interpola��o (somente quando IWR.GT.0).
    
    FILE* IW = fopen("param2.dat", "w");
	if (IW == NULL){
		printf("N�o foi possivel abrir o arquivo param2.dat");
		exit(0);
	}
	
	fprintf(IW, "#  Interpolation error = %.7E\n", ERRM );
	fprintf(IW, "# Normalising constant = %.7E\n", *CRITAN_.CNORM);
	fprintf(IW, " #     X           PDF(X)          A             B             C           error\n");

	
	for (int I =1; I <= *CRITA_.NPM1; I++){
		switch (PDF){
			case 1: //rndg3f2
				PDFE = fmax(rndg3f2_(CRITA_.XT[I-1]), ZEROT)* *CRITAN_.CNORM;
				break;
			case 2: //dcsel2
				PDFE = fmax(dcsel2_(CRITA_.XT[I-1]), ZEROT)* *CRITAN_.CNORM;
				break;	
			case 3: //dcsel2
				PDFE = fmax(graaf22_(CRITA_.XT[I-1]), ZEROT)* *CRITAN_.CNORM;
				break;
		}
	
		fprintf(IW, "%f  %.6E %f  %f  %f  %f\n", CRITA_.XT[I-1], PDFE, CRITA_.A[I-1], CRITA_.B[I-1], C[I-1], ERR[I-1] );
	}
	fclose(IW);
	
	
	FILE *IW1 = fopen("table.dat", "w");
	if (IW == NULL){
		printf("N�o foi possivel abrir o arquivo table2.dat");
		exit(0);
	}
	
	fprintf(IW1, "#  Interpolation error = %.7E\n", ERRM );
	fprintf(IW1, "# Normalising constant = %.7E\n", *CRITAN_.CNORM);
	fprintf(IW1, " #      X             PDF_ex          PDF_ap           err\n");
	for (int I =1; I <= *CRITA_.NPM1; I++){
		DX = (CRITA_.XT[I+1-1] - CRITA_.XT[I-1]) / (double)(NIP-1);
		for (int K = 1; K <= NIP; K = K+5){
			XTAU = CRITA_.XT[I-1] + (K-1)*DX;
			switch (PDF){
				case 1: //rndg3f2
					P1 = fmax(rndg3f2_(XTAU) * *CRITAN_.CNORM, ZEROT);
					break;
				case 2: //dcsel2
					P1 = fmax(dcsel2_(XTAU) * *CRITAN_.CNORM, ZEROT);
					break;
				case 3: //dcsel2
					P1 = fmax(graaf22_(XTAU) * *CRITAN_.CNORM, ZEROT);
					break;
			}
			//Interpolacao racional
			TAU = (XTAU - CRITA_.XT[I-1]) / (CRITA_.XT[I+1-1] - CRITA_.XT[I-1]);
			CON1 = 2.0e0*CRITA_.B[I-1]*TAU;
			CON2 = C[I-1] - CRITA_.A[I-1]*TAU;
			if (fabs(CON1) > 1.0e-10*fabs(CON2)){
				ETA = CON2*(1.0e0 - sqrt(1.0e0-2.0e0*TAU*CON1/pow(CON2,2)))/CON1;
			} else{
				ETA = TAU / CON2;
			}
			
			P2 = CRITA_.DPAC[I-1] * pow(1.0e0+ (CRITA_.A[I-1] + CRITA_.B[I-1]*ETA) * ETA, 2) /
				 ((1.0e0-CRITA_.B[I-1] * ETA *ETA) * C[I-1] * (CRITA_.XT[I+1-1] - CRITA_.XT[I-1]));
			fprintf(IW1, "%f   %.8E   %f   %f   \n", XTAU, P1, P2, (P1-P2)/P1);	 
		}
	}
	
	fclose(IW1);
	
	FILE* IW2 = fopen("limits.dat", "w");
	if (IW2 == NULL){
		printf("N�o foi possivel abrir o arquivo limits2.dat");
		exit(0);
	}
	
	fprintf(IW2, " #   I      PAC(ITL)         (I-1)/NPM1           I/NPM1            PAC(ITU)            PAC(I)\n");
	for (int I = 1; I <= *CRITA_.NPM1; I++){
		fprintf(IW2, "%5d  %f   %.11E   %f   %f   %f\n", I, CRITA_.PAC[CRITA_.IL[I-1] -1], (I-1)*BIN, I*BIN, CRITA_.PAC[CRITA_.IU[I-1]], CRITA_.PAC[I-1]);
		
		if (((CRITA_.PAC[CRITA_.IL[I-1]-1]) > (I-1)*BIN+EPS) || (CRITA_.PAC[CRITA_.IU[I-1]-1] < I*BIN-EPS)){
			fprintf(IW2, "WARNING: The first four values should be in creasing order.\n");
		}
	}
	 
	 fclose(IW2);
	 
//	 free(XS);
//	 free(YS);
//	 free(SUMI);
//	 free(ERR);
//	 free(C);

	 
	 
	
//	 printf("\nRITAI02\n");
	
}



 


void relaxr2_(FILE *IRD, FILE *IWR, int *INFO){


	
	/*
	
	Esta sub-rotina l� dados de relaxa��o at�mica, para um �nico elemento,
  a partir do arquivo de defini��o de material (unidade IRD) e inicializa o
  algoritmo para amostragem aleat�ria de cascatas de desexcita��o deste
  elemento. (Veja os coment�rios do cabe�alho na sub-rotina RELAX).
	
	*/
	static const int NTRAN = 2500;
	static const int NRX = 60000;
	
	int *JS0 = (int *) malloc(NTRAN * sizeof(int));
    int *JS1 = (int *) malloc(NTRAN * sizeof(int));
    int *JS2 = (int *) malloc(NTRAN * sizeof(int));
    int *KK = (int *) malloc(NTRAN * sizeof(int));
    int *IORD = (int *) malloc(NTRAN * sizeof(int));
    int *ISR = (int *) malloc(NTRAN * sizeof(int));
    double *PR = (double *) malloc(NTRAN * sizeof(double));
	double *ER = (double *) malloc(NTRAN * sizeof(double));
	double *WW = (double *) malloc(NTRAN * sizeof(double));
	double *FF = (double *) malloc(NTRAN * sizeof(double));

		
	
	int IQQ[30];
	double EE[30];
	double ALWR[30];
	double CP0P[30];
	double ALTIME;
	
	char LSHELL[30][3] = { {"K"},{"L1"},{"L2"},{"L3"},{"M1"},{"M2"},{"M3"},{"M4"},{"M5"},{"N1"},
       	   	   	   	   	   {"N2"},{"N3"},{"N4"},{"N5"},{"N6"},{"N7"},{"O1"},{"O2"},{"O3"},{"O4"},
      	  	  	  	  	   {"O5"},{"O6"},{"O7"},{"P1"},{"P2"},{"P3"},{"P4"},{"P5"},{"Q1"},{"X"} }; 
	
	
	
	
	char CH5[6];
	
	int IZ, NSHR, NT;
	
	*CRELAX_.MODER = 1; //RELAXE o modo de opera��o normal.
	
	//dados de transicao de entrada
	
	fgets(LINHA, sizeof(LINHA), IRD);
	extrairString(APOIO, LINHA, 16, 3);
	IZ = atoi(APOIO);
	extrairString(APOIO, LINHA, 37, 3);
	NSHR = atoi(APOIO);
	extrairString(APOIO, LINHA,63, 5);
	NT = atoi(APOIO);
	

	
	if (*INFO >= 2){
		fprintf(IWR, "\n *** RELAX:  Z =%3d, no. of shells = %3d, no. of transitions = %3d\n", IZ, NSHR, NT);
	}
	if (NT > NTRAN){
		printf("RELAXR. NTRAN needs to be increased.\n");
		exit(0);
	}
	if (*CRELAX_.NCUR + NT > NRX){
		fprintf(IWR, "Insufficient memory storage in RELAXR.\n");
		fprintf(IWR, "Increase the value of the parameter NRX to %d.\n",*CRELAX_.NCUR + NT );
		printf("Insufficient memory storage in RELAXR.\n");
		printf("Increase the value of the parameter NRX to %d.\n",*CRELAX_.NCUR + NT );
		exit(0);
	}
	
	if (*INFO >= 2){
		fprintf(IWR, "\n  i  Shell    f    Ui (eV)    Gamma(1/eV)  lifetime (s)    Ji(0)\n");
     	fprintf(IWR, "-----------------------------------------------------------------\n");
	}
	
	for (int IS = 1; IS <= NSHR; IS++){
		fgets(LINHA, sizeof(LINHA), IRD);
	   	extrairString(APOIO, LINHA, 1, 3);
	   	ISR[IS-1] = atoi(APOIO);
	   	extrairString(APOIO, LINHA, 5, 5);
	   	strcpy(CH5, APOIO);
	   	CH5[5] = '\0';
	   	extrairString(APOIO, LINHA, 11, 1);
	   	IQQ[IS-1] = atoi(APOIO);
	   	extrairString(APOIO, LINHA, 12, 13);
	   	EE[IS-1] = atof(APOIO);
	   	extrairString(APOIO, LINHA, 25, 13);
	   	ALWR[IS-1] = atof(APOIO);
	   	extrairString(APOIO, LINHA, 38, 13);
	   	CP0P[IS-1] = atof(APOIO);
	   	
	   	if (ALWR[IS-1] > 0.0)
	   		ALTIME = HBAR / ALWR[IS-1];
		else
			ALTIME = 0.0e0;
		
		if (*INFO >= 2)
			fprintf(IWR, "%3d %s %s %1d %.5E %.5E %.5E %.5E\n",ISR[IS-1], CH5, LSHELL[ISR[IS-1]-1], IQQ[IS-1], EE[IS-1], ALWR[IS-1],ALTIME, CP0P[IS-1]);	
	}
 
	if (NT > 0){
		if (*INFO >= 2){
			fprintf(IWR, "\n S0 S1 S2   Probability     Energy (eV)\n");
		    fprintf(IWR, " ------------------------------------------ \n");
		}
		for (int I = 1; I <= NT; I++){
			fgets(LINHA, sizeof(LINHA), IRD);
			
	 	  	extrairString(APOIO, LINHA, 1, 3);
		   	JS0[I-1] = atoi(APOIO);
		   	extrairString(APOIO, LINHA, 4, 3);
		    JS1[I-1] = atoi(APOIO);
		   	extrairString(APOIO, LINHA, 7, 3);
		   	JS2[I-1] = atoi(APOIO);
		   	extrairString(APOIO, LINHA, 10, 13);
		   	PR[I-1] = atof(APOIO);
		   	extrairString(APOIO, LINHA, 23, 13);
		   	ER[I-1] = atof(APOIO);	
		   	
		   	if (*INFO >= 2){
				   	fprintf(IWR, " %s %s %s %.5E %.5E\n",LSHELL[JS0[I-1]-1], LSHELL[JS1[I-1]-1], LSHELL[JS2[I-1]-1], PR[I-1],ER[I-1]);	
			}
			   
			if (PR[I-1] < 1.0e-35){
				if (*INFO > 2) 
					fprintf(IWR, " %s %s %s %.5E %.5E\n",LSHELL[JS0[I-1]-1], LSHELL[JS1[I-1]-1], LSHELL[JS2[I-1]-1], PR[I-1],ER[I-1]);	
				printf("RELAXR. Negative transition probability?\n");
				exit(0);
			}		
		}
	}

	
	//Verifique se os dados deste elemento j� foram carregados.

	
	if (CRELAX_.IFIRST[1-1][IZ-1] != 0){
		free(JS0);
	free(JS1);
	free(JS2);
	free(KK);
	free(IORD);
	free(ISR);
	free(PR);
	free(ER);
	free(WW);
	free(FF);
		
		return;
	}
		
 	CADATA_.NSHT[IZ-1] = NSHR;
 	
 	for (int IS = 1; IS <= NSHR; IS++){
	    CADATA_.IKS[IS-1][IZ-1] = ISR[IS-1];
 		CADATA_.IFI[ISR[IS-1]-1][IZ-1] = IQQ[IS-1];
 		CADATA_.EB[ISR[IS-1]-1][IZ-1] = EE[IS-1];
 		if (ALWR[IS-1] > 0.0e0)
 			CADATA_.ALW[ISR[IS-1]-1][IZ-1] = HBAR/ALWR[IS-1];
 		else
 			CADATA_.ALW[ISR[IS-1]-1][IZ-1] = 0.0e0;
 		CADATA_.CP0[ISR[IS-1]-1][IZ-1] = CP0P[IS-1];
	 }
	  
	 if (NT == 0){
		 for (int IS = 1; IS <= 16; IS++){
			CRELAX_.IFIRST[IS-1][IZ-1] = *CRELAX_.NCUR+1;
			CRELAX_.ILAST[IS-1][IZ-1] = *CRELAX_.NCUR+1;
	//A matriz IS0 cont�m os valores de alias.
			CRELAX_.IS0[*CRELAX_.NCUR+1-1] = *CRELAX_.NCUR+1;
			CRELAX_.P[*CRELAX_.NCUR+1-1] = 1.0e0;
			CRELAX_.F[*CRELAX_.NCUR+1-1] = 1.0e0;
			CRELAX_.ET[*CRELAX_.NCUR+1-1] = 0.0e0;
			CRELAX_.IS1[*CRELAX_.NCUR+1-1] = 1;
			CRELAX_.IS2[*CRELAX_.NCUR+1-1] = 1;
			*CRELAX_.NCUR = *CRELAX_.NCUR + 1;	
		 }
		 free(JS0);
	free(JS1);
	free(JS2);
	free(KK);
	free(IORD);
	free(ISR);
	free(PR);
	free(ER);
	free(WW);
	free(FF);
		 	
		 return;
	 }
	 //Walker's aliasing.
	 
	for (int IS = 1; IS <= 16; IS++){
	 	int N = 0;
	 	for (int J = 1; J <= NT; J++){
	 		if (JS0[J-1] == IS){
				 N = N + 1;
				 IORD[N-1] = J;
				 WW[N-1] = PR[J-1];	 
			 }
		 }
		 
		 if (N > 1){
		 
			 irnd02_(WW, FF, KK, &N);
			 
			 CRELAX_.IFIRST[IS-1][IZ-1] = *CRELAX_.NCUR +1;
			 CRELAX_.ILAST[IS-1][IZ-1] = *CRELAX_.NCUR + N;
			 for (int L = 1; L <= N; L++){
			 
			   	CRELAX_.P[*CRELAX_.NCUR+L-1] = WW[L-1];
				CRELAX_.F[*CRELAX_.NCUR+L-1] = FF[L-1];
				CRELAX_.ET[*CRELAX_.NCUR+L-1] = ER[IORD[L-1]-1];
				//The array IS0 contains the alias values.
				CRELAX_.IS0[*CRELAX_.NCUR+L-1] = CRELAX_.IFIRST[IS-1][IZ-1]+ KK[L-1] - 1;
				CRELAX_.IS1[*CRELAX_.NCUR+L-1] = JS1[IORD[L-1]-1];
				CRELAX_.IS2[*CRELAX_.NCUR+L-1] = JS2[IORD[L-1]-1];
				
				 
			 }
			 *CRELAX_.NCUR = *CRELAX_.NCUR + N;
		 }
		 else{
		 	*CRELAX_.NCUR = *CRELAX_.NCUR + 1;
		 	CRELAX_.IFIRST[IS-1][IZ-1] = *CRELAX_.NCUR;
			CRELAX_.ILAST[IS-1][IZ-1] = *CRELAX_.NCUR;
			CRELAX_.IS0[*CRELAX_.NCUR+1-1] = *CRELAX_.NCUR;;
			CRELAX_.P[*CRELAX_.NCUR-1] = 1.0e0;
			CRELAX_.F[*CRELAX_.NCUR-1] = 1.0e0;
			CRELAX_.ET[*CRELAX_.NCUR-1] = ER[1-1];
			CRELAX_.IS1[*CRELAX_.NCUR-1] = JS1[1-1];
			CRELAX_.IS2[*CRELAX_.NCUR-1] = JS2[1-1];	 
		 } 
	}
	
	//Verifique se as probabilidades de transi��o s�o reproduzidas corretamente.
	
	double TST = 0.0;
	
	for (int IS = 1; IS <= 16; IS++){
		double I0 = CRELAX_.IFIRST[IS-1][IZ-1];
		double IN = CRELAX_.ILAST[IS-1][IZ-1];
		double PT = 0.0;
		for (int I = I0; I <= IN; I++){
			PT = PT + CRELAX_.P[I-1];
		}
		
		for (int I = I0; I <= IN; I++){
			double PPI = 0.0;
			for (int J = I0; J<= IN; J++){
				if (CRELAX_.IS0[J-1] == I) 
					PPI = PPI + (1.0 - CRELAX_.F[J-1]);
			}
			PPI = (PPI + CRELAX_.F[I-1]) / (IN - I0 + 1);
			TST = fmax(TST, fabs(1.0 - PPI * PT/CRELAX_.P[I-1]));	
		}
	}
	if (TST > 1.0e-12) {
		printf("RELAXR. Rounding error is too large.");
		exit(0);
	}
	
	free(JS0);
	free(JS1);
	free(JS2);
	free(KK);
	free(IORD);
	free(ISR);
	free(PR);
	free(ER);
	free(WW);
	free(FF);
	
//	printf("\n\nRELAXR2\n\n");
	
		

}


void esiar2_(int *M, FILE *IRD, FILE *IWR, int *INFO){
	
	/*
	Esta sub-rotina l� se��es transversais para ioniza��o de camada interna por
	Impacto do el�tron dos elementos no material M e prepara simula��o Tabelas.
	*/

	static const int NDIN=850;
	static const int NRP=8000;
	
	double *E = (double *) malloc(NDIN * sizeof(double));
	double *X = (double *) malloc(NDIN * sizeof(double));
	double *Y = (double *) malloc(NDIN * sizeof(double));
	
	double XESIR[16][NDIN];
	int IZZ;
	int NSHR;
	int NDATA;	
	int IC;
	double XC;
	int J;
	int NSHA;
	int N;
	double DX;
	char CS5[16][6] = {{"CS-K "},{"CS-L1"},{"CS-L2"},{"CS-L3"},{"CS-M1"},{"CS-M2"},
	   	   	   	   	   {"CS-M3"},{"CS-M4"},{"CS-M5"},{"CS-N1"},{"CS-N2"},{"CS-N3"},
                       {"CS-N4"},{"CS-N5"},{"CS-N6"},{"CS-N7"}};
		
	//Leia as tabelas da se��o x do elemento
	
	for (int IEL = 1; IEL <= COMPOS_.NELEM[*M-1]; IEL++){
		fgets(LINHA, sizeof(LINHA), IRD);
	   	extrairString(APOIO, LINHA, 46, 3);
	    IZZ = atoi(APOIO);
	   	extrairString(APOIO, LINHA, 60, 3);
	   	NSHR = atoi(APOIO);
	   	extrairString(APOIO, LINHA, 73, 4);
	   	NDATA = atoi(APOIO);
	   	
	   	if (*INFO >= 2 )
	   		fprintf(IWR, "\n*** Electron impact ionisation cross sections, ,  IZ= %3d,  NSHELL =%3d,  NDATA=%4d\n", IZZ, NSHR, NDATA);
	   	if (IZZ != COMPOS_.IZ[IEL-1][*M-1]){
	   		printf("ESIaR. Corrupt material data file.\n");
	   		exit(0);
		}
		if (NDATA > NDIN){
	   		printf("ESIaR. Corrupt material data fileToo many data points.\n");
	   		exit(0);
		}
		if (NSHR > 16){
	   		printf("ESIaR. Too many shells.\n");
	   		exit(0);
		}
		
		for (int IE = 1; IE <= NDATA; IE++){
			fgets(LINHA, sizeof(LINHA), IRD);
            extrairString(APOIO, LINHA, 0, 12);
  	   	    E[IE-1] = atof(APOIO);   	
   	   	    for (int IS = 1; IS <= NSHR; IS++ ){
	   		   	extrairString(APOIO, LINHA, 12*IS, 12);
	   			XESIR[IS-1][IE-1] = atof(APOIO);
			}
		}
		
		//Remova os inv�lucros com energias de ioniza��o inferiores a 50 eV.
		
		if (NSHR > 1){
			NSHA = NSHR;
			for (int IS = NSHA; IS >= 1; IS--){
				if (CADATA_.EB[IS-1][IZZ-1] < 50.0e0)
					NSHR = NSHR-1;
				else goto L1;
			}
			if (NSHR < 1)
				NSHR = 1;
		}
L1:;	
		if (*INFO >= 2){
			fprintf(IWR, "\n   Energy        ");
			for (int IS = 1; IS <= NSHR; IS++){
				fprintf(IWR, "%s       ", CS5[IS-1]);
			}
			fprintf(IWR, "\n");
			for (int IE = 1; IE <= NDATA; IE++){
				double TCS = 0.0;
				for (int IS = 1; IS <= NSHR; IS++){
					TCS = TCS + XESIR[IS-1][IE-1];
				}
				fprintf(IWR, " %.5E", E[IE-1]);
				for (int IS = 1; IS <= NSHR; IS++){ 
					fprintf(IWR, " %.5E", XESIR[IS-1][IE-1]);
				}
				fprintf(IWR, " %.5E\n",TCS);	
			}
		}
		    
		CESI0_.NSESI[IZZ-1] = NSHR;
		if (CESI0_.IESIF[IZZ-1] == 0){
			CESI0_.IESIF[IZZ-1] = *CESI0_.NCURE + 1;
			if (*CESI0_.NCURE + NEGP > NRP){
				fprintf(IWR, "Insufficient memory storage in ESIaR.\n");
				fprintf(IWR, "Increase the value of the parameter NRP to %d.\n", (*CESI0_.NCURE + NEGP));
				printf( "Insufficient memory storage in ESIaR.\n");
				printf( "Increase the value of the parameter NRP to %d.\n", (*CESI0_.NCURE + NEGP));
				exit(0);
			}
			
			for (int IS = 1; IS <= NSHR; IS++){
			    N = 0;
				for (int I = 1; I <= NDATA; I++){
					if (XESIR[IS-1][I-1] > 1.0e-35){
						N = N+1;
						X[N-1] = log(E[I-1]);
						if (N > 1) 
							X[N-1] = max(X[N-1], X[N-1-1] + 1.0e-6);
						Y[N-1] = log(XESIR[IS-1][I-1]);
					}
				}
				if (N > 4){
					for (int I = 1; I <= NEGP; I++){
						IC = *CESI0_.NCURE + I;
						XC = CEGRID_.DLEMP[I-1];
						if (XC > X[1-1]){
							findi2_(X, &XC, &N, &J);		
							if (J == N) 
								J = N -1;
							DX = X[J+1-1] - X[J-1];
							if (DX > 1.0e-6){
								CESI0_.XESI[IS-1][IC-1] = Y[J-1] + (XC - X[J-1]) * ( Y[J+1-1] - Y[J-1]) / DX;
								
							} else{
								CESI0_.XESI[IS-1][IC-1] = (Y[J+1-1] + Y[J-1]) / 2.0e0;
								
							}
							
						} else{
							CESI0_.XESI[IS-1][IC-1] = -80.6e0;
						}
					}
				} else{
					for (int I = 1; I <= NEGP; I++){
						IC = *CESI0_.NCURE + I;
						CESI0_.XESI[IS-1][IC-1] = -80.6;	
					}
				}
				
			}
			*CESI0_.NCURE = *CESI0_.NCURE + NEGP;
			CESI0_.IESIL[IZZ-1] = *CESI0_.NCURE;	
		}
	}
	
	free(E);
	free(X);
	free(Y);
	
//	printf("\n\nESIAR2\n\n");
	
}


void findi2_(double *X, double *XC, int *N, int *I){
	
	/*
	Esta sub-rotina encontra o intervalo (X (I), X (I + 1)) que cont�m o
  valor XC usando pesquisa bin�ria.

  Entrada:
     X (I) (I = 1: N) ... pontos de grade (os valores de X devem ser crescentes
                      pedido).
     XC ............. ponto a ser localizado.
     N ............. n�mero de pontos da grade.
  Sa�da:
     I .............. �ndice de intervalo.
	*/
	
	if (*XC > X[*N-1]){
		*I = *N;
		return;
	}
	
	if (*XC < X[1-1]){
		*I = 1;
		return;
	}
	
	*I = 1;
	int I1 = *N;
L1:;
	int IT = (*I + I1)/2;
	if (*XC > X[IT-1]){
		*I = IT;
	} else{
		I1 = IT;
	}
	if ((I1 - *I) > 1)
		goto L1;
	
//	printf("\n\nFINDI2\n\n");	
	

}
	

void psiar2_(int *M, FILE *IRD, FILE *IWR, int *INFO){
	/*
	Esta sub-rotina l� se��es transversais para ioniza��o de camada interna por
   impacto p�sitron dos elementos no material M e prepara simula��o Tabelas.
	*/
	
	static const int NDIN=800;
	static const int NRP=8000;
	
		double *E = (double *) malloc(NDIN * sizeof(double));
	double *X = (double *) malloc(NDIN * sizeof(double));
	double *Y = (double *) malloc(NDIN * sizeof(double));
	
	double XPSIR[16][NDIN];
	int IZZ;
	int NSHR;
	int NDATA;	
	int IC;
	double XC;
	int J;
	int NSHA;
	double DX;
	
	
	char CS5[16][6] = {{"CS-K "},{"CS-L1"},{"CS-L2"},{"CS-L3"},{"CS-M1"},{"CS-M2"},
	   	   	   	   	   {"CS-M3"},{"CS-M4"},{"CS-M5"},{"CS-N1"},{"CS-N2"},{"CS-N3"},
                       {"CS-N4"},{"CS-N5"},{"CS-N6"},{"CS-N7"}};
	
	//Leia as tabelas da se��o x do elemento
	
	for (int IEL = 1; IEL <= COMPOS_.NELEM[*M-1]; IEL++){
		fgets(LINHA, sizeof(LINHA), IRD);
	   	extrairString(APOIO, LINHA, 46, 3);
	    IZZ = atoi(APOIO);
	   	extrairString(APOIO, LINHA, 60, 3);
	   	NSHR = atoi(APOIO);
	   	extrairString(APOIO, LINHA, 73, 4);
	   	NDATA = atoi(APOIO);
	   	
	   	if (*INFO >= 2 )
	   		fprintf(IWR, "\n*** Positron impact ionisation cross sections,  IZ= %3d,  NSHELL =%3d,  NDATA=%4d\n", IZZ, NSHR, NDATA);
	   	if (IZZ != COMPOS_.IZ[IEL-1][*M-1]){
	   		printf("PSIaR. Corrupt material data file.\n");
	   		exit(0);
		}
		if (NDATA > NDIN){
	   		printf("PSIaR. Corrupt material data file. Too many data points.\n");
	   		exit(0);
		}
		if (NSHR > 16){
	   		printf("PSIaR. Too many shells.\n");
	   		exit(0);
		}
		
		
		for (int IE = 1; IE <= NDATA; IE++){
			fgets(LINHA, sizeof(LINHA), IRD);
            extrairString(APOIO, LINHA, 0, 12);
  	   	    E[IE-1] = atof(APOIO);   	
   	   	    for (int IS = 1; IS <= NSHR; IS++ ){
	   		   	extrairString(APOIO, LINHA, 12*IS, 12);
	   			XPSIR[IS-1][IE-1] = atof(APOIO);
			}
		}
		
		
		//Remova os inv�lucros com energias de ioniza��o inferiores a 50 eV.
		
		if (NSHR > 1){
			NSHA = NSHR;
			for (int IS = NSHA; IS >= 1; IS--){
				if (CADATA_.EB[IS-1][IZZ-1] < 50.0e0)
					NSHR = NSHR-1;
				else goto L1;
			}
			if (NSHR < 1)
				NSHR = 1;
		}
L1:;	

		if (*INFO >= 2){
			fprintf(IWR, "\n   Energy        ");
			for (int IS = 1; IS <= NSHR; IS++){
				fprintf(IWR, "%s       ", CS5[IS-1]);
			}
			fprintf(IWR, "\n");
			for (int IE = 1; IE <= NDATA; IE++){
				double TCS = 0.0;
				for (int IS = 1; IS <= NSHR; IS++){
					TCS = TCS + XPSIR[IS-1][IE-1];
				}
				fprintf(IWR, " %.5E", E[IE-1]);
				for (int IS = 1; IS <= NSHR; IS++){ 
					fprintf(IWR, " %.5E", XPSIR[IS-1][IE-1]);
				}
				fprintf(IWR, " %.5E\n",TCS);	
			}
		}
		
		CPSI0_.NSPSI[IZZ-1] = NSHR;
		if (CPSI0_.IPSIF[IZZ-1] == 0){
			CPSI0_.IPSIF[IZZ-1] = *CPSI0_.NCURP + 1;
			if (*CPSI0_.NCURP + NEGP > NRP){
				fprintf(IWR, "Insufficient memory storage in PSIaR.\n");
				fprintf(IWR, "Increase the value of the parameter NRP to %d.\n", (*CPSI0_.NCURP + NEGP));
				printf( "Insufficient memory storage in PSIaR.\n");
				printf( "Increase the value of the parameter NRP to %d.\n", (*CPSI0_.NCURP + NEGP));
				exit(0);
			}
			
			for (int IS = 1; IS <= NSHR; IS++){
				int N = 0;
				for (int I = 1; I <= NDATA; I++){
					if (XPSIR[IS-1][I-1] > 1.0e-35){
						N = N+1;
						X[N-1] = log(E[I-1]);
						if (N > 1) 
							X[N-1] = max(X[N-1], X[N-1-1] + 1.0e-6);
						Y[N-1] = log(XPSIR[IS-1][I-1]);
					}
				}
				if (N > 4){
					for (int I = 1; I <= NEGP; I++){
						IC = *CPSI0_.NCURP + I;
						XC = CEGRID_.DLEMP[I-1];
						if (XC > X[1-1]){
							findi2_(X, &XC, &N, &J);		
							if (J == N) 
								J = N -1;
							DX = X[J+1-1] - X[J-1];
							if (DX > 1.0e-6){
								CPSI0_.XPSI[IS-1][IC-1] = Y[J-1] + (XC - X[J-1]) * (Y[J+1-1] - Y[J-1]) / DX;
							} else{
								CPSI0_.XPSI[IS-1][IC-1] = (Y[J+1-1] + Y[J-1]) / 2.0e0;	
							}
						} else{
							CPSI0_.XPSI[IS-1][IC-1] = -80.6;
						}
					}
				} else{
					for (int I = 1; I <= NEGP; I++){
						IC = *CPSI0_.NCURP + I;
						CPSI0_.XPSI[IS-1][IC-1] = -80.6;	
					}
				}
				
			}
			*CPSI0_.NCURP = *CPSI0_.NCURP + NEGP;
			CPSI0_.IPSIL[IZZ-1] = *CPSI0_.NCURP;	
		}
	}
	
		free(E);
	free(X);
	free(Y);
//	printf("\n\nPSIAR2\n\n");	
	return;	
	
}


void ebrar2_(double *WCRM, int *M, FILE *IRD, FILE *IWR, int *INFO){
	
 /*
	Esta sub-rotina l� a se��o transversal em escala de brems para el�trons
  no material M do arquivo de dados do material. Ele calcula restrito
  se��es transversais integradas e inicializa o algoritmo para simula��o
  ��o da emiss�o de brems por el�trons e p�sitrons.
 */	
 
//	ebrar_(0,0,0,0,0);

    const static int NEGP = 200;
    const static int NBE = 57;
    const static int NBW = 32;

    double *A = (double *) malloc(NEGP * sizeof(double));
	double *B = (double *) malloc(NEGP * sizeof(double));
	double *C = (double *) malloc(NEGP * sizeof(double));
	double *D = (double *) malloc(NEGP * sizeof(double));
	double *PAC = (double *) malloc(NEGP * sizeof(double));
	double *PDF = (double *) malloc(NEGP * sizeof(double));

 
//	double A[NEGP], B[NEGP], C[NEGP], D[NEGP], PAC[NEGP], PDF[NEGP];
	
	
	double WB0[NBW] = {1.0e-12, 0.025e0,0.05e0,0.075e0,0.1e0,0.15e0,0.2e0,0.25e0,
 	 	 	 	        0.3e0,0.35e0,0.4e0,0.45e0,0.5e0,0.55e0,0.6e0,0.65e0,0.7e0,
 	 	 	 	        0.75e0,0.8e0,0.85e0,0.9e0,0.925e0,0.95e0,0.97e0,0.99e0,
 	 	 	 	        0.995e0,0.999e0,0.9995e0,0.9999e0,0.99995e0,0.99999e0,1.0e0};

							
    double ZBR;
    double ELL;
    int J;
    int NBER;
    
    
							
							
 	//Lendo a tabela de se��o transversal em escala.
 	
	fgets(LINHA, sizeof(LINHA), IRD);
	extrairString(APOIO, LINHA, 45, 12);
	ZBR = atof(APOIO);
	extrairString(APOIO, LINHA, 67, 4);
	NBER = atoi(APOIO);
	
	if (*INFO >= 2) 
		fprintf(IWR, "\n*** Electron scaled bremss x-section,  ZEQ = %.5E,   NDATA =%4d\n", ZBR, NBER);
	
	if (NBER != NBE) {
		printf("EBRR. Inconsistent format.\n");
		exit(0);
	}
	CEBR_.ZBR2[*M-1] = ZBR*ZBR;
	
	
	for (int IE = 1; IE <= NBE; IE++){
		fgets(LINHA, sizeof(LINHA), IRD);
		
		extrairString(APOIO, LINHA, 0, 9);	
		CEBR01_.EBT[IE-1] = atof(APOIO);
		
		if (*INFO >= 2) 
			fprintf(IWR, " %.2E", CEBR01_.EBT[IE-1]);

		int IWA = -1;
	
		for (int I = 0; I <= 5; I++){
		
			extrairString(APOIO, LINHA, 9, 12);
			IWA++;
			CEBR01_.XS[IWA][IE-1] = atof(APOIO);
		
			extrairString(APOIO, LINHA, 21, 12);
			IWA++;
			CEBR01_.XS[IWA][IE-1] = atof(APOIO);
		
			extrairString(APOIO, LINHA, 33, 12);
			IWA++;
			CEBR01_.XS[IWA][IE-1] = atof(APOIO);
		
			extrairString(APOIO, LINHA, 45, 12);
			IWA++;
			CEBR01_.XS[IWA][IE-1] = atof(APOIO);
		
			extrairString(APOIO, LINHA, 57, 12);
			IWA++;
			CEBR01_.XS[IWA][IE-1] = atof(APOIO);
			
			fgets(LINHA, sizeof(LINHA), IRD);
		}

		extrairString(APOIO, LINHA, 9, 12);
		IWA++;
		CEBR01_.XS[IWA][IE-1] = atof(APOIO);
		
		extrairString(APOIO, LINHA, 21, 12);
		IWA++;
		CEBR01_.XS[IWA][IE-1] = atof(APOIO);
		
		extrairString(APOIO, LINHA, 69, 10);
		CEBR01_.TXS[IE-1] = atof(APOIO);	
		
		if (*INFO >= 2) {
			for (int IW = 1; IW <= NBW; IW++){
				fprintf(IWR, " %.5E", CEBR01_.XS[IW-1][IE-1]);
				if ((IW == 5) || (IW == 10) || (IW == 15) || (IW == 20) || (IW == 25) || (IW == 30))
					fprintf(IWR, "\n");
			}
			fprintf(IWR, "                                      %.3E \n", CEBR01_.TXS[IE-1]);
		}
		
		CEBR01_.X[IE-1] = log(CEBR01_.EBT[IE-1]);
	}
	
	
	//inicializacao das rotinas de calculo
	
		
	for (int I = 1; I <= NBW; I++){
		CEBR_.WB[I-1] = WB0[I-1];	
	}
	
	//Calcule a distribui��o de perda de energia em escala e amostragem Par�metros para as energias na grade de simula��o.
	
	//interpola��o em E
	
		
	for (int IW = 1; IW <= NBW; IW++){
		for (int IE = 1; IE<= NBE; IE++){
			CEBR01_.Y[IE-1] = log(CEBR01_.XS[IW-1][IE-1]); 
		}
		
		spline2_(CEBR01_.X, CEBR01_.Y, A, B, C, D, 0.0, 0.0, NBE);
		
		for (int I = 1; I <= NEGP; I++){
			
			ELL = CEGRID_.DLEMP[I-1];
			int NBEV = NBE;
			
			findi2_(CEBR01_.X, &ELL, &NBEV, &J);
			
			if (ELL > CEBR01_.X[1-1]){
			
				CEBR02_.P0[IW-1][I-1][*M-1] = exp(A[J-1] + ELL * (B[J-1] + ELL * (C[J-1] + ELL * D[J-1])));
			
			} else{
				double F1 = A[1-1] + CEBR01_.X[1-1] * (B[1-1] + CEBR01_.X[1-1] * (C[1-1] + CEBR01_.X[1-1] * D[1-1]));
				double FP1 = B[1-1] + CEBR01_.X[1-1] * (2.0 * C[1-1] + CEBR01_.X[1-1] * 3.0 * D[1-1]);
			
				CEBR02_.P0[IW-1][I-1][*M-1] = exp(F1 + FP1 * (ELL - CEBR01_.X[1-1]));
			
			}	
		}	
	}
	
	
	
	for (int IE = 1; IE <= NEGP; IE++ ){
		for (int IW = 1; IW <= NBW; IW++){
			PDF[IW-1] = CEBR02_.P0[IW-1][IE-1][*M-1];
			
		
		}
		
		int NBWV = NBW;
			
		rlpac2_(CEBR_.WB, PDF, PAC, &NBWV);
	
		
		
		for (int IW = 1; IW <= NBW; IW++){
			
			CEBR_.PDFB[IW-1][IE-1][*M-1] = PDF[IW-1];
			CEBR_.PACB[IW-1][IE-1][*M-1] = PAC[IW-1];
		}
			
			
		for (int IW = 1; IW <= NBW-1; IW++){	
			CEBR_.DPDFB[IW-1][IE-1][*M-1] = CEBR_.PDFB[IW+1-1][IE-1][*M-1] - CEBR_.PDFB[IW-1][IE-1][*M-1];
		}
			
		CEBR_.DPDFB[NBW-1][IE-1][*M-1] = 0.0e0;
		
		//A perda de energia em escala de corte � ligeiramente modificada para garantir que a rotina de amostragem EBR cobre o intervalo de perda de energia permitido.
		double XC;
		if (IE < NEGP)
			XC = *WCRM / CEGRID_.ET[IE+1-1];
		else
			XC = *WCRM / CEGRID_.ET[NEGP-1];
		
	//	int MOM = -1;
		
		
		if (XC < 1.0e0){
			CEBR_.PBCUT[IE-1][*M-1] = rlmom2_(CEBR_.WB, PDF, XC, NBW, -1);
			CEBR_.WBCUT[IE-1][*M-1] = XC;
		} else{
			//double WXC = 1.0;
			CEBR_.PBCUT[IE-1][*M-1] = rlmom2_(CEBR_.WB, PDF, 1.0, NBW, -1);
			CEBR_.WBCUT[IE-1][*M-1] = 1.0;
		}
	}
	
	free(A);
	free(B);
	free(C);
	free(D);
	free(PAC);
	free(PDF);
	
//	printf("\n\nEBRAR2\n\n");	


	
}


void spline2_(double *X, double *Y, double *A, double *B, double *C, double *D, double S1, double SN, int const N){
	
	/*
	interpola��o de spline c�bica de dados tabulados.

  Entrada:
     X (I) (I = 1: N) ... pontos de grade (os valores de X devem estar em aumento
                      pedido).
     S (I) (I = 1: N) ... ou valores de fun��o respondentes.
     S1, SN .......... derivadas de segundo em X (1) e X (N). O natural
                      spline corresponds to take S1 = SN = 0.
     N .............. n�mero de pontos da grade.
  Sa�da:
     A (I), B (I), (I), D (I) (I = 1: N) ... coeffiients spline.

  O polin�mio cubico de interpola��o no intervalo I-�simo, de X (I) a
  X (I + 1), �
               P (x) = A (I) + x * (B (I) + x * ((I) + x * D (I)))

  Refereno: M.J. Maron, 'Numerial Analysis: a Pratial Approah',
             MaMillan Publ. o., Nova York, 1982
             
             
             
	*/
	
	int K;
	double R;
	double SI;
    double H;
    double HI;
    double FNP;
    double SI1;
 	
	if (N < 4){
		printf("*** Error in SPLINE: interpolation cannot be performed with %4d points.\n", N);
		exit(0);
	}
	
	int N1 = N-1;
	int N2 = N-2;
	
	
	
	
	for (int I = 1; I <= N1; I++){

		A[I-1] = X[I+1-1] - X[I-1];
			
	
		if (A[I-1] < 1.0e-12 * fmax(fabs(X[I-1]), fabs(X[I+1-1]))){
		   	printf("*** Error in SPLINE: X values not in increasing order.");
		   	exit(0);
		}
		D[I-1] = (Y[I+1-1] - Y[I-1]) / A[I-1];
	}
	
	//Matriz de coeficiente sim�trica (aumentada).

	for (int I = 1; I <= N2; I++){
		B[I-1] = 2.0e0 * (A[I-1] + A[I+1-1]);
		K = N1 - I + 1;
		D[K-1] = 6.0e0*(D[K-1]-D[K-1-1]);
	}
	
	D[2-1] = D[2-1] - A[1-1] * S1;
	D[N1-1] = D[N1-1] -A[N1-1] * SN;
	
	//Solu��o de Gauss do sistema tridiagonal.
	
	
	for (int I = 2; I <= N2; I++){
		R = A[I-1] / B[I-1-1];
		B[I-1] = B[I-1] - R*A[I-1];
		D[I+1-1] = D[I+1-1] -R*D[I-1];	
	}
	
	//Os coeficientes SIGMA s�o armazenados na matriz D.	
	
	D[N1-1] = D[N1-1] / B[N2-1];
	
		
	for (int I = 2; I <= N2; I++){
		K = N1 - I + 1;
		D[K-1] = (D[K-1] - A[K-1] * D[K+1-1]) / B[K-1-1];		
	}
	D[N-1] = SN;
	
	//Coeficientes Spline
	
	SI1 = S1;

	for (int I = 1; I <= N1; I++){
		SI = SI1;
		SI1 = D[I+1-1];
		H = A[I-1];
		HI = 1.0e0/H;
		A[I-1] = (HI / 6.0e0) * (SI * pow(X[I+1-1],3) - SI1 * pow(X[I-1],3))
				 +HI * (Y[I-1] * X[I+1-1] - Y[I+1-1] * X[I-1])
				 +(H / 6.0e0) * (SI1 * X[I-1] - SI * X[I+1-1]);	 
		B[I-1] = (HI / 2.0e0) * (SI1 * pow(X[I-1],2) - SI * pow(X[I+1-1],2))
				 +HI * (Y[I+1-1] - Y[I-1]) + (H / 6.0e0) * (SI - SI1);
		C[I-1] = (HI / 2.0e0) * (SI * X[I+1-1] - SI1 * X[I-1]);
		D[I-1] = (HI / 6.0e0) * (SI1 - SI);					 
	}
//Spline c�bico natural para X.GT.X (N).

	FNP = B[N1-1] + X[N-1] * (2.0e0 * C[N1-1] + X[N-1] * 3.0e0 * D[N1-1]);
	A[N-1] = Y[N-1] - X[N-1] * FNP;
	B[N-1] = FNP;
	C[N-1] = 0.0e0;
	D[N-1] = 0.0e0;
	
//	printf("\n\nSPLINE2\n\n");	
	
	return;

}


void rlpac2_(double *X, double *PDF, double *PAC, int *NP){
	
	/*
	Fun��o de distribui��o cumulativa de PDF (X) / X, obtida a partir de
  interpola��o em uma tabela de valores PDF.
  A vari�vel independente X assume apenas valores positivos.

    X ..... array de valores da vari�vel (em ordem crescente).
    PDF ... valores PDF correspondentes.
    PAC ... fun��o de probabilidade cumulativa.
    NP .... n�mero de pontos na tabela.
	*/
	
	double EPS = 1.0e-35;
	double X1, X2, Y1, Y2, DX, DY, B, A, DS;
	
	PAC[1-1] = 0.0;
	for (int I = 2; I <= *NP; I++){
		X1 = fmax(X[I-1-1], EPS);
		Y1 = PDF[I-1-1];
		X2 = fmax(X[I-1], EPS);
		Y2 = PDF[I-1];
		DX = X2 - X1;
		DY = Y2 - Y1;
		B = DY / DX;
		A = Y1 - B * X1;
		DS = A * log(X2/X1) + B*(X2-X1);
		PAC[I-1] = PAC[I-1-1] + DS;
		
	}
}


double rlmom2_(double *X, double *FCT, double XC, int NP, int MOM){
	
	/*
	C�lculo da integral de (X ** MOM) * FCT (X) no intervalo de
  X (1) a XC, obtido por interpola��o linear sobre uma tabela de FCT.
  A vari�vel independente X assume apenas valores positivos.

    X ....... array de valores da vari�vel (em ordem crescente).
    FCT ..... valores FCT correspondentes.
    NP ...... n�mero de pontos na tabela.
    XC ...... limite superior da integral, X (1) .LE.XC.LE.X (NP).
    MOM ..... ordem de momento (GE.-1).
    RLMOM ... integral de (X ** MOM) * FCT (X) no intervalo de X (1)
              para XC.
	*/
	double EPS = 1.0e-35;
	
	double X1, X2, Y1, Y2, XTC, DX, DY, A, B, DS;
	
	if (MOM < -1){
		printf("RLMOM. Error code 0\n");
		exit(0);
	}
	
	if (NP < 2){
		printf("RLMOM. Error code 1.\n");
		exit(0);	
	}
	
	if (X[1-1] < 0.0e0) {
		printf("RLMOM. Error code 2\n");
		exit(0);
	}
	for (int I = 2; I <= NP; I++){
		if (X[I-1] < 0.0e0){
			printf("RLMOM. Error code 3\n");
			exit(0);	
		}
		if (X[I-1] < X[I-1-1]){
			printf("RLMOM. Error code 4\n");
			exit(0);	
		}	
	}
	
	double resultado = 0.0e0;
	
	if (XC < X[1-1])
		return resultado;
	

	
	int IEND = 0;
	double XT = fmin(XC, X[NP-1]);
	for (int I = 1; I <= NP - 1; I++){
	
		X1 = fmax(X[I-1], EPS);
		Y1 = FCT[I-1];
		X2 = fmax(X[I+1-1], EPS);
		Y2 = FCT[I+1-1];
		if (XT < X2) {
			XTC = XT;
			IEND = 1;
		} else{
			XTC = X2;
		}
		DX = X2 - X1;
		DY = Y2 - Y1;
		
			
		
		if (fabs(DX) > 1.0e-14 * fabs(DY)){
			B = DY / DX;
			A = Y1 - B * X1;
			if (MOM == -1){
			   	DS = A * log(XTC / X1) + B * (XTC-X1);
			  
			} else{
			   	DS = A * (pow(XTC, (MOM+1)) - pow(X1, (MOM+1))) / (MOM+1) +
			   		 B * (pow(XTC, (MOM+2)) - pow(X1, (MOM+2))) / (MOM+2);
			   		 
			}
		} else{
			DS = 0.5e0 * (Y1 + Y2) * (XTC - X1) * pow(XTC, MOM); 
		
		
		}
	
		resultado = resultado + DS;

		if (IEND != 0){
			return resultado;
		
		}
	}

	return resultado;
	
	
}


void braar2_(int *M, FILE *IRD, FILE *IWR, int *INFO){
	
	/*
	Esta sub-rotina l� os par�metros de distribui��o angular bremsstrahlung
	do material M do arquivo de dados do material. Ele tamb�m inicializa o
	Algoritmo para gera��o da dire��o inicial de bremsstrahlung	F�tons C.
	*/
	
	
	double E[6], RK[4], Q1R[4][6], Q2R[4][6], Q1[21][6], Q2[21][6], X[6], A[6], B[6], C[6], D[6];
	double REV = 5.10998928e5; //Energia de repouso de el�trons (eV)
	double TREV = 2.0 * REV;
	
	E[0] = 1.0e3;
	E[1] = 5.0e3;
	E[2] = 1.0e4;
	E[3] = 5.0e4;
	E[4] = 1.0e5;
	E[5] = 5.0e5;
	
	RK[0] = 0.0e0;
	RK[1] = 0.6e0;
	RK[2] = 0.8e0;
	RK[3] = 0.95e0;
	

	
	for (int IE = 1; IE <= 6; IE++){
		CBRANG_.BET[IE-1] = sqrt(E[IE-1] * (E[IE-1] + TREV)) / (E[IE-1] + REV);
	}
	
//Grade de energias reduzidas de f�tons.

	for (int IK = 1; IK <= 21; IK++){
		CBRANG_.BK[IK-1] = (IK-1) * 0.05e0;
	}
	
//Leia os par�metros de distribui��o angular do arquivo (unidade IRD).	


	double ZEQ;
	int NDATA;
	
	fgets(LINHA, sizeof(LINHA), IRD);
	
//	printf("\nlinha %s\n", LINHA);
	extrairString(APOIO, LINHA, 40, 12);
	ZEQ = atof(APOIO);
	extrairString(APOIO, LINHA, 62, 4);
	NDATA = atoi(APOIO);
	if (*INFO >= 2){
		fprintf(IWR, "\n*** Bremss angular distribution,  ZEQ = %.5E ,  NDATA =%4d \n",ZEQ, NDATA );
	}
	
	if (NDATA != 24){
		printf("BRAR. Inconsistent data.\n");
		exit(0);
	}
	
	CBRANG_.ZBEQ[*M-1] = fmin(fmax(ZEQ, 2.0e0), 92.0e0);
	
	for (int IE1 = 1; IE1 <= 6; IE1++){
		for (int IK1 = 1; IK1 <= 4; IK1++){
			fgets(LINHA, sizeof(LINHA), IRD);
	    	extrairString(APOIO, LINHA, 0, 2);
	    	int IEE = atoi(APOIO);
	  	   	extrairString(APOIO, LINHA, 2, 2);
	 	    int IKK = atoi(APOIO);
    		extrairString(APOIO, LINHA, 5, 10);
	    	double ER = atof(APOIO);
  			extrairString(APOIO, LINHA, 16, 10);
	    	double RKR = atof(APOIO);
  			extrairString(APOIO, LINHA, 27, 14);
	    	double Q1RR = atof(APOIO);
	    	extrairString(APOIO, LINHA, 42, 14);
	    	double Q2RR = atof(APOIO);
	    	
	    	if ((fabs(ER - E[IEE-1]) < 1.0e-6 ) && (fabs(RKR-RK[IKK-1]) < 1.0e-6)) {
		        Q1R[IKK-1][IEE-1] = Q1RR/ZEQ;
		        Q2R[IKK-1][IEE-1] = Q2RR;
					
			} else{
				printf("Corrupt data file (pdbrang.p08)\n");
				exit(0);
			}	
		}	
	}
	
	
	if (*INFO >= 2){
		for (int IE = 1; IE <= 6; IE++){
			for (int IK = 1; IK <= 4; IK++){
				fprintf(IWR, " %.3E  %.3E  %.8E  %.8E\n", E[IE-1], RK[IK-1], Q1R[IK-1][IE-1] * ZEQ, Q2R[IK-1][IE-1]);	
			}
		}
	}
	
	//Tabela expandida de par�metros de distribui��o.
	
	int J;
	int WW = 4;
	
	for (int IE = 1; IE <= 6; IE++){
		for (int IK = 1; IK <= 4; IK++){
			X[IK-1] = log(Q1R[IK-1][IE-1]);
		}
		spline2_(RK, X, A, B, C, D, 0.0, 0.0, 4);
		for (int IK = 1; IK <= 21; IK++){
			
			findi2_(RK, &CBRANG_.BK[IK-1], &WW, &J);
			Q1[IK-1][IE-1] = A[J-1] + CBRANG_.BK[IK-1]*(B[J-1] + CBRANG_.BK[IK-1] * (C[J-1] + CBRANG_.BK[IK-1] * D[J-1]));
		}
		
		for (int IK = 1; IK <= 4; IK++){
			X[IK-1] = Q2R[IK-1][IE-1];
		}
		spline2_(RK, X, A, B, C, D, 0.0, 0.0, 4);
		
		for (int IK = 1; IK <= 21; IK ++){
			findi2_(RK, &CBRANG_.BK[IK-1], &WW, &J);
			Q2[IK-1][IE-1] = A[J-1] + CBRANG_.BK[IK-1]*(B[J-1] + CBRANG_.BK[IK-1] * (C[J-1] + CBRANG_.BK[IK-1] * D[J-1]));
		}
	}
	
	
	//e interpola��es spline c�bicas naturais.

	
	for (int IK = 1; IK <= 21; IK++){
		for (int IE = 1; IE <= 6; IE++){
			X[IE-1] = Q1[IK-1][IE-1];
		}
		
		spline2_(CBRANG_.BET, X, A, B, C, D, 0.0, 0.0, 6);
		for (int IE = 1; IE <= 6; IE++){
			CBRANG_.BP1[1-1][IK-1][IE-1][*M-1] = A[IE-1];
			CBRANG_.BP1[2-1][IK-1][IE-1][*M-1] = B[IE-1];
			CBRANG_.BP1[3-1][IK-1][IE-1][*M-1] = C[IE-1];
			CBRANG_.BP1[4-1][IK-1][IE-1][*M-1] = D[IE-1];	
		}
		
		for (int IE = 1; IE <= 6; IE++){
			X[IE-1] = Q2[IK-1][IE-1];
			
		}
		
		spline2_(CBRANG_.BET, X, A, B, C, D, 0.0, 0.0, 6);
		
		for (int IE = 1; IE <= 6; IE++){
			CBRANG_.BP2[1-1][IK-1][IE-1][*M-1] = A[IE-1];
			CBRANG_.BP2[2-1][IK-1][IE-1][*M-1] = B[IE-1];
			CBRANG_.BP2[3-1][IK-1][IE-1][*M-1] = C[IE-1];
			CBRANG_.BP2[4-1][IK-1][IE-1][*M-1] = D[IE-1];	
		}	
	}
	
	//printf("\n\nBRAAR2\n\n");	

}


void einat2_(double &E, double &WCCM, double &XH0, double &XH1, double &XH2, double &XS0, double &XS1, double &XS2, double &XT1, double &XT2, double &DELTA, int *M){ //REVISTO
	
	/*
	Se��es transversais integradas para colis�es inel�sticas de el�trons de
  energia E no material M, restrita a perdas de energia maiores que e
  menos do que a energia de corte WCCM.

  Modelo Sternheimer-Liljequist GOS.

  Argumentos de sa�da:
    XH0 ... se��o transversal total para colls r�gidos. (cm ** 2).
    XH1 ... parando a se��o transversal para colls r�gidos. (eV * cm ** 2).
    XH2 ... se��o transversal extensa para colls r�gidos. (eV ** 2 * cm ** 2).
    XS0 ... se��o transversal total para colls macios. (cm ** 2).
    XS1 ... parando a se��o transversal para colls macios. (eV * cm ** 2)
    XS2 ... se��o transversal extensa para colls macios. (eV ** 2 * cm ** 2).
    XT1 ... 1� se��o transversal de transporte para soft colls. (cm ** 2).
    XT2 ... 2� se��o transversal de transporte para soft colls. (cm ** 2).
    DELTA ... Corre��o do efeito de densidade de Fermi.
	
	*/
	
	//constantes]
	
	static const int NO = 512;
	
	double H0, H1, H2, S0, S1, S2, R0, R1, R2, XT0;
	double WL2L;
	double WL2U;
	
	double REV = 5.10998928e5; //  ! Electron rest energy (eV)
	
	double GAM = 1.0e0 + E / REV;
	double GAM2 = GAM * GAM;
	double UK, WK;
	
	// Densidade efetiva
	
	//****  Sternheimer's resonance energy (WL2=L**2).
	
	double TST = COMPOS_.ZT[*M-1] / (GAM2 * CEIN_.OP2[*M-1]);

	
	double WL2 = 0.0e0;
	double FDEL = 0.0e0;

	for (int I = 1; I <= CEIN_.NOSC[*M-1]; I++){
		FDEL = FDEL + CEIN_.F[I-1][*M-1] / (pow(CEIN_.WRI[I-1][*M-1], 2) + WL2);
	}

	
	if (FDEL < TST){
		DELTA = 0.0e0;
		goto L3;
	}
	WL2 = pow(CEIN_.WRI[CEIN_.NOSC[*M-1]-1][*M-1], 2);
	
L1:;
	
	WL2 = WL2 + WL2;
	FDEL = 0.0e0;
	for (int I = 1; I <= CEIN_.NOSC[*M-1]; I++){
		FDEL = FDEL + CEIN_.F[I-1][*M-1] / (pow(CEIN_.WRI[I-1][*M-1], 2) + WL2);
	}
	
	if (FDEL > TST)
		goto L1;
	WL2L = 0.0e0;
	WL2U = WL2;
L2:;
	WL2 = 0.5e0 * (WL2L + WL2U);
	FDEL = 0.0e0;
	for (int I = 1; I <= CEIN_.NOSC[*M-1]; I++){
		FDEL = FDEL + CEIN_.F[I-1][*M-1] / (pow(CEIN_.WRI[I-1][*M-1], 2) + WL2);
	}
	
	if (FDEL > TST){
		WL2L = WL2;
	}else{
		WL2U = WL2;
	}
	
	if (WL2U - WL2L > (1.0e-12 * WL2))
		goto L2;
	
	// Density effect correction (delta).
	
	DELTA = 0.0e0;
	
	for (int I = 1; I <= CEIN_.NOSC[*M-1]; I++){
		
		DELTA = DELTA + CEIN_.F[I-1][*M-1] * log(1.0e0 + WL2 / pow(CEIN_.WRI[I-1][*M-1], 2));
	}
	DELTA = (DELTA / COMPOS_.ZT[*M-1]) - WL2 / (GAM2*CEIN_.OP2[*M-1]);

L3:;


//Se��es transversais do oscilador de concha.

	for (int I = 1; I <= CEIN_.NOSC[*M-1]; I++){
		CEIN00_.SEH0[I-1] = 0.0e0;
		CEIN00_.SEH1[I-1] = 0.0e0;
		CEIN00_.SEH2[I-1] = 0.0e0;
		CEIN00_.SES0[I-1] = 0.0e0;
		CEIN00_.SES1[I-1] = 0.0e0;
		CEIN00_.SES2[I-1] = 0.0e0;
		CEIN00_.SET0[I-1] = 0.0e0;
		CEIN00_.SET1[I-1] = 0.0e0;
		CEIN00_.SET2[I-1] = 0.0e0;
	}
	
	XH0=0.0e0;
	XH1=0.0e0;
	XH2=0.0e0;
	XS0=0.0e0;
	XS1=0.0e0;
	XS2=0.0e0;
	XT0=0.0e0;
	XT1=0.0e0;
	XT2=0.0e0;

	
	for (int K = 1; K <= CEIN_.NOSC[*M-1]; K++){
		UK = CEIN_.UI[K-1][*M-1];
		WK = CEIN_.WRI[K-1][*M-1];
		
		einat12_(E, UK, WK, DELTA, WCCM, H0, H1, H2, S0, S1, S2, R0, R1, R2);

		CEIN00_.SEH0[K-1] = CEIN_.F[K-1][*M-1] * H0;
		CEIN00_.SEH1[K-1] = CEIN_.F[K-1][*M-1] * H1;
		CEIN00_.SEH2[K-1] = CEIN_.F[K-1][*M-1] * H2;
		CEIN00_.SES0[K-1] = CEIN_.F[K-1][*M-1] * S0;
		CEIN00_.SES1[K-1] = CEIN_.F[K-1][*M-1] * S1;
		CEIN00_.SES2[K-1] = CEIN_.F[K-1][*M-1] * S2;
		CEIN00_.SET0[K-1] = CEIN_.F[K-1][*M-1] * R0;
		CEIN00_.SET1[K-1] = CEIN_.F[K-1][*M-1] * 2.0e0*R1;
		CEIN00_.SET2[K-1] = CEIN_.F[K-1][*M-1] * 6.0e0*(R1-R2);
		
		XH0= XH0 + CEIN00_.SEH0[K-1];
		XH1= XH1 + CEIN00_.SEH1[K-1];
		XH2= XH2 + CEIN00_.SEH2[K-1];
		XS0= XS0 + CEIN00_.SES0[K-1];
		XS1= XS1 + CEIN00_.SES1[K-1];
		XS2= XS2 + CEIN00_.SES2[K-1];
		XT0= XT0 + CEIN00_.SET0[K-1];
		XT1= XT1 + CEIN00_.SET1[K-1];
		XT2= XT2 + CEIN00_.SET2[K-1];	
	}
	
//		printf("\n\nEINAT2\n\n");	
	
		
	
}


void einat12_(double &E, double &UK, double &WK, double &DELTA, double &WCCM, double &H0, double &H1, double &H2, double &S0, double &S1, double &S2, double &R0, double &R1, double &R2){//REVISTO
	/*
	Se��es transversais integradas para colis�es inel�sticas de el�trons com
  um oscilador de camada �nica, restrito a perdas de energia maiores que,
  e menor do que, a perda de energia de corte WCCM.

  Modelo do oscilador Sternheimer-Liljequist.

  Argumentos de entrada:
    E ..... energia cin�tica (eV).
    UK .... energia de ioniza��o (eV).
    WK .... energia de resson�ncia (eV).
    DELTA ... Corre��o do efeito de densidade de Fermi.
    WCCM ... perda de energia de corte (eV).

  Argumentos de sa�da:
    H0 .... se��o transversal total para colls duros. (cm ** 2).
    H1 .... parando a se��o transversal para colls duros. (eV * cm ** 2).
    H2 .... se��o transversal de straggling para hard colls. (eV ** 2 * cm ** 2).
    S0 .... se��o transversal total para colls macios. (cm ** 2).
    S1 .... parando a se��o transversal para colls macios. (eV * cm ** 2).
    S2 .... se��o transversal extensa para colls macios. (eV ** 2 * cm ** 2).
    R0 .... se��o transversal total para colls macios. (cm ** 2).
    R1 .... 1� se��o transversal de transporte para soft colls. (cm ** 2).
    R2 .... 2� se��o transversal de transporte para soft colls. (cm ** 2).
	*/
	
	double REV = 5.10998928e5;
	double ELRAD= 2.8179403267e-13;
	double TREV= 2.0e0 * REV;
	double RTREV= 1.0e0 / TREV;
	double PI = 3.1415926535897932e0;
	double PIELR2= PI * ELRAD * ELRAD;
	
	double WTHR;
	double WM, WKP, QKP, WCMAX, WDMAX;
	double QM, SDL1, SDT1, CPPS, CPP, A, B, RMU1, BA;
	double F0, F1, F2, SD1, WL, WU;
	double CP2S, CP2, RMU2, CP3S, CP3, RMU3; 
	
	
	H0 = 0.0e0;
	H1 = 0.0e0;
	H2 = 0.0e0;
	S0 = 0.0e0;
	S1 = 0.0e0;
	S2 = 0.0e0;
	R0 = 0.0e0;
	R1 = 0.0e0;
	R2 = 0.0e0;
	
	if (UK > 1.0e-3){
		WTHR = UK;
	} else{
		WTHR = WK;
	}
	
	if (E < WTHR + 1.0e-6)
		return;
	
	
	//constantes
	
	*CEIN01_.EI = E;
	double GAM = 1.0e0 + E / REV;
	double GAM2 = GAM*GAM;
	double BETA2 = (GAM2 - 1.0e0) / GAM2;
	double CONST = PIELR2 * TREV / BETA2;
	
	*CEIN01_.CPS = E * (E+TREV);
	double CP = sqrt(*CEIN01_.CPS);
	*CEIN01_.AMOL = pow(E/(E+REV), 2);
	 
	 
	//Truque: A energia de resson�ncia e a energia de recuo de corte de
	//As camadas internas  s�o variadas para produzir um limite suave. 
	
	
	if (UK > 1.0e-3){
		WM = 3.0e0*WK - 2.0e0*UK;
		if (E > WM){
			WKP = WK;
			QKP = UK;
		} else{
			WKP = (E + 2.0e0 * UK) / 3.0e0;
			QKP = UK * (E/WM);
			WM = E;
			
		}
		*CEIN01_.EE = E + UK;
		WCMAX = 0.5e0 * *CEIN01_.EE;
		WDMAX = fmin(WM, WCMAX);
	}else{
		WM = E;
		WKP = WK;
		QKP = WK;
		*CEIN01_.EE = E;
		WCMAX = 0.5e0 * *CEIN01_.EE;
		WDMAX = WKP + 1.0e0;	
	}
	
	//Interacoes distantes
	
	SDL1 = 0.0e0;
	SDT1 = 0.0e0;
	if (WDMAX > (WTHR + 1.0e-6)){
		CPPS = (E - WKP) * (E - WKP+TREV);
		CPP = sqrt(CPPS);
		A = 4.0e0 * CP * CPP;
		B = pow(CP - CPP, 2);
		
		if (WKP > (1.0e-6 * E)){
			QM = sqrt(pow(CP - CPP, 2) + pow(REV, 2)) - REV;
		} else{
			QM = WKP * WKP / (BETA2 * TREV);
			QM = QM * (1.0e0 - QM*RTREV);
		}
		
		if (QM < QKP){
			SDL1 = log(QKP * (QM + TREV ) / (QM * (QKP + TREV )));
			SDT1=fmax(log(GAM2) - BETA2 - DELTA, 0.0e0);
			//Momentos de transporte distantes suaves das ordens 0-2.
			if (WCCM > WTHR){
				BA=B/A;
            	RMU1=(QKP*(QKP+TREV)-B)/A;
            	R0=log((RMU1+BA)/BA);
            	R1=RMU1-BA*R0;
            	R2=pow(BA,2)*R0+0.5e0*RMU1*(RMU1-2.0e0*BA);
            	R0=R0/WKP;
            	R1=R1/WKP;
            	R2=R2/WKP;
            	R0=R0+SDT1/WKP;
			}
		}
	}
	
	SD1=SDL1+SDT1;
	
	if (SD1 > 0.0e0){
		if (UK > 1.0e-3){
		//Excita��es internas (distribui��o triangular).	
			F0=1.0e0/pow(WM-UK, 2);
			F1=2.0e0*F0*SD1/WKP;
			
			if (WCCM < UK){
				WL=UK;
            	WU=WDMAX;
            	H0=F1*(WM*(WU-WL)-(pow(WU,2)- pow(WL,2))/2.0e0);
            	H1=F1*(WM*(pow(WU,2)- pow(WL,2))/2.0e0-(pow(WU,3)- pow(WL,3))/3.0e0);
            	H2=F1*(WM*(pow(WU,3)- pow(WL,3))/3.0e0-(pow(WU,4)- pow(WL,4))/4.0e0);
			} else{
				
				if (WCCM > WDMAX){
					WL=UK;
            		WU=WDMAX;
            		S0=F1*(WM*(WU-WL)-(pow(WU,2)- pow(WL,2))/2.0e0);
            		S1=F1*(WM*(pow(WU,2)- pow(WL,2))/2.0e0-(pow(WU,3)- pow(WL,3))/3.0e0);
            		S2=F1*(WM*(pow(WU,3)- pow(WL,3))/3.0e0-(pow(WU,4)- pow(WL,4))/4.0e0);
					
				} else{
					WL=WCCM;
            		WU=WDMAX;
            		H0=F1*(WM*(WU-WL)-(pow(WU,2)- pow(WL,2))/2.0e0);
            		H1=F1*(WM*(pow(WU,2)- pow(WL,2))/2.0e0-(pow(WU,3)- pow(WL,3))/3.0e0);
            		H2=F1*(WM*(pow(WU,3)- pow(WL,3))/3.0e0-(pow(WU,4)- pow(WL,4))/4.0e0);
            		WL=UK;
            		WU=WCCM;
            		S0=F1*(WM*(WU-WL)-(pow(WU,2)- pow(WL,2))/2.0e0);
            		S1=F1*(WM*(pow(WU,2)- pow(WL,2))/2.0e0-(pow(WU,3)- pow(WL,3))/3.0e0);
            		S2=F1*(WM*(pow(WU,3)- pow(WL,3))/3.0e0-(pow(WU,4)- pow(WL,4))/4.0e0);	
				}
				F2=F0*(2.0e0*WM*(WU-WL)-(pow(WU,2)- pow(WL,2)));
            	R0=F2*R0;
            	R1=F2*R1;
            	R2=F2*R2;
			}	
		} else{
			
			//Excita��es externas (oscilador delta).
			if (WCCM < WKP){
				H1=SD1;
            	H0=SD1/WKP;
            	H2=SD1*WKP;
				
			} else{
				S1=SD1;
            	S0=SD1/WKP;
  	 	 	 	S2=SD1*WKP;
			}	
		}
	}
	
	//Fechar colis�es (se��o transversal de Moller).
	
	if (WCMAX < (WTHR + 1.0e-6))
		goto L1;
	

	if (WCCM < WTHR){ //No soft interactions
		WL=WTHR;
        WU=WCMAX;
        H0=H0+(1.0e0/(*CEIN01_.EE-WU))-(1.0e0/(*CEIN01_.EE-WL))
         -(1.0e0/WU)+(1.0e0/WL)
         +(1.0e0-*CEIN01_.AMOL)*log(((*CEIN01_.EE-WU)*WL)/((*CEIN01_.EE-WL)*WU))/ *CEIN01_.EE
         +*CEIN01_.AMOL*(WU-WL)/pow(*CEIN01_.EE,2);
         
         H1=H1+log(WU/WL)+(*CEIN01_.EE/(*CEIN01_.EE-WU))-(*CEIN01_.EE/(*CEIN01_.EE-WL))
         +(2.0e0-*CEIN01_.AMOL)*log((*CEIN01_.EE-WU)/(*CEIN01_.EE-WL))
         +*CEIN01_.AMOL*(pow(WU,2)- pow(WL,2))/(2.0e0*pow(*CEIN01_.EE,2));
         
         H2=H2+(2.0e0-*CEIN01_.AMOL)*(WU-WL)+(WU*(2.0e0* *CEIN01_.EE-WU)/(*CEIN01_.EE-WU))
         -(WL*(2.0e0* *CEIN01_.EE-WL)/(*CEIN01_.EE-WL))
         +(3.0e0-*CEIN01_.AMOL)* *CEIN01_.EE*log((*CEIN01_.EE-WU)/(*CEIN01_.EE-WL))
		 +*CEIN01_.AMOL*(pow(WU,3)- pow(WL,3))/(3.0e0*pow(*CEIN01_.EE,2));
	} else{
		if (WCCM > WCMAX){
			WL=WTHR;
        	WU=WCMAX;
	        S0=S0+(1.0e0/(*CEIN01_.EE-WU))-(1.0e0/(*CEIN01_.EE-WL))
	         -(1.0e0/WU)+(1.0e0/WL)
	         +(1.0e0-*CEIN01_.AMOL)*log(((*CEIN01_.EE-WU)*WL)/((*CEIN01_.EE-WL)*WU))/ *CEIN01_.EE
	         +*CEIN01_.AMOL*(WU-WL)/pow(*CEIN01_.EE,2);
	         
	         S1=S1+log(WU/WL)+(*CEIN01_.EE/(*CEIN01_.EE-WU))-(*CEIN01_.EE/(*CEIN01_.EE-WL))
	         +(2.0e0-*CEIN01_.AMOL)*log((*CEIN01_.EE-WU)/(*CEIN01_.EE-WL))
	         +*CEIN01_.AMOL*(pow(WU,2)- pow(WL,2))/(2.0e0*pow(*CEIN01_.EE,2));
	         
	         S2=S2+(2.0e0-*CEIN01_.AMOL)*(WU-WL)+(WU*(2.0e0* *CEIN01_.EE-WU)/(*CEIN01_.EE-WU))
	         -(WL*(2.0e0* *CEIN01_.EE-WL)/(*CEIN01_.EE-WL))
	         +(3.0e0-*CEIN01_.AMOL)* *CEIN01_.EE*log((*CEIN01_.EE-WU)/(*CEIN01_.EE-WL))
			 +*CEIN01_.AMOL*(pow(WU,3)- pow(WL,3))/(3.0e0*pow(*CEIN01_.EE,2));
		} else{
			WL=WCCM;
	        WU=WCMAX;
    		H0=H0+(1.0e0/(*CEIN01_.EE-WU))-(1.0e0/(*CEIN01_.EE-WL))
	         -(1.0e0/WU)+(1.0e0/WL)
	         +(1.0e0-*CEIN01_.AMOL)*log(((*CEIN01_.EE-WU)*WL)/((*CEIN01_.EE-WL)*WU))/ *CEIN01_.EE
	         +*CEIN01_.AMOL*(WU-WL)/pow(*CEIN01_.EE,2);
	         
	         H1=H1+log(WU/WL)+(*CEIN01_.EE/(*CEIN01_.EE-WU))-(*CEIN01_.EE/(*CEIN01_.EE-WL))
	         +(2.0e0-*CEIN01_.AMOL)*log((*CEIN01_.EE-WU)/(*CEIN01_.EE-WL))
	         +*CEIN01_.AMOL*(pow(WU,2)- pow(WL,2))/(2.0e0*pow(*CEIN01_.EE,2));
	         
	         H2=H2+(2.0e0-*CEIN01_.AMOL)*(WU-WL)+(WU*(2.0e0**CEIN01_.EE-WU)/(*CEIN01_.EE-WU))
	         -(WL*(2.0e0* *CEIN01_.EE-WL)/(*CEIN01_.EE-WL))
	         +(3.0e0-*CEIN01_.AMOL)* *CEIN01_.EE*log((*CEIN01_.EE-WU)/(*CEIN01_.EE-WL))
			 +*CEIN01_.AMOL*(pow(WU,3)- pow(WL,3))/(3.0e0*pow(*CEIN01_.EE,2));
			 
		    WL=WTHR;
          	WU=WCCM;
	        S0=S0+(1.0e0/(*CEIN01_.EE-WU))-(1.0e0/(*CEIN01_.EE-WL))
	         -(1.0e0/WU)+(1.0e0/WL)
	         +(1.0e0-*CEIN01_.AMOL)*log(((*CEIN01_.EE-WU)*WL)/((*CEIN01_.EE-WL)*WU))/ *CEIN01_.EE
	         +*CEIN01_.AMOL*(WU-WL)/pow(*CEIN01_.EE,2);
	         
	         S1=S1+log(WU/WL)+(*CEIN01_.EE/(*CEIN01_.EE-WU))-(*CEIN01_.EE/(*CEIN01_.EE-WL))
	         +(2.0e0-*CEIN01_.AMOL)*log((*CEIN01_.EE-WU)/(*CEIN01_.EE-WL))
	         +*CEIN01_.AMOL*(pow(WU,2)- pow(WL,2))/(2.0e0*pow(*CEIN01_.EE,2));
	         
	         S2=S2+(2.0e0-*CEIN01_.AMOL)*(WU-WL)+(WU*(2.0e0* *CEIN01_.EE-WU)/(*CEIN01_.EE-WU))
	         -(WL*(2.0e0* *CEIN01_.EE-WL)/(*CEIN01_.EE-WL))
	         +(3.0e0-*CEIN01_.AMOL)* *CEIN01_.EE*log((*CEIN01_.EE-WU)/(*CEIN01_.EE-WL))
			 +*CEIN01_.AMOL*(pow(WU,3)- pow(WL,3))/(3.0e0*pow(*CEIN01_.EE,2));		
		}
		
		//Momentos de transporte de fechamento suave das ordens 0-2.
		
    	CP2S=(E-WL)*(E-WL+TREV);
        CP2=sqrt(CP2S);
        RMU2=(WL*(WL+TREV)-pow((CP-CP2),2))/(4.0e0*CP*CP2);
        CP3S=(E-WU)*(E-WU+TREV);
        CP3=sqrt(CP3S);
        RMU3=(WU*(WU+TREV)-pow((CP-CP3),2))/(4.0e0*CP*CP3);
        *CEIN01_.MOM = 0;
        int FUNCAO = 1;
        double TOL = 1.0e-7;
     //   printf("\n\nC++ RMU2: %.16e, RMU3: %.17f\n\n", RMU2, RMU3);
        
        R0= R0+ sumga2_(FUNCAO, RMU2,RMU3,TOL);
        *CEIN01_.MOM = 1;
        R1 = R1+ sumga2_(FUNCAO, RMU2,RMU3,TOL);
        *CEIN01_.MOM = 2;
        R2 = R2+ sumga2_(FUNCAO, RMU2,RMU3,TOL);
          
      //  exit(0);
	}
	
L1:;

//0.14877436420578918
//0.14877436420578921

    H0=CONST*H0;
	H1=CONST*H1;
	H2=CONST*H2;
	S0=CONST*S0;
	S1=CONST*S1;
	S2=CONST*S2;
	R0=CONST*R0;
	R1=CONST*R1;
	R2=CONST*R2;
	
	//	printf("\n\nEINAT12\n\n");	
	
}


double einads2_(double RMU){//OK
	
	/*
	Se��o transversal diferencial angular para soft close inel�stico colis�es de el�trons.
	*/
	
//	printf("\n\nRMU 2: %f\n\n", RMU);
	double resultado = 0.0e0;
	
	double REV = 5.10998928e5; 
	
	double AUX=2.0e0*RMU*(1.0e0-RMU);
	double DENOM= *CEIN01_.EI * AUX + REV;
	double W= *CEIN01_.CPS * AUX / DENOM;
	double DWDMU = *CEIN01_.CPS * REV * (2.0e0-4.0e0*RMU) / pow(DENOM,2);
	
	 resultado = (1.0e0 + pow( W / (*CEIN01_.EE - W),2) - (1.0e0 - *CEIN01_.AMOL)*(W/(*CEIN01_.EE - W))
           + *CEIN01_.AMOL * pow(W / *CEIN01_.EE, 2)) * DWDMU * pow(RMU,*CEIN01_.MOM) / pow(W,2);
	
//		   printf("\n\nresultado %f\n\n", resultado);
	return resultado;
	
//	printf("\n\nEINADS2\n\n");
}


double sumga2_(int &funcao, double &XL, double &XU, double &TOL){//OK
	
	/*
	Esta fun��o calcula o valor SUMGA da integral do
  (externo) fun��o FCT no intervalo (XL, XU) usando o ponto de 20
  M�todo da quadratura de Gauss com um esquema de bissec��o adaptativa.

  TOL � a toler�ncia, ou seja, erro relativo m�ximo permitido; deveria
  n�o pode ser inferior a 1.0D-13. Uma mensagem de aviso � escrita na unidade 6 quando
  a precis�o necess�ria n�o � atingida. O bloco comum CSUMGA pode ser
  usado para transferir o sinalizador de erro IERGA e o n�mero de
  valores de fun��o para o programa de chamada.
  
  C�digos de erro de sa�da:
     IERGA = 0, sem problemas, o c�lculo convergiu.
           = 1, muitos subintervalos abertos.
           = 2, muitas chamadas de fun��o.
           = 3, os subintervalos s�o muito estreitos.
  
  
	*/
	
//	printf("\n\nTOL: %.8f, Funcao: %d, XL: %.12f, XU: %.12f\n\n", TOL, funcao, XL, XU);
	
	static const int NP = 10;
	static const int NOIT = 128;
	
	int NP2 = 2 * NP;
	int NP4 = 4 * NP;
	int NOIT5 = NOIT / 5;
	int NCALLT = 100000;
	
	double X[NP], W[NP], XM[NP], XP[NP];
	double S[NOIT], SN[NOIT], XR[NOIT], XRN[NOIT];
	
	double TOL1, TOL2, TOL3, resultado, H, HH, X1, X2, SP, AHH, SUMR, SI, S1, S2, S12;
	int NOI, IDONE, NOIP;
	
	//F�rmula de integra��o de 20 pontos de Gauss. Abscissas.
	
	X[0] = 7.6526521133497334e-02;
	X[1] = 2.2778585114164508e-01;
	X[2] = 3.7370608871541956e-01;
	X[3] = 5.1086700195082710e-01;
	X[4] = 6.3605368072651503e-01;
	X[5] = 7.4633190646015079e-01;
	X[6] = 8.3911697182221882e-01;
	X[7] = 9.1223442825132591e-01;
	X[8] = 9.6397192727791379e-01;
	X[9] = 9.9312859918509492e-01;
	
	//Pesos
	W[0] = 1.5275338713072585e-01;
	W[1] = 1.4917298647260375e-01;
	W[2] = 1.4209610931838205e-01;
	W[3] = 1.3168863844917663e-01;
	W[4] = 1.1819453196151842e-01;
	W[5] = 1.0193011981724044e-01;
	W[6] = 8.3276741576704749e-02;
	W[7] = 6.2672048334109064e-02;
	W[8] = 4.0601429800386941e-02;
	W[9] = 1.7614007139152118e-02;
	
	for (int I = 1; I <= NP; I++){
		XM[I-1] = 1.0e0 - X[I-1];
		XP[I-1] = 1.0e0 + X[I-1];	
	}
	
	//Toler�ncias globais e parciais.
	
	TOL1=fmin(fmax(TOL,1.0e-13),1.0e-5);  // Global tolerance.
	TOL2=TOL1;  //! Effective tolerance.
	TOL3=1.0e-13;  //! Round-off protection.
	resultado = 0.0e0;
	*CSUMGA_.IERGA = 0;
	
	//Integra��o direta de XL a XU.	
	
	H = XU - XL;
    HH=0.5e0*H;
    X1 = XL;
    switch (funcao){
		case 1: // einads
			SP = W[1-1] * (einads2_(X1+XM[1-1] * HH)+ einads2_(X1+XP[1-1] * HH));
			break;
		case 2: // pinads	
			SP = W[1-1] * (pinads2_(X1+XM[1-1] * HH)+ pinads2_(X1+XP[1-1] * HH));
			break;
	}
	
	for (int J = 2; J <= NP; J++){
		switch (funcao){
			case 1: // einads
				SP = SP + W[J-1] * (einads2_( X1 + XM[J-1] * HH) + einads2_( X1 + XP[J-1] * HH));
				break;
		
			case 2: // einads
				SP = SP + W[J-1] * (pinads2_( X1 + XM[J-1] * HH) + pinads2_( X1 + XP[J-1] * HH));
				break;	
		
	    }	
	}
	
	S[1-1]= SP * HH;
	XR[1-1] = X1;
	*CSUMGA_.NCALL = NP2;
	NOI = 1;
	IDONE=1; //  ! To prevent a compilation warning.
	
	
	//Esquema de bissec��o adaptativa.
	
L1:;
	H=HH;  //! Subinterval length.
	HH=0.5e0*H;
	AHH=fabs(HH);
	if (TOL2 > 0.01e0 * TOL1)
		TOL2 = TOL2 * 0.5e0;
	
	SUMR=0.0e0;
	NOIP=NOI;
	NOI=0;
//	printf("\n\nNOIP: %d\n\n", NOIP);
	for (int I = 1; I <= NOIP; I++){
		SI=S[I-1];  //! Bisect the I-th open interval.
        X1=XR[I-1];
        if(AHH < fabs(X1) * TOL3)
			*CSUMGA_.IERGA = 3; // ! The interval is too narrow.
        switch (funcao){
			case 1: // einads
				SP = W[1-1] * (einads2_(X1+XM[1-1] * HH)+ einads2_(X1+XP[1-1] * HH));
				break;
			case 2: // pinads
				SP = W[1-1] * (pinads2_(X1+XM[1-1] * HH)+ pinads2_(X1+XP[1-1] * HH));
				break;
			}
	
		for (int J = 2; J <= NP; J++){
			switch (funcao){
				case 1: // einads
					SP = SP + W[J-1] * (einads2_( X1 + XM[J-1] * HH) + einads2_( X1 + XP[J-1] * HH));
					break;
				case 2: // pinads
					SP = SP + W[J-1] * (pinads2_( X1 + XM[J-1] * HH) + pinads2_( X1 + XP[J-1] * HH));
					break;
	  	    }	
	    }	
	    S1 = SP * HH;
	    X2 = X1 + H;
	    
	    if(AHH < fabs(X2) * TOL3)
			*CSUMGA_.IERGA = 3; // ! The interval is too narrow.
        switch (funcao){
			case 1: // einads
				SP = W[1-1] * (einads2_(X2 + XM[1-1] * HH)+ einads2_(X2 + XP[1-1] * HH));
				break;
			case 2: // pinads
				SP = W[1-1] * (pinads2_(X2 + XM[1-1] * HH)+ pinads2_(X2 + XP[1-1] * HH));
				break;
			}
	
		for (int J = 2; J <= NP; J++){
			switch (funcao){
				case 1: // einads
					SP = SP + W[J-1] * (einads2_(X2 + XM[J-1] * HH) + einads2_( X2 + XP[J-1] * HH));
					break;
				case 2: // pinads
					SP = SP + W[J-1] * (pinads2_(X2 + XM[J-1] * HH) + pinads2_( X2 + XP[J-1] * HH));
					break;
	  	    }	
	    }	
	    
	    S2= SP * HH;
        IDONE = I;
        *CSUMGA_.NCALL = *CSUMGA_.NCALL + NP4;
        S12 = S1 + S2; // ! Sum of integrals on the two subintervals.
        
        if( fabs(S12 - SI) < fmax( TOL2 * fabs(S12), 1.0e-35)) {
        	//A integral sobre o intervalo pai convergiu.
        	resultado = resultado + S12;
			
		} else{
			SUMR = SUMR + S12;
 	  	    NOI = NOI + 2;
            if(NOI < NOIT){
            	
            //Intervalos abertos da loja
            	SN[NOI-1-1] = S1;
            	XRN[NOI-1-1] = X1;
            	SN[NOI-1] = S2;
            	XRN[NOI-1] = X2;	
			} else{
				//Muitos intervalos abertos.
				*CSUMGA_.IERGA = 1;
            	goto L2;
			}
		}
		
		if (*CSUMGA_.NCALL > NCALLT) {
			//muitas chamdas da funcao FCT
			*CSUMGA_.IERGA =2;
			goto L2;	
		}
	}
	
	//An�lise de resultados parciais e controle de erros.

	if (*CSUMGA_.IERGA == 3) { //! Intervals are too narrow.
		
		if (NOI < NOIT5) {
			
			*CSUMGA_.IERGA = 0; //  ! The result is probably correct.
            
			resultado = resultado + SUMR;
            	
            return resultado;
		}
		goto L2;		
	}  
	
	if(*CSUMGA_.IERGA == 0) {
		if ((fabs(SUMR) < fmax(TOL1 *fabs(resultado + SUMR), 1.0e-35)) || (NOI == 0)){
	        return resultado;	
		} else{
			for (int I = 1; I <= NOI; I++){
				S[I-1] = SN[I-1];
            	XR[I-1] = XRN[I-1];	
			}
			goto L1;
		}
	}
	
	//Mensagem de advert�ncia (baixa precis�o).
	
	
	
L2:;
	
	if (IDONE < NOIP){
		for (int I = IDONE+1; I <= NOIP; I++){
			SUMR = SUMR + S[I-1];
		}
		NOI=NOI+(NOIP-IDONE);	
	}
	
	resultado = resultado + SUMR;
	
	if (CSGAWR_.ISGAW == 0)
		return resultado;
	
	printf("\n>>> SUMGA. Gauss adaptive-bisection quadrature.\n");
	printf("  XL = %.8e , XU = %.8e, TOL = %.8e\n", XL, XU, TOL);
	if (fabs(resultado) > 1.0e-35){
		double RERR = fabs(SUMR) / fabs(resultado);
		printf("resultado = %.15e , relative error = %.1e\n", resultado, RERR);
	} else{
		double AERR = fabs(SUMR);
		printf("resultado = %.15e , relative error = %.1e\n", resultado, AERR);
	}
	
	printf("NCALL = %6d, open subintervals =' %4d, H = %.3e\n", *CSUMGA_.NCALL, NOI, HH);
	if (*CSUMGA_.IERGA == 1){
		printf("IERGA = 1, too many open subintervals.\n");
	} else if (*CSUMGA_.IERGA == 2){
		printf("IERGA = 2, too many open subintervals.\n");
		
	} else if (*CSUMGA_.IERGA == 3){
		printf("IERGA = 3, too many open subintervals.\n");
	}
	
	printf("WARNING: the required accuracy has not been attained.\n");
	
//	printf("\n\nSUMGA2\n\n");
	exit(0);
	
	
	return resultado;
}


void pinat2_(double &E, double &WCCM, double &XH0, double &XH1, double &XH2, double &XS0, double &XS1, double &XS2, double &XT1, double &XT2, double &DELTA, int *M){

	
		/*
	Se��es transversais integradas para colis�es inel�sticas de positrons de
  energia E no material M, restrita a perdas de energia maiores que e
  menos do que a energia de corte WCCM.

  Modelo Sternheimer-Liljequist GOS.

  Argumentos de sa�da:
    XH0 ... se��o transversal total para colls r�gidos. (cm ** 2).
    XH1 ... parando a se��o transversal para colls r�gidos. (eV * cm ** 2).
    XH2 ... se��o transversal extensa para colls r�gidos. (eV ** 2 * cm ** 2).
    XS0 ... se��o transversal total para colls macios. (cm ** 2).
    XS1 ... parando a se��o transversal para colls macios. (eV * cm ** 2)
    XS2 ... se��o transversal extensa para colls macios. (eV ** 2 * cm ** 2).
    XT1 ... 1� se��o transversal de transporte para soft colls. (cm ** 2).
    XT2 ... 2� se��o transversal de transporte para soft colls. (cm ** 2).
    DELTA ... Corre��o do efeito de densidade de Fermi.
	
	*/
	double REV = 5.10998928e5; //  ! Electron rest energy (eV)
	double GAM, GAM2, FDEL, WL2, WL2L, WL2U, TST;	
	double H0, H1, H2, S0, S1, S2, R0, R1, R2, XT0;
	double UK, WK;
		
		
    //Constants
	GAM=1.0e0+E/REV;
    GAM2=GAM*GAM;
    
    //Density effect.
    
    //Sternheimer's resonance energy (WL2=L**2).

    TST=COMPOS_.ZT[*M-1]/(GAM2*CEIN_.OP2[*M-1]);
    WL2=0.0e0;
    FDEL=0.0e0;
    
    for (int I = 1; I <= CEIN_.NOSC[*M-1]; I++){
    	FDEL=FDEL+CEIN_.F[I-1][*M-1]/(pow(CEIN_.WRI[I-1][*M-1],2)+WL2);	
	}
	
	if (FDEL < TST){
		DELTA=0.0e0;
        goto L3;
	}
	
    WL2=pow(CEIN_.WRI[CEIN_.NOSC[*M-1]-1][*M-1],2);
L1:;
	WL2=WL2+WL2;
    FDEL=0.0e0;
    
    for (int I = 1; I <= CEIN_.NOSC[*M-1]; I++){
    	FDEL=FDEL+CEIN_.F[I-1][*M-1]/(pow(CEIN_.WRI[I-1][*M-1],2)+WL2);	
	}
	if (FDEL > TST)
		goto L1;
	
	WL2L=0.0e0;
    WL2U=WL2;
L2:;
    WL2=0.5e0*(WL2L+WL2U);
    FDEL=0.0e0;
    
    for (int I = 1; I <= CEIN_.NOSC[*M-1]; I++){
    	FDEL=FDEL+CEIN_.F[I-1][*M-1]/(pow(CEIN_.WRI[I-1][*M-1],2)+WL2);
	}
	
	if (FDEL > TST)
		WL2L=WL2;
	else
        WL2U=WL2;
    
    if ((WL2U-WL2L) > (1.0e-12*WL2)) 
    	goto L2;
    
    DELTA = 0.0e0;
    
    for (int I = 1; I <= CEIN_.NOSC[*M-1]; I++){
    	DELTA=DELTA+CEIN_.F[I-1][*M-1]*log(1.0e0+WL2/pow(CEIN_.WRI[I-1][*M-1], 2));
	}
	DELTA=(DELTA/COMPOS_.ZT[*M-1])-WL2/(GAM2*CEIN_.OP2[*M-1]);
	
L3:;
		
	for (int I = 1; I <=  CEIN_.NOSC[*M-1]; I++){
		CPIN00_.SPH0[I-1] = 0.0;
		CPIN00_.SPH1[I-1] = 0.0;
		CPIN00_.SPH2[I-1] = 0.0;
		CPIN00_.SPS0[I-1] = 0.0;
		CPIN00_.SPS1[I-1] = 0.0;
		CPIN00_.SPS2[I-1] = 0.0;
		CPIN00_.SPT0[I-1] = 0.0;
		CPIN00_.SPT1[I-1] = 0.0;
		CPIN00_.SPT2[I-1] = 0.0;	
	}
	
	XH0=0.0e0;
	XH1=0.0e0;
	XH2=0.0e0;
	XS0=0.0e0;
	XS1=0.0e0;
	XS2=0.0e0;
	XT0=0.0e0;
	XT1=0.0e0;
	XT2=0.0e0;
	
	for (int K = 1; K <= CEIN_.NOSC[*M-1]; K++){
		UK=CEIN_.UI[K-1][*M-1];
        WK=CEIN_.WRI[K-1][*M-1];	
        pinat12_(E, UK, WK, DELTA, WCCM, H0, H1, H2, S0, S1, S2, R0, R1, R2);
        
        CPIN00_.SPH0[K-1] = CEIN_.F[K-1][*M-1] * H0;
		CPIN00_.SPH1[K-1] = CEIN_.F[K-1][*M-1] * H1;
		CPIN00_.SPH2[K-1] = CEIN_.F[K-1][*M-1] * H2;
		CPIN00_.SPS0[K-1] = CEIN_.F[K-1][*M-1] * S0;
		CPIN00_.SPS1[K-1] = CEIN_.F[K-1][*M-1] * S1;
		CPIN00_.SPS2[K-1] = CEIN_.F[K-1][*M-1] * S2;
		CPIN00_.SPT0[K-1] = CEIN_.F[K-1][*M-1] * R0;
		CPIN00_.SPT1[K-1] = CEIN_.F[K-1][*M-1] * 2.0e0 *R1;
		CPIN00_.SPT2[K-1] = CEIN_.F[K-1][*M-1] * 6.0e0 * (R1-R2);
       // CPIN00_.SPS2[K-1] = 10;
	   // CPIN00_.SPT0[K-1] = 20;
	  	//CPIN00_.SPT1[K-1] = 30;
	  //	CPIN00_.SPT2[K-1] = 40;
		
		XH0= XH0 + CPIN00_.SPH0[K-1];
		XH1= XH1 + CPIN00_.SPH1[K-1];
		XH2= XH2 + CPIN00_.SPH2[K-1];
		XS0= XS0 + CPIN00_.SPS0[K-1];
		XS1= XS1 + CPIN00_.SPS1[K-1];
		XS2= XS2 + CPIN00_.SPS2[K-1];
		XT0= XT0 + CPIN00_.SPT0[K-1];
		XT1= XT1 + CPIN00_.SPT1[K-1];
		XT2= XT2 + CPIN00_.SPT2[K-1];
	}	
		
	
	//constantes]
	
	//pinat_(E, WCCM, XH0, XH1, XH2, XS0, XS1, XS2, XT1, XT2, DELTA, M);
//	return;
	
/*	static const int NO = 512;
	
    double H0, H1, H2, S0, S1, S2, R0, R1, R2, XT0;
	double WL2L;
	double WL2U;
	double UK, WK;
	
	double REV = 5.10998928e5; //  ! Electron rest energy (eV)
	
	double GAM = 1.0e0 + E / REV;
	double GAM2 = GAM * GAM;
	
	// Densidade efetiva
	
	//****  Sternheimer's resonance energy (WL2=L**2).
	
	double TST = COMPOS_.ZT[*M-1] / (GAM2 * CEIN_.OP2[*M-1]);
	double WL2 = 0.0e0;
	double FDEL = 0.0e0;
	for (int I = 1; I <= CEIN_.NOSC[*M-1]; I++){
		FDEL = FDEL + CEIN_.F[I-1][*M-1] / (pow(CEIN_.WRI[I-1][*M-1], 2) + WL2);
	}
	
	if (FDEL < TST){
		DELTA = 0.0e0;
		goto L3;
	}
	
	WL2 = pow(CEIN_.WRI[CEIN_.NOSC[*M-1]-1][*M-1], 2);
L1:;
	WL2 = WL2 + WL2;
	FDEL = 0.0e0;

	for (int I = 1; I <= CEIN_.NOSC[*M-1]; I++){
		FDEL = FDEL + CEIN_.F[I-1][*M-1] / (pow(CEIN_.WRI[I-1][*M-1], 2) + WL2);
	}
	
	if (FDEL > TST)
		goto L1;
	WL2L = 0.0e0;
	WL2U = WL2;
L2:;
	WL2 = 0.5e0 * (WL2L + WL2U);
	FDEL = 0.0e0;
	for (int I = 1; I <= CEIN_.NOSC[*M-1]; I++){
		FDEL = FDEL + CEIN_.F[I-1][*M-1] / (pow(CEIN_.WRI[I-1][*M-1], 2) + WL2);
	}
	
	if (FDEL > TST){
		WL2L = WL2;
	}else{
		WL2U = WL2;
	}
	
	if ((WL2U - WL2L) > (1.0e-12 * WL2))
		goto L2;
	
	// Density effect correction (delta).
	
	DELTA = 0.0e0;
	for (int I = 1; I <= CEIN_.NOSC[*M-1]; I++){
		DELTA = DELTA + CEIN_.F[I-1][*M-1] * log(1.0e0 + WL2 / pow(CEIN_.WRI[I-1][*M-1], 2));
	}
	DELTA = (DELTA / COMPOS_.ZT[*M-1]) - WL2 / (GAM2*CEIN_.OP2[*M-1]);

L3:;


//Se��es transversais do oscilador de concha.

	for (int I = 1; I <= CEIN_.NOSC[*M-1]; I++){
		CPIN00_.SPH0[I-1] = 0.0;
		CPIN00_.SPH1[I-1] = 0.0;
		CPIN00_.SPH2[I-1] = 0.0;
		CPIN00_.SPS0[I-1] = 0.0;
		CPIN00_.SPS1[I-1] = 0.0;
		CPIN00_.SPS2[I-1] = 0.0;
		CPIN00_.SPT0[I-1] = 0.0;
		CPIN00_.SPT1[I-1] = 0.0;
		CPIN00_.SPT2[I-1] = 0.0;
	}
	
	
	
	XH0=0.0e0;
	XH1=0.0e0;
	XH2=0.0e0;
	XS0=0.0e0;
	XS1=0.0e0;
	XS2=0.0e0;
	XT0=0.0e0;
	XT1=0.0e0;
	XT2=0.0e0;

	
	for (int K = 1; K <= CEIN_.NOSC[*M-1]; K++){
		UK = CEIN_.UI[K-1][*M-1];
		WK = CEIN_.WRI[K-1][*M-1];

		pinat12_(E, UK, WK, DELTA, WCCM, H0, H1, H2, S0, S1, S2, R0, R1, R2);

		CPIN00_.SPH0[K-1] = CEIN_.F[K-1][*M-1] * H0;
		CPIN00_.SPH1[K-1] = CEIN_.F[K-1][*M-1] * H1;
		CPIN00_.SPH2[K-1] = CEIN_.F[K-1][*M-1] * H2;
		CPIN00_.SPS0[K-1] = CEIN_.F[K-1][*M-1] * S0;
		CPIN00_.SPS1[K-1] = CEIN_.F[K-1][*M-1] * S1;
		CPIN00_.SPS2[K-1] = CEIN_.F[K-1][*M-1] * S2;
		//CPIN00_.SPT0[K-1] = CEIN_.F[K-1][*M-1] * R0;
		CPIN00_.SPT0[K-1] = 10e0;
		CPIN00_.SPT1[K-1] = CEIN_.F[K-1][*M-1] * 2.0e0 *R1;
		CPIN00_.SPT2[K-1] = CEIN_.F[K-1][*M-1] * 6.0e0 * (R1-R2);
		
		XH0= XH0 + CPIN00_.SPH0[K-1];
		XH1= XH1 + CPIN00_.SPH1[K-1];
		XH2= XH2 + CPIN00_.SPH2[K-1];
		XS0= XS0 + CPIN00_.SPS0[K-1];
		XS1= XS1 + CPIN00_.SPS1[K-1];
		XS2= XS2 + CPIN00_.SPS2[K-1];
		XT0= XT0 + CPIN00_.SPT0[K-1];
		XT1= XT1 + CPIN00_.SPT1[K-1];
		XT2= XT2 + CPIN00_.SPT2[K-1];	
	}	
	
//	printf("\n\nPINAT2\n\n");*/
			
}

void pinat12_(double &E, double &UK, double &WK, double &DELTA, double &WCCM, double &H0, double &H1, double &H2, double &S0, double &S1, double &S2, double &R0, double &R1, double &R2){
	
/*
	Se��es transversais integradas para colis�es inel�sticas de positrons com
  um oscilador de camada �nica, restrito a perdas de energia maiores que,
  e menor do que, a perda de energia de corte WCCM.

  Modelo do oscilador Sternheimer-Liljequist.

  Argumentos de entrada:
    E ..... energia cin�tica (eV).
    UK .... energia de ioniza��o (eV).
    WK .... energia de resson�ncia (eV).
    DELTA ... Corre��o do efeito de densidade de Fermi.
    WCCM ... perda de energia de corte (eV).

  Argumentos de sa�da:
    H0 .... se��o transversal total para colls duros. (cm ** 2).
    H1 .... parando a se��o transversal para colls duros. (eV * cm ** 2).
    H2 .... se��o transversal de straggling para hard colls. (eV ** 2 * cm ** 2).
    S0 .... se��o transversal total para colls macios. (cm ** 2).
    S1 .... parando a se��o transversal para colls macios. (eV * cm ** 2).
    S2 .... se��o transversal extensa para colls macios. (eV ** 2 * cm ** 2).
    R0 .... se��o transversal total para colls macios. (cm ** 2).
    R1 .... 1� se��o transversal de transporte para soft colls. (cm ** 2).
    R2 .... 2� se��o transversal de transporte para soft colls. (cm ** 2).
	*/
	
//	pinat1_(E, UK, WK,DELTA, WCCM, H0, H1, H2, S0, S1, S2, R0, R1, R2);
//	return;
	
	double REV = 5.10998928e5;
	double ELRAD= 2.8179403267e-13;
	double TREV= 2.0e0 * REV;
	double RTREV= 1.0e0 / TREV;
	double PIELR2= PI * ELRAD * ELRAD;
	
	double WTHR;
	double WM, WKP, QKP, WCMAX, WDMAX;
	double QM, SDL1, SDT1, CPPS, CPP, A, B, RMU1, BA;
	double F0, F1, F2, SD1, WL, WU;
	double CP2S, CP2, RMU2, CP3S, CP3, RMU3; 
	
	
	H0 = 0.0e0;
	H1 = 0.0e0;
	H2 = 0.0e0;
	S0 = 0.0e0;
	S1 = 0.0e0;
	S2 = 0.0e0;
	R0 = 0.0e0;
	R1 = 0.0e0;
	R2 = 0.0e0;
	
	if (UK > 1.0e-3){
		WTHR = UK;
	} else{
		WTHR = WK;
	}
	
	if (E < (WTHR + 1.0e-6))
		return;
	
	
	//constantes
	
	
	*CPIN01_.EI = E;
	double GAM = 1.0e0 + E / REV;
	double GAM2 = GAM*GAM;
	double BETA2 = (GAM2 - 1.0e0) / GAM2;
	double CONST = PIELR2 * TREV / BETA2;
	
	*CPIN01_.CPS = E * (E+TREV);
	double CP = sqrt(*CPIN01_.CPS);
	double AMOL = pow(E /(E+REV), 2);
	double G12 = pow(GAM + 1.0e0, 2);
	*CPIN01_.BHA1 = AMOL * (2.0e0 * G12 - 1.0e0) / (GAM2 - 1.0e0);
    *CPIN01_.BHA2 = AMOL * (3.0e0 + 1.0e0 / G12);
    *CPIN01_.BHA3 = AMOL * 2.0e0 * GAM * (GAM - 1.0e0) / G12;
    *CPIN01_.BHA4 = AMOL * pow(GAM - 1.0e0, 2) / G12;
	
	 
	 
	//Truque: A energia de resson�ncia e a energia de recuo de corte de
	//As camadas internas  s�o variadas para produzir um limite suave. 
	
	
	if (UK > 1.0e-3){
		WM = 3.0e0 * WK - 2.0e0 * UK;
		if (E > WM){
			WKP = WK;
			QKP = UK;
		} else{
			WKP = (E + 2.0e0 * UK) / 3.0e0;
			QKP = UK * (E/WM);
			WM = E;
			
		}
	//	*CEIN01_.EE = E + UK;
		WCMAX = E;
		WDMAX = fmin(WM, WCMAX);
	}else{
		WM = E;
		WKP = WK;
		QKP = WK;
		WCMAX = E;
		WDMAX = WKP + 1.0e0;	
	}
	
	//Interacoes distantes
	
	SDL1 = 0.0e0;
	SDT1 = 0.0e0;
	if (WDMAX > (WTHR + 1.0e-6)){
		CPPS = (E - WKP) * (E - WKP+TREV);
		CPP = sqrt(CPPS);
		A = 4.0e0 * CP * CPP;
		B = pow(CP - CPP, 2);
		
		if (WKP > (1.0e-6 * E)){
			QM = sqrt(pow(CP - CPP, 2) + pow(REV, 2)) - REV;
		} else{
			QM = WKP * WKP / (BETA2 * TREV);
			QM = QM * (1.0e0 - QM * RTREV);
		}
		
		if (QM < QKP){
			SDL1 = log(QKP * (QM + TREV ) / (QM * (QKP + TREV )));
			SDT1 = fmax(log(GAM2) - BETA2 - DELTA, 0.0e0);
			//Momentos de transporte distantes suaves das ordens 0-2.
			if (WCCM > WTHR){
				BA=B/A;
            	RMU1 = (QKP*(QKP+TREV)-B)/A;
            	R0 = log((RMU1+BA)/BA);
            	R1 = RMU1-BA*R0;
            	R2 = pow(BA,2)*R0+0.5e0*RMU1*(RMU1-2.0e0*BA);
            	R0 = R0/WKP;
            	R1 = R1/WKP;
            	R2 = R2/WKP;
            	R0 = R0+SDT1/WKP;
			}
		}
	}
	
	SD1=SDL1+SDT1;
	
	if (SD1 > 0.0e0){
		if (UK > 1.0e-3){
		//Excita��es internas (distribui��o triangular).	
			F0=1.0e0/pow(WM-UK, 2);
			F1=2.0e0 * F0 * SD1 / WKP;
			
			if (WCCM < UK){
				WL=UK;
            	WU=WDMAX;
            	H0=F1*(WM*(WU-WL)-(pow(WU,2)- pow(WL,2))/2.0e0);
            	H1=F1*(WM*(pow(WU,2)- pow(WL,2))/2.0e0-(pow(WU,3)- pow(WL,3))/3.0e0);
            	H2=F1*(WM*(pow(WU,3)- pow(WL,3))/3.0e0-(pow(WU,4)- pow(WL,4))/4.0e0);
			} else{
				
				if (WCCM > WDMAX){
					WL=UK;
            		WU=WDMAX;
            		S0=F1*(WM*(WU-WL)-(pow(WU,2)- pow(WL,2))/2.0e0);
            		S1=F1*(WM*(pow(WU,2)- pow(WL,2))/2.0e0-(pow(WU,3)- pow(WL,3))/3.0e0);
            		S2=F1*(WM*(pow(WU,3)- pow(WL,3))/3.0e0-(pow(WU,4)- pow(WL,4))/4.0e0);
					
				} else{
					WL=WCCM;
            		WU=WDMAX;
            		H0=F1*(WM*(WU-WL)-(pow(WU,2)- pow(WL,2))/2.0e0);
            		H1=F1*(WM*(pow(WU,2)- pow(WL,2))/2.0e0-(pow(WU,3)- pow(WL,3))/3.0e0);
            		H2=F1*(WM*(pow(WU,3)- pow(WL,3))/3.0e0-(pow(WU,4)- pow(WL,4))/4.0e0);
            		WL=UK;
            		WU=WCCM;
            		S0=F1*(WM*(WU-WL)-(pow(WU,2)- pow(WL,2))/2.0e0);
            		S1=F1*(WM*(pow(WU,2)- pow(WL,2))/2.0e0-(pow(WU,3)- pow(WL,3))/3.0e0);
            		S2=F1*(WM*(pow(WU,3)- pow(WL,3))/3.0e0-(pow(WU,4)- pow(WL,4))/4.0e0);	
				}
				F2=F0*(2.0e0*WM*(WU-WL)-(pow(WU,2)- pow(WL,2)));
            	R0=F2*R0;
            	R1=F2*R1;
            	R2=F2*R2;
			}	
		} else{
			
			//Excita��es externas (oscilador delta).
			if (WCCM < WKP){
				H1=SD1;
            	H0=SD1/WKP;
            	H2=SD1*WKP;
				
			} else{
				S1=SD1;
            	S0=SD1/WKP;
  	 	 	 	S2=SD1*WKP;
			}	
		}
	}
	
	//Fechar colis�es (se��o transversal de Moller).
	
	if (WCMAX < (WTHR + 1.0e-6))
		goto L1;
	
	if (WCCM < WTHR){ //Sem interacoes macias
		WL=WTHR;
        WU=WCMAX;
        H0=H0+(1.0e0/WL)-(1.0e0/WU)
		 -*CPIN01_.BHA1*log(WU/WL)/E
         +*CPIN01_.BHA2*(WU-WL)/pow(E,2)
		 -*CPIN01_.BHA3*(pow(WU,2)-pow(WL,2))/(2.0e0*pow(E,3))
         +*CPIN01_.BHA4*(pow(WU,3)-pow(WL,3))/(3.0e0*pow(E,4));
         
        H1=H1+log(WU/WL)
		 -*CPIN01_.BHA1*(WU-WL)/E
         +*CPIN01_.BHA2*(pow(WU,2)-pow(WL,2))/(2.0e0*pow(E,2))
         -*CPIN01_.BHA3*(pow(WU,3)-pow(WL,3))/(3.0e0*pow(E,3))
         +*CPIN01_.BHA4*(pow(WU,4)-pow(WL,4))/(4.0e0*pow(E,4));
         
        H2=H2+WU-WL
		 -*CPIN01_.BHA1*(pow(WU,2)-pow(WL,2))/(2.0e0*E)
         +*CPIN01_.BHA2*(pow(WU,3)-pow(WL,3))/(3.0e0*pow(E,2))
         -*CPIN01_.BHA3*(pow(WU,4)-pow(WL,4))/(4.0e0*pow(E,3))
         +*CPIN01_.BHA4*(pow(WU,5)-pow(WL,5))/(5.0e0*pow(E,4));
         
	} else{
		if (WCCM > WCMAX){
			WL=WTHR;
        	WU=WCMAX;
        	S0=S0+(1.0e0/WL)-(1.0e0/WU)
		     -*CPIN01_.BHA1*log(WU/WL)/E
         	 +*CPIN01_.BHA2*(WU-WL)/pow(E,2)
		 	 -*CPIN01_.BHA3*(pow(WU,2)-pow(WL,2))/(2.0e0*pow(E,3))
         	 +*CPIN01_.BHA4*(pow(WU,3)-pow(WL,3))/(3.0e0*pow(E,4));
         
		    S1=S1+log(WU/WL)
		     -*CPIN01_.BHA1*(WU-WL)/E
         	 +*CPIN01_.BHA2*(pow(WU,2)-pow(WL,2))/(2.0e0*pow(E,2))
         	 -*CPIN01_.BHA3*(pow(WU,3)-pow(WL,3))/(3.0e0*pow(E,3))
         	 +*CPIN01_.BHA4*(pow(WU,4)-pow(WL,4))/(4.0e0*pow(E,4));
         
        	S2=S2+WU-WL
		 	 -*CPIN01_.BHA1*(pow(WU,2)-pow(WL,2))/(2.0e0*E)
         	 +*CPIN01_.BHA2*(pow(WU,3)-pow(WL,3))/(3.0e0*pow(E,2))
         	 -*CPIN01_.BHA3*(pow(WU,4)-pow(WL,4))/(4.0e0*pow(E,3))
         	 +*CPIN01_.BHA4*(pow(WU,5)-pow(WL,5))/(5.0e0*pow(E,4));
	        
		} else{
			WL=WCCM;
	        WU=WCMAX;
	        H0=H0+(1.0e0/WL)-(1.0e0/WU)
		 	 -*CPIN01_.BHA1*log(WU/WL)/E
         	 +*CPIN01_.BHA2*(WU-WL)/pow(E,2)
		 	 -*CPIN01_.BHA3*(pow(WU,2)-pow(WL,2))/(2.0e0*pow(E,3))
         	 +*CPIN01_.BHA4*(pow(WU,3)-pow(WL,3))/(3.0e0*pow(E,4));
         
        	H1=H1+log(WU/WL)
		     -*CPIN01_.BHA1*(WU-WL)/E
         	 +*CPIN01_.BHA2*(pow(WU,2)-pow(WL,2))/(2.0e0*pow(E,2))
         	 -*CPIN01_.BHA3*(pow(WU,3)-pow(WL,3))/(3.0e0*pow(E,3))
         	 +*CPIN01_.BHA4*(pow(WU,4)-pow(WL,4))/(4.0e0*pow(E,4));
         
        	H2=H2+WU-WL
		 	 -*CPIN01_.BHA1*(pow(WU,2)-pow(WL,2))/(2.0e0*E)
         	 +*CPIN01_.BHA2*(pow(WU,3)-pow(WL,3))/(3.0e0*pow(E,2))
         	 -*CPIN01_.BHA3*(pow(WU,4)-pow(WL,4))/(4.0e0*pow(E,3))
         	 +*CPIN01_.BHA4*(pow(WU,5)-pow(WL,5))/(5.0e0*pow(E,4));
         	 
         	 WL=WTHR;
 	         WU=WCCM;
 	         
 	         S0=S0+(1.0e0/WL)-(1.0e0/WU)
		    -*CPIN01_.BHA1*log(WU/WL)/E
         	 +*CPIN01_.BHA2*(WU-WL)/pow(E,2)
		 	 -*CPIN01_.BHA3*(pow(WU,2)-pow(WL,2))/(2.0e0*pow(E,3))
         	 +*CPIN01_.BHA4*(pow(WU,3)-pow(WL,3))/(3.0e0*pow(E,4));
         
		    S1=S1+log(WU/WL)
			-*CPIN01_.BHA1*(WU-WL)/E
         	 +*CPIN01_.BHA2*(pow(WU,2)-pow(WL,2))/(2.0e0*pow(E,2))
         	 -*CPIN01_.BHA3*(pow(WU,3)-pow(WL,3))/(3.0e0*pow(E,3))
         	 +*CPIN01_.BHA4*(pow(WU,4)-pow(WL,4))/(4.0e0*pow(E,4));
         
        	S2=S2+WU-WL
		 	 -*CPIN01_.BHA1*(pow(WU,2)-pow(WL,2))/(2.0e0*E)
         	 +*CPIN01_.BHA2*(pow(WU,3)-pow(WL,3))/(3.0e0*pow(E,2))
         	 -*CPIN01_.BHA3*(pow(WU,4)-pow(WL,4))/(4.0e0*pow(E,3))
         	 +*CPIN01_.BHA4*(pow(WU,5)-pow(WL,5))/(5.0e0*pow(E,4));
 	         		
		}
		
		//Momentos de transporte de fechamento suave das ordens 0-2.
		
    	CP2S=(E-WL)*(E-WL+TREV);
        CP2=sqrt(CP2S);
        RMU2=(WL*(WL+TREV)-pow(CP-CP2,2))/(4.0e0*CP*CP2);
        if (WU < (E - 1.0e0)){
			CP3S=(E-WU)*(E-WU+TREV);
        	CP3=sqrt(CP3S);
        	RMU3=(WU*(WU+TREV)-pow(CP-CP3,2))/(4.0e0*CP*CP3);
		} else{
			RMU3 = 0.5e0;
			
		}
        *CPIN01_.MOM = 0;
        int FUNCAO = 2;
        double TOL = 1.0e-7;
        R0= R0+ sumga2_(FUNCAO, RMU2,RMU3,TOL);
        *CPIN01_.MOM = 1;
        R1 = R1+ sumga2_(FUNCAO, RMU2,RMU3,TOL);
        *CPIN01_.MOM = 2;
        R2 = R2+ sumga2_(FUNCAO, RMU2,RMU3,TOL);
	}
	
L1:;

    H0=CONST*H0;
	H1=CONST*H1;
	H2=CONST*H2;
	S0=CONST*S0;
	S1=CONST*S1;
	S2=CONST*S2;
	R0=CONST*R0;
	R1=CONST*R1;
	R2=CONST*R2;
	
//	printf("\n\nPINAT12\n\n");
//	exit(0);
	
}

double pinads2_(double RMU){
	
	/*
	Se��o transversal diferencial angular para soft close inel�stico colis�es de positrons.
	*/
	
	double REV = 5.10998928e5; 
	
	double resultado;
	
	double AUX=2.0e0*RMU*(1.0e0-RMU);
	double DENOM= *CPIN01_.EI * AUX + REV;
	double W= *CPIN01_.CPS * AUX / DENOM;
	double DWDMU = *CPIN01_.CPS * REV * (2.0e0-4.0e0*RMU) / pow(DENOM,2);
	double WE = W / *CPIN01_.EI;
	
//	printf("\n\nPINADS2\n\n");
	
	resultado = (1.0e0 - WE * ( *CPIN01_.BHA1 - WE *( *CPIN01_.BHA2 - WE *( *CPIN01_.BHA3 - WE * *CPIN01_.BHA4)))) * DWDMU * pow(RMU,*CPIN01_.MOM)/ pow(W,2);
	return resultado;
}


void ebrat2_(double &E, double &WCRM, double &XH0, double &XH1, double &XH2, double &XS1, double &XS2, int *M){
	
	/*
	Se��es transversais integradas para emiss�o de brems por el�trons de energia
  E no material M, restrito a perdas de energia maiores e menores
  do que a energia de corte WCRM.

  Argumentos de sa�da:
    XH0 ... se��o transversal total para emiss�o dura (cm ** 2).
    XH1 ... se��o transversal de parada para emiss�o forte (eV * cm ** 2).
    XH2 ... se��o transversal difusa para emiss�o dura (eV ** 2 * cm ** 2).
    XS1 ... se��o transversal de parada para emiss�o suave (eV * cm ** 2).
    XS2 ... se��o transversal difusa para emiss�o suave (eV ** 2 * cm ** 2).
	
	*/
	static const int NBE = 57;
	static const int NBW = 32;
	
	double REV = 5.10998928e5;  //Electron rest energy (eV)
	double TREV = 2.0e0 * REV;
	
	*CEGRID_.XEL= fmax(log(E), *CEGRID_.DLEMP1);
	*CEGRID_.XE= 1.0e0+(*CEGRID_.XEL - *CEGRID_.DLEMP1) * *CEGRID_.DLFC;
	*CEGRID_.KE= *CEGRID_.XE;
	*CEGRID_.XEK= *CEGRID_.XE - *CEGRID_.KE;
	
//Fator de se��o x global.

	double FACT = CEBR_.ZBR2[*M-1] * (pow( E + REV,2) / ( E * ( E + TREV))) * 1.0e-27;
	
	//Momentos da se��o x dos bremss em escala.
	
	double WCRE=WCRM/E;
	
	for (int IW = 1; IW <= NBW; IW++){
		CEBR01_.X[IW-1]=CEBR_.WB[IW-1];
        CEBR01_.Y[IW-1]=CEBR02_.P0[IW-1][*CEGRID_.KE-1][*M-1];	
      //  printf("C++ KE: %d, XE: %f", *CEGRID_.KE, *CEGRID_.XE);
      //  printf("\nC++ P0: %.15f\n", CEBR02_.P0[IW-1][*CEGRID_.KE-1][*M-1]);
	}
	
	double XH0A = rlmom2_(CEBR01_.X, CEBR01_.Y, CEBR01_.X[NBW-1], NBW, -1) - rlmom2_(CEBR01_.X,CEBR01_.Y, WCRE, NBW, -1);
	double XS1A = rlmom2_(CEBR01_.X, CEBR01_.Y, WCRE, NBW, 0);
	double XS2A = rlmom2_(CEBR01_.X, CEBR01_.Y, WCRE, NBW, 1);
	double XH1A = rlmom2_(CEBR01_.X, CEBR01_.Y, CEBR01_.X[NBW-1], NBW, 0) - XS1A;
	double XH2A = rlmom2_(CEBR01_.X, CEBR01_.Y, CEBR01_.X[NBW-1], NBW, 1) - XS2A;
    
    for (int IW = 1; IW <= NBW; IW++){
    	CEBR01_.Y[IW-1] = CEBR02_.P0[IW-1][min(*CEGRID_.KE+1, NEGP)-1][*M-1];
	}
	
	double XH0B = rlmom2_(CEBR01_.X,CEBR01_.Y, CEBR01_.X[NBW-1], NBW, -1) - rlmom2_(CEBR01_.X,CEBR01_.Y, WCRE, NBW, -1);
	double XS1B = rlmom2_(CEBR01_.X,CEBR01_.Y, WCRE, NBW, 0);
	double XS2B = rlmom2_(CEBR01_.X,CEBR01_.Y, WCRE, NBW, 1);

	double XH1B = rlmom2_(CEBR01_.X,CEBR01_.Y, CEBR01_.X[NBW-1], NBW, 0) - XS1B;
	double XH2B = rlmom2_(CEBR01_.X,CEBR01_.Y, CEBR01_.X[NBW-1], NBW, 1) - XS2B;
	
//		printf("\nWCRE %.15e, XS1A: %.15e, XS1B: %.15e, XS2A: %.15e, XS2B: %.15e, M: %d\n", WCRE, XS1A, XS1B, XS2A, XS2B, *M);
	
	XH0 = ((1.0e0 - *CEGRID_.XEK) * XH0A + *CEGRID_.XEK * XH0B) * FACT;
	XS1 = ((1.0e0 - *CEGRID_.XEK) * XS1A + *CEGRID_.XEK * XS1B) * FACT * E;
	XH1 = ((1.0e0 - *CEGRID_.XEK) * XH1A + *CEGRID_.XEK * XH1B) * FACT * E;
	XS2 = ((1.0e0 - *CEGRID_.XEK) * XS2A + *CEGRID_.XEK * XS2B) * FACT * E * E;
	XH2 = ((1.0e0 - *CEGRID_.XEK) * XH2A + *CEGRID_.XEK * XH2B) * FACT * E * E;
	
//	printf("\nXS1: %.15e, XS2: %.15e", XS1, XS2);
	
//	printf("\n\nEBRAT2\n\n");
	
}


void sinteg2_(double *X, double *A, double *B, double *C, double *D,  double &XL, double &XU, double &SUM, int &N){
	/*
	Calcula a integral de uma fun��o spline c�bica.

  Entrada:
     X (I) (I = 1: N) ... pontos de grade (os valores de X devem ser crescentes
                      pedido).
     A (I), B (I), (I), D (I) (I = 1: N) ... coeficientes spline.
     N ............. n�mero de pontos da grade.
     XL ............. limite inferior da integral.
     XU ............. limite superior da integral.
  Sa�da:
     SUM ............ valor do integral.
	*/
	
	double XLL, XUU, SIGN, X1, X2, SUMP;
	int IL, IU;
	
	//Defina os limites de integra��o em ordem crescente.
	if (XU > XL){
		XLL = XL;
		XUU = XU;
		SIGN = 1.0e0;
	}else{
		XLL=XU;
        XUU=XL;
        SIGN=-1.0e0;
	}
	
	//Verifique os limites integrais.
	if ((XLL < X[1-1]) || (XUU > X[N-1])) {
		printf("ntegral limits out of range. Stop.\n");
		exit(0);
	}
	
	//Encontre os intervalos envolvidos.
	
	SUM=0.0e0;
	findi2_(X,&XLL,&N,&IL);
	findi2_(X,&XUU,&N,&IU);
	
	if (IL == IU){
		//Apenas um �nico intervalo envolvido.
		X1=XLL;
        X2=XUU;
        SUM= X2*(A[IL-1]+X2*((B[IL-1]/2)+X2*((C[IL-1]/3)+X2*D[IL-1]/4)))
            -X1*(A[IL-1]+X1*((B[IL-1]/2)+X1*((C[IL-1]/3)+X1*D[IL-1]/4)));
	} else{
		//Contribui��es de v�rios intervalos.7
		X1=XLL;
        X2=X[IL+1-1];
        SUM= X2*(A[IL-1]+X2*((B[IL-1]/2)+X2*((C[IL-1]/3)+X2*D[IL-1]/4)))
            -X1*(A[IL-1]+X1*((B[IL-1]/2)+X1*((C[IL-1]/3)+X1*D[IL-1]/4)));
        IL=IL+1;
        
        for (int I = IL; IL <= IU; IL++){
  			X1=X[I-1];
            X2=X[I+1-1];
            if (I == IU) 
				X2=XUU;
            SUMP= X2*(A[I-1]+X2*((B[I-1]/2)+X2*((C[I-1]/3)+X2*D[I-1]/4)))
                 -X1*(A[I-1]+X1*((B[I-1]/2)+X1*((C[I-1]/3)+X1*D[I-1]/4)));
            SUM=SUM+SUMP;
		}
	}
	
	SUM=SIGN*SUM;
	
//	printf("\n\nSINTEG2\n\n");
		
}


double rmomx2_(double *X, double *PDF, double XD, double XU, int NP, int MOM){
	
	/*
	C�lculo de momentos de um pdf, PDF (X), obtido a partir de
  interpola��o log-log da tabela de entrada. A vari�vel independente X
  assume-se que assume apenas valores positivos.

     X ........ array de valores de vari�veis ??(em ordem crescente).
     PDF ...... valores de PDF correspondentes (n�o deve ser negativo).
     NP ....... n�mero de pontos na tabela.
     XD, XU ... limites do intervalo de integra��o.
     MAM� ... ordem do momento.
     RMOM = INTEGRAL (X ** N) * PDF (X) dX no intervalo (XD, XU).
	*/
	double EPS = 1.0e-12;
	double ZERO = 1.0e-35;
	
	double XLOW, XUP, XIL, XFL, YIL, YFL, X1, X2, DENOM, Y1, Y2, DXL, DYL, AP1, DSUM;
	int IL, IU;
	double resultado = 0.0;
	
	
	if (NP < 2){
		printf("RMOMX. Error code 1.\n");
		exit(0);
	}
	
	if ((X[1-1] < 0.0e0) || (PDF[1-1] < 0.0e0)){
		printf("X(1),PDF(1) = %f, %f    RMOMX. Error code 2.\n", X[1-1],PDF[1-1]);
		exit(0);
	}
	
	for (int I = 2; I<=NP; I++){
		if ((X[I-1] < 0.0e0) || (PDF[I-1] < 0.0e0)){
		   	printf("X(1),PDF(1) = %f, %f    RMOMX. Error code 2.\n", X[I-1],PDF[I-1]);
		   	exit(0);
	   	}
	   	if (X[I-1] < X[I-1-1]){
	   		printf("RMOMX. Error code 4.\n");
		   	exit(0);
			  
		} 	
	}
	
	XLOW=fmax(X[1-1],XD);
	if (XLOW < ZERO)
		XLOW=ZERO;
	XUP=fmin(X[NP-1],XU);
	
	if (XLOW >= XUP){
		printf("WARNING: XLOW is greater than XUP in RMOMX.\n");
		printf("XLOW = %f, XUP = %f\n", XLOW, XUP);
		resultado = 0.0e0;
		return resultado;	
	}
	
	IL=1;
	IU=NP-1;
	
	for (int I = 1; I <= NP-1; I++){
		if (X[I-1] < XLOW) 
			IL=I;
        if (X[I-1] < XUP)
			IU=I;
	}
	
	//Um �nico intervalo.

	if (IU == IL){
		XIL=log(fmax(X[IL-1],ZERO));
        XFL=log(X[IL+1-1]);
        YIL=log(fmax(PDF[IL-1],ZERO));
        YFL=log(fmax(PDF[IL+1-1],ZERO));
        X1=XLOW;
        X2=XUP;
        DENOM=XFL-XIL;
        if (fabs(DENOM) > ZERO){
        	Y1=exp(YIL+(YFL-YIL)*(log(X1)-XIL)/DENOM)*pow(X1,MOM);
 	        Y2=exp(YIL+(YFL-YIL)*(log(X2)-XIL)/DENOM)*pow(X2,MOM);	
		}else{
			Y1=exp(YIL)* pow(X1,MOM);
            Y2=exp(YIL)* pow(X2,MOM);
		}
		DXL=log(X2)-log(X1);
        DYL=log(fmax(Y2,ZERO))-log(fmax(Y1,ZERO));
		
		if (fabs(DXL) > EPS*fabs(DYL)){
			AP1=1.0e0+(DYL/DXL);
			if (fabs(AP1) > EPS){
				DSUM=(Y2*X2-Y1*X1)/AP1;
			}else{
				DSUM=Y1*X1*DXL;	
			}
		}else{
			DSUM=0.5e0*(Y1+Y2)*(X2-X1);
		}
		resultado = DSUM;
		return resultado;
	}
	
	//Multiplos intervalos
	
	XIL=log(fmax(X[IL-1],ZERO));
    XFL=log(X[IL+1-1]);
    YIL=log(fmax(PDF[IL-1],ZERO));
    YFL=log(fmax(PDF[IL+1-1],ZERO));
    X1=XLOW;
    DENOM=XFL-XIL;
    if (fabs(DENOM) > ZERO){
  		Y1=exp(YIL+(YFL-YIL)*(log(X1)-XIL)/DENOM)*pow(X1,MOM);
	} else{
		Y1=exp(YIL)* pow(X1,MOM);	
	}
	X2=X[IL+1-1];
	Y2=fmax(PDF[IL+1-1],ZERO)* pow(X2,MOM);
	DXL=log(X2)-log(X1);
    DYL=log(fmax(Y2,ZERO))-log(fmax(Y1,ZERO));
    if (fabs(DXL) > EPS*fabs(DYL)){
		AP1=1.0e0+(DYL/DXL);
		if (fabs(AP1) > EPS){
			DSUM=(Y2*X2-Y1*X1)/AP1;
		}else{
			DSUM=Y1*X1*DXL;	
		}
		
	}else{
		DSUM=0.5e0*(Y1+Y2)*(X2-X1);
	}
	resultado = DSUM;
	
	if (IU > IL+1){
		for (int I = IL+1; I <= IU-1; I++){
			X1=X[I-1];
            Y1=fmax(PDF[I-1],ZERO)* pow(X1,MOM);
            X2=X[I+1-1];
            Y2=fmax(PDF[I+1-1],ZERO)* pow(X2,MOM);
            DXL=log(X2)-log(X1);
    		DYL=log(fmax(Y2,ZERO))-log(fmax(Y1,ZERO));
    		if (fabs(DXL) > EPS*fabs(DYL)){
    			AP1=1.0e0+(DYL/DXL);
    			if (fabs(AP1) > EPS){
			   	   	DSUM=(Y2*X2-Y1*X1)/AP1;
			   	}else{
			   	   	DSUM=Y1*X1*DXL;	
		   	   	}
   		   	}else{
				DSUM=0.5e0*(Y1+Y2)*(X2-X1);
	   	   	}
	   	   	resultado = resultado + DSUM;		
		}	
	}
	
	X1=X[IU-1];
	Y1=fmax(PDF[IU-1],ZERO)*pow(X1,MOM);
	XIL=log(X[IU-1]);
	XFL=log(X[IU+1-1]);
	YIL=log(fmax(PDF[IU-1],ZERO));
	YFL=log(fmax(PDF[IU+1-1],ZERO));
	X2=XUP;
	DENOM=XFL-XIL;
	
	if (fabs(DENOM) > ZERO){
		Y2=exp(YIL+(YFL-YIL)*(log(X2)-XIL)/DENOM)*pow(X2,MOM);
	}else{
		Y2=exp(YIL)*pow(X2,MOM);
	}
	
	DXL=log(X2)-log(X1);
    DYL=log(fmax(Y2,ZERO))-log(fmax(Y1,ZERO));
    
    if (fabs(DXL) > EPS*fabs(DYL)){
		AP1=1.0e0+(DYL/DXL);
    	if (fabs(AP1) > EPS){
		   	DSUM=(Y2*X2-Y1*X1)/AP1;
		}else{
		   	DSUM=Y1*X1*DXL;	
		}
	}else{
		DSUM=0.5e0*(Y1+Y2)*(X2-X1);
	}
	resultado = resultado + DSUM;
	
//	printf("\n\nRMOMX2\n\n");
	
	
	return resultado;

	
}


void pbrat2_(double &E, double &WCRM, double &XH0, double &XH1, double &XH2, double &XS1, double &XS2, int *M){
	
	/*
	Se��es transversais integradas para emiss�o de brems por positrons de energia
  E no material M, restrito a perdas de energia maiores e menores
  do que a energia de corte WCRM.

  Argumentos de sa�da:
    XH0 ... se��o transversal total para emiss�o dura (cm ** 2).
    XH1 ... se��o transversal de parada para emiss�o forte (eV * cm ** 2).
    XH2 ... se��o transversal difusa para emiss�o dura (eV ** 2 * cm ** 2).
    XS1 ... se��o transversal de parada para emiss�o suave (eV * cm ** 2).
    XS2 ... se��o transversal difusa para emiss�o suave (eV ** 2 * cm ** 2).
	
	*/
	
	
	double REV = 5.10998928e5;  //Electron rest energy (eV)
	double TREV = 2.0e0 * REV;
	
	static const int NBE = 57;
	static const int NBW = 32;
	
	
	
	*CEGRID_.XEL= fmax(log(E), *CEGRID_.DLEMP1);
	*CEGRID_.XE= 1.0e0+(*CEGRID_.XEL - *CEGRID_.DLEMP1) * *CEGRID_.DLFC;
	*CEGRID_.KE= *CEGRID_.XE;
	*CEGRID_.XEK= *CEGRID_.XE - *CEGRID_.KE;
	
//Fator de se��o x global.

	double FACT = CEBR_.ZBR2[*M-1] * (pow( E + REV,2) / ( E * ( E + TREV))) * 1.0e-27;
	
	//Fator de corre��o de p�sitrons.
	double T=log(1.0e0+1.0e6*E/(REV * CEBR_.ZBR2[*M-1]));
	double FPOS=1.0e0-exp(-T*(1.2359e-1-T*(6.1274e-2-T*(3.1516e-2-T
         	    *(7.7446e-3-T*(1.0595e-3-T*(7.0568e-5-T*1.8080e-6)))))));
	FACT=FACT*FPOS;
	
	//Momentos da se��o x dos bremss em escala.
	
	double WCRE=WCRM/E;
	
	for (int IW = 1; IW <= NBW; IW++){
		CEBR01_.X[IW-1]=CEBR_.WB[IW-1];
        CEBR01_.Y[IW-1]=CEBR02_.P0[IW-1][*CEGRID_.KE-1][*M-1];	
	}
	
	double XH0A = rlmom2_(CEBR01_.X, CEBR01_.Y, CEBR01_.X[NBW-1], NBW, -1) - rlmom2_(CEBR01_.X,CEBR01_.Y, WCRE, NBW, -1);
	double XS1A = rlmom2_(CEBR01_.X, CEBR01_.Y, WCRE, NBW, 0);
	double XS2A = rlmom2_(CEBR01_.X, CEBR01_.Y, WCRE, NBW, 1);
	double XH1A = rlmom2_(CEBR01_.X, CEBR01_.Y, CEBR01_.X[NBW-1], NBW, 0) - XS1A;
	double XH2A = rlmom2_(CEBR01_.X, CEBR01_.Y, CEBR01_.X[NBW-1], NBW, 1) - XS2A;
    
    for (int IW = 1; IW <= NBW; IW++){
    	CEBR01_.Y[IW-1] = CEBR02_.P0[IW-1][min(*CEGRID_.KE+1, NEGP)-1][*M-1];
	}
	
	double XH0B = rlmom2_(CEBR01_.X,CEBR01_.Y, CEBR01_.X[NBW-1], NBW, -1) - rlmom2_(CEBR01_.X,CEBR01_.Y, WCRE, NBW, -1);
	double XS1B = rlmom2_(CEBR01_.X,CEBR01_.Y, WCRE, NBW, 0);
	double XS2B = rlmom2_(CEBR01_.X,CEBR01_.Y, WCRE, NBW, 1);
	double XH1B = rlmom2_(CEBR01_.X,CEBR01_.Y, CEBR01_.X[NBW-1], NBW, 0) - XS1B;
	double XH2B = rlmom2_(CEBR01_.X,CEBR01_.Y, CEBR01_.X[NBW-1], NBW, 1) - XS2B;
	
	XH0 = ((1.0e0 - *CEGRID_.XEK) * XH0A + *CEGRID_.XEK * XH0B) * FACT;
	XS1 = ((1.0e0 - *CEGRID_.XEK) * XS1A + *CEGRID_.XEK * XS1B) * FACT * E;
	XH1 = ((1.0e0 - *CEGRID_.XEK) * XH1A + *CEGRID_.XEK * XH1B) * FACT * E;
	XS2 = ((1.0e0 - *CEGRID_.XEK) * XS2A + *CEGRID_.XEK * XS2B) * FACT * E * E;
	XH2 = ((1.0e0 - *CEGRID_.XEK) * XH2A + *CEGRID_.XEK * XH2B) * FACT * E * E;
	
//	printf("\n\nPBRAT2\n\n");
	
}


void panat2_(double &E, double &TXS){
	
	/*
	Se��o transversal total (por el�tron) para aniquila��o de p�sitrons com
  energia cin�tica E. Calculada a partir da f�rmula dcs de Heitler para aniquila-
  com el�trons livres em repouso.

  Argumento de sa�da:
    XST ... se��o transversal de aniquila��o total (cm ** 2).
	*/
	
	double REV=5.10998928e5;
    double	ELRAD=2.8179403267e-13;
	double PI=3.1415926535897932e0;
	double PIELR2=PI*ELRAD*ELRAD;
	
	double GAM=1.0e0+fmax(E,1.0e0)/REV;
	double GAM2=GAM*GAM;
	double F2=GAM2-1.0e0;
	double F1=sqrt(F2);
	TXS=PIELR2*((GAM2+4.0e0*GAM+1.0e0)*log(GAM+F1)/F2
        -(GAM+3.0e0)/F1)/(GAM+1.0e0);
        
  // printf("\n\nPNAT2\n\n");
	
	
}


void eelar2_(int *M, FILE *IRD, FILE *IWR, int *INFO){
	/*
	Esta sub-rotina l� se��es transversais el�sticas para el�trons e positrons
	 no material M do arquivo de dados do material. Ele tamb�m inicializa
    o algoritmo para simula��o de espalhamento el�stico de el�trons e P�sitrons 
	*/
	
	//Lendo a tabela de se��o transversal de entrada.
	
	int NDATA, J;
	double EC, XS0, XS1, XS2, FPEL, FPT1, FPST, XS0H, AAA, BBB, RNDC, XS1S, XS2S, FP0, FP1, HMFP;
	
	fgets(LINHA, sizeof(LINHA), IRD);
	extrairString(APOIO, LINHA, 59, 4);
	NDATA = atoi(APOIO);
	
	if (*INFO >= 2){
		fprintf(IWR, "***\n Electron and positron elastic cross sections,   NDATA = %4d\n\n", NDATA);
		fprintf(IWR, "  Energy       CS0,e-      CS1,e-      CS2,e-      CS0,e+      CS1,e+      CS2,e+\n");
		fprintf(IWR, "   (eV)        (cm**2)     (cm**2)     (cm**2)     (cm**2)     (cm**2)     (cm**2)\n");
		fprintf(IWR, " -----------------------------------------------------------------------------------------\n");
	}
	
	
//	printf("\n\nNDATA: %d\n\n", NDATA);
	
	
	for (int I = 1; I <= NDATA; I++){
		fgets(LINHA, sizeof(LINHA), IRD);
		extrairString(APOIO, LINHA, 0, 10);
	    CEEL00_.EJT[I-1] = atof(APOIO);
	    extrairString(APOIO, LINHA, 10, 12);
	    CEEL00_.XE0[I-1] = atof(APOIO);
	    extrairString(APOIO, LINHA, 22, 12);
	    CEEL00_.XE1[I-1] = atof(APOIO);
	    extrairString(APOIO, LINHA, 34, 12);
	    CEEL00_.XE2[I-1] = atof(APOIO);
	    extrairString(APOIO, LINHA, 46, 12);
	    CEEL00_.XP0[I-1] = atof(APOIO);
	    extrairString(APOIO, LINHA, 58, 12);
	    CEEL00_.XP1[I-1] = atof(APOIO);
	    extrairString(APOIO, LINHA, 70, 12);
	    CEEL00_.XP2[I-1] = atof(APOIO);
	    
	    
	    if (*INFO >= 3){
			fprintf(IWR, "%.3E %.5E %.5E %.5E %.5E %.5E %.5E\n", CEEL00_.EJT[I-1],CEEL00_.XE0[I-1],CEEL00_.XE1[I-1],CEEL00_.XE2[I-1],CEEL00_.XP0[I-1],CEEL00_.XP1[I-1],CEEL00_.XP2[I-1] );
		}
		CEEL00_.EJTL[I-1]=log(CEEL00_.EJT[I-1]);
	}
	
	
	
	//Espalhamento el�stico de el�trons.
	
	for (int I = 1; I <= NDATA; I++){
		CEEL00_.FJL[I-1]=log(CEEL00_.XE0[I-1]);
	}
	spline2_(CEEL00_.EJTL,CEEL00_.FJL,CEEL00_.A,CEEL00_.B,CEEL00_.C,CEEL00_.D,0.0e0,0.0e0,NDATA);
//	double www1 = 0.0e0;
//	double www2 = 0.0e0;
	
//	spline_(CEEL00_.EJTL,CEEL00_.FJL,CEEL00_.A,CEEL00_.B,CEEL00_.C,CEEL00_.D,www1,www2,NDATA);
	for (int I = 1; I<= NEGP; I++){
		EC= CEGRID_.DLEMP[I-1];
		findi2_(CEEL00_.EJTL,&EC,&NDATA,&J);
	//	findi_(CEEL00_.EJTL,EC,NDATA,J);
		CEEL00_.XE0[I-1]=exp(CEEL00_.A[J-1]+EC*(CEEL00_.B[J-1]+EC*(CEEL00_.C[J-1]+EC*CEEL00_.D[J-1])));
	}
	
	for (int I = 1; I <= NDATA; I++){
		CEEL00_.FJL[I-1]=log(CEEL00_.XE1[I-1]);
	}
	spline2_(CEEL00_.EJTL,CEEL00_.FJL,CEEL00_.A,CEEL00_.B,CEEL00_.C,CEEL00_.D,0.0e0,0.0e0,NDATA);
	for (int I = 1; I<= NEGP; I++){
		EC= CEGRID_.DLEMP[I-1];
		findi2_(CEEL00_.EJTL,&EC,&NDATA,&J);
		CEEL00_.XE1[I-1]=exp(CEEL00_.A[J-1]+EC*(CEEL00_.B[J-1]+EC*(CEEL00_.C[J-1]+EC*CEEL00_.D[J-1])));
	}
	
	for (int I = 1; I <= NDATA; I++){
		CEEL00_.FJL[I-1]=log(CEEL00_.XE2[I-1]);
	}
	
	spline2_(CEEL00_.EJTL,CEEL00_.FJL,CEEL00_.A,CEEL00_.B,CEEL00_.C,CEEL00_.D,0.0e0,0.0e0,NDATA);
	for (int I = 1; I<= NEGP; I++){
		EC= CEGRID_.DLEMP[I-1];
		findi2_(CEEL00_.EJTL,&EC,&NDATA,&J);
		CEEL00_.XE2[I-1]=exp(CEEL00_.A[J-1]+EC*(CEEL00_.B[J-1]+EC*(CEEL00_.C[J-1]+EC*CEEL00_.D[J-1])));
	}
	
	
	
	for (int I = 1; I <= NEGP; I++){
		XS0=CEEL00_.XE0[I-1];
        XS1=CEEL00_.XE1[I-1];
        XS2=CEEL00_.XE2[I-1];
        FPEL=1.0e0/(XS0*COMPOS_.VMOL[*M-1]);
        FPT1=1.0e0/(XS1*COMPOS_.VMOL[*M-1]);
        FPST=CEGRID_.ET[I-1]/(CEIMFP_.CSTPE[I-1][*M-1]+CEIMFP_.RSTPE[I-1][*M-1]);
        XS0H=1.0e0/(COMPOS_.VMOL[*M-1] * fmax(FPEL, fmin(PENELOPE_mod_.C1[*M-1] * FPT1,PENELOPE_mod_.C2[*M-1]*FPST)));
		eela02_(XS0,XS1,XS2,XS0H,AAA,BBB,RNDC,XS1S,XS2S);
        CEIMFP_.SEHEL[I-1][*M-1]=XS0H*COMPOS_.VMOL[*M-1];
        CEIMFP_.RNDCE[I-1][*M-1]=RNDC;
        CEIMFP_.AE[I-1][*M-1]=AAA;
        CEIMFP_.BE[I-1][*M-1]=BBB;
        CEEL00_.T1E0[I-1]=XS1S;
        CEIMFP_.T1E[I-1][*M-1]=CEINTF_.T1EI[I-1]+XS1S*COMPOS_.VMOL[*M-1];
        CEEL00_.T2E0[I-1]=XS2S;
        CEIMFP_.T2E[I-1][*M-1]=CEINTF_.T2EI[I-1]+XS2S*COMPOS_.VMOL[*M-1];	
	}
	

	//Imprima tabelas de espalhamento el�stico de el�trons.
	
	if (*INFO >= 3){
		fprintf(IWR, "\n PENELOPE >>>  Elastic scattering of electrons\n");
		fprintf(IWR, "\n   E (eV)      MFP (mtu)   TMFP1 (mtu)  MFPh (mtu)        A           B           RNDC\n");
     	fprintf(IWR, " -----------------------------------------------------------------------------------------\n");
	}
	
	for (int I = 1; I <= NEGP; I++){
		FP0=COMPOS_.RHO[*M-1]/(CEEL00_.XE0[I-1]*COMPOS_.VMOL[*M-1]);
        FP1=COMPOS_.RHO[*M-1]/(CEEL00_.XE1[I-1]*COMPOS_.VMOL[*M-1]);
        HMFP=COMPOS_.RHO[*M-1]/CEIMFP_.SEHEL[I-1][*M-1];
        if (*INFO >= 3){
			fprintf(IWR, "%.5E %.5E %.5E %.5E %.5E %.5E %.5E\n", 
			CEGRID_.ET[I-1],FP0,FP1,HMFP,CEIMFP_.AE[I-1][*M-1],CEIMFP_.BE[I-1][*M-1],CEIMFP_.RNDCE[I-1][*M-1]);
		}
		CEIMFP_.SEHEL[I-1][*M-1]=log(CEIMFP_.SEHEL[I-1][*M-1]);
        CEIMFP_.AE[I-1][*M-1]=log(CEIMFP_.AE[I-1][*M-1]);
        
        //Os eventos de dispers�o suave s�o desligados quando T1E � muito pequeno.
		
		if (CEIMFP_.T1E[I-1][*M-1] > 1.0e-6*CEEL00_.XE1[I-1]*COMPOS_.VMOL[*M-1]){
			CEIMFP_.T1E[I-1][*M-1]=log(fmax(CEIMFP_.T1E[I-1][*M-1],1.0e-35));
			CEIMFP_.T2E[I-1][*M-1]=log(fmax(CEIMFP_.T2E[I-1][*M-1],1.0e-35));
		} else{
			CEIMFP_.T1E[I-1][*M-1] = -100.0e0;
			CEIMFP_.T2E[I-1][*M-1] = -100.0e0;
		}	
	}
	

	
	//Espalhamento el�stico de p�sitrons.
	
	for (int I = 1; I <= NDATA; I++){
		CEEL00_.FJL[I-1]=log(CEEL00_.XP0[I-1]);
	}
	spline2_(CEEL00_.EJTL,CEEL00_.FJL,CEEL00_.A,CEEL00_.B,CEEL00_.C,CEEL00_.D,0.0e0,0.0e0,NDATA);
	for (int I = 1; I<= NEGP; I++){
		EC= CEGRID_.DLEMP[I-1];
		findi2_(CEEL00_.EJTL,&EC,&NDATA,&J);
		CEEL00_.XP0[I-1]=exp(CEEL00_.A[J-1]+EC*(CEEL00_.B[J-1]+EC*(CEEL00_.C[J-1]+EC*CEEL00_.D[J-1])));
	}
	
	for (int I = 1; I <= NDATA; I++){
		CEEL00_.FJL[I-1]=log(CEEL00_.XP1[I-1]);
	}
	spline2_(CEEL00_.EJTL,CEEL00_.FJL,CEEL00_.A,CEEL00_.B,CEEL00_.C,CEEL00_.D,0.0e0,0.0e0,NDATA);
	for (int I = 1; I<= NEGP; I++){
		EC= CEGRID_.DLEMP[I-1];
		findi2_(CEEL00_.EJTL,&EC,&NDATA,&J);
		CEEL00_.XP1[I-1]=exp(CEEL00_.A[J-1]+EC*(CEEL00_.B[J-1]+EC*(CEEL00_.C[J-1]+EC*CEEL00_.D[J-1])));
	}
	
	for (int I = 1; I <= NDATA; I++){
		CEEL00_.FJL[I-1]=log(CEEL00_.XP2[I-1]);
	}
	spline2_(CEEL00_.EJTL,CEEL00_.FJL,CEEL00_.A,CEEL00_.B,CEEL00_.C,CEEL00_.D,0.0e0,0.0e0,NDATA);
	for (int I = 1; I<= NEGP; I++){
		EC= CEGRID_.DLEMP[I-1];
		findi2_(CEEL00_.EJTL,&EC,&NDATA,&J);
		CEEL00_.XP2[I-1]=exp(CEEL00_.A[J-1]+EC*(CEEL00_.B[J-1]+EC*(CEEL00_.C[J-1]+EC*CEEL00_.D[J-1])));
	}
	
	for (int I = 1; I <= NEGP; I++){
		XS0=CEEL00_.XP0[I-1];
        XS1=CEEL00_.XP1[I-1];
        XS2=CEEL00_.XP2[I-1];
        FPEL=1.0e0/(XS0*COMPOS_.VMOL[*M-1]);
        FPT1=1.0e0/(XS1*COMPOS_.VMOL[*M-1]);
        FPST=CEGRID_.ET[I-1]/(CPIMFP_.CSTPP[I-1][*M-1]+CPIMFP_.RSTPP[I-1][*M-1]);
        XS0H=1.0e0/(COMPOS_.VMOL[*M-1] * fmax(FPEL, fmin(PENELOPE_mod_.C1[*M-1] * FPT1,PENELOPE_mod_.C2[*M-1]*FPST)));
        eela02_(XS0,XS1,XS2,XS0H,AAA,BBB,RNDC,XS1S,XS2S);
        CPIMFP_.SPHEL[I-1][*M-1]=XS0H*COMPOS_.VMOL[*M-1];
        CPIMFP_.RNDCP[I-1][*M-1]=RNDC;
        CPIMFP_.AP[I-1][*M-1]=AAA;
        CPIMFP_.BP[I-1][*M-1]=BBB;
        CEEL00_.T1P0[I-1]=XS1S;
        CPIMFP_.T1P[I-1][*M-1]=CEINTF_.T1PI[I-1]+XS1S*COMPOS_.VMOL[*M-1];
        CEEL00_.T2P0[I-1]=XS2S;
        CPIMFP_.T2P[I-1][*M-1]=CEINTF_.T2PI[I-1]+XS2S*COMPOS_.VMOL[*M-1];	
	}
	
	//Imprimir tabelas de dispers�o el�stica de p�sitron.
	
	if (*INFO >= 3){
		fprintf(IWR, "\n PENELOPE >>>  Elastic scattering of positrons\n");
		fprintf(IWR, "\n   E (eV)      MFP (mtu)   TMFP1 (mtu)  MFPh (mtu)        A           B           RNDC\n");
  	    fprintf(IWR, " -----------------------------------------------------------------------------------------\n");

	}
	
	for (int I = 1; I <= NEGP; I++){
		FP0=COMPOS_.RHO[*M-1]/(CEEL00_.XP0[I-1]*COMPOS_.VMOL[*M-1]);
        FP1=COMPOS_.RHO[*M-1]/(CEEL00_.XP1[I-1]*COMPOS_.VMOL[*M-1]);
        HMFP=COMPOS_.RHO[*M-1]/CPIMFP_.SPHEL[I-1][*M-1];
        if (*INFO >= 3){
			fprintf(IWR, "%.5E %.5E %.5E %.5E %.5E %.5E %.5E\n", 
			CEGRID_.ET[I-1],FP0,FP1,HMFP,CPIMFP_.AP[I-1][*M-1],CPIMFP_.BP[I-1][*M-1],CPIMFP_.RNDCP[I-1][*M-1]);
		}
		CPIMFP_.SPHEL[I-1][*M-1]=log(CPIMFP_.SPHEL[I-1][*M-1]);
        CPIMFP_.AP[I-1][*M-1]=log(CPIMFP_.AP[I-1][*M-1]);
        
        //Os eventos de dispers�o suave s�o desligados quando T1E � muito pequeno.
		
		if (CPIMFP_.T1P[I-1][*M-1] > 1.0e-6*CEEL00_.XP1[I-1]*COMPOS_.VMOL[*M-1]){
			CPIMFP_.T1P[I-1][*M-1]=log(fmax(CPIMFP_.T1P[I-1][*M-1],1.0e-35));
			CPIMFP_.T2P[I-1][*M-1]=log(fmax(CPIMFP_.T2P[I-1][*M-1],1.0e-35));
		} else{
			CPIMFP_.T1P[I-1][*M-1] = -100.0e0;
			CPIMFP_.T2P[I-1][*M-1] = -100.0e0;
		}	
	}
	
//	printf("\n\nEELAR2\n\n");
}


void eela02_(double &XS0, double &XS1, double &XS2, double &XS0H, double &A, double &B, double &RNDC, double &XS1S, double &XS2S){
	
	/*
	Esta sub-rotina determina os par�metros do modelo MW para
  espalhamento el�stico de el�trons e p�sitrons e inicializa o
  algoritmo de simula��o mista (para part�culas com uma determinada energia).

  Argumentos de entrada:
    XS0 .... se��o x total (cm ** 2).
    XS1 .... 1a se��o x de transporte (cm ** 2).
    XS2 .... 2o transporte se��o x (cm ** 2).
    XS0H ... se��o x sugerida para eventos dif�ceis (cm ** 2).

  Valores de sa�da:
    A, B ... par�metros de distribui��o angular.
    RNDC ... probabilidade de corte.
    XS0H ... se��o x adotada para eventos dif�ceis (cm ** 2).
    XS1S ... 1a se��o x de transporte para eventos leves (cm ** 2).
    XS2S ... 2a se��o x de transporte para eventos leves (cm ** 2).
	
	*/	
	
	double TST, AU, AL, RMU1, RMU2, RMU1W, RMU2W, A1, B1, RND0, RMUC, WB, C1, C2, D1, D2, D3, F, RNDMB, RNDCM;
	
	
	if (XS0 < 0.0e0){
		printf("EELa0. Negative total cross section.\n");
		exit(0);
	}
	RMU1=fmin(XS1/(2.0e0*XS0),0.48e0);  //! Ensures numerical consistency.
	RMU2=fmin((3.0e0*XS1-XS2)/(6.0e0*XS0),0.32e0);


	if ((RMU1 < 0.0e0) || (RMU1 < RMU2)){
		printf("*** The arguments in subroutine EELa0 are inconsistent. XS0 = %f, XS1 = %f\n", XS0,XS1);
		exit(0);	
	}
	
	//Par�metro de triagem Wentzel.
	A=1.0e0;
L10:;
	A=A+A;
    TST=A*(A+1.0e0)*log((1.0e0+A)/A)-A-RMU1;
    if (TST < 0.0e0)
		goto L10;
    AU=A;
    AL=0.0e0;
L1:;
	A=0.5e0*(AL+AU);
	TST=A*(A+1.0e0)*log((1.0e0+A)/A)-A-RMU1;
	if (TST > 0.0e0){
		AU = A;
	} else{
		AL = A;
	}
	
	if (fabs(TST) > 1.0e-15)
		goto L1;
	
	/*Em altas energias, quando erros de truncamento nas tabelas de entrada
	s�o significativos, usamos espalhamento delta.*/
	if ((RMU2-pow(RMU1,2) < 1.0e-12) || (A < 1.0e-9)){
		B=1.0e0;
       	RNDC=1.0e0-XS0H/XS0;
       	if (RNDC < 1.0e-14)
			RNDC=0.0e0;
       	XS1S=XS1*RNDC;
       	XS2S=XS2*RNDC;
       	return;
	}
	
	RMU1W=A*(A+1.0e0)*log((1.0e0+A)/A)-A;
    RMU2W=A*(1.0e0-2.0e0*RMU1W);
    B=(RMU2W-RMU2)/(RMU2W-RMU1W*RMU1W);
	
	// CASO I
	
	if(B > 0.0e0) {
        RNDC=1.0e0-XS0H/XS0;
        if (RNDC < 1.0e-6) {
			RNDC=0.0e0;
            XS0H=XS0;
            XS1S=0.0e0;
            XS2S=0.0e0;
            return;
		}
		A1=A+1.0e0;
        B1=1.0e0-B;
        RND0=B1*A1*RMU1/(A+RMU1);
        RNDC=1.0e0-XS0H/XS0;
        if (RNDC < RND0){
        	RMUC=RNDC*A/(B1*A1-RNDC);
            XS1S=B1*A*A1*(log((A+RMUC)/A)-(RMUC/(A+RMUC)));
            XS2S=B1*(A*A1*pow(RMUC,2)/(A+RMUC))-2.0e0*A*XS1S;
		} else if (RNDC > RND0+B){
			RNDMB=RNDC-B;
            RMUC=RNDMB*A/(B1*A1-RNDMB);
            XS1S=B1*A*A1*(log((A+RMUC)/A)-(RMUC/(A+RMUC)));
            XS2S=B1*(A*A1*pow(RMUC,2)/(A+RMUC))-2.0e0*A*XS1S;
            XS1S=XS1S+B*RMU1;
            XS2S=XS2S+B*pow(RMU1,2);
		} else{
			RMUC=RMU1;
            WB=RNDC-RND0;
            XS1S=B1*A*A1*(log((A+RMUC)/A)-(RMUC/(A+RMUC)));
            XS2S=B1*(A*A1*pow(RMUC,2)/(A+RMUC))-2.0e0*A*XS1S;
            XS1S=XS1S+WB*RMU1;
            XS2S=XS2S+WB*pow(RMU1,2);	
		}
		XS2S=6.0e0*XS0*(XS1S-XS2S);
        XS1S=2.0e0*XS0*XS1S;
        return;
	}
	if (B > -1.0e-12){
		B=0.0e0;
        RNDC=1.0e0-XS0H/XS0;
        A1=A+1.0e0;
        RMUC=RNDC*A/(A1-RNDC);
        XS1S=A*A1*(log((A+RMUC)/A)-(RMUC/(A+RMUC)));
        XS2S=(A*A1*pow(RMUC,2)/(A+RMUC))-2.0e0*A*XS1S;
        XS2S=6.0e0*XS0*(XS1S-XS2S);
        XS1S=2.0e0*XS0*XS1S;
        return;
	}
	
	//Caso II
	C1=8.333333333333333e-1;
    C2=7.083333333333333e-1;
    D1=C2-RMU2;
    D2=C1-RMU1;
    D3=C2*RMU1-C1*RMU2;
    AL=1.0e-24;
    AU=A;
L2:;
	A=0.5e0*(AL+AU);
    RMU1W=A*(A+1.0e0)*log((1.0e0+A)/A)-A;
    RMU2W=A*(1.0e0-2.0e0*RMU1W);
    F=D1*RMU1W-D2*RMU2W-D3;
    if (F < 0.0e0){
		AL = A;
	} else{
		AU = A;
	}
	
	if (AU-AL > 1.0e-14*A)
		goto L2;
	
	B=(RMU1W-RMU1)/(C1-RMU1W);
    RNDC=1.0e0-XS0H/XS0;
    
    if (RNDC < 1.0e-10){
    	RNDC=0.0e0;
        XS0H=XS0;
        XS1S=0.0e0;
        XS2S=0.0e0;
        return;
	}
	
	A1=A+1.0e0;
    B1=1.0e0+B;
    RNDCM=B1*A1*0.5e0/(A+0.5e0);
    
    if (RNDC > RNDCM){
    	RNDC=RNDCM;
        XS0H=XS0*(1.0e0-RNDC);
	}
	
	RMUC=RNDC*A/(B1*A1-RNDC);
    XS1S=B1*A*A1*(log((A+RMUC)/A)-(RMUC/(A+RMUC)));
    XS2S=B1*(A*A1*pow(RMUC,2)/(A+RMUC))-2.0e0*A*XS1S;
    XS2S=6.0e0*XS0*(XS1S-XS2S);
    XS1S=2.0e0*XS0*XS1S;
    
  //  printf("\n\nEELA01\n\n");
    return;
}  


void eeldr2_(int *M, FILE *IRD, FILE *IWR, int *INFO){
	/*
	Esta sub-rotina l� se��es transversais el�sticas para el�trons e
	p�sitrons no material M do banco de dados de espalhamento el�stico. Isso tamb�m
   inicializa o algoritmo para simula��o de colis�es el�sticas de
    El�trons e p�sitrons
	*/
	
	static const int NE = 96;
	static const int NA = 606;
	static const int NM = 512;
	static const int NP = 128;
	
	
	double PI = 3.1415926535897932e0;
//	double NE = 96;
//	double NA = 606;
//	double NM = 512;
//	double NP = 128;
	double EGRD[16];
	FILE *F;
	
	double  FGRID, EV,   ETSIE, ECS0, ECS1, ECS2, TCS1, TCS2, TS0, TS1, TS2, TSTE,
			FPEL, FPT1, FPST, XS0H, RNDC, RU,RR, DPRO, CI, RMUC, FP0, FP1, HMFP, ERRM,
			XM0, XM0A, XM1, XM2  ; 

	int IE, IGRID,II, IELEC, NPP, NU, IEME, IEMP, J, K;
		
	EGRD[0] = 1.0e0;
	EGRD[1] = 1.25e0;
	EGRD[2] = 1.50e0;
	EGRD[3] = 1.75e0;
	EGRD[4] = 2.00e0;
	EGRD[5] = 2.50e0;
	EGRD[6] = 3.00e0;
	EGRD[7] = 3.50e0;
	EGRD[8] = 4.00e0;
	EGRD[9] = 4.50e0;
	EGRD[10] = 5.00e0;
	EGRD[11] = 6.00e0;
	EGRD[12] = 7.00e0;
	EGRD[13] = 8.00e0;
	EGRD[14] = 9.00e0;
	EGRD[15] = 1.00e1;
	
	//Pontos da malha de energia (em eV).
	
	IE=0;
	IGRID=10;
	FGRID=10.0e0;
	
L10:;

	IGRID=IGRID+1;
	EV=EGRD[IGRID-1]*FGRID;
	
	if (IGRID == 16){
		IGRID=1;
        FGRID=10.0e0*FGRID;
	}
	
	IE=IE+1;
	CDCSEP_.ETS[IE-1]=EV;
	CDCSEP_.ETL[IE-1]=log(CDCSEP_.ETS[IE-1]);
	
	if (IE < NE)
		goto L10;
	
	//Grade angular (TH em graus, XMU = (1.0D0-COS (TH)) / 2).
	
	II=1;
    CDCSEP_.TH[II-1]=0.0e0;
    CDCSEP_.THR[II-1]=CDCSEP_.TH[II-1]*PI/180.0e0;
    CDCSEP_.XMU[II-1]=(1.0e0-cos(CDCSEP_.THR[II-1]))/2.0e0;
    CDCSEP_.XMUL[II-1]=log(1.0e-35);
    II=2;
    CDCSEP_.TH[II-1]=1.0e-4;
    CDCSEP_.THR[II-1]=CDCSEP_.TH[II-1]*PI/180.0e0;
    CDCSEP_.XMU[II-1]=(1.0e0-cos(CDCSEP_.THR[II-1]))/2.0e0;
    CDCSEP_.XMUL[II-1]=log(fmax(CDCSEP_.XMU[II-1],1.0e-35));
       
L20:;
	II=II+1;
	if (CDCSEP_.TH[II-1-1] < 0.9999e-3)
		CDCSEP_.TH[II-1]=CDCSEP_.TH[II-1-1]+2.5e-5;
	else if (CDCSEP_.TH[II-1-1] < 0.9999e-2)
		CDCSEP_.TH[II-1]=CDCSEP_.TH[II-1-1]+2.5e-4;
	else if (CDCSEP_.TH[II-1-1] < 0.9999e-1)
		CDCSEP_.TH[II-1]=CDCSEP_.TH[II-1-1]+2.5e-3;
	else if (CDCSEP_.TH[II-1-1] < 0.9999e+0)
		CDCSEP_.TH[II-1]=CDCSEP_.TH[II-1-1]+2.5e-2;
	else if (CDCSEP_.TH[II-1-1] < 0.9999e+1)
		CDCSEP_.TH[II-1]=CDCSEP_.TH[II-1-1]+1.0e-1;
	else if (CDCSEP_.TH[II-1-1] < 2.4999e+1)
		CDCSEP_.TH[II-1]=CDCSEP_.TH[II-1-1]+2.5e-1;
	else
		CDCSEP_.TH[II-1]=CDCSEP_.TH[II-1-1]+5.0e-1;
	
	CDCSEP_.THR[II-1]=CDCSEP_.TH[II-1]*PI/180.0e0;
	CDCSEP_.XMU[II-1]=fmax((1.0e0-cos(CDCSEP_.THR[II-1]))/2.0e0,1.0e-35);
	CDCSEP_.XMUL[II-1]=log(fmax(CDCSEP_.XMU[II-1],1.0e-35));
	
	if (II < NA)
		goto L20;
	
	//Leia tabelas el�sticas DCS.
	
	if (*INFO >= 3){
		fprintf(IWR, "\n *** Electron elastic differential cross sections");
	}
	
	IELEC = -1;
	
	
	
	fgets(LINHA, sizeof(LINHA), IRD);
	
	
	for (int IE = 1; IE <= NE; IE++){
		
		fgets(LINHA, sizeof(LINHA), IRD);
		extrairString(APOIO, LINHA, 0, 3);
		IELEC = atoi(APOIO);
		extrairString(APOIO, LINHA, 3, 10);
		ETSIE = atof(APOIO);
		extrairString(APOIO, LINHA, 13, 12);
		CDCSEP_.ECS[IE-1] = atof(APOIO);
		extrairString(APOIO, LINHA, 25, 12);
		CDCSEP_.ETCS1[IE-1] = atof(APOIO);
		extrairString(APOIO, LINHA, 37, 12);
		CDCSEP_.ETCS2[IE-1] = atof(APOIO);
		
		if (*INFO >= 3){
			fprintf(IWR, "\n%3d %.3E %.5E %.5E %.5E",IELEC, CDCSEP_.ETS[IE-1], CDCSEP_.ECS[IE-1], CDCSEP_.ETCS1[IE-1], CDCSEP_.ETCS2[IE-1] );
		}
		
		if ((IELEC != -1) || (fabs(ETSIE - CDCSEP_.ETS[IE-1]) > 0.1e0)){
			fprintf(IWR, "Error reading electron elastic DCS data.\n");
			printf("Error reading electron elastic DCS data.\n");
			exit(0);
		}
		
		int JJ = 0;
		int coluna = 0;
		for (int K = 1; K<= NA; K++){
			if (JJ == 0){
				fgets(LINHA, sizeof(LINHA), IRD);
				JJ = 10;	
				coluna = 0;
				if (*INFO >= 3)
				   	fprintf(IWR, "\n");	
			}
			if (JJ > 0){
				extrairString(APOIO, LINHA, coluna, 12);
				CDCSEP_.EDCS[K-1][IE-1] = atof(APOIO);
				coluna = coluna + 12;
				JJ--;
				if (*INFO >= 3)
					fprintf(IWR, "%.5E ", CDCSEP_.EDCS[K-1][IE-1]);
			}	
		}
		
		//Teste de Consistencia
		
		
		for (int K = 1; K <= NA; K++){
			CDCSEP_.DCSI[K-1]=CDCSEP_.EDCS[K-1][IE-1];
		}
		ECS0=4.0e0*PI*rmomx2_(CDCSEP_.XMU,CDCSEP_.DCSI,0.0e0,1.0e0,NA,0);
        ECS1=4.0e0*PI*rmomx2_(CDCSEP_.XMU,CDCSEP_.DCSI,0.0e0,1.0e0,NA,1);
        ECS2=4.0e0*PI*rmomx2_(CDCSEP_.XMU,CDCSEP_.DCSI,0.0e0,1.0e0,NA,2);
        TCS1=2.0e0*ECS1;
        TCS2=6.0e0*(ECS1-ECS2);
        TS0=(ECS0-CDCSEP_.ECS[IE-1])/CDCSEP_.ECS[IE-1];
        TS1=(TCS1-CDCSEP_.ETCS1[IE-1])/CDCSEP_.ETCS1[IE-1];
        TS2=(TCS2-CDCSEP_.ETCS2[IE-1])/CDCSEP_.ETCS2[IE-1];
        TSTE=fmax(fabs(TS0), fmax(fabs(TS1),fabs(TS2)));
        
        if (TSTE > 1.0e-4){
			fprintf(IWR, " E= %.5E\n", CDCSEP_.ETS[IE-1] );
			fprintf(IWR, "Electron cross section data are corrupt\n");
			printf("Electron cross section data are corrupt\n");
			exit(0);
		}	
	}
	

	IELEC = +1;
	
	if (*INFO >= 3){
		fprintf(IWR, "\n\n*** Positron elastic differential cross sections");
	}

	fgets(LINHA, sizeof(LINHA), IRD);
	for (int IE = 1; IE <= NE; IE++){
		fgets(LINHA, sizeof(LINHA), IRD);
		extrairString(APOIO, LINHA, 0, 3);
		IELEC = atoi(APOIO);
		extrairString(APOIO, LINHA, 3, 10);
		ETSIE = atof(APOIO);
		extrairString(APOIO, LINHA, 13, 12);
		CDCSEP_.PCS[IE-1] = atof(APOIO);
		extrairString(APOIO, LINHA, 25, 12);
		CDCSEP_.PTCS1[IE-1] = atof(APOIO);
		extrairString(APOIO, LINHA, 37, 12);
		CDCSEP_.PTCS2[IE-1] = atof(APOIO);
		
		if (*INFO >= 3){
			fprintf(IWR, "\n%3d %.3E %.5E %.5E %.5E",IELEC, CDCSEP_.ETS[IE-1], CDCSEP_.PCS[IE-1], CDCSEP_.PTCS1[IE-1], CDCSEP_.PTCS2[IE-1] );
		}
		
		if ((IELEC != +1) || (fabs(ETSIE - CDCSEP_.ETS[IE-1]) > 0.1e0)){
			fprintf(IWR, "Error reading positron elastic DCS data.\n");
			printf("Error reading positron elastic DCS data.\n");
			exit(0);
		}
		
		int J = 0;
		int coluna = 0;
		for (int K = 1; K<= NA; K++){
			if (J == 0){
				fgets(LINHA, sizeof(LINHA), IRD);
				J = 10;	
				coluna = 0;
				if (*INFO >= 3)
				   	fprintf(IWR, "\n");	
			}
			if (J > 0){
				extrairString(APOIO, LINHA, coluna, 12);
				CDCSEP_.PDCS[K-1][IE-1] = atof(APOIO);
				coluna = coluna + 12;
				J--;
				if (*INFO >= 3)
					fprintf(IWR, " %.5E", CDCSEP_.PDCS[K-1][IE-1]);
			}	
		}
		
		//Teste de Consistencia
		
		for (int K = 1; K <= NA; K++){
			CDCSEP_.DCSI[K-1]=CDCSEP_.PDCS[K-1][IE-1];
		}
		ECS0=4.0e0*PI*rmomx2_(CDCSEP_.XMU,CDCSEP_.DCSI,0.0e0,1.0e0,NA,0);
        ECS1=4.0e0*PI*rmomx2_(CDCSEP_.XMU,CDCSEP_.DCSI,0.0e0,1.0e0,NA,1);
        ECS2=4.0e0*PI*rmomx2_(CDCSEP_.XMU,CDCSEP_.DCSI,0.0e0,1.0e0,NA,2);
        TCS1=2.0e0*ECS1;
        TCS2=6.0e0*(ECS1-ECS2);
        TS0=(ECS0-CDCSEP_.PCS[IE-1])/CDCSEP_.PCS[IE-1];
        TS1=(TCS1-CDCSEP_.PTCS1[IE-1])/CDCSEP_.PTCS1[IE-1];
        TS2=(TCS2-CDCSEP_.PTCS2[IE-1])/CDCSEP_.PTCS2[IE-1];
        TSTE=fmax(fabs(TS0), fmax(fabs(TS1),fabs(TS2)));
        
        if (TSTE > 1.0e-4){
			fprintf(IWR, " E= %12.5E\n", CDCSEP_.ETS[IE-1] );
			fprintf(IWR, "Positron cross section data are corrupt\n");
			printf("Positron cross section data are corrupt\n");
			exit(0);
		}	
	}
	

	
	NPP=NP;
	NU=-NPP/4;
	
	//Eletrons
	double W0 = 0.0e0;
	double W1 = 1.0e0; 
	IEME=0;
	for (int KE = 1; KE <= NEGP; KE++){
		if (CEGRID_.ET[KE-1] > 0.999999e8)
			goto L100;
			
		dcsel02_(CEGRID_.ET[KE-1], -1);
		int PDF = 2;
		
		ritai02_(PDF,W0,W1,NPP,NU,ERRM,F);
		
		for (int I = 1; I <= NP; I++){
			CEELDB_.XSE[*M-1][KE-1][I-1]=CRITA_.XTI[I-1];
            CEELDB_.PSE[*M-1][KE-1][I-1]=CRITA_.PACI[I-1];
            CEELDB_.ASE[*M-1][KE-1][I-1]=CRITA_.AI[I-1];
            CEELDB_.BSE[*M-1][KE-1][I-1]=CRITA_.BI[I-1];
            CEELDB_.ITLE[*M-1][KE-1][I-1]=CRITA_.ITLI[I-1];
            CEELDB_.ITUE[*M-1][KE-1][I-1]=CRITA_.ITUI[I-1];
		}
		
		ritam2_(0.0e0,1.0e0,XM0A,XM1,XM2);

		
		ECS0=*CDCSEP_.CSI;
        ECS1=*CDCSEP_.CSI*XM1/XM0A;
        ECS2=*CDCSEP_.CSI*XM2/XM0A;
        CEEL00_.XE0[KE-1]=ECS0;
        CEEL00_.XE1[KE-1]=2.0e0*ECS1;
        CEEL00_.XE2[KE-1]=6.0e0*(ECS1-ECS2);
        
        FPEL=1.0e0/(CEEL00_.XE0[KE-1]*COMPOS_.VMOL[*M-1]);
        FPT1=1.0e0/(CEEL00_.XE1[KE-1]*COMPOS_.VMOL[*M-1]);
        FPST=CEGRID_.ET[KE-1]/(CEIMFP_.CSTPE[KE-1][*M-1]+CEIMFP_.RSTPE[KE-1][*M-1]);
        XS0H=1.0e0/(COMPOS_.VMOL[*M-1]*fmax(FPEL,fmin(PENELOPE_mod_.C1[*M-1]*FPT1,PENELOPE_mod_.C2[*M-1]*FPST)));
        RNDC=fmax(1.0e0-XS0H/CEEL00_.XE0[KE-1],1.0e-10);
        
        if (RNDC < 1.0e-6)
        	RNDC=0.0e0;
        
        CELSEP_.RNDCED[KE-1][*M-1]=RNDC;
        
        RU=RNDC;
        II=1;
        J=NP;
L1:;      
        K=(II+J)/2;
        if (RU > CEELDB_.PSE[*M-1][KE-1][K-1])
        	II = K;
        else
        	J=K;
        
        if (J-II > 1)
        	goto L1;
        
        RR=RU-CEELDB_.PSE[*M-1][KE-1][II-1];
        DPRO=CEELDB_.PSE[*M-1][KE-1][II+1-1]-CEELDB_.PSE[*M-1][KE-1][II-1];
        
        if (DPRO < 1.0e-10)
        	RMUC=CEELDB_.XSE[*M-1][KE-1][II-1];
        else{
        	CI=(1.0e0+CEELDB_.ASE[*M-1][KE-1][II-1]+CEELDB_.BSE[*M-1][KE-1][II-1])*DPRO;
 	   	    RMUC=CEELDB_.XSE[*M-1][KE-1][II-1]+(CI*RR/(pow(DPRO,2)+(DPRO*CEELDB_.ASE[*M-1][KE-1][II-1]
               +CEELDB_.BSE[*M-1][KE-1][II-1]*RR)*RR))*(CEELDB_.XSE[*M-1][KE-1][II+1-1]-CEELDB_.XSE[*M-1][KE-1][II-1]);
		}
		
		/*
		Momentos do PDF no intervalo restrito (0, RMUC).
		Total e se��es transversais de transporte para intera��es suaves.
		*/
		
		ritam2_(0.0e0,RMUC,XM0,XM1,XM2);
		ECS1=*CDCSEP_.CSI*XM1/XM0A;
        ECS2=*CDCSEP_.CSI*XM2/XM0A;
        TCS1=2.0e0*ECS1;
        TCS2=6.0e0*(ECS1-ECS2);
        CEIMFP_.SEHEL[KE-1][*M-1]=XS0H*COMPOS_.VMOL[*M-1];
        CEEL00_.T1E0[KE-1]=TCS1;
        CEIMFP_.T1E[KE-1][*M-1]=CEINTF_.T1EI[KE-1]+TCS1*COMPOS_.VMOL[*M-1];
        CEEL00_.T2E0[KE-1]=TCS2;
        CEIMFP_.T2E[KE-1][*M-1]=CEINTF_.T2EI[KE-1]+TCS2*COMPOS_.VMOL[*M-1];
        IEME=KE;
	}
		
L100:;
	CELSEP_.EELMAX[*M-1]=fmin(CEGRID_.ET[IEME-1]-1.0e0,0.999999e8);
	
	//Imprima tabelas de espalhamento el�stico de el�trons.
	
	if (*INFO >= 3){
		fprintf(IWR, "\n\n PENELOPE >>>  Elastic scattering of electrons  (ELSEPA database)\n");
		fprintf(IWR, "\n   E (eV)       MFP (mtu)      TMFP1 (mtu)      MFPh (mtu)\n");
		fprintf(IWR, " ------------------------------------------------------------\n");
	}
	
	for (int I = 1; I <= IEME; I++){
		FP0=COMPOS_.RHO[*M-1]/(CEEL00_.XE0[I-1]*COMPOS_.VMOL[*M-1]);
        FP1=COMPOS_.RHO[*M-1]/(CEEL00_.XE1[I-1]*COMPOS_.VMOL[*M-1]);
        HMFP=COMPOS_.RHO[*M-1]/CEIMFP_.SEHEL[I-1][*M-1];
        
        if (*INFO >= 3){
			fprintf(IWR, " %.5E  %.5E  %.5E  %.5E\n", CEGRID_.ET[I-1],FP0,FP1,HMFP);
		}
		CEIMFP_.SEHEL[I-1][*M-1]=log(CEIMFP_.SEHEL[I-1][*M-1]);
		
		//Os eventos de dispers�o suave s�o desligados quando T1E � muito pequeno.
		
		if (CEIMFP_.T1E[I-1][*M-1] > 1.0e-6*CEEL00_.XE1[I-1]*COMPOS_.VMOL[*M-1]){
			CEIMFP_.T1E[I-1][*M-1]=log(fmax(CEIMFP_.T1E[I-1][*M-1],1.0e-35));
 	        CEIMFP_.T2E[I-1][*M-1]=log(fmax(CEIMFP_.T2E[I-1][*M-1],1.0e-35));	
		}else{
			CEIMFP_.T1E[I-1][*M-1]=-100.0e0;
 	        CEIMFP_.T2E[I-1][*M-1]=-100.0e0;
		}
	}
	
	//Positrons
	
	IEMP=0;
	for (int KE = 1; KE <= NEGP; KE++){
		if (CEGRID_.ET[KE-1] > 0.999999e8)
			goto L200;
		
		dcsel02_(CEGRID_.ET[KE-1], +1);
		int PDF = 2;
		ritai02_(PDF,W0,W1,NPP,NU,ERRM,F);
		
		for (int I = 1; I <= NP; I++){
			CPELDB_.XSP[*M-1][KE-1][I-1]=CRITA_.XTI[I-1];
            CPELDB_.PSP[*M-1][KE-1][I-1]=CRITA_.PACI[I-1];
            CPELDB_.ASP[*M-1][KE-1][I-1]=CRITA_.AI[I-1];
            CPELDB_.BSP[*M-1][KE-1][I-1]=CRITA_.BI[I-1];
            CPELDB_.ITLP[*M-1][KE-1][I-1]=CRITA_.ITLI[I-1];
            CPELDB_.ITUP[*M-1][KE-1][I-1]=CRITA_.ITUI[I-1];
		}
		
		ritam2_(0.0e0,1.0e0,XM0A,XM1,XM2);
		
		ECS0=*CDCSEP_.CSI;
        ECS1=*CDCSEP_.CSI*XM1/XM0A;
        ECS2=*CDCSEP_.CSI*XM2/XM0A;
        CEEL00_.XP0[KE-1]=ECS0;
        CEEL00_.XP1[KE-1]=2.0e0*ECS1;
        CEEL00_.XP2[KE-1]=6.0e0*(ECS1-ECS2);
        
        FPEL=1.0e0/(CEEL00_.XP0[KE-1]*COMPOS_.VMOL[*M-1]);
        FPT1=1.0e0/(CEEL00_.XP1[KE-1]*COMPOS_.VMOL[*M-1]);
        FPST=CEGRID_.ET[KE-1]/(CPIMFP_.CSTPP[KE-1][*M-1]+CPIMFP_.RSTPP[KE-1][*M-1]);
        XS0H=1.0e0/(COMPOS_.VMOL[*M-1]*fmax(FPEL,fmin(PENELOPE_mod_.C1[*M-1]*FPT1,PENELOPE_mod_.C2[*M-1]*FPST)));
        RNDC=fmax(1.0e0-XS0H/CEEL00_.XP0[KE-1],1.0e-10);
        
        if (RNDC < 1.0e-6)
        	RNDC=0.0e0;
        
        CELSEP_.RNDCPD[KE-1][*M-1]=RNDC;
        
        RU=RNDC;
        II=1;
        J=NP;
L2:;      
        K=(II+J)/2;
        if (RU > CPELDB_.PSP[*M-1][KE-1][K-1])
        	II = K;
        else
        	J=K;
        
        if (J-II > 1)
        	goto L2;
        
        RR=RU-CPELDB_.PSP[*M-1][KE-1][II-1];
        DPRO=CPELDB_.PSP[*M-1][KE-1][II+1-1]-CPELDB_.PSP[*M-1][KE-1][II-1];
        
        if (DPRO < 1.0e-10)
        	RMUC=CPELDB_.XSP[*M-1][KE-1][II-1];
        else{
        	CI=(1.0e0+CPELDB_.ASP[*M-1][KE-1][II-1]+CPELDB_.BSP[*M-1][KE-1][II-1])*DPRO;
 	   	    RMUC=CPELDB_.XSP[*M-1][KE-1][II-1]+(CI*RR/(pow(DPRO,2)+(DPRO*CPELDB_.ASP[*M-1][KE-1][II-1]
               +CPELDB_.BSP[*M-1][KE-1][II-1]*RR)*RR))*(CPELDB_.XSP[*M-1][KE-1][II+1-1]-CPELDB_.XSP[*M-1][KE-1][II-1]);
		}
		
		/*
		Momentos do PDF no intervalo restrito (0, RMUC).
		Total e se��es transversais de transporte para intera��es suaves.
		*/
		
		ritam2_(0.0e0,RMUC,XM0,XM1,XM2);
		ECS1=*CDCSEP_.CSI*XM1/XM0A;
        ECS2=*CDCSEP_.CSI*XM2/XM0A;
        TCS1=2.0e0*ECS1;
        TCS2=6.0e0*(ECS1-ECS2);
        CPIMFP_.SPHEL[KE-1][*M-1]=XS0H*COMPOS_.VMOL[*M-1];
        CEEL00_.T1P0[KE-1]=TCS1;
        CPIMFP_.T1P[KE-1][*M-1]=CEINTF_.T1PI[KE-1]+TCS1*COMPOS_.VMOL[*M-1];
        CEEL00_.T2P0[KE-1]=TCS2;
        CPIMFP_.T2P[KE-1][*M-1]=CEINTF_.T2PI[KE-1]+TCS2*COMPOS_.VMOL[*M-1];
        IEMP=KE;
	}
L200:;
	CELSEP_.PELMAX[*M-1]=fmin(CEGRID_.ET[IEMP-1]-1.0e0,0.999999e8);
	
	//Imprima tabelas de espalhamento el�stico de el�trons.
	
	if (*INFO >= 3){
		fprintf(IWR, "\n PENELOPE >>>  Elastic scattering of positrons (ELSEPA database)\n");
		fprintf(IWR, "\n   E (eV)       MFP (mtu)     TMFP1 (mtu)     MFPh (mtu)\n");
		fprintf(IWR, " -------------------------------------------------------------------------\n");
	}
	
	for (int I = 1; I <= IEMP; I++){
		FP0=COMPOS_.RHO[*M-1]/(CEEL00_.XP0[I-1]*COMPOS_.VMOL[*M-1]);
        FP1=COMPOS_.RHO[*M-1]/(CEEL00_.XP1[I-1]*COMPOS_.VMOL[*M-1]);
        HMFP=COMPOS_.RHO[*M-1]/CPIMFP_.SPHEL[I-1][*M-1];
        
        if (*INFO >= 3){
			fprintf(IWR, " %.5E  %.5E  %.5E  %.5E\n", CEGRID_.ET[I-1],FP0,FP1,HMFP);
		}
		CPIMFP_.SPHEL[I-1][*M-1]=log(CPIMFP_.SPHEL[I-1][*M-1]);
		
		//Os eventos de dispers�o suave s�o desligados quando T1E � muito pequeno.
		
		if (CPIMFP_.T1P[I-1][*M-1] > 1.0e-6*CEEL00_.XP1[I-1]*COMPOS_.VMOL[*M-1]){
			CPIMFP_.T1P[I-1][*M-1]=log(fmax(CPIMFP_.T1P[I-1][*M-1],1.0e-35));
 	        CPIMFP_.T2P[I-1][*M-1]=log(fmax(CPIMFP_.T2P[I-1][*M-1],1.0e-35));	
		}else{
			CPIMFP_.T1P[I-1][*M-1]=-100.0e0;
 	        CPIMFP_.T2P[I-1][*M-1]=-100.0e0;
		}
	}
	
//	printf("\n\nEELDR2\n\n");
	
}


void dcsel02_(double E, int IELEC){
	
	/*
	Essas sub-rotinas l�em o banco de dados ELSEPA para espalhamento el�stico
 de el�trons e p�sitrons por �tomos neutros, e gerar a mol�cula
 lar DCS de um composto para valores arbitr�rios da energia (de 50 eV
 a 100 MeV) e a deflex�o angular RMU = (1-COS (THETA)) / 2.
	*/
	
	/*
	inicializa o c�lculo de DCSs
 para el�trons (IELEC = -1) ou p�sitrons (IELEC = + 1) com energia cin�tica
 E (eV). Ele constr�i uma tabela de valores DCS para a grade padr�o de
 deflex�es angulares RMU usando interpola��o de spline c�bica log-log em E
 das tabelas preparadas pela sub-rotina ELINIT.
	*/
	
	static const int NE = 96;
	static const int NA = 606;
	
	double Y[NE];
	double A[NE];
	double B[NE];
	double C[NE];
	double D[NE];
	
	double EL;
	int JE;
	
	int WNE = NE;
	
	if ((E < 49.999e0) || (E > 1.0001e8)){
		printf(" Error in DCSEL0: Energy out of range.\n");
		exit(0);
	}
	
	EL=log(E);
	findi2_(CDCSEP_.ETL,&EL,&WNE,&JE);
	
	for (int IA = 1; IA <= NA; IA++){
		for (int IE = 1; IE <= NE; IE++){
			if (IELEC == -1)
				Y[IE-1]=log(CDCSEP_.EDCS[IA-1][IE-1]);
			else
				Y[IE-1]=log(CDCSEP_.PDCS[IA-1][IE-1]);	
		}
		spline2_(CDCSEP_.ETL,Y,A,B,C,D,0.0e0,0.0e0,NE);
		
		
		CDCSEP_.DCSIL[IA-1]=A[JE-1]+EL*(B[JE-1]+EL*(C[JE-1]+EL*D[JE-1]));
        CDCSEP_.DCSI[IA-1]=exp(CDCSEP_.DCSIL[IA-1]);
	}
	
	for (int IE = 1; IE <= NE; IE++){
		if (IELEC == -1)
			Y[IE-1]=log(CDCSEP_.ECS[IE-1]);
		else
			Y[IE-1]=log(CDCSEP_.PCS[IE-1]);
	}
	
	spline2_(CDCSEP_.ETL,Y,A,B,C,D,0.0e0,0.0e0,NE);
	*CDCSEP_.CSI=exp(A[JE-1]+EL*(B[JE-1]+EL*(C[JE-1]+EL*D[JE-1])));
	
	for (int IE = 1; IE <= NE; IE++){
		if (IELEC == -1)
			Y[IE-1]=log(CDCSEP_.ETCS1[IE-1]);
		else
			Y[IE-1]=log(CDCSEP_.PTCS1[IE-1]);	
	}
	
	spline2_(CDCSEP_.ETL,Y,A,B,C,D,0.0e0,0.0e0,NE);
	*CDCSEP_.TCS1I=exp(A[JE-1]+EL*(B[JE-1]+EL*(C[JE-1]+EL*D[JE-1])));
	
	for (int IE = 1; IE <= NE; IE++){
		if (IELEC == -1)
			Y[IE-1]=log(CDCSEP_.ETCS2[IE-1]);
		else
			Y[IE-1]=log(CDCSEP_.PTCS2[IE-1]);
	}
	
	spline2_(CDCSEP_.ETL,Y,A,B,C,D,0.0e0,0.0e0,NE);
	*CDCSEP_.TCS2I=exp(A[JE-1]+EL*(B[JE-1]+EL*(C[JE-1]+EL*D[JE-1])));
	
//	printf("\n\nDCSE102\n\n");
}


double dcsel2_(double RMU){
	
	/*
	Esta fun��o calcula o DCS em (cm ** 2 / sr) por log-log linear entre-
	Pola��o C em RMU = (1-cos (teta)) / 2. � inicializado por sub-rotina
  DCSEL0, que deve ser chamado antes de usar DCSEL.
	*/
	
	static const int NE = 96;
	static const int NA = 606;
	
	double RMUL, resultado;
	int I;
	int WNA = NA;
	RMUL=log(fmin(fmax(RMU,1.0e-35),0.999999999999e0));
	
	findi2_(CDCSEP_.XMUL,&RMUL,&WNA,&I);
	resultado = exp(CDCSEP_.DCSIL[I-1]+(CDCSEP_.DCSIL[I+1-1]-CDCSEP_.DCSIL[I-1])*
			   ((RMUL-CDCSEP_.XMUL[I-1])/(CDCSEP_.XMUL[I+1-1]-CDCSEP_.XMUL[I-1])));

			   
  //  printf("\n\nDCSE12\n\n");
    return resultado;

}


void ritam2_(double XD, double XU, double &XM0, double &XM1, double &XM2){
	/*
	C�lculo de momentos (restritos) de um pdf, PDF (X), obtido a partir de
sua aproxima��o RITA.

 XD, XU ... limites do intervalo de integra��o.
 XM0 ...... momento de 0� ordem.
 XM1 ...... momento de 1� ordem.
 XM2 ...... momento de 2� ordem.
	
	*/
	
	
	static const int NM = 512;
	static const int NIP = 51;
	
	double XS[NIP];
	double YS[NIP];
	double DX, ETA, TAU, X1, X2, CON1, CON2, CONS, SUM;
	
	XM0=0.0e0;
	XM1=0.0e0;
	XM2=0.0e0;

	for (int I = 1; I <= *CRITA_.NPM1; I++){
		if (CRITA_.XT[I+1-1] >= XD && CRITA_.XT[I-1] <= XU){
			X1=fmax(CRITA_.XT[I-1],XD);
            X2=fmin(CRITA_.XT[I+1-1],XU);
            DX=(X2-X1)/(NIP-1);
            
            for (int K = 1; K <= NIP; K++){
            	XS[K-1]=X1+(K-1)*DX;
            	//Valor do pdf RITA racional no ponto XS (K).
            	TAU=(XS[K-1]-CRITA_.XT[I-1])/(CRITA_.XT[I+1-1]-CRITA_.XT[I-1]);
            	CON1=2.0e0*CRITA_.B[I-1] * TAU;
            	CON2=1.0e0+CRITA_.B[I-1]+CRITA_.A[I-1]*(1.0e0 - TAU);
            	if (fabs(CON1) > 1.0e-10*fabs(CON2))
            		ETA=CON2*(1.0e0-sqrt(1.0e0-2.0e0*TAU*CON1/pow(CON2,2)))/CON1;
            	else
            		ETA=TAU/CON2;
            	
            	YS[K-1]=CRITA_.DPAC[I-1]*pow((1.0e0+(CRITA_.A[I-1]+CRITA_.B[I-1]*ETA)*ETA),2) / 
					  ((1.0e0-CRITA_.B[I-1]*ETA*ETA)*(1.0e0+CRITA_.A[I-1]+CRITA_.B[I-1])*(CRITA_.XT[I+1-1]-CRITA_.XT[I-1]));	
			}
			
			//Integra��o de Simpson.
			CONS=DX*3.3333333333333333e-1;
 	        SUM=0.0e0;
 	        
 	        for (int L =3; L <= NIP; L=L+2){
 	        	SUM=SUM+YS[L-2-1]+4.0e0*YS[L-1-1]+YS[L-1];	 
			}
			XM0=XM0+SUM*CONS;
			for (int K = 1; K <= NIP; K++){
				YS[K-1]=YS[K-1]*XS[K-1];	
			}
			SUM=0.0e0;
			
			for (int L =3; L <= NIP; L=L+2){
 	        	SUM=SUM+YS[L-2-1]+4.0e0*YS[L-1-1]+YS[L-1];	 
			}
			XM1=XM1+SUM*CONS;
			for (int K = 1; K <= NIP; K++){
				YS[K-1]=YS[K-1]*XS[K-1];	
			}
			SUM=0.0e0;
			
			for (int L =3; L <= NIP; L=L+2){
 	        	SUM=SUM+YS[L-2-1]+4.0e0*YS[L-1-1]+YS[L-1];	 
			}
			XM2=XM2+SUM*CONS;
		}
	}
	
//	printf("\n\nRITAM2\n\n");	
	
}


void graar2_(int *M, FILE *IRD, FILE *IWR, int *INFO){

	/*
	Esta sub-rotina lê o fator de forma molecular quadrado e o DCS
    para espalhamento Rayleigh de fótons no material M. Essas duas funções
    são tabulados usando as mesmas grades para todos os materiais.
    A amostragem aleatória do ângulo de espalhamento é realizada usando o
    algoritmo RITA.

	*/

	
	double REV=5.10998928e5;
	double RREV = 1.0e0/REV;
	
	static const int NQ = 250;
	static const int NEX = 1024;
	static const int NP = 150;
	static const int NIP = 51;
	static const int NM = 512;
	
	double *Q = (double *) malloc(NQ * sizeof(double));
	double *F = (double *) malloc(NQ * sizeof(double));
	double *FFI = (double *) malloc(NQ * sizeof(double));
	double *ER = (double *) malloc(NEX * sizeof(double));
	double *A = (double *) malloc(NQ * sizeof(double));
	double *B = (double *) malloc(NQ * sizeof(double));
	double *C = (double *) malloc(NQ * sizeof(double));
	double *D = (double *) malloc(NQ * sizeof(double));
	double *QI = (double *) malloc(NIP * sizeof(double));
	double *FUN = (double *) malloc(NIP * sizeof(double));
	double *SUM = (double *) malloc(NIP * sizeof(double));
	
	//double Q[NQ], F[NQ], FFI[NQ], ER[NEX], A[NQ], B[NQ], C[NQ], D[NQ], QI[NIP], FUN[NIP], SUM[NIP];
	
	int NQQ;
	double XS1, Q2MIN, ERRM, DQ, Q1, Q2, QM, Q2M, TAU, CON1, CI, CON2, ETAP, EE, J1, XSMAX;
	int J, NPT, NU, NPI, II, IT;
	
	fgets(LINHA, sizeof(LINHA), IRD);
	extrairString(APOIO, LINHA, 32, 3);
	NQQ = atoi(APOIO);
	extrairString(APOIO, LINHA, 43, 4);
	*CGRA01_.NE = atoi(APOIO);
	
	if (*INFO >= 2){
		fprintf(IWR, "\n*** Rayleigh scattering.  NQ = %3d, NE= %4d\n", NQ, *CGRA01_.NE);
	}
	
	for (int I = 1; I <= NQ; I++){
		fgets(LINHA, sizeof(LINHA), IRD);
		extrairString(APOIO, LINHA, 0, 9);
		Q[I-1] = atof(APOIO);
		extrairString(APOIO, LINHA, 9, 12);
		FFI[I-1] = atof(APOIO);
		if (I == 1)
			CGRA02_.FF0[*M-1]=pow(FFI[I-1],2);
		F[I-1]=log(pow(FFI[I-1],2));
        CGRA01_.FF[I-1][*M-1]=FFI[I-1];
	}
	
	for (int I = 1; I <= *CGRA01_.NE; I++){
		fgets(LINHA, sizeof(LINHA), IRD);
		extrairString(APOIO, LINHA, 0, 12);
		ER[I-1] = atof(APOIO);
		extrairString(APOIO, LINHA, 12, 12);
		CGRA01_.XSRA[I-1][*M-1] = atof(APOIO);
		CGRA01_.XSRA[I-1][*M-1]=log(CGRA01_.XSRA[I-1][*M-1]*COMPOS_.VMOL[*M-1]);	
	}
	
	if (ER[*CGRA01_.NE-1] < *CEGRID_.EU){
		XS1=CGRA01_.XSRA[*CGRA01_.NE-1-1][*M-1]+(CGRA01_.XSRA[*CGRA01_.NE-1][*M-1]-CGRA01_.XSRA[*CGRA01_.NE-1-1][*M-1])
     	     *(log(*CEGRID_.EU)-log(ER[*CGRA01_.NE-1-1]))/(log(ER[*CGRA01_.NE-1])-log(ER[*CGRA01_.NE-1-1]));
        CGRA01_.XSRA[*CGRA01_.NE-1][*M-1]=XS1;
        ER[*CGRA01_.NE-1]=*CEGRID_.EU;	
	}
	
	if (*M == 1){
		*CGRA02_.QQM=pow(Q[NQ-1],2);
		for (int I = 1; I <= NQ; I++){
			CGRA02_.QQ[I-1]=log(pow(Q[I-1],2));	
		}
		
		for (int I= 1; I <= *CGRA01_.NE; I++){
			CGRA01_.ERA[I-1]=log(ER[I-1]);
		}
		
		for (int I = 1; I <= NEGP; I++){
			findi2_(CGRA01_.ERA,&CEGRID_.DLEMP[I-1],CGRA01_.NE,&J);
			CGRA01_.IED[I-1]=J;
		}
		
		for (int I = 1; I <= NEGP-1; I++){
			CGRA01_.IEU[I-1]=min(CGRA01_.IED[I+1-1]+1,*CGRA01_.NE);
		}
		CGRA01_.IEU[NEGP-1]=min(CGRA01_.IED[NEGP-1]+1,*CGRA01_.NE);
	}
	
	spline2_(CGRA02_.QQ,F,A,B,C,D,0.0e0,0.0e0,NQ);
	for (int I = 1; I <= NQ; I++){
    	CGRA02_.AR[I-1][*M-1]=A[I-1];
        CGRA02_.BR[I-1][*M-1]=B[I-1];
        CGRA02_.CR[I-1][*M-1]=C[I-1];
        CGRA02_.DR[I-1][*M-1]=D[I-1];
	}
	
	/*
	Se��o transversal total nos pontos da grade de simula��o (ligeiramente
     aumentado para simplificar a interpola��o).
	*/
	
	EE=CEGRID_.DLEMP[1-1];
	findi2_(CGRA01_.ERA,&EE,CGRA01_.NE,&J);
	XS1=CGRA01_.XSRA[J-1][*M-1]+(CGRA01_.XSRA[J+1-1][*M-1]-CGRA01_.XSRA[J-1][*M-1])*(EE-CGRA01_.ERA[J-1])/(CGRA01_.ERA[J+1-1]-CGRA01_.ERA[J-1]);
	for (int IE = 1; IE <= NEGP-1; IE++){
		XSMAX=XS1;
        J1=J+1;
        EE=CEGRID_.DLEMP[IE+1-1];
        findi2_(CGRA01_.ERA,&EE,CGRA01_.NE,&J);
  		XS1=CGRA01_.XSRA[J-1][*M-1]+(CGRA01_.XSRA[J+1-1][*M-1]-CGRA01_.XSRA[J-1][*M-1])*(EE-CGRA01_.ERA[J-1])/(CGRA01_.ERA[J+1-1]-CGRA01_.ERA[J-1]);
        XSMAX=fmax(XSMAX,XS1);	
		
		if (J1 < J){
			for (int I = J1; I <= J; I++){
				XSMAX=fmax(XSMAX,CGRA01_.XSRA[I-1][*M-1]);		
			}
		}
		CGIMFP_.SGRA[IE-1][*M-1]=exp(XSMAX);	
	}
	CGIMFP_.SGRA[NEGP-1][*M-1] = CGIMFP_.SGRA[NEGP-1-1][*M-1]; 
	
	if (*INFO >= 2){
		fprintf(IWR, "\n   Q/me*c     Form factor\n");
		fprintf(IWR, " ----------------------------------\n");
	
		for (int I = 1; I <= NQ; I++){
			fprintf(IWR, " %.5E  %.5E\n", Q[I-1],FFI[I-1]);
		}
		fprintf(IWR,"\n   Energy       CS-Rayl\n    (eV)        (cm**2)\n");
		fprintf(IWR, " ----------------------------------\n");
		for (int I = 1; I <= *CGRA01_.NE; I++){
			fprintf(IWR, " %.5E  %.5E\n", ER[I-1],exp(CGRA01_.XSRA[I-1][*M-1])/COMPOS_.VMOL[*M-1]);
		}
		
	}
	
		
	//Inicializa��o do algoritmo RITA para amostragem aleat�ria do
	//Transfer�ncia de momento quadrado do fator de forma molecular ao quadrado.
	
	*CGRA00_.MM=*M;
    Q2MIN=0.0e0;
    *CGRA00_.Q2MAX=0.0e0;
    NPT=NP;
    NU=NPT/4;
    
    for (int I = 2; I <= NQ; I++){
		if (graaf22_(pow(Q[I-1],2)) > 1.0e-35)
			*CGRA00_.Q2MAX=pow(Q[I-1-1],2);		
	}
	int PDF = 3;
	
	ritai02_(PDF,Q2MIN,*CGRA00_.Q2MAX,NPT,NU,ERRM,0);

	
	NPI = *CRITA_.NPM1I+1;
	if (NPI != NP){
		printf("GRAaR. RITA initialisation error.\n");
		printf("The number of fixed grid points is %d\n", NPI);
		printf("The required number of grid points was %d\n", NP);
		exit(0);
	}
	
	if (ERRM > 1.0e-5){
		printf("GRAaR. RITA interpolation error is too large.\n");
		printf("The interpolation error is %f\n", ERRM);
		exit(0);
	}
	
	//Limite superior do intervalo X2 para as energias da rede PENELOPE.
	
	for (int IE = 1; IE <= NEGP; IE++){
		QM=2.0e0*CEGRID_.ET[IE-1]*RREV;
        Q2M=QM*QM;
        if (Q2M > CRITA_.QTI[1-1]){
        	if (Q2M < CRITA_.QTI[NP-1]){
        		II=1;
            	J=NPI;
L1:;			
				IT=(II+J)/2;
				if (Q2M > CRITA_.QTI[IT-1])	
					II=IT;
				else
					J=IT;
				if (J-II > 1)
					goto L1;
				
				
				Q1=CRITA_.QTI[II-1];
            	Q2=Q2M;
            	DQ=(Q2-Q1)/(NIP-1);
            		
            	for (int K =1; K <= NIP; K++){
            		QI[K-1]=Q1+(K-1)*DQ;
              	    TAU=(QI[K-1]-CRITA_.QTI[II-1])/(CRITA_.QTI[II+1-1]-CRITA_.QTI[II-1]);
              	    CON1=2.0e0*CRITA_.BI[II-1]*TAU;
              	    CI=1.0e0+CRITA_.AI[II-1]+CRITA_.BI[II-1];
              	    CON2=CI-CRITA_.AI[II-1]*TAU;
              	    if (fabs(CON1) > 1.0e-16*fabs(CON2))
              	    	ETAP=CON2*(1.0e0-sqrt(1.0e0-2.0e0*TAU*CON1/pow(CON2,2)))/CON1;
              	    else
              	    	ETAP=TAU/CON2;
              	    
              	    FUN[K-1]=CRITA_.DPACI[II-1]*pow(1.0e0+(CRITA_.AI[II-1]+CRITA_.BI[II-1]*ETAP)*ETAP,2)
                           /((1.0e0-CRITA_.BI[II-1]*ETAP*ETAP)*CI*(CRITA_.QTI[II+1-1]-CRITA_.QTI[II-1]));	
				}
				
				slag62_(DQ,FUN,SUM,NIP);
				CGRA03_.PMAX[*M-1][IE-1]=CRITA_.PACI[II-1]+SUM[NIP-1];									
			} else{
				CGRA03_.PMAX[*M-1][IE-1]=1.0e0;
			}
		}else{
			CGRA03_.PMAX[*M-1][IE-1]=CRITA_.PACI[1-1];
		}	
	}
	for (int I = 1; I <= NP; I++){
    	CGRA03_.QRA[*M-1][I-1]=CRITA_.QTI[I-1];
        CGRA03_.PRA[*M-1][I-1]=CRITA_.PACI[I-1];
        CGRA03_.DPRA[*M-1][I-1]=CRITA_.DPACI[I-1];
        CGRA03_.ARA[*M-1][I-1]=CRITA_.AI[I-1];
        CGRA03_.BRA[*M-1][I-1]=CRITA_.BI[I-1];
        CGRA03_.ITLRA[*M-1][I-1]=CRITA_.ITLI[I-1];
        CGRA03_.ITURA[*M-1][I-1]=CRITA_.ITUI[I-1];
	}	
	
	
	free(Q);
	free(F);
	free(FFI);
	free(ER);
	free(A);
	free(B);
	free(C);
	free(D);
	free(QI);
	free(FUN);
	free(SUM);	
	
//	printf("\n\nGRAAR2\n\n");
	
}

double graaf22_(double Q2){
	
	static const int NQ = 250;
	//Fator de forma molecularao quadrado, como uma fun��o de (Q * SL / REV) ** 2.
	double resultado, QL, F2;
	int I;
	int NQW = NQ;
	
	
	
	if (Q2 < 1.0e-9)
		resultado = CGRA02_.FF0[*CGRA00_.MM-1];
	else 
		if (Q2 > *CGRA02_.QQM)
			resultado = 0.0e0;
		else{
			QL=log(Q2);
        	findi2_(CGRA02_.QQ,&QL,&NQW,&I);
        	F2=CGRA02_.AR[I-1][*CGRA00_.MM-1]+QL*(CGRA02_.BR[I-1][*CGRA00_.MM-1]
               +QL*(CGRA02_.CR[I-1][*CGRA00_.MM-1]+QL*CGRA02_.DR[I-1][*CGRA00_.MM-1]));
        	resultado=exp(F2);
		}

   return resultado;
   
 //  printf("\n\nGRAAF22\n\n");
	
		
	
	
}


void slag62_(double H, double *Y, double *S, int N){
	
	/*
	Integra��o de Lagrange de seis pontos por partes de uma tabela uniformemente tabulada
  fun��o.

  Argumentos de entrada:
     H ............ espa�amento da grade.
     Y (I) (1: N) ... array de valores de fun��o (abscissas ordenadas).
     N ............ n�mero de pontos de dados.

  Argumento de sa�da:
     S (I) (1: N) ... matriz de valores integrais definidos como
                S (I) = INTEGRAL (Y) de X (1) a X (I) = X (1) + (I-1) * H.
	*/
	
	double HR, Y1, Y2, Y3, Y4, Y5, Y6;
	
	if (N < 6){
		printf("SLAG6: too few data points.\n");
	    exit(0);
	}
	HR=H/1440.0e0;
    Y1=0.0e0;
    Y2=Y[1-1];
    Y3=Y[2-1];
    Y4=Y[3-1];
    S[1-1]=0.0e0;
    S[2-1]=HR*(475*Y2+1427*Y3-798*Y4+482*Y[4-1]-173*Y[5-1]+27*Y[6-1]);
    S[3-1]=S[2-1]+HR*(-27*Y2+637*Y3+1022*Y4-258*Y[4-1]+77*Y[5-1]-11*Y[6-1]);
	
	for (int I = 4; I <= N-2; I++){
		Y1=Y2;
        Y2=Y3;
        Y3=Y4;
        Y4=Y[I-1];
        S[I-1]=S[I-1-1]+HR*(11*(Y1+Y[I+2-1])-93*(Y2+Y[I+1-1])+802*(Y3+Y4));
	}
	
	Y5=Y[N-1-1];
    Y6=Y[N-1];
    S[N-1-1]=S[N-2-1]+HR*(-27*Y6+637*Y5+1022*Y4-258*Y3+77*Y2-11*Y1);
    S[N-1]=S[N-1-1]+HR*(475*Y6+1427*Y5-798*Y4+482*Y3-173*Y2+27*Y1);
    
 //   printf("\n\nSLAG62\n\n");
	
	

}

void gphar2_(int *M, FILE *IRD, FILE *IWR, int *INFO){
	
	/*
	Esta sub-rotina l� se��es transversais fotoel�tricas dos elementos em
  material M e prepara tabelas de simula��o.

  NOTA: A matriz SGPH (M, IE) define uma fun��o constante por partes que
  � maior do que a se��o transversal fotoel�trica real. SGPH (M, IE) �
  definido como o maior valor da se��o x fotoel�trica no
  intervalo de energia de ET (IE) a ET (IE + 1). O caminho livre m�dio do f�ton �
  amostrados usando esta se��o transversal "aumentada" e, para compensar
  para isso, o f�ton sobrevive (ou seja, n�o � absorvido) com um problema
  capacidade de modo que o coeficiente de atenua��o fotoel�trico "exato"
  � reproduzido. Este truque permite que a sub-rotina JUMP ignore o
  exist�ncia de bordas de absor��o na se��o x fotoel�trica e
  realizar o rastreamento de f�tons mais r�pido.
	*/
    const static int NTP = 12000;
    const static int NDIM = 12000;
	double XGPHR[17][NDIM];
	
/*	double X1[NDIM]; 
	double Y1[NDIM];
	double X2[NDIM];
	double Y2[NDIM];*/

	double EG1, F1, EG2, DX, F2, FM;
    int ISH[17], IZZ, NSHR, NDATA, IC, NSHA, LST, N2, IST, N1, I1, I2;
    
 	double *X1; 
	double *Y1;
	double *X2;
	double *Y2;
	
	X1 = (double *) malloc(12000 * sizeof(double));
	X2 = (double *) malloc(12000 * sizeof(double));
	Y1 = (double *) malloc(12000 * sizeof(double));
	Y2 = (double *) malloc(12000 * sizeof(double));
	
	if ((X1 == NULL) || (X2 == NULL) || (Y1 == NULL) || (Y2 == NULL)){
		printf("\n\nmemoria insuficiente\n\n");
		
	}
	
//	double *X1 = new double[NDIM];
//	double *Y1 = new double[NDIM];
//	double *X2 = new double[NDIM];
//	double *Y2 = new double[NDIM];
	
   // char CS5[17][];
    
  /*  CS5[0]='total';
    CS5[1]='CS-K ';
    CS5[2]='CS-L1';
    CS5[3]='CS-L2';
    CS5[4]='CS-L3';
    CS5[5]='CS-M1';
    CS5[6]='CS-M2';
    CS5[7]='CS-M3';
    CS5[8]='CS-M4';
    CS5[9]='CS-M5';
    CS5[10]='CS-N1';
    CS5[11]='CS-N2';
    CS5[12]='CS-N3';
    CS5[13]='CS-N4';
    CS5[14]='CS-N5';
    CS5[15]='CS-N6';
    CS5[16]='CS-N7';*/
    
    char CS5[17][6] = {{"total"},{"CS-K "},{"CS-L1"},{"CS-L2"},{"CS-L3"},{"CS-M1"},{"CS-M2"},
	   	   	   	   	   {"CS-M3"},{"CS-M4"},{"CS-M5"},{"CS-N1"},{"CS-N2"},{"CS-N3"},
                       {"CS-N4"},{"CS-N5"},{"CS-N6"},{"CS-N7"}};
    
    //Leia as tabelas da se��o x do elemento
    
    for (int IEL = 1; IEL <= COMPOS_.NELEM[*M-1]; IEL++){
    	fgets(LINHA, sizeof(LINHA), IRD);
		extrairString(APOIO, LINHA, 40, 3);
        IZZ = atoi(APOIO);
        extrairString(APOIO, LINHA, 54, 3);
        NSHR = atoi(APOIO);
        extrairString(APOIO, LINHA, 67, 4);
        NDATA = atoi(APOIO);
        if (*INFO >= 2){
			fprintf(IWR, "\n*** Photoelectric cross sections,  IZ =%3d,  NSHELL =%3d, NDATA =%4d\n", IZZ,NSHR,NDATA);
		}	
		if (IZZ != COMPOS_.IZ[IEL-1][*M-1]){
			printf("GPHaR. Corrupt material data file.\n");
			exit(0);
			}
		if (NDATA > NDIM){
			printf("GPHaR. Too many data points.\n");
			exit(0);
		}	
		if (NSHR > 16){
			printf("GPHaR. Too many shells.\n");
			exit(0);
		}
		fgets(LINHA, sizeof(LINHA), IRD);
		int coluna = 13;
		for (int IS = 1; IS <= NSHR+1; IS++){
			extrairString(APOIO, LINHA, coluna, 2);
        	ISH[IS-1] = atoi(APOIO);
        	coluna = coluna +12;
		}
		
		for (int IE = 1; IE <= NDATA; IE++){
			fgets(LINHA, sizeof(LINHA), IRD);
			extrairString(APOIO, LINHA, 0, 12);
        	CGPH01_.ER[IE-1] = atof(APOIO);
			coluna = 12;
			for (int IS = 1; IS <= NSHR+1; IS++){
				extrairString(APOIO, LINHA, coluna, 12);
        		XGPHR[IS-1][IE-1] = atof(APOIO);
        		coluna = coluna +12;
			}
	    }
	    
	    //Remova os inv�lucros com energias de ioniza��o inferiores a 50 eV.
	    
	    if (NSHR > 1){
	    	NSHA=NSHR;
	    	for (int IS = NSHA+1; IS >=2; IS--){
				if (CADATA_.EB[ISH[IS-1]-1][IZZ-1] < 50.0e0)
					NSHR=NSHR-1;
				else
					goto L1;
			}
			if (NSHR < 0)
				NSHR=0;
		}
		
L1:;

		if (*INFO >= 2){
			fprintf(IWR, "\n   Energy       ");
			for (int IS = 1; IS <= NSHR+1; IS++){
				fprintf(IWR, " %s       ", CS5[ISH[IS-1]+1-1]);
			}
			fprintf(IWR,"\n");
			
			for (int IE = 1; IE <= NDATA; IE++){
				fprintf(IWR, " %.5E ", CGPH01_.ER[IE-1]);
				for (int IS = 1; IS <= NSHR+1; IS++){
				   	fprintf(IWR, "%.5E ",XGPHR[IS-1][IE-1] );
				}
				fprintf(IWR,"\n");
			}	
		}

        if (CGPH00_.NPHS[IZZ-1] == 0){
        	CGPH00_.IPHF[IZZ-1]=*CGPH00_.NCUR+1;
        	if ((*CGPH00_.NCUR+NDATA) > NTP){
        		fprintf(IWR, "Insufficient memory storage in GPHaR.\n");
        		fprintf(IWR, "Increase the value of the parameter NTP to %d.\n", *CGPH00_.NCUR+NDATA);
        		printf("Insufficient memory storage in GPHaR.\n");
        		printf("Increase the value of the parameter NTP to %d.\n", *CGPH00_.NCUR+NDATA);
        		exit(0);	
			}
		
			for (int IE = 1; IE <= NDATA; IE++){
				IC=*CGPH00_.NCUR+IE;
            	CGPH00_.EPH[IC-1]=log(CGPH01_.ER[IE-1]);
            	for (int IS = 1; IS <= NSHR+1; IS++){
            		CGPH00_.XPH[ISH[IS-1]+1-1][IC-1]=log(fmax(XGPHR[ISH[IS-1]+1-1][IE-1],1.0e-35));	
				}
				
			}
       	    *CGPH00_.NCUR=*CGPH00_.NCUR+NDATA;
            CGPH00_.IPHL[IZZ-1]=*CGPH00_.NCUR;
            CGPH00_.NPHS[IZZ-1]=NSHR;			
		}	
	}
	
	//Coeficiente de atenua��o fotoel�trica total.
	
	IZZ=COMPOS_.IZ[1-1][*M-1];
    IST=CGPH00_.IPHF[IZZ-1];
    
    LST=CGPH00_.IPHL[IZZ-1];
    *CGPH01_.NPHD=0;
    
   // printf("\n\nIZZ= %d, IST= %d, LST = %d\n\n", IZZ, IST, LST);
   
    for (int I = IST; I <= LST; I++){
    	*CGPH01_.NPHD=*CGPH01_.NPHD+1;
        CGPH01_.ER[*CGPH01_.NPHD-1]=exp(CGPH00_.EPH[I-1]);
        CGPH01_.XSR[*CGPH01_.NPHD-1]=COMPOS_.STF[1-1][*M-1]*exp(CGPH00_.XPH[1-1][I-1]);
	}
 
	if (COMPOS_.NELEM[*M-1] > 1){
		for (int IEL = 2; IEL <= COMPOS_.NELEM[*M-1]; IEL++){
			N1=*CGPH01_.NPHD;
		//	printf("\n\nN1= %d\n\n", N1);
			for (int I = 1; I <= N1; I++){
				X1[I-1]=CGPH01_.ER[I-1];
            	Y1[I-1]=CGPH01_.XSR[I-1];		
			}
  	   		IZZ=COMPOS_.IZ[IEL-1][*M-1];
    		IST=CGPH00_.IPHF[IZZ-1];
  	  	    LST=CGPH00_.IPHL[IZZ-1];
            N2=0;
            for (int I = IST; I <= LST; I++){
       	   	    N2=N2+1;
            	X2[N2-1]=exp(CGPH00_.EPH[I-1]);
            	Y2[N2-1]=COMPOS_.STF[IEL-1][*M-1]*exp(CGPH00_.XPH[1-1][I-1]);
			}
	
			merge22_(X1, Y1, X2, Y2, CGPH01_.ER, CGPH01_.XSR, N1, N2, *CGPH01_.NPHD);
			
		}
	}
	
    //Se��o transversal fotoel�trica total nos pontos da grade de simula��o (ligeiramente aumentado para simplificar a interpola��o).
    
    for (int I = 1; I <= *CGPH01_.NPHD; I++){
    	X1[I-1]=log(CGPH01_.ER[I-1]);
        Y1[I-1]=log(CGPH01_.XSR[I-1]*COMPOS_.VMOL[*M-1]);	
	}
	
	for (int IE = 1; IE <= NEGP-1; IE++){
		EG1=CEGRID_.DLEMP[IE-1];
        findi2_(X1,&EG1,CGPH01_.NPHD,&I1);
        if (I1 == *CGPH01_.NPHD) 
			I1= *CGPH01_.NPHD-1;
        DX=X1[I1+1-1]-X1[I1-1];
        if (DX > 1.0e-15)
        	F1=Y1[I1-1]+(Y1[I1+1-1]-Y1[I1-1])*(EG1-X1[I1-1])/DX;
        else
        	F1=Y1[I1-1];
        
        EG2=CEGRID_.DLEMP[IE+1-1];
        findi2_(X1,&EG2,CGPH01_.NPHD,&I2);
        if(I2 == *CGPH01_.NPHD)
        	I2= *CGPH01_.NPHD-1;
        DX=X1[I2+1-1]-X1[I2-1];
        if (DX > 1.0e-15)
        	F2=Y1[I2-1]+(Y1[I2+1-1]-Y1[I2-1])*(EG2-X1[I2-1])/DX;
        else
        	F2=Y1[I2-1];
        
        /*
		Para evitar a interpola��o das tabelas de coeficientes de atenua��o, n�s
		substituir o caminho livre m�dio inverso fotoel�trico (imfp) em cada
		Intervalo de energia  por seu limite superior. O aumento do imfp
   	    � interpretado como o imfp das intera��es delta
		*/
		
		FM=fmax(F1,F2);
		if (I1+1 <= I2){
			for (int I = I1+1; I <= I2; I++){
				FM=fmax(FM,Y1[I-1]);	
			}
		}
		
		CGIMFP_.SGPH[IE-1][*M-1]=exp(FM);	
	}
	
	CGIMFP_.SGPH[NEGP-1][*M-1]=CGIMFP_.SGPH[NEGP-1-1][*M-1];
	
	free(X1);
	free(X2);
	free(Y1);
	free(Y2);
	
//	printf("\n\nGPHAR2\n\n");
	
}

void merge22_(double *X1,double *Y1, double *X2, double *Y2, double *XM, double *YM, int &N1, int &N2, int &N){
	
	/*
	Esta sub-rotina mescla duas tabelas (X1,Y1), (X2,Y2) de duas fun��es
  para produzir uma tabela (XM,YM) da soma dessas fun��es, com abs-
  cissas em ordem crescente. As abcissas e os valores da fun��o s�o
  assumido como positivo. N1, N2 e N s�o os n�meros de pontos de grade
  nas tabelas de entrada e mescladas. Uma descontinuidade na fun��o �
  descrito dando duas vezes a abcissa. Interpola��o linear log-log
  � usado para interpolar as tabelas de entrada.
	
	*/
	
	const static int NP = 12000;
	double EPS = 1.0e-10;
	const static int NP2 = NP+NP;
	double *X = (double *) malloc(NP2 * sizeof(double));
	double *Y12 = (double *) malloc(NP2 * sizeof(double));
	double XC, TST1, TST2, TST12, TST3, TST4, TST34, TST, XMIN, YI1, YI2, SAVE ;
	int  IMIN, I1, J;
	
//	X1 = (double *) malloc(12000 * sizeof(double));
//	X2 =(double *) malloc(12000 * sizeof(double));
//	Y1 = (double *) malloc(12000 * sizeof(double));
//	Y2 =(double *) malloc(12000 * sizeof(double));
	
	
	
	if ((N1 > NP) || (N2 > NP)){
		printf("MERGE2. Increase the value of the parameter NP = %f\n", fmax(N1, N2));
		exit(0);
	}

	sort22_(X1,Y1,N1);
    sort22_(X2,Y2,N2);
    	
    for (int I1 = 1; I1 <= N1; I1++){
    	X[I1-1]=X1[I1-1];
	}
	
	N=N1;
	
	for (int I2 = 1; I2 <= N2; I2++){
		XC=X2[I2-1];
        findi2_(X1,&XC,&N1,&I1);
        
        if (I1 == N1)
        	I1=N1-1;
        
        TST1=fabs(XC-X1[I1-1]);
        TST2=fabs(XC-X1[I1+1-1]);
        TST12=fmin(TST1,TST2);
        if (I2 > 1)
        	TST3=fabs(XC-X2[I2-1-1]);
        else
        	TST3=1.0e0;
        if (I2 < N2)
        	TST4=fabs(XC-X2[I2+1-1]);
        else
        	TST4=1.0e0;
        
        TST34=fmin(TST3,TST4);
        TST=EPS*XC;
        if (TST34 > TST){
			if (TST12 > TST){
				N=N+1;
            	X[N-1]=XC;				
			}
		}else{
			N=N+1;
            X[N-1]=XC;
			
		}	
	}
	
	//Classifique e limpe a grade mesclada.
L1:;
    for (int I = 1; I <= N-1; I++){
    	IMIN=I;
        XMIN=X[I-1];
        for (int J = I+1; J <= N; J++){
			if (X[J-1] < XMIN){
				IMIN=J;
            	XMIN=X[J-1];
			}
		}
		SAVE=X[I-1];
        X[I-1]=X[IMIN-1];
        X[IMIN-1]=SAVE;	
	}	
	
	for (int I = 1; I <= N-2; I ++){
		if (X[I-1] > X[I+2-1]*(1.0e0-EPS)){
			X[I+1-1]=X[N-1];
            N=N-1;
			goto L1;
		}
	}
	
	for (int I =1; I <= N; I++){
		XC=X[I-1];
		if (I < N){
			if (X[I-1] > X[I+1-1]*(1.0e0-EPS))
				XC=X[I-1]*(1.0e0-EPS);
		}
		if (I > 1){
			if (X[I-1] < X[I-1-1]*(1.0e0+EPS))
				XC=X[I-1]*(1.0e0+EPS);
		}
		
		findi2_(X1,&XC,&N1,&J);
		if (J == N1)
			J = N1-1;
		if (X1[J+1-1] > X1[J-1]+EPS)
		    YI1=exp(log(Y1[J-1])+log(XC/X1[J-1])*log(Y1[J+1-1]/Y1[J-1])/log(X1[J+1-1]/X1[J-1]));
		else
			YI1=Y1[J-1];
		
		
		findi2_(X2,&XC,&N2,&J);
		if (J == N2)
			J = N2-1;
		if (X2[J+1-1] > X2[J-1]+EPS)
		    YI2=exp(log(Y2[J-1])+log(XC/X2[J-1])*log(Y2[J+1-1]/Y2[J-1])/log(X2[J+1-1]/X2[J-1]));
		else
			YI2=Y2[J-1];
		
		Y12[I-1]=YI1+YI2;
		if (Y12[I-1] < 1.0e-75)	
			Y12[I-1]=1.0e-75;		
	}
	
	if (N > NP){
		printf("MERGE2. Increase the value of the parameter NP. %d\n", N);
		exit(0);
	}
	
	for (int I = 1; I <= N; I++){
		XM[I-1]=X[I-1];
        YM[I-1]=Y12[I-1];
	}
	free(X);
	free(Y12);
	
//	printf("\n\nMERGE22\n\n");
	
}

void sort22_(double *X, double *Y, int N){
	
	/*
	Esta sub-rotina classifica uma tabela (X, Y) de uma fun��o com n pontos de dados.
 Uma descontinuidade da fun��o � descrita dando duas vezes o abscissa. Assume-se que a 
 fun��o � estritamente positiva (negativa Os valores C de Y s�o definidos como zero).
	*/
	const static int NP = 12000;
	
//	double IORDER[NP];
	double *IORDER = (double *) malloc(NP * sizeof(double));
	
	double  XMIN, ISAVE, SAVE;
	int I2, IMIN;
	
	if (N > NP){
		printf("SORT2. Increase the value of the parameter NP. %d\n", N);
		exit(0);
	}
	
	if (N == 1)
		return;
	
	for (int I = 1; I <= N; I++){
		IORDER[I-1]=I;
		if (Y[I-1] < 1.0e-75)
			Y[I-1]=1.0e-75;	
	}
	
	for (int I = 1; I <= N-1; I ++){
		IMIN=I;
        XMIN=X[I-1];
        for (int J = I+1; J <= N; J++){
			if (X[J-1] < XMIN){
				IMIN=J;
            	XMIN=X[J-1];
			}
		}
		SAVE=X[I-1];
        X[I-1]=X[IMIN-1];
        X[IMIN-1]=SAVE;
        SAVE=Y[I-1];
        Y[I-1]=Y[IMIN-1];
        Y[IMIN-1]=SAVE;
        ISAVE=IORDER[I-1];
        IORDER[I-1]=IORDER[IMIN-1];
        IORDER[IMIN-1]=ISAVE;
        if (I == 1)
        	goto L1;
        
        if ((IORDER[I-1] < IORDER[I-1-1]) && (fabs(X[I-1]-X[I-1-1]) < 1.0e-15)){
  			SAVE=X[I-1-1];
            X[I-1-1]=X[I-1];
            X[I-1]=SAVE;
            SAVE=Y[I-1-1];
            Y[I-1-1]=Y[I-1];
            Y[I-1]=SAVE;
            ISAVE=IORDER[I-1-1];
            IORDER[I-1-1]=IORDER[I-1];
            IORDER[I-1]=ISAVE;	
		}
L1:;
	}
	I2=N;
	if ((IORDER[I2-1] < IORDER[I2-1-1]) && (fabs(X[I2-1]-X[I2-1-1]) < 1.0e-15)){
		SAVE=X[I2-1-1];
        X[I2-1-1]=X[I2-1];
        X[I2-1]=SAVE;
        SAVE=Y[I2-1-1];
        Y[I2-1-1]=Y[I2-1];
        Y[I2-1]=SAVE;	
	}
	
	free(IORDER);
	
//	printf("\n\nSORT22\n\n");
		
}


void graati2_(double E, double &ECS, int *M){
	
/*
Se��o transversal total para espalhamento de f�tons Rayleigh (coerente), em
 cm ** 2. Interpolado a partir dos dados de entrada.
*/	

	const static int NQ = 250;
	const static int NEX = 1024;
	
	double XELN, XEN;
	int  KEN, II,IU, IT;
	
	
	XELN=log(E);
    XEN=1.0e0+(XELN-*CEGRID_.DLEMP1)* *CEGRID_.DLFC;
    KEN=XEN;
    if (KEN == 0) 
		KEN=1;
	
	//busca binaria
	II=CGRA01_.IED[KEN-1];
    IU=CGRA01_.IEU[KEN-1];
L1:;
	IT=(II+IU)/2;

    if (XELN > CGRA01_.ERA[IT-1])
    	II=IT;
    else
    	IU=IT;
    
    if (IU-II > 1)
    	goto L1;
    
    ECS=exp(CGRA01_.XSRA[II-1][*M-1]+(CGRA01_.XSRA[II+1-1][*M-1]-CGRA01_.XSRA[II-1][*M-1])*(XELN-CGRA01_.ERA[II-1])
		/(CGRA01_.ERA[II+1-1]-CGRA01_.ERA[II-1]))/COMPOS_.VMOL[*M-1];
		
    //printf("\n\nGRAATI2\n\n");
		
		
}

void gppa02_(int *M){
	/*
	
Inicializa��o do algoritmo de amostragem para o par el�tron-p�sitron
Produ��o de  por f�tons no cruzamento diferencial de material M. Bethe-Heitler
Se��o .
	
	*/
	
	double SL = 137.035999074e0; //Velocidade da luz (1 / alfa)
	
//N�mero at�mico efetivo.

	double FACT, ALZ, FC, A;
	int IZZ;
	
	FACT=0.0e0;
	for (int I = 1; I <= COMPOS_.NELEM[*M-1]; I++){
    	IZZ=COMPOS_.IZ[I-1][*M-1];
        FACT=FACT+IZZ*CADATA_.ATW[IZZ-1]*COMPOS_.STF[I-1][*M-1];
	}
	
	CGPP00_.ZEQPP[*M-1]=FACT/COMPOS_.AT[*M-1];
    IZZ=CGPP00_.ZEQPP[*M-1]+0.25e0;
    if (IZZ <= 0)
    	IZZ=1;
    if (IZZ > 99)
    	IZZ=99;
    	
	//Corre��o de DBM Coulomb.
	
	ALZ=CGPP00_.ZEQPP[*M-1]/SL;
    A=ALZ*ALZ;
    FC=A*(0.202059e0-A*(0.03693e0-A*(0.00835e0-A*(0.00201e0-A*
       (0.00049e0-A*(0.00012e0-A*0.00003e0)))))+1.0e0/(A+1.0e0));
       
    //Fun��es de triagem e corre��o de baixa energia.
    
	CGPP00_.BCB[*M-1]=2.0e0/CADATA_.RSCR[IZZ-1];
	CGPP00_.F0[1-1][*M-1]=4.0e0*log(CADATA_.RSCR[IZZ-1]);
	CGPP00_.F0[2-1][*M-1]=CGPP00_.F0[1-1][*M-1]-4.0e0*FC;
	
 //	printf("\n\nGPPA02\n\n");
}

void pmrdr2_(){

	//Lê o arquivo de entrada e inicializa PENELOPE e PENGEOM.


    char LINHA[100];
	char KWORD[100];
	char APOIO[100];
	char LIT[3];
	char *PCH;
	char PMFILE[MAXMAT][100];
	char PFILE[21];
	char PFILER[21];
	char SPCDIO[21];
	char SPCDEO[21];
	char PSFDIO[21];
	char SPCFSO[21];
	char SPCAGE[21];
	char BUFFER[121];
	char BUF2[121];
	char KWTITL[7] = "TITLE "; char KWKPAR[7] = "SKPAR "; char KWSENE[7] = "SENERG"; char KWSPEC[7] = "SPECTR"; 
    char KWSPOL[7] = "SGPOL "; char KWSPOS[7] = "SPOSIT"; char KWSBOX[7] = "SBOX  "; char KWSBOD[7] = "SBODY "; 
    char KWSCON[7] = "SCONE "; char KWSREC[7] = "SRECTA"; char KWPSFN[7] = "IPSFN "; char KWPSPL[7] = "IPSPLI"; 
    char KWRRSP[7] = "WGTWIN"; char KWEMAX[7] = "EPMAX "; char KWMATF[7] = "MFNAME"; char KWSIMP[7] = "MSIMPA"; 
    char KWGEOM[7] = "GEOMFN"; char KWGPAR[7] = "PARINP"; char KWSMAX[7] = "DSMAX "; char KWEABS[7] = "EABSB ";
    char KWIFOR[7] = "IFORCE"; char KBRSPL[7] = "IBRSPL"; char KXRSPL[7] = "IXRSPL"; char KWNBE [7] = "NBE   ";
    char KWNBAN[7] = "NBANGL"; char KWIDET[7] = "IMPDET"; char KWIPSF[7] = "IDPSF "; char KWISPC[7] = "IDSPC ";
    char KWIFLN[7] = "IDFLNC"; char KWDIAL[7] = "IDAGEL"; char KWDIAF[7] = "IDAGEF"; char KWIBOD[7] = "IDBODY";
    char KWIPAR[7] = "IDKPAR"; char KWEDET[7] = "ENDETC"; char KWESPC[7] = "EDSPC "; char KWEBOD[7] = "EDBODY";
    char KGRDXX[7] = "GRIDX "; char KGRDYY[7] = "GRIDY "; char KGRDZZ[7] = "GRIDZ "; char KGRDRR[7] = "GRIDR ";
    char KWRESU[7] = "RESUME"; char KWDUMP[7] = "DUMPTO"; char KWDMPP[7] = "DUMPP "; char KWRSEE[7] = "RSEED "; 
    char KWNSIM[7] = "NSIMSH"; char KWTIME[7] = "TIME  "; char KWCOMM[7] = "      "; 


    double PI=3.1415926535897932e0;
	double DE2RA=PI/180.0e0;
	double FSAFE=1.000000001e0;

	static const int NPINPM=500;
	double PARINP[NPINPM];
	static const int NSEM = 1000;
	static const int NB = 5000;

	double SALPHA, SPHI, STHETA, THETLD,THETUD, PHILD,PHIUD, EPMAXR, EAB1, EAB2, EAB3, EMIN, EMAX,EDIL,EDIU, AGEU, AGEL, EDEL, EDEU, XLD, XUD, YLD, YUD, ZLD, ZUD, SAVE,TIMEA, SHNA, CPUTA;

	int KBSMAX, ISEC, KB, ISOURC, NMATR, IHEAD, IP, NMATG, NBE, NBTH, NBPH, NDBOD,NDICH, NAGE,ITST, KPARD, NDECH, NBX, NBY, NBZ, IDOSE, NBR, IRESUM;

    int NPINP = 0;



	for (int I = 1; I <= NPINPM; I++){
		PARINP[I-1] = NPINPM*0.0e0;
	}

	FILE* IWR = fopen("penmain2.dat", "w");
	if (IWR == NULL){
		printf("Nao foi possivel abrir o arquivo penmain2.dat");
		exit(0);
	}

	FILE* IRD = fopen("entrada.in", "r");
	if (IRD == NULL){
		printf("Nao foi possivel abrir o arquivo entrada.in");
		exit(0);
	}
	for (int I = 1; I <= 3; I++){
		CNT0_.PRIM[I-1]=0.0e0;
        CNT0_.PRIM2[I-1]=0.0e0;
		for (int K = 1; K <= 3; K++){
			CNT0_.SEC[I-1][K-1]=0.0e0;
            CNT0_.SEC2[I-1][K-1]=0.0e0;
		}
	}
	

	
	for (int I = 1; I <= 2; I++){
		CNT0_.AVW[I-1]=0.0e0;
        CNT0_.AVW2[I-1]=0.0e0;
        CNT0_.AVA[I-1]=0.0e0;
        CNT0_.AVA2[I-1]=0.0e0;
        CNT0_.AVE[I-1]=0.0e0;
        CNT0_.AVE2[I-1]=0.0e0;
	}


	for (int I = 1; I <= NSEM; I++){
		CNT2_.SHIST[I-1]=0.0e0;
		for (int K = 1; K <= 3; K++){
			CNT3_.SEDS[I-1][K-1]=0.0e0;
            CNT3_.SEDS2[I-1][K-1]=0.0e0;
		}
	}

	for (int I = 1; I <= NB; I++){
		CNT1_.TDEBO[I-1]=0.0e0;
        CNT1_.TDEBO2[I-1]=0.0e0;
        CSPGEO_.EABSB[I-1][1-1]=50.0e0;
 		CSPGEO_.EABSB[I-1][2-1]=50.0e0;
        CSPGEO_.EABSB[I-1][3-1]=50.0e0;
	}

	*CSOUR0_.CTHL=0.0e0;
    *CSOUR0_.DCTH=0.0e0;
    *CSOUR0_.PHIL=0.0e0;
    *CSOUR0_.DPHI=0.0e0;

	for(int KB = 1; KB <= NB; KB++){
		CSOUR3_.IXSBOD[KB-1]=0;
	}

	//Inicialização do contador de tempo.
//	auto begin = std::chrono::high_resolution_clock::now();
	




	//criar a chamada de uma funcao time para contar o tempo

	//Lendo o arquivo de entrda .in

	fprintf(IWR, "\n\n   *************************************************************\n");
	fprintf(IWR, "   **   Program PENMAIN.  Input data and run-time messages.   **\n");
	fprintf(IWR, "   *************************************************************\n\n");
	
    // Criar função para gerar uma data em c++

	fprintf(IWR, "   Date and time: \n");



	fgets(LINHA, sizeof(LINHA), IRD);

    extrairString(KWORD, LINHA, 0, 6);

	PCH = strtok(LINHA, " ");
	PCH = strtok(NULL, " ");
	strcpy(CTITLE_.TITLE, PCH);
	

//	extrairString(TITLE, LINHA, 7, 65);
	
	if (!strcmp(KWORD, KWTITL)){
		fprintf(IWR, "\n   %s\n", CTITLE_.TITLE);
	} else{
		fprintf(IWR, "\nThe input file must begin with the TITLE line\n");
		printf("The input file must begin with the TITLE line\n");
		exit(0);
	}
	fprintf(IWR, "   ------------------------------------------------------------------------\n");

    *CSOUR0_.KPARP=1;
    *CSOUR4_.NPSF=0;
    *CSOUR4_.NPSN=0;
    *CNT2_.NSEB=1;
    *CSOUR4_.NSPLIT=1;
    *CSOUR4_.RLREAD=0.0e0;
    *CSOUR4_.RWGMIN=1.0e35;
    KBSMAX=0;
    ISEC=0;
    *CDUMP_.LDUMP=false;
    *CSOUR0_.LPSF=false;
    *CSOUR2_.LSPEC=false;
    *CSOUR0_.LGPOL=false;
    *CSOUR3_.LEXSRC=false;
    *CSOUR3_.LEXBD=false;
    *CSOUR0_.LSCONE=true;
    *CSOUR0_.JOBEND=0;

	//Descrição da Fonte

L11:;

	fgets(LINHA, sizeof(LINHA), IRD);
	extrairString(KWORD, LINHA, 0, 6);
	extrairString(BUFFER, LINHA, 7, 127);
	if (!strcmp(KWORD, KWCOMM))
		goto L11;
	if (!strcmp(KWORD, KWPSFN))
		goto L21;

	fprintf(IWR, "   >>>>>>>> Source definition.\n");

	if (!strcmp(KWORD, KWKPAR)){
		extrairString(APOIO, BUFFER, 0, 1);
		*CSOUR0_.KPARP = atoi(APOIO);
L12:;
		fgets(LINHA, sizeof(LINHA), IRD);
		extrairString(KWORD, LINHA, 0, 6);
		extrairString(BUFFER, LINHA, 7, 127);
		if (!strcmp(KWORD, KWCOMM))
			goto L12;
		if (!strcmp(KWORD, KWPSFN))
			goto L21;
	}

	if ((*CSOUR0_.KPARP < 0) || (*CSOUR0_.KPARP > 3)){
		fprintf(IWR, "'KPARP = %d\n", *CSOUR0_.KPARP);
		fprintf(IWR, "Incorrect particle type.\n");
		printf("Incorrect particle type.\n");
		exit(0);
	}

	if (*CSOUR0_.KPARP == 1)
		fprintf(IWR, "\n   Primary particles: electrons\n");
	if (*CSOUR0_.KPARP == 2)
		fprintf(IWR, "\n   Primary particles: photons\n");
	if (*CSOUR0_.KPARP == 3)
		fprintf(IWR, "\n   Primary particles: positrons\n");
	if (*CSOUR0_.KPARP == 0)
		fprintf(IWR,"\n   Primary particles: set by the user subroutine SOURCE\n");

	//Fonte Monoenergetica

	if (!strcmp(KWORD, KWSENE)){
		extrairString(APOIO, BUFFER, 7, 15);
		*CSOUR1_.E0 = atof(APOIO);
		fprintf(IWR, "Initial energy = %.6E", *CSOUR1_.E0);
L13:;
		fgets(LINHA, sizeof(LINHA), IRD);
		extrairString(KWORD, LINHA, 0, 6);
		extrairString(BUFFER, LINHA, 7, strlen(LINHA));
		if (!strcmp(KWORD, KWCOMM))
			goto L13;
		if (!strcmp(KWORD, KWPSFN))
			goto L21;	
	} else{
			//Espectro de Energia Continuo
		if (!strcmp(KWORD, KWSPEC)){
			*CSOUR2_.LSPEC=true;
            *CNT2_.NSEB=0;
L14:;
			*CNT2_.NSEB=*CNT2_.NSEB+1;
			if (*CNT2_.NSEB > NSEM){
				fprintf(IWR, "Source energy spectrum.\n");
				fprintf(IWR, "The number of energy bins is too large.\n");
				printf("The number of energy bins is too large.\n");
				exit(0);
			}

			PCH = strtok(BUFFER, " ");
			CSOUR2_.ESRC[*CNT2_.NSEB-1] = atof(PCH);
			PCH = strtok(NULL, " ");
			CSOUR2_.PSRC[*CNT2_.NSEB-1] = atof(PCH);
			CSOUR2_.PSRC[*CNT2_.NSEB-1]=fmax(CSOUR2_.PSRC[*CNT2_.NSEB-1],0.0e0);
L15:;
			fgets(LINHA, sizeof(LINHA), IRD);
			extrairString(KWORD, LINHA, 0, 6);
			extrairString(BUFFER, LINHA, 7, strlen(LINHA));
			if (!strcmp(KWORD, KWCOMM))
				goto L15;
			if (!strcmp(KWORD, KWSPEC))
				goto L14;
			if (!strcmp(KWORD, KWPSFN))
				goto L21;		
		} else{
			*CSOUR1_.E0=1.0e9;
			if (*CSOUR0_.KPARP != 0)
				fprintf(IWR, "%.6E", *CSOUR1_.E0);
		}
	}

	if (*CSOUR2_.LSPEC){
		if (*CNT2_.NSEB > 1){
			sort22_(CSOUR2_.ESRC, CSOUR2_.PSRC, *CNT2_.NSEB);
			fprintf(IWR, "\n   Spectrum:       I    E_low(eV)    E_high(eV)     P_sum(E)\n");
			fprintf(IWR, "                -------------------------------------------------------\n");
			for (int I = 1; I <= *CNT2_.NSEB-1; I++){
				fprintf(IWR, "                %4d  %.6E  %.6E  %.6E\n", I,CSOUR2_.ESRC[I-1],CSOUR2_.ESRC[I+1-1],CSOUR2_.PSRC[I-1]);
			}
			*CSOUR1_.E0=CSOUR2_.ESRC[*CNT2_.NSEB-1];
            *CNT2_.NSEB =*CNT2_.NSEB-1;
			irnd02_(CSOUR2_.PSRC,CSOUR2_.FSRC,CSOUR2_.IASRC,CNT2_.NSEB);

		}else{
			fprintf(IWR, "The source energy spectrum is not defined.\n");
			printf("The source energy spectrum is not defined.\n");
			exit(0);
		}
	}

	if (*CSOUR1_.E0 < 50.0e0){
		fprintf(IWR, "The initial energy E0 is too small\n");
		printf("The initial energy E0 is too small\n");
		exit(0);
	}

	*CSOUR1_.EPMAX=*CSOUR1_.E0;

 /*
 Os pósitrons eventualmente dão raios gama de aniquilação. O máximo
A energia C dos fótons de aniquilação é .lt. 1,21*(E0+me*c**2).
 */

	if (*CSOUR0_.KPARP == 3)
		*CSOUR1_.EPMAX=1.21e0*(*CSOUR1_.E0+5.12e5);


//Efeitos de polarização de fótons (somente para fótons primários).

	if (!strcmp(KWORD, KWSPOL)){
		PCH = strtok(BUFFER, " ");
		*CSOUR1_.SP10 = atof(PCH);
		PCH = strtok(NULL, " ");
		*CSOUR1_.SP20 = atof(PCH);
		PCH = strtok(NULL, " ");
		*CSOUR1_.SP30 = atof(PCH);
		*CSOUR0_.LGPOL=true;
L20:;
		fgets(LINHA, sizeof(LINHA), IRD);
		extrairString(KWORD, LINHA, 0, 6);
		extrairString(BUFFER, LINHA, 7, strlen(LINHA));
		if (!strcmp(KWORD, KWCOMM))
			goto L20;
		if (!strcmp(KWORD, KWPSFN))
			goto L21;

		fprintf(IWR, "   Polarised primary photons. Stokes Parameters:\n");
		fprintf(IWR, "P1 = %.6E (linear polarisation at 45 deg azimuth)\n", *CSOUR1_.SP10);
		fprintf(IWR, "P2 = %.6E (circular polarisation)\n", *CSOUR1_.SP20);
		fprintf(IWR, "P3 = %.6E (linear polarisation at zero azimuth)\n", *CSOUR1_.SP30);			
	}

	//Posição da fonte pontual.

	if (!strcmp(KWORD, KWSPOS)){
		PCH = strtok(BUFFER, " ");
		*CSOUR3_.SX0 = atof(PCH);
		PCH = strtok(NULL, " ");
		*CSOUR3_.SY0 = atof(PCH);
		PCH = strtok(NULL, " ");
		*CSOUR3_.SZ0 = atof(PCH);
L16:;
		fgets(LINHA, sizeof(LINHA), IRD);
		extrairString(KWORD, LINHA, 0, 6);
		extrairString(BUFFER, LINHA, 7, strlen(LINHA));
		if (!strcmp(KWORD, KWCOMM))
			goto L16;
		if (!strcmp(KWORD, KWPSFN))
			goto L21;
	} else{
		*CSOUR3_.SX0=0.0e0;
        *CSOUR3_.SY0=0.0e0;
        *CSOUR3_.SZ0=0.0e0;
	}
	fprintf(IWR, "\n\n   Coordinates of centre:     SX0 =  %.6E cm\n", *CSOUR3_.SX0);
    fprintf(IWR, "                              SY0 =  %.6E cm\n", *CSOUR3_.SY0);
	fprintf(IWR, "                              SZ0 =  %.6E cm\n", *CSOUR3_.SZ0);

	// Fonte Extendida
	if (!strcmp(KWORD, KWSBOX)){
		*CSOUR3_.LEXSRC=true;
		PCH = strtok(BUFFER, " ");
		*CSOUR3_.SSX = atof(PCH);
		PCH = strtok(NULL, " ");
		*CSOUR3_.SSY = atof(PCH);
		PCH = strtok(NULL, " ");
		*CSOUR3_.SSZ = atof(PCH);
		*CSOUR3_.SSX=fabs(*CSOUR3_.SSX);
        *CSOUR3_.SSY=fabs(*CSOUR3_.SSY);
        *CSOUR3_.SSZ=fabs(*CSOUR3_.SSZ);
L17:;
		fgets(LINHA, sizeof(LINHA), IRD);
		extrairString(KWORD, LINHA, 0, 6);
		extrairString(BUFFER, LINHA, 7, strlen(LINHA));
		if (!strcmp(KWORD, KWCOMM))
			goto L17;
		if (!strcmp(KWORD, KWPSFN))
			goto L21;
	} else{
		*CSOUR3_.LEXSRC=false;
		*CSOUR3_.SSX=0.0e0;
        *CSOUR3_.SSY=0.0e0;
        *CSOUR3_.SSZ=0.0e0;
	}
	if (*CSOUR3_.LEXSRC){
		fprintf(IWR, "\n\n   Source size:               SSX = %.6E cm\n", *CSOUR3_.SSX);
   		fprintf(IWR, "                              SSY =  %.6E cm\n", *CSOUR3_.SSY);
		fprintf(IWR, "                              SSZ =  %.6E cm\n", *CSOUR3_.SSZ);
	}

	//Corpos ativos de uma fonte estendida (etiquetas internas do PENGEOM).
L717:;
	if (!strcmp(KWORD, KWSBOD)){
		PCH = strtok(BUFFER, " ");
		KB = atoi(PCH);
		if ((KB < 1) || (KB > NB)){
			fprintf(IWR, "%s %s\n", KWORD, BUFFER);
			fprintf(IWR, "Incorrect body label\n");
			printf("Incorrect body label\n");
			exit(0);
		}
		fprintf(IWR,"   Active body = %4d\n", KB);
		CSOUR3_.IXSBOD[KB-1]=1;
        *CSOUR3_.LEXBD=true;
        KBSMAX=max(KBSMAX,KB);
L766:;
		fgets(LINHA, sizeof(LINHA), IRD);
		extrairString(KWORD, LINHA, 0, 6);
		extrairString(BUFFER, LINHA, 7, strlen(LINHA));
		if (!strcmp(KWORD, KWCOMM))
			goto L766;
		if (!strcmp(KWORD, KWSBOD))
			goto L717;
	}

	//Distribuição angular de partículas primárias.

	ISOURC=0;
L777:;
	if (!strcmp(KWORD, KWSCON)){
		*CSOUR0_.LSCONE = true;
		PCH = strtok(BUFFER, " ");
		STHETA = atof(PCH);
		PCH = strtok(NULL, " ");
		SPHI = atof(PCH);
		PCH = strtok(NULL, " ");
		SALPHA = atof(PCH);

		if ((STHETA < -1.0e-9) || (STHETA-180.0e0 > 1.0e-9)) {
			fprintf(IWR, "%s %s\n", KWORD,BUFFER);
			fprintf(IWR, "   THETA must be between 0 and 180 deg.\n");
			printf("THETA must be between 0 and 180 deg.\n");
			exit(0);
		}
		if ((SPHI < -1.0e-9) || (SPHI-360.0e0 > 1.0e-9)) {
			fprintf(IWR, "%s %s\n", KWORD,BUFFER);
			fprintf(IWR, "   PHI must be between 0 and 360 deg.\n");
			printf("PHI must be between 0 and 360 deg.\n");
			exit(0);
		}
		if ((SALPHA < -1.0e-9) || (SALPHA-180.0e0 > 1.0e-9)) {
			fprintf(IWR, "%s %s\n", KWORD,BUFFER);
			fprintf(IWR, "   ALPHA must be between 0 and 180 deg.\n");
			printf("ALPHA must be between 0 and 180 deg.\n");
			exit(0);
		}

		gcone02_(STHETA*DE2RA,SPHI*DE2RA,SALPHA*DE2RA);
		ISOURC=1;

L18:;
		fgets(LINHA, sizeof(LINHA), IRD);
		extrairString(KWORD, LINHA, 0, 6);
		extrairString(BUFFER, LINHA, 7, strlen(LINHA));
		if (!strcmp(KWORD, KWCOMM))
			goto L18;
		if (!strcmp(KWORD, KWPSFN))
			goto L721;
		goto L777;
	} else if (!strcmp(KWORD, KWSREC)){
			*CSOUR0_.LSCONE=false;
			PCH = strtok(BUFFER, " ");
			THETLD = atof(PCH);
			PCH = strtok(NULL, " ");
			THETUD = atof(PCH);
			PCH = strtok(NULL, " ");
			PHILD = atof(PCH);
			PCH = strtok(NULL, " ");
			PHIUD = atof(PCH);

			if ((fmin(THETLD,THETUD) < -1.0e-9) ||  (fmax(THETLD,THETUD)-180.0e0 > 1.0e-9)){
				fprintf(IWR, "%s %s\n", KWORD,BUFFER);
				fprintf(IWR,"   THETA must be between 0 and 180 deg.\n");
				printf("THETA must be between 0 and 180 deg.\n");
				exit(0);
			}

			if ((fmin(PHILD,PHIUD) < -1.0e-9) ||  (fmax(PHILD,PHIUD)-360.0e0 > 1.0e-9)){
				fprintf(IWR, "%s %s\n", KWORD,BUFFER);
				fprintf(IWR, "   PHI must be between 0 and 360 deg.\n");
				printf("PHI must be between 0 and 360 deg.\n");
				exit(0);
			}

			*CSOUR0_.CTHL=cos(THETLD*DE2RA);
        	*CSOUR0_.DCTH=cos(THETUD*DE2RA)-*CSOUR0_.CTHL;
        	*CSOUR0_.PHIL=PHILD*DE2RA;
        	*CSOUR0_.DPHI=PHIUD*DE2RA-*CSOUR0_.PHIL;
        	ISOURC=1;
L19:;
			fgets(LINHA, sizeof(LINHA), IRD);
			extrairString(KWORD, LINHA, 0, 6);
			extrairString(BUFFER, LINHA, 7, strlen(LINHA));
			if (!strcmp(KWORD, KWCOMM))
				goto L19;
			if (!strcmp(KWORD, KWPSFN))
				goto L721;
			goto L777;
		} else if (ISOURC == 0){
			*CSOUR0_.LSCONE=true;
        	STHETA=0.0e0;
        	SPHI=0.0e0;
        	SALPHA=0.0e0;
			gcone02_(STHETA*DE2RA,SPHI*DE2RA,SALPHA*DE2RA);
		}
	
L721:;
	if (*CSOUR0_.LSCONE){
		fprintf(IWR, "   *** Conical beam:\n   Beam axis direction:     THETA = %.6E deg\n", STHETA);
		fprintf(IWR, "                              PHI = %.6E deg\n", SPHI);
		fprintf(IWR, "   Beam aperture:           ALPHA = %.6E deg\n", SALPHA);
	} else{
		fprintf(IWR, "   *** Rectangular beam:\n   Angular window: THETA = (%.6E , %.6E) deg\n", THETLD,THETUD);
		fprintf(IWR, "                     PHI = (%.6E , %.6E) deg\n", PHILD,PHIUD);
	}

	//Variáveis ​​de estado de partículas lidas de um arquivo de espaço de fase


L21:;

	if (!strcmp(KWORD, KWPSFN)){
		if (*CSOUR0_.KPARP == 0) {
			fprintf(IWR, "   With KPARP=0 (subroutine SOURCE activated)\n");
			fprintf(IWR, "   we cannot read particles from a phase-space file.\n");
			printf("Inconsistent definition of the primary source.\n");
			exit(0);
		}
		*CSOUR4_.NPSF=*CSOUR4_.NPSF+1;
		if (*CSOUR4_.NPSF == 1){
			fprintf(IWR, "   ------------------------------------------------------------------------\n");
			fprintf(IWR, "   >>>>>>  Input phase-space files.\n");
		}
		if (*CSOUR4_.NPSF > NPSFM){
			fprintf(IWR, "   Too many phase-space files.\n");
			printf("Too many phase-space files.\n");
			exit(0);
		}

		PCH = strtok(BUFFER, " ");
		strcpy(CSOUR5_.PSFI[*CSOUR4_.NPSF-1], PCH);
		fprintf(IWR, "   Phase-space file # %4d : %s\n", *CSOUR4_.NPSF,CSOUR5_.PSFI[*CSOUR4_.NPSF-1]);
		*CSOUR0_.LPSF=true;
        *CSOUR2_.LSPEC=false;
L22:;
		fgets(LINHA, sizeof(LINHA), IRD);
		extrairString(KWORD, LINHA, 0, 6);
		extrairString(BUFFER, LINHA, 7, strlen(LINHA));
		if (!strcmp(KWORD, KWCOMM))
			goto L22;
		if (!strcmp(KWORD, KWPSFN))
			goto L21;

		if (!strcmp(KWORD, KWPSPL)){
			PCH = strtok(BUFFER, " ");
			*CSOUR4_.NSPLIT = atoi(PCH);
			if ((*CSOUR4_.NSPLIT > 1) && (*CSOUR4_.NSPLIT <= 1000)){
				fprintf(IWR, "   Particle splitting number = %3d\n", *CSOUR4_.NSPLIT);
			} else{
				*CSOUR4_.NSPLIT=min(1000,abs(*CSOUR4_.NSPLIT)); 
				fprintf(IWR, "   Particle splitting number = %3d (modified)\n", *CSOUR4_.NSPLIT);
			}
L23:;
			fgets(LINHA, sizeof(LINHA), IRD);
			extrairString(KWORD, LINHA, 0, 6);
			extrairString(BUFFER, LINHA, 7, strlen(LINHA));
			if (!strcmp(KWORD, KWCOMM))
				goto L23;
		}

		if (!strcmp(KWORD, KWRRSP)){
			PCH = strtok(BUFFER, " ");
			*CSOUR4_.WGMIN = atof(PCH);
			PCH = strtok(NULL, " ");
			*CSOUR4_.WGMAX = atof(PCH);
			*CSOUR4_.WGMIN=fabs(*CSOUR4_.WGMIN);
            *CSOUR4_.WGMAX=fmin(fabs(*CSOUR4_.WGMAX),1.0e10);
			if (*CSOUR4_.WGMIN > *CSOUR4_.WGMAX){
				fprintf(IWR, "WGMIN = %e\n", *CSOUR4_.WGMIN);
				fprintf(IWR, "WGMAX = %e\n", *CSOUR4_.WGMAX);
				fprintf(IWR, "   Inconsistent window end points.\n");
				printf("Inconsistent window end points.\n");
				exit(0);
			}
			fprintf(IWR, "Initial weight window = ( %.6E, %.6E)\n", *CSOUR4_.WGMIN, *CSOUR4_.WGMAX);
L24:;
			fgets(LINHA, sizeof(LINHA), IRD);
			extrairString(KWORD, LINHA, 0, 6);
			extrairString(BUFFER, LINHA, 7, strlen(LINHA));
			if (!strcmp(KWORD, KWCOMM))
				goto L24;
		}
		*CSOUR4_.RWGMIN=1.0e0 / *CSOUR4_.WGMIN;
	}

	if(*CSOUR0_.KPARP == 0){
		//chamada da funçao SOURCE que cria um arquiv de phases
		*CNT3_.NSDE=400;
        *CNT3_.DSDE=FSAFE * *CSOUR1_.EPMAX / (*CNT3_.NSDE);
        *CNT3_.RDSDE=1.0e0 / *CNT3_.DSDE;
	}

	//Energia maxima de particulas
	if (!strcmp(KWORD, KWEMAX)){
		PCH = strtok(BUFFER, " ");
		EPMAXR = atof(PCH);
		if (*CSOUR0_.KPARP != 0)
			*CSOUR1_.EPMAX=EPMAXR;
		fprintf(IWR, "\n   Maximum particle energy =  %.6E eV\n",*CSOUR1_.EPMAX);
L25:;		
		fgets(LINHA, sizeof(LINHA), IRD);
		extrairString(KWORD, LINHA, 0, 6);
		extrairString(BUFFER, LINHA, 7, strlen(LINHA));
		if (!strcmp(KWORD, KWCOMM))
			goto L25;	
	}else{
		if (*CSOUR0_.LPSF){
			*CSOUR1_.EPMAX=1.0e9;
			fprintf(IWR, "   Maximum particle energy =  %.6E\n", *CSOUR1_.EPMAX);
			fprintf(IWR,"   WARNING: You should have specified the maximum energy EPMAX.\n");
		}
		fprintf(IWR, "   Maximum particle energy =  %.6E\n", *CSOUR1_.EPMAX);
			
	}	

	//Dados de materiais e parâmetros de simulação.	
	fprintf(IWR, "\n\n   ---------------------------------------------------------------------------\n");
	fprintf(IWR, "   >>>>>>  Material data and simulation parameters.\n");

	//Parametros de Simulacao

	for (int M = 1; M <= MAXMAT; M++){
		PENELOPE_mod_.EABS[M-1][1-1]=0.010e0* *CSOUR1_.EPMAX;
        PENELOPE_mod_.EABS[M-1][2-1]=0.001e0* *CSOUR1_.EPMAX;
        PENELOPE_mod_.EABS[M-1][3-1]=0.010e0* *CSOUR1_.EPMAX;
        PENELOPE_mod_.C1[M-1]=0.10e0;
        PENELOPE_mod_.C2[M-1]=0.10e0;
        PENELOPE_mod_.WCC[M-1]=PENELOPE_mod_.EABS[M-1][1-1];
        PENELOPE_mod_.WCR[M-1]=PENELOPE_mod_.EABS[M-1][2-1];
	}

	for (int IB = 1; IB <= NB; IB++){
		CSPGEO_.DSMAX[IB-1]=1.0e20;
	}

	NMATR =0;
L31:;
	if (!strcmp(KWORD, KWMATF)){
		NMATR=NMATR+1;
		PCH = strtok(BUFFER, " ");
		strcpy(PMFILE[NMATR-1], PCH);
L32:;
		fgets(LINHA, sizeof(LINHA), IRD);
		extrairString(KWORD, LINHA, 0, 6);
		extrairString(BUFFER, LINHA, 7, strlen(LINHA));
		if (!strcmp(KWORD, KWCOMM))
			goto L32;
		if (!strcmp(KWORD, KWMATF))
			goto L31;
	}

	if (!strcmp(KWORD, KWSIMP)){
		PCH = strtok(BUFFER, " ");
		PENELOPE_mod_.EABS[NMATR-1][1-1] = atof(PCH);
		PCH = strtok(NULL, " ");
		PENELOPE_mod_.EABS[NMATR-1][2-1] = atof(PCH);
		PCH = strtok(NULL, " ");
		PENELOPE_mod_.EABS[NMATR-1][3-1] = atof(PCH);
		PCH = strtok(NULL, " ");
		PENELOPE_mod_.C1[NMATR-1] = atof(PCH);
		PCH = strtok(NULL, " ");
		PENELOPE_mod_.C2[NMATR-1] = atof(PCH);
		PCH = strtok(NULL, " ");
		PENELOPE_mod_.WCC[NMATR-1] = atof(PCH);
		PCH = strtok(NULL, " ");
		PENELOPE_mod_.WCR[NMATR-1] = atof(PCH);
L33:;
		fgets(LINHA, sizeof(LINHA), IRD);
		extrairString(KWORD, LINHA, 0, 6);
		extrairString(BUFFER, LINHA, 7, strlen(LINHA));
		if (!strcmp(KWORD, KWCOMM))
			goto L33;
		if (!strcmp(KWORD, KWMATF))
			goto L31;
	}

	if (NMATR == 0){
		fprintf(IWR, "%s %s", KWORD, BUFFER);
		fprintf(IWR, "You have to specify a material file (line MFNAME).\n" );
		printf("You have to specify a material file (line MFNAME).\n");
		exit(0);

	}

	if (NMATR > MAXMAT){
		fprintf(IWR, "Wrong number of materials.\n");
		fprintf(IWR, "NMAT = %4d is larger than MAXMAT = %4d",NMATR,MAXMAT );
		printf("Wrong number of materials.");
		exit(0);
	}

	for (int M = 1; M <= NMATR; M++){
		if (M == 1)
			strcpy(LIT, "st");
		if (M == 2)
			strcpy(LIT, "nd");
		if (M == 3)
			strcpy(LIT, "rd");
		if (M > 3)
			strcpy(LIT, "th");
		fprintf(IWR,"\n   **** %2d %s material\n", M, LIT);
		fprintf(IWR, "   Material data file: %s\n", PMFILE[M-1]);
		if (PENELOPE_mod_.EABS[M-1][1-1] < 5.0e1)
			PENELOPE_mod_.EABS[M-1][1-1]=5.0e1;
		if (PENELOPE_mod_.EABS[M-1][2-1] < 5.0e1)
			PENELOPE_mod_.EABS[M-1][2-1]=5.0e1;
		if (PENELOPE_mod_.EABS[M-1][3-1] < 5.0e1)
			PENELOPE_mod_.EABS[M-1][3-1]=5.0e1;
		fprintf(IWR, "   Electron absorption energy = %.6E\n", PENELOPE_mod_.EABS[M-1][1-1]);
		fprintf(IWR, "     Photon absorption energy = %.6E\n", PENELOPE_mod_.EABS[M-1][2-1]);
		fprintf(IWR, "   Positron absorption energy = %.6E\n", PENELOPE_mod_.EABS[M-1][3-1]);
		fprintf(IWR, "   Electron-positron simulation parameters:\n");
		fprintf(IWR, "    C1 = %.6E,       C2 = %.6E\n", PENELOPE_mod_.C1[M-1], PENELOPE_mod_.C2[M-1]);
		fprintf(IWR, "   Wcc = %.6E eV,   Wcr = %.6E eV\n", PENELOPE_mod_.WCC[M-1], PENELOPE_mod_.WCR[M-1]);

	}

	//Inicializando o PENELOPE

	printf("  Initialising PENELOPE 2...\n");

	FILE* MATERIAL2 = fopen("material2.dat", "w");
	if (MATERIAL2 == NULL){
		printf("Nao foi possivel abrir o arquivo material2.dat");
		exit(0);
	}
	int INFO = 0;

	peinit2_(CSOUR1_.EPMAX, NMATR, MATERIAL2, &INFO, PMFILE);

	fclose(MATERIAL2);

	

	// Definicão da Geometria

	printf("  Initialising PENGEOM 2...\n");
	if (!strcmp(KWORD, KWGEOM)){
		PCH = strtok(BUFFER, " ");
		strcpy(PFILE, PCH);
		fprintf(IWR, "\n\n   ---------------------------------------------------------------------------\n");
		fprintf(IWR, "   >>>>>>  Geometry definition.\n   PENGEOM''s geometry file: %s", PFILE);
		FILE* GEOMETRIA = fopen(PFILE, "r");
		if (GEOMETRIA == NULL){
			fprintf(IWR, "Nao foi possivel abrir o arquivo %s\n", PFILE);
			printf("Nao foi possivel abrir o arquivo %s\n", PFILE);
			exit(0);
		}

		NPINP=0;
        IHEAD=0;
L34:;
		fgets(LINHA, sizeof(LINHA), IRD);
		extrairString(KWORD, LINHA, 0, 6);
		extrairString(BUFFER, LINHA, 7, strlen(LINHA));
		if (!strcmp(KWORD, KWCOMM))
			goto L34;
		if (!strcmp(KWORD, KWGPAR)){
			PCH = strtok(BUFFER, " ");
			IP = atoi(PCH);
			if (IP < 1){
				fprintf(IWR, "IP = %4d\n", IP);
				printf("The PARINP index must be positive.\n");
				exit(0);
			}
			NPINP=max(NPINP,IP);
			if (NPINP > NPINPM){
				fprintf(IWR, "Too many modified parameters.\n");
				fprintf(IWR, "NPINP = %4d,  must be less than %4d\n", NPINP,NPINPM);
				printf("Too many modified parameters\n");
				exit(0);	
			}
			PCH = strtok(NULL, " ");
			PARINP[IP-1] = atof(PCH);
			if (IHEAD == 0){
				fprintf(IWR, "   Replaced parameters: PARINP(%4d) = %.6E", IP, PARINP[IP-1]);
				IHEAD=1;
			}else{
				fprintf(IWR, "                        PARINP(%4d) = %.6E", IP, PARINP[IP-1]);
			}
			goto L34;
		}
		FILE* GEOMETRY2 = fopen("geometry2.rep", "w");
		if (GEOMETRY2 == NULL){
			fprintf(IWR, "Nao foi possivel abrir o arquivo geometry2.rep\n");
			printf("Nao foi possivel abrir o arquivo geometry2.rep\n");
			exit(0);
		}

		geomin2_(PARINP, &NPINP, &NMATG, PENGEOM_mod_.NBODY, GEOMETRIA, GEOMETRY2);
		fclose(GEOMETRIA);
		fclose(GEOMETRY2);

	

		if (NMATG < 1){
			fprintf(IWR, "NMATG must be greater than 0.\n");
			printf("NMATG must be greater than 0.\n");
		}

		if (*PENGEOM_mod_.NBODY > NB){
			fprintf(IWR, "      Too many bodies.\n");
			printf("Too many bodies.\n");
			exit(0);
		}

		if (NMATG > *PENELOPE_mod_.NMAT){
			fprintf(IWR, "      Too many different materials.\n");
			printf("Too many different materials.\n");
			exit(0);
		}

		if (KBSMAX > *PENGEOM_mod_.NBODY){
			fprintf(IWR,"      KBSMAX = %4d\n", KBSMAX);
			fprintf(IWR,"      NBODY = %4d\n", *PENGEOM_mod_.NBODY);
			fprintf(IWR,"      Some source bodies are undefined. STOP.\n");
			printf("Some source bodies are undefined. STOP.\n");
		}
	} else{
		fprintf(IWR, "%s %s\n", KWORD,BUFFER);
		fprintf(IWR, "You have to specify a geometry file.\n");
		printf("You have to specify a geometry file\n");
	}



	//Comprimentos máximos de passos de elétrons e pósitrons.

	if (!strcmp(KWORD, KWSMAX)){
L35_1:;
		PCH = strtok(BUFFER, " ");
		KB = atoi(PCH);
		PCH = strtok(NULL, " ");
		CSPGEO_.DSMAX[KB-1] = atof(PCH);
		if ((KB < 1) || (KB > *PENGEOM_mod_.NBODY)){
			fprintf(IWR, "%s %s",KWORD,BUFFER );
			fprintf(IWR, "Incorrect body number.\n");
			printf("Incorrect body number.\n");
			exit(0);
		}

		if (CSPGEO_.DSMAX[KB-1] < 1.0e-7)
			CSPGEO_.DSMAX[KB-1]=1.0e20;
L35_2:;
		fgets(LINHA, sizeof(LINHA), IRD);
		extrairString(KWORD, LINHA, 0, 6);
		extrairString(BUFFER, LINHA, 7, strlen(LINHA));
		if (!strcmp(KWORD, KWCOMM))
			goto L35_2;
		if (!strcmp(KWORD, KWSMAX))
			goto L35_1;	
	}

     /*Energias de absorção local (útil para reduzir trabalho de simulação
       em regiões de menor interesse).*/
    
	for (int IB = 1; IB <= *PENGEOM_mod_.NBODY; IB++){
		int M=PENGEOM_mod_.MATER[IB-1];
		if (M > 0){
			CSPGEO_.EABSB[IB-1][1-1]=PENELOPE_mod_.EABS[M-1][1-1];
            CSPGEO_.EABSB[IB-1][2-1]=PENELOPE_mod_.EABS[M-1][2-1];
            CSPGEO_.EABSB[IB-1][3-1]=PENELOPE_mod_.EABS[M-1][3-1];
		}
	}
	

	if (!strcmp(KWORD, KWEABS)){
L36:;
		PCH = strtok(BUFFER, " ");
		KB = atoi(PCH);
		PCH = strtok(NULL, " ");
		EAB1 = atof(PCH);
		PCH = strtok(NULL, " ");
		EAB2 = atof(PCH);
		PCH = strtok(NULL, " ");
		EAB3 = atof(PCH);

		if ((KB < 1) || (KB > *PENGEOM_mod_.NBODY)){
			fprintf(IWR, "%s %s",KWORD,BUFFER );
			fprintf(IWR, "Incorrect body number.\n");
			printf("Incorrect body number.\n");
			exit(0);
		}

		if (PENGEOM_mod_.MATER[KB-1] > 0){
			CSPGEO_.EABSB[KB-1][1-1]=fmax(CSPGEO_.EABSB[KB-1][1-1], EAB1);
            CSPGEO_.EABSB[KB-1][2-1]=fmax(CSPGEO_.EABSB[KB-1][2-1], EAB2);
            CSPGEO_.EABSB[KB-1][3-1]=fmax(CSPGEO_.EABSB[KB-1][3-1], EAB3);
		}
L37:;
		fgets(LINHA, sizeof(LINHA), IRD);
		extrairString(KWORD, LINHA, 0, 6);
		extrairString(BUFFER, LINHA, 7, strlen(LINHA));
		if (!strcmp(KWORD, KWCOMM))
			goto L37;
		if (!strcmp(KWORD, KWEABS))
			goto L36;	
	}

	fprintf(IWR, "\n\n         Maximum allowed step lengths of electrons and positrons\n         and local absorption energies (non-void bodies)\n\n");
	fprintf(IWR, "   Body    DSMAX(IB)     EABSB(1,IB)    EABSB(2,IB)    EABSB(3,IB)\n");
	fprintf(IWR, "    IB        (cm)           (eV)           (eV)           (eV)\n");
	for (int IB = 1; IB <= *PENGEOM_mod_.NBODY; IB++){
		if (PENGEOM_mod_.MATER[IB-1] > 0)
		    fprintf(IWR, "   %4d  %.6E  %.6E  %.6E  %.6E\n", IB,CSPGEO_.DSMAX[IB-1],
				    CSPGEO_.EABSB[IB-1][1-1],CSPGEO_.EABSB[IB-1][2-1],CSPGEO_.EABSB[IB-1][3-1]);
	}

	//Reduções de Variancia (Não será implementado nessa versão do programa)

	/*Forçando Interação
	IFORCE : Ativa o forçamento de interações do tipo ICOL de partículas
           KPAR no corpo KB. FORCER é o fator forçante, que deve
           ser maior que a unidade. WLOW e WHIG são a parte inferior e superior
           limites da janela de peso onde o forçamento de interação é
           aplicado. Quando vários mecanismos de interação são forçados em
           mesmo corpo, a janela de peso efetivo é igual a
           a interseção das janelas para esses mecanismos.
             PADRÃO: sem forçar interação

           Se o caminho livre médio para interações reais do tipo ICOL é
           MFP, o programa simulará interações deste tipo
           (real ou forçado) com um caminho livre médio efetivo igual a
           MFP/FORCER.

           TRICK: um valor de entrada negativo de FORCER, -FN, é assumido como
           significa que uma partícula com energia E=EPMAX deve interagir,
           em média, +FN vezes no curso de sua desaceleração para
           repouso, para elétrons e pósitrons, ou ao longo de uma média livre
           caminho, para fótons. Isso é muito útil, por exemplo, para gerar
           espectros de raios-x de amostras em massa.

  		O efeito real do forçamento de interação na eficiência não é fácil
  		prever. 
	*/


 /*
     >>>>>>>> Divisão de Bremsstrahlung.

  IBRSPL : Ativa a divisão bremsstrahlung no corpo KB para elétrons
           e pósitrons com pesos na janela (WLOW,WHIG) onde
           força de interação é aplicada. O inteiro IBRSPL é o
           fator de divisão.
             PADRÃO: sem divisão de bremsstrahlung

           Observe que a divisão bremsstrahlung é aplicada em combinação
           com forçamento de interação e, consequentemente, é ativado
           apenas naqueles corpos onde o forçamento de interação está ativo.
 /*
    >>>>>>>> Divisão de raios-X.

  IXRSPL : Divisão de raios X característicos emitidos no corpo KB, de
           qualquer elemento. Cada raio x não dividido com ILB(2)=2 (ou seja, do
           segunda geração) quando extraído da pilha secundária
           é dividido em quanta IXRSPL. Os novos quanta, mais leves, são
           direções aleatórias atribuídas distribuídas isotropicamente.
             PADRÃO: sem divisão de raios-x

 */

	
	
   // Energia e distribuições angulares de partículas .

	fprintf(IWR, "\n   ------------------------------------------------------------------------\n");
    fprintf(IWR, "   >>>>>>  Energy and angular distributions of emerging particles.\n");
	if (!strcmp(KWORD, KWNBE)){
		PCH = strtok(BUFFER, " ");
		EMIN = atof(PCH);
		PCH = strtok(NULL, " ");
		EMAX = atof(PCH);
		PCH = strtok(NULL, " ");
		NBE = atoi(PCH);
L51:;
		fgets(LINHA, sizeof(LINHA), IRD);
		extrairString(KWORD, LINHA, 0, 6);
		extrairString(BUFFER, LINHA, 7, strlen(LINHA));
		if (!strcmp(KWORD, KWCOMM))
			goto L51;
	}else{
		EMIN=0.0e0;
        EMAX = *CSOUR1_.EPMAX;
        NBE=100;
	}

	if (EMIN < 1.0e0)
		EMIN=0.0e0;
	if (EMAX < 1.0e0)
		EMAX= *CSOUR1_.EPMAX;

	if (NBE < 0){
		EMIN=fmax(EMIN,1.0e0);
		fprintf(IWR, "   E:       NBE = %3d,   EMIN = %.6E eV,   EMAX = %.6E eV\n", NBE,EMIN,EMAX);
		fprintf(IWR, "            (logarithmic scale, bin width increases with E)\n");
	}else if (NBE > 0){
		fprintf(IWR, "   E:       NBE = %3d,   EMIN = %.6E eV,   EMAX = %.6E eV\n", NBE,EMIN,EMAX);
		fprintf(IWR, "            (linear scale, uniform bin width)\n");
	}
	if (NBE == 0){
		fprintf(IWR, "   NBE is equal to zero.\n");
		printf("NBE equal to 0\n");
		exit(0);
	}
	if (!strcmp(KWORD, KWNBAN)){
		PCH = strtok(BUFFER, " ");
		NBTH = atoi(PCH);
		PCH = strtok(NULL, " ");
		NBPH = atoi(PCH);
L52:;
		fgets(LINHA, sizeof(LINHA), IRD);
		extrairString(KWORD, LINHA, 0, 6);
		extrairString(BUFFER, LINHA, 7, strlen(LINHA));
		if (!strcmp(KWORD, KWCOMM))
			goto L52;
	} else{
		NBTH=90;
        NBPH=1;
	}
	if (NBTH == 0){
		fprintf(IWR,"   NBTH is equal to zero.\n");
		printf("NBTH equal to 0.\n");
		exit(0);
	}
	if (NBTH > 0)
		fprintf(IWR,"   Theta:  NBTH = %3d (linear scale)\n", NBTH);
	else if (NBTH < 0)
		fprintf(IWR,"   Theta:  NBTH = %3d (logarithmic scale)\n", NBTH);

	fprintf(IWR, "   Phi:    NBPH = %3d", NBPH);
	if (NBPH < 1){
		fprintf(IWR,"   Wrong number of PHI bins.\n");
		printf("Wrong number of PHI bins.\n");
		exit(0);
	}


	enang02_(EMIN,EMAX,NBE,NBTH,NBPH,IWR);
	


	//Detectores de impacto

	for (int KD = 1; KD <= NIDM; KD++){
		CNT4_.KKDI[1-1][KD-1]=0;
        CNT4_.KKDI[2-1][KD-1]=0;
        CNT4_.KKDI[3-1][KD-1]=0;
	}

	NDBOD=0;
    *CNT4_.NPSFO=0;
	*CNT4_.NID=0;
L61:;
	if (!strcmp(KWORD, KWIDET)){
		*CNT4_.NID=*CNT4_.NID+1;
        NDBOD=0;
		fprintf(IWR,"   ------------------------------------------------------------------------\n");
		fprintf(IWR, "   >>>>>>  Impact detector # %2d\n", *CNT4_.NID);
		PCH = strtok(BUFFER, " ");
		EDIL = atof(PCH);
		PCH = strtok(NULL, " ");
		EDIU = atof(PCH);
		PCH = strtok(NULL, " ");
		NDICH = atoi(PCH);
		PCH = strtok(NULL, " ");
		CNT4_.IPSF[*CNT4_.NID-1] = atoi(PCH);
		PCH = strtok(NULL, " ");
		CNT4_.IDCUT[*CNT4_.NID-1] = atoi(PCH);

		if (NDICH == 0){
			fprintf(IWR, "%s %s\n", KWORD,BUFFER);
			fprintf(IWR, "Incorrect number of energy bins.\n");
			printf("Incorrect number of energy bins.\n");
			exit(0);
		}

		if (EDIL < 50.0e0)
			EDIL=50.0e0;

		if (EDIU < 50.0e0)
			EDIU= *CSOUR1_.EPMAX;

		if (NDICH < 0){
			EDIL=max(EDIL,50.0e0);
			fprintf(IWR, "   Energy window = (%.5E, %.5E) eV\n",EDIL,EDIU );
			fprintf(IWR, "   Number of energy bins = %4d   (logarithmic scale)\n", abs(NDICH));
		}else{
			fprintf(IWR, "   Energy window = (%.5E, %.5E) eV\n",EDIL,EDIU );
			fprintf(IWR, "   Number of energy bins = %4d   (linear scale)\n", NDICH);
		}

		if (EDIU < EDIL+1.0e0){
			fprintf(IWR, "%s %s\n", KWORD,BUFFER );
			fprintf(IWR, "Incorrect energy limits.\n");
			printf("Incorrect energy limits.\n");
			exit(0);
		}

		if ((CNT4_.IPSF[*CNT4_.NID-1] < 0) || (CNT4_.IPSF[*CNT4_.NID-1] > 1)){
			fprintf(IWR, "%s %s\n", KWORD,BUFFER );
			fprintf(IWR, "Wrong IPSF value.\n");
			printf("Wrong IPSF value.\n");
			exit(0);
		}

		if ((CNT4_.IDCUT[*CNT4_.NID-1] < 0) || (CNT4_.IDCUT[*CNT4_.NID-1] > 2)){
			fprintf(IWR, "%s %s\n", KWORD,BUFFER );
			fprintf(IWR, "Wrong IDCUT value.\n");
			printf("Wrong IDCUT value.\n");
			exit(0);
		}

L62:;
		fgets(LINHA, sizeof(LINHA), IRD);
		extrairString(KWORD, LINHA, 0, 6);
		extrairString(BUFFER, LINHA, 7, strlen(LINHA));
		if (!strcmp(KWORD, KWCOMM))
			goto L62;

		if (!strcmp(KWORD, KWISPC)){
			if (*CNT4_.NID == 0){
				fprintf(IWR, "%s %s\n", KWORD,BUFFER );
				fprintf(IWR, "No impact detector has been defined yet.\n");
				printf("No impact detector has been defined yet.\n");
				exit(0);
			}
			PCH = strtok(BUFFER, " ");
			strcpy(SPCDIO, PCH);
			fprintf(IWR, "   Output energy spectrum: %s\n", SPCDIO);
L64:;
			fgets(LINHA, sizeof(LINHA), IRD);
			extrairString(KWORD, LINHA, 0, 6);
			extrairString(BUFFER, LINHA, 7, strlen(LINHA));
			if (!strcmp(KWORD, KWCOMM))
				goto L64;
		}else {
			sprintf(BUF2, "%d", 1000 + *CNT4_.NID);
			sprintf(SPCDIO, "spc-impdet-%c%c.dat", BUF2[3], BUF2[4] );
			fprintf(IWR, "   Output energy spectrum: %s\n", SPCDIO);	
		}

		if (!strcmp(KWORD, KWIPSF)){
			if (*CNT4_.NID == 0){
				fprintf(IWR, "%s %s\n", KWORD,BUFFER );
				fprintf(IWR, "No impact detector has been defined yet.\n");
				printf("No impact detector has been defined yet.\n");
				exit(0);
			}
			PCH = strtok(BUFFER, " ");
			strcpy(PSFDIO, PCH);
			if (CNT4_.IPSF[*CNT4_.NID-1] > 0){
				fprintf(IWR, "   Output phase-space file: %s\n", PSFDIO);
				*CNT4_.NPSFO=*CNT4_.NPSFO+1;
				if (*CNT4_.NPSFO > 1){
					fprintf(IWR, "You cannot generate more than one PSF in a single run.\n");
					printf("Only one PSF can be generated in a each run\n");
					exit(0);
				}
			}else{
				fprintf(IWR, "No phase-space file is generated.\n");
			}
L63:;
			fgets(LINHA, sizeof(LINHA), IRD);
			extrairString(KWORD, LINHA, 0, 6);
			extrairString(BUFFER, LINHA, 7, strlen(LINHA));
			if (!strcmp(KWORD, KWCOMM))
				goto L63;
		}else{
			if (CNT4_.IPSF[*CNT4_.NID-1] > 0){
				sprintf(BUF2, "%d", 1000 + *CNT4_.NID);
				sprintf(PSFDIO, "psf-impdet-%c%c.dat", BUF2[3], BUF2[4] );
				fprintf(IWR, "   Output phase-space file: %s\n", PSFDIO);
				*CNT4_.NPSFO=*CNT4_.NPSFO+1;
				if (*CNT4_.NPSFO > 1){
					fprintf(IWR, "You cannot generate more than one PSF in a single run.\n");
					printf("Only one PSF can be generated in a each run\n");
					exit(0);
				}
			}
		}

		if (abs(CNT4_.IDCUT[*CNT4_.NID-1]) == 0){
			fprintf(IWR, "   Detected particles are absorbed\n");
		}else{
			fprintf(IWR, "   Particles are transported through this detector\n");
		}

		strcpy(SPCFSO, "none");

		if (!strcmp(KWORD, KWIFLN)){
			if (*CNT4_.NID == 0){
				fprintf(IWR, "%s %s\n", KWORD,BUFFER );
				fprintf(IWR, "No impact detector has been defined yet.\n");
				printf("No impact detector has been defined yet.\n");
				exit(0);
			}
			PCH = strtok(BUFFER, " ");
			strcpy(SPCFSO, PCH);
L83:;
			fgets(LINHA, sizeof(LINHA), IRD);
			extrairString(KWORD, LINHA, 0, 6);
			extrairString(BUFFER, LINHA, 7, strlen(LINHA));
			if (!strcmp(KWORD, KWCOMM))
				goto L83;

			if (CNT4_.IDCUT[*CNT4_.NID-1] == 2){
				fprintf(IWR, "Output fluence distribution: %s", SPCFSO);
			}
		}else{
			if (CNT4_.IDCUT[*CNT4_.NID-1] == 2){
				sprintf(BUF2, "%d", 1000 + *CNT4_.NID);
				sprintf(SPCFSO, "fln-impdet-%c%c.dat", BUF2[3], BUF2[4] );
				fprintf(IWR, "   Output fluence distribution: %s\n", SPCFSO);
			}
		}

		AGEU=-1.0e6;
        AGEL=0.0e0;
        strcpy(SPCAGE,"none");
        NAGE=0;

		if (!strcmp(KWORD, KWDIAL)){
			if (*CNT4_.NID == 0){
				fprintf(IWR, "%s %s\n", KWORD,BUFFER );
				fprintf(IWR, "No impact detector has been defined yet.\n");
				printf("No impact detector has been defined yet.\n");
				exit(0);
			}
			PCH = strtok(BUFFER, " ");
			AGEL = atof(PCH);
			PCH = strtok(NULL, " ");
			AGEU = atof(PCH);
			PCH = strtok(NULL, " ");
			NAGE = atoi(PCH);

			if (NAGE == 0){
				fprintf(IWR, "%s %s\n", KWORD,BUFFER );
				fprintf(IWR, "Incorrect number of age bins.\n");
				printf("Incorrect number of age bins.\n");
				exit(0);
			}

			AGEL=fmax(AGEL,0.0e0);
			if (NAGE < 0){
				AGEL=fmax(AGEL,1.0e-20);
				fprintf(IWR, "   Particle age window = (%.5E, %.5E) seconds",AGEL,AGEU );
				fprintf(IWR, "   Number of age bins = %.4d (logarithmic scale)", abs(NAGE));
			}else{
				fprintf(IWR, "   Particle age window = (%.5E, %.5E) seconds",AGEL,AGEU );
				fprintf(IWR, "   Number of age bins = %.4d (linear scale)", NAGE);
			}

			if (AGEU < AGEL+1.0e-19){
				fprintf(IWR, "%s %s\n", KWORD,BUFFER );
				fprintf(IWR, "Incorrect age limits.\n");
				printf("Incorrect age limits.\n");
				exit(0);
			}
L84:;	
			fgets(LINHA, sizeof(LINHA), IRD);
			extrairString(KWORD, LINHA, 0, 6);
			extrairString(BUFFER, LINHA, 7, strlen(LINHA));
			if (!strcmp(KWORD, KWCOMM))
				goto L84;
		}

		if (!strcmp(KWORD, KWDIAF)){
			if (*CNT4_.NID == 0){
				fprintf(IWR, "%s %s\n", KWORD,BUFFER );
				fprintf(IWR, "No impact detector has been defined yet.\n");
				printf("No impact detector has been defined yet.\n");
				exit(0);
			}

			if (AGEL < 0.0e0){
				fprintf(IWR, "%s %s\n", KWORD,BUFFER );
				fprintf(IWR, "Undefined age distribution limits.\n");
				printf("Undefined age distribution limits.\n");
				exit(0);
			}
			PCH = strtok(BUFFER, " ");
			strcpy(SPCAGE, PCH);
			fprintf(IWR, "Output age distribution: %s\n", SPCAGE);
L85:;
			
			fgets(LINHA, sizeof(LINHA), IRD);
			extrairString(KWORD, LINHA, 0, 6);
			extrairString(BUFFER, LINHA, 7, strlen(LINHA));
			if (!strcmp(KWORD, KWCOMM))
				goto L85;

		}else{
			if (AGEU > 0.0e0){
				sprintf(BUF2, "%d", 1000 + *CNT4_.NID);
				sprintf(SPCAGE, "age-impdet-%c%c.dat", BUF2[3], BUF2[4] );
				fprintf(IWR, "  Output age distribution: %s\n", SPCAGE);
			}
		}

L65:;

		if (!strcmp(KWORD, KWIBOD)){
			if (*CNT4_.NID == 0){
				fprintf(IWR, "%s %s\n", KWORD,BUFFER );
				fprintf(IWR, "No impact detector has been defined yet.\n");
				printf("No impact detector has been defined yet.\n");
				exit(0);
			}
			PCH = strtok(BUFFER, " ");
			KB = atoi(PCH);

			
			if ((KB < 1) || (KB > *PENGEOM_mod_.NBODY)){
				fprintf(IWR, "%s %s",KWORD,BUFFER );
				fprintf(IWR, "Incorrect body number.\n");
				printf("Incorrect body number.\n");
				exit(0);
			}

			if (PENGEOM_mod_.KDET[KB-1] != 0){
				fprintf(IWR, "%s %s\n", KWORD,BUFFER );
				fprintf(IWR, "   A body cannot be part of two detectors.\n");
				printf("   A body cannot be part of two detectors.\n");
				exit(0);
			}

			if (PENGEOM_mod_.MATER[KB-1] != 0){
				fprintf(IWR, "%s %s\n", KWORD,BUFFER );
				fprintf(IWR, "   A void body cannot be part of two detectors.\n");
				printf("   A void body cannot be part of two detectors.\n");
				exit(0);
			}

			fprintf(IWR, "   Active body = %4d\n", KB);
			PENGEOM_mod_.KDET[KB-1]= *CNT4_.NID;
          	NDBOD=NDBOD+1;
L66:;
			fgets(LINHA, sizeof(LINHA), IRD);
			extrairString(KWORD, LINHA, 0, 6);
			extrairString(BUFFER, LINHA, 7, strlen(LINHA));
			if (!strcmp(KWORD, KWCOMM))
				goto L66;
			if (!strcmp(KWORD, KWIBOD))
				goto L65;
			if (!strcmp(KWORD, KWIDET)){
				if (NDBOD == 0){
					fprintf(IWR, "This detector has no active bodies.\n");
					printf("This detector has no active bodies.\n");
					exit(0);
				}
				ITST=fmax(CNT4_.KKDI[1-1][*CNT4_.NID-1], fmax(CNT4_.KKDI[2-1][*CNT4_.NID-1], CNT4_.KKDI[3-1][*CNT4_.NID-1]));
				if (ITST == 0){
					CNT4_.KKDI[1-1][*CNT4_.NID-1] = 1;
					CNT4_.KKDI[2-1][*CNT4_.NID-1] = 1;
					CNT4_.KKDI[3-1][*CNT4_.NID-1] = 1;
					fprintf(IWR, "   Detected particles = electrons, photons and positrons\n");
				}

				imdet02_(EDIL,EDIU,NDICH,AGEL,AGEU,NAGE,CNT4_.IDCUT[*CNT4_.NID-1],SPCDIO,SPCFSO,SPCAGE, *CNT4_.NID, IWR);
				
				goto L61;
			}
		}
L67:;

		if (!strcmp(KWORD, KWIPAR)){
			if (*CNT4_.NID == 0){
				fprintf(IWR, "%s %s\n", KWORD,BUFFER );
				fprintf(IWR, "No impact detector has been defined yet.\n");
				printf("No impact detector has been defined yet.\n");
				exit(0);
			}
			if (NDBOD == 0){
				fprintf(IWR, "This detector has no active bodies.\n");
				printf("This detector has no active bodies.\n");
				exit(0);
			}
			PCH = strtok(BUFFER, " ");
			KPARD = atoi(PCH);

			if (KPARD == 1){
				CNT4_.KKDI[1-1][*CNT4_.NID-1]=1;
				fprintf(IWR,"   Detected particles = electrons\n");
			}else if (KPARD == 2){
				CNT4_.KKDI[2-1][*CNT4_.NID-1]=1;
				fprintf(IWR,"   Detected particles = photons\n");
			}else if (KPARD == 3){
				CNT4_.KKDI[3-1][*CNT4_.NID-1]=1;
				fprintf(IWR,"   Detected particles = positrons\n");
			}

L68:;
			fgets(LINHA, sizeof(LINHA), IRD);
			extrairString(KWORD, LINHA, 0, 6);
			extrairString(BUFFER, LINHA, 7, strlen(LINHA));
			if (!strcmp(KWORD, KWCOMM))
				goto L68;
			if (!strcmp(KWORD, KWIPAR))
				goto L67;
			if (!strcmp(KWORD, KWIBOD))
				goto L65;
			if (!strcmp(KWORD, KWIDET)){
				if (NDBOD == 0){
					fprintf(IWR, "This detector has no active bodies.\n");
					printf("This detector has no active bodies.\n");
					exit(0);
				}
				ITST=max(CNT4_.KKDI[1-1][*CNT4_.NID-1],max(CNT4_.KKDI[2-1][*CNT4_.NID-1],CNT4_.KKDI[3-1][*CNT4_.NID-1]));
				if (ITST == 0){
					CNT4_.KKDI[1-1][*CNT4_.NID-1]=1;
					CNT4_.KKDI[2-1][*CNT4_.NID-1]=1;
					CNT4_.KKDI[3-1][*CNT4_.NID-1]=1;
					fprintf(IWR, "Detected particles = electrons, photons and positrons\n");
				}

				imdet02_(EDIL,EDIU,NDICH,AGEL,AGEU,NAGE,CNT4_.IDCUT[*CNT4_.NID-1],SPCDIO,SPCFSO,SPCAGE, *CNT4_.NID, IWR);
				
				goto L61;
			}
		}		
	}

	if (*CNT4_.NID > 0){
		if (NDBOD == 0){
			fprintf(IWR, "This detector has no active bodies.\n");
			printf("This detector has no active bodies.\n");
			exit(0);
		}

		ITST=max(CNT4_.KKDI[1-1][*CNT4_.NID-1],max(CNT4_.KKDI[2-1][*CNT4_.NID-1],CNT4_.KKDI[3-1][*CNT4_.NID-1]));
		if (ITST == 0){
			CNT4_.KKDI[1-1][*CNT4_.NID-1]=1;
			CNT4_.KKDI[2-1][*CNT4_.NID-1]=1;
			CNT4_.KKDI[3-1][*CNT4_.NID-1]=1;
			fprintf(IWR, "Detected particles = electrons, photons and positrons\n");
		}
		imdet02_(EDIL,EDIU,NDICH,AGEL,AGEU,NAGE,CNT4_.IDCUT[*CNT4_.NID-1],SPCDIO,SPCFSO,SPCAGE, *CNT4_.NID, IWR);
			
	}



	//Detectores de deposição de energia

	for (int KB = 1; KB <= *PENGEOM_mod_.NBODY; KB++){
		CNT5_.KBDE[KB-1]=0;
	}

	NDBOD=0;
    *CNT5_.NED=0;
L43:;
	if (!strcmp(KWORD, KWEDET)){
		if (*CNT5_.NED > 0){
			if (NDBOD == 0){
				fprintf(IWR, "This detector has no active bodies.\n");
				printf("This detector has no active bodies.\n");
				exit(0);
			}
		}
		*CNT5_.NED=*CNT5_.NED+1;
        NDBOD=0;
		fprintf(IWR,"   ---------------------------------------------------------------------------\n");
		fprintf(IWR, "   '>>>>>>  Energy-deposition detector # %2d", *CNT5_.NED);
		PCH = strtok(BUFFER, " ");
		EDEL = atof(PCH);
		PCH = strtok(NULL, " ");
		EDEU = atof(PCH);
		PCH = strtok(NULL, " ");
		NDECH = atoi(PCH);

		if (NDECH == 0){
			fprintf(IWR, "Incorrect number of energy bins.\n");
			printf("Incorrect number of energy bins.\n");
			exit(0);
		}

		if (EDEL < 1.0e0)
			EDEL=0.0e0;
		if (EDEU < 1.0e0)
			EDEU=*CSOUR1_.EPMAX;

		if (NDECH < 0){
			EDEL=max(EDEL,1.0e0);
			fprintf(IWR, "   'Energy window = (%.5E, %.5E) seconds",EDEL,EDEU );
			fprintf(IWR, "   Number of age bins = %.4d (logarithmic scale)", abs(NDECH));
		}else{
			fprintf(IWR, "   'Energy window = (%.5E, %.5E) seconds",EDEL,EDEU );
			fprintf(IWR, "   Number of age bins = %.4d (linear scale)", NDECH);
		}

		if (EDEU < EDEL+1.0e0){
			fprintf(IWR, "Incorrect energy limits.\n");
			printf("Incorrect energy limits.\n");
			exit(0);
		}

L44:;
		fgets(LINHA, sizeof(LINHA), IRD);
		extrairString(KWORD, LINHA, 0, 6);
		extrairString(BUFFER, LINHA, 7, strlen(LINHA));
		if (!strcmp(KWORD, KWCOMM))
			goto L44;
		if (!strcmp(KWORD, KWESPC)){
			PCH = strtok(BUFFER, " ");
			strcpy(SPCDEO, PCH);
			fprintf(IWR,"Output spectrum: %s\n",SPCDEO );
L45:;
			fgets(LINHA, sizeof(LINHA), IRD);
			extrairString(KWORD, LINHA, 0, 6);
			extrairString(BUFFER, LINHA, 7, strlen(LINHA));
			if (!strcmp(KWORD, KWCOMM))
				goto L45;
		}else{
			sprintf(BUF2, "%d", 1000 + *CNT5_.NED);
			sprintf(SPCDEO, "spc-impdet-%c%c.dat", BUF2[3], BUF2[4] );
			fprintf(IWR, "   Output energy spectrum: %s\n", SPCDEO);
		}

L46:;
		if (!strcmp(KWORD, KWEBOD)){
			PCH = strtok(BUFFER, " ");
			KB = atoi(PCH);
			if ((KB < 1) || (KB > *PENGEOM_mod_.NBODY)){
				fprintf(IWR, "%s %s\n", KWORD, BUFFER);
				fprintf(IWR, "Incorrect body label\n");
				printf("Incorrect body label\n");
				exit(0);
			}
			if (CNT5_.KBDE[KB-1] != 0){
				fprintf(IWR, "%s %s\n", KWORD, BUFFER);
				fprintf(IWR, "A body cannot be part of two detectors.\n");
				printf("A body cannot be part of two detectors.\n");
				exit(0);
			}

			if (PENGEOM_mod_.MATER[KB-1] == 0){
				fprintf(IWR, "%s %s\n", KWORD, BUFFER);
				fprintf(IWR, "A void body cannot be part of a detector.\n");
				printf("A void body cannot be part of a detector.\n");
				exit(0);
			}
			fprintf(IWR, "Active body = %4d", KB);
			if ((CFORCI_.LFORCE[1-1][KB-1]) || (CFORCI_.LFORCE[2-1][KB-1]) || (CFORCI_.LFORCE[3-1][KB-1])){
				fprintf(IWR, "   #  WARNING: Spectrum may be strongly  biased\n");
				fprintf(IWR, "              when interaction forcing is used!\n");
			}

			CNT5_.KBDE[KB-1]=*CNT5_.NED;
        	NDBOD=NDBOD+1;
		}
L47:;
		fgets(LINHA, sizeof(LINHA), IRD);
		extrairString(KWORD, LINHA, 0, 6);
		extrairString(BUFFER, LINHA, 7, strlen(LINHA));
		if (!strcmp(KWORD, KWCOMM))
			goto L47;
		if (!strcmp(KWORD, KWEBOD))
			goto L46;
		if (!strcmp(KWORD, KWEDET)){
			endet02_(EDEL,EDEU,NDECH,SPCDEO,*CNT5_.NED, IWR);
			goto L43;
		}
	}

	if (*CNT5_.NED > 0){
		if (NDBOD == 0){
			fprintf(IWR, "This detector has no active bodies.\n");
			printf("This detector has no active bodies.\n");
			exit(0);
		}
		endet02_(EDEL,EDEU,NDECH,SPCDEO,*CNT5_.NED, IWR);
	}

	//Distribuição de dose


    *CNT6_.LDOSEM=false;
    NBX=0;
    NBY=0;
    NBZ=0;

	if (!strcmp(KWORD, KGRDXX)){
		*CNT6_.LDOSEM=true;
        IDOSE=1;
		fprintf(IWR,"\n\n   ------------------------------------------------------------------------\n");
		fprintf(IWR, "   >>>>>>  Dose distribution in a box.\n");
		PCH = strtok(BUFFER, " ");
		XLD = atof(PCH);
		PCH = strtok(NULL, " ");
		XUD = atof(PCH);
		PCH = strtok(NULL, " ");
		NBX = atoi(PCH);

		if (XLD > XUD){
			SAVE=XLD;
            XLD=XUD;
            XUD=SAVE;
		}
		if (XUD < XLD+1.0e-6){
			fprintf(IWR, "%s %s", KWORD,BUFFER);
			fprintf(IWR," XU must be greater than XL+1.0E-6.\n");
			printf("XU must be greater than XL+1.0E-6.\n");
		}
		NBX=max(1,NBX);
		fprintf(IWR, "   XL = %.6E  cm,  XU = %.6E cm,  NDBX = %.3d\n", XLD,XUD,NBX);
L70:;
		fgets(LINHA, sizeof(LINHA), IRD);
		extrairString(KWORD, LINHA, 0, 6);
		extrairString(BUFFER, LINHA, 7, strlen(LINHA));
		if (!strcmp(KWORD, KWCOMM))
			goto L70;
	}


	if (!strcmp(KWORD, KGRDYY)){
		PCH = strtok(BUFFER, " ");
		YLD = atof(PCH);
		PCH = strtok(NULL, " ");
		YUD = atof(PCH);
		PCH = strtok(NULL, " ");
		NBY = atoi(PCH);

		if (YLD > YUD){
			SAVE=YLD;
            YLD=YUD;
            YUD=SAVE;
		}
		if (NBX < 1){
			fprintf(IWR, "%s %s", KWORD,BUFFER);
			fprintf(IWR," Incorrect keyword.\n");
			printf("Incorect keyword.\n");
			exit(0);
		}
		if (YUD < YLD+1.0e-6){
			fprintf(IWR, "%s %s", KWORD,BUFFER);
			fprintf(IWR," YU must be greater than YL+1.0E-6.\n");
			printf("YU must be greater than YL+1.0E-6.\n");
			exit(0);
		}
		NBY=max(1,NBY);
		fprintf(IWR, "   YL = %.6E  cm,  YU = %.6E cm,  NDBY = %.3d\n", YLD,YUD,NBY);
L71:;
		fgets(LINHA, sizeof(LINHA), IRD);
		extrairString(KWORD, LINHA, 0, 6);
		extrairString(BUFFER, LINHA, 7, strlen(LINHA));
		if (!strcmp(KWORD, KWCOMM))
			goto L71;
	}else{
		if (NBX > 0){
			fprintf(IWR, "%s %s", KWORD,BUFFER);
			fprintf(IWR," Incorrect keyword.\n");
			printf("Incorect keyword.\n");
			exit(0);
		}
	}

	if (!strcmp(KWORD, KGRDZZ)){
		PCH = strtok(BUFFER, " ");
		ZLD = atof(PCH);
		PCH = strtok(NULL, " ");
		ZUD = atof(PCH);
		PCH = strtok(NULL, " ");
		NBZ = atoi(PCH);

		if (ZLD > ZUD){
			SAVE=ZLD;
            ZLD=ZUD;
            ZUD=SAVE;
		}
		if (ZUD < ZLD+1.0e-6){
			fprintf(IWR, "%s %s", KWORD,BUFFER);
			fprintf(IWR," ZU must be greater than ZL+1.0E-6.\n");
			printf("ZU must be greater than ZL+1.0E-6.\n");
			exit(0);
		}
		NBZ=max(1,NBZ);
		if (NBX > 0){
			fprintf(IWR, "   ZL = %.6E  cm,  ZU = %.6E ' cm,  NDBZ = %.3d\n", ZLD,ZUD,NBZ);
		}else{
			*CNT6_.LDOSEM=true;
         	IDOSE=2;
			fprintf(IWR,"   ------------------------------------------------------------------------\n");
			fprintf(IWR, "   >>>>>>  Dose distribution in a cylinder.\n");
			fprintf(IWR, "   ZL = %.6E  cm,  ZU = %.6E cm,  NDBZ = %.3d\n", ZLD,ZUD,NBZ);
		}
L72:;
		fgets(LINHA, sizeof(LINHA), IRD);
		extrairString(KWORD, LINHA, 0, 6);
		extrairString(BUFFER, LINHA, 7, strlen(LINHA));
		if (!strcmp(KWORD, KWCOMM))
			goto L72;
	}else{
		if (NBX > 0){
			fprintf(IWR, "%s %s", KWORD,BUFFER);
			fprintf(IWR," Unrecognized keyword.\n");
			printf("Unrecognized keyword.\n");
			exit(0);
		}
	}

	if (!strcmp(KWORD, KGRDRR)){
		PCH = strtok(BUFFER, " ");
		XUD = atof(PCH);
		PCH = strtok(NULL, " ");
		NBR = atof(PCH);

		if (XUD < 1.0e-6){
			fprintf(IWR, "%s %s", KWORD,BUFFER);
			fprintf(IWR," RU must be greater than 1.0E-6.\n");
			printf("RU must be greater than 1.0E-6.\n");
			exit(0);
		}

		if ((NBX > 0) ||(NBY > 0)){
			fprintf(IWR, "%s %s", KWORD,BUFFER);
			fprintf(IWR," Incorrect keyword.\n");
			printf("Incorrect keyword.\n");
			exit(0);
		}

		NBX=max(1,NBR);
        XLD=0.0e0;

		if (NBZ > 0){
			fprintf(IWR, "   RU = %.6E cm,  NDBR = %4d\n", XUD,NBX);
			NBY=1;
		}else{
			*CNT6_.LDOSEM=true;
          	IDOSE=3;
			fprintf(IWR,"   ------------------------------------------------------------------------\n");
			fprintf(IWR, "   >>>>>>  Dose distribution in a sphere.\n");
			fprintf(IWR, "   RU = %.6E cm,  NDBR = %3d\n", XUD,NBX);
			NBY=1;
        	NBZ=1;
		}
L73:;

		fgets(LINHA, sizeof(LINHA), IRD);
		extrairString(KWORD, LINHA, 0, 6);
		extrairString(BUFFER, LINHA, 7, strlen(LINHA));
		if (!strcmp(KWORD, KWCOMM))
			goto L73;
	}

	if (*CNT6_.LDOSEM)
		dose02_(XLD,XUD,YLD,YUD,ZLD,ZUD,NBX,NBY,NBZ,IDOSE,IWR);



	//Caracteristicas do trabalho

	fprintf(IWR,"\n\n   ------------------------------------------------------------------------\n");
	fprintf(IWR, "   >>>>>>  Job characteristics.\n");

	IRESUM=0;

	if (!strcmp(KWORD, KWRESU)){
		PCH = strtok(BUFFER, " ");
		strcpy(PFILER, PCH);
		fprintf(IWR,"   Resume simulation from previous dump file: %s\n", PFILER);
		IRESUM=1;
L75:; 
		fgets(LINHA, sizeof(LINHA), IRD);
		extrairString(KWORD, LINHA, 0, 6);
		extrairString(BUFFER, LINHA, 7, strlen(LINHA));
		if (!strcmp(KWORD, KWCOMM))
			goto L75;
	}

	*CNTRL_.DUMPP=1.0e15;
	if (!strcmp(KWORD, KWDUMP)){
		PCH = strtok(BUFFER, " ");
		strcpy(CDUMP_.PFILED, PCH);
		fprintf(IWR,"   Write final counter values on the dump file: %s\n", CDUMP_.PFILED);
		*CDUMP_.LDUMP=true;
L76:;
		fgets(LINHA, sizeof(LINHA), IRD);
		extrairString(KWORD, LINHA, 0, 6);
		extrairString(BUFFER, LINHA, 7, strlen(LINHA));
		if (!strcmp(KWORD, KWCOMM))
			goto L76;	
	}

	if (!strcmp(KWORD, KWDMPP)){
		PCH = strtok(BUFFER, " ");
		*CNTRL_.DUMPP = atof(PCH);

		if (*CDUMP_.LDUMP){
			if (*CNTRL_.DUMPP < 15.0e0)
				*CNTRL_.DUMPP=15.0e0;
			if (*CNTRL_.DUMPP > 86400.0e0)
				*CNTRL_.DUMPP=86400.0e0;
			fprintf(IWR, "   Dumping period: DUMPP = %.6E\n", *CNTRL_.DUMPP);
		}
L77:;
		fgets(LINHA, sizeof(LINHA), IRD);
		extrairString(KWORD, LINHA, 0, 6);
		extrairString(BUFFER, LINHA, 7, strlen(LINHA));
		if (!strcmp(KWORD, KWCOMM))
			goto L77;	
	}

	if (!strcmp(KWORD, KWRSEE)){
		PCH = strtok(BUFFER, " ");
		*RSEED_.ISEED1 = atoi(PCH);
		PCH = strtok(NULL, " ");
		*RSEED_.ISEED2 = atoi(PCH);

		if (*RSEED_.ISEED1 < 0)
		    rand02_(-*RSEED_.ISEED1);
L79:;
		fgets(LINHA, sizeof(LINHA), IRD);
		extrairString(KWORD, LINHA, 0, 6);
		extrairString(BUFFER, LINHA, 7, strlen(LINHA));
		if (!strcmp(KWORD, KWCOMM))
			goto L79;		
	}else{
		*RSEED_.ISEED1=1;
		*RSEED_.ISEED2=1;
	}
	fprintf(IWR,"   Random-number generator seeds = %10d, %10d\n",*RSEED_.ISEED1, *RSEED_.ISEED2);


	if (!strcmp(KWORD, KWNSIM)){
		PCH = strtok(BUFFER, " ");
		*CNTRL_.DSHN = atof(PCH);
		if (*CNTRL_.DSHN < 1.0e0)
			*CNTRL_.DSHN=2.0e9;
L78:;
		fgets(LINHA, sizeof(LINHA), IRD);
		extrairString(KWORD, LINHA, 0, 6);
		extrairString(BUFFER, LINHA, 7, strlen(LINHA));
		if (!strcmp(KWORD, KWCOMM))
			goto L78;
	}else{
		*CNTRL_.DSHN=2.0e9;
	}
	fprintf( IWR,  "\n   Number of showers to be simulated = %.6E\n", *CNTRL_.DSHN);

	if (!strcmp(KWORD, KWTIME)){
		PCH = strtok(BUFFER, " ");
		TIMEA = atof(PCH);
	} else
		TIMEA = 2.0e9;

	if (TIMEA < 1.0e0)
		TIMEA=100.0e0;
	fprintf(IWR, "   Computation time available = %.6E sec\n", TIMEA);

//	auto end = std::chrono::high_resolution_clock::now();
  //  auto TSECIN = std::chrono::duration_cast<std::chrono::seconds>(end - begin);

	//double TSECIN;
	//timer2_(TSECIN); analisar como chamar uma função para controlar o tempo de execucao

//	*CNTRL_.TSECA=TIMEA+TSECIN.count();
  //  *CNTRL_.TSECAD=TSECIN.count();

	//end = clock();
	//double tempo = (double)(end - start) / CLOCKS_PER_SEC;
	
	

	//printf("TEMPO DE EXECUCAO: %f\n",*CNTRL_.TSECAD );
	

	fprintf(IWR, "\n   -----------------------------------------------------------------------------\n");

	//Se 'RESUME' estiver ativo, leia o gerado anteriormente


	SHNA=0.0e0;
    CPUTA=0.0e0;

    *CNTRL_.N=0;

    *CNT4_.RLAST=0.0e0;
    *CNT4_.RWRITE=0.0e0;

	char LINHADUMP[50];
	double aux;
	string line;




//	printf("PFILER: %s\n", PFILER);
	if (IRESUM == 1){

		FILE* IWDUMP = fopen("iwdump.dmp", "w");
		if (IWDUMP == NULL){
			printf("Não foi possível abrir o arquivo %s\n", "iwdump.dmp");
			exit(0);
		}

		ifstream inFILE;
	    string line;

	

		inFILE.open(PFILER);
		getline(inFILE, line);

		

	    char *sArray = (char *) malloc((line.length()+1) * sizeof(char));

		if (line.length() == 0){
		//	printf("\nLinha dump1.dmp zerada\n");
			
			goto L91;

		}
		
		strcpy(sArray, line.c_str());
		PCH = strtok(sArray, " ");
		SHNA = atof(PCH);
		PCH = strtok(NULL, " ");
        CPUTA = atof(PCH);
		printf("  Reading the DUMP file ...\n");
		fprintf(IWDUMP, "	%f		%f\n", SHNA, CPUTA);
		

        getline(inFILE, line);
		sArray = (char *) realloc(NULL, (line.length()+1) * sizeof(char));
		strcpy(sArray, line.c_str());
		PCH = strtok(sArray, " ");
		strcpy(CTITLE_.TITLE2, PCH);
		fprintf(IWDUMP, "%s\n", CTITLE_.TITLE2);
		if (strcmp(CTITLE_.TITLE2, CTITLE_.TITLE)){
			fprintf(IWR, "The dump file is corrupted (the TITLE does not match).\n");
			printf("The dump file is corrupted (the TITLE does not match).\n");
			exit(0);
		}

		getline(inFILE, line);
		sArray = (char *) realloc(NULL, (line.length()+1) * sizeof(char));
		strcpy(sArray, line.c_str());
		PCH = strtok(sArray, " ");
		*RSEED_.ISEED1 = atoi(PCH);
		PCH = strtok(NULL, " ");
		*RSEED_.ISEED1 = atoi(PCH);
		fprintf(IWDUMP, "	%d		%d\n", *RSEED_.ISEED1, *RSEED_.ISEED2);

		getline(inFILE, line);
		sArray = (char *) realloc(NULL, (line.length()+1) * sizeof(char));
		strcpy(sArray, line.c_str());
		PCH = strtok(sArray, " ");
		*CSOUR4_.NPSN = atoi(PCH);
		PCH = strtok(NULL, " ");
		*CSOUR4_.RLREAD = atof(PCH);
		fprintf(IWDUMP, "	%d		%f\n", *CSOUR4_.NPSN , *CSOUR4_.RLREAD);

		getline(inFILE, line);
		sArray = (char *) realloc(NULL, (line.length()+1) * sizeof(char));
		strcpy(sArray, line.c_str());
		PCH = strtok(sArray, " ");
		*CSOUR0_.KPARP = atoi(PCH);
		fprintf(IWDUMP, "	%d\n", *CSOUR0_.KPARP);

		if (*CSOUR0_.KPARP == 0){
			getline(inFILE, line);
			sArray = (char *) realloc(NULL, (line.length()+1) * sizeof(char));
			strcpy(sArray, line.c_str());
			PCH = strtok(sArray, " ");
			*CNT3_.NSDE = atoi(PCH);
			PCH = strtok(NULL, " ");
			*CNT3_.DSDE = atof(PCH);
			PCH = strtok(NULL, " ");
			*CNT3_.RDSDE = atof(PCH);
			fprintf(IWDUMP, "	%d		%.10f		%.10f\n", *CNT3_.NSDE, *CNT3_.DSDE,*CNT3_.RDSDE );


			getline(inFILE, line);
			sArray = (char *) realloc(NULL, (line.length()+1) * sizeof(char));
			strcpy(sArray, line.c_str());
			PCH = strtok(sArray, " ");
			for  (int I = 1; I <= *CNT3_.NSDE; I++){
				for (int K = 1; K <= 3; K++){
					CNT3_.SEDS[I-1][K-1] = atof(PCH);
					fprintf(IWDUMP, " %.10f		", CNT3_.SEDS[I-1][K-1] );
					PCH = strtok(NULL, " ");
				}
			}
			fprintf(IWDUMP, "\n");

			getline(inFILE, line);
			sArray = (char *) realloc(NULL, (line.length()+1) * sizeof(char));
			strcpy(sArray, line.c_str());
			PCH = strtok(sArray, " ");
			for  (int I = 1; I <= *CNT3_.NSDE; I++){
				for (int K = 1; K <= 3; K++){
					CNT3_.SEDS2[I-1][K-1] = atof(PCH);
					fprintf(IWDUMP, " %.10f		", CNT3_.SEDS2[I-1][K-1] );
					PCH = strtok(NULL, " ");
				}
			}
			fprintf(IWDUMP, "\n");
		} else{
			getline(inFILE, line);
			sArray = (char *) realloc(NULL, (line.length()+1) * sizeof(char));
			strcpy(sArray, line.c_str());
			PCH = strtok(sArray, " ");
			*CNT2_.NSEB = atoi(PCH);
			fprintf(IWDUMP, "	%d\n", *CNT2_.NSEB);

			getline(inFILE, line);
			sArray = (char *) realloc(NULL, (line.length()+1) * sizeof(char));
			strcpy(sArray, line.c_str());
			PCH = strtok(sArray, " ");
			for (int I = 1; I <= *CNT2_.NSEB+1; I++){
				CSOUR2_.ESRC[I-1] = atof(PCH);
				fprintf(IWDUMP, " %.10f		",  CSOUR2_.ESRC[I-1]);
				PCH = strtok(NULL, " ");
			}
			fprintf(IWDUMP, "\n");

			getline(inFILE, line);
			sArray = (char *) realloc(NULL, (line.length()+1) * sizeof(char));
			strcpy(sArray, line.c_str());
			PCH = strtok(sArray, " ");
			for (int I = 1; I <= *CNT2_.NSEB+1; I++){
				CSOUR2_.PSRC[I-1] = atof(PCH);
				fprintf(IWDUMP, " %.10E		",  CSOUR2_.PSRC[I-1]);
				PCH = strtok(NULL, " ");
			}
			fprintf(IWDUMP, "\n");

			getline(inFILE, line);
			sArray = (char *) realloc(NULL, (line.length()+1) * sizeof(char));
			strcpy(sArray, line.c_str());
			PCH = strtok(sArray, " ");
			for (int I = 1; I <= *CNT2_.NSEB+1; I++){
				CNT2_.SHIST[I-1] = atof(PCH);
				fprintf(IWDUMP, " %.10f		",  CNT2_.SHIST[I-1]);
				PCH = strtok(NULL, " ");
			}
			fprintf(IWDUMP, "\n");
		}

		getline(inFILE, line);
		sArray = (char *) realloc(NULL, (line.length()+1) * sizeof(char));
		strcpy(sArray, line.c_str());
		PCH = strtok(sArray, " ");
		for (int I = 1; I <= 3; I++){
			CNT0_.PRIM[I-1] = atof(PCH);
			fprintf(IWDUMP, " %.10f		",  CNT0_.PRIM[I-1]);
			PCH = strtok(NULL, " ");
		}
		fprintf(IWDUMP, "\n");

		getline(inFILE, line);
		sArray = (char *) realloc(NULL, (line.length()+1) * sizeof(char));
		strcpy(sArray, line.c_str());
		PCH = strtok(sArray, " ");
		for (int I = 1; I <= 3; I++){
			CNT0_.PRIM2[I-1] = atof(PCH);
			fprintf(IWDUMP, " %.10f		",  CNT0_.PRIM2[I-1]);
			PCH = strtok(NULL, " ");
		}
		fprintf(IWDUMP, "\n");

		getline(inFILE, line);
		sArray = (char *) realloc(NULL, (line.length()+1) * sizeof(char));
		strcpy(sArray, line.c_str());
		PCH = strtok(sArray, " ");
		for (int I = 1; I <= 3; I++){
			for (int K = 1; K <= 3; K++){
				CNT0_.SEC[I-1][K-1] = atof(PCH);
				fprintf(IWDUMP, " %.10f		",  CNT0_.SEC[I-1][K-1]);
				PCH = strtok(NULL, " ");
			}
		}
		fprintf(IWDUMP, "\n");

		getline(inFILE, line);
		sArray = (char *) realloc(NULL, (line.length()+1) * sizeof(char));
		strcpy(sArray, line.c_str());
		PCH = strtok(sArray, " ");
		for (int I = 1; I <= 3; I++){
			for (int K = 1; K <= 3; K++){
				CNT0_.SEC2[I-1][K-1] = atof(PCH);
				fprintf(IWDUMP, " %.10f		",  CNT0_.SEC2[I-1][K-1]);
				PCH = strtok(NULL, " ");
			}
		}
		fprintf(IWDUMP, "\n");

		getline(inFILE, line);
		sArray = (char *) realloc(NULL, (line.length()+1) * sizeof(char));
		strcpy(sArray, line.c_str());
		PCH = strtok(sArray, " ");
		for (int I = 1; I <= 2; I++){
			CNT0_.AVW[I-1] = atof(PCH);
			fprintf(IWDUMP, " %.10f		",  CNT0_.AVW[I-1]);
			PCH = strtok(NULL, " ");
		}
		fprintf(IWDUMP, "\n");

		getline(inFILE, line);
		sArray = (char *) realloc(NULL, (line.length()+1) * sizeof(char));
		strcpy(sArray, line.c_str());
		PCH = strtok(sArray, " ");
		for (int I = 1; I <= 2; I++){
			CNT0_.AVW2[I-1] = atof(PCH);
			fprintf(IWDUMP, " %.10f		",  CNT0_.AVW2[I-1]);
			PCH = strtok(NULL, " ");
		}
		fprintf(IWDUMP, "\n");

		getline(inFILE, line);
		sArray = (char *) realloc(NULL, (line.length()+1) * sizeof(char));
		strcpy(sArray, line.c_str());
		PCH = strtok(sArray, " ");
		for (int I = 1; I <= 2; I++){
			CNT0_.AVA[I-1] = atof(PCH);
			fprintf(IWDUMP, " %.10f		",  CNT0_.AVA[I-1]);
			PCH = strtok(NULL, " ");
		}
		fprintf(IWDUMP, "\n");

		getline(inFILE, line);
		sArray = (char *) realloc(NULL, (line.length()+1) * sizeof(char));
		strcpy(sArray, line.c_str());
		PCH = strtok(sArray, " ");
		for (int I = 1; I <= 2; I++){
			CNT0_.AVA2[I-1] = atof(PCH);
			fprintf(IWDUMP, " %.10f		",  CNT0_.AVA2[I-1]);
			PCH = strtok(NULL, " ");
		}
		fprintf(IWDUMP, "\n");

		getline(inFILE, line);
		sArray = (char *) realloc(NULL, (line.length()+1) * sizeof(char));
		strcpy(sArray, line.c_str());
		PCH = strtok(sArray, " ");
		for (int I = 1; I <= 2; I++){
			CNT0_.AVE[I-1] = atof(PCH);
			fprintf(IWDUMP, " %.10f		",  CNT0_.AVE[I-1]);
			PCH = strtok(NULL, " ");
		}
		fprintf(IWDUMP, "\n");

		getline(inFILE, line);
		sArray = (char *) realloc(NULL, (line.length()+1) * sizeof(char));
		strcpy(sArray, line.c_str());
		PCH = strtok(sArray, " ");
		for (int I = 1; I <= 2; I++){
			CNT0_.AVE2[I-1] = atof(PCH);
			fprintf(IWDUMP, " %.10f		",  CNT0_.AVE2[I-1]);
			PCH = strtok(NULL, " ");
		}
		fprintf(IWDUMP, "\n");

		getline(inFILE, line);
		sArray = (char *) realloc(NULL, (line.length()+1) * sizeof(char));
		strcpy(sArray, line.c_str());
		PCH = strtok(sArray, " ");
		*PENGEOM_mod_.NBODY = atoi(PCH);
		fprintf(IWDUMP, "		 %d\n",  *PENGEOM_mod_.NBODY);

		getline(inFILE, line);
		sArray = (char *) realloc(NULL, (line.length()+1) * sizeof(char));
		strcpy(sArray, line.c_str());
		PCH = strtok(sArray, " ");
		for (int I = 1; I <= *PENGEOM_mod_.NBODY; I++){
			CNT1_.TDEBO[I-1] = atof(PCH);
			fprintf(IWDUMP, " %.10f		",  CNT1_.TDEBO[I-1]);
			PCH = strtok(NULL, " ");
		}
		fprintf(IWDUMP, "\n");

		getline(inFILE, line);
		sArray = (char *) realloc(NULL, (line.length()+1) * sizeof(char));
		strcpy(sArray, line.c_str());
		PCH = strtok(sArray, " ");
		for (int I = 1; I <= *PENGEOM_mod_.NBODY; I++){
			CNT1_.TDEBO2[I-1] = atof(PCH);
			fprintf(IWDUMP, " %.10f		",  CNT1_.TDEBO2[I-1]);
			PCH = strtok(NULL, " ");
		}
		fprintf(IWDUMP, "\n");

		enangr2_(inFILE, IWDUMP); //Energia e distribuições angulares

		getline(inFILE, line);
		sArray = (char *) realloc(NULL, (line.length()+1) * sizeof(char));
		strcpy(sArray, line.c_str());
		PCH = strtok(sArray, " ");
		*CNT4_.NID = atoi(PCH);
		PCH = strtok(NULL, " ");
		*CNT5_.NED = atoi(PCH);
		PCH = strtok(NULL, " ");
		strcpy(APOIO, PCH);
		if (!strcmp("T", APOIO)){
			*CNT6_.LDOSEM = true;
		}
		else {
			*CNT6_.LDOSEM = false;
		}
		fprintf(IWDUMP, "		 %d			%d			%d\n",  *CNT4_.NID, *CNT5_.NED ,*CNT6_.LDOSEM);

		if (*CNT4_.NID > 0){
			getline(inFILE, line);
			sArray = (char *) realloc(NULL, (line.length()+1) * sizeof(char));
			strcpy(sArray, line.c_str());
			PCH = strtok(sArray, " ");
			*CNT4_.RLAST = atof(PCH);
			PCH = strtok(NULL, " ");
			*CNT4_.RWRITE = atof(PCH);
			imdetr2_(inFILE, IWDUMP); //Detectores de Impacto
		}

		if (*CNT5_.NED > 0){
			endetr2_(inFILE, IWDUMP); //Detectores de deposição de energia
		}

		if (*CNT6_.LDOSEM){
			doser2_(inFILE, IWDUMP);
		}
		inFILE.close();
		fclose(IWDUMP);

		fprintf(IWR, "   Simulation has been resumed from dump file: %s\n", PFILER);
		goto L92;


L90:;
		fprintf(IWR, "   The dump file is empty or corrupted.\n");
		printf("The dump file is empty or corrupted.\n");
		exit(0);
	
L91:;
		fprintf(IWR, "   WARNING: Could not resume from dump file...\n");
		inFILE.close();
		IRESUM = 0;
	}
L92:;
	*CNT4_.IPSFO=21;
	*CSOUR4_.IPSFI=20;



 //Na sequencia é feito a leitura dos dados referentes ao arquivo Phase-Espace. Porém, esta
 //versão nao contepla a utilização desse arquivo.

 //Inicializar  Constantes


	start = clock();
	double TSECIN;
	timer2_(TSECIN);

	*CNTRL_.TSECA=TIMEA+TSECIN;
    *CNTRL_.TSECAD=TSECIN;



	*CNTRL_.SHN=SHNA;  //Contador de simulações partículas do arquivo de despejo
    *CNTRL_.N=fmod(*CNTRL_.SHN,2.0e9)+0.5e0;
    *CNTRL_.TSIM=CPUTA;
   // end = clock();
	//*CNTRL_.CPUT0= (double)(end - start) / CLOCKS_PER_SEC;
	*CNTRL_.CPUT0=cputim2_();
	if (*CNTRL_.SHN > *CNTRL_.DSHN-0.5e0){
		fprintf(IWR, "  **** The simulation was already completed. \n");
		printf("  **** The simulation was already completed.\n");
		*CSOUR0_.JOBEND=3;
	}else{
		printf("   The simulation is started ...\n");
	}



//printf("\n\nimprimndo tabelas\n\n");
//tabelas_();
	//		printf("\n\nCaracteristicas do trabalh\n\n");
//		exit(0);

}


void gcone02_(double THETA, double PHI, double ALPHA){



	/*
	Esta sub-rotina define os parâmetros para amostragem de direções aleatórias
    uniformemente dentro de um cone com eixo na direção (THETA,PHI) e
    Abertura ALPHA (em rad).
	*/
	*CGCONE_.CPCT=cos(PHI)*cos(THETA);
    *CGCONE_.CPST=cos(PHI)*sin(THETA);
    *CGCONE_.SPCT=sin(PHI)*cos(THETA);
    *CGCONE_.SPST=sin(PHI)*sin(THETA);
    *CGCONE_.SPHI=sin(PHI);
    *CGCONE_.CPHI=cos(PHI);
    *CGCONE_.STHE=sin(THETA);
    *CGCONE_.CTHE=cos(THETA);
    *CGCONE_.CAPER=cos(ALPHA);	

}

void enang02_(double &EMIN, double &EMAX, int &NBE, int &NBTH, int &NBPH, FILE *IWR){


	/*
	Calcula energia e distribuições angulares de partículas emergentes,
    grava e carrega arquivos de despejo, acumula arquivos de despejo de diferentes
    é executado e grava os resultados.
	*/

	static const int NBEM=1500;
	static const int NBTHM=1800;
	static const int NBPHM=180;
	double PI = 3.1415926535897932e0;
	double TWOPI=2.0e0*PI;
	double RA2DE=180.0e0/PI;
	double DE2RA=PI/180.0e0;
	double FSAFE=1.000000001e0; //fator de segurança

	if (abs(NBE) > NBEM){
		fprintf(IWR, "ENANG: NBE is too large.\n");
		fprintf(IWR, "ENANG: Set the parameter NBEM equal to %.d\n", abs(NBE));
		printf("ENANG: NBE is too large.\n");
		exit(0);
	}

	if (abs(NBTH) > NBTHM){
		fprintf(IWR, "ENANG: NBTH is too large.\n");
		fprintf(IWR, "ENANG: Set the parameter NBTHM equal to %.d\n", abs(NBTH));
		printf("ENANG: NBTH is too large.\n");
		exit(0);
	}

	if (NBPH > NBPHM){
		fprintf(IWR, "ENANG: NBPH is too large.\n");
		fprintf(IWR, "ENANG: Set the parameter NBPHM equal to %.d\n", abs(NBPH));
		printf("ENANG: NBPH is too large.\n");
		exit(0);
	}

	if (NBE < 0){
		*CENANG_.LLE=1;
        *CENANG_.EL=log(EMIN);
        *CENANG_.EU=log(EMAX);
        *CENANG_.NE=-NBE;
	}else{
		*CENANG_.LLE=0;
        *CENANG_.EL=EMIN;
        *CENANG_.EU=EMAX;
        *CENANG_.NE=NBE;
	}

	*CENANG_.BSE=FSAFE*(*CENANG_.EU - *CENANG_.EL) / *CENANG_.NE;
    *CENANG_.RBSE=1.0e0/ *CENANG_.BSE;

	if (NBTH < 0){
		*CENANG_.LLTH=1;
        *CENANG_.THL=log(1.0e-2);
        *CENANG_.THU=log(180.0e0);
        *CENANG_.NTH=-NBTH;
	}else{
		*CENANG_.LLTH=0;
        *CENANG_.THL=0.0e0;
        *CENANG_.THU=180.0e0;
        *CENANG_.NTH=NBTH;
	}

	*CENANG_.BSTH=FSAFE*(*CENANG_.THU-*CENANG_.THL)/ *CENANG_.NTH;
    *CENANG_.RBSTH=1.0e0/ *CENANG_.BSTH;

	if (NBPH < 0)
		*CENANG_.NPH=-NBPH;
	else
		*CENANG_.NPH=NBPH;

	*CENANG_.BSPH=FSAFE*360.0e0/ *CENANG_.NPH;
    *CENANG_.RBSPH=1.0e0/ *CENANG_.BSPH;

	

	for (int I = 1; I <= 3; I++){
		for (int J = 1; J <= 2; J++){
			for (int K = 1; K <= NBEM; K++){
				CENANG_.PDE[K-1][J-1][I-1]=0.0e0;
                CENANG_.PDE2[K-1][J-1][I-1]=0.0e0;
                CENANG_.PDEP[K-1][J-1][I-1]=0.0e0;
                CENANG_.LPDE[K-1][J-1][I-1]=0;
			}
		}
	}



	for (int I = 1; I <= 3; I++){
		for (int J = 1; J <= NBTHM; J++){
			for (int K = 1; K <= NBPHM; K++){
				CENANG_.PDA[K-1][J-1][I-1]=0.0e0;
                CENANG_.PDA2[K-1][J-1][I-1]=0.0e0;
                CENANG_.PDAP[K-1][J-1][I-1]=0.0e0;
                CENANG_.LPDA[K-1][J-1][I-1]=0;
			}
		}
	}


}

void imdet02_(double &EMIN, double &EMAX, int &NBE, double &AGEMIN, double &AGEMAX, int &NBAGE, int &ICUT, char *FNSPC, char *FNFLU, char *FNAGE, int &ID, FILE *IWR){

	/*
	Calcula espectros de detectores de impacto, grava e carrega arquivos de despejo,
	 acumula arquivos de despejo de diferentes execuções e grava os resultados.
	*/

	double FSAFE=1.000000001e0;

	static const int NIDM=25;
	static const int NBEM=1000;


	static int NIDS = 0;

	if (ID <= NIDS){
		fprintf(IWR, " SIMDET: Detector already defined.\n");
		printf("SIMDET: Detector already defined.\n");
		exit(0);
	}

	NIDS=ID;
    *CIMDET_.NID=ID;

	if (*CIMDET_.NID > NIDM){
		fprintf(IWR, "   NID = %4d\n",*CIMDET_.NID );
		fprintf(IWR, "SIMDET: Too many detectors.\n");
		printf("SIMDET: Too many detectors.\n");
		exit(0);
	}

	strcpy(SPCDIO[ID-1], FNSPC);
    strcpy(SPCFLO[ID-1],FNFLU);
    strcpy(SPCAGE[ID-1],FNAGE);
    CIMDET_.IDCUT[ID-1]=ICUT;

	if (abs(NBE) > NBEM){
		fprintf(IWR, "SIMDET: NB is too large.\n");
		fprintf(IWR, "SIMDET: Set the parameter NBEM equal to %d\n",abs(NBE) );
		printf("SIMDET: NB is too large.\n");
		exit(0);
	}

	if (NBE < 0){
		CIMDET_.LLE[ID-1]=1;
        CIMDET_.EL[ID-1]=log(EMIN);
        CIMDET_.EU[ID-1]=log(EMAX);
        CIMDET_.NE[ID-1]=-NBE;
	}else{
		CIMDET_.LLE[ID-1]=0;
        CIMDET_.EL[ID-1]=EMIN;
        CIMDET_.EU[ID-1]=EMAX;
        CIMDET_.NE[ID-1]=NBE;
	}

	CIMDET_.BSE[ID-1]=FSAFE*(CIMDET_.EU[ID-1]-CIMDET_.EL[ID-1])/(CIMDET_.NE[ID-1]);
    CIMDET_.RBSE[ID-1]=1.0e0/CIMDET_.BSE[ID-1];

	for (int J = 1; J <= CIMDET_.NE[ID-1]+1; J++){
		if (CIMDET_.LLE[ID-1] == 1)
			CIMDET_.ET[J-1][ID-1]=exp(CIMDET_.EL[ID-1]+(J-1)*CIMDET_.BSE[ID-1]);
		else
			CIMDET_.ET[J-1][ID-1]=CIMDET_.EL[ID-1]+(J-1)*CIMDET_.BSE[ID-1];
	}

	CIMDET_.EDEP[ID-1]=0.0e0;
    CIMDET_.EDEP2[ID-1]=0.0e0;
    CIMDET_.EDEPP[ID-1]=0.0e0;
    CIMDET_.LEDEP[ID-1]=0;

	for (int J = 1;J <= NBEM; J++){
		CIMDET_.DIT[J-1][ID-1]=0.0e0;
        CIMDET_.DIT2[J-1][ID-1]=0.0e0;
        CIMDET_.DITP[J-1][ID-1]=0.0e0;
        CIMDET_.LDIT[J-1][ID-1]=0;
		for (int K = 1; K <= 3; K++){
			CIMDET_.DIP[K-1][J-1][ID-1]=0.0e0;
            CIMDET_.DIP2[K-1][J-1][ID-1]=0.0e0;
            CIMDET_.DIPP[K-1][J-1][ID-1]=0.0e0;
            CIMDET_.LDIP[K-1][J-1][ID-1]=0;
		}
	}

	for (int J = 1; J <= NBEM; J++){
		CIMDET_.FLT[J-1][ID-1]=0.0e0;
        CIMDET_.FLT2[J-1][ID-1]=0.0e0;
        CIMDET_.FLTP[J-1][ID-1]=0.0e0;
        CIMDET_.LFLT[J-1][ID-1]=0;
		for (int K = 1; K <= 3; K++){
			CIMDET_.FLP[K-1][J-1][ID-1]=0.0e0;
            CIMDET_.FLP2[K-1][J-1][ID-1]=0.0e0;
            CIMDET_.FLPP[K-1][J-1][ID-1]=0.0e0;
            CIMDET_.LFLP[K-1][J-1][ID-1]=0;
		}
	}

	if (NBAGE < 0){
		CIMDET_.LLAGE[ID-1]=1;
        CIMDET_.AGEL[ID-1]=log(AGEMIN);
        CIMDET_.AGEU[ID-1]=log(AGEMAX);
        CIMDET_.NAGE[ID-1]=-NBAGE;
	}else if (NBAGE > 0){
		CIMDET_.LLAGE[ID-1]=0;
        CIMDET_.AGEL[ID-1]=AGEMIN;
        CIMDET_.AGEU[ID-1]=AGEMAX;
        CIMDET_.NAGE[ID-1]=NBAGE;
	}else{
		CIMDET_.LLAGE[ID-1]=0;
        CIMDET_.NAGE[ID-1]=0;
        CIMDET_.AGEL[ID-1]=0.0e0;
        CIMDET_.AGEU[ID-1]=0.0e0;
        CIMDET_.BAGE[ID-1]=0.0e0;
        CIMDET_.RBAGE[ID-1]=0.0e0;
	}

	if (CIMDET_.NAGE[ID-1] > 0){
		*TRACK_mod_.LAGE=true;
        CIMDET_.BAGE[ID-1]=FSAFE*(CIMDET_.AGEU[ID-1]-CIMDET_.AGEL[ID-1])/(CIMDET_.NAGE[ID-1]);
        CIMDET_.RBAGE[ID-1]=1.0e0/CIMDET_.BAGE[ID-1];
		for (int J = 1; J <= NBEM; J++){
			CIMDET_.AGE[J-1][ID-1]=0.0e0;
            CIMDET_.AGE2[J-1][ID-1]=0.0e0;
            CIMDET_.AGEP[J-1][ID-1]=0.0e0;
            CIMDET_.LAGEA[J-1][ID-1]=0;
		}
	}


}

void endet02_(double &EMIN, double &EMAX, int &NB, char *FNSPC, int &ID, FILE *IWR){

	/*Registra espectros de detectores de deposição de energia, grava e carrega
	 despeja arquivos, acumula arquivos de despejo de diferentes execuções e grava
	Resultados .*/

	double FSAFE=1.000000001e0;

	static const int NIDM=25;
	static const int NBEM=1000;


	static int NIDS = 0;

	if (ID <= NIDS){
		fprintf(IWR, " SENDET: Detector already defined.\n");
		printf("SENDET: Detector already defined.\n");
		exit(0);
	}

	NIDS=ID;
    *CENDET_.NID=ID;

	if (*CENDET_.NID > NIDM){
		fprintf(IWR, "   NID = %4d\n",*CENDET_.NID );
		fprintf(IWR, "SENDET: Too many detectors.\n");
		printf("SENDET: Too many detectors.\n");
		exit(0);
	}

	strcpy(SPCDEO[ID-1], FNSPC);

	if (abs(NB) > NBEM){
		fprintf(IWR, "SENDET: NB is too large.\n");
		fprintf(IWR, "SENDET: Set the parameter NBEM equal to %d\n",abs(NB) );
		printf("SENDET: NB is too large.\n");
		exit(0);
	}

	if (NB < 0){
		CENDET_.LLE[ID-1]=1;
        CENDET_.EL[ID-1]=log(EMIN);
        CENDET_.EU[ID-1]=log(EMAX);
        CENDET_.NE[ID-1]=-NB;
	}else{
		CENDET_.LLE[ID-1]=0;
        CENDET_.EL[ID-1]=EMIN;
        CENDET_.EU[ID-1]=EMAX;
        CENDET_.NE[ID-1]=NB;
	}

	CENDET_.BSE[ID-1]=FSAFE*(CENDET_.EU[ID-1]-CENDET_.EL[ID-1])/(CENDET_.NE[ID-1]);
    CENDET_.RBSE[ID-1]=1.0e0/CENDET_.BSE[ID-1];
    CENDET_.EDEP[ID-1]=0.0e0;
    CENDET_.EDEP2[ID-1]=0.0e0;

	for (int J = 1; J <= NBEM; J++){
		CENDET_.DET[J-1][ID-1]=0.0e0;
	}

}

void dose02_(double &XL, double &XU,double &YL,double &YU,double &ZL,double &ZU, int &NBX, int &NBY, int &NBZ, int &IDOSE, FILE *IWR){

	/*
	Registra a distribuição de dose dentro da caixa de dose, grava e carrega
 	despeja arquivos, acumula arquivos de despejo de diferentes execuções e grava
	Resultados .
	*/

	static const int NDXM = 201;
	static const int NDYM = 201;
	static const int NDZM = 201;

	double FSAFE = 1.000000001e0;
	double PI = 3.1415926535897932e0;

	int NCS = 5;
	double RNCS=1.0e0/5.0e0;

	double DX, DY, DZ, VOXEL, FNORM, DDX, DDY, DDZ, UO, VO, DENT, DR, DDR, TMASS,VOLUM;

	//Inicialização: massas médias de voxel, contadores de dose.

	if ((IDOSE < 1) || (IDOSE > 3)){
		fprintf(IWR, "(IDOSE = %6d\n", IDOSE);
		fprintf(IWR, "IDOSE should be 1, 2, or 3.\n");
		printf("SDOSE: IDOSE should be 1, 2, or 3.\n");
		exit(0);
	}

	*CDOSE1_.KDOSE=IDOSE;

	if ((NBX < 0) || (NBX > NDXM)){
		fprintf(IWR, "NBX = %6d\n", NBX);
		fprintf(IWR, "NBX must be .GT.0. and .LT. %4d\n", NDXM);
		fprintf(IWR, "Increase the value of the parameter NDXM\n");
		printf("SDOSE: NBX must be .GT.0. and .LE.NDXM");
		exit(0);
	}

	if ((NBY < 0) || (NBY > NDYM)){
		fprintf(IWR, "NBY = %6d\n", NBY);
		fprintf(IWR, "NBY must be .GT.0. and .LT. %4d\n", NDYM);
		fprintf(IWR, "Increase the value of the parameter NDYM\n");
		printf("SDOSE: NBY must be .GT.0. and .LE.NDYM");
		exit(0);
	}

	if ((NBZ < 0) || (NBZ > NDZM)){
		fprintf(IWR, "NBZ = %6d\n", NBZ);
		fprintf(IWR, "NBZ must be .GT.0. and .LT. %4d\n", NDZM);
		fprintf(IWR, "Increase the value of the parameter NDZM\n");
		printf("SDOSE: NBZ must be .GT.0. and .LE.NDZM");
		exit(0);
	}

	CDOSE3_.DXL[1-1]=XL;
    CDOSE3_.DXU[1-1]=XU;
    CDOSE3_.DXL[2-1]=YL;
    CDOSE3_.DXU[2-1]=YU;
    CDOSE3_.DXL[3-1]=ZL;
    CDOSE3_.DXU[3-1]=ZU;

	if (*CDOSE1_.KDOSE == 1){
		CDOSE3_.NDB[1-1]=2*(NBX/2)+1;
        CDOSE3_.NDB[2-1]=2*(NBY/2)+1;
        CDOSE3_.NDB[3-1]=2*(NBZ/2)+1;
	}else if (*CDOSE1_.KDOSE == 2){
		CDOSE3_.DXL[1-1]=0.0e0;
		CDOSE3_.DXL[2-1]=0.0e0;
        CDOSE3_.DXU[2-1]=1.0e0;
        CDOSE3_.NDB[1-1]=NBX;
        CDOSE3_.NDB[2-1]=1;
        CDOSE3_.NDB[3-1]=NBZ;
	}else{
		CDOSE3_.DXL[1-1]=0.0e0;
        CDOSE3_.DXL[2-1]=0.0e0;
        CDOSE3_.DXU[2-1]=1.0e0;
        CDOSE3_.DXL[3-1]=0.0e0;
        CDOSE3_.DXU[3-1]=1.0e0;
        CDOSE3_.NDB[1-1]=NBX;
        CDOSE3_.NDB[2-1]=1;
        CDOSE3_.NDB[3-1]=1;
	}

	for (int I = 1; I <= 3; I++){
		CDOSE3_.BDOSE[I-1]=FSAFE*(CDOSE3_.DXU[I-1]-CDOSE3_.DXL[I-1])/(CDOSE3_.NDB[I-1]);
		if (fabs(CDOSE3_.BDOSE[I-1]) < 1.0e-35)
			CDOSE3_.RBDOSE[I-1]=1.0e35;
		else
		    CDOSE3_.RBDOSE[I-1]=1.0e0/CDOSE3_.BDOSE[I-1];
	}

	//Contadores de Dose

	for (int K = 1; K <= NDZM; K++){
		for (int J = 1; J <= NDYM; J++){
			for (int I =1; I <= NDXM; I++){
				CDOSE1_.DOSE[K-1][J-1][I-1]=0.0e0;
                CDOSE1_.DOSE2[K-1][J-1][I-1]=0.0e0;
                CDOSE1_.DOSEP[K-1][J-1][I-1]=0.0e0;
                CDOSE1_.LDOSE[K-1][J-1][I-1]=0;
                CDOSE4_.VMASS[K-1][J-1][I-1]=0.0e0;
			}
		}
		CDOSE2_.DDOSE[K-1]=0.0e0;
        CDOSE2_.DDOSE2[K-1]=0.0e0;
        CDOSE2_.DDOSEP[K-1]=0.0e0;
        CDOSE2_.LDDOSE[K-1]=0;
	}

	//Massas de Voxels

	//Dose em uma Caixa

	if (*CDOSE1_.KDOSE == 1){
		DX =CDOSE3_.BDOSE[1-1];
		DY=CDOSE3_.BDOSE[2-1];
		DZ=CDOSE3_.BDOSE[3-1];
        VOXEL =DX*DY*DZ;
        FNORM=VOXEL*pow(RNCS,3);
        DDX=RNCS*DX; 
		DDY=RNCS*DY; 
		DDZ=RNCS*DZ;
		for (int K = 1; K <= CDOSE3_.NDB[3-1]; K++){
			for (int J = 1; J <= CDOSE3_.NDB[2-1]; J++){
				for (int I = 1; I <= CDOSE3_.NDB[1-1]; I++){
					UO=1.0e0;
					VO=1.0e0;
                	*TRACK_mod_.X=CDOSE3_.DXL[1-1]+(I-1)*DX+0.5e0*DDX;
                	*TRACK_mod_.Y=CDOSE3_.DXL[2-1]+(J-1)*DY+0.5e0*DDY;
               	 	*TRACK_mod_.Z=CDOSE3_.DXL[3-1]+(K-1)*DZ+0.5e0*DDZ;
					locate2_();
					if (*TRACK_mod_.MAT == 0)
						DENT=0.0e0;
					else	
						DENT=PENELOPE_mod_.DEN[*TRACK_mod_.MAT-1];

					for (int J3 = 1; J3 <= NCS; J3++){
						for (int J2 = 1; J2 <= NCS; J2++){
							for (int J1 = 1; J1 <= NCS-1; J1++){
								*TRACK_mod_.U=UO;
								*TRACK_mod_.V=0.0e0;
								*TRACK_mod_.W=0.0e0;
                      			*TRACK_mod_.X= *TRACK_mod_.X + DDX * *TRACK_mod_.U;
							    locate2_();
								if (*TRACK_mod_.MAT > 0)
									DENT=DENT+PENELOPE_mod_.DEN[*TRACK_mod_.MAT-1];
							}
							UO=-UO;
							if (J2 < NCS){
								*TRACK_mod_.U=0.0e0;
								*TRACK_mod_.V=VO;
								*TRACK_mod_.W=0.0e0;
                      			*TRACK_mod_.Y=*TRACK_mod_.Y + DDY * *TRACK_mod_.V;
								locate2_();
								if (*TRACK_mod_.MAT > 0)
									DENT=DENT+PENELOPE_mod_.DEN[*TRACK_mod_.MAT-1];
							}
						}
						VO=-VO;
						if (J3 < NCS){
							*TRACK_mod_.U=0.0e0; 
							*TRACK_mod_.V=0.0e0; 
							*TRACK_mod_.W=1.0e0;
                    		*TRACK_mod_.Z=*TRACK_mod_.Z + DDZ * *TRACK_mod_.W;
							locate2_();
							if (*TRACK_mod_.MAT > 0)
								DENT=DENT+PENELOPE_mod_.DEN[*TRACK_mod_.MAT-1];
						}

					}
					CDOSE4_.VMASS[K-1][J-1][I-1]=DENT*FNORM;
				}
			}
		}

		for (int K = 1; K <= CDOSE3_.NDB[3-1]; K++){
			for (int J = 1; J <= CDOSE3_.NDB[2-1]; J++){
				for (int I = 1; I <= CDOSE3_.NDB[1-1]; I++){
					if (CDOSE4_.VMASS[K-1][J-1][I-1] > 1.0e-35)
						CDOSE4_.VMASS[K-1][J-1][I-1]=1.0e0/CDOSE4_.VMASS[K-1][J-1][I-1];
					else
						CDOSE4_.VMASS[K-1][J-1][I-1]=0.0e0;
				}
			}
		}


	} else if (*CDOSE1_.KDOSE == 2){ //Dose em um Cilindro
		int J = 1;
		DR=CDOSE3_.BDOSE[1-1];
		DZ=CDOSE3_.BDOSE[3-1];
        DDR=RNCS*DR; 
		DDZ=RNCS*DZ;
		*TRACK_mod_.U=1.0e0; 
		*TRACK_mod_.V=0.0e0; 
		*TRACK_mod_.W=0.0e0;
		for (int K = 1; K <= CDOSE3_.NDB[3-1]; K++){
			for (int I = 1; I <= CDOSE3_.NDB[1-1]; I++){
				TMASS=0.0e0;
				for (int J3 = 1; J3 <= NCS; J3++){
					for (int J1 = 1; J1 <= NCS; J1++){
						*TRACK_mod_.X=(I-1)*DR+((J1)-0.5e0)*DDR;
                  		*TRACK_mod_.Y=0.0e0;
                  		*TRACK_mod_.Z=CDOSE3_.DXL[3-1]+(K-1)*DZ+((J3)-0.5e0)*DDZ;
					    locate2_();
						if (*TRACK_mod_.MAT != 0){
							VOLUM=2.0e0*PI* *TRACK_mod_.X *DDR*DDZ;
                    		TMASS=TMASS+PENELOPE_mod_.DEN[*TRACK_mod_.MAT-1]*VOLUM;
						}
					}
				}
				CDOSE4_.VMASS[K-1][J-1][I-1]=TMASS;
				if (CDOSE4_.VMASS[K-1][J-1][I-1] > 1.0e-35)
					CDOSE4_.VMASS[K-1][J-1][I-1]=1.0e0/CDOSE4_.VMASS[K-1][J-1][I-1];
				else
					CDOSE4_.VMASS[K-1][J-1][I-1]=0.0e0;
			}
		}
	}else{//Dose em uma esfera.
		int J = 1;
		int K = 1;
		DR=CDOSE3_.BDOSE[1-1];
        DDR=RNCS*DR;
        *TRACK_mod_.U=0.0e0; 
		*TRACK_mod_.V=0.0e0; 
		*TRACK_mod_.W=1.0e0;
		for (int I = 1; I <= CDOSE3_.NDB[1-1]; I++){
			TMASS=0.0e0;
			for (int J1 = 1; J1 <= NCS; J1++){
				*TRACK_mod_.X=(I-1)*DR+((J1)-0.5e0)*DDR;
                *TRACK_mod_.Y=0.0e0;
                *TRACK_mod_.Z=0.0e0;
			    locate2_();
				if (*TRACK_mod_.MAT != 0){
					VOLUM=3.0e0*pow(*TRACK_mod_.X,2)+0.25e0*pow(DDR,2);
                    TMASS=TMASS+PENELOPE_mod_.DEN[*TRACK_mod_.MAT-1]*VOLUM;
				}
			}
			CDOSE4_.VMASS[K-1][J-1][I-1]=TMASS*(4.0e0*PI/3.0e0)*DDR;
			if (CDOSE4_.VMASS[K-1][J-1][I-1] > 1.0e-35)
					CDOSE4_.VMASS[K-1][J-1][I-1]=1.0e0/CDOSE4_.VMASS[K-1][J-1][I-1];
				else
					CDOSE4_.VMASS[K-1][J-1][I-1]=0.0e0;
		}
	} 


//	printf("\n\ndose02\n\n");
}

void rand02_(int N){

	/*
	Em cálculos paralelos, precisamos inicializar RAND() de tal forma
  que diferentes processadores produzem sequências verdadeiramente independentes de
  Números aleatórios. Isso pode ser feito alimentando cada processador
  com sementes iniciais que pertencem a uma única longa sequência gerada por
  RAND(), mas estão distantes o suficiente um do outro para evitar a sobreposição do
  subsequências geradas pelos diferentes processadores. A lista abaixo
  foi obtido executando um programa escrito por Andreu Badal e Josep
  Sempau [Comp. Física Comum. 175 (2006) 440-450]. Contém pares de
  sementes que pertencem a uma longa sequência e cuja separação é 10**14
  chamadas. Ou seja, se começarmos com o N-ésimo par de sementes, após 10**14
  chamadas para RAND() obtemos o par (N+1)-ésimo.

  Uma chamada para a sub-rotina presente RAND0(N) carrega o N-ésimo par de
  sementes na lista. Assim, em simulações paralelas, podemos inicializar
  RAND() nos diferentes processadores chamando RAND0 com diferentes
  valores do argumento N.
	*/

	int IS1[1002];
	int IS2[1002];

	  IS1[1   ] =          1; IS2[1   ] =          1;
      IS1[2   ] = 1088794366; IS2[2   ] =  722792456;
      IS1[3   ] = 1993751964; IS2[3   ] =  753089694;
      IS1[4   ] =  610005387; IS2[4   ] = 1748134360;
      IS1[5   ] =   27944595; IS2[5   ] =  774572312;
      IS1[6   ] = 1108394934; IS2[6   ] =  620441713;
      IS1[7   ] =  580582953; IS2[7   ] =   73753031;
      IS1[8   ] =  479095664; IS2[8   ] =  801971873;
      IS1[9   ] =  204282520; IS2[9   ] = 1624172377;
      IS1[10  ] =  327796803; IS2[10  ] =  915478773;
      IS1[11  ] =  918882992; IS2[11  ] =  858672133;
      IS1[12  ] =  210516912; IS2[12  ] =  158496513;
      IS1[13  ] = 1155963390; IS2[13  ] =  704005330;
      IS1[14  ] = 1929580647; IS2[14  ] =  324681933;
      IS1[15  ] = 2061216471; IS2[15  ] =  209323758;
      IS1[16  ] = 1963942692; IS2[16  ] = 1686825906;
      IS1[17  ] = 2124896501; IS2[17  ] =  460463303;
      IS1[18  ] =  233299547; IS2[18  ] = 2134948096;
      IS1[19  ] = 1675549170; IS2[19  ] = 1144329343;
      IS1[20  ] = 1421385813; IS2[20  ] =  672532778;
      IS1[21  ] = 2069007070; IS2[21  ] = 1309916099;
      IS1[22  ] = 1905029285; IS2[22  ] = 1817420725;
      IS1[23  ] = 1056478052; IS2[23  ] = 1943244114;
      IS1[24  ] = 1628262065; IS2[24  ] =  782316895;
      IS1[25  ] = 1368245324; IS2[25  ] = 2070604605;
      IS1[26  ] =  414218397; IS2[26  ] =  519652740;
      IS1[27  ] = 1032135553; IS2[27  ] = 1729805542;
      IS1[28  ] = 1709285619; IS2[28  ] = 1859376254;
      IS1[29  ] = 1363836736; IS2[29  ] =  433738184;
      IS1[30  ] = 1307886537; IS2[30  ] = 1841770585;
      IS1[31  ] =  944675654; IS2[31  ] = 1438406465;
      IS1[32  ] =  733913518; IS2[32  ] =  219173860;
      IS1[33  ] = 1994697060; IS2[33  ] = 1473267531;
      IS1[34  ] = 1816284687; IS2[34  ] = 1806014371;
      IS1[35  ] = 1204077982; IS2[35  ] =  313513580;
      IS1[36  ] = 1322304106; IS2[36  ] = 1896720184;
      IS1[37  ] = 1152681299; IS2[37  ] = 1347076753;
      IS1[38  ] = 1910575054; IS2[38  ] = 1881544152;
      IS1[39  ] = 1158165315; IS2[39  ] = 1769280665;
      IS1[40  ] = 1256652686; IS2[40  ] = 1678458141;
      IS1[41  ] =  149156960; IS2[41  ] =  257442270;
      IS1[42  ] = 1531563057; IS2[42  ] =  769379539;
      IS1[43  ] =   33227674; IS2[43  ] = 1940657859;
      IS1[44  ] =  563899375; IS2[44  ] = 1495802342;
      IS1[45  ] =  534312669; IS2[45  ] = 1627605675;
      IS1[46  ] =   20883767; IS2[46  ] = 1860146428;
      IS1[47  ] =  654784964; IS2[47  ] = 1449089950;
      IS1[48  ] = 1622565819; IS2[48  ] =  300004830;
      IS1[49  ] =  720217438; IS2[49  ] =   59015257;
      IS1[50  ] = 1162709863; IS2[50  ] = 1966621140;
      IS1[51  ] =  360537627; IS2[51  ] =  133123709;
      IS1[52  ] =  345935498; IS2[52  ] = 1389524816;
      IS1[53  ] =  876665055; IS2[53  ] = 1171092501;
      IS1[54  ] = 1998564166; IS2[54  ] = 1037986676;
      IS1[55  ] = 2066661466; IS2[55  ] =  366166383;
      IS1[56  ] = 1279875030; IS2[56  ] =  582766249;
      IS1[57  ] = 1955021266; IS2[57  ] =  116351171;
      IS1[58  ] =  195816083; IS2[58  ] =  933489046;
      IS1[59  ] = 2055228393; IS2[59  ] =  283595456;
      IS1[60  ] =   93961337; IS2[60  ] =  611643703;
      IS1[61  ] = 1446789139; IS2[61  ] = 1248992867;
      IS1[62  ] =  414524648; IS2[62  ] =  803732805;
      IS1[63  ] = 1327723283; IS2[63  ] =  201031458;
      IS1[64  ] =  535511919; IS2[64  ] =  938162313;
      IS1[65  ] =   86014814; IS2[65  ] =    7981320;
      IS1[66  ] = 1646051327; IS2[66  ] = 1690586644;
      IS1[67  ] = 1717255367; IS2[67  ] =  569003207;
      IS1[68  ] = 1794770573; IS2[68  ] = 1043545502;
      IS1[69  ] =  853934312; IS2[69  ] =  848156009;
      IS1[70  ] =  903676328; IS2[70  ] = 1753511888;
      IS1[71  ] =  888673974; IS2[71  ] = 2014364429;
      IS1[72  ] = 1763908015; IS2[72  ] =  835386162;
      IS1[73  ] =  415969825; IS2[73  ] =  751041626;
      IS1[74  ] = 2041881831; IS2[74  ] =  527275421;
      IS1[75  ] = 1223957624; IS2[75  ] =  282891294;
      IS1[76  ] =  216252620; IS2[76  ] = 1933221826;
      IS1[77  ] =  956133864; IS2[77  ] = 1506952100;
      IS1[78  ] = 1265794918; IS2[78  ] =  713084270;
      IS1[79  ] = 1474654961; IS2[79  ] =   88797131;
      IS1[80  ] =  122865758; IS2[80  ] = 2073278963;
      IS1[81  ] =     258943; IS2[81  ] =  664687714;
      IS1[82  ] = 1152463120; IS2[82  ] =  839536051;
      IS1[83  ] =  181367474; IS2[83  ] = 1038508765;
      IS1[84  ] =  618933039; IS2[84  ] =  131404490;
      IS1[85  ] = 1185139338; IS2[85  ] = 1987460903;
      IS1[86  ] = 2078683375; IS2[86  ] = 1391654211;
      IS1[87  ] = 1157287301; IS2[87  ] = 1870939725;
      IS1[88  ] =  490572205; IS2[88  ] = 1828672152;
      IS1[89  ] =  713452544; IS2[89  ] = 1880501392;
      IS1[90  ] = 1399731654; IS2[90  ] =  661442337;
      IS1[91  ] = 1434784182; IS2[91  ] = 1598489021;
      IS1[92  ] =  157980824; IS2[92  ] =  512999623;
      IS1[93  ] = 1631668015; IS2[93  ] = 2021483107;
      IS1[94  ] =  695840037; IS2[94  ] =  316143310;
      IS1[95  ] = 2011902133; IS2[95  ] =  460681770;
      IS1[96  ] = 1482456963; IS2[96  ] = 1832621179;
      IS1[97  ] =  983630146; IS2[97  ] = 1244860854;
      IS1[98  ] =  424488068; IS2[98  ] = 1894080738;
      IS1[99  ] =  592109479; IS2[99  ] =  254844002;
      IS1[100 ] =  698713089; IS2[100 ] =  852473215;
      IS1[101 ] =  698429770; IS2[101 ] = 1978724894;
      IS1[102 ] = 1991339714; IS2[102 ] =  590904772;
      IS1[103 ] = 1812612029; IS2[103 ] = 1455273894;
      IS1[104 ] =  878555690; IS2[104 ] =  849354664;
      IS1[105 ] = 1415741666; IS2[105 ] = 1989849407;
      IS1[106 ] =  740336773; IS2[106 ] = 2094756349;
      IS1[107 ] = 1357150877; IS2[107 ] = 1752614313;
      IS1[108 ] =  446288602; IS2[108 ] =  605474927;
      IS1[109 ] =  156511135; IS2[109 ] = 1112863746;
      IS1[110 ] = 1315731039; IS2[110 ] = 1420661365;
      IS1[111 ] = 1590179518; IS2[111 ] =  797341276;
      IS1[112 ] =  210183789; IS2[112 ] = 1297307595;
      IS1[113 ] =   94234820; IS2[113 ] = 1379679575;
      IS1[114 ] =  273023900; IS2[114 ] = 1595938728;
      IS1[115 ] =  868831745; IS2[115 ] = 1149726246;
      IS1[116 ] =  170384549; IS2[116 ] =  129431617;
      IS1[117 ] =   13185587; IS2[117 ] =  487985593;
      IS1[118 ] = 1362898234; IS2[118 ] =  342281893;
      IS1[119 ] =  574105532; IS2[119 ] =  245617051;
      IS1[120 ] = 1822103760; IS2[120 ] = 1034640984;
      IS1[121 ] =  658812725; IS2[121 ] = 1829379549;
      IS1[122 ] = 1005671726; IS2[122 ] =  183096918;
      IS1[123 ] = 1254435204; IS2[123 ] = 1402468728;
      IS1[124 ] = 1798478003; IS2[124 ] =  432519190;
      IS1[125 ] =  601935466; IS2[125 ] =  253536637;
      IS1[126 ] =  341666647; IS2[126 ] =  118329947;
      IS1[127 ] = 1713183513; IS2[127 ] = 2027318311;
      IS1[128 ] = 1196735493; IS2[128 ] = 1530795526;
      IS1[129 ] = 1348986425; IS2[129 ] =  668877863;
      IS1[130 ] =  532005772; IS2[130 ] =  203741901;
      IS1[131 ] = 1141813962; IS2[131 ] = 1863091192;
      IS1[132 ] = 1741278758; IS2[132 ] =  104388873;
      IS1[133 ] =   86859261; IS2[133 ] = 2120737529;
      IS1[134 ] =  726923420; IS2[134 ] = 1212110048;
      IS1[135 ] = 1058541527; IS2[135 ] =  166531451;
      IS1[136 ] = 2131549752; IS2[136 ] = 1447596029;
      IS1[137 ] = 2033958033; IS2[137 ] =  926546635;
      IS1[138 ] =  969574276; IS2[138 ] = 1513067587;
      IS1[139 ] = 1924153065; IS2[139 ] = 1224270071;
      IS1[140 ] = 1789201564; IS2[140 ] =  244925517;
      IS1[141 ] =  970004038; IS2[141 ] =  827424326;
      IS1[142 ] =  584997635; IS2[142 ] =  900020443;
      IS1[143 ] = 1121567821; IS2[143 ] = 1551288142;
      IS1[144 ] = 1901695144; IS2[144 ] =  802507892;
      IS1[145 ] = 1381949729; IS2[145 ] =  338664653;
      IS1[146 ] =  931183121; IS2[146 ] =  255723333;
      IS1[147 ] =  425040923; IS2[147 ] = 1183865313;
      IS1[148 ] = 2063648383; IS2[148 ] = 1931656776;
      IS1[149 ] =  372520795; IS2[149 ] = 1381463141;
      IS1[150 ] =  112959009; IS2[150 ] =  328713331;
      IS1[151 ] =  156172654; IS2[151 ] = 2028805919;
      IS1[152 ] = 1206630112; IS2[152 ] = 1317701868;
      IS1[153 ] = 1141325584; IS2[153 ] = 1635316096;
      IS1[154 ] = 1226401966; IS2[154 ] =  891065751;
      IS1[155 ] = 1044869640; IS2[155 ] = 1695983251;
      IS1[156 ] = 1468275435; IS2[156 ] =  827674970;
      IS1[157 ] =  791167482; IS2[157 ] =  645339068;
      IS1[158 ] =  464714547; IS2[158 ] = 1616229175;
      IS1[159 ] = 1084778504; IS2[159 ] =  563204166;
      IS1[160 ] =  258987949; IS2[160 ] = 1152030385;
      IS1[161 ] =  479486796; IS2[161 ] =  238097920;
      IS1[162 ] = 1499316991; IS2[162 ] = 1401468685;
      IS1[163 ] =  519780735; IS2[163 ] =  481196391;
      IS1[164 ] = 1679611087; IS2[164 ] =  604709095;
      IS1[165 ] = 1691996345; IS2[165 ] =  989109993;
      IS1[166 ] =  196557364; IS2[166 ] =  113178392;
      IS1[167 ] =  851784008; IS2[167 ] = 1222308139;
      IS1[168 ] =  143277176; IS2[168 ] = 1931319584;
      IS1[169 ] = 1973626791; IS2[169 ] = 1586075498;
      IS1[170 ] =  585402145; IS2[170 ] =  125579035;
      IS1[171 ] = 1926622811; IS2[171 ] = 1383430112;
      IS1[172 ] =  614117445; IS2[172 ] =  809143743;
      IS1[173 ] =  209722147; IS2[173 ] = 1485956780;
      IS1[174 ] =  445830939; IS2[174 ] = 1964280618;
      IS1[175 ] = 1346542997; IS2[175 ] =  748234912;
      IS1[176 ] =  576549607; IS2[176 ] =  554392763;
      IS1[177 ] = 1852906063; IS2[177 ] =  675586306;
      IS1[178 ] = 1415103532; IS2[178 ] =  839112213;
      IS1[179 ] =  868356749; IS2[179 ] = 1226343583;
      IS1[180 ] = 1373222177; IS2[180 ] = 1149205884;
      IS1[181 ] = 1024191502; IS2[181 ] =  938910203;
      IS1[182 ] =  463832557; IS2[182 ] =  441736082;
      IS1[183 ] =  599161815; IS2[183 ] =  880842978;
      IS1[184 ] =   27305702; IS2[184 ] = 1259875018;
      IS1[185 ] = 1622662871; IS2[185 ] = 1545669937;
      IS1[186 ] = 1314825492; IS2[186 ] = 1230485856;
      IS1[187 ] = 1919359439; IS2[187 ] =  649279764;
      IS1[188 ] =  776491967; IS2[188 ] =  940088497;
      IS1[189 ] =  153029369; IS2[189 ] =  604610332;
      IS1[190 ] = 1075161827; IS2[190 ] =  333444224;
      IS1[191 ] =  914109165; IS2[191 ] =  803254235;
      IS1[192 ] = 1944937918; IS2[192 ] = 1451340862;
      IS1[193 ] = 1738752454; IS2[193 ] =  499708706;
      IS1[194 ] =  321360177; IS2[194 ] = 1192481545;
      IS1[195 ] = 1078034966; IS2[195 ] =  318667790;
      IS1[196 ] = 1983953435; IS2[196 ] = 1865542330;
      IS1[197 ] = 1964072934; IS2[197 ] = 1092645796;
      IS1[198 ] = 1951113931; IS2[198 ] = 1549447078;
      IS1[199 ] = 1059124001; IS2[199 ] = 1047835649;
      IS1[200 ] = 1695266076; IS2[200 ] =  163616804;
      IS1[201 ] =  996688518; IS2[201 ] = 2093478795;
      IS1[202 ] = 1354445549; IS2[202 ] =  335951295;
      IS1[203 ] =  998773555; IS2[203 ] = 1219483032;
      IS1[204 ] =    4058649; IS2[204 ] =  559033528;
      IS1[205 ] =  470886335; IS2[205 ] = 1667952318;
      IS1[206 ] =  158745647; IS2[206 ] = 2012602568;
      IS1[207 ] =  461380034; IS2[207 ] =  799155945;
      IS1[208 ] =  105655873; IS2[208 ] = 1187620434;
      IS1[209 ] =  915491195; IS2[209 ] =  858595438;
      IS1[210 ] =  247536109; IS2[210 ] =  727545379;
      IS1[211 ] = 1343291889; IS2[211 ] = 1203423504;
      IS1[212 ] = 2001578188; IS2[212 ] =  990991500;
      IS1[213 ] = 1002265824; IS2[213 ] =  219121455;
      IS1[214 ] =  501880984; IS2[214 ] =  846802413;
      IS1[215 ] = 1411611767; IS2[215 ] =  443023724;
      IS1[216 ] = 1514586313; IS2[216 ] =  974179120;
      IS1[217 ] =  261179117; IS2[217 ] =   18434323;
      IS1[218 ] =   20916375; IS2[218 ] = 2145251247;
      IS1[219 ] = 1863207976; IS2[219 ] =  176589398;
      IS1[220 ] =  969221669; IS2[220 ] =  683394530;
      IS1[221 ] = 1783116228; IS2[221 ] =  282077422;
      IS1[222 ] = 1847032911; IS2[222 ] = 1518960264;
      IS1[223 ] =  754180209; IS2[223 ] =  423392320;
      IS1[224 ] = 1927259077; IS2[224 ] =  249198027;
      IS1[225 ] =  285155942; IS2[225 ] = 1861278512;
      IS1[226 ] = 1748963900; IS2[226 ] = 1181877087;
      IS1[227 ] =  691637076; IS2[227 ] = 1593810530;
      IS1[228 ] =  378323627; IS2[228 ] =  196639057;
      IS1[229 ] =  985536851; IS2[229 ] =  991565678;
      IS1[230 ] = 1218539427; IS2[230 ] = 1989132277;
      IS1[231 ] =  552569869; IS2[231 ] = 1861318300;
      IS1[232 ] =  300092931; IS2[232 ] =  550437007;
      IS1[233 ] = 1582757652; IS2[233 ] = 1690689155;
      IS1[234 ] =  482469061; IS2[234 ] =  126744526;
      IS1[235 ] =  225088150; IS2[235 ] = 1140456485;
      IS1[236 ] =  566037433; IS2[236 ] =  654854217;
      IS1[237 ] = 1644348452; IS2[237 ] =  737109295;
      IS1[238 ] =  337491116; IS2[238 ] = 2098977589;
      IS1[239 ] =  655153746; IS2[239 ] =   43299124;
      IS1[240 ] = 1499772543; IS2[240 ] = 1730234211;
      IS1[241 ] =  839738220; IS2[241 ] = 1673889598;
      IS1[242 ] =  528696432; IS2[242 ] = 1902853997;
      IS1[243 ] =   50552080; IS2[243 ] = 2075206178;
      IS1[244 ] =  130137340; IS2[244 ] = 1283599409;
      IS1[245 ] = 1013040075; IS2[245 ] =  676361106;
      IS1[246 ] = 1793063152; IS2[246 ] = 1860713192;
      IS1[247 ] = 1912078503; IS2[247 ] =  259429094;
      IS1[248 ] =  695750580; IS2[248 ] = 1364366801;
      IS1[249 ] =  851302736; IS2[249 ] = 1708079864;
      IS1[250 ] = 1365371254; IS2[250 ] =  227135534;
      IS1[251 ] = 1436544680; IS2[251 ] =  196014388;
      IS1[252 ] =  255720485; IS2[252 ] = 1188024965;
      IS1[253 ] =  517569477; IS2[253 ] =   63939330;
      IS1[254 ] =  379457723; IS2[254 ] =  417116556;
      IS1[255 ] = 1714565676; IS2[255 ] =  870984368;
      IS1[256 ] =  427585641; IS2[256 ] =  200189082;
      IS1[257 ] = 1982058823; IS2[257 ] = 1003464933;
      IS1[258 ] = 1738564860; IS2[258 ] = 1247323567;
      IS1[259 ] =  487708829; IS2[259 ] =  462209958;
      IS1[260 ] = 1646850245; IS2[260 ] =   61645060;
      IS1[261 ] = 1590006138; IS2[261 ] = 1754830438;
      IS1[262 ] = 1733095787; IS2[262 ] = 1907130822;
      IS1[263 ] = 1660369047; IS2[263 ] = 1665129257;
      IS1[264 ] = 1104503590; IS2[264 ] = 2005102976;
      IS1[265 ] =  557868773; IS2[265 ] = 1956949606;
      IS1[266 ] =  665813373; IS2[266 ] = 2050201793;
      IS1[267 ] = 1832810072; IS2[267 ] =  206523161;
      IS1[268 ] =  420890754; IS2[268 ] = 1367078059;
      IS1[269 ] =  517192491; IS2[269 ] = 1150925658;
      IS1[270 ] = 1564894415; IS2[270 ] = 2058874982;
      IS1[271 ] =   88748046; IS2[271 ] = 1574175136;
      IS1[272 ] =  214105914; IS2[272 ] = 1747446780;
      IS1[273 ] =  667844668; IS2[273 ] =  188322609;
      IS1[274 ] = 1127730224; IS2[274 ] =  995741669;
      IS1[275 ] =  367330131; IS2[275 ] =  816142314;
      IS1[276 ] = 1246444093; IS2[276 ] =  549605711;
      IS1[277 ] =  847974161; IS2[277 ] =  183325985;
      IS1[278 ] = 1721669301; IS2[278 ] =  479407779;
      IS1[279 ] = 1915327339; IS2[279 ] = 1088013018;
      IS1[280 ] = 1578470586; IS2[280 ] = 1065050626;
      IS1[281 ] =  940315134; IS2[281 ] = 1848214072;
      IS1[282 ] = 1473252673; IS2[282 ] =    2877464;
      IS1[283 ] = 1827676712; IS2[283 ] = 1664447670;
      IS1[284 ] =  839285700; IS2[284 ] = 1640026298;
      IS1[285 ] =  751020328; IS2[285 ] = 1588681006;
      IS1[286 ] = 1091191738; IS2[286 ] = 1790306835;
      IS1[287 ] =  177083683; IS2[287 ] = 1254397576;
      IS1[288 ] = 2055613182; IS2[288 ] =  939654007;
      IS1[289 ] = 1473470878; IS2[289 ] =  335189253;
      IS1[290 ] = 1800767926; IS2[290 ] =  290320395;
      IS1[291 ] = 1676463528; IS2[291 ] = 1603536943;
      IS1[292 ] = 1650288797; IS2[292 ] =  388374267;
      IS1[293 ] = 2035708356; IS2[293 ] =  518387205;
      IS1[294 ] = 1452341404; IS2[294 ] =  985322233;
      IS1[295 ] = 1661040488; IS2[295 ] =   68406361;
      IS1[296 ] =  895503595; IS2[296 ] = 1876688389;
      IS1[297 ] = 1951727335; IS2[297 ] =  185661402;
      IS1[298 ] =  195345739; IS2[298 ] = 1532771577;
      IS1[299 ] =  268247973; IS2[299 ] = 1395541411;
      IS1[300 ] = 1057595183; IS2[300 ] =  128171866;
      IS1[301 ] = 1466997063; IS2[301 ] = 1372373334;
      IS1[302 ] = 1588823091; IS2[302 ] = 1473252323;
      IS1[303 ] = 1825606993; IS2[303 ] =  398379605;
      IS1[304 ] =  836937167; IS2[304 ] = 2001897494;
      IS1[305 ] =   60259114; IS2[305 ] = 2113852924;
      IS1[306 ] = 1298040193; IS2[306 ] =  597799372;
      IS1[307 ] =  892667157; IS2[307 ] =   98544655;
      IS1[308 ] =  349769327; IS2[308 ] = 1119519495;
      IS1[309 ] = 1659599388; IS2[309 ] =  848436478;
      IS1[310 ] =  347450508; IS2[310 ] =  197988152;
      IS1[311 ] = 1564752903; IS2[311 ] =  303075472;
      IS1[312 ] =  271104778; IS2[312 ] = 2101408514;
      IS1[313 ] = 1071650412; IS2[313 ] =  557206316;
      IS1[314 ] = 2129784998; IS2[314 ] =  319291050;
      IS1[315 ] =  149442067; IS2[315 ] = 1161643665;
      IS1[316 ] = 1382871443; IS2[316 ] =  623775604;
      IS1[317 ] =  217752411; IS2[317 ] =  888310836;
      IS1[318 ] = 1265937666; IS2[318 ] =  276467372;
      IS1[319 ] =  569940604; IS2[319 ] = 1477486653;
      IS1[320 ] =  379102003; IS2[320 ] =  377940861;
      IS1[321 ] =  885729895; IS2[321 ] = 1974466429;
      IS1[322 ] =  163050126; IS2[322 ] =   23068033;
      IS1[323 ] = 1194934955; IS2[323 ] = 1962429405;
      IS1[324 ] = 1988644587; IS2[324 ] = 1176217709;
      IS1[325 ] =  223953149; IS2[325 ] =  165457549;
      IS1[326 ] =  173062795; IS2[326 ] = 1058081267;
      IS1[327 ] = 1975967738; IS2[327 ] =  793527446;
      IS1[328 ] =  903886181; IS2[328 ] =   72550470;
      IS1[329 ] = 1844109661; IS2[329 ] = 1278970903;
      IS1[330 ] =  224746454; IS2[330 ] = 2027692323;
      IS1[331 ] =  381257506; IS2[331 ] =  782649282;
      IS1[332 ] = 1521365813; IS2[332 ] = 1328897351;
      IS1[333 ] =  838883107; IS2[333 ] =  636334681;
      IS1[334 ] = 1958225287; IS2[334 ] = 1397155038;
      IS1[335 ] = 1880542285; IS2[335 ] = 2136705686;
      IS1[336 ] =  351096309; IS2[336 ] =   67624347;
      IS1[337 ] =  276181341; IS2[337 ] =  720002598;
      IS1[338 ] = 1483813475; IS2[338 ] = 1496522845;
      IS1[339 ] = 1721418406; IS2[339 ] =  298856548;
      IS1[340 ] = 1646984747; IS2[340 ] =  613517180;
      IS1[341 ] = 1115726648; IS2[341 ] = 1979141747;
      IS1[342 ] =  497311653; IS2[342 ] =  431235843;
      IS1[343 ] = 1847923343; IS2[343 ] = 1440287460;
      IS1[344 ] = 1612185030; IS2[344 ] = 1770007478;
      IS1[345 ] = 1436854618; IS2[345 ] = 2062850297;
      IS1[346 ] = 1289356410; IS2[346 ] =  773503574;
      IS1[347 ] = 1618534459; IS2[347 ] =  405022273;
      IS1[348 ] = 1923690772; IS2[348 ] = 1950961468;
      IS1[349 ] = 2005241207; IS2[349 ] =  854563799;
      IS1[350 ] = 1191066623; IS2[350 ] =  808720040;
      IS1[351 ] =  463667715; IS2[351 ] =  466536804;
      IS1[352 ] = 1792026494; IS2[352 ] = 1569987550;
      IS1[353 ] = 1145504660; IS2[353 ] = 1622132321;
      IS1[354 ] =  885012929; IS2[354 ] = 1092533602;
      IS1[355 ] = 1224139137; IS2[355 ] = 1840751652;
      IS1[356 ] = 1930672614; IS2[356 ] = 1637040668;
      IS1[357 ] = 1073822199; IS2[357 ] = 1783545033;
      IS1[358 ] =  921097169; IS2[358 ] = 1938213600;
      IS1[359 ] = 1425751390; IS2[359 ] = 1172781558;
      IS1[360 ] = 1214921245; IS2[360 ] =  825459365;
      IS1[361 ] = 2040403296; IS2[361 ] = 1533648867;
      IS1[362 ] =  896839067; IS2[362 ] = 1828625926;
      IS1[363 ] = 1663221476; IS2[363 ] =  623151978;
      IS1[364 ] = 1820439500; IS2[364 ] = 1190628682;
      IS1[365 ] = 1787778713; IS2[365 ] = 2139868031;
      IS1[366 ] = 1312411209; IS2[366 ] =  124016638;
      IS1[367 ] = 1511089893; IS2[367 ] = 1123651614;
      IS1[368 ] =  630843631; IS2[368 ] = 1212169916;
      IS1[369 ] = 1564923744; IS2[369 ] =  514797409;
      IS1[370 ] =  258126650; IS2[370 ] = 1007434416;
      IS1[371 ] = 1135521143; IS2[371 ] =  630960850;
      IS1[372 ] =  830276638; IS2[372 ] =  851724397;
      IS1[373 ] =  376035913; IS2[373 ] =  808391452;
      IS1[374 ] = 2023656286; IS2[374 ] =  465517081;
      IS1[375 ] = 1766701603; IS2[375 ] = 1993165647;
      IS1[376 ] = 1259750908; IS2[376 ] =  951063358;
      IS1[377 ] = 1644557611; IS2[377 ] = 1363160828;
      IS1[378 ] = 1583850975; IS2[378 ] = 1328161074;
      IS1[379 ] =  579027304; IS2[379 ] = 1626248155;
      IS1[380 ] = 1175583557; IS2[380 ] = 1137630999;
      IS1[381 ] =  253025339; IS2[381 ] =  222102409;
      IS1[382 ] = 1864136836; IS2[382 ] = 1013284156;
      IS1[383 ] =  447381646; IS2[383 ] =  720482175;
      IS1[384 ] = 1422107210; IS2[384 ] =  101344372;
      IS1[385 ] =   60187744; IS2[385 ] = 1036992166;
      IS1[386 ] =  736865928; IS2[386 ] = 1011413694;
      IS1[387 ] = 1223364420; IS2[387 ] = 1661421549;
      IS1[388 ] =  199571836; IS2[388 ] = 1149316001;
      IS1[389 ] =  118600702; IS2[389 ] =  498570418;
      IS1[390 ] = 1053021159; IS2[390 ] = 1200634496;
      IS1[391 ] = 1054745378; IS2[391 ] = 1927760137;
      IS1[392 ] = 1594338747; IS2[392 ] = 1093514040;
      IS1[393 ] =  694107219; IS2[393 ] =  541441173;
      IS1[394 ] = 1957890210; IS2[394 ] =  910785465;
      IS1[395 ] = 1470833484; IS2[395 ] =  426201828;
      IS1[396 ] =  585868751; IS2[396 ] = 1069807689;
      IS1[397 ] = 1134131445; IS2[397 ] = 1067457516;
      IS1[398 ] =  921089340; IS2[398 ] = 1038655815;
      IS1[399 ] =  616921523; IS2[399 ] = 1366192583;
      IS1[400 ] =   41330973; IS2[400 ] =  255560572;
      IS1[401 ] =  164189022; IS2[401 ] =   49014916;
      IS1[402 ] = 1161928238; IS2[402 ] =  209936365;
      IS1[403 ] = 2020361273; IS2[403 ] = 1950362287;
      IS1[404 ] =  511719771; IS2[404 ] =  325061593;
      IS1[405 ] =  454126255; IS2[405 ] = 1574510902;
      IS1[406 ] = 1859816564; IS2[406 ] = 1632823687;
      IS1[407 ] = 1432776991; IS2[407 ] =  227252761;
      IS1[408 ] =  765835313; IS2[408 ] = 2029746355;
      IS1[409 ] =  547912517; IS2[409 ] =  591050613;
      IS1[410 ] = 1863964385; IS2[410 ] =  712242677;
      IS1[411 ] = 1995501485; IS2[411 ] = 1312458862;
      IS1[412 ] = 1905501124; IS2[412 ] =  276433488;
      IS1[413 ] =  651011325; IS2[413 ] =  278589745;
      IS1[414 ] =  616628955; IS2[414 ] = 1231979682;
      IS1[415 ] =  625576690; IS2[415 ] =   76923407;
      IS1[416 ] =  214901782; IS2[416 ] =  956950803;
      IS1[417 ] = 1874851100; IS2[417 ] =  630418924;
      IS1[418 ] = 1796015654; IS2[418 ] = 1799191741;
      IS1[419 ] =  378100074; IS2[419 ] =  405443763;
      IS1[420 ] = 2115599125; IS2[420 ] = 1158325172;
      IS1[421 ] =   51748422; IS2[421 ] = 1108481669;
      IS1[422 ] = 1519507484; IS2[422 ] = 1078568949;
      IS1[423 ] =  231731663; IS2[423 ] = 1639105310;
      IS1[424 ] = 1813510342; IS2[424 ] = 1756686695;
      IS1[425 ] = 1180641209; IS2[425 ] =  817484587;
      IS1[426 ] =  547335020; IS2[426 ] = 1740348176;
      IS1[427 ] =  981294631; IS2[427 ] = 1645483561;
      IS1[428 ] =  267267642; IS2[428 ] =  525489921;
      IS1[429 ] =  975682868; IS2[429 ] = 1681199536;
      IS1[430 ] = 1070208555; IS2[430 ] =  572418479;
      IS1[431 ] = 1998085041; IS2[431 ] = 1713183034;
      IS1[432 ] =  173175676; IS2[432 ] =  486314861;
      IS1[433 ] = 1393518568; IS2[433 ] =   48960372;
      IS1[434 ] = 1341932065; IS2[434 ] = 1865938542;
      IS1[435 ] =  801752013; IS2[435 ] =  341065424;
      IS1[436 ] = 2132704165; IS2[436 ] =  494928752;
      IS1[437 ] =  422454854; IS2[437 ] =  743233701;
      IS1[438 ] = 1402307772; IS2[438 ] =  303431257;
      IS1[439 ] = 1557964062; IS2[439 ] = 1825819623;
      IS1[440 ] = 1894865913; IS2[440 ] = 1365259674;
      IS1[441 ] =  622187675; IS2[441 ] = 1865578472;
      IS1[442 ] = 1326023361; IS2[442 ] =  122041713;
      IS1[443 ] = 1730100218; IS2[443 ] =  365085301;
      IS1[444 ] =  243501496; IS2[444 ] = 1010792790;
      IS1[445 ] = 1899088012; IS2[445 ] =  162527744;
      IS1[446 ] =  421847337; IS2[446 ] = 2076217683;
      IS1[447 ] = 1760083121; IS2[447 ] =  891099538;
      IS1[448 ] = 2057638875; IS2[448 ] = 1503480695;
      IS1[449 ] =  195671618; IS2[449 ] = 2085854796;
      IS1[450 ] = 1810716138; IS2[450 ] =  477554292;
      IS1[451 ] = 1801885889; IS2[451 ] =  710751106;
      IS1[452 ] = 1154047452; IS2[452 ] = 1841903658;
      IS1[453 ] =  667706413; IS2[453 ] = 1964945942;
      IS1[454 ] = 2017974505; IS2[454 ] = 1663725788;
      IS1[455 ] =  884726865; IS2[455 ] =  522792338;
      IS1[456 ] =  926150544; IS2[456 ] = 1559103945;
      IS1[457 ] =  690922273; IS2[457 ] =   59346276;
      IS1[458 ] =  712157685; IS2[458 ] =  287197618;
      IS1[459 ] = 2059778138; IS2[459 ] =  815135778;
      IS1[460 ] = 1190461438; IS2[460 ] = 1905576318;
      IS1[461 ] = 1564333110; IS2[461 ] =  542814824;
      IS1[462 ] =  271883897; IS2[462 ] =  163552060;
      IS1[463 ] =  344754143; IS2[463 ] = 1429455140;
      IS1[464 ] =  668346479; IS2[464 ] =  907377254;
      IS1[465 ] =  906777901; IS2[465 ] =  973540088;
      IS1[466 ] =  211937991; IS2[466 ] = 1829082247;
      IS1[467 ] = 2057662804; IS2[467 ] =  466664141;
      IS1[468 ] =  685469316; IS2[468 ] =  801959481;
      IS1[469 ] = 1806627086; IS2[469 ] = 1933314854;
      IS1[470 ] =   27530138; IS2[470 ] = 1590842779;
      IS1[471 ] = 1972873114; IS2[471 ] = 1790212125;
      IS1[472 ] = 2028079049; IS2[472 ] =  909199739;
      IS1[473 ] =  568804232; IS2[473 ] =  180866254;
      IS1[474 ] = 1897611427; IS2[474 ] = 1205049755;
      IS1[475 ] = 1529982236; IS2[475 ] =  727044114;
      IS1[476 ] = 1930649184; IS2[476 ] = 1157145550;
      IS1[477 ] =  579071696; IS2[477 ] = 1002305204;
      IS1[478 ] = 1522526588; IS2[478 ] =  102634606;
      IS1[479 ] =  502645745; IS2[479 ] = 2002850332;
      IS1[480 ] = 1455547110; IS2[480 ] = 1863533555;
      IS1[481 ] = 1865827035; IS2[481 ] =   87808690;
      IS1[482 ] = 2146621882; IS2[482 ] = 2037357581;
      IS1[483 ] = 1348397757; IS2[483 ] =  252539429;
      IS1[484 ] =   74374264; IS2[484 ] =  584457262;
      IS1[485 ] =  909905711; IS2[485 ] =  641791254;
      IS1[486 ] =  406627724; IS2[486 ] = 1970801280;
      IS1[487 ] = 1869585944; IS2[487 ] = 2086692085;
      IS1[488 ] =  471312487; IS2[488 ] =   58607088;
      IS1[489 ] =  314332810; IS2[489 ] = 1762002696;
      IS1[490 ] =  714059427; IS2[490 ] =  892062884;
      IS1[491 ] =   88551984; IS2[491 ] =  707506711;
      IS1[492 ] = 1764182800; IS2[492 ] = 1566043351;
      IS1[493 ] = 1660801101; IS2[493 ] =   91743653;
      IS1[494 ] = 2053618389; IS2[494 ] =   87543326;
      IS1[495 ] = 1872198017; IS2[495 ] = 1419845282;
      IS1[496 ] = 1718051970; IS2[496 ] =   16608354;
      IS1[497 ] =  106783453; IS2[497 ] = 1579552005;
      IS1[498 ] =  198639753; IS2[498 ] = 1260841572;
      IS1[499 ] =  444341049; IS2[499 ] =  185823881;
      IS1[500 ] =  822270701; IS2[500 ] =  703588888;
      IS1[501 ] = 1741688389; IS2[501 ] = 1199341216;
      IS1[502 ] = 1740532989; IS2[502 ] =  506564970;
      IS1[503 ] =  306294863; IS2[503 ] = 1989939786;
      IS1[504 ] = 1282800170; IS2[504 ] =  909139593;
      IS1[505 ] =  869069844; IS2[505 ] =  759737034;
      IS1[506 ] = 1099376549; IS2[506 ] =  592166139;
      IS1[507 ] =  530526621; IS2[507 ] =  372525993;
      IS1[508 ] =  132237605; IS2[508 ] = 1235725904;
      IS1[509 ] = 1250152263; IS2[509 ] =  734059529;
      IS1[510 ] =  997470030; IS2[510 ] =  853534414;
      IS1[511 ] = 1214905199; IS2[511 ] = 1227201813;
      IS1[512 ] = 1024791465; IS2[512 ] = 1264083624;
      IS1[513 ] =  215462734; IS2[513 ] = 1323361591;
      IS1[514 ] = 1746861828; IS2[514 ] = 1171508714;
      IS1[515 ] =   98049234; IS2[515 ] =      76692;
      IS1[516 ] = 1941536515; IS2[516 ] = 1557540564;
      IS1[517 ] = 1062043665; IS2[517 ] = 1536279542;
      IS1[518 ] = 1933220889; IS2[518 ] =  531737550;
      IS1[519 ] = 1920928412; IS2[519 ] = 1961452165;
      IS1[520 ] =  924668593; IS2[520 ] = 1126181753;
      IS1[521 ] =  483780576; IS2[521 ] = 1943663885;
      IS1[522 ] = 1172795790; IS2[522 ] =  902336756;
      IS1[523 ] =   57254329; IS2[523 ] =  548344687;
      IS1[524 ] =  742985485; IS2[524 ] =   75812010;
      IS1[525 ] = 1138516383; IS2[525 ] =  704793701;
      IS1[526 ] =  516069233; IS2[526 ] =  658536656;
      IS1[527 ] =  767025613; IS2[527 ] = 1696634101;
      IS1[528 ] = 1183876758; IS2[528 ] =  436794231;
      IS1[529 ] =   68594352; IS2[529 ] = 1019241011;
      IS1[530 ] = 1445194529; IS2[530 ] = 1652427192;
      IS1[531 ] = 1984022317; IS2[531 ] =  634620520;
      IS1[532 ] = 1581637534; IS2[532 ] =  715105076;
      IS1[533 ] = 1852560766; IS2[533 ] = 1849389822;
      IS1[534 ] = 1786797677; IS2[534 ] = 1775016593;
      IS1[535 ] =  282102855; IS2[535 ] =  812059288;
      IS1[536 ] = 2120900530; IS2[536 ] = 1367713004;
      IS1[537 ] =  483020346; IS2[537 ] =  224667086;
      IS1[538 ] = 1686388582; IS2[538 ] = 1056110078;
      IS1[539 ] = 1080296513; IS2[539 ] = 1000944206;
      IS1[540 ] = 1128089199; IS2[540 ] =  211017438;
      IS1[541 ] = 1873948292; IS2[541 ] = 1459653839;
      IS1[542 ] = 1478295042; IS2[542 ] = 2091011370;
      IS1[543 ] =  780336939; IS2[543 ] = 1878440217;
      IS1[544 ] =   89453090; IS2[544 ] =  496618994;
      IS1[545 ] =   84723786; IS2[545 ] =  193890549;
      IS1[546 ] =  916751048; IS2[546 ] =  529107147;
      IS1[547 ] =  937673116; IS2[547 ] =  141932466;
      IS1[548 ] = 1522160349; IS2[548 ] =  617355429;
      IS1[549 ] = 1323199052; IS2[549 ] =  759342156;
      IS1[550 ] = 1238578537; IS2[550 ] = 1328836664;
      IS1[551 ] = 1112923558; IS2[551 ] = 1026465383;
      IS1[552 ] =  505271372; IS2[552 ] = 1384592778;
      IS1[553 ] =  566594858; IS2[553 ] =  934194365;
      IS1[554 ] = 1039240942; IS2[554 ] = 2009330113;
      IS1[555 ] = 1832500056; IS2[555 ] = 1947683833;
      IS1[556 ] =  360636801; IS2[556 ] = 1001734064;
      IS1[557 ] =  669874416; IS2[557 ] = 1595554733;
      IS1[558 ] = 1832363129; IS2[558 ] = 1805004882;
      IS1[559 ] = 1913361231; IS2[559 ] = 1861860225;
      IS1[560 ] = 1042646163; IS2[560 ] = 1027660606;
      IS1[561 ] = 1997920265; IS2[561 ] = 1340547150;
      IS1[562 ] =  212027369; IS2[562 ] = 1980768873;
      IS1[563 ] =  955366244; IS2[563 ] = 1250568151;
      IS1[564 ] =   62449968; IS2[564 ] = 2054390312;
      IS1[565 ] =  407925365; IS2[565 ] =  387411782;
      IS1[566 ] =  136158279; IS2[566 ] =  868480095;
      IS1[567 ] =  383555850; IS2[567 ] = 1375373714;
      IS1[568 ] = 1578435951; IS2[568 ] =  212038261;
      IS1[569 ] =  358815004; IS2[569 ] = 1539319712;
      IS1[570 ] =  207947798; IS2[570 ] =   84668320;
      IS1[571 ] = 1234390761; IS2[571 ] =  410483487;
      IS1[572 ] = 1493324788; IS2[572 ] =  487843369;
      IS1[573 ] = 2056136989; IS2[573 ] = 1938329879;
      IS1[574 ] = 1636698515; IS2[574 ] =  698986119;
      IS1[575 ] =  304252344; IS2[575 ] = 1411426968;
      IS1[576 ] =  541389504; IS2[576 ] =  660279563;
      IS1[577 ] = 1976796554; IS2[577 ] =  721255515;
      IS1[578 ] = 1994743666; IS2[578 ] =  260308500;
      IS1[579 ] = 1129912793; IS2[579 ] =   70624725;
      IS1[580 ] =  394341156; IS2[580 ] =  160687023;
      IS1[581 ] =  698480071; IS2[581 ] = 1131283564;
      IS1[582 ] =   15953128; IS2[582 ] = 1327808851;
      IS1[583 ] = 2047761093; IS2[583 ] = 1655969917;
      IS1[584 ] = 1514379033; IS2[584 ] =  628302318;
      IS1[585 ] =  455080996; IS2[585 ] = 2123515005;
      IS1[586 ] = 1145649301; IS2[586 ] = 1563293737;
      IS1[587 ] = 1631296493; IS2[587 ] = 1232136613;
      IS1[588 ] =  376739480; IS2[588 ] =  694184162;
      IS1[589 ] =   62700700; IS2[589 ] = 1621822750;
      IS1[590 ] = 1443922028; IS2[590 ] =  394412632;
      IS1[591 ] = 2034833661; IS2[591 ] =   56369217;
      IS1[592 ] =  179925148; IS2[592 ] =   46595906;
      IS1[593 ] =  943763422; IS2[593 ] =  359712423;
      IS1[594 ] =  294551736; IS2[594 ] = 2091885037;
      IS1[595 ] = 1612117876; IS2[595 ] = 1416720025;
      IS1[596 ] = 2060353278; IS2[596 ] = 1999037873;
      IS1[597 ] = 2027668315; IS2[597 ] = 1597386667;
      IS1[598 ] =  565740086; IS2[598 ] = 1681038972;
      IS1[599 ] =   85508641; IS2[599 ] =  422362053;
      IS1[600 ] =  632136951; IS2[600 ] =  967303111;
      IS1[601 ] =  287527126; IS2[601 ] = 2006085515;
      IS1[602 ] =  872188325; IS2[602 ] =  973826090;
      IS1[603 ] = 1528678858; IS2[603 ] =  870128621;
      IS1[604 ] = 2133000311; IS2[604 ] = 1630341425;
      IS1[605 ] =    9267403; IS2[605 ] =   82683218;
      IS1[606 ] = 1055725918; IS2[606 ] = 1071217436;
      IS1[607 ] = 1601015878; IS2[607 ] =  460901436;
      IS1[608 ] = 1173244981; IS2[608 ] =  575155810;
      IS1[609 ] =  190740363; IS2[609 ] =  209424492;
      IS1[610 ] =  373750904; IS2[610 ] = 1037445515;
      IS1[611 ] = 1903087315; IS2[611 ] =  199145625;
      IS1[612 ] =  935216143; IS2[612 ] =  767883538;
      IS1[613 ] =  616643835; IS2[613 ] =  249787085;
      IS1[614 ] = 1269743498; IS2[614 ] = 2038689023;
      IS1[615 ] = 1906186820; IS2[615 ] =  455180313;
      IS1[616 ] = 1341988859; IS2[616 ] = 1474088526;
      IS1[617 ] =  999778032; IS2[617 ] = 1005749219;
      IS1[618 ] =  620954154; IS2[618 ] = 1241433020;
      IS1[619 ] = 1863565816; IS2[619 ] = 1065421906;
      IS1[620 ] = 1497283145; IS2[620 ] = 2115805116;
      IS1[621 ] =  406679876; IS2[621 ] =  510949186;
      IS1[622 ] =  912988730; IS2[622 ] =   49564907;
      IS1[623 ] = 1764585681; IS2[623 ] =  313681775;
      IS1[624 ] =  493773352; IS2[624 ] =  791156315;
      IS1[625 ] = 2086722153; IS2[625 ] =  107352467;
      IS1[626 ] = 1286084681; IS2[626 ] =  208200070;
      IS1[627 ] =  808109356; IS2[627 ] = 1731413771;
      IS1[628 ] = 1486614524; IS2[628 ] = 2064085170;
      IS1[629 ] = 2129645386; IS2[629 ] = 1284378691;
      IS1[630 ] = 1014423030; IS2[630 ] =  701300786;
      IS1[631 ] =   60658846; IS2[631 ] =  397903983;
      IS1[632 ] = 1048804021; IS2[632 ] = 1593351979;
      IS1[633 ] = 1842774497; IS2[633 ] = 1285982663;
      IS1[634 ] =  854481213; IS2[634 ] = 1322235277;
      IS1[635 ] =  798953202; IS2[635 ] = 1538431839;
      IS1[636 ] = 1052424288; IS2[636 ] = 1822372594;
      IS1[637 ] =  601068089; IS2[637 ] = 1425648861;
      IS1[638 ] = 1849885612; IS2[638 ] =  733138526;
      IS1[639 ] = 1269451977; IS2[639 ] =  917070258;
      IS1[640 ] = 2144361786; IS2[640 ] = 1842934549;
      IS1[641 ] = 1496984691; IS2[641 ] = 1297819612;
      IS1[642 ] =  438239309; IS2[642 ] = 1148023460;
      IS1[643 ] = 1737313181; IS2[643 ] =  645117283;
      IS1[644 ] =  430834434; IS2[644 ] =  284660368;
      IS1[645 ] = 1143775914; IS2[645 ] = 1381744399;
      IS1[646 ] = 1995585326; IS2[646 ] = 2120820043;
      IS1[647 ] =  135161363; IS2[647 ] = 1799867404;
      IS1[648 ] =  836364692; IS2[648 ] = 1029908703;
      IS1[649 ] = 1757227577; IS2[649 ] =  546266239;
      IS1[650 ] =  623798852; IS2[650 ] =  585377965;
      IS1[651 ] = 1098341577; IS2[651 ] =  636341909;
      IS1[652 ] = 1520291052; IS2[652 ] =  913917239;
      IS1[653 ] =  975535163; IS2[653 ] = 1598597453;
      IS1[654 ] = 1447444469; IS2[654 ] = 1937942110;
      IS1[655 ] = 1432720174; IS2[655 ] =  839452541;
      IS1[656 ] = 1295341632; IS2[656 ] = 2103887297;
      IS1[657 ] = 1300522179; IS2[657 ] =  809881664;
      IS1[658 ] = 1425114463; IS2[658 ] = 1208521323;
      IS1[659 ] = 1257800427; IS2[659 ] = 1205340370;
      IS1[660 ] = 1186021021; IS2[660 ] = 1115454768;
      IS1[661 ] = 1164381967; IS2[661 ] = 1695441674;
      IS1[662 ] = 1339932055; IS2[662 ] =  626668376;
      IS1[663 ] = 1054295865; IS2[663 ] = 1077719907;
      IS1[664 ] =  908887630; IS2[664 ] =  375160191;
      IS1[665 ] = 1172794729; IS2[665 ] = 1568832201;
      IS1[666 ] =  192588897; IS2[666 ] =  917870514;
      IS1[667 ] =  643461236; IS2[667 ] =  484049433;
      IS1[668 ] =  315353239; IS2[668 ] =  397852714;
      IS1[669 ] =  931623820; IS2[669 ] = 1720458459;
      IS1[670 ] = 1580993163; IS2[670 ] =  734491598;
      IS1[671 ] = 1513146206; IS2[671 ] = 1441389702;
      IS1[672 ] =  701699279; IS2[672 ] = 1253582219;
      IS1[673 ] = 1597280672; IS2[673 ] =  920294785;
      IS1[674 ] = 2083951540; IS2[674 ] = 1899768161;
      IS1[675 ] = 1604320624; IS2[675 ] =  154338943;
      IS1[676 ] = 1832407686; IS2[676 ] = 2039753469;
      IS1[677 ] = 1522755360; IS2[677 ] = 1558890156;
      IS1[678 ] = 1696355490; IS2[678 ] = 1445912335;
      IS1[679 ] = 1308491933; IS2[679 ] = 1267480279;
      IS1[680 ] = 1946363807; IS2[680 ] = 1896823905;
      IS1[681 ] = 2028942171; IS2[681 ] = 1457946439;
      IS1[682 ] =  510086891; IS2[682 ] =  536540300;
      IS1[683 ] = 1486809993; IS2[683 ] = 1547406458;
      IS1[684 ] = 1317061925; IS2[684 ] = 1591791104;
      IS1[685 ] =  989995158; IS2[685 ] = 1604821909;
      IS1[686 ] =  326298337; IS2[686 ] = 1171640652;
      IS1[687 ] =  152659804; IS2[687 ] =  495837027;
      IS1[688 ] = 1222941036; IS2[688 ] =  796199205;
      IS1[689 ] =  907350872; IS2[689 ] =  727966425;
      IS1[690 ] = 1686084314; IS2[690 ] = 1613446594;
      IS1[691 ] = 1642955746; IS2[691 ] = 1880495478;
      IS1[692 ] =  247284465; IS2[692 ] = 1706304962;
      IS1[693 ] = 1611723103; IS2[693 ] = 1706608231;
      IS1[694 ] =  719689499; IS2[694 ] =   31477369;
      IS1[695 ] = 1079226399; IS2[695 ] = 1782920006;
      IS1[696 ] = 1420885629; IS2[696 ] = 1072981519;
      IS1[697 ] = 2084453400; IS2[697 ] =  224386433;
      IS1[698 ] = 1047203160; IS2[698 ] =  614309249;
      IS1[699 ] =  586617884; IS2[699 ] = 1582667003;
      IS1[700 ] =   11984589; IS2[700 ] = 1558266095;
      IS1[701 ] =  755652237; IS2[701 ] =  866088075;
      IS1[702 ] = 2017712409; IS2[702 ] = 1954340696;
      IS1[703 ] =  589844984; IS2[703 ] =  314030536;
      IS1[704 ] =  129229479; IS2[704 ] =  272111716;
      IS1[705 ] =  870276471; IS2[705 ] =  465006906;
      IS1[706 ] = 1612019958; IS2[706 ] =   21334935;
      IS1[707 ] = 1315108425; IS2[707 ] = 1818825397;
      IS1[708 ] = 1258602567; IS2[708 ] = 1066619326;
      IS1[709 ] =  166077102; IS2[709 ] =  909070060;
      IS1[710 ] = 1848198256; IS2[710 ] =  385880783;
      IS1[711 ] =  488779996; IS2[711 ] =  241674453;
      IS1[712 ] = 1098298571; IS2[712 ] = 1348905710;
      IS1[713 ] =  561394508; IS2[713 ] = 1782802528;
      IS1[714 ] =  254791280; IS2[714 ] =  354432011;
      IS1[715 ] = 1214976755; IS2[715 ] =  774386703;
      IS1[716 ] =   90779321; IS2[716 ] = 1418378337;
      IS1[717 ] = 1355457939; IS2[717 ] = 1690663694;
      IS1[718 ] = 1541432462; IS2[718 ] = 1040751740;
      IS1[719 ] =  394291901; IS2[719 ] = 1603027222;
      IS1[720 ] = 1239001540; IS2[720 ] =  146841931;
      IS1[721 ] =   38818735; IS2[721 ] = 1839009879;
      IS1[722 ] = 1674954341; IS2[722 ] = 1920405940;
      IS1[723 ] = 1640316191; IS2[723 ] = 1265350958;
      IS1[724 ] = 1937228975; IS2[724 ] =  996533450;
      IS1[725 ] =  423424194; IS2[725 ] =  561329945;
      IS1[726 ] =  323756417; IS2[726 ] = 1682627390;
      IS1[727 ] = 1549193098; IS2[727 ] =  805029685;
      IS1[728 ] = 1360079132; IS2[728 ] =  925188637;
      IS1[729 ] =   78241293; IS2[729 ] = 1040585429;
      IS1[730 ] = 1074975265; IS2[730 ] =  797619830;
      IS1[731 ] = 1582340080; IS2[731 ] =  721022974;
      IS1[732 ] =  629043128; IS2[732 ] =  610470736;
      IS1[733 ] = 1657680582; IS2[733 ] =  963797322;
      IS1[734 ] = 1011918751; IS2[734 ] = 1560942165;
      IS1[735 ] =   42122891; IS2[735 ] = 1481369897;
      IS1[736 ] = 2043026443; IS2[736 ] = 1863587933;
      IS1[737 ] = 1332181389; IS2[737 ] =  854812560;
      IS1[738 ] = 1986410252; IS2[738 ] = 1040318983;
      IS1[739 ] =  521386266; IS2[739 ] = 1950110774;
      IS1[740 ] =  216896173; IS2[740 ] = 1685734611;
      IS1[741 ] = 1867435681; IS2[741 ] = 2057062478;
      IS1[742 ] = 1401811081; IS2[742 ] = 1014131666;
      IS1[743 ] = 1538740757; IS2[743 ] =  620335187;
      IS1[744 ] =  923455231; IS2[744 ] = 1901856320;
      IS1[745 ] =  597449802; IS2[745 ] =  706565272;
      IS1[746 ] =  683455885; IS2[746 ] = 1033766701;
      IS1[747 ] =  828147504; IS2[747 ] = 1580010438;
      IS1[748 ] =  451283302; IS2[748 ] =  781324118;
      IS1[749 ] = 2048256218; IS2[749 ] = 1184873148;
      IS1[750 ] = 1583574204; IS2[750 ] = 1032841150;
      IS1[751 ] = 1652016656; IS2[751 ] =   33056864;
      IS1[752 ] = 1532136667; IS2[752 ] = 1581149947;
      IS1[753 ] = 1462299459; IS2[753 ] =  687082954;
      IS1[754 ] =  698944454; IS2[754 ] = 2122885802;
      IS1[755 ] = 1195045208; IS2[755 ] = 1678424394;
      IS1[756 ] =  502707485; IS2[756 ] = 1444358879;
      IS1[757 ] =  941731361; IS2[757 ] = 1717503286;
      IS1[758 ] =  206769672; IS2[758 ] = 1166864868;
      IS1[759 ] = 1312487368; IS2[759 ] = 1171921906;
      IS1[760 ] =   70908405; IS2[760 ] = 1544257314;
      IS1[761 ] =  666380866; IS2[761 ] = 1849868712;
      IS1[762 ] =  305790335; IS2[762 ] = 2037569416;
      IS1[763 ] =  334326322; IS2[763 ] = 1721074287;
      IS1[764 ] = 2101405245; IS2[764 ] = 1240524239;
      IS1[765 ] =  760831995; IS2[765 ] =  351651496;
      IS1[766 ] = 2076975119; IS2[766 ] =  732446407;
      IS1[767 ] = 1457683031; IS2[767 ] =  930496443;
      IS1[768 ] = 1947001405; IS2[768 ] =  590201248;
      IS1[769 ] =  275187592; IS2[769 ] =   62024761;
      IS1[770 ] = 1289777261; IS2[770 ] =  542257293;
      IS1[771 ] =  363521237; IS2[771 ] =  517555072;
      IS1[772 ] = 1063822524; IS2[772 ] =   90991909;
      IS1[773 ] =  312343146; IS2[773 ] = 1445115042;
      IS1[774 ] = 2081825980; IS2[774 ] = 1071980321;
      IS1[775 ] = 1790415941; IS2[775 ] =  818819165;
      IS1[776 ] =   28581357; IS2[776 ] =  877661732;
      IS1[777 ] = 1803064654; IS2[777 ] = 1992656309;
      IS1[778 ] = 1030876307; IS2[778 ] =  513683199;
      IS1[779 ] =  492391370; IS2[779 ] = 1206877439;
      IS1[780 ] =  211680892; IS2[780 ] = 1783884173;
      IS1[781 ] =  744879183; IS2[781 ] =  984195787;
      IS1[782 ] =  490863602; IS2[782 ] = 1663478249;
      IS1[783 ] =  756240663; IS2[783 ] =  927897638;
      IS1[784 ] =  980102031; IS2[784 ] = 1228241271;
      IS1[785 ] = 1517579622; IS2[785 ] = 1565288529;
      IS1[786 ] = 2018942706; IS2[786 ] = 1486937165;
      IS1[787 ] =  914892050; IS2[787 ] = 1011671153;
      IS1[788 ] = 2071502238; IS2[788 ] =  910410508;
      IS1[789 ] = 1144674111; IS2[789 ] = 1035198034;
      IS1[790 ] =  551571043; IS2[790 ] = 1704515839;
      IS1[791 ] = 2067235260; IS2[791 ] =   62628367;
      IS1[792 ] = 1183655237; IS2[792 ] = 1825596188;
      IS1[793 ] = 1720738448; IS2[793 ] = 1426908311;
      IS1[794 ] = 1428394554; IS2[794 ] = 1331528227;
      IS1[795 ] = 1236406902; IS2[795 ] =  739342228;
      IS1[796 ] =  314637105; IS2[796 ] = 1635162390;
      IS1[797 ] = 1631561757; IS2[797 ] = 1259987681;
      IS1[798 ] = 1113570671; IS2[798 ] =  941650185;
      IS1[799 ] = 1316684934; IS2[799 ] = 1653069887;
      IS1[800 ] = 1026427146; IS2[800 ] =  713191356;
      IS1[801 ] = 1616622265; IS2[801 ] = 1445073589;
      IS1[802 ] = 1184120268; IS2[802 ] =  844684601;
      IS1[803 ] = 1916487469; IS2[803 ] =  957905046;
      IS1[804 ] = 1196992255; IS2[804 ] =  143852508;
      IS1[805 ] =  810274414; IS2[805 ] =  670905422;
      IS1[806 ] =  994554837; IS2[806 ] = 1752024033;
      IS1[807 ] =  631323451; IS2[807 ] = 1938843572;
      IS1[808 ] = 2108796165; IS2[808 ] =  686849224;
      IS1[809 ] =  675199957; IS2[809 ] = 2066177454;
      IS1[810 ] =  500608835; IS2[810 ] =  501044809;
      IS1[811 ] =  318482166; IS2[811 ] =  929306814;
      IS1[812 ] = 1833923717; IS2[812 ] =  224005423;
      IS1[813 ] = 1473405260; IS2[813 ] =  584253050;
      IS1[814 ] = 1922717185; IS2[814 ] =  725230049;
      IS1[815 ] =  358747736; IS2[815 ] = 1894346138;
      IS1[816 ] = 1262935388; IS2[816 ] = 1123083929;
      IS1[817 ] =  170752885; IS2[817 ] =  282349087;
      IS1[818 ] = 1766873876; IS2[818 ] = 1639448540;
      IS1[819 ] = 1327238154; IS2[819 ] = 2086656898;
      IS1[820 ] =  608102000; IS2[820 ] = 1953835572;
      IS1[821 ] = 1933730221; IS2[821 ] =  592600179;
      IS1[822 ] =   20093493; IS2[822 ] = 1802818520;
      IS1[823 ] = 2047560831; IS2[823 ] = 1690864475;
      IS1[824 ] = 2120624346; IS2[824 ] = 1399722254;
      IS1[825 ] =  728200766; IS2[825 ] = 1251271247;
      IS1[826 ] = 1996802725; IS2[826 ] = 1182594334;
      IS1[827 ] = 1732977781; IS2[827 ] =  494629971;
      IS1[828 ] = 1333989141; IS2[828 ] = 1463491202;
      IS1[829 ] =  414434960; IS2[829 ] =  427148665;
      IS1[830 ] = 2058685774; IS2[830 ] = 1258425844;
      IS1[831 ] = 1570688571; IS2[831 ] = 1718768035;
      IS1[832 ] = 1361600109; IS2[832 ] =  732095097;
      IS1[833 ] = 1464532988; IS2[833 ] = 1592327040;
      IS1[834 ] = 1154029608; IS2[834 ] =  457075509;
      IS1[835 ] =  504833970; IS2[835 ] =  282459181;
      IS1[836 ] =  720633547; IS2[836 ] = 1754749459;
      IS1[837 ] = 1542782084; IS2[837 ] =  556876143;
      IS1[838 ] =  498983980; IS2[838 ] = 1448363633;
      IS1[839 ] =  815335207; IS2[839 ] = 1938426616;
      IS1[840 ] =  250652265; IS2[840 ] = 1560814150;
      IS1[841 ] =   61330462; IS2[841 ] = 1822327770;
      IS1[842 ] =  252780969; IS2[842 ] = 2058641830;
      IS1[843 ] = 2103405990; IS2[843 ] =  532243551;
      IS1[844 ] = 1908261717; IS2[844 ] = 2064263729;
      IS1[845 ] = 1977211665; IS2[845 ] =  777733094;
      IS1[846 ] =  237237934; IS2[846 ] =  528493150;
      IS1[847 ] =   73678412; IS2[847 ] =  219112977;
      IS1[848 ] = 1963696531; IS2[848 ] = 1929981191;
      IS1[849 ] =  101047790; IS2[849 ] =  216743219;
      IS1[850 ] =    2093800; IS2[850 ] =   74011539;
      IS1[851 ] =  632655512; IS2[851 ] =  198607329;
      IS1[852 ] =  941788307; IS2[852 ] = 1227252584;
      IS1[853 ] =  545304972; IS2[853 ] = 1963545088;
      IS1[854 ] =   55853502; IS2[854 ] =  498296470;
      IS1[855 ] = 1891044982; IS2[855 ] =  212219604;
      IS1[856 ] =  713967553; IS2[856 ] = 1094926756;
      IS1[857 ] =  126818203; IS2[857 ] =  569771356;
      IS1[858 ] = 1701715475; IS2[858 ] =  840368587;
      IS1[859 ] =  904419474; IS2[859 ] = 1160456593;
      IS1[860 ] =  914755144; IS2[860 ] = 1602166631;
      IS1[861 ] =  719105598; IS2[861 ] =  768277780;
      IS1[862 ] =  436565842; IS2[862 ] = 1526045329;
      IS1[863 ] =  642339928; IS2[863 ] = 1308502226;
      IS1[864 ] =  952737893; IS2[864 ] =  889261161;
      IS1[865 ] =  750349739; IS2[865 ] =  419661629;
      IS1[866 ] = 1240092349; IS2[866 ] =  385969468;
      IS1[867 ] =  947883679; IS2[867 ] =  858658062;
      IS1[868 ] =  412164731; IS2[868 ] =  227225801;
      IS1[869 ] = 1077025161; IS2[869 ] = 1809495121;
      IS1[870 ] =  623779545; IS2[870 ] = 1748437918;
      IS1[871 ] = 1462115422; IS2[871 ] = 1828054930;
      IS1[872 ] = 1793988879; IS2[872 ] =  971499218;
      IS1[873 ] =  401387846; IS2[873 ] = 2063529218;
      IS1[874 ] =  939523117; IS2[874 ] =  466501459;
      IS1[875 ] = 1325434731; IS2[875 ] =  933144734;
      IS1[876 ] = 1844466921; IS2[876 ] = 1723628496;
      IS1[877 ] =  612243172; IS2[877 ] =  675935430;
      IS1[878 ] =  122322491; IS2[878 ] =  700754464;
      IS1[879 ] = 1118985067; IS2[879 ] =   77681872;
      IS1[880 ] =  673808420; IS2[880 ] = 2027510724;
      IS1[881 ] = 1900746742; IS2[881 ] =  875745816;
      IS1[882 ] = 1951431584; IS2[882 ] = 1326338561;
      IS1[883 ] = 1185595160; IS2[883 ] = 1508083812;
      IS1[884 ] =   20387986; IS2[884 ] =    2502650;
      IS1[885 ] = 1448512176; IS2[885 ] =  406078533;
      IS1[886 ] =  341831642; IS2[886 ] =  444907341;
      IS1[887 ] = 1749622481; IS2[887 ] = 1458991053;
      IS1[888 ] =  179921081; IS2[888 ] = 1671951076;
      IS1[889 ] =  928183806; IS2[889 ] = 1991458904;
      IS1[890 ] =  603248036; IS2[890 ] = 1824888100;
      IS1[891 ] = 1063724212; IS2[891 ] = 1890874257;
      IS1[892 ] =   78830689; IS2[892 ] = 2043989239;
      IS1[893 ] = 1753470474; IS2[893 ] =  830139537;
      IS1[894 ] = 1566789535; IS2[894 ] = 1388524135;
      IS1[895 ] = 1518518357; IS2[895 ] = 1787113559;
      IS1[896 ] =  891266992; IS2[896 ] = 1954215738;
      IS1[897 ] =  725221948; IS2[897 ] =  471108830;
      IS1[898 ] =  396459481; IS2[898 ] =  560156443;
      IS1[899 ] =   54855828; IS2[899 ] = 1250479705;
      IS1[900 ] = 1671412588; IS2[900 ] =  238648368;
      IS1[901 ] = 1457198819; IS2[901 ] = 1108923041;
      IS1[902 ] = 1864168313; IS2[902 ] = 2034120136;
      IS1[903 ] =  737458311; IS2[903 ] =  420578061;
      IS1[904 ] =  693032926; IS2[904 ] = 1415068309;
      IS1[905 ] =  549217560; IS2[905 ] =  285822649;
      IS1[906 ] = 2028105476; IS2[906 ] =  257057932;
      IS1[907 ] =    5253877; IS2[907 ] =  467436652;
      IS1[908 ] =  124062287; IS2[908 ] =  913845906;
      IS1[909 ] =  983237044; IS2[909 ] = 1573260196;
      IS1[910 ] = 1688115577; IS2[910 ] = 1088025993;
      IS1[911 ] =  456193194; IS2[911 ] = 1237163793;
      IS1[912 ] =  818167884; IS2[912 ] = 2137528272;
      IS1[913 ] =     670825; IS2[913 ] =  326538226;
      IS1[914 ] = 1256025768; IS2[914 ] = 1517888550;
      IS1[915 ] = 1603246774; IS2[915 ] =  955989622;
      IS1[916 ] = 1723321062; IS2[916 ] =  756662277;
      IS1[917 ] =  548919448; IS2[917 ] =  670807456;
      IS1[918 ] =  765198119; IS2[918 ] = 1636394764;
      IS1[919 ] = 1940460545; IS2[919 ] = 1546053813;
      IS1[920 ] = 1253731346; IS2[920 ] =  462240916;
      IS1[921 ] =  452873281; IS2[921 ] = 1640963727;
      IS1[922 ] =  563455527; IS2[922 ] =  678314347;
      IS1[923 ] =  296152006; IS2[923 ] = 1921509503;
      IS1[924 ] = 1488389520; IS2[924 ] = 2075130919;
      IS1[925 ] = 1416431702; IS2[925 ] =  400649975;
      IS1[926 ] =  187540584; IS2[926 ] =  280946768;
      IS1[927 ] =  265064824; IS2[927 ] =  292386889;
      IS1[928 ] =   17812467; IS2[928 ] =  158666141;
      IS1[929 ] =  678154378; IS2[929 ] =  273032591;
      IS1[930 ] = 1008995524; IS2[930 ] =  225822851;
      IS1[931 ] =  932640361; IS2[931 ] = 1794357721;
      IS1[932 ] = 1256165221; IS2[932 ] = 1362164825;
      IS1[933 ] = 1566130220; IS2[933 ] =  880937875;
      IS1[934 ] = 1719055144; IS2[934 ] = 1475807990;
      IS1[935 ] = 1511255203; IS2[935 ] = 1398326134;
      IS1[936 ] =   40137809; IS2[936 ] =  844274026;
      IS1[937 ] = 1514777596; IS2[937 ] = 1176190656;
      IS1[938 ] =  862983829; IS2[938 ] = 1444976675;
      IS1[939 ] = 1419376580; IS2[939 ] =  897595798;
      IS1[940 ] = 2004253892; IS2[940 ] =  394099144;
      IS1[941 ] = 1208598747; IS2[941 ] =  710801376;
      IS1[942 ] = 2134068686; IS2[942 ] = 1199555698;
      IS1[943 ] =  383571065; IS2[943 ] = 1799022351;
      IS1[944 ] = 1896509659; IS2[944 ] = 1062349110;
      IS1[945 ] = 2089619889; IS2[945 ] =  103925910;
      IS1[946 ] =   15951717; IS2[946 ] = 1115129379;
      IS1[947 ] = 1209658212; IS2[947 ] = 1067867972;
      IS1[948 ] = 1533825359; IS2[948 ] =  707403901;
      IS1[949 ] =  878388702; IS2[949 ] =  641697588;
      IS1[950 ] =  370529890; IS2[950 ] =  306771059;
      IS1[951 ] = 1054078483; IS2[951 ] = 1574620834;
      IS1[952 ] = 1512913863; IS2[952 ] =  619850280;
      IS1[953 ] =  516041141; IS2[953 ] = 1250978720;
      IS1[954 ] =  964083750; IS2[954 ] = 1794316764;
      IS1[955 ] = 1222510673; IS2[955 ] = 1010199648;
      IS1[956 ] =  716095488; IS2[956 ] =  890400554;
      IS1[957 ] =  257132284; IS2[957 ] =  233348130;
      IS1[958 ] = 1317716326; IS2[958 ] = 2012396977;
      IS1[959 ] = 1982982843; IS2[959 ] =  584095052;
      IS1[960 ] = 1376656999; IS2[960 ] = 1982442382;
      IS1[961 ] = 1158326173; IS2[961 ] =  841561919;
      IS1[962 ] =  223831123; IS2[962 ] =  115318433;
      IS1[963 ] = 1612316526; IS2[963 ] =  190139923;
      IS1[964 ] =  811660944; IS2[964 ] = 1411875518;
      IS1[965 ] =  461792825; IS2[965 ] = 1500187934;
      IS1[966 ] = 2080858235; IS2[966 ] = 2112319914;
      IS1[967 ] =   34295036; IS2[967 ] =  281625837;
      IS1[968 ] = 1796061661; IS2[968 ] = 1731981711;
      IS1[969 ] =  338906098; IS2[969 ] =  474926566;
      IS1[970 ] = 1049296891; IS2[970 ] = 1071577019;
      IS1[971 ] = 1254386847; IS2[971 ] =  867276511;
      IS1[972 ] =  823835412; IS2[972 ] = 1169979512;
      IS1[973 ] =   12773203; IS2[973 ] = 1151341486;
      IS1[974 ] =  145388856; IS2[974 ] = 1918711308;
      IS1[975 ] = 1163564225; IS2[975 ] =  588792593;
      IS1[976 ] = 1585160972; IS2[976 ] = 1831535360;
      IS1[977 ] =  151747303; IS2[977 ] =  893561329;
      IS1[978 ] = 2103142931; IS2[978 ] = 1524302572;
      IS1[979 ] = 1076384122; IS2[979 ] = 1181095863;
      IS1[980 ] = 1047454590; IS2[980 ] =  212580845;
      IS1[981 ] = 1391900713; IS2[981 ] = 1595457237;
      IS1[982 ] = 1980182019; IS2[982 ] = 2099452891;
      IS1[983 ] = 1134036187; IS2[983 ] = 1086465811;
      IS1[984 ] = 1561015123; IS2[984 ] = 2087543080;
      IS1[985 ] = 1091162968; IS2[985 ] = 2039625632;
      IS1[986 ] =  905907344; IS2[986 ] = 1707901257;
      IS1[987 ] = 1097576632; IS2[987 ] =  399477627;
      IS1[988 ] =  854163724; IS2[988 ] =  727734495;
      IS1[989 ] =  993622338; IS2[989 ] = 1208219252;
      IS1[990 ] = 1137639335; IS2[990 ] = 1202540324;
      IS1[991 ] = 1405884001; IS2[991 ] = 1769921362;
      IS1[992 ] = 1002962410; IS2[992 ] =  726546017;
      IS1[993 ] = 1906752935; IS2[993 ] = 1700064870;
      IS1[994 ] =  547343311; IS2[994 ] =  947713211;
      IS1[995 ] =  154484285; IS2[995 ] = 1863597615;
      IS1[996 ] = 1283816755; IS2[996 ] =  282974211;
      IS1[997 ] = 1206555620; IS2[997 ] = 1748596686;
      IS1[998 ] =  830620896; IS2[998 ] =  922833376;
      IS1[999 ] =  470833195; IS2[999 ] = 2082830087;
      IS1[1000] = 1275774316; IS2[1000] =  515177342;
      IS1[1001] =  731077242; IS2[1001] = 2126489340;

	  if ((N > 0) && (N < 1002)){
		  *RSEED_.ISEED1=IS1[N];
       	  *RSEED_.ISEED2=IS2[N];
	  }else{
		  *RSEED_.ISEED1=1;
       	  *RSEED_.ISEED2=1;
	  }
}
	  
void enangr2_(ifstream &IRD, FILE *IWDUMP){

	string line;
	char *PCH;

	getline(IRD, line);
	char *sArray = (char *) malloc((line.length()+1) * sizeof(char));
	strcpy(sArray, line.c_str());
	PCH = strtok(sArray, " ");
	*CENANG_.NE = atoi(PCH);
	PCH = strtok(NULL, " ");
	*CENANG_.EL = atof(PCH);
	PCH = strtok(NULL, " ");
	*CENANG_.EU = atof(PCH);
	PCH = strtok(NULL, " ");
	*CENANG_.BSE = atof(PCH);
	PCH = strtok(NULL, " ");
	*CENANG_.RBSE = atof(PCH);
	PCH = strtok(NULL, " ");
	*CENANG_.LLE = atoi(PCH);
	fprintf(IWDUMP, "		 %d		%.10f		%.10f		%.10f		%.10E		%d\n", *CENANG_.NE, *CENANG_.EL, *CENANG_.EU, *CENANG_.BSE, *CENANG_.RBSE, *CENANG_.LLE);

	getline(IRD, line);
	sArray = (char *) realloc(NULL, (line.length()+1) * sizeof(char));
	strcpy(sArray, line.c_str());
	PCH = strtok(sArray, " ");
	for (int K = 1; K <= *CENANG_.NE; K++){
		for (int J = 1; J <= 2; J ++){
			for (int I = 1; I <= 3; I++){
				CENANG_.PDE[I-1][J-1][K-1] = atof(PCH);
				fprintf(IWDUMP, "		%f", CENANG_.PDE[I-1][J-1][K-1]);
				PCH = strtok(NULL, " ");
			}
		}
	}
	fprintf(IWDUMP, "\n");

	getline(IRD, line);
	sArray = (char *) realloc(NULL, (line.length()+1) * sizeof(char));
	strcpy(sArray, line.c_str());
	PCH = strtok(sArray, " ");
	for (int K = 1; K <= *CENANG_.NE; K++){
		for (int J = 1; J <= 2; J ++){
			for (int I = 1; I <= 3; I++){
				CENANG_.PDE2[I-1][J-1][K-1] = atof(PCH);
				fprintf(IWDUMP, "		%f", CENANG_.PDE2[I-1][J-1][K-1]);
				PCH = strtok(NULL, " ");
			}
		}
	}
	fprintf(IWDUMP, "\n");

	getline(IRD, line);
	sArray = (char *) realloc(NULL, (line.length()+1) * sizeof(char));
	strcpy(sArray, line.c_str());
	PCH = strtok(sArray, " ");
	*CENANG_.NTH = atoi(PCH);
	PCH = strtok(NULL, " ");
	*CENANG_.THL = atof(PCH);
	PCH = strtok(NULL, " ");
	*CENANG_.THU = atof(PCH);
	PCH = strtok(NULL, " ");
	*CENANG_.BSTH = atof(PCH);
	PCH = strtok(NULL, " ");
	*CENANG_.RBSTH = atof(PCH);
	PCH = strtok(NULL, " ");
	*CENANG_.NPH = atoi(PCH);
	PCH = strtok(NULL, " ");
	*CENANG_.BSPH = atof(PCH);
	PCH = strtok(NULL, " ");
	*CENANG_.RBSPH = atof(PCH);
	PCH = strtok(NULL, " ");
	*CENANG_.LLTH = atoi(PCH);
	fprintf(IWDUMP, "		 %d		%.10f		%.10f		%.10f		%.10f		%d		%.10f		%.10f		%d\n", 
		    *CENANG_.NTH, *CENANG_.THL, *CENANG_.THU, *CENANG_.BSTH, *CENANG_.RBSTH, *CENANG_.NPH, *CENANG_.BSPH, *CENANG_.RBSPH, *CENANG_.LLTH);


	getline(IRD, line);
	sArray = (char *) realloc(NULL, (line.length()+1) * sizeof(char));
	strcpy(sArray, line.c_str());
	PCH = strtok(sArray, " ");
	for (int K = 1; K <= *CENANG_.NPH; K++){
		for (int J = 1; J <= *CENANG_.NTH; J ++){
			for (int I = 1; I <= 3; I++){
				CENANG_.PDA[I-1][J-1][K-1] = atof(PCH);
				fprintf(IWDUMP, "		%f", CENANG_.PDA[I-1][J-1][K-1]);
				PCH = strtok(NULL, " ");
			}
		}
	}
	fprintf(IWDUMP, "\n");

	getline(IRD, line);
	sArray = (char *) realloc(NULL, (line.length()+1) * sizeof(char));
	strcpy(sArray, line.c_str());
	PCH = strtok(sArray, " ");
	for (int K = 1; K <= *CENANG_.NPH; K++){
		for (int J = 1; J <= *CENANG_.NTH; J ++){
			for (int I = 1; I <= 3; I++){
				CENANG_.PDA2[I-1][J-1][K-1] = atof(PCH);
				fprintf(IWDUMP, "		%f", CENANG_.PDA2[I-1][J-1][K-1]);
				PCH = strtok(NULL, " ");
			}
		}
	}
	fprintf(IWDUMP, "\n");

	for (int I = 1; I <= 3; I++){
		for (int J = 1; J <=2; J++){
			for (int K = 1; K <= *CENANG_.NE; K++){
				CENANG_.PDEP[K-1][J-1][I-1]=0.0e0;
              	CENANG_.LPDE[K-1][J-1][I-1]=0;
			}
		}
	}

	for (int I = 1; I <= 3; I++){
		for (int J = 1; J <= *CENANG_.NTH; J++){
			for (int K = 1; K <= *CENANG_.NPH; K++){
				CENANG_.PDAP[K-1][J-1][I-1]=0.0e0;
             	CENANG_.LPDA[K-1][J-1][I-1]=0;
			}
		}
	}

}

void imdetr2_(ifstream &IRD, FILE *IWDUMP){

	char *PCH;
	string line;


	getline(IRD, line);
	char *sArray = (char *) malloc((line.length()+1) * sizeof(char));
	strcpy(sArray, line.c_str());
	PCH = strtok(sArray, " ");
	*CIMDET_.NID = atoi(PCH);
	fprintf(IWDUMP, "	%d\n", *CIMDET_.NID);

	getline(IRD, line);
	sArray = (char *) realloc(NULL, (line.length()+1) * sizeof(char));
	strcpy(sArray, line.c_str());
	PCH = strtok(sArray, " ");
	for (int I = 1; I <= *CIMDET_.NID; I++){
		CIMDET_.EDEP[I-1] = atof(PCH);
		fprintf(IWDUMP, "		%f", CIMDET_.EDEP[I-1]);
		PCH = strtok(NULL, " ");
	}
	fprintf(IWDUMP, "\n");

	getline(IRD, line);
	sArray = (char *) realloc(NULL, (line.length()+1) * sizeof(char));
	strcpy(sArray, line.c_str());
	PCH = strtok(sArray, " ");
	for (int I = 1; I <= *CIMDET_.NID; I++){
		CIMDET_.EDEP2[I-1] = atof(PCH);
		fprintf(IWDUMP, "		%f", CIMDET_.EDEP2[I-1]);
		PCH = strtok(NULL, " ");
	}
	fprintf(IWDUMP, "\n");

	getline(IRD, line);
	sArray = (char *) realloc(NULL, (line.length()+1) * sizeof(char));
	strcpy(sArray, line.c_str());
	PCH = strtok(sArray, " ");
	for (int I = 1; I <= *CIMDET_.NID; I++){
		CIMDET_.IDCUT[I-1] = atoi(PCH);
		fprintf(IWDUMP, "		%d", CIMDET_.IDCUT[I-1]);
		PCH = strtok(NULL, " ");
	}
	fprintf(IWDUMP, "\n");

	for (int I = 1; I <= *CIMDET_.NID; I++){
		getline(IRD, line);
		sArray = (char *) realloc(NULL, (line.length()+1) * sizeof(char));
		strcpy(sArray, line.c_str());
		PCH = strtok(sArray, " ");
		strcpy(SPCDIO[I-1], PCH);
		fprintf(IWDUMP, "	%s\n", SPCDIO[I-1]);

		getline(IRD, line);
		char *sArray = (char *) malloc((line.length()+1) * sizeof(char));
		strcpy(sArray, line.c_str());
		PCH = strtok(sArray, " ");
		CIMDET_.NE[I-1] = atoi(PCH);
		PCH = strtok(NULL, " ");
		CIMDET_.EL[I-1] = atof(PCH);
		PCH = strtok(NULL, " ");
		CIMDET_.EU[I-1] = atof(PCH);
		PCH = strtok(NULL, " ");
		CIMDET_.BSE[I-1] = atof(PCH);
		PCH = strtok(NULL, " ");
		CIMDET_.RBSE[I-1] = atof(PCH);
		PCH = strtok(NULL, " ");
		*CIMDET_.LLE = atoi(PCH);
		fprintf(IWDUMP, "		 %d		%.10f		%.10f		%.10f		%.10E		%d\n", CIMDET_.NE[I-1], CIMDET_.EL[I-1], CIMDET_.EU[I-1], CIMDET_.BSE[I-1], CIMDET_.RBSE[I-1], CIMDET_.LLE[I-1]);

		getline(IRD, line);
		sArray = (char *) realloc(NULL, (line.length()+1) * sizeof(char));
		strcpy(sArray, line.c_str());
		PCH = strtok(sArray, " ");
		for (int J = 1; I <= CIMDET_.NE[I-1]; I++){
			CIMDET_.DIT[J-1][I-1] = atof(PCH);
			fprintf(IWDUMP, "		%.10f", CIMDET_.DIT[J-1][I-1]);
			PCH = strtok(NULL, " ");
		}
		fprintf(IWDUMP, "\n");

		getline(IRD, line);
		sArray = (char *) realloc(NULL, (line.length()+1) * sizeof(char));
		strcpy(sArray, line.c_str());
		PCH = strtok(sArray, " ");
		for (int J = 1; I <= CIMDET_.NE[I-1]; I++){
			CIMDET_.DIT2[J-1][I-1] = atof(PCH);
			fprintf(IWDUMP, "		%.10f", CIMDET_.DIT2[J-1][I-1]);
			PCH = strtok(NULL, " ");
		}
		fprintf(IWDUMP, "\n");

		getline(IRD, line);
		sArray = (char *) realloc(NULL, (line.length()+1) * sizeof(char));
		strcpy(sArray, line.c_str());
		PCH = strtok(sArray, " ");
		for (int J = 1; I <= CIMDET_.NE[I-1]; I++){
			for (int K = 1; K <= 3; K++){
				CIMDET_.DIP[K-1][J-1][I-1] = atof(PCH);
				fprintf(IWDUMP, "		%.10f", CIMDET_.DIP[K-1][J-1][I-1]);
				PCH = strtok(NULL, " ");
			}
		}
		fprintf(IWDUMP, "\n");

		getline(IRD, line);
		sArray = (char *) realloc(NULL, (line.length()+1) * sizeof(char));
		strcpy(sArray, line.c_str());
		PCH = strtok(sArray, " ");
		for (int J = 1; I <= CIMDET_.NE[I-1]; I++){
			for (int K = 1; K <= 3; K++){
				CIMDET_.DIP2[K-1][J-1][I-1] = atof(PCH);
				fprintf(IWDUMP, "		%.10f", CIMDET_.DIP2[K-1][J-1][I-1]);
				PCH = strtok(NULL, " ");
			}
		}
		fprintf(IWDUMP, "\n");

		if (CIMDET_.IDCUT[I-1] == 2){
			getline(IRD, line);
			sArray = (char *) realloc(NULL, (line.length()+1) * sizeof(char));
			strcpy(sArray, line.c_str());
			PCH = strtok(sArray, " ");
			strcpy(SPCFLO[I-1], PCH);
			fprintf(IWDUMP, "	%s\n", SPCFLO[I-1]);

			getline(IRD, line);
			sArray = (char *) realloc(NULL, (line.length()+1) * sizeof(char));
			strcpy(sArray, line.c_str());
			PCH = strtok(sArray, " ");
			for (int J = 1; I <= CIMDET_.NE[I-1]; I++){
				CIMDET_.FLT[J-1][I-1] = atof(PCH);
				fprintf(IWDUMP, "		%.10f", CIMDET_.FLT[J-1][I-1]);
				PCH = strtok(NULL, " ");
			}
			fprintf(IWDUMP, "\n");

			getline(IRD, line);
			sArray = (char *) realloc(NULL, (line.length()+1) * sizeof(char));
			strcpy(sArray, line.c_str());
			PCH = strtok(sArray, " ");
			for (int J = 1; I <= CIMDET_.NE[I-1]; I++){
				CIMDET_.FLT2[J-1][I-1] = atof(PCH);
				fprintf(IWDUMP, "		%.10f", CIMDET_.FLT2[J-1][I-1]);
				PCH = strtok(NULL, " ");
			}
			fprintf(IWDUMP, "\n");

			getline(IRD, line);
			sArray = (char *) realloc(NULL, (line.length()+1) * sizeof(char));
			strcpy(sArray, line.c_str());
			PCH = strtok(sArray, " ");
			for (int J = 1; I <= CIMDET_.NE[I-1]; I++){
				for (int K = 1; K <= 3; K++){
					CIMDET_.FLP[K-1][J-1][I-1] = atof(PCH);
					fprintf(IWDUMP, "		%.10f", CIMDET_.FLP[K-1][J-1][I-1]);
					PCH = strtok(NULL, " ");
				}
			}
			fprintf(IWDUMP, "\n");

			getline(IRD, line);
			sArray = (char *) realloc(NULL, (line.length()+1) * sizeof(char));
			strcpy(sArray, line.c_str());
			PCH = strtok(sArray, " ");
			for (int J = 1; I <= CIMDET_.NE[I-1]; I++){
				for (int K = 1; K <= 3; K++){
					CIMDET_.FLP2[K-1][J-1][I-1] = atof(PCH);
					fprintf(IWDUMP, "		%.10f", CIMDET_.FLP2[K-1][J-1][I-1]);
					PCH = strtok(NULL, " ");
				}
			}
			fprintf(IWDUMP, "\n");
		}
		

		getline(IRD, line);
		sArray = (char *) realloc(NULL, (line.length()+1) * sizeof(char));
		strcpy(sArray, line.c_str());
		PCH = strtok(sArray, " ");
		CIMDET_.NAGE[I-1] = atoi(PCH);
		PCH = strtok(NULL, " ");
		CIMDET_.AGEL[I-1] = atof(PCH);
		PCH = strtok(NULL, " ");
		CIMDET_.AGEU[I-1] = atof(PCH);
		PCH = strtok(NULL, " ");
		CIMDET_.BAGE[I-1] = atof(PCH);
		PCH = strtok(NULL, " ");
		CIMDET_.RBAGE[I-1] = atof(PCH);
		PCH = strtok(NULL, " ");
		*CIMDET_.LLAGE = atoi(PCH);
		fprintf(IWDUMP, "		 %d		%.10f		%.10f		%.10f		%.10E		%d\n", CIMDET_.NAGE[I-1], CIMDET_.AGEL[I-1], CIMDET_.AGEU[I-1], CIMDET_.BAGE[I-1], CIMDET_.RBAGE[I-1], CIMDET_.LLAGE[I-1]);

		if (CIMDET_.NAGE[I-1] > 0){
			getline(IRD, line);
			sArray = (char *) realloc(NULL, (line.length()+1) * sizeof(char));
			strcpy(sArray, line.c_str());
			PCH = strtok(sArray, " ");
			strcpy(SPCAGE[I-1], PCH);
			fprintf(IWDUMP, "	%s\n", SPCAGE[I-1]);

			getline(IRD, line);
			sArray = (char *) realloc(NULL, (line.length()+1) * sizeof(char));
			strcpy(sArray, line.c_str());
			PCH = strtok(sArray, " ");
			for (int J = 1; I <= CIMDET_.NAGE[I-1]; I++){
				CIMDET_.AGE[J-1][I-1] = atof(PCH);
				fprintf(IWDUMP, "		%.10f", CIMDET_.AGE[J-1][I-1]);
				PCH = strtok(NULL, " ");
			}
			fprintf(IWDUMP, "\n");

			getline(IRD, line);
			sArray = (char *) realloc(NULL, (line.length()+1) * sizeof(char));
			strcpy(sArray, line.c_str());
			PCH = strtok(sArray, " ");
			for (int J = 1; I <= CIMDET_.NAGE[I-1]; I++){
				CIMDET_.AGE2[J-1][I-1] = atof(PCH);
				fprintf(IWDUMP, "		%.10f", CIMDET_.AGE2[J-1][I-1]);
				PCH = strtok(NULL, " ");
			}
			fprintf(IWDUMP, "\n");
		}
	}
	for (int I = 1; I <= *CIMDET_.NID; I++){
		CIMDET_.EDEPP[I-1]=0.0e0;
        CIMDET_.LEDEP[I-1]=0;
		for (int J = 1; CIMDET_.NE[I-1]; J++){
			CIMDET_.DITP[J-1][I-1]=0.0e0;
            CIMDET_.LDIT[J-1][I-1]=0;
			for (int K = 1; K <= 3; K++){
				CIMDET_.DIPP[K-1][J-1][I-1]=0.0e0;
              	CIMDET_.LDIP[K-1][J-1][I-1]=0;
			}
		}
	}

	for (int I = 1; *CIMDET_.NID; I++){
		for (int J = 1; CIMDET_.NE[I-1]; J++){
			CIMDET_.FLTP[J-1][I-1]=0.0e0;
            CIMDET_.LFLT[J-1][I-1]=0;
			for (int K = 1; K <= 3; K++){
				 CIMDET_.FLPP[K-1][J-1][I-1]=0.0e0;
                 CIMDET_.LFLP[K-1][J-1][I-1]=0;
			}
		}
	}

	for (int I = 1; *CIMDET_.NID; I++){
		for (int J = 1; CIMDET_.NAGE[I-1]; J++){
			CIMDET_.AGEP[J-1][I-1]=0.0e0;
            CIMDET_.LAGEA[J-1][I-1]=0;
		}
	}
}

void endetr2_(ifstream &IRD, FILE *IWDUMP){

	char *PCH;
	string line;


	getline(IRD, line);
	char *sArray = (char *) malloc((line.length()+1) * sizeof(char));
	strcpy(sArray, line.c_str());
	PCH = strtok(sArray, " ");
	*CENDET_.NID = atoi(PCH);
	fprintf(IWDUMP, "	%d\n", *CENDET_.NID);

	getline(IRD, line);
	sArray = (char *) realloc(NULL, (line.length()+1) * sizeof(char));
	strcpy(sArray, line.c_str());
	PCH = strtok(sArray, " ");
	for (int I = 1; I <= *CENDET_.NID; I++){
		CENDET_.EDEP[I-1] = atof(PCH);
		fprintf(IWDUMP, "		%f", CENDET_.EDEP[I-1]);
		PCH = strtok(NULL, " ");
	}
	fprintf(IWDUMP, "\n");

	getline(IRD, line);
	sArray = (char *) realloc(NULL, (line.length()+1) * sizeof(char));
	strcpy(sArray, line.c_str());
	PCH = strtok(sArray, " ");
	for (int I = 1; I <= *CENDET_.NID; I++){
		CENDET_.EDEP2[I-1] = atof(PCH);
		fprintf(IWDUMP, "		%f", CENDET_.EDEP2[I-1]);
		PCH = strtok(NULL, " ");
	}
	fprintf(IWDUMP, "\n");

	for (int I = 1; I <= *CENDET_.NID; I++){
		getline(IRD, line);
		sArray = (char *) realloc(NULL, (line.length()+1) * sizeof(char));
		strcpy(sArray, line.c_str());
		PCH = strtok(sArray, " ");
		strcpy(SPCDEO[I-1], PCH);
		fprintf(IWDUMP, "	%s\n", SPCDEO[I-1]);

		getline(IRD, line);
		char *sArray = (char *) malloc((line.length()+1) * sizeof(char));
		strcpy(sArray, line.c_str());
		PCH = strtok(sArray, " ");
		CENDET_.NE[I-1] = atoi(PCH);
		PCH = strtok(NULL, " ");
		CENDET_.EL[I-1] = atof(PCH);
		PCH = strtok(NULL, " ");
		CENDET_.EU[I-1] = atof(PCH);
		PCH = strtok(NULL, " ");
		CENDET_.BSE[I-1] = atof(PCH);
		PCH = strtok(NULL, " ");
		CENDET_.RBSE[I-1] = atof(PCH);
		PCH = strtok(NULL, " ");
		*CENDET_.LLE = atoi(PCH);
		fprintf(IWDUMP, "		 %d		%.10f		%.10f		%.10f		%.10E		%d\n", CENDET_.NE[I-1], CENDET_.EL[I-1], CENDET_.EU[I-1], CENDET_.BSE[I-1], CENDET_.RBSE[I-1], CENDET_.LLE[I-1]);

		getline(IRD, line);
		sArray = (char *) realloc(NULL, (line.length()+1) * sizeof(char));
		strcpy(sArray, line.c_str());
		PCH = strtok(sArray, " ");
		for (int J = 1; I <= CENDET_.NE[I-1]; I++){
			CENDET_.DET[J-1][I-1] = atof(PCH);
			fprintf(IWDUMP, "		%.10f", CENDET_.DET[J-1][I-1]);
			PCH = strtok(NULL, " ");
		}
		fprintf(IWDUMP, "\n");
	}
}

void doser2_(ifstream &IRD, FILE *IWDUMP){

	char *PCH;
	string line;


	getline(IRD, line);
	char *sArray = (char *) malloc((line.length()+1) * sizeof(char));
	strcpy(sArray, line.c_str());
	PCH = strtok(sArray, " ");
	for (int I = 1; I <= 3; I++){
		CDOSE3_.NDB[I-1] = atoi(PCH);
		fprintf(IWDUMP, "		%d", CDOSE3_.NDB[I-1]);
		PCH = strtok(NULL, " ");
	}
	fprintf(IWDUMP, "\n");

	getline(IRD, line);
	sArray = (char *) realloc(NULL, (line.length()+1) * sizeof(char));
	strcpy(sArray, line.c_str());
	PCH = strtok(sArray, " ");
	for (int I = 1; I <= 3; I++){
		CDOSE3_.DXL[I-1] = atof(PCH);
		fprintf(IWDUMP, "		%.10f", CDOSE3_.DXL[I-1]);
		PCH = strtok(NULL, " ");
	}
	fprintf(IWDUMP, "\n");

	getline(IRD, line);
	sArray = (char *) realloc(NULL, (line.length()+1) * sizeof(char));
	strcpy(sArray, line.c_str());
	PCH = strtok(sArray, " ");
	for (int I = 1; I <= 3; I++){
		CDOSE3_.DXU[I-1] = atof(PCH);
		fprintf(IWDUMP, "		%.10f", CDOSE3_.DXU[I-1]);
		PCH = strtok(NULL, " ");
	}
	fprintf(IWDUMP, "\n");

	getline(IRD, line);
	sArray = (char *) realloc(NULL, (line.length()+1) * sizeof(char));
	strcpy(sArray, line.c_str());
	PCH = strtok(sArray, " ");
	for (int I = 1; I <= 3; I++){
		CDOSE3_.BDOSE[I-1] = atof(PCH);
		fprintf(IWDUMP, "		%.10f", CDOSE3_.BDOSE[I-1]);
		PCH = strtok(NULL, " ");
	}
	for (int I = 1; I <= 3; I++){
		CDOSE3_.RBDOSE[I-1] = atof(PCH);
		fprintf(IWDUMP, "		%.10f", CDOSE3_.RBDOSE[I-1]);
		PCH = strtok(NULL, " ");
	}
	fprintf(IWDUMP, "\n");

	getline(IRD, line);
	sArray = (char *) realloc(NULL, (line.length()+1) * sizeof(char));
	strcpy(sArray, line.c_str());
	PCH = strtok(sArray, " ");
	*CDOSE1_.KDOSE = atoi(PCH);
	fprintf(IWDUMP, "		%.d\n", *CDOSE1_.KDOSE);

	getline(IRD, line);
	sArray = (char *) realloc(NULL, (line.length()+1) * sizeof(char));
	strcpy(sArray, line.c_str());
	PCH = strtok(sArray, " ");
	for (int I3 = 1; I3 <= CDOSE3_.NDB[3-1]; I3++){
		for (int I2 = 1; I2 <= CDOSE3_.NDB[2-1]; I2++){
			for (int I1 = 1; I1 <= CDOSE3_.NDB[1-1]; I1++){
				CDOSE4_.VMASS[I3-1][I2-1][I1-1] = atof(PCH);
				fprintf(IWDUMP, "		%.10f", CDOSE4_.VMASS[I3-1][I2-1][I1-1]);
				PCH = strtok(NULL, " ");
			}
		}
	}
	fprintf(IWDUMP, "\n");

	getline(IRD, line);
	sArray = (char *) realloc(NULL, (line.length()+1) * sizeof(char));
	strcpy(sArray, line.c_str());
	PCH = strtok(sArray, " ");
	for (int I3 = 1; I3 <= CDOSE3_.NDB[3-1]; I3++){
		for (int I2 = 1; I2 <= CDOSE3_.NDB[2-1]; I2++){
			for (int I1 = 1; I1 <= CDOSE3_.NDB[1-1]; I1++){
				CDOSE1_.DOSE[I3-1][I2-1][I1-1] = atof(PCH);
				fprintf(IWDUMP, "		%.10f", CDOSE1_.DOSE[I3-1][I2-1][I1-1]);
				PCH = strtok(NULL, " ");
			}
		}
	}
	fprintf(IWDUMP, "\n");

	getline(IRD, line);
	sArray = (char *) realloc(NULL, (line.length()+1) * sizeof(char));
	strcpy(sArray, line.c_str());
	PCH = strtok(sArray, " ");
	for (int I3 = 1; I3 <= CDOSE3_.NDB[3-1]; I3++){
		for (int I2 = 1; I2 <= CDOSE3_.NDB[2-1]; I2++){
			for (int I1 = 1; I1 <= CDOSE3_.NDB[1-1]; I1++){
				CDOSE1_.DOSE2[I3-1][I2-1][I1-1] = atof(PCH);
				fprintf(IWDUMP, "		%.10f", CDOSE1_.DOSE2[I3-1][I2-1][I1-1]);
				PCH = strtok(NULL, " ");
			}
		}
	}
	fprintf(IWDUMP, "\n");

	getline(IRD, line);
	sArray = (char *) realloc(NULL, (line.length()+1) * sizeof(char));
	strcpy(sArray, line.c_str());
	PCH = strtok(sArray, " ");
	for (int I3 = 1; I3 <= CDOSE3_.NDB[3-1]; I3++){
		CDOSE2_.DDOSE[I3-1] = atof(PCH);
		fprintf(IWDUMP, "		%.10f", CDOSE2_.DDOSE[I3-1]);
		PCH = strtok(NULL, " ");
	}
	fprintf(IWDUMP, "\n");

	getline(IRD, line);
	sArray = (char *) realloc(NULL, (line.length()+1) * sizeof(char));
	strcpy(sArray, line.c_str());
	PCH = strtok(sArray, " ");
	for (int I3 = 1; I3 <= CDOSE3_.NDB[3-1]; I3++){
		CDOSE2_.DDOSE2[I3-1] = atof(PCH);
		fprintf(IWDUMP, "		%.10f", CDOSE2_.DDOSE2[I3-1]);
		PCH = strtok(NULL, " ");
	}
	fprintf(IWDUMP, "\n");

}




void shower2_(){

	//Simula uma nova particula e registra as quantidades relevantes.
	//printf("shower2\n");

/*		if (imprimiu==0){
		printf("\n\nSHOWER2\n\n");
		imprimiu++;
	}*/

	bool LINTF;

	double REV=5.10998928e5;
	double TREV=2.0e0*REV;
	double PI=3.1415926535897932e0;
	double TWOPI=2.0e0*PI;

	int IEXIT, METAST, NTRIAL, K, KEn, IBODYL, NCROSS, IDET, NSHJ, MATL, ICOL, LEFT;
	double RN, RNF, UV, PHI, DSEF, DS, DEP, XL, YL, ZL, DECSD, DSEFR, XD, YD, ZD, DE, WS, US, VS, SDTS, DF;

	//A simulação da particula começa aqui.

	//Contadores de particulas primarias

L101:;
	for (int I = 1; I <= 3; I++){
		CNT0_.DPRIM[I-1]=0.0e0;
		for (int K = 1; K <= 3; K++){
			CNT0_.DSEC[I-1][K-1]=0.0e0;
		}
	}

	for (int I = 1; I <= 2; I++){
		CNT0_.DAVW[I-1]=0.0e0;
        CNT0_.DAVA[I-1]=0.0e0;
        CNT0_.DAVE[I-1]=0.0e0;
	}

	for (int KB = 1; KB <= *PENGEOM_mod_.NBODY; KB++){
		CNT1_.DEBO[KB-1]=0.0e0; //Energias depositadas nos diversos corpos KB
	}

	IEXIT=0;
    METAST=0;

	cleans2_(); //Limpa a pilha secundaria

	//if (*TRACK_MOD_.LAGE)
	//	PAGE0; //define o tempo de voo igual a 0
	//Não será realizado a contabilização do tempo individual de vida de uma particula nesta versão do programa


	//Definindo o estado inicial da partícula primária.

L201:;

	if (*CSOUR0_.KPARP == 0){ //não será implementado o tratamento para eletrons como particula primaria nesta versão do programas.
		printf("Nao sera implementado simulação com particula primaria sendo eletron, apenas fotons\n");

	}else{
		//Fonte Externa
		*CNTRL_.SHN= *CNTRL_.SHN+1.0e0;
        *CNTRL_.N= *CNTRL_.N+1;

		if (*CNTRL_.N > 2000000000)
			*CNTRL_.N= *CNTRL_.N-2000000000;

		*TRACK_mod_.KPAR=*CSOUR0_.KPARP;
        *TRACK_mod_.WGHT=1.0e0;

		//Posição inicial da particula
		if (*CSOUR3_.LEXSRC){
			if (*CSOUR3_.LEXBD){
				NTRIAL=0;
L301:;			
				*TRACK_mod_.X= *CSOUR3_.SX0+(rand2_(1.0e0)-0.5e0)* *CSOUR3_.SSX;
            	*TRACK_mod_.Y= *CSOUR3_.SY0+(rand2_(2.0e0)-0.5e0)* *CSOUR3_.SSY;
            	*TRACK_mod_.Z= *CSOUR3_.SZ0+(rand2_(3.0e0)-0.5e0)* *CSOUR3_.SSZ;
            	locate2_();
				NTRIAL=NTRIAL+1;
				if (NTRIAL > 200){
					printf("   WARNING: the sampling of initial positions may be very inefficient.");
					//imrpimir em um arquivo tbm
				}

				if (CSOUR3_.IXSBOD[*TRACK_mod_.IBODY-1] == 0)
					goto L301;
			}else{
				*TRACK_mod_.X= *CSOUR3_.SX0+(rand2_(1.0e0)-0.5e0)* *CSOUR3_.SSX;
            	*TRACK_mod_.Y= *CSOUR3_.SY0+(rand2_(2.0e0)-0.5e0)* *CSOUR3_.SSY;
            	*TRACK_mod_.Z= *CSOUR3_.SZ0+(rand2_(3.0e0)-0.5e0)* *CSOUR3_.SSZ;
			}
		} else{
			*TRACK_mod_.X= *CSOUR3_.SX0;
            *TRACK_mod_.Y= *CSOUR3_.SY0;
            *TRACK_mod_.Z= *CSOUR3_.SZ0;
		}

		//Direção Inicial
		if (*CSOUR0_.LSCONE){
			gcone2_(*TRACK_mod_.U, *TRACK_mod_.V, *TRACK_mod_.W); //Feixe Conico
			
		}else{ //Feixe Retangular

			*TRACK_mod_.W=*CSOUR0_.CTHL+rand2_(4.0e0)* *CSOUR0_.DCTH; 
            UV=sqrt(1.0e0-*TRACK_mod_.W * *TRACK_mod_.W);
            PHI=*CSOUR0_.PHIL+rand2_(5.0e0)* *CSOUR0_.DPHI;
            *TRACK_mod_.U=UV*cos(PHI);
            *TRACK_mod_.V=UV*sin(PHI);
		}
		//Energia Inicial
		if (*CSOUR2_.LSPEC){
			RN=rand2_(6.0e0)* *CNT2_.NSEB + 1;
            K=RN; //Espectro contínuo. E amostrado pelo método de Walker.
            RNF=RN-K;
			if (RNF > CSOUR2_.FSRC[K-1]){
				
				KEn=CSOUR2_.IASRC[K-1];
		
			}else{
				KEn=K;
			}
			*TRACK_mod_.E=CSOUR2_.ESRC[KEn-1]+rand2_(7.0e0)*(CSOUR2_.ESRC[KEn+1-1]-CSOUR2_.ESRC[KEn-1]);
          	CNT2_.SHIST[KEn-1]=CNT2_.SHIST[KEn-1]+1.0e0;
		}else{
			*TRACK_mod_.E= *CSOUR1_.E0; //Fonte MonoEnergetica.
         	 CNT2_.SHIST[1-1]=CNT2_.SHIST[1-1]+1.0e0;
		}

		TRACK_mod_.ILB[1-1]=1;  //Identifica partículas primárias.
        TRACK_mod_.ILB[2-1]=0;
        TRACK_mod_.ILB[3-1]=0;
        TRACK_mod_.ILB[4-1]=0;
        TRACK_mod_.ILB[5-1]=0;

		if (*TRACK_mod_.KPAR == 2){
			if (*CSOUR0_.LGPOL){
				*TRACK_mod_.IPOL=1;  //Polarizacao de Fotons
           	 	*TRACK_mod_.SP1=*CSOUR1_.SP10;
            	*TRACK_mod_.SP2=*CSOUR1_.SP20;
            	*TRACK_mod_.SP3=*CSOUR1_.SP30;
			}
			else{
				*TRACK_mod_.IPOL=0;
			}
		}else{
			*TRACK_mod_.IPOL=0;
		}
	}

	//Verifique se a trajetória cruza o sistema de materiais.
	//implementacao da simulacao

	//A partir daqui entra as implmentações em GPU.
L302:;

	locate2_();
	if (*TRACK_mod_.MAT == 0){
		IBODYL=*TRACK_mod_.IBODY;	
		DS = 1.0e30;
		step2_(&DS,&DSEF,&NCROSS);
		/*if (*TRACK_mod_.LAGE)
			DPAGE(DSEF,DSTOT)*/ //Funcao para contabilizar o tempo de vida de uma particula, não sera implementada
		
		if (*TRACK_mod_.MAT == 0){//A particula não entrou no sistema
			if (*TRACK_mod_.W > 0){
				IEXIT = 1; //Rotula partículas ascendentes emergentes.
			}else{
				IEXIT=2; //Rotula partículas descendentes emergentes.
			}
			goto L104;
		} 

		//Detetores de Impacto

		IDET=PENGEOM_mod_.KDET[*TRACK_mod_.IBODY-1];
		if (IDET != 0){
			if ((PENGEOM_mod_.KDET[IBODYL-1] != IDET)  && (CNT4_.KKDI[*TRACK_mod_.KPAR-1][IDET-1] == 1)){
				//Esse trecho faz gravação no arquivo Phase-Space que não será implementado nesta versao
				/*if (CNT4_.IPSF[IDET-1] == 1){
					NSHJ=*CNTRL_.SHN - *CNT4_.RLAST;
					wrpsf2_(CNT4_.IPSFO,NSHJ,0);
					*CNT4_.RWRITE=*CNT4_.RWRITE+1.0e0;
					*CNT4_.RLAST=*CNTRL_.SHN;*/

				simdet2_(*CNTRL_.N,IDET);

				if (CNT4_.IDCUT[IDET-1] == 0){
					CNT1_.DEBO[*TRACK_mod_.IBODY-1]=CNT1_.DEBO[*TRACK_mod_.IBODY-1]+*TRACK_mod_.E * *TRACK_mod_.WGHT;
					IEXIT=3;
					goto L104;
				}
			}
		}

	}

	//Aniquiliação de positron quando a energia da particula é muito pequena

	if (*TRACK_mod_.E < CSPGEO_.EABSB[*TRACK_mod_.IBODY-1][*TRACK_mod_.KPAR-1]){ //energia é muito baixa
		DEP= *TRACK_mod_.E * *TRACK_mod_.WGHT;
		if ((*TRACK_mod_.KPAR == 3) && (*TRACK_mod_.E > 1.0e-6)){ //aniquilação de positrion
			panar2_(CSPGEO_.EABSB[*TRACK_mod_.IBODY-1][2-1]);
          	DEP= DEP + TREV * *TRACK_mod_.WGHT;
		}
		*TRACK_mod_.E=0.0e0;

		CNT1_.DEBO[*TRACK_mod_.IBODY-1]=CNT1_.DEBO[*TRACK_mod_.IBODY-1]+DEP;
		if (*CNT6_.LDOSEM){
			sdose2_(DEP,*TRACK_mod_.X,*TRACK_mod_.Y,*TRACK_mod_.Z,*TRACK_mod_.MAT,*CNTRL_.N);
		}
		IEXIT=3;
		goto L104;

	}

	//Simulação da historia da particula inicia aqui.

	//Divisão de particulas e roleta russa . Apenas para particulas lidas de um arquivoco Phase-Space.
	//Essa versão não irá implementar reduções de variancia e leitura de arquivo Phase-Space.

L102:;

	if (*TRACK_mod_.E < CSPGEO_.EABSB[*TRACK_mod_.IBODY-1][*TRACK_mod_.KPAR-1]){ //energia é muito baixa
		DEP= *TRACK_mod_.E * *TRACK_mod_.WGHT;
		if ((*TRACK_mod_.KPAR == 3) && (*TRACK_mod_.E > 1.0e-6)){ //aniquilação de positrion
			panar2_(CSPGEO_.EABSB[*TRACK_mod_.IBODY-1][2-1]);
          	DEP= DEP + TREV* *TRACK_mod_.WGHT;
		}
		*TRACK_mod_.E=0.0e0;

		//A energia é depositada localmente no material.
		CNT1_.DEBO[*TRACK_mod_.IBODY-1]=CNT1_.DEBO[*TRACK_mod_.IBODY-1]+DEP;
		if (*CNT6_.LDOSEM){ //Partícula dentro da caixa de dose.
			sdose2_(DEP,*TRACK_mod_.X,*TRACK_mod_.Y,*TRACK_mod_.Z,*TRACK_mod_.MAT,*CNTRL_.N);
		}
		IEXIT=3; //Marca partículas absorvidas.
		goto L104; //saida
	}

	start2_(); //Inicia a simulação no meio atual.

	//Comprimento do caminho livre para o próximo evento de interação.
L103:;
	IBODYL=*TRACK_mod_.IBODY; 
	MATL=*TRACK_mod_.MAT; 
	XL=*TRACK_mod_.X; 
	YL=*TRACK_mod_.Y; 
	ZL=*TRACK_mod_.Z;

	/*if ((CFORCI_.LFORCE[KPAR-1][IBODY-1]) && ((WGHT > WLOW[KPAR-1][IBODY-1]) && (WGHT < WHIG[KPAR-1][IBODY-1]))){
		//JUMPF(DSMAX(IBODY),DS)  //Força de Interação
        //LINTF=true;
		//Esse trecho nao sera implementado pois a versão nao irá tratar redução de variancia.
	}else`{
		jumtp2(DSMAX(IBODY),DS)  ! Analogue simulation.
        LINTF=false
	}*/

	jump2_(CSPGEO_.DSMAX[*TRACK_mod_.IBODY-1], DS);
    LINTF=false;

	step2_(&DS,&DSEF,&NCROSS); //Determina a posição final do passo.

	//Distribuição de energia de fluência.
	IDET=PENGEOM_mod_.KDET[IBODYL-1];
	if (IDET != 0){
		if (CNT4_.IDCUT[IDET-1] == 2){
			if (CNT4_.KKDI[*TRACK_mod_.KPAR-1][IDET-1] == 1){
				if (*TRACK_mod_.KPAR == 2){
					fimdet2_(*CNTRL_.N,IDET,DSEF);
				}else{
					DECSD=*PENELOPE_mod_.SSOFT*DSEF;
					if (DECSD > 1.0e-12){ // A distribuição pode se estender
						fimdes2_(*CNTRL_.N,IDET,*PENELOPE_mod_.E0STEP,DECSD,DSEF); //abaixo do EABS.
					}else{
						fimdet2_(*CNTRL_.N,IDET,DSEF);
					}
				}
			}
		}
	}

	//A partícula cruzou uma interface.

	if (NCROSS > 0){
		//Correção da perda de energia suave (CSDA).
		if (*TRACK_mod_.KPAR != 2){
			*TRACK_mod_.E=*PENELOPE_mod_.E0STEP-*PENELOPE_mod_.SSOFT*DSEF;
			if (*CJUMP1_.MHINGE == 0){
				DEP=*PENELOPE_mod_.SSOFT*DSEF* *TRACK_mod_.WGHT;
				if (*CNT6_.LDOSEM){
					DSEFR=rand2_(8.0e0)*DSEF;
              		XD=XL+*TRACK_mod_.U*DSEFR;
              		YD=YL+*TRACK_mod_.V*DSEFR;
              		ZD=ZL+*TRACK_mod_.W*DSEFR;
					sdose2_(DEP,XD,YD,ZD,MATL,*CNTRL_.N);
				}
			}else{
				DEP=-*PENELOPE_mod_.SSOFT*(DS-DSEF)* *TRACK_mod_.WGHT;
				if (*CNT6_.LDOSEM){
					sdose2_(DEP,XL,YL,ZL,MATL,*CNTRL_.N);
				}
			}
			CNT1_.DEBO[IBODYL-1]=CNT1_.DEBO[IBODYL-1]+DEP;
		}
		//Verifique se a partícula está fora do invólucro.
		if (*TRACK_mod_.MAT == 0){ //A partícula está fora do recinto.
			if (*TRACK_mod_.W > 0.0e0){
				IEXIT=1; //Marca partículas ascendentes emergentes.
			}else{
				IEXIT=2; //Rotula partículas descendentes emergentes.
			}
			goto L104; //Saida
		}

		//Detectores de Impacto

		IDET=PENGEOM_mod_.KDET[*TRACK_mod_.IBODY-1];
		if (IDET != 0){
			if ((PENGEOM_mod_.KDET[IBODYL-1] != IDET)  && (CNT4_.KKDI[*TRACK_mod_.KPAR-1][IDET-1] == 1)){
				//Esse trecho faz gravação no arquivo Phase-Space que não será implementado nesta versao
				/*if (CNT4_.IPSF[IDET-1] == 1){
					NSHJ=*CNTRL_.SHN - *CNT4_.RLAST;
					wrpsf2_(CNT4_.IPSFO,NSHJ,0);
					*CNT4_.RWRITE=*CNT4_.RWRITE+1.0e0;
					*CNT4_.RLAST=*CNTRL_.SHN;*/

				simdet2_(*CNTRL_.N,IDET);

				if (CNT4_.IDCUT[IDET-1] == 0){
					CNT1_.DEBO[*TRACK_mod_.IBODY-1]=CNT1_.DEBO[*TRACK_mod_.IBODY-1]+*TRACK_mod_.E * *TRACK_mod_.WGHT;
					IEXIT=3;
					goto L104;
				}
			}
		}
		goto L102;
	}

	//Simulação dos eventos de interação da particula com a materia.
	if (LINTF){
		// knockf(DE, ICOL) // Não será implementado redução de variancia nesta versão
	}else{
		knock2_(DE, ICOL);
	}

	if (*TRACK_mod_.E < CSPGEO_.EABSB[*TRACK_mod_.IBODY-1][*TRACK_mod_.KPAR-1]){ //A partícula foi absorvida.
		DE=DE+ *TRACK_mod_.E;
		if ((*TRACK_mod_.KPAR == 3) && (*TRACK_mod_.E > 1.0e-6)){ //Aniquilação de positron
			panar2_(CSPGEO_.EABSB[*TRACK_mod_.IBODY-1][2-1]); // Quando Absorvida
			DE=DE+TREV;
		}
		*TRACK_mod_.E=0.0e0;
	}

	DEP=DE* *TRACK_mod_.WGHT;
    CNT1_.DEBO[*TRACK_mod_.IBODY-1]=CNT1_.DEBO[*TRACK_mod_.IBODY-1]+DEP;

	if (*CNT6_.LDOSEM)
		sdose2_(DEP,*TRACK_mod_.X,*TRACK_mod_.Y,*TRACK_mod_.Z,*TRACK_mod_.MAT,*CNTRL_.N);

	if (*TRACK_mod_.E < CSPGEO_.EABSB[*TRACK_mod_.IBODY-1][*TRACK_mod_.KPAR-1]){ ////A partícula foi absorvida.
		IEXIT=3; //Marca partículas absorvidas.
		goto L104; //saida
	}

	goto L103;

	//A simulação da particula termina aqui.

L104:;

	//Incrementar contadores de partículas.
	if (TRACK_mod_.ILB[1-1] == 1){
		CNT0_.DPRIM[IEXIT-1]=CNT0_.DPRIM[IEXIT-1]+ *TRACK_mod_.WGHT;
		/*if (*CSOUR0_.LPSF){  não sera implementado redução de variancia como divisão de particulas ou roleta russa
			if (NSPL1 > 1)
				CNT0_.DPRIM[IEXIT-1]=CNT0_.DPRIM[IEXIT-1]+*TRACK_mod_.WGHT*(NSPL1-1);
		}
		*/

		if (IEXIT < 3){
			CNT0_.DAVW[IEXIT-1]=CNT0_.DAVW[IEXIT-1]+*TRACK_mod_.W* *TRACK_mod_.WGHT;
            CNT0_.DAVA[IEXIT-1]=CNT0_.DAVA[IEXIT-1]+acos(*TRACK_mod_.W)* *TRACK_mod_.WGHT;
            CNT0_.DAVE[IEXIT-1]=CNT0_.DAVE[IEXIT-1]+*TRACK_mod_.E* *TRACK_mod_.WGHT;
		}

	}else{
		CNT0_.DSEC[IEXIT-1][*TRACK_mod_.KPAR-1]=CNT0_.DSEC[IEXIT-1][*TRACK_mod_.KPAR-1]+ *TRACK_mod_.WGHT;
	}

	if (IEXIT < 3){
		tenang2_(IEXIT,*CNTRL_.N);  //Energia e distribuições angulares.
	}


	//Particulas Secundarias

L202:;


	secpar2_(LEFT);
	if (LEFT > 0){
		if (TRACK_mod_.ILB[1-1] == 1){ //Fonte da particula primaria
			KEn=*TRACK_mod_.E* *CNT3_.RDSDE+1.0e0;
          	CNT3_.SEDS[KEn-1][*TRACK_mod_.KPAR-1]=CNT3_.SEDS[KEn-1][*TRACK_mod_.KPAR-1]+*TRACK_mod_.WGHT;
          	CNT3_.SEDS2[KEn-1][*TRACK_mod_.KPAR-1]=CNT3_.SEDS2[KEn-1][*TRACK_mod_.KPAR-1]+pow(*TRACK_mod_.WGHT,2);
			//if (*TRACK_mod_.LAGE) //nao irá contabilizar a idade da particula
			//	page02_();
			goto L302; //A energia não é removida do local.
		}
		if (*TRACK_mod_.E > CSPGEO_.EABSB[*TRACK_mod_.IBODY-1][*TRACK_mod_.KPAR-1]){
			//Divisoes de raios X
			if (*TRACK_mod_.KPAR == 2){
				if (TRACK_mod_.ILB[4-1] > 0){ //caracteristica de raio X
					if ((TRACK_mod_.ILB[1-1] == 2) && (TRACK_mod_.ILB[3-1] < 9)){
						if (CXRSPL_.LXRSPL[*TRACK_mod_.IBODY-1]){
							*TRACK_mod_.WGHT=*TRACK_mod_.WGHT/(CXRSPL_.IXRSPL[*TRACK_mod_.IBODY-1]);
							CXRSPL_.ILBA[1-1]=TRACK_mod_.ILB[1-1];
							CXRSPL_.ILBA[2-1]=TRACK_mod_.ILB[2-1];
							CXRSPL_.ILBA[3-1]=9;
							CXRSPL_.ILBA[4-1]=TRACK_mod_.ILB[4-1];
							CXRSPL_.ILBA[5-1]=TRACK_mod_.ILB[5-1];
							for (int I = 2; I <= CXRSPL_.IXRSPL[*TRACK_mod_.IBODY-1]; I++){
								WS=-1.0e0+2.0e0*rand2_(9.0e0);
								SDTS=sqrt(1.0e0-WS*WS);
								DF=TWOPI*rand2_(10.0e0);
								US=cos(DF)*SDTS;
								VS=sin(DF)*SDTS;
								stores2_(*TRACK_mod_.E,*TRACK_mod_.X,*TRACK_mod_.Y,*TRACK_mod_.Z,US,VS,WS,*TRACK_mod_.WGHT,*TRACK_mod_.KPAR,CXRSPL_.ILBA,wIPOLI);
							}
						}
					}

				}
			}
			DEP=*TRACK_mod_.E* *TRACK_mod_.WGHT;
          	CNT1_.DEBO[*TRACK_mod_.IBODY-1]=CNT1_.DEBO[*TRACK_mod_.IBODY-1]-DEP; 
			if (*CNT6_.LDOSEM){
				double wDEP = -DEP;
				sdose2_(wDEP,*TRACK_mod_.X,*TRACK_mod_.Y,*TRACK_mod_.Z,*TRACK_mod_.MAT,*CNTRL_.N);
			}
		}else{
			goto L202;
		}
		goto L102;
	}


	/*if (*CSOUR0_.LPSF){ Não sera implementado o Phase-Space
		if (ISEC == 1)
			goto L201;

	}*/

	/*Energias depositadas em diferentes corpos e detectores.
	Calculando os espectros dos detectores de deposição de energia.*/

	if (*CNT5_.NED > 0){
		for (int KD = 1; KD <= *CNT5_.NED; KD++){
			CNT5_.DEDE[KD-1]=0.0e0;
		}
		for (int KB=1; KB <= *PENGEOM_mod_.NBODY; KB++){
			IDET=CNT5_.KBDE[KB-1];
			if (IDET != 0)
				CNT5_.DEDE[IDET-1]=CNT5_.DEDE[IDET-1]+CNT1_.DEBO[KB-1];

		}

		for (int IDET=1; IDET <= *CNT5_.NED; IDET++){
			sendet2_(CNT5_.DEDE[IDET-1],IDET);
		}
	}

	
	for (int KB = 1; KB <= *PENGEOM_mod_.NBODY; KB++){
		CNT1_.TDEBO[KB-1]=CNT1_.TDEBO[KB-1]+CNT1_.DEBO[KB-1];
        CNT1_.TDEBO2[KB-1]=CNT1_.TDEBO2[KB-1]+pow(CNT1_.DEBO[KB-1],2);
	//	printf("TDBEBO: %f\n", CNT1_.TDEBO[KB-1]);
	//	printf("TDBEBO2: %f\n\n", CNT1_.TDEBO2[KB-1]);

	}

	//Contadores de estado final
	for (int I=1; I <= 3; I++){
		CNT0_.PRIM[I-1]=CNT0_.PRIM[I-1]+CNT0_.DPRIM[I-1];
        CNT0_.PRIM2[I-1]=CNT0_.PRIM2[I-1]+pow(CNT0_.DPRIM[I-1],2);
		for (int K =1; K <= 3; K++){
			CNT0_.SEC[I-1][K-1]=CNT0_.SEC[I-1][K-1]+CNT0_.DSEC[I-1][K-1];
          	CNT0_.SEC2[I-1][K-1]=CNT0_.SEC2[I-1][K-1]+pow(CNT0_.DSEC[I-1][K-1],2);
		}
	}

	for (int I=1; I <= 2; I++){
		CNT0_.AVW[I-1]=CNT0_.AVW[I-1]+CNT0_.DAVW[I-1];
        CNT0_.AVW2[I-1]=CNT0_.AVW2[I-1]+pow(CNT0_.DAVW[I-1],2);
        CNT0_.AVA[I-1]=CNT0_.AVA[I-1]+CNT0_.DAVA[I-1];
        CNT0_.AVA2[I-1]=CNT0_.AVA2[I-1]+pow(CNT0_.DAVA[I-1],2);
        CNT0_.AVE[I-1]=CNT0_.AVE[I-1]+CNT0_.DAVE[I-1];
        CNT0_.AVE2[I-1]=CNT0_.AVE2[I-1]+pow(CNT0_.DAVE[I-1],2);
	}

}

void cleans2_(){

	//Esta sub-rotina inicializa a pilha secundária. Deve ser chamada antes de iniciar a simulação de cada pista primária.

	*SECST_.NSEC=0;

/*	if (imprimiu==0){
		printf("\n\ncleans2\n\n");
		imprimiu++;
	}*/

}

double rand2_(double DUMMY){ //gerador de numeros aleatorios
	/*
	Esta é uma versão adaptada da sub-rotina RANECU escrita por F. James
	(Comput. Phys. Commun. 60 (1990) 329-344), que foi modificado para
	dá um único número aleatório em cada chamada.

	As 'sementes' ISEED1 e ISEED2 devem ser inicializadas no programa principal
	e transferido através do bloco comum nomeado /RSEED/.

	Alguns compiladores incorporam um gerador de números aleatórios intrínseco com
	o mesmo nome (mas com diferentes listas de argumentos). Para evitar conflitos,
	é aconselhável declarar RAND como uma função externa em todas as sub-
	programas em que o chamam.
	*/
	/*if (imprimiu == 0){
		printf("\n\nrand2\n\n");
		imprimiu++;
	}*/
	double USCALE=1.0e0/2.147483563e9;
	int I1, I2, IZ;

	I1=*RSEED_.ISEED1/53668;
    *RSEED_.ISEED1=40014*(*RSEED_.ISEED1-I1*53668)-I1*12211;
	if (*RSEED_.ISEED1 < 0)
		*RSEED_.ISEED1=*RSEED_.ISEED1+2147483563;

	I2=*RSEED_.ISEED2/52774;
    *RSEED_.ISEED2=40692*(*RSEED_.ISEED2-I2*52774)-I2*3791;
    if (*RSEED_.ISEED2 < 0)
		 *RSEED_.ISEED2=*RSEED_.ISEED2+2147483399;

	IZ=*RSEED_.ISEED1-*RSEED_.ISEED2;
    if (IZ < 1) 
		IZ=IZ+2147483562;

	return IZ*USCALE;

}

void gcone2_(double &UF, double &VF, double &WF){

	/*
	Esta sub-rotina amostra uma direção aleatória uniformemente dentro de um cone
	 com eixo central na direção (THETA,PHI) e abertura ALPHA.
	Os parâmetros são inicializados chamando a sub-rotina GCONE0.
	*/

	double PI=3.1415926535897932e0;
	double TWOPI=2.0e0*PI;

	double WT,DF, SUV, UT, VT, DXY, DXYZ, FNORM;

	//Defina uma direção relativa ao eixo z.
	WT=*CGCONE_.CAPER+(1.0e0-*CGCONE_.CAPER)*rand2_(1.0e0);
    DF=TWOPI*rand2_(2.0e0);
    SUV=sqrt(1.0e0-WT*WT);
    UT=SUV*cos(DF);
    VT=SUV*sin(DF);


	//Rotacao para a direção do eixo do feixe

	UF=*CGCONE_.CPCT*UT-*CGCONE_.SPHI*VT+*CGCONE_.CPST*WT;
    VF=*CGCONE_.SPCT*UT+*CGCONE_.CPHI*VT+*CGCONE_.SPST*WT;
    WF=-*CGCONE_.STHE*UT+*CGCONE_.CTHE*WT;

	//Normalizaçao
	DXY=UF*UF+VF*VF;
    DXYZ=DXY+WF*WF;
	if (fabs(DXYZ-1.0e0) > 1.0e-14){
		FNORM=1.0e0/sqrt(DXYZ);
        UF=FNORM*UF;
        VF=FNORM*VF;
        WF=FNORM*WF;
	}
	/*if (imprimiu==0){
		printf("\n\ngcone2\n\n");
		imprimiu++;
	}*/
}

void simdet2_(int &N, int &ID){

	/*Calcula espectros de detectores de impacto, grava e carrega arquivos de despejo,
    acumula arquivos de despejo de diferentes execuções e grava os resultados.*/

	double FSAFE=1.000000001e0;
	static const int NIDM=25;
	static const int NBEM=1000;

	int IE, IT;

	//Espectro de energia das partículas que entram.

	if (CIMDET_.LLE[ID-1] == 1){
		IE=1.0e0+(log(*TRACK_mod_.E)-CIMDET_.EL[ID-1])*CIMDET_.RBSE[ID-1];
	}else{
		IE=1.0e0+(*TRACK_mod_.E-CIMDET_.EL[ID-1])*CIMDET_.RBSE[ID-1];
	}

	if (N != CIMDET_.LEDEP[ID-1]){
		CIMDET_.EDEP[ID-1]=CIMDET_.EDEP[ID-1]+CIMDET_.EDEPP[ID-1];
        CIMDET_.EDEP2[ID-1]=CIMDET_.EDEP2[ID-1]+pow(CIMDET_.EDEPP[ID-1],2);
        CIMDET_.EDEPP[ID-1]= *TRACK_mod_.E * *TRACK_mod_.WGHT;
        CIMDET_.LEDEP[ID-1]=N;
	}else{
		CIMDET_.EDEPP[ID-1]=CIMDET_.EDEPP[ID-1]+*TRACK_mod_.E * *TRACK_mod_.WGHT;
	}

	if ((IE > 0) && (IE <= CIMDET_.NE[ID-1])){
		if (N != CIMDET_.LDIT[IE-1][ID-1]){
			CIMDET_.DIT[IE-1][ID-1]=CIMDET_.DIT[IE-1][ID-1]+CIMDET_.DITP[IE-1][ID-1];
            CIMDET_.DIT2[IE-1][ID-1]=CIMDET_.DIT2[IE-1][ID-1]+pow(CIMDET_.DITP[IE-1][ID-1],2);
            CIMDET_.DITP[IE-1][ID-1]=*TRACK_mod_.WGHT;
            CIMDET_.LDIT[IE-1][ID-1]=N;
		}else{
			CIMDET_.DITP[IE-1][ID-1]=CIMDET_.DITP[IE-1][ID-1]+*TRACK_mod_.WGHT;
		}

		if (N != CIMDET_.LDIP[*TRACK_mod_.KPAR-1][IE-1][ID-1]){
			CIMDET_.DIP[*TRACK_mod_.KPAR-1][IE-1][ID-1]=CIMDET_.DIP[*TRACK_mod_.KPAR-1][IE-1][ID-1]+CIMDET_.DIPP[*TRACK_mod_.KPAR-1][IE-1][ID-1];
            CIMDET_.DIP2[*TRACK_mod_.KPAR-1][IE-1][ID-1]=CIMDET_.DIP2[*TRACK_mod_.KPAR-1][IE-1][ID-1]+pow(CIMDET_.DIPP[*TRACK_mod_.KPAR-1][IE-1][ID-1],2);
            CIMDET_.DIPP[*TRACK_mod_.KPAR-1][IE-1][ID-1]=*TRACK_mod_.WGHT;
            CIMDET_.LDIP[*TRACK_mod_.KPAR-1][IE-1][ID-1]=N;
		}else{
			CIMDET_.DIPP[*TRACK_mod_.KPAR-1][IE-1][ID-1]=CIMDET_.DIPP[*TRACK_mod_.KPAR-1][IE-1][ID-1]+*TRACK_mod_.WGHT;
		}

		//Distribuição de Idades das particulas
		if (CIMDET_.NAGE[ID-1] > 0){
			if (CIMDET_.LLAGE[ID-1] == 1){
				IT=1.0e0+(log(*TRACK_mod_.PAGE)-CIMDET_.AGEL[ID-1])*CIMDET_.RBAGE[ID-1];
			}else{
				IT=1.0e0+(*TRACK_mod_.PAGE-CIMDET_.AGEL[ID-1])*CIMDET_.RBAGE[ID-1];
			}
			if ((IT > 0) && (IT <= CIMDET_.NAGE[ID-1])){
				if (N != CIMDET_.LAGEA[IT-1][ID-1]) {
					CIMDET_.AGE[IT-1][ID-1]=CIMDET_.AGE[IT-1][ID-1]+CIMDET_.AGEP[IT-1][ID-1];
                	CIMDET_.AGE2[IT-1][ID-1]=CIMDET_.AGE2[IT-1][ID-1]+pow(CIMDET_.AGEP[IT-1][ID-1],2);
                	CIMDET_.AGEP[IT-1][ID-1]=*TRACK_mod_.WGHT;
                	CIMDET_.LAGEA[IT-1][ID-1]=N;
				}else{
					CIMDET_.AGEP[IT-1][ID-1]=CIMDET_.AGEP[IT-1][ID-1]+*TRACK_mod_.WGHT;
				}
			}
		}
	}


}

void panar2_(double &ECUT){

	/*Simulação da aniquilação de pósitrons em repouso. quanta de aniquilação
    são armazenados na pilha secundária somente quando ECUT é menor que REV.*/

	double REV=5.10998928e5;
	double TREV=2.0e0*REV;
	double PI=3.1415926535897932e0;
	double TWOPI=PI+PI;

	double US, VS, WS, CDT1, DF;

	int ILBA[5];

	if (REV < ECUT)
		return;

	US=*TRACK_mod_.U;
    VS=*TRACK_mod_.V;
    WS=*TRACK_mod_.W;
    CDT1=-1.0e0+2.0e0*rand2_(1.0e0);
    DF=TWOPI*rand2_(2.0e0);

	direct2_(CDT1,DF,US,VS,WS);

	ILBA[1-1]=TRACK_mod_.ILB[1-1]+1;
    ILBA[2-1]=3;
    ILBA[3-1]=6;
    ILBA[4-1]=0;
    ILBA[5-1]=TRACK_mod_.ILB[5-1];
	int WKPARP = 2;
    stores2_(REV,*TRACK_mod_.X,*TRACK_mod_.Y,*TRACK_mod_.Z,US,VS,WS,*TRACK_mod_.WGHT,WKPARP,ILBA,wIPOLI);
	US = -US; VS = -VS; WS = -WS;
    stores2_(REV,*TRACK_mod_.X,*TRACK_mod_.Y,*TRACK_mod_.Z,US,VS,WS,*TRACK_mod_.WGHT,WKPARP,ILBA,wIPOLI);

/*	if (imprimiu==0){
		printf("\n\npanar2\n\n");
		imprimiu++;
	}*/


}

void direct2_(double &CDT, double &DF, double &U, double &V, double &W){

	/*
	Esta sub-rotina calcula os novos cossenos de direção da partícula
	Velocidade após uma colisão com determinado espalhamento polar e azimutal
	ângulos.

	Entrada: U,V,W ... cossenos da direção inicial.
	CDT ..... cosseno do ângulo de dispersão polar.
	DF ...... ângulo de espalhamento azimutal (rad).

	Saída: U,V,W ... novos cossenos de direção.
	CDT e DF permanecem inalterados
		
	*/



	double UV, UVW, FNORM, SDT, SDTSDF, SDTCDF, SUV, UN, VN;

	//Garante a Normalizacao

	UV=U*U+V*V;
    UVW=UV+W*W;

	if (fabs(UVW-1.0e0) > 1.0e-13){
		FNORM=1.0e0/sqrt(UVW);
        U=FNORM*U;
        V=FNORM*V;
        W=FNORM*W;
        UV=U*U+V*V;
	}

	//Calculando a nova direção

	if ((1.0e0-fabs(CDT)) > 1.0e-8){
		SDT=sqrt(1.0e0-CDT*CDT);
	}else{
		 SDT=sqrt(2.0e0*(1.0e0-fabs(CDT)));
	}

	if (SDT < 1.0e-13){
		if (CDT < 0.0e0){
			U=-U;
          	V=-V;
          	W=-W;
		}
	}else{
		SDTSDF=SDT*sin(DF);
        SDTCDF=SDT*cos(DF);
		if (UV > 1.0e-26){
			SUV=sqrt(UV);
            UN=U/SUV;
            VN=V/SUV;
            U=U*CDT+(UN*W*SDTCDF-VN*SDTSDF);
            V=V*CDT+(VN*W*SDTCDF+UN*SDTSDF);
            W=W*CDT-SUV*SDTCDF;
		}else{
			if (W > 0.0e0){
				U=SDTCDF;
           		V=SDTSDF;
           		W=CDT;
			}else{
				U=-SDTCDF;
            	V=-SDTSDF;
           	 	W=-CDT;
			}
		}
	}
/*	if (imprimiu==0){
		printf("\n\nDIRECT2\n\n");
		imprimiu++;
	}*/

}

void stores2_(double &EI, double &XI, double &YI, double &ZI, double &UI, double &VI, double &WI, double &WGHTI, int &KPARI, int *ILBI, int &IPOLI){

/*
	Esta sub-rotina armazena o estado inicial de uma nova partícula secundária
	na pilha secundária. Os valores de entrada são:
	EI ........... energia inicial.
	XI, YI, ZI ... coordenadas da posição inicial.
	UI, VI, WI ... cossenos de direção inicial.
	WGHTI ........ peso (=1 na simulação analógica).
	KPARI ........ tipo de partícula (1: elétron, 2: fóton,
	3: pósitron).
	ILBI(5) ...... rótulos de partículas.
	IPOLI ........ sinalizador de polarização.

	O parâmetro NMS fixa o tamanho da pilha secundária (ou seja, 
	número máximo de partículas que podem ser armazenadas). Se este número for
	excedido, uma mensagem de aviso é impressa na unidade 26. Quando a memória
	O armazenamento de é esgotado, cada nova partícula secundária é armazenada no
	posição do elétron secundário menos energético ou fóton já
	produzido, que é assim descartado.
*/

	/*if (imprimiu==0){
		printf("\n\nSTORES2\n\n");
		imprimiu++;
	}*/
	double EME, EMG;
	int IE, IG, IS;

	if (*SECST_.NSEC < NMS){
		*SECST_.NSEC=*SECST_.NSEC+1;
        IS=*SECST_.NSEC;	
	}else{
		if (*CERSEC_.IERSEC == 0){
			printf("   *** WARNING: (STORES) not enough storage for secondaries.\n EABS(KPAR,MAT) or the parameter NMS should be enlarged\n");
			*CERSEC_.IERSEC = 1;
		}
		*SECST_.NSEC=NMS;
        EME=1.0e35;
        EMG=1.0e35;
        IE=0;
        IG=0;
		for (int I =1; I <= NMS; I ++){
			if (SECST_.KS[I-1] == 1){
				if (SECST_.ES[I-1] < EME){
					EME=SECST_.ES[I-1] ;
            		IE=I;
				}
			}else if (SECST_.KS[I-1] == 2){
				if (SECST_.ES[I-1] < EMG){
					EMG=SECST_.ES[I-1];
            		IG=I;
				}
			}
		}

		if (IE > 0){
			IS=IE;
		}else if (IG > 0){
			IS=IG;
		}else{
			printf("   *** Not enough storage for secondary positrons./n       JOB ABORTED.");
			IS=0;
			exit(0);
		}
	}

	SECST_.ES[IS-1]=EI;
    SECST_.XS[IS-1]=XI;
    SECST_.YS[IS-1]=YI;
    SECST_.ZS[IS-1]=ZI;
    SECST_.US[IS-1]=UI;
    SECST_.VS[IS-1]=VI;
    SECST_.WS[IS-1]=WI;
    SECST_.WGHTS[IS-1]=WGHTI;
    SECST_.KS[IS-1]=KPARI;
    SECST_.IBODYS[IS-1]=*TRACK_mod_.IBODY;
    SECST_.MS[IS-1]=*TRACK_mod_.MAT;
    SECST_.ILBS[IS-1][1-1]=ILBI[1-1];
    SECST_.ILBS[IS-1][2-1]=ILBI[2-1];
    SECST_.ILBS[IS-1][3-1]=ILBI[3-1];
    SECST_.ILBS[IS-1][4-1]=ILBI[4-1];
    SECST_.ILBS[IS-1][5-1]=ILBI[5-1];

	if (IPOLI == 1){
		SECST_.SP1S[IS-1]=*TRACK_mod_.SP1;
        SECST_.SP2S[IS-1]=*TRACK_mod_.SP2;
        SECST_.SP3S[IS-1]=*TRACK_mod_.SP3;
        SECST_.IPOLS[IS-1]=IPOLI;
	}else{
		SECST_.SP1S[IS-1]=0.e0;
        SECST_.SP2S[IS-1]=0.e0;
        SECST_.SP3S[IS-1]=0.e0;
        SECST_.IPOLS[IS-1]=0;
	}

	SECST_.PAGES[IS-1]=*TRACK_mod_.PAGE;


}

void sdose2_(double &DEP, double &XD, double &YD, double &ZD, int &MATC, int &N){

/*
	Registra a distribuição de dose dentro da caixa de dose, grava e carrega
	despeja arquivos, acumula arquivos de despejo de diferentes execuções e grava
	Resultados .
*/
	
	double PI=3.1415926535897932e0;
	double FSAFE=1.000000001e0;
	double RNCS=1.0e0/5.0e0;
	int NCS = 5;
	static const int NDXM=201;
	static const int NDYM=201;
	static const int NDZM=201;

	int I1, I2, I3;
	double RD;

	//Distribuição de Dose
	//printf("\nKDOSE: %d\n", *CDOSE1_.KDOSE);
	if (*CDOSE1_.KDOSE == 1){ //Caixa
		if ((ZD > CDOSE3_.DXL[3-1]) && (ZD < CDOSE3_.DXU[3-1])){
			I3=1.0e0+(ZD-CDOSE3_.DXL[3-1])*CDOSE3_.RBDOSE[3-1];
			if (N != CDOSE2_.LDDOSE[I3-1]){
				CDOSE2_.DDOSE[I3-1]=CDOSE2_.DDOSE[I3-1]+CDOSE2_.DDOSEP[I3-1];
            	CDOSE2_.DDOSE2[I3-1]=CDOSE2_.DDOSE2[I3-1]+pow(CDOSE2_.DDOSEP[I3-1],2);
            	CDOSE2_.DDOSEP[I3-1]=DEP* PENELOPE_mod_.RDEN[MATC-1];
            	CDOSE2_.LDDOSE[I3-1]=N;
			}else{
				CDOSE2_.DDOSEP[I3-1]=CDOSE2_.DDOSEP[I3-1]+DEP*PENELOPE_mod_.RDEN[MATC-1];
			}

			if (((XD > CDOSE3_.DXL[1-1])  && (XD< CDOSE3_.DXU[1-1])) && ((YD > CDOSE3_.DXL[2-1]) && (YD < CDOSE3_.DXU[2-1]))){
				I1=1.0e0+(XD-CDOSE3_.DXL[1-1])*CDOSE3_.RBDOSE[1-1];
           		I2=1.0e0+(YD-CDOSE3_.DXL[2-1])*CDOSE3_.RBDOSE[2-1];
				if (N != CDOSE1_.LDOSE[I3-1][I2-1][I1-1]){
					CDOSE1_.DOSE[I3-1][I2-1][I1-1]=CDOSE1_.DOSE[I3-1][I2-1][I1-1]+CDOSE1_.DOSEP[I3-1][I2-1][I1-1];
              		CDOSE1_.DOSE2[I3-1][I2-1][I1-1]=CDOSE1_.DOSE2[I3-1][I2-1][I1-1]+pow(CDOSE1_.DOSEP[I3-1][I2-1][I1-1],2);
              		CDOSE1_.DOSEP[I3-1][I2-1][I1-1]=DEP;
              		CDOSE1_.LDOSE[I3-1][I2-1][I1-1]=N;
				}else{
					CDOSE1_.DOSEP[I3-1][I2-1][I1-1]=CDOSE1_.DOSEP[I3-1][I2-1][I1-1]+DEP;
				}
			}
		}
	}else if (*CDOSE1_.KDOSE == 2){ //Cilindro
		if ((ZD > CDOSE3_.DXL[3-1]) && (ZD < CDOSE3_.DXU[3-1])){
			I3=1.0e0+(ZD-CDOSE3_.DXL[3-1])*CDOSE3_.RBDOSE[3-1];
			if (N != CDOSE2_.LDDOSE[I3-1]){
				CDOSE2_.DDOSE[I3-1]=CDOSE2_.DDOSE[I3-1]+CDOSE2_.DDOSEP[I3-1];
				CDOSE2_.DDOSE2[I3-1]=CDOSE2_.DDOSE2[I3-1]+pow(CDOSE2_.DDOSEP[I3-1],2);
				CDOSE2_.DDOSEP[I3-1]=DEP* PENELOPE_mod_.RDEN[MATC-1];
				CDOSE2_.LDDOSE[I3-1]=N;
			}else{
				CDOSE2_.DDOSEP[I3-1]=CDOSE2_.DDOSEP[I3-1]+DEP*PENELOPE_mod_.RDEN[MATC-1];
			}

			RD=sqrt(XD*XD+YD*YD);
			if (RD < CDOSE3_.DXU[1-1]){
				I1=1.0e0+RD*CDOSE3_.RBDOSE[1-1];
				I2=1;
				if (N != CDOSE1_.LDOSE[I3-1][I2-1][I1-1]){
					CDOSE1_.DOSE[I3-1][I2-1][I1-1]=CDOSE1_.DOSE[I3-1][I2-1][I1-1]+CDOSE1_.DOSEP[I3-1][I2-1][I1-1];
					CDOSE1_.DOSE2[I3-1][I2-1][I1-1]=CDOSE1_.DOSE2[I3-1][I2-1][I1-1]+pow(CDOSE1_.DOSEP[I3-1][I2-1][I1-1],2);
					CDOSE1_.DOSEP[I3-1][I2-1][I1-1]=DEP;
					CDOSE1_.LDOSE[I3-1][I2-1][I1-1]=N;
				}else{
					CDOSE1_.DOSEP[I3-1][I2-1][I1-1]=CDOSE1_.DOSEP[I3-1][I2-1][I1-1]+DEP;
				}
			}
		}
	}else{ //Esfera
		RD=sqrt(XD*XD+YD*YD+ZD*ZD);
		if (RD < CDOSE3_.DXU[1-1]){
			I1=1.0e0+RD*CDOSE3_.RBDOSE[1-1];
          	I2=1;
          	I3=1;
			if (N != CDOSE1_.LDOSE[I3-1][I2-1][I1-1]){
				CDOSE1_.DOSE[I3-1][I2-1][I1-1]=CDOSE1_.DOSE[I3-1][I2-1][I1-1]+CDOSE1_.DOSEP[I3-1][I2-1][I1-1];
				CDOSE1_.DOSE2[I3-1][I2-1][I1-1]=CDOSE1_.DOSE2[I3-1][I2-1][I1-1]+pow(CDOSE1_.DOSEP[I3-1][I2-1][I1-1],2);
				CDOSE1_.DOSEP[I3-1][I2-1][I1-1]=DEP;
				CDOSE1_.LDOSE[I3-1][I2-1][I1-1]=N;
			}else{
				CDOSE1_.DOSEP[I3-1][I2-1][I1-1]=CDOSE1_.DOSEP[I3-1][I2-1][I1-1]+DEP;
			}
		}
	}

/*		if (imprimiu==0){
		printf("\n\nsdose2\n\n");
		imprimiu++;
	}*/


}

void start2_(){

/*
	Esta sub-rotina força o próximo evento a ser um evento soft artificial.
	Deve ser chamado quando uma nova trilha de partículas (primária ou secundária) é
	é iniciado e quando cruza uma interface.
*/

	if ((*TRACK_mod_.E < *CEGRID_.EMIN) || (*TRACK_mod_.E > 0.99999999e0 * *CEGRID_.EU)){
		printf("   *** Energy out of range. KPAR = %2d, E = %.5E eV\n", *TRACK_mod_.KPAR,*TRACK_mod_.E);
		for (int J = 1; J <= 5; J++){
			printf("       ILB%d = %d\n", J, TRACK_mod_.ILB[J-1]);
		}
		printf("       EMIN = %.5E eV, EMAX = %.5E eV\n", *CEGRID_.EL,*CEGRID_.EU);
		printf("       Check the values of EABS(KPAR,M) and EMAX.\n");
		exit(0);
	}
	*CJUMP1_.MHINGE=0;
    *CJUMP1_.ELAST1=*TRACK_mod_.E+1.0e30;
    *CJUMP1_.ELAST2=*CJUMP1_.ELAST1;

	/*	if (imprimiu==0){
		printf("\n\nSTART2\n\n");
		imprimiu++;
	}*/

}

void jump2_(double &DSMAX, double &DS){
	/*
	Cálculo do caminho livre do ponto de partida até a posição
	do próximo evento e das probabilidades de ocorrência de diferentes
	eventos .

	Argumentos :
	DSMAX .... comprimento máximo permitido do passo (entrada),
	DS ....... comprimento do segmento (saída).

	Saída , através do módulo PENELOPE_mod:
	E0STEP ... energia no início do segmento,
	DESOFT ... perda de energia devido a interações suaves ao longo da etapa,
	SSOFT .... poder de parada devido a interações suaves,
	= DESOFT/passo_comprimento.
	*/

	/*	if (imprimiu==0){
		printf("\n\nJUMP2\n\n");
		imprimiu++;
	}*/

    double DSMAXP, DSMC, EDE0, VDE0, FSEDE, FSVDE, EDEM, VDEM, W21, ELOWER, XE1, XEK1, STLWR, EDE, VDE, SIGMA;
	double RU, EDE2, VDE3, PNULL;
	int KE1;

	if (*TRACK_mod_.KPAR == 1){ //eletrons
		if (*CJUMP1_.MHINGE == 1){
			if (*TRACK_mod_.E < *CJUMP1_.ELAST1){
				*CEGRID_.XEL=log(*TRACK_mod_.E);
            	*CEGRID_.XE=1.0e0+(*CEGRID_.XEL-*CEGRID_.DLEMP1)* *CEGRID_.DLFC;
            	*CEGRID_.KE=*CEGRID_.XE;
            	*CEGRID_.XEK=*CEGRID_.XE-*CEGRID_.KE;
				int wvar = 1;
            	eimfp2_(wvar);
            	*CJUMP1_.ELAST1=*TRACK_mod_.E;
			}
			DS=*CJUMP0_.DSR;
			return;
		}
		*PENELOPE_mod_.E0STEP=*TRACK_mod_.E;
		if (*TRACK_mod_.E < *CJUMP1_.ELAST2){
			*CEGRID_.XEL=log(*TRACK_mod_.E);
            *CEGRID_.XE=1.0e0+(*CEGRID_.XEL-*CEGRID_.DLEMP1)* *CEGRID_.DLFC;
            *CEGRID_.KE=*CEGRID_.XE;
            *CEGRID_.XEK=*CEGRID_.XE-*CEGRID_.KE;
			int wvar = 2;
            eimfp2_(wvar);
            *CJUMP1_.ELAST2=*TRACK_mod_.E;
            *CJUMP1_.ELAST1=*TRACK_mod_.E;
		}

		//Caminho livre médio rígido inverso (probabilidade de interação por unidade comprimento do caminho).
		*CJUMP0_.ST=CJUMP0_.P[2-1]+CJUMP0_.P[3-1]+CJUMP0_.P[4-1]+CJUMP0_.P[5-1]+CJUMP0_.P[8-1];
        DSMAXP=DSMAX;

		/*
		Interações de parada suave.
 		KSOFTI=1, parada suave está ativa,
 		KSOFTI=0, a parada suave não está ativa.
		*/

		if (*CJUMP0_.W1 > 1.0e-20){
			*CJUMP1_.KSOFTI=1;
			/*
			O comprimento máximo do passo, DSMAXP, é determinado em termos do
			valor DSMAX de entrada (que é especificado pelo usuário) e a média
			caminho livre para interações difíceis (1/ST).
			*/
			DSMC=4.0e0 / *CJUMP0_.ST;
			if (DSMAXP > DSMC){
				DSMAXP=DSMC;
			}else if (DSMAXP < 1.0e-8){
				DSMAXP=DSMC;
			}

			//O valor de DSMAXP é randomizado para eliminar artefatos de dose no final da primeira etapa.
			DSMAXP=(0.5e0+rand2_(1.0e0)*0.5e0)*DSMAXP;

			//Limite superior para a probabilidade de interação ao longo da etapa (incluindo straggling de energia suave).

			EDE0=*CJUMP0_.W1*DSMAXP;
            VDE0=*CJUMP0_.W2*DSMAXP;
            FSEDE=fmax(1.0e0-CEIMFP_.DW1EL[*CEGRID_.KE-1][*TRACK_mod_.MAT-1]*EDE0,0.75e0);
            FSVDE=fmax(1.0e0-CEIMFP_.DW2EL[*CEGRID_.KE-1][*TRACK_mod_.MAT-1]*EDE0,0.75e0);
            EDEM=EDE0*FSEDE;
            VDEM=VDE0*FSVDE;
            W21=VDEM/EDEM;

			if (EDEM > 9.0e0*W21){
				ELOWER=fmax(*TRACK_mod_.E-(EDEM+3.0e0*sqrt(VDEM)),*CEGRID_.EMIN);
			}else if (EDEM > 3.0e0*W21){
				ELOWER=fmax(*TRACK_mod_.E-(EDEM+sqrt(3.0e0*VDEM)),*CEGRID_.EMIN);
			}else{
				ELOWER=fmax(*TRACK_mod_.E-1.5e0*(EDEM+W21),*CEGRID_.EMIN);
			}

			XE1=1.0e0+(log(ELOWER)-*CEGRID_.DLEMP1)* *CEGRID_.DLFC;
            KE1=XE1;
            XEK1=XE1-KE1;
            STLWR=exp(CEIMFP_.SETOT[KE1-1][*TRACK_mod_.MAT-1]+(CEIMFP_.SETOT[KE1+1-1][*TRACK_mod_.MAT-1]-CEIMFP_.SETOT[KE1-1][*TRACK_mod_.MAT-1])*XEK1);
            *CJUMP0_.ST=fmax(*CJUMP0_.ST,STLWR);
		}else{
			*CJUMP1_.KSOFTI=0;
            *PENELOPE_mod_.DESOFT=0.0e0;
            *PENELOPE_mod_.SSOFT=0.0e0;
		}

		/*
		Dispersão elástica suave.
		KSOFTE=1, dispersão suave está ativa,
		KSOFTE=0, dispersão suave não está ativa.
		*/

		if (*CJUMP0_.T1 > 1.0e-20){
			*CJUMP1_.KSOFTE=1;
		}else{
			*CJUMP1_.KSOFTE=0;
		}

		/*
		Interações delta.
		KDELTA=0, segue-se uma interação difícil,
		KDELTA=1, segue-se uma interação delta.
		*/

		*CJUMP0_.DST=-log(rand2_(2.0e0))/ *CJUMP0_.ST;
		if (*CJUMP0_.DST< DSMAXP){
			*CJUMP1_.KDELTA=0;
		}else{
			*CJUMP0_.DST=DSMAXP;
            *CJUMP1_.KDELTA=1;
		}

		if (*CJUMP1_.KSOFTE+*CJUMP1_.KSOFTI == 0){
			*CJUMP1_.MHINGE=1;
         	DS=*CJUMP0_.DST;
		}else{
			DS=*CJUMP0_.DST*rand2_(3.0e0);
         	*CJUMP0_.DSR=*CJUMP0_.DST-DS;
			 if (*CJUMP1_.KSOFTI == 1){
				 if (*CJUMP0_.DST < 1.0e-8){
					*PENELOPE_mod_.SSOFT=*CJUMP0_.W1;
              		*PENELOPE_mod_.DESOFT=*PENELOPE_mod_.SSOFT * *CJUMP0_.DST;
				 }else{
					EDE0=*CJUMP0_.W1* *CJUMP0_.DST;
              		VDE0=*CJUMP0_.W2* *CJUMP0_.DST;
              		FSEDE=fmax(1.0e0-CEIMFP_.DW1EL[*CEGRID_.KE-1][*TRACK_mod_.MAT-1]*EDE0,0.75e0);
             		FSVDE=fmax(1.0e0-CEIMFP_.DW2EL[*CEGRID_.KE-1][*TRACK_mod_.MAT-1]*EDE0,0.75e0);
             		EDE=EDE0*FSEDE;
              		VDE=VDE0*FSVDE;

					//Geração de valores aleatórios DE com média EDE e variância VDE.
					SIGMA=sqrt(VDE);
					if (SIGMA < 0.333333333e0*EDE){
						//Distribuição gaussiana truncada.
						*PENELOPE_mod_.DESOFT=EDE+rndg32_()*SIGMA;
					}else{
						RU=rand2_(4.0e0);
                		EDE2=EDE*EDE;
                		VDE3=3.0e0*VDE;
						if (EDE2 < VDE3){
							PNULL=(VDE3-EDE2)/(VDE3+3.0e0*EDE2);
							if (RU < PNULL){
								*PENELOPE_mod_.DESOFT=0.0e0;
                    			*PENELOPE_mod_.SSOFT=0.0e0;
								if (*CJUMP1_.KSOFTE == 0){
                      				*CJUMP1_.MHINGE=1;
                      				DS=*CJUMP0_.DST;
								}else{
									*CJUMP1_.KSOFTI=0;
								}
								return;
							}else{
								//Distribuição Uniforme
								*PENELOPE_mod_.DESOFT=1.5e0*(EDE+VDE/EDE)*(RU-PNULL)/(1.0e0-PNULL);
							}
						}else{
							*PENELOPE_mod_.DESOFT=EDE+(2.0e0*RU-1.0e0)*sqrt(VDE3);
						}
					}
					*PENELOPE_mod_.SSOFT=*PENELOPE_mod_.DESOFT/ *CJUMP0_.DST;
				 }
			 }
		}
		return;
	}else if (*TRACK_mod_.KPAR == 3){ //Positrons

		if (*CJUMP1_.MHINGE == 1){
			if (*TRACK_mod_.E < *CJUMP1_.ELAST1){
				*CEGRID_.XEL=log(*TRACK_mod_.E);
            	*CEGRID_.XE=1.0e0+(*CEGRID_.XEL-*CEGRID_.DLEMP1)* *CEGRID_.DLFC;
            	*CEGRID_.KE=*CEGRID_.XE;
            	*CEGRID_.XEK=*CEGRID_.XE-*CEGRID_.KE;
				int wvar = 1;
            	pimfp2_(wvar);
            	*CJUMP1_.ELAST1=*TRACK_mod_.E;
			}
			DS=*CJUMP0_.DSR;
            return;
		}

		*PENELOPE_mod_.E0STEP=*TRACK_mod_.E;
		if (*TRACK_mod_.E < *CJUMP1_.ELAST2){
			*CEGRID_.XEL=log(*TRACK_mod_.E);
          	*CEGRID_.XE=1.0e0+(*CEGRID_.XEL-*CEGRID_.DLEMP1)* *CEGRID_.DLFC;
          	*CEGRID_.KE=*CEGRID_.XE;
          	*CEGRID_.XEK=*CEGRID_.XE-*CEGRID_.KE;
         	int wvar = 2;
            pimfp2_(wvar);
         	*CJUMP1_.ELAST2=*TRACK_mod_.E;
          	*CJUMP1_.ELAST1=*TRACK_mod_.E;
		}

		//Caminho livre médio rígido inverso (probabilidade de interação por unidade comprimento do caminho).
		*CJUMP0_.ST=CJUMP0_.P[2-1]+CJUMP0_.P[3-1]+CJUMP0_.P[4-1]+CJUMP0_.P[5-1]+CJUMP0_.P[6-1]+CJUMP0_.P[8-1];
        DSMAXP=DSMAX;

		/*
		Interações de parada suave.
 		KSOFTI=1, parada suave está ativa,
 		KSOFTI=0, a parada suave não está ativa.
		*/

		if (*CJUMP0_.W1 > 1.0e-20){
			*CJUMP1_.KSOFTI=1;
			/*
			O comprimento máximo do passo, DSMAXP, é determinado em termos do
			valor DSMAX de entrada (que é especificado pelo usuário) e a média
			caminho livre para interações difíceis (1/ST).
			*/
			DSMC=4.0e0 / *CJUMP0_.ST;
			if (DSMAXP > DSMC){
				DSMAXP=DSMC;
			}else if (DSMAXP < 1.0e-8){
				DSMAXP=DSMC;
			}

			//O valor de DSMAXP é randomizado para eliminar artefatos de dose no final da primeira etapa.
			DSMAXP=(0.5e0+rand2_(1.0e0)*0.5e0)*DSMAXP;

			//Limite superior para a probabilidade de interação ao longo da etapa (incluindo straggling de energia suave).

			EDE0=*CJUMP0_.W1*DSMAXP;
            VDE0=*CJUMP0_.W2*DSMAXP;
            FSEDE=fmax(1.0e0-CPIMFP_.DW1PL[*CEGRID_.KE-1][*TRACK_mod_.MAT-1]*EDE0,0.75e0);
            FSVDE=fmax(1.0e0-CPIMFP_.DW2PL[*CEGRID_.KE-1][*TRACK_mod_.MAT-1]*EDE0,0.75e0);
            EDEM=EDE0*FSEDE;
            VDEM=VDE0*FSVDE;
            W21=VDEM/EDEM;

			if (EDEM > 9.0e0*W21){
				ELOWER=fmax(*TRACK_mod_.E-(EDEM+3.0e0*sqrt(VDEM)),*CEGRID_.EMIN);
			}else if (EDEM > 3.0e0*W21){
				ELOWER=fmax(*TRACK_mod_.E-(EDEM+sqrt(3.0e0*VDEM)),*CEGRID_.EMIN);
			}else{
				ELOWER=fmax(*TRACK_mod_.E-1.5e0*(EDEM+W21),*CEGRID_.EMIN);
			}

			XE1=1.0e0+(log(ELOWER)-*CEGRID_.DLEMP1)* *CEGRID_.DLFC;
            KE1=XE1;
            XEK1=XE1-KE1;
            STLWR=exp(CPIMFP_.SPTOT[KE1-1][*TRACK_mod_.MAT-1]+(CPIMFP_.SPTOT[KE1+1-1][*TRACK_mod_.MAT-1]-CPIMFP_.SPTOT[KE1-1][*TRACK_mod_.MAT-1])*XEK1);
            *CJUMP0_.ST=fmax(*CJUMP0_.ST,STLWR);
		}else{
			*CJUMP1_.KSOFTI=0;
            *PENELOPE_mod_.DESOFT=0.0e0;
            *PENELOPE_mod_.SSOFT=0.0e0;
		}

		/*
		Dispersão elástica suave.
		KSOFTE=1, dispersão suave está ativa,
		KSOFTE=0, dispersão suave não está ativa.
		*/

		if (*CJUMP0_.T1 > 1.0e-20){
			*CJUMP1_.KSOFTE=1;
		}else{
			*CJUMP1_.KSOFTE=0;
		}

		/*
		Interações delta.
		KDELTA=0, segue-se uma interação difícil,
		KDELTA=1, segue-se uma interação delta.
		*/

		*CJUMP0_.DST=-log(rand2_(2.0e0))/ *CJUMP0_.ST;
		if (*CJUMP0_.DST< DSMAXP){
			*CJUMP1_.KDELTA=0;
		}else{
			*CJUMP0_.DST=DSMAXP;
            *CJUMP1_.KDELTA=1;
		}

		if (*CJUMP1_.KSOFTE+*CJUMP1_.KSOFTI == 0){
			*CJUMP1_.MHINGE=1;
         	DS=*CJUMP0_.DST;
		}else{
			DS=*CJUMP0_.DST*rand2_(3.0e0);
         	*CJUMP0_.DSR=*CJUMP0_.DST-DS;
			 if (*CJUMP1_.KSOFTI == 1){
				 if (*CJUMP0_.DST < 1.0e-8){
					*PENELOPE_mod_.SSOFT=*CJUMP0_.W1;
              		*PENELOPE_mod_.DESOFT=*PENELOPE_mod_.SSOFT * *CJUMP0_.DST;
				 }else{
					EDE0=*CJUMP0_.W1* *CJUMP0_.DST;
              		VDE0=*CJUMP0_.W2* *CJUMP0_.DST;
              		FSEDE=fmax(1.0e0-CPIMFP_.DW1PL[*CEGRID_.KE-1][*TRACK_mod_.MAT-1]*EDE0,0.75e0);
             		FSVDE=fmax(1.0e0-CPIMFP_.DW2PL[*CEGRID_.KE-1][*TRACK_mod_.MAT-1]*EDE0,0.75e0);
             		EDE=EDE0*FSEDE;
              		VDE=VDE0*FSVDE;

					//Geração de valores aleatórios DE com média EDE e variância VDE.
					SIGMA=sqrt(VDE);
					if (SIGMA < 0.333333333e0*EDE){
						//Distribuição gaussiana truncada.
						*PENELOPE_mod_.DESOFT=EDE+rndg32_()*SIGMA;
					}else{
						RU=rand2_(4.0e0);
                		EDE2=EDE*EDE;
                		VDE3=3.0e0*VDE;
						if (EDE2 < VDE3){
							PNULL=(VDE3-EDE2)/(VDE3+3.0e0*EDE2);
							if (RU < PNULL){
								*PENELOPE_mod_.DESOFT=0.0e0;
                    			*PENELOPE_mod_.SSOFT=0.0e0;
								if (*CJUMP1_.KSOFTE == 0){
                      				*CJUMP1_.MHINGE=1;
                      				DS=*CJUMP0_.DST;
								}else{
									*CJUMP1_.KSOFTI=0;
								}
								return;
							}else{
								//Distribuição Uniforme
								*PENELOPE_mod_.DESOFT=1.5e0*(EDE+VDE/EDE)*(RU-PNULL)/(1.0e0-PNULL);
							}
						}else{
							*PENELOPE_mod_.DESOFT=EDE+(2.0e0*RU-1.0e0)*sqrt(VDE3);
						}
					}
					*PENELOPE_mod_.SSOFT=*PENELOPE_mod_.DESOFT/ *CJUMP0_.DST;
				 }
			 }
		}
		return;
	}else{ //Fotons

		if (*TRACK_mod_.E < *CJUMP1_.ELAST1){
			*CEGRID_.XEL=log(*TRACK_mod_.E);
			*CEGRID_.XE=1.0e0+(*CEGRID_.XEL-*CEGRID_.DLEMP1)* *CEGRID_.DLFC;
			*CEGRID_.KE=*CEGRID_.XE;
			*CEGRID_.XEK=*CEGRID_.XE-*CEGRID_.KE;
			gimfp2_();
			*CJUMP1_.ELAST1=*TRACK_mod_.E;
			*CJUMP0_.ST=CJUMP0_.P[1-1]+CJUMP0_.P[2-1]+CJUMP0_.P[3-1]+CJUMP0_.P[4-1]+CJUMP0_.P[8-1];
		}

		DS=-log(rand2_(1.0e0))/ *CJUMP0_.ST;
	}

/*	if (imprimiu==0){
		printf("\n\nJUMP2\n\n");
		imprimiu++;
	}*/

}

double rndg32_(){

	/*
	Esta função entrega valores aleatórios no intervalo (-3.0,3.0)
 	amostrado de uma distribuição gaussiana truncada que tem média zero e
	variação da unidade. A amostragem é realizada pelo método RITA.
	*/

	/*if (imprimiu==0){
		printf("\n\nRNDG32\n\n");
		imprimiu++;
	}*/

	static const int NR = 128;
	// Selection of the interval (Walker's aliasing).

	double RN, TST, RR, D, resultado;
	int K, I;

	RN=rand2_(1.0e0)* *CRNDG3_.NPM1+1.0e0;
    K=int(RN);
    TST=RN-K;
	if (TST < CRNDG3_.F[K-1]){
		I=K;
        RR=TST;
        D=CRNDG3_.F[K-1];
	}else{
		I=CRNDG3_.KA[K-1];
        RR=TST-CRNDG3_.F[K-1];
        D=1.0e0-CRNDG3_.F[K-1];
	}

	//Amostragem da distribuição cumulativa inversa racional.
	if (RR > 1.0e-12){
		resultado = CRNDG3_.X[I-1]+((1.0e0+CRNDG3_.A[I-1]+CRNDG3_.B[I-1])*D*RR/(D*D+(CRNDG3_.A[I-1]*D+CRNDG3_.B[I-1]*RR)*RR))*(CRNDG3_.X[I+1-1]-CRNDG3_.X[I-1]);
	}else{
		resultado=CRNDG3_.X[I-1]+rand2_(2.0e0)*(CRNDG3_.X[I+1-1]-CRNDG3_.X[I-1]);
	}

	return resultado;

}

void pimfp2_(int &IEND){
	/*
	Esta sub-rotina calcula os caminhos livres médios inversos para
	ações de pósitrons com a energia atual no material M.
	*/

   /*if (imprimiu==0){
		printf("\n\nPIMFP2\n\n");
		imprimiu++;
	}*/

	CJUMP0_.P[2-1]=exp(CPIMFP_.SPHEL[*CEGRID_.KE-1][*TRACK_mod_.MAT-1]+(CPIMFP_.SPHEL[*CEGRID_.KE+1-1][*TRACK_mod_.MAT-1]-CPIMFP_.SPHEL[*CEGRID_.KE-1][*TRACK_mod_.MAT-1])* *CEGRID_.XEK);
    CJUMP0_.P[3-1]=exp(CPIMFP_.SPHIN[*CEGRID_.KE-1][*TRACK_mod_.MAT-1]+(CPIMFP_.SPHIN[*CEGRID_.KE+1-1][*TRACK_mod_.MAT-1]-CPIMFP_.SPHIN[*CEGRID_.KE-1][*TRACK_mod_.MAT-1])* *CEGRID_.XEK);
    CJUMP0_.P[4-1]=exp(CPIMFP_.SPHBR[*CEGRID_.KE-1][*TRACK_mod_.MAT-1]+(CPIMFP_.SPHBR[*CEGRID_.KE+1-1][*TRACK_mod_.MAT-1]-CPIMFP_.SPHBR[*CEGRID_.KE-1][*TRACK_mod_.MAT-1])* *CEGRID_.XEK);
    CJUMP0_.P[5-1]=exp(CPIMFP_.SPISI[*CEGRID_.KE-1][*TRACK_mod_.MAT-1]+(CPIMFP_.SPISI[*CEGRID_.KE+1-1][*TRACK_mod_.MAT-1]-CPIMFP_.SPISI[*CEGRID_.KE-1][*TRACK_mod_.MAT-1])* *CEGRID_.XEK);
    CJUMP0_.P[6-1]=exp(CPIMFP_.SPAN[*CEGRID_.KE-1][*TRACK_mod_.MAT-1]+(CPIMFP_.SPAN[*CEGRID_.KE+1-1][*TRACK_mod_.MAT-1]-CPIMFP_.SPAN[*CEGRID_.KE-1][*TRACK_mod_.MAT-1])* *CEGRID_.XEK);
    CJUMP0_.P[8-1]=0.0e0;

	if (IEND == 1){
		return;
	}

	if (CPIMFP_.W1P[*CEGRID_.KE+1-1][*TRACK_mod_.MAT-1] > -78.3e0){
		*CJUMP0_.W1=exp(CPIMFP_.W1P[*CEGRID_.KE-1][*TRACK_mod_.MAT-1]+(CPIMFP_.W1P[*CEGRID_.KE+1-1][*TRACK_mod_.MAT-1]-CPIMFP_.W1P[*CEGRID_.KE-1][*TRACK_mod_.MAT-1])* *CEGRID_.XEK);
        *CJUMP0_.W2=exp(CPIMFP_.W2P[*CEGRID_.KE-1][*TRACK_mod_.MAT-1]+(CPIMFP_.W2P[*CEGRID_.KE+1-1][*TRACK_mod_.MAT-1]-CPIMFP_.W2P[*CEGRID_.KE-1][*TRACK_mod_.MAT-1])* *CEGRID_.XEK);
	}else{
		*CJUMP0_.W1=0.0e0;
		*CJUMP0_.W2=0.0e0;
	}

	if (CPIMFP_.T1P[*CEGRID_.KE+1-1][*TRACK_mod_.MAT-1] > -78.3e0){
		*CJUMP0_.T1=exp(CPIMFP_.T1P[*CEGRID_.KE-1][*TRACK_mod_.MAT-1]+(CPIMFP_.T1P[*CEGRID_.KE+1-1][*TRACK_mod_.MAT-1]-CPIMFP_.T1P[*CEGRID_.KE-1][*TRACK_mod_.MAT-1])* *CEGRID_.XEK);
        *CJUMP0_.T2=exp(CPIMFP_.T2P[*CEGRID_.KE-1][*TRACK_mod_.MAT-1]+(CPIMFP_.T2P[*CEGRID_.KE+1-1][*TRACK_mod_.MAT-1]-CPIMFP_.T2P[*CEGRID_.KE-1][*TRACK_mod_.MAT-1])* *CEGRID_.XEK);
	}else{
		*CJUMP0_.T1=0.0e0;
		*CJUMP0_.T2=0.0e0;
	}
}

void eimfp2_(int &IEND){
	/*
	Esta sub-rotina calcula os caminhos livres médios inversos para
	ações de eletrons com a energia atual no material M.
	*/
	/*if (imprimiu==0){
		printf("\n\nEIMFP2\n\n");
		imprimiu++;
	}*/

	CJUMP0_.P[2-1]=exp(CEIMFP_.SEHEL[*CEGRID_.KE-1][*TRACK_mod_.MAT-1]+(CEIMFP_.SEHEL[*CEGRID_.KE+1-1][*TRACK_mod_.MAT-1]-CEIMFP_.SEHEL[*CEGRID_.KE-1][*TRACK_mod_.MAT-1])* *CEGRID_.XEK);
    CJUMP0_.P[3-1]=exp(CEIMFP_.SEHIN[*CEGRID_.KE-1][*TRACK_mod_.MAT-1]+(CEIMFP_.SEHIN[*CEGRID_.KE+1-1][*TRACK_mod_.MAT-1]-CEIMFP_.SEHIN[*CEGRID_.KE-1][*TRACK_mod_.MAT-1])* *CEGRID_.XEK);
    CJUMP0_.P[4-1]=exp(CEIMFP_.SEHBR[*CEGRID_.KE-1][*TRACK_mod_.MAT-1]+(CEIMFP_.SEHBR[*CEGRID_.KE+1-1][*TRACK_mod_.MAT-1]-CEIMFP_.SEHBR[*CEGRID_.KE-1][*TRACK_mod_.MAT-1])* *CEGRID_.XEK);
    CJUMP0_.P[5-1]=exp(CEIMFP_.SEISI[*CEGRID_.KE-1][*TRACK_mod_.MAT-1]+(CEIMFP_.SEISI[*CEGRID_.KE+1-1][*TRACK_mod_.MAT-1]-CEIMFP_.SEISI[*CEGRID_.KE-1][*TRACK_mod_.MAT-1])* *CEGRID_.XEK);
    CJUMP0_.P[8-1]=0.0e0;

	if (IEND == 1){
		return;
	}

	if (CEIMFP_.W1E[*CEGRID_.KE+1-1][*TRACK_mod_.MAT-1] > -78.3e0){
		*CJUMP0_.W1=exp(CEIMFP_.W1E[*CEGRID_.KE-1][*TRACK_mod_.MAT-1]+(CEIMFP_.W1E[*CEGRID_.KE+1-1][*TRACK_mod_.MAT-1]-CEIMFP_.W1E[*CEGRID_.KE-1][*TRACK_mod_.MAT-1])* *CEGRID_.XEK);
        *CJUMP0_.W2=exp(CEIMFP_.W2E[*CEGRID_.KE-1][*TRACK_mod_.MAT-1]+(CEIMFP_.W2E[*CEGRID_.KE+1-1][*TRACK_mod_.MAT-1]-CEIMFP_.W2E[*CEGRID_.KE-1][*TRACK_mod_.MAT-1])* *CEGRID_.XEK);
	}else{
		*CJUMP0_.W1=0.0e0;
		*CJUMP0_.W2=0.0e0;
	}

	if (CEIMFP_.T1E[*CEGRID_.KE+1-1][*TRACK_mod_.MAT-1] > -78.3e0){
		*CJUMP0_.T1=exp(CEIMFP_.T1E[*CEGRID_.KE-1][*TRACK_mod_.MAT-1]+(CEIMFP_.T1E[*CEGRID_.KE+1-1][*TRACK_mod_.MAT-1]-CEIMFP_.T1E[*CEGRID_.KE-1][*TRACK_mod_.MAT-1])* *CEGRID_.XEK);
        *CJUMP0_.T2=exp(CEIMFP_.T2E[*CEGRID_.KE-1][*TRACK_mod_.MAT-1]+(CEIMFP_.T2E[*CEGRID_.KE+1-1][*TRACK_mod_.MAT-1]-CEIMFP_.T2E[*CEGRID_.KE-1][*TRACK_mod_.MAT-1])* *CEGRID_.XEK);
	}else{
		*CJUMP0_.T1=0.0e0;
		*CJUMP0_.T2=0.0e0;
	}
}

void gimfp2_(){
	/*
	Esta sub-rotina calcula os caminhos livres médios inversos para interações
	de fótons com a energia atual no material M.
	*/

	/*if (imprimiu==0){
		printf("\n\nGIMFP2\n\n");
		imprimiu++;
	}*/

	CJUMP0_.P[1-1]=CGIMFP_.SGRA[*CEGRID_.KE-1][*TRACK_mod_.MAT-1];
    CJUMP0_.P[2-1]=exp(CGIMFP_.SGCO[*CEGRID_.KE-1][*TRACK_mod_.MAT-1]+(CGIMFP_.SGCO[*CEGRID_.KE+1-1][*TRACK_mod_.MAT-1]-CGIMFP_.SGCO[*CEGRID_.KE-1][*TRACK_mod_.MAT-1])* *CEGRID_.XEK);
    CJUMP0_.P[3-1]=CGIMFP_.SGPH[*CEGRID_.KE-1][*TRACK_mod_.MAT-1];

	if (*TRACK_mod_.E < 1.023e6){
		CJUMP0_.P[4-1]=0.0e0;
	}else{
		CJUMP0_.P[4-1]=exp(CGIMFP_.SGPP[*CEGRID_.KE-1][*TRACK_mod_.MAT-1]+(CGIMFP_.SGPP[*CEGRID_.KE+1-1][*TRACK_mod_.MAT-1]-CGIMFP_.SGPP[*CEGRID_.KE-1][*TRACK_mod_.MAT-1])* *CEGRID_.XEK);
	}
	CJUMP0_.P[8-1]=0.0e0;
}

void fimdet2_(int &N, int &ID, double &DSEF){

	/*
	Distribuição de fluência de partículas dentro do
	detector. Apenas colisões discretas.
	
	*/
	/*if (imprimiu==0){
		printf("\n\nFIMDET2\n\n");
		imprimiu++;
	}*/
	
	double FSAFE=1.000000001e0;
	static const int NIDM=25;
	static const int NBEM=1000;

	int IE, IT;

	if (CIMDET_.LLE[ID-1] == 1){
		IE=1.0e0+(log(*TRACK_mod_.E)-CIMDET_.EL[ID-1])*CIMDET_.RBSE[ID-1];
	}else{
		IE=1.0e0+(*TRACK_mod_.E-CIMDET_.EL[ID-1])*CIMDET_.RBSE[ID-1];
	}

	if ((IE > 0) && (IE <= CIMDET_.NE[ID-1])){
		if (N != CIMDET_.LFLT[IE-1][ID-1]){
			CIMDET_.FLT[IE-1][ID-1]=CIMDET_.FLT[IE-1][ID-1]+CIMDET_.FLTP[IE-1][ID-1];
            CIMDET_.FLT2[IE-1][ID-1]=CIMDET_.FLT2[IE-1][ID-1]+pow(CIMDET_.FLTP[IE-1][ID-1],2);
            CIMDET_.FLTP[IE-1][ID-1]=*TRACK_mod_.WGHT * DSEF;
            CIMDET_.LFLT[IE-1][ID-1]=N;
		}else{
			CIMDET_.FLTP[IE-1][ID-1]=CIMDET_.FLTP[IE-1][ID-1]+*TRACK_mod_.WGHT*DSEF;
		}

		if (N != CIMDET_.LFLP[*TRACK_mod_.KPAR-1][IE-1][ID-1]){
			CIMDET_.FLP[*TRACK_mod_.KPAR-1][IE-1][ID-1]=CIMDET_.FLP[*TRACK_mod_.KPAR-1][IE-1][ID-1]+CIMDET_.FLPP[*TRACK_mod_.KPAR-1][IE-1][ID-1];
            CIMDET_.FLP2[*TRACK_mod_.KPAR-1][IE-1][ID-1]=CIMDET_.FLP2[*TRACK_mod_.KPAR-1][IE-1][ID-1]+pow(CIMDET_.FLPP[*TRACK_mod_.KPAR-1][IE-1][ID-1],2);
            CIMDET_.FLPP[*TRACK_mod_.KPAR-1][IE-1][ID-1]=*TRACK_mod_.WGHT*DSEF;
            CIMDET_.LFLP[*TRACK_mod_.KPAR-1][IE-1][ID-1]=N;
		}else{
			CIMDET_.FLPP[*TRACK_mod_.KPAR-1][IE-1][ID-1]=CIMDET_.FLPP[*TRACK_mod_.KPAR-1][IE-1][ID-1]+*TRACK_mod_.WGHT*DSEF;
		}
	}


}

void fimdes2_(int &N, int &ID, double &EI, double &DECSD, double &DSEF ){

	/*
	Distribuição de fluência de partículas dentro do
	detector C. Desaceleração contínua
	*/

/*	if (imprimiu==0){
		printf("\n\nFINDES2\n\n");
		imprimiu++;
	}*/

	double FSAFE=1.000000001e0;
	static const int NIDM=25;
	static const int NBEM=1000;

	int IE, IT, IEI, IEF;

	double EIC, EF, FACT, EA, EB, TLBIN;

	if (EI < CIMDET_.EL[ID-1])
		return;
	EIC=fmin(EI,CIMDET_.EU[ID-1]);
    EF=fmax(EI-DECSD,CIMDET_.EL[ID-1]);

	if (EF > EIC)
		return;

	if (CIMDET_.LLE[ID-1] == 1){
		IEI=1.0e0+(log(EIC)-CIMDET_.EL[ID-1])*CIMDET_.RBSE[ID-1];
        IEF=1.0e0+(log(EF)-CIMDET_.EL[ID-1])*CIMDET_.RBSE[ID-1];
	}else{
		IEI=1.0e0+(EIC-CIMDET_.EL[ID-1])*CIMDET_.RBSE[ID-1];
        IEF=1.0e0+(EF-CIMDET_.EL[ID-1])*CIMDET_.RBSE[ID-1];
	}

	FACT=DSEF/DECSD;

	for (int IE = IEF; IE <= IEI; IE++){
		EA=fmax(EF,CIMDET_.ET[IE-1][ID-1]);
        EB=fmin(EIC,CIMDET_.ET[IE+1-1][ID-1]);
        TLBIN=(EB-EA)*FACT;
		if (N != CIMDET_.LFLT[IE-1][ID-1]){
			CIMDET_.FLT[IE-1][ID-1]=CIMDET_.FLT[IE-1][ID-1]+CIMDET_.FLTP[IE-1][ID-1];
            CIMDET_.FLT2[IE-1][ID-1]=CIMDET_.FLT2[IE-1][ID-1]+pow(CIMDET_.FLTP[IE-1][ID-1],2);
            CIMDET_.FLTP[IE-1][ID-1]=*TRACK_mod_.WGHT * TLBIN;
            CIMDET_.LFLT[IE-1][ID-1]=N;
		}else{
			CIMDET_.FLTP[IE-1][ID-1]=CIMDET_.FLTP[IE-1][ID-1]+*TRACK_mod_.WGHT*TLBIN;
		}

		if (N != CIMDET_.LFLP[*TRACK_mod_.KPAR-1][IE-1][ID-1]){
			CIMDET_.FLP[*TRACK_mod_.KPAR-1][IE-1][ID-1]=CIMDET_.FLP[*TRACK_mod_.KPAR-1][IE-1][ID-1]+CIMDET_.FLPP[*TRACK_mod_.KPAR-1][IE-1][ID-1];
            CIMDET_.FLP2[*TRACK_mod_.KPAR-1][IE-1][ID-1]=CIMDET_.FLP2[*TRACK_mod_.KPAR-1][IE-1][ID-1]+pow(CIMDET_.FLPP[*TRACK_mod_.KPAR-1][IE-1][ID-1],2);
            CIMDET_.FLPP[*TRACK_mod_.KPAR-1][IE-1][ID-1]=*TRACK_mod_.WGHT*TLBIN;
            CIMDET_.LFLP[*TRACK_mod_.KPAR-1][IE-1][ID-1]=N;
		}else{
			CIMDET_.FLPP[*TRACK_mod_.KPAR-1][IE-1][ID-1]=CIMDET_.FLPP[*TRACK_mod_.KPAR-1][IE-1][ID-1]+*TRACK_mod_.WGHT*TLBIN;
		}
	}

}

void knock2_(double &DE, int &ICOL){

	/*
	Simulação de dobradiças aleatórias e eventos de interação difícil.

	Argumentos de saída:
	DE ..... energia depositada pela partícula no material. Isto é
	geralmente igual à diferença entre as energias
	antes e depois da interação.
	ICOL ... tipo de interação sofrida pela partícula.
	*/

/*	if (imprimiu==0){
		printf("\n\nKNOCK2\n\n");
		imprimiu++;
	}*/

	double PI=3.1415926535897932e0;
	double TWOPI=PI+PI; 
	double REV=5.10998928e5;
	double RREV=1.0e0/REV;
	double TREV=2.0e0*REV;

	double EMU1, EMU2, PNUM, PDEN, PMU0, PA, RND, CDT, DF, STNOW, STS, SS, TRNDC, TA, TB, RMU, DELTA, ES, EP, CDTS;
	double DFS, US, VS, WS, ECDT, CONS, EE, CDTE, CDTP, E1, CDT1, E2, CDT2;
	int IOSC, IZA, ISA, IEFF;

	static const int NO=512;
	static const int NOCO=512;
	static const int NBW=32;

	if (*TRACK_mod_.KPAR == 1)
		goto L1000;
	else if (*TRACK_mod_.KPAR == 2)
		goto L2000;
	else if (*TRACK_mod_.KPAR == 3)
		goto L3000;
	else{
		printf("   KNOCK: Incorrect particle type.\n");
		exit(0);
	}

L1000:;
	if (*CJUMP1_.MHINGE == 1)
		goto L1100;

	//Eletrons 
	//Dobradiça, evento suave artificial (ICOL=1).

	ICOL=1;
    *CJUMP1_.MHINGE=1;

	//Perda de Energia

	if (*CJUMP1_.KSOFTI == 1){
		DE=*PENELOPE_mod_.DESOFT;
        *TRACK_mod_.E=*TRACK_mod_.E-DE;
		if (*TRACK_mod_.E < PENELOPE_mod_.EABS[*TRACK_mod_.MAT-1][1-1]){
			DE=*PENELOPE_mod_.E0STEP;
            *TRACK_mod_.E=0.0e0;
			return;
		}
		*PENELOPE_mod_.E0STEP=*PENELOPE_mod_.E0STEP-*PENELOPE_mod_.SSOFT*(*CJUMP0_.DST-*CJUMP0_.DSR);
		if (*CJUMP1_.KSOFTE == 0)
			return;
		*CEGRID_.XEL=log(*PENELOPE_mod_.E0STEP);
        *CEGRID_.XE=1.0e0+(*CEGRID_.XEL-*CEGRID_.DLEMP1)* *CEGRID_.DLFC;
        *CEGRID_.KE=*CEGRID_.XE;
        *CEGRID_.XEK=*CEGRID_.XE-*CEGRID_.KE;
	}else{
		DE=0.0e0;
	}
	

	//Deflexão Angular
	if (CEIMFP_.T1E[*CEGRID_.KE+1-1][*TRACK_mod_.MAT-1] > -78.3e0){
		*CJUMP0_.T1=exp(CEIMFP_.T1E[*CEGRID_.KE-1][*TRACK_mod_.MAT-1]+(CEIMFP_.T1E[*CEGRID_.KE+1-1][*TRACK_mod_.MAT-1]-CEIMFP_.T1E[*CEGRID_.KE-1][*TRACK_mod_.MAT-1])* *CEGRID_.XEK);
        *CJUMP0_.T2=exp(CEIMFP_.T2E[*CEGRID_.KE-1][*TRACK_mod_.MAT-1]+(CEIMFP_.T2E[*CEGRID_.KE+1-1][*TRACK_mod_.MAT-1]-CEIMFP_.T2E[*CEGRID_.KE-1][*TRACK_mod_.MAT-1])* *CEGRID_.XEK);
	}else{
		*CJUMP0_.T1=0.0e0;
        *CJUMP0_.T2=0.0e0;
	}
	if (*CJUMP0_.T1 < 1.0e-20)
	    return;
	//1º e 2º momentos da distribuição angular.
	EMU1=0.5e0*(1.0e0-exp(- *CJUMP0_.DST* *CJUMP0_.T1));
    EMU2=EMU1-(1.0e0-exp(- *CJUMP0_.DST* *CJUMP0_.T2))/6.0e0;
	//Amostragem de um histograma de duas barras com esses momentos.
	PNUM=2.0e0*EMU1-3.0e0*EMU2;
    PDEN=1.0e0-2.0e0*EMU1;
    PMU0=PNUM/PDEN;
    PA=PDEN+PMU0;
    RND=rand2_(2.0e0);

	if (RND < PA){
		CDT=1.0e0-2.0e0*PMU0*(RND/PA);
	}else{
		CDT=1.0e0-2.0e0*(PMU0+(1.0e0-PMU0)*((RND-PA)/(1.0e0-PA)));
	}
	DF=TWOPI*rand2_(3.0e0);
	direct2_(CDT,DF,*TRACK_mod_.U,*TRACK_mod_.V,*TRACK_mod_.W);
	return;

	//Evento duro

L1100:;

	*CJUMP1_.MHINGE=0;
	//Uma interação delta (ICOL=7) ocorre quando o máximo comprimento de passo permitido é excedido.
	if (*CJUMP1_.KDELTA == 1){
		ICOL=7;
        DE=0.0e0;
		return;
	}

	//Amostragem aleatória do tipo de interação.
	STNOW=CJUMP0_.P[2-1]+CJUMP0_.P[3-1]+CJUMP0_.P[4-1]+CJUMP0_.P[5-1]+CJUMP0_.P[8-1];
    STS=fmax(STNOW,*CJUMP0_.ST)*rand2_(4.0e0);
    SS=CJUMP0_.P[2-1];
	if (SS > STS)
		goto L1200;

	SS=SS+CJUMP0_.P[3-1];
	if (SS > STS)
		goto L1300;

	SS=SS+CJUMP0_.P[4-1];
	if (SS > STS)
		goto L1400;

	SS=SS+CJUMP0_.P[5-1];
	if (SS > STS)
		goto L1500;

	SS=SS+CJUMP0_.P[8-1];
	if (SS > STS)
		goto L1800;

	/*
	Uma interação delta (ICOL=7) pode ocorrer quando o total
	A probabilidade de interação por unidade de comprimento do caminho, ST, é maior que STNOW.
	*/

	ICOL=7;
    DE=0.0e0;
	return;

	//Colisão Elastica Dura ICOL=2

L1200:;
	ICOL=2;
	if (*TRACK_mod_.E >= CELSEP_.EELMAX[*TRACK_mod_.MAT-1]){
		TRNDC=CEIMFP_.RNDCE[*CEGRID_.KE-1][*TRACK_mod_.MAT-1]+(CEIMFP_.RNDCE[*CEGRID_.KE+1-1][*TRACK_mod_.MAT-1]-CEIMFP_.RNDCE[*CEGRID_.KE-1][*TRACK_mod_.MAT-1])* *CEGRID_.XEK;
        TA=exp(CEIMFP_.AE[*CEGRID_.KE-1][*TRACK_mod_.MAT-1]+(CEIMFP_.AE[*CEGRID_.KE+1-1][*TRACK_mod_.MAT-1]-CEIMFP_.AE[*CEGRID_.KE-1][*TRACK_mod_.MAT-1])* *CEGRID_.XEK);
        TB=CEIMFP_.BE[*CEGRID_.KE-1][*TRACK_mod_.MAT-1]+(CEIMFP_.BE[*CEGRID_.KE+1-1][*TRACK_mod_.MAT-1]-CEIMFP_.BE[*CEGRID_.KE-1][*TRACK_mod_.MAT-1])* *CEGRID_.XEK;
        eela2_(TA,TB,TRNDC,RMU);
	} else{
		//Implementacao do modelo alternativo utilzando  ELSEPA database
		TRNDC=CELSEP_.RNDCED[*CEGRID_.KE-1][*TRACK_mod_.MAT-1]+(CELSEP_.RNDCED[*CEGRID_.KE+1-1][*TRACK_mod_.MAT-1]-CELSEP_.RNDCED[*CEGRID_.KE-1][*TRACK_mod_.MAT-1])* *CEGRID_.XEK;
        eeld2_(TRNDC,RMU);
	}

	CDT=1.0e0-(RMU+RMU);
    DF=TWOPI*rand2_(5.0e0);
	direct2_(CDT,DF,*TRACK_mod_.U,*TRACK_mod_.V,*TRACK_mod_.W);
	DE=0.0e0;
	return;

	//Colisão Dura inelastica (ICOL=3)

L1300:;
	ICOL=3;
    DELTA=CEIMFP_.DEL[*CEGRID_.KE-1][*TRACK_mod_.MAT-1]+(CEIMFP_.DEL[*CEGRID_.KE+1-1][*TRACK_mod_.MAT-1]-CEIMFP_.DEL[*CEGRID_.KE-1][*TRACK_mod_.MAT-1])* *CEGRID_.XEK;
	eina2_(*TRACK_mod_.E,DELTA,DE,EP,CDT,ES,CDTS,*TRACK_mod_.MAT,IOSC);

	//Ângulos de espalhamento (elétron primário).
	DF=TWOPI*rand2_(6.0e0);
	//Raio Delta
	if (ES > PENELOPE_mod_.EABS[*TRACK_mod_.MAT-1][1-1]){
		DFS=DF+PI;
        US=*TRACK_mod_.U;
        VS=*TRACK_mod_.V;
        WS=*TRACK_mod_.W;
        direct2_(CDTS,DFS,US,VS,WS);
        CHIST_.ILBA[1-1]=TRACK_mod_.ILB[1-1]+1;
        CHIST_.ILBA[2-1]=*TRACK_mod_.KPAR;
        CHIST_.ILBA[3-1]=ICOL;
        CHIST_.ILBA[4-1]=0;
        CHIST_.ILBA[5-1]=TRACK_mod_.ILB[5-1];
		int wvar = 1;
        stores2_(ES,*TRACK_mod_.X,*TRACK_mod_.Y,*TRACK_mod_.Z,US,VS,WS,*TRACK_mod_.WGHT,wvar,CHIST_.ILBA,wIPOLI);
	}
	//Nova energia e direção.
	if (EP > PENELOPE_mod_.EABS[*TRACK_mod_.MAT-1][1-1]){
		*TRACK_mod_.E=EP;
        direct2_(CDT,DF,*TRACK_mod_.U,*TRACK_mod_.V,*TRACK_mod_.W);
	}else{
		DE=*TRACK_mod_.E;
        *TRACK_mod_.E=0.0e0;
	}

	return;

	//Emissão bremsstrahlung dura (ICOL=4).
L1400:;
	ICOL=4;
	ebra2_(*TRACK_mod_.E, DE, *TRACK_mod_.MAT);

	//fóton de Bremsstrahlung.

	if (DE > PENELOPE_mod_.EABS[*TRACK_mod_.MAT-1][2-1]){
		ebraa2_(*TRACK_mod_.E,DE,CDTS,*TRACK_mod_.MAT);
        DFS=TWOPI*rand2_(7.0e0);
        US=*TRACK_mod_.U;
        VS=*TRACK_mod_.V;
        WS=*TRACK_mod_.W;
        direct2_(CDTS,DFS,US,VS,WS);
        CHIST_.ILBA[1-1]=TRACK_mod_.ILB[1-1]+1;
        CHIST_.ILBA[2-1]=*TRACK_mod_.KPAR;
        CHIST_.ILBA[3-1]=ICOL;
        CHIST_.ILBA[4-1]=0;
        CHIST_.ILBA[5-1]=TRACK_mod_.ILB[5-1];
		int wvar = 2;
        stores2_(DE,*TRACK_mod_.X,*TRACK_mod_.Y,*TRACK_mod_.Z,US,VS,WS,*TRACK_mod_.WGHT,wvar,CHIST_.ILBA,wIPOLI);
	}
	//Nova energia
	*TRACK_mod_.E=*TRACK_mod_.E-DE;
	if (*TRACK_mod_.E < PENELOPE_mod_.EABS[*TRACK_mod_.MAT-1][1-1]){
		DE=*TRACK_mod_.E+DE;
        *TRACK_mod_.E=0.0e0;
	}
	return;

	//Ionização de uma casca interna (ICOL=5).

L1500:;
	ICOL=5;
	DELTA=CEIMFP_.DEL[*CEGRID_.KE-1][*TRACK_mod_.MAT-1]+(CEIMFP_.DEL[*CEGRID_.KE+1-1][*TRACK_mod_.MAT-1]-CEIMFP_.DEL[*CEGRID_.KE-1][*TRACK_mod_.MAT-1])* *CEGRID_.XEK;
	esia2_(*TRACK_mod_.E,DELTA,DE,EP,CDT,ES,CDTS,*TRACK_mod_.MAT,IZA,ISA);
	//relaxamento atomico
	if (IZA > 2){
		CHIST_.ILBA[3-1]=ICOL;
        relax2_(IZA,ISA);
	}

	//Ângulos de espalhamento (elétron primário).
	DF=TWOPI*rand2_(8.0e0);
	//raio delta
	if (ES > PENELOPE_mod_.EABS[*TRACK_mod_.MAT-1][1-1]){
		DFS=DF+PI;
        US=*TRACK_mod_.U;
        VS=*TRACK_mod_.V;
        WS=*TRACK_mod_.W;
        direct2_(CDTS,DFS,US,VS,WS);
        CHIST_.ILBA[1-1]=TRACK_mod_.ILB[1-1]+1;
        CHIST_.ILBA[2-1]=*TRACK_mod_.KPAR;
        CHIST_.ILBA[3-1]=ICOL;
        CHIST_.ILBA[4-1]=0;
        CHIST_.ILBA[5-1]=TRACK_mod_.ILB[5-1];
		int wvar = 1;
        stores2_(ES,*TRACK_mod_.X,*TRACK_mod_.Y,*TRACK_mod_.Z,US,VS,WS,*TRACK_mod_.WGHT,wvar,CHIST_.ILBA,wIPOLI);
	}

	//Nova energia e direção
	if (EP > PENELOPE_mod_.EABS[*TRACK_mod_.MAT-1][1-1]){
		*TRACK_mod_.E=EP;
        direct2_(CDT,DF,*TRACK_mod_.U,*TRACK_mod_.V,*TRACK_mod_.W);
	}else{
		DE=*TRACK_mod_.E;
        *TRACK_mod_.E=0.0e0;
	}
	return;

	//Mecanismo fictício auxiliar (ICOL=8).

L1800:;
	ICOL=8;
	DE=0.0e0;
	eaux2_();
	return;

	//Fotons KPAR=2

L2000:;

	STS=*CJUMP0_.ST*rand2_(1.0e0);
	SS=CJUMP0_.P[1-1];
	if (SS > STS)
		goto L2100;

	SS=SS+CJUMP0_.P[2-1];
	if (SS > STS)
		goto L2200;

	SS=SS+CJUMP0_.P[3-1];
	if (SS > STS)
		goto L2300;

	SS=SS+CJUMP0_.P[4-1];
	if (SS > STS)
		goto L2400;

	SS=SS+CJUMP0_.P[8-1];
	if (SS > STS)
		goto L2800;


	//Espalhamento Rayleigh ICOL=1

L2100:;
	DE=0.0e0;
	graa2_(*TRACK_mod_.E,CDT,IEFF,*TRACK_mod_.MAT);

	/*
	Interação delta. Introduzido para corrigir o uso de um
	limite superior do coeficiente de atenuação de Rayleigh.
	*/

	if (IEFF == 0){
		ICOL=7;
		return;
	}
	ICOL=1;
	if (*TRACK_mod_.IPOL == 1){
		double wCONS = 0.0e0;
		dirpol2_(CDT,DF,wCONS,*TRACK_mod_.SP1,*TRACK_mod_.SP2,*TRACK_mod_.SP3,*TRACK_mod_.U,*TRACK_mod_.V,*TRACK_mod_.W);
	}else{
		DF=TWOPI*rand2_(2.0e0);
		direct2_(CDT,DF,*TRACK_mod_.U,*TRACK_mod_.V,*TRACK_mod_.W);
	}

	TRACK_mod_.ILB[1-1]=TRACK_mod_.ILB[1-1]+1;
    TRACK_mod_.ILB[2-1]=*TRACK_mod_.KPAR;
    TRACK_mod_.ILB[3-1]=ICOL;
	return;

	//Espalhamento Compton ICOL=2

L2200:;
	ICOL=2;
	gcoa2_(*TRACK_mod_.E,DE,EP,CDT,ES,CDTS,*TRACK_mod_.MAT,IZA,ISA);
	US=*TRACK_mod_.U;
    VS=*TRACK_mod_.V;
    WS=*TRACK_mod_.W;
    DF=-1.0e0;
	if ((IZA > 0) && (ISA < 17)){
		CHIST_.ILBA[3-1]=ICOL;
		relax2_(IZA,ISA);
	}

	//Nova direção e energia
	if (EP > PENELOPE_mod_.EABS[*TRACK_mod_.MAT-1][2-1]){
		if (*TRACK_mod_.IPOL == 1){
			ECDT=*TRACK_mod_.E*RREV*(1.0e0-CDT);
          	CONS=ECDT*ECDT/(1.0e0+ECDT);
         	dirpol2_(CDT,DF,CONS,*TRACK_mod_.SP1,*TRACK_mod_.SP2,*TRACK_mod_.SP3,*TRACK_mod_.U,*TRACK_mod_.V,*TRACK_mod_.W);
		}else{
			DF=TWOPI*rand2_(3.0e0);
			direct2_(CDT,DF,*TRACK_mod_.U,*TRACK_mod_.V,*TRACK_mod_.W);
		}
		*TRACK_mod_.E=EP;
	}else{
		DE=*TRACK_mod_.E;
        *TRACK_mod_.E=0.0e0;
	}

	//Electron Compton - particula secundaria
	if (ES > PENELOPE_mod_.EABS[*TRACK_mod_.MAT-1][1-1]){
		if (DF < -0.5e0)
			DF=TWOPI*rand2_(4.0e0);
		DFS=DF+PI;
        direct2_(CDTS,DFS,US,VS,WS);
        CHIST_.ILBA[1-1]=TRACK_mod_.ILB[1-1]+1;
        CHIST_.ILBA[2-1]=*TRACK_mod_.KPAR;
        CHIST_.ILBA[3-1]=ICOL;
        CHIST_.ILBA[4-1]=0;
        CHIST_.ILBA[5-1]=TRACK_mod_.ILB[5-1];
		int wKPARP = 1;
        stores2_(ES,*TRACK_mod_.X,*TRACK_mod_.Y,*TRACK_mod_.Z,US,VS,WS,*TRACK_mod_.WGHT,wKPARP,CHIST_.ILBA,wIPOLI);
	}
	TRACK_mod_.ILB[1-1]=TRACK_mod_.ILB[1-1]+1;
    TRACK_mod_.ILB[2-1]=*TRACK_mod_.KPAR;
    TRACK_mod_.ILB[3-1]=ICOL;
	return;

	//Absorção Fotoeletrica ICOL=3

L2300:;
	ICOL=3;
	gpha2_(ES,IZA,ISA);
	/*
	Interação delta. Introduzido para corrigir o uso de um
	limite superior do coeficiente de atenuação fotoelétrica.
	*/

	if (IZA == 0){
		ICOL=7;
        DE=0.0e0;
        return;
	}

	if (ES > PENELOPE_mod_.EABS[*TRACK_mod_.MAT-1][1-1]){
		sauter2_(ES,CDTS);
        DFS=TWOPI*rand2_(5.0e0);
        US=*TRACK_mod_.U;
        VS=*TRACK_mod_.V;
        WS=*TRACK_mod_.W;
		direct2_(CDTS,DFS,US,VS,WS);
		CHIST_.ILBA[1-1]=TRACK_mod_.ILB[1-1]+1;
        CHIST_.ILBA[2-1]=*TRACK_mod_.KPAR;
        CHIST_.ILBA[3-1]=ICOL;
        CHIST_.ILBA[4-1]=0;
        CHIST_.ILBA[5-1]=TRACK_mod_.ILB[5-1];
		int wKPARP = 1;
        stores2_(ES,*TRACK_mod_.X,*TRACK_mod_.Y,*TRACK_mod_.Z,US,VS,WS,*TRACK_mod_.WGHT,wKPARP,CHIST_.ILBA,wIPOLI);
	}
	if (ISA < 17){
		CHIST_.ILBA[3-1]=ICOL;
        relax2_(IZA,ISA);
	}
	DE=*TRACK_mod_.E;
    *TRACK_mod_.E=0.0e0;
	return;


	//Produção de pares eletron-positron ICOL=4

L2400:;

	ICOL=4;
	gppa2_(EE,CDTE,EP,CDTP,IZA,ISA);
	DE=*TRACK_mod_.E;
	//Eletron
	if (EE >  PENELOPE_mod_.EABS[*TRACK_mod_.MAT-1][1-1]){
		DF=TWOPI*rand2_(6.0e0);
        US=*TRACK_mod_.U;
        VS=*TRACK_mod_.V;
        WS=*TRACK_mod_.W;
		direct2_(CDTE,DF,US,VS,WS);
		CHIST_.ILBA[1-1]=TRACK_mod_.ILB[1-1]+1;
        CHIST_.ILBA[2-1]=*TRACK_mod_.KPAR;
        CHIST_.ILBA[3-1]=ICOL;
        CHIST_.ILBA[4-1]=0;
        CHIST_.ILBA[5-1]=TRACK_mod_.ILB[5-1];
		int wKPARP = 1;
        stores2_(EE,*TRACK_mod_.X,*TRACK_mod_.Y,*TRACK_mod_.Z,US,VS,WS,*TRACK_mod_.WGHT,wKPARP,CHIST_.ILBA,wIPOLI);
	}
	//Positron
	if (EP >  PENELOPE_mod_.EABS[*TRACK_mod_.MAT-1][3-1]){
		DF=TWOPI*rand2_(7.0e0);
        US=*TRACK_mod_.U;
        VS=*TRACK_mod_.V;
        WS=*TRACK_mod_.W;
		direct2_(CDTP,DF,US,VS,WS);
		CHIST_.ILBA[1-1]=TRACK_mod_.ILB[1-1]+1;
        CHIST_.ILBA[2-1]=*TRACK_mod_.KPAR;
        CHIST_.ILBA[3-1]=ICOL;
        CHIST_.ILBA[4-1]=0;
        CHIST_.ILBA[5-1]=TRACK_mod_.ILB[5-1];
		int wKPARP = 3;
        stores2_(EP,*TRACK_mod_.X,*TRACK_mod_.Y,*TRACK_mod_.Z,US,VS,WS,*TRACK_mod_.WGHT,wKPARP,CHIST_.ILBA,wIPOLI);
		//O pósitron carrega uma energia 'latente' de 1022 keV.
		DE=DE-TREV;
	}else{
		panar2_(PENELOPE_mod_.EABS[*TRACK_mod_.MAT-1][2-1]);
	}
	*TRACK_mod_.E=0.0e0;

	//Relaxamento atômico após a produção de tripletos.
	if (ISA < 17){
		CHIST_.ILBA[3-1]=ICOL;
        relax2_(IZA,ISA);
	}
	return;

	//Mecanismo fictício auxiliar (ICOL=8).

L2800:;

	ICOL=8;
	DE=0.0e0;
	gaux2_();
	return;

	//Positrons KPAR=3


L3000:;	

	if (*CJUMP1_.MHINGE == 1)
		goto L3100;


	//Dobradiça, evento suave artificial (ICOL=1).

	ICOL=1;
    *CJUMP1_.MHINGE=1;

	//Perda de Energia

	if (*CJUMP1_.KSOFTI == 1){
		DE=*PENELOPE_mod_.DESOFT;
        *TRACK_mod_.E=*TRACK_mod_.E-DE;
		if (*TRACK_mod_.E < PENELOPE_mod_.EABS[*TRACK_mod_.MAT-1][3-1]){
			panar2_(PENELOPE_mod_.EABS[*TRACK_mod_.MAT-1][2-1]); //Aniquilação em repouso.
			DE=*PENELOPE_mod_.E0STEP+TREV;
            *TRACK_mod_.E=0.0e0;
			return;
		}
		*PENELOPE_mod_.E0STEP=*PENELOPE_mod_.E0STEP-*PENELOPE_mod_.SSOFT*(*CJUMP0_.DST-*CJUMP0_.DSR);
		if (*CJUMP1_.KSOFTE == 0)
			return;

		*CEGRID_.XEL=log(*PENELOPE_mod_.E0STEP);
        *CEGRID_.XE=1.0e0+(*CEGRID_.XEL-*CEGRID_.DLEMP1)* *CEGRID_.DLFC;
        *CEGRID_.KE=*CEGRID_.XE;
        *CEGRID_.XEK=*CEGRID_.XE-*CEGRID_.KE;
	}else{
		DE=0.0e0;
	}
	

	//Deflexão Angular
	if (CPIMFP_.T1P[*CEGRID_.KE+1-1][*TRACK_mod_.MAT-1] > -78.3e0){
		*CJUMP0_.T1=exp(CPIMFP_.T1P[*CEGRID_.KE-1][*TRACK_mod_.MAT-1]+(CPIMFP_.T1P[*CEGRID_.KE+1-1][*TRACK_mod_.MAT-1]-CPIMFP_.T1P[*CEGRID_.KE-1][*TRACK_mod_.MAT-1])* *CEGRID_.XEK);
        *CJUMP0_.T2=exp(CPIMFP_.T2P[*CEGRID_.KE-1][*TRACK_mod_.MAT-1]+(CPIMFP_.T2P[*CEGRID_.KE+1-1][*TRACK_mod_.MAT-1]-CPIMFP_.T2P[*CEGRID_.KE-1][*TRACK_mod_.MAT-1])* *CEGRID_.XEK);
	}else{
		*CJUMP0_.T1=0.0e0;
        *CJUMP0_.T2=0.0e0;
	}
	if (*CJUMP0_.T1 < 1.0e-20)
	    return;
	//1º e 2º momentos da distribuição angular.
	EMU1=0.5e0*(1.0e0-exp(- *CJUMP0_.DST* *CJUMP0_.T1));
    EMU2=EMU1-(1.0e0-exp(- *CJUMP0_.DST* *CJUMP0_.T2))/6.0e0;
	//Amostragem de um histograma de duas barras com esses momentos.
	PNUM=2.0e0*EMU1-3.0e0*EMU2;
    PDEN=1.0e0-2.0e0*EMU1;
    PMU0=PNUM/PDEN;
    PA=PDEN+PMU0;
    RND=rand2_(2.0e0);

	if (RND < PA){
		CDT=1.0e0-2.0e0*PMU0*(RND/PA);
	}else{
		CDT=1.0e0-2.0e0*(PMU0+(1.0e0-PMU0)*((RND-PA)/(1.0e0-PA)));
	}
	DF=TWOPI*rand2_(3.0e0);
	direct2_(CDT,DF,*TRACK_mod_.U,*TRACK_mod_.V,*TRACK_mod_.W);
	return;

	//Evento duro

L3100:;

	*CJUMP1_.MHINGE=0;
	//Uma interação delta (ICOL=7) ocorre quando o máximo comprimento de passo permitido é excedido.
	if (*CJUMP1_.KDELTA == 1){
		ICOL=7;
        DE=0.0e0;
		return;
	}

	//Amostragem aleatória do tipo de interação.
	STNOW=CJUMP0_.P[2-1]+CJUMP0_.P[3-1]+CJUMP0_.P[4-1]+CJUMP0_.P[5-1]+CJUMP0_.P[6-1]+CJUMP0_.P[8-1];
    STS=fmax(STNOW,*CJUMP0_.ST)*rand2_(4.0e0);
    SS=CJUMP0_.P[2-1];
	if (SS > STS)
		goto L3200;

	SS=SS+CJUMP0_.P[3-1];
	if (SS > STS)
		goto L3300;

	SS=SS+CJUMP0_.P[4-1];
	if (SS > STS)
		goto L3400;

	SS=SS+CJUMP0_.P[5-1];
	if (SS > STS)
		goto L3500;

	SS=SS+CJUMP0_.P[6-1];
	if (SS > STS)
		goto L3600;

	SS=SS+CJUMP0_.P[8-1];
	if (SS > STS)
		goto L3800;

	/*
	Uma interação delta (ICOL=7) pode ocorrer quando o total
	A probabilidade de interação por unidade de comprimento do caminho, ST, é maior que STNOW.
	*/

	ICOL=7;
    DE=0.0e0;
	return;

	//Colisão Elastica Dura ICOL=2

L3200:;
	ICOL=2;
	if (*TRACK_mod_.E >= CELSEP_.PELMAX[*TRACK_mod_.MAT-1]){
		TRNDC=CPIMFP_.RNDCP[*CEGRID_.KE-1][*TRACK_mod_.MAT-1]+(CPIMFP_.RNDCP[*CEGRID_.KE+1-1][*TRACK_mod_.MAT-1]-CPIMFP_.RNDCP[*CEGRID_.KE-1][*TRACK_mod_.MAT-1])* *CEGRID_.XEK;
        TA=exp(CPIMFP_.AP[*CEGRID_.KE-1][*TRACK_mod_.MAT-1]+(CPIMFP_.AP[*CEGRID_.KE+1-1][*TRACK_mod_.MAT-1]-CPIMFP_.AP[*CEGRID_.KE-1][*TRACK_mod_.MAT-1])* *CEGRID_.XEK);
        TB=CPIMFP_.BP[*CEGRID_.KE-1][*TRACK_mod_.MAT-1]+(CPIMFP_.BP[*CEGRID_.KE+1-1][*TRACK_mod_.MAT-1]-CPIMFP_.BP[*CEGRID_.KE-1][*TRACK_mod_.MAT-1])* *CEGRID_.XEK;
        eela2_(TA,TB,TRNDC,RMU);
	} else{
		//Implementacao do modelo alternativo utilzando  ELSEPA database
		TRNDC=CELSEP_.RNDCPD[*CEGRID_.KE-1][*TRACK_mod_.MAT-1]+(CELSEP_.RNDCPD[*CEGRID_.KE+1-1][*TRACK_mod_.MAT-1]-CELSEP_.RNDCPD[*CEGRID_.KE-1][*TRACK_mod_.MAT-1])* *CEGRID_.XEK;
        peld2_(TRNDC,RMU);
	}

	CDT=1.0e0-(RMU+RMU);
    DF=TWOPI*rand2_(5.0e0);
	direct2_(CDT,DF,*TRACK_mod_.U,*TRACK_mod_.V,*TRACK_mod_.W);
	DE=0.0e0;
	return;

	//Colisão Dura inelastica (ICOL=3)

L3300:;
	ICOL=3;
    DELTA=CEIMFP_.DEL[*CEGRID_.KE-1][*TRACK_mod_.MAT-1]+(CEIMFP_.DEL[*CEGRID_.KE+1-1][*TRACK_mod_.MAT-1]-CEIMFP_.DEL[*CEGRID_.KE-1][*TRACK_mod_.MAT-1])* *CEGRID_.XEK;
	pina2_(*TRACK_mod_.E,DELTA,DE,EP,CDT,ES,CDTS,*TRACK_mod_.MAT,IOSC);

	//Ângulos de espalhamento (positron primário).
	DF=TWOPI*rand2_(6.0e0);
	//Raio Delta
	if (ES > PENELOPE_mod_.EABS[*TRACK_mod_.MAT-1][1-1]){
		DFS=DF+PI;
        US=*TRACK_mod_.U;
        VS=*TRACK_mod_.V;
        WS=*TRACK_mod_.W;
        direct2_(CDTS,DFS,US,VS,WS);
        CHIST_.ILBA[1-1]=TRACK_mod_.ILB[1-1]+1;
        CHIST_.ILBA[2-1]=*TRACK_mod_.KPAR;
        CHIST_.ILBA[3-1]=ICOL;
        CHIST_.ILBA[4-1]=0;
        CHIST_.ILBA[5-1]=TRACK_mod_.ILB[5-1];
		int wvar = 1;
        stores2_(ES,*TRACK_mod_.X,*TRACK_mod_.Y,*TRACK_mod_.Z,US,VS,WS,*TRACK_mod_.WGHT,wvar,CHIST_.ILBA,wIPOLI);
	}
	//Nova energia e direção.
	if (EP > PENELOPE_mod_.EABS[*TRACK_mod_.MAT-1][3-1]){
		*TRACK_mod_.E=EP;
        direct2_(CDT,DF,*TRACK_mod_.U,*TRACK_mod_.V,*TRACK_mod_.W);
	}else{
		panar2_(PENELOPE_mod_.EABS[*TRACK_mod_.MAT-1][2-1]);
		DE=*TRACK_mod_.E+TREV;
        *TRACK_mod_.E=0.0e0;
	}

	return;

	//Emissão bremsstrahlung dura (ICOL=4).
L3400:;
	ICOL=4;
	ebra2_(*TRACK_mod_.E, DE, *TRACK_mod_.MAT);

	//fóton de Bremsstrahlung.

	if (DE > PENELOPE_mod_.EABS[*TRACK_mod_.MAT-1][2-1]){
		ebraa2_(*TRACK_mod_.E,DE,CDTS,*TRACK_mod_.MAT);
        DFS=TWOPI*rand2_(7.0e0);
        US=*TRACK_mod_.U;
        VS=*TRACK_mod_.V;
        WS=*TRACK_mod_.W;
        direct2_(CDTS,DFS,US,VS,WS);
        CHIST_.ILBA[1-1]=TRACK_mod_.ILB[1-1]+1;
        CHIST_.ILBA[2-1]=*TRACK_mod_.KPAR;
        CHIST_.ILBA[3-1]=ICOL;
        CHIST_.ILBA[4-1]=0;
        CHIST_.ILBA[5-1]=TRACK_mod_.ILB[5-1];
		int wvar = 2;
        stores2_(DE,*TRACK_mod_.X,*TRACK_mod_.Y,*TRACK_mod_.Z,US,VS,WS,*TRACK_mod_.WGHT,wvar,CHIST_.ILBA,wIPOLI);
	}
	//Nova energia
	*TRACK_mod_.E=*TRACK_mod_.E-DE;
	if (*TRACK_mod_.E < PENELOPE_mod_.EABS[*TRACK_mod_.MAT-1][3-1]){ 
		panar2_(PENELOPE_mod_.EABS[*TRACK_mod_.MAT-1][2-1]);//Aniquilação em repouso.
		DE=*TRACK_mod_.E+DE+TREV;
        *TRACK_mod_.E=0.0e0;
	}
	return;

	//Ionização de uma casca interna (ICOL=5).

L3500:;
	ICOL=5;
	DELTA=CEIMFP_.DEL[*CEGRID_.KE-1][*TRACK_mod_.MAT-1]+(CEIMFP_.DEL[*CEGRID_.KE+1-1][*TRACK_mod_.MAT-1]-CEIMFP_.DEL[*CEGRID_.KE-1][*TRACK_mod_.MAT-1])* *CEGRID_.XEK;
	psia2_(*TRACK_mod_.E,DELTA,DE,EP,CDT,ES,CDTS,*TRACK_mod_.MAT,IZA,ISA);
	//relaxamento atomico
	if (IZA > 2){
		CHIST_.ILBA[3-1]=ICOL;
        relax2_(IZA,ISA);
	}

	//Ângulos de espalhamento (elétron primário).
	DF=TWOPI*rand2_(8.0e0);
	//raio delta
	if (ES > PENELOPE_mod_.EABS[*TRACK_mod_.MAT-1][1-1]){
		DFS=DF+PI;
        US=*TRACK_mod_.U;
        VS=*TRACK_mod_.V;
        WS=*TRACK_mod_.W;
        direct2_(CDTS,DFS,US,VS,WS);
        CHIST_.ILBA[1-1]=TRACK_mod_.ILB[1-1]+1;
        CHIST_.ILBA[2-1]=*TRACK_mod_.KPAR;
        CHIST_.ILBA[3-1]=ICOL;
        CHIST_.ILBA[4-1]=0;
        CHIST_.ILBA[5-1]=TRACK_mod_.ILB[5-1];
		int wvar = 1;
        stores2_(ES,*TRACK_mod_.X,*TRACK_mod_.Y,*TRACK_mod_.Z,US,VS,WS,*TRACK_mod_.WGHT,wvar,CHIST_.ILBA,wIPOLI);
	}

	//Nova energia e direção
	if (EP > PENELOPE_mod_.EABS[*TRACK_mod_.MAT-1][3-1]){
		*TRACK_mod_.E=EP;
        direct2_(CDT,DF,*TRACK_mod_.U,*TRACK_mod_.V,*TRACK_mod_.W);
	}else{
		panar2_(PENELOPE_mod_.EABS[*TRACK_mod_.MAT-1][2-1]);
		DE=*TRACK_mod_.E+TREV;
        *TRACK_mod_.E=0.0e0;
	}
	return;


L3600:; //Aniquilação de positron em voo

	ICOL=6;
	pana2_(*TRACK_mod_.E,E1,CDT1,E2,CDT2,*TRACK_mod_.MAT);
	DF=TWOPI*rand2_(9.0e0);
	if (E1 > PENELOPE_mod_.EABS[*TRACK_mod_.MAT-1][2-1]){
		US=*TRACK_mod_.U;
        VS=*TRACK_mod_.V;
        WS=*TRACK_mod_.W;
        direct2_(CDT1,DF,US,VS,WS);
        CHIST_.ILBA[1-1]=TRACK_mod_.ILB[1-1]+1;
        CHIST_.ILBA[2-1]=*TRACK_mod_.KPAR;
        CHIST_.ILBA[3-1]=ICOL;
        CHIST_.ILBA[4-1]=0;
        CHIST_.ILBA[5-1]=TRACK_mod_.ILB[5-1];
		int wvar = 2;
        stores2_(E1,*TRACK_mod_.X,*TRACK_mod_.Y,*TRACK_mod_.Z,US,VS,WS,*TRACK_mod_.WGHT,wvar,CHIST_.ILBA,wIPOLI);
	}
	if (E2 > PENELOPE_mod_.EABS[*TRACK_mod_.MAT-1][2-1]){
		DF=DF+PI;
		US=*TRACK_mod_.U;
        VS=*TRACK_mod_.V;
        WS=*TRACK_mod_.W;
        direct2_(CDT2,DF,US,VS,WS);
        CHIST_.ILBA[1-1]=TRACK_mod_.ILB[1-1]+1;
        CHIST_.ILBA[2-1]=*TRACK_mod_.KPAR;
        CHIST_.ILBA[3-1]=ICOL;
        CHIST_.ILBA[4-1]=0;
        CHIST_.ILBA[5-1]=TRACK_mod_.ILB[5-1];
		int wvar = 2;
        stores2_(E2,*TRACK_mod_.X,*TRACK_mod_.Y,*TRACK_mod_.Z,US,VS,WS,*TRACK_mod_.WGHT,wvar,CHIST_.ILBA,wIPOLI);
	}
	DE=*TRACK_mod_.E+TREV;
    *TRACK_mod_.E=0.0e0;
	return;


	//Mecanismo fictício auxiliar (ICOL=8).

L3800:;
	ICOL=8;
	DE=0.0e0;
	paux2_();
	return;

}

void eela2_(double &A, double &B, double &RNDC, double &RMU){

	/*
	Simulação de eventos elásticos duros. Modelo Wentzel modificado.

	Argumentos de entrada:
	A, B ... parâmetros de distribuição angular.
	RNDC ... probabilidade de corte.
	Valores de saída :
	RMU .... deflexão angular, =(1-CDT)/2.
	*/

/*	if (imprimiu==0){
		printf("\n\nEELA2\n\n");
		imprimiu++;
	}*/

	double A1, B1, RMUAV, RND0, RND, RNDMB, BB, RMUC, PW, RNDRC;

	A1=A+1.0e0;

	if (B >= 0.0e0){
		//Caso I

		RMUAV=A*A1*log(A1/A)-A;
        B1=1.0e0-B;
        RND0=B1*A1*RMUAV/(A+RMUAV);
        RND=RNDC+rand2_(1.0e0)*(1.0e0-RNDC);

		if (RND < RND0){
			RMU=RND*A/(B1*A1-RND);
		}else if (RND > RND0+B){
			RNDMB=RND-B;
          	RMU=RNDMB*A/(B1*A1-RNDMB);
		}else{
			RMU=RMUAV;
		}
	}else{
		//Caso II
		BB=-B;
        B1=1.0e0-BB;
        RMUC=RNDC*A/(B1*A1-RNDC);
        PW=B1*A*(1.0e0-RMUC)/(A+RMUC);
		if (rand2_(2.0e0)*(BB+PW) < BB){
			RMU=0.5e0*(1.0e0+sqrt(rand2_(3.0e0)));
		}else{
			RNDRC=rand2_(3.0e0)*(1.0e0-RMUC);
          	RMU=(A*RNDRC+A1*RMUC)/(A1-RNDRC);
		}
	}

}

void eeld2_(double &RNDC, double &RMU){

	/*
	Simulação de eventos elásticos rígidos de elétrons. Seções transversais de
 	a base de dados numérica ELSEPA.

	Valor do argumento C:
 	RNDC ... valor de corte do número aleatório uniforme
 	(apenas eventos difíceis são simulados).
 	RMU .... deflexão angular amostrada, =(1-CDT)/2.
	
	*/

	/*if (imprimiu==0){
		printf("\n\nEELD2\n\n");
		imprimiu++;
	}*/

	static const int NP=128;
	static const int NPM1=NP-1;

	double PK, RU, RR, PP, XX, AA, BB, D;
	int ITN, JE, I, J, K;

	//Ponto da rede de energia


	PK=(*CEGRID_.XEL-CEGRID_.DLEMP[*CEGRID_.KE-1])* *CEGRID_.DLFC;

	if (rand2_(1.0e0) < PK){
		JE=*CEGRID_.KE+1;
	}else{
		JE=*CEGRID_.KE;
	}

	//Ponto
	RU=RNDC+rand2_(2.0e0)*(1.0e0-RNDC);

	//Selection of the interval (binary search in a restricted interval).
	ITN=RU*NPM1+1;
    I=CEELDB_.ITLE[*TRACK_mod_.MAT-1][JE-1][ITN-1];
    J=CEELDB_.ITUE[*TRACK_mod_.MAT-1][JE-1][ITN-1];
	if ((J-I) < 2)
		goto L2;
L1:;
	K=(I+J)/2;
	if (RU > CEELDB_.PSE[*TRACK_mod_.MAT-1][JE-1][K-1]){
		I=K;
	}else{
		J=K;
	}

	if ((J-I) > 1)
		goto L1;

	//Amostragem da distribuição cumulativa inversa racional.
L2:;
	PP=CEELDB_.PSE[*TRACK_mod_.MAT-1][JE-1][I-1];
    RR=RU-PP;
	if (RR > 1.0e-16){
		XX=CEELDB_.XSE[*TRACK_mod_.MAT-1][JE-1][I-1];
        AA=CEELDB_.ASE[*TRACK_mod_.MAT-1][JE-1][I-1];
        BB=CEELDB_.BSE[*TRACK_mod_.MAT-1][JE-1][I-1];
        D=CEELDB_.PSE[*TRACK_mod_.MAT-1][JE-1][I+1-1]-PP;
        RMU=XX+((1.0e0+AA+BB)*D*RR/(D*D+(AA*D+BB*RR)*RR))*(CEELDB_.XSE[*TRACK_mod_.MAT-1][JE-1][I+1-1]-XX);
	}else{
		RMU=CEELDB_.XSE[*TRACK_mod_.MAT-1][JE-1][I-1];
	}
}

void eina2_(double &E, double &DELTA, double &DE, double &EP, double &CDT, double &ES, double &CDTS, int &M, int &IOSC){

	/*
	Amostragem aleatória de colisões inelásticas duras de elétrons.

	Modelo C Sternheimer-Liljequist GOS.

	Argumentos de entrada:
	E ....... energia do elétron (eV).
	M ....... material onde os elétrons se propagam.
	DELTA ... Correção do efeito de densidade de Fermi.
	Argumentos de saída:
	DE ...... perda de energia (eV).
	EP ...... energia do elétron espalhado (eV).
	CDT ..... cosseno do ângulo de dispersão polar.
	ES ...... energia do elétron secundário emitido (eV).
	CDTS .... cosseno polar de direção do elétron secundário.
	IOSC .... índice do oscilador que foi 'ionizado'.
	
	*/

	/*if (imprimiu==0){
		printf("\n\nEINA2\n\n");
		imprimiu++;
	}*/

	double REV=5.10998928e5;
	bool LDIST;
	double RREV=1.0e0/REV;
	double TREV=2.0e0*REV;
	double RTREV=1.0e0/TREV;

	static const int NO=512;

	double WCCM, PK, TST, UK, WK,WTHR, WM, WKP, QKP, EE, WCMAX, WDMAX, RB, GAM, GAM2, BETA2, AMOL, CPS, CP;
	double CPPS, CPP, QM, RWKP, XHDL, XHDT, F0, RCL, RL1, RRL1, XHC, XHTOT, TS1, A, ARCL, FB, RK, RK2, RKF, PHI;
	double QS, Q, QTREV;
	int IO, JO, IT, JE;

	WCCM=PENELOPE_mod_.WCC[M-1];

	if (WCCM > E){
		DE=0.0e0;
        EP=E;
        CDT=1.0e0;
        ES=0.0e0;
        CDTS=0.0e0;
        IOSC=NO;
        return;
	}
	//Ponto da rede de energia
	PK=(*CEGRID_.XEL-CEGRID_.DLEMP[*CEGRID_.KE-1])* *CEGRID_.DLFC;
	if (rand2_(1.0e0) < PK){
		JE=*CEGRID_.KE+1;
	}else{
		JE=*CEGRID_.KE;
	}

	//Seleção do oscilador ativo.
	TST=rand2_(2.0e0);
	//Busca binaria
	IO=1;
    JO=CEINAC_.NEIN[M-1]+1;
L1:;
    IT=(IO+JO)/2;
	if (TST > CEINAC_.EINAC[IT-1][JE-1][M-1]){
		IO=IT;
	}else{
		JO=IT;
	}
	if (JO-IO > 1)
		goto L1;
	IOSC=CEINAC_.IEIN[IO-1][M-1];
    UK=CEIN_.UI[IOSC-1][M-1];
    WK=CEIN_.WRI[IOSC-1][M-1];
	
	if (UK > 1.0e-3){
		WTHR=fmax(WCCM,UK);
	}else{
		WTHR=fmax(WCCM,WK);
	}

	if (E < WTHR+1.0e-6){
		DE=0.0e0;
        EP=E;
        CDT=1.0e0;
        ES=0.0e0;
        CDTS=0.0e0;
        IOSC=NO;
        return;
	}

	/*
	Truque: A energia de ressonância e a energia de recuo de corte de
	conchas internas são variadas para produzir um limiar suave.
	*/

	LDIST=true;
	if (UK > 1.0e-3){
		WM=3.0e0*WK-2.0e0*UK;
		if (E > WM){
			WKP=WK;
            QKP=UK;
		}else{
			WKP=(E+2.0e0*UK)/3.0e0;
            QKP=UK*(E/WM);
            WM=E;
		}
		if (WCCM > WM)
		    LDIST=false;

		EE=E+UK;
        WCMAX=0.5e0*EE;
        WDMAX=fmin(WM,WCMAX);
        if (WTHR > WDMAX) 
			LDIST=false;
	}else{
		if (WCCM > WK)
			LDIST=false;
		WKP=WK;
        QKP=WK;
        WM=E;
        EE=E;
        WCMAX=0.5e0*EE;
        WDMAX=WKP+1.0e0;
	}

	//Constantes

	RB=E+TREV;
    GAM=1.0e0+E*RREV;
    GAM2=GAM*GAM;
    BETA2=(GAM2-1.0e0)/GAM2;
    AMOL=pow(((GAM-1.0e0)/GAM),2);
    CPS=E*RB;
    CP=sqrt(CPS);

	//Seções transversais parciais do oscilador ativo.
	//Excitalçoes Distantes
	if (LDIST){
		CPPS=(E-WKP)*(E-WKP+TREV);
        CPP=sqrt(CPPS);
		if (WKP > 1.0e-6*E){
			QM=sqrt(pow((CP-CPP),2)+REV*REV)-REV;
		}else{
			QM=pow(WKP,2)/(BETA2*TREV);
            QM=QM*(1.0e0-QM*RTREV);
		}
		if (QM < QKP){
			RWKP=1.0e0/WKP;
         	XHDL=log(QKP*(QM+TREV)/(QM*(QKP+TREV)))*RWKP;
          	XHDT=fmax(log(GAM2)-BETA2-DELTA,0.0e0)*RWKP;
			if (UK > 1.0e-3){
				F0=(WDMAX-WTHR)*(WM+WM-WDMAX-WTHR)/pow((WM-UK),2);
             	XHDL=F0*XHDL;
             	XHDT=F0*XHDT;
			}
		}else{
			XHDL=0.0e0;
            XHDT=0.0e0;
		}
	}else{
		QM=0.0e0;    //Definido para evitar avisos de compilação.
        CPP=0.0e0;   
        CPPS=0.0e0;  
        XHDL=0.0e0;
        XHDT=0.0e0;
	}

	//Colisoes Fechadas
	RCL=WTHR/EE;

	if (RCL < 0.5e0){
		RL1=1.0e0-RCL;
        RRL1=1.0e0/RL1;
        XHC=(AMOL*(0.5e0-RCL)+1.0e0/RCL-RRL1+(1.0e0-AMOL)*log(RCL*RRL1))/EE;
	}else{
		XHC=0.0e0;
	}

	XHTOT=XHC+XHDL+XHDT;
	if (XHTOT < 1.0e-35){
		DE=0.0e0;
        EP=E;
        CDT=1.0e0;
        ES=0.0e0;
        CDTS=0.0e0;
        IOSC=NO;
        return;
	}

	//Amostragem de variáveis ​​de estado final.

	TST=rand2_(3.0e0)*XHTOT;

	//Colisão fechada dura

	TS1=XHC;
	if (TST < TS1){
		A=5.0e0*AMOL;
        ARCL=A*0.5e0*RCL;
L2:;
		FB=(1.0e0+ARCL)*rand2_(4.0e0);
		if (FB < 1.0e0){
			RK=RCL/(1.0e0-FB*(1.0e0-(RCL+RCL)));
		}else{
			RK=RCL+(FB-1.0e0)*(0.5e0-RCL)/ARCL;
		}
		RK2=RK*RK;
        RKF=RK/(1.0e0-RK);
        PHI=1.0e0+pow(RKF,2)-RKF+AMOL*(RK2+RKF);
		if (rand2_(5.0e0)*(1.0e0+A*RK2) > PHI)
			goto L2;
		//Energia e ângulo de espalhamento (elétron primário).
		DE=RK*EE;
        EP=E-DE;
        CDT=sqrt(EP*RB/(E*(RB-DE)));
		//Energia e ângulo de emissão do raio delta.
		if (CEIN_.KS[IOSC-1][M-1] < 17){
			if (UK > CECUTR_.ECUTR[M-1]){
				ES=DE-UK;
			}else{
				ES=DE;
			}
		}else{
			ES=DE;
		}
		CDTS=sqrt(DE*RB/(E*(DE+TREV)));
		return;
	}

	//Interação longitudinal dura distante.
	TS1=TS1+XHDL;
	if (UK > 1.0e-3){
		DE=WM-sqrt(pow((WM-WTHR),2)-rand2_(7.0e0)*(WDMAX-WTHR)*(WM+WM-WDMAX-WTHR));
	}else{
		DE=WKP;
	}
	EP=E-DE;
	if (TST < TS1){
		QS=QM/(1.0e0+QM*RTREV);
        Q=QS/(pow(((QS/QKP)*(1.0e0+QKP*RTREV)), rand2_(6.0e0))-(QS*RTREV));
        QTREV=Q*(Q+TREV);
        CDT=(CPPS+CPS-QTREV)/(2.0e0*CP*CPP);
		if (CDT > 1.0e0){
			CDT=1.0e0;
		}
		//Energia e ângulo de emissão do raio delta.
		if (CEIN_.KS[IOSC-1][M-1] < 17){
			ES=DE-UK; //Apenas conchas internas.
		}else{
			ES=DE;
		}
		CDTS=0.5e0*(WKP*(E+RB-WKP)+QTREV)/sqrt(CPS*QTREV);
		if (CDTS > 1.0e0)
		    CDTS=1.0e0;
		return;
	}

	//Interação transversal distante difícil.
	CDT=1.0e0;
	//Energia e ângulo de emissão do raio delta.
	if (CEIN_.KS[IOSC-1][M-1] < 17){
		if (UK >  CECUTR_.ECUTR[M-1]){
			ES=DE-UK; //Apenas conchas internas.
		}else{
			ES=DE;
		}
	}else{
		ES=DE;
	}
	CDTS=1.0e0;

}

void ebra2_(double &E, double &W, int &M){

 /*
	Simulação da emissão de bremsstrahlung por elétrons ou pósitrons em
	material M.
 */

	/*if (imprimiu==0){
		printf("\n\nEBRA2\n\n");
		imprimiu++;
	}*/


	static const int NBW=32;
	int IE, I, J, K;
	double PT, W1, W2, DW, B, A, PMAX;


	if (PENELOPE_mod_.WCR[M-1] > E){
		W=0.0e0;
        return;
	}

	//Seleção do ponto da rede de energia.

	if (rand2_(1.0e0) < *CEGRID_.XEK){
		IE=*CEGRID_.KE+1;
	}else{
		IE=*CEGRID_.KE;
	}
	//Ponteiro
L1:;
	PT=CEBR_.PBCUT[IE-1][M-1]+rand2_(2.0e0)*(CEBR_.PACB[NBW-1][IE-1][M-1]-CEBR_.PBCUT[IE-1][M-1]);

	//Pesquisa binária do intervalo W.
	I=1;
    J=NBW;
L2:;
    K=(I+J)/2;
	if (PT > CEBR_.PACB[K-1][IE-1][M-1]){
		I=K;
	}else{
		J=K;
	}

	if ((J-I) > 1){
	//	printf("loop ebra2 gotoL2\n");
	//	printf("PT= %f, PACB= %f, I= %d, J= %d, K= %d\n\n", PT,  CEBR_.PACB[K-1][IE-1][M-1], I, J, K );
		goto L2;

	}

	//Amostragem da energia do fóton (método de rejeição).
	W1=CEBR_.WB[I-1];
    W2=CEBR_.WB[I+1-1];
    DW=W2-W1;
    B=CEBR_.DPDFB[I-1][IE-1][M-1]/DW;
    A=CEBR_.PDFB[I-1][IE-1][M-1]-B*W1;
	if (W1 < CEBR_.WBCUT[IE-1][M-1]){
		W1=CEBR_.WBCUT[IE-1][M-1];
	}
	if (W2 < W1){
		printf(" **** WARNING: EBR. Conflicting end-point values.\n");
		W=W1;
		return;
	}
	PMAX=fmax(A+B*W1,A+B*W2);
L3:;
	W=W1*pow((W2/W1),rand2_(3.0e0));
	if ((rand2_(4.0e0)*PMAX) > (A+B*W)){
	//	printf("loop ebra2 gotoL3\n");
		goto L3;
	}
	W=W*E;
	if (W < PENELOPE_mod_.WCR[M-1]){
	//	printf("loop ebra2 gotoL1\n");
		goto L1;
		

	}

}

void ebraa2_(double &E, double &DE, double &CDT, int &M){
	/*
		Amostragem aleatória da direção inicial dos fótons bremss, relativa
	na direção do projétil.
	Ajuste numérico/interpolação de funções de forma de onda parcial dadas por
	Kissel, Quarles e Pratt; ANDT 28(1993)381.

	Parâmetros de entrada:
	M ..... material onde o projétil se move.
	E ..... energia cinética do projétil.
	DE .... energia do fóton emitido.
	Parâmetro de saída:
	CDT ... cosseno do ângulo de emissão polar.
	*/

	/*if (imprimiu==0){
		printf("\n\nEBRAA2\n\n");
		imprimiu++;
	}*/

	double REV=5.10998928e5;
	double TREV=2.0e0*REV;

	double BETA, RK, P10, P11, P1, P20, P21, P2, BETAP;
	int IE, IE1, IET, IK;

	//Parametros de distribuição

	BETA=sqrt(E*(E+TREV))/(E+REV);

	//Uma distribuição dipolar pura é usada para E>500 keV.

	if (E > 500.0e3){
		CDT=2.0e0*rand2_(1.0e0)-1.0e0;
		if (rand2_(2.0e0) > 0.75e0){
			if (CDT < 0.0e0){
				 CDT=-pow((-CDT),0.333333333333333e0);
			}else{
				CDT=pow(CDT,0.333333333333333e0);
			}
		}
		CDT=(CDT+BETA)/(1.0e0+BETA*CDT);
		return;
	}

	if (BETA > CBRANG_.BET[6-1]){
		IE=6;
		goto L20;
	}
	if (BETA < CBRANG_.BET[1-1]){
		IE=1;
		goto L20;
	}
	IE=1;
    IE1=6;
L10:;
    IET=(IE+IE1)/2;
	if (BETA > CBRANG_.BET[IET-1]){
		IE=IET;
	}else{
		IE1=IET;
	}
	if ((IE1-IE) > 1)
		goto L10;
L20:;

	RK=1.0e0+20.0e0*DE/E;
    IK=min(int(RK),20);

	P10=CBRANG_.BP1[1-1][IK-1][IE-1][M-1]+BETA*(CBRANG_.BP1[2-1][IK-1][IE-1][M-1]
     	+BETA*(CBRANG_.BP1[3-1][IK-1][IE-1][M-1]+BETA*CBRANG_.BP1[4-1][IK-1][IE-1][M-1]));
    P11=CBRANG_.BP1[1-1][IK+1-1][IE-1][M-1]+BETA*(CBRANG_.BP1[2-1][IK+1-1][IE-1][M-1]
    	+BETA*(CBRANG_.BP1[3-1][IK+1-1][IE-1][M-1]+BETA*CBRANG_.BP1[4-1][IK+1-1][IE-1][M-1]));
    P1=P10+(RK-IK)*(P11-P10);

    P20=CBRANG_.BP2[1-1][IK-1][IE-1][M-1]+BETA*(CBRANG_.BP2[2-1][IK-1][IE-1][M-1]
		+BETA*(CBRANG_.BP2[3-1][IK-1][IE-1][M-1]+BETA*CBRANG_.BP2[4-1][IK-1][IE-1][M-1]));
    P21=CBRANG_.BP2[1-1][IK+1-1][IE-1][M-1]+BETA*(CBRANG_.BP2[2-1][IK+1-1][IE-1][M-1]
    	+BETA*(CBRANG_.BP2[3-1][IK+1-1][IE-1][M-1]+BETA*CBRANG_.BP2[4-1][IK+1-1][IE-1][M-1]));
    P2=P20+(RK-IK)*(P21-P20);

	//Amostragem das distribuições dipolo transformadas por Lorentz.

	P1=fmin(exp(P1)/BETA,1.0e0);
    BETAP=fmin(fmax(BETA*(1.0e0+P2/BETA),0.0e0),0.999999999e0);

	if (rand2_(3.0e0) < P1){
L1:;
		CDT=2.0e0*rand2_(4.0e0)-1.0e0;
		if ((2.0e0*rand2_(5.0e0)) > (1.0e0+CDT*CDT))
			goto L1;

	}else{
L2:;
		CDT=2.0e0*rand2_(4.0e0)-1.0e0;
		if (rand2_(5.0e0) > 1.0e0-CDT*CDT)
			goto L2;
	}
	CDT=(CDT+BETAP)/(1.0e0+BETAP*CDT);

}

void esia2_(double &E, double &DELTA, double &DE, double &EP, double &CDT, double &ES, double &CDTS, int &M, int &IZZ, int &ISH){

/*
	Amostragem aleatória de ionização de camada interna por impacto de elétrons.

	Modelo Sternheimer-Liljequist GOS.

	Argumentos de entrada:
	E ....... energia do elétron (eV).
	M ....... material onde os elétrons se propagam.
	DELTA ... Correção do efeito de densidade de Fermi.
	Argumentos de saída:
	DE ...... perda de energia (eV).
	EP ...... energia do elétron espalhado (eV).
	CDT ..... cosseno do ângulo de dispersão polar.
	ES ...... energia do elétron secundário emitido (eV).
	CDTS .... cosseno polar de direção do elétron secundário.
	IZZ ..... número atômico do elemento onde a ionização ocorreu.
	ISH ..... camada de elétrons atômicos que foi ionizada.
*/

	/*if (imprimiu==0){
		printf("\n\nESIA2\n\n");
		imprimiu++;
	}*/

	double REV=5.10998928e5;
	double RREV=1.0e0/REV;

	double TREV=2.0e0*REV;
	double RTREV=1.0e0/TREV;

	static const int NO=512;

	double WCCM, PK, TST, UK, WK,WTHR, WM, WKP, QKP, EE, WCMAX, WDMAX, RB, GAM, GAM2, BETA2, AMOL, CPS, CP;
	double CPPS, CPP, QM, RWKP, XHDL, XHDT, F0, RCL, RL1, RRL1, XHC, XHTOT, TS1, A, ARCL, FB, RK, RK2, RKF, PHI;
	double QS, Q, QTREV;
	int IO, JO, IT, IOSC, JE;

	//Ponto da rede de energia

	PK=(*CEGRID_.XEL-CEGRID_.DLEMP[*CEGRID_.KE-1])* *CEGRID_.DLFC;
	if (rand2_(1.0e0) < PK){
		JE=*CEGRID_.KE+1;
	}else{
		JE=*CEGRID_.KE;
	}

	//Seleção do oscilador ativo.

	TST=rand2_(2.0e0);

	//Busca binaria
	IO=1;
    JO=CESIAC_.NESI[M-1]+1;
L1:;
    IT=(IO+JO)/2;
	if (TST > CESIAC_.ESIAC[IT-1][JE-1][M-1]){
		IO=IT;
	}else{
		JO=IT;
	}
	if ((JO-IO) > 1)
		goto L1;
	
	IOSC=CESIAC_.IESI[IO-1][M-1];
    IZZ=CEIN_.KZ[IOSC-1][M-1];
    ISH=CEIN_.KS[IOSC-1][M-1];
    UK=CEIN_.UI[IOSC-1][M-1];
    WK=CEIN_.WRI[IOSC-1][M-1];

	if (UK > 1.0e-3){
		WTHR=UK;
	}else{
		WTHR=WK;
	}

	if (E < (WTHR+1.0e-6)){
		DE=UK;
        EP=E-DE;
        CDT=1.0e0;
        ES=0.0e0;
        CDTS=0.0e0;
        return;
	}

	/*
	Truque: A energia de ressonância e a energia de recuo de corte de
	 conchas internas são variadas para produzir um limiar suave.
	*/

 	WM=3.0e0*WK-2.0e0*UK;
	if (E > WM){
		WKP=WK;
        QKP=UK;
	}else{
		WKP=(E+2.0e0*UK)/3.0e0;
        QKP=UK*(E/WM);
        WM=E;
	}
	EE=E+UK;
    WCMAX=0.5e0*EE;
    WDMAX=fmin(WM,WCMAX);


	//Constantes

	RB=E+TREV;
    GAM=1.0e0+E*RREV;
    GAM2=GAM*GAM;
    BETA2=(GAM2-1.0e0)/GAM2;
    AMOL=pow(((GAM-1.0e0)/GAM),2);
    CPS=E*RB;
    CP=sqrt(CPS);

	//Seções transversais parciais do oscilador ativo.

	//Excitações distantes.

	CPPS=(E-WKP)*(E-WKP+TREV);
    CPP=sqrt(CPPS);

	if (WKP > (1.0e-6*E)){
		QM=sqrt(pow((CP-CPP),2)+REV*REV)-REV;
	}else{
		QM=pow(WKP,2)/(BETA2*TREV);
        QM=QM*(1.0e0-QM*RTREV);
	}

	if (QM < QKP){
		RWKP=1.0e0/WKP;
        XHDL=log(QKP*(QM+TREV)/(QM*(QKP+TREV)))*RWKP;
        XHDT=fmax(log(GAM2)-BETA2-DELTA,0.0e0)*RWKP;
        F0=(WDMAX-WTHR)*(WM+WM-WDMAX-WTHR)/pow((WM-UK),2);
        XHDL=F0*XHDL;
        XHDT=F0*XHDT;
	}else{
		XHDL=0.0e0;
        XHDT=0.0e0;
	}

	//Colisoes Proximas
	RCL=WTHR/EE;
    RL1=1.0e0-RCL;
    RRL1=1.0e0/RL1;
    XHC=(AMOL*(0.5e0-RCL)+1.0e0/RCL-RRL1+(1.0e0-AMOL)*log(RCL*RRL1))/EE;

    XHTOT=XHC+XHDL+XHDT;

	if (XHTOT < 1.0e-35){
		DE=UK;
        EP=E-DE;
        CDT=1.0e0;
        ES=0.0e0;
        CDTS=0.0e0;
        return;
	}

	//Amostragem de variáveis ​​de estado final.

	TST=rand2_(3.0e0)*XHTOT;

	//Colisão aproximada dura.

	TS1=XHC;
	if (TST < TS1){
		A=5.0e0*AMOL;
        ARCL=A*0.5e0*RCL;
L2:;
		FB=(1.0e0+ARCL)*rand2_(4.0e0);
		if (FB < 1.0e0){
			RK=RCL/(1.0e0-FB*(1.0e0-(RCL+RCL)));
		}else{
			RK=RCL+(FB-1.0e0)*(0.5e0-RCL)/ARCL;
		}
		RK2=RK*RK;
		RKF=RK/(1.0e0-RK);
		PHI=1.0e0+pow(RKF,2)-RKF+AMOL*(RK2+RKF);
		if (rand2_(5.0e0)*(1.0e0+A*RK2) > PHI)
			goto L2;
		//Energia e ângulo de espalhamento (elétron primário).
		DE=RK*EE;
		EP=E-DE;
		CDT=sqrt(EP*RB/(E*(RB-DE)));
		//Energia e ângulo de emissão do raio delta.
		ES=DE-UK;
		CDTS=sqrt(DE*RB/(E*(DE+TREV)));
		return;
	}

	//Interação longitudinal dura distante.
	TS1=TS1+XHDL;
    DE=WM-sqrt(pow((WM-WTHR),2)-rand2_(7.0e0)*(WDMAX-WTHR)*(WM+WM-WDMAX-WTHR));
    EP=E-DE;

	if (TST < TS1){
		QS=QM/(1.0e0+QM*RTREV);
        Q=QS/(pow(((QS/QKP)*(1.0e0+QKP*RTREV)),rand2_(6.0e0))-(QS*RTREV));
        QTREV=Q*(Q+TREV);
        CDT=(CPPS+CPS-QTREV)/(2.0e0*CP*CPP);
		if (CDT > 1.0e0)
			CDT=1.0e0;

		//Energia e ângulo de emissão do raio delta.
		ES=DE-UK;
        CDTS=0.5e0*(WKP*(E+RB-WKP)+QTREV)/sqrt(CPS*QTREV);
		if (CDTS >1.0e0)
			CDTS=1.0e0;
		return;
	}

	//Interação transversal distante difícil.

	CDT=1.0e0;

	//Energia e ângulo de emissão do raio delta.

	ES=DE-UK;
    CDTS=1.0e0;

}

void relax2_(int &IZ, int &IS){
    /*
	Esta sub-rotina simula o relaxamento de um átomo ionizado
	o elemento IZ com uma vaga no shell IS (o shell K ou um L
	subcamada C ou M). Essa vacância inicial é preenchida por elétrons de
	camadas externas através de transições radiativas e não radiativas, que
	pode produzir vagas adicionais.

	Usamos a seguinte notação para designar as transições possíveis:
	* Radiativo: IS0-IS1 (um elétron da camada IS1 preenche o
	vacância na casca IS0, deixando um buraco na casca IS1).
	* Não radiativo: IS0-IS1-IS2 (um elétron da camada IS1 preenche
	a vacância na camada IS0, e a energia liberada é tomada
	afastado por um elétron na camada IS2; este processo deixa dois
	, nas carcaças IS1 e IS2).
	A cascata de desexcitação (ou seja, o conjunto de transições que ocorrem para
	uma determinada vaga inicial) é amostrado a partir das probabilidades de transição
	contido na Biblioteca de Dados Atômicos Avaliados de Livermore (EADL). O
	A energia C da radiação emitida em cada transição é lida do
	Base de dados C PENELOPE.

    A simulação da cascata de desexcitação é descontinuada
	quando as conchas K, L e M estiverem cheias ou quando não houver
	energia suficiente para produzir radiação 'ativa' (com energia maior que
	EAB). A energia de excitação do íon residual é assumida como sendo
	depositado localmente. Desconsideramos a emissão e transporte de soft
	raios-x e elétrons lentos, cujas energias são menores do que a ligação
	energia C da camada N1 do elemento mais pesado no meio. este
	estabelece um limite inferior para o intervalo de energia que pode ser coberto pelo
	programa de simulação C de forma consistente.

	Os dados de desexcitação para os elementos carregados são armazenados no
	Bloco /CRELAX/, em uma forma projetada para minimizar a quantidade de memória
	e para facilitar a amostragem aleatória. As quantidades em comum
	bloco são os seguintes:
	IFIRST(99,16) ... dados de desexcitação para uma vaga no shell IS de
	o elemento IZ começa na posição K=IFIRST(IZ,IS) no
	matrizes de armazenamento. Os valores permitidos para IS são de 1 a 16 (K shell
	subcamadas L, M e N).
	ILAST(99,16) ... os dados de desexcitação para uma vaga no shell
	IS do elemento IZ termina na posição K=ILAST(IZ,IS) no
	matrizes de armazenamento.
	IS1(K), IS2(K) ... shells que estão ativos na transição (veja o
	código do rótulo do shell abaixo). Para transições radiativas, IS2(K)=0.
	P(K) ... probabilidade relativa para a transição IS-IS1(K)-IS2(K).
	ET(K) ... energia da partícula secundária emitida na transição.
	F(K), IAL(K) ... valores de corte e alias (método de amostragem de Walker).

  ---------------------------------------------------------------------
  Label code IS for electron shells:
      1 = K  (1s1/2),     11 = N2 (4p1/2),     21 = O5 (5d5/2),
      2 = L1 (2s1/2),     12 = N3 (4p3/2),     22 = O6 (5f5/2),
      3 = L2 (2p1/2),     13 = N4 (4d3/2),     23 = O7 (5f7/2),
      4 = L3 (2p3/2),     14 = N5 (4d5/2),     24 = P1 (6s1/2),
      5 = M1 (3s1/2),     15 = N6 (4f5/2),     25 = P2 (6p1/2),
      6 = M2 (3p1/2),     16 = N7 (4f7/2),     26 = P3 (6p3/2),
      7 = M3 (3p3/2),     17 = O1 (5s1/2),     27 = P4 (6d3/2),
      8 = M4 (3d3/2),     18 = O2 (5p1/2),     28 = P5 (6d5/2),
      9 = M5 (3d5/2),     19 = O3 (5p3/2),     29 = Q1 (7s1/2),
     10 = N1 (4s1/2),     20 = O4 (5d3/2),     30 = outer shells.
  ---------------------------------------------------------------------
    */

   /*	if (imprimiu==0){
		printf("\n\nRELAX2\n\n");
		imprimiu++;
	}*/

   double PI=3.1415926535897932e0;
   double TWOPI=PI+PI;

   double PTIM[256];
   double ISV[256];

   static const int NRX=60000;

   int NV, ISP, KF, KL, K1, NVIS, IS1K, IS2K, KPARS;
   double PAGE0, RN, TST, WS, US, VS, DF, SDTS;

   //Inicializacao

    if ((IZ < 3) || (IS > 16)){
	   return; 
    }

    if (CADATA_.EB[IS-1][IZ-1] < CECUTR_.ECUTR[*TRACK_mod_.MAT-1]) ////Se a energia de ionização da casca for menor que ECUTR, a cascata não é seguido.
		return;

	NV=1;
    ISV[1-1]=IS;
    PAGE0=*TRACK_mod_.PAGE;
    PTIM[1-1]=*TRACK_mod_.PAGE;

    //Proxima transição

L1:;
	ISP=ISV[NV-1];
    *TRACK_mod_.PAGE=PTIM[NV-1];
    KF=CRELAX_.IFIRST[ISP-1][IZ-1];
    KL=CRELAX_.ILAST[ISP-1][IZ-1];
    NV=NV-1;

	if (KL > KF){
		//Algoritmo de amostragem de Walker.
		RN=rand2_(1.0e0)*(KL-KF+1);
        K1=int(RN);
        TST=RN-K1;
		if (TST > CRELAX_.F[KF+K1-1]){
			*CRELAX_.KS=CRELAX_.IS0[KF+K1-1];
		}else{
			*CRELAX_.KS=KF+K1;
		}
	}else{
		*CRELAX_.KS=KF;
	}

	/*
	Se MODER=0, o controle é retornado ao programa chamador após
	determinando a primeira transição, KS. Útil para testar o aleatório
	amostragem. Para operação normal, podemos comentar o seguinte
	declaração .
	*/

	if (*CRELAX_.MODER == 0)
		return;
		//Se LAGE=true, a idade da partícula é registrada.
	if (*TRACK_mod_.LAGE){
		NVIS=1;
		if (NV > 1) {//Várias vagas no shell ativo?
			for (int ISC = 1; ISC <= NV; ISC++){
				if (ISV[ISC-1] == ISP)
				    NVIS=NVIS+1;
			}
		}
		*TRACK_mod_.PAGE=*TRACK_mod_.PAGE-(CADATA_.ALW[ISP-1][IZ-1]/(NVIS))*log(rand2_(2.0e0));
	}

	//radiação fluorescente

	IS1K=CRELAX_.IS1[*CRELAX_.KS-1];
    IS2K=CRELAX_.IS2[*CRELAX_.KS-1];
	if (IS2K == 0){
		KPARS=2;
		if (IS1K < 17){
			if (CADATA_.EB[IS1K-1][IZ-1] > CECUTR_.ECUTR[*TRACK_mod_.MAT-1]){
				NV=NV+1;
            	ISV[NV-1]=IS1K;
            	PTIM[NV-1]=*TRACK_mod_.PAGE;
			}
		}
	}else{
		KPARS=1;
		if (IS1K < 17){
			if (CADATA_.EB[IS1K-1][IZ-1] > CECUTR_.ECUTR[*TRACK_mod_.MAT-1]){
				NV=NV+1;
            	ISV[NV-1]=IS1K;
            	PTIM[NV-1]=*TRACK_mod_.PAGE;
			}
		}
		if (IS2K < 17){
			if (CADATA_.EB[IS2K-1][IZ-1] > CECUTR_.ECUTR[*TRACK_mod_.MAT-1]){
				NV=NV+1;
            	ISV[NV-1]=IS2K;
            	PTIM[NV-1]=*TRACK_mod_.PAGE;
			}
		}
	}

	/*	
	A partícula emitida é armazenada na pilha secundária quando
	sua energia ET(K) é maior que EABS.
	*/

	if (CRELAX_.ET[*CRELAX_.KS-1] > PENELOPE_mod_.EABS[*TRACK_mod_.MAT-1][KPARS-1]){
		//Direção inicial (isotrópica).
		WS=-1.0e0+2.0e0*rand2_(2.0e0);
        SDTS=sqrt(1.0e0-WS*WS);
        DF=TWOPI*rand2_(3.0e0);
        US=cos(DF)*SDTS;
        VS=sin(DF)*SDTS;
        CHIST_.ILBA[1-1]=TRACK_mod_.ILB[1-1]+1;
        CHIST_.ILBA[2-1]=*TRACK_mod_.KPAR;
        CHIST_.ILBA[4-1]=IZ*1000000+ISP*10000+IS1K*100+IS2K;
        CHIST_.ILBA[5-1]=TRACK_mod_.ILB[5-1];
        stores2_(CRELAX_.ET[*CRELAX_.KS-1],*TRACK_mod_.X,*TRACK_mod_.Y,*TRACK_mod_.Z,US,VS,WS,*TRACK_mod_.WGHT,KPARS,CHIST_.ILBA,wIPOLI);
	}

	//Existem vagas não preenchidas nas conchas internas?

	if (NV > 0)
		goto L1;
	*TRACK_mod_.PAGE=PAGE0;


}

void eaux2_(){
	/*	
	Mecanismo de interação auxiliar para elétrons, definível pelo usuário.
	Geralmente não está ativo.
	*/
	printf("Warning: Subroutine EAUX has been entered.\n");
}

void graa2_(double &E, double &CDT, int &IEFF, int &M){

	/*if (imprimiu==0){
		printf("\n\nGRAA2\n\n");
		imprimiu++;
	}*/

	//Amostragem aleatória de espalhamento coerente (Rayleigh)

	double REV=5.10998928e5;
	double RREV=1.0e0/REV;
	static const int NQ=250;
	static const int NEX=1024;
	static const int NP=150;
	static const int NPM1=NP-1;

	int II, IU, IT, ITN, I, J, K;
	double XSE, QMAX, Q2MAX, G, RU, RR, D, XX;

	//Busca Binaria

	II=CGRA01_.IED[*CEGRID_.KE-1];
    IU=CGRA01_.IEU[*CEGRID_.KE-1];
L1:;
    IT=(II+IU)/2;

	if (*CEGRID_.XEL > CGRA01_.ERA[IT-1])
		II=IT;
	else
		IU=IT;
	if (IU-II > 1)
		goto L1;
	
	XSE=exp(CGRA01_.XSRA[II-1][M-1]+(CGRA01_.XSRA[II+1-1][M-1]-CGRA01_.XSRA[II-1][M-1])*(*CEGRID_.XEL-CGRA01_.ERA[II-1])/(CGRA01_.ERA[II+1-1]-CGRA01_.ERA[II-1]));
	
	if (rand2_(1.0e0)*CGIMFP_.SGRA[*CEGRID_.KE-1][M-1] > XSE){
		IEFF=0;
    	CDT=1.0e0;
		return;
	}

	IEFF=1;
    QMAX=2.0e0*E*RREV;

	if (QMAX < 1.0e-10){
L2:;
		CDT=1.0e0-2.0e0*rand2_(1.0e0);
        G=0.5e0*(1.0e0+CDT*CDT);	
		if (rand2_(2.0e0) > G)
			goto L2;
		return;
	}
	Q2MAX=fmin(QMAX*QMAX,CGRA03_.QRA[M-1][NP-1]);

L3:;
	RU=rand2_(3.0e0)*CGRA03_.PMAX[M-1][*CEGRID_.KE+1-1];

	//Seleção do intervalo (busca binária dentro de limites pré-calculados).

	ITN=RU*NPM1+1;
    I=CGRA03_.ITLRA[M-1][ITN-1];
    J=CGRA03_.ITURA[M-1][ITN-1];

	if ((J-I) < 2)
		goto L5;
L4:;
	K=(I+J)/2;
	if (RU > CGRA03_.PRA[M-1][K-1])
		I=K;
	else
		J=K;
	if ((J-I) > 1)
		goto L4;
	
	//Amostragem da distribuição cumulativa inversa racional.

L5:;
	RR=RU-CGRA03_.PRA[M-1][I-1];
	if (RR > 1.0e-16){
		D=CGRA03_.DPRA[M-1][I-1];
        XX=CGRA03_.QRA[M-1][I-1]+((1.0e0+CGRA03_.ARA[M-1][I-1]+CGRA03_.BRA[M-1][I-1])*D*RR/(D*D+(CGRA03_.ARA[M-1][I-1]*D+CGRA03_.BRA[M-1][I-1]*RR)*RR))
			*(CGRA03_.QRA[M-1][I+1-1]-CGRA03_.QRA[M-1][I-1]);
	}else{
		XX=CGRA03_.QRA[M-1][I-1];
	}
	if (XX > Q2MAX)
		goto L3;
	CDT=1.0e0-2.0e0*XX/Q2MAX;
	//Rejeição
	G=0.5e0*(1.0e0+CDT*CDT);
	if (rand2_(4.0e0) > G)
		goto L3;



}

void dirpol2_(double &CDT, double &DF, double &CONS, double &SP1, double &SP2, double &SP3, double &U, double &V, double &W){

	/*
	Esta sub-rotina calcula os cossenos de direção _e_ os Stokes
	Parâmetros de um fóton polarizado após espalhamento com um dado polar
	ângulo.

	Entrada: U,V,W ... cossenos da direção inicial.
	SP1,SP2,SP3 ... parâmetros Stokes iniciais.
	CDT ..... cosseno do ângulo de dispersão polar.
	CONS .... constante no PDF do ângulo azimutal.
	Saída: U,V,W ... novos cossenos de direção.
	SP1,SP2,SP3 ... novos parâmetros Stokes.
	DF ...... ângulo de dispersão azimutal.
	CDT e CONS permanecem inalterados.
	*/

	/*if (imprimiu==0){
		printf("\n\nDIRPOL2\n\n");
		imprimiu++;
	}*/

	double PI=3.1415926535897932e0;
	double TWOPI=2.0e0*PI;

	double SP1P, RSP0, CDT2, CDT21, PHA, PHB, SP0MAX, SDF, CDF, S2DF, C2DF, SP3P, SP0P, UV, UVW, FNORM, SDT, SDTSDF, SDTCDF,SUV, UN, VN;

	//Amostragem do ângulo de espalhamento azimutal.

	CDT2=CDT*CDT;
    CDT21=CDT2+1.0e0;
    PHA=CDT21+CONS;
    PHB=1.0e0-CDT2;
    SP0MAX=PHA+PHB*sqrt(SP1*SP1+SP3*SP3+1.0e-35);
L1:;
	DF=rand2_(1.0e0)*TWOPI;
    SDF=sin(DF);
    CDF=cos(DF);
    S2DF=2.0e0*SDF*CDF;
    C2DF=CDF*CDF-SDF*SDF;
    SP3P=S2DF*SP1+C2DF*SP3; //Parâmetro Stokes com novo zero azimute.
    SP0P=PHA-PHB*SP3P;

	if (rand2_(2.0e0)*SP0MAX > SP0P)
		goto L1;

	//Calcular novos parâmetros Stokes
 	SP1P=C2DF*SP1-S2DF*SP3; //Parâmetro Stokes com novo zero azimute.
    RSP0=1.0e0/SP0P;
    SP1=2.0e0*CDT*SP1P*RSP0;
    SP2=(2.0e0+CONS)*CDT*SP2*RSP0;
    SP3=(CDT21*SP3P-PHB)*RSP0;

    //Garanta a normalizacao
	UV=U*U+V*V;
    UVW=UV+W*W;
	if (fabs(UVW-1.0e0) > 1.0e-13){
		FNORM=1.0e0/sqrt(UVW);
        U=FNORM*U;
        V=FNORM*V;
        W=FNORM*W;
        UV=U*U+V*V;
	}

	//Calcula a nova direção
	if (1.0e0-fabs(CDT) > 1.0e-8)
		SDT=sqrt(PHB);
	else
		SDT=sqrt(2.0e0*(1.0e0-fabs(CDT)));

	if (SDT < 1.0e-13){
		if (CDT < 0.0e0){
			U=-U;
          	V=-V;
          	W=-W;
		}
	}else{
		SDTSDF=SDT*SDF;
        SDTCDF=SDT*CDF;
		if (UV > 1.0e-26){
			SUV=sqrt(UV);
          	UN=U/SUV;
          	VN=V/SUV;
         	U=U*CDT+(UN*W*SDTCDF-VN*SDTSDF);
          	V=V*CDT+(VN*W*SDTCDF+UN*SDTSDF);
          	W=W*CDT-SUV*SDTCDF;
		}else{
			if (W > 0.0e0){
				U=SDTCDF;
            	V=SDTSDF;
            	W=CDT;
			}else{
				U=-SDTCDF;
            	V=-SDTSDF;
            	W=-CDT;
			}
		}
	}
}

void gcoa2_(double &E, double &DE, double &EP, double &CDT, double &ES, double &CDTS, int &M, int &IZZ, int &ISH){

 /*
	Amostragem aleatória de espalhamento incoerente (Compton) de fótons.

	Argumentos de entrada:
	E ..... energia do fóton incidente (eV).
	M ..... material onde os fótons se propagam.
	Argumento de saída:
	DE .... perda de energia (eV).
	EP .... energia do fóton espalhado (eV).
	CDT ... cosseno do ângulo de espalhamento polar.
	ES .... energia do elétron emitido (eV).
	CDTS .. cosseno polar de direção do elétron.
	IZZ ... número atômico do átomo onde ocorreu a dispersão.
	ISH ... camada de elétrons atômica que foi ionizada.
 */

	/*if (imprimiu==0){
		printf("\n\nGCOA2\n\n");
		imprimiu++;
	}*/

	double REV=5.10998928e5;
	double RREV=1.0e0/REV;
	double D2=1.4142135623731e0;
	double D1=1.0e0/D2;
	double D12=0.5e0;

	static const int NOCO=512;

	double RN[NOCO];
	double PAC[NOCO];

	double EK, EK2, EKS, EK1, TAUMIN, TAUM2, A1, A2, S0, AUX, PZOMC, RNI, TAU, CDT1, S, TST, JO, A, XQC, AF;
	double FPZ, FPZMAX, T, B1, B2, Q2; 
	int I2, ISHELL, I3;

	EK=E*RREV;
    EK2=EK+EK+1.0e0;
    EKS=EK*EK;
    EK1=EKS-EK2-1.0e0;
    TAUMIN=1.0e0/EK2;
    TAUM2=TAUMIN*TAUMIN;
    A1=log(EK2);
    A2=A1+2.0e0*EK*(1.0e0+EK)*TAUM2;

	if (E > 5.0e6)
		goto L4;

	//Função de dispersão incoerente para theta=PI

	S0=0.0e0;
	for (int I = 1; I <= CGCO_.NOSCCO[M-1]; I++){
		if (CGCO_.UICO[I-1][M-1] < E){
			AUX=E*(E-CGCO_.UICO[I-1][M-1])*2.0e0;
            PZOMC=CGCO_.FJ0[I-1][M-1]*(AUX-REV*CGCO_.UICO[I-1][M-1])
     			/(REV*sqrt(AUX+AUX+pow(CGCO_.UICO[I-1][M-1],2)));
			if (PZOMC > 0.0e0)
				RNI=1.0e0-0.5e0*exp(D12-pow((D1+D2*PZOMC),2));
			else
				RNI=0.5e0*exp(D12-pow((D1-D2*PZOMC),2));
			S0=S0+CGCO_.FCO[I-1][M-1]*RNI;
		}
	}

	//Amostragem de tau.

L1:;
	if (rand2_(1.0e0)*A2 < A1)
		TAU=pow(TAUMIN,rand2_(2.0e0));
	else
		TAU=sqrt(1.0e0+rand2_(3.0e0)*(TAUM2-1.0e0));
	CDT1=(1.0e0-TAU)/(EK*TAU);

	//Função de dispersão incoerente.
	S=0.0e0;
	for (int I = 1;I <= CGCO_.NOSCCO[M-1]; I++ ){
		if (CGCO_.UICO[I-1][M-1] < E){
			AUX=E*(E-CGCO_.UICO[I-1][M-1])*CDT1;
			PZOMC=CGCO_.FJ0[I-1][M-1]*(AUX-REV*CGCO_.UICO[I-1][M-1])
				/(REV*sqrt(AUX+AUX+pow(CGCO_.UICO[I-1][M-1],2)));
			if (PZOMC > 0.0e0)
				RN[I-1]=1.0e0-0.5e0*exp(D12-pow((D1+D2*PZOMC),2));
			else
				RN[I-1]=0.5e0*exp(D12-pow((D1-D2*PZOMC),2));
			S=S+CGCO_.FCO[I-1][M-1]*RN[I-1];
			PAC[I-1]=S;
		}else{
			PAC[I-1]=S;
		}
	}

	//Funcao de Rejeição
	TST=S*(1.0e0+TAU*(EK1+TAU*(EK2+TAU*EKS)))/(EKS*TAU*(1.0e0+TAU*TAU));
	if (rand2_(4.0e0)*S0 > TST)
		goto L1;
	CDT=1.0e0-CDT1;

	//Escudo do elétron alvo.

L2:
	TST=S*rand2_(5.0e0);
	//Busca Binaria
	if (TST < PAC[1-1]){
		ISHELL=1;
	}else{
		ISHELL=1;
        JO=CGCO_.NOSCCO[M-1]+1;
L3:;
    	I2=(ISHELL+JO)/2;
		if (TST > PAC[I2-1])
			ISHELL=I2;
		else
			JO=I2;
		if (JO-ISHELL > 1)
			goto L3;
		ISHELL=ISHELL+1;
	}

	//Momento projetado do elétron alvo.
	A=rand2_(6.0e0)*RN[ISHELL-1];
	if (A < 0.5e0)
	 	 PZOMC=(D1-sqrt(D12-log(A+A)))/(D2*CGCO_.FJ0[ISHELL-1][M-1]);
	else
		PZOMC=(sqrt(D12-log(2.0e0-A-A))-D1)/(D2*CGCO_.FJ0[ISHELL-1][M-1]);
	if (PZOMC < -1.0e0)
		goto L2;

	//Rejeição F(EP).
	XQC=1.0e0+TAU*(TAU-2.0e0*CDT);
    AF=sqrt(XQC)*(1.0e0+TAU*(TAU-CDT)/XQC);

	if (AF > 0.0e0)
		FPZMAX=1.0e0+AF*0.2e0;
	else	
		FPZMAX=1.0e0-AF*0.2e0;

	FPZ=1.0e0+AF*fmax(fmin(PZOMC,0.2e0),-0.2e0);
	if (rand2_(7.0e0)*FPZMAX > FPZ)
		goto L2;

	//Energia do fóton espalhado
	T=pow(PZOMC,2);
    B1=1.0e0-T*TAU*TAU;
    B2=1.0e0-T*TAU*CDT;
	if (PZOMC > 0.0e0)
		EP=E*(TAU/B1)*(B2+sqrt(fabs(B2*B2-B1*(1.0e0-T))));
	else
		EP=E*(TAU/B1)*(B2-sqrt(fabs(B2*B2-B1*(1.0e0-T))));
	goto L6;

	//Sem alargamento Doppler para E maior que 5 MeV.

L4:;
	if (rand2_(8.0e0)*A2 < A1)
		TAU=pow(TAUMIN,rand2_(9.0e0));
	else
		TAU=sqrt(1.0e0+rand2_(10.0e0)*(TAUM2-1.0e0));

	//Funcao de Rejeição
	TST=(1.0e0+TAU*(EK1+TAU*(EK2+TAU*EKS)))/(EKS*TAU*(1.0e0+TAU*TAU));
	if (rand2_(11.0e0) > TST)
		goto L4;
	EP=TAU*E;
    CDT=1.0e0-(1.0e0-TAU)/(EK*TAU);

	//Camada eletrônica alvo.
	TST=rand2_(12.0e0);

	//busca Binaria
	if (TST < CGCO_.PTRSH[1-1][M-1]){
		ISHELL=1;
	}else{
		ISHELL=1;
        JO=CGCO_.NOSCCO[M-1]+1;
L5:;
    	I3=(ISHELL+JO)/2;
		if (TST > CGCO_.PTRSH[I3-1][M-1])
			ISHELL=I3;
		else
			JO=I3;
		if (JO-ISHELL > 1)
			goto L5;
		ISHELL=ISHELL+1;
	}

	if (EP > (E-CGCO_.UICO[ISHELL-1][M-1]))
		goto L4;

L6:;
	DE=E-EP;
	if (CGCO_.KSCO[ISHELL-1][M-1] < 17){
		if (CGCO_.UICO[ISHELL-1][M-1] > CECUTR_.ECUTR[M-1]){
			ES=DE-CGCO_.UICO[ISHELL-1][M-1];
		}else{
			ES=DE;
		}
	}else{
		ES=DE;
	}

	Q2=E*E+EP*(EP-2.0e0*E*CDT);
	if (Q2 > 1.0e-12)
		CDTS=(E-EP*CDT)/sqrt(Q2);
	else
		CDTS=1.0e0;

	IZZ=CGCO_.KZCO[ISHELL-1][M-1];
    ISH=CGCO_.KSCO[ISHELL-1][M-1];


}

void  gpha2_(double &ES, int &IZZ, int &ISH){

 /*
	Simulação de absorção fotoelétrica no material M.

	Argumentos de saída:
	ES .... energia cinética do fotoelétron.
	IZZ ... número atômico do átomo onde ocorreu a absorção.
	ISH ... camada de elétrons atômica que foi ionizada.

	NOTA: O JUMP usa uma seção transversal fotoelétrica ligeiramente maior
	do que seu valor 'verdadeiro'. Para corrigir isso, o fóton pode
	'sobrevive' a um evento fotoelétrico. A sobrevivência do fóton é sinalizada por
	configuração IZZ=0, ISH=0, ES=0.0D0 (a energia E do fóton é mantida
	inalterado.
 */

	/*if (imprimiu==0){
		printf("\n\nGPHA2\n\n");
		imprimiu++;
	}*/

	static const int NTP=12000;
	double ACP[35];
	double IP[35];

	int I, IU, IT, IELAC, J;
	double PTOT, DEE, PCSL, TST, PIS, EBB;

	//Coeficientes de atenuação parcial.

	PTOT=0.0e0;
	for (int IEL=1; IEL <= COMPOS_.NELEM[*TRACK_mod_.MAT-1]; IEL++){
		IZZ=COMPOS_.IZ[IEL-1][*TRACK_mod_.MAT-1];
       //Busca Binaria
        I=CGPH00_.IPHF[IZZ-1];
        IU=CGPH00_.IPHL[IZZ-1];
L1:;
    	IT=(I+IU)/2;
		if (*CEGRID_.XEL > CGPH00_.EPH[IT-1])
			I=IT;
		else
			IU=IT;
		if (IU-I > 1)
			goto L1;

		IP[IEL-1]=I;
        DEE=CGPH00_.EPH[I+1-1]-CGPH00_.EPH[I-1];
		if (DEE > 1.0e-15){
			PCSL=CGPH00_.XPH[1-1][I-1]+(CGPH00_.XPH[1-1][I+1-1]-CGPH00_.XPH[1-1][I-1])*(*CEGRID_.XEL-CGPH00_.EPH[I-1])/DEE;
		}else{
			PCSL=CGPH00_.XPH[1-1][I-1];
		}
		PTOT=PTOT+COMPOS_.STF[IEL-1][*TRACK_mod_.MAT-1]*exp(PCSL);
        ACP[IEL-1]=PTOT;
	}
	if (PTOT*COMPOS_.VMOL[*TRACK_mod_.MAT-1] > CGIMFP_.SGPH[*CEGRID_.KE-1][*TRACK_mod_.MAT-1]){
		printf("WARNING: SGPH is less than the actual mac.\n");
	}

	//Faça uma amostra do elemento ativo.
	TST=rand2_(1.0e0)*CGIMFP_.SGPH[*CEGRID_.KE-1][*TRACK_mod_.MAT-1]/COMPOS_.VMOL[*TRACK_mod_.MAT-1];
	for (int IEL = 1; IEL <= COMPOS_.NELEM[*TRACK_mod_.MAT-1]; IEL++){
		if (ACP[IEL-1] > TST){
			IELAC=IEL;
          	IZZ=COMPOS_.IZ[IEL-1][*TRACK_mod_.MAT-1];
          	goto L2;
		}
	}

	//nteração delta. Introduzido para corrigir o uso de um limite superior do coeficiente de atenuação fotoelétrica.
	IZZ=0;
    ISH=0;
    ES=0.0e0;
    return;

L2:;
	//Seleção do shell ativo.
	I=IP[IELAC-1];
    DEE=CGPH00_.EPH[I+1-1]-CGPH00_.EPH[I-1];
    PIS=0.0e0;

	if (DEE > 1.0e-15){
		PTOT=exp(CGPH00_.XPH[1-1][I-1]+(CGPH00_.XPH[1-1][I+1-1]-CGPH00_.XPH[1-1][I-1])*(*CEGRID_.XEL-CGPH00_.EPH[I-1])/DEE);
        TST=rand2_(2.0e0)*PTOT;
	
		for (int IS=1; IS <= CGPH00_.NPHS[IZZ-1]; IS++){
			J=IS+1;
          	PCSL=CGPH00_.XPH[J-1][I-1]+(CGPH00_.XPH[J-1][I+1-1]-CGPH00_.XPH[J-1][I-1])*(*CEGRID_.XEL-CGPH00_.EPH[I-1])/DEE;
         	PIS=PIS+exp(PCSL);
			 if (PIS > TST){
				 ISH=IS;
				 goto L3;
			 }
		}
	}else{
		PTOT=exp(CGPH00_.XPH[1-1][I-1]);
        TST=rand2_(2.0e0)*PTOT;
		for (int IS=1; IS <= CGPH00_.NPHS[IZZ-1]; IS++){
			PIS=PIS+exp(CGPH00_.XPH[IS+1-1][I-1]);
			if (PIS > TST){
				ISH=IS;
				goto L3;
			}
		}
	}
	ISH=17;

	//Emissão do FotoEletron
L3:;
	if (ISH < 17){
		EBB=CADATA_.EB[ISH-1][IZZ-1];
		if (EBB > CECUTR_.ECUTR[*TRACK_mod_.MAT-1]){
			ES=*TRACK_mod_.E-EBB;
		}else{
			ES=*TRACK_mod_.E;
          	ISH=17;
		}
	}else{
		ES=*TRACK_mod_.E;
	}

}

void sauter2_(double &ES, double &CDTS){

	/*
	Amostragem aleatória da direção inicial dos fotoelétrons do
	Distribuição Sauter.
	*/

	/*if (imprimiu==0){
		printf("\n\nSAUTER2\n\n");
		imprimiu++;
	}*/

	double REV=5.10998928e5;

	double GAM, GAM2, BETA, AC, A1, A2, GTMAX, RU, TSAM, GTR;

	if (ES > 1.0e9){
		CDTS=1.0e0;
		return;
	}
	GAM=1.0e0+ES/REV;
    GAM2=GAM*GAM;
    BETA=sqrt((GAM2-1.0e0)/GAM2);
    AC=1.0e0/BETA-1.0e0;
    A1=0.5e0*BETA*GAM*(GAM-1.0e0)*(GAM-2.0e0);
    A2=AC+2.0e0;
    GTMAX=2.0e0*(A1+1.0e0/AC);
L1:;
	RU=rand2_(1.0e0);
    TSAM=2.0e0*AC*(2.0e0*RU+A2*sqrt(RU))/(A2*A2-4.0e0*RU);
    GTR=(2.0e0-TSAM)*(A1+1.0e0/(AC+TSAM));
	if (rand2_(2.0e0)*GTMAX > GTR)
		goto L1;
	CDTS=1.0e0-TSAM;

}

void gppa2_(double &EE, double &CDTE, double &EP, double &CDTP, int &IZZ, int &ISH){

 /*
	Amostragem aleatória do par elétron-pósitron e produção de tripletos por
	fótons. Seção transversal diferencial Bethe-Heitler.

	Valores de saída :
	EE ..... energia cinética do elétron.
	CDTE ... direção polar cosseno do elétron.
	EP ..... energia cinética do pósitron.
	CDTP ... cosseno de direção polar do pósitron.
	IZZ ... número atômico do átomo onde ocorreu a absorção.
	ISH .... camada de elétrons atômicos que foi ionizada.
 */

	/*if (imprimiu==0){
		printf("\n\nGPPA2\n\n");
		imprimiu++;
	}*/

	double REV=5.10998928e5;
	double SL=137.035999074e0;
	double  TREV=2.0e0*REV;

	static const int NOCO=512;

	double EKI, EPS, ALZ, T, F00, G0, BMIN, G1MIN, G2MIN, XR, A1, A2, P1, RU2M1, B, G1, G2, TRIPL, TST;
	int I, ISHELL, JO;

	EKI=REV/ *TRACK_mod_.E;

	if (*TRACK_mod_.E < 1.1e6){
		EPS=EKI+(1.0e0-2.0e0*EKI)*rand2_(1.0e0);
		goto L3;
	}

	//Low-energy and Coulomb corrections.
	ALZ=CGPP00_.ZEQPP[*TRACK_mod_.MAT-1]/SL;
    T=sqrt(2.0e0*EKI);
    F00=(-1.774e0-1.210e1*ALZ+1.118e1*ALZ*ALZ)*T
        +(8.523e0+7.326e1*ALZ-4.441e1*ALZ*ALZ)*pow(T,2)
        -(1.352e1+1.211e2*ALZ-9.641e1*ALZ*ALZ)*pow(T,3)
        +(8.946e0+6.205e1*ALZ-6.341e1*ALZ*ALZ)*pow(T,4);
    G0=CGPP00_.F0[2-1][*TRACK_mod_.MAT-1]+F00;
    BMIN=4.0e0*EKI/CGPP00_.BCB[*TRACK_mod_.MAT-1];
	schiff2_(BMIN,G1,G2);
	G1MIN=G1+G0;
    G2MIN=G2+G0;
    XR=0.5e0-EKI;
    A1=6.666666666666666e-1*G1MIN*pow(XR,2);
    P1=A1/(A1+G2MIN);

	//Amostragem aleatória de EPS.

L1:;
	if (rand2_(2.020) > P1)
		goto L2;
	RU2M1=2.0e0*rand2_(3.0e0)-1.0e0;
	if (RU2M1 < 0.0e0)
		EPS=0.5e0-XR*pow(fabs(RU2M1),3.333333333333333e-1);
	else
		EPS=0.5e0+XR*pow(RU2M1,3.333333333333333e-1);
	B=EKI/(CGPP00_.BCB[*TRACK_mod_.MAT-1]*EPS*(1.0e0-EPS));
	schiff2_(B,G1,G2);
	G1=fmax(G1+G0,0.0e0);
	if (rand2_(4.0e0)*G1MIN > G1)
		goto L1;
	goto L3;
L2:;
	EPS=EKI+2.0e0*XR*rand2_(5.0e0);
    B=EKI/(CGPP00_.BCB[*TRACK_mod_.MAT-1]*EPS*(1.0e0-EPS));
    schiff2_(B,G1,G2);
    G2=fmax(G2+G0,0.0e0);
	if (rand2_(6.0e0)*G2MIN > G2)
		goto L1;
L3:;
	//Eletron
	EE=EPS**TRACK_mod_.E-REV;
    CDTE=2.0e0*rand2_(7.0e0)-1.0e0;
    A1=EE+REV;
    A2=sqrt(EE*(EE+TREV));
    CDTE=(CDTE*A1+A2)/(A1+CDTE*A2);
	//Positron
	EP=(1.0e0-EPS)* *TRACK_mod_.E-REV;
    CDTP=2.0e0*rand2_(8.0e0)-1.0e0;
    A1=EP+REV;
    A2=sqrt(EP*(EP+TREV));
    CDTP=(CDTP*A1+A2)/(A1+CDTP*A2);

	//Produção de Tripletos
	TRIPL=CGPP01_.TRIP[*CEGRID_.KE-1][*TRACK_mod_.MAT-1]+(CGPP01_.TRIP[*CEGRID_.KE+1-1][*TRACK_mod_.MAT-1]-CGPP01_.TRIP[*CEGRID_.KE-1][*TRACK_mod_.MAT-1])* *CEGRID_.XEK;
    IZZ=0;
    ISH=30;
	if (TRIPL < 1.0e-5)
		return;
	if (rand2_(9.0e0) > TRIPL)
		return;
	TST=rand2_(10.0e0);
	//Busca binaria
	if (TST < CGCO_.PTRSH[1-1][*TRACK_mod_.MAT-1]){
		ISHELL=1;
	}else{
		ISHELL=1;
        JO=CGCO_.NOSCCO[*TRACK_mod_.MAT-1]+1;
L4:;
    	I=(ISHELL+JO)/2;
		if (TST > CGCO_.PTRSH[I-1][*TRACK_mod_.MAT-1]){
			ISHELL=I;

		}else{
			JO=I;
		}
		if (JO-ISHELL > 1)
			goto L4;
		ISHELL=ISHELL+1;
	}
	IZZ=CGCO_.KZCO[ISHELL-1][*TRACK_mod_.MAT-1];
    ISH=CGCO_.KSCO[ISHELL-1][*TRACK_mod_.MAT-1];


}

void schiff2_(double &B, double &G1, double &G2){

	/*
	Funções de triagem F1(B) e F2(B) no diferencial Bethe-Heitler
	seção transversal para produção de pares.
	*/

	/*if (imprimiu==0){
		printf("\n\nSCHIFF2\n\n");
		imprimiu++;
	}*/

	double PI=3.1415926535897932e0;
	double TWOPI=PI+PI;
	double B1, F1, F2, A0, B2;

	B2=B*B;
    F1=2.0e0-2.0e0*log(1.0e0+B2);
    F2=F1-6.666666666666666e-1;
	if (B < 1.0e-10){
		F1=F1-TWOPI*B;
	}else{
		A0=4.0e0*B*atan2(1.0e0,B);
        F1=F1-A0;
        F2=F2+2.0e0*B2*(4.0e0-A0-3.0e0*log((1.0e0+B2)/B2));
	}
	G1=0.5e0*(3.0e0*F1-F2);
    G2=0.25e0*(3.0e0*F1+F2);

}

void gaux2_(){
	/*	
	Mecanismo de interação auxiliar para fotons, definível pelo usuário.
	Geralmente não está ativo.
	*/
	printf("Warning: Subroutine GAUX has been entered.\n");
}

void peld2_(double &RNDC, double &RMU){

	/*
	Simulação de eventos elásticos rígidos de positrons. Seções transversais de
 	a base de dados numérica ELSEPA.

	Valor do argumento :
 	RNDC ... valor de corte do número aleatório uniforme
 	(apenas eventos difíceis são simulados).
 	RMU .... deflexão angular amostrada, =(1-CDT)/2.
	
	*/

	/*if (imprimiu==0){
		printf("\n\nPELD2\n\n");
		imprimiu++;
	}*/

	static const int NP=128;
	static const int NPM1=NP-1;

	double PK, RU, RR, PP, XX, AA, BB, D;
	int ITN, JE, I, J, K;

	//Ponto da rede de energia


	PK=(*CEGRID_.XEL-CEGRID_.DLEMP[*CEGRID_.KE-1])* *CEGRID_.DLFC;

	if (rand2_(1.0e0) < PK){
		JE=*CEGRID_.KE+1;
	}else{
		JE=*CEGRID_.KE;
	}

	//Ponto
	RU=RNDC+rand2_(2.0e0)*(1.0e0-RNDC);

	//Selection of the interval (binary search in a restricted interval).
	ITN=RU*NPM1+1;
    I=CPELDB_.ITLP[*TRACK_mod_.MAT-1][JE-1][ITN-1];
    J=CPELDB_.ITUP[*TRACK_mod_.MAT-1][JE-1][ITN-1];
	if ((J-I) < 2)
		goto L2;
L1:;
	K=(I+J)/2;
	if (RU > CPELDB_.PSP[*TRACK_mod_.MAT-1][JE-1][K-1]){
		I=K;
	}else{
		J=K;
	}

	if ((J-I) > 1)
		goto L1;

	//Amostragem da distribuição cumulativa inversa racional.
L2:;
	PP=CPELDB_.PSP[*TRACK_mod_.MAT-1][JE-1][I-1];
    RR=RU-PP;
	if (RR > 1.0e-16){
		XX=CPELDB_.XSP[*TRACK_mod_.MAT-1][JE-1][I-1];
        AA=CPELDB_.ASP[*TRACK_mod_.MAT-1][JE-1][I-1];
        BB=CPELDB_.BSP[*TRACK_mod_.MAT-1][JE-1][I-1];
        D=CPELDB_.PSP[*TRACK_mod_.MAT-1][JE-1][I+1-1]-PP;
        RMU=XX+((1.0e0+AA+BB)*D*RR/(D*D+(AA*D+BB*RR)*RR))*(CPELDB_.XSP[*TRACK_mod_.MAT-1][JE-1][I+1-1]-XX);
	}else{
		RMU=CPELDB_.XSP[*TRACK_mod_.MAT-1][JE-1][I-1];
	}
}

void pina2_(double &E, double &DELTA, double &DE, double &EP, double &CDT, double &ES, double &CDTS, int &M, int &IOSC){

	/*
	Amostragem aleatória de colisões inelásticas duras de positrons.

	Modelo C Sternheimer-Liljequist GOS.

	Argumentos de entrada:
	E ....... energia do elétron (eV).
	M ....... material onde os elétrons se propagam.
	DELTA ... Correção do efeito de densidade de Fermi.
	Argumentos de saída:
	DE ...... perda de energia (eV).
	EP ...... energia do elétron espalhado (eV).
	CDT ..... cosseno do ângulo de dispersão polar.
	ES ...... energia do elétron secundário emitido (eV).
	CDTS .... cosseno polar de direção do elétron secundário.
	IOSC .... índice do oscilador que foi 'ionizado'.
	
	*/

	/*if (imprimiu==0){
		printf("\n\nPINA2\n\n");
		imprimiu++;
	}*/

	double REV=5.10998928e5;
	bool LDIST;
	double RREV=1.0e0/REV;
	double TREV=2.0e0*REV;
	double RTREV=1.0e0/TREV;

	static const int NO=512;

	double WCCM, PK, TST, UK, WK,WTHR, WM, WKP, QKP, WCMAX, WDMAX, RB, GAM, GAM2, BETA2, AMOL, CPS, CP;
	double CPPS, CPP, QM, RWKP, XHDL, XHDT, F0, RCL, RL1, RRL1, XHC, XHTOT, TS1, A, ARCL, FB, RK, RK2, RKF, PHI;
	double QS, Q, QTREV, G12, BHA1, BHA2, BHA3, BHA4;
	int IO, JO, IT, JE;

	WCCM=PENELOPE_mod_.WCC[M-1];

	if (WCCM > E){
		DE=0.0e0;
        EP=E;
        CDT=1.0e0;
        ES=0.0e0;
        CDTS=0.0e0;
        IOSC=NO;
        return;
	}
	//Ponto da rede de energia
	PK=(*CEGRID_.XEL-CEGRID_.DLEMP[*CEGRID_.KE-1])* *CEGRID_.DLFC;
	if (rand2_(1.0e0) < PK){
		JE=*CEGRID_.KE+1;
	}else{
		JE=*CEGRID_.KE;
	}

	//Seleção do oscilador ativo.
	TST=rand2_(2.0e0);
	//Busca binaria
	IO=1;
    JO=CPINAC_.NPIN[M-1]+1;
L1:;
    IT=(IO+JO)/2;
	if (TST > CPINAC_.PINAC[IT-1][JE-1][M-1]){
		IO=IT;
	}else{
		JO=IT;
	}
	if (JO-IO > 1)
		goto L1;
	IOSC=CPINAC_.IPIN[IO-1][M-1];
    UK=CEIN_.UI[IOSC-1][M-1];
    WK=CEIN_.WRI[IOSC-1][M-1];
	
	if (UK > 1.0e-3){
		WTHR=fmax(WCCM,UK);
	}else{
		WTHR=fmax(WCCM,WK);
	}

	if (E < (WTHR+1.0e-6)){
		DE=0.0e0;
        EP=E;
        CDT=1.0e0;
        ES=0.0e0;
        CDTS=0.0e0;
        IOSC=NO;
        return;
	}

	/*
	Truque: A energia de ressonância e a energia de recuo de corte de
	conchas internas são variadas para produzir um limiar suave.
	*/

	LDIST=true;
	if (UK > 1.0e-3){
		WM=3.0e0*WK-2.0e0*UK;
		if (E > WM){
			WKP=WK;
            QKP=UK;
		}else{
			WKP=(E+2.0e0*UK)/3.0e0;
            QKP=UK*(E/WM);
            WM=E;
		}
		if (WCCM > WM)
		    LDIST=false;
        WCMAX=E;
        WDMAX=fmin(WM,WCMAX);
        if (WTHR > WDMAX) 
			LDIST=false;
	}else{
		if (WCCM > WK)
			LDIST=false;
		WKP=WK;
        QKP=WK;
        WM=E;
        WCMAX=E;
        WDMAX=WKP+1.0e0;
	}

	//Constantes

	RB=E+TREV;
    GAM=1.0e0+E*RREV;
    GAM2=GAM*GAM;
    BETA2=(GAM2-1.0e0)/GAM2;
	G12=pow((GAM+1.0e0),2);
    AMOL=pow(((GAM-1.0e0)/GAM),2);
	BHA1=AMOL*(2.0e0*G12-1.0e0)/(GAM2-1.0e0);
    BHA2=AMOL*(3.0e0+1.0e0/G12);
    BHA3=AMOL*2.0e0*GAM*(GAM-1.0e0)/G12;
    BHA4=AMOL*pow((GAM-1.0e0),2)/G12;
    CPS=E*RB;
    CP=sqrt(CPS);

	//Seções transversais parciais do oscilador ativo.
	//Excitalçoes Distantes
	if (LDIST){
		CPPS=(E-WKP)*(E-WKP+TREV);
        CPP=sqrt(CPPS);
		if (WKP > 1.0e-6*E){
			QM=sqrt(pow((CP-CPP),2)+REV*REV)-REV;
		}else{
			QM=pow(WKP,2)/(BETA2*TREV);
            QM=QM*(1.0e0-QM*RTREV);
		}
		if (QM < QKP){
			RWKP=1.0e0/WKP;
         	XHDL=log(QKP*(QM+TREV)/(QM*(QKP+TREV)))*RWKP;
          	XHDT=fmax(log(GAM2)-BETA2-DELTA,0.0e0)*RWKP;
			if (UK > 1.0e-3){
				F0=(WDMAX-WTHR)*(WM+WM-WDMAX-WTHR)/pow((WM-UK),2);
             	XHDL=F0*XHDL;
             	XHDT=F0*XHDT;
			}
		}else{
			XHDL=0.0e0;
            XHDT=0.0e0;
		}
	}else{
		QM=0.0e0;    //Definido para evitar avisos de compilação.
        CPP=0.0e0;   
        CPPS=0.0e0;  
        XHDL=0.0e0;
        XHDT=0.0e0;
	}

	//Colisoes Fechadas
	RCL=WTHR/E;
	RL1=1.0e0-RCL;
      XHC=((1.0e0/RCL-1.0e0)+BHA1*log(RCL)+BHA2*RL1
     	 +(BHA3/2.0e0)*(pow(RCL,2)-1.0e0)+(BHA4/3.0e0)*(1.0e0-pow(RCL,3)))/E;

	XHTOT=XHC+XHDL+XHDT;
	if (XHTOT < 1.0e-35){
		DE=0.0e0;
        EP=E;
        CDT=1.0e0;
        ES=0.0e0;
        CDTS=0.0e0;
        IOSC=NO;
        return;
	}

	//Amostragem de variáveis ​​de estado final.

	TST=rand2_(3.0e0)*XHTOT;

	//Colisão fechada dura

	TS1=XHC;
	if (TST < TS1){
L2:;
		RK=RCL/(1.0e0-rand2_(4.0e0)*RL1);
        PHI=1.0e0-RK*(BHA1-RK*(BHA2-RK*(BHA3-BHA4*RK)));
		if (rand2_(5.0e0) > PHI)
			goto L2;
		//Energia e ângulo de espalhamento (positron primário).
		DE=RK*E;
        EP=E-DE;
        CDT=sqrt(EP*RB/(E*(RB-DE)));
		//Energia e ângulo de emissão do raio delta.
		if (CEIN_.KS[IOSC-1][M-1] < 17){
			if (UK > CECUTR_.ECUTR[M-1]){
				ES=DE-UK;
			}else{
				ES=DE;
			}
		}else{
			ES=DE;
		}
		CDTS=sqrt(DE*RB/(E*(DE+TREV)));
		return;
	}

	//Interação longitudinal dura distante.
	TS1=TS1+XHDL;
	if (UK > 1.0e-3){
		DE=WM-sqrt(pow((WM-WTHR),2)-rand2_(7.0e0)*(WDMAX-WTHR)*(WM+WM-WDMAX-WTHR));
	}else{
		DE=WKP;
	}
	EP=E-DE;
	if (TST < TS1){
		QS=QM/(1.0e0+QM*RTREV);
        Q=QS/(pow(((QS/QKP)*(1.0e0+QKP*RTREV)), rand2_(6.0e0))-(QS*RTREV));
        QTREV=Q*(Q+TREV);
        CDT=(CPPS+CPS-QTREV)/(2.0e0*CP*CPP);
		if (CDT > 1.0e0){
			CDT=1.0e0;
		}
		//Energia e ângulo de emissão do raio delta.
		if (CEIN_.KS[IOSC-1][M-1] < 17){
			ES=DE-UK; //Apenas conchas internas.
		}else{
			ES=DE;
		}
		CDTS=0.5e0*(WKP*(E+RB-WKP)+QTREV)/sqrt(CPS*QTREV);
		if (CDTS > 1.0e0)
		    CDTS=1.0e0;
		return;
	}

	//Interação transversal distante difícil.
	CDT=1.0e0;
	//Energia e ângulo de emissão do raio delta.
	if (CEIN_.KS[IOSC-1][M-1] < 17){
		if (UK >  CECUTR_.ECUTR[M-1]){
			ES=DE-UK; //Apenas conchas internas.
		}else{
			ES=DE;
		}
	}else{
		ES=DE;
	}
	CDTS=1.0e0;

}

void psia2_(double &E, double &DELTA, double &DE, double &EP, double &CDT, double &ES, double &CDTS, int &M, int &IZZ, int &ISH){

/*
	Amostragem aleatória de ionização de camada interna por impacto de positrons.

	Modelo Sternheimer-Liljequist GOS.

	Argumentos de entrada:
	E ....... energia do elétron (eV).
	M ....... material onde os elétrons se propagam.
	DELTA ... Correção do efeito de densidade de Fermi.
	Argumentos de saída:
	DE ...... perda de energia (eV).
	EP ...... energia do elétron espalhado (eV).
	CDT ..... cosseno do ângulo de dispersão polar.
	ES ...... energia do elétron secundário emitido (eV).
	CDTS .... cosseno polar de direção do elétron secundário.
	IZZ ..... número atômico do elemento onde a ionização ocorreu.
	ISH ..... camada de elétrons atômicos que foi ionizada.
*/

	/*if (imprimiu==0){
		printf("\n\nPSIA2\n\n");
		imprimiu++;
	}*/


	double REV=5.10998928e5;
	double RREV=1.0e0/REV;

	double TREV=2.0e0*REV;
	double RTREV=1.0e0/TREV;

	static const int NO=512;

	double WCCM, PK, TST, UK, WK,WTHR, WM, WKP, QKP, WCMAX, WDMAX, RB, GAM, GAM2, BETA2, AMOL, CPS, CP;
	double CPPS, CPP, QM, RWKP, XHDL, XHDT, F0, RCL, RL1, RRL1, XHC, XHTOT, TS1, A, ARCL, FB, RK, RK2, RKF, PHI;
	double QS, Q, QTREV, G12, BHA1, BHA2, BHA3, BHA4;
	int IO, JO, IT, IOSC, JE;

	//Ponto da rede de energia

	PK=(*CEGRID_.XEL-CEGRID_.DLEMP[*CEGRID_.KE-1])* *CEGRID_.DLFC;
	if (rand2_(1.0e0) < PK){
		JE=*CEGRID_.KE+1;
	}else{
		JE=*CEGRID_.KE;
	}

	//Seleção do oscilador ativo.

	TST=rand2_(2.0e0);

	//Busca binaria
	IO=1;
    JO=CPSIAC_.NPSI[M-1]+1;
L1:;
    IT=(IO+JO)/2;
	if (TST > CPSIAC_.PSIAC[IT-1][JE-1][M-1]){
		IO=IT;
	}else{
		JO=IT;
	}
	if ((JO-IO) > 1)
		goto L1;
	
	IOSC=CPSIAC_.IPSI[IO-1][M-1];
    IZZ=CEIN_.KZ[IOSC-1][M-1];
    ISH=CEIN_.KS[IOSC-1][M-1];
    UK=CEIN_.UI[IOSC-1][M-1];
    WK=CEIN_.WRI[IOSC-1][M-1];

	if (UK > 1.0e-3){
		WTHR=UK;
	}else{
		WTHR=WK;
	}

	if (E < WTHR+1.0e-6){
		DE=UK;
        EP=E-DE;
        CDT=1.0e0;
        ES=0.0e0;
        CDTS=0.0e0;
        return;
	}

	/*
	Truque: A energia de ressonância e a energia de recuo de corte de
	 conchas internas são variadas para produzir um limiar suave.
	*/

 	WM=3.0e0*WK-2.0e0*UK;
	if (E > WM){
		WKP=WK;
        QKP=UK;
	}else{
		WKP=(E+2.0e0*UK)/3.0e0;
        QKP=UK*(E/WM);
        WM=E;
	}
	
    WCMAX=E;
    WDMAX=fmin(WM,WCMAX);


	//Constantes

	RB=E+TREV;
    GAM=1.0e0+E*RREV;
    GAM2=GAM*GAM;
    BETA2=(GAM2-1.0e0)/GAM2;
	G12=pow((GAM+1.0e0),2);
    AMOL=pow(((GAM-1.0e0)/GAM),2);
	BHA1=AMOL*(2.0e0*G12-1.0e0)/(GAM2-1.0e0);
    BHA2=AMOL*(3.0e0+1.0e0/G12);
    BHA3=AMOL*2.0e0*GAM*(GAM-1.0e0)/G12;
    BHA4=AMOL*pow((GAM-1.0e0),2)/G12;
    CPS=E*RB;
    CP=sqrt(CPS);

	//Seções transversais parciais do oscilador ativo.

	//Excitações distantes.

	CPPS=(E-WKP)*(E-WKP+TREV);
    CPP=sqrt(CPPS);

	if (WKP > 1.0e-6*E){
		QM=sqrt(pow((CP-CPP),2)+REV*REV)-REV;
	}else{
		QM=pow(WKP,2)/(BETA2*TREV);
        QM=QM*(1.0e0-QM*RTREV);
	}

	if (QM < QKP){
		RWKP=1.0e0/WKP;
        XHDL=log(QKP*(QM+TREV)/(QM*(QKP+TREV)))*RWKP;
        XHDT=fmax(log(GAM2)-BETA2-DELTA,0.0e0)*RWKP;
        F0=(WDMAX-WTHR)*(WM+WM-WDMAX-WTHR)/pow((WM-UK),2);
        XHDL=F0*XHDL;
        XHDT=F0*XHDT;
	}else{
		XHDL=0.0e0;
        XHDT=0.0e0;
	}

	//Colisoes Proximas
	RCL=WTHR/E;
    RL1=1.0e0-RCL;
    RRL1=1.0e0/RL1;
    XHC=((1.0e0/RCL-1.0e0)+BHA1*log(RCL)+BHA2*RL1
       +(BHA3/2.0e0)*(pow(RCL,2)-1.0e0)+(BHA4/3.0e0)*(1.0e0-pow(RCL,3)))/E;

    XHTOT=XHC+XHDL+XHDT;

	if (XHTOT < 1.0e-35){
		DE=UK;
        EP=E-DE;
        CDT=1.0e0;
        ES=0.0e0;
        CDTS=0.0e0;
        return;
	}

	//Amostragem de variáveis ​​de estado final.

	TST=rand2_(3.0e0)*XHTOT;

	//Colisão aproximada dura.

	TS1=XHC;
	if (TST < TS1){
L2:;
		RK=RCL/(1.0e0-rand2_(4.0e0)*RL1);
        PHI=1.0e0-RK*(BHA1-RK*(BHA2-RK*(BHA3-BHA4*RK)));
		if (rand2_(5.0e0) > PHI)
			goto L2;
		//Energia e ângulo de espalhamento (elétron primário).
		DE=RK*E;
		EP=E-DE;
		CDT=sqrt(EP*RB/(E*(RB-DE)));
		//Energia e ângulo de emissão do raio delta.
		ES=DE-UK;
		CDTS=sqrt(DE*RB/(E*(DE+TREV)));
		return;
	}

	//Interação longitudinal dura distante.
	TS1=TS1+XHDL;
    DE=WM-sqrt(pow((WM-WTHR),2)-rand2_(7.0e0)*(WDMAX-WTHR)*(WM+WM-WDMAX-WTHR));
    EP=E-DE;

	if (TST < TS1){
		QS=QM/(1.0e0+QM*RTREV);
        Q=QS/(pow(((QS/QKP)*(1.0e0+QKP*RTREV)),rand2_(6.0e0))-(QS*RTREV));
        QTREV=Q*(Q+TREV);
        CDT=(CPPS+CPS-QTREV)/(2.0e0*CP*CPP);
		if (CDT > 1.0e0)
			CDT=1.0e0;

		//Energia e ângulo de emissão do raio delta.
		ES=DE-UK;
        CDTS=0.5e0*(WKP*(E+RB-WKP)+QTREV)/sqrt(CPS*QTREV);
		if (CDTS >1.0e0)
			CDTS=1.0e0;
		return;
	}

	//Interação transversal distante difícil.

	CDT=1.0e0;

	//Energia e ângulo de emissão do raio delta.

	ES=DE-UK;
    CDTS=1.0e0;

}

void pana2_(double &E, double &E1, double &CDT1, double &E2, double &CDT2, int &M){

	/*
	Simulação de aniquilação de pósitrons (em repouso ou em voo) em
	material M. Ei e CDTi são as energias e os cossenos da direção polar
	dos dois fótons de aniquilação.
	*/

	/*if (imprimiu==0){
		printf("\n\nPANA2\n\n");
		imprimiu++;
	}*/

	double REV=5.10998928e5;
	double TREV=2.0e0*REV;

	double GAM, GAM21, ANI, CHIMIN, RCHI, GT0, CHI, GREJ, DET, CHIP;

	//Pósitrons lentos (assumidos em repouso).

	if (E < PENELOPE_mod_.EABS[M-1][3-1]){
		E1=0.5e0*(E+TREV);
        E2=E1;
        CDT1=-1.0e0+2.0e0*rand2_(1.0e0);
        CDT2=-CDT1;

	}else{
		/*
		Aniquilação em vôo (dois fótons com energia e direções
		determinado a partir da dcs e conservação energia-momento).
		*/
		GAM=1.0e0+fmax(E,1.0e0)/REV;
        GAM21=sqrt(GAM*GAM-1.0e0);
        ANI=1.0e0+GAM;
        CHIMIN=1.0e0/(ANI+GAM21);
        RCHI=(1.0e0-CHIMIN)/CHIMIN;
        GT0=ANI*ANI-2.0e0;
L1:;
        CHI=CHIMIN*pow(RCHI,rand2_(2.0e0));
        GREJ=ANI*ANI*(1.0e0-CHI)+GAM+GAM-1.0e0/CHI;
		if (rand2_(3.0e0)*GT0 > GREJ)
			goto L1;
		
		DET=E+TREV;
        E1=CHI*DET;
        CDT1=(ANI-1.0e0/CHI)/GAM21;
        CHIP=1.0e0-CHI;
        E2=DET-E1;
        CDT2=(ANI-1.0e0/CHIP)/GAM21;
	}

}

void paux2_(){
	/*	
	Mecanismo de interação auxiliar para fotons, definível pelo usuário.
	Geralmente não está ativo.
	*/
	printf("Warning: Subroutine PAUX has been entered.\n");
}

void tenang2_(int &IEXIT, int &N){

	/*
	Calcula energia e distribuições angulares de partículas emergentes,
	grava e carrega arquivos de despejo, acumula arquivos de despejo de diferentes
	é executado e grava os resultados.
	*/

	/*if (imprimiu==0){
		printf("\n\nTENANG2\n\n");
		imprimiu++;
	}*/

	double PI=3.1415926535897932e0;
	double TWOPI=2.0e0*PI;
	double RA2DE=180.0e0/PI;
	double DE2RA=PI/180.0e0;
	double FSAFE=1.000000001e0;
	static const int NBEM=1500; //
	static const int NBTHM=1800; //numero maximo de caixas de energia.
	static const int NBPHM=180; //numero maximo de caixas angulares

	//Pontue as contribuições de uma nova partícula.

	//Distribuição de energia de partículas emergentes.

	int KEn, KTH, KPH;
	double THETA, PHI;

	if (*CENANG_.LLE == 1){
		KEn=1.0e0+(log(*TRACK_mod_.E)-*CENANG_.EL)* *CENANG_.RBSE;

	}else{
		KEn=1.0e0+(*TRACK_mod_.E-*CENANG_.EL)* *CENANG_.RBSE;
	}
	if ((KEn > 0) && (KEn <= *CENANG_.NE)){
		if (N != CENANG_.LPDE[KEn-1][IEXIT-1][*TRACK_mod_.KPAR-1]){
			CENANG_.PDE[KEn-1][IEXIT-1][*TRACK_mod_.KPAR-1]=CENANG_.PDE[KEn-1][IEXIT-1][*TRACK_mod_.KPAR-1]+CENANG_.PDEP[KEn-1][IEXIT-1][*TRACK_mod_.KPAR-1];
			CENANG_.PDE2[KEn-1][IEXIT-1][*TRACK_mod_.KPAR-1]=CENANG_.PDE2[KEn-1][IEXIT-1][*TRACK_mod_.KPAR-1]+pow(CENANG_.PDEP[KEn-1][IEXIT-1][*TRACK_mod_.KPAR-1],2);
			CENANG_.PDEP[KEn-1][IEXIT-1][*TRACK_mod_.KPAR-1]=*TRACK_mod_.WGHT;
			CENANG_.LPDE[KEn-1][IEXIT-1][*TRACK_mod_.KPAR-1]=N;	
		}else{
			CENANG_.PDEP[KEn-1][IEXIT-1][*TRACK_mod_.KPAR-1]=CENANG_.PDEP[KEn-1][IEXIT-1][*TRACK_mod_.KPAR-1]+*TRACK_mod_.WGHT;
		}
	}

	//Distribuição angular de partículas emergentes.
	THETA=acos(*TRACK_mod_.W);
	if (*CENANG_.LLTH == 1){
		KTH=1.0e0+(log(fmax(THETA,1.0e-12)*RA2DE)-*CENANG_.THL)* *CENANG_.RBSTH;
		if (KTH < 1)
			KTH=1;
	}else{
		KTH=1.0e0+THETA*RA2DE* *CENANG_.RBSTH;
	}
	if (fabs(*TRACK_mod_.U) > 1.0e-16){
		PHI=atan2(*TRACK_mod_.V,*TRACK_mod_.U);
	}else if (fabs(*TRACK_mod_.V) > 1.0e-16){
		PHI=atan2(*TRACK_mod_.V,*TRACK_mod_.U);
	}else{
		PHI=0.0e0;
	}
	if (PHI < 0.0e0)
		PHI=TWOPI+PHI;
	KPH=1.0e0+PHI*RA2DE* *CENANG_.RBSPH;
	if (N != CENANG_.LPDA[KPH-1][KTH-1][*TRACK_mod_.KPAR-1]){
		CENANG_.PDA[KPH-1][KTH-1][*TRACK_mod_.KPAR-1]=CENANG_.PDA[KPH-1][KTH-1][*TRACK_mod_.KPAR-1]+CENANG_.PDAP[KPH-1][KTH-1][*TRACK_mod_.KPAR-1];
		CENANG_.PDA2[KPH-1][KTH-1][*TRACK_mod_.KPAR-1]=CENANG_.PDA2[KPH-1][KTH-1][*TRACK_mod_.KPAR-1]+pow(CENANG_.PDAP[KPH-1][KTH-1][*TRACK_mod_.KPAR-1],2);
		CENANG_.PDAP[KPH-1][KTH-1][*TRACK_mod_.KPAR-1]=*TRACK_mod_.WGHT;
		CENANG_.LPDA[KPH-1][KTH-1][*TRACK_mod_.KPAR-1]=N;	
	}else{
		CENANG_.PDAP[KPH-1][KTH-1][*TRACK_mod_.KPAR-1]=CENANG_.PDAP[KPH-1][KTH-1][*TRACK_mod_.KPAR-1]+*TRACK_mod_.WGHT;
	}

}

void sendet2_(double &ED, int &ID){

   /* if (imprimiu==0){
		printf("\n\nSENDET2\n\n");
		imprimiu++;
	}*/

	/*
	Espectro de energia depositada.
	ED inclui o peso da partícula.
	*/
	int IE;


	CENDET_.EDEP[ID-1]=CENDET_.EDEP[ID-1]+ED;
    CENDET_.EDEP2[ID-1]=CENDET_.EDEP2[ID-1]+pow(ED,2);
	if (ED > 1.0e-5){
		if (CENDET_.LLE[ID-1] == 1){
			IE=1.0e0+(log(ED)-CENDET_.EL[ID-1])*CENDET_.RBSE[ID-1];
		}else{
			IE=1.0e0+(ED-CENDET_.EL[ID-1])*CENDET_.RBSE[ID-1];
		}
		if ((IE > 0) && (IE <= CENDET_.NE[ID-1]))
			CENDET_.DET[IE-1][ID-1]=CENDET_.DET[IE-1][ID-1]+1.0e0;
	}

}

void secpar2_(int &LEFT){
	/*
	Esta sub-rotina entrega o estado inicial de uma partícula secundária
	produzido durante a simulação anterior do chuveiro. Esta partícula
	é removido da pilha secundária, de modo que será perdido se um novo
	A chamada  para SECPAR é realizada antes de simular sua trajetória até
	o fim.

	LEFT é o número de partículas na pilha secundária na chamada
	hora . Quando LEFT=0, a simulação do chuveiro foi concluída.
	*/

 	/*if (imprimiu==0){
		printf("\n\nSECPAR2\n\n");
		imprimiu++;
	}*/

	if (*SECST_.NSEC > 0){
		LEFT=*SECST_.NSEC;
        *TRACK_mod_.E=SECST_.ES[*SECST_.NSEC-1];
        *TRACK_mod_.X=SECST_.XS[*SECST_.NSEC-1];
        *TRACK_mod_.Y=SECST_.YS[*SECST_.NSEC-1];
        *TRACK_mod_.Z=SECST_.ZS[*SECST_.NSEC-1];
        *TRACK_mod_.U=SECST_.US[*SECST_.NSEC-1];
        *TRACK_mod_.V=SECST_.VS[*SECST_.NSEC-1];
        *TRACK_mod_.W=SECST_.WS[*SECST_.NSEC-1];
        *TRACK_mod_.WGHT=SECST_.WGHTS[*SECST_.NSEC-1];
        *TRACK_mod_.KPAR=SECST_.KS[*SECST_.NSEC-1];
        *TRACK_mod_.IBODY=SECST_.IBODYS[*SECST_.NSEC-1];
        *TRACK_mod_.MAT=SECST_.MS[*SECST_.NSEC-1];
        *TRACK_mod_.IPOL=SECST_.IPOLS[*SECST_.NSEC-1];
        *TRACK_mod_.SP1=SECST_.SP1S[*SECST_.NSEC-1];
        *TRACK_mod_.SP2=SECST_.SP2S[*SECST_.NSEC-1];
        *TRACK_mod_.SP3=SECST_.SP3S[*SECST_.NSEC-1];
        *TRACK_mod_.PAGE=SECST_.PAGES[*SECST_.NSEC-1];
		for (int I = 1; I <= 5; I++){
			TRACK_mod_.ILB[I-1]=SECST_.ILBS[*SECST_.NSEC-1][I-1];
		}
		*SECST_.NSEC=*SECST_.NSEC-1;
	}else{
		LEFT=0;
	}
}


void pmwrt2_(int ICLOSE){
	/*
	Calcula médias e grava resultados em arquivos de saída.

	ICLOSE é um sinalizador de fechamento, que é usado somente quando o espaço de fase
	O arquivo de um detector de impacto está sendo gerado.
	-- Quando ICLOSE .GT. 0, as partículas restantes no buffer são
	transferido para o psf e a unidade de saída do psf é fechada.
	-- Se ICLOSE .LT. 0, as partículas são movidas para o psf, mas o psf
	A unidade de saída C permanece aberta
		*/

	double PI=3.1415926535897932e0;
	double RA2DE=180.0e0/PI;

	double FT, ERR1, ERR2, ERR, FB, FA, DF, QER, QAV, EFFIC, PTOT, YAV, YERR, EINTL, FACT, FNT;

	double WSEC[3][3];
	double WSEC2[3][3];
	double WAVE2[2], WAVE[2], WAVW2[2],WAVW[2],WAVA2[2],WAVA[2];

	//Se 'DUMPTO' estiver ativo, grave os contadores em um arquivo de despejo.

	if (*CDUMP_.LDUMP){
		FILE* IWR = fopen(CDUMP_.PFILED, "w");
		if (IWR == NULL){
			printf("Não foi possível abrir o arquivo %s", CDUMP_.PFILED);
			exit(0);
		}

		fprintf(IWR, "	%.16f		%.16f\n", *CNTRL_.SHN,*CNTRL_.TSIM);
		fprintf(IWR, "%s\n",  CTITLE_.TITLE);
		fprintf(IWR, "	%d		%d\n", *RSEED_.ISEED1, *RSEED_.ISEED2);
		fprintf(IWR, "	%d		%.16f\n", *CSOUR4_.NPSN,*CSOUR4_.RLREAD);
		fprintf(IWR, "	%d\n", *CSOUR0_.KPARP);
		if (*CSOUR0_.KPARP == 0){
			fprintf(IWR, "	%d		%16.f		%16.f\n", *CNT3_.NSDE,*CNT3_.DSDE,*CNT3_.RDSDE);
			for (int I = 1; I <= *CNT3_.NSDE; I++){
				for (int K = 1; K <= 3; K++){
					fprintf(IWR, "	%11.f	",CNT3_.SEDS[I-1][K-1]);
				}
			}
			fprintf(IWR, "\n");

			for (int I = 1; I <= *CNT3_.NSDE; I++){
				for (int K = 1; K <= 3; K++){
					fprintf(IWR, "	%16.f	",CNT3_.SEDS2[I-1][K-1]);
				}
			}
			fprintf(IWR, "\n");
		}else{
			fprintf(IWR, "	%d	\n", *CNT2_.NSEB);
			for (int I = 1; I <= *CNT2_.NSEB+1; I++){
				fprintf(IWR, "	%.16f	", CSOUR2_.ESRC[I-1]);
			}
			fprintf(IWR, "\n");

			for (int I = 1; I <= *CNT2_.NSEB+1; I++){
				fprintf(IWR, "	%.16f	", CSOUR2_.PSRC[I-1]);
			}
			fprintf(IWR, "\n");

			for (int I = 1; I <= *CNT2_.NSEB+1; I++){
				fprintf(IWR, "	%.16f	", CNT2_.SHIST[I-1]);
			}
			fprintf(IWR, "\n");
		}
		for (int I = 1; I <= 3; I++){
			fprintf(IWR, "	%.16f	", CNT0_.PRIM[I-1]);
		}
		fprintf(IWR, "\n");

		for (int I = 1; I <= 3; I++){
			fprintf(IWR, "	%.16f	", CNT0_.PRIM2[I-1]);
		}
		fprintf(IWR, "\n");

		for (int I = 1; I <= 3; I++){
			for (int K = 1; K <= 3; K++){
				fprintf(IWR, "	%.16f	", CNT0_.SEC[I-1][K-1]);
			}
		}
		fprintf(IWR, "\n");

		for (int I = 1; I <= 3; I++){
			for (int K = 1; K <= 3; K++){
				fprintf(IWR, "	%.16f	", CNT0_.SEC2[I-1][K-1]);
			}
		}
		fprintf(IWR, "\n");

		for (int I = 1; I <= 2; I++){
			fprintf(IWR, "	%.16f	", CNT0_.AVW[I-1]);
		}
		fprintf(IWR, "\n");

		for (int I = 1; I <= 2; I++){
			fprintf(IWR, "	%.16f	", CNT0_.AVW2[I-1]);
		}
		fprintf(IWR, "\n");

		for (int I = 1; I <= 2; I++){
			fprintf(IWR, "	%.16f	", CNT0_.AVA[I-1]);
		}
		fprintf(IWR, "\n");

		for (int I = 1; I <= 2; I++){
			fprintf(IWR, "	%.16f	", CNT0_.AVA2[I-1]);
		}
		fprintf(IWR, "\n");

		for (int I = 1; I <= 2; I++){
			fprintf(IWR, "	%.16f	", CNT0_.AVE[I-1]);
		}
		fprintf(IWR, "\n");

		for (int I = 1; I <= 2; I++){
			fprintf(IWR, "	%.16f	", CNT0_.AVE2[I-1]);
		}
		fprintf(IWR, "\n");

		fprintf(IWR, " 	%d	\n", *PENGEOM_mod_.NBODY);

		for (int I = 1; I <= *PENGEOM_mod_.NBODY; I++){
			fprintf(IWR, "	%.16f	", CNT1_.TDEBO[I-1]);
		}
		fprintf(IWR, "\n");

		for (int I = 1; I <= *PENGEOM_mod_.NBODY; I++){
			fprintf(IWR, "	%.16f	", CNT1_.TDEBO2[I-1]);
		}
		fprintf(IWR, "\n");

		//Energia e distribuições angulares.
		enangd2_(IWR);
		fprintf(IWR, "	%d	%d	%d	\n", *CNT4_.NID,*CNT5_.NED,*CNT6_.LDOSEM);
		if (*CNT4_.NID > 0){
			fprintf(IWR, "	%.16f	%.16f\n", *CNT4_.RLAST,*CNT4_.RWRITE);
			imdetd2_(IWR); //Detectores de Impacto
		}

		if (*CNT5_.NED > 0){
			endetd2_(IWR); //Deposito de energia nos detectores
		}

		if (*CNT6_.LDOSEM){
			dosed2_(IWR);
		}

		fprintf(IWR, "\n   *** END ***\n");
		fclose(IWR);
	}

	//Escrevendo os resultados da Simulação

	// IEXIT: 1=upbound, 2=downbound, 3=absorbed.

	FILE* IWR2 = fopen("penmain-res2.dat", "w");
	if (IWR2 == NULL){
		printf("Não foi possível abrir o arquivo penmain-res2.dat");
		exit(0);
	}
	fprintf(IWR2,"\n\n   ***********************************\n");
	fprintf(IWR2,"   **   Program PENMAIN2. Results.  **\n");
	fprintf(IWR2,"   ***********************************\n\n");

	fprintf(IWR2,"   Date and time: %s\n\n", CDATE_.DATE23);
	fprintf(IWR2, "   %s\n\n", CTITLE_.TITLE);

	fprintf(IWR2, "   Simulation time ......................... %.6E sec\n", *CNTRL_.TSIM);
	double TAVS=*CNTRL_.SHN / *CNTRL_.TSIM;
	fprintf(IWR2, "   Simulation speed ........................ %.6E showers/sec\n\n\n", TAVS);

	fprintf(IWR2, "   Simulated primary particles ............. %.6E sec\n\n", *CNTRL_.SHN);
	if (*CSOUR0_.KPARP == 1)
		fprintf(IWR2,"   Primary particles: electrons\n\n");
	if (*CSOUR0_.KPARP == 2)
		fprintf(IWR2,"   Primary particles: photons\n\n");
	if (*CSOUR0_.KPARP == 3)
		fprintf(IWR2,"   Primary particles: positrons\n\n");
	if (*CSOUR0_.KPARP == 0)
		fprintf(IWR2,"   Primary particles: set by the user subroutine SOURCE\n");

	fprintf(IWR2, "   Upbound primary particles ............... %.6E\n", CNT0_.PRIM[1-1]);
	fprintf(IWR2, "   Downbound primary particles ............. %.6E\n", CNT0_.PRIM[2-1]);
	fprintf(IWR2, "   Absorbed primary particles .............. %.6E\n\n", CNT0_.PRIM[3-1]);

	FNT=1.0e0 / *CNTRL_.SHN;
	if (*CSOUR0_.KPARP != 0){
		FT=(CNT0_.PRIM[1-1]+CNT0_.SEC[1-1][*CSOUR0_.KPARP-1])*FNT;
        ERR1=3.0e0*FNT*sqrt(fabs(CNT0_.PRIM2[1-1]-pow(CNT0_.PRIM[1-1],2)*FNT));
        ERR2=3.0e0*FNT*sqrt(fabs(CNT0_.SEC2[1-1][*CSOUR0_.KPARP-1]-pow(CNT0_.SEC[1-1][*CSOUR0_.KPARP-1],2)*FNT));
        ERR=ERR1+ERR2;
		fprintf(IWR2, "   Upbound fraction ................... %.6E +- %.1E\n", FT, ERR);

		FB=(CNT0_.PRIM[2-1]+CNT0_.SEC[2-1][*CSOUR0_.KPARP-1])*FNT;
        ERR1=3.0e0*FNT*sqrt(fabs(CNT0_.PRIM2[2-1]-pow(CNT0_.PRIM[2-1],2)*FNT));
        ERR2=3.0e0*FNT*sqrt(fabs(CNT0_.SEC2[2-1][*CSOUR0_.KPARP-1]-pow(CNT0_.SEC[2-1][*CSOUR0_.KPARP-1],2)*FNT));
        ERR=ERR1+ERR2;
		fprintf(IWR2, "   Downbound fraction ................. %.6E +- %.1E\n", FB, ERR);

		FA=CNT0_.PRIM[3-1]*FNT;
        ERR=3.0e0*FNT*sqrt(fabs(CNT0_.PRIM2[3-1]-pow(CNT0_.PRIM[3-1],2)*FNT));
		fprintf(IWR2, "   Absorption fraction ................ %.6E +- %.1E\n\n", FA, ERR);
	}

	for (int K =1; K <= 3; K++){
		for (int I = 1; I <= 3; I++){
			WSEC2[I-1][K-1]=3.0e0*FNT*sqrt(fabs(CNT0_.SEC2[I-1][K-1]-pow(CNT0_.SEC[I-1][K-1],2)*FNT));
         	WSEC[I-1][K-1]=CNT0_.SEC[I-1][K-1]*FNT;
		}
	}

	fprintf(IWR2,"   Secondary-particle generation probabilities:\n" );
	fprintf(IWR2,"                   ----------------------------------------------\n" );
	fprintf(IWR2, "                   |  electrons    |   photons     |  positrons  |\n");
	fprintf(IWR2, "   --------------------------------------------------------------\n");
	fprintf(IWR2, "   |   upbound     | %.6E | %.6E | %.6E |\n", WSEC[1-1][1-1],WSEC[1-1][2-1],WSEC[1-1][3-1]);
	fprintf(IWR2, "   |               |  +- %.1E  |  +- %.1E  |  +- %.1E  |\n",  WSEC2[1-1][1-1],WSEC2[1-1][2-1],WSEC2[1-1][3-1]);
    fprintf(IWR2, "   --------------------------------------------------------------\n");
	fprintf(IWR2, "   |   downbound   | %.6E | %.6E | %.6E |\n", WSEC[2-1][1-1],WSEC[2-1][2-1],WSEC[2-1][3-1]);
	fprintf(IWR2, "   |               |  +- %.1E  |  +- %.1E  |  +- %.1E  |\n",  WSEC2[2-1][1-1],WSEC2[2-1][2-1],WSEC2[2-1][3-1]);
    fprintf(IWR2, "   --------------------------------------------------------------\n");
	fprintf(IWR2, "   |   absorbed    | %.6E | %.6E | %.6E |\n", WSEC[3-1][1-1],WSEC[3-1][2-1],WSEC[3-1][3-1]);
	fprintf(IWR2, "   |               |  +- %.1E  |  +- %.1E  |  +- %.1E  |\n",  WSEC2[3-1][1-1],WSEC2[3-1][2-1],WSEC2[3-1][3-1]);
    fprintf(IWR2, "   --------------------------------------------------------------\n\n");

	for (int I = 1; I <=2; I++){
		DF=1.0e0/fmax(CNT0_.PRIM[I-1],1.0e0);
        WAVE2[I-1]=3.0e0*DF*sqrt(fabs(CNT0_.AVE2[I-1]-pow(CNT0_.AVE[I-1],2)*DF));
        WAVE[I-1]=CNT0_.AVE[I-1]*DF;
        WAVW2[I-1]=3.0e0*DF*sqrt(fabs(CNT0_.AVW2[I-1]-pow(CNT0_.AVW[I-1],2)*DF));
        WAVW[I-1]=CNT0_.AVW[I-1]*DF;
        WAVA2[I-1]=3.0e0*DF*RA2DE*sqrt(fabs(CNT0_.AVA2[I-1]-pow(CNT0_.AVA[I-1],2)*DF));
        WAVA[I-1]=CNT0_.AVA[I-1]*RA2DE*DF;
	}

	fprintf(IWR2, "   Average final energy:\n");
	fprintf(IWR2, "      Upbound primary particles ....... %.6E +- %.1E eV\n", WAVE[1-1],WAVE2[1-1]);
	fprintf(IWR2, "      Downbound primary particles ..... %.6E +- %.1E eV\n\n", WAVE[2-1],WAVE2[2-1]);

	fprintf(IWR2, "   Mean value of the polar cosine of the exit direction:\n");
	fprintf(IWR2, "      Upbound primary particles ....... %.6E +- %.1E\n", WAVW[1-1],WAVW2[1-1]);
	fprintf(IWR2, "      Downbound primary particles ..... %.6E +- %.1E\n\n", WAVW[2-1],WAVW2[2-1]);

	fprintf(IWR2, "   Mean value of the polar angle of the exit direction:\n");
	fprintf(IWR2, "      Upbound primary particles ....... %.6E +- %.1E deg\n", WAVA[1-1],WAVA2[1-1]);
	fprintf(IWR2, "      Downbound primary particles ..... %.6E +- %.1E deg\n\n", WAVA[2-1],WAVA2[2-1]);

	//Energias médias depositadas nos corpos..

	DF=1.0e0 / *CNTRL_.SHN;

	fprintf(IWR2, "Average deposited energies (bodies):\n");
	for (int KB = 1; KB <= *PENGEOM_mod_.NBODY; KB++){
		if (PENGEOM_mod_.MATER[KB-1] != 0){
			QER=3.0e0*DF*sqrt(fabs(CNT1_.TDEBO2[KB-1]-pow(CNT1_.TDEBO[KB-1],2)*DF));
          	QAV=CNT1_.TDEBO[KB-1]*DF;
			  if (QER > 1.0e-10*fabs(QAV)){
				  EFFIC=pow(QAV,2)/(pow((QER/3.0e0),2)* *CNTRL_.TSIM);
			  }else{
				  EFFIC=0.0e0;
			  }
			  fprintf(IWR2, "      Body %4d ...... %.6E +- %.1E eV    (effic. = %.2E)\n", KB,QAV,QER,EFFIC);
		}
	}

	//Saída de detectores de impacto.
//	if (*CNT4_.NID > 0)
		//imdetw2_(CNTRL_.SHN,CNTRL_.TSIM,IWR); Não será implementado detectores de impacto

	//Saída de detectores de deposição de energia.
	//if (*CNT5_.NED > 0)
		//endetw2_(CNTRL_.SHN,CNTRL_.TSIM,IWR); //nao será implementado dectores de deposicao de energia


	//Espectro de energia da fonte (conforme definido no PENMAIN).
	if (*CSOUR2_.LSPEC){
		FILE* IWR3 = fopen("psource2.dat", "w");
		if (IWR3 == NULL){
			printf("Não foi possível abrir o arquivo psource2.dat");
			exit(0);
		}

		fprintf(IWR3, " #  Results from PENMAIN. \n");
		fprintf(IWR3, " #  Source energy spectrum.\n");
		fprintf(IWR3, " #  1st column: E (eV). 2nd column: spectrum (1/eV).\n");
		fprintf(IWR3, " #  3rd and 4th columns: simul. pdf limits (3SD, 1/eV).\n");
		
		PTOT=0.0e0;
		fprintf(IWR3, " %.6E %.6E %.6E %.6E\n", CSOUR2_.ESRC[1-1],PTOT,PTOT,PTOT);
		for (int KEn = 1; KEn <= *CNT2_.NSEB; KEn++){
			PTOT=PTOT+CSOUR2_.PSRC[KEn-1];
		}
		for (int KEn = 1; KEn <= *CNT2_.NSEB; KEn++){
			YAV=CNT2_.SHIST[KEn-1]*DF;
          	YERR=3.0e0*sqrt(fabs(YAV*(1.0e0-YAV)*DF));
          	EINTL=CSOUR2_.ESRC[KEn+1-1]-CSOUR2_.ESRC[KEn-1];
			  if (EINTL > 1.0e-15){
				  FACT=1.0e0/EINTL;
			  }else{
				  FACT=1.0e15;
			  }
			  fprintf(IWR3," %.6E %.6E %.6E %.6E\n", CSOUR2_.ESRC[KEn-1], CSOUR2_.PSRC[KEn-1]*FACT/PTOT, (YAV-YERR)*FACT,(YAV+YERR)*FACT);
			  fprintf(IWR3," %.6E %.6E %.6E %.6E\n", CSOUR2_.ESRC[KEn+1-1], CSOUR2_.PSRC[KEn-1]*FACT/PTOT, (YAV-YERR)*FACT,(YAV+YERR)*FACT);
		}
		PTOT=0.0e0;
		fprintf(IWR3, " %.6E %.6E %.6E %.6E\n", CSOUR2_.ESRC[*CNT2_.NSEB-1],PTOT,PTOT,PTOT);
		fclose(IWR3);
	}

	/*Espectro de energia da fonte (resultado da simulação).
	Não sera implementado pois parte da não especificação de uma particula primaria KPARP=0*/

	//Energia e distribuições angulares de partículas emergentes.
	enangw2_(*CNTRL_.SHN);

	//Distribuição de dose

	if (*CNT6_.LDOSEM){
		dosew2_(*CNTRL_.SHN,*CNTRL_.TSIM, IWR2);
	}

	fprintf(IWR2, "\n   Last random seeds = %d , %d\n\n", *RSEED_.ISEED1, *RSEED_.ISEED2);
	fprintf(IWR2,"  ------------------------------------------------------------------------\n");


}


void enangd2_(FILE *IWR){

	//Transferir contadores parciais para contadores globais.

	for (int KP = 1; KP <= 3; KP++){
		for (int IEX =1; IEX <=2; IEX++){
			for (int KEn = 1; KEn <= *CENANG_.NE; KEn++){
				CENANG_.PDE[KEn-1][IEX-1][KP-1]=CENANG_.PDE[KEn-1][IEX-1][KP-1]+CENANG_.PDEP[KEn-1][IEX-1][KP-1];
              	CENANG_.PDE2[KEn-1][IEX-1][KP-1]=CENANG_.PDE2[KEn-1][IEX-1][KP-1]+pow(CENANG_.PDEP[KEn-1][IEX-1][KP-1],2);
              	CENANG_.PDEP[KEn-1][IEX-1][KP-1]=0.0e0;
              	CENANG_.LPDE[KEn-1][IEX-1][KP-1]=0;
			}
		}
	}

	for (int KP = 1; KP <= 3; KP++){
		for (int KTH =1; KTH <= *CENANG_.NTH; KTH++){
			for (int KPH = 1; KPH <= *CENANG_.NPH; KPH++){
				CENANG_.PDA[KPH-1][KTH-1][KP-1]=CENANG_.PDA[KPH-1][KTH-1][KP-1]+CENANG_.PDAP[KPH-1][KTH-1][KP-1];
              	CENANG_.PDA2[KPH-1][KTH-1][KP-1]=CENANG_.PDA2[KPH-1][KTH-1][KP-1]+pow(CENANG_.PDAP[KPH-1][KTH-1][KP-1],2);
              	CENANG_.PDAP[KPH-1][KTH-1][KP-1]=0.0e0;
              	CENANG_.LPDA[KPH-1][KTH-1][KP-1]=0;
			}
		}
	}

	//Escreva parâmetros e contadores.

	fprintf(IWR, "	%d		%.15f	%.15f	%.15f	%.15f 	%d\n", *CENANG_.NE, *CENANG_.EL,*CENANG_.EU,*CENANG_.BSE,*CENANG_.RBSE,*CENANG_.LLE);
	for (int I = 1; I <= 3; I++){
		for (int J =1; J <= 2; J++){
			for (int K = 1; K <= NE; K++){
				fprintf(IWR,"	%.16f    ", CENANG_.PDE[K-1][J-1][I-1]);
			}
		}
	}
	fprintf(IWR,"\n");

	for (int I = 1; I <= 3; I++){
		for (int J =1; J <= 2; J++){
			for (int K = 1; K <= *CENANG_.NE; K++){
				fprintf(IWR,"	%.16f    ", CENANG_.PDE2[K-1][J-1][I-1]);
			}
		}
	}
	fprintf(IWR,"\n");

	fprintf(IWR, "	%d		%.15f	%.15f	%.15f	%.15f 	%d    %.15f	%.15f 	%d	\n", *CENANG_.NTH,*CENANG_.THL,*CENANG_.THU,*CENANG_.BSTH,*CENANG_.RBSTH,*CENANG_.NPH,*CENANG_.BSPH,*CENANG_.RBSPH,*CENANG_.LLTH);
	
	for (int I = 1; I <= 3; I++){
		for (int J =1; J <= *CENANG_.NTH; J++){
			for (int K = 1; K <= *CENANG_.NPH; K++){
				fprintf(IWR,"	%.16f    ", CENANG_.PDA[K-1][J-1][I-1]);
			}
		}
	}
	fprintf(IWR,"\n");

	for (int I = 1; I <= 3; I++){
		for (int J =1; J <= *CENANG_.NTH; J++){
			for (int K = 1; K <= *CENANG_.NPH; K++){
				fprintf(IWR,"	%.16f    ", CENANG_.PDA2[K-1][J-1][I-1]);
			}
		}
	}
	fprintf(IWR,"\n");


}

void imdetd2_(FILE *IWR){
	//Transferir contadores parciais para contadores globais.
	for (int I = 1; I <= *CIMDET_.NID; I++){
		CIMDET_.EDEP[I-1]=CIMDET_.EDEP[I-1]+CIMDET_.EDEPP[I-1];
        CIMDET_.EDEP2[I-1]=CIMDET_.EDEP2[I-1]+pow(CIMDET_.EDEPP[I-1],2);
        CIMDET_.EDEPP[I-1]=0.0e0;
        CIMDET_.LEDEP[I-1]=0;
		for (int J = 1; J <= CIMDET_.NE[I-1]; J++){
			CIMDET_.DIT[J-1][I-1]=CIMDET_.DIT[J-1][I-1]+CIMDET_.DITP[J-1][I-1];
            CIMDET_.DIT2[J-1][I-1]=CIMDET_.DIT2[J-1][I-1]+pow(CIMDET_.DITP[J-1][I-1],2);
            CIMDET_.DITP[J-1][I-1]=0.0e0;
            CIMDET_.LDIT[J-1][I-1]=0;
			for (int K = 1; K <= 3; K++){
				CIMDET_.DIP[K-1][J-1][I-1]=CIMDET_.DIP[K-1][J-1][I-1]+CIMDET_.DIPP[K-1][J-1][I-1];
              	CIMDET_.DIP2[K-1][J-1][I-1]=CIMDET_.DIP2[K-1][J-1][I-1]+pow(CIMDET_.DIPP[K-1][J-1][I-1],2);
              	CIMDET_.DIPP[K-1][J-1][I-1]=0.0e0;
             	CIMDET_.LDIP[K-1][J-1][I-1]=0;
			}
		}
	}

	for (int I =1; I <= *CIMDET_.NID; I++){
		for (int J = 1; J <= CIMDET_.NE[I-1]; J++){
			CIMDET_.FLT[J-1][I-1]=CIMDET_.FLT[J-1][I-1]+CIMDET_.FLTP[J-1][I-1];
            CIMDET_.FLT2[J-1][I-1]=CIMDET_.FLT2[J-1][I-1]+pow(CIMDET_.FLTP[J-1][I-1],2);
            CIMDET_.FLTP[J-1][I-1]=0.0e0;
            CIMDET_.LFLT[J-1][I-1]=0;
			for (int K = 1; K <= 3; K++){
				CIMDET_.FLP[K-1][J-1][I-1]=CIMDET_.FLP[K-1][J-1][I-1]+CIMDET_.FLPP[K-1][J-1][I-1];
              	CIMDET_.FLP2[K-1][J-1][I-1]=CIMDET_.FLP2[K-1][J-1][I-1]+pow(CIMDET_.FLPP[K-1][J-1][I-1],2);
              	CIMDET_.FLPP[K-1][J-1][I-1]=0.0e0;
              	CIMDET_.LFLP[K-1][J-1][I-1]=0;
			}
		}
	}

	for (int I = 1; I <= *CIMDET_.NID; I++){
		for (int J = 1; J <= CIMDET_.NAGE[I-1]; J++){
			CIMDET_.AGE[J-1][I-1]=CIMDET_.AGE[J-1][I-1]+CIMDET_.AGEP[J-1][I-1];
            CIMDET_.AGE2[J-1][I-1]=CIMDET_.AGE2[J-1][I-1]+pow(CIMDET_.AGEP[J-1][I-1],2);
            CIMDET_.AGEP[J-1][I-1]=0.0e0;
            CIMDET_.LAGEA[J-1][I-1]=0;
		}
	}

	//Escrevendo paramentros e contadores

	fprintf(IWR, "	%d	\n", *CIMDET_.NID);
	for (int I = 1; I <= *CIMDET_.NID; I++){
		fprintf(IWR, "	%.16f	", CIMDET_.EDEP[I-1]);
	}
	fprintf(IWR, "\n");

	for (int I = 1; I <= *CIMDET_.NID; I++){
		fprintf(IWR, "	%.16f	", CIMDET_.EDEP2[I-1]);
	}
	fprintf(IWR, "\n");

	for (int I = 1; I <= *CIMDET_.NID; I++){
		fprintf(IWR, "	%d	", CIMDET_.IDCUT[I-1]);
	}
	fprintf(IWR, "\n");

	for (int I = 1; I <= *CIMDET_.NID; I++){
		fprintf(IWR, "	%s	\n", SPCDIO[I-1]);
		fprintf(IWR, "	%d	%.16f	%.16f	%.16f	%.16f	%.d\n", CIMDET_.NE[I-1],CIMDET_.EL[I-1],CIMDET_.EU[I-1],CIMDET_.BSE[I-1],CIMDET_.RBSE[I-1],CIMDET_.LLE[I-1]);
		for (int J=1; J <= CIMDET_.NE[I-1]; J++){
			fprintf(IWR, "	%.16f	", CIMDET_.DIT[I-1][J-1]);
		}
		fprintf(IWR,"\n");

		for (int J=1; J <= CIMDET_.NE[I-1]; J++){
			fprintf(IWR, "	%.16f	", CIMDET_.DIT2[I-1][J-1]);
		}
		fprintf(IWR,"\n");

		for (int J=1; J <= CIMDET_.NE[I-1]; J++){
			for (int K = 1; K <= 3; K++){
				fprintf(IWR, "	%.16f	", CIMDET_.DIP[K-1][I-1][J-1]);
			}
		}
		fprintf(IWR,"\n");

		for (int J=1; J <= CIMDET_.NE[I-1]; J++){
			for (int K = 1; K <= 3; K++){
				fprintf(IWR, "	%.16f	", CIMDET_.DIP2[K-1][I-1][J-1]);
			}
		}
		fprintf(IWR,"\n");

		if (CIMDET_.IDCUT[I-1] == 2){
			fprintf(IWR, " 	%s	\n", SPCFLO[I-1]);
			for (int J=1; J <= CIMDET_.NE[I-1]; J++){
				fprintf(IWR, "	%.16f	", CIMDET_.FLT[I-1][J-1]);
			}
			fprintf(IWR,"\n");

			for (int J=1; J <= CIMDET_.NE[I-1]; J++){
				fprintf(IWR, "	%.16f	", CIMDET_.FLT2[I-1][J-1]);
			}
			fprintf(IWR,"\n");

			for (int J=1; J <= CIMDET_.NE[I-1]; J++){
				for (int K = 1; K <= 3; K++){
					fprintf(IWR, "	%.16f	", CIMDET_.FLP[K-1][I-1][J-1]);
				}
			}
			fprintf(IWR,"\n");

			for (int J=1; J <= CIMDET_.NE[I-1]; J++){
				for (int K = 1; K <= 3; K++){
					fprintf(IWR, "	%.16f	", CIMDET_.FLP2[K-1][I-1][J-1]);
				}
			}
			fprintf(IWR,"\n");
		}

		fprintf(IWR, "		%d		%.16f		%.16f		%.16f		 %.16f		 %d		\n", CIMDET_.NAGE[I-1],CIMDET_.AGEL[I-1],CIMDET_.AGEU[I-1],CIMDET_.BAGE[I-1],CIMDET_.RBAGE[I-1],CIMDET_.LLAGE[I-1]);
		if (CIMDET_.NAGE[I-1] > 0){
			fprintf(IWR, "		%s	\n",SPCAGE[I-1] );
			
			for (int J=1; J <= CIMDET_.NAGE[I-1]; J++){
				fprintf(IWR, "	%.16f	", CIMDET_.AGE[I-1][J-1]);
			}
			fprintf(IWR,"\n");

			for (int J=1; J <= CIMDET_.NAGE[I-1]; J++){
				fprintf(IWR, "	%.16f	", CIMDET_.AGE2[I-1][J-1]);
			}
			fprintf(IWR,"\n");
		}
	}

}

void endetd2_(FILE *IWR){
	fprintf(IWR," 	%d	\n", *CENDET_.NID);
	
	for (int I = 1; I <= *CENDET_.NID; I++){
			fprintf(IWR, "	%.16f	", CENDET_.EDEP[I-1]);
		}
	fprintf(IWR, "\n");

	for (int I = 1; I <= *CENDET_.NID; I++){
			fprintf(IWR, "	%.16f	", CENDET_.EDEP2[I-1]);
		}
	fprintf(IWR, "\n");

	for (int I =1; I <= *CENDET_.NID; I++){
		fprintf(IWR, "		%s	\n", SPCDEO[I-1]);
		fprintf(IWR, "	%d		%.16f	%.16f		%.16f		%.16f		%d\n", CENDET_.NE[I-1],CENDET_.EL[I-1],CENDET_.EU[I-1],CENDET_.BSE[I-1],CENDET_.RBSE[I-1],CENDET_.LLE[I-1]);
		for (int J = 1; J <= CENDET_.NE[I-1]; J++){
			fprintf(IWR,"	%.16f	", CENDET_.DET[J-1][I-1]);
		}
		fprintf(IWR,"\n");
	}
}

void dosed2_(FILE *IWR){
	//Transferir contadores parciais para contadores globais.

	for (int I3 =1; I3 <= CDOSE3_.NDB[3-1]; I3++){
		for (int I2 =1; I2 <= CDOSE3_.NDB[2-1]; I2++){
			for (int I1 =1; I1 <= CDOSE3_.NDB[1-1]; I1++){
				CDOSE1_.DOSE[I3-1][I2-1][I1-1]=CDOSE1_.DOSE[I3-1][I2-1][I1-1]+CDOSE1_.DOSEP[I3-1][I2-1][I1-1];
          		CDOSE1_.DOSE2[I3-1][I2-1][I1-1]=CDOSE1_.DOSE2[I3-1][I2-1][I1-1]+pow(CDOSE1_.DOSEP[I3-1][I2-1][I1-1],2);
          		CDOSE1_.DOSEP[I3-1][I2-1][I1-1]=0.0e0;
          		CDOSE1_.LDOSE[I3-1][I2-1][I1-1]=0;
			}
		}
		CDOSE2_.DDOSE[I3-1]=CDOSE2_.DDOSE[I3-1]+CDOSE2_.DDOSEP[I3-1];
        CDOSE2_.DDOSE2[I3-1]=CDOSE2_.DDOSE2[I3-1]+pow(CDOSE2_.DDOSEP[I3-1],2);
        CDOSE2_.DDOSEP[I3-1]=0.0e0;
        CDOSE2_.LDDOSE[I3-1]=0;
	}

	// Write parameters and counters.

	for (int I =1; I <= 3; I++){
		fprintf(IWR, "	%d	", CDOSE3_.NDB[I-1]);
	}
	fprintf(IWR,"\n");

	for (int I =1; I <= 3; I++){
		fprintf(IWR, "	%.16f	", CDOSE3_.DXL[I-1]);
	}
	fprintf(IWR,"\n");

	for (int I =1; I <= 3; I++){
		fprintf(IWR, "	%.16f	", CDOSE3_.DXU[I-1]);
	}
	fprintf(IWR,"\n");

	for (int I =1; I <= 3; I++){
		fprintf(IWR, "	%.16f	", CDOSE3_.BDOSE[I-1]);
	}
	for (int I =1; I <= 3; I++){
		fprintf(IWR, "	%.16f	", CDOSE3_.RBDOSE[I-1]);
	}
	fprintf(IWR, "\n");

	fprintf(IWR, "	%d	\n", *CDOSE1_.KDOSE);

	for (int I3 =1; I3 <= CDOSE3_.NDB[3-1]; I3++){
		for (int I2 =1; I2 <= CDOSE3_.NDB[2-1]; I2++){
			for (int I1 =1; I1 <= CDOSE3_.NDB[1-1]; I1++){
				fprintf(IWR, "	%.16f	", CDOSE4_.VMASS[I3-1][I2-1][I1-1]);
			}
		}
	}
	fprintf(IWR,"\n");

	for (int I3 =1; I3 <= CDOSE3_.NDB[3-1]; I3++){
		for (int I2 =1; I2 <= CDOSE3_.NDB[2-1]; I2++){
			for (int I1 =1; I1 <= CDOSE3_.NDB[1-1]; I1++){
				fprintf(IWR, "	%.16f	", CDOSE1_.DOSE[I3-1][I2-1][I1-1]);
			}
		}
	}
	fprintf(IWR,"\n");

	for (int I3 =1; I3 <= CDOSE3_.NDB[3-1]; I3++){
		for (int I2 =1; I2 <= CDOSE3_.NDB[2-1]; I2++){
			for (int I1 =1; I1 <= CDOSE3_.NDB[1-1]; I1++){
				fprintf(IWR, "	%.16f	", CDOSE1_.DOSE2[I3-1][I2-1][I1-1]);
			}
		}
	}
	fprintf(IWR,"\n");

	for (int I3 =1; I3 <= CDOSE3_.NDB[3-1]; I3++){
		fprintf(IWR, "	%.16f	", CDOSE2_.DDOSE[I3-1]);
	}
	fprintf(IWR, "\n");

	for (int I3 =1; I3 <= CDOSE3_.NDB[3-1]; I3++){
		fprintf(IWR, "	%.16f	", CDOSE2_.DDOSE2[I3-1]);
	}
	fprintf(IWR, "\n");

}

void enangw2_(double &SHN){


	double PI=3.1415926535897932e0;
	double TWOPI=2.0e0*PI;
	double RA2DE=180.0e0/PI;
	double DE2RA=PI/180.0e0;
	double FSAFE=1.000000001e0;
	static const int NBEM=1500; //
	static const int NBTHM=1800; //numero maximo de caixas de energia.
	static const int NBPHM=180; //numero maximo de caixas angulares

	double DF, XLOW, XUPP, XX, BINS, YERR1, YAV1, YERR2, YAV2, YERR3, YAV3, DSANG, YY, S1, S2;

	//Transferir contadores parciais para contadores globais.
	for (int KP = 1; KP <= 3; KP++){
		for (int IEX = 1; IEX <=2; IEX++){
			for (int KEn = 1; KEn <= *CENANG_.NE; KEn++){
				CENANG_.PDE[KEn-1][IEX-1][KP-1]=CENANG_.PDE[KEn-1][IEX-1][KP-1]+CENANG_.PDEP[KEn-1][IEX-1][KP-1];
              	CENANG_.PDE2[KEn-1][IEX-1][KP-1]=CENANG_.PDE2[KEn-1][IEX-1][KP-1]+pow(CENANG_.PDEP[KEn-1][IEX-1][KP-1],2);
              	CENANG_.PDEP[KEn-1][IEX-1][KP-1]=0.0e0;
              	CENANG_.LPDE[KEn-1][IEX-1][KP-1]=0;
			}
		}
	}

	for (int KP =1; KP <= 3; KP++){
		for (int KTH=1; KTH<= *CENANG_.NTH; KTH++){
			for (int KPH=1; KPH <= *CENANG_.NPH; KPH++){
				CENANG_.PDA[KPH-1][KTH-1][KP-1]=CENANG_.PDA[KPH-1][KTH-1][KP-1]+CENANG_.PDAP[KPH-1][KTH-1][KP-1];
              	CENANG_.PDA2[KPH-1][KTH-1][KP-1]=CENANG_.PDA2[KPH-1][KTH-1][KP-1]+pow(CENANG_.PDAP[KPH-1][KTH-1][KP-1],2);
              	CENANG_.PDAP[KPH-1][KTH-1][KP-1]=0.0e0;
              	CENANG_.LPDA[KPH-1][KTH-1][KP-1]=0;
			}
		}
	}

	DF=1.0e0/ SHN;

	//Distribuições de energia de partículas emergentes.
	//Partículas ascendentes.

	FILE* IWR = fopen("energy-up2.dat", "w");
	if (IWR == NULL){
		printf("Não foi possível abrir o arquivo energy-up.dat");
		exit(0);
	}

	fprintf(IWR, " #  Results from PENMAIN.\n");
	fprintf(IWR, " #  Energy distributions of upbound particles.\n");
	fprintf(IWR, " #  1st column: E (eV).\n");
	fprintf(IWR, " #  2nd and 3rd columns: PDF and STU for electrons.\n");
	fprintf(IWR, " #  4th and 5th columns: PDF and STU for photons.\n");
	fprintf(IWR, " #  6th and 7th columns: PDF and STU for positrons.\n");
	fprintf(IWR, " #    PDFs and STUs in units of 1/(eV*primary_particle).\n");

	for (int KEn = 1; KEn <= *CENANG_.NE; KEn++){
		if (*CENANG_.LLE == 1){
			XLOW=exp(*CENANG_.EL+(KEn-1)* *CENANG_.BSE);
            XUPP=exp(*CENANG_.EL+KEn* *CENANG_.BSE);
            XX=0.5e0*(XUPP+XLOW);
            BINS=XUPP-XLOW;
		}else{
			XX=*CENANG_.EL+(KEn-0.5e0)* *CENANG_.BSE;
            BINS=*CENANG_.BSE;
		}
		YERR1=3.0e0*sqrt(fabs(CENANG_.PDE2[KEn-1][1-1][1-1]-pow(CENANG_.PDE[KEn-1][1-1][1-1],2)*DF));
        YAV1=CENANG_.PDE[KEn-1][1-1][1-1]*DF/BINS;
        YERR1=YERR1*DF/BINS;
        YERR2=3.0e0*sqrt(fabs(CENANG_.PDE2[KEn-1][1-1][2-1]-pow(CENANG_.PDE[KEn-1][1-1][2-1],2)*DF));
        YAV2=CENANG_.PDE[KEn-1][1-1][2-1]*DF/BINS;
        YERR2=YERR2*DF/BINS;
        YERR3=3.0e0*sqrt(fabs(CENANG_.PDE2[KEn-1][1-1][3-1]-pow(CENANG_.PDE[KEn-1][1-1][3-1],2)*DF));
        YAV3=CENANG_.PDE[KEn-1][1-1][3-1]*DF/BINS;
        YERR3=YERR3*DF/BINS;
		fprintf(IWR,"  %.6E  %.6E  %.2E  %.6E  %.2E  %.6E  %.2E\n", XX,fmax(YAV1,1.0e-35),YERR1, fmax(YAV2,1.0e-35),YERR2,fmax(YAV3,1.0e-35),YERR3);
	}
	fclose(IWR);

	//Partículas descendentes.

	FILE* IWR2 = fopen("energy-down2.dat", "w");
	if (IWR2 == NULL){
		printf("Não foi possível abrir o arquivo energy-down2.dat");
		exit(0);
	}

	fprintf(IWR2, " #  Results from PENMAIN.\n");
	fprintf(IWR2, " #  Energy distributions of downbound particles.\n");
	fprintf(IWR2, " #  1st column: E (eV).\n");
	fprintf(IWR2, " #  2nd and 3rd columns: PDF and STU for electrons.\n");
	fprintf(IWR2, " #  4th and 5th columns: PDF and STU for photons.\n");
	fprintf(IWR2, " #  6th and 7th columns: PDF and STU for positrons.\n");
	fprintf(IWR2, " #    PDFs and STUs in units of 1/(eV*primary_particle).\n");

	for (int KEn = 1; KEn <= *CENANG_.NE; KEn++){
		if (*CENANG_.LLE == 1){
			XLOW=exp(*CENANG_.EL+(KEn-1)* *CENANG_.BSE);
            XUPP=exp(*CENANG_.EL+KEn* *CENANG_.BSE);
            XX=0.5e0*(XUPP+XLOW);
            BINS=XUPP-XLOW;
		}else{
			XX=*CENANG_.EL+(KEn-0.5e0)* *CENANG_.BSE;
            BINS=*CENANG_.BSE;
		}
		YERR1=3.0e0*sqrt(fabs(CENANG_.PDE2[KEn-1][2-1][1-1]-pow(CENANG_.PDE[KEn-1][2-1][1-1],2)*DF));
        YAV1=CENANG_.PDE[KEn-1][2-1][1-1]*DF/BINS;
        YERR1=YERR1*DF/BINS;
        YERR2=3.0e0*sqrt(fabs(CENANG_.PDE2[KEn-1][2-1][2-1]-pow(CENANG_.PDE[KEn-1][2-1][2-1],2)*DF));
        YAV2=CENANG_.PDE[KEn-1][2-1][2-1]*DF/BINS;
        YERR2=YERR2*DF/BINS;
        YERR3=3.0e0*sqrt(fabs(CENANG_.PDE2[KEn-1][2-1][3-1]-pow(CENANG_.PDE[KEn-1][2-1][3-1],2)*DF));
        YAV3=CENANG_.PDE[KEn-1][2-1][3-1]*DF/BINS;
        YERR3=YERR3*DF/BINS;
		fprintf(IWR2,"  %.6E  %.6E  %.2E  %.6E  %.2E  %.6E  %.2E\n", XX,fmax(YAV1,1.0e-35),YERR1, fmax(YAV2,1.0e-35),YERR2,fmax(YAV3,1.0e-35),YERR3);
	}
	fclose(IWR2);

	//Distribuições angulares de partículas emergentes.

	if (*CENANG_.NPH > 1){
		FILE* IWR3 = fopen("angle2.dat", "w");
		if (IWR3 == NULL){
			printf("Não foi possível abrir o arquivo angle2.dat");
			exit(0);
		}

		fprintf(IWR3, " #  Results from PENMAIN.\n");
		fprintf(IWR3, " #  Angular distributions of emerging particles.\n");
		fprintf(IWR3, " #  1st and 2nd columns: THETA and PHI (deg).\n");
		fprintf(IWR3, " #  3rd and 4th columns: PDF and STU for electrons.\n");
		fprintf(IWR3, " #  5th and 6th columns: PDF and STU for photons.\n");
		fprintf(IWR3, " #  7th and 8th columns: PDF and STU for positrons.\n");
		fprintf(IWR3, " #    PDFs and STUs in units of 1/(sr*primary_particle).\n");

		for (int KTH = 1; KTH <= *CENANG_.NTH; KTH++){
			if (*CENANG_.LLTH == 1){
				XLOW=exp(*CENANG_.THL+(KTH-1)* *CENANG_.BSTH);
				XUPP=exp(*CENANG_.THL+KTH* *CENANG_.BSTH);
			}else{
				XLOW=*CENANG_.THL+(KTH-1)* *CENANG_.BSTH;
				XUPP=*CENANG_.THL+KTH* *CENANG_.BSTH;
			}
			XX=0.5e0*(XUPP+XLOW);
			BINS=XUPP-XLOW;
			DSANG=(cos(XLOW* DE2RA)-cos(XUPP* DE2RA))*(*CENANG_.BSPH* DE2RA);
			for (int L =1; L<= *CENANG_.NPH; L++){
				YY=(L-0.5e0)* *CENANG_.BSPH;

				YERR1=3.0e0*sqrt(fabs(CENANG_.PDA2[L-1][KTH-1][1-1]-pow(CENANG_.PDA[L-1][KTH-1][1-1],2)*DF));
				YAV1=CENANG_.PDA[L-1][KTH-1][1-1]*DF/DSANG;
				YERR1=YERR1*DF/DSANG;
				YERR2=3.0e0*sqrt(fabs(CENANG_.PDA2[L-1][KTH-1][2-1]-pow(CENANG_.PDA[L-1][KTH-1][2-1],2)*DF));
				YAV2=CENANG_.PDA[L-1][KTH-1][1-1]*DF/DSANG;
				YERR2=YERR2*DF/DSANG;
				YERR3=3.0e0*sqrt(fabs(CENANG_.PDE2[L-1][KTH-1][3-1]-pow(CENANG_.PDE[L-1][KTH-1][3-1],2)*DF));
				YAV3=CENANG_.PDE[L-1][KTH-1][3-1]*DF/DSANG;
				YERR3=YERR3*DF/DSANG;
				fprintf(IWR3,"  %.6E  %.6E  %.6E  %.2E  %.6E  %.2E  %.6E  %.2E\n", XX,YY,fmax(YAV1,1.0e-35),YERR1, fmax(YAV2,1.0e-35),YERR2,fmax(YAV3,1.0e-35),YERR3);
			}
			fprintf(IWR3,"\n");
			}
		fclose(IWR3);
	}

	FILE* IWR4 = fopen("polar-angle2.dat", "w");
	if (IWR4 == NULL){
		printf("Não foi possível abrir o arquivo polar-angle2.dat");
		exit(0);
	}

	fprintf(IWR4, " #  Results from PENMAIN.\n");
	fprintf(IWR4, " #   Angular distributions of emerging particles.\n");
	fprintf(IWR4, " #  1st column: E (eV).\n");
	fprintf(IWR4, " #  2nd and 3rd columns: PDF and STU for electrons.\n");
	fprintf(IWR4, " #  4th and 5th columns: PDF and STU for photons.\n");
	fprintf(IWR4, " #  6th and 7th columns: PDF and STU for positrons.\n");
	fprintf(IWR4, " #    PDFs and STUs in units of 1/(sr*primary_particle).\n");

	for (int KTH=1; KTH <= *CENANG_.NTH; KTH++){
		if (*CENANG_.LLTH == 1){
			XLOW=exp(*CENANG_.THL+(KTH-1)* *CENANG_.BSTH);
			XUPP=exp(*CENANG_.THL+KTH* *CENANG_.BSTH);
		}else{
			XLOW=*CENANG_.THL+(KTH-1)* *CENANG_.BSTH;
			XUPP=*CENANG_.THL+KTH* *CENANG_.BSTH;
		}
		XX=0.5e0*(XUPP+XLOW);
		BINS=XUPP-XLOW;
		DSANG=(cos(XLOW* DE2RA)-cos(XUPP* DE2RA))*TWOPI;

		S1=0.0e0;
        S2=0.0e0;
		for (int L = 1; L <= *CENANG_.NPH; L++){
			S1=S1+CENANG_.PDA[L-1][KTH-1][1-1];
            S2=S2+CENANG_.PDA2[L-1][KTH-1][1-1];
		}
		YERR1=3.0e0*sqrt(fabs(S2-pow(S1,2)*DF));
        YAV1=S1*DF/DSANG;
        YERR1=YERR1*DF/DSANG;

		S1=0.0e0;
        S2=0.0e0;
		for (int L = 1; L <= *CENANG_.NPH; L++){
			S1=S1+CENANG_.PDA[L-1][KTH-1][2-1];
            S2=S2+CENANG_.PDA2[L-1][KTH-1][2-1];
		}
		YERR2=3.0e0*sqrt(fabs(S2-pow(S1,2)*DF));
        YAV2=S1*DF/DSANG;
        YERR2=YERR2*DF/DSANG;

		S1=0.0e0;
        S2=0.0e0;
		for (int L = 1; L <= *CENANG_.NPH; L++){
			S1=S1+CENANG_.PDA[L-1][KTH-1][3-1];
            S2=S2+CENANG_.PDA2[L-1][KTH-1][3-1];
		}
		YERR3=3.0e0*sqrt(fabs(S2-pow(S1,2)*DF));
        YAV3=S1*DF/DSANG;
        YERR3=YERR3*DF/DSANG;
		fprintf(IWR4,"  %.6E  %.6E  %.2E  %.6E  %.2E  %.6E  %.2E\n", XX,fmax(YAV1,1.0e-35),YERR1, fmax(YAV2,1.0e-35),YERR2,fmax(YAV3,1.0e-35),YERR3);
	}
	fclose(IWR4);

}

void dosew2_(double &SHN, double &TSIM, FILE *IWR){

	double DF, DMAX,  QAV, QAV2, QER, EFFIC, ZZ, YAV, YAV2, YERR, XX, YY, XYZ, RR, ZAV, ZAV2, ZERR;
	int I1M, I2M, I3M, I1C, I2C, I3C;

	DF=1.0e0/SHN;

	//Transferir contadores parciais para contadores globais.

	for (int I3 =1; I3 <= CDOSE3_.NDB[3-1]; I3++){
		for (int I2 =1; I2 <= CDOSE3_.NDB[2-1]; I2++){
			for (int I1 =1; I1 <= CDOSE3_.NDB[1-1]; I1++){
				CDOSE1_.DOSE[I3-1][I2-1][I1-1]=CDOSE1_.DOSE[I3-1][I2-1][I1-1]+CDOSE1_.DOSEP[I3-1][I2-1][I1-1];
          		CDOSE1_.DOSE2[I3-1][I2-1][I1-1]=CDOSE1_.DOSE2[I3-1][I2-1][I1-1]+pow(CDOSE1_.DOSEP[I3-1][I2-1][I1-1],2);
          		CDOSE1_.DOSEP[I3-1][I2-1][I1-1]=0.0e0;
          		CDOSE1_.LDOSE[I3-1][I2-1][I1-1]=0;
			}
		}
		CDOSE2_.DDOSE[I3-1]=CDOSE2_.DDOSE[I3-1]+CDOSE2_.DDOSEP[I3-1];
        CDOSE2_.DDOSE2[I3-1]=CDOSE2_.DDOSE2[I3-1]+pow(CDOSE2_.DDOSEP[I3-1],2);
        CDOSE2_.DDOSEP[I3-1]=0.0e0;
        CDOSE2_.LDDOSE[I3-1]=0;
	}

	DMAX=0.0e0;
    I1M=1;
    I2M=1;
    I3M=1;

	for (int I1 =1; I1 <= CDOSE3_.NDB[1-1]; I1++){
		for (int I2 =1; I2 <= CDOSE3_.NDB[2-1]; I2++){
			for (int I3 =1; I3 <= CDOSE3_.NDB[3-1]; I3++){
				if (CDOSE1_.DOSE[I3-1][I2-1][I1-1]*CDOSE4_.VMASS[I3-1][I2-1][I1-1] > DMAX){
					I1M=I1;
					I2M=I2;
					I3M=I3;
					DMAX=CDOSE1_.DOSE[I3-1][I2-1][I1-1]*CDOSE4_.VMASS[I3-1][I2-1][I1-1];
				}
			}
		}
	}

	QAV=CDOSE1_.DOSE[I3M-1][I2M-1][I1M-1];
    QAV2=CDOSE1_.DOSE2[I3M-1][I2M-1][I1M-1];
    QER=3.0E0*sqrt(fabs(QAV2-pow(QAV,2)*DF));
    QAV=QAV*DF*CDOSE4_.VMASS[I3M-1][I2M-1][I1M-1];
    QER=QER*DF*CDOSE4_.VMASS[I3M-1][I2M-1][I1M-1];

	if (QER > 1.0e-10*fabs(QAV)){
		EFFIC=pow(QAV,2)/(pow((QER/3.0e0),2)*TSIM);
	}else{
		 EFFIC=0.0e0;
	}

	fprintf(IWR, "\n      Maximum dose ... %.6E +- %.1E eV/g  (effic. = %.2E)\n", QAV,QER,EFFIC);

	if (*CDOSE1_.KDOSE == 1){ //Caixa

		//distribuição de dose em profundidade
		FILE* IWR2 = fopen("depth-dose2.dat", "w");
		if (IWR2 == NULL){
			printf("Nao foi possivel abrir o arquivo depth-dose2.dat");
			exit(0);
		}

		fprintf(IWR2," #  Results from PENMAIN. Depth-dose distribution.\n");
		fprintf(IWR2," #  (integrated over X and Y within the volume of the material system).\n");
		fprintf(IWR2," #  1st column: z coordinate (cm).\n");
		fprintf(IWR2," #  2nd column: depth-dose (eV/(g/cm**2)).\n");
		fprintf(IWR2," #  3rd column: statistical uncertainty (3 sigma).\n");
		fprintf(IWR2," #  NOTE: The calculated dose distribution is correct only when the\n");
		fprintf(IWR2," #         Z bins have uniform mass density.\n #\n");

		for (int I3 = 1; I3 <= CDOSE3_.NDB[3-1]; I3++){
			ZZ=CDOSE3_.DXL[3-1]+(I3-0.5e0)*CDOSE3_.BDOSE[3-1];
			YAV=CDOSE2_.DDOSE[I3-1];
			YAV2=CDOSE2_.DDOSE2[I3-1];
			YERR=3.0e0*sqrt(fabs(YAV2-pow(YAV,2)*DF));
			YAV=YAV*DF*CDOSE3_.RBDOSE[3-1];
			YERR=YERR*DF*CDOSE3_.RBDOSE[3-1];
			fprintf(IWR2, " %.6E  %.6E  %.2E\n", ZZ,fmax(YAV,1.0e-35),fmax(YERR,1.0e-35));
		}
		fclose(IWR2);

		//Mapa de dose em 3D
		FILE* IWR3 = fopen("3d-dose-map2.dat", "w");
		if (IWR3 == NULL){
			printf("Nao foi possivel abrir o arquivo 3d-dose-map2.dat");
			exit(0);
		}

		fprintf(IWR3, " #  Results from PENMAIN. 3D dose distribution.\n");
		fprintf(IWR3, " #  Dose-map box:  XL = %.6E cm,  XU = %.6E cm\n", CDOSE3_.DXL[1-1],CDOSE3_.DXU[1-1]);
		fprintf(IWR3, " #                 YL = %.6E cm,  YU = %.6E cm\n", CDOSE3_.DXL[2-1],CDOSE3_.DXU[2-1]);
		fprintf(IWR3, " #                 ZL = %.6E cm,  ZU = %.6E cm\n", CDOSE3_.DXL[3-1],CDOSE3_.DXU[3-1]);
		fprintf(IWR3, " #  Numbers of bins:     NBX = %d, NBY = %d, NBZ = %d\n #\n", CDOSE3_.NDB[1-1],CDOSE3_.NDB[2-1],CDOSE3_.NDB[3-1]);
		fprintf(IWR3, " #  columns 1 to 3: coordinates X,Y,Z of the bin  centres.\n");
		fprintf(IWR3, " #  4th column: dose (eV/g).\n");
		fprintf(IWR3, " #  5th column: statistical uncertainty (3 sigma).\n");
		fprintf(IWR3, " #  columns 6 to 8: bin indices IX,IY,IZ.\n");

		for (int I3 = 1; I3 <= CDOSE3_.NDB[3-1]; I3++){
			ZZ=CDOSE3_.DXL[3-1]+(I3-0.5e0)*CDOSE3_.BDOSE[3-1];
			for (int I1=1; I1 <= CDOSE3_.NDB[1-1]; I1++){
				XX=CDOSE3_.DXL[1-1]+(I1-0.5e0)*CDOSE3_.BDOSE[1-1];
				for (int I2=1; I2 <= CDOSE3_.NDB[2-1]; I2++){
					YY=CDOSE3_.DXL[2-1]+(I2-0.5e0)*CDOSE3_.BDOSE[2-1];
					YAV=CDOSE1_.DOSE[I3-1][I2-1][I1-1];
					YAV2=CDOSE1_.DOSE2[I3-1][I2-1][I1-1];
					YERR=3.0e0*sqrt(fabs(YAV2-pow(YAV,2)*DF));
					YAV=YAV*DF*CDOSE4_.VMASS[I3-1][I2-1][I1-1];
					YERR=YERR*DF*CDOSE4_.VMASS[I3-1][I2-1][I1-1];
					fprintf(IWR3, " %.3E  %.3E  %.3E  %.6E  %.2E %4d  %4d  %4d\n", XX,YY,ZZ,fmax(YAV,1.0e-35),fmax(YERR,1.0e-35),I1,I2,I3);
				}
				fprintf(IWR3, "	\n");
			}
			fprintf(IWR3, "	\n");
		}
		fclose(IWR3);

		//Distribuições de dose nos eixos centrais.
		I1C=(CDOSE3_.NDB[1-1]/2)+1;
        I2C=(CDOSE3_.NDB[2-1]/2)+1;
        I3C=(CDOSE3_.NDB[3-1]/2)+1;

		if (CDOSE3_.NDB[1-1] > 1){
			FILE* IWR4 = fopen("x-dose2.dat", "w");
			if (IWR4 == NULL){
				printf("Nao foi possivel abrir o arquivo x-dose2.dat");
				exit(0);
			}
			fprintf(IWR4," #  Results from PENMAIN.\n");
			fprintf(IWR4," #  Dose distribution along the central X axis.\n");
			fprintf(IWR4," #  1st column: x (cm).\n");
			fprintf(IWR4," #  2nd column: dose (eV/g).\n");
			fprintf(IWR4," #  3rd column: statistical uncertainty (3 sigma).\n");

			for (int I1 = 1; I1 <= CDOSE3_.NDB[1-1]; I1++){
				XYZ=CDOSE3_.DXL[1-1]+(I1-0.5e0)*CDOSE3_.BDOSE[1-1];
				YAV=CDOSE1_.DOSE[I3C-1][I2C-1][I1-1];
				YAV2=CDOSE1_.DOSE2[I3C-1][I2C-1][I1-1];
				YERR=3.0e0*sqrt(fabs(YAV2-pow(YAV,2)*DF));
				YAV=YAV*DF*CDOSE4_.VMASS[I3C-1][I2C-1][I1-1];
				YERR=YERR*DF*CDOSE4_.VMASS[I3C-1][I2C-1][I1-1];
				fprintf(IWR4, " %.6E  %.6E  %.2E\n", XYZ,fmax(YAV,1.0e-35),fmax(YERR,1.0e-35));
			}
			fclose(IWR4);
		}

		if (CDOSE3_.NDB[2-1] > 1){
			FILE* IWR5 = fopen("y-dose2.dat", "w");
			if (IWR5 == NULL){
				printf("Nao foi possivel abrir o arquivo y-dose2.dat");
				exit(0);
			}
			fprintf(IWR5," #  Results from PENMAIN.\n");
			fprintf(IWR5," #  Dose distribution along the central Y axis.\n");
			fprintf(IWR5," #  1st column: y (cm).\n");
			fprintf(IWR5," #  2nd column: dose (eV/g).\n");
			fprintf(IWR5," #  3rd column: statistical uncertainty (3 sigma).\n");
			fprintf(IWR5," #  NOTE: The calculated dose distribution is correct only when the\n");
			fprintf(IWR5," #         Z bins have uniform mass density.\n #\n");

			for (int I2 = 1; I2 <= CDOSE3_.NDB[2-1]; I2++){
				XYZ=CDOSE3_.DXL[2-1]+(I2-0.5e0)*CDOSE3_.BDOSE[2-1];
				YAV=CDOSE1_.DOSE[I3C-1][I2-1][I1C-1];
				YAV2=CDOSE1_.DOSE2[I3C-1][I2-1][I1C-1];
				YERR=3.0e0*sqrt(fabs(YAV2-pow(YAV,2)*DF));
				YAV=YAV*DF*CDOSE4_.VMASS[I3C-1][I2-1][I1C-1];
				YERR=YERR*DF*CDOSE4_.VMASS[I3C-1][I2-1][I1C-1];
				fprintf(IWR5, " %.6E  %.6E  %.2E\n", XYZ,fmax(YAV,1.0e-35),fmax(YERR,1.0e-35));
			}
			fclose(IWR5);
		}

		if (CDOSE3_.NDB[2-1] > 1){
			FILE* IWR6 = fopen("z-dose2.dat", "w");
			if (IWR6 == NULL){
				printf("Nao foi possivel abrir o arquivo z-dose2.dat");
				exit(0);
			}
			fprintf(IWR6," #  Results from PENMAIN.\n");
			fprintf(IWR6," #  Dose distribution along the central Z axis.\n");
			fprintf(IWR6," #  1st column: z (cm).\n");
			fprintf(IWR6," #  2nd column: dose (eV/g).\n");
			fprintf(IWR6," #  3rd column: statistical uncertainty (3 sigma).\n");

			for (int I3 = 1; I3 <= CDOSE3_.NDB[3-1]; I3++){
				XYZ=CDOSE3_.DXL[3-1]+(I3-0.5e0)*CDOSE3_.BDOSE[3-1];
				YAV=CDOSE1_.DOSE[I3-1][I2C-1][I1C-1];
				YAV2=CDOSE1_.DOSE2[I3-1][I2C-1][I1C-1];
				YERR=3.0e0*sqrt(fabs(YAV2-pow(YAV,2)*DF));
				YAV=YAV*DF*CDOSE4_.VMASS[I3-1][I2C-1][I1C-1];
				YERR=YERR*DF*CDOSE4_.VMASS[I3-1][I2C-1][I1C-1];
				fprintf(IWR6, " %.6E  %.6E  %.2E\n", XYZ,fmax(YAV,1.0e-35),fmax(YERR,1.0e-35));
			}
			fclose(IWR6);
		}

	} else if (*CDOSE1_.KDOSE == 2){ //Cilindro

	//Mapa de dose 2D
		FILE* IWR7 = fopen("2d-dose-map.dat", "w");
		if (IWR7 == NULL){
			printf("Nao foi possivel abrir o arquivo 2d-dose-map.dat");
			exit(0);
		}
		fprintf(IWR7, " #  Results from PENMAIN. Dose distribution.\n");
		fprintf(IWR7, " #  Dose-map cylinder:         RU = %.6E cm\n", CDOSE3_.DXU[1-1]);
		fprintf(IWR7, " #     ZL = %.6E cm,  ZU = %.6E cm\n", CDOSE3_.DXL[3-1], CDOSE3_.DXU[3-1]);
		fprintf(IWR7, " #  Numbers of bins:     NBR = %d, NBZ = %d\n", CDOSE3_.NDB[1-1],CDOSE3_.NDB[3-1]);
		fprintf(IWR7, " #  columns 1 and 2: coordinates R,Z of the bin  centres\n");
		fprintf(IWR7, " #  3rd column: dose (eV/g).\n");
		fprintf(IWR7, " #  4th column: statistical uncertainty (3 sigma).\n");

		for (int I1=1; I1 <= CDOSE3_.NDB[1-1]; I1++){
			RR=(I1-0.5e0)*CDOSE3_.BDOSE[1-1];
			for (int I3=1; I3<= CDOSE3_.NDB[3-1]; I3++){
				ZZ=CDOSE3_.DXL[3-1]+(I3-0.5e0)*CDOSE3_.BDOSE[3-1];
				YAV=CDOSE1_.DOSE[I3-1][1-1][I1-1];
				YAV2=CDOSE1_.DOSE2[I3-1][1-1][I1-1];
				YERR=3.0e0*sqrt(fabs(YAV2-pow(YAV,2)*DF));
				YAV=YAV*DF*CDOSE4_.VMASS[I3-1][1-1][I1-1];
				YERR=YERR*DF*CDOSE4_.VMASS[I3-1][1-1][I1-1];
				fprintf(IWR7, " %.6E  %.6E  %.6E  %.2E\n", RR,ZZ,fmax(YAV,1.0e-35),fmax(YERR,1.0e-35));
			}
			fprintf(IWR7, "	\n");
		}
		fclose(IWR7);

		FILE* IWR8 = fopen("depth-dose2.dat", "w");
		if (IWR8 == NULL){
			printf("Nao foi possivel abrir o arquivo depth-dose2.dat");
			exit(0);
		}

		fprintf(IWR8," #  Results from PENMAIN. Depth-dose distribution.\n");
		fprintf(IWR8," #  (integrated over X and Y within the volume of the material system).\n");
		fprintf(IWR8," #  1st column: z coordinate (cm).\n");
		fprintf(IWR8," #  2nd column: depth-dose (eV/(g/cm**2)).\n");
		fprintf(IWR8," #  3rd column: statistical uncertainty (3 sigma).\n");
		fprintf(IWR8," #  NOTE: The calculated dose distribution is correct only when the\n");
		fprintf(IWR8," #         Z bins have uniform mass density.\n #\n");

		for (int I3 = 1; I3 <= CDOSE3_.NDB[3-1]; I3++){
			ZZ=CDOSE3_.DXL[3-1]+(I3-0.5e0)*CDOSE3_.BDOSE[3-1];
			YAV=CDOSE2_.DDOSE[I3-1];
			YAV2=CDOSE2_.DDOSE2[I3-1];
			YERR=3.0e0*sqrt(fabs(YAV2-pow(YAV,2)*DF));
			YAV=YAV*DF*CDOSE3_.RBDOSE[3-1];
			YERR=YERR*DF*CDOSE3_.RBDOSE[3-1];
			fprintf(IWR8, " %.6E  %.6E  %.2E\n", ZZ,fmax(YAV,1.0e-35),fmax(YERR,1.0e-35));
		}
		fclose(IWR8);

		FILE* IWR9 = fopen("z-dose2.dat", "w");
			if (IWR9 == NULL){
				printf("Nao foi possivel abrir o arquivo z-dose2.dat");
				exit(0);
			}
			fprintf(IWR9," #  Results from PENMAIN.\n");
			fprintf(IWR9," #  Dose distribution along the central Z axis.\n");
			fprintf(IWR9," #  1st column: z (cm).\n");
			fprintf(IWR9," #  2nd column: dose (eV/g).\n");
			fprintf(IWR9," #  3rd column: statistical uncertainty (3 sigma).\n");

			for (int I3 = 1; I3 <= CDOSE3_.NDB[3-1]; I3++){
				XYZ=CDOSE3_.DXL[3-1]+(I3-0.5e0)*CDOSE3_.BDOSE[3-1];
				YAV=CDOSE1_.DOSE[I3-1][I2C-1][I1C-1];
				YAV2=CDOSE1_.DOSE2[I3-1][I2C-1][I1C-1];
				YERR=3.0e0*sqrt(fabs(YAV2-pow(YAV,2)*DF));
				YAV=YAV*DF*CDOSE4_.VMASS[I3-1][I2C-1][I1C-1];
				YERR=YERR*DF*CDOSE4_.VMASS[I3-1][I2C-1][I1C-1];
				fprintf(IWR9, " %.6E  %.6E  %.2E\n", XYZ,fmax(YAV,1.0e-35),fmax(YERR,1.0e-35));
			}
			fclose(IWR9);
	}else{ //Esfera
		//Distribuição radial da dose.
		FILE* IWR10 = fopen("radial-dose2.dat", "w");
		if (IWR10 == NULL){
			printf("Nao foi possivel abrir o arquivo radial-dose2.dat");
			exit(0);
		}

		fprintf(IWR10, " #  Results from PENMAIN. Dose distribution.\n");
		fprintf(IWR10, " #  Dose-map sphere:         RU = %.6E cm\n", CDOSE3_.DXU[1-1]);
		fprintf(IWR10, " #  Numbers of bins:     NBR = %d\n", CDOSE3_.NDB[1-1]);
		fprintf(IWR10, " #  column 1: radius R of the bin centres.\n");
		fprintf(IWR10, " #  2nd column: absorbed dose (eV/g).\n");
		fprintf(IWR10, " #  3rd column: statistical uncertainty (3 sigma).\n");
		fprintf(IWR10, " #  4th column: deposited energy per unit radius (eV/cm).\n");
		fprintf(IWR10, " #  5th column: statistical uncertainty (3 sigma).\n");

		for (int I1=1; I1 <= CDOSE3_.NDB[1-1]; I1++){
			RR=(I1-0.5e0)*CDOSE3_.BDOSE[1-1];
			YAV=CDOSE1_.DOSE[1-1][1-1][I1-1];
			YAV2=CDOSE1_.DOSE2[1-1][1-1][I1-1];
			YERR=3.0e0*sqrt(fabs(YAV2-pow(YAV,2)*DF));
			YAV=YAV*DF*CDOSE4_.VMASS[1-1][1-1][I1-1];
			YERR=YERR*DF*CDOSE4_.VMASS[1-1][1-1][I1-1];

			ZAV=CDOSE1_.DOSE[1-1][1-1][I1-1];
			ZAV2=CDOSE1_.DOSE2[1-1][1-1][I1-1];
			ZERR=3.0e0*sqrt(fabs(ZAV2-pow(ZAV,2)*DF));
			ZAV=ZAV*DF*CDOSE3_.RBDOSE[1-1];
			ZERR=ZERR*DF*CDOSE3_.RBDOSE[1-1];
			fprintf(IWR10," %.6E  %.6E  %.2E  %.6E  %.2E\n", RR,YAV,YERR,ZAV,ZERR);
		}
		fclose(IWR10);
	}
}


int main(){

	inicializarStructs();
	

	//Leia os arquivos de entrada e inicialize os pacotes de simulação.
	pmrdr2_();
	

	if (*CSOUR0_.JOBEND != 0)
		goto L103;

L101:;
	//Simulação de uma nova ducha e pontuação.
	shower2_();
	if (*CSOUR0_.JOBEND != 0)
		goto L102;

	timer2_(*CNTRL_.TSEC);

	//Terminar a simulação após o tempo previsto ou após completar Chuveiros DSHN.
	if ((*CNTRL_.TSEC < *CNTRL_.TSECA) && (*CNTRL_.SHN < *CNTRL_.DSHN)){
		//Escreva os resultados parciais após cada período de despejo.
		if (*CDUMP_.LDUMP){
			if (*CNTRL_.TSEC-*CNTRL_.TSECAD > *CNTRL_.DUMPP){
				*CNTRL_.TSIM=*CNTRL_.TSIM+cputim2_()-*CNTRL_.CPUT0;
                pmwrt2_(-1);
				printf("Number of simulated showers = %.6E\n", *CNTRL_.SHN);
                *CNTRL_.TSECAD=*CNTRL_.TSEC;
                *CNTRL_.CPUT0=cputim2_();
				goto L101;
			}
		}
		goto L101;
	}

L102:;
	*CNTRL_.TSIM=*CNTRL_.TSIM+cputim2_()-*CNTRL_.CPUT0;
L103:;
	pmwrt2_(1);
	printf("Number of simulated showers = %.6E\n", *CNTRL_.SHN);
	printf("*** END ***\n");

	memoryFree();

	//system("PAUSE");   
	return 0;
}

void timer2_(double &SEC){

	SEC = (double)(clock() - start) / CLOCKS_PER_SEC;

}

double cputim2_(){

	return (double)(clock() - start) / CLOCKS_PER_SEC;
}

void inicializarStructs(){


/*se o seu compilador suporta arrays de comprimento variável ou se colsé uma constante de tempo de compilação, 
você nem precisa calcular os deslocamentos por conta própria; se você usar
 int (*mat)[cols] = malloc(rows * sizeof *mat), você pode acessar os elementos via mat[i][j]


int (*mat)[col];
mat=(int (*)[col])malloc(sizeof(*mat)*row);*/


	//printf("\ninicializarStructs\n");
	//QSURF
	QSURF_.AXX = (double *) malloc(NS*sizeof(double));
	QSURF_.AXY = (double *) malloc(NS*sizeof(double));
	QSURF_.AXZ =(double *) malloc(NS*sizeof(double));
	QSURF_.AYY = (double *) malloc(NS*sizeof(double));
	QSURF_.AYZ = (double *) malloc(NS*sizeof(double));
	QSURF_.AZZ =(double *) malloc(NS*sizeof(double));
	QSURF_.AX = (double *) malloc(NS*sizeof(double));
	QSURF_.AY =(double *) malloc(NS*sizeof(double));
	QSURF_.AZ = (double *) malloc(NS*sizeof(double));
	QSURF_.A0 = (double *) malloc(NS*sizeof(double));
	QSURF_.NSURF = (int *) malloc(sizeof(int));
	QSURF_.KPLANE =(int *) malloc(NS*sizeof(int));		

	//QTREE
	QTREE_.NBODYS = (int *) malloc(sizeof(int));
	QTREE_.KMOTH = (int *) malloc(NB*sizeof(int));
	QTREE_.KDGHT = (int (*)[NB])malloc(NXG*NB*sizeof(int));
	QTREE_.KSURF = (int (*)[NB])malloc(NXG*NB*sizeof(int));
	QTREE_.KFLAG = (int (*)[NB])malloc(NXG*NB*sizeof(int));
	QTREE_.KSP = (int *) malloc(NS*sizeof(int));
	QTREE_.NWARN = (int *) malloc(sizeof(int));

	//TRACK_mod
	TRACK_mod_.E = (double *) malloc(sizeof(double));
	TRACK_mod_.X = (double *) malloc(sizeof(double));
	TRACK_mod_.Y = (double *) malloc(sizeof(double));
	TRACK_mod_.Z = (double *) malloc(sizeof(double));
	TRACK_mod_.U = (double *) malloc(sizeof(double));
	TRACK_mod_.V = (double *) malloc(sizeof(double));
	TRACK_mod_.W = (double *) malloc(sizeof(double));
	TRACK_mod_.WGHT = (double *) malloc(sizeof(double));
	TRACK_mod_.SP1 = (double *) malloc(sizeof(double));
	TRACK_mod_.SP2 = (double *) malloc(sizeof(double));
	TRACK_mod_.SP3 = (double *) malloc(sizeof(double));
	TRACK_mod_.PAGE = (double *) malloc(sizeof(double));
	TRACK_mod_.KPAR = (int *) malloc(sizeof(int));
	TRACK_mod_.IBODY = (int *) malloc(sizeof(int));
	TRACK_mod_.MAT = (int *) malloc(sizeof(int));
	TRACK_mod_.ILB = (int *) malloc(5*sizeof(int));
	TRACK_mod_.IPOL = (int *) malloc(sizeof(int));
	TRACK_mod_.LAGE = (bool *) malloc(sizeof(bool));

	*TRACK_mod_.IPOL = 0;
	*TRACK_mod_.LAGE = false;
	*TRACK_mod_.PAGE = 0.0e0;


	//PENELOPE_mod
	PENELOPE_mod_.EABS =  (double (*)[3])malloc(MAXMAT*3*sizeof(double));
	PENELOPE_mod_.C1 = (double *) malloc(MAXMAT*sizeof(double));
	PENELOPE_mod_.C2 = (double *) malloc(MAXMAT*sizeof(double));
	PENELOPE_mod_.WCC = (double *) malloc(MAXMAT*sizeof(double));
	PENELOPE_mod_.WCR = (double *) malloc(MAXMAT*sizeof(double));
	PENELOPE_mod_.DEN = (double *) malloc(MAXMAT*sizeof(double));
	PENELOPE_mod_.RDEN = (double *) malloc(MAXMAT*sizeof(double));
	PENELOPE_mod_.E0STEP = (double *) malloc(sizeof(double));
	PENELOPE_mod_.DESOFT = (double *) malloc(sizeof(double));
	PENELOPE_mod_.SSOFT = (double *) malloc(sizeof(double));
	PENELOPE_mod_.NMS = (int *) malloc(sizeof(int));
	PENELOPE_mod_.NEGP = (int *) malloc(sizeof(int));
	PENELOPE_mod_.NMAT = (int *) malloc(sizeof(int));


	for (int I = 1; I <= MAXMAT; I++){
		for (int J = 1; J <= 3; J++){
			PENELOPE_mod_.EABS[I-1][J-1] = 50.0e0;
		}
		PENELOPE_mod_.C1[I-1] = 0.01e0;
		PENELOPE_mod_.C2[I-1] = 0.01e0;
		PENELOPE_mod_.WCC[I-1] = 1.0e2;
		PENELOPE_mod_.WCR[I-1] = 1.0e2;
		PENELOPE_mod_.DEN[I-1] = 1.0e0;
		PENELOPE_mod_.RDEN[I-1] = 1.0e0;


	}
	*PENELOPE_mod_.NEGP = 200;
	*PENELOPE_mod_.NMS = 1000;




	//PENGEOM_mod
	PENGEOM_mod_.BALIAS = (char (*)[5]) malloc(NB*5*sizeof(char));
	PENGEOM_mod_.DSTOT = (double *) malloc(sizeof(double));
	PENGEOM_mod_.MATER = (int *) malloc(NB*sizeof(int));
	PENGEOM_mod_.KDET = (int *) malloc(NB*sizeof(int));
	PENGEOM_mod_.KSLAST = (int *) malloc(sizeof(int));
	PENGEOM_mod_.NBODY = (int *) malloc(sizeof(int));
	PENGEOM_mod_.LVERB = (bool *) malloc(sizeof(bool));

	*PENGEOM_mod_.LVERB = false;

	for (int I = 1; I <= NB; I++){
		PENGEOM_mod_.KDET[I-1] = 0; 
		PENGEOM_mod_.MATER[I-1] = 0;
		strcpy(PENGEOM_mod_.BALIAS[I-1], "    ");
	}

	//QBODY
	QBODY_.KBODY = (int (*)[NB]) malloc(NXG*NB*sizeof(int));
	QBODY_.KBOMO = (int *) malloc(NB*sizeof(int));

	//CECUTR
	CECUTR_.ECUTR = (double *) malloc(MAXMAT*sizeof(double));

	//CSGAWR
	CSGAWR_.ISGAW = (int *) malloc(sizeof(int));

	//CERSEC
	CERSEC_.IERSEC = (int *) malloc(sizeof(int));

	//CEGRID
	CEGRID_.EMIN = (double *) malloc(sizeof(double));
	CEGRID_.EL = (double *) malloc(sizeof(double));
	CEGRID_.EU = (double *) malloc(sizeof(double));
	CEGRID_.ET = (double *) malloc(NEGP*sizeof(double));
    CEGRID_.DLEMP = (double *) malloc(NEGP*sizeof(double));
	CEGRID_.DLEMP1 = (double *) malloc(sizeof(double));
	CEGRID_.DLFC = (double *) malloc(sizeof(double));
	CEGRID_.XEL = (double *) malloc(sizeof(double));
	CEGRID_.XE = (double *) malloc(sizeof(double));
	CEGRID_.XEK = (double *) malloc(sizeof(double));
	CEGRID_.KE = (int *) malloc(sizeof(int));

	//CESI0
	CESI0_.XESI = (double (*)[NRP]) malloc(16*NRP*sizeof(double));
	CESI0_.IESIF = (int *) malloc(99*sizeof(int));
	CESI0_.IESIL = (int *) malloc(99*sizeof(int));
	CESI0_.NSESI = (int *) malloc(99*sizeof(int));
	CESI0_.NCURE = (int *) malloc(sizeof(int));

	//CPSI0
	CPSI0_.XPSI = (double (*)[NRP]) malloc(16*NRP*sizeof(double));
	CPSI0_.IPSIF = (int *) malloc(99*sizeof(int));
	CPSI0_.IPSIL = (int *) malloc(99*sizeof(int));
	CPSI0_.NSPSI =(int *) malloc(99*sizeof(int));
	CPSI0_.NCURP = (int *) malloc(sizeof(int));

	//CADATA_
	CADATA_.ATW = (double *) malloc(99*sizeof(double));
	CADATA_.EPX = (double *) malloc(99*sizeof(double));
	CADATA_.RSCR = (double *) malloc(99*sizeof(double));
	CADATA_.ETA = (double *) malloc(99*sizeof(double));
	CADATA_.EB = (double (*)[99]) malloc(30*99*sizeof(double));
	CADATA_.ALW = (double (*)[99]) malloc(30*99*sizeof(double));
	CADATA_.CP0 = (double (*)[99]) malloc(30*99*sizeof(double));	
	CADATA_.IFI = (int (*)[99]) malloc(30*99*sizeof(int));
	CADATA_.IKS = (int (*)[99]) malloc(30*99*sizeof(int));
	CADATA_.NSHT = (int *) malloc(99*sizeof(int));
	CADATA_.LASYMB = (char (*)[2]) malloc(99*2*sizeof(char));


	
	double ATW[] = {1.0079e0,4.0026e0,6.9410e0,9.0122e0,1.0811e1,
					1.2011e1,1.4007e1,1.5999e1,1.8998e1,2.0179e1,2.2990e1,
					2.4305e1,2.6982e1,2.8086e1,3.0974e1,3.2066e1,3.5453e1,
					3.9948e1,3.9098e1,4.0078e1,4.4956e1,4.7880e1,5.0942e1,
					5.1996e1,5.4938e1,5.5847e1,5.8933e1,5.8690e1,6.3546e1,
					6.5390e1,6.9723e1,7.2610e1,7.4922e1,7.8960e1,7.9904e1,
					8.3800e1,8.5468e1,8.7620e1,8.8906e1,9.1224e1,9.2906e1,
					9.5940e1,9.7907e1,1.0107e2,1.0291e2,1.0642e2,1.0787e2,
					1.1241e2,1.1482e2,1.1871e2,1.2175e2,1.2760e2,1.2690e2,
					1.3129e2,1.3291e2,1.3733e2,1.3891e2,1.4012e2,1.4091e2,
					1.4424e2,1.4491e2,1.5036e2,1.5196e2,1.5725e2,1.5893e2,
					1.6250e2,1.6493e2,1.6726e2,1.6893e2,1.7304e2,1.7497e2,
					1.7849e2,1.8095e2,1.8385e2,1.8621e2,1.9020e2,1.9222e2,
					1.9508e2,1.9697e2,2.0059e2,2.0438e2,2.0720e2,2.0898e2,
					2.0898e2,2.0999e2,2.2202e2,2.2302e2,2.2603e2,2.2703e2,
					2.3204e2,2.3104e2,2.3803e2,2.3705e2,2.3905e2,2.4306e2,
					2.4707e2,2.4707e2,2.5108e2,2.5208e2};

	for (int I = 1; I <= 99; I++){
		CADATA_.ATW[I-1] = ATW[I-1];
	}

	double EPX[] = {19.2e0, 41.8e0, 40.0e0, 63.7e0, 76.0e0, 81.0e0,
					82.0e0, 95.0e0,115.0e0,137.0e0,149.0e0,156.0e0,166.0e0,
					173.0e0,173.0e0,180.0e0,174.0e0,188.0e0,190.0e0,191.0e0,
					216.0e0,233.0e0,245.0e0,257.0e0,272.0e0,286.0e0,297.0e0,
					311.0e0,322.0e0,330.0e0,334.0e0,350.0e0,347.0e0,348.0e0,
					343.0e0,352.0e0,363.0e0,366.0e0,379.0e0,393.0e0,417.0e0,
					424.0e0,428.0e0,441.0e0,449.0e0,470.0e0,470.0e0,469.0e0,
					488.0e0,488.0e0,487.0e0,485.0e0,491.0e0,482.0e0,488.0e0,
					491.0e0,501.0e0,523.0e0,535.0e0,546.0e0,560.0e0,574.0e0,
					580.0e0,591.0e0,614.0e0,628.0e0,650.0e0,658.0e0,674.0e0,
					684.0e0,694.0e0,705.0e0,718.0e0,727.0e0,736.0e0,746.0e0,
					757.0e0,790.0e0,790.0e0,800.0e0,810.0e0,823.0e0,823.0e0,
					830.0e0,825.0e0,794.0e0,827.0e0,826.0e0,841.0e0,847.0e0,
					878.0e0,890.0e0,902.0e0,921.0e0,934.0e0,939.0e0,952.0e0,
					966.0e0,980.0e0};

	for (int I = 1; I <= 99; I++){
		CADATA_.EPX[I-1] = EPX[I-1];
	}

	double RSCR[] = {1.2281e2,7.3167e1,6.9228e1,6.7301e1,6.4696e1,
					6.1228e1,5.7524e1,5.4033e1,5.0787e1,4.7851e1,4.6373e1,
					4.5401e1,4.4503e1,4.3815e1,4.3074e1,4.2321e1,4.1586e1,
					4.0953e1,4.0524e1,4.0256e1,3.9756e1,3.9144e1,3.8462e1,
					3.7778e1,3.7174e1,3.6663e1,3.5986e1,3.5317e1,3.4688e1,
					3.4197e1,3.3786e1,3.3422e1,3.3068e1,3.2740e1,3.2438e1,
					3.2143e1,3.1884e1,3.1622e1,3.1438e1,3.1142e1,3.0950e1,
					3.0758e1,3.0561e1,3.0285e1,3.0097e1,2.9832e1,2.9581e1,
					2.9411e1,2.9247e1,2.9085e1,2.8930e1,2.8721e1,2.8580e1,
					2.8442e1,2.8312e1,2.8139e1,2.7973e1,2.7819e1,2.7675e1,
					2.7496e1,2.7285e1,2.7093e1,2.6911e1,2.6705e1,2.6516e1,
					2.6304e1,2.6108e1,2.5929e1,2.5730e1,2.5577e1,2.5403e1,
					2.5245e1,2.5100e1,2.4941e1,2.4790e1,2.4655e1,2.4506e1,
					2.4391e1,2.4262e1,2.4145e1,2.4039e1,2.3922e1,2.3813e1,
					2.3712e1,2.3621e1,2.3523e1,2.3430e1,2.3331e1,2.3238e1,
					2.3139e1,2.3048e1,2.2967e1,2.2833e1,2.2694e1,2.2624e1,
					2.2545e1,2.2446e1,2.2358e1,2.2264e1};

	for (int I = 1; I <= 99; I++){
		CADATA_.RSCR[I-1] = RSCR[I-1];
	}

	double ETA[] = {1.1570e0,1.1690e0,1.2190e0,1.2010e0,1.1890e0,
        1.1740e0,1.1760e0,1.1690e0,1.1630e0,1.1570e0,1.1740e0,
        1.1830e0,1.1860e0,1.1840e0,1.1800e0,1.1780e0,1.1750e0,
        1.1700e0,1.1800e0,1.1870e0,1.1840e0,1.1800e0,1.1770e0,
        1.1660e0,1.1690e0,1.1660e0,1.1640e0,1.1620e0,1.1540e0,
        1.1560e0,1.1570e0,1.1580e0,1.1570e0,1.1580e0,1.1580e0,
        1.1580e0,1.1660e0,1.1730e0,1.1740e0,1.1750e0,1.1700e0,
        1.1690e0,1.1720e0,1.1690e0,1.1680e0,1.1640e0,1.1670e0,
        1.1700e0,1.1720e0,1.1740e0,1.1750e0,1.1780e0,1.1790e0,
        1.1800e0,1.1870e0,1.1940e0,1.1970e0,1.1960e0,1.1940e0,
        1.1940e0,1.1940e0,1.1940e0,1.1940e0,1.1960e0,1.1970e0,
        1.1960e0,1.1970e0,1.1970e0,1.1980e0,1.1980e0,1.2000e0,
        1.2010e0,1.2020e0,1.2040e0,1.2050e0,1.2060e0,1.2080e0,
        1.2070e0,1.2080e0,1.2120e0,1.2150e0,1.2180e0,1.2210e0,
        1.2240e0,1.2270e0,1.2300e0,1.2370e0,1.2430e0,1.2470e0,
        1.2500e0,1.2510e0,1.2520e0,1.2550e0,1.2560e0,1.2570e0,
        1.2590e0,1.2620e0,1.2620e0,1.2650e0};

	for (int I = 1; I <= 99; I++){
		CADATA_.ETA[I-1] = ETA[I-1];
	}


	//CGPH00_
	CGPH00_.EPH = (double *) malloc(NTP*sizeof(double));
	CGPH00_.XPH = (double (*)[NTP]) malloc(17*NTP*sizeof(double));
	CGPH00_.IPHF = (int *) malloc(99*sizeof(int));
	CGPH00_.IPHL = (int *) malloc(99*sizeof(int));
	CGPH00_.NPHS =(int *) malloc(99*sizeof(int));
	CGPH00_.NCUR = (int *) malloc(sizeof(int));

	//CRELAX_
	CRELAX_.P = (double *) malloc(NRX*sizeof(double));
	CRELAX_.ET = (double *) malloc(NRX*sizeof(double));
	CRELAX_.F = (double *) malloc(NRX*sizeof(double));
	CRELAX_.IS0 = (int *) malloc(NRX*sizeof(int));
	CRELAX_.IS1 =(int *) malloc(NRX*sizeof(int));
	CRELAX_.IS2 = (int *) malloc(NRX*sizeof(int));
	CRELAX_.IFIRST = (int (*)[99]) malloc(16*99*sizeof(int));
	CRELAX_.ILAST = (int (*)[99]) malloc(16*99*sizeof(int));
	CRELAX_.NCUR = (int *) malloc(sizeof(int));
	CRELAX_.KS = (int *) malloc(sizeof(int));
	CRELAX_.MODER = (int *) malloc(sizeof(int));

	//CRITAA_
	CRITAA_.XA = (double *) malloc(NM*sizeof(double));
	CRITAA_.AA = (double *) malloc(NM*sizeof(double));
	CRITAA_.BA = (double *) malloc(NM*sizeof(double));; 
	CRITAA_.FA = (double *) malloc(NM*sizeof(double));
	CRITAA_.IA = (int *) malloc(NM*sizeof(int));
	CRITAA_.NPM1A =  (int *) malloc(sizeof(int));

	//CRNDG3_
	CRNDG3_.X = (double *) malloc(NR*sizeof(double));
	CRNDG3_.A = (double *) malloc(NR*sizeof(double));
	CRNDG3_.B = (double *) malloc(NR*sizeof(double));
	CRNDG3_.F = (double *) malloc(NR*sizeof(double));
	CRNDG3_.KA =(int *) malloc(NR*sizeof(int));
	CRNDG3_.NPM1 = (int *) malloc(sizeof(int));

	//CRITA_
	CRITA_.XT = (double *) malloc(NM*sizeof(double));
	CRITA_.PAC =(double *) malloc(NM*sizeof(double));
	CRITA_.DPAC = (double *) malloc(NM*sizeof(double));
	CRITA_.A = (double *) malloc(NM*sizeof(double));
	CRITA_.B = (double *) malloc(NM*sizeof(double));
	CRITA_.IL = (int *) malloc(NM*sizeof(int));
	CRITA_.IU = (int *) malloc(NM*sizeof(int));
	CRITA_.NPM1 = (int *) malloc(sizeof(int));
	
	CRITA_.QTI = CRITA_.XT;
	CRITA_.PACI = CRITA_.PAC;
	CRITA_.DPACI = CRITA_.DPAC;
	CRITA_.AI = CRITA_.A;
	CRITA_.BI = CRITA_.B; 
	CRITA_.ITLI = CRITA_.IL;
	CRITA_.ITUI = CRITA_.IU;
	CRITA_.NPM1I = CRITA_.NPM1;
	
	CRITA_.XTI = CRITA_.XT;

	//CRITAN_
	CRITAN_.CNORM = (double *) malloc(sizeof(double));

	//COMPOS_
	COMPOS_.STF = (double (*)[MAXMAT]) malloc(30*MAXMAT*sizeof(double));
	COMPOS_.ZT = (double *) malloc(MAXMAT*sizeof(double));
	COMPOS_.AT = (double *) malloc(MAXMAT*sizeof(double));
	COMPOS_.RHO = (double *) malloc(MAXMAT*sizeof(double));
	COMPOS_.VMOL = (double *) malloc(MAXMAT*sizeof(double));
	COMPOS_.IZ = (int (*)[MAXMAT]) malloc(30*MAXMAT*sizeof(int));
	COMPOS_.NELEM = (int *) malloc(MAXMAT*sizeof(int));

	//CRANGE_
	CRANGE_.RANGE = (double (*)[MAXMAT][3])malloc(NEGP*MAXMAT*3*sizeof(double));
	CRANGE_.RANGEL = (double (*)[MAXMAT][3])malloc(NEGP*MAXMAT*3*sizeof(double));

	//CEIN_
	CEIN_.EXPOT = (double *) malloc(MAXMAT*sizeof(double));
	CEIN_.OP2 = (double *) malloc(MAXMAT*sizeof(double));
	CEIN_.F = (double (*)[MAXMAT]) malloc(NO*MAXMAT*sizeof(double));
	CEIN_.UI = (double (*)[MAXMAT]) malloc(NO*MAXMAT*sizeof(double));
	CEIN_.WRI = (double (*)[MAXMAT]) malloc(NO*MAXMAT*sizeof(double));
	CEIN_.KZ = (int (*)[MAXMAT]) malloc(NO*MAXMAT*sizeof(int));
	CEIN_.KS = (int (*)[MAXMAT]) malloc(NO*MAXMAT*sizeof(int));
	CEIN_.NOSC = (int *) malloc(MAXMAT*sizeof(int));

	//CEINTF_
	CEINTF_.T1EI = (double *) malloc(NEGP*sizeof(double));
	CEINTF_.T2EI = (double *) malloc(NEGP*sizeof(double));
	CEINTF_.T1PI = (double *) malloc(NEGP*sizeof(double));
	CEINTF_.T2PI = (double *) malloc(NEGP*sizeof(double));

	//CEIN00_
	CEIN00_.SEH0 = (double *) malloc(NO*sizeof(double));
	CEIN00_.SEH1 = (double *) malloc(NO*sizeof(double));
	CEIN00_.SEH2 = (double *) malloc(NO*sizeof(double));
	CEIN00_.SES0 = (double *) malloc(NO*sizeof(double));
	CEIN00_.SES1 = (double *) malloc(NO*sizeof(double));
	CEIN00_.SES2 = (double *) malloc(NO*sizeof(double));
	CEIN00_.SET0 = (double *) malloc(NO*sizeof(double));
	CEIN00_.SET1 = (double *) malloc(NO*sizeof(double));
	CEIN00_.SET2 = (double *) malloc(NO*sizeof(double));

	//CPIN00_
	CPIN00_.SPH0 = (double *) malloc(NO*sizeof(double));
	CPIN00_.SPH1 = (double *) malloc(NO*sizeof(double)); 
	CPIN00_.SPH2 = (double *) malloc(NO*sizeof(double));
	CPIN00_.SPS0 = (double *) malloc(NO*sizeof(double));
	CPIN00_.SPS1 = (double *) malloc(NO*sizeof(double));
	CPIN00_.SPS2 = (double *) malloc(NO*sizeof(double));
	CPIN00_.SPT0 = (double *) malloc(NO*sizeof(double));
	CPIN00_.SPT1 = (double *) malloc(NO*sizeof(double));
	CPIN00_.SPT2 = (double *) malloc(NO*sizeof(double));

	//CEINAC_
	CEINAC_.EINAC = (double (*)[NEGP][MAXMAT])malloc(NO*NEGP*MAXMAT*sizeof(double));
	CEINAC_.IEIN = (int (*)[MAXMAT])malloc(NO*MAXMAT*sizeof(int));
	CEINAC_.NEIN = (int *)malloc(MAXMAT*sizeof(int));

	//CESIAC_
	CESIAC_.ESIAC = (double (*)[NEGP][MAXMAT])malloc(NO*NEGP*MAXMAT*sizeof(double));
	CESIAC_.IESI = (int (*)[MAXMAT])malloc(NO*MAXMAT*sizeof(int));
	CESIAC_.NESI = (int *)malloc(MAXMAT*sizeof(int));

	//CESIN_
	CESIN_.XSEIN = (double (*)[NEGP])malloc(NO*NEGP*sizeof(double));
	CESIN_.XSESI = (double (*)[NEGP])malloc(NO*NEGP*sizeof(double));
	CESIN_.ISIE = (int *)malloc(NO*sizeof(int));

	//CPINAC_
	CPINAC_.PINAC = (double (*)[NEGP][MAXMAT])malloc(NO*NEGP*MAXMAT*sizeof(double));
	CPINAC_.IPIN = (int (*)[MAXMAT])malloc(NO*MAXMAT*sizeof(int));
	CPINAC_.NPIN = (int *)malloc(MAXMAT*sizeof(int));

	//CPSIAC_
	CPSIAC_.PSIAC = (double (*)[NEGP][MAXMAT])malloc(NO*NEGP*MAXMAT*sizeof(double));
	CPSIAC_.IPSI = (int (*)[MAXMAT])malloc(NO*MAXMAT*sizeof(int));
	CPSIAC_.NPSI = (int *)malloc(MAXMAT*sizeof(int));

	//CPSIN_
	CPSIN_.XSPIN = (double (*)[NEGP])malloc(NO*NEGP*sizeof(double));
	CPSIN_.XSPSI = (double (*)[NEGP])malloc(NO*NEGP*sizeof(double));
	CPSIN_.ISIP = (int *)malloc(NO*sizeof(int));

	//CGCO_
	CGCO_.FCO   =  (double (*)[MAXMAT])malloc(NOCO*MAXMAT*sizeof(double));
	CGCO_.UICO  =  (double (*)[MAXMAT])malloc(NOCO*MAXMAT*sizeof(double));
	CGCO_.FJ0   =  (double (*)[MAXMAT])malloc(NOCO*MAXMAT*sizeof(double));
	CGCO_.PTRSH =  (double (*)[MAXMAT])malloc(NOCO*MAXMAT*sizeof(double));
	CGCO_.KZCO  =  (int (*)[MAXMAT])malloc(NOCO*MAXMAT*sizeof(int));
	CGCO_.KSCO  =  (int (*)[MAXMAT])malloc(NOCO*MAXMAT*sizeof(int));
	CGCO_.NOSCCO = (int *)malloc(MAXMAT*sizeof(int));

	//CEIMFP_
	CEIMFP_.SEHEL = (double (*)[MAXMAT])malloc(NEGP*MAXMAT*sizeof(double));
	CEIMFP_.SEHIN = (double (*)[MAXMAT])malloc(NEGP*MAXMAT*sizeof(double));
	CEIMFP_.SEISI = (double (*)[MAXMAT])malloc(NEGP*MAXMAT*sizeof(double));
	CEIMFP_.SEHBR = (double (*)[MAXMAT])malloc(NEGP*MAXMAT*sizeof(double));
	CEIMFP_.SEAUX = (double (*)[MAXMAT])malloc(NEGP*MAXMAT*sizeof(double));
    CEIMFP_.SETOT =(double (*)[MAXMAT])malloc(NEGP*MAXMAT*sizeof(double));
	CEIMFP_.CSTPE = (double (*)[MAXMAT])malloc(NEGP*MAXMAT*sizeof(double));
	CEIMFP_.RSTPE = (double (*)[MAXMAT])malloc(NEGP*MAXMAT*sizeof(double));
	CEIMFP_.DEL =   (double (*)[MAXMAT])malloc(NEGP*MAXMAT*sizeof(double));
	CEIMFP_.W1E =   (double (*)[MAXMAT])malloc(NEGP*MAXMAT*sizeof(double));
    CEIMFP_.W2E =   (double (*)[MAXMAT])malloc(NEGP*MAXMAT*sizeof(double));
	CEIMFP_.DW1EL = (double (*)[MAXMAT])malloc(NEGP*MAXMAT*sizeof(double));
	CEIMFP_.DW2EL = (double (*)[MAXMAT])malloc(NEGP*MAXMAT*sizeof(double));
	CEIMFP_.RNDCE = (double (*)[MAXMAT])malloc(NEGP*MAXMAT*sizeof(double));
	CEIMFP_.AE =    (double (*)[MAXMAT])malloc(NEGP*MAXMAT*sizeof(double));
	CEIMFP_.BE =    (double (*)[MAXMAT])malloc(NEGP*MAXMAT*sizeof(double));
    CEIMFP_.T1E =   (double (*)[MAXMAT])malloc(NEGP*MAXMAT*sizeof(double));
	CEIMFP_.T2E =   (double (*)[MAXMAT])malloc(NEGP*MAXMAT*sizeof(double));

	//CLAS1E_
	CLAS1E_.TSTPE = (double (*)[MAXMAT])malloc(NEGP*MAXMAT*sizeof(double));
	CLAS1E_.TSTRE = (double (*)[MAXMAT])malloc(NEGP*MAXMAT*sizeof(double));
	CLAS1E_.TRL1E = (double (*)[MAXMAT])malloc(NEGP*MAXMAT*sizeof(double));
	CLAS1E_.TRL2E = (double (*)[MAXMAT])malloc(NEGP*MAXMAT*sizeof(double));

	//CPIMFP_
	CPIMFP_.SPHEL = (double (*)[MAXMAT])malloc(NEGP*MAXMAT*sizeof(double));
	CPIMFP_.SPHIN = (double (*)[MAXMAT])malloc(NEGP*MAXMAT*sizeof(double));
	CPIMFP_.SPISI = (double (*)[MAXMAT])malloc(NEGP*MAXMAT*sizeof(double));
	CPIMFP_.SPHBR = (double (*)[MAXMAT])malloc(NEGP*MAXMAT*sizeof(double));
	CPIMFP_.SPAN =  (double (*)[MAXMAT])malloc(NEGP*MAXMAT*sizeof(double));
	CPIMFP_.SPAUX = (double (*)[MAXMAT])malloc(NEGP*MAXMAT*sizeof(double));
	CPIMFP_.SPTOT = (double (*)[MAXMAT])malloc(NEGP*MAXMAT*sizeof(double));
	CPIMFP_.CSTPP = (double (*)[MAXMAT])malloc(NEGP*MAXMAT*sizeof(double));
	CPIMFP_.RSTPP = (double (*)[MAXMAT])malloc(NEGP*MAXMAT*sizeof(double));
	CPIMFP_.W1P =   (double (*)[MAXMAT])malloc(NEGP*MAXMAT*sizeof(double));
	CPIMFP_.W2P =   (double (*)[MAXMAT])malloc(NEGP*MAXMAT*sizeof(double));
	CPIMFP_.DW1PL = (double (*)[MAXMAT])malloc(NEGP*MAXMAT*sizeof(double));
	CPIMFP_.DW2PL = (double (*)[MAXMAT])malloc(NEGP*MAXMAT*sizeof(double));
	CPIMFP_.RNDCP = (double (*)[MAXMAT])malloc(NEGP*MAXMAT*sizeof(double));
	CPIMFP_.AP =    (double (*)[MAXMAT])malloc(NEGP*MAXMAT*sizeof(double));
	CPIMFP_.BP =    (double (*)[MAXMAT])malloc(NEGP*MAXMAT*sizeof(double));
	CPIMFP_.T1P =   (double (*)[MAXMAT])malloc(NEGP*MAXMAT*sizeof(double));
	CPIMFP_.T2P =   (double (*)[MAXMAT])malloc(NEGP*MAXMAT*sizeof(double));

	//CLAS1P_
	CLAS1P_.TSTPP = (double (*)[MAXMAT])malloc(NEGP*MAXMAT*sizeof(double)); 
	CLAS1P_.TSTRP = (double (*)[MAXMAT])malloc(NEGP*MAXMAT*sizeof(double)); 
	CLAS1P_.TRL1P = (double (*)[MAXMAT])malloc(NEGP*MAXMAT*sizeof(double)); 
	CLAS1P_.TRL2P = (double (*)[MAXMAT])malloc(NEGP*MAXMAT*sizeof(double));

	//CEEL00_.
	CEEL00_.EJT  = (double *)malloc(NEGP*sizeof(double)); 
	CEEL00_.XE0  = (double *)malloc(NEGP*sizeof(double));  
	CEEL00_.XE1  =(double *)malloc(NEGP*sizeof(double)); 
	CEEL00_.XE2  = (double *)malloc(NEGP*sizeof(double)); 
	CEEL00_.XP0  = (double *)malloc(NEGP*sizeof(double)); 
	CEEL00_.XP1  = (double *)malloc(NEGP*sizeof(double));  
	CEEL00_.XP2  = (double *)malloc(NEGP*sizeof(double));  
	CEEL00_.T1E0 = (double *)malloc(NEGP*sizeof(double)); 
	CEEL00_.T2E0 = (double *)malloc(NEGP*sizeof(double)); 
	CEEL00_.T1P0 = (double *)malloc(NEGP*sizeof(double)); 
	CEEL00_.T2P0 = (double *)malloc(NEGP*sizeof(double)); 
	CEEL00_.EJTL = (double *)malloc(NEGP*sizeof(double)); 
	CEEL00_.FJL  = (double *)malloc(NEGP*sizeof(double)); 
	CEEL00_.A    = (double *)malloc(NEGP*sizeof(double));    
	CEEL00_.B    = (double *)malloc(NEGP*sizeof(double));    
	CEEL00_.C    = (double *)malloc(NEGP*sizeof(double));    
	CEEL00_.D =    (double *)malloc(NEGP*sizeof(double)); 

	//CBRYLD_
	CBRYLD_.EBRY = (double (*)[MAXMAT])malloc(NEGP*MAXMAT*sizeof(double)); 
	CBRYLD_.PBRY = (double (*)[MAXMAT])malloc(NEGP*MAXMAT*sizeof(double)); 

	//CGIMFP_
	CGIMFP_.SGRA = (double (*)[MAXMAT])malloc(NEGP*MAXMAT*sizeof(double)); 
	CGIMFP_.SGCO =  (double (*)[MAXMAT])malloc(NEGP*MAXMAT*sizeof(double)); 
	CGIMFP_.SGPH =  (double (*)[MAXMAT])malloc(NEGP*MAXMAT*sizeof(double)); 
	CGIMFP_.SGPP =  (double (*)[MAXMAT])malloc(NEGP*MAXMAT*sizeof(double)); 
	CGIMFP_.SGAUX = (double (*)[MAXMAT])malloc(NEGP*MAXMAT*sizeof(double)); 

	//CGPH01_
	CGPH01_.ER = (double *)malloc(NDIM*sizeof(double)); 
	CGPH01_.XSR = (double *)malloc(NDIM*sizeof(double)); 
	CGPH01_.NPHD = (int *)malloc(sizeof(int)); 

	//CGPP01_
	CGPP01_.TRIP = (double (*)[MAXMAT])malloc(NEGP*MAXMAT*sizeof(double));

	//CEBR_
	CEBR_.WB = (double *)malloc(NBW*sizeof(double)); 
	CEBR_.PBCUT = (double (*)[MAXMAT])malloc(NEGP*MAXMAT*sizeof(double)); 
	CEBR_.WBCUT = (double (*)[MAXMAT])malloc(NEGP*MAXMAT*sizeof(double)); 
	CEBR_.PDFB = (double (*)[NEGP][MAXMAT])malloc(NBW*NEGP*MAXMAT*sizeof(double)); 
	CEBR_.DPDFB = (double (*)[NEGP][MAXMAT])malloc(NBW*NEGP*MAXMAT*sizeof(double)); 
	CEBR_.PACB = (double (*)[NEGP][MAXMAT])malloc(NBW*NEGP*MAXMAT*sizeof(double));
	CEBR_.ZBR2 = (double *)malloc(MAXMAT*sizeof(double)); 

	//CEBR01_
	CEBR01_.EBT = (double *)malloc(NBE*sizeof(double)); 
	CEBR01_.XS = (double (*)[NBE])malloc(NBE*NBW*sizeof(double)); 
	CEBR01_.TXS = (double *)malloc(NBE*sizeof(double));
	CEBR01_.X = (double *)malloc(NBE*sizeof(double));
	CEBR01_.Y = (double *)malloc(NBE*sizeof(double));

	//CEBR02_
	CEBR02_.P0 = (double (*)[NEGP][MAXMAT])malloc(NBW*NEGP*MAXMAT*sizeof(double));

	//CBRANG_
	CBRANG_.BET = (double *)malloc(6*sizeof(double)); 
	CBRANG_.BK = (double *)malloc(21*sizeof(double)); 
	CBRANG_.BP1 = (double (*)[21][6][MAXMAT])malloc(4*21*6*MAXMAT*sizeof(double));
	CBRANG_.BP2 = (double (*)[21][6][MAXMAT])malloc(4*21*6*MAXMAT*sizeof(double));
	CBRANG_.ZBEQ = (double *)malloc(MAXMAT*sizeof(double)); 

	//CEIN01_
	CEIN01_.EI = (double *)malloc(sizeof(double));
	CEIN01_.EE = (double *)malloc(sizeof(double));
	CEIN01_.CPS = (double *)malloc(sizeof(double)); 
	CEIN01_.AMOL = (double *)malloc(sizeof(double)); 
	CEIN01_.MOM = (int *)malloc(sizeof(int));

	//CSUMGA_
	CSUMGA_.IERGA = (int *)malloc(sizeof(int));
	CSUMGA_.NCALL = (int *)malloc(sizeof(int));

	//CPIN01_
	CPIN01_.EI = (double *)malloc(sizeof(double));
	CPIN01_.CPS = (double *)malloc(sizeof(double));
	CPIN01_.BHA1 = (double *)malloc(sizeof(double));
	CPIN01_.BHA2 = (double *)malloc(sizeof(double));
	CPIN01_.BHA3 = (double *)malloc(sizeof(double));
	CPIN01_.BHA4 = (double *)malloc(sizeof(double));
	CPIN01_.MOM = (int *)malloc(sizeof(int));

	//CDCSEP_
	CDCSEP_.ETS = (double *)malloc(NE*sizeof(double));
	CDCSEP_.ETL = (double *)malloc(NE*sizeof(double));
	CDCSEP_.TH = (double *)malloc(NA*sizeof(double));
	CDCSEP_.THR = (double *)malloc(NA*sizeof(double));
	CDCSEP_.XMU = (double *)malloc(NA*sizeof(double));
	CDCSEP_.XMUL = (double *)malloc(NA*sizeof(double));
	CDCSEP_.ECS = (double *)malloc(NE*sizeof(double)); 
	CDCSEP_.ETCS1 = (double *)malloc(NE*sizeof(double));
	CDCSEP_.ETCS2 = (double *)malloc(NE*sizeof(double));
	CDCSEP_.EDCS = (double (*)[NE])malloc(NA*NE*sizeof(double));
	CDCSEP_.PCS = (double *)malloc(NE*sizeof(double));
	CDCSEP_.PTCS1 = (double *)malloc(NE*sizeof(double));
	CDCSEP_.PTCS2 = (double *)malloc(NE*sizeof(double));
	CDCSEP_.PDCS = (double (*)[NE])malloc(NA*NE*sizeof(double));
	CDCSEP_.DCSI =  (double *)malloc(NA*sizeof(double));
	CDCSEP_.DCSIL = (double *)malloc(NA*sizeof(double));
	CDCSEP_.CSI = (double *)malloc(sizeof(double));
	CDCSEP_.TCS1I = (double *)malloc(sizeof(double));
	CDCSEP_.TCS2I =  (double *)malloc(sizeof(double));

	//CEELDB_
	CEELDB_.XSE = (double (*)[NEGP][NP])malloc(MAXMAT*NEGP*NP*sizeof(double));
	CEELDB_.PSE = (double (*)[NEGP][NP])malloc(MAXMAT*NEGP*NP*sizeof(double));
	CEELDB_.ASE = (double (*)[NEGP][NP])malloc(MAXMAT*NEGP*NP*sizeof(double));
	CEELDB_.BSE = (double (*)[NEGP][NP])malloc(MAXMAT*NEGP*NP*sizeof(double));
	CEELDB_.ITLE = (int (*)[NEGP][NP])malloc(MAXMAT*NEGP*NP*sizeof(int));
	CEELDB_.ITUE = (int (*)[NEGP][NP])malloc(MAXMAT*NEGP*NP*sizeof(int));

	//CPELDB_
	CPELDB_.XSP = (double (*)[NEGP][NP])malloc(MAXMAT*NEGP*NP*sizeof(double));
	CPELDB_.PSP =(double (*)[NEGP][NP])malloc(MAXMAT*NEGP*NP*sizeof(double));
	CPELDB_.ASP =(double (*)[NEGP][NP])malloc(MAXMAT*NEGP*NP*sizeof(double));
	CPELDB_.BSP = (double (*)[NEGP][NP])malloc(MAXMAT*NEGP*NP*sizeof(double));
	CPELDB_.ITLP = (int (*)[NEGP][NP])malloc(MAXMAT*NEGP*NP*sizeof(int));
	CPELDB_.ITUP = (int (*)[NEGP][NP])malloc(MAXMAT*NEGP*NP*sizeof(int));

	//CELSEP_
	CELSEP_.EELMAX = (double *)malloc(MAXMAT*sizeof(double));
	CELSEP_.PELMAX = (double *)malloc(MAXMAT*sizeof(double));
	CELSEP_.RNDCED = (double (*)[MAXMAT])malloc(NEGP*MAXMAT*sizeof(double));
	CELSEP_.RNDCPD = (double (*)[MAXMAT])malloc(NEGP*MAXMAT*sizeof(double));

	//CGRA00_
	CGRA00_.FACTE =	(double *)malloc(sizeof(double));
	CGRA00_.Q2MAX = (double *)malloc(sizeof(double));
	CGRA00_.MM = (int *)malloc(sizeof(int));
	CGRA00_.MOM = (int *)malloc(sizeof(int));

	//CGRA01_
	CGRA01_.FF =(double (*)[MAXMAT])malloc(NQ*MAXMAT*sizeof(double));
	CGRA01_.ERA = (double *)malloc(NEX*sizeof(double));
	CGRA01_.XSRA = (double (*)[MAXMAT])malloc(NEX*MAXMAT*sizeof(double));
	CGRA01_.IED = (int *)malloc(NEGP*sizeof(int));
	CGRA01_.IEU = (int *)malloc(NEGP*sizeof(int));
	CGRA01_.NE = (int *)malloc(sizeof(int));

	//CGRA02_
	CGRA02_.QQ = (double *)malloc(NQ*sizeof(double));
	CGRA02_.AR =(double (*)[MAXMAT])malloc(NQ*MAXMAT*sizeof(double));
	CGRA02_.BR = (double (*)[MAXMAT])malloc(NQ*MAXMAT*sizeof(double));
	CGRA02_.CR = (double (*)[MAXMAT])malloc(NQ*MAXMAT*sizeof(double));
	CGRA02_.DR = (double (*)[MAXMAT])malloc(NQ*MAXMAT*sizeof(double));
	CGRA02_.FF0 = (double *)malloc(MAXMAT*sizeof(double));
	CGRA02_.QQM = (double *)malloc(sizeof(double));

	//
	CGRA03_.QRA = (double (*)[NP2])malloc(MAXMAT*NP2*sizeof(double));
	CGRA03_.PRA = (double (*)[NP2])malloc(MAXMAT*NP2*sizeof(double));
	CGRA03_.DPRA = (double (*)[NP2])malloc(MAXMAT*NP2*sizeof(double));
	CGRA03_.ARA = (double (*)[NP2])malloc(MAXMAT*NP2*sizeof(double));
	CGRA03_.BRA = (double (*)[NP2])malloc(MAXMAT*NP2*sizeof(double));
	CGRA03_.PMAX = (double (*)[NEGP])malloc(MAXMAT*NEGP*sizeof(double));
	CGRA03_.ITLRA = (int (*)[NP2])malloc(MAXMAT*NP2*sizeof(int));
	CGRA03_.ITURA = (int (*)[NP2])malloc(MAXMAT*NP2*sizeof(int));

	//CGPP00_
	CGPP00_.ZEQPP = (double *)malloc(MAXMAT*sizeof(double));
	CGPP00_.F0 = (double (*)[MAXMAT])malloc(2*MAXMAT*sizeof(double));
	CGPP00_.BCB = (double *)malloc(MAXMAT*sizeof(double));

	//RSEED_
	RSEED_.ISEED1 = (int *)malloc(sizeof(int));
	RSEED_.ISEED2 = (int *)malloc(sizeof(int));

	//CTITLE_
	CTITLE_.TITLE = (char *)malloc(100*sizeof(char));
	CTITLE_.TITLE2 = (char *)malloc(100*sizeof(char));

	//CDATE_
	CDATE_.DATE23 = (char *)malloc(50*sizeof(char));

	//CSPGEO_
	CSPGEO_.DSMAX = (double *)malloc(NB*sizeof(double)); 
	CSPGEO_.EABSB = (double (*)[3])malloc(3*NB*sizeof(double));

	//CFORCI_
	CFORCI_.WLOW = (double (*)[NBV])malloc(NBV*3*sizeof(double));
	CFORCI_.WHIG = (double (*)[NBV])malloc(NBV*3*sizeof(double));
	CFORCI_.LFORCE = (bool (*)[NBV])malloc(NBV*3*sizeof(bool));

	//CXRSPL_
	CXRSPL_.IXRSPL = (int *)malloc(NBV*sizeof(int));
	CXRSPL_.ILBA = (int *)malloc(5*sizeof(int));
	CXRSPL_.LXRSPL = (bool *)malloc(NBV*sizeof(bool));

	//CSOUR0_
	CSOUR0_.CTHL =   (double *)malloc(sizeof(double)); 
    CSOUR0_.DCTH =   (double *)malloc(sizeof(double)); 
    CSOUR0_.PHIL =   (double *)malloc(sizeof(double)); 
    CSOUR0_.DPHI =   (double *)malloc(sizeof(double)); 
    CSOUR0_.KPARP =  (int *)malloc(sizeof(int)); 
    CSOUR0_.JOBEND = (int *)malloc(sizeof(int)); 
    CSOUR0_.LSCONE = (bool *)malloc(sizeof(bool));
    CSOUR0_.LGPOL =  (bool *)malloc(sizeof(bool));
    CSOUR0_.LPSF =   (bool *)malloc(sizeof(bool));

	//CSOUR1_
	CSOUR1_.E0    = (double *)malloc(sizeof(double));   
	CSOUR1_.EPMAX = (double *)malloc(sizeof(double));
	CSOUR1_.SP10  = (double *)malloc(sizeof(double));
	CSOUR1_.SP20  =(double *)malloc(sizeof(double));
	CSOUR1_.SP30  =(double *)malloc(sizeof(double));

	// CSOUR2_
	CSOUR2_.ESRC = (double *)malloc(NSEM*sizeof(double));   
	CSOUR2_.PSRC = (double *)malloc(NSEM*sizeof(double)); 
	CSOUR2_.IASRC = (int *)malloc(NSEM*sizeof(int)); 
	CSOUR2_.FSRC = (double *)malloc(NSEM*sizeof(double)); 
	CSOUR2_.LSPEC = (bool *)malloc(sizeof(bool)); 

	//CSOUR3_
	CSOUR3_.SX0 = (double *)malloc(sizeof(double)); 
	CSOUR3_.SY0 = (double *)malloc(sizeof(double)); 
	CSOUR3_.SZ0 = (double *)malloc(sizeof(double)); 
	CSOUR3_.SSX = (double *)malloc(sizeof(double)); 
	CSOUR3_.SSY = (double *)malloc(sizeof(double)); 
	CSOUR3_.SSZ = (double *)malloc(sizeof(double));
	CSOUR3_.IXSBOD = (int *)malloc(NB*sizeof(int));
	CSOUR3_.LEXSRC = (bool *)malloc(sizeof(bool));
	CSOUR3_.LEXBD = (bool *)malloc(sizeof(bool));

	//CSOUR4_
	CSOUR4_.WGMIN =   (double *)malloc(sizeof(double));
	CSOUR4_.RWGMIN =	(double *)malloc(sizeof(double));
	CSOUR4_.WGMAX =	 (double *)malloc(sizeof(double));;
	CSOUR4_.RLREAD =	(double *)malloc(sizeof(double));
	CSOUR4_.IPSFI =	(int *)malloc(sizeof(int));
	CSOUR4_.NPSF =	(int *)malloc(sizeof(int));
	CSOUR4_.NPSN =	(int *)malloc(sizeof(int));
	CSOUR4_.NSPLIT =	(int *)malloc(sizeof(int));
	CSOUR4_.KODEPS =	(int *)malloc(sizeof(int));

	//CSOUR5_
	CSOUR5_.PSFI = (char (*)[20])malloc(NPSFM*20*sizeof(char));

	//CNT0
	CNT0_.PRIM =  (double *) malloc(3*sizeof(double));
	CNT0_.PRIM2	= (double *)malloc(3*sizeof(double));
	CNT0_.DPRIM	= (double *)malloc(3*sizeof(double));
	CNT0_.SEC	= (double (*)[3])malloc(3*3*sizeof(double));
	CNT0_.SEC2	= (double (*)[3])malloc(3*3*sizeof(double));
	CNT0_.DSEC	= (double (*)[3])malloc(3*3*sizeof(double));
	CNT0_.AVW	= (double *)malloc(2*sizeof(double));
	CNT0_.AVW2	= (double *)malloc(2*sizeof(double));
	CNT0_.DAVW	= (double *)malloc(2*sizeof(double));
	CNT0_.AVA	= (double *)malloc(2*sizeof(double));
	CNT0_.AVA2	= (double *)malloc(2*sizeof(double));
	CNT0_.DAVA	= (double *)malloc(2*sizeof(double));
	CNT0_.AVE	= (double *)malloc(2*sizeof(double));
	CNT0_.AVE2	= (double *)malloc(2*sizeof(double));
	CNT0_.DAVE	= (double *)malloc(2*sizeof(double));

	//CNT1_
	CNT1_.TDEBO = (double *)malloc(NB*sizeof(double));
	CNT1_.TDEBO2 = (double *)malloc(NB*sizeof(double));
	CNT1_.DEBO = (double *)malloc(NB*sizeof(double));

	//CNT2_
	CNT2_.SHIST = (double *)malloc(NSEM*sizeof(double)); 
	CNT2_.NSEB = (int *)malloc(sizeof(int));

	//CNT3_
	CNT3_.SEDS = (double (*)[3])malloc(NSEM*3*sizeof(double)); 
	CNT3_.SEDS2 = (double (*)[3])malloc(NSEM*3*sizeof(double)); 
	CNT3_.DSDE = (double *)malloc(sizeof(double)); 
	CNT3_.RDSDE = (double *)malloc(sizeof(double)); 
	CNT3_.NSDE =  (int *)malloc(sizeof(int)); 

	//CNT4_
	CNT4_.RLAST	= (double *)malloc(sizeof(double)); 
	CNT4_.RWRITE =	(double *)malloc(sizeof(double)); 
	CNT4_.IDCUT	= (int *)malloc(NIDM*sizeof(int)); 
	CNT4_.KKDI	= (int (*)[NIDM])malloc(NIDM*3*sizeof(int));
	CNT4_.IPSF =	(int *)malloc(sizeof(int)); 
	CNT4_.NID =	(int *)malloc(NIDM*sizeof(int)); 
	CNT4_.NPSFO	 = (int *)malloc(NIDM*sizeof(int)); 
	CNT4_.IPSFO	= (int *)malloc(NIDM*sizeof(int)); 

	//CNT5_
	CNT5_.DEDE = (double *)malloc(NIDM*sizeof(double)); 
	CNT5_.KBDE = (int *)malloc(NB*sizeof(int)); 
	CNT5_.NED = (int *)malloc(sizeof(int)); 

	//CNT6_
	CNT6_.LDOSEM =(bool *)malloc(sizeof(bool)); 

	//CDUMP_
	CDUMP_.LDUMP = (bool *)malloc(sizeof(bool)); 
	CDUMP_.PFILED = (char *)malloc(50*sizeof(char)); 

	//CNTRL_
	CNTRL_.TSIM	= (double *)malloc(sizeof(double)); 
	CNTRL_.TSEC =	(double *)malloc(sizeof(double)); 
	CNTRL_.TSECA =	(double *)malloc(sizeof(double)); 
	CNTRL_.TSECAD =(double *)malloc(sizeof(double)); 
	CNTRL_.CPUT0 =	(double *)malloc(sizeof(double)); 
	CNTRL_.DUMPP =	(double *)malloc(sizeof(double)); 
	CNTRL_.DSHN =	(double *)malloc(sizeof(double)); 
	CNTRL_.SHN =	(double *)malloc(sizeof(double)); 
	CNTRL_.N =	(int *)malloc(sizeof(int)); 

	//CGCONE_
	CGCONE_.CPCT = (double *)malloc(sizeof(double)); 
	CGCONE_.CPST = (double *)malloc(sizeof(double)); 
	CGCONE_.SPCT = (double *)malloc(sizeof(double)); 
	CGCONE_.SPST = (double *)malloc(sizeof(double)); 
	CGCONE_.SPHI = (double *)malloc(sizeof(double)); 
	CGCONE_.CPHI = (double *)malloc(sizeof(double)); 
	CGCONE_.STHE = (double *)malloc(sizeof(double)); 
	CGCONE_.CTHE = (double *)malloc(sizeof(double)); 
	CGCONE_.CAPER = (double *)malloc(sizeof(double)); 

	//CENANG_
	CENANG_.EL =	(double *)malloc(sizeof(double)); 
	CENANG_.EU =(double *)malloc(sizeof(double)); 
	CENANG_.THL =	(double *)malloc(sizeof(double)); 
	CENANG_.THU	= (double *)malloc(sizeof(double)); 
	CENANG_.BSE	= (double *)malloc(sizeof(double)); 
	CENANG_.RBSE	= (double *)malloc(sizeof(double)); 
	CENANG_.BSTH	= (double *)malloc(sizeof(double)); 
	CENANG_.RBSTH	= (double *)malloc(sizeof(double)); 
	CENANG_.BSPH	= (double *)malloc(sizeof(double)); 
	CENANG_.RBSPH	= (double *)malloc(sizeof(double)); 
	CENANG_.PDE    = (double (*)[2][3])malloc(NBEM*2*3*sizeof(double)); 
	CENANG_.PDE2	= (double (*)[2][3])malloc(NBEM*2*3*sizeof(double)); 
	CENANG_.PDEP	= (double (*)[2][3])malloc(NBEM*2*3*sizeof(double)); 
	CENANG_.PDA    = (double (*)[NBTHM][3])malloc(NBPHM*NBTHM*3*sizeof(double)); 
	CENANG_.PDA2	= (double (*)[NBTHM][3])malloc(NBPHM*NBTHM*3*sizeof(double));
	CENANG_.PDAP	= (double (*)[NBTHM][3])malloc(NBPHM*NBTHM*3*sizeof(double));
	CENANG_.LPDE    = (bool (*)[2][3])malloc(NBEM*2*3*sizeof(bool)); 
	CENANG_.LPDA	= (bool (*)[NBTHM][3])malloc(NBPHM*NBTHM*3*sizeof(bool));
	CENANG_.NE	= (int *)malloc(sizeof(int)); 
	CENANG_.NTH	= (int *)malloc(sizeof(int)); 
	CENANG_.NPH	= (int *)malloc(sizeof(int)); 
	CENANG_.LLE	= (bool *)malloc(sizeof(bool)); 
	CENANG_.LLTH	= (bool *)malloc(sizeof(bool)); 

	//CIMDET_
	CIMDET_.EL =	(double *)malloc(NIDM*sizeof(double)); 	
	CIMDET_.EU =	(double *)malloc(NIDM*sizeof(double));
	CIMDET_.BSE =	(double *)malloc(NIDM*sizeof(double));
	CIMDET_.RBSE =	(double *)malloc(NIDM*sizeof(double));
	CIMDET_.ET =	(double (*)[NIDM])malloc(NIDM*(NBEM2+1)*sizeof(double)); 
	CIMDET_.EDEP =	(double *)malloc(NIDM*sizeof(double));
	CIMDET_.EDEP2 =	(double *)malloc(NIDM*sizeof(double));
	CIMDET_.EDEPP =	(double *)malloc(NIDM*sizeof(double));
	CIMDET_.DIT =	(double (*)[NIDM])malloc(NBEM2*NIDM*sizeof(double)); 
	CIMDET_.DIT2 =	(double (*)[NIDM])malloc(NBEM2*NIDM*sizeof(double));  
	CIMDET_.DITP =	(double (*)[NIDM])malloc(NBEM2*NIDM*sizeof(double));   
	CIMDET_.DIP =	(double (*)[NBEM2][NIDM])malloc(3*NBEM2*NIDM*sizeof(double)); 
	CIMDET_.DIP2 =	(double (*)[NBEM2][NIDM])malloc(3*NBEM2*NIDM*sizeof(double));  
	CIMDET_.DIPP =	(double (*)[NBEM2][NIDM])malloc(3*NBEM2*NIDM*sizeof(double));  
	CIMDET_.FLT =	(double (*)[NIDM])malloc(NBEM2*NIDM*sizeof(double)); 
	CIMDET_.FLT2 =	(double (*)[NIDM])malloc(NBEM2*NIDM*sizeof(double)); 
	CIMDET_.FLTP =	(double (*)[NIDM])malloc(NBEM2*NIDM*sizeof(double)); 
	CIMDET_.FLP =	(double (*)[NBEM2][NIDM])malloc(3*NBEM2*NIDM*sizeof(double)); 
	CIMDET_.FLP2 =	(double (*)[NBEM2][NIDM])malloc(3*NBEM2*NIDM*sizeof(double)); 
	CIMDET_.FLPP =	(double (*)[NBEM2][NIDM])malloc(3*NBEM2*NIDM*sizeof(double));  
	CIMDET_.AGEL =	(double *)malloc(NIDM*sizeof(double)); 
	CIMDET_.AGEU =	(double *)malloc(NIDM*sizeof(double));
	CIMDET_.BAGE =	(double *)malloc(NIDM*sizeof(double));
	CIMDET_.RBAGE =	(double *)malloc(NIDM*sizeof(double)); 
	CIMDET_.AGE =	(double (*)[NIDM])malloc(NBEM2*NIDM*sizeof(double)); 
	CIMDET_.AGE2 =	(double (*)[NIDM])malloc(NBEM2*NIDM*sizeof(double)); 
	CIMDET_.AGEP =	(double (*)[NIDM])malloc(NBEM2*NIDM*sizeof(double)); 
	CIMDET_.LEDEP =	(bool *)malloc(NIDM*sizeof(bool));  
	CIMDET_.LDIT =	(bool (*)[NIDM])malloc(NBEM2*NIDM*sizeof(bool));  
	CIMDET_.LDIP =	(bool (*)[NBEM2][NIDM])malloc(3*NBEM2*NIDM*sizeof(bool));  
	CIMDET_.LFLT =	(bool (*)[NIDM])malloc(NBEM2*NIDM*sizeof(bool));  
	CIMDET_.LFLP =	(bool (*)[NBEM2][NIDM])malloc(3*NBEM2*NIDM*sizeof(bool));  
	CIMDET_.LAGEA =	(bool (*)[NIDM])malloc(NBEM2*NIDM*sizeof(bool));   
	CIMDET_.IDCUT =	(int *)malloc(NIDM*sizeof(int)); 
	CIMDET_.NE =	(int *)malloc(NIDM*sizeof(int));  
	CIMDET_.LLE =	(bool *)malloc(NIDM*sizeof(bool)); 
	CIMDET_.LLAGE =	(bool *)malloc(NIDM*sizeof(bool)); 
	CIMDET_.NAGE =	(int *)malloc(NIDM*sizeof(int));  
	CIMDET_.NID =	(int *)malloc(sizeof(int)); 

	//CENDET
	CENDET_.EL = (double *)malloc(NIDM*sizeof(double)); 
	CENDET_.EU = (double *)malloc(NIDM*sizeof(double)); 
	CENDET_.BSE = (double *)malloc(NIDM*sizeof(double)); 
	CENDET_.RBSE = (double *)malloc(NIDM*sizeof(double)); 
	CENDET_.EDEP = (double *)malloc(NIDM*sizeof(double)); 
	CENDET_.EDEP2 = (double *)malloc(NIDM*sizeof(double)); 
	CENDET_.DET = (double (*)[NIDM])malloc(NBEM2*NIDM*sizeof(double)); 
	CENDET_.NE = (int *)malloc(NIDM*sizeof(int)); 
	CENDET_.NID = (int *)malloc(sizeof(int)); 
	CENDET_.LLE = (bool *)malloc(NIDM*sizeof(bool)); 

	//CDOSE1_
	CDOSE1_.DOSE = (double (*)[NDYM][NDXM])malloc(NDZM*NDYM*NDXM*sizeof(double));  
	CDOSE1_.DOSE2 = (double (*)[NDYM][NDXM])malloc(NDZM*NDYM*NDXM*sizeof(double)); 
	CDOSE1_.DOSEP = (double (*)[NDYM][NDXM])malloc(NDZM*NDYM*NDXM*sizeof(double)); 
	CDOSE1_.LDOSE = (int (*)[NDYM][NDXM])malloc(NDZM*NDYM*NDXM*sizeof(int)); 
	CDOSE1_.KDOSE = (int *)malloc(sizeof(int));

	//CDOSE2_
	CDOSE2_.DDOSE = (double *)malloc(NDZM*sizeof(double)); 
	CDOSE2_.DDOSE2 = (double *)malloc(NDZM*sizeof(double)); 
	CDOSE2_.DDOSEP = (double *)malloc(NDZM*sizeof(double)); 
	CDOSE2_.LDDOSE = (int *)malloc(NDZM*sizeof(int)); 

	//CDOSE3_
	CDOSE3_.DXL = (double *)malloc(3*sizeof(double)); 
	CDOSE3_.DXU = (double *)malloc(3*sizeof(double));
	CDOSE3_.BDOSE = (double *)malloc(3*sizeof(double));
	CDOSE3_.RBDOSE = (double *)malloc(3*sizeof(double));
	CDOSE3_.NDB = (int *)malloc(3*sizeof(int));

	//CDOSE4_
	CDOSE4_.VMASS = (double (*)[NDYM][NDXM])malloc(NDZM*NDYM*NDXM*sizeof(double)); 

	//CJUMP1_
	CJUMP1_.ELAST1 = (double *)malloc(sizeof(double)); 
	CJUMP1_.ELAST2 = (double *)malloc(sizeof(double)); 
	CJUMP1_.MHINGE = (int *)malloc(sizeof(int)); 
	CJUMP1_.KSOFTE = (int *)malloc(sizeof(int)); 
	CJUMP1_.KSOFTI = (int *)malloc(sizeof(int)); 
	CJUMP1_.KDELTA = (int *)malloc(sizeof(int)); 

	//SECST_
	SECST_.ES = (double *)malloc(NMS*sizeof(double)); 
	SECST_.XS = (double *)malloc(NMS*sizeof(double)); 
	SECST_.YS = (double *)malloc(NMS*sizeof(double)); 
	SECST_.ZS = (double *)malloc(NMS*sizeof(double)); 
	SECST_.US = (double *)malloc(NMS*sizeof(double)); 
	SECST_.VS = (double *)malloc(NMS*sizeof(double)); 
	SECST_.WS = (double *)malloc(NMS*sizeof(double)); 
	SECST_.WGHTS = (double *)malloc(NMS*sizeof(double)); 
	SECST_.SP1S = (double *)malloc(NMS*sizeof(double)); 
	SECST_.SP2S = (double *)malloc(NMS*sizeof(double)); 
	SECST_.SP3S = (double *)malloc(NMS*sizeof(double)); 
	SECST_.PAGES = (double *)malloc(NMS*sizeof(double)); 
    SECST_.KS = (int *)malloc(NMS*sizeof(int)); 
	SECST_.IBODYS =(int *)malloc(NMS*sizeof(int)); 
	SECST_.MS = (int *)malloc(NMS*sizeof(int)); 
	SECST_.ILBS = (int (*)[5])malloc(NMS*5*sizeof(int)); 
	SECST_.IPOLS = (int *)malloc(NMS*sizeof(int)); 
	SECST_.NSEC = (int *)malloc(sizeof(int)); 

	//CJUMP0_
	CJUMP0_.P = (double *)malloc(8*sizeof(double)); 
	CJUMP0_.ST = (double *)malloc(sizeof(double)); 
	CJUMP0_.DST = (double *)malloc(sizeof(double)); 
	CJUMP0_.DSR =(double *)malloc(sizeof(double)); 
	CJUMP0_.W1 = (double *)malloc(sizeof(double)); 
	CJUMP0_.W2 = (double *)malloc(sizeof(double)); 
	CJUMP0_.T1 = (double *)malloc(sizeof(double)); 
	CJUMP0_.T2 = (double *)malloc(sizeof(double)); 

	//CHIST_
	CHIST_.ILBA = (int *)malloc(5*sizeof(int)); 


}

void memoryFree(){

	free(QSURF_.AXX );
	free(QSURF_.AXY );
	free(QSURF_.AXZ );
	free(QSURF_.AYY );
	free(QSURF_.AYZ );
	free(QSURF_.AZZ );
	free(QSURF_.AX );
	free(QSURF_.AY );
	free(QSURF_.AZ );
	free(QSURF_.A0 );
	free(QSURF_.NSURF );
	free(QSURF_.KPLANE	);	

	//QTREE);
	free(QTREE_.NBODYS);
	free(QTREE_.KMOTH );
	free(QTREE_.KDGHT   );
	free(QTREE_.KSURF );
	free(QTREE_.KFLAG );
	free(QTREE_.KSP);
	free(QTREE_.NWARN );

	//TRACK_mod);
	free(TRACK_mod_.E  );
	free(TRACK_mod_.X );
	free(TRACK_mod_.Y  );
	free(TRACK_mod_.Z  );
	free(TRACK_mod_.U  );
	free(TRACK_mod_.V  );
	free(TRACK_mod_.W );
	free(TRACK_mod_.WGHT);
	free(TRACK_mod_.SP1 );
	free(TRACK_mod_.SP2 );
	free(TRACK_mod_.SP3 );
	free(TRACK_mod_.PAGE);
	free(TRACK_mod_.KPAR);
	free(TRACK_mod_.IBODY);
	free(TRACK_mod_.MAT );
	free(TRACK_mod_.ILB );
	free(TRACK_mod_.IPOL );
	free(TRACK_mod_.LAGE );

	//PENELOPE_mod);
	free(PENELOPE_mod_.EABS );
	free(PENELOPE_mod_.C1); 
	free(PENELOPE_mod_.C2); 
	free(PENELOPE_mod_.WCC);
	free(PENELOPE_mod_.WCR);
	free(PENELOPE_mod_.DEN);
	free(PENELOPE_mod_.RDEN);
	free(PENELOPE_mod_.E0STEP);
	free(PENELOPE_mod_.DESOFT);
	free(PENELOPE_mod_.SSOFT); 
	free(PENELOPE_mod_.NMS); 
	free(PENELOPE_mod_.NEGP);
	free(PENELOPE_mod_.NMAT);

	//PENGEOM_mod);
	free(PENGEOM_mod_.BALIAS);
	free(PENGEOM_mod_.DSTOT); 
	free(PENGEOM_mod_.MATER); 
	free(PENGEOM_mod_.KDET); 
	free(PENGEOM_mod_.KSLAST);
	free(PENGEOM_mod_.NBODY); 
	free(PENGEOM_mod_.LVERB); 

	//QBODY);
	free(QBODY_.KBODY );
	free(QBODY_.KBOMO );

	//CECUTR);
	free(CECUTR_.ECUTR);

	//CSGAWR);
	free(CSGAWR_.ISGAW);

	//CERSEC);
	free(CERSEC_.IERSEC);

	//CEGRID);
	free(CEGRID_.EMIN );
	free(CEGRID_.EL  );
	free(CEGRID_.EU  );
	free(CEGRID_.ET );
    free(CEGRID_.DLEMP);
	free(CEGRID_.DLEMP1);
	free(CEGRID_.DLFC );
	free(CEGRID_.XEL );
	free(CEGRID_.XE );
	free(CEGRID_.XEK  );
	free(CEGRID_.KE );

	//CESI0);
	free(CESI0_.XESI );
	free(CESI0_.IESIF );
	free(CESI0_.IESIL );
	free(CESI0_.NSESI );
	free(CESI0_.NCURE );

	//CPSI0);
	free(CPSI0_.XPSI );
	free(CPSI0_.IPSIF );
	free(CPSI0_.IPSIL );
	free(CPSI0_.NSPSI );
	free(CPSI0_.NCURP );

	//CADATA_);
	free(CADATA_.ATW );
	free(CADATA_.EPX );
	free(CADATA_.RSCR);
	free(CADATA_.ETA );
	free(CADATA_.EB);
	free(CADATA_.ALW);
	free(CADATA_.CP0 	);
	free(CADATA_.IFI );
	free(CADATA_.IKS );
	free(CADATA_.NSHT);
	free(CADATA_.LASYMB);


	//CGPH00_);
	free(CGPH00_.EPH );
	free(CGPH00_.XPH );
	free(CGPH00_.IPHF);
	free(CGPH00_.IPHL);
	free(CGPH00_.NPHS);
	free(CGPH00_.NCUR);

	//CRELAX_);
	free(CRELAX_.P );
	free(CRELAX_.ET );
	free(CRELAX_.F  );
	free(CRELAX_.IS0 );
	free(CRELAX_.IS1 );
	free(CRELAX_.IS2 );
	free(CRELAX_.IFIRST);
	free(CRELAX_.ILAST );
	free(CRELAX_.NCUR );
	free(CRELAX_.KS );
	free(CRELAX_.MODER );

	//CRITAA_);
	free(CRITAA_.XA );
	free(CRITAA_.AA );
	free(CRITAA_.BA );
	free(CRITAA_.FA );
	free(CRITAA_.IA );
	free(CRITAA_.NPM1A );

	//CRNDG3_);
	free(CRNDG3_.X );
	free(CRNDG3_.A );
	free(CRNDG3_.B );
	free(CRNDG3_.F );
	free(CRNDG3_.KA);
	free(CRNDG3_.NPM1);

	//CRITA_
	free(CRITA_.XT );
	free(CRITA_.PAC );
	free(CRITA_.DPAC);
	free(CRITA_.A  );
	free(CRITA_.B  );
	free(CRITA_.IL );
	free(CRITA_.IU );
	free(CRITA_.NPM1);
	
	/*CRITA_.QTI = CRIT);A_.XT;
	CRITA_.PACI = CRITA);_.PAC;
	CRITA_.DPACI = CRIT);A_.DPAC;
	CRITA_.AI = CRITA_.);A;
	CRITA_.BI = CRITA_.);B; 
	CRITA_.ITLI = CRITA);_.IL;
	CRITA_.ITUI = CRITA);_.IU;
	CRITA_.NPM1I = CRIT);A_.NPM1;
	);
	CRITA_.XTI = CRITA_);.XT;*/

	//CRITAN_);
	free(CRITAN_.CNORM);

	//COMPOS_);
	free(COMPOS_.STF );
	free(COMPOS_.ZT);
	free(COMPOS_.AT);
	free(COMPOS_.RHO );
	free(COMPOS_.VMOL);
	free(COMPOS_.IZ);
	free(COMPOS_.NELEM);

	//CRANGE_);
	free(CRANGE_.RANGE );
	free(CRANGE_.RANGEL);

	//CEIN_);
	free(CEIN_.EXPOT);
	free(CEIN_.OP2 );
	free(CEIN_.F );
	free(CEIN_.UI );
	free(CEIN_.WRI);
	free(CEIN_.KZ );
	free(CEIN_.KS );
	free(CEIN_.NOSC );

	//CEINTF_);
	free(CEINTF_.T1EI );
	free(CEINTF_.T2EI );
	free(CEINTF_.T1PI );
	free(CEINTF_.T2PI );

	//CEIN00_);
	free(CEIN00_.SEH0 );
	free(CEIN00_.SEH1 );
	free(CEIN00_.SEH2 );
	free(CEIN00_.SES0 );
	free(CEIN00_.SES1 );
	free(CEIN00_.SES2 );
	free(CEIN00_.SET0 );
	free(CEIN00_.SET1 );
	free(CEIN00_.SET2 );

	//CPIN00_);
	free(CPIN00_.SPH0 );
	free(CPIN00_.SPH1 );
	free(CPIN00_.SPH2 );
	free(CPIN00_.SPS0 );
	free(CPIN00_.SPS1 );
	free(CPIN00_.SPS2 );
	free(CPIN00_.SPT0 );
	free(CPIN00_.SPT1 );
	free(CPIN00_.SPT2 );

	//CEINAC_);
	free(CEINAC_.EINAC);
	free(CEINAC_.IEIN );
	free(CEINAC_.NEIN );

	//CESIAC_);
	free(CESIAC_.ESIAC);
	free(CESIAC_.IESI );
	free(CESIAC_.NESI );

	//CESIN_);
	free(CESIN_.XSEIN );
	free(CESIN_.XSESI );
	free(CESIN_.ISIE );

	//CPINAC_);
	free(CPINAC_.PINAC);
	free(CPINAC_.IPIN );
	free(CPINAC_.NPIN );

	//CPSIAC_);
	free(CPSIAC_.PSIAC );
	free(CPSIAC_.IPSI );
	free(CPSIAC_.NPSI );

	//CPSIN_);
	free(CPSIN_.XSPIN );
	free(CPSIN_.XSPSI );
	free(CPSIN_.ISIP );

	//CGCO_);
	free(CGCO_.FCO    );
	free(CGCO_.UICO  );
	free(CGCO_.FJ0   );
	free(CGCO_.PTRSH  );
	free(CGCO_.KZCO   );
	free(CGCO_.KSCO   );
	free(CGCO_.NOSCCO );

	//CEIMFP_);
	free(CEIMFP_.SEHEL );
	free(CEIMFP_.SEHIN );
	free(CEIMFP_.SEISI );
	free(CEIMFP_.SEHBR );
	free(CEIMFP_.SEAUX );
    free(CEIMFP_.SETOT );
	free(CEIMFP_.CSTPE );
	free(CEIMFP_.RSTPE );
	free(CEIMFP_.DEL  );
	free(CEIMFP_.W1E );
    free(CEIMFP_.W2E );
	free(CEIMFP_.DW1EL );
	free(CEIMFP_.DW2EL );
	free(CEIMFP_.RNDCE );
	free(CEIMFP_.AE   );
	free(CEIMFP_.BE   );
    free(CEIMFP_.T1E  );
	free(CEIMFP_.T2E );

	//CLAS1E_);
	free(CLAS1E_.TSTPE );
	free(CLAS1E_.TSTRE );
	free(CLAS1E_.TRL1E );
	free(CLAS1E_.TRL2E );

	//CPIMFP_);
	free(CPIMFP_.SPHEL );
	free(CPIMFP_.SPHIN );
	free(CPIMFP_.SPISI );
	free(CPIMFP_.SPHBR );
	free(CPIMFP_.SPAN );
	free(CPIMFP_.SPAUX );
	free(CPIMFP_.SPTOT );
	free(CPIMFP_.CSTPP );
	free(CPIMFP_.RSTPP );
	free(CPIMFP_.W1P  );
	free(CPIMFP_.W2P  );
	free(CPIMFP_.DW1PL );
	free(CPIMFP_.DW2PL );
	free(CPIMFP_.RNDCP );
	free(CPIMFP_.AP   );
	free(CPIMFP_.BP   );
	free(CPIMFP_.T1P );
	free(CPIMFP_.T2P  );

	//CLAS1P_);
	free(CLAS1P_.TSTPP );
	free(CLAS1P_.TSTRP );
	free(CLAS1P_.TRL1P );
	free(CLAS1P_.TRL2P );

	//CEEL00_.);
	free(CEEL00_.EJT  );
	free(CEEL00_.XE0  );
	free(CEEL00_.XE1  );
	free(CEEL00_.XE2  );
	free(CEEL00_.XP0  );
	free(CEEL00_.XP1  );
	free(CEEL00_.XP2  );
	free(CEEL00_.T1E0 );
	free(CEEL00_.T2E0 );
	free(CEEL00_.T1P0 );
	free(CEEL00_.T2P0 );
	free(CEEL00_.EJTL );
	free(CEEL00_.FJL  );
	free(CEEL00_.A     );
	free(CEEL00_.B     );
	free(CEEL00_.C     );
	free(CEEL00_.D  );

	//CBRYLD_);
	free(CBRYLD_.EBRY );
	free(CBRYLD_.PBRY );

	//CGIMFP_);
	free(CGIMFP_.SGRA );
	free(CGIMFP_.SGCO );
	free(CGIMFP_.SGPH );
	free(CGIMFP_.SGPP );
	free(CGIMFP_.SGAUX);

	//CGPH01_);
	free(CGPH01_.ER );
	free(CGPH01_.XSR );
	free(CGPH01_.NPHD);

	//CGPP01_);
	free(CGPP01_.TRIP );

	//CEBR_);
	free(CEBR_.WB  );
	free(CEBR_.PBCUT );
	free(CEBR_.WBCUT );
	free(CEBR_.PDFB  );
	free(CEBR_.DPDFB);
	free(CEBR_.PACB );
	free(CEBR_.ZBR2 );

	//CEBR01_);
	free(CEBR01_.EBT );
	free(CEBR01_.XS );
	free(CEBR01_.TXS );
	free(CEBR01_.X );
	free(CEBR01_.Y );

	//CEBR02_);
	free(CEBR02_.P0 );

	//CBRANG_);
	free(CBRANG_.BET );
	free(CBRANG_.BK  );
	free(CBRANG_.BP1 );
	free(CBRANG_.BP2 );
	free(CBRANG_.ZBEQ );

	//CEIN01_);
	free(CEIN01_.EI  );
	free(CEIN01_.EE  );
	free(CEIN01_.CPS );
	free(CEIN01_.AMOL );
	free(CEIN01_.MOM );

	//CSUMGA_);
	free(CSUMGA_.IERGA);
	free(CSUMGA_.NCALL);

	//CPIN01_);
	free(CPIN01_.EI );
	free(CPIN01_.CPS );
	free(CPIN01_.BHA1 );
	free(CPIN01_.BHA2 );
	free(CPIN01_.BHA3 );
	free(CPIN01_.BHA4 );
	free(CPIN01_.MOM );

	//CDCSEP_);
	free(CDCSEP_.ETS );
	free(CDCSEP_.ETL   );
	free(CDCSEP_.TH );
	free(CDCSEP_.THR );
	free(CDCSEP_.XMU );
	free(CDCSEP_.XMUL );
	free(CDCSEP_.ECS );
	free(CDCSEP_.ETCS1);
	free(CDCSEP_.ETCS2);
	free(CDCSEP_.EDCS );
	free(CDCSEP_.PCS );
	free(CDCSEP_.PTCS1);
	free(CDCSEP_.PTCS2);
	free(CDCSEP_.PDCS );
	free(CDCSEP_.DCSI );
	free(CDCSEP_.DCSIL);
	free(CDCSEP_.CSI );
	free(CDCSEP_.TCS1I);
	free(CDCSEP_.TCS2I);

	//CEELDB_);
	free(CEELDB_.XSE );
	free(CEELDB_.PSE );
	free(CEELDB_.ASE );
	free(CEELDB_.BSE );
	free(CEELDB_.ITLE );
	free(CEELDB_.ITUE );

	//CPELDB_);
	free(CPELDB_.XSP );
	free(CPELDB_.PSP );
	free(CPELDB_.ASP );
	free(CPELDB_.BSP );
	free(CPELDB_.ITLP );
	free(CPELDB_.ITUP );

	//CELSEP_);
	free(CELSEP_.EELMAX); 
	free(CELSEP_.PELMAX); 
	free(CELSEP_.RNDCED); 
	free(CELSEP_.RNDCPD); 

	//CGRA00_);
	free(CGRA00_.FACTE );
	free(CGRA00_.Q2MAX );
	free(CGRA00_.MM );
	free(CGRA00_.MOM );

	//CGRA01_);
	free(CGRA01_.FF);
	free(CGRA01_.ERA );
	free(CGRA01_.XSRA );
	free(CGRA01_.IED );
	free(CGRA01_.IEU );
	free(CGRA01_.NE );

	//CGRA02_);
	free(CGRA02_.QQ );
	free(CGRA02_.AR );
	free(CGRA02_.BR );
	free(CGRA02_.CR );
	free(CGRA02_.DR );
	free(CGRA02_.FF0 );
	free(CGRA02_.QQM );

	//);
	free(CGRA03_.QRA );
	free(CGRA03_.PRA );
	free(CGRA03_.DPRA);
	free(CGRA03_.ARA );
	free(CGRA03_.BRA );
	free(CGRA03_.PMAX);
	free(CGRA03_.ITLRA );
	free(CGRA03_.ITURA );

	//CGPP00_);
	free(CGPP00_.ZEQPP );
	free(CGPP00_.F0 );
	free(CGPP00_.BCB );

	//RSEED_);
	free(RSEED_.ISEED1 );
	free(RSEED_.ISEED2 );

	//CTITLE_);
	free(CTITLE_.TITLE );
	free(CTITLE_.TITLE2);

	//CDATE_);
	free(CDATE_.DATE23 );

	//CSPGEO_);
	free(CSPGEO_.DSMAX );
	free(CSPGEO_.EABSB );

	//CFORCI_);
	free(CFORCI_.WLOW );
	free(CFORCI_.WHIG );
	free(CFORCI_.LFORCE);

	//CXRSPL_);
	free(CXRSPL_.IXRSPL);
	free(CXRSPL_.ILBA );
	free(CXRSPL_.LXRSPL);

	//CSOUR0_);
	free(CSOUR0_.CTHL );
    free(CSOUR0_.DCTH );
    free(CSOUR0_.PHIL );
    free(CSOUR0_.DPHI );
    free(CSOUR0_.KPARP );
    free(CSOUR0_.JOBEND);
    free(CSOUR0_.LSCONE);
    free(CSOUR0_.LGPOL );
    free(CSOUR0_.LPSF );

	//CSOUR1_);
	free(CSOUR1_.E0    );
	free(CSOUR1_.EPMAX );
	free(CSOUR1_.SP10  );
	free(CSOUR1_.SP20  );
	free(CSOUR1_.SP30  );

	// CSOUR2_);
	free(CSOUR2_.ESRC  );
	free(CSOUR2_.PSRC );
	free(CSOUR2_.IASRC );
	free(CSOUR2_.FSRC );
	free(CSOUR2_.LSPEC );

	//CSOUR3_);
	free(CSOUR3_.SX0  );
	free(CSOUR3_.SY0  );
	free(CSOUR3_.SZ0  );
	free(CSOUR3_.SSX  );
	free(CSOUR3_.SSY );
	free(CSOUR3_.SSZ );
	free(CSOUR3_.IXSBOD);
	free(CSOUR3_.LEXSRC);
	free(CSOUR3_.LEXBD );

	//CSOUR4_);
	free(CSOUR4_.WGMIN );
	free(CSOUR4_.RWGMIN); 
	free(CSOUR4_.WGMAX );
	free(CSOUR4_.RLREAD); 
	free(CSOUR4_.IPSFI );
	free(CSOUR4_.NPSF );
	free(CSOUR4_.NPSN );
	free(CSOUR4_.NSPLIT); 
	free(CSOUR4_.KODEPS); 

	//CSOUR5_);
	free(CSOUR5_.PSFI  );

	//CNT0);
	free(CNT0_.PRIM );
	free(CNT0_.PRIM2);
	free(CNT0_.DPRIM);
	free(CNT0_.SEC	);
	free(CNT0_.SEC2);
	free(CNT0_.DSEC);
	free(CNT0_.AVW	);
	free(CNT0_.AVW2);
	free(CNT0_.DAVW);
	free(CNT0_.AVA	);
	free(CNT0_.AVA2);
	free(CNT0_.DAVA);
	free(CNT0_.AVE	);
	free(CNT0_.AVE2);
	free(CNT0_.DAVE);

	//CNT1_);
	free(CNT1_.TDEBO );
	free(CNT1_.TDEBO2);
	free(CNT1_.DEBO );

	//CNT2_);
	free(CNT2_.SHIST );
	free(CNT2_.NSEB );

	//CNT3_);
	free(CNT3_.SEDS );
	free(CNT3_.SEDS2  );
	free(CNT3_.DSDE );
	free(CNT3_.RDSDE );
	free(CNT3_.NSDE );

	//CNT4_);
	free(CNT4_.RLAST	);
	free(CNT4_.RWRITE);
	free(CNT4_.IDCUT	);
	free(CNT4_.KKDI);
	free(CNT4_.IPSF);
	free(CNT4_.NID 	);
	free(CNT4_.NPSFO	);
	free(CNT4_.IPSFO	);

	//CNT5_);
	free(CNT5_.DEDE );
	free(CNT5_.KBDE );
	free(CNT5_.NED );

	//CNT6_);
	free(CNT6_.LDOSEM );

	//CDUMP_);
	free(CDUMP_.LDUMP );
	free(CDUMP_.PFILED );

	//CNTRL_
	free(CNTRL_.TSIM	);
	free(CNTRL_.TSEC );
	free(CNTRL_.TSECA );
	free(CNTRL_.TSECAD );
	free(CNTRL_.CPUT0 );
	free(CNTRL_.DUMPP );
	free(CNTRL_.DSHN );
	free(CNTRL_.SHN );
	free(CNTRL_.N );

	//CGCONE_);
	free(CGCONE_.CPCT );
	free(CGCONE_.CPST );
	free(CGCONE_.SPCT );
	free(CGCONE_.SPST );
	free(CGCONE_.SPHI );
	free(CGCONE_.CPHI );
	free(CGCONE_.STHE );
	free(CGCONE_.CTHE );
	free(CGCONE_.CAPER);

	//CENANG_);
	free(CENANG_.EL );
	free(CENANG_.EU );
	free(CENANG_.THL );
	free(CENANG_.THU	);
	free(CENANG_.BSE	);
	free(CENANG_.RBSE);
	free(CENANG_.BSTH);
	free(CENANG_.RBSTH);
	free(CENANG_.BSPH);
	free(CENANG_.RBSPH);
	free(CENANG_.PDE    );
	free(CENANG_.PDE2);
	free(CENANG_.PDEP);
	free(CENANG_.PDA    );
	free(CENANG_.PDA2);
	free(CENANG_.PDAP);
	free(CENANG_.LPDE   );
	free(CENANG_.LPDA);
	free(CENANG_.NE);
	free(CENANG_.NTH	);
	free(CENANG_.NPH);
	free(CENANG_.LLE	);
	free(CENANG_.LLTH      );

	//CIMDET_
	free(CIMDET_.EL); 
	free(CIMDET_.EU); 
	free(CIMDET_.BSE); 
	free(CIMDET_.RBSE); 
	free(CIMDET_.ET); 
	free(CIMDET_.EDEP); 
	free(CIMDET_.EDEP2);
	free(CIMDET_.EDEPP); 
	free(CIMDET_.DIT); 
	free(CIMDET_.DIT2); 
	free(CIMDET_.DITP); 
	free(CIMDET_.DIP); 
	free(CIMDET_.DIP2); 
	free(CIMDET_.DIPP); 
	free(CIMDET_.FLT); 
	free(CIMDET_.FLT2); 
	free(CIMDET_.FLTP); 
	free(CIMDET_.FLP); 
	free(CIMDET_.FLP2); 
	free(CIMDET_.FLPP); 
	free(CIMDET_.AGEL); 
	free(CIMDET_.AGEU); 
	free(CIMDET_.BAGE); 
	free(CIMDET_.RBAGE); 
	free(CIMDET_.AGE); 
	free(CIMDET_.AGE2); 
	free(CIMDET_.AGEP); 
	free(CIMDET_.LEDEP); 
	free(CIMDET_.LDIT); 
	free(CIMDET_.LDIP); 
	free(CIMDET_.LFLT); 
	free(CIMDET_.LFLP); 
	free(CIMDET_.LAGEA); 
	free(CIMDET_.IDCUT); 
	free(CIMDET_.NE); 
	free(CIMDET_.LLE); 
	free(CIMDET_.LLAGE); 
	free(CIMDET_.NAGE); 
	free(CIMDET_.NID); 

	//CENDET);
	free(CENDET_.EL  );
	free(CENDET_.EU  );
	free(CENDET_.BSE );
	free(CENDET_.RBSE );
	free(CENDET_.EDEP );
	free(CENDET_.EDEP2);
	free(CENDET_.DET );
	free(CENDET_.NE  );
	free(CENDET_.NID );
	free(CENDET_.LLE );

	//CDOSE1_);
	free(CDOSE1_.DOSE );
	free(CDOSE1_.DOSE2);
	free(CDOSE1_.DOSEP);
	free(CDOSE1_.LDOSE);
	free(CDOSE1_.KDOSE);

	//CDOSE2_);
	free(CDOSE2_.DDOSE);
	free(CDOSE2_.DDOSE2 );
	free(CDOSE2_.DDOSEP );
	free(CDOSE2_.LDDOSE );

	//CDOSE3_);
	free(CDOSE3_.DXL );
	free(CDOSE3_.DXU );
	free(CDOSE3_.BDOSE    );
	free(CDOSE3_.RBDOSE );
	free(CDOSE3_.NDB );

	//CDOSE4_);
	free(CDOSE4_.VMASS );
	//CJUMP1_);
	free(CJUMP1_.ELAST1  );
	free(CJUMP1_.ELAST2 );
	free(CJUMP1_.MHINGE );
	free(CJUMP1_.KSOFTE );
	free(CJUMP1_.KSOFTI );
	free(CJUMP1_.KDELTA );

	//SECST_);
	free(SECST_.ES );
	free(SECST_.XS );
	free(SECST_.YS );
	free(SECST_.ZS );
	free(SECST_.US );
	free(SECST_.VS );
	free(SECST_.WS );
	free(SECST_.WGHTS );
	free(SECST_.SP1S );
	free(SECST_.SP2S );
	free(SECST_.SP3S );
	free(SECST_.PAGES );
    free(SECST_.KS );
	free(SECST_.IBODYS);
	free(SECST_.MS );
	free(SECST_.ILBS );
	free(SECST_.IPOLS );
	free(SECST_.NSEC );

	//CJUMP0_);
	free(CJUMP0_.P );
	free(CJUMP0_.ST );
	free(CJUMP0_.DST);
	free(CJUMP0_.DSR);
	free(CJUMP0_.W1 );
	free(CJUMP0_.W2 );
	free(CJUMP0_.T1 );
	free(CJUMP0_.T2 ); 

	//CHIST_
	free(CHIST_.ILBA); 




}



