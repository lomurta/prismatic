
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <cmath>
#include <math.h>
#include <algorithm>
#include <stdint.h>

using namespace std;
 //***

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

char LINHA[200];
char APOIO[200];




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
} CEGRID; //Rede de energia e constantes de interpola��o para a energia atual.


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
	int (*LFORCE)[NBV];

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
	double *PSFI;

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
	int *IDCUT, (*KKDI)[NDIM], *IPSF, *NID, *NPSFO, *IPSFO;

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

void transfcforci_(double (*WLOW)[NBV], double(*WHIG)[NBV], int (*LFORCE)[NBV]);

void transfcxrspl_(int *IXRSPL, int *ILBA, bool *LXRSPL);

void transfcsour0_(double *CTHL, double *DCTH, double *PHIL, double *DPHI, int *KPARP, int *JOBEND, bool *LSCONE, bool *LGPOL, bool *LPSF);

void transfcsour1_(double *E0, double *EPMAX, double *SP10, double *SP20, double *SP30);

void transfcsour2_(double *ESRC, double *PSRC, int *IASRC, double *FSRC, bool *LSPEC);

void transfcsour3_(double *SX0, double *SY0, double *SZ0, double *SSX, double *SSY, double *SSZ, int *IXSBOD, bool *LEXSRC, bool *LEXBD);

void transfcsour4_(double *WGMIN, double *RWGMIN, double *WGMAX, double *RLREAD, int *IPSFI, int *NPSF, int *NPSN, int *NSPLIT, int *KODEPS);

void transfcsour5_(double *PSFI);

void transfcnt0_(double *PRIM , double *PRIM2, double *DPRIM, double (*SEC)[3], double (*SEC2)[3], double (*DSEC)[3], 
			    double *AVW, double *AVW2, double *DAVW, double *AVA, double *AVA2, double *DAVA, double *AVE, double *AVE2, double *DAVE);

void transfcnt1_(double *TDEBO, double *TDEBO2, double *DEBO);

void transfcnt2_(double *SHIST, int *NSEB);

void transfcnt3_(double (*SEDS)[3], double (*SEDS2)[3], double *DSDE, double *RDSDE, int *NSDE);

void transfcnt4_(double *RLAST, double *RWRITE, int *IDCUT, int (*KKDI)[NDIM], int *IPSF, int *NID, int *NPSFO, int *IPSFO);

void transfcnt5_(double *DEDE, int *KBDE, int *NED);

void transfcnt6_(bool *LDOSEM);

void transfcdump_(bool *LDUMP, char *PFILED);

void transfcntrl_(double *TSIM, double *TSEC, double *TSECA, double *TSECAD, double *CPUT0, double *DUMPP, double *DSHN, double *SHN, int *N);











void einat_(double &E, double &WCCM, double &XH0, double &XH1, double &XH2, double &XS0, double &XS1, double &XS2, double &XT1, double &XT2, double &DELTA, int *M); //REVISTO


void ebrar_(double EBT, int M, int IRD, int IWR, int INFO);


void locate2_();
    
void fsurf2_(int &KS, double &A, double &B, double &C);

void step2_(double *KS, double *DSEF, int *NCROSS);

void stepsi2_(int &KB, double *S, int *IS, int &NSC);

void steplb2_(int &KB, int &IERR);

void geomin2_(double *PARINP, int *NPINP, int *NMATG, int *NBOD, int *IRD, int *IWR);

void rotshf2_(double &OMEGA, double &THETA, double &PHI, double &DX, double &DY, double &DZ, double &AXX, double &AXY, double &AXZ, double &AYY, double &AYZ, double &AZZ, double &AX, double &AY, double &AZ, double &A0);

void peinit2_(double *EMAX, int *NMATER, int *IWR, int *INFO, char (*PMFILE)[20]);

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

void pmrdr2();

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

void transfcforci_(double (*WLOW)[NBV], double(*WHIG)[NBV], int (*LFORCE)[NBV]){
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

void transfcsour5_(double *PSFI){
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

void transfcnt4_(double *RLAST, double *RWRITE, int *IDCUT, int (*KKDI)[NDIM], int *IPSF, int *NID, int *NPSFO, int *IPSFO){
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








  
  
  
  
  
      
void fsurf2_(int& KS, double& A, double& B, double& C){
	
	
	/*
	Calcula os par�metros da fun��o mestre da superf�cie KS e o raio (X, Y, Z) + S * (U, V, W).
	*/
     KS = KS - 1;
     
    if (QSURF_.KPLANE[KS] == 0) {
        A = *TRACK_mod_.U * (QSURF_.AXX[KS] * *TRACK_mod_.U + QSURF_.AXY[KS] * *TRACK_mod_.V + QSURF_.AXZ[KS] * *TRACK_mod_.W) +
            *TRACK_mod_.V * (QSURF_.AYY[KS] * *TRACK_mod_.V + QSURF_.AYZ[KS] * *TRACK_mod_.W) + *TRACK_mod_.W * QSURF_.AZZ[KS] * *TRACK_mod_.W;
		double XXX = QSURF_.AXX[KS] * *TRACK_mod_.X + QSURF_.AXY[KS] * *TRACK_mod_.Y + QSURF_.AXZ[KS] * *TRACK_mod_.Z + QSURF_.AX[KS];
        double YYY = QSURF_.AYY[KS] * *TRACK_mod_.Y + QSURF_.AYZ[KS] * *TRACK_mod_.Z + QSURF_.AY[KS];
        double ZZZ = QSURF_.AZZ[KS] * *TRACK_mod_.Z + QSURF_.AZ[KS];

        B = *TRACK_mod_.U * (QSURF_.AXX[KS] * *TRACK_mod_.X + XXX) + *TRACK_mod_.V * (QSURF_.AXY[KS] * *TRACK_mod_.X + QSURF_.AYY[KS] * *TRACK_mod_.Y + YYY) +
            *TRACK_mod_.W * (QSURF_.AXZ[KS] * *TRACK_mod_.X + QSURF_.AYZ[KS] * *TRACK_mod_.Y + QSURF_.AZZ[KS] * *TRACK_mod_.Z + ZZZ);

        C = *TRACK_mod_.X * XXX + *TRACK_mod_.Y * YYY + *TRACK_mod_.Z * ZZZ + QSURF_.A0[KS];
    }
    else{
        A = 0.0e0;
    	B = *TRACK_mod_.U * QSURF_.AX[KS] + *TRACK_mod_.V * QSURF_.AY[KS] + *TRACK_mod_.W * QSURF_.AZ[KS];
    	C = *TRACK_mod_.X * QSURF_.AX[KS] + *TRACK_mod_.Y * QSURF_.AY[KS] + *TRACK_mod_.Z * QSURF_.AZ[KS] + QSURF_.A0[KS];
	} 
     
     KS = KS + 1;
}

void locate2_(){
	
	
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
											
    for (int i = 1; i <= *QSURF_.NSURF; ++i) { 
	    QTREE_.KSP[i-1] = 0;
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
	
	
	
	
	
	*DSEF = 0.0e0;
	*PENGEOM_mod_.DSTOT = 0.0e0;
	*NCROSS = 0;
	*PENGEOM_mod_.KSLAST = 0;
	double DSRES;
	int KB1;
	double S[NS2M];
	int IS[NS2M];
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
	for (int i = 1; i <= *QSURF_.NSURF; i++){
	 	QTREE_.KSP[i - 1]  = 0; //Ponteiros laterais das superfícies avaliadas.
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
		stepsi2_(KB1, S, IS, NSC);
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
					for (int i = 1; i <= NSC; i++){
						S[i-1] = S[i-1] - DSP;
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
					stepsi2_(KB1, S, IS, NSC);
					goto L100;
				} else{
					//A particula entrou em um corpo material
					if (*TRACK_mod_.MAT != 0){
						*NCROSS = 1;
						return;
					} else{
						KB1 = *TRACK_mod_.IBODY;
						stepsi2_(KB1, S, IS, NSC);
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
	stepsi2_(KB1, S, IS, NSC);
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
	if (PENGEOM_mod_.KDET[*TRACK_mod_.IBODY - 1] != PENGEOM_mod_.KDET[IBODYL -1] || *TRACK_mod_.MAT != MAT0){
		*NCROSS = 1;
		*DSEF = 0.0e0;
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
			for (int i = 1; i <= NSC; i++){
				S[i-1] = S[i-1] - DSP;
			}
		}
		steplb2_(KB1, IERR);
		
		L201:;
		KB1 = *TRACK_mod_.IBODY;
		if (IERR == -1){
			// A particula entrou em um submodulo
			stepsi2_(KB1, S, IS, NSC);
			steplb2_(KB1, IERR);
			goto L201;
		} 
		else if (IERR == 1){
			//A partícula deixa o corpo ou módulo.
			if (*TRACK_mod_.IBODY <= *QTREE_.NBODYS){
				stepsi2_(KB1, S, IS, NSC);
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
				return;
			}
			L202:;
			stepsi2_(KB1, S, IS, NSC);
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
	return;
}


void stepsi2_(int &KB, double *S, int *IS, int &NSC){
	/*Calcula as interseções da trajetória com o limite
	Superfícies do corpo KB. Os cruzamentos são adicionados à lista e
	classificados em ordem decrescente. 
	Esta sub-rotina funciona apenas quando chamada de dentro da sub-rotina STEP.*/
	
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
		
		/*As interseções com uma determinada superfície são calculadas apenas uma vez.
		O ponteiro lateral de uma superfície deve ser alterado cada vez que o a superfície está cruzada.*/
		
		KFL = QTREE_.KFLAG[KSS - 1 ][KB -1];
		if (KFL > 4)
			goto L100;
		KS = QTREE_.KSURF[KSS - 1][KB - 1];
		if (QTREE_.KSP[KS - 1] != 0) 
			goto L100;
		fsurf2_(KS, A, B, C);
		ABSA = fabs(A);
		ABSB = fabs(B);
		
		// Plano, unica raiz
		if(ABSA < 1e-36){
			if (ABSB > 0.0e0){
				if (C < -FUZZL){
					QTREE_.KSP[KS -1] = 1;
				} 
				else if (C > FUZZL){
						QTREE_.KSP [KS -1] = 2;
				} else{
					if (B < 0.0e0){
						QTREE_.KSP [KS -1] = 1;
					} else{
						QTREE_.KSP [KS -1] = 2;
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
	
	if (NSC > 1){
		for (int KI = 1; KI <= NSC - 1; KI++){
			SMAX = S[KI -1];
			KMAX = KI;
			for ( int KJ = KI +1; KJ <= NSC; KJ++){
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
	
	return;
		
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
				if ((KF < 3) && (QTREE_.KSP[KS -1] != KF))
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
	return;
}



void geomin2_(double *PARINP, int *NPINP, int *NMATG, int *NBOD, int *IRD, int *IWR){
	
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
	
	FILE* IW = fopen("geometry-res.geo", "w");
	if (IW == NULL){
		printf("Não foi possível abrir o arquivo de geometria-res");
		exit(0);
	}
	
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
	C1[0] = CA[NINCL];
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
			VALUE = PARINP[ICHPAR];
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
	if (ICHPAR > 0)
		if (ICHPAR <= *NPINP){
			VALUE = PARINP[ICHPAR-1];
			ICHPAR = -ICHPAR; //Desativa a opção de alteração de parâmetro.
	} else{
		fprintf(IW, "%s(%f,%d,%s)\n", LKEYW, VALUE, ICHPAR, LANGLE);
		//fputs("linha 539", IW);
		fputs("*** NPINP is too small (check PARINP).", IW);
		exit(0);
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
		fprintf(IW, "%s(%$d)\n", LKEYW, KB);
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
			VALUE = PARINP[ICHPAR];
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
			VALUE = PARINP[ICHPAR];
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
					for ( int K = 1; K <=NSB; K++){
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
							QTREE_.KSURF[NXG-1][KB-1] = NSB -1;
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
			for (int J = 1; J <= QTREE_.KSURF[I-1][KB-1]; J++){
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
	for (int i = 1; i <= NXG; i++){ 
		if (QTREE_.KSURF[i-1][KB-1] != 0){
			fprintf(IW, "%5d", QTREE_.KSURF[i-1][KB-1]);
			if (i == 15)
			fprintf(IW, "\n       ");
		}
	}
	fprintf(IW, "\n");
	
}

void imprimirKFLAG(FILE* IW, int &KB){
	fprintf(IW, "KFLAG =");
	for (int i = 1; i <= NXG; i++){
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
	for (int i = 1; i <= NXG; i++){
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
	for (int i = 1; i <= NXG; i++){
		if (QTREE_.KDGHT[i-1][KB-1] != 0){
			fprintf(IW, "%5d", QTREE_.KDGHT[i-1][KB-1]);
			if (i == 15)
			   	fprintf(IW, "\n       ");
		}
	}
	fprintf(IW, "\n");
	
}

void peinit2_(double *EMAX, int *NMATER, int *IW, int *INFO, char (*PMFILE)[20]){
	
	
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
	char APOIO[80];
	
	
	FILE* IWR = fopen("material2.dat", "w");
	if (IWR == NULL){
		printf("N�o foi possivel abrir o arquivo material2.dat");
		exit(0);
	}
	fprintf(IWR, "\n **********************************\n");
	fprintf(IWR, " **   PENELOPE  (version 2014)   **\n");
	fprintf(IWR, " **********************************\n");
	
	*CSGAWR_.ISGAW = 0;
	*CERSEC_.IERSEC = 0;
	
	
	*PENELOPE_mod_.NMAT = *NMATER;
	for (int M = 1; M <= *PENELOPE_mod_.NMAT; M++){
		if (PENELOPE_mod_.EABS[M-1][1-1] < 49.999e0){
			fprintf(IWR, "EABS(1, %2d) = %.4E eV \n ERROR: electron absorption energy cannot be less than 50 eV\n", M, PENELOPE_mod_.EABS[M-1, 1-1]);
			printf("Electron absorption energy less than 50 eV.");
			exit(0);
		}
		EABS0[M-1][1-1] = PENELOPE_mod_.EABS[M-1][1-1];
		
		if (PENELOPE_mod_.EABS[M-1][2-1] < 49.999e0){
			fprintf(IWR, "EABS(2, %2d) = %.4E eV \n ERROR: photon absorption energy cannot be less than 50 eV\n", M, PENELOPE_mod_.EABS[M-1, 2-1]);
			printf("Photon absorption energy less than 50 eV.");
			exit(0);
		}
		EABS0[M-1][2-1] = PENELOPE_mod_.EABS[M-1][2-1];
		
		if (PENELOPE_mod_.EABS[M-1][3-1] < 49.999e0){
			fprintf(IWR, "EABS(3, %2d) = %.4E eV \n ERROR: positron absorption energy cannot be less than 50 eV.\n", M, PENELOPE_mod_.EABS[M-1, 3-1]);
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
	relax0_();
	relax02_(); //Inicializa rotinas de relaxamento at�mico.
	rndg30_();
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
		
		fprintf(IWR, "\nMaterial data file: %s\n\n", APOIO);
		
		FILE* IRD = fopen(APOIO, "r");
	
		if (IRD == NULL){
			fprintf(IWR, "N�o foi possivel abrir o arquivo %s\n", APOIO);
  	   		printf("N�o foi possivel abrir o arquivo %s\n", APOIO);
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
	
	printf("\n\nPEINIT\n\n");
	
	
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
	*CEGRID_.DLFC = log(*CEGRID_.EU / *CEGRID_.EL) / (*PENELOPE_mod_.NEGP - 1);
	*CEGRID_.DLEMP1 = log(*CEGRID_.EL);
	CEGRID_.DLEMP[1-1] = *CEGRID_.DLEMP1;
	CEGRID_.ET[1-1] = *CEGRID_.EL;
	
	for (int I = 2; I <= *PENELOPE_mod_.NEGP; I++){
		CEGRID_.DLEMP[I-1] = CEGRID_.DLEMP[I-1-1] + *CEGRID_.DLFC;
		CEGRID_.ET[I-1] = exp(CEGRID_.DLEMP[I-1]);
	} 
	
	*CEGRID_.DLFC = 1.0e0 / *CEGRID_.DLFC;
	
	printf("\n\nEGRID2\n\n");
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
	printf("\n\nESIA02\n\n");
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
	
	printf("\n\nPSIA02\n\n");
	
	
	
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
	printf("\nGPHA02\n");
		
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
	
	printf("\nRELAX02\n");
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
	
	printf("ERRM: %2.12f\n", ERRM );
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
	
	printf("\nRNDG302\n");
	
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
	
	printf("\nRITA02\n");
	
}


void pematr2_(int *M, FILE *IRD, FILE *IWR, int *INFO){	
	
	/*
	
	Esta sub-rotina l� o arquivo de defini��o do material M (unidade IRD) 
	e inicializa as rotinas de simula��o para este material. 
	A informa��o � escrita na unidade IWR, a quantidade de informa��o escrita 
	� determinada pelo valor de INFO.
	
	*/
	
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
	
    static const int NO = 512;
    
	
	
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
	
	char LNAME[] = " PENELOPE (v. 2014)  Material data file ...............";
	
	fgets(LINHA, sizeof(LINHA), IRD);
	
	extrairString(NAME, LINHA, 0, 55);
	if (strcmp(NAME, LNAME)){
		fprintf(IWR, "I/O error. Corrupt material data file.\n");
		fprintf(IWR, "     The first line is: %s\n", NAME);
		fprintf(IWR, "     ... and should be: %s\n", LNAME);
		printf("I/O error. Corrupt material data file.\n");
		printf("     The first line is: %s\n\n", NAME);
		printf("     ... and should be: %s\n", LNAME);
		printf("PEMATR. Corrupt material data file\n", LNAME);
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
            	for (int IEL = 1; IEL <= COMPOS_.NELEM[*M-1]; IEL++){
    		   		if (IZZ == COMPOS_.IZ[IEL-1][*M-1])
            		   	STFI[NS-1] = COMPOS_.STF[IEL-1][*M-1];
				}   		
			} else{
				NI=NI+1;
            	CPINAC_.IPIN[NI-1][*M-1] = KO;
            	CPSIN_.ISIP[KO-1] = -NI;
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
		
		//CEIMFP_.DEL[I-1][*M-1] = DELTA;
        
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
		
    
	
	

	printf("FIM PMATR\n");
//	 exit(0);
}

void irnd02_(double *W, double *F, int *K, int *N){ //OK
	
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
	FACT = (*N) * (*CRITAN_.CNORM);
//	printf("\n\nN: %d, CNORM: %f, valor FACT: %f", *N, *CRITAN_.CNORM, FACT);
	for (int I = 1; I <= *N; I++){
		K[I-1] = I;
		F[I-1] = W[I-1] * FACT;	
	}
	
	if (*N == 1){
		return;
	}
	
	for (int I = 1; I <= *N - 1; I++){
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
	
	
//	double XS[NIP];
//	double YS[NIP];
//	double SUMI[NIP];
//	double ERR[NM];
//	double C[NM];
	
	double *XS = (double *) malloc(NIP * sizeof(double));
	double *YS = (double *) malloc(NIP * sizeof(double));
	double *SUMI = (double *) malloc(NIP * sizeof(double));
	double *ERR = (double *) malloc(NM * sizeof(double));
	double *C = (double *) malloc(NM * sizeof(double));
	
	
	
	
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
		printf("Error in RITAI0: XLOW must be larger than XHIGH. XLOW=E%.6, XHIGH = E%.6\n", XLOW, XHIGH);
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
        	printf("Error in RITAI0: XLOW and NU are negative. XLOW= E%.7, NU=%11d", XLOW, NU);
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
    
    FILE* IW = fopen("param.dat", "w");
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
				PDFE = fmax(rndg3f2_(XS[I-1]), ZEROT)* *CRITAN_.CNORM;
				break;
			case 2: //dcsel2
				PDFE = fmax(dcsel2_(XS[I-1]), ZEROT)* *CRITAN_.CNORM;
				break;	
			case 3: //dcsel2
				PDFE = fmax(graaf22_(XS[I-1]), ZEROT)* *CRITAN_.CNORM;
				break;
		}
		fprintf(IW, "%f  %.6E %f  %f  %f  %d\n", CRITA_.XT[I-1], PDFE, CRITA_.A[I-1], CRITA_.B[I-1], C[I-1], ERR[I-1] );
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
			fprintf(IW1, "%f   %.8E   %f   %f   %f\n", XTAU, P1, P2, (P1-P2)/P1);	 
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
	 
	 free(XS);
	 free(YS);
	 free(SUMI);
	 free(ERR);
	 free(C);

	 
	 
	
//	 printf("\nRITAI02\n");
	
}

int teste = 1;

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
		printf("TST %.6E\n", TST);
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
	
//	printf("\n\nBRAAR2\n\n");	

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
	
	DELTA = 0.0;
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
		printf("X(1),PDF(1) = %d, %d    RMOMX. Error code 2.\n", X[1-1],PDF[1-1]);
		exit(0);
	}
	
	for (int I = 2; I<=NP; I++){
		if ((X[I-1] < 0.0e0) || (PDF[I-1] < 0.0e0)){
		   	printf("X(1),PDF(1) = %d, %d    RMOMX. Error code 2.\n", X[I-1],PDF[I-1]);
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
		printf("XLOW = %d, XUP = %d\n", XLOW, XUP);
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
		printf("*** The arguments in subroutine EELa0 are inconsistent. XS0 = %d, XS1 = %d\n", XS0,XS1);
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
	This subroutine reads the squared molecular form factor and the DCS
    for Rayleigh scattering of photons in material M. These two functions
    are tabulated using the same grids for all materials.
    The random sampling of the scattering angle is performed using the
    RITA algorithm.
	
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
		printf("MERGE2. Increase the value of the parameter NP = %d\n", fmax(N1, N2));
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

void pmrdr2(){

	//Lê o arquivo de entrada e inicializa PENELOPE e PENGEOM.


	char PMFILE[MAXMAT][21];
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







}















































