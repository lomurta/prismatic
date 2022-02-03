C
C  ************  Common blocks of the PENMAIN subroutines.
C  >>>>>>>>>>>>>>>>>>>>>>>>  DO NOT MODIFY!  <<<<<<<<<<<<<<<<<<<<<<<<<<<
C
      LOGICAL LDUMP,LSPEC,LGPOL,LPSF,LEXSRC,LEXBD,LSCONE,LFORCE,
     1   LXRSPL,LDOSEM
      CHARACTER PSFI*20,PFILED*20
      CHARACTER DATE23*23
      CHARACTER TITLE*65,TITLE2*65
C
      COMMON/CTITLE/TITLE,TITLE2
      COMMON/CDATE/DATE23
C
      COMMON/CSPGEO/DSMAX(NB),EABSB(3,NB)
C  ****  Interaction forcing, weight windows.
      COMMON/CFORCI/WLOW(NBV,3),WHIG(NBV,3),LFORCE(NBV,3)
C  ****  X-ray splitting.
      COMMON/CXRSPL/IXRSPL(NBV),ILBA(5),LXRSPL(NBV)
C
C  ****  Source definition.
C
C  ----  Primary particles.
      COMMON/CSOUR0/CTHL,DCTH,PHIL,DPHI,KPARP,JOBEND,LSCONE,LGPOL,LPSF
      COMMON/CSOUR1/E0,EPMAX,SP10,SP20,SP30
C  ----  Energy spectrum.
      PARAMETER (NSEM=1000)
      COMMON/CSOUR2/ESRC(NSEM),PSRC(NSEM),IASRC(NSEM),FSRC(NSEM),LSPEC
C  ----  Extended source.
      COMMON/CSOUR3/SX0,SY0,SZ0,SSX,SSY,SSZ,IXSBOD(NB),LEXSRC,LEXBD
C
C  ----  Input phase-space file.
      PARAMETER (NPSFM=100)
      COMMON/CSOUR4/WGMIN,RWGMIN,WGMAX,RLREAD,IPSFI,NPSF,NPSN,NSPLIT,
     1  KODEPS
      COMMON/CSOUR5/PSFI(NPSFM)
C
C  ****  Discrete counters.
C
      COMMON/CNT0/
     1  PRIM(3),PRIM2(3),DPRIM(3),    ! Numbers of IEXIT particles.
     1  SEC(3,3),SEC2(3,3),DSEC(3,3), ! Generated secondary particles.
     1  AVW(2),AVW2(2),DAVW(2),       ! Final polar director cosine.
     1  AVA(2),AVA2(2),DAVA(2),       ! Final polar angle.
     1  AVE(2),AVE2(2),DAVE(2)        ! Final energy.
C  ----  Deposited energies in various bodies.
      COMMON/CNT1/TDEBO(NB),TDEBO2(NB),DEBO(NB)
C
C  ****  Continuous distributions.
C
C  ----  Energy spectrum of the source.
      COMMON/CNT2/SHIST(NSEM),NSEB  ! Definition.
      COMMON/CNT3/SEDS(3,NSEM),SEDS2(3,NSEM),DSDE,RDSDE,NSDE  ! Check.
C
C  ----  Detectors (up to NIDM different detectors).
      PARAMETER (NIDM=25)
      COMMON/CNT4/RLAST,RWRITE,IDCUT(NIDM),KKDI(NIDM,3),IPSF(NIDM),
     1  NID,NPSFO,IPSFO
      COMMON/CNT5/DEDE(NIDM),KBDE(NB),NED
      COMMON/CNT6/LDOSEM
C
C  ****  Job details.
C
C  ----  Dump file.
      COMMON/CDUMP/LDUMP,PFILED
C  ----  Time control and shower counters.
      COMMON/CNTRL/TSIM,TSEC,TSECA,TSECAD,CPUT0,DUMPP,DSHN,SHN,N
