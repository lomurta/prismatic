C
C    This program adds the values of the counters in a number of dump
C  files generated from independent runs of PENMAIN. It allows the user
C  to run the same problem on several processors, using different seeds
C  of the random number generator, and then combine the results to
C  produce a single set of output files, with accumulated statistics.
C  It works as an effective "poor man's" parallelisation device.
C
C    The program reads the list of relative paths and filenames of the
C  dump files from a text file (standard input, unit 5, one filename in
C  each line). All dump files must correspond to the same problem, that
C  is, they must result from PENMAIN simulations with input files that
C  differ only in the seeds of the random number generator. The
C  responsibility of ensuring that the dump files do correspond to the
C  same problem rests with the user.
C
C                                        Francesc Salvat. August, 2014.
C


C  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      MODULE TRACK_mod   ! Particle TRACK variables.
C  ----  Energy, position, direction, and weight.
      DOUBLE PRECISION :: E,X,Y,Z,U,V,W,WGHT
C  ----  Particle type, current body (geometry parameter), and material.
      INTEGER*4 :: KPAR,IBODY,MAT
C  ----  Particle history flags.
      INTEGER*4, DIMENSION (5) :: ILB
C  ****  Photon polarisation
C  ----  Polarised photons if IPOL=1, otherwise unpolarised photons.
      INTEGER*4 :: IPOL=0
C  ----  Stokes parameters.
      DOUBLE PRECISION :: SP1,SP2,SP3
C  ****  The particle age (time elapsed since the start of the shower)
C        is recorded when LAGE=.TRUE.
      LOGICAL :: LAGE =.FALSE.
      DOUBLE PRECISION :: PAGE=0.0D0
      END MODULE TRACK_mod
C  <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


C  *********************************************************************
C                       MAIN PROGRAM
C  *********************************************************************
C
      USE TRACK_mod
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
      CHARACTER PFILER*40
      DIMENSION ISED1(100),ISED2(100)
C
C
      LOGICAL LDUMP,LSPEC,LGPOL,LPSF,LEXSRC,LEXBD,LSCONE,LFORCE,
     1   LXRSPL,LDOSEM
      CHARACTER PSFI*20,PFILED*20
      CHARACTER DATE23*23
      CHARACTER TITLE*65,TITLE2*65
C
      COMMON/CTITLE/TITLE,TITLE2
      COMMON/CDATE/DATE23
      COMMON/SUMAUX/NBODY
C
C  ****  Geometry.
      PARAMETER (NS=10000,NB=5000,NXG=250)
      COMMON/QBODY/KBODY(NB,NXG),KBOMO(NB)
      COMMON/QTREE/NBODYS,KMOTH(NB),KDGHT(NB,NXG),KSURF(NB,NXG),
     1  KFLAG(NB,NXG),KSP(NS),NWARN
      COMMON/QMAT/MATER(NB)
C
C  ****  Interaction forcing parameters.
      COMMON/CFORCE/FORCE(NB,3,8),IBRSPL(NB)
      COMMON/CFORCI/WLOW(NB,3),WHIG(NB,3),LFORCE(NB,3)
C  ****  X-ray splitting.
      COMMON/CXRSPL/IXRSPL(NB),ILBA(5),LXRSPL(NB)
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
     1  AVW(2),AVW2(2),DAVW(2),        ! Final polar director cosine.
     1  AVA(2),AVA2(2),DAVA(2),        ! Final polar angle.
     1  AVE(2),AVE2(2),DAVE(2)         ! Final energy.
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
      COMMON/CDUMP/LDUMP,PFILED,ISEED1,ISEED2
C  ----  Time control and shower counters.
      COMMON/CNTRL/TSIM,TSEC,TSECA,TSECAD,CPUT0,DUMPP,DSHN,SHN,N
C
C
      DIMENSION SEDSP(3,NSEM),AUX1(NSEM)
C
      DATE23='runs at different times'
C
C  ************  Read the first dump file.
C
      IDMP=1
      READ(5,'(A)',ERR=300,END=300) PFILER
      OPEN(9,FILE=PFILER)
      WRITE(6,'(A,A)') '  Reading the dump file ',PFILER
C
      READ (9,*,ERR=300,END=300) SHN,CPUT
      READ (9,'(A65)',ERR=300) TITLE
      READ (9,*,ERR=300) ISED1(IDMP),ISED2(IDMP)
      READ (9,*,ERR=300) NPSN,RLREAD
      READ (9,*,ERR=300) KPARP
      IF(KPARP.EQ.0) THEN
        READ (9,*,ERR=300) NSDE,DSDE,RDSDE
        READ (9,*,ERR=300) ((SEDS(K,I),I=1,NSDE),K=1,3)
        READ (9,*,ERR=300) ((SEDS2(K,I),I=1,NSDE),K=1,3)
      ELSE
        READ (9,*,ERR=300) NSEB
        READ (9,*,ERR=300) (ESRC(I),I=1,NSEB+1)
        READ (9,*,ERR=300) (PSRC(I),I=1,NSEB+1)
        READ (9,*,ERR=300) (SHIST(I),I=1,NSEB+1)
      ENDIF
      READ (9,*,ERR=300) (PRIM(I),I=1,3)
      READ (9,*,ERR=300) (PRIM2(I),I=1,3)
      READ (9,*,ERR=300) ((SEC(K,I),I=1,3),K=1,3)
      READ (9,*,ERR=300) ((SEC2(K,I),I=1,3),K=1,3)
      READ (9,*,ERR=300) (AVW(I),I=1,2)
      READ (9,*,ERR=300) (AVW2(I),I=1,2)
      READ (9,*,ERR=300) (AVA(I),I=1,2)
      READ (9,*,ERR=300) (AVA2(I),I=1,2)
      READ (9,*,ERR=300) (AVE(I),I=1,2)
      READ (9,*,ERR=300) (AVE2(I),I=1,2)
      READ (9,*,ERR=300) NBODY
      READ (9,*,ERR=300) (TDEBO(I),I=1,NBODY)
      READ (9,*,ERR=300) (TDEBO2(I),I=1,NBODY)
      CALL ENANGR(9)  ! Energy and angular distributions.
      READ (9,*,ERR=300) NID,NED,LDOSEM
      IF(NID.GT.0) THEN
        READ (9,*,ERR=300) RLAST,RWRITE
        CALL IMDETR(9)  ! Impact detectors.
      ENDIF
      IF(NED.GT.0) THEN
        CALL ENDETR(9)  ! Energy-deposition detectors.
      ENDIF
      IF(LDOSEM) THEN
        CALL DOSER(9)  ! Dose distribution.
      ENDIF
      CLOSE(9)
C
      DO KB=1,NBODY
        MATER(KB)=1
      ENDDO
C
C  ************  Read the second and following dump files.
C
  100 CONTINUE
      READ(5,'(A)',ERR=100,END=200) PFILER
      IDMP=IDMP+1
      OPEN(9,FILE=PFILER)
      WRITE(6,'(A,A)') '  Reading the dump file ',PFILER
C
      READ (9,*,ERR=100,END=100) SHNA,CPUTA
      SHN=SHN+SHNA
      CPUT=CPUT+CPUTA
      READ (9,'(A65)',ERR=100) TITLE2
      IF(TITLE2.NE.TITLE) THEN
        WRITE(6,*)
     1    'The dump file is corrupted (the TITLE does not match).'
        STOP 'The dump file is corrupted (the TITLE does not match).'
      ENDIF
      READ (9,*,ERR=100) ISED1(IDMP),ISED2(IDMP)
      READ (9,*,ERR=100) NPSN,RLREAD
      READ (9,*,ERR=100) KPARP
      IF(KPARP.EQ.0) THEN
        READ (9,*,ERR=100) NSDE,DSDE,RDSDE
        READ (9,*,ERR=100) ((SEDSP(K,I),I=1,NSDE),K=1,3)
        DO I=1,NSDE
          DO K=1,3
            SEDS(K,I)=SEDS(K,I)+SEDSP(K,I)
          ENDDO
        ENDDO
        READ (9,*,ERR=100) ((SEDSP(K,I),I=1,NSDE),K=1,3)
        DO I=1,NSDE
          DO K=1,3
            SEDS2(K,I)=SEDS2(K,I)+SEDSP(K,I)
          ENDDO
        ENDDO
      ELSE
        READ (9,*,ERR=100) NSEB
        READ (9,*,ERR=300) (ESRC(I),I=1,NSEB+1)
        READ (9,*,ERR=300) (PSRC(I),I=1,NSEB+1)
        READ (9,*,ERR=100) (AUX1(I),I=1,NSEB+1)
        DO I=1,NSEB
          SHIST(I)=SHIST(I)+AUX1(I)
        ENDDO
      ENDIF
C
      READ (9,*,ERR=100) (DPRIM(I),I=1,3)
      DO I=1,3
        PRIM(I)=PRIM(I)+DPRIM(I)
      ENDDO
C
      READ (9,*,ERR=100) (DPRIM(I),I=1,3)
      DO I=1,3
        PRIM2(I)=PRIM2(I)+DPRIM(I)
      ENDDO
C
      READ (9,*,ERR=100) ((DSEC(K,I),I=1,3),K=1,3)
      DO K=1,3
        DO I=1,3
          SEC(K,I)=SEC(K,I)+DSEC(K,I)
        ENDDO
      ENDDO
C
      READ (9,*,ERR=100) ((DSEC(K,I),I=1,3),K=1,3)
      DO K=1,3
        DO I=1,3
          SEC2(K,I)=SEC2(K,I)+DSEC(K,I)
        ENDDO
      ENDDO
C
      READ (9,*,ERR=100) (DAVW(I),I=1,2)
      DO I=1,2
        AVW(I)=AVW(I)+DAVW(I)
      ENDDO
C
      READ (9,*,ERR=100) (DAVW(I),I=1,2)
      DO I=1,2
        AVW2(I)=AVW2(I)+DAVW(I)
      ENDDO
C
      READ (9,*,ERR=100) (DAVA(I),I=1,2)
      DO I=1,2
        AVA(I)=AVA(I)+DAVA(I)
      ENDDO
C
      READ (9,*,ERR=100) (DAVA(I),I=1,2)
      DO I=1,2
        AVA2(I)=AVA2(I)+DAVA(I)
      ENDDO
C
      READ (9,*,ERR=100) (DAVE(I),I=1,2)
      DO I=1,2
        AVE(I)=AVE(I)+DAVE(I)
      ENDDO
C
      READ (9,*,ERR=100) (DAVE(I),I=1,2)
      DO I=1,2
        AVE2(I)=AVE2(I)+DAVE(I)
      ENDDO
C
      READ (9,*,ERR=100) NBODY
      READ (9,*,ERR=100) (DEBO(I),I=1,NBODY)
      DO I=1,NBODY
        TDEBO(I)=TDEBO(I)+DEBO(I)
      ENDDO
C
      READ (9,*,ERR=100) (DEBO(I),I=1,NBODY)
      DO I=1,NBODY
        TDEBO2(I)=TDEBO2(I)+DEBO(I)
      ENDDO
C
      CALL ENANGA(9)  ! Energy and angular distributions.
      READ (9,*,ERR=100) NID,NED,LDOSEM
      IF(NID.GT.0) THEN
        READ (9,*,ERR=100) RLAST,RWRITE
        CALL IMDETA(9)  ! Impact detectors.
      ENDIF
      IF(NED.GT.0) THEN
        CALL ENDETA(9)  ! Energy-deposition detectors.
      ENDIF
      IF(LDOSEM) THEN
        CALL DOSEA(9)  ! Dose distribution.
      ENDIF
      CLOSE(9)
C
      GO TO 100
  200 CONTINUE
C
      ISEED1=ISED1(IDMP)
      ISEED2=ISED2(IDMP)
      LDUMP=.TRUE.
      NPSFO=0
      TSIM=CPUT
      PFILED='all-dump.dat'
      CALL PMWRT
      WRITE(6,3040) SHN
 3040 FORMAT(2X,'Number of simulated showers =',1P,E14.7)
      WRITE(6,'(2X,''***  END  ***'')')
      STOP
C
  300 CONTINUE
      WRITE(6,'(2X,''***  Inconsistent dump file format.'')')
      STOP
      END


C  *********************************************************************
C                       SUBROUTINE PMWRT
C  *********************************************************************
      SUBROUTINE PMWRT
C
C  Computes averages and writes results in output files.
C
      USE TRACK_mod
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
C
      LOGICAL LDUMP,LSPEC,LGPOL,LPSF,LEXSRC,LEXBD,LSCONE,LFORCE,
     1   LXRSPL,LDOSEM
      CHARACTER PSFI*20,PFILED*20
      CHARACTER DATE23*23
      CHARACTER TITLE*65,TITLE2*65
C
      PARAMETER (PI=3.1415926535897932D0,RA2DE=180.0D0/PI)
C
      COMMON/CTITLE/TITLE,TITLE2
      COMMON/CDATE/DATE23
      COMMON/SUMAUX/NBODY
C
C  ****  Geometry.
      PARAMETER (NS=10000,NB=5000,NXG=250)
      COMMON/QBODY/KBODY(NB,NXG),KBOMO(NB)
      COMMON/QTREE/NBODYS,KMOTH(NB),KDGHT(NB,NXG),KSURF(NB,NXG),
     1  KFLAG(NB,NXG),KSP(NS),NWARN
      COMMON/QMAT/MATER(NB)
C
C  ****  Interaction forcing parameters.
      COMMON/CFORCE/FORCE(NB,3,8),IBRSPL(NB)
      COMMON/CFORCI/WLOW(NB,3),WHIG(NB,3),LFORCE(NB,3)
C  ****  X-ray splitting.
      COMMON/CXRSPL/IXRSPL(NB),ILBA(5),LXRSPL(NB)
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
      COMMON/CDUMP/LDUMP,PFILED,ISEED1,ISEED2
C  ----  Time control and shower counters.
      COMMON/CNTRL/TSIM,TSEC,TSECA,TSECAD,CPUT0,DUMPP,DSHN,SHN,N
C
      DIMENSION WSEC(3,3),WSEC2(3,3)
      DIMENSION WAVE2(2),WAVE(2),WAVW2(2),WAVW(2),WAVA2(2),WAVA(2)
C
C  ************  If 'DUMPTO' is active, write counters in a dump file.
C
      IF(LDUMP) THEN
        OPEN(9,FILE=PFILED)
        WRITE(9,*) SHN,TSIM
        WRITE(9,'(A65)') TITLE
        WRITE(9,*) ISEED1,ISEED2
C       WRITE(9,*) NPSN,RLREAD
        WRITE(9,*) 0,0.0D0
        WRITE(9,*) KPARP
        IF(KPARP.EQ.0) THEN
          WRITE(9,*) NSDE,DSDE,RDSDE
          WRITE(9,*) ((SEDS(K,I),I=1,NSDE),K=1,3)
          WRITE(9,*) ((SEDS2(K,I),I=1,NSDE),K=1,3)
        ELSE
          WRITE(9,*) NSEB
          WRITE(9,*) (ESRC(I),I=1,NSEB+1)
          WRITE(9,*) (PSRC(I),I=1,NSEB+1)
          WRITE(9,*) (SHIST(I),I=1,NSEB+1)
        ENDIF
        WRITE(9,*) (PRIM(I),I=1,3)
        WRITE(9,*) (PRIM2(I),I=1,3)
        WRITE(9,*) ((SEC(K,I),I=1,3),K=1,3)
        WRITE(9,*) ((SEC2(K,I),I=1,3),K=1,3)
        WRITE(9,*) (AVW(I),I=1,2)
        WRITE(9,*) (AVW2(I),I=1,2)
        WRITE(9,*) (AVA(I),I=1,2)
        WRITE(9,*) (AVA2(I),I=1,2)
        WRITE(9,*) (AVE(I),I=1,2)
        WRITE(9,*) (AVE2(I),I=1,2)
        WRITE(9,*) NBODY
        WRITE(9,*) (TDEBO(I),I=1,NBODY)
        WRITE(9,*) (TDEBO2(I),I=1,NBODY)
        CALL ENANGD(9)  ! Energy and angular distributions.
        WRITE(9,*) NID,NED,LDOSEM
        IF(NID.GT.0) THEN
          WRITE(9,*) RLAST,RWRITE
          CALL IMDETD(9)  ! Impact detectors.
        ENDIF
        IF(NED.GT.0) THEN
          CALL ENDETD(9)  ! Energy-deposition detectors.
        ENDIF
        IF(LDOSEM) THEN
          CALL DOSED(9)  ! Dose distribution.
        ENDIF
        WRITE(9,'(/3X,''*** END ***'')')
        CLOSE(9)
      ENDIF
C
C  ------------------------  Write simulation results.
C
C     IEXIT: 1=upbound, 2=downbound, 3=absorbed.
C
      OPEN(27,FILE='penmain-res.dat')
      WRITE(27,3000)
 3000 FORMAT(//3X,35('*')/3X,'**   Program PENMAIN. Results.   **',
     1  /3X,35('*'))
C
      WRITE(27,1001) DATE23
 1001 FORMAT(/3X,'Date and time: ',A23)
      WRITE(27,'(/3X,A65)') TITLE
C
      WRITE(27,3001) TSIM
 3001 FORMAT(/3X,'Simulation time ......................... ',
     1  1P,E13.6,' sec')
      TAVS=SHN/TSIM
      WRITE(27,3002) TAVS
 3002 FORMAT(3X,'Simulation speed ........................ ',
     1  1P,E13.6,' showers/sec')
      WRITE(27,3003) SHN
 3003 FORMAT(//3X,'Simulated primary particles ............. ',
     1  1P,E13.6)
C
      IF(KPARP.EQ.1) WRITE(27,1110)
 1110 FORMAT(/3X,'Primary particles: electrons')
      IF(KPARP.EQ.2) WRITE(27,1111)
 1111 FORMAT(/3X,'Primary particles: photons')
      IF(KPARP.EQ.3) WRITE(27,1112)
 1112 FORMAT(/3X,'Primary particles: positrons')
      IF(KPARP.EQ.0) WRITE(27,1113)
 1113 FORMAT(/3X,'Primary particles: set by the user subroutine SOURCE')
C
      WRITE(27,3004) PRIM(1)
 3004 FORMAT(/3X,'Upbound primary particles ............... ',
     1  1P,E13.6)
      WRITE(27,3005) PRIM(2)
 3005 FORMAT(3X,'Downbound primary particles ............. ',
     1  1P,E13.6)
      WRITE(27,3006) PRIM(3)
 3006 FORMAT(3X,'Absorbed primary particles .............. ',
     1  1P,E13.6)
C
      FNT=1.0D0/SHN
      IF(KPARP.NE.0) THEN
        FT=(PRIM(1)+SEC(KPARP,1))*FNT
        ERR1=3.0D0*FNT*SQRT(ABS(PRIM2(1)-PRIM(1)**2*FNT))
        ERR2=3.0D0*FNT*SQRT(ABS(SEC2(KPARP,1)-SEC(KPARP,1)**2*FNT))
        ERR=ERR1+ERR2
        WRITE(27,3007) FT,ERR
 3007   FORMAT(/3X,'Upbound fraction ................... ',
     1    1P,E13.6,' +-',E8.1)
        FB=(PRIM(2)+SEC(KPARP,2))*FNT
        ERR1=3.0D0*FNT*SQRT(ABS(PRIM2(2)-PRIM(2)**2*FNT))
        ERR2=3.0D0*FNT*SQRT(ABS(SEC2(KPARP,2)-SEC(KPARP,2)**2*FNT))
        ERR=ERR1+ERR2
        WRITE(27,3008) FB,ERR
 3008   FORMAT(3X,'Downbound fraction ................. ',
     1    1P,E13.6,' +-',E8.1)
        FA=PRIM(3)*FNT
        ERR=3.0D0*FNT*SQRT(ABS(PRIM2(3)-PRIM(3)**2*FNT))
        WRITE(27,3009) FA,ERR
 3009   FORMAT(3X,'Absorption fraction ................ ',
     1    1P,E13.6,' +-',E8.1)
      ENDIF
C
      DO K=1,3
        DO I=1,3
          WSEC2(K,I)=3.0D0*FNT*SQRT(ABS(SEC2(K,I)-SEC(K,I)**2*FNT))
          WSEC(K,I)=SEC(K,I)*FNT
        ENDDO
      ENDDO
      WRITE(27,3010)
     1  WSEC(1,1),WSEC(2,1),WSEC(3,1),WSEC2(1,1),WSEC2(2,1),WSEC2(3,1),
     1  WSEC(1,2),WSEC(2,2),WSEC(3,2),WSEC2(1,2),WSEC2(2,2),WSEC2(3,2),
     1  WSEC(1,3),WSEC(2,3),WSEC(3,3),WSEC2(1,3),WSEC2(2,3),WSEC2(3,3)
 3010 FORMAT(/3X,'Secondary-particle generation probabilities:',
     1  /19X,46('-'),
     1  /19X,'|  electrons   |   photons    |  positrons   |',1P,
     1  /3X,62('-')/3X,'|   upbound     |',3(E13.6,1X,'|'),
     1  /3X,'|               |',3('  +-',E8.1,2X,'|'),
     1  /3X,62('-')/3X,'|   downbound   |',3(E13.6,1X,'|'),
     1  /3X,'|               |',3('  +-',E8.1,2X,'|'),
     1  /3X,62('-')/3X,'|   absorbed    |',3(E13.6,1X,'|'),
     1  /3X,'|               |',3('  +-',E8.1,2X,'|'),
     1  /3X,62('-'))
C
      DO I=1,2
        DF=1.0D0/MAX(PRIM(I),1.0D0)
        WAVE2(I)=3.0D0*DF*SQRT(ABS(AVE2(I)-AVE(I)**2*DF))
        WAVE(I)=AVE(I)*DF
        WAVW2(I)=3.0D0*DF*SQRT(ABS(AVW2(I)-AVW(I)**2*DF))
        WAVW(I)=AVW(I)*DF
        WAVA2(I)=3.0D0*DF*RA2DE*SQRT(ABS(AVA2(I)-AVA(I)**2*DF))
        WAVA(I)=AVA(I)*RA2DE*DF
      ENDDO
C
      WRITE(27,3026)
 3026 FORMAT(/3X,'Average final energy:')
      WRITE(27,3027) WAVE(1),WAVE2(1)
 3027 FORMAT(6X,'Upbound primary particles ....... ',1P,E13.6,
     1  ' +-',E8.1,' eV')
      WRITE(27,3028) WAVE(2),WAVE2(2)
 3028 FORMAT(6X,'Downbound primary particles ..... ',1P,E13.6,
     1  ' +-',E8.1,' eV')
C
      WRITE(27,3033)
 3033 FORMAT(/3X,'Mean value of the polar cosine of the exit',
     1  ' direction:')
      WRITE(27,3034) WAVW(1),WAVW2(1)
 3034 FORMAT(6X,'Upbound primary particles ....... ',1P,E13.6,
     1  ' +-',E8.1)
      WRITE(27,3035) WAVW(2),WAVW2(2)
 3035 FORMAT(6X,'Downbound primary particles ..... ',1P,E13.6,
     1  ' +-',E8.1)
C
      WRITE(27,3036)
 3036 FORMAT(/3X,'Mean value of the polar angle of the exit dir',
     1  'ection:')
      WRITE(27,3037) WAVA(1),WAVA2(1)
 3037 FORMAT(6X,'Upbound primary particles ....... ',1P,E13.6,
     1  ' +-',E8.1,' deg')
      WRITE(27,3038) WAVA(2),WAVA2(2)
 3038 FORMAT(6X,'Downbound primary particles ..... ',1P,E13.6,
     1  ' +-',E8.1,' deg')
C
C  ****  Average energies deposited in bodies..
C
      DF=1.0D0/SHN
      WRITE(27,3011)
 3011 FORMAT(/3X,'Average deposited energies (bodies):')
      DO KB=1,NBODY
        IF(MATER(KB).NE.0) THEN
          QER=3.0D0*DF*SQRT(ABS(TDEBO2(KB)-TDEBO(KB)**2*DF))
          QAV=TDEBO(KB)*DF
          IF(QER.GT.1.0D-10*ABS(QAV)) THEN
            EFFIC=QAV**2/((QER/3.0D0)**2*TSIM)
          ELSE
            EFFIC=0.0D0
          ENDIF
          WRITE(27,3012) KB,QAV,QER,EFFIC
        ENDIF
      ENDDO
 3012 FORMAT(6X,'Body ',I4, ' ...... ',1P,E13.6,' +-',E8.1,' eV',4X,
     1  '(effic. =',E9.2,')')
C
C  ****  Output of impact detectors.
C
      IF(NID.GT.0) THEN
        CALL IMDETW(SHN,TSIM,27)
      ENDIF
C
C  ****  Output of energy-deposition detectors.
C
      IF(NED.GT.0) THEN
        CALL ENDETW(SHN,TSIM,27)
      ENDIF
C
C  ************  Energy spectrum of the source (as defined in PENMAIN).
C
      IF(NSEB.GT.1) THEN
        OPEN(9,FILE='psource.dat')
        WRITE(9,9000)
 9000 FORMAT(
     1  1X,'#  Results from PENMAIN. ',
     1 /1X,'#  Source energy spectrum.',
     1 /1X,'#  1st column: E (eV). 2nd column: spectrum (1/eV).',
     1 /1X,'#  3rd and 4th columns: simul. pdf limits (3SD, 1/eV).')
        PTOT=0.0D0
        WRITE(9,'(1P,4E14.6)') ESRC(1),PTOT,PTOT,PTOT
        DO KEn=1,NSEB
          PTOT=PTOT+PSRC(KEn)
        ENDDO
        DO KEn=1,NSEB
          YAV=SHIST(KEn)*DF
          YERR=3.0D0*SQRT(ABS(YAV*(1.0D0-YAV)*DF))
          EINTL=ESRC(KEn+1)-ESRC(KEn)
          IF(EINTL.GT.1.0D-15) THEN
            FACT=1.0D0/EINTL
          ELSE
            FACT=1.0D15
          ENDIF
          WRITE(9,'(1P,4E14.6)') ESRC(KEn),PSRC(KEn)*FACT/PTOT,
     1      (YAV-YERR)*FACT,(YAV+YERR)*FACT
          WRITE(9,'(1P,4E14.6)') ESRC(KEn+1),PSRC(KEn)*FACT/PTOT,
     1      (YAV-YERR)*FACT,(YAV+YERR)*FACT
        ENDDO
        PTOT=0.0D0
        WRITE(9,'(1P,4E14.6)') ESRC(NSEB+1),PTOT,PTOT,PTOT
        CLOSE(9)
      ENDIF
C
C  ************  Energy spectrum of the source (simulation result).
C
      IF(KPARP.EQ.0) THEN
        OPEN(9,FILE='usource.dat')
        WRITE(9,9010)
 9010 FORMAT(
     1  1X,'#  Results from PENMAIN. User-defined source.',
     1 /1X,'#  Energy spectra of particles from the source.',
     1 /1X,'#  1st column: particle energy (eV).',
     1 /1X,'#  2nd-3rd columns: electron spectrum and STU (3 sigma).',
     1 /1X,'#  4th-5th columns: photon spectrum and STU (3 sigma).',
     1 /1X,'#  6th-7th columns: positron spectrum and STU (3 sigma).')
        DO KEn=1,NSDE
          XX=(KEn-0.5D0)*DSDE
C
          YERR1=3.0D0*SQRT(ABS(SEDS2(1,KEn)-SEDS(1,KEn)**2*DF))
          YAV1=MAX(SEDS(1,KEn)*DF*RDSDE,1.0D-35)
          YERR1=MAX(YERR1*DF*RDSDE,1.0D-35)
C
          YERR2=3.0D0*SQRT(ABS(SEDS2(2,KEn)-SEDS(2,KEn)**2*DF))
          YAV2=MAX(SEDS(2,KEn)*DF*RDSDE,1.0D-35)
          YERR2=MAX(YERR2*DF*RDSDE,1.0D-35)
C
          YERR3=3.0D0*SQRT(ABS(SEDS2(3,KEn)-SEDS(3,KEn)**2*DF))
          YAV3=MAX(SEDS(3,KEn)*DF*RDSDE,1.0D-35)
          YERR3=MAX(YERR3*DF*RDSDE,1.0D-35)
C
          WRITE(9,'(1P,E14.6,3(E14.6,E10.2))')
     1      XX,YAV1,YERR1,YAV2,YERR2,YAV3,YERR3
        ENDDO
        CLOSE(9)
      ENDIF
C
C  ************  Energy and angular distributions of emerging particles.
C
      CALL ENANGW(SHN)  ! Energy and angular distributions.
C
C  ************  Dose distributions.
C
      IF(LDOSEM) THEN
        CALL DOSEW(SHN,TSIM)
      ENDIF
C
      WRITE(27,3030) ISEED1,ISEED2
 3030 FORMAT(/3X,'Last random seeds = ',I10,' , ',I10)
      WRITE(27,'(/3X,72(''-''))')
      CLOSE(27)
C
      RETURN
      END


C  *********************************************************************
C                       SUBROUTINE TENANG
C  *********************************************************************
      SUBROUTINE TENANG(IEXIT,N)
C
C  Tallies energy and angular distributions of emerging particles,
C  writes and loads dump files, accumulates dump files from different
C  runs, and writes results.
C
      USE TRACK_mod
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
      PARAMETER (PI=3.1415926535897932D0, TWOPI=2.0D0*PI,
     1  RA2DE=180.0D0/PI, DE2RA=PI/180.0D0)
      PARAMETER (FSAFE=1.000000001D0)  ! Safety factor.
      PARAMETER (NBEM=1500)  ! Max. no. of energy bins.
      PARAMETER (NBTHM=1800,NBPHM=180)  ! Max. nos of angular bins.
      COMMON/CENANG/EL,EU,THL,THU,BSE,RBSE,BSTH,RBSTH,BSPH,RBSPH,
     1  PDE(3,2,NBEM),PDE2(3,2,NBEM),PDEP(3,2,NBEM),
     2  PDA(3,NBTHM,NBPHM),PDA2(3,NBTHM,NBPHM),PDAP(3,NBTHM,NBPHM),
     4  LPDE(3,2,NBEM),LPDA(3,NBTHM,NBPHM),NE,NTH,NPH,LLE,LLTH
C
C  ************  Score contributions of a new particle.
C
C  ****  Energy distribution of emerging particles.
      IF(LLE.EQ.1) THEN
        KEn=1.0D0+(LOG(E)-EL)*RBSE
      ELSE
        KEn=1.0D0+(E-EL)*RBSE
      ENDIF
      IF(KEn.GT.0.AND.KEn.LE.NE) THEN
        IF(N.NE.LPDE(KPAR,IEXIT,KEn)) THEN
          PDE(KPAR,IEXIT,KEn)=PDE(KPAR,IEXIT,KEn)+PDEP(KPAR,IEXIT,KEn)
          PDE2(KPAR,IEXIT,KEn)=
     1      PDE2(KPAR,IEXIT,KEn)+PDEP(KPAR,IEXIT,KEn)**2
          PDEP(KPAR,IEXIT,KEn)=WGHT
          LPDE(KPAR,IEXIT,KEn)=N
        ELSE
          PDEP(KPAR,IEXIT,KEn)=PDEP(KPAR,IEXIT,KEn)+WGHT
        ENDIF
      ENDIF
C  ****  Angular distribution of emerging particles.
      THETA=ACOS(W)
      IF(LLTH.EQ.1) THEN
        KTH=1.0D0+(LOG(MAX(THETA,1.0D-12)*RA2DE)-THL)*RBSTH
        IF(KTH.LT.1) KTH=1
      ELSE
        KTH=1.0D0+THETA*RA2DE*RBSTH
      ENDIF
      IF(ABS(U).GT.1.0D-16) THEN  ! Azimuthal bin number corrected.
         PHI=ATAN2(V,U)
      ELSE IF(ABS(V).GT.1.0D-16) THEN
         PHI=ATAN2(V,U)
      ELSE
         PHI=0.0D0
      ENDIF
      IF(PHI.LT.0.0D0) PHI=TWOPI+PHI
      KPH=1.0D0+PHI*RA2DE*RBSPH
      IF(N.NE.LPDA(KPAR,KTH,KPH)) THEN
        PDA(KPAR,KTH,KPH)=PDA(KPAR,KTH,KPH)+PDAP(KPAR,KTH,KPH)
        PDA2(KPAR,KTH,KPH)=PDA2(KPAR,KTH,KPH)+PDAP(KPAR,KTH,KPH)**2
        PDAP(KPAR,KTH,KPH)=WGHT
        LPDA(KPAR,KTH,KPH)=N
      ELSE
        PDAP(KPAR,KTH,KPH)=PDAP(KPAR,KTH,KPH)+WGHT
      ENDIF
      RETURN
C
C  ************  Initialise the histograms.
C  Logarithmic (NB.LT.0) or uniform (NB.GT.0) scale.
C
      ENTRY ENANG0(EMIN,EMAX,NBE,NBTH,NBPH)  !<<<<<<<<<<<<<<<<<<<<<<<<<<
C
        IF (ABS(NBE).GT.NBEM) THEN
          WRITE(26,*) 'ENANG: NBE is too large.'
          WRITE(26,*) 'ENANG: Set the parameter NBEM equal to ',
     1      ABS(NBE)
          STOP 'ENANG: NBE is too large.'
        ENDIF
C
        IF (ABS(NBTH).GT.NBTHM) THEN
          WRITE(26,*) 'ENANG: NBTH is too large.'
          WRITE(26,*) 'ENANG: Set the parameter NBTHM equal to ',
     1      ABS(NBTH)
          STOP 'ENANG: NBTH is too large.'
        ENDIF
C
        IF (NBPH.GT.NBPHM) THEN
          WRITE(26,*) 'ENANG: NBPH is too large.'
          WRITE(26,*) 'ENANG: Set the parameter NBPHM equal to ',NBPH
          STOP 'ENANG: NBPH is too large.'
        ENDIF
C
        IF(NBE.LT.0) THEN
          LLE=1
          EL=LOG(EMIN)
          EU=LOG(EMAX)
          NE=-NBE
        ELSE
          LLE=0
          EL=EMIN
          EU=EMAX
          NE=NBE
        ENDIF
        BSE=FSAFE*(EU-EL)/DBLE(NE)
        RBSE=1.0D0/BSE
C
        IF(NBTH.LT.0) THEN
          LLTH=1
          THL=LOG(1.0D-2)
          THU=LOG(180.0D0)
          NTH=-NBTH
        ELSE
          LLTH=0
          THL=0.0D0
          THU=180.0D0
          NTH=NBTH
        ENDIF
        BSTH=FSAFE*(THU-THL)/DBLE(NTH)
        RBSTH=1.0D0/BSTH
C
        IF(NBPH.LT.0) THEN
          NPH=-NBPH
        ELSE
          NPH=NBPH
        ENDIF
        BSPH=FSAFE*360.0D0/DBLE(NPH)
        RBSPH=1.0D0/BSPH
C
        DO I=1,3
          DO J=1,2
            DO K=1,NBEM
              PDE(I,J,K)=0.0D0
              PDE2(I,J,K)=0.0D0
              PDEP(I,J,K)=0.0D0
              LPDE(I,J,K)=0
            ENDDO
          ENDDO
        ENDDO
C
        DO I=1,3
          DO J=1,NBTHM
            DO K=1,NBPHM
              PDA(I,J,K)=0.0D0
              PDA2(I,J,K)=0.0D0
              PDAP(I,J,K)=0.0D0
              LPDA(I,J,K)=0
            ENDDO
          ENDDO
        ENDDO
      RETURN
C
C  ************  DUMP. Write the histograms.
C
      ENTRY ENANGD(IWR)  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
C
C  ****  Transfer partial counters to global counters.
        DO KP=1,3
          DO IEX=1,2
            DO KEn=1,NE
              PDE(KP,IEX,KEn)=PDE(KP,IEX,KEn)+PDEP(KP,IEX,KEn)
              PDE2(KP,IEX,KEn)=PDE2(KP,IEX,KEn)+PDEP(KP,IEX,KEn)**2
              PDEP(KP,IEX,KEn)=0.0D0
              LPDE(KP,IEX,KEn)=0
            ENDDO
          ENDDO
        ENDDO
C
        DO KP=1,3
          DO KTH=1,NTH
            DO KPH=1,NPH
              PDA(KP,KTH,KPH)=PDA(KP,KTH,KPH)+PDAP(KP,KTH,KPH)
              PDA2(KP,KTH,KPH)=PDA2(KP,KTH,KPH)+PDAP(KP,KTH,KPH)**2
              PDAP(KP,KTH,KPH)=0.0D0
              LPDA(KP,KTH,KPH)=0
            ENDDO
          ENDDO
        ENDDO
C  ****  Write parameters and counters.
        WRITE(IWR,*) NE,EL,EU,BSE,RBSE,LLE
        WRITE(IWR,*) (((PDE(I,J,K),K=1,NE),J=1,2),I=1,3)
        WRITE(IWR,*) (((PDE2(I,J,K),K=1,NE),J=1,2),I=1,3)
        WRITE(IWR,*) NTH,THL,THU,BSTH,RBSTH,NPH,BSPH,RBSPH,LLTH
        WRITE(IWR,*) (((PDA(I,J,K),K=1,NPH),J=1,NTH),I=1,3)
        WRITE(IWR,*) (((PDA2(I,J,K),K=1,NPH),J=1,NTH),I=1,3)
      RETURN
C
C  ************  RESUME. Read histograms from a previous run.
C
      ENTRY ENANGR(IRD)  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
C
        READ(IRD,*) NE,EL,EU,BSE,RBSE,LLE
        READ(IRD,*) (((PDE(I,J,K),K=1,NE),J=1,2),I=1,3)
        READ(IRD,*) (((PDE2(I,J,K),K=1,NE),J=1,2),I=1,3)
        READ(IRD,*) NTH,THL,THU,BSTH,RBSTH,NPH,BSPH,RBSPH,LLTH
        READ(IRD,*) (((PDA(I,J,K),K=1,NPH),J=1,NTH),I=1,3)
        READ(IRD,*) (((PDA2(I,J,K),K=1,NPH),J=1,NTH),I=1,3)
C
        DO I=1,3
          DO J=1,2
            DO K=1,NE
              PDEP(I,J,K)=0.0D0
              LPDE(I,J,K)=0
            ENDDO
          ENDDO
        ENDDO
C
        DO I=1,3
          DO J=1,NTH
            DO K=1,NPH
              PDAP(I,J,K)=0.0D0
              LPDA(I,J,K)=0
            ENDDO
          ENDDO
        ENDDO
      RETURN
C
C  ************  ADD results from previous runs.
C  Use ENANGR to load histograms from one of the runs and then call
C  ENANGA to accumulate results from other runs.
C
      ENTRY ENANGA(IRD)  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
C
        READ(IRD,*) NE,EL,EU,BSE,RBSE,LLE
        READ(IRD,*) (((PDEP(I,J,K),K=1,NE),J=1,2),I=1,3)
        DO I=1,3
          DO J=1,2
            DO K=1,NE
              PDE(I,J,K)=PDE(I,J,K)+PDEP(I,J,K)
            ENDDO
          ENDDO
        ENDDO
        READ(IRD,*) (((PDEP(I,J,K),K=1,NE),J=1,2),I=1,3)
        DO I=1,3
          DO J=1,2
            DO K=1,NE
              PDE2(I,J,K)=PDE2(I,J,K)+PDEP(I,J,K)
              PDEP(I,J,K)=0.0D0
            ENDDO
          ENDDO
        ENDDO
        READ(IRD,*) NTH,THL,THU,BSTH,RBSTH,NPH,BSPH,RBSPH,LLTH
        READ(IRD,*) (((PDAP(I,J,K),K=1,NPH),J=1,NTH),I=1,3)
        DO I=1,3
          DO J=1,NTH
            DO K=1,NPH
              PDA(I,J,K)=PDA(I,J,K)+PDAP(I,J,K)
            ENDDO
          ENDDO
        ENDDO
        READ(IRD,*) (((PDAP(I,J,K),K=1,NPH),J=1,NTH),I=1,3)
        DO I=1,3
          DO J=1,NTH
            DO K=1,NPH
              PDA2(I,J,K)=PDA2(I,J,K)+PDAP(I,J,K)
              PDAP(I,J,K)=0.0D0
            ENDDO
          ENDDO
        ENDDO
      RETURN
C
C  ************  Write simulation results in separate output files.
C
      ENTRY ENANGW(SHN)  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
C
C  ****  Transfer partial counters to global counters.
        DO KP=1,3
          DO IEX=1,2
            DO KEn=1,NE
              PDE(KP,IEX,KEn)=PDE(KP,IEX,KEn)+PDEP(KP,IEX,KEn)
              PDE2(KP,IEX,KEn)=PDE2(KP,IEX,KEn)+PDEP(KP,IEX,KEn)**2
              PDEP(KP,IEX,KEn)=0.0D0
              LPDE(KP,IEX,KEn)=0
            ENDDO
          ENDDO
        ENDDO
C
        DO KP=1,3
          DO KTH=1,NTH
            DO KPH=1,NPH
              PDA(KP,KTH,KPH)=PDA(KP,KTH,KPH)+PDAP(KP,KTH,KPH)
              PDA2(KP,KTH,KPH)=PDA2(KP,KTH,KPH)+PDAP(KP,KTH,KPH)**2
              PDAP(KP,KTH,KPH)=0.0D0
              LPDA(KP,KTH,KPH)=0
            ENDDO
          ENDDO
        ENDDO
C
        DF=1.0D0/SHN
C
C  ****  Energy distributions of emerging particles.
C
C  ----  Upbound electrons.
        OPEN(9,FILE='energy-up.dat')
        WRITE(9,9110)
 9110   FORMAT(
     1    1X,'#  Results from PENMAIN.',
     1   /1X,'#  Energy distributions of upbound particles.',
     1   /1X,'#  1st column: E (eV).',
     1   /1X,'#  2nd and 3rd columns: PDF and STU for electrons.',
     1   /1X,'#  4th and 5th columns: PDF and STU for photons.',
     1   /1X,'#  6th and 7th columns: PDF and STU for positrons.',
     1   /1X,'#    PDFs and STUs in units of 1/(eV*primary_particle).')
        DO KEn=1,NE
          IF(LLE.EQ.1) THEN
            XLOW=EXP(EL+(KEn-1)*BSE)
            XUPP=EXP(EL+KEn*BSE)
            XX=0.5D0*(XUPP+XLOW)
            BINS=XUPP-XLOW
          ELSE
            XX=EL+(KEn-0.5D0)*BSE
            BINS=BSE
          ENDIF
          YERR1=3.0D0*SQRT(ABS(PDE2(1,1,KEn)-PDE(1,1,KEn)**2*DF))
          YAV1=PDE(1,1,KEn)*DF/BINS
          YERR1=YERR1*DF/BINS
          YERR2=3.0D0*SQRT(ABS(PDE2(2,1,KEn)-PDE(2,1,KEn)**2*DF))
          YAV2=PDE(2,1,KEn)*DF/BINS
          YERR2=YERR2*DF/BINS
          YERR3=3.0D0*SQRT(ABS(PDE2(3,1,KEn)-PDE(3,1,KEn)**2*DF))
          YAV3=PDE(3,1,KEn)*DF/BINS
          YERR3=YERR3*DF/BINS
          WRITE(9,'(1P,E14.6,3(E14.6,E10.2))')
     1      XX,MAX(YAV1,1.0D-35),YERR1,
     1      MAX(YAV2,1.0D-35),YERR2,MAX(YAV3,1.0D-35),YERR3
        ENDDO
        CLOSE(9)
C  ----  Downbound particles.
        OPEN(9,FILE='energy-down.dat')
        WRITE(9,9120)
 9120   FORMAT(
     1    1X,'#  Results from PENMAIN.',
     1   /1X,'#  Energy distributions of downbound particles.',
     1   /1X,'#  1st column: E (eV).',
     1   /1X,'#  2nd and 3rd columns: PDF and STU for electrons.',
     1   /1X,'#  4th and 5th columns: PDF and STU for photons.',
     1   /1X,'#  6th and 7th columns: PDF and STU for positrons.',
     1   /1X,'#    PDFs and STUs in units of 1/(eV*primary_particle).')
        DO KEn=1,NE
          IF(LLE.EQ.1) THEN
            XLOW=EXP(EL+(KEn-1)*BSE)
            XUPP=EXP(EL+KEn*BSE)
            XX=0.5D0*(XUPP+XLOW)
            BINS=XUPP-XLOW
          ELSE
            XX=EL+(KEn-0.5D0)*BSE
            BINS=BSE
          ENDIF
          YERR1=3.0D0*SQRT(ABS(PDE2(1,2,KEn)-PDE(1,2,KEn)**2*DF))
          YAV1=PDE(1,2,KEn)*DF/BINS
          YERR1=YERR1*DF/BINS
          YERR2=3.0D0*SQRT(ABS(PDE2(2,2,KEn)-PDE(2,2,KEn)**2*DF))
          YAV2=PDE(2,2,KEn)*DF/BINS
          YERR2=YERR2*DF/BINS
          YERR3=3.0D0*SQRT(ABS(PDE2(3,2,KEn)-PDE(3,2,KEn)**2*DF))
          YAV3=PDE(3,2,KEn)*DF/BINS
          YERR3=YERR3*DF/BINS
          WRITE(9,'(1P,E14.6,3(E14.6,E10.2))')
     1      XX,MAX(YAV1,1.0D-35),YERR1,
     1      MAX(YAV2,1.0D-35),YERR2,MAX(YAV3,1.0D-35),YERR3
        ENDDO
        CLOSE(9)
C
C  ************  Angular distributions of emerging particles.
C
        IF(NPH.GT.1) THEN
          OPEN(9,FILE='angle.dat')
          WRITE(9,9210)
 9210     FORMAT(
     1    1X,'#  Results from PENMAIN.',
     1     /1X,'#  Angular distributions of emerging particles.',
     1     /1X,'#  1st and 2nd columns: THETA and PHI (deg).',
     1     /1X,'#  3rd and 4th columns: PDF and STU for electrons.',
     1     /1X,'#  5th and 6th columns: PDF and STU for photons.',
     1     /1X,'#  7th and 8th columns: PDF and STU for positrons.',
     1     /1X,'#  PDFs and STUs in units of 1/(sr*primary_particle).')
          DO KTH=1,NTH
            IF(LLTH.EQ.1) THEN
              XLOW=EXP(THL+(KTH-1)*BSTH)
              XUPP=EXP(THL+KTH*BSTH)
            ELSE
              XLOW=THL+(KTH-1)*BSTH
              XUPP=THL+KTH*BSTH
            ENDIF
            XX=0.5D0*(XUPP+XLOW)
            BINS=XUPP-XLOW
            DSANG=(COS(XLOW*DE2RA)-COS(XUPP*DE2RA))*(BSPH*DE2RA)
            DO L=1,NPH
              YY=(L-0.5D0)*BSPH
              YERR1=3.0D0*SQRT(ABS(PDA2(1,KTH,L)-PDA(1,KTH,L)**2*DF))
              YAV1=PDA(1,KTH,L)*DF/DSANG
              YERR1=YERR1*DF/DSANG
              YERR2=3.0D0*SQRT(ABS(PDA2(2,KTH,L)-PDA(2,KTH,L)**2*DF))
              YAV2=PDA(2,KTH,L)*DF/DSANG
              YERR2=YERR2*DF/DSANG
              YERR3=3.0D0*SQRT(ABS(PDA2(3,KTH,L)-PDA(3,KTH,L)**2*DF))
              YAV3=PDA(3,KTH,L)*DF/DSANG
              YERR3=YERR3*DF/DSANG
              WRITE(9,'(1P,2E14.6,3(E14.6,E10.2))')
     1          XX,YY,MAX(YAV1,1.0D-35),YERR1,
     1          MAX(YAV2,1.0D-35),YERR2,MAX(YAV3,1.0D-35),YERR3
            ENDDO
            WRITE(9,*) '   '
          ENDDO
          CLOSE(9)
        ENDIF
C
        OPEN(9,FILE='polar-angle.dat')
        WRITE(9,9310)
 9310   FORMAT(
     1    1X,'#  Results from PENMAIN.',
     1   /1X,'#  Angular distributions of emerging particles.',
     1   /1X,'#  1st column: THETA (deg).',
     1   /1X,'#  2nd and 3rd columns: PDF and STU for electrons.',
     1   /1X,'#  4th and 5th columns: PDF and STU for photons.',
     1   /1X,'#  6th and 7th columns: PDF and STU for positrons',
     1   /1X,'#  PDFs and STUs in units of 1/(sr*primary_particle).')
        DO KTH=1,NTH
          IF(LLTH.EQ.1) THEN
            XLOW=EXP(THL+(KTH-1)*BSTH)
            XUPP=EXP(THL+KTH*BSTH)
          ELSE
            XLOW=THL+(KTH-1)*BSTH
            XUPP=THL+KTH*BSTH
          ENDIF
          XX=0.5D0*(XUPP+XLOW)
          BINS=XUPP-XLOW
          DSANG=(COS(XLOW*DE2RA)-COS(XUPP*DE2RA))*TWOPI
C
          S1=0.0D0
          S2=0.0D0
          DO L=1,NPH
            S1=S1+PDA(1,KTH,L)
            S2=S2+PDA2(1,KTH,L)
          ENDDO
          YERR1=3.0D0*SQRT(ABS(S2-S1**2*DF))
          YAV1=S1*DF/DSANG
          YERR1=YERR1*DF/DSANG
C
          S1=0.0D0
          S2=0.0D0
          DO L=1,NPH
            S1=S1+PDA(2,KTH,L)
            S2=S2+PDA2(2,KTH,L)
          ENDDO
          YERR2=3.0D0*SQRT(ABS(S2-S1**2*DF))
          YAV2=S1*DF/DSANG
          YERR2=YERR2*DF/DSANG
C
          S1=0.0D0
          S2=0.0D0
          DO L=1,NPH
            S1=S1+PDA(3,KTH,L)
            S2=S2+PDA2(3,KTH,L)
          ENDDO
          YERR3=3.0D0*SQRT(ABS(S2-S1**2*DF))
          YAV3=S1*DF/DSANG
          YERR3=YERR3*DF/DSANG
C
          WRITE(9,'(1P,E14.6,3(E14.6,E10.2))')
     1      XX,MAX(YAV1,1.0D-35),YERR1,
     1      MAX(YAV2,1.0D-35),YERR2,MAX(YAV3,1.0D-35),YERR3
        ENDDO
        CLOSE(9)
C
      RETURN
      END


C  *********************************************************************
C                       SUBROUTINE SIMDET
C  *********************************************************************
      SUBROUTINE SIMDET(N,ID)
C
C  Tallies spectra from impact detectors, writes and loads dump files,
C  accumulates dump files from different runs, and writes results.
C
      USE TRACK_mod
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
      PARAMETER (FSAFE=1.000000001D0)  ! Safety factor.
      CHARACTER FNSPC*20,FNFLU*20,FNAGE*20
      PARAMETER (NIDM=25,NBEM=1000)
      CHARACTER SPCDIO(NIDM)*20,SPCFLO(NIDM)*20,SPCAGE(NIDM)*20
      SAVE SPCDIO,SPCFLO,SPCAGE
C
      COMMON/CIMDET/EL(NIDM),EU(NIDM),BSE(NIDM),RBSE(NIDM),
     1  ET(NIDM,NBEM+1),EDEP(NIDM),EDEP2(NIDM),EDEPP(NIDM),
     1  DIT(NIDM,NBEM),DIT2(NIDM,NBEM),DITP(NIDM,NBEM),
     1  DIP(NIDM,NBEM,3),DIP2(NIDM,NBEM,3),DIPP(NIDM,NBEM,3),
     1  FLT(NIDM,NBEM),FLT2(NIDM,NBEM),FLTP(NIDM,NBEM),
     1  FLP(NIDM,NBEM,3),FLP2(NIDM,NBEM,3),FLPP(NIDM,NBEM,3),
     1  AGEL(NIDM),AGEU(NIDM),BAGE(NIDM),RBAGE(NIDM),AGE(NIDM,NBEM),
     1  AGE2(NIDM,NBEM),AGEP(NIDM,NBEM),LEDEP(NIDM),LDIT(NIDM,NBEM),
     1  LDIP(NIDM,NBEM,3),LFLT(NIDM,NBEM),LFLP(NIDM,NBEM,3),
     1  LAGEA(NIDM,NBEM),IDCUT(NIDM),NE(NIDM),LLE(NIDM),LLAGE(NIDM),
     1  NAGE(NIDM),NID
C
      DATA NIDS/0/
      SAVE NIDS
C
C  ************  Energy spectrum of entering particles.
C
        IF(LLE(ID).EQ.1) THEN
          IE=1.0D0+(LOG(E)-EL(ID))*RBSE(ID)
        ELSE
          IE=1.0D0+(E-EL(ID))*RBSE(ID)
        ENDIF
C
        IF(N.NE.LEDEP(ID)) THEN
          EDEP(ID)=EDEP(ID)+EDEPP(ID)
          EDEP2(ID)=EDEP2(ID)+EDEPP(ID)**2
          EDEPP(ID)=E*WGHT
          LEDEP(ID)=N
        ELSE
          EDEPP(ID)=EDEPP(ID)+E*WGHT
        ENDIF
C
        IF(IE.GT.0.AND.IE.LE.NE(ID)) THEN
          IF(N.NE.LDIT(ID,IE)) THEN
            DIT(ID,IE)=DIT(ID,IE)+DITP(ID,IE)
            DIT2(ID,IE)=DIT2(ID,IE)+DITP(ID,IE)**2
            DITP(ID,IE)=WGHT
            LDIT(ID,IE)=N
          ELSE
            DITP(ID,IE)=DITP(ID,IE)+WGHT
          ENDIF
C
          IF(N.NE.LDIP(ID,IE,KPAR)) THEN
            DIP(ID,IE,KPAR)=DIP(ID,IE,KPAR)+DIPP(ID,IE,KPAR)
            DIP2(ID,IE,KPAR)=
     1        DIP2(ID,IE,KPAR)+DIPP(ID,IE,KPAR)**2
            DIPP(ID,IE,KPAR)=WGHT
            LDIP(ID,IE,KPAR)=N
          ELSE
            DIPP(ID,IE,KPAR)=DIPP(ID,IE,KPAR)+WGHT
          ENDIF
C  ****  Age distribution.
          IF(NAGE(ID).GT.0) THEN
            IF(LLAGE(ID).EQ.1) THEN
              IT=1.0D0+(LOG(PAGE)-AGEL(ID))*RBAGE(ID)
            ELSE
              IT=1.0D0+(PAGE-AGEL(ID))*RBAGE(ID)
            ENDIF
            IF(IT.GT.0.AND.IT.LE.NAGE(ID)) THEN
              IF(N.NE.LAGEA(ID,IT)) THEN
                AGE(ID,IT)=AGE(ID,IT)+AGEP(ID,IT)
                AGE2(ID,IT)=AGE2(ID,IT)+AGEP(ID,IT)**2
                AGEP(ID,IT)=WGHT
                LAGEA(ID,IT)=N
              ELSE
                AGEP(ID,IT)=AGEP(ID,IT)+WGHT
              ENDIF
            ENDIF
          ENDIF
        ENDIF
      RETURN
C
C  ************  Fluence distribution of particles within the
C                detector. Discrete collisions only.
C
      ENTRY FIMDET(N,ID,DSEF)  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
C
        IF(LLE(ID).EQ.1) THEN
          IE=1.0D0+(LOG(E)-EL(ID))*RBSE(ID)
        ELSE
          IE=1.0D0+(E-EL(ID))*RBSE(ID)
        ENDIF
C
        IF(IE.GT.0.AND.IE.LE.NE(ID)) THEN
          IF(N.NE.LFLT(ID,IE)) THEN
            FLT(ID,IE)=FLT(ID,IE)+FLTP(ID,IE)
            FLT2(ID,IE)=FLT2(ID,IE)+FLTP(ID,IE)**2
            FLTP(ID,IE)=WGHT*DSEF
            LFLT(ID,IE)=N
          ELSE
            FLTP(ID,IE)=FLTP(ID,IE)+WGHT*DSEF
          ENDIF
          IF(N.NE.LFLP(ID,IE,KPAR)) THEN
            FLP(ID,IE,KPAR)=FLP(ID,IE,KPAR)+FLPP(ID,IE,KPAR)
            FLP2(ID,IE,KPAR)=
     1          FLP2(ID,IE,KPAR)+FLPP(ID,IE,KPAR)**2
            FLPP(ID,IE,KPAR)=WGHT*DSEF
            LFLP(ID,IE,KPAR)=N
          ELSE
            FLPP(ID,IE,KPAR)=FLPP(ID,IE,KPAR)+WGHT*DSEF
          ENDIF
        ENDIF
      RETURN
C
C  ************  Fluence distribution of particles within the
C                detector. Continuous slowing down
C
      ENTRY FIMDES(N,ID,EI,DECSD,DSEF)  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
C
        IF(EI.LT.EL(ID)) RETURN
        EIC=MIN(EI,EU(ID))
        EF=MAX(EI-DECSD,EL(ID))
        IF(EF.GT.EIC) RETURN
        IF(LLE(ID).EQ.1) THEN
          IEI=1.0D0+(LOG(EIC)-EL(ID))*RBSE(ID)
          IEF=1.0D0+(LOG(EF)-EL(ID))*RBSE(ID)
        ELSE
          IEI=1.0D0+(EIC-EL(ID))*RBSE(ID)
          IEF=1.0D0+(EF-EL(ID))*RBSE(ID)
        ENDIF
C
        FACT=DSEF/DECSD
        DO IE=IEF,IEI
          EA=MAX(EF,ET(ID,IE))
          EB=MIN(EIC,ET(ID,IE+1))
          TLBIN=(EB-EA)*FACT
          IF(N.NE.LFLT(ID,IE)) THEN
            FLT(ID,IE)=FLT(ID,IE)+FLTP(ID,IE)
            FLT2(ID,IE)=FLT2(ID,IE)+FLTP(ID,IE)**2
            FLTP(ID,IE)=WGHT*TLBIN
            LFLT(ID,IE)=N
          ELSE
            FLTP(ID,IE)=FLTP(ID,IE)+WGHT*TLBIN
          ENDIF
          IF(N.NE.LFLP(ID,IE,KPAR)) THEN
            FLP(ID,IE,KPAR)=FLP(ID,IE,KPAR)+FLPP(ID,IE,KPAR)
            FLP2(ID,IE,KPAR)=
     1          FLP2(ID,IE,KPAR)+FLPP(ID,IE,KPAR)**2
            FLPP(ID,IE,KPAR)=WGHT*TLBIN
            LFLP(ID,IE,KPAR)=N
          ELSE
            FLPP(ID,IE,KPAR)=FLPP(ID,IE,KPAR)+WGHT*TLBIN
          ENDIF
        ENDDO
      RETURN
C
C  ************  Initialise a new detector.
C  Logarithmic (NB.LT.0) or uniform (NB.GT.0) scale.
C
      ENTRY IMDET0(EMIN,EMAX,NBE,AGEMIN,AGEMAX,NBAGE,  !<<<<<<<<<<<<<<<<
     1        ICUT,FNSPC,FNFLU,FNAGE,ID)
C
        IF(ID.LE.NIDS) THEN
          WRITE(26,*) 'SIMDET: Detector already defined.',ID
          STOP 'SIMDET: Detector cannot be defined.'
        ENDIF
        NIDS=ID
        NID=ID
C
        IF(NID.GT.NIDM) THEN
          WRITE(26,'(3X,''NID = '',I4)') NID
          WRITE(26,*) 'SIMDET: Too many detectors.'
          STOP 'SIMDET: Too many detectors.'
        ENDIF
C
        SPCDIO(ID)=FNSPC
        SPCFLO(ID)=FNFLU
        SPCAGE(ID)=FNAGE
        IDCUT(ID)=ICUT
C
        IF (ABS(NBE).GT.NBEM) THEN
          WRITE(26,*) 'SIMDET: NB is too large.'
          WRITE(26,*) 'SIMDET: Set the parameter NBEM equal to ',
     1      ABS(NBE)
          STOP 'SIMDET: NB is too large.'
        ENDIF
C
        IF(NBE.LT.0) THEN
          LLE(ID)=1
          EL(ID)=LOG(EMIN)
          EU(ID)=LOG(EMAX)
          NE(ID)=-NBE
        ELSE
          LLE(ID)=0
          EL(ID)=EMIN
          EU(ID)=EMAX
          NE(ID)=NBE
        ENDIF
        BSE(ID)=FSAFE*(EU(ID)-EL(ID))/DBLE(NE(ID))
        RBSE(ID)=1.0D0/BSE(ID)
C
        DO J=1,NE(ID)+1
          IF(LLE(ID).EQ.1) THEN
            ET(ID,J)=EXP(EL(ID)+(J-1)*BSE(ID))
          ELSE
            ET(ID,J)=EL(ID)+(J-1)*BSE(ID)
          ENDIF
        ENDDO
C
        EDEP(ID)=0.0D0
        EDEP2(ID)=0.0D0
        EDEPP(ID)=0.0D0
        LEDEP(ID)=0
C
        DO J=1,NBEM
          DIT(ID,J)=0.0D0
          DIT2(ID,J)=0.0D0
          DITP(ID,J)=0.0D0
          LDIT(ID,J)=0
          DO K=1,3
            DIP(ID,J,K)=0.0D0
            DIP2(ID,J,K)=0.0D0
            DIPP(ID,J,K)=0.0D0
            LDIP(ID,J,K)=0
          ENDDO
        ENDDO
C
        DO J=1,NBEM
          FLT(ID,J)=0.0D0
          FLT2(ID,J)=0.0D0
          FLTP(ID,J)=0.0D0
          LFLT(ID,J)=0
          DO K=1,3
            FLP(ID,J,K)=0.0D0
            FLP2(ID,J,K)=0.0D0
            FLPP(ID,J,K)=0.0D0
            LFLP(ID,J,K)=0
          ENDDO
        ENDDO
C
        IF(NBAGE.LT.0) THEN
          LLAGE(ID)=1
          AGEL(ID)=LOG(AGEMIN)
          AGEU(ID)=LOG(AGEMAX)
          NAGE(ID)=-NBAGE
        ELSE IF(NBAGE.GT.0) THEN
          LLAGE(ID)=0
          AGEL(ID)=AGEMIN
          AGEU(ID)=AGEMAX
          NAGE(ID)=NBAGE
        ELSE
          LLAGE(ID)=0
          NAGE(ID)=0
          AGEL(ID)=0.0D0
          AGEU(ID)=0.0D0
          BAGE(ID)=0.0D0
          RBAGE(ID)=0.0D0
        ENDIF
C
        IF(NAGE(ID).GT.0) THEN
          LAGE=.TRUE.
          BAGE(ID)=FSAFE*(AGEU(ID)-AGEL(ID))/DBLE(NAGE(ID))
          RBAGE(ID)=1.0D0/BAGE(ID)
          DO J=1,NBEM
            AGE(ID,J)=0.0D0
            AGE2(ID,J)=0.0D0
            AGEP(ID,J)=0.0D0
            LAGEA(ID,J)=0
          ENDDO
        ENDIF
      RETURN
C
C  ************  DUMP. Write the histograms.
C
      ENTRY IMDETD(IWR)  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
C
C  ****  Transfer partial counters to global counters.
        DO I=1,NID
          EDEP(I)=EDEP(I)+EDEPP(I)
          EDEP2(I)=EDEP2(I)+EDEPP(I)**2
          EDEPP(I)=0.0D0
          LEDEP(I)=0
          DO J=1,NE(I)
            DIT(I,J)=DIT(I,J)+DITP(I,J)
            DIT2(I,J)=DIT2(I,J)+DITP(I,J)**2
            DITP(I,J)=0.0D0
            LDIT(I,J)=0
            DO K=1,3
              DIP(I,J,K)=DIP(I,J,K)+DIPP(I,J,K)
              DIP2(I,J,K)=DIP2(I,J,K)+DIPP(I,J,K)**2
              DIPP(I,J,K)=0.0D0
              LDIP(I,J,K)=0
            ENDDO
          ENDDO
        ENDDO
C
        DO I=1,NID
          DO J=1,NE(I)
            FLT(I,J)=FLT(I,J)+FLTP(I,J)
            FLT2(I,J)=FLT2(I,J)+FLTP(I,J)**2
            FLTP(I,J)=0.0D0
            LFLT(I,J)=0
            DO K=1,3
              FLP(I,J,K)=FLP(I,J,K)+FLPP(I,J,K)
              FLP2(I,J,K)=FLP2(I,J,K)+FLPP(I,J,K)**2
              FLPP(I,J,K)=0.0D0
              LFLP(I,J,K)=0
            ENDDO
          ENDDO
        ENDDO
C
        DO I=1,NID
          DO J=1,NAGE(I)
            AGE(I,J)=AGE(I,J)+AGEP(I,J)
            AGE2(I,J)=AGE2(I,J)+AGEP(I,J)**2
            AGEP(I,J)=0.0D0
            LAGEA(I,J)=0
          ENDDO
        ENDDO
C  ****  Write parameters and counters.
        WRITE(IWR,*) NID
        WRITE(IWR,*) (EDEP(I),I=1,NID)
        WRITE(IWR,*) (EDEP2(I),I=1,NID)
        WRITE(IWR,*) (IDCUT(I),I=1,NID)
        DO I=1,NID
          WRITE(IWR,'(A)') SPCDIO(I)
          WRITE(IWR,*) NE(I),EL(I),EU(I),BSE(I),RBSE(I),LLE(I)
          WRITE(IWR,*) (DIT(I,J),J=1,NE(I))
          WRITE(IWR,*) (DIT2(I,J),J=1,NE(I))
          WRITE(IWR,*) ((DIP(I,J,K),J=1,NE(I)),K=1,3)
          WRITE(IWR,*) ((DIP2(I,J,K),J=1,NE(I)),K=1,3)
          IF(IDCUT(I).EQ.2) THEN
            WRITE(IWR,'(A)') SPCFLO(I)
            WRITE(IWR,*) (FLT(I,J),J=1,NE(I))
            WRITE(IWR,*) (FLT2(I,J),J=1,NE(I))
            WRITE(IWR,*) ((FLP(I,J,K),J=1,NE(I)),K=1,3)
            WRITE(IWR,*) ((FLP2(I,J,K),J=1,NE(I)),K=1,3)
          ENDIF
          WRITE(IWR,*)
     1      NAGE(I),AGEL(I),AGEU(I),BAGE(I),RBAGE(I),LLAGE(I)
          IF(NAGE(I).GT.0) THEN
            WRITE(IWR,'(A)') SPCAGE(I)
            WRITE(IWR,*) (AGE(I,J),J=1,NAGE(I))
            WRITE(IWR,*) (AGE2(I,J),J=1,NAGE(I))
          ENDIF
        ENDDO
      RETURN
C
C  ************  RESUME. Read histograms from a previous run.
C
      ENTRY IMDETR(IRD)  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
C
        READ(IRD,*) NID
        READ(IRD,*) (EDEP(I),I=1,NID)
        READ(IRD,*) (EDEP2(I),I=1,NID)
        READ(IRD,*) (IDCUT(I),I=1,NID)
        DO I=1,NID
          READ(IRD,'(A)') SPCDIO(I)
          READ(IRD,*) NE(I),EL(I),EU(I),BSE(I),RBSE(I),LLE(I)
          READ(IRD,*) (DIT(I,J),J=1,NE(I))
          READ(IRD,*) (DIT2(I,J),J=1,NE(I))
          READ(IRD,*) ((DIP(I,J,K),J=1,NE(I)),K=1,3)
          READ(IRD,*) ((DIP2(I,J,K),J=1,NE(I)),K=1,3)
          IF(IDCUT(I).EQ.2) THEN
            READ(IRD,'(A)') SPCFLO(I)
            READ(IRD,*) (FLT(I,J),J=1,NE(I))
            READ(IRD,*) (FLT2(I,J),J=1,NE(I))
            READ(IRD,*) ((FLP(I,J,K),J=1,NE(I)),K=1,3)
            READ(IRD,*) ((FLP2(I,J,K),J=1,NE(I)),K=1,3)
          ENDIF
          READ(IRD,*)
     1      NAGE(I),AGEL(I),AGEU(I),BAGE(I),RBAGE(I),LLAGE(I)
          IF(NAGE(I).GT.0) THEN
            READ(IRD,'(A)') SPCAGE(I)
            READ(IRD,*) (AGE(I,J),J=1,NAGE(I))
            READ(IRD,*) (AGE2(I,J),J=1,NAGE(I))
          ENDIF
        ENDDO
C
        DO I=1,NID
          EDEPP(I)=0.0D0
          LEDEP(I)=0
          DO J=1,NE(I)
            DITP(I,J)=0.0D0
            LDIT(I,J)=0
            DO K=1,3
              DIPP(I,J,K)=0.0D0
              LDIP(I,J,K)=0
            ENDDO
          ENDDO
        ENDDO
C
        DO I=1,NID
          DO J=1,NE(I)
            FLTP(I,J)=0.0D0
            LFLT(I,J)=0
            DO K=1,3
              FLPP(I,J,K)=0.0D0
              LFLP(I,J,K)=0
            ENDDO
          ENDDO
        ENDDO
C
        DO I=1,NID
          DO J=1,NAGE(I)
            AGEP(I,J)=0.0D0
            LAGEA(I,J)=0
          ENDDO
        ENDDO
      RETURN
C
C  ************  ADD results from previous runs.
C  Use IMDETR to load histograms from one of the runs and then call
C  IMDETA to accumulate results from other runs.
C
      ENTRY IMDETA(IRD)  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
C
        READ(IRD,*) NID
        READ(IRD,*) (EDEPP(I),I=1,NID)
        DO I=1,NID
          EDEP(I)=EDEP(I)+EDEPP(I)
        ENDDO
        READ(IRD,*) (EDEPP(I),I=1,NID)
        DO I=1,NID
          EDEP2(I)=EDEP2(I)+EDEPP(I)
          EDEPP(I)=0.0D0
        ENDDO
        READ(IRD,*) (IDCUT(I),I=1,NID)
        DO I=1,NID
          READ(IRD,'(A)') SPCDIO(I)
          READ(IRD,*) NE(I),EL(I),EU(I),BSE(I),RBSE(I),LLE(I)
          READ(IRD,*) (DITP(I,J),J=1,NE(I))
          DO J=1,NE(I)
            DIT(I,J)=DIT(I,J)+DITP(I,J)
          ENDDO
          READ(IRD,*) (DITP(I,J),J=1,NE(I))
          DO J=1,NE(I)
            DIT2(I,J)=DIT2(I,J)+DITP(I,J)
            DITP(I,J)=0.0D0
          ENDDO
          READ(IRD,*) ((DIPP(I,J,K),J=1,NE(I)),K=1,3)
          DO J=1,NE(I)
            DO K=1,3
              DIP(I,J,K)=DIP(I,J,K)+DIPP(I,J,K)
            ENDDO
          ENDDO
          READ(IRD,*) ((DIPP(I,J,K),J=1,NE(I)),K=1,3)
          DO J=1,NE(I)
            DO K=1,3
              DIP2(I,J,K)=DIP2(I,J,K)+DIPP(I,J,K)
              DIPP(I,J,K)=0.0D0
            ENDDO
          ENDDO
          IF(IDCUT(I).EQ.2) THEN
            READ(IRD,'(A)') SPCFLO(I)
            READ(IRD,*) (FLTP(I,J),J=1,NE(I))
            DO J=1,NE(I)
              FLT(I,J)=FLT(I,J)+FLTP(I,J)
            ENDDO
            READ(IRD,*) (FLTP(I,J),J=1,NE(I))
            DO J=1,NE(I)
              FLT2(I,J)=FLT2(I,J)+FLTP(I,J)
              FLTP(I,J)=0.0D0
            ENDDO
            READ(IRD,*) ((FLPP(I,J,K),J=1,NE(I)),K=1,3)
            DO J=1,NE(I)
              DO K=1,3
                FLP(I,J,K)=FLP(I,J,K)+FLPP(I,J,K)
              ENDDO
            ENDDO
            READ(IRD,*) ((FLPP(I,J,K),J=1,NE(I)),K=1,3)
            DO J=1,NE(I)
              DO K=1,3
                FLP2(I,J,K)=FLP2(I,J,K)+FLPP(I,J,K)
                FLPP(I,J,K)=0.0D0
              ENDDO
            ENDDO
          ENDIF
          READ(IRD,*)
     1      NAGE(I),AGEL(I),AGEU(I),BAGE(I),RBAGE(I),LLAGE(I)
          IF(NAGE(I).GT.0) THEN
            READ(IRD,'(A)') SPCAGE(I)
            READ(IRD,*) (AGEP(I,J),J=1,NAGE(I))
            DO J=1,NAGE(I)
              AGE(I,J)=AGE(I,J)+AGEP(I,J)
            ENDDO
            READ(IRD,*) (AGEP(I,J),J=1,NAGE(I))
            DO J=1,NAGE(I)
              AGE2(I,J)=AGE2(I,J)+AGEP(I,J)
              AGEP(I,J)=0.0D0
            ENDDO
          ENDIF
        ENDDO
      RETURN
C
C  ************  Write simulation results in separate output files.
C
      ENTRY IMDETW(SHN,TSIM,IWR)  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
C
      IF(NID.EQ.0) RETURN
C  ****  Transfer partial counters to global counters.
        DO I=1,NID
          EDEP(I)=EDEP(I)+EDEPP(I)
          EDEP2(I)=EDEP2(I)+EDEPP(I)**2
          EDEPP(I)=0.0D0
          LEDEP(I)=0
          DO J=1,NE(I)
            DIT(I,J)=DIT(I,J)+DITP(I,J)
            DIT2(I,J)=DIT2(I,J)+DITP(I,J)**2
            DITP(I,J)=0.0D0
            LDIT(I,J)=0
            DO K=1,3
              DIP(I,J,K)=DIP(I,J,K)+DIPP(I,J,K)
              DIP2(I,J,K)=DIP2(I,J,K)+DIPP(I,J,K)**2
              DIPP(I,J,K)=0.0D0
              LDIP(I,J,K)=0
            ENDDO
          ENDDO
        ENDDO
C
        DO I=1,NID
          DO J=1,NE(I)
            FLT(I,J)=FLT(I,J)+FLTP(I,J)
            FLT2(I,J)=FLT2(I,J)+FLTP(I,J)**2
            FLTP(I,J)=0.0D0
            LFLT(I,J)=0
            DO K=1,3
              FLP(I,J,K)=FLP(I,J,K)+FLPP(I,J,K)
              FLP2(I,J,K)=FLP2(I,J,K)+FLPP(I,J,K)**2
              FLPP(I,J,K)=0.0D0
              LFLP(I,J,K)=0
            ENDDO
          ENDDO
        ENDDO
C
        DO I=1,NID
          DO J=1,NAGE(I)
            AGE(I,J)=AGE(I,J)+AGEP(I,J)
            AGE2(I,J)=AGE2(I,J)+AGEP(I,J)**2
            AGEP(I,J)=0.0D0
            LAGEA(I,J)=0
          ENDDO
        ENDDO
C
      DF=1.0D0/SHN
C
C  ****  Average energies 'collected' by impact detectors.
C
      WRITE(IWR,3013)
 3013 FORMAT(/3X,'Average incoming energies (impact detectors):')
      DO I=1,NID
        QER=3.0D0*DF*SQRT(ABS(EDEP2(I)-EDEP(I)**2*DF))
        QAV=EDEP(I)*DF
        IF(QER.GT.1.0D-10*ABS(QAV)) THEN
          EFFIC=QAV**2/((QER/3.0D0)**2*TSIM)
        ELSE
          EFFIC=0.0D0
        ENDIF
        WRITE(IWR,3014) I,QAV,QER,EFFIC
      ENDDO
 3014 FORMAT(6X,'Detector #',I2,' ... ',1P,E13.6,' +-',E8.1,' eV',4X,
     1  '(effic. =',E9.2,')')
C
C  ****  Detector efficiencies.
C
      WRITE(IWR,3016)
 3016 FORMAT(/3X,'Spectral areas (impact detectors):')
      DO I=1,NID
        SUM=0.0D0
        ESUM=0.0D0
        DO J=1,NE(I)
          YERR=3.0D0*SQRT(ABS(DIT2(I,J)-DIT(I,J)**2*DF))
          YAV=MAX(DIT(I,J)*DF,1.0D-35)
          YERR=MAX(YERR*DF,1.0D-35)
          SUM=SUM+YAV
          ESUM=ESUM+YERR
        ENDDO
        WRITE(IWR,3017) I,SUM,ESUM
      ENDDO
 3017 FORMAT(6X,'Detector #',I2,' ... ',1P,E13.6,' +-',E8.1)
C
C  ************  Spectra from impact detectors.
C
      DO I=1,NID
        OPEN(9,FILE=SPCDIO(I))
        WRITE(9,9310) I
 9310 FORMAT(
     1  1X,'#  Results from PENMAIN. Output from impact detector #',I3,
     1 /1X,'#  Energy spectra of incident particles.',
     1 /1X,'#  1st column: particle energy (eV).',
     1 /1X,'#  2nd column: probability density (1/(eV*particle)).',
     1 /1X,'#  3rd column: statistical uncertainty, STU (3 sigma).',
     1 /1X,'#  4th-5th columns: electron spectrum and STU (3 sigma).',
     1 /1X,'#  6th-7th columns: photon spectrum and STU (3 sigma).',
     1 /1X,'#  8th-9th columns: positron spectrum and STU (3 sigma).')
        DO J=1,NE(I)
          IF(LLE(I).EQ.1) THEN
            XLOW=EXP(EL(I)+(J-1)*BSE(I))
            XUPP=EXP(EL(I)+J*BSE(I))
            XX=0.5D0*(XUPP+XLOW)
            BINS=XUPP-XLOW
          ELSE
            XX=EL(I)+(J-0.5D0)*BSE(I)
            BINS=BSE(I)
          ENDIF
          YERR=3.0D0*SQRT(ABS(DIT2(I,J)-DIT(I,J)**2*DF))
          YAV=MAX(DIT(I,J)*DF/BINS,1.0D-35)
          YERR=MAX(YERR*DF/BINS,1.0D-35)
C
          YERR1=3.0D0*SQRT(ABS(DIP2(I,J,1)-DIP(I,J,1)**2*DF))
          YAV1=MAX(DIP(I,J,1)*DF/BINS,1.0D-35)
          YERR1=MAX(YERR1*DF/BINS,1.0D-35)
C
          YERR2=3.0D0*SQRT(ABS(DIP2(I,J,2)-DIP(I,J,2)**2*DF))
          YAV2=MAX(DIP(I,J,2)*DF/BINS,1.0D-35)
          YERR2=MAX(YERR2*DF/BINS,1.0D-35)
C
          YERR3=3.0D0*SQRT(ABS(DIP2(I,J,3)-DIP(I,J,3)**2*DF))
          YAV3=MAX(DIP(I,J,3)*DF/BINS,1.0D-35)
          YERR3=MAX(YERR3*DF/BINS,1.0D-35)
          WRITE(9,'(1P,E14.6,5(E14.6,E10.2))')
     1      XX,YAV,YERR,YAV1,YERR1,YAV2,YERR2,YAV3,YERR3
        ENDDO
        CLOSE(9)
      ENDDO
C
C  ************  Distributions of fluence with respect to energy
C  (integrated over the volume of the impact detectors).
C
      DO I=1,NID
        IF(IDCUT(I).EQ.2) THEN
        OPEN(9,FILE=SPCFLO(I))
        WRITE(9,9320) I
 9320 FORMAT(
     1  1X,'#  Results from PENMAIN. Output from impact detector #',I3,
     1 /1X,'#  Fluences integrated over the volume of the detector ',
     1   '(in cm/(eV*particle)).',
     1 /1X,'#  WARNING: fluence of charged particles may extend below',
     1  ' EABS (CSDA tail).',
     1 /1X,'#  1st column: particle energy (eV).',
     1 /1X,'#  2nd column: energy distribution of total fluence.',
     1 /1X,'#  3rd column: statistical uncertainty, STU (3 sigma).',
     1 /1X,'#  4th-5th columns: electron fluence and STU (3 sigma).',
     1 /1X,'#  6th-7th columns: photon fluence and STU (3 sigma).',
     1 /1X,'#  8th-9th columns: positron fluence and STU (3 sigma).')
        DO J=1,NE(I)
          IF(LLE(I).EQ.1) THEN
            XLOW=EXP(EL(I)+(J-1)*BSE(I))
            XUPP=EXP(EL(I)+J*BSE(I))
            XX=0.5D0*(XUPP+XLOW)
            BINS=XUPP-XLOW
          ELSE
            XX=EL(I)+(J-0.5D0)*BSE(I)
            BINS=BSE(I)
          ENDIF
          YERR=3.0D0*SQRT(ABS(FLT2(I,J)-FLT(I,J)**2*DF))
          YAV=MAX(FLT(I,J)*DF/BINS,1.0D-35)
          YERR=MAX(YERR*DF/BINS,1.0D-35)
C
          YERR1=3.0D0*SQRT(ABS(FLP2(I,J,1)-FLP(I,J,1)**2*DF))
          YAV1=MAX(FLP(I,J,1)*DF/BINS,1.0D-35)
          YERR1=MAX(YERR1*DF/BINS,1.0D-35)
C
          YERR2=3.0D0*SQRT(ABS(FLP2(I,J,2)-FLP(I,J,2)**2*DF))
          YAV2=MAX(FLP(I,J,2)*DF/BINS,1.0D-35)
          YERR2=MAX(YERR2*DF/BINS,1.0D-35)
C
          YERR3=3.0D0*SQRT(ABS(FLP2(I,J,3)-FLP(I,J,3)**2*DF))
          YAV3=MAX(FLP(I,J,3)*DF/BINS,1.0D-35)
          YERR3=MAX(YERR3*DF/BINS,1.0D-35)
          WRITE(9,'(1P,E14.6,5(E14.6,E10.2))')
     1      XX,YAV,YERR,YAV1,YERR1,YAV2,YERR2,YAV3,YERR3
        ENDDO
        CLOSE(9)
        ENDIF
      ENDDO
C
C  ************  Time-of-flight distributions from impact detectors.
C
      DO I=1,NID
        IF(NAGE(I).GT.0) THEN
          OPEN(9,FILE=SPCAGE(I))
          WRITE(9,9330) I
 9330 FORMAT(
     1  1X,'#  Results from PENMAIN. Output from impact detector #',I3,
     1 /1X,'#  Age distribution of incident particles.',
     1 /1X,'#  1st column: age (seconds).',
     1 /1X,'#  2nd column: probability density (1/seconds).',
     1 /1X,'#  3rd column: statistical uncertainty, STU (3 sigma).')
          DO J=1,NAGE(I)
            IF(LLAGE(I).EQ.1) THEN
              XLOW=EXP(AGEL(I)+(J-1)*BAGE(I))
              XUPP=EXP(AGEL(I)+J*BAGE(I))
              XX=0.5D0*(XUPP+XLOW)
              BINS=XUPP-XLOW
            ELSE
              XX=AGEL(I)+(J-0.5D0)*BAGE(I)
              BINS=BAGE(I)
            ENDIF
            YERR=3.0D0*SQRT(ABS(AGE2(I,J)-AGE(I,J)**2*DF))
            YAV=MAX(AGE(I,J)*DF/BINS,1.0D-35)
            YERR=MAX(YERR*DF/BINS,1.0D-35)
            WRITE(9,'(1P,E14.6,5(E14.6,E10.2))') XX,YAV,YERR
          ENDDO
          CLOSE(9)
        ENDIF
      ENDDO
C
      RETURN
      END


C  *********************************************************************
C                       SUBROUTINE SENDET
C  *********************************************************************
      SUBROUTINE SENDET(ED,ID)
C
C  Tallies spectra from energy-deposition detectors, writes and loads
C  dump files, accumulates dump files from different runs, and writes
C  results.
C
      USE TRACK_mod
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
      PARAMETER (FSAFE=1.000000001D0)  ! Safety factor.
C
      CHARACTER FNSPC*20
      PARAMETER (NIDM=25,NBEM=1000)
      CHARACTER SPCDEO(NIDM)*20
      SAVE SPCDEO
C
      DIMENSION AUX1(NIDM),AUX2(NIDM,NBEM)
      COMMON/CENDET/EL(NIDM),EU(NIDM),BSE(NIDM),RBSE(NIDM),
     1  EDEP(NIDM),EDEP2(NIDM),DET(NIDM,NBEM),
     1  NE(NIDM),LLE(NIDM),NID
C
      DATA NIDS/0/
      SAVE NIDS
C
C  ************  Deposited energy spectrum.
C                ED includes the particle weight.
C
        EDEP(ID)=EDEP(ID)+ED
        EDEP2(ID)=EDEP2(ID)+ED**2
C
        IF(ED.GT.1.0D-5) THEN
          IF(LLE(ID).EQ.1) THEN
            IE=1.0D0+(LOG(ED)-EL(ID))*RBSE(ID)
          ELSE
            IE=1.0D0+(ED-EL(ID))*RBSE(ID)
          ENDIF
          IF(IE.GT.0.AND.IE.LE.NE(ID)) THEN
            DET(ID,IE)=DET(ID,IE)+1.0D0
          ENDIF
        ENDIF
      RETURN
C
C  ************  Initialise a new detector.
C  Logarithmic (NB.LT.0) or uniform (NB.GT.0) scale.
C
      ENTRY ENDET0(EMIN,EMAX,NB,FNSPC,ID)  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<
C
        IF(ID.LE.NIDS) THEN
          WRITE(26,*) 'SENDET: Detector already defined.',ID
          STOP 'SENDET: Detector cannot be defined.'
        ENDIF
        NIDS=ID
        NID=ID
C
        IF(NID.GT.NIDM) THEN
          WRITE(26,'(3X,''NID = '',I4)') NID
          WRITE(26,*) 'SENDET: Too many detectors.'
          STOP 'SENDET: Too many detectors.'
        ENDIF
C
        SPCDEO(ID)=FNSPC
C
        IF (ABS(NB).GT.NBEM) THEN
          WRITE(26,*) 'SENDET: NB is too large.'
          WRITE(26,*) 'SENDET: Set the parameter NBEM equal to ',
     1      ABS(NB)
          STOP 'SENDET: NB is too large.'
        ENDIF
C
        IF(NB.LT.0) THEN
          LLE(ID)=1
          EL(ID)=LOG(EMIN)
          EU(ID)=LOG(EMAX)
          NE(ID)=-NB
        ELSE
          LLE(ID)=0
          EL(ID)=EMIN
          EU(ID)=EMAX
          NE(ID)=NB
        ENDIF
        BSE(ID)=FSAFE*(EU(ID)-EL(ID))/DBLE(NE(ID))
        RBSE(ID)=1.0D0/BSE(ID)
C
        EDEP(ID)=0.0D0
        EDEP2(ID)=0.0D0
C
        DO J=1,NBEM
          DET(ID,J)=0.0D0
        ENDDO
      RETURN
C
C  ************  DUMP. Write the histograms.
C
      ENTRY ENDETD(IWR)  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
C
        WRITE(IWR,*) NID
        WRITE(IWR,*) (EDEP(I),I=1,NID)
        WRITE(IWR,*) (EDEP2(I),I=1,NID)
        DO I=1,NID
          WRITE(IWR,'(A)') SPCDEO(I)
          WRITE(IWR,*) NE(I),EL(I),EU(I),BSE(I),RBSE(I),LLE(I)
          WRITE(IWR,*) (DET(I,J),J=1,NE(I))
        ENDDO
      RETURN
C
C  ************  RESUME. Read histograms from a previous run.
C
      ENTRY ENDETR(IRD)  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
C
        READ(IRD,*) NID
        READ(IRD,*) (EDEP(I),I=1,NID)
        READ(IRD,*) (EDEP2(I),I=1,NID)
        DO I=1,NID
          READ(IRD,'(A)') SPCDEO(I)
          READ(IRD,*) NE(I),EL(I),EU(I),BSE(I),RBSE(I),LLE(I)
          READ(IRD,*) (DET(I,J),J=1,NE(I))
        ENDDO
      RETURN
C
C  ************  ADD results from previous runs.
C  Use ENDETR to load histograms from one of the runs and then call
C  ENDETA to accumulate results from other runs.
C
      ENTRY ENDETA(IRD)  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
C
        READ(IRD,*) NID
        READ(IRD,*) (AUX1(I),I=1,NID)
        DO I=1,NID
          EDEP(I)=EDEP(I)+AUX1(I)
        ENDDO
        READ(IRD,*) (AUX1(I),I=1,NID)
        DO I=1,NID
          EDEP2(I)=EDEP2(I)+AUX1(I)
        ENDDO
        DO I=1,NID
          READ(IRD,'(A)') SPCDEO(I)
          READ(IRD,*) NE(I),EL(I),EU(I),BSE(I),RBSE(I),LLE(I)
          READ(IRD,*) (AUX2(I,J),J=1,NE(I))
          DO J=1,NE(I)
            DET(I,J)=DET(I,J)+AUX2(I,J)
          ENDDO
        ENDDO
      RETURN
C
C  ************  Write simulation results in separate output files.
C
      ENTRY ENDETW(SHN,TSIM,IWR)  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
C
      IF(NID.EQ.0) RETURN
      DF=1.0D0/SHN
C
C  ****  Average deposited energies in energy-deposition detectors.
C
      WRITE(IWR,3015)
 3015 FORMAT(/3X,'Average deposited energies (energy detectors):')
      DO I=1,NID
        QER=3.0D0*DF*SQRT(ABS(EDEP2(I)-EDEP(I)**2*DF))
        QAV=EDEP(I)*DF
        IF(QER.GT.1.0D-10*ABS(QAV)) THEN
          EFFIC=QAV**2/((QER/3.0D0)**2*TSIM)
        ELSE
          EFFIC=0.0D0
        ENDIF
        WRITE(IWR,3014) I,QAV,QER,EFFIC
      ENDDO
 3014 FORMAT(6X,'Detector #',I2,' ... ',1P,E13.6,' +-',E8.1,' eV',4X,
     1  '(effic. =',E9.2,')')
C
C  ****  Detector efficiencies.
C
      WRITE(IWR,3018)
 3018 FORMAT(/3X,'Spectral areas (energy-deposition detectors):')
      DO I=1,NID
        SUM=0.0D0
        ESUM=0.0D0
        DO J=1,NE(I)
          YAV=DET(I,J)*DF
          YERR=3.0D0*SQRT(ABS(YAV*(1.0D0-YAV))*DF)
          SUM=SUM+YAV
          ESUM=ESUM+YERR
        ENDDO
        WRITE(IWR,3017) I,SUM,ESUM
      ENDDO
 3017 FORMAT(6X,'Detector #',I2,' ... ',1P,E13.6,' +-',E8.1)
C
C  ************  Spectra from energy-deposition detectors.
C
      DO I=1,NID
        OPEN(9,FILE=SPCDEO(I))
        WRITE(9,9310) I
 9310 FORMAT(
     1  1X,'#  Results from PENMAIN. Output from energy-deposition',
     1  ' detector #',I3,
     1 /1X,'#  Deposited energy spectrum.',
     1 /1X,'#  WARNING: May be strongly biased if interaction ',
     1  'forcing is used!',
     1 /1X,'#  1st column: deposited energy (eV).',
     1 /1X,'#  2nd column: probability density (1/(eV*particle)).',
     1 /1X,'#  3rd column: statistical uncertainty (3 sigma).')
        DO J=1,NE(I)
          IF(LLE(I).EQ.1) THEN
            XLOW=EXP(EL(I)+(J-1)*BSE(I))
            XUPP=EXP(EL(I)+J*BSE(I))
            XX=0.5D0*(XUPP+XLOW)
            BINS=XUPP-XLOW
          ELSE
            XX=EL(I)+(J-0.5D0)*BSE(I)
            BINS=BSE(I)
          ENDIF
          YAV=DET(I,J)*DF
          YERR=3.0D0*SQRT(ABS(YAV*(1.0D0-YAV))*DF)
          YAV=YAV/BINS
          YERR=YERR/BINS
          WRITE(9,'(1P,2E14.6,E10.2)') XX,MAX(YAV,1.0D-35),YERR
        ENDDO
        CLOSE(9)
      ENDDO
C
      RETURN
      END


C  *********************************************************************
C                       SUBROUTINE SDOSE
C  *********************************************************************
      SUBROUTINE SDOSE(DEP,XD,YD,ZD,MATC,N)
C
C  Tallies the dose distribution within the dose box, writes and loads
C  dump files, accumulates dump files from different runs, and writes
C  results.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
      PARAMETER (PI=3.1415926535897932D0)
      PARAMETER (FSAFE=1.000000001D0)  ! Safety factor.
C  ----  Dose distribution (up to NDXM*NDYM*NDZM volume bins).
      PARAMETER (NDXM=201,NDYM=201,NDZM=201)
      COMMON/CDOSE1/DOSE(NDXM,NDYM,NDZM),DOSE2(NDXM,NDYM,NDZM),
     1  DOSEP(NDXM,NDYM,NDZM),LDOSE(NDXM,NDYM,NDZM),KDOSE
      COMMON/CDOSE2/DDOSE(NDZM),DDOSE2(NDZM),DDOSEP(NDZM),LDDOSE(NDZM)
      COMMON/CDOSE3/DXL(3),DXU(3),BDOSE(3),RBDOSE(3),NDB(3)
C
      PARAMETER (NCS=5,RNCS=1.0D0/5.0D0)
      COMMON/CDOSE4/VMASS(NDXM,NDYM,NDZM)
C
C  ************  Dose distribution.
C
      IF(KDOSE.EQ.1) THEN  ! Box
        IF(ZD.GT.DXL(3).AND.ZD.LT.DXU(3)) THEN
          I3=1.0D0+(ZD-DXL(3))*RBDOSE(3)
          IF(N.NE.LDDOSE(I3)) THEN
            DDOSE(I3)=DDOSE(I3)+DDOSEP(I3)
            DDOSE2(I3)=DDOSE2(I3)+DDOSEP(I3)**2
            DDOSEP(I3)=DEP*DBLE(MATC)  ! Modified.
            LDDOSE(I3)=N
          ELSE
            DDOSEP(I3)=DDOSEP(I3)+DEP*DBLE(MATC)  ! Modified.
          ENDIF
C
          IF((XD.GT.DXL(1).AND.XD.LT.DXU(1)).AND.
     1       (YD.GT.DXL(2).AND.YD.LT.DXU(2))) THEN
            I1=1.0D0+(XD-DXL(1))*RBDOSE(1)
            I2=1.0D0+(YD-DXL(2))*RBDOSE(2)
            IF(N.NE.LDOSE(I1,I2,I3)) THEN
              DOSE(I1,I2,I3)=DOSE(I1,I2,I3)+DOSEP(I1,I2,I3)
              DOSE2(I1,I2,I3)=DOSE2(I1,I2,I3)+DOSEP(I1,I2,I3)**2
              DOSEP(I1,I2,I3)=DEP
              LDOSE(I1,I2,I3)=N
            ELSE
              DOSEP(I1,I2,I3)=DOSEP(I1,I2,I3)+DEP
            ENDIF
          ENDIF
        ENDIF
C
C
      ELSE IF(KDOSE.EQ.2) THEN  ! Cylinder.
        IF(ZD.GT.DXL(3).AND.ZD.LT.DXU(3)) THEN
          I3=1.0D0+(ZD-DXL(3))*RBDOSE(3)
          IF(N.NE.LDDOSE(I3)) THEN
            DDOSE(I3)=DDOSE(I3)+DDOSEP(I3)
            DDOSE2(I3)=DDOSE2(I3)+DDOSEP(I3)**2
            DDOSEP(I3)=DEP*DBLE(MATC)  ! Modified.
            LDDOSE(I3)=N
          ELSE
            DDOSEP(I3)=DDOSEP(I3)+DEP*DBLE(MATC)  ! Modified.
          ENDIF
C
          RD=SQRT(XD*XD+YD*YD)
          IF(RD.LT.DXU(1)) THEN
            I1=1.0D0+RD*RBDOSE(1)
            I2=1
            IF(N.NE.LDOSE(I1,I2,I3)) THEN
              DOSE(I1,I2,I3)=DOSE(I1,I2,I3)+DOSEP(I1,I2,I3)
              DOSE2(I1,I2,I3)=DOSE2(I1,I2,I3)+DOSEP(I1,I2,I3)**2
              DOSEP(I1,I2,I3)=DEP
              LDOSE(I1,I2,I3)=N
            ELSE
              DOSEP(I1,I2,I3)=DOSEP(I1,I2,I3)+DEP
            ENDIF
          ENDIF
        ENDIF
C
C
      ELSE  ! Sphere.
        RD=SQRT(XD*XD+YD*YD+ZD*ZD)
        IF(RD.LT.DXU(1)) THEN
          I1=1.0D0+RD*RBDOSE(1)
          I2=1
          I3=1
          IF(N.NE.LDOSE(I1,I2,I3)) THEN
            DOSE(I1,I2,I3)=DOSE(I1,I2,I3)+DOSEP(I1,I2,I3)
            DOSE2(I1,I2,I3)=DOSE2(I1,I2,I3)+DOSEP(I1,I2,I3)**2
            DOSEP(I1,I2,I3)=DEP
            LDOSE(I1,I2,I3)=N
          ELSE
            DOSEP(I1,I2,I3)=DOSEP(I1,I2,I3)+DEP
          ENDIF
        ENDIF
      ENDIF
      RETURN
C
C  ************  Initialisation:  voxel average masses, dose counters.
C
      ENTRY DOSE0(XL,XU,YL,YU,ZL,ZU,NBX,NBY,NBZ,IDOSE)  !<<<<<<<<<<<<<<<
C
        IF(IDOSE.LT.1.OR.IDOSE.GT.3) THEN
          WRITE(26,'(''IDOSE ='',I6)') IDOSE
          WRITE(26,'(''IDOSE should be 1, 2, or 3.'')')
          STOP 'SDOSE: IDOSE should be 1, 2, or 3.'
        ENDIF
        KDOSE=IDOSE
        IF(NBX.LT.0.OR.NBX.GT.NDXM) THEN
          WRITE(26,'(''NBX ='',I6)') NBX
          WRITE(26,'(''NBX must be .GT.0. and .LT.'',I4)') NDXM
          WRITE(26,*) 'Increase the value of the parameter NDXM.'
          STOP 'SDOSE: NBX must be .GT.0. and .LE.NDXM'
        ENDIF
        IF(NBY.LT.0.OR.NBY.GT.NDYM) THEN
          WRITE(26,'(''NBY ='',I6)') NBY
          WRITE(26,'(''NBY must be .GT.0. and .LT.'',I4)') NDYM
          WRITE(26,*) 'Increase the value of the parameter NDYM.'
          STOP 'SDOSE: NBY must be .GT.0. and .LE.NDYM'
        ENDIF
        IF(NBZ.LT.0.OR.NBZ.GT.NDZM) THEN
          WRITE(26,'(''NBZ ='',I6)') NBZ
          WRITE(26,'(''NBZ must be .GT.0. and .LT.'',I4)') NDZM
          WRITE(26,*) 'Increase the value of the parameter NDZM.'
          STOP 'SDOSE: NBZ must be .GT.0. and .LE.NDZM'
        ENDIF
C
        DXL(1)=XL
        DXU(1)=XU
        DXL(2)=YL
        DXU(2)=YU
        DXL(3)=ZL
        DXU(3)=ZU
        IF(KDOSE.EQ.1) THEN
          NDB(1)=2*(NBX/2)+1
          NDB(2)=2*(NBY/2)+1
          NDB(3)=2*(NBZ/2)+1
        ELSE IF(KDOSE.EQ.2) THEN
          DXL(1)=0.0D0
          DXL(2)=0.0D0
          DXU(2)=1.0D0
          NDB(1)=NBX
          NDB(2)=1
          NDB(3)=NBZ
        ELSE
          DXL(1)=0.0D0
          DXL(2)=0.0D0
          DXU(2)=1.0D0
          DXL(3)=0.0D0
          DXU(3)=1.0D0
          NDB(1)=NBX
          NDB(2)=1
          NDB(3)=1
        ENDIF
C
        DO I=1,3
          BDOSE(I)=FSAFE*(DXU(I)-DXL(I))/DBLE(NDB(I))
          IF(ABS(BDOSE(I)).LT.1.0D-35) THEN
            RBDOSE(I)=1.0D35
          ELSE
            RBDOSE(I)=1.0D0/BDOSE(I)
          ENDIF
        ENDDO
C
C  ****  Dose counters.
C
        DO K=1,NDZM
          DO J=1,NDYM
            DO I=1,NDXM
              DOSE(I,J,K)=0.0D0
              DOSE2(I,J,K)=0.0D0
              DOSEP(I,J,K)=0.0D0
              LDOSE(I,J,K)=0
              VMASS(I,J,K)=0.0D0
            ENDDO
          ENDDO
          DDOSE(K)=0.0D0
          DDOSE2(K)=0.0D0
          DDOSEP(K)=0.0D0
          LDDOSE(K)=0
        ENDDO
      RETURN
C
C  ************  DUMP. Write the histograms.
C
      ENTRY DOSED(IWR)  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
C
C  ************  Transfer partial counters to global counters.
        DO I3=1,NDB(3)
        DO I2=1,NDB(2)
        DO I1=1,NDB(1)
          DOSE(I1,I2,I3)=DOSE(I1,I2,I3)+DOSEP(I1,I2,I3)
          DOSE2(I1,I2,I3)=DOSE2(I1,I2,I3)+DOSEP(I1,I2,I3)**2
          DOSEP(I1,I2,I3)=0.0D0
          LDOSE(I1,I2,I3)=0
        ENDDO
        ENDDO
          DDOSE(I3)=DDOSE(I3)+DDOSEP(I3)
          DDOSE2(I3)=DDOSE2(I3)+DDOSEP(I3)**2
          DDOSEP(I3)=0.0D0
          LDDOSE(I3)=0
        ENDDO
C  ****  Write parameters and counters.
        WRITE(IWR,*) (NDB(I),I=1,3)
        WRITE(IWR,*) (DXL(I),I=1,3)
        WRITE(IWR,*) (DXU(I),I=1,3)
        WRITE(IWR,*) (BDOSE(I),I=1,3),(RBDOSE(I),I=1,3)
        WRITE(IWR,*) KDOSE
        WRITE(IWR,*)
     1    (((VMASS(I1,I2,I3),I3=1,NDB(3)),I2=1,NDB(2)),I1=1,NDB(1))
        WRITE(IWR,*)
     1    (((DOSE(I1,I2,I3),I3=1,NDB(3)),I2=1,NDB(2)),I1=1,NDB(1))
        WRITE(IWR,*)
     1    (((DOSE2(I1,I2,I3),I3=1,NDB(3)),I2=1,NDB(2)),I1=1,NDB(1))
        WRITE(IWR,*) (DDOSE(I3),I3=1,NDB(3))
        WRITE(IWR,*) (DDOSE2(I3),I3=1,NDB(3))
      RETURN
C
C  ************  RESUME. Read histograms from a previous run.
C
      ENTRY DOSER(IRD)  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
C
        READ(IRD,*) (NDB(I),I=1,3)
        READ(IRD,*) (DXL(I),I=1,3)
        READ(IRD,*) (DXU(I),I=1,3)
        READ(IRD,*) (BDOSE(I),I=1,3),(RBDOSE(I),I=1,3)
        READ(IRD,*) KDOSE
        READ(IRD,*)
     1    (((VMASS(I1,I2,I3),I3=1,NDB(3)),I2=1,NDB(2)),I1=1,NDB(1))
        READ(IRD,*)
     1    (((DOSE(I1,I2,I3),I3=1,NDB(3)),I2=1,NDB(2)),I1=1,NDB(1))
        READ(IRD,*)
     1    (((DOSE2(I1,I2,I3),I3=1,NDB(3)),I2=1,NDB(2)),I1=1,NDB(1))
        READ(IRD,*) (DDOSE(I3),I3=1,NDB(3))
        READ(IRD,*) (DDOSE2(I3),I3=1,NDB(3))
      RETURN
C
C  ************  ADD results from previous runs.
C  Use DOSER to load histograms from one of the runs and then call
C  DOSEA to accumulate results from other runs.
C
      ENTRY DOSEA(IRD)  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
C
        READ(IRD,*) (NDB(I),I=1,3)
        READ(IRD,*) (DXL(I),I=1,3)
        READ(IRD,*) (DXU(I),I=1,3)
        READ(IRD,*) (BDOSE(I),I=1,3),(RBDOSE(I),I=1,3)
        READ(IRD,*) KDOSE
        READ(IRD,*)
     1    (((VMASS(I1,I2,I3),I3=1,NDB(3)),I2=1,NDB(2)),I1=1,NDB(1))
        READ(IRD,*)
     1    (((DOSEP(I1,I2,I3),I3=1,NDB(3)),I2=1,NDB(2)),I1=1,NDB(1))
        DO I1=1,NDB(1)
          DO I2=1,NDB(2)
            DO I3=1,NDB(3)
              DOSE(I1,I2,I3)=DOSE(I1,I2,I3)+DOSEP(I1,I2,I3)
            ENDDO
          ENDDO
        ENDDO
        READ(IRD,*)
     1    (((DOSEP(I1,I2,I3),I3=1,NDB(3)),I2=1,NDB(2)),I1=1,NDB(1))
        DO I1=1,NDB(1)
          DO I2=1,NDB(2)
            DO I3=1,NDB(3)
              DOSE2(I1,I2,I3)=DOSE2(I1,I2,I3)+DOSEP(I1,I2,I3)
              DOSEP(I1,I2,I3)=0.0D0
            ENDDO
          ENDDO
        ENDDO
        READ(IRD,*) (DDOSEP(I3),I3=1,NDB(3))
        DO I3=1,NDB(3)
          DDOSE(I3)=DDOSE(I3)+DDOSEP(I3)
        ENDDO
        READ(IRD,*) (DDOSEP(I3),I3=1,NDB(3))
        DO I3=1,NDB(3)
          DDOSE2(I3)=DDOSE2(I3)+DDOSEP(I3)
          DDOSEP(I3)=0.0D0
        ENDDO
      RETURN
C
C  ************  Write simulation results in output files.
C
      ENTRY DOSEW(SHN,TSIM)  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        DF=1.0D0/SHN
C
C  ************  Transfer partial counters to global counters.
        DO I3=1,NDB(3)
        DO I2=1,NDB(2)
        DO I1=1,NDB(1)
          DOSE(I1,I2,I3)=DOSE(I1,I2,I3)+DOSEP(I1,I2,I3)
          DOSE2(I1,I2,I3)=DOSE2(I1,I2,I3)+DOSEP(I1,I2,I3)**2
          DOSEP(I1,I2,I3)=0.0D0
          LDOSE(I1,I2,I3)=0
        ENDDO
        ENDDO
          DDOSE(I3)=DDOSE(I3)+DDOSEP(I3)
          DDOSE2(I3)=DDOSE2(I3)+DDOSEP(I3)**2
          DDOSEP(I3)=0.0D0
          LDDOSE(I3)=0
        ENDDO
C
        DMAX=0.0D0
        I1M=1
        I2M=1
        I3M=1
        DO I1=1,NDB(1)
          DO I2=1,NDB(2)
            DO I3=1,NDB(3)
              IF(DOSE(I1,I2,I3)*VMASS(I1,I2,I3).GT.DMAX) THEN
                I1M=I1
                I2M=I2
                I3M=I3
                DMAX=DOSE(I1,I2,I3)*VMASS(I1,I2,I3)
              ENDIF
            ENDDO
          ENDDO
        ENDDO
C
        QAV=DOSE(I1M,I2M,I3M)
        QAV2=DOSE2(I1M,I2M,I3M)
        QER=3.0D0*SQRT(ABS(QAV2-QAV**2*DF))
        QAV=QAV*DF*VMASS(I1M,I2M,I3M)
        QER=QER*DF*VMASS(I1M,I2M,I3M)
        IF(QER.GT.1.0D-10*ABS(QAV)) THEN
          EFFIC=QAV**2/((QER/3.0D0)**2*TSIM)
        ELSE
          EFFIC=0.0D0
        ENDIF
        WRITE(27,3020) QAV,QER,EFFIC
 3020   FORMAT(/6X,'Maximum dose ... ',1P,E13.6,' +-',E8.1,' eV/g',2X,
     1    '(effic. =',E9.2,')')
C
      IF(KDOSE.EQ.1) THEN  ! Box.
C
C  ****  Depth-dose distribution.
C
        OPEN(9,FILE='depth-dose.dat')
        WRITE(9,9410)
 9410   FORMAT(
     1     1X,'#  Results from PENMAIN. Depth-dose distribution.',
     1    /1X,'#  (integrated over X and Y within the volume of the ',
     1      'material system).',
     1    /1X,'#  1st column: z coordinate (cm).',
     1    /1X,'#  2nd column: depth-dose (eV/(g/cm**2)).',
     1    /1X,'#  3rd column: statistical uncertainty (3 sigma).')
        WRITE(9,9411)
 9411   FORMAT(1X,'#',
     1    /1X,'#  NOTE: The calculated dose distribution is correct',
     1       ' only when the ',
     1    /1X,'#',8X,'Z bins have uniform mass density.',/1X,'#')
        DO I3=1,NDB(3)
          ZZ=DXL(3)+(I3-0.5D0)*BDOSE(3)
          YAV=DDOSE(I3)
          YAV2=DDOSE2(I3)
          YERR=3.0D0*SQRT(ABS(YAV2-YAV**2*DF))
          YAV=YAV*DF*RBDOSE(3)
          YERR=YERR*DF*RBDOSE(3)
          WRITE(9,'(1P,E14.6,E14.6,E10.2)')
     1      ZZ,MAX(YAV,1.0D-35),MAX(YERR,1.0D-35)
        ENDDO
        CLOSE(9)
C
C  ****  3D dose map.
C
        OPEN(9,FILE='3d-dose-map.dat')
        WRITE(9,9420)
 9420   FORMAT(1X,'#  Results from PENMAIN. 3D dose distribution.')
        WRITE(9,9421) DXL(1),DXU(1)
 9421   FORMAT(1X,'#  Dose-map box:  XL = ',1P,E13.6,
     1    ' cm,  XU = ',E13.6,' cm')
        WRITE(9,9422) DXL(2),DXU(2)
 9422   FORMAT(1X,'#',17X,'YL = ',1P,E13.6,' cm,  YU = ',E13.6,' cm')
        WRITE(9,9423) DXL(3),DXU(3)
 9423   FORMAT(1X,'#',17X,'ZL = ',1P,E13.6,' cm,  ZU = ',E13.6,' cm')
        WRITE(9,9424) NDB(1),NDB(2),NDB(3)
 9424   FORMAT(1X,'#  Numbers of bins:     NBX =',I4,', NBY =',I4,
     1        ', NBZ =',I4,/1X,'#')
        WRITE(9,9425)
 9425   FORMAT(1X,'#  columns 1 to 3: coordinates X,Y,Z of the bin',
     1    ' centres.',/1X,'#  4th column: dose (eV/g).',
     1    /1X,'#  5th column: statistical uncertainty (3 sigma).',
     1    /1X,'#  columns 6 to 8: bin indices IX,IY,IZ.')
        DO I3=1,NDB(3)
          ZZ=DXL(3)+(I3-0.5D0)*BDOSE(3)
          DO I1=1,NDB(1)
            XX=DXL(1)+(I1-0.5D0)*BDOSE(1)
            DO I2=1,NDB(2)
              YY=DXL(2)+(I2-0.5D0)*BDOSE(2)
              YAV=DOSE(I1,I2,I3)
              YAV2=DOSE2(I1,I2,I3)
              YERR=3.0D0*SQRT(ABS(YAV2-YAV**2*DF))
              YAV=YAV*DF*VMASS(I1,I2,I3)
              YERR=YERR*DF*VMASS(I1,I2,I3)
              WRITE(9,9426) XX,YY,ZZ,
     1          MAX(YAV,1.0D-35),MAX(YERR,1.0D-35),I1,I2,I3
            ENDDO
            WRITE(9,*) '   '
          ENDDO
          WRITE(9,*) '   '
        ENDDO
 9426   FORMAT(1P,3E11.3,E14.6,E10.2,3I4)
        CLOSE(9)
C
C  ****  Dose distributions at the central axes.
C
        I1C=(NDB(1)/2)+1
        I2C=(NDB(2)/2)+1
        I3C=(NDB(3)/2)+1
C
        IF(NDB(1).GT.1) THEN
          OPEN(9,FILE='x-dose.dat')
            WRITE(9,9440)
 9440     FORMAT(
     1       1X,'#  Results from PENMAIN.',
     1      /1X,'#  Dose distribution along the central X axis.',
     1      /1X,'#  1st column: x (cm).',
     1      /1X,'#  2nd column: dose (eV/g).',
     1      /1X,'#  3rd column: statistical uncertainty (3 sigma).')
          DO I1=1,NDB(1)
            XYZ=DXL(1)+(I1-0.5D0)*BDOSE(1)
            YAV=DOSE(I1,I2C,I3C)
            YAV2=DOSE2(I1,I2C,I3C)
            YERR=3.0D0*SQRT(ABS(YAV2-YAV**2*DF))
            YAV=YAV*DF*VMASS(I1,I2C,I3C)
            YERR=YERR*DF*VMASS(I1,I2C,I3C)
            WRITE(9,'(1P,E14.6,E14.6,E10.2)')
     1        XYZ,MAX(YAV,1.0D-35),MAX(YERR,1.0D-35)
          ENDDO
          CLOSE(9)
        ENDIF
C
        IF(NDB(2).GT.1) THEN
          OPEN(9,FILE='y-dose.dat')
            WRITE(9,9450)
 9450     FORMAT(
     1       1X,'#  Results from PENMAIN.',
     1      /1X,'#  Dose distribution along the central Y axis.',
     1      /1X,'#  1st column: y (cm).',
     1      /1X,'#  2nd column: dose (eV/g).',
     1      /1X,'#  3rd column: statistical uncertainty (3 sigma).')
          WRITE(9,9411)
          DO I2=1,NDB(2)
            XYZ=DXL(2)+(I2-0.5D0)*BDOSE(2)
            YAV=DOSE(I1C,I2,I3C)
            YAV2=DOSE2(I1C,I2,I3C)
            YERR=3.0D0*SQRT(ABS(YAV2-YAV**2*DF))
            YAV=YAV*DF*VMASS(I1C,I2,I3C)
            YERR=YERR*DF*VMASS(I1C,I2,I3C)
            WRITE(9,'(1P,E14.6,E14.6,E10.2)')
     1        XYZ,MAX(YAV,1.0D-35),MAX(YERR,1.0D-35)
          ENDDO
          CLOSE(9)
        ENDIF
C
        IF(NDB(3).GT.1) THEN
          OPEN(9,FILE='z-dose.dat')
            WRITE(9,9460)
 9460     FORMAT(
     1       1X,'#  Results from PENMAIN.',
     1      /1X,'#  Dose distribution along the central Z axis.',
     1      /1X,'#  1st column: z (cm).',
     1      /1X,'#  2nd column: dose (eV/g).',
     1      /1X,'#  3rd column: statistical uncertainty (3 sigma).')
          DO I3=1,NDB(3)
            XYZ=DXL(3)+(I3-0.5D0)*BDOSE(3)
            YAV=DOSE(I1C,I2C,I3)
            YAV2=DOSE2(I1C,I2C,I3)
            YERR=3.0D0*SQRT(ABS(YAV2-YAV**2*DF))
            YAV=YAV*DF*VMASS(I1C,I2C,I3)
            YERR=YERR*DF*VMASS(I1C,I2C,I3)
            WRITE(9,'(1P,E14.6,E14.6,E10.2)')
     1        XYZ,MAX(YAV,1.0D-35),MAX(YERR,1.0D-35)
          ENDDO
          CLOSE(9)
        ENDIF
C
C
      ELSE IF(KDOSE.EQ.2) THEN  ! Cylinder.
C
C  ****  2D dose map.
C
        OPEN(9,FILE='2d-dose-map.dat')
        WRITE(9,9520)
 9520   FORMAT(1X,'#  Results from PENMAIN. Dose distribution.')
        WRITE(9,9521) DXU(1)
 9521   FORMAT(1X,'#  Dose-map cylinder:         RU = ',1P,E13.6,' cm')
        WRITE(9,9523) DXL(3),DXU(3)
 9523   FORMAT(1X,'#',5X,'ZL = ',1P,E13.6,' cm,  ZU = ',E13.6,' cm')
        WRITE(9,9524) NDB(1),NDB(3)
 9524   FORMAT(1X,'#  Numbers of bins:     NBR =',I4,', NBZ =',I4,
     1    /1X,'#')
        WRITE(9,9525)
 9525   FORMAT(1X,'#  columns 1 and 2: coordinates R,Z of the bin',
     1    ' centres',/1X,'#  3rd column: dose (eV/g).',
     1    /1X,'#  4th column: statistical uncertainty (3 sigma).')
        DO I1=1,NDB(1)
          RR=(I1-0.5D0)*BDOSE(1)
          DO I3=1,NDB(3)
            ZZ=DXL(3)+(I3-0.5D0)*BDOSE(3)
            YAV=DOSE(I1,1,I3)
            YAV2=DOSE2(I1,1,I3)
            YERR=3.0D0*SQRT(ABS(YAV2-YAV**2*DF))
            YAV=YAV*DF*VMASS(I1,1,I3)
            YERR=YERR*DF*VMASS(I1,1,I3)
            WRITE(9,9526) RR,ZZ,MAX(YAV,1.0D-35),MAX(YERR,1.0D-35)
          ENDDO
          WRITE(9,*) '   '
        ENDDO
 9526   FORMAT(1P,2E14.6,E14.6,E10.2)
        CLOSE(9)
C
        OPEN(9,FILE='depth-dose.dat')
        WRITE(9,9410)
        WRITE(9,9411)
        DO I3=1,NDB(3)
          ZZ=DXL(3)+(I3-0.5D0)*BDOSE(3)
          YAV=DDOSE(I3)
          YAV2=DDOSE2(I3)
          YERR=3.0D0*SQRT(ABS(YAV2-YAV**2*DF))
          YAV=YAV*DF*RBDOSE(3)
          YERR=YERR*DF*RBDOSE(3)
          WRITE(9,'(1P,E14.6,E14.6,E10.2)')
     1      ZZ,MAX(YAV,1.0D-35),MAX(YERR,1.0D-35)
        ENDDO
        CLOSE(9)
C
        OPEN(9,FILE='z-dose.dat')
        WRITE(9,9460)
        DO I3=1,NDB(3)
          XYZ=DXL(3)+(I3-0.5D0)*BDOSE(3)
          YAV=DOSE(1,1,I3)
          YAV2=DOSE2(1,1,I3)
          YERR=3.0D0*SQRT(ABS(YAV2-YAV**2*DF))
          YAV=YAV*DF*VMASS(1,1,I3)
          YERR=YERR*DF*VMASS(1,1,I3)
          WRITE(9,'(1P,E14.6,E14.6,E10.2)')
     1      XYZ,MAX(YAV,1.0D-35),MAX(YERR,1.0D-35)
        ENDDO
        CLOSE(9)
C
C
      ELSE  ! Sphere.
C
C  ****  Radial dose distribution.
C
        OPEN(9,FILE='radial-dose.dat')
        WRITE(9,9620)
 9620   FORMAT(1X,'#  Results from PENMAIN. Dose distribution.')
        WRITE(9,9621) DXU(1)
 9621   FORMAT(1X,'#  Dose-map sphere:  RU = ',1P,E13.6,' cm')
        WRITE(9,9624) NDB(1)
 9624   FORMAT(1X,'#  Number of bins:  NBR =',I4,/1X,'#')
        WRITE(9,9625)
 9625   FORMAT(1X,'#  column 1: radius R of the bin centres.',
     1    /1X,'#  2nd column: absorbed dose (eV/g).',
     1    /1X,'#  3rd column: statistical uncertainty (3 sigma).',
     1  /1X,'#  4th column: deposited energy per unit radius (eV/cm).',
     1    /1X,'#  5th column: statistical uncertainty (3 sigma).')
        DO I1=1,NDB(1)
          RR=(I1-0.5D0)*BDOSE(1)
          YAV=DOSE(I1,1,1)
          YAV2=DOSE2(I1,1,1)
          YERR=3.0D0*SQRT(ABS(YAV2-YAV**2*DF))
          YAV=YAV*DF*VMASS(I1,1,1)
          YERR=YERR*DF*VMASS(I1,1,1)
          ZAV=DOSE(I1,1,1)
          ZAV2=DOSE2(I1,1,1)
          ZERR=3.0D0*SQRT(ABS(ZAV2-ZAV**2*DF))
          ZAV=ZAV*DF*RBDOSE(1)
          ZERR=ZERR*DF*RBDOSE(1)
          WRITE(9,9626) RR,YAV,YERR,ZAV,ZERR
        ENDDO
 9626   FORMAT(1P,E14.6,2(E14.6,E10.2))
        CLOSE(9)
      ENDIF
        RETURN
      END
