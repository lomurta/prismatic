C  *********************************************************************
C                       SUBROUTINE SOURCE
C  *********************************************************************
      SUBROUTINE SOURCE(METAST)
C
C  >>>>>> Example of source simulation routine.
C
C  The program PENMAIN allows the use of arbitrary sources of primary
C  radiations. When KPARP (type of primary particle) is set to 0, the
C  initial states of primary particles are assumed to be delivered by
C  the present subroutine, which may be edited by the user to define and
C  simulate the primary radiation source under consideration.
C
C  It is assumed that in the first call to subroutine SOURCE (when the
C  value of the parameter IDEF is equal to 0) we define _all_ the source
C  parameters that are required and have not been set previously by
C  PENMAIN. Definitions in this first call supersede any information
C  that has been read previously by PENMAIN. The value of EPMAX (highest
C  energy of the primary particles) should be defined in the first call
C  to SOURCE; this parameter is needed for the initialisation of
C  PENELOPE.
C
C  Each subsequent call to SOURCE sets the initial states of a primary
C  particle, or of a set of particles released simultaneously in a
C  single primary event (e.g., in a cascade de-excitation of a radio-
C  active nucleus). The state variables E, X,Y,Z, U,V,W, WGHT, KPAR and
C  ILB(1:5) of primary particles are stored in memory. One of the
C  particles is stored directly in common TRACK, the others (if any) are
C  sent to the secondary stack (by calling subroutine STORES of the
C  PENELOPE package). The particles sent to the secondary stack must be
C  marked with ILB(1)=1, to be treated as proper primary particles by
C  subroutine SHOWER (their initial energies will not be subtracted
C  from the locally absorbed energy).
C
C  In the case of polarised photons, the polarization state (Stokes
C  parameters) and the IPOL flag are loaded in common STOKES. Note that
C  subroutine STORES gets information through the common blocks TRACK
C  and STOKES. Hence, the physics variables in these commons (i.e., all
C  quantities except IBODY and MAT) should be defined before calling
C  STORES.
C
C  A radioactive nucleus may reach a metaestable state (i.e., a state
C  with a life time longer than the detector resolution time) in the
C  course of its de-excitation cascade. Such states are treated as
C  effective halts of the decay cascade, and radiations emitted between
C  successive halts are considered to be released simultaneously. When a
C  metastable state is reached, the output value of METAST is set to 1
C  and control is returned to the main program. At the next call to
C  subroutine SOURCE, the simulation of the cascade continues down to
C  the following metastable state or to the ground level. The output
C  value METAST = 0 indicates that the nucleus has reached its ground
C  level. This operation scheme allows simulating the response of real
C  detectors, which may give several output counts as the result of the
C  decay of a single nucleus through a metastable state.
C
C  Input/output:
C    METAST ... output flag:
C               = 0, the de-excitation cascade has been completed.
C               = 1, a metastable state has been found.
C
      USE TRACK_mod
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
      LOGICAL LSPEC,LGPOL,LPSF,LEXSRC,LEXBD,LSCONE
*     PARAMETER (PI=3.1415926535897932D0, TWOPI=2.0D0*PI)
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
      PARAMETER (NB=5000)
      COMMON/CSOUR3/SX0,SY0,SZ0,SSX,SSY,SSZ,IXSBOD(NB),LEXSRC,LEXBD
C
      DIMENSION ILBS(5)
      EXTERNAL RAND
C  ****  Save the position of the decaying nucleus.
      SAVE X0,Y0,Z0
C
C  First call only.  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
C
C  ************  Initialisation. Definition of the radiation source.
C  Save IDEF and all parameters needed to set the initial states of
C  primary particles.
C
      SAVE E1,E2,TAV1,IDEF
      DATA IDEF/0/
C
C  ****  Initialisation: Define here the source characteristics.
C
      IF(IDEF.EQ.0) THEN  ! Example: 60-Co source, two gamma rays.
        E1=1.17324D6  ! Energy of 1st transition.
        E2=1.33251D6  ! Energy of 2nd transition.
        TAV1=0.713D-12  ! Mean life of the intermediate level.
        EPMAX=MAX(E1,E2)+1.0D0  ! Note: for positrons add 5.12D5.
        IDEF=1
        RETURN
      ENDIF
C  <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
C
C
C  ************  Set the initial states of particles released in a new
C                nuclear decay.
C
      IF(METAST.NE.0) GO TO 200  ! Nucleus in a metastable level.
C
C  ----  Positition of the decaying nucleus (same as in PENMAIN).
C
      IF(LEXSRC) THEN  ! Extended source.
        IF(LEXBD) THEN  ! Active bodies have been defined.
          NTRIAL=0
   10     CONTINUE
          X=SX0+(RAND(1.0D0)-0.5D0)*SSX
          Y=SY0+(RAND(2.0D0)-0.5D0)*SSY
          Z=SZ0+(RAND(3.0D0)-0.5D0)*SSZ
          CALL LOCATE
          NTRIAL=NTRIAL+1
          IF(NTRIAL.GT.100) THEN
            WRITE(26,'(3X,''WARNING: the sampling of initial '',
     1        ''positions may be very inefficient.'')')
            WRITE(6,'(3X,''WARNING: the sampling of initial '',
     1        ''positions may be very inefficient.'')')
          ENDIF
          IF(IXSBOD(IBODY).EQ.0) GO TO 10
        ELSE
          X=SX0+(RAND(1.0D0)-0.5D0)*SSX
          Y=SY0+(RAND(2.0D0)-0.5D0)*SSY
          Z=SZ0+(RAND(3.0D0)-0.5D0)*SSZ
        ENDIF
      ELSE  ! Point source.
        X=SX0
        Y=SY0
        Z=SZ0
      ENDIF
C
      X0=X
      Y0=Y
      Z0=Z
      PAGE=0.0D0
C
C  ****  First gamma ray. Loaded in common track.
C
      KPAR=2
      E=E1
C  ----  Initial direction.
      IF(LSCONE) THEN
        CALL GCONE(U,V,W)  ! Conical beam.
      ELSE
        W=CTHL+RAND(4.0D0)*DCTH  ! Rectangular beam.
        UV=SQRT(1.0D0-W*W)
        PHI=PHIL+RAND(5.0D0)*DPHI
        U=UV*COS(PHI)
        V=UV*SIN(PHI)
      ENDIF
C
      WGHT=1.0D0
      ILB(1)=1
      ILB(2)=0
      ILB(3)=0
      ILB(4)=0
      ILB(5)=0
      IPOL=0
C
C  ****  Second gamma ray (stored in the secondary stack).
C
      KPARS=2
      ES=E2
C  ----  Initial direction (same as in PENMAIN). .
      IF(LSCONE) THEN
        CALL GCONE(US,VS,WS)  ! Conical beam.
      ELSE
        WS=CTHL+RAND(4.0D0)*DCTH  ! Rectangular beam.
        UV=SQRT(1.0D0-WS*WS)
        PHI=PHIL+RAND(5.0D0)*DPHI
        US=UV*COS(PHI)
        VS=UV*SIN(PHI)
      ENDIF
C
      WGHTS=1.0D0
      ILBS(1)=1  ! Identifies primary particles.
      ILBS(2)=0
      ILBS(3)=0
      ILBS(4)=0
      ILBS(5)=0
      IPOLS=0
C
      IF(LAGE) THEN
        PAGE0=PAGE
        PAGE=PAGE-TAV1*LOG(RAND(2.0D0))
        CALL STORES(ES,X,Y,Z,US,VS,WS,WGHTS,KPARS,ILBS,IPOLS)
        PAGE=PAGE0
      ELSE
        CALL STORES(ES,X,Y,Z,US,VS,WS,WGHTS,KPARS,ILBS,IPOLS)
      ENDIF
C
      METAST=0
      RETURN
C
C  ************  De-excitation continues from a metastable level.
C
  200 CONTINUE
C  ----  Initial position (that of the decaying nucleus).
      X=X0
      Y=Y0
      Z=Z0
C  ....  add more particles, if needed.
C
      RETURN
      END
