C  ****  Files included to simplify compilation. Note that the packages
C  'penelope.f', 'pengeom.f', and 'penvared.f' contain modules; they
C  must be compiled first (in the indicated order), to allow other
C  program units to access these modules.
      INCLUDE 'penelope.f'
      INCLUDE 'pengeom.f'
      INCLUDE 'penvared.f'
      INCLUDE 'rita.f'
      INCLUDE 'timer.f'
      INCLUDE 'source.f'
	  
	  


C  Common blocks are defined in the included file 'pmcomms.f'.

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C         PPPPP   EEEEEE  N    N  M    M    AA    IIII  N    N         C
C         P    P  E       NN   N  MM  MM   A  A    II   NN   N         C
C         P    P  E       N N  N  M MM M  A    A   II   N N  N         C
C         PPPPP   EEEE    N  N N  M    M  AAAAAA   II   N  N N         C
C         P       E       N   NN  M    M  A    A   II   N   NN         C
C         P       EEEEEE  N    N  M    M  A    A  IIII  N    N         C
C                                                                      C
C                                                   (version 2014).    C
C                                                                      C
C  This program performs Monte Carlo simulation of electron-photon     C
C  showers in material structures described with the constructive      C
C  quadric geometry package 'PENGEOM.F'.                               C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C  PENELOPE/PENGEOM (version 2014)                                     C
C  Copyright (c) 2001-2014                                             C
C  Universitat de Barcelona                                            C
C                                                                      C
C  Permission to use, copy, modify, distribute and sell this software  C
C  and its documentation for any purpose is hereby granted without     C
C  fee, provided that the above copyright notice appears in all        C
C  copies and that both that copyright notice and this permission      C
C  notice appear in all supporting documentation. The Universitat de   C
C  Barcelona makes no representations about the suitability of this    C
C  software for any purpose. It is provided 'as is' without express    C
C  or implied warranty.                                                C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C  In the default mode, the program assumes that primary particles of a
C  given type, KPARP, are emitted from a point or an extended source,
C  either with fixed energy SE0 or with a specified (histogram-like)
C  energy spectrum. The initial direction of the primary particles is
C  sampled uniformly within a circle of the unit sphere (conical beam),
C  or within a 'rectangular' window on the unit sphere (rectangular
C  beam).
C
C  Alternatively, the program can read the initial state variables of
C  'primary' particles from a pre-calculated phase-space file. This
C  option may be used to split the simulation of complex problems into
C  several consecutive stages. The program also admits radiation sources
C  defined by the user, although this requires some editing work.
C
C  >>>>>>>> NOTE: All energies and lengths are given in eV and cm,
C                 respectively.
C
C
C  ************  Structure of the input data file.
C
C  Each line in the input data file consists of a 6-character keyword
C  (columns 1-6) followed either by numerical data (in free format) or
C  by a character string, which start at the 8th column. Keywords are
C  explicitly used/verified by the program (which is case sensitive!).
C  Notice also that the order of the data lines is important. The
C  keyword '______' (6 blanks) indicates comment lines, these can be
C  placed anywhere in the input file. The program ignores any text
C  following the first blank after the last numerical datum, or after
C  the character string, in each line (thus, in the table given below,
C  the comments in square brackets are ignored by the program). Lines
C  with some keywords (e.g., 'SPECTR', 'IPSFN') can appear an arbitrary
C  number of times, limited only by the allocated amount of memory.
C
C  The program assigns default values to many input variables; lines
C  that declare default values may be eliminated from the input file.
C
C
C  The structure of the input file is the following,
C
C  ....+....1....+....2....+....3....+....4....+....5....+....6....+....
C  TITLE  Title of the job, up to 65 characters.
C         . (the dot prevents editors from removing trailing blanks)
C         >>>>>>>> Source definition.
C  SKPAR  KPARP    [Primary particles: 1=electron, 2=photon, 3=positron]
C         KPARP=0 activates a user-defined SOURCE model.
C  SENERG SE0              [Initial energy (monoenergetic sources only)]
C  SPECTR Ei,Pi                 [E bin: lower-end and total probability]
C  SGPOL  SP1,SP2,SP3          [Stokes parameters for polarised photons]
C  SPOSIT SX0,SY0,SZ0                        [Coordinates of the source]
C  SBOX   SSX,SSY,SSZ                            [Source box dimensions]
C  SBODY  KB                [Active source body; one line for each body]
C  SCONE  THETA,PHI,ALPHA                  [Conical beam; angles in deg]
C  SRECTA THETAL,THETAU,PHIL,PHIU      [Rectangular beam; angles in deg]
C         .
C         >>>>>>>> Input phase-space file (psf).
C  IPSFN  psf-filename.ext         [Input psf name, up to 20 characters]
C  IPSPLI NSPLIT                                      [Splitting number]
C  WGTWIN WGMIN,WGMAX         [Weight window, RR & spl of psf particles]
C  EPMAX  EPMAX                 [Maximum energy of particles in the psf]
C         .
C         >>>>>>>> Material data and simulation parameters.
C                  Up to MAXMAT materials; 2 lines for each material.
C  MFNAME mat-filename.ext               [Material file, up to 20 chars]
C  MSIMPA EABS(1:3),C1,C2,WCC,WCR              [EABS(1:3),C1,C2,WCC,WCR]
C         .
C         >>>>>>>> Geometry and local simulation parameters.
C  GEOMFN geo-filename.ext               [Geometry file, up to 20 chars]
C  PARINP IP,PARINP(IP)                          [Replacement parameter]
C  DSMAX  KB,DSMAX(KB)              [KB, maximum step length in body KB]
C  EABSB  KB,EABSB(1:3,KB)   [KB, local absorption energies, EABSB(1:3)]
C         .
C         >>>>>>>> Interaction forcing.
C  IFORCE KB,KPAR,ICOL,FORCER,WLOW,WHIG  [KB,KPAR,ICOL,FORCER,WLOW,WHIG]
C         .
C         >>>>>>>> Bremsstrahlung splitting.
C  IBRSPL KB,IBRSPL                                [KB,splitting factor]
C         .
C         >>>>>>>> X-ray splitting.
C  IXRSPL KB,IXRSPL                                [KB,splitting factor]
C         .
C         >>>>>>>> Emerging particles. Energy and angular distributions.
C  NBE    EL,EU,NBE                      [Energy window and no. of bins]
C  NBANGL NBTH,NBPH           [No. of bins for the angles THETA and PHI]
C         .
C         >>>>>>>> Impact detectors (up to 25 different detectors).
C         IPSF=0; no psf is created.
C         IPSF=1; a psf is created (for only one detector).
C         IDCUT=0; tracking is discontinued at the detector entrance.
C         IDCUT=1; the detector does not affect the tracking.
C         IDCUT=2; the detector does not affect tracking, the energy
C                  distribution of particle fluence (integrated over the
C                  volume of the detector) is calculated.
C  IMPDET EL,EU,NBE,IPSF,IDCUT      [E-window, no. of bins, IPSF, IDCUT]
C  IDSPC  spc-impdet-##.dat               [Spectrum file name, 20 chars]
C  IDPSF  psf-impdet-##.dat            [Phase-space file name, 20 chars]
C  IDFLNC fln-impdet-##.dat       [Fluence spectrum file name, 20 chars]
C  IDAGEL AGEL,AGEU,NAGE                  [Age interval and no. of bins]
C  IDAGEF age-impdet-##.dat       [Age-distribution file name, 20 chars]
C  IDBODY KB                                               [Active body]
C  IDKPAR KPAR                              [Type of detected particles]
C         .
C         >>>>>>>> Energy-deposition detectors (up to 25).
C  ENDETC EL,EU,NBE                      [Energy window and no. of bins]
C  EDSPC  spc-enddet-##.dat        [Output spectrum file name, 20 chars]
C  EDBODY KB                                               [Active body]
C         .
C         >>>>>>>> Absorbed dose distribution.
C  GRIDX  XL,XU,NDBX         [X coords of the box vertices, no. of bins]
C  GRIDY  YL,YU,NDBY         [Y coords of the box vertices, no. of bins]
C  GRIDZ  ZL,ZU,NDBZ         [Z coords of the box vertices, no. of bins]
C  GRIDR  RU,NDBR               [Radius of the dose volume, no. of bins]
C         .
C         >>>>>>>> Job properties.
C  RESUME dump1.dmp               [Resume from this dump file, 20 chars]
C  DUMPTO dump2.dmp                  [Generate this dump file, 20 chars]
C  DUMPP  DUMPP                                 [Dumping period, in sec]
C         .
C  RSEED  ISEED1,ISEED2           [Seeds of the random-number generator]
C  NSIMSH DSHN                     [Desired number of simulated showers]
C  TIME   TIMEA                       [Allotted simulation time, in sec]
C  END                                  [Ends the reading of input data]
C  ....+....1....+....2....+....3....+....4....+....5....+....6....+....
C
C
C  The following listing describes the function of each of the keywords,
C  the accompanying data and their default values. For clarity, blanks
C  in keywords are indicated as '_'.
C
C  TITLE_ : Title of the job (up to 65 characters).
C             DEFAULT: none (the input file must start with this line)
C
C           The TITLE string is used to mark dump files. To prevent the
C           improper use of wrong resuming files, change the title each
C           time you modify basic parameters of your problem. The code
C           will then be able to identify the inconsistency and to print
C           an error message before stopping.
C
C  >>>>>>>> Source definition.
C
C  SKPAR_ : Type of primary particle KPARP (1=electrons, 2=photons or
C           3=positrons).
C             DEFAULT: KPARP=1
C           If KPARP=0, the initial states of primary particles are
C           set by subroutine SOURCE, to be provided by the user. An
C           example of that subroutine, corresponding to a 60-Co source
C           (two gamma rays in each nuclear deexcitation), is included
C           in the PENMAIN package (file 'source.f').
C
C  SENERG : For a monoenergetic source, initial energy SE0 of primary
C           particles.
C             DEFAULT: SE0=1.0E6
C
C  SPECTR : For a source with continuous (stepwise constant) spectrum,
C           each 'SPECTR' line gives the lower end-point of an energy
C           bin of the source spectrum (Ei) and the associated relative
C           probability (Pi), integrated over the bin. Up to NSEM=1000
C           lines, in arbitrary order. The upper end of the spectrum is
C           defined by entering a line with Ei equal to the upper energy
C           end point and with a negative Pi value.
C             DEFAULT: none
C
C  SGPOL_ : This line activates the simulation of polarisation effects
C           in the scattering of photons (electrons and positrons are
C           assumed to be unpolarised). SP1, SP2, SP3 are the Stokes
C           parameters of primary photons, which define the degrees of
C           linear polarisation at 45 deg azimuth, of circular
C           polarisation, and of linear polarisation at zero azimuth,
C           respectively. It is assumed that secondary photons are
C           emitted with null polarisation (SP1=SP2=SP3=0).
C             DEFAULT: none
C
C  SPOSIT : Coordinates of the source centre.
C             DEFAULT: SX0=SY0=0.0, SZ0=0.0D0
C
C  SBOX__ : Extended source box. The source has uniform activity within
C           the volume of a right prism centred at the point (SX0,SY0,
C           SZ0) and whose sides have lengths SSX, SSY and SSZ.
C             DEFAULT: SSX=SSY=SSZ=0.0
C
C  In the case of a extended source, the active volume can be restricted
C  to that of a body or a set of bodies, which must be defined as parts
C  of the geometry. The activity of the source is assumed to be uniform
C  within the volume of the intersection of the active bodies and the
C  source box. Note that the initial coordinates of primary particles
C  are sampled by the rejection method; the sampling efficiency is equal
C  to the fraction of the source box volume that is occupied by active
C  bodies.
C
C  To define each active source body, add the following line:
C
C  SBODY_ : Active source body (PENGEOM sequential body label). One
C           line for each body.
C             DEFAULT: None
C           The program stops if the source box has not been defined
C           previously.
C
C  The initial direction of primary particles is sampled uniformly
C  within a circle on the unit sphere (conical beam), or within a
C  'rectangular' window on the unit sphere (rectangular beam).
C
C  SCONE_ : Conical source beam. Polar and azimuthal angles of the
C           beam axis direction, THETA and PHI, and angular aperture,
C           ALPHA, in deg.
C             DEFAULTS: THETA=0.0, PHI=0.0, ALPHA=0.0
C
C           The case ALPHA=0 defines a monodirectional source, and ALPHA
C           =180 deg corresponds to an isotropic source.
C
C  SRECTA : Rectangular source beam. Limiting polar and azimuthal angles
C           of the source beam window, (THETAL,THETAU)x(PHIL,PHIU), in
C           deg.
C             DEFAULTS: THETAL=0.0, THETAU=0.0, PHIL=0.0, PHIU=0.0
C
C           The case THETAL=THETAU, PHIL=PHIU defines a monodirectional
C           source. To define an isotropic source, set THETAL=0, THETAU=
C           180, PHIL=0 and PHIU=360.
C
C           --> Notice that the default source is a pencil beam that
C           moves upwards along the Z-axis.
C
C  >>>>>>>> Input phase-space file.
C
C  The initial state variables of primary particles can be read directly
C  from a set of pre-calculated phase-space files (psf). When this
C  option is active, previous definitions about the source are ignored.
C  Photons from the psf's are assumed to be unpolarised.
C
C  IPSFN_ : Name of an input psf (up to 20 characters).
C             DEFAULT: none
C           Up to 100 psf's may be declared. They are read sequentially.
C
C  The input psf is in ASCII format. Each line defines the initial state
C  of a particle; it contains the following quantities in free format
C  (and in the order they are listed here):
C    -- KPAR, type of particle (1, electron; 2, photon; 3, positron).
C    -- E, energy.
C    -- X,Y,Z, position coordinates.
C    -- U,V,W, direction cosines.
C    -- WGHT, weight.
C    -- ILB(1),ILB(2),ILB(3),ILB(4), a set of indices that provide
C           information on how the particle was generated (see the file
C           'manual.txt').
C    -- NSHI, incremental shower number (difference between the shower
C           numbers of the present particle and the one preceding it
C           in the psf).
C  Phase-space files can be generated by running PENMAIN using an impact
C  detector with the flag IPSF=1 (see below).
C
C  Because of the limited size of the psf's, the results of analogue
C  simulations tend to be 'too noisy'. This can be partially corrected
C  by splitting the particles from the psf.
C
C  IPSPLI : Splitting number. Each particle in the psf's will be split
C           into NSPLIT equivalent particles, with weights equal to
C           WGHT/NSPLIT.
C             DEFAULT: NSPLIT=1 (no splitting)
C
C  --> WARNING: Notice that there is a 'latent' uncertainty in the psf,
C  which sets a limit to the accuracy that can be attained by using
C  large splitting numbers.
C
C  WGTWIN : Weight window, (WGMIN,WGMAX). Particles in the phase-space
C           file that have initial weights WGHT less than WGMIN will be
C           subjected to Russian roulette, and those with WGHT larger
C           than WGMAX will be split. Note that the weight window has
C           preference over the splitting option, i.e., a particle will
C           be split into NSPLIT or less particles only if the latter
C           have weights larger than WGMIN.
C             DEFAULTS: WGMIN=1.0E-35, WGMAX=1.0E35  (no action)
C
C  EPMAX_ : Maximum energy (in eV) of particles in the psf's.
C           EPMAX is the upper limit of the energy interval covered by
C           the simulation lookup tables. To minimize interpolation
C           errors, EPMAX should not be much larger than the maximum
C           energy actually occurring during the simulation.
C
C           When the initial state variables of particles are read from
C           a psf, this parameter is required to initialise PENELOPE and
C           is critical; the code crashes if it finds a particle that
C           has energy larger than EPMAX.
C             DEFAULT: EPMAX=1.0E9 (interpolation is not optimal)
C
C  >>>>>>>> Material data and simulation parameters.
C
C  Each material is defined by introducing the following _two_ lines;
C
C  MFNAME : Name of a PENELOPE input material data file (up to 20
C           characters). This file must be generated in advance by
C           running the program MATERIAL.
C             DEFAULT: none
C
C  MSIMPA : Set of simulation parameters for this material; absorption
C           energies, EABS(1:3,M), elastic scattering parameters, C1(M)
C           and C2(M), and cutoff energy losses for inelastic collisions
C           and bremsstrahlung emission, WCC(M) and WCR(M).
C             DEFAULTS: EABS(1,M)=EABS(3,M)=0.01*EPMAX,
C                       EABS(2,M)=0.001*EPMAX
C                       C1(M)=C2(M)=0.1, WCC=EABS(1,M), WCR=EABS(2,M)
C           EPMAX is the upper limit of the energy interval covered by
C           the simulation lookup tables.
C
C  Note that we must declare a separate material data file name and a
C  set of simulation parameters for each material. The label (material
C  number) assigned by PENELOPE to each material is determined by the
C  ordering of the material data files in the PENMAIN input file. That
C  is, the first, second, ... materials are assigned the labels 1, 2,
C  ... These labels are also used in the geometry definition file.
C
C  The original programs in the distribution package allow up to 10
C  materials. This number can be increased by changing the value of the
C  parameter MAXMAT in the original source files.
C
C  >>>>>>>> Geometry and local simulation parameters.
C
C  GEOMFN : PENGEOM geometry definition file name (a string of up to
C           20 characters).
C             DEFAULT: none.
C
C           --> The geometry definition file can be debugged/visualised
C           with the viewers GVIEW2D and GVIEW3D (operable only under
C           Windows).
C
C           The bodies in the material structure are normally identified
C           by the sequential labels assigned by PENGEOM. For complex
C           geometries, however, it may be more practical to employ user
C           labels, i.e., the four-character strings that identify the
C           body in the geometry definition file. In PENMAIN (and only
C           in the parts of the code that follow the definition of the
C           geometry), a body can be specified by giving either its
C           PENGEOM numerical label or its user label enclosed in a
C           pair of apostrophes (e.g., 'BOD1'). However, bodies that
C           result from the cloning of modules (as well as those defined
C           in an INCLUDEd geometry file) do not have a user label and
C           only the PENGEOM numerical label is acceptable.
C
C  PARINP : The values of certain parameters of the geometry definition
C           may be defined from the main program by means of the array
C           PARINP (an input argument of the GEOMIN subroutine). The
C           entered PARINP(IP) value replaces the parameter values that
C           are marked with the index IP in the geometry definition
C           file.
C             DEFAULT:  none
C
C  DSMAX_ : Maximum step length DSMAX(KB) of electrons and positrons in
C           body KB. This parameter is important only for thin bodies;
C           it should be given a value of the order of one tenth of the
C           body thickness or less.
C             DEFAULT: DSMAX=1.0E20 (no step length control)
C
C  EABSB_ : Local absorption energies EABSB(KPAR,KB) of particles of
C           type KPAR in body KB. These values must be larger than
C           EABS(KPAR,M), where M is the material of body KB. When the
C           particle is moving within body KB, the absorption energy
C           EABS(KPAR,M) is temporarily set equal to EABSB(KPAR,KB).
C           Thus, the simulation of the particle history is discontinued
C           when the energy becomes less than EABSB(KPAR,KB). This
C           feature can be used, e.g., to reduce the simulation work in
C           regions of lesser interest.
C             DEFAULTS: EABSB(KPAR,KB)=EABS(KPAR,M)  (no action)
C
C  >>>>>>>> Interaction forcing.
C
C  IFORCE : Activates forcing of interactions of type ICOL of particles
C           KPAR in body KB. FORCER is the forcing factor, which must
C           be larger than unity. WLOW and WHIG are the lower and upper
C           limits of the weight window where interaction forcing is
C           applied. When several interaction mechanisms are forced in
C           the same body, the effective weight window is set equal to
C           the intersection of the windows for these mechanisms.
C             DEFAULT: no interaction forcing
C
C           If the mean free path for real interactions of type ICOL is
C           MFP, the program will simulate interactions of this type
C           (real or forced) with an effective mean free path equal to
C           MFP/FORCER.
C
C           TRICK: a negative input value of FORCER, -FN, is assumed to
C           mean that a particle with energy E=EPMAX should interact,
C           on average, +FN times in the course of its slowing down to
C           rest, for electrons and positrons, or along a mean free
C           path, for photons. This is very useful, e.g., to generate
C           x-ray spectra from bulk samples.
C
C  The real effect of interaction forcing on the efficiency is not easy
C  to predict. Please, do tentative runs with different FORCER values
C  and check the efficiency gain (or loss!).
C
C  >>>>>>>> Bremsstrahlung splitting.
C
C  IBRSPL : Activates bremsstrahlung splitting in body KB for electrons
C           and positrons with weights in the window (WLOW,WHIG) where
C           interaction forcing is applied. The integer IBRSPL is the
C           splitting factor.
C             DEFAULT: no bremsstrahlung splitting
C
C           Note that bremsstrahlung splitting is applied in combination
C           with interaction forcing and, consequently, it is activated
C           only in those bodies where interaction forcing is active.
C
C  >>>>>>>> X-ray splitting.
C
C  IXRSPL : Splitting of characteristic x rays emitted in body KB, from
C           any element. Each unsplit x ray with ILB(2)=2 (i.e., of the
C           second generation) when extracted from the secondary stack
C           is split into IXRSPL quanta. The new, lighter, quanta are
C           assigned random directions distributed isotropically.
C             DEFAULT: no x-ray splitting
C
C  >>>>>>>> Energy and angular distributions of emerging particles.
C
C  NBE___ : Limits EL and EU of the interval where energy distributions
C           of emerging particles are tallied. Number of energy bins
C           (.LE. 1000).
C             DEFAULT: EL=0.0, EU=EPMAX, NBE=500
C
C           NBE is the number of bins in the output energy distribution
C           (.LE. 1000). If NBE is positive, energy bins have uniform
C           width, DE=(EU-EL)/NBE. When NBE is negative, the bin width
C           increases geometrically with the energy, i.e., the energy
C           bins have uniform width on a logarithmic scale.
C
C  NBANGL : Numbers of bins for the polar angle THETA and the azimuthal
C           angle PHI, respectively, NBTH and NBPH (.LE. 3600 and 180,
C           respectively).
C             DEFAULT: NBTH=90, NBPH=1 (azimuthal average)
C
C           If NBTH is positive, angular bins have uniform width,
C           DTH=180./NBTHE. When NBTH is negative, the bin width
C           increases geometrically with THETA, i.e., the bins have
C           uniform width on a logarithmic scale.
C
C           NOTE: In the output files, the terms 'upbound' and
C           'downbound' are used to denote particles that leave the
C           material system moving upwards (W>0) and downwards (W<0),
C           respectively.
C
C  >>>>>>>> Impact detectors.
C
C  Each impact detector consists of a set of active bodies, which must
C  have been defined as parts of the geometry. The output spectrum is
C  the energy distribution of particles that entered any of the active
C  bodies coming from a body that is not active (i.e. that is not part
C  of the detector). Notice that a detected particle can re-enter the
C  detector volume and, consequently, be 'counted' several times (except
C  when the flag IDCUT is set equal to 0, see below).
C
C  Active bodies cannot be void, because the geometry routines would not
C  stop particles at their limiting surfaces. In case you need to define
C  detectors outside the material system, fill them with an arbitrary
C  material of very small density to avoid perturbing the transport
C  process.
C
C  To define each impact detector, insert the following block of lines;
C
C  IMPDET : Starts the definition of a new detector. Up to 25 different
C           detectors can be considered.
C           EL and EU are the lower and upper limits of the energy
C             window covered by the impact detector.
C           NBE is the number of bins in the output energy spectrum of
C             the detector (.LE. 1000). If NBE is positive, energy bins
C             have uniform width, DE=(EU-EL)/NBE. When NBE is negative,
C             the bin width increases geometrically with the energy,
C             i.e., the energy bins have uniform width on a logarithmic
C             scale.
C
C           The integer flag IPSF serves to activate the creation of a
C           phase-space file (psf), which contains the state variables
C           of all particles that enter the detector. Use this option
C           with care, because psf's may grow very fast.
C           IPSF=0; no psf is created.
C           IPSF=1; the psf is created. Only one PSF can be created in
C             each simulation run.
C
C           The integer flag IDCUT allows discontinuing the tracking of
C           particles that enter the detector.
C           IDCUT=0; the simulation of a particle is discontinued when
C             it enters the detector (useful to stop the simulation of
C             particles recorded in a psf).
C           IDCUT=1; the presence of the detector does not affect the
C             tracking of particles.
C           IDCUT=2; the presence of the detector does not affect the
C             tracking of particles. The distribution of particle
C             fluence with respect to energy (integrated over the volume
C             of the detector) is tallied. The calculated distribution
C             has dimensions of length/energy.
C
C             DEFAULTS: None
C
C  IDPSF_ : Name of the output phase-space file (up to 20 characters).
C             DEFAULT: 'psf-impdet-##.dat'
C
C  IDSPC_ : Name of the output energy spectrum file (up to 20
C           characters).
C             DEFAULT: 'spc-impdet-##.dat'
C
C  IDFLNC : Name of the output file with the energy distribution of
C           particle fluence (20 characters). This file is generated
C           only when IDCUT=2.
C             DEFAULT: 'fln-impdet-##.dat'
C
C  IDAGEL : Activates the evaluation of the age of particles, defined as
C           the time elapsed since the start of the primary particle
C           that originated the shower. The program generates the age
C           distribution of detected particles, i.e., particles of the
C           types declared in lines IDKPAR (see below) that enter the
C           detector with energy in the window (EL,EU). The distribution
C           is tallied for ages in the interval between AGEL and AGEU
C           (both in seconds), which is partitioned into NAGE bins. If
C           NAGE is positive, the age bins have uniform width. When
C           NAGE is negative, the width of age bins is uniform on a
C           logarithmic scale.
C
C             DEFAULTS: NAGE=100, AGEL=0.0, AGEU must always be
C                       specified
C
C  IDAGEF : Name of the output age distribution file (up to 20
C           characters)
C             DEFAULT: 'age-impdet-##.dat'
C
C  IDBODY : Active body of the detector. One line for each active body.
C             DEFAULT: None
C           --> Notice that a body cannot be part of more than one
C           impact detector.
C
C  IDKPAR : Type of particle that is detected (1=electrons, 2=photons or
C           3=positrons). One line for each type.
C
C           The detector has no effect for particles that are not
C           detected. This feature can be used, e.g., to make a body or
C           a set of bodies opaque to particles of a certain type.
C
C             DEFAULT: All particles are detected
C
C  >>>>>>>> Energy-deposition detectors.
C
C  Each energy-deposition detector consists of a set of active bodies,
C  which must have been defined as parts of the geometry. The output
C  spectrum is the distribution of absorbed energy (per primary shower)
C  in the active bodies.
C
C           *** WARNING: The energy-deposition spectrum may be strongly
C           biased when interaction forcing is applied, even outside the
C           detector bodies.
C
C  To define each energy-deposition detector insert the following block
C  of lines;
C
C  ENDETC : Starts the definition of a new energy-deposition detector.
C           Up to 25 different detectors can be considered.
C           EL and EU are the lower and upper limits of the energy
C             window covered by the detector.
C           NBE is the number of bins in the output energy spectrum of
C             the detector (.LE. 1000). If NBE is positive, energy bins
C             have uniform width, DE=(EU-EL)/NBE. When NBE is negative,
C             the bin width increases geometrically with the energy,
C             i.e., the energy bins have uniform width on a logarithmic
C             scale.
C
C             DEFAULTS: None
C
C  EDSPC_ : Name of the output spectrum file (up to 20 characters).
C             DEFAULT: 'spc-enddet-##.dat'
C
C  EDBODY : Active body KB of the detector. One line for each active
C           body.
C             DEFAULT: None
C           --> Notice that a body cannot be part of more than one
C           energy-deposition detector.
C
C  >>>>>>>> Dose map.
C
C  Generally, the program can calculate the dose distribution inside a
C  parallelepiped (dose box) whose edges are parallel to the axes of the
C  laboratory frame. The dose box is defined by giving the coordinates
C  of its vertices. The dose is tallied using a uniform orthogonal grid
C  with NDBX, NDBY and NDBZ bins (= voxels) along the directions of
C  the respective coordinate axes. These numbers should be odd, to make
C  sure that each 'central' axis (i.e., the line that join the centres
C  of two opposite faces of the box) goes through the centres of a row
C  of voxels.
C
C  GRIDX_ : X-coordinates of the vertices of the dose box and number of
C           bins in the X direction.
C             DEFAULT: None
C  GRIDY_ : Y-coordinates of the vertices of the dose box and number of
C           bins in the Y direction.
C             DEFAULT: None
C  GRIDZ_ : Z-coordinates of the vertices of the dose box and number of
C           bins in the Z direction.
C             DEFAULT: None
C
C  The efficiency of the dose map calculation can be increased by taking
C  advantage of possible symmetries of the system (source and geometry).
C  In problems with axial symmetry about the Z axis, it is advantageous
C  to tally the dose distribution in the volume of a cylinder of radius
C  RU, about the Z axis, limited by the planes Z=ZL and Z=ZU. For
C  systems with spherical symmetry about the origin of coordinates, it
C  is most convenient to consider the radial dose distribution in a
C  sphere of radius RU. The generation of these symmetric dose maps is
C  activated by entering the line
C
C  GRIDR_ : Radius RU of the dose zone and number of radial bins.
C             DEFAULT: None
C
C  Specifically, when the input file contains only the lines GRIDZ_ and
C  GRIDR_, the program assumes that the dose distribution is axially
C  symmetric and generates a cylindrical map. When the input file has
C  only the line GRIDR_, spherical symmetry is assumed and the radial
C  distribution of absorbed dose is tallied.
C
C  The different types of dose maps are mutually exclusive. Notice that
C  when the assumed symmetry does not hold, the program may not be able
C  to evaluate the masses of voxels correctly.
C
C  --> The meshes defined here to calculate dose distributions can be
C  used to tally other 3D distributions (e.g. the space distribution of
C  inner-shell ionization events, used in electron-probe microanalysis).
C  This, however, requires to edit the present source file.
C
C  >>>>>>>> Job properties.
C
C  RESUME : The program will read the dump file named `dump1.dmp' (up to
C           20 characters) and resume the simulation from the point
C           where it was left. Use this option very, _VERY_ carefully.
C           Make sure that the input data file is fully consistent with
C           the one used to generate the dump file.
C             DEFAULT: off
C
C  DUMPTO : Generate a dump file named 'dump2.dmp' (up to 20 characters)
C           after completing the simulation run. This allows the
C           simulation to be resumed later on to improve statistics.
C             DEFAULT: off
C
C           NOTE: If the file 'dump2.dmp' already exists, it is
C           overwritten.
C
C  DUMPP_ : When the DUMPTO option is activated, simulation results are
C           written in the output files every DUMPP seconds. This option
C           is useful to check the progress of long simulations. It also
C           allows the program to be run with a long execution time and
C           to be stopped when the required statistical uncertainty has
C           been reached.
C             DEFAULT: DUMPP=1.0E15
C
C  RSEED_ : Seeds of the random-number generator. When ISEED1 is equal
C           to a negative integer, -N, the seeds are set by calling
C           subroutine RAND0(N) with the input argument equal to N. This
C           ensures that sequences of random numbers used in different
C           runs of the program (with different values of N) are truly
C           independent.
C             DEFAULT: ISEED1=1, ISEED2=1
C
C  NSIMSH : Desired number of simulated showers.
C             DEFAULT: DSHN=2.0E9
C
C  TIME__ : Allotted simulation time, in sec.
C             DEFAULT: TIMEA=2.0E9
C
C  END___ : Ends the reading of the input file. This line is needed only
C           when the TIME__ line is missing.
C
C  The program is aborted when an incorrect input datum is found. The
C  conflicting quantity usually appears in the last line of the output
C  file 'penmain.dat'. If the trouble is with arrays having dimensions
C  smaller than required, the program indicates how the problem can be
C  solved (this usually requires editing the source file, be careful).
C
C  The clock subroutine (TIMER) may have to be adapted to your specific
C  computer-compiler configuration; standard Fortran 77 does not provide
C  timing tools. However, the routines in module TIMER.F do work for
C  many Fortran compilers.
C
C  ************  Generating the executable PENMAIN and running it.
C
C  To generate the executable binary file PENMAIN.EXE, compile and link
C  the Fortran source files PENMAIN.F, PENELOPE.F, PENGEOM.F, PENVARED.F
C  and TIMER.F. For example, if you are using the GNU Fortran compiler
C  under Windows, place these five files in the same directory, open a
C  command window and from that directory enter the command
C    'gfortran -Wall -Os PENMAIN.F -o PENMAIN.EXE'
C  (The same, with file names in lowercase, should work under Linux).
C
C  To run PENMAIN, you have to generate an input data file, let's call
C  it PENMAIN.IN, and the corresponding geometry definition and material
C  data files. Place these files and the binary file PEMAIN.EXE in the
C  same directory and, from there, issue the command
C    'PENMAIN.EXE < PENMAIN.IN'
C
C  The calculated distributions are written in separate files, whose
C  names have the extension '.dat'. These files are in a format suited
C  for direct visualisation with GNUPLOT (version 4.2).
C
C  *********************************************************************
C                       MAIN PROGRAM
C  *********************************************************************
C
      USE PENELOPE_mod
      USE TRACK_mod
      USE PENGEOM_mod
      USE PENVARED_mod
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
	  PARAMETER (NRP=8000)
	  CHARACTER*2 LASYMB
	  PARAMETER (NTP=12000)
	  PARAMETER (NRX=60000)
	  PARAMETER (NM=512)
	  PARAMETER (NR=128)
	  PARAMETER (NBE=57, NBW=32)
	  PARAMETER (NQ=250,NEX=1024)
	  COMMON/QSURF/AXX(NS),AXY(NS),AXZ(NS),AYY(NS),AYZ(NS),AZZ(NS),
     1    AX(NS),AY(NS),AZ(NS),A0(NS),NSURF,KPLANE(NS)
      COMMON/QBODY/KBODY(NB,NXG),KBOMO(NB)
      COMMON/QTREE/NBODYS,KMOTH(NB),KDGHT(NB,NXG),KSURF(NB,NXG),
     1  KFLAG(NB,NXG),KSP(NS),NWARN
	 
	  COMMON/CECUTR/ECUTR(MAXMAT)
	  COMMON/CSGAWR/ISGAW  ! Controls warning messages from SUMGA.
      COMMON/CERSEC/IERSEC
	  
	  COMMON/CEGRID/EMIN,EL,EU,ET(NEGP),DLEMP(NEGP),DLEMP1,DLFC,
     1  XEL,XE,XEK,KE
	  COMMON/CESI0/XESI(NRP,16),IESIF(99),IESIL(99),NSESI(99),NCURE
	  COMMON/CPSI0/XPSI(NRP,16),IPSIF(99),IPSIL(99),NSPSI(99),NCURP
	
      COMMON/CADATA/ATW(99),EPX(99),RSCR(99),ETA(99),EB(99,30),
     1  ALW(99,30),CP0(99,30),IFI(99,30),IKS(99,30),NSHT(99),LASYMB(99)
	 
      COMMON/CGPH00/EPH(NTP),XPH(NTP,17),IPHF(99),IPHL(99),NPHS(99),NCUR
c     COMMON/CRELAX/P(NRX),ET(NRX),F(NRX),IS0(NRX),IS1(NRX),IS2(NRX),
c     1              IFIRST(99,16),ILAST(99,16),NCUR,KS,MODER

	  COMMON/CRITAA/XA(NM),AA(NM),BA(NM),FA(NM),IA(NM),NPM1A
      COMMON/CRITA/XT(NM),PAC(NM),DPAC(NM),A(NM),B(NM),IL(NM),IU(NM),
     1  NPM1
	  COMMON/CRITAN/CNORM
	  
	  
	  
C  ****  Composition data.
      COMMON/COMPOS/STF(MAXMAT,30),ZT(MAXMAT),AT(MAXMAT),RHO(MAXMAT),
     1  VMOL(MAXMAT),IZ(MAXMAT,30),NELEM(MAXMAT)

      COMMON/CRANGE/RANGE(3,MAXMAT,NEGP),RANGEL(3,MAXMAT,NEGP)
C  ****  E/P inelastic collisions.
      PARAMETER (NO=512)
      COMMON/CEIN/EXPOT(MAXMAT),OP2(MAXMAT),F(MAXMAT,NO),UI(MAXMAT,NO),
     1  WRI(MAXMAT,NO),KZ(MAXMAT,NO),KS(MAXMAT,NO),NOSC(MAXMAT)
      COMMON/CEINTF/T1EI(NEGP),T2EI(NEGP),T1PI(NEGP),T2PI(NEGP)
C  ****  Partial cross sections of individual shells/oscillators.
      COMMON/CEIN00/SEH0(NO),SEH1(NO),SEH2(NO),SES0(NO),SES1(NO),
     1              SES2(NO),SET0(NO),SET1(NO),SET2(NO)
      COMMON/CPIN00/SPH0(NO),SPH1(NO),SPH2(NO),SPS0(NO),SPS1(NO),
     1              SPS2(NO),SPT0(NO),SPT1(NO),SPT2(NO)
C  ****  Inner-shell ionisation by electron and positron impact.
C  ****  
      COMMON/CEINAC/EINAC(MAXMAT,NEGP,NO),IEIN(MAXMAT,NO),NEIN(MAXMAT)
      COMMON/CESIAC/ESIAC(MAXMAT,NEGP,NO),IESI(MAXMAT,NO),NESI(MAXMAT)
      COMMON/CESIN/XSEIN(NEGP,NO),XSESI(NEGP,NO),ISIE(NO)
C  ****  Positron inelastic coll. and inner-shell ionisation tables.
      COMMON/CPINAC/PINAC(MAXMAT,NEGP,NO),IPIN(MAXMAT,NO),NPIN(MAXMAT)
      COMMON/CPSIAC/PSIAC(MAXMAT,NEGP,NO),IPSI(MAXMAT,NO),NPSI(MAXMAT)
      COMMON/CPSIN/XSPIN(NEGP,NO),XSPSI(NEGP,NO),ISIP(NO)
	  
	  PARAMETER (NOCO=512)
      COMMON/CGCO/FCO(MAXMAT,NOCO),UICO(MAXMAT,NOCO),FJ0(MAXMAT,NOCO),
     1  PTRSH(MAXMAT,NOCO),KZCO(MAXMAT,NOCO),KSCO(MAXMAT,NOCO),
     2  NOSCCO(MAXMAT)
C  ****  Electron simulation tables.
      COMMON/CEIMFP/SEHEL(MAXMAT,NEGP),SEHIN(MAXMAT,NEGP),
     1  SEISI(MAXMAT,NEGP),SEHBR(MAXMAT,NEGP),SEAUX(MAXMAT,NEGP),
     2  SETOT(MAXMAT,NEGP),CSTPE(MAXMAT,NEGP),RSTPE(MAXMAT,NEGP),
     3  DEL(MAXMAT,NEGP),W1E(MAXMAT,NEGP),W2E(MAXMAT,NEGP),
     4  DW1EL(MAXMAT,NEGP),DW2EL(MAXMAT,NEGP),
     5  RNDCE(MAXMAT,NEGP),AE(MAXMAT,NEGP),BE(MAXMAT,NEGP),
     6  T1E(MAXMAT,NEGP),T2E(MAXMAT,NEGP)
      COMMON/CLAS1E/TSTPE(MAXMAT,NEGP),TSTRE(MAXMAT,NEGP),
     1  TRL1E(MAXMAT,NEGP),TRL2E(MAXMAT,NEGP)
C  ****  Positron simulation tables.
      COMMON/CPIMFP/SPHEL(MAXMAT,NEGP),SPHIN(MAXMAT,NEGP),
     1  SPISI(MAXMAT,NEGP),SPHBR(MAXMAT,NEGP),SPAN(MAXMAT,NEGP),
     2  SPAUX(MAXMAT,NEGP),SPTOT(MAXMAT,NEGP),CSTPP(MAXMAT,NEGP),
     3  RSTPP(MAXMAT,NEGP),W1P(MAXMAT,NEGP),W2P(MAXMAT,NEGP),
     4  DW1PL(MAXMAT,NEGP),DW2PL(MAXMAT,NEGP),
     5  RNDCP(MAXMAT,NEGP),AP(MAXMAT,NEGP),BP(MAXMAT,NEGP),
     6  T1P(MAXMAT,NEGP),T2P(MAXMAT,NEGP)
      COMMON/CLAS1P/TSTPP(MAXMAT,NEGP),TSTRP(MAXMAT,NEGP),
     1  TRL1P(MAXMAT,NEGP),TRL2P(MAXMAT,NEGP)
C  ****  Elastic scattering of electrons and positrons.
	  COMMON/CEEL00/EJT(NEGP),XE0(NEGP),XE1(NEGP),XE2(NEGP),XP0(NEGP),
     1  XP1(NEGP),XP2(NEGP),T1E0(NEGP),T2E0(NEGP),T1P0(NEGP),T2P0(NEGP),
     2  EJTL(NEGP),FJL(NEGP),A2(NEGP),B2(NEGP),C(NEGP),D(NEGP)



C  ****  Electron and positron radiative yields.
      COMMON/CBRYLD/EBRY(MAXMAT,NEGP),PBRY(MAXMAT,NEGP)
C  ****  Photon simulation tables.
      COMMON/CGIMFP/SGRA(MAXMAT,NEGP),SGCO(MAXMAT,NEGP),
     1  SGPH(MAXMAT,NEGP),SGPP(MAXMAT,NEGP),SGAUX(MAXMAT,NEGP)
      PARAMETER (NDIM=12000)
      COMMON/CGPH01/ER(NDIM),XSR(NDIM),NPHD
      COMMON/CGPP01/TRIP(MAXMAT,NEGP)
	  
	  
	  
	  COMMON/CEBR/WB(NBW),PBCUT(MAXMAT,NEGP),WBCUT(MAXMAT,NEGP),
     1  PDFB(MAXMAT,NEGP,NBW),DPDFB(MAXMAT,NEGP,NBW),
     2  PACB(MAXMAT,NEGP,NBW),ZBR2(MAXMAT)
	 
	  COMMON/CEBR01/EBT(NBE),XS(NBE,NBW),TXS(NBE),X1(NBE),Y1(NBE)
      COMMON/CEBR02/P0(MAXMAT,NEGP,NBW)
	  
	  COMMON/CBRANG/BET(6),BK(21),BP1(MAXMAT,6,21,4),
     1              BP2(MAXMAT,6,21,4),ZBEQ(MAXMAT)
	 
	 COMMON/CEIN01/EI,EE,CPS,AMOL,MOM
	 
	 COMMON/CSUMGA/IERGA,NCALL
	 
      PARAMETER (NE=96,NA=606)
      COMMON/CDCSEP/ETS(NE),ETL(NE),TH(NA),THR(NA),XMU(NA),XMUL(NA),
     1       ECS(NE),ETCS1(NE),ETCS2(NE),EDCS(NE,NA),
     2       PCS(NE),PTCS1(NE),PTCS2(NE),PDCS(NE,NA),
     3       DCSI(NA),DCSIL(NA),CSI,TCS1I,TCS2I
	 
      PARAMETER (NP=128)
      COMMON/CEELDB/XSE(NP,NEGP,MAXMAT),PSE(NP,NEGP,MAXMAT),
     1              ASE(NP,NEGP,MAXMAT),BSE(NP,NEGP,MAXMAT),
     2              ITLE(NP,NEGP,MAXMAT),ITUE(NP,NEGP,MAXMAT)
      COMMON/CPELDB/XSP(NP,NEGP,MAXMAT),PSP(NP,NEGP,MAXMAT),
     1              ASP(NP,NEGP,MAXMAT),BSP(NP,NEGP,MAXMAT),
     2              ITLP(NP,NEGP,MAXMAT),ITUP(NP,NEGP,MAXMAT)
      COMMON/CELSEP/EELMAX(MAXMAT),PELMAX(MAXMAT),
     1              RNDCEd(MAXMAT,NEGP),RNDCPd(MAXMAT,NEGP)
	 
	 
	  COMMON/CGRA00/FACTE,Q2MAX,MM,MOM2
      COMMON/CGRA01/FF(MAXMAT,NQ),ERA(NEX),XSRA(MAXMAT,NEX),
     1    IED(NEGP),IEU(NEGP),NE2
      COMMON/CGRA02/QQ(NQ),AR(MAXMAT,NQ),BR(MAXMAT,NQ),CR(MAXMAT,NQ),
     1              DR(MAXMAT,NQ),FF0(MAXMAT),QQM
      PARAMETER (NP2=150)
      COMMON/CGRA03/QRA(NP2,MAXMAT),PRA(NP2,MAXMAT),DPRA(NP2,MAXMAT),
     1  ARA(NP2,MAXMAT),BRA(NP2,MAXMAT),PMAX(NEGP,MAXMAT),
     2  ITLRA(NP2,MAXMAT),ITURA(NP2,MAXMAT)

      COMMON/CGPP00/ZEQPP(MAXMAT),F0(MAXMAT,2),BCB(MAXMAT)
	  
	  COMMON/CPIN01/EI3,CPS3,BHA1,BHA2,BHA3,BHA4,MOM3
	  
	 
C
      INCLUDE 'pmcomms.f'

	  
	 CALL transfQTREE(NBODYS,KMOTH,KDGHT,KSURF, KFLAG,KSP, NWARN)

	 CALL transfQSURF(AXX,AXY,AXZ,AYY,AYZ,AZZ,
     1    AX,AY,AZ,A0,NSURF,KPLANE)
	 
	 CALL transfTRACK_mod(E ,X ,Y,Z,U,V,W,WGHT,
     1    SP1,SP2,SP3,PAGE,KPAR,IBODY,MAT,ILB,IPOL,LAGE)
	 
	 CALL transfPENELOPE_mod(EABS,C1,C2,WCC,WCR,DEN,
     1    RDEN,E0STEP,DESOFT,SSOFT,NMS,NEGP, NMAT)
	 
	 CALL transfPENGEOM_mod(BALIAS,DSTOT,MATER,KDET,
     1    KSLAST,NBODY,LVERB)
	 
	 CALL transfQBODY(KBODY, KBOMO)	
	 
	 CALL transfCECUTR(ECUTR)
	 
	 CALL transfCSGAWR(ISGAW)
	 
	 CALL transfCERSEC(IESERC)
	 
	 CALL transfCEGRID(EMIN ,EL ,EU,ET,DLEMP,DLEMP1,
     1    DLFC,XEL,XE,XEK,KE)
	 
	 CALL transfCESI0(XESI, IESIF, IESIL, NSESI, NCURE)
	 
	 CALL transfCPSI0(XPSI, IPSIF, IPSIL, NSPSI, NCURP)
	 
	 CALL transfCADATA(ATW,EPX,RSCR,ETA,EB,ALW,CP0,
     1    IFI, IKS, NSHT, LASYMB)
	 
	 CALL transfCGPH00(IPHF, IPHL, NPHS, NCUR, EPH, XPH)

c	 CALL transfcrelax(P,ET,F,IS0,IS1,IS2,IFIRST,ILAST,NCUR,KS,MODER)
	
	 CALL transfCRITAA(XA, AA, BA, FA, IA, NPM1A)

	 CALL transfCRITA(XT, PAC, DPAC, A, B, IL, IU, NPM1)	
	
C	 CALL RELAX0
	 
	 CALL transfCRITAN(CNORM)
	 
	 CALL transfCOMPOS(STF, ZT, AT, RHO, VMOL, IZ, NELEM)
	 
	 CALL transfCRANGE(RANGE, RANGEL)
	 
	 CALL transfCEIN(EXPOT, OP2, F, UI, WRI, KZ, KS, NOSC)
	 
	 CALL transfCEINTF(T1EI, T2EI, T1PI, T2PI)
	 
	 CALL transfCEIN00(SEH0,SEH1,SEH2,SES0,SES1,SES2,SET0,SET1,SET2)
	 
	 CALL transfCPIN00(SPH0,SPH1,SPH2,SPS0,SPS1,SPS2,SPT0,SPT1,SPT2)
	 
	 CALL transfCEINAC(EINAC, IEIN, NEIN)
	 
	 CALL transfCESIN(XSEIN, XSESI, ISIE)
	 
	 CALL transfCESIAC(ESIAC, IESI, NESI);
	 
	 CALL transfCPINAC(PINAC, IPIN, NPIN)
	 
	 CALL transfCPSIAC(PSIAC, IPSI, NPSI)
	 
	 CALL transfCPSIN(XSPIN, XSPSI, ISIP)
	 
	 CALL transfCGCO(FCO,UICO ,FJ0,
     1  PTRSH ,KZCO ,KSCO ,
     2  NOSCCO)
	 
	 CALL transfCEIMFP(SEHEL,SEHIN,
     1  SEISI,SEHBR,SEAUX,
     2  SETOT,CSTPE,RSTPE,
     3  DEL,W1E,W2E,
     4  DW1EL,DW2EL,
     5  RNDCE,AE,BE,
     6  T1E,T2E)
	 
	 CALL transfCLAS1E(TSTPE,TSTRE,
     1  TRL1E,TRL2E)
	 
	 CALL transfCPIMFP(SPHEL,SPHIN,
     1  SPISI,SPHBR,SPAN,
     2  SPAUX,SPTOT,CSTPP,
     3  RSTPP,W1P,W2P,
     4  DW1PL,DW2PL,
     5  RNDCP,AP,BP,
     6  T1P,T2P)
	 
	 CALL transfCLAS1P(TSTPP,TSTRP,
     1  TRL1P,TRL2P)
	 
	 CALL transfCBRYLD(EBRY, PBRY)
	 
	 CALL transfCGIMFP(SGRA, SGCO, SGPH, SGPP, SGAUX)
	 
	 CALL transfCGPH01(ER, XSR, NPHD)
	 
	 CALL transfCGPP01(TRIP)
	 
	 CALL transfCEBR(WB, PBCUT, WBCUT, PDFB, DPDFB, PACB, ZBR2)
	 
	 CALL transfCEBR01(EBT,XS,TXS,X1,Y1)
	 
	 CALL transfCEBR02(P0)
	 
	 CALL transfCEEL00(EJT,XE0,XE1,XE2,XP0,
     1  XP1,XP2,T1E0,T2E0,T1P0,T2P0,
     2  EJTL,FJL,A2,B2,C,D)
	 
	 CALL transfcbrang(BET, BK, BP1, BP2, ZBEQ)
	 
	 CALL transfcein01(EI,EE,CPS,AMOL,MOM)
	 
	 CALL transfcsumga(IERGA,NCALL)
	 
c	 CALL PINaT1(0.e0,0.e0,0.e0,0.e0,0.e0,0.e0,0.e0,0.e0,
c	 1	0.e0,0.e0,0.e0,0.e0,0.e0,0.e0)

c	 CALL PINaT1(0.e0,0.e0,0.e0,0.e0,0.e0,0.e0,0.e0,0.e0,
c    1  0.e0,0.e0,0.e0,0.e0,0.e0,0.e0)
	 
C	 CALL PINaT1(E,UK,WK,DELTA,WCCM,H0,H1,H2,S0,S1,S2,R0,R1,R2)
	 
	 
	 CALL transfcdcsep(ETS,ETL,TH,THR,XMU,XMUL,ECS,ETCS1,
     1  ETCS2,EDCS,PCS,PTCS1,PTCS2,PDCS,
     2  DCSI,DCSIL,CSI,TCS1I,TCS2I)
	 
	 CALL transfcEELDB(XSE,PSE,ASE,BSE,ITLE,ITUE)
	 
	 CALL transfcPELDB(XSP,PSP,ASP,BSP,ITLP,ITUP)
	 
	 CALL transfcELSEP(EELMAX,PELMAX,RNDCEd,RNDCPd)
	 
	 CALL transfcgra00(FACTE,Q2MAX,MM,MOM2)
	 
	 CALL transfcgra01(FF,ERA,XSRA,IED,IEU,NE2)
	 
	 CALL transfcgra02(QQ,AR,BR,CR,DR,FF0,QQM)
	 
	 CALL transfcgra03(QRA,PRA,DPRA,
     1  ARA,BRA,PMAX,
     2  ITLRA,ITURA)
	 
	 CALL transfcgpp00(ZEQPP, F0, BCB)
	 
	 CALL transfCPIN01(EI3,CPS3,BHA1,BHA2,BHA3,BHA4,MOM3)
	 

	 
C
C  ****  Read input files and initialise the simulation packages.
C



      CALL PMRDR
	  
	  CALL TABELAS
	  
 
	  
	  
	  
	  
      IF(JOBEND.NE.0) GO TO 103 ! The simulation was already completed.
	  
	 
	  
	  
	  
C
C  ****  Simulation of a new shower and scoring.
C
  101 CONTINUE
      CALL SHOWER
      IF(JOBEND.NE.0) GO TO 102  ! The simulation is completed.
C
C  ****  End the simulation after the allotted time or after completing
C        DSHN showers.
C

      CALL TIMER(TSEC)
      IF(TSEC.LT.TSECA.AND.SHN.LT.DSHN) THEN
C  ****  Write partial results after each dumping period.
        IF(LDUMP) THEN
          IF(TSEC-TSECAD.GT.DUMPP) THEN
            TSIM=TSIM+CPUTIM()-CPUT0
            CALL PMWRT(-1)
            WRITE(6,1001) SHN
 1001       FORMAT(3X,'Number of simulated showers =',1P,E14.7)
            CALL TIMER(TSEC)
            TSECAD=TSEC
            CPUT0=CPUTIM()
            GO TO 101
          ENDIF
        ENDIF
        GO TO 101
      ENDIF
C
  102 CONTINUE
      TSIM=TSIM+CPUTIM()-CPUT0
  103 CONTINUE
      CALL PMWRT(1)
      WRITE(6,1002) SHN
 1002 FORMAT(3X,'Number of simulated showers =',1P,E14.7,
     1  /3X,'*** END ***')
C
      WRITE(26,'(/3X,''***  END  ***'')')
      WRITE(26,'(/3X,72(''-''))')
      CLOSE(26)
C
      STOP
      END


C  *********************************************************************
C                       SUBROUTINE PMRDR
C  *********************************************************************
      SUBROUTINE PMRDR
C
C  Reads the input file and initialises PENELOPE and PENGEOM.
C
C
      USE PENELOPE_mod
      USE TRACK_mod
      USE PENGEOM_mod
      USE PENVARED_mod
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
      CHARACTER LIT*2,WORD*4
      CHARACTER*20 PMFILE,PFILE,PFILER
      CHARACTER*20 SPCDIO,SPCDEO,PSFDIO,SPCFSO,SPCAGE
      CHARACTER*120 BUFFER,BUF2
C
      CHARACTER*6 KWORD,
     1  KWTITL,KWKPAR,KWSENE,KWSPEC,  KWSPOL,KWSPOS,KWSBOX,KWSBOD,
     1  KWSCON,KWSREC,KWPSFN,KWPSPL,  KWRRSP,KWEMAX,KWMATF,KWSIMP,
     1  KWGEOM,KWGPAR,KWSMAX,KWEABS,  KWIFOR,KBRSPL,KXRSPL,KWNBE,
     1  KWNBAN,KWIDET,KWIPSF,KWISPC,  KWIFLN,KWDIAL,KWDIAF,KWIBOD,
     1  KWIPAR,KWEDET,KWESPC,KWEBOD,  KGRDXX,KGRDYY,KGRDZZ,KGRDRR,
     1  KWRESU,KWDUMP,KWDMPP,KWRSEE,  KWNSIM,KWTIME,KWCOMM
      PARAMETER (
     1  KWTITL='TITLE ',KWKPAR='SKPAR ',KWSENE='SENERG',KWSPEC='SPECTR',
     1  KWSPOL='SGPOL ',KWSPOS='SPOSIT',KWSBOX='SBOX  ',KWSBOD='SBODY ',
     1  KWSCON='SCONE ',KWSREC='SRECTA',KWPSFN='IPSFN ',KWPSPL='IPSPLI',
     1  KWRRSP='WGTWIN',KWEMAX='EPMAX ',KWMATF='MFNAME',KWSIMP='MSIMPA',
     1  KWGEOM='GEOMFN',KWGPAR='PARINP',KWSMAX='DSMAX ',KWEABS='EABSB ',
     1  KWIFOR='IFORCE',KBRSPL='IBRSPL',KXRSPL='IXRSPL',KWNBE ='NBE   ',
     1  KWNBAN='NBANGL',KWIDET='IMPDET',KWIPSF='IDPSF ',KWISPC='IDSPC ',
     1  KWIFLN='IDFLNC',KWDIAL='IDAGEL',KWDIAF='IDAGEF',KWIBOD='IDBODY',
     1  KWIPAR='IDKPAR',KWEDET='ENDETC',KWESPC='EDSPC ',KWEBOD='EDBODY',
     1  KGRDXX='GRIDX ',KGRDYY='GRIDY ',KGRDZZ='GRIDZ ',KGRDRR='GRIDR ',
     1  KWRESU='RESUME',KWDUMP='DUMPTO',KWDMPP='DUMPP ',KWRSEE='RSEED ',
     1  KWNSIM='NSIMSH',KWTIME='TIME  ',KWCOMM='      ')
C
      PARAMETER (PI=3.1415926535897932D0, DE2RA=PI/180.0D0)
      PARAMETER (FSAFE=1.000000001D0)
C
C  ----  Seeds of the random number generator.
      COMMON/RSEED/ISEED1,ISEED2
C
      INCLUDE 'pmcomms.f'
C
      PARAMETER (NPINPM=500)
      DIMENSION PARINP(NPINPM)
      DATA PARINP/NPINPM*0.0D0/
      DIMENSION PMFILE(MAXMAT)
C
      OPEN(26,FILE='penmain.dat')  ! Global output/message file.
C
      DO I=1,3
        PRIM(I)=0.0D0
        PRIM2(I)=0.0D0
        DO K=1,3
          SEC(K,I)=0.0D0
          SEC2(K,I)=0.0D0
        ENDDO
      ENDDO
C
      DO I=1,2
        AVW(I)=0.0D0
        AVW2(I)=0.0D0
        AVA(I)=0.0D0
        AVA2(I)=0.0D0
        AVE(I)=0.0D0
        AVE2(I)=0.0D0
      ENDDO
C
      DO I=1,NSEM
        SHIST(I)=0.0D0
        DO K=1,3
          SEDS(K,I)=0.0D0
          SEDS2(K,I)=0.0D0
        ENDDO
      ENDDO
C
      DO I=1,NB
        TDEBO(I)=0.0D0
        TDEBO2(I)=0.0D0
        EABSB(1,I)=50.0D0
        EABSB(2,I)=50.0D0
        EABSB(3,I)=50.0D0
      ENDDO

      CTHL=0.0D0
      DCTH=0.0D0
      PHIL=0.0D0
      DPHI=0.0D0
      DO KB=1,NB
        IXSBOD(KB)=0
      ENDDO
C
C  ****  Time counter initiation.
C
      CALL TIME0
C
C  ------------------------  Read input data file.
C
      WRITE(26,1000)
 1000 FORMAT(//3X,61('*'),/3X,'**   Program PENMAIN. ',
     1 ' Input data and run-time messages.   **',/3X,61('*'))
C
      CALL PDATET(DATE23)
      WRITE(26,1001) DATE23
 1001 FORMAT(/3X,'Date and time: ',A23)
C
C  ****  Title.
C
      READ(5,'(A6,1X,A65)') KWORD,TITLE
      IF(KWORD.EQ.KWTITL) THEN
        WRITE(26,'(/3X,A65)') TITLE
      ELSE
        WRITE(26,*) 'The input file must begin with the TITLE line.'
        STOP 'The input file must begin with the TITLE line.'
      ENDIF
C
      KPARP=1
      NPSF=0
      NPSN=0
      NSEB=1
      NSPLIT=1
      RLREAD=0.0D0
      RWGMIN=1.0D35
      KBSMAX=0
      ISEC=0
      LDUMP=.FALSE.
      LPSF=.FALSE.
      LSPEC=.FALSE.
      LGPOL=.FALSE.
      LEXSRC=.FALSE.
      LEXBD=.FALSE.
      LSCONE=.TRUE.
      JOBEND=0
C
C  ************  Source description.
C
   11 CONTINUE
      READ(5,'(A6,1X,A120)') KWORD,BUFFER
      IF(KWORD.EQ.KWCOMM) GO TO 11
      IF(KWORD.EQ.KWPSFN) GO TO 21
C
      WRITE(26,1100)
 1100 FORMAT(/3X,72('-'),/3X,'>>>>>>  Source description.')
C
      IF(KWORD.EQ.KWKPAR) THEN
        READ(BUFFER,*) KPARP
   12   CONTINUE
        READ(5,'(A6,1X,A120)') KWORD,BUFFER
        IF(KWORD.EQ.KWCOMM) GO TO 12
        IF(KWORD.EQ.KWPSFN) GO TO 21
      ENDIF
      IF(KPARP.LT.0.OR.KPARP.GT.3) THEN
        WRITE(26,*) 'KPARP =',KPARP
        WRITE(26,*) 'Incorrect particle type.'
        STOP 'Incorrect particle type.'
      ENDIF
      IF(KPARP.EQ.1) WRITE(26,1110)
 1110 FORMAT(/3X,'Primary particles: electrons')
      IF(KPARP.EQ.2) WRITE(26,1111)
 1111 FORMAT(/3X,'Primary particles: photons')
      IF(KPARP.EQ.3) WRITE(26,1112)
 1112 FORMAT(/3X,'Primary particles: positrons')
      IF(KPARP.EQ.0) WRITE(26,1113)
 1113 FORMAT(/3X,'Primary particles: set by the user subroutine SOURCE')
C
C  ****  Monoenergetic source.
C
      IF(KWORD.EQ.KWSENE) THEN
        READ(BUFFER,*) E0
        WRITE(26,1120) E0
 1120   FORMAT(3X,'Initial energy = ',1P,E13.6,' eV')
   13   CONTINUE
        READ(5,'(A6,1X,A120)') KWORD,BUFFER
        IF(KWORD.EQ.KWCOMM) GO TO 13
        IF(KWORD.EQ.KWPSFN) GO TO 21
C
C  ****  Continuous energy spectrum.
C
      ELSE IF(KWORD.EQ.KWSPEC) THEN
        LSPEC=.TRUE.
        NSEB=0
   14   CONTINUE
        NSEB=NSEB+1
        IF(NSEB.GT.NSEM) THEN
          WRITE(26,*) 'Source energy spectrum.'
          WRITE(26,*) 'The number of energy bins is too large.'
          STOP 'The number of energy bins is too large.'
        ENDIF
        READ(BUFFER,*) ESRC(NSEB),PSRC(NSEB)
        PSRC(NSEB)=MAX(PSRC(NSEB),0.0D0)
   15   CONTINUE
        READ(5,'(A6,1X,A120)') KWORD,BUFFER
        IF(KWORD.EQ.KWCOMM) GO TO 15
        IF(KWORD.EQ.KWSPEC) GO TO 14
        IF(KWORD.EQ.KWPSFN) GO TO 21
      ELSE
        E0=1.0D9
        IF(KPARP.NE.0) WRITE(26,1120) E0
      ENDIF
      IF(LSPEC) THEN
        IF(NSEB.GT.1) THEN
          CALL SORT2(ESRC,PSRC,NSEB)
          WRITE(26,1121)
 1121     FORMAT(/3X,'Spectrum:',7X,'I',4X,'E_low(eV)',4x,'E_high(eV)',
     1      5X,'P_sum(E)',/16X,45('-'))
          DO I=1,NSEB-1
            WRITE(26,'(16X,I4,1P,3E14.6)') I,ESRC(I),ESRC(I+1),PSRC(I)
          ENDDO
          WRITE(26,*) '  '
          E0=ESRC(NSEB)
          NSEB=NSEB-1
          CALL IRND0(PSRC,FSRC,IASRC,NSEB)
        ELSE
          WRITE(26,*) 'The source energy spectrum is not defined.'
          STOP 'The source energy spectrum is not defined.'
        ENDIF
      ENDIF
      IF(E0.LT.50.0D0) THEN
        WRITE(26,*) 'The initial energy E0 is too small.'
        STOP 'The initial energy E0 is too small.'
      ENDIF
      EPMAX=E0
C  ----  Positrons eventually give annihilation gamma rays. The maximum
C        energy of annihilation photons is .lt. 1.21*(E0+me*c**2).
      IF(KPARP.EQ.3) EPMAX=1.21D0*(E0+5.12D5)
C
C  ****  Photon polarisation effects (only for primary photons).
C
      IF(KWORD.EQ.KWSPOL) THEN
        READ(BUFFER,*) SP10,SP20,SP30
        LGPOL=.TRUE.
   20   CONTINUE
        READ(5,'(A6,1X,A120)') KWORD,BUFFER
        IF(KWORD.EQ.KWCOMM) GO TO 20
        IF(KWORD.EQ.KWPSFN) GO TO 21
        WRITE(26,1131) SP10,SP20,SP30
 1131   FORMAT(/3X,'Polarised primary photons. Stokes Parameters:',
     1    /6X,'P1 = ',E13.6,' (linear polarisation at 45 deg azimuth)',
     1    /6X,'P2 = ',E13.6,' (circular polarisation)',
     1    /6X,'P3 = ',E13.6,' (linear polarisation at zero azimuth)')
      ENDIF
C
C  ****  Position of the point source.
C
      IF(KWORD.EQ.KWSPOS) THEN
        READ(BUFFER,*) SX0,SY0,SZ0
   16   CONTINUE
        READ(5,'(A6,1X,A120)') KWORD,BUFFER
        IF(KWORD.EQ.KWCOMM) GO TO 16
        IF(KWORD.EQ.KWPSFN) GO TO 21
      ELSE
        SX0=0.0D0
        SY0=0.0D0
        SZ0=0.0D0
      ENDIF
      WRITE(26,1132) SX0,SY0,SZ0
 1132 FORMAT(/3X,'Coordinates of centre:     SX0 = ',1P,E13.6,
     1  ' cm',/30X,'SY0 = ',E13.6,' cm',/30X,'SZ0 = ',E13.6,' cm')
C
C  ****  Extended source.
C
      IF(KWORD.EQ.KWSBOX) THEN
        LEXSRC=.TRUE.
        READ(BUFFER,*) SSX,SSY,SSZ
        SSX=ABS(SSX)
        SSY=ABS(SSY)
        SSZ=ABS(SSZ)
   17   CONTINUE
        READ(5,'(A6,1X,A120)') KWORD,BUFFER
        IF(KWORD.EQ.KWCOMM) GO TO 17
        IF(KWORD.EQ.KWPSFN) GO TO 21
      ELSE
        LEXSRC=.FALSE.
        SSX=0.0D0
        SSY=0.0D0
        SSZ=0.0D0
      ENDIF
      IF(LEXSRC) WRITE(26,1135) SSX,SSY,SSZ
 1135 FORMAT(3X,'Source size:',15X,'SSX = ',1P,E13.6,' cm',
     1  /30X,'SSY = ',E13.6,' cm',/30X,'SSZ = ',E13.6,' cm')
C
C  ****  Active bodies of an extended source (PENGEOM internal labels).
C
  717 CONTINUE
      IF(KWORD.EQ.KWSBOD) THEN
        READ(BUFFER,*) KB
        IF(KB.LT.1.OR.KB.GT.NB) THEN
          WRITE(26,'(A6,1X,A120)') KWORD,BUFFER
          WRITE(26,*) 'Incorrect body label.'
          STOP 'Incorrect body label.'
        ENDIF
        WRITE(26,7620) KB
 7620   FORMAT(3X,'Active body = ',I4)
        IXSBOD(KB)=1
        LEXBD=.TRUE.
        KBSMAX=MAX(KBSMAX,KB)
  766   CONTINUE
        READ(5,'(A6,1X,A120)') KWORD,BUFFER
        IF(KWORD.EQ.KWCOMM) GO TO 766
        IF(KWORD.EQ.KWSBOD) GO TO 717
      ENDIF
C
C  ****  Angular distribution of primary particles.
C
      ISOURC=0
  777 CONTINUE
      IF(KWORD.EQ.KWSCON) THEN
        LSCONE=.TRUE.
        READ(BUFFER,*) STHETA,SPHI,SALPHA
        IF(STHETA.LT.-1.0D-9.OR.STHETA-180.0D0.GT.1.0D-9) THEN
          WRITE(26,'(A6,1X,A120)') KWORD,BUFFER
          WRITE(26,'(/3X,''THETA must be between 0 and 180 deg.'')')
          STOP 'THETA must be between 0 and 180 deg'
        ENDIF
        IF(SPHI.LT.-1.0D-9.OR.SPHI-360.0D0.GT.1.0D-9) THEN
          WRITE(26,'(A6,1X,A120)') KWORD,BUFFER
          WRITE(26,'(/3X,''PHI must be between 0 and 360 deg.'')')
          STOP 'PHI must be between 0 and 360 deg'
        ENDIF
        IF(SALPHA.LT.-1.0D-9.OR.SALPHA-180.0D0.GT.1.0D-9) THEN
          WRITE(26,'(A6,1X,A120)') KWORD,BUFFER
          WRITE(26,'(/3X,''ALPHA must be between 0 and 180 deg.'')')
          STOP 'ALPHA must be between 0 and 180 deg'
        ENDIF
        CALL GCONE0(STHETA*DE2RA,SPHI*DE2RA,SALPHA*DE2RA)
        ISOURC=1
   18   CONTINUE
        READ(5,'(A6,1X,A120)') KWORD,BUFFER
        IF(KWORD.EQ.KWCOMM) GO TO 18
        IF(KWORD.EQ.KWPSFN) GO TO 721
        GO TO 777
      ELSE IF(KWORD.EQ.KWSREC) THEN
        LSCONE=.FALSE.
        READ(BUFFER,*) THETLD,THETUD,PHILD,PHIUD
        IF(MIN(THETLD,THETUD).LT.-1.0D-9.OR.
     1     MAX(THETLD,THETUD)-180.0D0.GT.1.0D-9) THEN
          WRITE(26,'(A6,1X,A120)') KWORD,BUFFER
          WRITE(26,'(/3X,''THETA must be between 0 and 180 deg.'')')
          STOP 'THETA must be between 0 and 180 deg'
        ENDIF
        IF(MIN(PHILD,PHIUD).LT.-1.0D-9.OR.
     1     MAX(PHILD,PHIUD)-360.0D0.GT.1.0D-9) THEN
          WRITE(26,'(A6,1X,A120)') KWORD,BUFFER
          WRITE(26,'(/3X,''PHI must be between 0 and 360 deg.'')')
          STOP 'PHI must be between 0 and 360 deg'
        ENDIF
        CTHL=COS(THETLD*DE2RA)
        DCTH=COS(THETUD*DE2RA)-CTHL
        PHIL=PHILD*DE2RA
        DPHI=PHIUD*DE2RA-PHIL
        ISOURC=1
   19   CONTINUE
        READ(5,'(A6,1X,A120)') KWORD,BUFFER
        IF(KWORD.EQ.KWCOMM) GO TO 19
        IF(KWORD.EQ.KWPSFN) GO TO 721
        GO TO 777
      ELSE IF(ISOURC.EQ.0) THEN
        LSCONE=.TRUE.
        STHETA=0.0D0
        SPHI=0.0D0
        SALPHA=0.0D0
        CALL GCONE0(STHETA*DE2RA,SPHI*DE2RA,SALPHA*DE2RA)
      ENDIF
C
  721 CONTINUE
      IF(LSCONE) THEN
        WRITE(26,1133) STHETA,SPHI
 1133   FORMAT(/3X,'*** Conical beam:'
     1    /3X,'Beam axis direction:     THETA = ',1P,E13.6,' deg',
     1    /30X,'PHI = ',E13.6,' deg')
        WRITE(26,1134) SALPHA
 1134   FORMAT(3X,'Beam aperture:',11X,'ALPHA = ',1P,E13.6,' deg')
      ELSE
        WRITE(26,1733) THETLD,THETUD,PHILD,PHIUD
 1733   FORMAT(/3X,'*** Rectangular beam:'
     1    /3X,'Angular window: THETA = (',1P,E13.6,',',E13.6,') deg',
     1    /21X,'PHI = (',E13.6,',',E13.6,') deg')
      ENDIF
C
C  ************  Particle state variables read from a phase-space file.
C
   21 CONTINUE
      IF(KWORD.EQ.KWPSFN) THEN
        IF(KPARP.EQ.0) THEN
          WRITE(26,'(/3X,''With KPARP=0 (subroutine SOURCE activat'',
     1      ''ed),'')')
          WRITE(26,'(3X,''we cannot read particles from a phase-'',
     1      ''space file.'')')
          STOP 'Inconsistent definition of the primary source.'
        ENDIF
        NPSF=NPSF+1
        IF(NPSF.EQ.1) THEN
          WRITE(26,1200)
 1200     FORMAT(/3X,72('-'),/3X,'>>>>>>  Input phase-space files.'/)
        ENDIF
        IF(NPSF.GT.NPSFM) THEN
          WRITE(26,'(/3X,''Too many phase-space files.'')')
          STOP 'Too many phase-space files'
        ENDIF
        READ(BUFFER,'(A20)') PSFI(NPSF)
        WRITE(26,1201) NPSF,PSFI(NPSF)
 1201   FORMAT(3X,'Phase-space file #',I4,': ',A20)
        LPSF=.TRUE.
        LSPEC=.FALSE.
   22   CONTINUE
        READ(5,'(A6,1X,A120)') KWORD,BUFFER
        IF(KWORD.EQ.KWCOMM) GO TO 22
        IF(KWORD.EQ.KWPSFN) GO TO 21
C
        IF(KWORD.EQ.KWPSPL) THEN
          READ(BUFFER,*) NSPLIT
          IF(NSPLIT.GT.1.AND.NSPLIT.LE.1000) THEN
            WRITE(26,1210) NSPLIT
 1210       FORMAT(/3X,'Particle splitting number = ',I3)
          ELSE
            NSPLIT=MIN(1000,ABS(NSPLIT))
            WRITE(26,1211) NSPLIT
 1211       FORMAT(/3X,'Particle splitting number = ',I3,
     1        ' (modified)')
          ENDIF
   23     CONTINUE
          READ(5,'(A6,1X,A120)') KWORD,BUFFER
          IF(KWORD.EQ.KWCOMM) GO TO 23
        ENDIF
C
        IF(KWORD.EQ.KWRRSP) THEN
          READ(BUFFER,*) WGMIN,WGMAX
          WGMIN=ABS(WGMIN)
          WGMAX=MIN(ABS(WGMAX),1.0D10)
          IF(WGMIN.GT.WGMAX) THEN
            WRITE(26,*) 'WGMIN =',WGMIN
            WRITE(26,*) 'WGMAX =',WGMAX
            WRITE(26,'(/3X,''Inconsistent window end points.'')')
            STOP 'Inconsistent window end points.'
          ENDIF
          WRITE(26,1291) WGMIN,WGMAX
 1291     FORMAT(/3X,'Initial weight window = (',1P,E13.6,',',
     1      E13.6,')')
   24     CONTINUE
          READ(5,'(A6,1X,A120)') KWORD,BUFFER
          IF(KWORD.EQ.KWCOMM) GO TO 24
        ENDIF
        RWGMIN=1.0D0/WGMIN
      ENDIF
C
      IF(KPARP.EQ.0) THEN
        CALL SOURCE(METAST)
        NSDE=400
        DSDE=FSAFE*EPMAX/DBLE(NSDE)
        RDSDE=1.0D0/DSDE
      ENDIF
C
C  ****  Maximum particle energy.
C
      IF(KWORD.EQ.KWEMAX) THEN
        READ(BUFFER,*) EPMAXR
        IF(KPARP.NE.0) EPMAX=EPMAXR
        WRITE(26,1220) EPMAX
 1220   FORMAT(/3X,'Maximum particle energy = ',1P,E13.6,' eV')
   25   CONTINUE
        READ(5,'(A6,1X,A120)') KWORD,BUFFER
        IF(KWORD.EQ.KWCOMM) GO TO 25
      ELSE
        IF(LPSF) THEN
          EPMAX=1.0D9
          WRITE(26,1220) EPMAX
          WRITE(26,'(3X,''WARNING: You should have specified the '',
     1      ''maximum energy EPMAX.'')')
        ENDIF
        WRITE(26,1220) EPMAX
      ENDIF
C
C  ************  Material data and simulation parameters.
C
      WRITE(26,1300)
 1300 FORMAT(/3X,72('-'),/
     1  3X,'>>>>>>  Material data and simulation parameters.')
C
C  ****  Simulation parameters.
C
      DO M=1,MAXMAT
        EABS(1,M)=0.010D0*EPMAX
        EABS(2,M)=0.001D0*EPMAX
        EABS(3,M)=0.010D0*EPMAX
        C1(M)=0.10D0
        C2(M)=0.10D0
        WCC(M)=EABS(1,M)
        WCR(M)=EABS(2,M)
      ENDDO
      DO IB=1,NB
        DSMAX(IB)=1.0D20
      ENDDO
C
      NMATR=0
   31 CONTINUE
      IF(KWORD.EQ.KWMATF) THEN
        NMATR=NMATR+1
        READ(BUFFER,'(A20)') PMFILE(NMATR)
   32   CONTINUE
        READ(5,'(A6,1X,A120)') KWORD,BUFFER
        IF(KWORD.EQ.KWCOMM) GO TO 32
        IF(KWORD.EQ.KWMATF) GO TO 31
      ENDIF
C
      IF(KWORD.EQ.KWSIMP) THEN
        READ(BUFFER,*) EABS(1,NMATR),EABS(2,NMATR),EABS(3,NMATR),
     1    C1(NMATR),C2(NMATR),WCC(NMATR),WCR(NMATR)
   33   CONTINUE
        READ(5,'(A6,1X,A120)') KWORD,BUFFER
        IF(KWORD.EQ.KWCOMM) GO TO 33
        IF(KWORD.EQ.KWMATF) GO TO 31
      ENDIF
C
      IF(NMATR.EQ.0) THEN
        WRITE(26,'(A6,1X,A120)') KWORD,BUFFER
        WRITE(26,*) 'You have to specify a material file (line MFNAME).'
        STOP 'You have to specify a material file (line MFNAME).'
      ENDIF
      IF(NMATR.GT.MAXMAT) THEN
        WRITE(26,*) 'Wrong number of materials.'
        WRITE(26,'(''NMAT ='',I4,'' is larger than MAXMAT ='',I4)')
     1    NMATR,MAXMAT
        STOP 'Wrong number of materials.'
      ENDIF
C
      DO M=1,NMATR
        IF(M.EQ.1) LIT='st'
        IF(M.EQ.2) LIT='nd'
        IF(M.EQ.3) LIT='rd'
        IF(M.GT.3) LIT='th'
        WRITE(26,1320) M,LIT
 1320   FORMAT(/3X,'**** ',I2,A2,' material')
        WRITE(26,1325) PMFILE(M)
 1325   FORMAT(3X,'Material data file: ',A)
        IF(EABS(1,M).LT.5.0D1) EABS(1,M)=5.0D1
        IF(EABS(2,M).LT.5.0D1) EABS(2,M)=5.0D1
        IF(EABS(3,M).LT.5.0D1) EABS(3,M)=5.0D1
        WRITE(26,1321) EABS(1,M)
 1321   FORMAT(3X,'Electron absorption energy = ',1P,E13.6,' eV')
        WRITE(26,1322) EABS(2,M)
 1322   FORMAT(3X,'  Photon absorption energy = ',1P,E13.6,' eV')
        WRITE(26,1323) EABS(3,M)
 1323   FORMAT(3X,'Positron absorption energy = ',1P,E13.6,' eV')
        WRITE(26,1324) C1(M),C2(M),WCC(M),WCR(M)
 1324   FORMAT(3X,'Electron-positron simulation parameters:',
     1    /4X,'C1 =',1P,E13.6,',      C2 =',E13.6,/3X,'Wcc =',E13.6,
     1    ' eV,  Wcr =',E13.6,' eV')
      ENDDO
C
C  ****  Initialisation of PENELOPE.
C
      WRITE(6,*) '  Initialising PENELOPE ...'
      IWR=16
      OPEN(IWR,FILE='material.dat')
        INFO=3
        CALL PEINIT(EPMAX,NMATR,IWR,INFO,PMFILE)
      CLOSE(IWR)
C
C  ************  Geometry definition.
C
      WRITE(6,*) '  Initialising PENGEOM ...'
      IF(KWORD.EQ.KWGEOM) THEN
        READ(BUFFER,'(A20)') PFILE
        WRITE(26,1340) PFILE
 1340   FORMAT(/3X,72('-'),/3X,'>>>>>>  Geometry definition.',
     1    /3X,'PENGEOM''s geometry file: ',A20)
        OPEN(15,FILE=PFILE,IOSTAT=KODE)
        IF(KODE.NE.0) THEN
          WRITE(26,'(''File '',A20,'' could not be opened.'')') PFILE
          STOP 'Geometry file could not be opened.'
        ENDIF
C
        NPINP=0
        IHEAD=0
   34   CONTINUE
        READ(5,'(A6,1X,A120)') KWORD,BUFFER
        IF(KWORD.EQ.KWCOMM) GO TO 34
        IF(KWORD.EQ.KWGPAR) THEN
          READ(BUFFER,*) IP
          IF(IP.LT.1) THEN
            WRITE(26,'(''IP = '',I4)') IP
            STOP 'The PARINP index must be positive.'
          ENDIF
          NPINP=MAX(NPINP,IP)
          IF(NPINP.GT.NPINPM) THEN
            WRITE(26,'(''Too many modified parameters.'')')
            WRITE(26,'(''NPINP = '',I4,'' must be less than'',I4)')
     1        NPINP,NPINPM
            STOP 'Too many modified parameters.'
          ENDIF
          READ(BUFFER,*) I,PARINP(IP)
          IF(IHEAD.EQ.0) THEN
            WRITE(26,1341) IP,PARINP(IP)
 1341       FORMAT(/3X,'Replaced parameters: PARINP(',I4,')= ',
     1        1P,E13.6)
            IHEAD=1
          ELSE
            WRITE(26,1342) IP,PARINP(IP)
 1342       FORMAT(24X,'PARINP(',I4,')= ',1P,E13.6)
          ENDIF
          GO TO 34
        ENDIF
C
        OPEN(16,FILE='geometry.rep')
        CALL GEOMIN(PARINP,NPINP,NMATG,NBODY,15,16)
        CLOSE(15)
        CLOSE(16)
        IF(NMATG.LT.1) THEN
          WRITE(26,*) 'NMATG must be greater than 0.'
          STOP 'NMATG must be greater than 0.'
        ENDIF
C
        IF(NBODY.GT.NB) THEN
          WRITE(26,'(/6X,''Too many bodies.'')')
          STOP 'Too many bodies.'
        ENDIF
C
        IF(NMATG.GT.NMAT) THEN
          WRITE(26,'(/6X,''Too many different materials.'')')
          STOP 'Too many different materials.'
        ENDIF
C
        IF(KBSMAX.GT.NBODY) THEN
          WRITE(26,'(/6X,''KBSMAX = '',I4)') KBSMAX
          WRITE(26,'(6X,'' NBODY = '',I4)') NBODY
          WRITE(26,'(6X,''Some source bodies are undefined. STOP.'')')
          STOP 'Some source bodies are undefined.'
        ENDIF
      ELSE
        WRITE(26,'(A6,1X,A120)') KWORD,BUFFER
        WRITE(26,*) 'You have to specify a geometry file.'
        STOP 'You have to specify a geometry file.'
      ENDIF
C
C  ****  Maximum step lengths of electrons and positrons.
C
      IF(KWORD.EQ.KWSMAX) THEN
        CALL FWORD(BUFFER,BUF2,WORD,LENGTH)
        IF(LENGTH.EQ.4) THEN
          DO IB=1,NBODY
            IF(WORD.EQ.BALIAS(IB)) KB=IB
          ENDDO
          READ(BUF2,*) DSMAX(KB)
        ELSE
          READ(BUFFER,*) KB,DSMAX(KB)
        ENDIF
        IF(KB.LT.1.OR.KB.GT.NBODY) THEN
          WRITE(26,'(A6,1X,A120)') KWORD,BUFFER
          WRITE(26,*) 'Incorrect body number.'
          STOP 'Incorrect body number.'
        ENDIF
        IF(DSMAX(KB).LT.1.0D-7) DSMAX(KB)=1.0D20
   35   CONTINUE
        READ(5,'(A6,1X,A120)') KWORD,BUFFER
        IF(KWORD.EQ.KWCOMM) GO TO 35
        IF(KWORD.EQ.KWSMAX) THEN
          CALL FWORD(BUFFER,BUF2,WORD,LENGTH)
          IF(LENGTH.EQ.4) THEN
            DO IB=1,NBODY
              IF(WORD.EQ.BALIAS(IB)) KB=IB
            ENDDO
            READ(BUF2,*) DSMAX(KB)
          ELSE
            READ(BUFFER,*) KB,DSMAX(KB)
          ENDIF
          IF(KB.LT.1.OR.KB.GT.NBODY) THEN
            WRITE(26,'(A6,1X,A120)') KWORD,BUFFER
            WRITE(26,*) 'Incorrect body number.'
            STOP 'Incorrect body number.'
          ENDIF
          IF(DSMAX(KB).LT.1.0D-7) DSMAX(KB)=1.0D20
          GO TO 35
        ENDIF
      ENDIF
C
C  ****  Local absorption energies (useful to reduce simulation work
C        in regions of lesser interest).
C
      DO IB=1,NBODY
        M=MATER(IB)
        IF(M.GT.0) THEN
          EABSB(1,IB)=EABS(1,M)
          EABSB(2,IB)=EABS(2,M)
          EABSB(3,IB)=EABS(3,M)
        ENDIF
      ENDDO
C
      IF(KWORD.EQ.KWEABS) THEN
   36   CONTINUE
        CALL FWORD(BUFFER,BUF2,WORD,LENGTH)
        IF(LENGTH.EQ.4) THEN
          DO IB=1,NBODY
            IF(WORD.EQ.BALIAS(IB)) KB=IB
          ENDDO
          READ(BUF2,*) EAB1,EAB2,EAB3
        ELSE
          READ(BUFFER,*) KB,EAB1,EAB2,EAB3
        ENDIF
        IF(KB.LT.1.OR.KB.GT.NBODY) THEN
          WRITE(26,'(A6,1X,A120)') KWORD,BUFFER
          WRITE(26,*) 'Incorrect body number.'
          STOP 'Incorrect body number.'
        ENDIF
        IF(MATER(KB).GT.0) THEN
          EABSB(1,KB)=MAX(EABSB(1,KB),EAB1)
          EABSB(2,KB)=MAX(EABSB(2,KB),EAB2)
          EABSB(3,KB)=MAX(EABSB(3,KB),EAB3)
        ENDIF
   37   CONTINUE
        READ(5,'(A6,1X,A120)') KWORD,BUFFER
        IF(KWORD.EQ.KWCOMM) GO TO 37
        IF(KWORD.EQ.KWEABS) GO TO 36
      ENDIF
C
      WRITE(26,1350)
 1350 FORMAT(/9X,'Maximum allowed step lengths of',
     1  ' electrons and positrons',/9X,'and local absorption ',
     2  'energies (non-void bodies).',
     3  //3X,'Body',3X,'DSMAX(IB)',4X,'EABSB(1,IB)',3X,
     3  'EABSB(2,IB)',3X,'EABSB(3,IB)',/4X,'IB',7X,'(cm)',10X,'(eV)',
     4  10X,'(eV)',10X,'(eV)')
      DO IB=1,NBODY
        IF(MATER(IB).GT.0) WRITE(26,1351) IB,DSMAX(IB),EABSB(1,IB),
     1    EABSB(2,IB),EABSB(3,IB)
 1351   FORMAT(3X,I4,1P,4E14.6)
      ENDDO
C
C  ************  Variance reduction.
C
      DO KB=1,NBV
        DO ICOLi=1,8
          DO KPARi=1,3
            FORCE(KB,KPARi,ICOLi)=1.0D0
          ENDDO
        ENDDO
        DO KPARi=1,3
          LFORCE(KB,KPARi)=.FALSE.
          WLOW(KB,KPARi)=0.0D0
          WHIG(KB,KPARi)=1.0D6
        ENDDO
        IBRSPL(KB)=1
      ENDDO
C
      DO I=1,NBV
        IXRSPL(I)=1
        LXRSPL(I)=.FALSE.
      ENDDO
C
C  ****  Interaction forcing.
C
      IF(KWORD.EQ.KWIFOR) THEN
        WRITE(26,1400)
 1400   FORMAT(/3X,72('-'),/
     1    3X,'>>>>>>  Interaction forcing: FORCE(IBODY,KPAR,ICOL)')
   41   CONTINUE
        CALL FWORD(BUFFER,BUF2,WORD,LENGTH)
        IF(LENGTH.EQ.4) THEN
          DO IB=1,NBODY
            IF(WORD.EQ.BALIAS(IB)) KB=IB
          ENDDO
          READ(BUF2,*) KPAR,ICOL,FORCER,WWLOW,WWHIG
        ELSE
          READ(BUFFER,*) KB,KPAR,ICOL,FORCER,WWLOW,WWHIG
        ENDIF
C  ****  Negative FORCER values are re-interpreted, as described in the
C        heading comments above.
        IF(FORCER.LT.-1.0D-6) THEN
          MM=MATER(KB)
          AVNCL=AVNCOL(E0,KPAR,MM,ICOL)
          IF(AVNCL.GT.1.0D-8) THEN
            FORCER=MAX(ABS(FORCER)/AVNCL,1.0D0)
          ELSE
            FORCER=MAX(ABS(FORCER),1.0D0)
          ENDIF
        ELSE
          FORCER=MAX(FORCER,1.0D0)
        ENDIF
        IF(WWLOW.LT.1.0D-6) WWLOW=1.0D-6
        IF(WWHIG.GT.1.0D6) WWHIG=1.0D6
        IF(KB.LT.1.OR.KB.GT.NBODY) THEN
          WRITE(26,'(A6,1X,A120)') KWORD,BUFFER
          WRITE(26,*) 'Incorrect KB value.'
          STOP 'Incorrect KB value.'
        ENDIF
        IF(KPAR.LT.1.OR.KPAR.GT.3) THEN
          WRITE(26,'(A6,1X,A120)') KWORD,BUFFER
          WRITE(26,*) 'Incorrect value of KPAR.'
          STOP 'Incorrect value of KPAR.'
        ENDIF
        IF(ICOL.LT.1.OR.ICOL.GT.8) THEN
          WRITE(26,'(A6,1X,A120)') KWORD,BUFFER
          WRITE(26,*) 'Incorrect value of ICOL.'
          STOP 'Incorrect value of ICOL.'
        ENDIF
        WLOW(KB,KPAR)=MAX(WLOW(KB,KPAR),WWLOW)*0.9999999999D0
        WHIG(KB,KPAR)=MIN(WHIG(KB,KPAR),WWHIG)*1.0000000001D0
        IF(WLOW(KB,KPAR).GT.WHIG(KB,KPAR)) THEN
          WRITE(26,'(A6,1X,A120)') KWORD,BUFFER
          WRITE(26,*) 'Incorrect weight window limits.'
          STOP 'Incorrect weight window limits.'
        ENDIF
        IF(FORCER.GT.1.0001D0) THEN
          LFORCE(KB,KPAR)=.TRUE.
          FORCE(KB,KPAR,ICOL)=FORCER
          WRITE(26,1410) KB,KPAR,ICOL,FORCER,WLOW(KB,KPAR),WHIG(KB,KPAR)
 1410     FORMAT(3X,'FORCE(',I4,',',I1,',',I1,') =',1P,E13.6,
     1      ',  weight window = (',E9.2,',',E9.2,')')
        ENDIF
   42   CONTINUE
        READ(5,'(A6,1X,A120)') KWORD,BUFFER
        IF(KWORD.EQ.KWCOMM) GO TO 42
        IF(KWORD.EQ.KWIFOR) GO TO 41
      ENDIF
C
C  ****  Bremsstrahlung splitting.
C
      IF(KWORD.EQ.KBRSPL) THEN
        WRITE(26,1420)
 1420   FORMAT(/3X,72('-'),/
     1    3X,'>>>>>>  Bremsstrahlung splitting')
   26   CONTINUE
        CALL FWORD(BUFFER,BUF2,WORD,LENGTH)
        IF(LENGTH.EQ.4) THEN
          DO IB=1,NBODY
            IF(WORD.EQ.BALIAS(IB)) KB=IB
          ENDDO
          READ(BUF2,*) IBR
        ELSE
          READ(BUFFER,*) KB,IBR
        ENDIF
        IF(KB.LT.1.OR.KB.GT.NBODY) THEN
          WRITE(26,'(A6,1X,A120)') KWORD,BUFFER
          WRITE(26,*) 'Incorrect KB value.'
          STOP 'Incorrect KB value.'
        ENDIF
        IF(IBR.LT.1) THEN
          WRITE(26,'(A6,1X,A120)') KWORD,BUFFER
          WRITE(26,*) 'Incorrect value of IBRSPL.'
          STOP 'Incorrect value of IBRSPL.'
        ENDIF
        IF(LFORCE(KB,KPAR)) THEN
          IBRSPL(KB)=IBR
        ELSE
          WRITE(26,'(A6,1X,A120)') KWORD,BUFFER
          WRITE(26,*) 'Interaction forcing unactive in this body.'
          STOP 'Interaction forcing unactive in this body.'
        ENDIF
        WRITE(26,1421) KB,IBRSPL(KB)
 1421   FORMAT(3X,'IBODY = ',I4,',  IBRSPL =',I4)
   27   CONTINUE
        READ(5,'(A6,1X,A120)') KWORD,BUFFER
        IF(KWORD.EQ.KWCOMM) GO TO 27
        IF(KWORD.EQ.KBRSPL) GO TO 26
      ENDIF
C
C  ****  X-ray splitting.
C
      IF(KWORD.EQ.KXRSPL) THEN
        WRITE(26,1422)
 1422   FORMAT(/3X,72('-'),/
     1    3X,'>>>>>>  X-ray splitting')
   57   CONTINUE
        CALL FWORD(BUFFER,BUF2,WORD,LENGTH)
        IF(LENGTH.EQ.4) THEN
          DO IB=1,NBODY
            IF(WORD.EQ.BALIAS(IB)) KB=IB
          ENDDO
          READ(BUF2,*) IXR
        ELSE
          READ(BUFFER,*) KB,IXR
        ENDIF
        IF(KB.LT.1.OR.KB.GT.NBODY) THEN
          WRITE(26,'(A6,1X,A120)') KWORD,BUFFER
          WRITE(26,*) 'Incorrect KB value.'
          STOP 'Incorrect KB value.'
        ENDIF
        IF(IXR.LT.1) THEN
          WRITE(26,'(A6,1X,A120)') KWORD,BUFFER
          WRITE(26,*) 'Incorrect value of IXRSPL.'
          STOP 'Incorrect value of IXRSPL.'
        ENDIF
        IXRSPL(KB)=IXR
        IF(IXR.GT.1) LXRSPL(KB)=.TRUE.
        IF(LXRSPL(KB)) WRITE(26,1423) KB,IXRSPL(KB)
 1423   FORMAT(3X,'IBODY = ',I4,',  IXRSPL =',I4)
   58   CONTINUE
        READ(5,'(A6,1X,A120)') KWORD,BUFFER
        IF(KWORD.EQ.KWCOMM) GO TO 58
        IF(KWORD.EQ.KXRSPL) GO TO 57
      ENDIF
C
C  ************  Energy and angular distributions of emerging
C                particles.
C
      WRITE(26,1500)
 1500 FORMAT(/3X,72('-'),/
     1  3X,'>>>>>>  Energy and angular distributions of emerging',
     1  ' particles.')
C
      IF(KWORD.EQ.KWNBE) THEN
        READ(BUFFER,*) EMIN,EMAX,NBE
   51   CONTINUE
        READ(5,'(A6,1X,A120)') KWORD,BUFFER
        IF(KWORD.EQ.KWCOMM) GO TO 51
      ELSE
        EMIN=0.0D0
        EMAX=EPMAX
        NBE=100
      ENDIF
      IF(EMIN.LT.1.0D0) EMIN=0.0D0
      IF(EMAX.LT.1.0D0) EMAX=EPMAX
 1510 FORMAT(3X,'E:       NBE = ',I3,
     1  ',  EMIN =',1P,E13.6,' eV,  EMAX =',E13.6,' eV')
      IF(NBE.LT.0) THEN
        EMIN=MAX(EMIN,1.0D0)
        WRITE(26,1510) NBE,EMIN,EMAX
        WRITE(26,1511)
 1511   FORMAT(12X,'(logarithmic scale, bin width increases with E)')
      ELSE IF(NBE.GT.0) THEN
        WRITE(26,1510) NBE,EMIN,EMAX
        WRITE(26,1512)
 1512   FORMAT(12X,'(linear scale, uniform bin width)')
      ENDIF
      IF (NBE.EQ.0) THEN
        WRITE(26,*) 'NBE is equal to zero.'
        STOP 'NBE equal to 0.'
      ENDIF
C
      IF(KWORD.EQ.KWNBAN) THEN
        READ(BUFFER,*) NBTH,NBPH
   52   CONTINUE
        READ(5,'(A6,1X,A120)') KWORD,BUFFER
        IF(KWORD.EQ.KWCOMM) GO TO 52
      ELSE
        NBTH=90
        NBPH=1
      ENDIF
      IF (NBTH.EQ.0) THEN
        WRITE(26,*) 'NBTH is equal to zero.'
        STOP 'NBTH equal to 0.'
      ENDIF
      IF(NBTH.GT.0) THEN
        WRITE(26,1520) NBTH
 1520   FORMAT(3X,'Theta:  NBTH = ',I3,' (linear scale)')
      ELSE
        WRITE(26,1521) NBTH
 1521   FORMAT(3X,'Theta:  NBTH = ',I3,' (logarithmic scale)')
      ENDIF
      WRITE(26,1530) NBPH
 1530 FORMAT(3X,'Phi:    NBPH = ',I3)
      IF(NBPH.LT.1) THEN
        WRITE(26,*) 'Wrong number of PHI bins.'
        STOP 'Wrong number of PHI bins.'
      ENDIF
C
      CALL ENANG0(EMIN,EMAX,NBE,NBTH,NBPH)  ! Initialisation.
C
C  ************  Impact detectors.
C
      DO KD=1,NIDM
        KKDI(KD,1)=0
        KKDI(KD,2)=0
        KKDI(KD,3)=0
      ENDDO
      NDBOD=0
      NPSFO=0
C
      NID=0
   61 CONTINUE
      IF(KWORD.EQ.KWIDET) THEN
        NID=NID+1
        NDBOD=0
        WRITE(26,1600) NID
 1600   FORMAT(/3X,72('-'),/
     1    3X,'>>>>>>  Impact detector #', I2)
        READ(BUFFER,*) EDIL,EDIU,NDICH,IPSF(NID),IDCUT(NID)
C
        IF(NDICH.EQ.0) THEN
          WRITE(26,'(A6,1X,A120)') KWORD,BUFFER
          WRITE(26,*) 'Incorrect number of energy bins.'
          STOP 'Incorrect number of energy bins.'
        ENDIF
C
        IF(EDIL.LT.50.0D0) EDIL=50.0D0
        IF(EDIU.LT.50.1D0) EDIU=EPMAX
C
        IF(NDICH.LT.0) THEN
          EDIL=MAX(EDIL,50.0D0)
          WRITE(26,1614) EDIL,EDIU,ABS(NDICH)
 1614     FORMAT(3X,'Energy window = (',1P,E12.5,',',E12.5,') eV',
     1      /3X,'Number of energy bins = ',I4,'  (logarithmic scale)')
        ELSE
          WRITE(26,1610) EDIL,EDIU,NDICH
 1610     FORMAT(3X,'Energy window = (',1P,E12.5,',',E12.5,') eV',
     1      /3X,'Number of energy bins = ',I4,'  (linear scale)')
        ENDIF
C
        IF(EDIU.LT.EDIL+1.0D0) THEN
          WRITE(26,'(A6,1X,A120)') KWORD,BUFFER
          WRITE(26,*) 'Incorrect energy limits.'
          STOP 'Incorrect energy limits.'
        ENDIF
C
        IF(IPSF(NID).LT.0.OR.IPSF(NID).GT.1) THEN
          WRITE(26,'(A6,1X,A120)') KWORD,BUFFER
          WRITE(26,*) 'Wrong IPSF value.'
          STOP 'Wrong IPSF value.'
        ENDIF
C
        IF(IDCUT(NID).LT.0.OR.IDCUT(NID).GT.2) THEN
          WRITE(26,'(A6,1X,A120)') KWORD,BUFFER
          WRITE(26,*) 'Wrong IDCUT value.'
          STOP 'Wrong IDCUT value.'
        ENDIF
C
   62   CONTINUE
        READ(5,'(A6,1X,A120)') KWORD,BUFFER
        IF(KWORD.EQ.KWCOMM) GO TO 62
C
        IF(KWORD.EQ.KWISPC) THEN
          IF(NID.EQ.0) THEN
            WRITE(26,'(A6,1X,A120)') KWORD,BUFFER
            WRITE(26,*) 'No impact detector has been defined yet.'
            STOP 'No impact detector has been defined yet.'
          ENDIF
          READ(BUFFER,'(A20)') SPCDIO
          WRITE(26,1612) SPCDIO
 1612     FORMAT(3X,'Output energy spectrum: ',A20)
   64     CONTINUE
          READ(5,'(A6,1X,A120)') KWORD,BUFFER
          IF(KWORD.EQ.KWCOMM) GO TO 64
        ELSE
          WRITE(BUF2,'(I5)') 1000+NID
          SPCDIO='spc-impdet-'//BUF2(4:5)//'.dat'
          WRITE(26,1612) SPCDIO
        ENDIF
C
        IF(KWORD.EQ.KWIPSF) THEN
          IF(NID.EQ.0) THEN
            WRITE(26,'(A6,1X,A120)') KWORD,BUFFER
            WRITE(26,*) 'No impact detector has been defined yet.'
            STOP 'No impact detector has been defined yet.'
          ENDIF
          READ(BUFFER,'(A20)') PSFDIO
          IF(IPSF(NID).GT.0) THEN
            WRITE(26,1611) PSFDIO
 1611       FORMAT(3X,'Output phase-space file: ',A20)
            NPSFO=NPSFO+1
            IF(NPSFO.GT.1) THEN
              WRITE(26,*) 'You cannot generate more than one PSF in ',
     1          'a single run.'
              STOP 'Only one PSF can be generated in a each run.'
            ENDIF
          ELSE
            WRITE(26,1613)
 1613       FORMAT(3X,'No phase-space file is generated.')
          ENDIF
   63     CONTINUE
          READ(5,'(A6,1X,A120)') KWORD,BUFFER
          IF(KWORD.EQ.KWCOMM) GO TO 63
        ELSE
          IF(IPSF(NID).GT.0) THEN
            WRITE(BUF2,'(I5)') 1000+NID
            PSFDIO='psf-impdet-'//BUF2(4:5)//'.dat'
            WRITE(26,1611) PSFDIO
            NPSFO=NPSFO+1
            IF(NPSFO.GT.1) THEN
              WRITE(26,*) 'You cannot generate more than one PSF in ',
     1          'a single run.'
              STOP 'Only one PSF can be generated in each run.'
            ENDIF
          ENDIF
        ENDIF
C
        IF(ABS(IDCUT(NID)).EQ.0) THEN
          WRITE(26,1601)
 1601     FORMAT(3X,'Detected particles are absorbed')
        ELSE
          WRITE(26,1602)
 1602     FORMAT(3X,'Particles are transported through this detector')
        ENDIF
C
        SPCFSO='none'
        IF(KWORD.EQ.KWIFLN) THEN
          IF(NID.EQ.0) THEN
            WRITE(26,'(A6,1X,A120)') KWORD,BUFFER
            WRITE(26,*) 'No impact detector has been defined yet.'
            STOP 'No impact detector has been defined yet.'
          ENDIF
          READ(BUFFER,'(A20)') SPCFSO
   83     CONTINUE
          READ(5,'(A6,1X,A120)') KWORD,BUFFER
          IF(KWORD.EQ.KWCOMM) GO TO 83
          IF(IDCUT(NID).EQ.2) THEN
            WRITE(26,1662) SPCFSO
 1662       FORMAT(3X,'Output fluence distribution: ',A20)
          ENDIF
        ELSE
          IF(IDCUT(NID).EQ.2) THEN
            WRITE(BUF2,'(I5)') 1000+NID
            SPCFSO='fln-impdet-'//BUF2(4:5)//'.dat'
            WRITE(26,1662) SPCFSO
          ENDIF
        ENDIF
C
        AGEU=-1.0D6
        AGEL=0.0D0
        SPCAGE='none'
        NAGE=0
        IF(KWORD.EQ.KWDIAL) THEN
          IF(NID.EQ.0) THEN
            WRITE(26,'(A6,1X,A120)') KWORD,BUFFER
            WRITE(26,*) 'No impact detector has been defined yet.'
            STOP 'No impact detector has been defined yet.'
          ENDIF
          READ(BUFFER,*) AGEL,AGEU,NAGE
          IF(NAGE.EQ.0) THEN
            WRITE(26,'(A6,1X,A120)') KWORD,BUFFER
            WRITE(26,*) 'Incorrect number of age bins.'
            STOP 'Incorrect number of age bins.'
          ENDIF
C
          AGEL=MAX(AGEL,0.0D0)
          IF(NAGE.LT.0) THEN
            AGEL=MAX(AGEL,1.0D-20)
            WRITE(26,1671) AGEL,AGEU,ABS(NAGE)
 1671       FORMAT(3X,'Particle age window = (',1P,E12.5,',',E12.5,
     1        ') seconds',/3X,'Number of age bins = ',I4,
     1        '  (logarithmic scale)')
          ELSE
            WRITE(26,1672) AGEL,AGEU,NAGE
 1672       FORMAT(3X,'Particle age window = (',1P,E12.5,',',E12.5,
     1        ') seconds',/3X,'Number of age bins = ',I4,
     1        '  (linear scale)')
          ENDIF
          IF(AGEU.LT.AGEL+1.0D-19) THEN
            WRITE(26,'(A6,1X,A120)') KWORD,BUFFER
            WRITE(26,*) 'Incorrect age limits.'
            STOP 'Incorrect age limits.'
          ENDIF
C
   84     CONTINUE
          READ(5,'(A6,1X,A120)') KWORD,BUFFER
          IF(KWORD.EQ.KWCOMM) GO TO 84
        ENDIF
C
        IF(KWORD.EQ.KWDIAF) THEN
          IF(NID.EQ.0) THEN
            WRITE(26,'(A6,1X,A120)') KWORD,BUFFER
            WRITE(26,*) 'No impact detector has been defined yet.'
            STOP 'No impact detector has been defined yet.'
          ENDIF
          IF(AGEL.LT.0.0D0) THEN
            WRITE(26,'(A6,1X,A120)') KWORD,BUFFER
            WRITE(26,*) 'Undefined age distribution limits.'
            STOP 'Undefined age distribution limits.'
          ENDIF
          READ(BUFFER,'(A20)') SPCAGE
          WRITE(26,1673) SPCAGE
 1673     FORMAT(3X,'Output age distribution: ',A20)
   85     CONTINUE
          READ(5,'(A6,1X,A120)') KWORD,BUFFER
          IF(KWORD.EQ.KWCOMM) GO TO 85
        ELSE
          IF(AGEU.GT.0.0D0) THEN
            WRITE(BUF2,'(I5)') 1000+NID
            SPCAGE='age-impdet-'//BUF2(4:5)//'.dat'
            WRITE(26,1673) SPCAGE
          ENDIF
        ENDIF
C
   65   CONTINUE
        IF(KWORD.EQ.KWIBOD) THEN
          IF(NID.EQ.0) THEN
            WRITE(26,'(A6,1X,A120)') KWORD,BUFFER
            WRITE(26,*) 'No impact detector has been defined yet.'
            STOP 'No impact detector has been defined yet.'
          ENDIF
          CALL FWORD(BUFFER,BUF2,WORD,LENGTH)
          IF(LENGTH.EQ.4) THEN
            DO IB=1,NBODY
              IF(WORD.EQ.BALIAS(IB)) KB=IB
            ENDDO
          ELSE
            READ(BUFFER,*) KB
          ENDIF
          IF(KB.LT.1.OR.KB.GT.NBODY) THEN
            WRITE(26,'(A6,1X,A120)') KWORD,BUFFER
            WRITE(26,*) 'Incorrect body label.'
            STOP 'Incorrect body label.'
          ENDIF
          IF(KDET(KB).NE.0) THEN
            WRITE(26,'(A6,1X,A120)') KWORD,BUFFER
            WRITE(26,*) 'A body cannot be part of two detectors.'
            STOP 'A body cannot be part of two detectors.'
          ENDIF
          IF(MATER(KB).EQ.0) THEN
            WRITE(26,'(A6,1X,A120)') KWORD,BUFFER
            WRITE(26,*) 'A void body cannot be part of a detector.'
            STOP 'A void body cannot be part of a detectors.'
          ENDIF
          WRITE(26,1620) KB
 1620     FORMAT(3X,'Active body = ',I4)
          KDET(KB)=NID
          NDBOD=NDBOD+1
   66     CONTINUE
          READ(5,'(A6,1X,A120)') KWORD,BUFFER
          IF(KWORD.EQ.KWCOMM) GO TO 66
          IF(KWORD.EQ.KWIBOD) GO TO 65
          IF(KWORD.EQ.KWIDET) THEN
            IF(NDBOD.EQ.0) THEN
              WRITE(26,*) 'This detector has no active bodies.'
              STOP 'This detector has no active bodies.'
            ENDIF
            ITST=MAX(KKDI(NID,1),KKDI(NID,2),KKDI(NID,3))
            IF(ITST.EQ.0) THEN
              KKDI(NID,1)=1
              KKDI(NID,2)=1
              KKDI(NID,3)=1
              WRITE(26,1630)
 1630         FORMAT(3X,'Detected particles = electrons, photons and ',
     1          'positrons')
            ENDIF
            CALL IMDET0(EDIL,EDIU,NDICH,AGEL,AGEU,NAGE,IDCUT(NID),
     1        SPCDIO,SPCFSO,SPCAGE,NID)
            GO TO 61
          ENDIF
        ENDIF
C
   67   CONTINUE
        IF(KWORD.EQ.KWIPAR) THEN
          IF(NID.EQ.0) THEN
            WRITE(26,'(A6,1X,A120)') KWORD,BUFFER
            WRITE(26,*) 'No impact detector has been defined yet.'
            STOP 'No impact detector has been defined yet.'
          ENDIF
          IF(NDBOD.EQ.0) THEN
            WRITE(26,*) 'This detector has no active bodies.'
            STOP 'This detector has no active bodies.'
          ENDIF
C
          READ(BUFFER,*) KPARD
          IF(KPARD.EQ.1) THEN
            KKDI(NID,1)=1
            WRITE(26,1631)
 1631       FORMAT(3X,'Detected particles = electrons')
          ELSE IF(KPARD.EQ.2) THEN
            KKDI(NID,2)=1
            WRITE(26,1632)
 1632       FORMAT(3X,'Detected particles = photons')
          ELSE IF(KPARD.EQ.3) THEN
            KKDI(NID,3)=1
            WRITE(26,1633)
 1633       FORMAT(3X,'Detected particles = positrons')
          ENDIF
   68     CONTINUE
          READ(5,'(A6,1X,A120)') KWORD,BUFFER
          IF(KWORD.EQ.KWCOMM) GO TO 68
          IF(KWORD.EQ.KWIPAR) GO TO 67
          IF(KWORD.EQ.KWIBOD) GO TO 65
          IF(KWORD.EQ.KWIDET) THEN
            IF(NDBOD.EQ.0) THEN
              WRITE(26,*) 'This detector has no active bodies.'
              STOP 'This detector has no active bodies.'
            ENDIF
            ITST=MAX(KKDI(NID,1),KKDI(NID,2),KKDI(NID,3))
            IF(ITST.EQ.0) THEN
              KKDI(NID,1)=1
              KKDI(NID,2)=1
              KKDI(NID,3)=1
              WRITE(26,1630)
            ENDIF
            CALL IMDET0(EDIL,EDIU,NDICH,AGEL,AGEU,NAGE,IDCUT(NID),
     1        SPCDIO,SPCFSO,SPCAGE,NID)
            GO TO 61
          ENDIF
        ENDIF
      ENDIF
C
      IF(NID.GT.0) THEN
        IF(NDBOD.EQ.0) THEN
          WRITE(26,*) 'This detector has no active bodies.'
          STOP 'This detector has no active bodies.'
        ENDIF
        ITST=MAX(KKDI(NID,1),KKDI(NID,2),KKDI(NID,3))
        IF(ITST.EQ.0) THEN
          KKDI(NID,1)=1
          KKDI(NID,2)=1
          KKDI(NID,3)=1
          WRITE(26,1630)
        ENDIF
        CALL IMDET0(EDIL,EDIU,NDICH,AGEL,AGEU,NAGE,IDCUT(NID),
     1    SPCDIO,SPCFSO,SPCAGE,NID)
      ENDIF
C
C  ************  Energy-deposition detectors.
C
      DO KB=1,NBODY
        KBDE(KB)=0
      ENDDO
      NDBOD=0
C
      NED=0
   43 CONTINUE
      IF(KWORD.EQ.KWEDET) THEN
        IF(NED.GT.0) THEN
          IF(NDBOD.EQ.0) THEN
            WRITE(26,*) 'This detector has no active bodies.'
            STOP 'This detector has no active bodies.'
          ENDIF
        ENDIF
        NED=NED+1
        NDBOD=0
        WRITE(26,1650) NED
 1650   FORMAT(/3X,72('-'),/
     1    3X,'>>>>>>  Energy-deposition detector #', I2)
        READ(BUFFER,*) EDEL,EDEU,NDECH
C
        IF(NDECH.EQ.0) THEN
          WRITE(26,'(A6,1X,A120)') KWORD,BUFFER
          WRITE(26,*) 'Incorrect number of energy bins.'
          STOP 'Incorrect number of energy bins.'
        ENDIF
C
        IF(EDEL.LT.1.0D0) EDEL=0.0D0
        IF(EDEU.LT.1.0D0) EDEU=EPMAX
C
        IF(NDECH.LT.0) THEN
          EDEL=MAX(EDEL,1.0D0)
          WRITE(26,1684) EDEL,EDEU,ABS(NDECH)
 1684     FORMAT(3X,'Energy window = (',1P,E12.5,',',E12.5,') eV',
     1      /3X,'Number of energy bins = ',I4,'  (logarithmic scale)')
        ELSE
          WRITE(26,1685) EDEL,EDEU,NDECH
 1685     FORMAT(3X,'Energy window = (',1P,E12.5,',',E12.5,') eV',
     1      /3X,'Number of energy bins = ',I4,'  (linear scale)')
        ENDIF
C
        IF(EDEU.LT.EDEL+1.0D0) THEN
          WRITE(26,'(A6,1X,A120)') KWORD,BUFFER
          WRITE(26,*) 'Incorrect energy limits.'
          STOP 'Incorrect energy limits.'
        ENDIF
C
   44   CONTINUE
        READ(5,'(A6,1X,A120)') KWORD,BUFFER
        IF(KWORD.EQ.KWCOMM) GO TO 44
        IF(KWORD.EQ.KWESPC) THEN
          READ(BUFFER,'(A20)') SPCDEO
          WRITE(26,1651) SPCDEO
 1651     FORMAT(3X,'Output spectrum: ',A20)
   45     CONTINUE
          READ(5,'(A6,1X,A120)') KWORD,BUFFER
          IF(KWORD.EQ.KWCOMM) GO TO 45
        ELSE
          WRITE(BUF2,'(I5)') 1000+NED
          SPCDEO='spc-enddet-'//BUF2(4:5)//'.dat'
          WRITE(26,1651) SPCDEO
        ENDIF
C
   46   CONTINUE
        IF(KWORD.EQ.KWEBOD) THEN
          CALL FWORD(BUFFER,BUF2,WORD,LENGTH)
          IF(LENGTH.EQ.4) THEN
            DO IB=1,NBODY
              IF(WORD.EQ.BALIAS(IB)) KB=IB
            ENDDO
          ELSE
            READ(BUFFER,*) KB
          ENDIF
          IF(KB.LT.1.OR.KB.GT.NBODY) THEN
            WRITE(26,'(A6,1X,A120)') KWORD,BUFFER
            WRITE(26,*) 'Incorrect body label.'
            STOP 'Incorrect body label.'
          ENDIF
          IF(KBDE(KB).NE.0) THEN
            WRITE(26,'(A6,1X,A120)') KWORD,BUFFER
            WRITE(26,*) 'A body cannot be part of two detectors.'
            STOP 'A body cannot be part of two detectors.'
          ENDIF
          IF(MATER(KB).EQ.0) THEN
            WRITE(26,'(A6,1X,A120)') KWORD,BUFFER
            WRITE(26,*) 'A void body cannot be part of a detector.'
            STOP 'A void body cannot be part of a detectors.'
          ENDIF
          WRITE(26,1652) KB
 1652     FORMAT(3X,'Active body = ',I4)
          IF(LFORCE(KB,1).OR.LFORCE(KB,2).OR.LFORCE(KB,3)) THEN
            WRITE(26,'(3X,''#  WARNING: Spectrum may be strongly'',
     1        '' biased'',/15X,''when interaction forcing is used!'')')
          ENDIF
          KBDE(KB)=NED
          NDBOD=NDBOD+1
        ENDIF
   47   CONTINUE
        READ(5,'(A6,1X,A120)') KWORD,BUFFER
        IF(KWORD.EQ.KWCOMM) GO TO 47
        IF(KWORD.EQ.KWEBOD) GO TO 46
        IF(KWORD.EQ.KWEDET) THEN
          CALL ENDET0(EDEL,EDEU,NDECH,SPCDEO,NED)
          GO TO 43
        ENDIF
      ENDIF
C
      IF(NED.GT.0) THEN
        IF(NDBOD.EQ.0) THEN
          WRITE(26,*) 'This detector has no active bodies.'
          STOP 'This detector has no active bodies.'
        ENDIF
        CALL ENDET0(EDEL,EDEU,NDECH,SPCDEO,NED)
      ENDIF
C
C  ************  Dose distribution.
C
      LDOSEM=.FALSE.
      NBX=0
      NBY=0
      NBZ=0
C
      IF(KWORD.EQ.KGRDXX) THEN
        LDOSEM=.TRUE.
        IDOSE=1
        WRITE(26,1700)
 1700   FORMAT(/3X,72('-'),/3X,'>>>>>>  Dose distribution in a box.')
        READ(BUFFER,*) XLD,XUD,NBX
        IF(XLD.GT.XUD) THEN
          SAVE=XLD
          XLD=XUD
          XUD=SAVE
        ENDIF
        IF(XUD.LT.XLD+1.0D-6) THEN
          WRITE(26,'(A6,1X,A120)') KWORD,BUFFER
          WRITE(26,*) 'XU must be greater than XL+1.0E-6.'
          STOP 'XU must be greater than XL+1.0E-6.'
        ENDIF
        NBX=MAX(1,NBX)
        WRITE(26,1710) XLD,XUD,NBX
 1710   FORMAT(3X,'XL = ',1P,E13.6,' cm,  XU = ',E13.6,
     1    ' cm,  NDBX =',I4)
   70   CONTINUE
        READ(5,'(A6,1X,A120)') KWORD,BUFFER
        IF(KWORD.EQ.KWCOMM) GO TO 70
      ENDIF
C
      IF(KWORD.EQ.KGRDYY) THEN
        READ(BUFFER,*) YLD,YUD,NBY
        IF(YLD.GT.YUD) THEN
          SAVE=YLD
          YLD=YUD
          YUD=SAVE
        ENDIF
        IF(NBX.LT.1) THEN
          WRITE(26,'(A6,1X,A120)') KWORD,BUFFER
          WRITE(26,*) 'Incorrect keyword.'
          STOP 'Incorect keyword.'
        ENDIF
        IF(YUD.LT.YLD+1.0D-6) THEN
          WRITE(26,'(A6,1X,A120)') KWORD,BUFFER
          WRITE(26,*) 'YU must be greater than YL+1.0E-6.'
          STOP 'YU must be greater than YL+1.0E-6.'
        ENDIF
        NBY=MAX(1,NBY)
        WRITE(26,1711) YLD,YUD,NBY
 1711   FORMAT(3X,'YL = ',1P,E13.6,' cm,  YU = ',E13.6,
     1    ' cm,  NDBY =',I4)
   71   CONTINUE
        READ(5,'(A6,1X,A120)') KWORD,BUFFER
        IF(KWORD.EQ.KWCOMM) GO TO 71
      ELSE
        IF(NBX.GT.0) THEN
          WRITE(26,'(A6,1X,A120)') KWORD,BUFFER
          WRITE(26,*) 'Incorrect keyword.'
          STOP 'Incorect keyword.'
        ENDIF
      ENDIF
C
      IF(KWORD.EQ.KGRDZZ) THEN
        READ(BUFFER,*) ZLD,ZUD,NBZ
        IF(ZLD.GT.ZUD) THEN
          SAVE=ZLD
          ZLD=ZUD
          ZUD=SAVE
        ENDIF
        IF(ZUD.LT.ZLD+1.0D-6) THEN
          WRITE(26,'(A6,1X,A120)') KWORD,BUFFER
          WRITE(26,*) 'ZU must be greater than ZL+1.0E-6.'
          STOP 'ZU must be greater than ZL+1.0E-6.'
        ENDIF
        NBZ=MAX(1,NBZ)
        IF(NBX.GT.0) THEN
          WRITE(26,1712) ZLD,ZUD,NBZ
 1712     FORMAT(3X,'ZL = ',1P,E13.6,' cm,  ZU = ',E13.6,
     1      ' cm,  NDBZ =',I4)
        ELSE
          LDOSEM=.TRUE.
          IDOSE=2
          WRITE(26,1701)
 1701     FORMAT(/3X,72('-'),/3X,'>>>>>>  Dose distribution in a',
     1      ' cylinder.')
          WRITE(26,1712) ZLD,ZUD,NBZ
        ENDIF
   72   CONTINUE
        READ(5,'(A6,1X,A120)') KWORD,BUFFER
        IF(KWORD.EQ.KWCOMM) GO TO 72
      ELSE
        IF(NBX.GT.0) THEN
          WRITE(26,'(A6,1X,A120)') KWORD,BUFFER
          WRITE(26,*) 'Unrecognized keyword.'
          STOP 'Unrecognized keyword.'
        ENDIF
      ENDIF
C
      IF(KWORD.EQ.KGRDRR) THEN
        READ(BUFFER,*) XUD,NBR
        IF(XUD.LT.1.0D-6) THEN
          WRITE(26,'(A6,1X,A120)') KWORD,BUFFER
          WRITE(26,*) 'RU must be greater than 1.0E-6.'
          STOP 'RU must be greater than 1.0E-6.'
        ENDIF
        IF(NBX.GT.0.OR.NBY.GT.0) THEN
          WRITE(26,'(A6,1X,A120)') KWORD,BUFFER
          WRITE(26,*) 'Incorrect keyword.'
          STOP 'Incorect keyword.'
        ENDIF
        NBX=MAX(1,NBR)
        XLD=0.0D0
        IF(NBZ.GT.0) THEN
          WRITE(26,1713) XUD,NBX
 1713     FORMAT(3X,'RU = ',1P,E13.6,' cm,  NDBR =',I4)
          NBY=1
        ELSE
          LDOSEM=.TRUE.
          IDOSE=3
          WRITE(26,1702)
 1702     FORMAT(/3X,72('-'),/3X,'>>>>>>  Dose distribution in a',
     1      ' sphere.')
          WRITE(26,1713) XUD,NBX
          NBY=1
          NBZ=1
        ENDIF
   73   CONTINUE
        READ(5,'(A6,1X,A120)') KWORD,BUFFER
        IF(KWORD.EQ.KWCOMM) GO TO 73
      ENDIF
C
      IF(LDOSEM) THEN
        CALL DOSE0(XLD,XUD,YLD,YUD,ZLD,ZUD,NBX,NBY,NBZ,IDOSE)
      ENDIF
C
C  ************  Job characteristics.
C
      WRITE(26,1800)
 1800 FORMAT(/3X,72('-'),/
     1  3X,'>>>>>>  Job characteristics.')
C
      IRESUM=0
      IF(KWORD.EQ.KWRESU) THEN
        READ(BUFFER,'(A20)') PFILER
        WRITE(26,1810) PFILER
 1810   FORMAT(3X,'Resume simulation from previous dump file: ',A20)
        IRESUM=1
   75   CONTINUE
        READ(5,'(A6,1X,A120)') KWORD,BUFFER
        IF(KWORD.EQ.KWCOMM) GO TO 75
      ENDIF
C
      DUMPP=1.0D15
      IF(KWORD.EQ.KWDUMP) THEN
        READ(BUFFER,'(A20)') PFILED
        WRITE(26,1820) PFILED
 1820   FORMAT(3X,'Write final counter values on the dump file: ',A20)
        LDUMP=.TRUE.
   76   CONTINUE
        READ(5,'(A6,1X,A120)') KWORD,BUFFER
        IF(KWORD.EQ.KWCOMM) GO TO 76
      ENDIF
C
      IF(KWORD.EQ.KWDMPP) THEN
        READ(BUFFER,*) DUMPP
        IF(LDUMP) THEN
          IF(DUMPP.LT.15.0D0) DUMPP=15.0D0
          IF(DUMPP.GT.86400.0D0) DUMPP=86400.0D0
          WRITE(26,1830) DUMPP
 1830     FORMAT(3X,'Dumping period: DUMPP =',1P,E13.6)
        ENDIF
   77   CONTINUE
        READ(5,'(A6,1X,A120)') KWORD,BUFFER
        IF(KWORD.EQ.KWCOMM) GO TO 77
      ENDIF
C
      IF(KWORD.EQ.KWRSEE) THEN
        READ(BUFFER,*) ISEED1,ISEED2
        IF(ISEED1.LT.0) CALL RAND0(-ISEED1)
   79   CONTINUE
        READ(5,'(A6,1X,A120)') KWORD,BUFFER
        IF(KWORD.EQ.KWCOMM) GO TO 79
      ELSE
        ISEED1=1
        ISEED2=1
      ENDIF
      WRITE(26,1850) ISEED1,ISEED2
 1850 FORMAT(3X,'Random-number generator seeds = ',I10,', ',I10)
C
      IF(KWORD.EQ.KWNSIM) THEN
        READ(BUFFER,*) DSHN
        IF(DSHN.LT.1.0D0) DSHN=2.0D9
   78   CONTINUE
        READ(5,'(A6,1X,A120)') KWORD,BUFFER
        IF(KWORD.EQ.KWCOMM) GO TO 78
      ELSE
        DSHN=2.0D9
      ENDIF
      WRITE(26,1840) DSHN
 1840 FORMAT(/3X,'Number of showers to be simulated =',1P,E13.6)
C
      IF(KWORD.EQ.KWTIME) THEN
        READ(BUFFER,*) TIMEA
      ELSE
        TIMEA=2.0D9
      ENDIF
      IF(TIMEA.LT.1.0D0) TIMEA=100.0D0
      WRITE(26,1860) TIMEA
 1860 FORMAT(3X,'Computation time available = ',1P,E13.6,' sec')
C
      CALL TIMER(TSECIN)
      TSECA=TIMEA+TSECIN
      TSECAD=TSECIN
      WRITE(26,1870)
 1870 FORMAT(/3X,72('-'))
C
C  ************  If 'RESUME' is active, read previously generated
C                counters...
C
      SHNA=0.0D0
      CPUTA=0.0D0
      N=0
C
      RLAST=0.0D0
      RWRITE=0.0D0
C
      IF(IRESUM.EQ.1) THEN
        OPEN(9,FILE=PFILER)
        READ (9,*,ERR=91,END=91) SHNA,CPUTA
        WRITE(6,*) '  Reading the DUMP file ...'
        READ (9,'(A65)',ERR=90) TITLE2
        IF(TITLE2.NE.TITLE) THEN
          WRITE(26,*)
     1      'The dump file is corrupted (the TITLE does not match).'
          STOP 'The dump file is corrupted (the TITLE does not match).'
        ENDIF
        READ (9,*,ERR=90) ISEED1,ISEED2
        READ (9,*,ERR=90) NPSN,RLREAD
        READ (9,*,ERR=90) KPARP
        IF(KPARP.EQ.0) THEN
          READ (9,*,ERR=90) NSDE,DSDE,RDSDE
          READ (9,*,ERR=90) ((SEDS(K,I),I=1,NSDE),K=1,3)
          READ (9,*,ERR=90) ((SEDS2(K,I),I=1,NSDE),K=1,3)
        ELSE
          READ (9,*,ERR=90) NSEB
          READ (9,*,ERR=90) (ESRC(I),I=1,NSEB+1)
          READ (9,*,ERR=90) (PSRC(I),I=1,NSEB+1)
          READ (9,*,ERR=90) (SHIST(I),I=1,NSEB+1)
        ENDIF
        READ (9,*,ERR=90) (PRIM(I),I=1,3)
        READ (9,*,ERR=90) (PRIM2(I),I=1,3)
        READ (9,*,ERR=90) ((SEC(K,I),I=1,3),K=1,3)
        READ (9,*,ERR=90) ((SEC2(K,I),I=1,3),K=1,3)
        READ (9,*,ERR=90) (AVW(I),I=1,2)
        READ (9,*,ERR=90) (AVW2(I),I=1,2)
        READ (9,*,ERR=90) (AVA(I),I=1,2)
        READ (9,*,ERR=90) (AVA2(I),I=1,2)
        READ (9,*,ERR=90) (AVE(I),I=1,2)
        READ (9,*,ERR=90) (AVE2(I),I=1,2)
        READ (9,*,ERR=90) NBODY
        READ (9,*,ERR=90) (TDEBO(I),I=1,NBODY)
        READ (9,*,ERR=90) (TDEBO2(I),I=1,NBODY)
        CALL ENANGR(9)  ! Energy and angular distributions.
        READ (9,*,ERR=90) NID,NED,LDOSEM
        IF(NID.GT.0) THEN
          READ (9,*,ERR=90) RLAST,RWRITE
          CALL IMDETR(9)  ! Impact detectors.
        ENDIF
        IF(NED.GT.0) THEN
          CALL ENDETR(9)  ! Energy-deposition detectors.
        ENDIF
        IF(LDOSEM) THEN
          CALL DOSER(9)  ! Dose distribution.
        ENDIF
        CLOSE(9)
        WRITE(26,1880) PFILER
 1880   FORMAT(3X,'Simulation has been resumed from dump file: ',A20)
        GO TO 92
   90   CONTINUE
        WRITE(26,*) 'The dump file is empty or corrupted.'
        STOP 'The dump file is empty or corrupted.'
   91   CONTINUE
        WRITE(26,1890)
 1890   FORMAT(3X,'WARNING: Could not resume from dump file...')
        CLOSE(9)
        IRESUM=0
      ENDIF
   92 CONTINUE
C
      IPSFO=21
      IF(NID.GT.0) THEN
        DO ID=1,NID
          IF(IPSF(ID).GT.0) THEN
            OPEN(IPSFO,FILE=PSFDIO,IOSTAT=KODE)
            IF(KODE.NE.0) THEN
              WRITE(26,'(''File '',A20,'' could not be opened.'')')
     1          PSFDIO
              STOP 'File could not be opened.'
            ENDIF
            RWR=0.0D0
            IF(RWRITE.GT.0.0D0) THEN
   93         CONTINUE
              KODEPS=0
              CALL RDPSF(IPSFO,NSHI,ISEC,KODEPS)
              IF(KODEPS.EQ.-1) THEN
                GO TO 94
              ELSE
                RWR=RWR+1.0D0
                IF(RWR.LT.RWRITE-0.5D0) GO TO 93
                GO TO 94
              ENDIF
            ENDIF
   94       CONTINUE
            IF(RWR.LT.0.5D0) THEN
              WRITE(IPSFO,1901) ID
 1901         FORMAT(1X,'#  Results from PENMAIN. Phase-space fi',
     1          'le of detector no.',I3,/1X,'#')
              WRITE(IPSFO,1902)
 1902         FORMAT(1X,'#/KPAR',2X,'E',12X,'X',12X,'Y',12X,'Z',12X,
     1          'U',12X,'V',12X,'W',11X,'WGHT',5X,'ILB(1:4)',7X,'NSHI',
     1          /1X,'#',125('-'))
            ENDIF
          ENDIF
        ENDDO
      ENDIF
C
      IPSFI=20
      IF(LPSF) THEN
        IF(IRESUM.EQ.1) THEN
          IF(NPSN.GT.NPSF) THEN
            WRITE(6,*) '  **** The simulation was already completed.'
            WRITE(26,*) '  **** The simulation was already completed.'
            SHN=SHNA
            TSIM=CPUTA
            CPUT0=CPUTIM()
            JOBEND=1
            RETURN
          ENDIF
          IF(NPSN.GT.1) THEN
            DO I=1,NPSN-1
              WRITE(6,1903) PSFI(I)
              WRITE(26,1903) PSFI(I)
 1903         FORMAT(/3X,'+ The phase-space file ',A20,
     1          ' was read in previous runs.')
            ENDDO
          ENDIF
   95     CONTINUE
          OPEN(IPSFI,FILE=PSFI(NPSN))
          KODEPS=0
          CALL RDPSF(IPSFI,NSHI,ISEC,KODEPS)
          WRITE(6,1904) PSFI(NPSN)
          WRITE(26,1904) PSFI(NPSN)
 1904     FORMAT(/3X,'+ The phase-space file ',A20,' is opened.')
          IF(KODEPS.EQ.-1) THEN
            WRITE(26,*) ' STOP. The file is empty or corrupted.'
            STOP 'The file is empty or corrupted.'
          ENDIF
          IF(RLREAD.GT.0.5D0) THEN
            RI=0.0D0
   96       CONTINUE
            RI=RI+1.0D0
            CALL RDPSF(IPSFI,NSHI,ISEC,KODEPS)
            IF(KODEPS.EQ.-1) THEN
              NPSN=NPSN+1
              IF(NPSN.GT.NPSF) THEN
                WRITE(6,*)
     1            '  **** The simulation was already completed.'
                WRITE(26,*)
     1            '  **** The simulation was already completed.'
                SHN=SHNA
                TSIM=CPUTA
                CPUT0=CPUTIM()
                JOBEND=2
                RETURN
              ELSE
                WRITE(6,1905) PSFI(NPSN)
                WRITE(26,1905) PSFI(NPSN)
 1905           FORMAT(/3X,'+ The phase-space file ',A20,
     1            ' was completed in the last run.')
                CLOSE(IPSFI)
                RLREAD=0.0D0
                GO TO 95
              ENDIF
            ENDIF
            IF(RI.LT.RLREAD-0.5D0) GO TO 96
          ELSE
            RLREAD=0.0D0
          ENDIF
        ELSE
          NPSN=1
          RLREAD=0.0D0
          WRITE(6,1904) PSFI(NPSN)
          WRITE(26,1904) PSFI(NPSN)
          OPEN(IPSFI,FILE=PSFI(NPSN),IOSTAT=KODE)
          IF(KODE.NE.0) THEN
            WRITE(26,'(''File '',A20,'' could not be opened.'')')
     1        PSFI(NPSN)
            STOP 'File could not be opened.'
          ENDIF
          KODEPS=0
          CALL RDPSF(IPSFI,NSHI,ISEC,KODEPS)
          IF(KODEPS.EQ.-1) THEN
            WRITE(26,*) ' STOP. The file is empty or corrupted.'
            STOP 'The file is empty or corrupted.'
          ENDIF
        ENDIF
      ENDIF
C
C  ************  Initialise constants.
C
      SHN=SHNA  ! Shower counter, particles from the dump file.
      N=MOD(SHN,2.0D9)+0.5D0
      TSIM=CPUTA
      CPUT0=CPUTIM()
      IF(SHN.GT.DSHN-0.5D0) THEN
        WRITE(6,*) '  **** The simulation was already completed.'
        WRITE(26,*) '  **** The simulation was already completed.'
        JOBEND=3
      ELSE
        WRITE(6,*) '  The simulation is started ...'
      ENDIF
C
      RETURN
      END


C  *********************************************************************
C                       SUBROUTINE SHOWER
C  *********************************************************************
      SUBROUTINE SHOWER
C
C  Simulates a new shower and keeps score of relevant quantities.
C
      USE PENELOPE_mod
      USE TRACK_mod
      USE PENGEOM_mod
      USE PENVARED_mod
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
      LOGICAL LINTF
C
C  ----  Seeds of the random number generator.
      COMMON/RSEED/ISEED1,ISEED2
C
      PARAMETER (REV=5.10998928D5)  ! Electron rest energy (eV)
      PARAMETER (TREV=2.0D0*REV)
C  ****  Current state in class II simulation.
C        MHINGE=0 (1) before (after) the hinge.
      COMMON/CJUMP1/ELAST1,ELAST2,MHINGE,KSOFTE,KSOFTI,KDELTA
C
      PARAMETER (PI=3.1415926535897932D0, TWOPI=2.0D0*PI)
      INCLUDE 'pmcomms.f'
      EXTERNAL RAND
C
C  ------------------------  Shower simulation starts here.
C
C  ************  Primary particle counters.
C
  101 CONTINUE
c     write(6,*) ISEED1,ISEED2
      DO I=1,3
        DPRIM(I)=0.0D0
        DO K=1,3
          DSEC(K,I)=0.0D0
        ENDDO
      ENDDO
C
      DO I=1,2
        DAVW(I)=0.0D0
        DAVA(I)=0.0D0
        DAVE(I)=0.0D0
      ENDDO
C
      DO KB=1,NBODY
        DEBO(KB)=0.0D0  ! Energies deposited in the various bodies KB.
      ENDDO
C
      IEXIT=0
      METAST=0
C
      CALL CLEANS  ! Cleans the secondary stack.
      IF(LAGE) CALL PAGE0  ! Sets the time of flight equal to 0.
C
C  **********  Set the initial state of the primary particle.
C
  201 CONTINUE
      IF(KPARP.EQ.0) THEN
C  ****  User-defined source.
        CALL SOURCE(METAST)
        SHN=SHN+1.0D0
        N=N+1
        IF(N.GT.2000000000) N=N-2000000000
        KEn=E*RDSDE+1.0D0
        SEDS(KPAR,KEn)=SEDS(KPAR,KEn)+WGHT
        SEDS2(KPAR,KEn)=SEDS2(KPAR,KEn)+WGHT**2
c     write(6,'(''n,kpar,gen,x,y,z,w,e,ibody='',3i3,1p,5e11.3,0p,i3)')
c    1    MOD(N,100),KPAR,ILB(1),X,Y,Z,W,E,IBODY
      ELSE IF(LPSF) THEN
C  ****  Phase-space file.
        CALL RDPSF(IPSFI,NSHI,ISEC,KODEPS)
        IF(KODEPS.EQ.-1) THEN
c       write(35,'(i2,1p,8e13.5,i3,i2,i2,i9,i9,2i3)')
c    1    KPAR,E,X,Y,Z,U,V,W,WGHT,ILB(1),ILB(2),ILB(3),ILB(4),
c    1    NSHI,ISEC,KODEPS
          CLOSE(IPSFI)
          IF(NPSN.EQ.NPSF) THEN
            JOBEND=4
            RETURN
          ENDIF
          NPSN=NPSN+1
          RLREAD=0.0D0
          OPEN(IPSFI,FILE=PSFI(NPSN))
          WRITE(6,3040) SHN
          WRITE(6,1904) PSFI(NPSN)
          WRITE(26,1904) PSFI(NPSN)
          KODEPS=0
          CALL RDPSF(IPSFI,NSHI,ISEC,KODEPS)
c       write(35,'(i2,1p,8e13.5,i3,i2,i2,i9,i9,2i3)')
c    1    KPAR,E,X,Y,Z,U,V,W,WGHT,ILB(1),ILB(2),ILB(3),ILB(4),
c    1    NSHI,ISEC,KODEPS
          GO TO 101
        ENDIF
 3040   FORMAT(2X,'Number of simulated showers =',1P,E14.7)
 1904   FORMAT(/3X,'+ The phase-space file ',A20,' is opened.')
C
        RLREAD=RLREAD+1.0D0
        SHN=SHN+DBLE(NSHI)
        N=N+NSHI
        IF(N.GT.2000000000) N=N-2000000000
        ILB(5)=0
        IPOL=0  ! Particles in the phase-space file are unpolarised.
      ELSE
C  ****  External source.
        SHN=SHN+1.0D0
        N=N+1
        IF(N.GT.2000000000) N=N-2000000000
        KPAR=KPARP
        WGHT=1.0D0
C  ----  Initial position ...
        IF(LEXSRC) THEN
          IF(LEXBD) THEN
            NTRIAL=0
  301       CONTINUE
            X=SX0+(RAND(1.0D0)-0.5D0)*SSX
            Y=SY0+(RAND(2.0D0)-0.5D0)*SSY
            Z=SZ0+(RAND(3.0D0)-0.5D0)*SSZ
            CALL locate
            NTRIAL=NTRIAL+1
            IF(NTRIAL.GT.200) THEN
              WRITE(26,'(3X,''WARNING: the sampling of initial '',
     1          ''positions may be very inefficient.'')')
              WRITE(6,'(3X,''WARNING: the sampling of initial '',
     1          ''positions may be very inefficient.'')')
            ENDIF
            IF(IXSBOD(IBODY).EQ.0) GO TO 301
          ELSE
            X=SX0+(RAND(1.0D0)-0.5D0)*SSX
            Y=SY0+(RAND(2.0D0)-0.5D0)*SSY
            Z=SZ0+(RAND(3.0D0)-0.5D0)*SSZ
          ENDIF
        ELSE
          X=SX0
          Y=SY0
          Z=SZ0
        ENDIF
C  ----  Initial direction ...
        IF(LSCONE) THEN
          CALL GCONE(U,V,W)  ! Conical beam.
        ELSE
          W=CTHL+RAND(4.0D0)*DCTH  ! Rectangular beam.
          UV=SQRT(1.0D0-W*W)
          PHI=PHIL+RAND(5.0D0)*DPHI
          U=UV*COS(PHI)
          V=UV*SIN(PHI)
        ENDIF
C  ----  Initial energy ...
        IF(LSPEC) THEN
          RN=RAND(6.0D0)*NSEB+1
          K=INT(RN) ! Continuous spectrum. E sampled by Walker's method.
          RNF=RN-K
          IF(RNF.GT.FSRC(K)) THEN
            KEn=IASRC(K)
          ELSE
            KEn=K
          ENDIF
          E=ESRC(KEn)+RAND(7.0D0)*(ESRC(KEn+1)-ESRC(KEn))
          SHIST(KEn)=SHIST(KEn)+1.0D0
        ELSE
          E=E0      ! Monoenergetic source.
          SHIST(1)=SHIST(1)+1.0D0
        ENDIF
        ILB(1)=1  ! Identifies primary particles.
        ILB(2)=0
        ILB(3)=0
        ILB(4)=0
        ILB(5)=0
        IF(KPAR.EQ.2) THEN
          IF(LGPOL) THEN
            IPOL=1  ! Photon polarisation.
            SP1=SP10
            SP2=SP20
            SP3=SP30
          ELSE
            IPOL=0
          ENDIF
        ELSE
          IPOL=0
        ENDIF
      ENDIF
C
C  ****  Check if the trajectory intersects the material system.
C
  302 CONTINUE
      CALL locate
c     write(6,'(''n,kpar,gen,x,y,z,w,e,ibody='',3i3,1p,5e11.3,0p,i3)')
c    1    MOD(N,100),KPAR,ILB(1),X,Y,Z,W,E,IBODY
C
      IF(MAT.EQ.0) THEN
        IBODYL=IBODY
        CALL STEP(1.0D30,DSEF,NCROSS)
c     write(6,'(''n,kpar,gen,x,y,z,w,e,ibody='',3i3,1p,5e11.3,0p,i3)')
c    1    MOD(N,100),KPAR,ILB(1),X,Y,Z,W,E,IBODY
        IF(LAGE) CALL DPAGE(DSEF,DSTOT)
        IF(MAT.EQ.0) THEN  ! The particle does not enter the system.
          IF(W.GT.0) THEN
            IEXIT=1        ! Labels emerging upbound particles.
          ELSE
            IEXIT=2        ! Labels emerging downbound particles.
          ENDIF
          GO TO 104        ! Exit.
        ENDIF

C  ----  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
C  ----  Impact detectors.
        IDET=KDET(IBODY)
        IF(IDET.NE.0) THEN
          IF(KDET(IBODYL).NE.IDET.AND.KKDI(IDET,KPAR).EQ.1) THEN
C
            IF(IPSF(IDET).EQ.1) THEN
              NSHJ=SHN-RLAST
              CALL WRPSF(IPSFO,NSHJ,0)
              RWRITE=RWRITE+1.0D0
              RLAST=SHN
            ENDIF
C
            CALL SIMDET(N,IDET)
C
            IF(IDCUT(IDET).EQ.0) THEN
              DEBO(IBODY)=DEBO(IBODY)+E*WGHT
              IEXIT=3
              GO TO 104
            ENDIF
          ENDIF
        ENDIF
C  ----  <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

      ENDIF
C
      IF(E.LT.EABSB(KPAR,IBODY)) THEN  ! The energy is too low.
        DEP=E*WGHT
        IF(KPAR.EQ.3.AND.E.GT.1.0D-6) THEN  ! Positrons annihilate
          CALL PANaR(EABSB(2,IBODY))  ! when absorbed.
          DEP=DEP+TREV*WGHT
        ENDIF
        E=0.0D0
C
        DEBO(IBODY)=DEBO(IBODY)+DEP
        IF(LDOSEM) THEN
          CALL SDOSE(DEP,X,Y,Z,MAT,N)
        ENDIF
        IEXIT=3
        GO TO 104
      ENDIF
C  ---------------------------------------------------------------------
C  ------------------------  Track simulation begins here.
C

C  >>>>>>>>>>>>  Particle splitting and Russian roulette  >>>>>>>>>>>>>>
C  ...  only for particles read from phase-space files.
      IF(LPSF) THEN
C  ----  Russian roulette (weight window).
        IF(WGHT.LT.WGMIN) THEN
          PSURV=WGHT*RWGMIN
          CALL VRR(PSURV)
          IF(E.LT.1.0D0) GO TO 201
        ENDIF
C  ----  Splitting.
        NSPL1=MAX(MIN(NSPLIT,INT(MIN(1.0D4,WGHT*RWGMIN))),1)
        IF(NSPL1.GT.1) THEN
          CALL VSPLIT(NSPL1)
C  ----  Energy is locally deposited in the material.
          DEP=(NSPL1-1)*E*WGHT
          DEBO(IBODY)=DEBO(IBODY)+DEP
          IF(LDOSEM) THEN  ! Particle inside the dose box.
            CALL SDOSE(DEP,X,Y,Z,MAT,N)
          ENDIF
        ENDIF
      ENDIF
C  <<<<<<<<<<<<  Particle splitting and Russian roulette  <<<<<<<<<<<<<<

  102 CONTINUE
C
C  ************  The particle energy is less than EABSB(KPAR,IBODY).
C
      IF(E.LT.EABSB(KPAR,IBODY)) THEN  ! The particle is absorbed.
        DEP=E*WGHT
        IF(KPAR.EQ.3.AND.E.GT.1.0D-6) THEN  ! Positrons annihilate
          CALL PANaR(EABSB(2,IBODY))  ! when absorbed.
          DEP=DEP+TREV*WGHT
        ENDIF
        E=0.0D0
C  ----  Energy is locally deposited in the material.
        DEBO(IBODY)=DEBO(IBODY)+DEP
        IF(LDOSEM) THEN  ! Particle inside the dose box.
          CALL SDOSE(DEP,X,Y,Z,MAT,N)
        ENDIF
        IEXIT=3                     ! Labels absorbed particles.
        GO TO 104                   ! Exit.
      ENDIF
C
      CALL START           ! Starts simulation in current medium.
C
C  ----  Free path length to the next interaction event.
C
  103 CONTINUE
      IBODYL=IBODY; MATL=MAT; XL=X; YL=Y; ZL=Z
c     write(6,'(''n,kpar,gen,x,y,z,w,e,ibody='',3i3,1p,5e11.3,0p,i3)')
c    1    MOD(N,100),KPAR,ILB(1),X,Y,Z,W,E,IBODY
      IF(LFORCE(IBODY,KPAR).AND.((WGHT.GT.WLOW(IBODY,KPAR)).AND.
     1  (WGHT.LT.WHIG(IBODY,KPAR)))) THEN
        CALL JUMPF(DSMAX(IBODY),DS)  ! Interaction forcing.
        LINTF=.TRUE.
      ELSE
        CALL JUMP(DSMAX(IBODY),DS)  ! Analogue simulation.
        LINTF=.FALSE.
      ENDIF
      CALL STEP(DS,DSEF,NCROSS)  ! Determines step end position.
      IF(LAGE) CALL DPAGE(DSEF,DSTOT)
C
C  ----  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
C  ----  Energy distribution of fluence.
      IDET=KDET(IBODYL)
      IF(IDET.NE.0) THEN
        IF(IDCUT(IDET).EQ.2) THEN
          IF(KKDI(IDET,KPAR).EQ.1) THEN
            IF(KPAR.EQ.2) THEN
              CALL FIMDET(N,IDET,DSEF)
            ELSE
              DECSD=SSOFT*DSEF
              IF(DECSD.GT.1.0D-12) THEN  ! The distribution may extend
                CALL FIMDES(N,IDET,E0STEP,DECSD,DSEF)  ! below EABS.
              ELSE
                CALL FIMDET(N,IDET,DSEF)
              ENDIF
            ENDIF
          ENDIF
        ENDIF
      ENDIF
C  ----  <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

C
C  ----  The particle has crossed an interface.
C
      IF(NCROSS.GT.0) THEN
C  ----  Correction to the soft energy loss (CSDA).
        IF(KPAR.NE.2) THEN
          E=E0STEP-SSOFT*DSEF
          IF(MHINGE.EQ.0) THEN
            DEP=SSOFT*DSEF*WGHT
            IF(LDOSEM) THEN
              DSEFR=RAND(8.0D0)*DSEF
              XD=XL+U*DSEFR
              YD=YL+V*DSEFR
              ZD=ZL+W*DSEFR
              CALL SDOSE(DEP,XD,YD,ZD,MATL,N)
            ENDIF
          ELSE
            DEP=-SSOFT*(DS-DSEF)*WGHT
            IF(LDOSEM) CALL SDOSE(DEP,XL,YL,ZL,MATL,N)
          ENDIF
          DEBO(IBODYL)=DEBO(IBODYL)+DEP
        ENDIF
C  ----  Check whether the particle is outside the enclosure.
        IF(MAT.EQ.0) THEN  ! The particle is outside the enclosure.
          IF(W.GT.0.0D0) THEN
            IEXIT=1        ! Labels emerging upbound particles.
          ELSE
            IEXIT=2        ! Labels emerging downbound particles.
          ENDIF
          GO TO 104        ! Exit.
        ENDIF

C  ----  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
C  ----  Impact detectors.
        IDET=KDET(IBODY)
        IF(IDET.NE.0) THEN
          IF(KDET(IBODYL).NE.IDET.AND.KKDI(IDET,KPAR).EQ.1) THEN
C
            IF(IPSF(IDET).EQ.1) THEN
              NSHJ=SHN-RLAST
              CALL WRPSF(IPSFO,NSHJ,0)
              RWRITE=RWRITE+1.0D0
              RLAST=SHN
            ENDIF
C
            CALL SIMDET(N,IDET)
C
            IF(IDCUT(IDET).EQ.0) THEN
              DEBO(IBODY)=DEBO(IBODY)+E*WGHT
              IEXIT=3
              GO TO 104
            ENDIF
          ENDIF
        ENDIF
C  ----  <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

        GO TO 102
      ENDIF
C
C  ----  Simulate next event.
C
      IF(LINTF) THEN
        CALL KNOCKF(DE,ICOL)  ! Interaction forcing is active.
      ELSE
        CALL KNOCK(DE,ICOL)  ! Analogue simulation.
      ENDIF
C
      IF(E.LT.EABSB(KPAR,IBODY)) THEN  ! The particle has been absorbed.
        DE=DE+E
        IF(KPAR.EQ.3.AND.E.GT.1.0D-6) THEN  ! Positrons annihilate
          CALL PANaR(EABSB(2,IBODY))  ! when absorbed.
          DE=DE+TREV
        ENDIF
        E=0.0D0
      ENDIF
C
      DEP=DE*WGHT
      DEBO(IBODY)=DEBO(IBODY)+DEP
      IF(LDOSEM) CALL SDOSE(DEP,X,Y,Z,MAT,N)
C
      IF(E.LT.EABSB(KPAR,IBODY)) THEN  ! The particle has been absorbed.
        IEXIT=3            ! Labels absorbed particles.
        GO TO 104          ! Exit.
      ENDIF
C
      GO TO 103
C  ------------------------  The simulation of the track ends here.
C  ---------------------------------------------------------------------
  104 CONTINUE
C
C  ************  Increment particle counters.
C
      IF(ILB(1).EQ.1) THEN
        DPRIM(IEXIT)=DPRIM(IEXIT)+WGHT
        IF(LPSF) THEN
          IF(NSPL1.GT.1) THEN
            DPRIM(IEXIT)=DPRIM(IEXIT)+WGHT*(NSPL1-1)
          ENDIF
        ENDIF
        IF(IEXIT.LT.3) THEN
          DAVW(IEXIT)=DAVW(IEXIT)+W*WGHT
          DAVA(IEXIT)=DAVA(IEXIT)+ACOS(W)*WGHT
          DAVE(IEXIT)=DAVE(IEXIT)+E*WGHT
        ENDIF
      ELSE
        DSEC(KPAR,IEXIT)=DSEC(KPAR,IEXIT)+WGHT
      ENDIF
C
      IF(IEXIT.LT.3) THEN
        CALL TENANG(IEXIT,N)  ! Energy and angular distributions.
      ENDIF
C
C  ************  Any secondary left?
C
  202 CONTINUE
      CALL SECPAR(LEFT)
      IF(LEFT.GT.0) THEN
c     write(6,'(/''new secondary'')')
c     write(6,'(''n,kpar,ilb(1:4),ibody,mat='',i6,4i4,i10,2i4)')
c    1    N,KPAR,ILB(1),ILB(2),ILB(3),ILB(4),IBODY,MAT
c     write(6,'(''  ener,x,y,z='',1p,5e14.6)') E,X,Y,Z
c     write(6,'(''  wght,u,v,w='',1p,5e14.6)') WGHT,U,V,W
        IF(ILB(1).EQ.1) THEN  ! Primary particle from SOURCE.
          KEn=E*RDSDE+1.0D0
          SEDS(KPAR,KEn)=SEDS(KPAR,KEn)+WGHT
          SEDS2(KPAR,KEn)=SEDS2(KPAR,KEn)+WGHT**2
          IF(LAGE) CALL PAGE0
          GO TO 302  ! Energy is not removed from the site.
        ENDIF
        IF(E.GT.EABSB(KPAR,IBODY)) THEN
C  ****  X-ray splitting.
          IF(KPAR.EQ.2) THEN
            IF(ILB(4).GT.0) THEN  ! Characteristic x rays.
              IF(ILB(1).EQ.2.AND.ILB(3).LT.9) THEN  ! Unsplit, 2nd gen.
                IF(LXRSPL(IBODY)) THEN
                  WGHT=WGHT/DBLE(IXRSPL(IBODY))
                  ILBA(1)=ILB(1)
                  ILBA(2)=ILB(2)
                  ILBA(3)=9
                  ILBA(4)=ILB(4)
                  ILBA(5)=ILB(5)
                  DO I=2,IXRSPL(IBODY)
                    WS=-1.0D0+2.0D0*RAND(9.0D0)
                    SDTS=SQRT(1.0D0-WS*WS)
                    DF=TWOPI*RAND(10.0D0)
                    US=COS(DF)*SDTS
                    VS=SIN(DF)*SDTS
                    CALL STORES(E,X,Y,Z,US,VS,WS,WGHT,KPAR,ILBA,0)
                  ENDDO
                ENDIF
              ENDIF
            ENDIF
          ENDIF
C
          DEP=E*WGHT
          DEBO(IBODY)=DEBO(IBODY)-DEP  ! Energy is removed.
          IF(LDOSEM) THEN  ! Particle inside the dose box.
            CALL SDOSE(-DEP,X,Y,Z,MAT,N)
          ENDIF
        ELSE
          GO TO 202
        ENDIF

C  >>>>>>>>>>>>>>>>>>>>>>>  Russian roulette  >>>>>>>>>>>>>>>>>>>>>>>>>>
C  ****  Russian roulette for photons moving downstream.
cr      IF(KPAR.EQ.2.AND.W.LT.0.0D0) THEN
cr        IF(WGHT.LT.0.1D0) THEN
cr          PSURV=0.25D0
cr          CALL VRR(PSURV)
cr          IF(E.LT.EABSB(KPAR,IBODY)) GO TO 202
cr        ENDIF
cr      ENDIF
C  <<<<<<<<<<<<<<<<<<<<<<<  Russian roulette  <<<<<<<<<<<<<<<<<<<<<<<<<<

        GO TO 102
      ENDIF
C
      IF(LPSF) THEN
        IF(ISEC.EQ.1) GO TO 201
      ENDIF
C
C  ----  Energies deposited in different bodies and detectors.
C
C  ----  Tallying the spectra from energy-deposition detectors.
      IF(NED.GT.0) THEN
        DO KD=1,NED
          DEDE(KD)=0.0D0
        ENDDO
        DO KB=1,NBODY
          IDET=KBDE(KB)
          IF(IDET.NE.0) THEN
            DEDE(IDET)=DEDE(IDET)+DEBO(KB)
          ENDIF
        ENDDO
C
        DO IDET=1,NED
          CALL SENDET(DEDE(IDET),IDET)
        ENDDO
      ENDIF
C
      DO KB=1,NBODY
        TDEBO(KB)=TDEBO(KB)+DEBO(KB)
        TDEBO2(KB)=TDEBO2(KB)+DEBO(KB)**2
      ENDDO
C  --  Final state counters.
      DO I=1,3
        PRIM(I)=PRIM(I)+DPRIM(I)
        PRIM2(I)=PRIM2(I)+DPRIM(I)**2
        DO K=1,3
          SEC(K,I)=SEC(K,I)+DSEC(K,I)
          SEC2(K,I)=SEC2(K,I)+DSEC(K,I)**2
        ENDDO
      ENDDO
C
      DO I=1,2
        AVW(I)=AVW(I)+DAVW(I)
        AVW2(I)=AVW2(I)+DAVW(I)**2
        AVA(I)=AVA(I)+DAVA(I)
        AVA2(I)=AVA2(I)+DAVA(I)**2
        AVE(I)=AVE(I)+DAVE(I)
        AVE2(I)=AVE2(I)+DAVE(I)**2
      ENDDO
C
C  ------------------------  The simulation of the shower ends here.
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C  ****  Source particle from decay of metastable states (mean life
C  longer than detector resolution time).
C
      IF(METAST.NE.0) THEN
        DO I=1,3
          DPRIM(I)=0.0D0
          DO K=1,3
            DSEC(K,I)=0.0D0
          ENDDO
        ENDDO
C
        DO I=1,2
          DAVW(I)=0.0D0
          DAVA(I)=0.0D0
          DAVE(I)=0.0D0
        ENDDO
C
        DO KB=1,NBODY
          DEBO(KB)=0.0D0  ! Energies deposited in the various bodies KB.
        ENDDO
C
        IEXIT=0
C
*       CALL CLEANS          ! Cleans the secondary stack.
        CALL SOURCE(METAST)
        ILB(1)=1
        KEn=E*RDSDE+1.0D0
        SEDS(KPAR,KEn)=SEDS(KPAR,KEn)+WGHT
        SEDS2(KPAR,KEn)=SEDS2(KPAR,KEn)+WGHT**2
        GO TO 302
      ENDIF
C
      RETURN
      END


C  *********************************************************************
C                       SUBROUTINE PMWRT
C  *********************************************************************
      SUBROUTINE PMWRT(ICLOSE)
C
C  Computes averages and writes results in output files.
C
C  ICLOSE is a closing flag, which is used only when the phase-space
C  file of an impact detector is being generated.
C  -- When ICLOSE .GT. 0, particles remaining in the buffer are
C     transferred to the psf, and the psf output unit is closed.
C  -- If ICLOSE .LT. 0, particles are moved to the psf, but the psf
C     output unit remains open.
C
      USE PENELOPE_mod
      USE TRACK_mod
      USE PENGEOM_mod
      USE PENVARED_mod
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
      PARAMETER (PI=3.1415926535897932D0, RA2DE=180.0D0/PI)
C
C  ----  Seeds of the random number generator.
      COMMON/RSEED/ISEED1,ISEED2
C
      INCLUDE 'pmcomms.f'
C
      DIMENSION WSEC(3,3),WSEC2(3,3)
      DIMENSION WAVE2(2),WAVE(2),WAVW2(2),WAVW(2),WAVA2(2),WAVA(2)
C
      IF(NPSFO.GT.0) THEN
        CALL WRPSF(IPSFO,0,ICLOSE)
      ENDIF
C
C  ************  If 'DUMPTO' is active, write counters in a dump file.
C
      IF(LDUMP) THEN
        OPEN(9,FILE=PFILED)
        WRITE(9,*) SHN,TSIM
        WRITE(9,'(A65)') TITLE
        WRITE(9,*) ISEED1,ISEED2
        WRITE(9,*) NPSN,RLREAD
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
      IF(LSPEC) THEN
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
C                       SUBROUTINE GCONE
C  *********************************************************************
      SUBROUTINE GCONE(UF,VF,WF)
C
C  This subroutine samples a random direction uniformly within a cone
C  with central axis in the direction (THETA,PHI) and aperture ALPHA.
C  Parameters are initialised by calling subroutine GCONE0.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
      PARAMETER (PI=3.1415926535897932D0, TWOPI=2.0D0*PI)
C  ****  Parameters for sampling directions within a cone.
      COMMON/CGCONE/CPCT,CPST,SPCT,SPST,SPHI,CPHI,STHE,CTHE,CAPER
C
      EXTERNAL RAND
C  ****  Define a direction relative to the z-axis.
      WT=CAPER+(1.0D0-CAPER)*RAND(1.0D0)
      DF=TWOPI*RAND(2.0D0)
      SUV=SQRT(1.0D0-WT*WT)
      UT=SUV*COS(DF)
      VT=SUV*SIN(DF)
C  **** Rotate to the beam axis direction.
      UF=CPCT*UT-SPHI*VT+CPST*WT
      VF=SPCT*UT+CPHI*VT+SPST*WT
      WF=-STHE*UT+CTHE*WT
C  ****  Ensure normalisation.
      DXY=UF*UF+VF*VF
      DXYZ=DXY+WF*WF
      IF(ABS(DXYZ-1.0D0).GT.1.0D-14) THEN
        FNORM=1.0D0/SQRT(DXYZ)
        UF=FNORM*UF
        VF=FNORM*VF
        WF=FNORM*WF
      ENDIF
      RETURN
      END
C  *********************************************************************
C                       SUBROUTINE GCONE0
C  *********************************************************************
      SUBROUTINE GCONE0(THETA,PHI,ALPHA)
C
C  This subroutine defines the parameters for sampling random directions
C  uniformly within a cone with axis in the direction (THETA,PHI) and
C  aperture ALPHA (in rad).
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
C  ****  Parameters for sampling directions within a cone.
      COMMON/CGCONE/CPCT,CPST,SPCT,SPST,SPHI,CPHI,STHE,CTHE,CAPER
C
      CPCT=COS(PHI)*COS(THETA)
      CPST=COS(PHI)*SIN(THETA)
      SPCT=SIN(PHI)*COS(THETA)
      SPST=SIN(PHI)*SIN(THETA)
      SPHI=SIN(PHI)
      CPHI=COS(PHI)
      STHE=SIN(THETA)
      CTHE=COS(THETA)
      CAPER=COS(ALPHA)
      RETURN
      END


C  *********************************************************************
C                       SUBROUTINE RDPSF
C  *********************************************************************
      SUBROUTINE RDPSF(IUNIT,NSHI,ISEC,KODEPS)
C
C  This subroutine reads the phase-space file (psf). Blank lines and
C  lines starting with the pound sign (#) are considered as comments and
C  are ignored.
C
C  The reading of a new psf is initiated by calling subroutine RDPSF
C  with KODEPS=0. If the psf contains some particles, the oputput value
C  of KODEPS is set equal to 1. In subsequent calls to subroutine RDPSF,
C  the output value KODEPS=1 indicates that a valid record has been read
C  and that the state variables of a particle have been loaded in common
C  TRACK. The output value of ISEC is 1 if the psf contains more
C  particles from the current shower, and 0 otherwise. That is, when
C  ISEC=0 the next particle in the psf belongs to a new shower. The
C  output value KODEPS=-1 indicates that all particle records in the psf
C  have been read (or that the psf was empty). Note that when KODEPS=-1,
C  the contents of module TRACK_mod remains unchanged.
C
C  To reduce disk-access time, particle records are read in bunches and
C  stored in a buffer. The buffer size (number of particles in a bunch)
C  is defined by the parameter NRBUFF.
C
C  To verify that the input psfs are read correctly, uncomment the lines
C  with 'CALL WRPSF' (which write the read particle records in an output
C  file, UNIT=13). You should also add a sentence 'CALL WRPSF(13,0,1)'
C  at the end of the main program, which forces the writing of particle
C  records that may remain in the buffer and closes the output file.
C  This file should contain all the records read from the psf files.
C
      USE TRACK_mod
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
C
      PARAMETER (NRBUFF=1000)  ! Buffer size.
      CHARACTER PSFR(NRBUFF)*136
      SAVE PSFR
C
      SAVE E0,X0,Y0,Z0,U0,V0,W0,WGHT0,KPAR0,ILB10,ILB20,ILB30,ILB40,
     1  NSHI0,ILAST,IRD,IEOF
C
C  ****  Initialisation (first call after opening the psf).
C
      IF(KODEPS.EQ.0) THEN
        NSHI=0
        ISEC=0
        IEOF=0
        DO I=1,NRBUFF
         PSFR(I)='#'  ! Identifies empty records.
        ENDDO
        READ(IUNIT,'(A)',END=10,ERR=10) (PSFR(I),I=1,NRBUFF)
        GO TO 11
   10   CONTINUE
        IEOF=1
   11   CONTINUE
C
        IRD=0
   20   CONTINUE
        IF(IRD.GE.NRBUFF) THEN
          IF(IEOF.EQ.1) THEN
            KODEPS=-1  ! End of file.
            RETURN
          ENDIF
          DO I=1,NRBUFF
           PSFR(I)='#'  ! Identifies empty records.
          ENDDO
          READ(IUNIT,'(A)',END=21,ERR=21) (PSFR(I),I=1,NRBUFF)
          GO TO 22
   21     CONTINUE
          IEOF=1
   22     CONTINUE
          IRD=0
        ENDIF
        IRD=IRD+1
        READ(PSFR(IRD),*,ERR=20,END=20) KPAR0,E0,X0,Y0,Z0,U0,V0,W0,
     1    WGHT0,ILB10,ILB20,ILB30,ILB40,NSHI0
        ILAST=0
        KODEPS=1
        RETURN
      ENDIF
C
C  ****  The last particle in the psf has already been loaded.
C
      IF(ILAST.EQ.1) THEN
        NSHI=0
        ISEC=0
        KODEPS=-1  ! End of file.
C       CALL WRPSF(13,NSHI,-1)  ! Uncomment to check psf reading.
        RETURN
      ENDIF
C
C  ****  A read particle is loaded.
C
      E=E0
      X=X0
      Y=Y0
      Z=Z0
      U=U0
      V=V0
      W=W0
      WGHT=WGHT0
      KPAR=KPAR0
      ILB(1)=ILB10
      ILB(2)=ILB20
      ILB(3)=ILB30
      ILB(4)=ILB40
      NSHI=NSHI0
C     CALL WRPSF(13,NSHI,0)  ! Uncomment to check psf reading.
      KODEPS=1
C
C  ****  Read a new particle and keep its state variables in memory.
C        This is needed to identify possible secondary particles that
C        belong to the same shower.
C
   30 CONTINUE
      IF(IRD.GE.NRBUFF) THEN
        IF(IEOF.EQ.1) THEN
          ILAST=1
          ISEC=0
          RETURN
        ENDIF
        DO I=1,NRBUFF
         PSFR(I)='#'  ! Identifies empty records.
        ENDDO
        READ(IUNIT,'(A)',END=31,ERR=31) (PSFR(I),I=1,NRBUFF)
        GO TO 32
   31   CONTINUE
        IEOF=1
   32   CONTINUE
        IRD=0
      ENDIF
      IRD=IRD+1
      READ(PSFR(IRD),*,ERR=30,END=30) KPAR0,E0,X0,Y0,Z0,U0,V0,W0,
     1    WGHT0,ILB10,ILB20,ILB30,ILB40,NSHI0
C
      IF(NSHI0.EQ.0) THEN
        ISEC=1
      ELSE
        ISEC=0
      ENDIF
      RETURN
      END
C  *********************************************************************
C                       SUBROUTINE WRPSF
C  *********************************************************************
      SUBROUTINE WRPSF(IPSFO,NSHI,ICLOSE)
C
C  This subroutine writes particle records in the output phase-space
C  file (psf). To reduce disk-access time, particle records are
C  collected in a buffer and moved to the psf when the buffer is full.
C  The buffer size (number of records) is set by the parameter NRBUFF.
C
C  Particles that remain in the buffer at the end of a simulation run
C  must be transferred to the psf by calling subroutine WRPSF with
C  ICLOSE=1.
C
C  Input parameters:
C     IPSFO ..... output unit.
C     NSHI ...... incremental shower number.
C     ICLOSE .... closing flag, acts only when its value is not 0.
C                 When ICLOSE .GT. 0, particles remaining in the buffer
C                 are transferred to the psf, and the output unit IPSFO
C                 is closed.
C                 If ICLOSE .LT. 0, particles in the buffer are moved to
C                 the psf, but the output unit remains open.
C
      USE TRACK_mod
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
      CHARACTER*10 LC10A,LC10B
C
      PARAMETER (NRBUFF=1000)  ! Buffer size.
      CHARACTER PSFR(NRBUFF)*136
      SAVE PSFR
C
      DATA NREC/0/
      SAVE NREC
C
      IF(ICLOSE.NE.0) THEN
        IF(NREC.GT.0) THEN
          DO I=1,NREC
            NCHAR=INDEX(PSFR(I),'*')-1
            WRITE(IPSFO,'(A)') PSFR(I)(1:NCHAR)
          ENDDO
        ENDIF
        IF(ICLOSE.GT.0) CLOSE(IPSFO)
        NREC=0
        RETURN
      ENDIF
C
      NREC=NREC+1
      CALL N2CH10(ILB(4),LC10A,NDIGA)
      CALL N2CH10(NSHI,LC10B,NDIGB)
      WRITE(PSFR(NREC),'(I2,1P,8E13.5,I3,I2,I2,1X,A,1X,A,A)')
     1  KPAR,E,X,Y,Z,U,V,W,WGHT,ILB(1),ILB(2),ILB(3),
     1  LC10A(1:NDIGA),LC10B(1:NDIGB),'*'
C
      IF(NREC.EQ.NRBUFF) THEN
        WRITE(IPSFO,'(A)') (PSFR(I)(1:INDEX(PSFR(I),'*')-1),I=1,NREC)
        NREC=0
      ENDIF
C
      RETURN
      END
C  *********************************************************************
C                       SUBROUTINE N2CH10
C  *********************************************************************
      SUBROUTINE N2CH10(N,L,NDIG)
C
C  This subroutine writes an integer number N in a 10-character string
C  L. The number is written at the left, followed by unused blanks. NDIG
C  is the number of decimal digits of N.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
      CHARACTER*10 L,LT
C
      WRITE(LT,'(I10)') N
      DO I=1,10
        IF(LT(I:I).NE.' ') THEN
          IT=I-1
          GO TO 1
        ENDIF
      ENDDO
      IT=9
    1 CONTINUE
      L=LT(IT+1:10)
      NDIG=10-IT
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
C  ----  Upbound particles.
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
     1      1X,'#  Results from PENMAIN.',
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
      USE PENELOPE_mod
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
      USE PENELOPE_mod
      USE TRACK_mod
	  
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
            DDOSEP(I3)=DEP*RDEN(MATC)
            LDDOSE(I3)=N
          ELSE
            DDOSEP(I3)=DDOSEP(I3)+DEP*RDEN(MATC)
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
            DDOSEP(I3)=DEP*RDEN(MATC)
            LDDOSE(I3)=N
          ELSE
            DDOSEP(I3)=DDOSEP(I3)+DEP*RDEN(MATC)
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
C
C  ****  Voxel masses.
C
C  ----  Dose in a box.
        IF(KDOSE.EQ.1) THEN
          DX=BDOSE(1); DY=BDOSE(2); DZ=BDOSE(3)
          VOXEL=DX*DY*DZ
          FNORM=VOXEL*RNCS**3
          DDX=RNCS*DX; DDY=RNCS*DY; DDZ=RNCS*DZ
          DO K=1,NDB(3)
            DO J=1,NDB(2)
              DO I=1,NDB(1)
                UO=1.0D0; VO=1.0D0
                X=DXL(1)+(I-1)*DX+0.5D0*DDX
                Y=DXL(2)+(J-1)*DY+0.5D0*DDY
                Z=DXL(3)+(K-1)*DZ+0.5D0*DDZ
                CALL LOCATE
                IF(MAT.EQ.0) THEN
                  DENT=0.0D0
                ELSE
                  DENT=DEN(MAT)
                ENDIF
                DO J3=1,NCS
                  DO J2=1,NCS
                    DO J1=1,NCS-1
                      U=UO; V=0.0D0; W=0.0D0
                      X=X+DDX*U
                      CALL LOCATE
                      IF(MAT.GT.0) DENT=DENT+DEN(MAT)
                    ENDDO
                    UO=-UO
                    IF(J2.LT.NCS) THEN
                      U=0.0D0; V=VO; W=0.0D0
                      Y=Y+DDY*V
                      CALL LOCATE
                      IF(MAT.GT.0) DENT=DENT+DEN(MAT)
                    ENDIF
                  ENDDO
                  VO=-VO
                  IF(J3.LT.NCS) THEN
                    U=0.0D0; V=0.0D0; W=1.0D0
                    Z=Z+DDZ*W
                    CALL LOCATE
                    IF(MAT.GT.0) DENT=DENT+DEN(MAT)
                  ENDIF
                ENDDO
                VMASS(I,J,K)=DENT*FNORM
              ENDDO
            ENDDO
          ENDDO
C
          DO K=1,NDB(3)
            DO J=1,NDB(2)
              DO I=1,NDB(1)
                IF(VMASS(I,J,K).GT.1.0D-35) THEN
                  VMASS(I,J,K)=1.0D0/VMASS(I,J,K)
                ELSE
                  VMASS(I,J,K)=0.0D0
                ENDIF
              ENDDO
            ENDDO
          ENDDO
C
C  ----  Dose in a cylinder.
        ELSE IF(KDOSE.EQ.2) THEN
          J=1
          DR=BDOSE(1); DZ=BDOSE(3)
          DDR=RNCS*DR; DDZ=RNCS*DZ
          U=1.0D0; V=0.0D0; W=0.0D0
          DO K=1,NDB(3)
            DO I=1,NDB(1)
              TMASS=0.0D0
              DO J3=1,NCS
                DO J1=1,NCS
                  X=(I-1)*DR+(DBLE(J1)-0.5D0)*DDR
                  Y=0.0D0
                  Z=DXL(3)+(K-1)*DZ+(DBLE(J3)-0.5D0)*DDZ
                  CALL LOCATE
                  IF(MAT.NE.0) THEN
                    VOLUM=2.0D0*PI*X*DDR*DDZ
                    TMASS=TMASS+DEN(MAT)*VOLUM
                  ENDIF
                ENDDO
              ENDDO
              VMASS(I,J,K)=TMASS
              IF(VMASS(I,J,K).GT.1.0D-35) THEN
                VMASS(I,J,K)=1.0D0/VMASS(I,J,K)
              ELSE
                VMASS(I,J,K)=0.0D0
              ENDIF
            ENDDO
          ENDDO
C  ----  Dose in a sphere.
        ELSE
          J=1
          K=1
          DR=BDOSE(1)
          DDR=RNCS*DR
          U=0.0D0; V=0.0D0; W=1.0D0
          DO I=1,NDB(1)
            TMASS=0.0D0
            DO J1=1,NCS
              X=(I-1)*DR+(DBLE(J1)-0.5D0)*DDR
              Y=0.0D0
              Z=0.0D0
              CALL LOCATE
              IF(MAT.NE.0) THEN
                VOLUM=3.0D0*X**2+0.25D0*DDR**2
                TMASS=TMASS+DEN(MAT)*VOLUM
              ENDIF
            ENDDO
            VMASS(I,J,K)=TMASS*(4.0D0*PI/3.0D0)*DDR
            IF(VMASS(I,J,K).GT.1.0D-35) THEN
              VMASS(I,J,K)=1.0D0/VMASS(I,J,K)
            ELSE
              VMASS(I,J,K)=0.0D0
            ENDIF
          ENDDO
        ENDIF
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


C  *********************************************************************
C                       SUBROUTINE DPAGE
C  *********************************************************************
      SUBROUTINE DPAGE(DSEF,DSTOT)
C
C  This subroutine keeps track of the age of particles, measured from
C  the start of the primary particle that originated the shower. At
C  output, the variable PAGE in the module TRACK_mod contains the age
C  of the particle at the end of the last step (in seconds).
C
C  The input parameters DSEF and DSTOT are, respectively, the length of
C  the step in the 'original' material and the total length of the step,
C  including segments in void volumes. When using PENGEOM, DSTOT is
C  given by subroutine STEP through module PENGEOM_mod.
C
C  Usage:
C  1) CALL PAGE0 just after starting each primary particle.
C  2) CALL DPAGE(DSEF,DSTOT) at the end of each free flight.
C
      USE PENELOPE_mod
      USE TRACK_mod
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
      PARAMETER (SLCM=2.99792458D10)  ! Speed of light (cm/s)
      PARAMETER (RSLCM=1.0D0/SLCM)  ! Reciprocal of SLCM (1/cm).
      PARAMETER (REV=5.10998928D5)  ! Electron rest energy (eV)
      PARAMETER (TREV=2.0D0*REV)
C
C  ****  Updating the particle's age after each free flight.
C
      IF(KPAR.EQ.2) THEN
        IF(MAT.EQ.0) THEN
          PAGE=PAGE+DSEF*RSLCM
        ELSE
          PAGE=PAGE+DSTOT*RSLCM
        ENDIF
      ELSE
        IF(SSOFT.GT.1.0D-16) THEN
          E0=E0STEP
          SQE=SQRT(E0*(E0+TREV))
          E1=E0-SSOFT*DSEF  ! Energy at the end of the step, S=DSEF.
          IF(E0.LT.1.000001D0*E1) THEN
            PAGE=PAGE+DSTOT*RSLCM*(E0+REV)/SQE
          ELSE
            IF(E1.LT.50.0D0) E1=50.0D0
            SQE1=SQRT(E1*(E1+TREV))
            PAGE=PAGE+RSLCM*(SQE-SQE1)/SSOFT
            IF(DSTOT.GT.DSEF.AND.MAT.NE.0)
     1        PAGE=PAGE+(DSTOT-DSEF)*RSLCM*(E1+REV)/SQE1
          ENDIF
        ELSE
          SQE=SQRT(E*(E+TREV))
          IF(MAT.EQ.0) THEN
            PAGE=PAGE+DSEF*RSLCM*(E+REV)/SQE
          ELSE
            PAGE=PAGE+DSTOT*RSLCM*(E+REV)/SQE
          ENDIF
        ENDIF
      ENDIF
      RETURN
C
C  ****  PAGE is initialised when the primary particle starts moving.
C
      ENTRY PAGE0
      PAGE=0.0D0
      SSOFT=0.0D0  ! The particle may start in vacuum.
      RETURN
      END


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                       SUBROUTINE FWORD
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE FWORD(STRING,TAIL,WORD,LENGTH)
C
C  This subroutine extracts the First WORD (substring between enclosing
C  apostrophes) of the string STRING, if it contains one. LENGTH is the
C  length of the WORD, excluding the apostrophes. TAIL is the remainder
C  of STRING after removing the word and leading blanks, commas, and
C  apostrophes.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
      CHARACTER*(*) STRING,TAIL,WORD
      CHARACTER TMP*132
C
      WORD=' '
      LENGTH=0
      I1=INDEX(STRING,'''')
      IF(I1.EQ.0) RETURN
      TMP=STRING(I1+1:)
      I2=INDEX(TMP,'''')
      LENGTH=I2-1
      IF(LENGTH.LT.1) RETURN
      WORD=TMP(1:LENGTH)
      TAIL=TMP(LENGTH+1:)
C
C  ****  Remove leading blanks, commas and apostrophes.
C
    1 CONTINUE
      IF(TAIL(1:1).EQ.''''.OR.TAIL(1:1).EQ.' '.OR.TAIL(1:1).EQ.',') THEN
        TMP=TAIL(2:)
        TAIL=TMP
        GO TO 1
      ENDIF
C
      RETURN
      END
	  
	  
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                       SUBROUTINE TABELAS
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC	  
      SUBROUTINE TABELAS
      USE PENELOPE_mod
      USE TRACK_mod
      USE PENGEOM_mod
      USE PENVARED_mod
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
	  PARAMETER (NRP=8000)
	  CHARACTER*2 LASYMB
	  PARAMETER (NTP=12000)
	  PARAMETER (NRX=60000)
	  PARAMETER (NM=512)
	  PARAMETER (NR=128)
	  PARAMETER (NBE=57, NBW=32)
	  PARAMETER (NQ=250,NEX=1024)
	  COMMON/QSURF/AXX(NS),AXY(NS),AXZ(NS),AYY(NS),AYZ(NS),AZZ(NS),
     1    AX(NS),AY(NS),AZ(NS),A0(NS),NSURF,KPLANE(NS)
      COMMON/QBODY/KBODY(NB,NXG),KBOMO(NB)
      COMMON/QTREE/NBODYS,KMOTH(NB),KDGHT(NB,NXG),KSURF(NB,NXG),
     1  KFLAG(NB,NXG),KSP(NS),NWARN
	 
	  COMMON/CECUTR/ECUTR(MAXMAT)
	  COMMON/CSGAWR/ISGAW  ! Controls warning messages from SUMGA.
      COMMON/CERSEC/IERSEC
	  
	  COMMON/CEGRID/EMIN,EL,EU,ET1(NEGP),DLEMP(NEGP),DLEMP1,DLFC,
     1  XEL,XE,XEK,KE
	  COMMON/CESI0/XESI(NRP,16),IESIF(99),IESIL(99),NSESI(99),NCURE
	  COMMON/CPSI0/XPSI(NRP,16),IPSIF(99),IPSIL(99),NSPSI(99),NCURP
	
      COMMON/CADATA/ATW(99),EPX(99),RSCR(99),ETA(99),EB(99,30),
     1  ALW(99,30),CP0(99,30),IFI(99,30),IKS(99,30),NSHT(99),LASYMB(99)
	 
      COMMON/CGPH00/EPH(NTP),XPH(NTP,17),IPHF(99),IPHL(99),NPHS(99),NCUR
     COMMON/CRELAX/P(NRX),ET(NRX),F(NRX),IS0(NRX),IS1(NRX),
     1  IS2(NRX),IFIRST(99,16),ILAST(99,16),NCUR,KS,MODER

	  COMMON/CRITAA/XA(NM),AA(NM),BA(NM),FA(NM),IA(NM),NPM1A
	  COMMON/CRITAN/CNORM
	  
	  
	  
C  ****  Composition data.
      COMMON/COMPOS/STF(MAXMAT,30),ZT(MAXMAT),AT(MAXMAT),RHO(MAXMAT),
     1  VMOL(MAXMAT),IZ(MAXMAT,30),NELEM(MAXMAT)

      COMMON/CRANGE/RANGE(3,MAXMAT,NEGP),RANGEL(3,MAXMAT,NEGP)
C  ****  E/P inelastic collisions.
      PARAMETER (NO=512)
      COMMON/CEIN/EXPOT(MAXMAT),OP2(MAXMAT),F1(MAXMAT,NO),UI(MAXMAT,NO),
     1  WRI(MAXMAT,NO),KZ(MAXMAT,NO),KS1(MAXMAT,NO),NOSC(MAXMAT)
      COMMON/CEINTF/T1EI(NEGP),T2EI(NEGP),T1PI(NEGP),T2PI(NEGP)
C  ****  Partial cross sections of individual shells/oscillators.
      COMMON/CEIN00/SEH0(NO),SEH1(NO),SEH2(NO),SES0(NO),SES1(NO),
     1              SES2(NO),SET0(NO),SET1(NO),SET2(NO)
      COMMON/CPIN00/SPH0(NO),SPH1(NO),SPH2(NO),SPS0(NO),SPS1(NO),
     1              SPS2(NO),SPT0(NO),SPT1(NO),SPT2(NO)
C  ****  Inner-shell ionisation by electron and positron impact.
C  ****  
      COMMON/CEINAC/EINAC(MAXMAT,NEGP,NO),IEIN(MAXMAT,NO),NEIN(MAXMAT)
      COMMON/CESIAC/ESIAC(MAXMAT,NEGP,NO),IESI(MAXMAT,NO),NESI(MAXMAT)
      COMMON/CESIN/XSEIN(NEGP,NO),XSESI(NEGP,NO),ISIE(NO)
C  ****  Positron inelastic coll. and inner-shell ionisation tables.
      COMMON/CPINAC/PINAC(MAXMAT,NEGP,NO),IPIN(MAXMAT,NO),NPIN(MAXMAT)
      COMMON/CPSIAC/PSIAC(MAXMAT,NEGP,NO),IPSI(MAXMAT,NO),NPSI(MAXMAT)
      COMMON/CPSIN/XSPIN(NEGP,NO),XSPSI(NEGP,NO),ISIP(NO)
	  
	  PARAMETER (NOCO=512)
      COMMON/CGCO/FCO(MAXMAT,NOCO),UICO(MAXMAT,NOCO),FJ0(MAXMAT,NOCO),
     1  PTRSH(MAXMAT,NOCO),KZCO(MAXMAT,NOCO),KSCO(MAXMAT,NOCO),
     2  NOSCCO(MAXMAT)
C  ****  Electron simulation tables.
      COMMON/CEIMFP/SEHEL(MAXMAT,NEGP),SEHIN(MAXMAT,NEGP),
     1  SEISI(MAXMAT,NEGP),SEHBR(MAXMAT,NEGP),SEAUX(MAXMAT,NEGP),
     2  SETOT(MAXMAT,NEGP),CSTPE(MAXMAT,NEGP),RSTPE(MAXMAT,NEGP),
     3  DEL(MAXMAT,NEGP),W1E(MAXMAT,NEGP),W2E(MAXMAT,NEGP),
     4  DW1EL(MAXMAT,NEGP),DW2EL(MAXMAT,NEGP),
     5  RNDCE(MAXMAT,NEGP),AE(MAXMAT,NEGP),BE(MAXMAT,NEGP),
     6  T1E(MAXMAT,NEGP),T2E(MAXMAT,NEGP)
      COMMON/CLAS1E/TSTPE(MAXMAT,NEGP),TSTRE(MAXMAT,NEGP),
     1  TRL1E(MAXMAT,NEGP),TRL2E(MAXMAT,NEGP)
C  ****  Positron simulation tables.
      COMMON/CPIMFP/SPHEL(MAXMAT,NEGP),SPHIN(MAXMAT,NEGP),
     1  SPISI(MAXMAT,NEGP),SPHBR(MAXMAT,NEGP),SPAN(MAXMAT,NEGP),
     2  SPAUX(MAXMAT,NEGP),SPTOT(MAXMAT,NEGP),CSTPP(MAXMAT,NEGP),
     3  RSTPP(MAXMAT,NEGP),W1P(MAXMAT,NEGP),W2P(MAXMAT,NEGP),
     4  DW1PL(MAXMAT,NEGP),DW2PL(MAXMAT,NEGP),
     5  RNDCP(MAXMAT,NEGP),AP(MAXMAT,NEGP),BP(MAXMAT,NEGP),
     6  T1P(MAXMAT,NEGP),T2P(MAXMAT,NEGP)
      COMMON/CLAS1P/TSTPP(MAXMAT,NEGP),TSTRP(MAXMAT,NEGP),
     1  TRL1P(MAXMAT,NEGP),TRL2P(MAXMAT,NEGP)
C  ****  Elastic scattering of electrons and positrons.
C  ****  Electron and positron radiative yields.
      COMMON/CBRYLD/EBRY(MAXMAT,NEGP),PBRY(MAXMAT,NEGP)
C  ****  Photon simulation tables.
      COMMON/CGIMFP/SGRA(MAXMAT,NEGP),SGCO(MAXMAT,NEGP),
     1  SGPH(MAXMAT,NEGP),SGPP(MAXMAT,NEGP),SGAUX(MAXMAT,NEGP)
      PARAMETER (NDIM=12000)
      COMMON/CGPH01/ER(NDIM),XSR(NDIM),NPHD
      COMMON/CGPP01/TRIP(MAXMAT,NEGP)
	  
	  COMMON/CEBR/WB(NBW),PBCUT(MAXMAT,NEGP),WBCUT(MAXMAT,NEGP),
     1  PDFB(MAXMAT,NEGP,NBW),DPDFB(MAXMAT,NEGP,NBW),
     2  PACB(MAXMAT,NEGP,NBW),ZBR2(MAXMAT)
      COMMON/CEBR02/P0(MAXMAT,NEGP,NBW)
	  
	  COMMON/CBRANG/BET(6),BK(21),BP1(MAXMAT,6,21,4),
     1              BP2(MAXMAT,6,21,4),ZBEQ(MAXMAT)
	 
	  COMMON/CEIN01/EI,EE,CPS,AMOL,MOM
	 
	  COMMON/CSUMGA/IERGA,NCALL
	 
	  COMMON/CPIN01/EI3,CPS3,BHA1,BHA2,BHA3,BHA4,MOM3
	 
      PARAMETER (NE=96,NA=606)
      COMMON/CDCSEP/ETS(NE),ETL(NE),TH(NA),THR(NA),XMU(NA),XMUL(NA),
     1       ECS(NE),ETCS1(NE),ETCS2(NE),EDCS(NE,NA),
     2       PCS(NE),PTCS1(NE),PTCS2(NE),PDCS(NE,NA),
     3       DCSI(NA),DCSIL(NA),CSI,TCS1I,TCS2I
	 
      PARAMETER (NP=128)
      COMMON/CEELDB/XSE(NP,NEGP,MAXMAT),PSE(NP,NEGP,MAXMAT),
     1              ASE(NP,NEGP,MAXMAT),BSE(NP,NEGP,MAXMAT),
     2              ITLE(NP,NEGP,MAXMAT),ITUE(NP,NEGP,MAXMAT)
      COMMON/CPELDB/XSP(NP,NEGP,MAXMAT),PSP(NP,NEGP,MAXMAT),
     1              ASP(NP,NEGP,MAXMAT),BSP(NP,NEGP,MAXMAT),
     2              ITLP(NP,NEGP,MAXMAT),ITUP(NP,NEGP,MAXMAT)
      COMMON/CELSEP/EELMAX(MAXMAT),PELMAX(MAXMAT),
     1              RNDCEd(MAXMAT,NEGP),RNDCPd(MAXMAT,NEGP)
	 
	 
	  COMMON/CGRA00/FACTE,Q2MAX,MM,MOM2
      COMMON/CGRA01/FF(MAXMAT,NQ),ERA(NEX),XSRA(MAXMAT,NEX),
     1    IED(NEGP),IEU(NEGP),NE2
      COMMON/CGRA02/QQ(NQ),AR(MAXMAT,NQ),BR(MAXMAT,NQ),CR(MAXMAT,NQ),
     1              DR(MAXMAT,NQ),FF0(MAXMAT),QQM
      PARAMETER (NP2=150)
      COMMON/CGRA03/QRA(NP2,MAXMAT),PRA(NP2,MAXMAT),DPRA(NP2,MAXMAT),
     1  ARA(NP2,MAXMAT),BRA(NP2,MAXMAT),PMAX(NEGP,MAXMAT),
     2  ITLRA(NP2,MAXMAT),ITURA(NP2,MAXMAT)

      COMMON/CGPP00/ZEQPP(MAXMAT),F0(MAXMAT,2),BCB(MAXMAT)
	  
	  COMMON/CRITA/XT(NM),PAC(NM),DPAC(NM),A1(NM),B1(NM),IL(NM),
     1  IU(NM),NPM11
      COMMON/CRNDG3/XX(NR),A(NR),B(NR),F2(NR),KA(NR),NPM1
	  
	  COMMON/CEBR01/EBT(NBE),XS(NBE,NBW),TXS(NBE),X1(NBE),Y1(NBE)
	  
	  COMMON/CEEL00/EJT(NEGP),XE0(NEGP),XE1(NEGP),XE2(NEGP),XP0(NEGP),
     1  XP1(NEGP),XP2(NEGP),T1E0(NEGP),T2E0(NEGP),T1P0(NEGP),T2P0(NEGP),
     2  EJTL(NEGP),FJL(NEGP),A2(NEGP),B2(NEGP),C(NEGP),D(NEGP)
	  
	  
      INCLUDE 'pmcomms.f'
	  
	  OPEN(IWR,FILE='PENELOPE_mod.txt')
	  write(IWR,*) 'PENELOPE_mod'
	  WRITE(IWR,'(10E14.5)') EABS,C1,C2,WCC,
     1    WCR,DEN,RDEN,E0STEP,DESOFT,SSOFT
	  WRITE(IWR,'(3I5)') NMS,NEGP,NMAT
	  CLOSE(IWR)
	  
	  OPEN(IWR,FILE='PENGEOM_mod.txt')
	  write(IWR,*) 'PENGEOM_mod'
	  WRITE(IWR,'(E14.5)') DSTOT
	  WRITE(IWR,'(4I5)') MATER,KDET,KSLAST,NBODY
	  WRITE(IWR,*) BALIAS
	  WRITE(IWR,*) LVERB
	  CLOSE(IWR)
	  
	  OPEN(IWR,FILE='TRACK_mod.txt')
	  write(IWR,*) 'TRACK_mod'
	  WRITE(IWR,'(12E14.5)') E ,X ,Y,Z,U,V,W,WGHT,
     1    SP1,SP2,SP3,PAGE
	  WRITE(IWR,'(5I5)') KPAR,IBODY,MAT,ILB,IPOL
	  WRITE(IWR,*) LAGE
	  CLOSE(IWR)
	 
	  OPEN(IWR,FILE='QSURF.txt')
	  write(IWR,*) 'QSURF'
	  DO I=1,NS
	    WRITE(IWR,'(I5,10E14.5,2I4)') I,AXX(I),AXY(I),AXZ(I),AYY(I),
     1      AYZ(I),AZZ(I),AX(I),AY(I),AZ(I),A0(I),NSURF,KPLANE(I)
      ENDDO
	  CLOSE(IWR)
	  
	  OPEN(IWR,FILE='QTREE.txt')
	  write(IWR,*) 'QTREE'
	  WRITE(IWR,'(7I5)') NBODYS,KMOTH,KDGHT,KSURF,KFLAG,KSP,NWARN
	  CLOSE(IWR)
	  
	  OPEN(IWR,FILE='QBODY1.txt')
	  write(IWR,*) 'QBODY1'
	  WRITE(IWR,'(10I5)') KBODY, KBOMO
	  CLOSE(IWR)

	  OPEN(IWR,FILE='CECUTR.txt')
	  write(IWR,*) 'ECUTR'
	  WRITE(IWR,'(10E14.5)') ECUTR
	  CLOSE(IWR)
	 
	  OPEN(IWR,FILE='CSGAWR.txt')
	  write(IWR,*) 'CSGAWR'
	  WRITE(IWR,'(I5)') ISGAW
	  CLOSE(IWR)
	  
	  OPEN(IWR,FILE='CERSEC.txt')
	  write(IWR,*) 'CERSEC'
	  WRITE(IWR,'(I5)') IERSEC
	  CLOSE(IWR)
	  
	  OPEN(IWR,FILE='CEGRID.txt')
	  write(IWR,*) 'CEGRID'
	  WRITE(IWR,'(10E14.5)') EMIN,EL,EU,ET1,DLEMP,DLEMP1,
     1  DLFC,XEL,XE,XEK
	  WRITE(IWR,'(I5)') KE
	  CLOSE(IWR)
	  
	  OPEN(IWR,FILE='CESI0.txt')
	  write(IWR,*) 'CESI0'
	  WRITE(IWR,'(E14.5)') XESI
	  WRITE(IWR,'(4I5)') IESIF,IESIL,NSESI,NCURE
	  CLOSE(IWR)
	  
	  OPEN(IWR,FILE='CPSI0.txt')
	  write(IWR,*) 'CPSI0'
	  WRITE(IWR,'(E14.5)') XPSI
	  WRITE(IWR,'(4I5)') IPSIF,IPSIL,NSPSI,NCURP
	  CLOSE(IWR)
	  
	  OPEN(IWR,FILE='CADATA.txt')
	  write(IWR,*) 'CADATA'
	  WRITE(IWR,'(7E14.5)') ATW,EPX,RSCR,ETA,EB,
     1  ALW,CP0
	  WRITE(IWR,'(3I5)') IFI,IKS,NSHT
	  WRITE(IWR,*) LASYMB
	  CLOSE(IWR)
	  
	  OPEN(IWR,FILE='CGPH00.txt')
	  write(IWR,*) 'CGPH00'
	  WRITE(IWR,'(2E14.5)') EPH,XPH
	  WRITE(IWR,'(4I5)') IPHF,IPHL,NPHS,NCUR
	  CLOSE(IWR)
	  
	  OPEN(IWR,FILE='CRELAX.txt')
	  write(IWR,*) 'CRELAX'
	  WRITE(IWR,'(3E14.5)') P,ET,F
	  WRITE(IWR,'(8I5)') IS0,IS1,IS2,IFIRST,ILAST,NCUR,KS,MODER
	  CLOSE(IWR)
	  
	  OPEN(IWR,FILE='CRITAA.txt')
	  write(IWR,*) 'CRITAA'
	  WRITE(IWR,'(4E14.5)') XA,AA,BA,FA
	  WRITE(IWR,'(2I5)') IA,NPM1A
	  CLOSE(IWR)
	  
	  OPEN(IWR,FILE='CRNDG3.txt')
	  write(IWR,*) 'CRNDG3'
	  WRITE(IWR,'(4E14.5)') XX,A,B,F2
	  WRITE(IWR,'(2I5)') KA,NPM1
	  CLOSE(IWR)
	  
	  OPEN(IWR,FILE='CRITA.txt')
	  write(IWR,*) 'CRITA'
	  WRITE(IWR,'(5E14.5)') XT,PAC,DPAC,A1,B1
	  WRITE(IWR,'(3I5)') IL,IU,NPM11
	  CLOSE(IWR)
	  
	  OPEN(IWR,FILE='CRITAN.txt')
	  write(IWR,*) 'CRITAN'
	  WRITE(IWR,'(E14.5)') CNORM
	  CLOSE(IWR)
	 
	  OPEN(IWR,FILE='COMPOS.txt')
	  write(IWR,*) 'COMPOS'
	  WRITE(IWR,'(5E14.5)') STF,ZT,AT,RHO,VMOL
	  WRITE(IWR,'(2I5)') IZ,NELEM
	  CLOSE(IWR)
	  
	  OPEN(IWR,FILE='CRANGE.txt')
	  write(IWR,*) 'CRANGE'
	  WRITE(IWR,'(7E14.5)') RANGE,RANGEL 
	  CLOSE(IWR)
	  
	  OPEN(IWR,FILE='CEIN.txt')
	  write(IWR,*) 'CEIN'
	  WRITE(IWR,'(5E14.5)') EXPOT,OP2,F,UI,WRI
	  WRITE(IWR,'(3I5)') KZ,KS,NOSC
	  CLOSE(IWR)
	 
	  OPEN(IWR,FILE='CEINTF.txt')
	  write(IWR,*) 'CEINTF'
	  WRITE(IWR,'(4E14.5)') T1EI,T2EI,T1PI,T2PI
	  CLOSE(IWR)	 

	  OPEN(IWR,FILE='CEIN00.txt')
	  write(IWR,*) 'CEIN00'
	  WRITE(IWR,'(9E14.5)') SEH0,SEH1,SEH2,SES0,SES1,
     1              SES2,SET0,SET1,SET2
	  CLOSE(IWR)	  

	  OPEN(IWR,FILE='CPIN00.txt')
	  write(IWR,*) 'CPIN00'
	  WRITE(IWR,'(9E14.5)') SPH0,SPH1,SPH2,SPS0,SPS1,
     1              SPS2,SPT0,SPT1,SPT2
	  CLOSE(IWR)	  
	
C  ****  Inner-shell ionisation by electron and positron impact.

      OPEN(IWR,FILE='CEINAC.txt')
	  write(IWR,*) 'CEINAC'
	  WRITE(IWR,'(5E14.5)') EINAC
	  WRITE(IWR,'(4I5)') IEIN,NEIN
	  CLOSE(IWR)
	  
	  OPEN(IWR,FILE='CESIAC.txt')
	  write(IWR,*) 'CESIAC'
	  WRITE(IWR,'(5E14.5)') ESIAC
	  WRITE(IWR,'(4I5)') IESI,NESI
	  CLOSE(IWR)

C  ****  Electron inelastic coll. and inner-shell ionisation tables.
	  OPEN(IWR,FILE='CESIN.txt')
	  write(IWR,*) 'CESIN'
	  WRITE(IWR,'(10E14.5)') XSEIN,XSESI
	  WRITE(IWR,'(10I5)') ISIE
	  CLOSE(IWR)
	  
C  ****  Positron inelastic coll. and inner-shell ionisation tables.
	  OPEN(IWR,FILE='CPINAC.txt')
	  write(IWR,*) 'CPINAC'
	  WRITE(IWR,'(5E14.5)') PINAC 
	  WRITE(IWR,'(4I5)') IPIN,NPIN
	  CLOSE(IWR)
	  
	  OPEN(IWR,FILE='CPSIAC.txt')
	  write(IWR,*) 'CPSIAC'
	  WRITE(IWR,'(5E14.5)') PSIAC 
	  WRITE(IWR,'(4I5)') IPSI,NPSI
	  CLOSE(IWR)
	  
	  OPEN(IWR,FILE='CPSIN.txt')
	  write(IWR,*) 'CPSIN'
	  WRITE(IWR,'(6E14.5)') XSPIN,XSPSI
	  WRITE(IWR,'(4I5)') ISIP
	  CLOSE(IWR)
	  
C  ****  Compton scattering.


	  OPEN(IWR,FILE='CGCO.txt')
	  write(IWR,*) 'CGCO'
	  WRITE(IWR,'(10E14.5)') FCO,UICO,FJ0,PTRSH,OSCCO
	  WRITE(IWR,'(10I5)') KZCO,KSCO
	  CLOSE(IWR)

C  ****  Electron simulation tables.

	  OPEN(IWR,FILE='CEIMFP.txt')
	  write(IWR,*) 'CEIMFP'
	  WRITE(IWR,'(18E14.5)') SEHEL,SEHIN,SEISI,SEHBR,SEAUX,
     1  SETOT,CSTPE,RSTPE,DEL,W1E,W2E,DW1EL,DW2EL,
     2  RNDCE,AE,BE,T1E,T2E
	  CLOSE(IWR)
	  
	  OPEN(IWR,FILE='CLAS1E.txt')
	  write(IWR,*) 'CLAS1E'
	  WRITE(IWR,'(4E14.5)') TSTPE,TSTRE,TRL1E,TRL2E
	  CLOSE(IWR)


C  ****  Positron simulation tables.

	  OPEN(IWR,FILE='CPIMFP.txt')
	  write(IWR,*) 'CPIMFP'
	  WRITE(IWR,'(18E14.5)') SPHEL,SPHIN,SPISI,SPHBR,SPAN,
     1  SPAUX,SPTOT,CSTPP,RSTPP,W1P,W2P,DW1PL,DW2PL,
     2  RNDCP,AP,BP,T1P,T2P
	  CLOSE(IWR)
	  
	  OPEN(IWR,FILE='CLAS1P.txt')
	  write(IWR,*) 'CLAS1P'
	  WRITE(IWR,'(4E14.5)') TSTPP,TSTRP,TRL1P,TRL2P
	  CLOSE(IWR)

C  ****  Elastic scattering of electrons and positrons.
 
	  OPEN(IWR,FILE='CEEL00.txt')
	  write(IWR,*) 'CEEL00'
	  WRITE(IWR,'(17E14.5)') EJT,XE0,XE1,XE2,XP0,
     1  XP1,XP2,T1E0,T2E0,T1P0,T2P0,
     2  EJTL,FJL,A2,B2,C,D
	  CLOSE(IWR)
	  

	  

C  ****  Electron and positron radiative yields.
      OPEN(IWR,FILE='CBRYLD.txt')
	  write(IWR,*) 'CBRYLD'
	  WRITE(IWR,'(4E14.5)') EBRY,PBRY
	  CLOSE(IWR)
 
C  ****  Photon simulation tables.
 
      OPEN(IWR,FILE='CGIMFP.txt')
	  write(IWR,*) 'CGIMFP'
	  WRITE(IWR,'(10E14.5)') SGRA,SGCO,SGPH,SGPP,SGAUX
	  CLOSE(IWR)
	  
	  OPEN(IWR,FILE='CGPH01.txt')
	  write(IWR,*) 'CGPH01'
	  WRITE(IWR,'(8E14.5)') ER,XSR
	  WRITE(IWR,'(10I5)') NPHD
	  CLOSE(IWR)
	  
	  OPEN(IWR,FILE='CGPP01.txt')
	  write(IWR,*) 'CGPP01'
	  WRITE(IWR,'(8E14.5)') TRIP
	  CLOSE(IWR)
	  
	  OPEN(IWR,FILE='CEBR.txt')
	  write(IWR,*) 'CEBR'
	  WRITE(IWR,'(7E14.5)') WB,PBCUT,WBCUT,
     1  PDFB,DPDFB,PACB,ZBR2
	  CLOSE(IWR)
	  
	  OPEN(IWR,FILE='CEBR01.txt')
	  write(IWR,*) 'CEBR01'
	  WRITE(IWR,'(7E14.5)') EBT,XS,TXS,X1,Y1
	  CLOSE(IWR)
	  
	  OPEN(IWR,FILE='CEBR02.txt')
	  write(IWR,*) 'CEBR02'
	  WRITE(IWR,'(7E14.5)') P0
	  CLOSE(IWR)
	  
	  OPEN(IWR,FILE='CBRANG.txt')
	  write(IWR,*) 'CBRANG'
	  WRITE(IWR,'(8E14.5)') BET,BK,BP1,BP2,ZBEQ
	  CLOSE(IWR)
	  
	  OPEN(IWR,FILE='CEIN01.txt')
	  write(IWR,*) 'CEIN01'
	  WRITE(IWR,'(8E14.5)') EI,EE,CPS,AMOL
	  WRITE(IWR,'(10I5)') MOM 
	  CLOSE(IWR)
	  
	  OPEN(IWR,FILE='CSUMGA.txt')
	  write(IWR,*) 'CSUMGA'
	  WRITE(IWR,'(10I5)') IERGA,NCALL
	  CLOSE(IWR)
	  
	  OPEN(IWR,FILE='CPIN01.txt')
	  write(IWR,*) 'CPIN01'
	  WRITE(IWR,'(8E14.5)') EI3,CPS3,BHA1,BHA2,BHA3,BHA4
	  WRITE(IWR,'(10I5)') MOM3
	  CLOSE(IWR)
	  
	  
	  OPEN(IWR,FILE='CDCSEP.txt')
	  write(IWR,*) 'CDCSEP'
	  WRITE(IWR,'(19E14.5)') ETS,ETL,TH,THR,XMU,XMUL,
     1       ECS,ETCS1,ETCS2,EDCS,
     2       PCS,PTCS1,PTCS2,PDCS,
     3       DCSI,DCSIL,CSI,TCS1I,TCS2I
	  CLOSE(IWR)
	  
	  OPEN(IWR,FILE='CEELDB.txt')
	  write(IWR,*) 'CEELDB'
	  WRITE(IWR,'(12E14.5)') XSE,PSE,ASE,BSE
	  WRITE(IWR,'(10I5)') ITLE,ITUE
	  CLOSE(IWR)
	  
	  OPEN(IWR,FILE='CPELDB.txt')
	  write(IWR,*) 'CPELDB'
	  WRITE(IWR,'(12E14.5)') XSP,PSP,ASP,BSP
	  WRITE(IWR,'(10I5)') ITLP,ITUP
	  CLOSE(IWR)
	  
	  OPEN(IWR,FILE='CELSEP.txt')
	  write(IWR,*) 'CELSEP'
	  WRITE(IWR,'(8E14.5)') EELMAX,PELMAX,RNDCEd,RNDCPd
	  CLOSE(IWR)
	  
	  OPEN(IWR,FILE='CGRA00.txt')
	  write(IWR,*) 'CGRA00'
	  WRITE(IWR,'(4E14.5)') FACTE,Q2MAX
	  WRITE(IWR,'(10I5)') MM,MOM2
	  CLOSE(IWR)
	  
	  OPEN(IWR,FILE='CGRA01.txt')
	  write(IWR,*) 'CGRA01'
	  WRITE(IWR,'(6E14.5)') FF,ERA,XSRA
	  WRITE(IWR,'(6I5)') IED,IEU,NE2
	  CLOSE(IWR)
	  
	  OPEN(IWR,FILE='CGRA02.txt')
	  write(IWR,*) 'CGRA02'
	  WRITE(IWR,'(7E14.5)') QQ,AR,BR,CR,DR,FF0,QQM
	  CLOSE(IWR)
	  
	  OPEN(IWR,FILE='CGRA03.txt')
	  write(IWR,*) 'CGRA03'
	  WRITE(IWR,'(7E14.5)') QRA,PRA,DPRA,ARA,BRA,PMAX
	  WRITE(IWR,'(6I5)') ITLRA,ITURA
	  CLOSE(IWR)
	  
	  OPEN(IWR,FILE='CGPP00.txt')
	  write(IWR,*) 'CGPP00'
	  WRITE(IWR,'(7E14.5)') ZEQPP,F0,BCB
	  CLOSE(IWR)
	 

      RETURN
      END
	  
	  
	  
