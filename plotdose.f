C
C     PROGRAM PLOTDOSE
C
C  This program reads the '3d-dose-map.dat' output file from PENMAIN and
C  generates a map of the absorbed dose on a plane perpendicular to
C  one of the coordinate axes. The program generates the GNUPLOT
C  scripts 'plotdose.gnu' (which displays the dose map and a twin image
C  of the materials in the view plane), and 'plotdose2.gnu' (which
C  displays the dose map as a two-dimensional surface).
C
C  The executable binary file 'plotdose.exe' must be placed in the same
C  directory as the files '3d-dose-map.dat' and 'geometry.rep' generated
C  by PENMAIN.
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


      INCLUDE 'pengeom.f'

C  *********************************************************************
C                       MAIN PROGRAM
C  *********************************************************************
C
      USE TRACK_mod
      USE PENGEOM_mod
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
      CHARACTER LINE*100,LAX*1,NUMA*10,NUMB*10
      DIMENSION DXL(3),DXU(3),NDB(3)
      DIMENSION BDOSE(3),RBDOSE(3)
      PARAMETER (NP=250)
      DIMENSION DOSE(NP,NP,NP),EDOSE(NP,NP,NP)
      DIMENSION PARINP(20)
C
      DO I=1,NP
      DO J=1,NP
      DO K=1,NP
        DOSE(I,J,K)=0.0D0
        EDOSE(I,J,K)=0.0D0
      ENDDO
      ENDDO
      ENDDO
C
      OPEN(15,FILE='geometry.rep')
      OPEN(16,FILE='geometry2.rep')
      NPINP=0
      CALL GEOMIN(PARINP,NPINP,NMATG,NBODY,15,16)
      CLOSE(15)
      CLOSE(16)
C
      OPEN(9,FILE='3d-dose-map.dat')
      OPEN(19,FILE='plotdose.dat')
      READ(9,'(A)') LINE
      WRITE(19,9420)
 9420 FORMAT(1X,'#  Results from PENMAIN. 3D dose distribution.')
      READ(9,'(24X,E13.6,11X,E13.6)') DXL(1),DXU(1)
      WRITE(19,9421) DXL(1),DXU(1)
 9421 FORMAT(1X,'#  Dose-map box:  XL = ',1P,E13.6,
     1  ' cm,  XU = ',E13.6,' cm')
      READ(9,'(24X,E13.6,11X,E13.6)') DXL(2),DXU(2)
      WRITE(19,9422) DXL(2),DXU(2)
 9422 FORMAT(1X,'#',17X,'YL = ',1P,E13.6,' cm,  YU = ',E13.6,' cm')
      READ(9,'(24X,E13.6,11X,E13.6)') DXL(3),DXU(3)
      WRITE(19,9423) DXL(3),DXU(3)
 9423 FORMAT(1X,'#',17X,'ZL = ',1P,E13.6,' cm,  ZU = ',E13.6,' cm')
      READ(9,'(30X,I4,7X,I4,7X,I4)')  NDB(1),NDB(2),NDB(3)
      WRITE(19,9424) NDB(1),NDB(2),NDB(3)
 9424 FORMAT(1X,'#  Numbers of bins:     NBX =',I4,', NBY =',I4,
     1      ', NBZ =',I4,/1X,'#')
      READ(9,'(A)') LINE
      READ(9,'(A)') LINE
      READ(9,'(A)') LINE
      READ(9,'(A)') LINE
      READ(9,'(A)') LINE
C
      FSAFE=1.000000001D0
      DO I=1,3
        BDOSE(I)=FSAFE*(DXU(I)-DXL(I))/DBLE(NDB(I))
        RBDOSE(I)=1.0D0/BDOSE(I)
      ENDDO
C
      DO I3=1,NDB(3)
        DO I1=1,NDB(1)
          DO I2=1,NDB(2)
            READ(9,*) X,Y,Z,DOSE(I1,I2,I3),EDOSE(I1,I2,I3)
          ENDDO
          READ(9,'(A)') LINE
        ENDDO
        READ(9,'(A)') LINE
      ENDDO
C
 1    CONTINUE
      WRITE(6,*) 'Select the axis (X, Y or Z) that is perpendicular ',
     1  'to the plane...'
      READ(5,'(A)') LAX
C
      IF(LAX.EQ.'x'.OR.LAX.EQ.'X') THEN
        WRITE(6,*) 'Enter the X coordinate of the plane...'
        READ(5,*) PCOOR
        WRITE(19,9431)
 9431   FORMAT(1X,'#  Columns 1 to 3: X,Y,Z coordinates of the',
     1    ' bin centres (cm).',/1X,'#  4th column: dose (eV/g).',
     1    /1X,'#  Columns 5 to 7, voxel indices IX,IY,IZ.',
     1    /1X,'#  Column 8, material MAT.',/)
        ICASO=1
        U=1.0D0
        V=0.0D0
        W=0.0D0
        X=PCOOR
C        I1=1+(PCOOR-DXL(1))*RBDOSE(1)
		DO I1=1,NDB(1)+1
          X=DXL(1)+(I1-0.5D0)*BDOSE(1)
          XW=X-0.5D0*BDOSE(1)
          X=MIN(MAX(X,DXL(1)+1.0D-9),DXU(1)-1.0D-9)
			DO I2=1,NDB(2)+1
				Y=DXL(2)+(I2-0.5D0)*BDOSE(2)
				YW=Y-0.5D0*BDOSE(2)
				Y=MIN(MAX(Y,DXL(2)+1.0D-9),DXU(2)-1.0D-9)
				DO I3=1,NDB(3)+1
					Z=DXL(3)+(I3-0.5D0)*BDOSE(3)
					ZW=Z-0.5D0*BDOSE(3)
					Z=MIN(MAX(Z,DXL(3)+1.0D-9),DXU(3)-1.0D-9)
					CALL LOCATE
					WRITE(19,9440) XW,YW,ZW,DOSE(I1,I2,I3),I1,I2,I3,IBODY
				ENDDO
				WRITE(19,*) '   '
			ENDDO	
			WRITE(19,*) '   '
        ENDDO
      ELSE IF(LAX.EQ.'y'.OR.LAX.EQ.'Y') THEN
        WRITE(6,*) 'Enter the Y coordinate of the plane...'
        READ(5,*) PCOOR
        WRITE(19,9432)
 9432   FORMAT(1X,'#  Columns 1 to 3: Y,X,Z coordinates of the',
     1    ' bin centres (cm).',/1X,'#  4th column: dose (eV/g).',
     1    /1X,'#  Columns 5 to 7, voxel indices IX,IY,IZ.',
     1    /1X,'#  Column 8, material MAT.',/)
        ICASO=2
        U=0.0D0
        V=1.0D0
        W=0.0D0
        Y=PCOOR
        I2=1+(PCOOR-DXL(2))*RBDOSE(2)
        DO I1=1,NDB(1)+1
          X=DXL(1)+(I1-0.5D0)*BDOSE(1)
          XW=X-0.5D0*BDOSE(1)
          X=MIN(MAX(X,DXL(1)+1.0D-9),DXU(1)-1.0D-9)
          DO I3=1,NDB(3)+1
            Z=DXL(3)+(I3-0.5D0)*BDOSE(3)
            ZW=Z-0.5D0*BDOSE(3)
            Z=MIN(MAX(Z,DXL(3)+1.0D-9),DXU(3)-1.0D-9)
            CALL LOCATE
            WRITE(19,9440) PCOOR,XW,ZW,DOSE(I1,I2,I3),I1,I2,I3,IBODY
          ENDDO
          WRITE(19,*) '   '
        ENDDO
      ELSE IF(LAX.EQ.'z'.OR.LAX.EQ.'Z') THEN
        WRITE(6,*) 'Enter the Z coordinate of the plane...'
        READ(5,*) PCOOR
        WRITE(19,9433)
 9433   FORMAT(1X,'#  Columns 1 to 3: Z,X,Y coordinates of the',
     1    ' bin centres (cm).',/1X,'#  4th column: dose (eV/g).',
     1    /1X,'#  Columns 5 to 7, voxel indices IX,IY,IZ.',
     1    /1X,'#  Column 8, material MAT.',/)
        ICASO=3
        U=0.0D0
        V=0.0D0
        W=1.0D0
        Z=PCOOR
        I3=1+(PCOOR-DXL(3))*RBDOSE(3)
        DO I1=1,NDB(1)+1
          X=DXL(1)+(I1-0.5D0)*BDOSE(1)
          XW=X-0.5D0*BDOSE(1)
          X=MIN(MAX(X,DXL(1)+1.0D-9),DXU(1)-1.0D-9)
          DO I2=1,NDB(2)+1
            Y=DXL(2)+(I2-0.5D0)*BDOSE(2)
            YW=Y-0.5D0*BDOSE(2)
            Y=MIN(MAX(Y,DXL(2)+1.0D-9),DXU(2)-1.0D-9)
            CALL LOCATE
            WRITE(19,9440) PCOOR,XW,YW,DOSE(I1,I2,I3),I1,I2,I3,IBODY
          ENDDO
          WRITE(19,*) '   '
        ENDDO
      ELSE
        GO TO 1
      ENDIF
 9440 FORMAT(1P,4E11.3,4I4)
      CLOSE(9)
      CLOSE(19)
C
C  ****  GNUPLOT script for displaying the material and dose maps.
C
      OPEN(9,FILE='plotdose.gnu')
      WRITE(9,'(A)') '# Gnuplot MS-Windows 32 bit version 4.2'
      WRITE(9,'(A)') '# Plots results from ''plotdose.f'' '
      WRITE(9,'(A)') '  '
      WRITE(9,'(A)') 'set zero 1.0e-60; unset mouse'
      WRITE(9,'(A)') 'set style line 100 lt 5 lw 0.25'
      WRITE(9,'(A)') '# set pm3d solid hidden3d 100'
      WRITE(9,'(A)') 'set pm3d corners2color c1 map'
      WRITE(9,'(A)') 'set xlabel offset 0,1'
      WRITE(9,'(A)') 'set ylabel offset .25,0'
      WRITE(9,'(A)') '  '
C
      IF(ICASO.EQ.1) THEN
        WRITE(NUMA,'(1P,E10.3)') PCOOR
        WRITE(9,'(3A)') 'set label ''Plane X =',NUMA,
     1    ' cm'' at screen 0.5,0.92 center'
        WRITE(9,'(A)') 'set xlabel ''Y (cm)'' offset 0,+0.75'
        WRITE(9,'(A)') 'set ylabel ''Z (cm)'' offset -1,0'
        WRITE(NUMA,'(1P,E10.3)') DXL(2)
        WRITE(NUMB,'(1P,E10.3)') DXU(2)
        WRITE(9,'(5A)') 'set xrange [',NUMA,':',NUMB,']'
        WRITE(NUMA,'(1P,E10.3)') DXL(3)
        WRITE(NUMB,'(1P,E10.3)') DXU(3)
        WRITE(9,'(5A)') 'set yrange [',NUMA,':',NUMB,']'
        WRITE(9,'(A)') 'set zrange [0.:*]'
        WRITE(9,'(A)') 'set cbrange [0.:*]'
      ELSE IF(ICASO.EQ.2) THEN
        WRITE(NUMA,'(1P,E10.3)') PCOOR
        WRITE(9,'(3A)') 'set label ''Plane Y =',NUMA,
     1    ' cm'' at screen 0.5,0.92 center'
        WRITE(9,'(A)') 'set xlabel ''X (cm)'' offset 0,+0.75'
        WRITE(9,'(A)') 'set ylabel ''Z (cm)'' offset -1,0'
        WRITE(NUMA,'(1P,E10.3)') DXL(1)
        WRITE(NUMB,'(1P,E10.3)') DXU(1)
        WRITE(9,'(5A)') 'set xrange [',NUMA,':',NUMB,']'
        WRITE(NUMA,'(1P,E10.3)') DXL(3)
        WRITE(NUMB,'(1P,E10.3)') DXU(3)
        WRITE(9,'(5A)') 'set yrange [',NUMA,':',NUMB,']'
        WRITE(9,'(A)') 'set zrange [0.:*]'
        WRITE(9,'(A)') 'set cbrange [0.:*]'
      ELSE
        WRITE(NUMA,'(1P,E10.3)') PCOOR
        WRITE(9,'(3A)') 'set label ''Plane Z =',NUMA,
     1    ' cm'' at screen 0.5,0.92 center'
        WRITE(9,'(A)') 'set xlabel ''X (cm)'' offset 0,+0.75'
        WRITE(9,'(A)') 'set ylabel ''Y (cm)'' offset -1,0'
        WRITE(NUMA,'(1P,E10.3)') DXL(1)
        WRITE(NUMB,'(1P,E10.3)') DXU(1)
        WRITE(9,'(5A)') 'set xrange [',NUMA,':',NUMB,']'
        WRITE(NUMA,'(1P,E10.3)') DXL(2)
        WRITE(NUMB,'(1P,E10.3)') DXU(2)
        WRITE(9,'(5A)') 'set yrange [',NUMA,':',NUMB,']'
        WRITE(9,'(A)') 'set zrange [0.:*]'
        WRITE(9,'(A)') 'set cbrange [0.:*]'
      ENDIF
C
      WRITE(9,'(A)') '  '
      WRITE(9,'(2A)') '# To generate an eps file with the plot,',
     1  ' uncomment the following line.'
      WRITE(9,'(2A)') '# set terminal postscript eps color; ',
     1  'set output ''mat2d.eps'' '
      WRITE(9,'(A)') 'set title ''Materials'' '
      WRITE(9,'(A)') 'set size ratio -1 '
      WRITE(9,'(A)') 'set xtics offset 0,+0.5'
      WRITE(9,'(A)') 'set ytics offset +0.25,0'
      WRITE(9,'(2A)') 'set cbtics offset -0.5,0 autofreq 1 ',
     1  'format ''%.0f'' '
      WRITE(9,'(A)') 'splot ''plotdose.dat'' u 2:3:8 notitle w pm3d'
      WRITE(9,'(A)') 'pause -1 ''Press enter to continue'' '
C
      WRITE(9,'(A)') '  '
      WRITE(9,'(2A)') '# To generate an eps file with the plot,',
     1  ' uncomment the following line.'
      WRITE(9,'(2A)') '# set terminal postscript eps color; ',
     1  'set output ''dose2d.eps'' '
      WRITE(9,'(A)') 'set title ''Log10(1+Dose/(eV/g))'' '
      WRITE(9,'(A)') 'set size ratio -1 '
      WRITE(9,'(A)') 'set xtics offset 0,+0.5'
      WRITE(9,'(A)') 'set ytics offset +0.25,0'
      WRITE(9,'(2A)') 'set cbtics offset -0.25,0 format ''%.1te%+-3T'' '
      WRITE(9,'(2A)') 'splot ''plotdose.dat'' u 2:3:(log10(1.0+$4))',
     1  ' notitle w pm3d'
      WRITE(9,'(A)') 'pause -1 ''Press enter to continue'' '
C
      WRITE(9,'(A)') '  '
      WRITE(9,'(2A)') '# To generate an eps file with the plot,',
     1  ' uncomment the following line.'
      WRITE(9,'(2A)') '# set terminal postscript eps color; ',
     1  'set output ''matdose2d.eps'' '
      WRITE(9,'(A)') 'set multiplot'
      WRITE(9,'(A)') '  '
      WRITE(9,'(A)') 'set size ratio -1 0.5,1'
      WRITE(9,'(A)') 'set origin 0,0'
      WRITE(9,'(A)') 'set title ''Materials'' '
      WRITE(9,'(A)') 'set xtics offset 0,+0.5'
      WRITE(9,'(A)') 'set ytics offset +0.25,0'
      WRITE(9,'(2A)') 'set cbtics offset -0.5,0 autofreq 1 ',
     1  'format ''%.0f'' '
      WRITE(9,'(A)') 'splot ''plotdose.dat'' u 2:3:8 notitle w pm3d'
      WRITE(9,'(A)') '  '
      WRITE(9,'(A)') 'set size ratio -1 0.5,1'
      WRITE(9,'(A)') 'set origin 0.47,0'
      WRITE(9,'(A)') 'unset ylabel'
      WRITE(9,'(A)') 'set title ''Log10(1+Dose/(eV/g))'' '
      WRITE(9,'(A)') 'set xtics offset 0,+0.5'
      WRITE(9,'(A)') 'set ytics offset +0.25,0'
      WRITE(9,'(2A)') 'set cbtics offset -0.25,0 format ''%.1te%+-3T'' '
      WRITE(9,'(2A)') 'splot ''plotdose.dat'' u 2:3:(log10(1.0+$4))',
     1  ' notitle w pm3d'
      WRITE(9,'(A)') '  '
      WRITE(9,'(A)') 'unset multiplot'
      WRITE(9,'(A)') '  '
      WRITE(9,'(A)') 'pause -1 ''Press enter to continue'' '
      CLOSE(9)
C
C  ****  GNUPLOT script for displaying the dose map as a surface.
C
      OPEN(9,FILE='plotdose2.gnu')
      WRITE(9,'(A)') '# Gnuplot MS-Windows 32 bit version 4.2'
      WRITE(9,'(A)') '# Plots results from ''plotdose.f'' '
      WRITE(9,'(A)') '  '
      WRITE(9,'(A)') '# White-to-black color palette'
      WRITE(9,'(A)') 'set palette rgbformulae -30,-31,-32'
      WRITE(9,'(A)') 'set style line 10 linecolor rgb ''black'' lw 1'
      WRITE(9,'(A)') 'set pm3d hidden3d 10'
      WRITE(9,'(A)') '  '
      WRITE(9,'(A)') 'set zero 1.0e-60'
      WRITE(9,'(A)') 'set ticslevel 0.05'
      WRITE(9,'(A)') 'unset colorbox'
      WRITE(9,'(A)') 'set format z ''%.1te%+-3T'' '
      WRITE(9,'(A)') 'set xtics offset 0,0'
      WRITE(9,'(A)') 'set ytics offset 0,-0.3'
      WRITE(9,'(A)') 'set ztics offset 0.5,0'
      WRITE(9,'(A)') '  '
C
      IF(ICASO.EQ.1) THEN
        WRITE(NUMA,'(1P,E10.3)') PCOOR
        WRITE(9,'(3A)')
     1   'set title ''2D dose distribution (eV/g). Plane X =',
     1   NUMA,' cm'' '
        WRITE(9,'(A)') 'set xlabel ''Y (cm)'' offset 0,0'
        WRITE(9,'(A)') 'set ylabel ''Z (cm)'' offset 0,0'
        WRITE(NUMA,'(1P,E10.3)') DXL(2)
        WRITE(NUMB,'(1P,E10.3)') DXU(2)
        WRITE(9,'(5A)') 'set xrange [',NUMA,':',NUMB,']'
        WRITE(NUMA,'(1P,E10.3)') DXL(3)
        WRITE(NUMB,'(1P,E10.3)') DXU(3)
        WRITE(9,'(5A)') 'set yrange [',NUMA,':',NUMB,']'
        WRITE(9,'(A)') 'set zrange [0.:*]'
      ELSE IF(ICASO.EQ.2) THEN
        WRITE(NUMA,'(1P,E10.3)') PCOOR
        WRITE(9,'(3A)')
     1   'set title ''2D dose distribution (eV/g). Plane Y =',
     1   NUMA,' cm'' '
        WRITE(9,'(A)') 'set xlabel ''X (cm)'' offset 0,0'
        WRITE(9,'(A)') 'set ylabel ''Z (cm)'' offset 0,0'
        WRITE(NUMA,'(1P,E10.3)') DXL(1)
        WRITE(NUMB,'(1P,E10.3)') DXU(1)
        WRITE(9,'(5A)') 'set xrange [',NUMA,':',NUMB,']'
        WRITE(NUMA,'(1P,E10.3)') DXL(3)
        WRITE(NUMB,'(1P,E10.3)') DXU(3)
        WRITE(9,'(5A)') 'set yrange [',NUMA,':',NUMB,']'
        WRITE(9,'(A)') 'set zrange [0.:*]'
      ELSE
        WRITE(NUMA,'(1P,E10.3)') PCOOR
        WRITE(9,'(3A)')
     1   'set title ''2D dose distribution (eV/g). Plane Z =',
     1   NUMA,' cm'' '
        WRITE(9,'(A)') 'set xlabel ''X (cm)'' offset 0,0'
        WRITE(9,'(A)') 'set ylabel ''Y (cm)'' offset 0,0'
        WRITE(NUMA,'(1P,E10.3)') DXL(1)
        WRITE(NUMB,'(1P,E10.3)') DXU(1)
        WRITE(9,'(5A)') 'set xrange [',NUMA,':',NUMB,']'
        WRITE(NUMA,'(1P,E10.3)') DXL(2)
        WRITE(NUMB,'(1P,E10.3)') DXU(2)
        WRITE(9,'(5A)') 'set yrange [',NUMA,':',NUMB,']'
        WRITE(9,'(A)') 'set zrange [0.:*]'
      ENDIF
C
      WRITE(9,'(A)') '  '
      WRITE(9,'(2A)') '# To generate an eps file with the plot,',
     1  ' uncomment the following line.'
      WRITE(9,'(2A,/)') '# set terminal postscript eps color; ',
     1  'set output ''dose2dsurf.eps'' '
      WRITE(9,'(A)') 'set size ratio -1 '
      WRITE(9,'(A)') 'set view 65,120'
      WRITE(9,'(A)') 'splot ''plotdose.dat'' u 2:3:4 notitle w pm3d'
      WRITE(9,'(A)') 'pause -1 ''Press enter to continue'' '
C
      WRITE(9,'(A)') '  '
      WRITE(NUMA,'(1P,E10.3)') PCOOR
      IF(ICASO.EQ.1) THEN
        WRITE(9,'(4A)')
     1   'set title ''2D dose distribution, Log10(1+Dose/(eV/g)).',
     1   ' Plane X =', NUMA,' cm'' '
      ELSE IF(ICASO.EQ.2) THEN
        WRITE(9,'(4A)')
     1   'set title ''2D dose distribution, Log10(1+Dose/(eV/g)).',
     1   ' Plane Y =', NUMA,' cm'' '
      ELSE
        WRITE(9,'(4A)')
     1   'set title ''2D dose distribution, Log10(1+Dose/(eV/g)).',
     1   ' Plane Z =', NUMA,' cm'' '
      ENDIF
      WRITE(9,'(A)')
     1  'splot ''plotdose.dat'' u 2:3:(log10(1.0+$4)) notitle w pm3d'
      WRITE(9,'(A)') 'pause -1 ''Press enter to continue'' '
C
      CLOSE(9)
      END
