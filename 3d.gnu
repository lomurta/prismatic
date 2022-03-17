set palette defined ( 0 '#000000',\
                      1 '#000fff',\
                      2 '#0090ff',\
                      3 '#0fffee',\
                      4 '#90ff70',\
                      5 '#ffee00',\
                      6 '#ff7000',\
                      7 '#ee0000',\
                      8 '#7f0000')
                      
                      
set xlabel 'X (cm)'
set ylabel 'Y (cm)' 
set zlabel 'Z (cm)' 
set xrange [-2.500E+01: 2.500E+01]
set yrange [-1.750E+01: 1.750E+01]
set zrange [-1.000E+01: 1.000E+01]
set cbrange [0.:*]


set xtics 3
set ytics 3
set ztics 3


set view equal xyz
set ticslevel 0.0
set pm3d
splot "plotdose.dat" using 1:2:3:((($8 == 1 || $8 == 9 || $8 == 10 || $8 == 11 || $8 == 12) )? log10(1.0 +$4): NaN) notitle w pm3d

pause -1 'Press para continuar'
