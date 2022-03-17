set palette defined ( 0 '#000000',\
                      1 '#000fff',\
                      2 '#0090ff',\
                      3 '#0fffee',\
                      4 '#90ff70',\
                      5 '#ffee00',\
                      6 '#ff7000',\
                      7 '#ee0000',\
                      8 '#7f0000')
                      
                      
set xlabel 'Y (cm)'
set ylabel 'Z (cm)' 

set xrange [-1.750E+01: 1.750E+01]
set yrange [-1.000E+01: 1.000E+01]
set cbrange [0.:*]

set xtics 2
set ytics 2

set view equal xy
set ticslevel 0.0
set view map

splot "plotdose.dat" u 2:3:(($5 == 26)? (log10(1.0 +$4)) : NaN) notitle w pm3d

pause -1 'Press para continuar'
