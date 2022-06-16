set datafile separator '\t'

set term qt 0
set xrange  [0:4]


#set title 'Eigenfunctions 3d'
plot 'psi3d.csv' using 1:2 w lines notitle
pause -1