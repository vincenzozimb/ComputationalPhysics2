set datafile separator '\t'
set term qt 0


plot "data.csv" using 1:2 w lines notitle
pause -1

set title "Free electronic density for N = 20"
plot "density_free.csv" using 1:2 w lines notitle
pause -1

set title "Interacting electronic density of N = 20"
plot "density.csv" using 1:2 w lines notitle
pause -1
