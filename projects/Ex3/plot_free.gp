set datafile separator '\t'
set term qt 0

plot "density_free.csv" using 1:2 w lines title "free electronic density for N = 8"
pause -1

plot "density.csv" using 1:2 w lines title "interacting electronic density of N = 8"
pause -1