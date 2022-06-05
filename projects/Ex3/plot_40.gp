set datafile separator '\t'
set term qt 0

# set title "Free electronic density for N = 40"
# plot "density_free.csv" using 1:2 w lines notitle
# pause -1

set title "Interacting electronic density of N = 40"
plot "density40.csv" using 1:2 w lines notitle
pause -1