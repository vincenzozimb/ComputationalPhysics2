set datafile separator '\t'
set term qt 0

set title "Interacting electronic densities"
plot "density8.csv" using 1:2 w lines title "N = 8" ,\
    "density20.csv" using 1:2 w lines title "N = 20",\
    "density40.csv" using 1:2 w lines title "N = 40"
pause -1
