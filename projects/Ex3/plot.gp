set datafile separator '\t'
set term qt 0

set xrange [0:20]
set title "Free electronic densities" font "Helvetica,16"
plot "density_free8.csv" using 1:2 w lines title "N = 8" ,\
    "density_free20.csv" using 1:2 w lines title "N = 20",\
    "density_free40.csv" using 1:2 w lines title "N = 40"
pause -1

set xrange [0:20]
set title "Interacting electronic densities" font "Helvetica,16"
plot "density8.csv" using 1:2 w lines title "N = 8" ,\
    "density20.csv" using 1:2 w lines title "N = 20",\
    "density40.csv" using 1:2 w lines title "N = 40"
pause -1
