set datafile separator '\t'
set term qt 0

set xrange [0.0:15]

plot 's.csv' using 1:2 w lines title "0s" ,\
     'p.csv' using 1:2 w lines title "0p"

pause -1
    
plot 'density.csv' using 1:2 w lines title "electron density"

pause -1