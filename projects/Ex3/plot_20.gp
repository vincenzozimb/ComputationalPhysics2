set datafile separator '\t'
set term qt 0

set xrange [0.0:15]

plot 's.csv' using 1:2 w lines title "0s" ,\
     's.csv' using 1:3 w lines title "1s" ,\
     'p.csv' using 1:2 w lines title "0p" ,\
     'd.csv' using 1:2 w lines title "0d"

pause -1
    
plot 'density.csv' using 1:2 w lines title "electron density"

pause -1