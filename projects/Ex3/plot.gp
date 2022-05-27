set datafile separator '\t'
set term qt 0

set xrange [0.8:30]

plot 'wf.csv' using 1:2 w lines notitle #,\
    #'potential.csv' using 1:2 w lines notitle 

pause -1