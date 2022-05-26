set datafile separator '\t'
set term qt 0

# set xrange [0:4]

plot 'potential.csv' using 1:2 w lines notitle ,\
    'wf.csv' using 1:2 w lines notitle 

pause -1