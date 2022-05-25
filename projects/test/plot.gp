set datafile separator '\t'
set term qt 0

# set xrange [-3:3]


plot 'potential.csv' using 1:2 w lines ,\
    'wf.csv' using 1:2 w lines notitle
pause -1