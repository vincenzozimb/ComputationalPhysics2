set datafile separator '\t'

set term qt 0

plot 'data.csv' using 1:2 w lines notitle
pause -1