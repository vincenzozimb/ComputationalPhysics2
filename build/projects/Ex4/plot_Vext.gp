set datafile separator '\t'
set term qt 0

set title "External potential V_{ext}"
plot "V_ext.csv" using 1:2 w lines notitle
pause -1
