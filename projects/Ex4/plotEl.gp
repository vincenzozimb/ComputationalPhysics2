set datafile separator '\t'
set term qt 0

set grid
plot "dataEnEl.csv" u 1:2 w l
set xlabel "eta" font ",15"
set ylabel "E_{tot} [eta]" font ",15"

pause -1
