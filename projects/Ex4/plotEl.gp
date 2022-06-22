set datafile separator '\t'
set term qt 0

plot "dataVext.csv" u 1:2 w l
set xlabel "r" font ",15"
set ylabel "V_{ext}(r)" font ",15"

pause -1
