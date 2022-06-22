set datafile separator '\t'
set term qt 0

plot "dataMC.csv" using 1:2 w l t "Monte Carlo", "dataEx.csv" using 1:2 w l t "Exact"
set xlabel "\alpha"
set ylabel "E_{tot}[\alpha]"
set title "Total energy"
pause -1
