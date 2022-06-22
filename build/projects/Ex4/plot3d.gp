set datafile separator '\t'
set term qt 0

plot "dataMC3d.csv" using 1:2  t "Monte Carlo", "dataEx3d.csv" using 1:2 w l t "Exact"
set xlabel "{/Symbol a}" font ",15"
set ylabel "E_{tot}[{/Symbol a}]" font ",15"
pause -1
