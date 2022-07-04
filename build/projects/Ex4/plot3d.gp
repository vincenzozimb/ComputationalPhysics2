set datafile separator '\t'
set term qt 0

set grid
plot "dataMC3d.csv" using 1:2  t "Monte Carlo", "dataEx3d.csv" using 1:2 w l t "Exact"
set xlabel "{/Symbol a}" font ",15"
set ylabel "E_{tot}[{/Symbol a}]" font ",15"
pause -1

set grid
plot "dataMC3d.csv" u 1:3 t "Monte Carlo" , "dataEx3d.csv" using 1:3 w l t "Exact"
#set xlabel "{/Symbol a}" font ",15"
#set ylabel "Variance[{/Symbol a}]" font ",15"

pause -1
