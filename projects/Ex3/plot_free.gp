set datafile separator '\t'
set term qt 0

plot "data.csv" using 1:2 w points title "final value"
pause -1

# xrange [:,5]
plot "wf.csv" using 1:2 w lines title "wavefunction"
pause -1