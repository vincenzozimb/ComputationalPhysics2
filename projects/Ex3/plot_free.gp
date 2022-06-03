set datafile separator '\t'
set term qt 0

# plot "data.csv" using 1:2 w dots notitle
# pause -1

plot "density.csv" using 1:2 w lines title "electronic density for N = 8"
pause -1