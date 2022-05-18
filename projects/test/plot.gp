set datafile separator '\t'
set term qt 0

# set xlabel "E^*"
# set ylabel "{/Symbol D}" offset 1,0 rotate by 0
# set title 'Cosh potential'

plot 'data.csv' using 1:2 w lines notitle
pause -1