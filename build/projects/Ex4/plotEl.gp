set datafile separator '\t'
set term qt 0

set grid
plot "dataEnEl.csv" u 1:2
#plot "data_mEn.csv" u 1:3

pause -1
