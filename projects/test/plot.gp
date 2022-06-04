set datafile separator '\t'
set term qt 0

plot 'data.csv' using 1:2 w lines title "0s" ,\
    'data.csv' using 1:3 w lines title "1s" ,\
    'data.csv' using 1:4 w lines title "2s" ,\
    'data.csv' using 1:5 w lines title "3s"
pause -1