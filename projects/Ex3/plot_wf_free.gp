set datafile separator '\t'

set term qt 0


#set title 'Eigenfunctions 3d'
plot 'wf1_free.csv' using 1:2 w lines t '0s'  , 'wf2_free.csv' using 1:2 w lines t '0p'
pause -1