set xlabel 'time'
set ylabel 'vacf'
set title 'Velocity autocorrelation function <v(t)xv(0)>'
set nokey
set grid
set xrange[0:*]
plot 'vacf.dat' us 1:2 notitle with li lw 3
pause -1 "Hit return to continue"
set title 'Integral <v(t)xv(0)>'
plot 'vacf.dat' us 1:3 notitle with li lw 3
pause -1 "Hit return to continue"
