set xlabel 'distance r'
set ylabel 'g(r)'
set title 'Radial distribution function'
set nokey
set grid
set xrange[0:*]
plot 'rdf.dat' us 1:2 notitle with li lw 3
pause -1 "Hit return to continue"
