set xlabel 'time'
set ylabel 'msd <r^2>'
set title 'Mean square displacement'
set nokey
set grid
set xrange[0:*]
plot 'msd.dat' us 1:2 notitle with li lw 3
pause -1 "Hit return to continue"
set ylabel 'D'
set title 'Diffusion coefficient'
plot 'msd.dat' us 1:3 notitle with li lw 3
pause -1 "Hit return to continue"
