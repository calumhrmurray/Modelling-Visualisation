set xrange [0:1000]
set yrange [0:1]
plot for [col=2:4] 'plot.dat' using 1:col with lines
pause 0.033
reread
