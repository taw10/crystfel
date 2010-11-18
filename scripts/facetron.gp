set xlabel "Calculated partiality"
set ylabel "Observed Partiality"
set xrange [0:1]
set yrange [0:1]
unset key
set size square
plot "iteration-1.dat" using 4:5 w p ps 1 pt 7 lc -1
