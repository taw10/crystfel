set ylabel "Rsplit (%)"
set y2label "CChalf (fraction)"
set ytics nomirror
set yrange [0:90]
set y2range [0:1]
set key bottom right
set y2tics

set xtics ("100" 0.1000, "50" 0.2000, "10" 1.000, "8" 1.250, "6" 1.667, "5" 2.000, "4" 2.500, "3.5" 2.857, "3" 3.333, "2.5" 4.000, "2" 5.000, "1.8" 5.556, "1.6" 6.250, "1.4" 7.143, "1.2" 8.333, "1" 10.00, "0.9" 11.11, "0.8" 12.50)
set xlabel "Resolution d (= lambda/2sin(theta)) / Angstrom"

plot "rsplit.dat" using 1:2 w l lw 3 axis x1y1 title "Rsplit / %"
replot "cchalf.dat" using 1:2 w l lw 3 axis x1y2 title "CChalf"
