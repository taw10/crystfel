set ylabel "R ( sum(|I1-kI2|) / sum(I1) ) (%) [blue]"
set ytics nomirror
set yrange [0:120]
#set y2range [0:50000]
#set xrange [0.1:1.15]
unset key
#set y2tics
set xlabel "Resolution (one over d) (1/nm)"
plot "shells.dat" using 1:2 w l lw 3 lc 3 axis x1y1
