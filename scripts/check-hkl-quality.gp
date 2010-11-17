set ylabel "Redundancy [blue]"
set y2label "Post-merging SNR [red]"
set ytics nomirror
#set yrange [0:100]
#set y2range [0:10]
#set xrange [0.1:1.15]
unset key
set y2tics
set xlabel "Resolution (one over d) (1/nm)"
plot "shells.dat" using 1:6 w l lw 3 lc 3 axis x1y1, \
     "shells.dat" using 1:7 w l lw 3 lc 1 axis x1y2
