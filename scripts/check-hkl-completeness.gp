set ylabel "Completeness [blue]"
set y2label "Reflections measured [red] Possible [pink]"
set ytics nomirror
set yrange [0:100]
#set y2range [0:10]
#set xrange [0.1:1.15]
unset key
set y2tics
set xlabel "Resolution (one over d) (1/nm)"
plot "shells.dat" using 1:2 w l lw 3 lc 1 axis x1y2, \
     "shells.dat" using 1:3 w l lw 3 lc 4 axis x1y2, \
     "shells.dat" using 1:4 w l lw 3 lc 3 axis x1y1
