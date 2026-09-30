a = 1.0
gauss(x)=a*exp(-x*x/2.0)
set xlabel "Sample quantile"
set ylabel "Theoretical quantile"
set y2label "Density"
set ytics nomirror
set y2tics
set key top left
plot [-5:5] "normal.plot" using 3:5 w lp title "Theoretical quantile"
replot "normal.plot" using 3:4 w fsteps axis x1y2 title "Sample density"
replot [-2:2] x title "x"
fit [-0.2:0.2] gauss(x) 'normal.plot' using 2:4 via a
replot gauss(x) w l axis x1y2
