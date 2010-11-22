set xlabel "Pattern number"
set ylabel "Mean intensity deviation"
unset key
plot "g-iteration-1.dat" using 1:2 w impulses lw 1 lc -1
