set xtics -0.2,0.2
set ytics -0.2,0.1

set format x "%g"
set format y "%g"
set xlabel ""
set ylabel ""
set xrange [-0.2:0.2]
set yrange [-0.3:0.05]
set xzeroaxis
set yzeroaxis

plot "C:/Users/hyoshida/Desktop/floquetic/phiN1.dat" u 1:2 notitle w l lc "black"
replot "C:/Users/hyoshida/Desktop/floquetic/phiN1.dat" u 3:4 notitle w l lc "black"
replot "C:/Users/hyoshida/Desktop/floquetic/simN1.dat" u 1:2 notitle lc "black" pt 4
replot "C:/Users/hyoshida/Desktop/floquetic/simN1.dat" u 1:3 notitle lc "black" pt 4

set output "C:/Users/hyoshida/Desktop/floquetic/phiN1.png"
# set term postscript eps enhanced
# set output "C:/Users/hyoshida/Desktop/floquetic/phiN11.eps"

set terminal png
set term png enhanced
replot
reset
