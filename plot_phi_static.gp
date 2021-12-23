set xtics -0.2,0.2
set ytics -0.2,0.1

set format x "%g"
set format y "%g"
set xlabel ""
set ylabel ""
set xrange [-0.5:0.5]
set yrange [-0.3:0.05]
set xzeroaxis
set yzeroaxis

plot "C:/Users/hyoshida/Desktop/floquetic/phi_211223_N2_2.dat" u 2:3 title "N=2" w l lc "black"
replot "C:/Users/hyoshida/Desktop/floquetic/phi_211223_N1_2.dat" u 2:3 title "N=1" w l lc "blue"
replot "C:/Users/hyoshida/Desktop/floquetic/simN1.dat" u 1:3 notitle lc "black" pt 4

set output "C:/Users/hyoshida/Desktop/floquetic/phi_static_compare.png"
# set term postscript eps enhanced
# set output "C:/Users/hyoshida/Desktop/floquetic/phiN11.eps"

set terminal png
set term png enhanced
replot
reset
