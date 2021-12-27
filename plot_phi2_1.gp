# J-phi(J) plot for simulation and zero calculation
set term png enhanced
set output "C:/Users/hyoshida/Desktop/floquetic/phi_211227_1.png"
# set term postscript eps enhanced
# set output "C:/Users/hyoshida/Desktop/floquetic/phi_211227_1.eps"
set multiplot

set xtics -0.1,0.1,0.1
set ytics -0.05,0.025

# d=0 with N=2
set lmargin screen 0.1
set rmargin screen 0.9
set bmargin screen 0.51
set tmargin screen 0.9
set format x ""
set format y "%g"
set xlabel ""
set ylabel ""
set arrow from 0,-0.05 to 0,0 nohead lw 1 lc rgb "gray"
plot [-0.2:0.2][-0.05:0.0]"C:/Users/hyoshida/Desktop/floquetic/phi_211227.dat" u ($1):($2) notitle w l lc "black"
replot "C:/Users/hyoshida/Desktop/floquetic/phi_211227.dat" u ($3):($4) notitle w l lc "blue" dt "-"
replot "C:/Users/hyoshida/Desktop/floquetic/sim_211227.dat" u 1:2 notitle lc "black" pt 2
# replot "C:/Users/hyoshida/Desktop/floquetic/sim_211227_N2_1.dat" u 1:2 every 1 notitle lc "black" pt 4

# d=0 with N=4
set lmargin screen 0.1
set rmargin screen 0.9
set bmargin screen 0.1
set tmargin screen 0.49
set format x "%g"
set format y "%g"
set xlabel ""
set ylabel ""
set arrow from 0,-0.05 to 0,0 nohead lw 1 lc rgb "gray"
plot [-0.2:0.2][-0.05:0.0]"C:/Users/hyoshida/Desktop/floquetic/phi_211227.dat" u ($5):($6) notitle w l lc "black"
replot "C:/Users/hyoshida/Desktop/floquetic/phi_211227.dat" u ($7):($8) notitle w l lc "blue" dt "-"
replot "C:/Users/hyoshida/Desktop/floquetic/sim_211227.dat" u 1:3 notitle lc "black" pt 2

unset multiplot

# set term png enhanced
replot
reset
