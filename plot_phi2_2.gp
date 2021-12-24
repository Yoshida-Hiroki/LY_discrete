# J-phi(J) plot for simulation and zero calculation
# set term png enhanced
# set output "C:/Users/hyoshida/Desktop/floquetic/phi_211224_2.png"
set term postscript eps enhanced
set output "C:/Users/hyoshida/Desktop/floquetic/phi_211224_2.eps"
set multiplot

set xtics -0.1,0.1,0.1
set ytics -0.05,0.025

# g=0 with N=2
set lmargin screen 0.1
set rmargin screen 0.9
set bmargin screen 0.51
set tmargin screen 0.9
set format x ""
set format y "%g"
set xlabel ""
set ylabel ""
set arrow from 0,-0.05 to 0,0 nohead lw 2 lc rgb "black"
plot [-0.15:0.15][-0.05:0.0]"C:/Users/hyoshida/Desktop/floquetic/phi_211224.dat" u ($9):($10) notitle w l lc "black"
replot "C:/Users/hyoshida/Desktop/floquetic/phi_211224.dat" u ($11):($12) notitle w l lc "black" dt "-"
replot "C:/Users/hyoshida/Desktop/floquetic/sim_211224.dat" u 1:4 notitle lc "black" pt 4

# g=0 with N=4
set lmargin screen 0.1
set rmargin screen 0.9
set bmargin screen 0.1
set tmargin screen 0.49
set format x "%g"
set format y "%g"
set xlabel ""
set ylabel ""
set arrow from 0,-0.05 to 0,0 nohead lw 2 lc rgb "black"
plot [-0.15:0.15][-0.05:0.0]"C:/Users/hyoshida/Desktop/floquetic/phi_211224.dat" u ($13):($14) notitle w l lc "black"
replot "C:/Users/hyoshida/Desktop/floquetic/phi_211224.dat" u ($15):($16) notitle w l lc "black" dt "-"
replot "C:/Users/hyoshida/Desktop/floquetic/sim_211224.dat" u 1:5 notitle lc "black" pt 4

unset multiplot

# set term png enhanced
replot
reset
