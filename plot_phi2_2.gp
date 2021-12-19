# J-phi(J) plot for simulation and zero calculation
# set term png enhanced
# set output "C:/Users/hyoshida/Desktop/floquetic/phi_211219_2.png"
set term postscript eps enhanced
set output "C:/Users/hyoshida/Desktop/floquetic/phi_211219_2.eps"
set multiplot

set xtics -0.2,0.2
set ytics -0.2,0.1

# g=0 with N=2
set lmargin screen 0.1
set rmargin screen 0.9
set bmargin screen 0.53
set tmargin screen 0.9
set format x ""
set format y "%g"
set xlabel ""
set ylabel ""
set arrow from 0,-0.2 to 0,0 nohead lw 2 lc rgb "black"
plot [-0.3:0.3][-0.2:0.001]"C:/Users/hyoshida/Desktop/floquetic/phi_211210.dat" u ($5):($6) notitle w l lc "black"
replot "C:/Users/hyoshida/Desktop/floquetic/phi_211210.dat" u ($7):($8) notitle w l lc "black" dt "-"
replot "C:/Users/hyoshida/Desktop/floquetic/sim_211210.dat" u 1:3 notitle lc "black" pt 4

# g=0 with N=4
set lmargin screen 0.1
set rmargin screen 0.9
set bmargin screen 0.1
set tmargin screen 0.47
set format x "%g"
set format y "%g"
set xlabel ""
set ylabel ""
set arrow from 0,-0.2 to 0,0 nohead lw 2 lc rgb "black"
plot [-0.3:0.3][-0.2:0.001]"C:/Users/hyoshida/Desktop/floquetic/phi_211210.dat" u ($13):($14) notitle w l lc "black"
replot "C:/Users/hyoshida/Desktop/floquetic/phi_211210.dat" u ($15):($16) notitle w l lc "black" dt (10,20)
replot "C:/Users/hyoshida/Desktop/floquetic/sim_211210.dat" u 1:5 notitle lc "black" pt 4

unset multiplot

# set term png enhanced
replot
