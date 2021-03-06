# J-phi(J) plot for simulation and zero calculation
set term png enhanced
set output "C:/Users/hyoshida/Desktop/floquetic/phi_211219_1.png"
# set term postscript eps enhanced
# set output "C:/Users/hyoshida/Desktop/floquetic/phi_211219_1.eps"
set multiplot

set xtics -0.2,0.2
set ytics -0.2,0.1

# g=0 with N=2 (left top)
set lmargin screen 0.1
set rmargin screen 0.49
set tmargin screen 0.9
set bmargin screen 0.52
set format x ""
set format y "%g"
set xlabel ""
set ylabel ""
set arrow from 0,-0.2 to 0,0 nohead lw 2 lc rgb "black"
plot [-0.3:0.3][-0.2:0.001]"C:/Users/hyoshida/Desktop/floquetic/phi_211210.dat" u ($5):($6) notitle w l lc "black"
replot "C:/Users/hyoshida/Desktop/floquetic/phi_211210.dat" u ($7):($8) notitle w l lc "black" dt "-"
replot "C:/Users/hyoshida/Desktop/floquetic/sim_211210.dat" u 1:3 notitle lc "black" pt 4

# d=0 with N=2 (left bottom)
set lmargin screen 0.1
set rmargin screen 0.49
set tmargin screen 0.48
set bmargin screen 0.1
set format x "%g"
set format y "%g"
set xlabel ""
set ylabel ""
set arrow from 0,-0.2 to 0,0 nohead lw 2 lc rgb "black"
plot [-0.3:0.3][-0.2:0.001]"C:/Users/hyoshida/Desktop/floquetic/phi_211210.dat" u ($1):($2) notitle w l lc "black"
replot "C:/Users/hyoshida/Desktop/floquetic/phi_211210.dat" u ($3):($4) notitle w l lc "black" dt "-"
replot "C:/Users/hyoshida/Desktop/floquetic/sim_211210.dat" u 1:2 notitle lc "black" pt 4

# g=0 with N=4 (right top)
set lmargin screen 0.51
set rmargin screen 0.9
set tmargin screen 0.9
set bmargin screen 0.52
set format x ""
set format y ""
set xlabel ""
set ylabel ""
set arrow from 0,-0.2 to 0,0 nohead lw 2 lc rgb "black"
plot [-0.3:0.3][-0.2:0.001]"C:/Users/hyoshida/Desktop/floquetic/phi_211210.dat" u ($13):($14) notitle w l lc "black"
replot "C:/Users/hyoshida/Desktop/floquetic/phi_211210.dat" u ($15):($16) notitle w l lc "black" dt (10,20)
replot "C:/Users/hyoshida/Desktop/floquetic/sim_211210.dat" u 1:5 notitle lc "black" pt 4

# d=0 with N=4 (right bottom)
set lmargin screen 0.51
set rmargin screen 0.9
set tmargin screen 0.48
set bmargin screen 0.1
set format x "%g"
set format y ""
set xlabel ""
set ylabel ""
set arrow from 0,-0.2 to 0,0 nohead lw 2 lc rgb "black"
plot [-0.3:0.3][-0.2:0.001]"C:/Users/hyoshida/Desktop/floquetic/phi_211210.dat" u ($9):($10) notitle w l lc "black"
replot "C:/Users/hyoshida/Desktop/floquetic/phi_211210.dat" u ($11):($12) notitle w l lc "black" dt "-"
replot "C:/Users/hyoshida/Desktop/floquetic/sim_211210.dat" u 1:4 notitle lc "black" pt 4

unset multiplot

# set term png enhanced
replot
