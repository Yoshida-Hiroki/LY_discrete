# J-phi(J) plot for simulation and zero calculation
# set term png enhanced
# set output "C:/Users/hyoshida/Desktop/floquetic/phi_211227_1.png"
set term postscript eps enhanced
set output "C:/Users/hyoshida/Desktop/floquetic/phi_211227_1.eps"
set multiplot

set xtics -0.1,0.1,0.1
set ytics -0.05,0.025

# d=0 with N=2
set lmargin screen 0.1
set rmargin screen 0.9
set bmargin screen 0.51
set tmargin screen 0.9
set xtics ("" -0.15 1,"" -0.1 0, "" -0.05 1,"" 0 0, "" 0.05 1, "" 0.1 0,"" 0.15 1)
set ytics ("-0.1" -0.1 0, "" -0.075 1, "-0.05" -0.05 0, "" -0.025 1, "0" 0 0)
set xlabel ""
set ylabel ""
set arrow from 0,-0.1 to 0,0.01 nohead lw 1
set arrow from -0.15,0 to 0.15,0 nohead lw 1
plot [-0.15:0.15][-0.1:0.01]"C:/Users/hyoshida/Desktop/floquetic/phi_211227.dat" u ($1):($2) notitle w l lc "black" lw 3
replot "C:/Users/hyoshida/Desktop/floquetic/phi_211227.dat" u ($3):($4) notitle w l lc "black" dt "-" lw 3
replot "C:/Users/hyoshida/Desktop/floquetic/sim_211225.dat" u 1:2 every 20 notitle lc "black" pt 2 ps 2

# d=0 with N=4
set lmargin screen 0.1
set rmargin screen 0.9
set bmargin screen 0.1
set tmargin screen 0.49
set xtics ("" -0.15 1,"-0.1" -0.1 0, "" -0.05 1,"0" 0 0, "" 0.05 1, "0.1" 0.1 0,"" 0.15 1)
set ytics ("-0.1" -0.1 0, "" -0.075 1, "-0.05" -0.05 0, "" -0.025 1, "0" 0 0)
set xlabel ""
set ylabel ""
set arrow from 0,-0.1 to 0,0.01 nohead lw 1
set arrow from -0.15,0 to 0.15,0 nohead lw 1
plot [-0.15:0.15][-0.1:0.01]"C:/Users/hyoshida/Desktop/floquetic/phi_211227.dat" u ($5):($6) notitle w l lc "black" lw 3
replot "C:/Users/hyoshida/Desktop/floquetic/phi_211227.dat" u ($7):($8) notitle w l lc "black" dt "-" lw 3
replot "C:/Users/hyoshida/Desktop/floquetic/sim_211225.dat" u 1:3 every 20 notitle lc "black" pt 2 ps 2

unset multiplot

# set term png enhanced
replot
reset
