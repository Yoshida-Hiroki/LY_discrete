################# rho plot ####################
# set output "C:/Users/hyoshida/Desktop/figure/rho_1.png"
# set terminal png
set output "C:/Users/hyoshida/Desktop/figure/rho_N1_1.eps"
set term postscript eps enhanced
set multiplot

set lmargin screen 0.1
set rmargin screen 0.9
set bmargin screen 0.51
set tmargin screen 0.9

set xzeroaxis
set yzeroaxis

# Jg = 0
set arrow from -1.5,0.0 to -1.5,5 nohead lw 1
set arrow from -0.66666666666,0.0 to -0.66666666666,5 nohead lw 1

set xtics -2,0.2,0
set xlabel ''
set ylabel ''
set ytics 0,1,5
plot [-2:0][0:5]'C:/Users/hyoshida/Desktop/floquetic/rho_211225_N1_1.dat' u ($1):($2) notitle w l lc "black" lw 1
reset

set lmargin screen 0.1
set rmargin screen 0.9
set bmargin screen 0.1
set tmargin screen 0.49

set xzeroaxis
set yzeroaxis

# Jg = 0
set arrow from -1.5,0.0 to -1.5,5 nohead lw 1
set arrow from -0.25,0.0 to -0.25,5 nohead lw 1

set xtics -2,0.2,0
set xlabel ''
set ylabel ''
set ytics 0,1,5
plot [-2:0][0:5]'C:/Users/hyoshida/Desktop/floquetic/rho_211225_N1_2.dat' u ($1):($2) notitle w l lc "black" lw 1

reset

unset multiplot
