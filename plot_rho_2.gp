################# rho plot ####################
# set output "C:/Users/hyoshida/Desktop/floquetic/rho_211224_2.png"
# set terminal png
set output "C:/Users/hyoshida/Desktop/floquetic/rho_211224_2.eps"
set term postscript eps enhanced
set multiplot

set lmargin screen 0.1
set rmargin screen 0.9
set bmargin screen 0.51
set tmargin screen 0.9

set xzeroaxis
set yzeroaxis

# Jd = 0
set arrow from 13.36210652,0.0 to 13.36210652,10 nohead lw 1
set arrow from 3,0.0 to 3,10 nohead lw 1
set arrow from 0.33333,0.0 to 0.33333,10 nohead lw 1
set arrow from 0.0748385,0.0 to 0.0748385,10 nohead lw 1

set logscale x
# set format y "%1.1e"
set xtics ("" 1000,"" 100 1,"" 10 1,"" 1,"" 0.1 1,"" 0.01 1,"" 1.0e-03,"" 1.0e-04)
set xlabel ''
set ylabel ''
set ytics ("0" 0 0, "" 2 1,"4" 4 0, "" 6 1, "8" 8 0)
plot [3.0e+03:5.0e-04][0:10]'C:/Users/hyoshida/Desktop/floquetic/rho_211224_N2_2.dat' u (abs($1)):($2) notitle w l lc "black" lw 5
reset

set lmargin screen 0.1
set rmargin screen 0.9
set bmargin screen 0.1
set tmargin screen 0.49

set xzeroaxis
set yzeroaxis

# Jd = 0
set arrow from 6.85734322e+02,0.0 to 6.85734322e+02,50 nohead lw 1
set arrow from 2.692611e+01,0.0 to 2.692611e+01,50 nohead lw 1
set arrow from 6.26131355,0.0 to 6.26131355,50 nohead lw 1
set arrow from 4.29807545,0.0 to 4.29807545,50 nohead lw 1
set arrow from 4.42992457e-01,0.0 to 4.42992457e-01,50 nohead lw 1
set arrow from 3.80722489e-01,0.0 to 3.80722489e-01,50 nohead lw 1
set arrow from 1.55794600e-02,0.0 to 1.55794600e-02,50 nohead lw 1
set arrow from 7.65902956e-04,0.0 to 7.65902956e-04,50 nohead lw 1

set logscale x
set xtics ("-10^{%L}" 1000,"" 100 1,"" 10 1,"-1" 1,"" 0.1 1,"" 0.01 1,"-10^{%L}" 1.0e-03)
set xlabel ''
set ylabel ''
set ytics ("0" 0 0, "" 10 1,"20" 20 0, "" 30 1, "40" 40 0)
plot [3.0e+03:5.0e-04][0:50]'C:/Users/hyoshida/Desktop/floquetic/rho_211224_N4_2.dat' u (abs($1)):($2) notitle w l lc "black" lw 5

reset

unset multiplot
