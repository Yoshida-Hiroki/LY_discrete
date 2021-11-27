################# rho plot ####################
set output "C:/Users/hyoshida/Desktop/floquetic/rho_211127_N2_2.png"
set terminal png
set multiplot

set size 1,1
set origin 0,0
set xzeroaxis
set yzeroaxis
# set logscale xy
# set format y "%1.1e"
set xlabel 'z'
set ylabel '{/Symbol-Oblique r}(z)'
plot [-55:0][-0.1:3]'C:/Users/hyoshida/Desktop/floquetic/rho_211127_N2_2.dat' u ($1):($2) title "rho"
replot 'C:/Users/hyoshida/Desktop/floquetic/rho_211127_N2_2.dat' u ($1):($3) title "rho_d"
replot 'C:/Users/hyoshida/Desktop/floquetic/rho_211127_N2_2.dat' u ($1):($4) title "rho_g"

set size 0.6,0.6
set origin 0.2,0.4
plot [-0.34:0][-0.1:11]'C:/Users/hyoshida/Desktop/floquetic/rho_211127_N2_2.dat' u ($1):($2) notitle
replot 'C:/Users/hyoshida/Desktop/floquetic/rho_211127_N2_2.dat' u ($1):($3) notitle
replot 'C:/Users/hyoshida/Desktop/floquetic/rho_211127_N2_2.dat' u ($1):($4) notitle

unset multiplot
reset

set term png enhanced
replot
# unset logscale xy
# unset format
