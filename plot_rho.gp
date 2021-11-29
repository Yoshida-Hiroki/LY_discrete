################# rho plot ####################
set output "C:/Users/hyoshida/Desktop/floquetic/rho_211129_N2_3.png"
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
plot [-10:0][-0.1:3]'C:/Users/hyoshida/Desktop/floquetic/rho_211129_N2_3.dat' u ($1):($2) title "rho"
replot 'C:/Users/hyoshida/Desktop/floquetic/rho_211129_N2_3.dat' u ($1):($3) title "rho_d"
replot 'C:/Users/hyoshida/Desktop/floquetic/rho_211129_N2_3.dat' u ($1):($4) title "rho_g"

set size 0.5,0.5
set origin 0.22,0.5
plot [-0.34:0][-0.1:11]'C:/Users/hyoshida/Desktop/floquetic/rho_211129_N2_3.dat' u ($1):($2) notitle
replot 'C:/Users/hyoshida/Desktop/floquetic/rho_211129_N2_3.dat' u ($1):($3) notitle
replot 'C:/Users/hyoshida/Desktop/floquetic/rho_211129_N2_3.dat' u ($1):($4) notitle

unset multiplot
reset

set term png enhanced
replot
# unset logscale xy
# unset format
