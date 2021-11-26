################# rho plot ####################
set xrange [-1.0e+3:-1.0e-2]
set yrange [-0.1:10]
set xzeroaxis
set yzeroaxis
# set logscale xy
# set format y "%1.1e"
set xlabel 'z'
set ylabel '{/Symbol-Oblique r}(z)'
plot 'C:/Users/hyoshida/Desktop/floquetic/rho_211126_2.dat' u ($1):($2) title "rho"
# replot 'C:/Users/hyoshida/Desktop/floquetic/rho_211126_2.dat' u (abs($1)):($3) title "rho_d"
# replot 'C:/Users/hyoshida/Desktop/floquetic/rho_211126_2.dat' u (abs($1)):($4) title "rho_g"
set output "C:/Users/hyoshida/Desktop/floquetic/rho_211126_2_1.png"
set terminal png
set term png enhanced
replot
unset logscale xy
unset format
