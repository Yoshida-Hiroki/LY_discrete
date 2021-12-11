################# rho plot ####################
set output "C:/Users/hyoshida/Desktop/floquetic/rho_211211_N2_1.png"
set terminal png

set xzeroaxis
set yzeroaxis
# set logscale xy
# set format y "%1.1e"
set xlabel 'z'
set ylabel '{/Symbol-Oblique r}(z)'
plot [-15:0][0:10]'C:/Users/hyoshida/Desktop/floquetic/rho_211211_N2_1.dat' u ($1):($2) title "rho"
