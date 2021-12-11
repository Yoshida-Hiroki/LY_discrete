################# rho plot ####################
set output "C:/Users/hyoshida/Desktop/floquetic/rho_211211_N2_2.png"
set terminal png

set xzeroaxis
set yzeroaxis
set logscale x
# set format y "%1.1e"
set xrange [] reverse
set xtics ("-10^{%L}" 100,"-10^{%L}" 10,"-10^{%L}" 1,"-10^{%L}" 0.1,"-10^{%L}" 0.01)
set xlabel 'z'
set ylabel '{/Symbol-Oblique r}(z)'
plot [:100][0:10]'C:/Users/hyoshida/Desktop/floquetic/rho_211211_N2_2.dat' u (abs($1)):($2) title "rho"

reset
