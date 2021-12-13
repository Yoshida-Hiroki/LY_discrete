################# rho plot ####################
set output "C:/Users/hyoshida/Desktop/floquetic/rho_211213_N2_1.png"
set terminal png

set xzeroaxis
set yzeroaxis
set arrow from 13.36210652,0.0 to 13.36210652,40 nohead lw 2 lc rgb "gray"
set arrow from 3,0.0 to 3,40 nohead lw 2 lc rgb "gray"
set arrow from 0.33333,0.0 to 0.33333,40 nohead lw 2 lc rgb "gray"
set arrow from 0.0748385,0.0 to 0.0748385,40 nohead lw 2 lc rgb "gray"
set logscale x
# set format y "%1.1e"
set xtics ("-10^{%L}" 10,"-1" 1,"-10^{%L}" 0.1,"-10^{%L}" 0.01)
set xlabel 'z'
set ylabel '{/Symbol-Oblique r}(z)'
plot [3.0e+01:5.0e-02][:40]'C:/Users/hyoshida/Desktop/floquetic/rho_211213_N2_1.dat' u (abs($1)):($2) notitle w l lc "black" lw 5

reset
