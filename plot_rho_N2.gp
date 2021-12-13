################# rho plot ####################
set output "C:/Users/hyoshida/Desktop/floquetic/rho_211213_N2_2.png"
set terminal png

set xzeroaxis
set yzeroaxis

# Jd = 0
# set arrow from 13.36210652,0.0 to 13.36210652,40 nohead lw 2 lc rgb "gray"
# set arrow from 3,0.0 to 3,40 nohead lw 2 lc rgb "gray"
# set arrow from 0.33333,0.0 to 0.33333,40 nohead lw 2 lc rgb "gray"
# set arrow from 0.0748385,0.0 to 0.0748385,40 nohead lw 2 lc rgb "gray"

# Jg = 0
set arrow from 4.35428941e+01,0.0 to 4.35428941e+01,40 nohead lw 2 lc rgb "gray"
set arrow from 8.79927236,0.0 to 8.79927236,40 nohead lw 2 lc rgb "gray"
set arrow from 5.01142277e-01,0.0 to 5.01142277e-01,40 nohead lw 2 lc rgb "gray"
set arrow from 5.20804715e-03,0.0 to 5.20804715e-03,40 nohead lw 2 lc rgb "gray"

set logscale x
# set format y "%1.1e"
set xtics ("-10^{%L}" 1000,"-10^{%L}" 100,"-10^{%L}" 10,"-1" 1,"-10^{%L}" 0.1,"-10^{%L}" 0.01,"-10^{%L}" 1.0e-03,"-10^{%L}" 1.0e-04)
set xlabel 'z'
set ylabel '{/Symbol-Oblique r}(z)'
plot [1.0e+02:1.0e-03][:40]'C:/Users/hyoshida/Desktop/floquetic/rho_211213_N2_2.dat' u (abs($1)):($2) notitle w l lc "black" lw 5

reset
