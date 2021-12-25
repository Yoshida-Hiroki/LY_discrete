################# rho plot ####################
set output "C:/Users/hyoshida/Desktop/floquetic/rho_211224_1.png"
set terminal png
# set output "C:/Users/hyoshida/Desktop/floquetic/rho_211219_1.eps"
# set term postscript eps enhanced
set multiplot

set lmargin screen 0.1
set rmargin screen 0.9
set bmargin screen 0.51
set tmargin screen 0.9

set xzeroaxis
set yzeroaxis

# Jg = 0
set arrow from 4.35428941e+01,0.0 to 4.35428941e+01,10 nohead lw 1 lc rgb "gray"
set arrow from 8.79927236,0.0 to 8.79927236,10 nohead lw 1 lc rgb "gray"
set arrow from 5.01142277e-01,0.0 to 5.01142277e-01,10 nohead lw 1 lc rgb "gray"
set arrow from 5.20804715e-03,0.0 to 5.20804715e-03,10 nohead lw 1 lc rgb "gray"

set logscale x
# set format y "%1.1e"
set xtics ("" 1000,"" 100,"" 10,"" 1,"" 0.1,"" 0.01,"" 1.0e-03,"" 1.0e-04)
set xlabel ''
set ylabel '{/Symbol-Oblique r}(z)'
plot [3.0e+03:5.0e-04][:10]'C:/Users/hyoshida/Desktop/floquetic/rho_211224_N2_1.dat' u (abs($1)):($2) notitle w l lc "black" lw 8
reset

set lmargin screen 0.1
set rmargin screen 0.9
set bmargin screen 0.1
set tmargin screen 0.49

set xzeroaxis
set yzeroaxis

# Jg = 0
set arrow from 5.55797520e+02,0.0 to 5.55797520e+02,40 nohead lw 1 lc rgb "gray"
set arrow from 2.86302669e+02,0.0 to 2.86302669e+02,40 nohead lw 1 lc rgb "gray"
set arrow from 9.10222864,0.0 to 9.10222864,40 nohead lw 1 lc rgb "gray"
set arrow from 6.50095807,0.0 to 6.50095807,40 nohead lw 1 lc rgb "gray"
set arrow from 5.82454763e-01,0.0 to 5.82454763e-01,40 nohead lw 1 lc rgb "gray"
set arrow from 4.28196908e-01,0.0 to 4.28196908e-01,40 nohead lw 1 lc rgb "gray"
set arrow from 4.04276486e-02,0.0 to 4.04276486e-02,40 nohead lw 1 lc rgb "gray"
set arrow from 1.06374282e-02,0.0 to 1.06374282e-02,40 nohead lw 1 lc rgb "gray"

set logscale x
# set format y "%1.1e"
# set xrange [] reverse
set xtics ("-10^{%L}" 1000,"-10^{%L}" 100,"-10^{%L}" 10,"-1" 1,"-10^{%L}" 0.1,"-10^{%L}" 0.01,"-10^{%L}" 1.0e-03)
set xlabel 'z'
set ylabel '{/Symbol-Oblique r}(z)'
plot [3.0e+03:5.0e-04][:40]'C:/Users/hyoshida/Desktop/floquetic/rho_211224_N4_1.dat' u (abs($1)):($2) notitle w l lc "black" lw 8

reset

unset multiplot
