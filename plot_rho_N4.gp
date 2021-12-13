################# rho plot ####################
set output "C:/Users/hyoshida/Desktop/floquetic/rho_211213_N4_2.png"
set terminal png

set xzeroaxis
set yzeroaxis
# Jd = 0
# set arrow from 6.85734322e+02,0.0 to 6.85734322e+02,400 nohead lw 2 lc rgb "gray"
# set arrow from 2.692611e+01,0.0 to 2.692611e+01,400 nohead lw 2 lc rgb "gray"
# set arrow from 6.26131355,0.0 to 6.26131355,400 nohead lw 2 lc rgb "gray"
# set arrow from 4.29807545,0.0 to 4.29807545,400 nohead lw 2 lc rgb "gray"
# set arrow from 4.42992457e-01,0.0 to 4.42992457e-01,400 nohead lw 2 lc rgb "gray"
# set arrow from 3.80722489e-01,0.0 to 3.80722489e-01,400 nohead lw 2 lc rgb "gray"
# set arrow from 1.55794600e-02,0.0 to 1.55794600e-02,400 nohead lw 2 lc rgb "gray"
# set arrow from 7.65902956e-04,0.0 to 7.65902956e-04,400 nohead lw 2 lc rgb "gray"

# Jg = 0
set arrow from 5.55797520e+02,0.0 to 5.55797520e+02,200 nohead lw 2 lc rgb "gray"
set arrow from 2.86302669e+02,0.0 to 2.86302669e+02,200 nohead lw 2 lc rgb "gray"
set arrow from 9.10222864,0.0 to 9.10222864,200 nohead lw 2 lc rgb "gray"
set arrow from 6.50095807,0.0 to 6.50095807,200 nohead lw 2 lc rgb "gray"
set arrow from 5.82454763e-01,0.0 to 5.82454763e-01,200 nohead lw 2 lc rgb "gray"
set arrow from 4.28196908e-01,0.0 to 4.28196908e-01,200 nohead lw 2 lc rgb "gray"
set arrow from 4.04276486e-02,0.0 to 4.04276486e-02,200 nohead lw 2 lc rgb "gray"
set arrow from 1.06374282e-02,0.0 to 1.06374282e-02,200 nohead lw 2 lc rgb "gray"

set logscale x
# set format y "%1.1e"
# set xrange [] reverse
set xtics ("-10^{%L}" 1000,"-10^{%L}" 100,"-10^{%L}" 10,"-1" 1,"-10^{%L}" 0.1,"-10^{%L}" 0.01,"-10^{%L}" 1.0e-03,"-10^{%L}" 1.0e-04)
set xlabel 'z'
set ylabel '{/Symbol-Oblique r}(z)'
plot [1.0e+03:1.0e-02][:200]'C:/Users/hyoshida/Desktop/floquetic/rho_211213_N4_2.dat' u (abs($1)):($2) notitle w l lc "black" lw 5

reset
