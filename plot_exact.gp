# J-phi(J) plot for simulation and zero calculation
set xrange [-0.05:0.05]
set yrange [-0.005:0.0001]
set xzeroaxis
set xlabel 'J'
set ylabel '{/Symbol-Oblique j}(J)'

plot 'C:/Users/hyoshida/Desktop/timedep/sim_211105_1.dat' title "sim N=100,M=100,iter=10^5"

replot for [i=1:10] 'C:/Users/hyoshida/Desktop/timedep/exact_zeros_211105_1_current.dat' u 2*i-1:2*i title sprintf("t=%d",i) w line

# replot 'C:/Users/hyoshida/Desktop/timedep/phi_211028_1_N5.dat' using 2:3 title "N=5" with line
# replot 'C:/Users/hyoshida/Desktop/timedep/phi_211028_1_N10.dat' using 2:3 title "N=10" with line
# replot 'C:/Users/hyoshida/Desktop/timedep/phi_211028_1_N50.dat' using 2:3 title "N=50" with line
# replot 'C:/Users/hyoshida/Desktop/timedep/phi_211028_1_N100.dat' using 2:3 title "N=100" with line

# replot 'C:/Users/hyoshida/Desktop/timedep/phi_211028_1_2points_2_N100000.dat' using 2:3 title "2 points" with line

set output "C:/Users/hyoshida/Desktop/timedep/phi_211105_1.png"
set terminal png
set term png enhanced
replot
