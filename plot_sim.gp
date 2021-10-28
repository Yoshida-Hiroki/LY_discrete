# J-phi(J) plot for simulation and zero calculation
set xrange [-0.25:0]
set yrange [-0.05:0.001]
set xzeroaxis
set xlabel 'J'
set ylabel '{/Symbol-Oblique j}(J)'
plot 'C:/Users/hyoshida/Desktop/timedep/sim211028_1_N_1.dat' title "sim N=100,M=1,iter=10^5"
replot 'C:/Users/hyoshida/Desktop/timedep/sim211028_1_N_2.dat' title "sim N=100,M=10,iter=10^5"
replot 'C:/Users/hyoshida/Desktop/timedep/sim211028_1_N_3.dat' title "sim N=100,M=50,iter=10^5"
replot 'C:/Users/hyoshida/Desktop/timedep/sim211028_1_N_4.dat' title "sim N=100,M=100,iter=10^5"
set output "C:/Users/hyoshida/Desktop/timedep/sim_211028_1_N.png"
set terminal png
set term png enhanced
replot
