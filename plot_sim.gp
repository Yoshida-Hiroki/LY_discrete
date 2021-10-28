# J-phi(J) plot for simulation and zero calculation
set xrange [-0.25:0]
set yrange [-0.05:0.001]
set xzeroaxis
set xlabel 'J'
set ylabel '{/Symbol-Oblique j}(J)'
plot 'C:/Users/hyoshida/Desktop/timedep/sim211028_1_M_1.dat' title "sim N=5,M=100,iter=10^5"
replot 'C:/Users/hyoshida/Desktop/timedep/sim211028_1_M_2.dat' title "sim N=10,M=100,iter=10^5"
replot 'C:/Users/hyoshida/Desktop/timedep/sim211028_1_M_3.dat' title "sim N=50,M=100,iter=10^5"
replot 'C:/Users/hyoshida/Desktop/timedep/sim211028_1_M_4.dat' title "sim N=100,M=100,iter=10^5"
set output "C:/Users/hyoshida/Desktop/timedep/sim_211028_1_M.png"
set terminal png
set term png enhanced
replot
