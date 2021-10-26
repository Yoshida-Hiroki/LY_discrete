# J-phi(J) plot for simulation and zero calculation
set xrange [-0.25:-0.1]
set yrange [-0.05:0.01]
set xzeroaxis
set xlabel 'J'
set ylabel '{/Symbol-Oblique j}(J)'
plot 'C:/Users/hyoshida/Desktop/timedep/sim_211026_10.dat' title "sim N=100,M=1,iter=10^5"
replot 'C:/Users/hyoshida/Desktop/timedep/sim_211026_11.dat' title "sim N=100,M=10,iter=10^5"
replot 'C:/Users/hyoshida/Desktop/timedep/sim_211026_12.dat' title "sim N=100,M=50,iter=10^5"
replot 'C:/Users/hyoshida/Desktop/timedep/sim_211026_13.dat' title "sim N=100,M=100,iter=10^5"
#replot 'C:/Users/hyoshida/Desktop/timedep/phi_211026_14.dat' using 2:3 title "-100" with line
#replot 'C:/Users/hyoshida/Desktop/timedep/phi_211026_15.dat' using 2:3 title "-1000" with line
#replot 'C:/Users/hyoshida/Desktop/timedep/phi_211026_16.dat' using 2:3 title "-10000" with line
set output "C:/Users/hyoshida/Desktop/timedep/211026_5.png"
set terminal png
set term png enhanced
replot

################# rho plot ####################
set xrange [-30:0]
set yrange [-0.1:2]
set xzeroaxis
set yzeroaxis
set xlabel 'z'
set ylabel '{/Symbol-Oblique r}(z)'
plot 'C:/Users/hyoshida/Desktop/timedep/rho_211026_10.dat'
