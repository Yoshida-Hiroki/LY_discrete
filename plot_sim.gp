# J-phi(J) plot for simulation and zero calculation
set xrange [-0.21:-0.13]
set yrange [-0.01:0.001]
set xzeroaxis
set xlabel 'J'
set ylabel '{/Symbol-Oblique j}(J)'
plot 'C:/Users/hyoshida/Desktop/timedep/sim211028_1_MN_1.dat' with points ps 0.5 title "sim N=5,M=200,iter=10^5"
replot 'C:/Users/hyoshida/Desktop/timedep/sim211028_1_MN_2.dat' with points ps 0.5 title "sim N=10,M=100,iter=10^5"
replot 'C:/Users/hyoshida/Desktop/timedep/sim211028_1_MN_3.dat' with points ps 0.5 title "sim N= 50,M=20,iter=10^5"
replot 'C:/Users/hyoshida/Desktop/timedep/sim211028_1_MN_4.dat' with points ps 0.5 title "sim N= 100,M=10,iter=10^5"
set output "C:/Users/hyoshida/Desktop/timedep/sim_211028_1_MN.png"
set terminal png
set term png enhanced
replot
