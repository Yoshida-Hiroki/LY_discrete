# J-phi(J) plot for simulation and zero calculation
set xrange [-2:0.2]
set yrange [-0.1:0.0001]
set xzeroaxis
set xlabel 'J'
set ylabel '{/Symbol-Oblique j}(J)'
set key right bottom
plot 'C:/Users/hyoshida/Desktop/floquetic/sim_211118_1.dat' title "sim rep = 2*10^3,iter=10^4"
replot 'C:/Users/hyoshida/Desktop/floquetic/phi_211118_1.dat' using 2:3 title "" with line

set output "C:/Users/hyoshida/Desktop/floquetic/phi_211118_1.png"
set terminal png
set term png enhanced
replot
