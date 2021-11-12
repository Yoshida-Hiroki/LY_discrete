# J-phi(J) plot for simulation and zero calculation
set xrange [-0.6:-0.3]
set yrange [-0.005:0.0001]
set xzeroaxis
set xlabel 'J'
set ylabel '{/Symbol-Oblique j}(J)'
# plot 'C:/Users/hyoshida/Desktop/floquetic/sim_211109_2.dat' title "sim N=100,M=100,iter=10^5"
plot 'C:/Users/hyoshida/Desktop/floquetic/phi_211111_1.dat' using 1:2 title "" with line

set output "C:/Users/hyoshida/Desktop/floquetic/phi_211111_1.png"
set terminal png
set term png enhanced
replot
