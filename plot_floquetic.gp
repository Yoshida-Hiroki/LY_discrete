# J-phi(J) plot for simulation and zero calculation
set xrange [-0.4:0.4]
set yrange [-0.3:0.0001]
set xzeroaxis
set xlabel 'J'
set ylabel '{/Symbol-Oblique j}(J)'
set key right bottom
plot 'C:/Users/hyoshida/Desktop/floquetic/sim_211203_N4_2.dat' title "sim rep = 2*10^3,iter=10^4"
replot 'C:/Users/hyoshida/Desktop/floquetic/phi_211203_N4_2.dat' using 2:3 title "{/Symbol-Oblique j}_U" with line
replot 'C:/Users/hyoshida/Desktop/floquetic/phi_211203_N4_2.dat' using 4:5 title "{/Symbol-Oblique j}_{{ad}}" with line

set output "C:/Users/hyoshida/Desktop/floquetic/phi_211203_N4_2.png"
set terminal png
set term png enhanced
replot
