# J-phi(J) plot for simulation and zero calculation
set xrange [-0.3:0.3]
set yrange [-0.2:0.0001]
set xzeroaxis
set xlabel 'J'
set ylabel '{/Symbol-Oblique j}(J)'
set key right bottom
plot 'C:/Users/hyoshida/Desktop/floquetic/sim_211210_N2_1.dat' title "sim rep = 2*10^3,iter=10^6"
replot 'C:/Users/hyoshida/Desktop/floquetic/phi_211210_N2_1.dat' using 2:3 title "{/Symbol-Oblique j}_U" with line
replot 'C:/Users/hyoshida/Desktop/floquetic/phi_211210_N2_1.dat' using 4:5 title "{/Symbol-Oblique j}_{{ad}}" with line

set output "C:/Users/hyoshida/Desktop/floquetic/phi_211210_N2_1.png"
set terminal png
set term png enhanced
replot
