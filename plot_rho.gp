################# rho plot ####################
set xrange [-6:0]
set yrange [-0.01:2]
set xzeroaxis
set yzeroaxis
set xlabel 'z'
set ylabel '{/Symbol-Oblique r}(z)'
plot 'C:/Users/hyoshida/Desktop/floquetic/rho_211118_1.dat'
set output "C:/Users/hyoshida/Desktop/floquetic/rho_211118_1.png"
set terminal png
set term png enhanced
replot
