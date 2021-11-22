################# rho plot ####################
set xrange [-5:0]
set yrange [-0.01:2000]
set xzeroaxis
set yzeroaxis
set xlabel 'z'
set ylabel '{/Symbol-Oblique r}(z)'
plot 'C:/Users/hyoshida/Desktop/floquetic/rho_211121_3_2.dat'
set output "C:/Users/hyoshida/Desktop/floquetic/rho_211121_3_2.png"
set terminal png
set term png enhanced
replot
