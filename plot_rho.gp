################# rho plot ####################
set xrange [0.001:100]
set yrange [:10000]
set xzeroaxis
set yzeroaxis
set xlabel '-z'
set ylabel '{/Symbol-Oblique r}(z)'
plot 'C:/Users/hyoshida/Desktop/floquetic/rho_211121_3.dat' using (abs($1)):($2)
set logscale y 10
set logscale x 10
set output "C:/Users/hyoshida/Desktop/floquetic/rho_211121_3_3.png"
set terminal png
set term png enhanced
replot

unset logscale xy
