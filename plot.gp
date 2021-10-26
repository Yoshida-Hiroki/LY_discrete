#plot [][-0.05:]"C:/Users/hyoshida/Desktop/timedep/time_density_211024.dat" using 1:($2>0 ? $2: 1/0) with line
#replot "C:/Users/hyoshida/Desktop/timedep/time_density_211024.dat" using 1:($3>0 ? $3: 1/0) with line

#plot [-0.2:0.2][-0.3:0.05]'C:/Users/hyoshida/Desktop/timedep/time_simulation_211024.dat'

set xrange [-0.25:-0.1]
set yrange [-0.05:0.01]
set xzeroaxis
set xlabel 'J'
set ylabel '{/Symbol-Oblique j}(J)'
plot 'C:/Users/hyoshida/Desktop/timedep/time_simulation_211026_10.dat' title "sim N=100,M=1,iter=10^5"
replot 'C:/Users/hyoshida/Desktop/timedep/time_simulation_211026_11.dat' title "sim N=100,M=10,iter=10^5"
replot 'C:/Users/hyoshida/Desktop/timedep/time_simulation_211026_12.dat' title "sim N=100,M=50,iter=10^5"
replot 'C:/Users/hyoshida/Desktop/timedep/time_simulation_211026_13.dat' title "sim N=100,M=100,iter=10^5"
#replot 'C:/Users/hyoshida/Desktop/timedep/time_density_211026_14.dat' using 2:3 title "-100" with line
#replot 'C:/Users/hyoshida/Desktop/timedep/time_density_211026_15.dat' using 2:3 title "-1000" with line
#replot 'C:/Users/hyoshida/Desktop/timedep/time_density_211026_16.dat' using 2:3 title "-10000" with line
set output "C:/Users/hyoshida/Desktop/timedep/211026_5.png"
set terminal png
set term png enhanced
replot

#plot 'C:/Users/hyoshida/Desktop/timedep/time_density_211026_5.dat' using 1:2 with points lc 2
