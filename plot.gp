# J-phi(J) plot for simulation and zero calculation
set xrange [-0.25:-0.1]
set yrange [-0.01:0.01]
set xzeroaxis
set xlabel 'J'
set ylabel '{/Symbol-Oblique j}(J)'
plot 'C:/Users/hyoshida/Desktop/timedep/sim_211028_1.dat' title "sim N=100,M=100,iter=10^5"
# replot 'C:/Users/hyoshida/Desktop/timedep/sim_211028_11.dat' title "sim N=100,M=10,iter=10^5"
# replot 'C:/Users/hyoshida/Desktop/timedep/sim_211028_12.dat' title "sim N=100,M=50,iter=10^5"
# replot 'C:/Users/hyoshida/Desktop/timedep/sim_211028_13.dat' title "sim N=100,M=100,iter=10^5" */
replot 'C:/Users/hyoshida/Desktop/timedep/phi_211028_1_N5.dat' using 2:3 title "N=5" with line
replot 'C:/Users/hyoshida/Desktop/timedep/phi_211028_1_N10.dat' using 2:3 title "N=10" with line
replot 'C:/Users/hyoshida/Desktop/timedep/phi_211028_1_N50.dat' using 2:3 title "N=50" with line
replot 'C:/Users/hyoshida/Desktop/timedep/phi_211028_1_N100.dat' using 2:3 title "N=100" with line
set output "C:/Users/hyoshida/Desktop/timedep/phi_211028_1.png"
set terminal png
set term png enhanced
replot
#
# ################# rho plot ####################
# set xrange [-30:0]
# set yrange [-0.1:2]
# set xzeroaxis
# set yzeroaxis
# set xlabel 'z'
# set ylabel '{/Symbol-Oblique r}(z)'
# plot 'C:/Users/hyoshida/Desktop/timedep/rho_211027_10.dat'
# set output "C:/Users/hyoshida/Desktop/timedep/rho_211027_5.png"
# set terminal png
# set term png enhanced
# replot

################## f(z) plot ####################
# set label 1 at graph 0.8,0.1 sprintf("x1= %f",x1)
# set label 2 at graph 0.8,0.05 sprintf("x2= %f",x2)


set xzeroaxis
set yzeroaxis
set xlabel 'z'
set ylabel '<(a+b)/2*sqrt((z-z1)(z-z2)/(z(1-z1)(1-z2)))'

set yrange [0:*]
stats 'C:/Users/hyoshida/Desktop/timedep/f_211028_1_N5.dat' u 3 every ::0::0 nooutput
set yrange [0:*]
x1_N5 = STATS_mean
stats 'C:/Users/hyoshida/Desktop/timedep/f_211028_1_N5.dat' using 4 every ::0::0 nooutput
set yrange [0:*]
x2_N5 = STATS_mean
plot 'C:/Users/hyoshida/Desktop/timedep/f_211028_1_N5.dat' using 1:2 with points pointsize 2 title sprintf("N=5, x1= %f, x2 = %f", x1_N5,x2_N5)

stats 'C:/Users/hyoshida/Desktop/timedep/f_211028_1_N10.dat' u 3 every ::0::0 nooutput
set yrange [0:*]
x1_N10 = STATS_mean
stats 'C:/Users/hyoshida/Desktop/timedep/f_211028_1_N10.dat' using 4 every ::0::0 nooutput
set yrange [0:*]
x2_N10 = STATS_mean
replot 'C:/Users/hyoshida/Desktop/timedep/f_211028_1_N10.dat' using 1:2 with points pointsize 2 title sprintf("N=10, x1= %f, x2 = %f", x1_N10,x2_N10)

stats 'C:/Users/hyoshida/Desktop/timedep/f_211028_1_N50.dat' u 3 every ::0::0 nooutput
set yrange [0:*]
x1_N50 = STATS_mean
stats 'C:/Users/hyoshida/Desktop/timedep/f_211028_1_N50.dat' using 4 every ::0::0 nooutput
set yrange [0:*]
x2_N50 = STATS_mean
replot 'C:/Users/hyoshida/Desktop/timedep/f_211028_1_N50.dat' using 1:2 with points pointsize 2 title sprintf("N=50, x1= %f, x2 = %f", x1_N50,x2_N50)

stats 'C:/Users/hyoshida/Desktop/timedep/f_211028_1_N100.dat' u 3 every ::0::0 nooutput
set yrange [0:*]
x1_N100 = STATS_mean
stats 'C:/Users/hyoshida/Desktop/timedep/f_211028_1_N100.dat' using 4 every ::0::0 nooutput
set yrange [0:*]
x2_N100 = STATS_mean
replot 'C:/Users/hyoshida/Desktop/timedep/f_211028_1_N100.dat' using 1:2 with points pointsize 2 title sprintf("N=100, x1= %f, x2 = %f", x1_N100,x2_N100)

set xrange [-30:0]
set yrange [0:0.5]

set output "C:/Users/hyoshida/Desktop/timedep/f_211028_1.png"
set terminal png
set term png enhanced
replot
