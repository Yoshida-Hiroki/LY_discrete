# # J-phi(J) plot for simulation and zero calculation
# set xrange [-0.25:-0.1]
# set yrange [-0.05:0.01]
# set xzeroaxis
# set xlabel 'J'
# set ylabel '{/Symbol-Oblique j}(J)'
# plot 'C:/Users/hyoshida/Desktop/timedep/sim_211027_10.dat' title "sim N=100,M=1,iter=10^5"
# replot 'C:/Users/hyoshida/Desktop/timedep/sim_211027_11.dat' title "sim N=100,M=10,iter=10^5"
# replot 'C:/Users/hyoshida/Desktop/timedep/sim_211027_12.dat' title "sim N=100,M=50,iter=10^5"
# replot 'C:/Users/hyoshida/Desktop/timedep/sim_211027_13.dat' title "sim N=100,M=100,iter=10^5" */
# #replot 'C:/Users/hyoshida/Desktop/timedep/phi_211027_14.dat' using 2:3 title "-100" with line
# #replot 'C:/Users/hyoshida/Desktop/timedep/phi_211027_15.dat' using 2:3 title "-1000" with line
# #replot 'C:/Users/hyoshida/Desktop/timedep/phi_211027_16.dat' using 2:3 title "-10000" with line
# set output "C:/Users/hyoshida/Desktop/timedep/phi_211027_5.png"
# set terminal png
# set term png enhanced
# replot
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
stats 'C:/Users/hyoshida/Desktop/timedep/f_211027_1.dat' u 3 every ::0::0 nooutput
set yrange [0:*]
x1 = STATS_mean
stats 'C:/Users/hyoshida/Desktop/timedep/f_211027_1.dat' using 4 every ::0::0 nooutput
set yrange [0:*]
x2 = STATS_mean
set label 1 at graph 0.8,0.1 sprintf("x1= %f",x1)
set label 2 at graph 0.8,0.05 sprintf("x2= %f",x2)

set xrange [-30:0]
set yrange [0:0.15]
set xzeroaxis
set yzeroaxis
set xlabel 'z'
set ylabel '<(a+b)/2*sqrt((z-z1)(z-z2)/(z(1-z1)(1-z2)))'
plot 'C:/Users/hyoshida/Desktop/timedep/f_211027_1.dat' using 1:2 title ""
set output "C:/Users/hyoshida/Desktop/timedep/f_211027_1.png"
set terminal png
set term png enhanced
replot
