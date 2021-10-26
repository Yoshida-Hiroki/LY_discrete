#plot [][-0.05:]"C:/Users/hyoshida/Desktop/timedep/time_density_211024.dat" using 1:($2>0 ? $2: 1/0) with line
#replot "C:/Users/hyoshida/Desktop/timedep/time_density_211024.dat" using 1:($3>0 ? $3: 1/0) with line

#plot [-0.2:0.2][-0.3:0.05]'C:/Users/hyoshida/Desktop/timedep/time_simulation_211024.dat'

plot [-0.45:0][-0.3:0.05]'C:/Users/hyoshida/Desktop/timedep/time_simulation_211026.dat'
replot 'C:/Users/hyoshida/Desktop/timedep/time_density_211026.dat' using 2:3 with points lc 2
