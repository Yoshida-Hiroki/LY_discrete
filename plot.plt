plot [][-0.05:]"C:/Users/hyoshida/Desktop/timedep/time_density_211024.dat" using 1:($2>0 ? $2: 1/0) with line
replot "C:/Users/hyoshida/Desktop/timedep/time_density_211024.dat" using 1:($3>0 ? $3: 1/0) with line
