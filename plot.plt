plot "C:/Users/hyoshida/Desktop/den.dat" using 1:($2>0 ? $2: 1/0) with line
replot "C:/Users/hyoshida/Desktop/den.dat" using 1:($3>0 ? $3: 1/0) with line
