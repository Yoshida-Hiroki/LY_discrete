set xtics -0.2,0.2
set ytics -0.2,0.1

set xtics ("-0.2" -0.2 0,"" -0.15 1,"-0.1" -0.1 0,"" -0.05 1,"0" 0 0,"" 0.05 1, "0.1" 0.1 0, "" 0.15 1,"0.2" 0.2 0)
set ytics ("-0.3" -0.3 0,"" -0.25 1, "-0.2" -0.2 0,"" -0.15 1,"-0.1" -0.1 0,"" -0.05 1, "0" 0 0,"" 0.05 1)
set xlabel ""
set ylabel ""
set xrange [-0.2:0.2]
set yrange [-0.3:0.05]
set arrow from -0.2,0.0 to 0.2,0 nohead lw 0.5
set arrow from 0,-0.3 to 0,0.05 nohead lw 0.5

plot "C:/Users/hyoshida/Desktop/floquetic/phi_211225_N1_1.dat" u 2:3 notitle w l lc "black" dt 1 lw 3
replot "C:/Users/hyoshida/Desktop/floquetic/phi_211225_N1_2.dat" u 2:3 notitle w l lc "black" dt 1 lw 3
replot "C:/Users/hyoshida/Desktop/floquetic/sim_211225_N1_1.dat" u 1:2 every 10 notitle lc "black" pt 2 ps 2
replot "C:/Users/hyoshida/Desktop/floquetic/sim_211225_N1_2.dat" u 1:2 every 10 notitle lc "black" pt 4 ps 2

# set output "C:/Users/hyoshida/Desktop/floquetic/phi_211225_N1.png"
set term postscript eps enhanced
set output "C:/Users/hyoshida/Desktop/floquetic/phi_211225_N1.eps"

# set terminal png
# set term png enhanced
replot
reset
