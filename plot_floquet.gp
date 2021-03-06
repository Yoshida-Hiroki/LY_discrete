set term png size 800,600
set output "C:/Users/hyoshida/Desktop/floquet.png"
# set output "C:/Users/hyoshida/Desktop/floquet.eps"
# set term postscript enhanced eps size 8.1,5.0
set datafile separator ","
set multiplot

set lmargin screen 0.02
set rmargin screen 0.98
set bmargin screen 0.1
set tmargin screen 0.9

set xlabel ''
set ylabel ''
set xtics ("0.0" 0 0,"" 0.05, "0.1" 0.1 0, "" 0.15, "0.2" 0.2 0)
set ytics ("-2" -2 0,"" -1.5,"-1" -1 0,"" -0.5, "0" 0 0)

plot [0:0.25][-2:0.25]'C:/Users/hyoshida/Desktop/210903_zeros.csv' u 1:2 notitle w l lc "black",\
'C:/Users/hyoshida/Desktop/210903_zeros.csv' u 1:3 notitle w l lc "black",\
'C:/Users/hyoshida/Desktop/210903_zeros.csv' u 1:4 notitle pt 4 lc "red" pointsize 4,\
'C:/Users/hyoshida/Desktop/210903_zeros.csv' u 1:5 notitle pt 4 lc "red" pointsize 4,\
'C:/Users/hyoshida/Desktop/210903_zeros.csv' u 1:7 notitle pt 2 lc "blue" pointsize 4,\
'C:/Users/hyoshida/Desktop/210903_zeros.csv' u 1:8 notitle pt 2 lc "blue" pointsize 4,\
'C:/Users/hyoshida/Desktop/210903_zeros.csv' u 1:9 notitle pt 2 lc "blue" pointsize 4

reset
set datafile separator ","

set lmargin screen 0.65
set rmargin screen 0.92
set bmargin screen 0.52
set tmargin screen 0.78

set xlabel ''
set ylabel ''
set xtics ("0.0" 0 0,"" 0.05, "0.1" 0.1 0, "" 0.15, "0.2" 0.2 0)
set ytics ("0" 0 0,"" 20000, "40000" 40000 0, "" 60000, "80000" 80000 0)

plot [0:0.25][0:90000]'C:/Users/hyoshida/Desktop/210903_zeros.csv' u 1:6 every 2 notitle pt 2 lc "blue" pointsize 4

reset

# set datafile separator ","
#
# set lmargin screen 0.57
# set rmargin screen 0.95
# set bmargin screen 0.2
# set tmargin screen 0.8
#
# set xlabel ''
# set ylabel ''
# set xtics ("0.0" 0 0,"" 0.05, "0.1" 0.1 0, "" 0.15, "0.2" 0.2 0)
# set ytics ("-0.05" -0.05 0,"" -0.045,"-0.04" -0.04 0,"" 0.-0.035, "-0.03" -0.03 0)
#
# plot [0:0.25][-0.05:-0.03] 'C:/Users/hyoshida/Desktop/floquetic/210903_currents.csv' u 1:2 notitle w l lw 1 lc "black",\
# 'C:/Users/hyoshida/Desktop/floquetic/210903_currents.csv' u 1:3 notitle w l dt 1 lc "black",\
# 'C:/Users/hyoshida/Desktop/floquetic/210903_currents.csv' u 1:4 notitle w l dt 2 lc "black",\
# 'C:/Users/hyoshida/Desktop/floquetic/210903_currents.csv' u 1:5 notitle w l dt 3 lc "black"
#
# reset
unset multiplot
