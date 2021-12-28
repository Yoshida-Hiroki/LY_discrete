set output "C:/Users/hyoshida/Desktop/floquetic/protocol_N.eps"
set term postscript eps enhanced
set multiplot

set lmargin screen 0.1
set rmargin screen 0.49
set bmargin screen 0.1
set tmargin screen 0.9

set xlabel ''
set ylabel ''
set xtics("0" 0 0,"" pi/2 1,"pi" pi 0,"" 3*pi/2 1, "2pi" 2*pi 0)
set ytics("0.0" 0 0,"" 0.2 1,"0.4" 0.4 0, "" 0.6 1, "0.8" 0.8 0)

# r(x) = 0.5
# phi_a(x) = pi*0.75+pi*0.2*cos(x)
# phi_b(x) = pi*0.5+0.4*pi*sin(x)
#
# a_L(x) = 0.5*(1+r(x))*sin(phi_a(x)/2)**2
# a_R(x) = 0.5*(1+r(x))*cos(phi_a(x)/2)**2
# b_L(x) = 0.5*(1-r(x))*sin(phi_b(x)/2)**2
# b_R(x) = 0.5*(1-r(x))*cos(phi_b(x)/2)**2
a_L(x) = 3/8.0*(1+cos(pi*0.25-pi*0.2*cos(x)))
a_R(x) = 3/8.0*(1-cos(pi*0.25-pi*0.2*cos(x)))
b_L(x) = 1/8.0*(1+sin(2*pi*0.2*sin(x)))
b_R(x) = 1/8.0*(1-sin(2*pi*0.2*sin(x)))

plot [0:2*pi][0:0.9] a_L(x) notitle dt 1 lc "black" lw 5
replot a_R(x) notitle dt 2 lc "black" lw 5
replot b_L(x) notitle dt 3 lc "black" lw 5
replot b_R(x) notitle dt 4 lc "black" lw 5

reset

set lmargin screen 0.51
set rmargin screen 0.9
set bmargin screen 0.1
set tmargin screen 0.9

set xlabel ''
set ylabel ''
set xtics("0" 0 0,"" pi/2 1,"pi" pi 0,"" 3*pi/2 1, "2pi" 2*pi 0)
set ytics("" 0 0,"" 0.2 1,"" 0.4 0, "" 0.6 1, "" 0.8 0)

# r(x) = 0.5+0.4*sin(x)
# phi_a(x) = pi*0.5+pi*0.2*cos(x)
#
# a_L(x) = 0.5*(1+r(x))*sin(phi_a(x)/2)**2
# a_R(x) = 0.5*(1+r(x))*cos(phi_a(x)/2)**2
# b_L(x) = 0.5*(1-r(x))*sin(phi_a(x)/2)**2
# b_R(x) = 0.5*(1-r(x))*cos(phi_a(x)/2)**2
a_L(x) = (0.75+0.2*sin(x))*(1+sin(pi*0.2*cos(x)))*0.5
a_R(x) = (0.75+0.2*sin(x))*(1-sin(pi*0.2*cos(x)))*0.5
b_L(x) = (0.25-0.2*sin(x))*(1+sin(pi*0.2*cos(x)))*0.5
b_R(x) = (0.25-0.2*sin(x))*(1-sin(pi*0.2*cos(x)))*0.5

plot [0:2*pi][0:0.9] a_L(x) notitle dt 1 lc "black" lw 5
replot a_R(x) notitle dt 2 lc "black" lw 5
replot b_L(x) notitle dt 3 lc "black" lw 5
replot b_R(x) notitle dt 4 lc "black" lw 5

reset
unset multiplot
