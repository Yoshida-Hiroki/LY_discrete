set output "C:/Users/hyoshida/Desktop/figure/protocol.eps"
set term postscript eps enhanced
set multiplot

set lmargin screen 0.05
set rmargin screen 0.45
set bmargin screen 0.2
set tmargin screen 0.8

set xlabel ''
set ylabel ''
set xtics 0,1,6
set ytics 0,0.1,0.9

r(x) = pi*0.25*sin(x)
phi_a(x) = pi*0.5+pi*0.25*cos(x)
phi_b(x) = pi*0.5+pi*0.25*sin(x)

a_L(x) = 0.5*(1+r(x))*sin(phi_a(x)/2)**2
a_R(x) = 0.5*(1+r(x))*cos(phi_a(x)/2)**2
b_L(x) = 0.5*(1-r(x))*sin(phi_b(x)/2)**2
b_R(x) = 0.5*(1-r(x))*cos(phi_b(x)/2)**2
# a_L(x) = 3/8.0*(1+cos(pi*0.25-pi*0.2*cos(x)))
# a_R(x) = 3/8.0*(1-cos(pi*0.25-pi*0.2*cos(x)))
# b_L(x) = 1/8.0*(1+sin(2*pi*0.2*sin(x)))
# b_R(x) = 1/8.0*(1-sin(2*pi*0.2*sin(x)))

plot [0:2*pi][0:0.9] a_L(x)+a_R(x) notitle dt 1 lc "black" lw 1,\
a_L(x) notitle dt 2 lc "black" lw 1,\
a_R(x) notitle dt 3 lc "black" lw 1

reset

set lmargin screen 0.55
set rmargin screen 0.95
set bmargin screen 0.2
set tmargin screen 0.8

set xlabel ''
set ylabel ''
set xtics 0,1,6
set ytics 0,0.1,0.9

plot [0:2*pi][0:0.9] b_L(x)+b_R(x) notitle dt 1 lc "black" lw 1,\
b_L(x) notitle dt 2 lc "black" lw 1,\
b_R(x) notitle dt 3 lc "black" lw 1

reset
unset multiplot
