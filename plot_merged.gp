# J-phi(J) plot for simulation and zero calculation
set output "C:/Users/hyoshida/Desktop/floquetic/phi_211210.png"
set terminal png
set multiplot

# d=0 with N=2 (left top)
set lmargin screen 0.1
set rmargin screen 0.49
set tmargin screen 0.9
set bmargin screen 0.51
plot [-0.3:0.3][-0.2:0.001]'C:/Users/hyoshida/Desktop/floquetic/phi_211210.dat' u ($1):($2) notitle
replot 'C:/Users/hyoshida/Desktop/floquetic/phi_211210.dat' u ($3):($4) notitle

# g=0 with N=2 (left bottom)
set lmargin screen 0.1
set rmargin screen 0.49
set tmargin screen 0.49
set bmargin screen 0.1
plot [-0.3:0.3][-0.2:0.001]'C:/Users/hyoshida/Desktop/floquetic/phi_211210.dat' u ($5):($6) notitle
replot 'C:/Users/hyoshida/Desktop/floquetic/phi_211210.dat' u ($7):($8) notitle

# d=0 with N=4 (right top)
set lmargin screen 0.51
set rmargin screen 0.9
set tmargin screen 0.9
set bmargin screen 0.51
plot [-0.3:0.3][-0.2:0.001]'C:/Users/hyoshida/Desktop/floquetic/phi_211210.dat' u ($9):($10) notitle
replot 'C:/Users/hyoshida/Desktop/floquetic/phi_211210.dat' u ($11):($12) notitle

# g=0 with N=4 (right bottom)
set lmargin screen 0.51
set rmargin screen 0.9
set tmargin screen 0.49
set bmargin screen 0.1
plot [-0.3:0.3][-0.2:0.001]'C:/Users/hyoshida/Desktop/floquetic/phi_211210.dat' u ($13):($14) notitle
replot 'C:/Users/hyoshida/Desktop/floquetic/phi_211210.dat' u ($15):($16) notitle

unset multiplot

set term png enhanced
replot
