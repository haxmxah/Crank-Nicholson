set terminal png
set output "P10-20-21-ﬁg1.png"
set term png size 1100, 600


set title  "Perfil de temperatures per l'estat estacionari" 
set title font "Times New Roman,20"
set xlabel font "Times New Roman,12"
set ylabel font "Times New Roman,12"
set xlabel "X[cm]" 
set ylabel "T[ºC]"
!set key outside
set xrange [0:50]
set yrange [0:300]


plot "check1.dat" i 0 w l lw 4 t"β=0.00004 s^{-1}", "check1.dat" i 1 w l lw 4 t"β=0.0003 s^{-1}", "check1.dat" i 2 w l lw 4 t"β=0.0025 s^{-1}"