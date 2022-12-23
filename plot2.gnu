set terminal png
set output "P10-20-21-ﬁg2.png"
set term png size 1100, 600


set title  "Evolució de les temperatures per x = 6, 42 i 126 cm" 
set title font "Times New Roman,20"
set xlabel font "Times New Roman,12"
set ylabel font "Times New Roman,12"
set xlabel "X[cm]" 
set ylabel "T[ºC]"
!set key outside

set key top left
plot "check2.dat" u 1:3 w l lw 3 t"x_p = 6 cm", "check2.dat" u 1:4 w l lw 3t"x_p = 42 cm","check2.dat" u 1:5 w l lw 3 t"x_p = 126 cm"