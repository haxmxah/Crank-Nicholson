set terminal png
set output "P10-20-21-ﬁg3.png"
set term png size 1100, 600


set title  "Evolució de les temperatures mitjanes per α_{au} i α_{fe} i β=0.00015 s^{-1}"  
set title font "Times New Roman,20"
set xlabel font "Times New Roman,12"
set ylabel font "Times New Roman,12"
set xlabel "X[cm]" 
set ylabel "T[ºC]"
!set key outside

set key top left
plot "check3.dat" i 0 u 1:2 w l lw 5 t"α_{fe}","check3.dat" i 1 u 1:2 w l lw 5 t"α_{au}"