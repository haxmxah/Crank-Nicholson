set term gif size 1000,600 animate  delay 10 loop 0 optimize
set output "P10-1819-fig5.gif"

# leeremos el numero de bloques de manera automatica
datafile ="check4.dat"
stats datafile
numblo=STATS_blocks

# pone el titulo
set title "T(x,t)"
set key top left

# fijamos los rangos del canvas, para que no se reajusten
set xrange[0:150]
set yrange[0:280]

#fijamos los títulos
set xlabel "x(cm)"
set ylabel "T(x,t) (grados)"

#bucle

do for[i=0:numblo-9:100]{
#escribe la etiqueta del tiempo
set label 2 sprintf('Time: %20.3f (sec)',i*0.002) at 10,120 right front font 'Verdana,12'

#dibuja
plot datafile index i u 2:3 w l lw 4 t"T(x,t)","" index 0 u 2:3 t"T_0(x)" w l lw 3

#borra la etiqueta
unset label 2
}
