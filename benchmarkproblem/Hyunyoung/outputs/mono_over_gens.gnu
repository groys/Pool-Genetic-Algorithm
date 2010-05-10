set term postscript eps enhanced color "Helvetica" 24


set output "monoplot.eps"
#set title "f1 for 500 Gen Pop size 50/thread"
set ylabel "fitness"
#set logscale y
#set yrange [1.0e-6: 0.1]
set yrange [-4000:-1500]
#set ylabel 3,0
#set logscale x 2
set xrange [0:900]
#set xtics ("1" 1, "2" 2, "4" 3, "8" 4, "16" 5, "32" 6, "64" 7)
set xlabel "generation number" 
set pointsize 1.5
#set nokey
#set key 800, 80 

plot \
"mono_over_gens" using 1:3 title "2 threads" with lines lt 1 lw 4,\
"mono_over_gens" using 1:4 title "4 threads" with lines lt 2 lw 4,\
"mono_over_gens" using 1:5 title "8 threads" with lines lt 3 lw 4,\
"mono_over_gens" using 1:6 title "16 threads" with lines lt 4 lw 4,\
"mono_over_gens" using 1:7 title "32 threads" with lines lt 7 lw 4
