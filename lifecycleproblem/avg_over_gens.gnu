# Pool GA - avg best value seen in each of 60 generations. 

set term postscript eps "Helvetica" 24

set output "pop50avg-each-gen.eps"
set title "Avg best value seen in each of 60 gens with pop size 50 per thread"
set ylabel "fitness"
set yrange [0.3:0.70]
#set ylabel 3,0
#set logscale x 2
set xrange [0:61]
#set xtics ("1" 1, "2" 2, "4" 3, "8" 4, "16" 5, "32" 6, "64" 7)
set xlabel "generation number" 
#set nokey
set key 55, 0.50
plot \
"avg_over_gens" using 1:2 title "1 thread" with linespoints,\
"avg_over_gens" using 1:3 title "2 threads" with linespoints,\
"avg_over_gens" using 1:4 title "4 threads" with linespoints,\
"avg_over_gens" using 1:5 title "8 threads" with linespoints,\
"avg_over_gens" using 1:6 title "16 threads" with linespoints,\
"avg_over_gens" using 1:7 title "32 threads" with linespoints,\
"avg_over_gens" using 1:8 title "64 threads" with linespoints
