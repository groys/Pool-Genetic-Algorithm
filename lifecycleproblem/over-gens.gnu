# Pool GA - best value seen over 30 generations. 

set term postscript eps "Helvetica" 24

set output "over-gens-pop1000.eps"
set title "Best value seen over 30 gens with pop size 1000 per thread"
set ylabel "fitness"
set yrange [0.79:0.96]
#set ylabel 3,0
#set logscale x 2
set xrange [0:31]
#set xtics ("1" 1, "2" 2, "4" 3, "8" 4, "16" 5, "32" 6, "64" 7)
set xlabel "generation number" 
#set nokey
set key 25, 0.94
plot \
"over_gens" using 1:2 title "1 thread" with linespoints,\
"over_gens" using 1:3 title "2 threads" with linespoints,\
"over_gens" using 1:4 title "4 threads" with linespoints,\
"over_gens" using 1:5 title "8 threads" with linespoints,\
"over_gens" using 1:6 title "16 threads" with linespoints,\
"over_gens" using 1:7 title "32 threads" with linespoints,\
"over_gens" using 1:8 title "64 threads" with linespoints
#"k15" using 1:2 title "k = 15" with linespoints,\
#"k17" using 1:2 title "k = 17" with linespoints
