# Pool GA - (monotonic) best value seen over 60 generations. 

#set term postscript eps "Helvetica" 24
set term postscript eps enhanced color "Helvetica" 24

set output "technophile-constant-pool.eps"
#set title "Technophile 100 gens,pop size 50/thread"
set ylabel "fitness"
#set autoscale
#set yrange [0.50:0.70]
#set ylabel 3,0
#set logscale x 2
#jset xrange [0:101]
#set xtics ("1" 1, "2" 2, "4" 3, "8" 4, "16" 5, "32" 6, "64" 7)
set xlabel "generation number" 
#set nokey
set key 75, 0.56
#set key 75, 0.81
#"mono_over_gens" using 1:2 title "1 thread" with linespoints,\
#"mono_over_gens" using 1:3 title "2 threads" with linespoints,\
#"mono_over_gens" using 1:4 title "4 threads" with linespoints,\
#"mono_over_gens" using 1:5 title "8 threads" with linespoints,\
#"mono_over_gens" using 1:6 title "16 threads" with linespoints,\
#"mono_over_gens" using 1:7 title "32 threads" with linespoints,\
#"mono_over_gens" using 1:8 title "64 threads" with linespoints
#"mono_over_gens" using 1:2 title "SGA" with lines lw 4,\

plot \
"mono_over_gens" using 1:2 title "SGA" with lines lw 4,\
"mono_over_gens" using 1:3 title "PGA 2 thread" with lines lw 4,\
"mono_over_gens" using 1:4 title "PGA 4 threads" with lines lw 4,\
"mono_over_gens" using 1:5 title "PGA 8 threads" with lines lw 4,\
"mono_over_gens" using 1:7 title "PGA 32 threads" with lines lw 4
