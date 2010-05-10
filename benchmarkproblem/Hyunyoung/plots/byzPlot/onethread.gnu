# Pool GA Plot data for a particular number of threads(e.g. 8 with varying percentage of failures 

set term postscript eps enhanced color "Helvetica" 24


set output "BAE-f1-8Threads.eps"
#set title "BAE-f2-900Gen-8Threads-Popsize=100/thread"
set ylabel "fitness"
set logscale y
#set yrange [1.0e-5: 0.1]
set yrange [10: 10000]
set xrange [0:1000]
#set xtics ("1" 1, "2" 2, "4" 3, "8" 4, "16" 5, "32" 6, "64" 7)
set xlabel "generation number" 
set pointsize 1.5
#set nokey
#set key 800, 80 
#"data8thread" using 1:2 title "1 thread" with linespoints,\
#"data8thread" using 1:3 title "2 threads" with linespoints,\
#"data8thread" using 1:4 title "4 threads" with linespoints,\
#"data8thread" using 1:5 title "8 threads" with linespoints,\
#"data8thread" using 1:6 title "16 threads" with linespoints,\
#"data8thread" using 1:7 title "32 threads" with linespoints,\
#"data8thread" using 1:8 title "64 threads" with linespoints

plot \
"f1data8thread" using 1:2 title "f=0.33" with lines lt 1 lw 4,\
"f1data8thread" using 1:3 title "f=0.40" with lines lt 2 lw 4,\
"f1data8thread" using 1:4 title "f=0.50" with lines lt 3 lw 4,\
"f1data8thread" using 1:5 title "f=0.60" with lines lt 4 lw 4,\
"f1data8thread" using 1:6 title "f=0.70" with lines lt 5 lw 4,\
"f1data8thread" using 1:7 title "f=0.80" with lines lt 7 lw 4
