# Pool GA Plot data for a particular number of threads(e.g. 8 with varying percentage of failures 

set term postscript eps enhanced color "Helvetica" 24


set output "BAE-f3-8Threads-50perthread.eps"
#set title "BAE-f3-900Gen-8Threads-Popsize=100/thread"
set ylabel "fitness"
set logscale y
#set yrange [1.0e-5: 0.1]
#set yrange [10: 10000]
set xrange [0:900]
#set xtics ("1" 1, "2" 2, "4" 3, "8" 4, "16" 5, "32" 6, "64" 7)
set xlabel "generation number" 
set pointsize 1.5
#set nokey
#set key 800, 80 

plot \
"f3data8thread" using 1:2 title "f=0.33" with lines lt 1 lw 4,\
"f3data8thread" using 1:3 title "f=0.40" with lines lt 2 lw 4,\
"f3data8thread" using 1:4 title "f=0.50" with lines lt 3 lw 4,\
"f3data8thread" using 1:5 title "f=0.60" with lines lt 4 lw 4,\
"f3data8thread" using 1:6 title "f=0.70" with lines lt 5 lw 4,\
"f3data8thread" using 1:7 title "f=0.80" with lines lt 7 lw 4
