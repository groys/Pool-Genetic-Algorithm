set term postscript eps enhanced color "Helvetica" 24

#set linestyle 1 lt 1 lw 50

#set output "f1-nofail-init-16-snap.eps"
set output "f1-nofail-final-4-snap.eps"
#set title "f1 for 500 Gen Pop size 50/thread"
set ylabel "Frequency"
#set logscale y
#set yrange [1.0e-6: 0.1]
#set yrange [0.00001:100]
#set xtics ("1" 1, "2" 2, "4" 3, "4" 4, "8" 5, "32" 6, "64" 7)
set xtics rotate by -90
#set tics scale 100000
set logscale x 
set xlabel "Interval" 
#set format x "%1.3f"
set nokey

plot \
"f1-nofail-snapshot/snapshot_final_4" using 1:2:xticlabels(1) title "f1 nofail failure" with imp lt 1 lw 30
