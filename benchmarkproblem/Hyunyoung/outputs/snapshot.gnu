set term postscript eps enhanced color "Helvetica" 24

#set linestyle 1 lt 1 lw 50

set output "syncf3-nofail-init-8-snap.eps"
#set output "syncf1-nofail-final-8-snap.eps"
#set title "f1 for 500 Gen Pop size 50/thread"
set ylabel "Frequency"
#set logscale y
#set yrange [0.00001:100]
#set tics scale 100000
set logscale x 
set xlabel "Interval" 
set xrange [1.0e-5: 1.0e9]
#set xtics ("1e-7" 1e-7, "1e-6" 1e-6, "1e-5" 1e-5, "1e-4" 1e-4, "1e-3" 1e-3, "1e-2" 1e-2, "1e-1" 1e-1, "1e0" 1e0, "1e1" 1e1, "1e2" 1e2, "1e3" 1e3, "1e3" 1e3, "1e3" 1e3, "1e3" 1e3, "1e3" 1e3, "1e3" 1e3, "1e3" 1e3)
set xtics ("1e-5" 1e-5, "1e-4" 1e-4, "1e-3" 1e-3, "1e-2" 1e-2, "1e-1" 1e-1, "1e0" 1e0, "1e1" 1e1, "1e2" 1e2, "1e3" 1e3, "1e4" 1e4, "1e5" 1e5, "1e6" 1e6, "1e7" 1e7, "1e8" 1e8, "1e9" 1e9)

set xtics rotate by -90
#set format x "%1.3f"
set nokey

plot \
"snapshot/snapshot_init_8" using 1:2 title "f3 no failure" with imp lt 1 lw 15
#"snapshot/snapshot_final_8" using 1:2 title "f3 no failure" with imp lt 1 lw 15
