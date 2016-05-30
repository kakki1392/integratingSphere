set term pdfcairo enhanced font ",12" size 5,3

#set style line 1 lt rgb "bdd7e7" lw 2 pt 1
#set style line 2 lt rgb "#6baed6" lw 2 pt 6
#set style line 3 lt rgb "#3182bd" lw 2 pt 2
#set style line 4 lt rgb "#08519c" lw 2 pt 9

set style line 1 linetype 1 linecolor rgb "#332288" linewidth 4 pointtype 1
set style line 2 linetype 1 linecolor rgb "#88CCEE" linewidth 4 pointtype 2
set style line 3 linetype 1 linecolor rgb "#999933" linewidth 4 pointtype 3
set style line 4 linetype 1 linecolor rgb "#AA4499" linewidth 4 pointtype 4

set output "test.pdf"
set xrange [0:3.5]
set xtics nomirror
set ytics nomirror
set border 3
plot sin(x+0.1) w lines ls 1, sin(x+0.2) w lines ls 2, sin(x+0.3) w lines ls 3, sin(x+0.4) w lines ls 4
unset output
