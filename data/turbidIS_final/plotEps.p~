reset
set term pdfcairo enhanced font ",12" size 5,3

#set style line 1 lt rgb "bdd7e7" lw 2 pt 1
#set style line 2 lt rgb "#6baed6" lw 2 pt 6
#set style line 3 lt rgb "#3182bd" lw 2 pt 2
#set style line 4 lt rgb "#08519c" lw 2 pt 9

set style line 1 linetype 1 linecolor rgb "#332288" linewidth 3 pointtype 4 pointsize 0.5
set style line 2 linetype 1 linecolor rgb "#88CCEE" linewidth 3 pointtype 5 pointsize 0.5
set style line 3 linetype 1 linecolor rgb "#999933" linewidth 3 pointtype 6 pointsize 0.5
set style line 4 linetype 1 linecolor rgb "#AA4499" linewidth 3 pointtype 7 pointsize 0.5

set output "absorbingEps_c.pdf"
set xtics nomirror
set ytics nomirror
set border 3 lw 2
set key top right spacing "1.3"

#set xrange[0:0.2]

#set xlabel "NA"
#set ylabel " "

plot "full_rho90.txt" u 1:(100*$2) w lines ls 1 t ' ', "full_rho95.txt" u 1:(100*$2) w lines ls 2 t ' ', "full_rho99.txt" u 1:(100*$2) w lines ls 3 t ' ' 
unset output
