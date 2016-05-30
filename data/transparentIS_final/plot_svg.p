set terminal svg size 500,300 dynamic enhanced fname 'arial'  fsize 14
set output "bilde.svg"

set xtics nomirror
set ytics nomirror


plot "105x105_rho90.txt" u 1:2 w linespoints pt 6 ps 0.5 t ' ', "105x105_rho93.txt" u 1:2 w linespoints pt 6 ps 0.3 t ' ', "105x105_rho96.txt" u 1:2 w linespoints pt 2 ps 0.8 t ' '
unset output
