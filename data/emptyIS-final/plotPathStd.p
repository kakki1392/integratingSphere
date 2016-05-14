reset

set key top center
set xrange [0.3:1]
plot "rho95-combination3.txt" u 1:8 t "95", "rho96-combination3.txt" u 1:8 t "96", "rho97-combination3.txt" u 1:8 t "97", "rho98-combination3.txt" u 1:8 t "98"
