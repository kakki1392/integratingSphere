reset
set term pdfcairo enhanced
set output "105coreNA0.22.pdf"

set ylabel("Collection efficiency [%]")
set xlabel("Sphere radius [mm]")

set title "Source: 105{/Symbol m}m core, NA = 0.22. Pick-up: 105{/Symbol m}m core, NA = 0.22"
set xtics nomirror
set ytics nomirror

set xrange[0.24:0.53]

plot "air_0.dat" u 1:(100*$2) lt 7 lc -1 ps 0.5 t '{/Symbol r} = 0.95', "air_1.dat" u 1:(100*$2) lt 7 lc 1 ps 0.5 t '{/Symbol r} = 0.96', "air_2.dat" u 1:(100*$2) lt 7 lc 2 ps 0.5 t '{/Symbol r} = 0.97', "air_3.dat" u 1:(100*$2) lt 7 lc 3 ps 0.5 t '{/Symbol r} = 0.98', "air_4.dat" u 1:(100*$2) lt 7 lc 4 ps 0.5 t '{/Symbol r} = 0.99'

unset output