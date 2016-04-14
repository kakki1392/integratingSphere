reset
set term pdfcairo enhanced
set output "200coreNA0.5.pdf"

set ylabel("Collection efficiency [%]")
set xlabel("Sphere radius [mm]")

set title "Source: 105{/Symbol m}m core, NA = 0.22. Pick-up: 200{/Symbol m}m core, NA = 0.5"
set xtics nomirror
set ytics nomirror

set xrange[0.24:0.53]

plot "air_10.dat" u 1:(100*$2) lt 7 lc -1 ps 0.5 t '{/Symbol r} = 0.95', "air_11.dat" u 1:(100*$2) lt 7 lc 1 ps 0.5 t '{/Symbol r} = 0.96', "air_12.dat" u 1:(100*$2) lt 7 lc 2 ps 0.5 t '{/Symbol r} = 0.97', "air_13.dat" u 1:(100*$2) lt 7 lc 3 ps 0.5 t '{/Symbol r} = 0.98', "air_14.dat" u 1:(100*$2) lt 7 lc 4 ps 0.5 t '{/Symbol r} = 0.99'

unset output
