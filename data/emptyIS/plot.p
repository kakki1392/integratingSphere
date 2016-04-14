reset
set term pdfcairo enhanced
set output "105core.pdf"

set ylabel("Collection efficiency")
set xlabel("Sphere radius [mm]")
plot "air_0.dat" u 1:2 lt 7 lc -1 ps 0.5 t '{/Symbol r} = 0.95', "air_1.dat" u 1:2 lt 7 lc 1 ps 0.5 t '{/Symbol r} = 0.96', "air_2.dat" u 1:2 lt 7 lc 2 ps 0.5 t '{/Symbol r} = 0.97', "air_3.dat" u 1:2 lt 7 lc 3 ps 0.5 t '{/Symbol r} = 0.98', "air_4.dat" u 1:2 lt 7 lc 4 ps 0.5 t '{/Symbol r} = 0.99'

unset output
