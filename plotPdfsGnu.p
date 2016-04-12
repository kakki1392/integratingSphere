set term pdfcairo enhanced size 6in,4in font ",12"

set output "pathLengthDistributionsR1.0.pdf"
set title('Distribution of path lengths, empty sphere. R = 1.0 mm')
set xtics nomirror
set ytics nomirror
set xlabel('Total path length L [mm]')
set ylabel('Kernel estimate of f(L)')
plot "pdfs.txt" u 1:2 w lines t'{/Symbol r} = 0.95', "" u 1:3 w lines  t'{/Symbol r} = 0.96', "" u 1:4 w lines  t'{/Symbol r} = 0.97', "" u 1:5 w lines  t'{/Symbol r} = 0.98', "" u 1:6 w lines  t'{/Symbol r} = 0.99'
unset output
