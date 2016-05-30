set terminal pngcairo enhanced size 1440,900 crop font "Palatino,27"
set output "bilde.png"
plot sin(x), cos(x), x
unset output
