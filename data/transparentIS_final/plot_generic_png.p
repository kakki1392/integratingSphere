set term pngcairo size 1920,1080 enhanced font 'Verdana,12'
set output "test.png"
plot sin(x)
unset output
