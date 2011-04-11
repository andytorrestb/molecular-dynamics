set xlabel 'Time'
set ylabel 'Energy'
set logscale y
plot 'tmp/temp-energy.dat' using 1:4 with lines title 'Total energy'
