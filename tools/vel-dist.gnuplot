set terminal postscript

binwidth = 0.2
bin(x,width) = width*floor(x/width) + binwidth/2

set xlabel 'Velocity'
set ylabel 'Frequency'
set boxwidth binwidth
plot '-' using (bin($1,binwidth)):(1.0) smooth freq with boxes
