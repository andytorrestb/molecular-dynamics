binwidth = 0.1
bin(x,width) = width*floor(x/width) + binwidth/2

set xlabel 'Velocity'
set xrange [0:*]
set ylabel 'Frequency'
set yrange [0:*]
set boxwidth binwidth
plot '-' using (bin($1,binwidth)):(1.0) smooth freq with boxes
