set xlabel 'Time'
set ylabel 'Position'
plot 'tmp/lj-pos.dat' using 1:2 title 'Atom 1' with lines,\
     'tmp/lj-pos.dat' using 1:3 title 'Atom 2' with lines
