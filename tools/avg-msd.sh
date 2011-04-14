#!/bin/sh
DENSITY=0.9
./md --avg-msd $DENSITY 0.2 1 0.05 > tmp/avg-msd.dat
gnuplot -persist tools/avg-msd.gnuplot
