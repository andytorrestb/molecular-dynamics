#!/bin/sh
DENSITY=0.5
./md --avg-msd $DENSITY 1 100 5 > tmp/avg-msd.dat
gnuplot -persist tools/avg-msd.gnuplot
