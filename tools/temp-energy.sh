#!/bin/sh
./md --temp-energy 30.0 $1 > tmp/temp-energy.dat
gnuplot -persist tools/temp-energy.gnuplot
