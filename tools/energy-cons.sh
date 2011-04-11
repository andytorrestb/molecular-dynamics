#!/bin/sh
./md --lj-pos 5000 > tmp/energy-cons.dat
gnuplot -persist tools/energy-cons.gnuplot
