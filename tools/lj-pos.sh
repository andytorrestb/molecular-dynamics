#!/bin/sh
./md --lj-pos 5000 > tmp/lj-pos.dat
gnuplot -persist tools/lj-pos.gnuplot
