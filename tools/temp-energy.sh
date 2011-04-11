#!/bin/sh
L=40.0
TEMP=10.0
./md --temp-energy $L $TEMP > tmp/temp-energy.dat
gnuplot -persist tools/temp-energy.gnuplot
