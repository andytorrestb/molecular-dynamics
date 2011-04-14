#!/bin/sh
L=36.0 # 2 unit cells => density = 32/36 = 0.9
TEMP=0.7
./md --temp-energy $L $TEMP > tmp/temp-energy.dat
gnuplot -persist tools/temp-energy.gnuplot
