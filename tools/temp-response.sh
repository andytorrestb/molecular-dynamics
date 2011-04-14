#!/bin/sh
./md --temp-response 0.9 0.1 1 0.025 > tmp/temp-response.dat
gnuplot -persist tools/temp-response.gnuplot
