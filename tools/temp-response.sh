#!/bin/sh
./md --temp-response 0.5 0.1 0.9 0.1 > tmp/temp-response.dat
./md --temp-response 0.5 1 10 1 >> tmp/temp-response.dat
./md --temp-response 0.5 15 100 5 >> tmp/temp-response.dat
gnuplot -persist tools/temp-response.gnuplot
