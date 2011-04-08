#!/bin/sh
(cat tools/vel-dist.gnuplot; ./md --vel-dist 100000 10.0) | gnuplot > tmp/vel-dist.png
eog tmp/vel-dist.png
