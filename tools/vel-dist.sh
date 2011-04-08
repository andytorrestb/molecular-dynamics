#!/bin/sh
(cat tools/vel-dist.gnuplot; ./md --vel-dist 100000 10.0) | gnuplot -persist
