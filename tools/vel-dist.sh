#!/bin/sh
(cat tools/vel-dist.gnuplot; ./md --vel-dist 100000) | gnuplot
