#!/bin/sh
(cat tools/fcc-pos.gnuplot; ./md --fcc-pos 3) | gnuplot -persist
