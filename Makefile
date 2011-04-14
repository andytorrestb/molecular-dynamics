all:
	gcc -o md -Wall -lm src/*c

report:
	cd doc; pdflatex report.tex

clean:
	rm -f md
	rm -f doc/report.aux doc/report.log doc/report.pdf
	rm -f report.aux report.log

.PHONY: all clean
