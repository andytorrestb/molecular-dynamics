all:
	gcc -o md -Wall -lm src/*c

clean:
	rm -f md

.PHONY: all clean
