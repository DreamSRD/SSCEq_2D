HEADER = -I/usr/local/include
LIBB = -L/usr/local/lib
LIBRA = -lfftw3_threads -lfftw3 -lm
SOURCES = SSCEq_2D_full.c
CFLAGS = -fopenmp

all:
	gcc -std=c99 $(CFLAGS) $(SOURCES) $(HEADER) $(LIBB) $(LIBRA) -o SSCEq_2D_full
clean:
	rm -rf *.o SSCEq_2D_full
