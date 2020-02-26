#
# Makefile for ChebQR
#

include ../make.inc

all: libchebqr.so

libchebqr.so: $(LIBCHEBQR_SOURCES) ../chebqr.h utils.o
	$(FC) -shared -o libchebqr.so $(FFLAGS) $(LIBCHEBQR_SOURCES) utils.o $(LIBS)

utils.o: utils.c
	$(CC) -o $@ -c $<

clean:
	rm -f libchebqr.so utils.o
