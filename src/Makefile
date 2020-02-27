#
# Makefile for ChebQR
#

include ../make.inc

LIBCHEBQR_SOURCES = \
	chebqr.f90 \
	fastfastqr.f90 \
	fastfastqr_double_noaed.f90 

all: libchebqr.so

libchebqr.so: $(LIBCHEBQR_SOURCES) ../chebqr.h utils.o
	$(FC) -shared -o libchebqr.so $(FFLAGS) $(LIBCHEBQR_SOURCES) utils.o $(LIBS)

utils.o: utils.c
	$(CC) -o $@ -c $<

clean:
	rm -f libchebqr.so utils.o
