#
# Makefile for ChebQR
#

include make.inc

all:
	make -C src/ all
	make -C examples/ all

clean:
	make -C examples/ clean
	make -C src/
