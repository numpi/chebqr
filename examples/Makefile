include ../make.inc

all: example1

example1: example1.f90 ../src/libchebqr.so
	$(FC) -o $@ $(FFLAGS) $< $(LIBS) -L.. -lchebqr '-Wl,-rpath=$$ORIGIN/../src'

clean:
	rm -f example1
