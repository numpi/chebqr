include ../make.inc

all: example1 example2 example3

example%: example%.f90 ../src/libchebqr.so
	$(FC) -o $@ $(FFLAGS) $< $(LIBS) -L../src/ -lchebqr '-Wl,-rpath=$$ORIGIN/../src'

clean:
	rm -f example1 example2
