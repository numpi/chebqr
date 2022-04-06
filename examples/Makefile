FC=gfortran

all: simple_poly

simple_poly: simple_poly.f90 ../src/utils.f90 ../src/single_shift.f90
	$(FC) -o simple_poly simple_poly.f90 ../src/utils.f90 ../src/single_shift.f90 -lblas -llapack

clean:
	rm -f simple_poly
