FC=gfortran
SOURCE_FILES=../src/single_shift.f90 ../src/double_shift.f90 ../src/utils.f90

all: simple_poly test_speed simple_poly_real

simple_poly: simple_poly.f90 ../src/utils.f90 ../src/single_shift.f90
	$(FC) -o simple_poly simple_poly.f90 ${SOURCE_FILES} -lblas -llapack

simple_poly_real: simple_poly_real.f90 ../src/utils.f90 ../src/double_shift.f90
	$(FC) -Wall -fbounds-check -g -O0 -o simple_poly_real simple_poly_real.f90 ${SOURCE_FILES}  -lblas -llapack

test_speed: test_speed.f90 ../src/utils.f90 ../src/single_shift.f90
	$(FC) -o test_speed test_speed.f90 ${SOURCE_FILES}  -lblas -llapack

clean:
	rm -f simple_poly test_speed simple_poly_real
