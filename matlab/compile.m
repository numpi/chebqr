mex FFLAGS="$FFLAGS -O0 -g -fopenmp" LDFLAGS="$LDFLAGS -O0 -g -fopenmp" cqr_eig.F90 ../src/single_shift.f90 ../src/double_shift.f90 -lmwlapack -lmwblas
