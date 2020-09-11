mex FFLAGS="$FFLAGS -fopenmp" LDFLAGS="$LDFLAGS -fopenmp" cqr_eig.F90 ../src/single_shift.f90 -lmwlapack -lmwblas
