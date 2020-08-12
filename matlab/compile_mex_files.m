mex cqr_qr_ss_aed.F90 ../src/fastfastqr.f90 -llapack -lblas
mex FFLAGS="$FFLAGS -fopenmp" LDFLAGS="$LDFLAGS -fopenmp" cqr_qr_ss_aed_par.F90 ../src/fastqr_par.f90 -llapack -lblas
