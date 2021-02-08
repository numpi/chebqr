function [x, its] = cqr_eig(d, beta, u, v, np)
%CQR_EIG Compute the eigenvalues of a Hermitian-plus-rank-one matrix. 
%
% X = CQR_EIG(D, BETA, U, V) computes the eigenvalues of a
% Hermitian-plus-rank-one matrix A such that: 
% 
%  - The diagonal entries of A are stored in the vector D; 
%  - The subdiagonal entries of A are stored in BETA; 
%  - The matrix A - U*V', where U,V are column vectors, is Hermitian. 
%
% In particular, D,U,V need to be of length N, and BETA of length N-1.
%
% X = CQR_EIG(D, BETA, U, V, NP) sets the number of processors for the QR
% iteration to NP; if not specified, NP is automatically chosen as the
% number of cores available on the system (as detected by OpenMP).
%
% [X, BW] = CQR_EIG(...) also returns an upper bound for the normwise
% backward error on the polynomial for the computed roots. 
%
% [X, BW, NITS] = CQR_EIG(...) returns the number of sweeps performed
% by the QR iteration. 
%
% This function is conly implemented as a MEX file, which can be compiled
% running the compile() function in this folder. 

error('CQR_EIG is only available as a MEX file; please run the command ''compile''');

end

