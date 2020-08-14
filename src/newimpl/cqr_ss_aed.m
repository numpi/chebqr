function [d, b, u, v, imin, imax, s, ns] = cqr_ss_aed(d, b, u, v, imin, imax)
%CQR_SS_AED 

k = 6;
s = zeros(k, 1);

% Compute the schur form of the trailing k x k block
Q = eye(k, 1); Q(1) = b(imax - k);

[d(imax-k+1:imax), b(imax-k+1:imax-1), u(imax-k+1:imax), v(imax-k+1:imax), ~, Q] = ...
    cqr_ss_qr1(d(imax-k+1:imax), b(imax-k+1:imax-1), u(imax-k+1:imax), v(imax-k+1:imax), false, true, Q);

% Count the number of eigenvalues to deflate
di = abs(d(imax));
for j = k : -1 : 1
    dim1 = di;
    di   = abs(d(imax-j));
    if abs(Q(j)) > ( dim1 + di ) * eps
        break;
    end
end

% Save the shifts for later use
ns = j;
s(1:ns) = d(imax - k + 1 : imax - k + j);

% Adjust imax
imax = imax - k + j;

% Restore the upper Hessenberg form of the rest of the matrix
[d, b, u, v] = cqr_ss_aed_restore_hess(d, b, u, v, Q, imax-j+1, imax);


end

