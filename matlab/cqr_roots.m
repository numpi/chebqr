function x = cqr_roots(f, parallel)
%CQR_ROOTS Find roots of a smooth function over [-1, 1].

if ~exist('parallel', 'var')
    parallel = false;
end

% Construct the polynomial interpolant that well-approximate f(x) over the
% interval [-1, 1]. 
p = cqr_interp(f);

% Build the structured representation of the colleague matrix
[d, beta, u, v] = cqr_colleague(p);

% Find all the roots of the polynomial using the structured QR iteration
if parallel
    x = cqr_qr_ss_aed_par(d, beta, u, v);
else
    x = cqr_qr_ss_aed(d, beta, u, v);
end

% Filter the roots (and make them real)
x = cqr_filter_roots(p, x);

% Sort them (now they are real, so this makes sense)
x = sort(x);

end

