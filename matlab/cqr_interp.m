function p = cqr_interp(f, tol)
%CQR_INTERP Interpolate f with a Chebyshev polynomial over [-1, 1]

% At the moment we rely on Chebfun, if available
pp = chebfun(f);
p  = chebcoeffs(pp);
p  = p(end:-1:1);

end

