function n = cqr_newton(p, x)
%CQR_NEWTON Evaluate the Newton correction at x

cp = chebfun(p(end:-1:1), 'coeffs');
cp1 = diff(cp);

n = cp(x) ./ cp1(x);

end

