function [be_fastqr, be_eig, be_eignb, deg, p] = test_be_f(f, mp_required)
%TEST_BE_F 

g=chebfun(f);
p=chebcoeffs(g);
p=p(end:-1:1)';

[d,beta,u,v] = cqr_colleague(p);
n=length(d) + 1;

[e] = cqr_eig(d,beta,u,v, 1);

% Dense QR
H = diag(ones(1, n-2), 1) + diag(ones(1,n-2), -1); H = H / 2;
H(end-1,end) = 1 / sqrt(2);
H(end,end-1) = 1 / sqrt(2);
e2 = eig(H + u*v);

% Dense QR without balancing
e2nb = eig(H + u*v, 'nobalance');

digts = 16;
if mp_required
    digts = 32; % 32 digits are usually enough to give accurate results.
end

be_fastqr = backward_error_pol(p(end:-1:1), e, digts);
be_eig    = backward_error_pol(p(end:-1:1), e2, digts);
be_eignb  = backward_error_pol(p(end:-1:1), e2nb, digts);
deg = length(p) - 1;

end

