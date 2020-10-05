function [d, beta, u, v] = cqr_colleague(p)
%CQR_COLLEAGUE Colleague linearization for a polynomial.

% Decrease the degree, if necessary
ll = 0;
while p(ll+1) == 0
    ll = ll + 1;
end
p = p(ll+1:end);

n = length(p) - 1;
p = p(2:end) ./ p(1);

d = zeros(n, 1);
beta = ones(n-1, 1) / 2; beta(end) = 1 / sqrt(2);
v = -p; v(end) = v(end) * sqrt(2); v = v / 2;
u = eye(n, 1);
d(1) = v(1);

end

